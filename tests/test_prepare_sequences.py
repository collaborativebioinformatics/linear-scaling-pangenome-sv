"""Tests for pipeline/prepare/prepare_sequences.py.

The faidx layer is stubbed because samtools is not guaranteed in CI. What
is under test here is the orchestration: interval arithmetic, PanSN
naming, strand handling, the length sanity gates, and the dual-filename
output that keeps build_baseline.sh working across branches.
"""
import csv
import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pipeline.prepare.faidx_utils as FU  # noqa: E402
import pipeline.prepare.prepare_sequences as PS  # noqa: E402


class _FakeFasta:
    """Stands in for samtools faidx over an in-memory genome."""

    def __init__(self):
        self.genome = {}

    def _key(self, p):
        return os.path.abspath(p)

    def add(self, path, contigs):
        self.genome[self._key(path)] = contigs
        with open(path, "w") as f:
            for c, s in contigs.items():
                f.write(">%s\n%s\n" % (c, s))

    def ensure_faidx(self, p):
        with open(p + ".fai", "w") as f:
            for c, s in self.genome.get(self._key(p), {}).items():
                f.write("%s\t%d\t0\t60\t61\n" % (c, len(s)))

    def extract(self, p, contig, s, e, strand="+"):
        seq = self.genome[self._key(p)][contig][s:e]
        return FU._revcomp(seq) if strand == "-" else seq

    def length(self, p, contig):
        return len(self.genome[self._key(p)][contig])


def _build_workspace(start=1000, end=2000, ref_len=5000, hap_len=4000,
                     src_start=500, src_end=1480, seed=7):
    """Create a temp project tree and return (dir, fake, ref_sequence)."""
    d = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(d)
    for sub in ("config", "work/reference", "work/downloads",
                "results/preparation"):
        os.makedirs(sub, exist_ok=True)

    fake = _FakeFasta()
    rnd = random.Random(seed)
    refseq = "".join(rnd.choice("ACGT") for _ in range(ref_len))
    fake.add("work/reference/GRCh38_chr21.fa", {"chr21": refseq})

    with open("config/pipeline.yaml", "w") as f:
        f.write("target:\n  chromosome: chr21\n  start: %d\n  end: %d\n"
                % (start, end))

    rows = []
    for sm, hp in [("HG00673", "1"), ("HG00673", "2"),
                   ("HG00733", "1"), ("HG00733", "2")]:
        an = "%s_%s_hprc_r2_v1.0.1" % (sm, "pat" if hp == "1" else "mat")
        ap = "work/downloads/%s.fa" % an
        fake.add(ap, {"h1tg000001l":
                      "".join(rnd.choice("ACGT") for _ in range(hap_len))})
        rows.append(dict(
            sample=sm, haplotype=hp,
            haplotype_label="paternal" if hp == "1" else "maternal",
            assembly_name=an, source_contig="h1tg000001l",
            source_start=str(src_start), source_end=str(src_end),
            strand="+" if hp == "1" else "-", status="mapped"))
    with open("results/preparation/sequence_mapping.tsv", "w",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    PS.ensure_faidx = fake.ensure_faidx
    PS.faidx_extract = fake.extract
    PS.get_contig_length = fake.length
    sys.argv = ["prepare_sequences.py"]
    return d, fake, refseq, cwd


def _records(path):
    """[(name, sequence)] from a FASTA."""
    out = []
    for blk in open(path).read().split(">")[1:]:
        lines = blk.split("\n")
        out.append((lines[0].strip(), "".join(lines[1:]).strip()))
    return out


def test_writes_both_stamped_and_stable_filenames():
    """REGRESSION: branches disagree on build_baseline.sh's default input.

    One expects chr21_multi.fa, another chr21_20000000_21000000_multi.fa.
    Writing both means `make baseline` works either way.
    """
    d, _, _, cwd = _build_workspace()
    try:
        PS.main()
        assert os.path.exists("results/preparation/chr21_1000_2000_multi.fa")
        assert os.path.exists("results/preparation/chr21_multi.fa")
    finally:
        os.chdir(cwd)


def test_paths_use_pansn_naming():
    d, _, _, cwd = _build_workspace()
    try:
        PS.main()
        names = [n for n, _ in
                 _records("results/preparation/chr21_1000_2000_multi.fa")]
        assert names == ["GRCh38#0#chr21", "HG00673#1#chr21",
                         "HG00673#2#chr21", "HG00733#1#chr21",
                         "HG00733#2#chr21"]
    finally:
        os.chdir(cwd)


def test_reference_interval_is_exact():
    """0-based half-open [start, end) — guards against an off-by-one."""
    d, _, refseq, cwd = _build_workspace(start=1000, end=2000)
    try:
        PS.main()
        recs = _records("results/preparation/chr21_1000_2000_multi.fa")
        assert recs[0][1] == refseq[1000:2000]
        assert len(recs[0][1]) == 1000
    finally:
        os.chdir(cwd)


def test_reverse_strand_is_reverse_complemented():
    d, fake, _, cwd = _build_workspace()
    try:
        PS.main()
        recs = _records("results/preparation/chr21_1000_2000_multi.fa")
        hap2 = dict(recs)["HG00673#2#chr21"]
        src = fake.genome[fake._key(
            "work/downloads/HG00673_mat_hprc_r2_v1.0.1.fa")]["h1tg000001l"]
        assert hap2 == FU._revcomp(src[500:1480])
    finally:
        os.chdir(cwd)


def test_reference_is_not_slurped_into_memory():
    """The reference must be pulled via faidx, like the haplotypes.

    The previous implementation read the whole reference FASTA into a
    Python string and sliced it, contradicting its own docstring and
    blowing up on a whole-genome reference.
    """
    d, fake, _, cwd = _build_workspace()
    calls = []
    orig = fake.extract

    def spy(p, contig, s, e, strand="+"):
        calls.append(os.path.basename(p))
        return orig(p, contig, s, e, strand)

    PS.faidx_extract = spy
    try:
        PS.main()
        assert "GRCh38_chr21.fa" in calls, \
            "reference was not extracted through faidx"
    finally:
        os.chdir(cwd)


def test_oversized_haplotype_is_rejected():
    """A haplotype far longer than the reference interval means bad mapping."""
    d, _, _, cwd = _build_workspace(start=1000, end=1100,
                                    src_start=0, src_end=3000)
    try:
        try:
            PS.main()
        except SystemExit as e:
            assert e.code != 0
        else:
            raise AssertionError("expected a non-zero exit on bad mapping")
    finally:
        os.chdir(cwd)


def test_path_count_follows_manifest_not_a_hardcoded_five():
    d, _, _, cwd = _build_workspace()
    try:
        rows = list(csv.DictReader(
            open("results/preparation/sequence_mapping.tsv"), delimiter="\t"))
        with open("results/preparation/sequence_mapping.tsv", "w",
                  newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            for r in rows[:2]:
                w.writerow(r)
        PS.main()
        recs = _records("results/preparation/chr21_1000_2000_multi.fa")
        assert len(recs) == 3  # 2 haplotypes + reference
    finally:
        os.chdir(cwd)
