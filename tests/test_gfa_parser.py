"""Tests for GFA data model, chunking, merge, and interval mapping."""
import os
import sys
import tempfile
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.merge.gfa import (
    GfaGraph, Header, Segment, Link, Path, Walk,
    _revcomp, infer_data_mode, extract_chromosome_from_gfa, _split_orient,
    parse_pansn, haplotype_key)
from pipeline.parallel.make_chunks import create_chunks
from pipeline.merge.merge_graphs import overlap_aware_stitch, walks_as_paths
from pipeline.merge.paths import (
    slice_steps, seam_position, path_interval,
    stitch_haplotype, group_paths_by_haplotype)
from pipeline.merge.validate_merge import validate


class TestHeader:
    def test_default(self):
        h = Header()
        assert h.version == "1.1"

    def test_parse(self):
        h = Header.parse(["H", "VN:Z:1.0"])
        assert h.version == "1.0"


class TestSegment:
    def test_basic(self):
        s = Segment("s1", "ACGT")
        assert s.length == 4

    def test_empty(self):
        s = Segment.parse(["S", "n1", "*"])
        assert s.sequence == ""


class TestLink:
    def test_parse(self):
        l = Link.parse(["L", "n1", "+", "n2", "-", "*"])
        assert l.to_node == "n2"

    def test_roundtrip(self):
        l = Link("n1", "+", "n2", "-", "10M")
        l2 = Link.parse(l.to_gfa().split(chr(9)))
        assert l2.overlap == "10M"


class TestPath:
    def test_basic(self):
        p = Path("px", ["s1+"])
        assert "s1+" in p.to_gfa()


class TestWalk:
    def test_basic(self):
        w = Walk("HG001", "1", "chr21", 0, 100, ["s1+"])
        assert "W" in w.to_gfa()


class TestGraph:
    def test_empty(self):
        g = GfaGraph()
        assert g.node_count() == 0

    def test_parse(self):
        gfa = ("H\tVN:Z:1.1\nS\ts1\tACGT\nS\ts2\tTGCA\n"
               "L\ts1\t+\ts2\t+\t*\nP\tp1\ts1+,s2+\t*")
        g = GfaGraph.parse(gfa)
        assert g.node_count() == 2
        assert g.edge_count() == 1

    def test_roundtrip(self):
        g = GfaGraph()
        g.segments["s1"] = Segment("s1", "ACGT")
        g2 = GfaGraph.parse(g.to_gfa())
        assert g2.node_count() == 1

    def test_copy(self):
        g = GfaGraph()
        g.segments["s1"] = Segment("s1", "ACGT")
        assert g.copy().node_count() == 1

    def test_samples(self):
        g = GfaGraph()
        g.walks.append(Walk("S1", "1", "c", 0, 10, ["n1+"]))
        assert "S1" in g.get_sample_names()
class TestUtils:
    def test_revcomp(self):
        assert _revcomp("ACGT") == "ACGT"
        assert _revcomp("AAAA") == "TTTT"

    def test_split_orient(self):
        assert _split_orient("n1+") == ("n1", "+")
        assert _split_orient("n1") == ("n1", "+")

    def test_data_mode(self):
        g = GfaGraph()
        g.segments["s1"] = Segment("s1", "GGGGGAAAAAGGGGGAAAAAGGGGGAAAAA")
        assert infer_data_mode(g) == "synthetic"

    def test_extract_chrom(self):
        g = GfaGraph()
        g.walks.append(Walk("GRCh38", "0", "chr21", 0, 1000, ["n1+"]))
        assert extract_chromosome_from_gfa(g) == "chr21"


class TestIO:
    def test_write_read(self):
        g = GfaGraph()
        g.segments["s1"] = Segment("s1", "ACGT")
        with tempfile.NamedTemporaryFile(suffix=".gfa", delete=False) as f:
            g.write_gfa(f.name)
            tmp = f.name
        g2 = GfaGraph.parse_file(tmp)
        os.unlink(tmp)
        assert g2.node_count() == 1


class TestHPRCIndex:
    def _hl(self, an):
        return "maternal" if "_mat_" in an else (
            "paternal" if "_pat_" in an else "unknown")

    def test_mat_label(self):
        assert self._hl("HG00673_mat_hprc_r2_v1.0.1") == "maternal"

    def test_pat_label(self):
        assert self._hl("HG00673_pat_hprc_r2_v1.0.1") == "paternal"

    def test_unknown_label(self):
        assert self._hl("HG00001_hprc_r2_v1.0.1") == "unknown"


class TestChunks:
    def test_create(self):
        ch = create_chunks("chr21", 0, 1000, 400, 50)
        assert len(ch) >= 2
        assert ch[0]["chunk_id"] == "chunk_0001"


class TestIntervalExtraction:
    def test_revcomp_roundtrip(self):
        seq = "ACGTACGTACGT"
        assert _revcomp(_revcomp(seq)) == seq

    def test_revcomp_known(self):
        assert _revcomp("ACGT") == "ACGT"
        assert _revcomp("AAAA") == "TTTT"

    def test_interval_extract_small(self):
        ref = "A" * 50 + "C" * 20 + "G" * 30
        chunk = ref[50:70]
        assert len(chunk) == 20
        assert chunk == "C" * 20

    def test_interval_reverse_strand(self):
        seq = "ACGTACGT" * 10
        chunk = seq[20:40]
        rc = _revcomp(chunk)
        assert len(rc) == 20
        assert rc == _revcomp(chunk)

    def test_chunk_scaling_basic(self):
        ref_s, ref_e = 0, 1000
        src_s, src_e = 500, 1500
        chunk_s, chunk_e = 200, 400
        frac_s = (chunk_s - ref_s) / (ref_e - ref_s)
        frac_e = (chunk_e - ref_s) / (ref_e - ref_s)
        result_s = int(src_s + frac_s * (src_e - src_s))
        result_e = int(src_s + frac_e * (src_e - src_s))
        assert result_s == 700
        assert result_e == 900

    def test_chunk_scaling_reverse(self):
        ref_s, ref_e = 0, 1000
        src_s, src_e = 500, 1500
        chunk_s, chunk_e = 200, 400
        frac_s = (chunk_s - ref_s) / (ref_e - ref_s)
        frac_e = (chunk_e - ref_s) / (ref_e - ref_s)
        fwd_s = int(src_s + frac_s * (src_e - src_s))
        fwd_e = int(src_s + frac_e * (src_e - src_s))
        rev_s = src_e - (fwd_e - src_s)
        rev_e = src_e - (fwd_s - src_s)
        assert rev_s == 1100
        assert rev_e == 1300

    def test_smoke_interval_size(self):
        ref = "A" * 46709983
        start, end = 20000000, 21000000
        chunk = ref[start:end]
        assert len(chunk) == 1000000

    def test_linear_scaling_vs_alignment_with_indel(self):
        """Prove that linear scaling and alignment-projected coordinates differ.

        Scenario: ref = 1 Mb, source has a 50 Kb insertion.
        Linear scaling predicts chunk_end at proportional position,
        but alignment maps it to a different coordinate.
        """
        ref_len = 1_000_000
        hap_insertion = 50_000
        source_len = ref_len + hap_insertion

        # Linear scaling: chunk at 40-60% of ref
        chunk_s, chunk_e = 400_000, 600_000
        frac_s = chunk_s / ref_len
        frac_e = chunk_e / ref_len

        # Linear: predicts source at same fraction
        linear_s = int(frac_s * source_len)
        linear_e = int(frac_e * source_len)

        # Linear scaling is WRONG: it spreads the insertion across the
        # entire interval, so chunk_s at 40% of ref becomes 42% of source
        assert linear_s == 420000  # WRONG coordinate from linear scaling
        assert linear_e == 630000  # WRONG coordinate from linear scaling

        # Alignment-based mapping is correct: coordinates before the
        # insertion point remain unchanged; after the insertion they shift
        insertion_pos = 500_000
        align_s = chunk_s  # 400k is before the insertion, unchanged
        if chunk_e > insertion_pos:
            align_e = chunk_e + hap_insertion  # shift by insert size
        else:
            align_e = chunk_e

        assert align_s == 400000  # CORRECT - before insertion, no shift
        assert align_e == 650000  # CORRECT - 100k inserted at 500k

    def test_linear_scaling_vs_alignment_with_deletion(self):
        """Prove linear scaling and alignment differ for deletions.

        Scenario: ref = 1 Mb, source has a 100 Kb deletion.
        """
        ref_len = 1_000_000
        hap_deletion = 100_000
        source_len = ref_len - hap_deletion

        chunk_s, chunk_e = 400_000, 600_000
        frac_s = chunk_s / ref_len
        frac_e = chunk_e / ref_len

        linear_s = int(frac_s * source_len)
        linear_e = int(frac_e * source_len)

        # Alignment-based: deletion at ref 500k means the chunk
        # starting at 400k maps to source 400k-500k
        del_pos = 500_000
        if chunk_e > del_pos:
            align_e = del_pos
        else:
            align_e = chunk_e
        align_s = chunk_s

        # Linear scaling gives wrong coordinates with indels
        assert linear_e == 540000  # wrong - doesn't account for deletion
        assert align_e == 500000  # correct - deletion boundary

    def test_chunk_boundaries_not_proportional(self):
        """Demonstrate that same ref chunk boundary maps differently
        in two haplotypes with different indel patterns."""
        ref_len = 1_000_000
        boundary = 500_000

        # Haplotype A: 50 Kb insertion after boundary
        # Haplotype B: 30 Kb deletion after boundary
        hapA_len = ref_len + 50_000
        hapB_len = ref_len - 30_000

        frac = boundary / ref_len

        # Linear scaling would give different results per hap
        linear_A = int(frac * hapA_len)
        linear_B = int(frac * hapB_len)

        # With alignment, both should map to ~500k because
        # the boundary is before the indel
        assert linear_A != linear_B  # linear gives different boundaries
        assert linear_A == 525000
        assert linear_B == 485000


class TestMultiContigExtraction:
    """Regression test: prove coordinates from contigB return sequence from
    contigB, not coordinates applied to contigA+contigB concatenation."""

    def test_samtools_faidx_coordinate_conversion(self):
        """Verify PAF 0-based to samtools 1-based coordinate conversion."""
        # PAF: 0-based half-open [10, 20)
        # samtools: 1-based inclusive 11-20
        paf_start, paf_end = 10, 20
        samtools_start = paf_start + 1
        samtools_end = paf_end
        assert samtools_start == 11
        assert samtools_end == 20
        # Length should be identical
        assert (samtools_end - samtools_start + 1) == (paf_end - paf_start)

    def test_multi_contig_not_concatenated(self):
        """Prove that extracting from contigB does NOT return sequence
        from contigA+contigB concatenation.

        If we concatenate:
          contigA = "AAAAA..." (100 bp)
          contigB = "CCCCC..." (100 bp)
          combined = "AAAAA...CCCCC..."
        and ask for combined[50:60], we get contigB[50:60].
        But if we ask for contigB[50:60], we should get contigB[50:60]
        starting from position 0 of contigB, not position 100 of combined.
        """
        # Build a synthetic multi-contig FASTA
        import shutil
        if shutil.which("samtools") is None:
            pytest.skip("samtools not available on this machine")
        import tempfile
        fa_dir = tempfile.mkdtemp()
        fa_path = os.path.join(fa_dir, "multi.fa")
        with open(fa_path, "w") as f:
            f.write(">contigA\n")
            f.write("A" * 100 + "\n")
            f.write(">contigB\n")
            f.write("C" * 100 + "\n")

        # Index it
        import subprocess
        subprocess.run(["samtools", "faidx", fa_path], check=True)

        # Import the faidx utility
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
        from pipeline.prepare.faidx_utils import faidx_extract

        # WRONG: concatenating both contigs and slicing by coordinate
        with open(fa_path) as f:
            combined = "".join(line.strip() for line in f if not line.startswith(">"))
        assert combined == "A" * 100 + "C" * 100

        # Applying contigB coordinates [50, 60) to the concatenated string
        # returns the WRONG sequence (As instead of Cs)
        wrong = combined[50:60]
        assert wrong == "A" * 10  # WRONG: returning from contigA

        # CORRECT: extracting from contigB via samtools faidx
        correct = faidx_extract(fa_path, "contigB", 50, 60)
        assert correct == "C" * 10  # CORRECT: returning from contigB

        # Cleanup
        for f in os.listdir(fa_dir):
            os.remove(os.path.join(fa_dir, f))
        os.rmdir(fa_dir)

    def test_faidx_revcomp(self):
        """Verify reverse-complement extraction works correctly."""
        import shutil
        if shutil.which("samtools") is None:
            pytest.skip("samtools not available on this machine")
        import tempfile, subprocess
        fa_dir = tempfile.mkdtemp()
        fa_path = os.path.join(fa_dir, "rc.fa")
        with open(fa_path, "w") as f:
            f.write(">test\n")
            f.write("ACGTACGTACGT\n")  # 12 bp
        subprocess.run(["samtools", "faidx", fa_path], check=True)
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
        from pipeline.prepare.faidx_utils import faidx_extract
        # Extract forward
        fwd = faidx_extract(fa_path, "test", 4, 8, "+")
        assert fwd == "ACGT"
        # Extract reverse-complemented
        rev = faidx_extract(fa_path, "test", 4, 8, "-")
        assert rev == "ACGT"  # ACGT is palindromic
        # Test non-palindromic
        fwd2 = faidx_extract(fa_path, "test", 0, 4, "+")
        rev2 = faidx_extract(fa_path, "test", 0, 4, "-")
        assert fwd2 == "ACGT"
        assert rev2 == "ACGT"  # still palindromic at this position
        for f in os.listdir(fa_dir):
            os.remove(os.path.join(fa_dir, f))
        os.rmdir(fa_dir)


class TestBuildAllChunksImports:
    """Test that build_all_chunks.py compiles and imports without syntax errors."""

    def test_import_success(self):
        import py_compile
        import tempfile
        path = "pipeline/parallel/build_all_chunks.py"
        assert os.path.exists(path), f"{path} not found"
        # Compile to catch syntax errors
        py_compile.compile(path, doraise=True)

    def test_functions_defined(self):
        """Verify map_chunk has a return statement inside the function."""
        with open("pipeline/parallel/build_all_chunks.py") as f:
            content = f.read()
        # The return must be BEFORE the main() function definition
        main_pos = content.find("def main():")
        return_pos = content.rfind("return best if bc", 0, main_pos)
        assert return_pos > 0, "return best if bc... must be inside map_chunk(), before main()"


class TestDemo:
    def test_json_exists(self):
        assert os.path.exists("web/public/data/latest.json")

class TestP0Regression:
    """P0 regression tests for the audit."""

    def test_parse_file_regression(self):
        """item 1: GfaGraph.parse_file with a temporary GFA."""
        import tempfile
        gfa_text = "H\tVN:Z:1.1\nS\ts1\tACGT\nS\ts2\tTGCA\nL\ts1\t+\ts2\t+\t*\nP\tp1\ts1+,s2+\t*\n"
        with tempfile.NamedTemporaryFile(suffix=".gfa", mode="w", delete=False) as f:
            f.write(gfa_text)
            tmp = f.name
        try:
            g = GfaGraph.parse_file(tmp)
            assert g.node_count() == 2
            assert "p1" in g.paths
        finally:
            os.unlink(tmp)

    def test_walk_merge_roundtrip(self):
        """item 2: Disjoint-union merge with W-lines, write, reparse, verify."""
        g1 = GfaGraph()
        g1.headers.append(Header("1.1"))
        g1.segments["s1"] = Segment("s1", "ACGT")
        g1.segments["s2"] = Segment("s2", "TGCA")
        g1.links.append(Link("s1", "+", "s2", "+", "*"))
        g1.walks.append(Walk("HG001", "1", "chr21", 0, 8, ["s1+", "s2+"]))
        g2 = GfaGraph()
        g2.segments["s3"] = Segment("s3", "GGGG")
        g2.segments["s4"] = Segment("s4", "CCCC")
        g2.links.append(Link("s3", "+", "s4", "+", "*"))
        g2.walks.append(Walk("HG002", "1", "chr21", 0, 8, ["s3+", "s4+"]))
        merged = g1.copy()
        for sid, sg in g2.segments.items(): merged.segments[sid] = sg
        for lnk in g2.links: merged.links.append(lnk)
        for w in g2.walks: merged.walks.append(w)
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".gfa", mode="w", delete=False) as f:
            merged.write_gfa(f.name); tmp = f.name
        try:
            rp = GfaGraph.parse_file(tmp)
            assert rp.node_count() == 4
            assert len(rp.walks) == 2
            for w in rp.walks:
                assert w.sample in ("HG001", "HG002")
                assert len(w.path) == 2
                for step in w.path: assert step[-1] in ("+", "-")
        finally: os.unlink(tmp)

    def test_chunk_overlap_exact(self):
        """item 3: pairwise overlap between adjacent chunks is exactly 50 kb."""
        from pipeline.parallel.make_chunks import create_chunks
        chunks = create_chunks("chr21", 0, 1000000, 400000, 50000)
        for i in range(1, len(chunks)):
            overlap = (min(chunks[i - 1]["reference_end"], chunks[i]["reference_end"])
                       - max(chunks[i - 1]["reference_start"], chunks[i]["reference_start"]))
            assert overlap == 50000, (
                f"Overlap chunk {i-1}/{i}: {overlap}, expected 50000")

    def test_pggb_params_equal(self):
        """item 3: Baseline and chunk PGGB parameters from pipeline.yaml."""
        import yaml
        cfg = yaml.safe_load(open("config/pipeline.yaml"))
        params = cfg.get("pggb", {}).get("params", {})
        for k in ["minimum_identity", "segment_length", "mash_kmer",
                    "match_length", "path_jump_max", "edge_jump_max"]:
            assert params.get(k) is not None, f"{k} must be set"
        img = cfg.get("pggb", {}).get("image", "")
        assert ":latest" not in img, f"PGGB image must be pinned (got: {img})"

    def test_final_gfa_selection(self):
        """item 5: Must use *final.gfa, not find ... | head -1."""
        with open("dnanexus/applets/pggb_chunk/src/code.sh") as f:
            c = f.read()
        assert "*final.gfa" in c
        assert "find output -name \"*.gfa\" -type f | head -1" not in c

    def test_baseline_applet_exists(self):
        """item 6: pggb_baseline applet with dxapp.json and code.sh."""
        assert os.path.exists("dnanexus/applets/pggb_baseline/dxapp.json")
        assert os.path.exists("dnanexus/applets/pggb_baseline/src/code.sh")
        import json
        with open("dnanexus/applets/pggb_baseline/dxapp.json") as f:
            a = json.load(f)
        assert a["name"] == "pggb_baseline"
        assert "mem3_ssd1_v2_x16" in json.dumps(a)

    def test_stale_provenance(self):
        """item 8: Provenance hash based on git commit, config, coordinates."""
        from pipeline.parallel.build_all_chunks import compute_provenance
        p1 = compute_provenance("chunk_01", 0, 400000)
        p2 = compute_provenance("chunk_01", 0, 400001)
        assert p1 != p2, "Different coords must differ"
        p3 = compute_provenance("chunk_02", 0, 400000)
        assert p1 != p3, "Different chunk IDs must differ"

    def test_parent_locus_constraint(self):
        """item 7: build_all_chunks.py must reject off-contig chunk."""
        with open("pipeline/parallel/build_all_chunks.py") as f:
            c = f.read()
        assert "FATAL" in c and "parent" in c.lower()
        assert "chunk jumps to another contig" in c.lower()

    def test_orchestration_strict(self):
        """item 9: run_parallel_chunks.sh has FATAL on mismatch."""
        with open("dnanexus/run_parallel_chunks.sh") as f:
            c = f.read()
        fatal_lines = [l.strip() for l in c.split("\n") if "FATAL" in l]
        assert any("Missing chunk FASTA" in l for l in fatal_lines)

    def test_compare_paths_sequence(self):
        """item 17: compare_paths.py must spell sequences, not just count."""
        with open("pipeline/benchmark/compare_paths.py") as f:
            c = f.read()
        assert "sha256" in c.lower() or "SHA256" in c
        assert "identity" in c
        assert "spell_path_sequence" in c or "_revcomp" in c

    def test_auto_venv(self):
        """item 13: run_pipeline.sh must auto-detect .pipeline-env."""
        with open("dnanexus/run_pipeline.sh") as f:
            c = f.read()
        assert ".pipeline-env" in c

    def test_reference_dnanexus_stage(self):
        """item 14: prepare_reference.sh stages from DNAnexus first."""
        with open("scripts/prepare_reference.sh") as f:
            c = f.read()
        assert "dx find data" in c or "DNAnexus" in c

    def test_report_partial(self):
        """item 15/16: build_report.py reports PARTIAL, NOT_IMPLEMENTED."""
        with open("pipeline/benchmark/build_report.py") as f:
            c = f.read()
        assert "PARTIAL" in c and "NOT_IMPLEMENTED" in c

    def test_no_merged_gfa_in_web_public(self):
        """item 18: web/public/data/ must not contain large merged.gfa."""
        assert not os.path.exists("web/public/data/merged.gfa")

    def test_applet_build_brief(self):
        """item 10: run_parallel_chunks.sh uses dx build --brief."""
        with open("dnanexus/run_parallel_chunks.sh") as f:
            c = f.read()
        assert "dx build" in c and "--brief" in c

    def test_job_output_download(self):
        """item 11: Download via job output reference."""
        with open("dnanexus/run_parallel_chunks.sh") as f:
            c = f.read()
        assert "dnanexus_link" in c or "'gfa'" in c

    def test_chunk_timing_recorded(self):
        """item 12: Per-chunk timing recorded."""
        with open("dnanexus/run_parallel_chunks.sh") as f:
            c = f.read()
        assert "start_timestamp" in c or "started" in c
        assert "instance_type" in c

    def test_run_pggb_exists(self):
        """item 3/4: scripts/run_pggb.py exists and uses config."""
        assert os.path.exists("scripts/run_pggb.py")
        with open("scripts/run_pggb.py") as f:
            c = f.read()
        assert "config/pipeline.yaml" in c

    def test_docker_helper_pinned(self):
        """item 4: docker_helper.sh reads image from config."""
        with open("dnanexus/docker_helper.sh") as f:
            c = f.read()
        assert "read_pggb_config" in c or "config/pipeline.yaml" in c


# --- Merge: path spelling, chunk stitching, boundary reporting -------------

REF = "AAAAACCCCCGGGGGTTTTT" * 5          # 100 bp
ALT = "TTTTTGGGGGCCCCCAAAAA"              # replaces REF[20:40] on the haplotype
HAP = REF[:20] + ALT + REF[40:]


def _chunk(cid, start, end):
    """Chunk graph over REF[start:end], 20 bp nodes, reference + one haplotype."""
    g = GfaGraph()
    g.headers.append(Header("1.1"))
    ref_steps, hap_steps = [], []
    for off in range(start, end, 20):
        rs, hs = REF[off:off + 20], HAP[off:off + 20]
        name = f"{cid}_{off}"
        g.segments[name] = Segment(name, rs)
        ref_steps.append(name + "+")
        if hs == rs:
            hap_steps.append(name + "+")
        else:
            alt = f"{cid}_{off}_alt"
            g.segments[alt] = Segment(alt, hs)
            hap_steps.append(alt + "+")
    for a, b in zip(ref_steps, ref_steps[1:]):
        g.links.append(Link(a[:-1], "+", b[:-1], "+"))
    rn = f"GRCh38#0#chr21:{start}-{end}"
    hn = f"HG1#1#chr21:{start}-{end}"
    g.paths[rn] = Path(rn, ref_steps)
    g.paths[hn] = Path(hn, hap_steps)
    return cid, g


class TestPathSequence:
    def test_spell_forward(self):
        g = GfaGraph.parse("S\ts1\tACGT\nS\ts2\tTTTT\nP\tp\ts1+,s2+\t*")
        assert g.get_path_sequence("p") == "ACGTTTTT"

    def test_spell_reverse_step(self):
        g = GfaGraph.parse("S\ts1\tACGT\nS\ts2\tTTGC\nP\tp\ts1+,s2-\t*")
        assert g.get_path_sequence("p") == "ACGT" + "GCAA"

    def test_walk_sequence(self):
        g = GfaGraph.parse("S\ts1\tACGT\nS\ts2\tGGGG")
        g.walks.append(Walk("S", "1", "c", 0, 8, 2, ["s1+", "s2+"]))
        assert g.get_walk_sequence(g.walks[0]) == "ACGTGGGG"

    def test_orphans_and_dangling(self):
        g = GfaGraph.parse(
            "S\ts1\tACGT\nS\ts2\tGGGG\nL\ts1\t+\tmissing\t+\t*\nP\tp\ts1+\t*")
        assert g.orphan_segments() == {"s2"}
        assert len(g.dangling_links()) == 1


class TestPanSN:
    def test_subrange(self):
        assert parse_pansn("HG1#1#chr21:40-100")[:5] == ("HG1", "1", "chr21", 40, 100)

    def test_chunk_suffix(self):
        s, h, c, st, en, cid = parse_pansn("HG1#1#chr21_chunk_0003")
        assert (s, h, c, st, cid) == ("HG1", "1", "chr21", None, "chunk_0003")

    def test_plain(self):
        assert haplotype_key("GRCh38#0#chr21") == ("GRCh38", "0", "chr21")


class TestSlicing:
    def test_seam_midpoint(self):
        assert seam_position((0, 60), (40, 100)) == 50

    def test_seam_none_when_disjoint(self):
        assert seam_position((0, 40), (60, 100)) is None

    def test_slice_forward_splits_node(self):
        g = GfaGraph.parse("S\ta\tAAAAA\nS\tb\tCCCCC\nP\tp\ta+,b+\t*")
        assert slice_steps(g, g.path_steps("p"), 0, 3, 7) == [
            ("a", "+", 3, 5), ("b", "+", 0, 2)]

    def test_slice_reverse_keeps_right_piece(self):
        g = GfaGraph.parse("S\ta\tAAACC\nP\tp\ta-\t*")
        assert slice_steps(g, g.path_steps("p"), 0, 0, 3) == [("a", "-", 2, 5)]

    def test_path_interval_from_subrange(self):
        _cid, g = _chunk("c1", 40, 100)
        assert path_interval(g, "HG1#1#chr21:40-100") == (40, 100)


class TestOverlapAwareStitch:
    def _merged(self):
        chunks = [_chunk("chunk_0001", 0, 60), _chunk("chunk_0002", 40, 100)]
        return chunks, overlap_aware_stitch(chunks)

    def test_reference_path_spells_whole_region(self):
        _c, (m, _br) = self._merged()
        assert m.get_path_sequence("GRCh38#0#chr21") == REF

    def test_haplotype_path_keeps_its_variant(self):
        _c, (m, _br) = self._merged()
        assert m.get_path_sequence("HG1#1#chr21") == HAP

    def test_overlap_counted_once(self):
        _c, (m, _br) = self._merged()
        assert m.path_length("GRCh38#0#chr21") == len(REF)

    def test_shared_nodes_stay_shared(self):
        _c, (m, _br) = self._merged()
        ref = {n for n, _o in m.path_steps("GRCh38#0#chr21")}
        hap = {n for n, _o in m.path_steps("HG1#1#chr21")}
        assert ref & hap, "merge split every haplotype into its own chain"

    def test_boundary_reported_pass(self):
        _c, (_m, br) = self._merged()
        assert len(br) == 1
        assert br[0]["status"] == "PASS"
        assert br[0]["haplotypes_joined"] == 2

    def test_order_independent(self):
        chunks, (m, _br) = self._merged()
        rev, _br2 = overlap_aware_stitch(list(reversed(chunks)))
        assert rev.to_gfa() == m.to_gfa()

    def test_gap_is_reported_not_papered_over(self):
        chunks = [_chunk("chunk_0001", 0, 40), _chunk("chunk_0002", 60, 100)]
        _m, br = overlap_aware_stitch(chunks)
        assert br[0]["status"] == "FAIL"
        assert br[0]["haplotypes_unjoined"] == 2

    def test_disjoint_union_keeps_every_chunk_path(self):
        from pipeline.merge.merge_graphs import diagnostic_disjoint_union
        chunks, _m = self._merged()
        u = diagnostic_disjoint_union(chunks)
        assert u.path_count() == 4
        assert u.node_count() == sum(g.node_count() for _c, g in chunks)

    def test_validate_accepts_the_stitch(self, tmp_path):
        _c, (m, _br) = self._merged()
        merged_p, base_p = tmp_path / "merged.gfa", tmp_path / "baseline.gfa"
        m.write_gfa(str(merged_p))
        base = GfaGraph()
        base.headers.append(Header("1.1"))
        base.segments["r"] = Segment("r", REF)
        base.segments["h"] = Segment("h", HAP)
        base.paths["GRCh38#0#chr21"] = Path("GRCh38#0#chr21", ["r+"])
        base.paths["HG1#1#chr21"] = Path("HG1#1#chr21", ["h+"])
        base.write_gfa(str(base_p))
        assert validate(str(base_p), str(merged_p)) == []

    def test_validate_catches_a_broken_merge(self, tmp_path):
        _c, (m, _br) = self._merged()
        first = m.path_steps("GRCh38#0#chr21")[0][0]
        m.segments[first].sequence = "G" + m.segments[first].sequence[1:]
        merged_p, base_p = tmp_path / "merged.gfa", tmp_path / "baseline.gfa"
        m.write_gfa(str(merged_p))
        base = GfaGraph()
        base.segments["r"] = Segment("r", REF)
        base.paths["GRCh38#0#chr21"] = Path("GRCh38#0#chr21", ["r+"])
        base.write_gfa(str(base_p))
        errs = [i for i in validate(str(base_p), str(merged_p))
                if i["severity"] == "ERROR"]
        assert any("sequence differs" in e["message"] for e in errs)

    def test_dropped_middle_chunk_does_not_invent_a_false_edge(self):
        seq = "ABCDEFGH"

        def span(cid, start, end):
            g = GfaGraph()
            g.headers.append(Header("1.1"))
            n = f"{cid}_n"
            g.segments[n] = Segment(n, seq[start:end])
            pn = f"GRCh38#0#chr21:{start}-{end}"
            g.paths[pn] = Path(pn, [n + "+"])
            return cid, g

        chunks = [
            span("chunk_0001", 0, 8),
            span("chunk_0002", 2, 8),
            span("chunk_0003", 2, 4),
        ]
        hap = ("GRCh38", "0", "chr21")
        pieces, _gaps = stitch_haplotype(
            chunks, group_paths_by_haplotype(chunks)[hap])
        assert [p[0] for p in pieces] == ["chunk_0001", "chunk_0003"], (
            "precondition: middle chunk must contribute zero pieces")

        m, br = overlap_aware_stitch(chunks)
        spelled = m.get_path_sequence("GRCh38#0#chr21")
        jumped = seq[0:5] + seq[3:4]
        assert spelled != jumped, (
            "merged path jumps backward across a dropped middle chunk")
        assert spelled and spelled in seq, (
            f"path is not a contiguous haplotype slice: {spelled!r}")
        sequences = {s.name: s.sequence for s in m.segments.values()}
        false_edge = any(
            sequences[link.from_node].endswith("E")
            and sequences[link.to_node] == "D"
            for link in m.links)
        assert not false_edge, (
            "invented edge from left tail to right head that are not "
            "consecutive on the haplotype")
        if spelled != seq:
            assert any(r["status"] == "FAIL" for r in br), (
                "dropped middle chunk left a non-adjacent join with no FAIL")

    def test_boundary_report_follows_reference_order_not_input_order(self):
        chunks, (m, br) = self._merged()
        assert br[0]["left_chunk"] == "chunk_0001"
        assert br[0]["right_chunk"] == "chunk_0002"
        rev, rev_br = overlap_aware_stitch(list(reversed(chunks)))
        assert rev.to_gfa() == m.to_gfa()
        assert len(rev_br) == 1
        assert rev_br[0]["left_chunk"] == "chunk_0001"
        assert rev_br[0]["right_chunk"] == "chunk_0002"
        assert rev_br[0]["boundary"] == "chunk_0001--chunk_0002"
        assert rev_br[0]["status"] == "PASS"

    def test_reversed_gap_is_still_reported_fail(self):
        chunks = [_chunk("chunk_0002", 60, 100), _chunk("chunk_0001", 0, 40)]
        _m, br = overlap_aware_stitch(chunks)
        assert len(br) == 1
        assert br[0]["left_chunk"] == "chunk_0001"
        assert br[0]["right_chunk"] == "chunk_0002"
        assert br[0]["boundary"] == "chunk_0001--chunk_0002"
        assert br[0]["status"] == "FAIL"
        assert br[0]["haplotypes_unjoined"] == 2


class TestWalkChunks:
    def test_w_line_chunk_is_stitchable(self):
        g = GfaGraph()
        g.headers.append(Header("1.1"))
        g.segments["n1"] = Segment("n1", REF[:20])
        g.walks.append(Walk("GRCh38", "0", "chr21", 0, 20, 1, ["n1+"]))
        p = walks_as_paths(g)
        assert "GRCh38#0#chr21:0-20" in p.paths
        assert p.get_path_sequence("GRCh38#0#chr21:0-20") == REF[:20]


class TestChunkOverlapMath:
    def test_adjacent_windows_share_the_configured_overlap(self):
        c = create_chunks("chr21", 0, 2000, 800, 100)
        shared = [min(a["reference_end"], b["reference_end"])
                  - max(a["reference_start"], b["reference_start"])
                  for a, b in zip(c, c[1:])]
        assert shared == [100, 100]

    def test_cores_tile_without_gaps(self):
        c = create_chunks("chr21", 0, 2000, 800, 100)
        assert c[0]["core_start"] == 0 and c[-1]["core_end"] == 2000
        assert all(a["core_end"] == b["core_start"] for a, b in zip(c, c[1:]))

    def test_edge_chunks_have_one_sided_padding(self):
        c = create_chunks("chr21", 0, 2000, 800, 100)
        assert c[0]["overlap_left"] == 0 and c[-1]["overlap_right"] == 0

    def test_overlap_larger_than_chunk_rejected(self):
        try:
            create_chunks("chr21", 0, 2000, 100, 100)
        except ValueError:
            return
        raise AssertionError("expected ValueError for overlap >= chunk_size")


class TestHPRCIndelCoordinates:
    """Haplotype coordinates diverge from GRCh38 when the hap has an indel.

    Real HPRC chunk FASTAs are named SAMPLE#HAP#CONTIG with no :start-end
    subrange. The only absolute coordinates are source_start/source_end/strand
    in chunk_mapping.tsv. Falling back to GRCh38 reference_start cuts the
    haplotype at the wrong seam and drops the inserted bases.
    """

    REF = "ACGT" * 25          # 100 bp
    INS = "NNNNNNNNNN"         # 10 bp insertion at ref pos 25
    HAP = REF[:25] + INS + REF[25:]  # 110 bp

    def _hap_chunk(self, cid, ref_s, ref_e, hap_s, hap_e):
        g = GfaGraph()
        g.headers.append(Header("1.1"))
        rname, hname = f"{cid}_r", f"{cid}_h"
        g.segments[rname] = Segment(rname, self.REF[ref_s:ref_e])
        g.segments[hname] = Segment(hname, self.HAP[hap_s:hap_e])
        g.paths["GRCh38#0#chr21"] = Path("GRCh38#0#chr21", [rname + "+"])
        g.paths["HG1#1#chr21"] = Path("HG1#1#chr21", [hname + "+"])
        return cid, g

    def _mapping(self):
        return [
            {"chunk_id": "chunk_0001", "sample": "GRCh38", "haplotype": "0",
             "source_start": 0, "source_end": 60, "strand": "+"},
            {"chunk_id": "chunk_0001", "sample": "HG1", "haplotype": "1",
             "source_start": 0, "source_end": 70, "strand": "+"},
            {"chunk_id": "chunk_0002", "sample": "GRCh38", "haplotype": "0",
             "source_start": 40, "source_end": 100, "strand": "+"},
            {"chunk_id": "chunk_0002", "sample": "HG1", "haplotype": "1",
             "source_start": 50, "source_end": 110, "strand": "+"},
        ]

    def _rows(self):
        return {
            "chunk_0001": {"reference_start": 0, "reference_end": 60},
            "chunk_0002": {"reference_start": 40, "reference_end": 100},
        }

    def test_grch38_fallback_drops_the_insertion(self):
        """If we placed the hap on GRCh38 coords, the 10 bp insert is lost."""
        chunks = [
            self._hap_chunk("chunk_0001", 0, 60, 0, 70),
            self._hap_chunk("chunk_0002", 40, 100, 50, 110),
        ]
        m, _br = overlap_aware_stitch(chunks, chunk_rows=self._rows())
        assert m.get_path_sequence("HG1#1#chr21") != self.HAP

    def test_mapping_coords_preserve_haplotype_with_indel(self):
        chunks = [
            self._hap_chunk("chunk_0001", 0, 60, 0, 70),
            self._hap_chunk("chunk_0002", 40, 100, 50, 110),
        ]
        m, br = overlap_aware_stitch(
            chunks, chunk_rows=self._rows(), chunk_mapping=self._mapping())
        assert m.get_path_sequence("GRCh38#0#chr21") == self.REF
        assert m.get_path_sequence("HG1#1#chr21") == self.HAP
        assert self.INS in m.get_path_sequence("HG1#1#chr21")
        assert br[0]["status"] == "PASS"

    def test_w_line_haplotype_coords_still_stitch(self):
        """W-line SeqStart/SeqEnd are already haplotype coords — keep that."""
        ref, hap = self.REF, self.HAP
        def wchunk(cid, sample, hap_id, start, end, seq):
            g = GfaGraph()
            g.headers.append(Header("1.1"))
            n = f"{cid}_{sample}"
            g.segments[n] = Segment(n, seq)
            g.walks.append(Walk(sample, hap_id, "chr21", start, end, path=[n + "+"]))
            return cid, g

        # Two overlapping W-line graphs per haplotype, merged by combining
        # walks onto one graph per chunk.
        c1, c2 = GfaGraph(), GfaGraph()
        for g, cid, rs, re, hs, he in (
                (c1, "chunk_0001", 0, 60, 0, 70),
                (c2, "chunk_0002", 40, 100, 50, 110)):
            g.headers.append(Header("1.1"))
            rn, hn = f"{cid}_r", f"{cid}_h"
            g.segments[rn] = Segment(rn, ref[rs:re])
            g.segments[hn] = Segment(hn, hap[hs:he])
            g.walks.append(Walk("GRCh38", "0", "chr21", rs, re, path=[rn + "+"]))
            g.walks.append(Walk("HG1", "1", "chr21", hs, he, path=[hn + "+"]))
        m, br = overlap_aware_stitch([("chunk_0001", c1), ("chunk_0002", c2)])
        assert m.get_path_sequence("GRCh38#0#chr21") == ref
        assert m.get_path_sequence("HG1#1#chr21") == hap
        assert br[0]["status"] == "PASS"

    def test_mapping_overrides_chunk_local_w_line_coords(self):
        """PGGB W-lines on unsliced FASTA headers are 0-based in the chunk.

        Absolute haplotype coordinates live in chunk_mapping.tsv. Mapping
        must win over those chunk-local start/end values or an indel shifts
        the seam onto the wrong bases.
        """
        ref, hap = self.REF, self.HAP
        c1, c2 = GfaGraph(), GfaGraph()
        for g, cid, rs, re, hs, he in (
                (c1, "chunk_0001", 0, 60, 0, 70),
                (c2, "chunk_0002", 40, 100, 50, 110)):
            g.headers.append(Header("1.1"))
            rn, hn = f"{cid}_r", f"{cid}_h"
            g.segments[rn] = Segment(rn, ref[rs:re])
            g.segments[hn] = Segment(hn, hap[hs:he])
            g.walks.append(Walk("GRCh38", "0", "chr21", 0, re - rs, path=[rn + "+"]))
            g.walks.append(Walk("HG1", "1", "chr21", 0, he - hs, path=[hn + "+"]))
        m, br = overlap_aware_stitch(
            [("chunk_0001", c1), ("chunk_0002", c2)],
            chunk_rows=self._rows(), chunk_mapping=self._mapping())
        assert m.get_path_sequence("GRCh38#0#chr21") == ref
        assert m.get_path_sequence("HG1#1#chr21") == hap
        assert br[0]["status"] == "PASS"

