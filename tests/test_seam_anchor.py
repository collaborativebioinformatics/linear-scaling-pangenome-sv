"""Seam cut: all haplotypes are cut at one homologous anchor node.

The original cut put each haplotype at the midpoint of its OWN overlap. Overlap
edges come from separate mappings of the reference window, so they differ by
tens of bp between haplotypes. When the cut lands inside a node every haplotype
shares, each one slices that node at a different offset and ends up on its own
copy of identical sequence. These tests rebuild that situation on a toy locus.
"""
import os
import random
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.merge.gfa import GfaGraph
from pipeline.merge.merge_graphs import overlap_aware_stitch
from pipeline.merge.paths import anchor_cuts

random.seed(7)
REF = "".join(random.choice("ACGT") for _ in range(600))

# haplotype coordinate = reference coordinate + shift (contig offset plus an
# upstream indel)
SHIFT = {"S1": 1000, "S2": 995, "S3": 1003}
# per-haplotype reference window that each chunk's mapping returned. The edges
# differ between haplotypes: this is the mapping fuzz that moves the midpoints.
WIN_A = {"S1": (0, 260), "S2": (0, 266), "S3": (0, 255)}
WIN_B = {"S1": (200, 500), "S2": (207, 506), "S3": (194, 494)}
WIN_C = {"S1": (440, 600), "S2": (447, 600), "S3": (434, 600)}


def _chunk(prefix, windows, extra_breaks=(), private=(), snp=None):
    """A chunk graph over REF, split at every window edge and extra break.

    `private` haplotypes get their own copy of every node (same sequence,
    different segment names), so they share no node with anyone.
    `snp` = (hap, position, base) gives that haplotype a one-base alt allele.
    """
    breaks = sorted({e for w in windows.values() for e in w} | set(extra_breaks))
    pieces = list(zip(breaks, breaks[1:]))
    lines = ["H\tVN:Z:1.1"]
    names = {}
    for lo, hi in pieces:
        names[(None, lo)] = f"{prefix}{lo}"
        lines.append(f"S\t{prefix}{lo}\t{REF[lo:hi]}")
    for hap in private:
        for lo, hi in pieces:
            names[(hap, lo)] = f"{prefix}{hap}_{lo}"
            lines.append(f"S\t{prefix}{hap}_{lo}\t{REF[lo:hi]}")
    if snp:
        _h, pos, base = snp
        lines.append(f"S\t{prefix}alt\t{base}")
    for hap, (ws, we) in windows.items():
        steps = []
        for lo, hi in pieces:
            if lo >= ws and hi <= we:
                if snp and hap == snp[0] and lo == snp[1]:
                    steps.append(f"{prefix}alt+")
                else:
                    steps.append(names[(hap if hap in private else None, lo)] + "+")
        sh = SHIFT[hap]
        lines.append(f"P\t{hap}#1#chr1:{ws + sh}-{we + sh}\t{','.join(steps)}\t*")
    return GfaGraph.parse("\n".join(lines) + "\n")


def _two(private=(), snp=None):
    b_breaks = (350,) if not snp else (350, snp[1], snp[1] + 1)
    return [("chunk_0001", _chunk("a", WIN_A, private=private)),
            ("chunk_0002", _chunk("b", WIN_B, b_breaks, private, snp))]


def _three():
    return [("chunk_0001", _chunk("a", WIN_A)),
            ("chunk_0002", _chunk("b", WIN_B, (350,))),
            ("chunk_0003", _chunk("c", WIN_C, (520,)))]


def _seq(m, hap):
    return m.get_path_sequence(f"{hap}#1#chr1")


def _expected(hap, end, snp=None):
    s = REF[:end]
    if snp and hap == snp[0]:
        s = s[:snp[1]] + snp[2] + s[snp[1] + 1:]
    return s


class TestAnchorSeam:
    def test_every_haplotype_still_spells_its_sequence(self):
        # S1, S2 and S3 reach 500, 506 and 494 in chunk B
        for mode in ("anchor", "midpoint"):
            m, _br = overlap_aware_stitch(_two(), seam_mode=mode)
            assert _seq(m, "S1") == REF[:500], mode
            assert _seq(m, "S2") == REF[:506], mode
            assert _seq(m, "S3") == REF[:494], mode

    def test_midpoint_cut_duplicates_shared_sequence(self):
        """The original behaviour: one copy of identical sequence per haplotype."""
        m, _br = overlap_aware_stitch(_two(), seam_mode="midpoint")
        assert m.total_sequence_bp() > 506

    def test_anchor_cut_keeps_shared_sequence_shared(self):
        m, br = overlap_aware_stitch(_two(), seam_mode="anchor")
        # the union of covered reference bases is [0, 506), each stored once
        # except the 12 bp tail only S2 reaches, which is its own node anyway
        assert m.total_sequence_bp() == 506
        assert br[0]["seam_mode"] == "anchor"
        assert br[0]["anchor_bp"] > 0
        assert br[0]["status"] == "PASS"

    def test_anchor_cut_adds_no_branches_when_nothing_differs(self):
        """All haplotypes carry identical sequence, so before the tail the graph
        is a line: no node has two different successors except at the ends."""
        m, _br = overlap_aware_stitch(_two(), seam_mode="anchor")
        out = {}
        for l in m.links:
            out.setdefault((l.from_node, l.from_orient), set()).add(l.to_node)
        branching = [k for k, v in out.items() if len(v) > 1]
        # the only branches are where haplotypes stop at different places
        assert len(branching) <= 2

    def test_snp_after_the_seam_is_the_only_bubble(self):
        snp = ("S3", 300, "A" if REF[300] != "A" else "C")
        m, br = overlap_aware_stitch(_two(snp=snp), seam_mode="anchor")
        assert m.total_sequence_bp() == 506 + 1
        assert _seq(m, "S3") == _expected("S3", 494, snp)
        assert _seq(m, "S1") == REF[:500]
        mid, _ = overlap_aware_stitch(_two(snp=snp), seam_mode="midpoint")
        assert mid.total_sequence_bp() > m.total_sequence_bp()

    def test_three_chunks_two_boundaries(self):
        m, br = overlap_aware_stitch(_three(), seam_mode="anchor")
        assert [b["seam_mode"] for b in br] == ["anchor", "anchor"]
        assert _seq(m, "S1") == REF[:600]
        assert _seq(m, "S2") == REF[:600]
        assert _seq(m, "S3") == REF[:600]
        # S1/S2/S3 differ only in how far each chunk window reaches, so the
        # merged graph stores the 600 reference bases once
        assert m.total_sequence_bp() == 600
        mid, _ = overlap_aware_stitch(_three(), seam_mode="midpoint")
        assert mid.total_sequence_bp() > 600

    def test_no_shared_node_falls_back_to_midpoint_and_says_so(self):
        m, br = overlap_aware_stitch(_two(private=("S3",)), seam_mode="anchor")
        assert br[0]["seam_mode"] == "midpoint"
        assert br[0]["anchor_bp"] == 0
        assert _seq(m, "S1") == REF[:500]
        assert _seq(m, "S3") == REF[:494]

    def test_default_is_anchor(self):
        _m, br = overlap_aware_stitch(_two())
        assert br[0]["seam_mode"] == "anchor"

    def test_bad_mode_is_rejected(self):
        with pytest.raises(ValueError):
            overlap_aware_stitch(_two(), seam_mode="nope")


class TestAnchorCuts:
    def _members(self):
        g = _chunk("a", WIN_A)
        out = []
        for hap, (ws, we) in WIN_A.items():
            b0 = WIN_B[hap][0]
            sh = SHIFT[hap]
            out.append((hap, f"{hap}#1#chr1:{ws + sh}-{we + sh}",
                        (ws + sh, we + sh), (b0 + sh, we + sh)))
        return g, out

    def test_cut_is_the_same_reference_base_for_everyone(self):
        g, members = self._members()
        node, cuts = anchor_cuts(g, members)
        ref_cut = {hap: cuts[hap] - SHIFT[hap] for hap in cuts}
        assert len(set(ref_cut.values())) == 1, ref_cut

    def test_candidate_must_end_inside_every_overlap(self):
        g, members = self._members()
        # shrink S3's overlap so it ends before the only shared node does
        key, pn, iv, (lo, hi) = members[2]
        members[2] = (key, pn, iv, (lo, lo + 5))
        assert anchor_cuts(g, members) is None

    def test_min_anchor_bp_filters_short_nodes(self):
        g, members = self._members()
        assert anchor_cuts(g, members, min_anchor_bp=10_000) is None
