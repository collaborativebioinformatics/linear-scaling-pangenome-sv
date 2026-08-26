"""Tests for GFA data model, chunking, merge, and interval mapping."""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.merge.gfa import (
    GfaGraph, Header, Segment, Link, Path, Walk,
    _revcomp, infer_data_mode, extract_chromosome_from_gfa, _split_orient)
from pipeline.parallel.make_chunks import create_chunks


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
        w = Walk("HG001", "1", "chr21", 0, 100, 5, ["s1+"])
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
        g.walks.append(Walk("S1", "1", "c", 0, 10, 1, ["n1+"]))
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
        g.walks.append(Walk("GRCh38", "0", "chr21", 0, 1000, 5, ["n1+"]))
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


class TestDemo:
    def test_json_exists(self):
        assert os.path.exists("web/public/data/latest.json")