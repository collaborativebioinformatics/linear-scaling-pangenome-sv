"""GFA data model for pangenome graph manipulation.

Core types: Header, Segment, Link, Path, Walk, GfaGraph.
Supports GFA v1 (P-lines) and GFA v1.1 (W-lines for walks).
Zero external dependencies.
"""

from __future__ import annotations
import hashlib
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

TAB = chr(9)
NL = chr(10)


def _parse_tags(fields, start=3):
    return {f.split(":")[0]: f for f in fields[start:] if ":" in f}


def _split_orient(seg_str):
    if seg_str and seg_str[-1] in ("+", "-"):
        return seg_str[:-1], seg_str[-1]
    return seg_str, "+"



def _walk_to_segments(walk_str):
    """Parse an official GFA 1.1 W-line walk: >n1<n2>n3 -> [n1+, n2-, n3+]."""
    if not walk_str:
        return []
    segs = []
    i = 0
    while i < len(walk_str):
        if walk_str[i] == ">":
            orient = "+"
            i += 1
        elif walk_str[i] == "<":
            orient = "-"
            i += 1
        else:
            # Tolerate a legacy comma-separated walk (pre-spec output).
            return [s for s in walk_str.split(",") if s]
        j = i
        while j < len(walk_str) and walk_str[j] not in "><":
            j += 1
        name = walk_str[i:j]
        if name:
            segs.append(name + orient)
        i = j
    return segs


def _segments_to_walk(segments):
    """Render [n1+, n2-] back to the official >n1<n2 form."""
    out = []
    for seg in segments:
        name, orient = _split_orient(seg)
        out.append(f">{name}" if orient == "+" else f"<{name}")
    return "".join(out)


@dataclass
class Header:
    version: str = "1.1"
    metadata: Dict[str, str] = field(default_factory=dict)

    def to_gfa(self):
        parts = ["H", f"VN:Z:{self.version}"]
        for k, v in self.metadata.items():
            parts.append(f"{k}:Z:{v}")
        return TAB.join(parts)

    @classmethod
    def parse(cls, fields):
        h = cls()
        for field in fields[1:]:
            if field.startswith("VN:Z:"):
                h.version = field[5:]
            elif ":" in field:
                tag, _, val = field.partition(":Z:")
                h.metadata[tag or field.partition(":")[0]] = val or field
        return h


@dataclass
class Segment:
    name: str
    sequence: str
    length: int = 0
    tags: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        if self.length == 0:
            self.length = len(self.sequence)

    def to_gfa(self):
        seq = self.sequence if self.sequence else "*"
        parts = ["S", self.name, seq]
        for k, v in self.tags.items():
            parts.append(f"{k}:Z:{v}")
        return TAB.join(parts)

    @classmethod
    def parse(cls, fields):
        return cls(
            name=fields[1],
            sequence=fields[2] if fields[2] != "*" else "",
            tags=_parse_tags(fields, 3),
        )

    @property
    def sha256(self):
        return hashlib.sha256(self.sequence.encode()).hexdigest()[:16]


@dataclass
class Link:
    from_node: str
    from_orient: str
    to_node: str
    to_orient: str
    overlap: str = "*"
    tags: Dict[str, str] = field(default_factory=dict)

    def to_gfa(self):
        parts = ["L", self.from_node, self.from_orient,
                 self.to_node, self.to_orient, self.overlap]
        for k, v in self.tags.items():
            parts.append(f"{k}:Z:{v}")
        return TAB.join(parts)

    @classmethod
    def parse(cls, fields):
        return cls(
            from_node=fields[1],
            from_orient=fields[2],
            to_node=fields[3],
            to_orient=fields[4],
            overlap=fields[5] if len(fields) > 5 else "*",
            tags=_parse_tags(fields, 6),
        )


@dataclass
class Path:
    path_name: str
    segment_names: List[str]
    overlaps: List[str] = field(default_factory=list)
    tags: Dict[str, str] = field(default_factory=dict)

    def to_gfa(self):
        seg_str = ",".join(self.segment_names)
        ov_str = ",".join(self.overlaps) if self.overlaps else "*"
        parts = ["P", self.path_name, seg_str, ov_str]
        for k, v in self.tags.items():
            parts.append(f"{k}:Z:{v}")
        return TAB.join(parts)

    @classmethod
    def parse(cls, fields):
        seg_raw = fields[2] if fields[2] != "*" else ""
        ov_raw = fields[3] if len(fields) > 3 and fields[3] != "*" else ""
        return cls(
            path_name=fields[1],
            segment_names=seg_raw.split(",") if seg_raw else [],
            overlaps=ov_raw.split(",") if ov_raw else [],
            tags=_parse_tags(fields, 4),
        )


@dataclass(init=False, repr=True, eq=True)
class Walk:
    """Official GFA 1.1 W-line.

    Spec:  W  SampleId  HapIndex  SeqId  SeqStart  SeqEnd  Walk
    The walk field uses > / < orientation notation, NOT a comma list, and
    there is NO step_count column.

    The constructor accepts BOTH call shapes so existing code keeps working:
        Walk(sample, hap, contig, start, end, path)              # spec
        Walk(sample, hap, contig, start, end, step_count, path)  # legacy
    In the legacy form step_count is accepted and ignored; the authoritative
    value is always len(path), exposed as the .step_count property.
    """
    sample: str
    haplotype: str
    contig: str
    start: int
    end: int
    path: List[str]
    tags: Dict[str, str] = field(default_factory=dict)

    def __init__(self, sample, haplotype, contig, start, end,
                 *args, tags=None, path=None, step_count=None):
        self.sample = sample
        self.haplotype = haplotype
        self.contig = contig
        self.start = start
        self.end = end
        resolved_path = path
        if args:
            # Legacy positional form put step_count before path. Detect by
            # type: an int in slot 6 is step_count, a list is the path.
            if isinstance(args[0], int) and not isinstance(args[0], bool):
                if len(args) > 1:
                    resolved_path = args[1]
            else:
                resolved_path = args[0]
            if len(args) > 2 and tags is None:
                tags = args[2]
        self.path = list(resolved_path) if resolved_path else []
        self.tags = tags if tags is not None else {}

    @property
    def step_count(self) -> int:
        """Derived, never stored — always consistent with the actual path."""
        return len(self.path)

    def to_gfa(self):
        parts = ["W", self.sample, self.haplotype, self.contig,
                 str(self.start), str(self.end), _segments_to_walk(self.path)]
        for k, v in self.tags.items():
            parts.append(f"{k}:Z:{v}")
        return TAB.join(parts)

    @classmethod
    def parse(cls, fields):
        wr = fields[6] if len(fields) > 6 and fields[6] != "*" else ""
        return cls(
            fields[1], fields[2], fields[3],
            int(fields[4]), int(fields[5]),
            path=_walk_to_segments(wr),
            tags=_parse_tags(fields, 7),
        )


class GfaGraph:
    def __init__(self):
        self.headers: List[Header] = []
        self.segments: Dict[str, Segment] = {}
        self.links: List[Link] = []
        self.paths: Dict[str, Path] = {}
        self.walks: List[Walk] = []
        self.source: Optional[str] = None

    @classmethod
    def parse(cls, gfa_str, source=None):
        graph = cls()
        graph.source = source
        dispatch = {
            "H": lambda f: graph.headers.append(Header.parse(f)),
            "S": lambda f: graph.segments.__setitem__(f[1], Segment.parse(f)),
            "L": lambda f: graph.links.append(Link.parse(f)),
            "P": lambda f: graph.paths.__setitem__(f[1], Path.parse(f)),
            "W": lambda f: graph.walks.append(Walk.parse(f)),
        }
        for ln, line in enumerate(gfa_str.strip().split(NL), 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split(TAB)
            if not fields:
                continue
            parser = dispatch.get(fields[0])
            if parser:
                try:
                    parser(fields)
                except Exception as e:
                    raise ValueError(f"Line {ln}: {e}") from e
        return graph

    @classmethod
    def parse_file(cls, path):
        with open(path) as f:
            return cls.parse(f.read(), source=path)

    def to_gfa(self):
        lines = []
        lines.extend(h.to_gfa() for h in self.headers)
        lines.extend(s.to_gfa() for s in self.segments.values())
        lines.extend(l.to_gfa() for l in self.links)
        lines.extend(p.to_gfa() for p in self.paths.values())
        lines.extend(w.to_gfa() for w in self.walks)
        return NL.join(lines)

    def write_gfa(self, path):
        with open(path, "w") as f:
            f.write(self.to_gfa() + NL)

    def node_count(self): return len(self.segments)
    def edge_count(self): return len(self.links)
    def path_count(self): return len(self.paths)
    def walk_count(self): return len(self.walks)
    def total_sequence_bp(self): return sum(s.length for s in self.segments.values())

    def path_steps(self, path_name):
        """[(segment_name, orient)] for a P-line path."""
        p = self.paths[path_name]
        return [_split_orient(s) for s in p.segment_names]

    def walk_steps(self, walk):
        """[(segment_name, orient)] for a W-line walk."""
        return [_split_orient(s) for s in walk.path]

    def segment_sequence(self, name, orient="+"):
        seg = self.segments.get(name)
        if seg is None:
            raise KeyError(f"segment not in graph: {name}")
        return seg.sequence if orient == "+" else _revcomp(seg.sequence)

    def get_path_sequence(self, path_name):
        """Spell out a P-line path. Raises KeyError on a dangling step."""
        return "".join(self.segment_sequence(n, o)
                       for n, o in self.path_steps(path_name))

    def get_walk_sequence(self, walk):
        """Spell out a W-line walk. Raises KeyError on a dangling step."""
        return "".join(self.segment_sequence(n, o)
                       for n, o in self.walk_steps(walk))

    def path_length(self, path_name):
        return sum(self.segments[n].length for n, _ in self.path_steps(path_name))

    def used_segments(self):
        """Names of segments touched by any path or walk."""
        used = set()
        for pn in self.paths:
            used.update(n for n, _ in self.path_steps(pn))
        for w in self.walks:
            used.update(n for n, _ in self.walk_steps(w))
        return used

    def orphan_segments(self):
        """Segments on no path and no walk."""
        return set(self.segments) - self.used_segments()

    def dangling_links(self):
        """Links whose endpoints are not segments in this graph."""
        return [l for l in self.links
                if l.from_node not in self.segments or l.to_node not in self.segments]

    def get_sample_names(self):
        samples = {w.sample for w in self.walks}
        for p in self.paths.values():
            n = p.path_name.split("#")[0] if "#" in p.path_name else p.path_name
            samples.add(n)
        return samples

    def copy(self):
        g = GfaGraph()
        g.headers = [Header(h.version, dict(h.metadata)) for h in self.headers]
        g.segments = {k: Segment(v.name, v.sequence, v.length, dict(v.tags))
                       for k, v in self.segments.items()}
        g.links = [Link(l.from_node, l.from_orient, l.to_node,
                        l.to_orient, l.overlap, dict(l.tags)) for l in self.links]
        g.paths = {k: Path(v.path_name, list(v.segment_names),
                           list(v.overlaps), dict(v.tags))
                   for k, v in self.paths.items()}
        g.walks = [Walk(w.sample, w.haplotype, w.contig, w.start, w.end,
                        w.step_count, list(w.path), dict(w.tags)) for w in self.walks]
        g.source = self.source
        return g


def _revcomp(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G",
            "a": "t", "t": "a", "g": "c", "c": "g",
            "N": "N", "n": "n"}
    return "".join(comp.get(c, c) for c in reversed(seq))


def extract_chromosome_from_gfa(graph, ref_name="GRCh38"):
    for w in graph.walks:
        if ref_name in w.sample:
            return w.contig
    for p in graph.paths:
        if ref_name in p:
            parts = p.split("#")
            if len(parts) >= 3: return parts[2]
    return None


def infer_data_mode(graph):
    for seg in graph.segments.values():
        s = seg.sequence.upper()
        if not s or len(s) <= 20: continue
        if len(set(s)) <= 3: return "synthetic"
        gc = s.count("G") + s.count("C")
        if gc / len(s) < 0.2 or gc / len(s) > 0.8: return "synthetic"
    return "real"


_SUBRANGE = re.compile(r"^(.*):(\d+)-(\d+)$")
_CHUNK_SUFFIX = re.compile(r"_chunk_\d+$")


def parse_pansn(name):
    """Split a PanSN path name into its parts.

    Accepts `SAMPLE#HAP#CONTIG`, an optional `:start-end` subrange (what PGGB
    emits when the input FASTA was sliced), and the `_chunk_NNNN` suffix the
    synthetic demo appends. Returns
    (sample, haplotype, contig, start, end, chunk_id); start/end are None when
    the name carries no subrange, chunk_id is None when it carries no suffix.
    """
    chunk_id = None
    m = _CHUNK_SUFFIX.search(name)
    if m:
        chunk_id = m.group(0)[1:]
        name = name[:m.start()]
    start = end = None
    m = _SUBRANGE.match(name)
    if m:
        name, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    parts = name.split("#")
    if len(parts) >= 3:
        sample, hap, contig = parts[0], parts[1], "#".join(parts[2:])
    elif len(parts) == 2:
        sample, hap, contig = parts[0], parts[1], parts[1]
    else:
        sample, hap, contig = name, "0", name
    return sample, hap, contig, start, end, chunk_id


def haplotype_key(name):
    """The (sample, haplotype, contig) a chunk-local path belongs to."""
    sample, hap, contig, _s, _e, _c = parse_pansn(name)
    return sample, hap, contig
