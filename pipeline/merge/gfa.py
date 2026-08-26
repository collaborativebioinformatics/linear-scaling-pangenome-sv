"""GFA data model for pangenome graph manipulation.
Supports GFA v1 (P-lines) and GFA v1.1 (W-lines for walks).
Official W format: W SampleId HapIndex SeqId SeqStart SeqEnd Walk
Walk string: >node1<node2>node3 (no step_count column)
Zero external dependencies.
"""
from __future__ import annotations
import hashlib
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set
TAB = chr(9)
NL = chr(10)

def _parse_tags(fields, start=3):
    return {f.split(":")[0]: f for f in fields[start:] if ":" in f}

def _split_orient(seg_str):
    if seg_str and seg_str[-1] in ("+", "-"): return seg_str[:-1], seg_str[-1]
    return seg_str, "+"

def _walk_to_segments(walk_str):
    """Parse official W-line: >node1<node2>node3 -> [node1+, node2-]"""
    if not walk_str: return []
    segs = []; i = 0
    while i < len(walk_str):
        if walk_str[i] == ">": orient = "+"; i += 1
        elif walk_str[i] == "<": orient = "-"; i += 1
        else:
            name = ""
            while i < len(walk_str) and walk_str[i] not in (">", "<"): name += walk_str[i]; i += 1
            if name: segs.append(name + orient)
    return segs

def _segments_to_walk(segments):
    result = []
    for seg in segments:
        name, orient = _split_orient(seg)
        result.append(f">{name}" if orient == "+" else f"<{name}")
    return "".join(result)

@dataclass
class Header:
    version: str = "1.1"
    metadata: Dict[str, str] = field(default_factory=dict)
    def to_gfa(self): parts = ["H",f"VN:Z:{self.version}"] + [f"{k}:Z:{v}" for k,v in self.metadata.items()]; return TAB.join(parts)
    @classmethod
    def parse(cls, fields):
        h = cls()
        for field in fields[1:]:
            if field.startswith("VN:Z:"): h.version = field[5:]
            elif ":" in field: tag,_,val = field.partition(":Z:"); h.metadata[tag or field.partition(":")[0]] = val or field
        return h

@dataclass
class Segment:
    name: str; sequence: str; length: int = 0
    tags: Dict[str, str] = field(default_factory=dict)
    def __post_init__(self):
        if self.length == 0: self.length = len(self.sequence)
    def to_gfa(self): seq = self.sequence if self.sequence else "*"; parts = ["S",self.name,seq] + [f"{k}:Z:{v}" for k,v in self.tags.items()]; return TAB.join(parts)
    @classmethod
    def parse(cls, fields): return cls(name=fields[1], sequence=fields[2] if fields[2]!="*" else "", tags=_parse_tags(fields,3))
    @property
    def sha256(self): return hashlib.sha256(self.sequence.encode()).hexdigest()[:16]

@dataclass
class Link:
    from_node: str; from_orient: str; to_node: str; to_orient: str
    overlap: str = "*"; tags: Dict[str, str] = field(default_factory=dict)
    def to_gfa(self): parts = ["L",self.from_node,self.from_orient,self.to_node,self.to_orient,self.overlap] + [f"{k}:Z:{v}" for k,v in self.tags.items()]; return TAB.join(parts)
    @classmethod
    def parse(cls, fields): return cls(from_node=fields[1], from_orient=fields[2], to_node=fields[3], to_orient=fields[4], overlap=fields[5] if len(fields)>5 else "*", tags=_parse_tags(fields,6))

@dataclass
class Path:
    path_name: str; segment_names: List[str]; overlaps: List[str] = field(default_factory=list); tags: Dict[str, str] = field(default_factory=dict)
    def to_gfa(self): seg_str=",".join(self.segment_names); ov_str=",".join(self.overlaps) if self.overlaps else "*"; parts=["P",self.path_name,seg_str,ov_str] + [f"{k}:Z:{v}" for k,v in self.tags.items()]; return TAB.join(parts)
    @classmethod
    def parse(cls, fields): sr = fields[2] if fields[2]!="*" else ""; ov = fields[3] if len(fields)>3 and fields[3]!="*" else ""; return cls(path_name=fields[1], segment_names=sr.split(",") if sr else [], overlaps=ov.split(",") if ov else [], tags=_parse_tags(fields,4))

@dataclass
class Walk:
    """Official GFA 1.1 W-line: W SampleId HapIndex SeqId SeqStart SeqEnd Walk (no step_count)"""
    sample: str; haplotype: str; contig: str; start: int; end: int; path: List[str]
    tags: Dict[str, str] = field(default_factory=dict)
    def to_gfa(self): ws = _segments_to_walk(self.path); parts = ["W",self.sample,self.haplotype,self.contig,str(self.start),str(self.end),ws] + [f"{k}:Z:{v}" for k,v in self.tags.items()]; return TAB.join(parts)
    @classmethod
    def parse(cls, fields): wr = fields[6] if len(fields)>6 and fields[6]!="*" else ""; return cls(sample=fields[1], haplotype=fields[2], contig=fields[3], start=int(fields[4]), end=int(fields[5]), path=_walk_to_segments(wr), tags=_parse_tags(fields,7))

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
        graph = cls(); graph.source = source
        dispatch = {"H":lambda f: graph.headers.append(Header.parse(f)), "S":lambda f: graph.segments.__setitem__(f[1],Segment.parse(f)), "L":lambda f: graph.links.append(Link.parse(f)), "P":lambda f: graph.paths.__setitem__(f[1],Path.parse(f)), "W":lambda f: graph.walks.append(Walk.parse(f))}
        for ln,line in enumerate(gfa_str.strip().split(NL),1):
            line=line.strip(); fields=line.split(TAB)
            if not line or line.startswith("#") or not fields: continue
            parser = dispatch.get(fields[0])
            if parser:
                try: parser(fields)
                except Exception as e: raise ValueError(f"Line {ln}: {e}") from e
        return graph
    @classmethod
    @classmethod
    def parse_file(cls, path):
        with open(path) as f: return cls.parse(f.read(), source=path)
    def to_gfa(self): lines = []; lines.extend(h.to_gfa() for h in self.headers); lines.extend(s.to_gfa() for s in self.segments.values()); lines.extend(l.to_gfa() for l in self.links); lines.extend(p.to_gfa() for p in self.paths.values()); lines.extend(w.to_gfa() for w in self.walks); return NL.join(lines)
    def write_gfa(self, path):
        with open(path,"w") as f: f.write(self.to_gfa()+NL)
    def node_count(self): return len(self.segments)
    def edge_count(self): return len(self.links)
    def path_count(self): return len(self.paths)
    def walk_count(self): return len(self.walks)
    def total_sequence_bp(self): return sum(s.length for s in self.segments.values())
    def get_sample_names(self): samples={w.sample for w in self.walks}; [samples.add(p.path_name.split("#")[0] if "#" in p.path_name else p.path_name) for p in self.paths.values()]; return samples
    def copy(self):
        g = GfaGraph()
        g.headers = [Header(h.version,dict(h.metadata)) for h in self.headers]
        g.segments = {k:Segment(v.name,v.sequence,v.length,dict(v.tags)) for k,v in self.segments.items()}
        g.links = [Link(l.from_node,l.from_orient,l.to_node,l.to_orient,l.overlap,dict(l.tags)) for l in self.links]
        g.paths = {k:Path(v.path_name,list(v.segment_names),list(v.overlaps),dict(v.tags)) for k,v in self.paths.items()}
        g.walks = [Walk(w.sample,w.haplotype,w.contig,w.start,w.end,list(w.path),dict(w.tags)) for w in self.walks]
        g.source = self.source; return g

def _revcomp(seq): c = {"A":"T","T":"A","G":"C","C":"G","a":"t","t":"a","g":"c","c":"g","N":"N","n":"n"}; return "".join(c.get(x,x) for x in reversed(seq))

def extract_chromosome_from_gfa(graph, ref_name="GRCh38"):
    for w in graph.walks:
        if ref_name in w.sample: return w.contig
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
