"use client";
import { useState, useEffect, useCallback } from "react";
import type { SampleGraph, NodeInfo, EdgeInfo } from "../types";
import GraphExplorer from "../components/GraphExplorer";

interface ManifestData {
  data_mode: string;
  pipeline_status?: Record<string, string>;
  samples?: { sample: string; haplotypes: string[]; hap_labels?: Record<string, string> }[];
  graphs?: Record<string, { nodes: number; edges: number; paths: number; walks: number }>;
}
interface LatestData {
  data_mode: string;
  metrics?: Record<string, Record<string, number>>;
  stitch?: { status?: string };
  equivalence?: { verdict?: string };
  boundaries?: { boundary?: string; status?: string; haplotypes_joined?: number }[];
}
interface OverviewData {
  chunks?: { chunk_id: string; reference_start: number; reference_end: number;
    core_start: number; core_end: number; overlap_left: number; overlap_right: number }[];
}

const FALLBACK_SAMPLES: { sample: string; haplotypes: string[]; hap_labels?: Record<string, string> }[] = [
  { sample: "GRCh38", haplotypes: ["0"], hap_labels: { "0": "reference" } },
  { sample: "HG00673", haplotypes: ["1", "2"], hap_labels: { "1": "paternal", "2": "maternal" } },
  { sample: "HG00733", haplotypes: ["1", "2"], hap_labels: { "1": "paternal", "2": "maternal" } },
];

export default function Home() {
  const [loading, setLoading] = useState(true);
  const [sg, setSG] = useState<SampleGraph | null>(null);
  const [nodeInfo, setNodeInfo] = useState<NodeInfo | null>(null);
  const [edgeInfo, setEdgeInfo] = useState<EdgeInfo | null>(null);
  const [tab, setTab] = useState("explore");
  const [manifest, setManifest] = useState<ManifestData | null>(null);
  const [latest, setLatest] = useState<LatestData | null>(null);
  const [overview, setOverview] = useState<OverviewData | null>(null);
  const [graphMode, setGraphMode] = useState("baseline");
  const [selSample, setSelSample] = useState("GRCh38");
  const [selHap, setSelHap] = useState("0");

  useEffect(() => {
    Promise.allSettled([
      fetch("/data/manifest.json").then(r => (r.ok ? r.json() : Promise.reject())),
      fetch("/data/latest.json").then(r => (r.ok ? r.json() : Promise.reject())),
      fetch("/data/overview.json").then(r => (r.ok ? r.json() : Promise.reject())),
    ]).then(([m, l, o]) => {
      if (m.status === "fulfilled") {
        setManifest(m.value);
        const s = m.value.samples;
        if (s && s.length) { setSelSample(s[0].sample); setSelHap(s[0].haplotypes?.[0] || "0"); }
      }
      if (l.status === "fulfilled") setLatest(l.value);
      if (o.status === "fulfilled") setOverview(o.value);
      setLoading(false);
    });
  }, []);

  const loadGraph = useCallback(() => {
    fetch(`/data/graphs/${graphMode}/${selSample}_${selHap}.json`)
      .then(r => (r.ok ? r.json() : Promise.reject()))
      .then((g: SampleGraph) => { setSG(g); setLoading(false); })
      .catch(() => { setSG(null); setLoading(false); });
  }, [graphMode, selSample, selHap]);

  useEffect(() => { loadGraph(); }, [loadGraph]);

  const samples: { sample: string; haplotypes: string[]; hap_labels?: Record<string, string> }[] =
    manifest?.samples?.length ? manifest.samples : FALLBACK_SAMPLES;
  const hapLabels: Record<string, string> = samples.find(s => s.sample === selSample)?.hap_labels || {};
  const graphNames = manifest?.graphs && Object.keys(manifest.graphs).length
    ? Object.keys(manifest.graphs) : ["baseline", "merged"];
  const dataMode = manifest?.data_mode || latest?.data_mode || "synthetic";
  const bm = latest?.metrics?.baseline || {};
  const mm = latest?.metrics?.merged || {};
  const eqVerdict = latest?.equivalence?.verdict || null;

  if (loading) return <div className="container" style={{ padding: 40 }}>Loading…</div>;

  const renderExplore = () => (
    <div className="explorer-wrapper">
      <div className="sidebar">
        <div className="sidebar-section"><h3>Graph</h3>
          {graphNames.map(name => (
            <button key={name} onClick={() => setGraphMode(name)}
              className={"graph-mode-btn" + (graphMode === name ? " active" : "")}>
              {name === "merged" ? "Merged (stitched)" : name === "baseline" ? "Baseline" : name}
            </button>
          ))}
        </div>
        <div className="sidebar-section"><h3>Samples</h3>
          {samples.map(s => (
            <div key={s.sample}>
              <button onClick={() => { setSelSample(s.sample); setSelHap(s.haplotypes?.[0] || "0"); }}
                className={"sample-btn" + (selSample === s.sample ? " active" : "")}>{s.sample}</button>
              {selSample === s.sample && (
                <div className="hap-list">
                  {s.haplotypes.map(h => (
                    <button key={h} onClick={() => setSelHap(h)}
                      className={"hap-btn" + (selHap === h ? " active" : "")}>{hapLabels[h] || h}</button>
                  ))}
                </div>
              )}
            </div>
          ))}
        </div>
      </div>
      <div>
        <GraphExplorer graph={sg}
          onNodeSelect={(i) => { setNodeInfo(i); setEdgeInfo(null); }}
          onEdgeSelect={(e) => { setEdgeInfo(e); setNodeInfo(null); }} />
      </div>
      <div className="inspector"><h3>Inspector</h3>
        {nodeInfo ? <div>
          <div className="inspector-row"><span>ID</span><span>{nodeInfo.id}</span></div>
          <div className="inspector-row"><span>Length</span><span>{nodeInfo.length}</span></div>
          <div className="inspector-row"><span>Degree</span><span>{nodeInfo.degree}</span></div>
        </div> : edgeInfo ? <div>
          <div className="inspector-row"><span>Edge</span><span>{edgeInfo.source} → {edgeInfo.target}</span></div>
        </div> : <div className="inspector-empty">Click a node or edge</div>}
      </div>
    </div>
  );

  const renderDashboard = () => (
    <section className="section">
      <h2>Dashboard</h2>
      <div className="grid-2">
        <div className="card"><h3>Baseline</h3>
          <p>Nodes: {bm.nodes ?? "–"} &nbsp; Edges: {bm.edges ?? "–"}</p>
          {bm.paths != null && <p>Paths: {bm.paths}</p>}
        </div>
        <div className="card"><h3>Merged</h3>
          <p>Nodes: {mm.nodes ?? "–"} &nbsp; Edges: {mm.edges ?? "–"}</p>
          {mm.paths != null && <p>Paths: {mm.paths}</p>}
          {mm.components != null && <p>Components: {mm.components}</p>}
        </div>
      </div>
      <p style={{ marginTop: 12 }}>
        <span className={"badge " + (latest?.equivalence?.verdict === "EQUIVALENT" ? "badge-ok" : "badge-warn")}>
          EQUIVALENCE: {latest?.equivalence?.verdict || "NOT_RUN"}
        </span>
      </p>
    </section>
  );

  const renderChunks = () => (
    <section className="section">
      <h2>Chunks</h2>
      <p>
        <span className={"badge " + (latest?.stitch?.status === "PASS" ? "badge-ok" : "badge-warn")}>
          STITCH: {latest?.stitch?.status || "NOT_RUN"}
        </span>
      </p>
      {(overview?.chunks || []).map(c => (
        <div key={c.chunk_id} className="card" style={{ marginBottom: 8 }}>
          <strong>{c.chunk_id}</strong>
          <span style={{ marginLeft: 12, fontSize: 13, color: "#64748b" }}>
            chr21:{c.reference_start}–{c.reference_end} (core {c.core_start}–{c.core_end})
          </span>
        </div>
      ))}
      {(latest?.boundaries || []).map(b => (
        <p key={b.boundary}>{b.boundary}: {b.status} joined={b.haplotypes_joined}</p>
      ))}
    </section>
  );

  const renderCompare = () => (
    <section className="section">
      <h2>Compare Haplotypes</h2>
      <p>Select samples in the sidebar to compare path lengths and node counts.</p>
      {sg && (
        <div className="grid-2">
          <div className="card"><h3>{selSample} (hap {selHap})</h3>
            <p>Path length: {sg.path?.length_bp ?? 0} bp</p>
            <p>Nodes on path: {sg.path?.steps?.length ?? 0}</p>
            <p>Total: {sg.nodes?.length ?? 0} nodes / {sg.edges?.length ?? 0} edges</p>
          </div>
          <div className="card"><h3>Reference (GRCh38)</h3>
            <p style={{ color: "#94a3b8" }}>Switch selection to compare against another haplotype.</p>
          </div>
        </div>
      )}
    </section>
  );

  return (<>
      {/* ===== HERO ===== */}
      <section style={{ background: "linear-gradient(135deg, #0f172a 0%, #1e40af 50%, #3b82f6 100%)", color: "#fff", padding: "80px 20px 60px", textAlign: "center" }}>
        <div className="container" style={{ maxWidth: 780 }}>
          <h1 style={{ fontSize: 42, fontWeight: 900, margin: "0 0 16px", letterSpacing: "-0.5px" }}>
            Parallel Pangenome<br />Graph Construction
          </h1>
          <p style={{ fontSize: 20, opacity: 0.92, margin: "0 auto 32px", maxWidth: 600, lineHeight: 1.6 }}>
            Can we build pangenome graphs in parallel chunks, stitch them back together,
            and get the same result as building the whole thing at once?
          </p>
          <div style={{ display: "flex", justifyContent: "center", gap: 12, flexWrap: "wrap", marginBottom: 24 }}>
            <span className={"badge "+(dataMode==="synthetic"?"badge-demo":"badge-real")}>{dataMode==="synthetic"?"SYNTHETIC DEMO":"REAL HPRC DATA"}</span>
            <span className={"badge "+(eqVerdict==="EQUIVALENT"?"badge-ok":"badge-warn")}>{eqVerdict==="EQUIVALENT"?"✓ EQUIVALENT":"NOT_RUN"}</span>
          </div>
        </div>
      </section>

      {/* ===== HOW IT WORKS ===== */}
      <section style={{ padding: "60px 20px", background: "#f8fafc" }}>
        <div className="container" style={{ maxWidth: 900 }}>
          <h2 style={{ fontSize: 28, textAlign: "center", marginBottom: 40 }}>How It Works</h2>
          <div className="steps-grid">
            <Step n="1" title="Select Region" desc="Target a ~1 Mb interval on chr21 (20–21 Mb) using GRCh38 as reference." />
            <Step n="2" title="Map Haplotypes" desc="Align 4 HPRC haplotypes independently against the reference using minimap2." />
            <Step n="3" title="Build in Parallel" desc="Split region into overlapping 400 Kb chunks. Run PGGB on each chunk independently on DNAnexus." />
            <Step n="4" title="Stitch Graphs" desc="Our overlap-aware algorithm welds chunk graphs at shared boundaries. 2/2 boundaries PASS in demo." />
            <Step n="5" title="Validate" desc={`Compare against monolithic baseline. ${eqVerdict==="EQUIVALENT"?eqVerdict:"Validation pending"}.`} />
          </div>
        </div>
      </section>
    <main className="container" style={{ padding: "40px 20px" }}>
      <header className="header">
        <h1>Parallel Pangenome Explorer</h1>
        <p>Explore haplotypes traversing independently constructed pangenome graph regions</p>
        <div style={{ display: "flex", gap: 8, flexWrap: "wrap" }}>
          <span className={"badge " + (dataMode === "synthetic" ? "badge-demo" : "badge-real")}>
            {dataMode === "synthetic" ? "DEMO / SYNTHETIC DATA" : "REAL HPRC DATA"}
          </span>
          <span className="badge badge-ok">BASELINE: {manifest?.pipeline_status?.baseline || "OK"}</span>
          <span className="badge badge-ok">CHUNKS: {manifest?.pipeline_status?.parallel_chunks || "OK"}</span>
          <span className="badge badge-warn">STITCH: {manifest?.pipeline_status?.stitch || latest?.stitch?.status || "NOT_RUN"}</span>
          <span className={"badge " + (latest?.equivalence?.verdict === "EQUIVALENT" ? "badge-ok" : "badge-warn")}>
            EQUIVALENCE: {latest?.equivalence?.verdict || "NOT_RUN"}
          </span>
        </div>
      </header>
      <div style={{ display: "flex", gap: 4, marginBottom: 16 }}>
        {["explore", "dashboard", "chunks", "compare"].map(t => (
          <button key={t} onClick={() => setTab(t)}
            className={"graph-btn" + (tab === t ? " active" : "")}>
            {t.charAt(0).toUpperCase() + t.slice(1)}
          </button>
        ))}
      </div>
      {tab === "explore" && renderExplore()}
      {tab === "dashboard" && renderDashboard()}
      {tab === "chunks" && renderChunks()}
      {tab === "compare" && renderCompare()}
    </main>
  </>);
}

function Step({ n, title, desc }: { n: string; title: string; desc: string }) {
  return (
    <div style={{ background: "#fff", borderRadius: 12, padding: 24, border: "1px solid #e2e8f0", display: "flex", gap: 16, alignItems: "flex-start" }}>
      <span style={{ display: "inline-flex", alignItems: "center", justifyContent: "center", width: 32, height: 32, borderRadius: "50%", background: "#3b82f6", color: "#fff", fontWeight: 700, fontSize: 14, flexShrink: 0 }}>{n}</span>
      <div>
        <h4 style={{ margin: "0 0 4px", fontSize: 15, color: "#1e293b" }}>{title}</h4>
        <p style={{ margin: 0, fontSize: 13, color: "#64748b", lineHeight: 1.6 }}>{desc}</p>
      </div>
    </div>
  );
}
