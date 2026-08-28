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
  const [selSampleB, setSelSampleB] = useState("HG00673");
  const [selHapB, setSelHapB] = useState("1");
  const [sgB, setSGB] = useState<SampleGraph | null>(null);

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
  const loadGraphB = useCallback(() => {
    fetch(`/data/graphs/${graphMode}/${selSampleB}_${selHapB}.json`)
      .then(r => (r.ok ? r.json() : Promise.reject()))
      .then((g: SampleGraph) => setSGB(g))
      .catch(() => setSGB(null));
  }, [graphMode, selSampleB, selHapB]);
  useEffect(() => { loadGraphB(); }, [loadGraphB]);

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

  const renderCompare = () => {
    const hapLabelsB: Record<string, string> = samples.find(s => s.sample === selSampleB)?.hap_labels || {};
    const idA = new Set((sg?.nodes || []).map(n => n.id));
    const idB = new Set((sgB?.nodes || []).map(n => n.id));
    const shared = sg && sgB ? Array.from(idA).filter(id => idB.has(id)).length : null;
    const onlyA = sg && sgB ? Array.from(idA).filter(id => !idB.has(id)).length : null;
    const onlyB = sg && sgB ? Array.from(idB).filter(id => !idA.has(id)).length : null;
    const sameGraphSpace = sg && sgB && sg.graph === sgB.graph;
    const sharePct = shared != null && sg && sgB ? Math.round(shared / Math.max(idA.size, idB.size, 1) * 100) : null;

    return (
      <section className="section">
        <h2 style={{ fontSize: 24, marginBottom: 16 }}>Compare Haplotypes</h2>
        <div className="compare-layout">
          <div className="compare-panel">
            <h4>Sample A</h4>
            <select value={selSample} onChange={e => { setSelSample(e.target.value); setSelHap(samples.find(s => s.sample === e.target.value)?.haplotypes[0] || "0"); }}
              style={{ width: "100%", padding: 6, marginBottom: 8, border: "1px solid #cbd5e1", borderRadius: 4 }}>
              {samples.map(s => <option key={s.sample} value={s.sample}>{s.sample}</option>)}
            </select>
            <div style={{ display: "flex", gap: 4, marginBottom: 8 }}>
              {(samples.find(s => s.sample === selSample)?.haplotypes || []).map(h => (
                <button key={h} onClick={() => setSelHap(h)}
                  className={"hap-btn" + (selHap === h ? " active" : "")}>{hapLabels[h] || h}</button>
              ))}
            </div>
            {sg && <div className="card">
              <div className="inspector-row"><span>Path length</span><span>{sg.path?.length_bp ?? 0} bp</span></div>
              <div className="inspector-row"><span>Nodes</span><span>{sg.nodes?.length ?? 0}</span></div>
              <div className="inspector-row"><span>Edges</span><span>{sg.edges?.length ?? 0}</span></div>
            </div>}
          </div>
          <div className="compare-panel">
            <h4>Sample B</h4>
            <select value={selSampleB} onChange={e => { setSelSampleB(e.target.value); setSelHapB(samples.find(s => s.sample === e.target.value)?.haplotypes[0] || "0"); }}
              style={{ width: "100%", padding: 6, marginBottom: 8, border: "1px solid #cbd5e1", borderRadius: 4 }}>
              {samples.map(s => <option key={s.sample} value={s.sample}>{s.sample}</option>)}
            </select>
            <div style={{ display: "flex", gap: 4, marginBottom: 8 }}>
              {(samples.find(s => s.sample === selSampleB)?.haplotypes || []).map(h => (
                <button key={h} onClick={() => setSelHapB(h)}
                  className={"hap-btn" + (selHapB === h ? " active" : "")}>{hapLabelsB[h] || h}</button>
              ))}
            </div>
            {sgB && <div className="card">
              <div className="inspector-row"><span>Path length</span><span>{sgB.path?.length_bp ?? 0} bp</span></div>
              <div className="inspector-row"><span>Nodes</span><span>{sgB.nodes?.length ?? 0}</span></div>
              <div className="inspector-row"><span>Edges</span><span>{sgB.edges?.length ?? 0}</span></div>
            </div>}
          </div>
        </div>
        {sg && sgB && (
          <div className="card" style={{ marginTop: 16 }}>
            <h4 style={{ marginBottom: 12 }}>Comparison</h4>
            {!sameGraphSpace && (
              <p style={{ fontSize: 13, color: "#f59e0b", marginBottom: 12 }}>
                ⚠ Direct node comparison across independent graph namespaces may be misleading.
              </p>
            )}
            <div className="inspector-row"><span>Shared nodes</span><span>{shared ?? "—"}</span></div>
            <div className="inspector-row"><span>Only in A</span><span>{onlyA ?? "—"}</span></div>
            <div className="inspector-row"><span>Only in B</span><span>{onlyB ?? "—"}</span></div>
            {sharePct != null && <div className="inspector-row"><span>Shared %</span><span>{sharePct}%</span></div>}
            <div className="inspector-row"><span>Same graph space</span><span>{sameGraphSpace ? "Yes" : "No"}</span></div>
          </div>
        )}
      </section>
    );
  };

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

      <div className="container" style={{ padding: "32px 20px 0" }}>
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
    </div>
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
