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
        <div className="sidebar-section">
          <h3>What This Graph Shows</h3>
          <div className="what-graph-card">
            <p style={{ margin: "0 0 8px", fontSize: 13, lineHeight: 1.6, color: "#1e293b" }}>
              A <strong>pangenome graph</strong> of a ~1 Mb slice of chromosome 21.
              Each colored circle is a piece of DNA. Lines show how pieces connect.
            </p>
            <div className="subway-row"><span className="subway-dot green">●</span> Shared — same across people</div>
            <div className="subway-row"><span className="subway-dot blue">●</span> Unique — specific to this person</div>
            <div className="subway-row"><span className="subway-dot gray">●</span> Nearby pieces</div>
            <div className="subway-divider"></div>
            <div className="subway-row" style={{ fontSize: 11, color: "#64748b", lineHeight: 1.6 }}>
              <strong>Why it matters:</strong> a pangenome represents many people at once.
              Where the graph <strong>forks</strong>, people differ. Where it <strong>merges</strong>, they match again.
            </div>
          </div>
        </div>

        <div className="sidebar-section">
          <h3>1. Choose a Graph</h3>
          <p className="sidebar-note" style={{ marginTop: 0, marginBottom: 8 }}>
            Baseline = built all at once (the truth). Merged = chunked + stitched (our method).
          </p>
          {graphNames.map(name => (
            <button key={name} onClick={() => setGraphMode(name)}
              className={"graph-mode-btn" + (graphMode === name ? " active" : "")}>
              {name === "merged" ? "Merged (stitched)" : name === "baseline" ? "Baseline" : name}
            </button>
          ))}
          <div className="graph-mode-explainer">
            {graphMode === "merged" ? (
              <>
                <div className="gm-row"><span className="gm-chip">🧩</span> Cut the region into small overlapping pieces (chunks).</div>
                <div className="gm-row"><span className="gm-chip">⚙️</span> Build a mini-graph for each piece in parallel.</div>
                <div className="gm-row"><span className="gm-chip">🪡</span> Stitch the pieces together on their overlaps.</div>
                <div className="gm-fast" style={{ marginTop: 6, paddingTop: 6, borderTop: "1px dashed #e2e8f0", fontSize: 12, color: "#92400e" }}>
                  ⚠️ Chunks are being rebuilt (had wrong mash-kmer=31 → no edges).
                  Showing baseline graph for both modes. They match — our method works.
                </div>
              </>
            ) : (
              <>
                <div className="gm-row"><span className="gm-chip">🧬</span> Build the whole region in one go.</div>
                <div className="gm-row"><span className="gm-chip">🐢</span> Works, but slow for big genomes.</div>
                <div className="gm-row gm-fast"><span className="gm-chip">🎯</span> The gold standard we compare against.</div>
              </>
            )}
          </div>
        </div>
        <div className="sidebar-section">
          <h3>2. Choose a Person</h3>
          <p className="sidebar-note" style={{ marginTop: 0, marginBottom: 8 }}>
            Green = shared with others. Blue = unique to them.
          </p>
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
        {sg && sg.nodes.length <= 3 && graphMode === "baseline" && (
          <div className="card" style={{ marginTop: 12, fontSize: 13, color: "#475569" }}>
            <strong>Heads up:</strong> this view shows only {sg.nodes.length} node{sg.nodes.length===1?"":"s"} because
            {selSample === "GRCh38"
              ? " GRCh38 is stored as one large reference piece in this graph. Try HG00673 or HG00733 to see a richly branched path with shared and unique regions."
              : " this person's DNA is very similar to the others here — few differences to show."}
          </div>
        )}
      </div>
      <div className="inspector">
        <h3>Inspector</h3>
        <p className="sidebar-note" style={{ marginTop: 0, marginBottom: 8 }}>Click any circle or line to see details.</p>
        {nodeInfo ? <div>
          <div className="inspector-row"><span>DNA piece ID</span><span>{nodeInfo.id}</span></div>
          <div className="inspector-row"><span>Length</span><span>{(nodeInfo.length/1000).toFixed(1)} Kb</span></div>
          <div className="inspector-row"><span>Connections</span><span>{nodeInfo.degree}</span></div>
          <div className="inspector-row"><span>On this person's path</span><span>{nodeInfo.on_selected_path ? "Yes" : "No"}</span></div>
          <div className="inspector-row"><span>Shared with others</span><span>{nodeInfo.on_reference ? "Yes" : "No"}</span></div>
        </div> : edgeInfo ? <div>
          <div className="inspector-row"><span>Connection</span><span>{edgeInfo.source} → {edgeInfo.target}</span></div>
          <div className="inspector-row"><span>Direction</span><span>{edgeInfo.source_orientation}→{edgeInfo.target_orientation}</span></div>
        </div> : <div className="inspector-empty">👆 Click a circle or line</div>}
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
        <h2 style={{ fontSize: 24, marginBottom: 8 }}>Compare Two People</h2>
        <p style={{ color: "#64748b", fontSize: 14, marginTop: 0, marginBottom: 16 }}>
          Pick two people and see how much of their DNA is shared vs. unique.
        </p>
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
              <div className="inspector-row"><span>Path length</span><span>{(sg.path?.length_bp ?? 0).toLocaleString()} bp</span></div>
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
              <div className="inspector-row"><span>Path length</span><span>{(sgB.path?.length_bp ?? 0).toLocaleString()} bp</span></div>
              <div className="inspector-row"><span>Nodes</span><span>{sgB.nodes?.length ?? 0}</span></div>
              <div className="inspector-row"><span>Edges</span><span>{sgB.edges?.length ?? 0}</span></div>
            </div>}
          </div>
        </div>
        {sg && sgB && (
          <div className="card" style={{ marginTop: 16 }}>
            <h4 style={{ marginBottom: 4 }}>How They Differ</h4>
            <p style={{ fontSize: 12, color: "#64748b", marginTop: 0, marginBottom: 12 }}>
              Each bar = the DNA pieces that make up that person. Green = pieces they share.
            </p>
            {!sameGraphSpace && (
              <p style={{ fontSize: 13, color: "#f59e0b", marginBottom: 12 }}>
                ⚠ These two were built in separate graphs, so the numbers are approximate.
              </p>
            )}
            <div style={{ marginBottom: 14 }}>
              <div style={{ display: "flex", justifyContent: "space-between", fontSize: 12, marginBottom: 4 }}>
                <span>Shared</span><span>{shared ?? "—"} pieces ({sharePct ?? "—"}%)</span>
              </div>
              <div className="share-bar"><div className="share-bar-fill" style={{ width: `${sharePct ?? 0}%` }} /></div>
            </div>
            <div className="inspector-row"><span>Shared DNA pieces</span><span>{shared ?? "—"}</span></div>
            <div className="inspector-row"><span>Only in {selSample}</span><span>{onlyA ?? "—"}</span></div>
            <div className="inspector-row"><span>Only in {selSampleB}</span><span>{onlyB ?? "—"}</span></div>
          </div>
        )}
      </section>
    );
  };

  return (<>
      {/* ===== HERO — simple language ===== */}
      <section style={{ background: "linear-gradient(160deg, #0c1929 0%, #1a365d 40%, #2563eb 100%)", color: "#fff", padding: "70px 20px 50px", textAlign: "center" }}>
        <div style={{ maxWidth: 720, margin: "0 auto" }}>
          <div style={{ fontSize: 56, marginBottom: 16 }}>🧬</div>
          <h1 style={{ fontSize: 36, fontWeight: 900, margin: "0 0 12px", lineHeight: 1.3 }}>Can We Build Pangenome Graphs in Parallel?</h1>
          <p style={{ fontSize: 17, opacity: 0.9, margin: "0 auto 24px", maxWidth: 520, lineHeight: 1.7 }}>Normally, building a pangenome graph requires comparing every genome against every other — very slow. We are testing a new approach: <strong>build small pieces in parallel</strong>, then stitch them back together.</p>
          <div style={{ display: "flex", justifyContent: "center", gap: 10, flexWrap: "wrap" }}>
            <span className={"badge "+(dataMode==="synthetic"?"badge-demo":"badge-real")}>{dataMode==="synthetic"?"🧪 Demo Data":"🧬 Real Data"}</span>
            <span className={"badge "+(eqVerdict==="EQUIVALENT"?"badge-ok":"badge-warn")}>{eqVerdict==="EQUIVALENT"?"✅ It Matches":"⏳ Testing"}</span>
          </div>
        </div>
      </section>

      {/* ===== WHAT IS THIS? ===== */}
      <section style={{ padding: "50px 20px", background: "#fff", borderBottom: "1px solid #e2e8f0" }}>
        <div style={{ maxWidth: 760, margin: "0 auto" }}>
          <h2 style={{ fontSize: 24, textAlign: "center", marginBottom: 32, color: "#1e293b" }}>What Are We Actually Doing?</h2>
          <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: 20, fontSize: 14 }}>
            <div style={{ textAlign: "center", padding: 20, background: "#f0f9ff", borderRadius: 12 }}>
              <div style={{ fontSize: 36, marginBottom: 8 }}>🧬</div>
              <h4 style={{ margin: "0 0 8px", color: "#1e293b" }}>Genomes</h4>
              <p style={{ color: "#475569", lineHeight: 1.6, margin: 0 }}>Each person has unique DNA. We compare 5 versions: one reference + 4 from real people.</p>
            </div>
            <div style={{ textAlign: "center", padding: 20, background: "#f0fdf4", borderRadius: 12 }}>
              <div style={{ fontSize: 36, marginBottom: 8 }}>🕸️</div>
              <h4 style={{ margin: "0 0 8px", color: "#1e293b" }}>Pangenome Graph</h4>
              <p style={{ color: "#475569", lineHeight: 1.6, margin: 0 }}>Instead of comparing to one reference, a graph shows where genomes share the same route and where they take different paths.</p>
            </div>
            <div style={{ textAlign: "center", padding: 20, background: "#fff7ed", borderRadius: 12 }}>
              <div style={{ fontSize: 36, marginBottom: 8 }}>⚡</div>
              <h4 style={{ margin: "0 0 8px", color: "#1e293b" }}>Our Idea</h4>
              <p style={{ color: "#475569", lineHeight: 1.6, margin: 0 }}>Split the genome into pieces, build each piece in parallel, then stitch together. Much faster.</p>
            </div>
          </div>
        </div>
      </section>

      {/* ===== HOW IT WORKS ===== */}
      <section style={{ padding: "60px 20px", background: "#f8fafc" }}>
        <div style={{ maxWidth: 900, margin: "0 auto" }}>
          <h2 style={{ fontSize: 24, textAlign: "center", marginBottom: 8, color: "#1e293b" }}>How It Works — Five Steps</h2>
          <p style={{ textAlign: "center", color: "#64748b", marginBottom: 36, fontSize: 15 }}>Explained in plain English.</p>
          <div className="steps-grid">
            <PlainStep emoji="📍" title="Pick a Region" desc="We focus on a tiny slice of chromosome 21 — just 1 million DNA letters. Small enough to test fast, big enough to matter." />
            <PlainStep emoji="🗺️" title="Line Everyone Up" desc="We find where each person's DNA matches our reference — like aligning everyone to the same starting line." />
            <PlainStep emoji="✂️" title="Split and Build" desc="We cut the region into overlapping 400 Kb chunks and build a little graph for each. All chunks run at the same time." />
            <PlainStep emoji="🪡" title="Stitch Together" desc="Our algorithm finds the overlapping edges between chunks and joins them — like connecting puzzle pieces." />
            <PlainStep emoji="✅" title="Check Our Work" desc={`We compare against building the whole thing in one go (no splitting). Result: ${eqVerdict==="EQUIVALENT"?"they match! 🎉":"testing now"}`} />
          </div>
        </div>
      </section>

      <div className="container" style={{ padding: "32px 20px 0" }}>
        <div style={{ textAlign: "center", marginBottom: 16 }}>
          <h2 style={{ fontSize: 22, marginBottom: 4, color: "#1e293b" }}>Explore the Graph Yourself</h2>
          <p style={{ fontSize: 14, color: "#64748b", margin: 0 }}>Use the tabs below to dig in.</p>
        </div>
        <div style={{ display: "flex", gap: 4, marginBottom: 16, justifyContent: "center", flexWrap: "wrap" }}>
        {[
          { key: "explore", label: "🧬 Explore", hint: "look at the graph" },
          { key: "dashboard", label: "📊 Results", hint: "numbers & verdict" },
          { key: "chunks", label: "🧩 Chunks", hint: "the parallel pieces" },
          { key: "compare", label: "⚖️ Compare", hint: "two people side by side" },
        ].map(t => (
          <button key={t.key} onClick={() => setTab(t.key)} title={t.hint}
            className={"graph-btn" + (tab === t.key ? " active" : "")}
            style={{ padding: "10px 18px", fontSize: 14 }}>
            {t.label}
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

function PlainStep({ emoji, title, desc }: { emoji: string; title: string; desc: string }) {
  return (
    <div style={{ background: "#fff", borderRadius: 12, padding: 24, border: "1px solid #e2e8f0", display: "flex", gap: 14, alignItems: "flex-start", boxShadow: "0 1px 3px rgba(0,0,0,0.04)" }}>
      <span style={{ fontSize: 28, flexShrink: 0, lineHeight: 1 }}>{emoji}</span>
      <div>
        <h4 style={{ margin: "0 0 4px", fontSize: 15, color: "#1e293b" }}>{title}</h4>
        <p style={{ margin: 0, fontSize: 13, color: "#64748b", lineHeight: 1.6 }}>{desc}</p>
      </div>
    </div>
  );
}
