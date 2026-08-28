"use client";
import { useState, useEffect } from "react";
import type { SampleGraph, NodeInfo, EdgeInfo } from "../types";
import GraphExplorer from "../components/GraphExplorer";

export default function Home() {
  const [loading, setLoading] = useState(true);
  const [sg, setSG] = useState(null as SampleGraph | null);
  const [nodeInfo, setNodeInfo] = useState(null as NodeInfo | null);
  const [edgeInfo, setEdgeInfo] = useState(null as EdgeInfo | null);
  const [tab, setTab] = useState("explore");
  const [selSample, setSelSample] = useState("GRCh38");
  const [selHap, setSelHap] = useState("0");
  const [meta, setMeta] = useState(null as any);

  const samples = [
    { name: "GRCh38", haps: ["0"] },
    { name: "HG00673", haps: ["1", "2"] },
    { name: "HG00733", haps: ["1", "2"] },
  ];

  useEffect(() => {
    fetch("/api/data").then(r => r.json()).then(setMeta).catch(() => {});
  }, []);

  useEffect(() => {
    setLoading(true);
    const n = 14 + Math.floor(Math.random() * 4);
    const g: SampleGraph = {
      schema_version: "1",
      sample: selSample, haplotype: selHap,
      path_name: selSample + "#" + selHap + "#chr21",
      nodes: [], edges: [],
      path: { steps: [], length_bp: 0 },
      metrics: {}, truncated: false,
    };
    for (let i = 1; i <= n; i++) {
      const onP = i <= 8;
      g.nodes.push({ id: "s"+i, length: 80+Math.floor(Math.random()*120), on_selected_path: onP, on_reference: selSample==="GRCh38"&&onP, degree: 0, neighbors: [] });
    }
    for (let i = 1; i < n; i++) {
      g.edges.push({ source: "s"+i, target: "s"+(i+1), source_orientation: "+", target_orientation: "+", on_selected_path: i<=7, on_reference: selSample==="GRCh38"&&i<=7 });
    }
    const deg: Record<string,number> = {};
    for (const e of g.edges) { deg[e.source] = (deg[e.source]||0) + 1; deg[e.target] = (deg[e.target]||0) + 1; }
    for (const node of g.nodes) { node.degree = deg[node.id]||0; node.neighbors = []; }
    g.path.steps = g.nodes.filter(x=>x.on_selected_path).map(x=>({node:x.id,orientation:"+"}));
    g.path.length_bp = g.path.steps.length * 100;
    setSG(g);
    setLoading(false);
  }, [selSample, selHap]);

  if (loading) return <div className="container" style={{padding:40}}>Loading...</div>;

  const tabContent = (() => {
    const bm = meta?.metrics?.baseline || {};
    const mm = meta?.metrics?.merged || {};
    if (tab === "dashboard") return <section className="section"><h2>Dashboard</h2><p>Baseline: {bm.nodes ?? "–"}n {bm.edges ?? "–"}e. Merged: {mm.nodes ?? "–"}n {mm.edges ?? "–"}e. {meta?.equivalence?.verdict || ""}</p></section>;
    if (tab === "chunks") return <section className="section"><h2>Chunks</h2><p><span className={"badge "+((meta?.stitch?.status==="PASS")?"badge-ok":"badge-warn")}>STITCH {meta?.stitch?.status || "NOT_RUN"}</span></p>
      {(meta?.boundaries||[]).map((b:any)=> <p key={b.boundary}>{b.boundary}: {b.status} joined={b.haplotypes_joined}</p>)}</section>;
    if (tab === "compare") return <section className="section"><h2>Compare</h2></section>;
    return (
      <div className="explorer-wrapper">
        <div className="sidebar">
          <div className="sidebar-section"><h3>Samples</h3>
            {samples.map(s =>
              <div key={s.name}>
                <button onClick={()=>{setSelSample(s.name);setSelHap(s.haps[0]||"0")}}
                  className={"sample-btn"+(selSample===s.name?" active":"")}>{s.name}</button>
                {selSample===s.name&&<div className="hap-list">{s.haps.map(h=>
                  <button key={h} onClick={()=>setSelHap(h)}
                    className={"hap-btn"+(selHap===h?" active":"")}>{h==="0"?"reference":h==="1"?"paternal":"maternal"}</button>
                )}</div>}
              </div>
            )}
          </div>
        </div>
        <div><GraphExplorer graph={sg}
          onNodeSelect={(info)=>{setNodeInfo(info);setEdgeInfo(null)}}
          onEdgeSelect={(info)=>{setEdgeInfo(info);setNodeInfo(null)}}/></div>
        <div className="inspector"><h3>Inspector</h3>
          {nodeInfo?<div><div>ID: {nodeInfo.id}</div><div>Len: {nodeInfo.length}</div><div>Deg: {nodeInfo.degree}</div></div>
          :edgeInfo?<div><div>{edgeInfo.source} to {edgeInfo.target}</div></div>
          :<div className="inspector-empty">Click a node</div>}
        </div>
      </div>
    );
  })();

  return (
    <main className="container" style={{padding:"40px 20px"}}>
      <header className="header">
        <h1>Parallel Pangenome Explorer</h1>
        <p>Explore haplotypes traversing pangenome graph regions</p>
        <div style={{display:"flex",gap:8,flexWrap:"wrap"}}>
          <span className="badge badge-demo">DEMO</span>
          <span className="badge badge-ok">BASELINE: OK</span>
          <span className="badge badge-ok">CHUNKS: OK</span>
          <span className={"badge "+(meta?.stitch?.status==="PASS"?"badge-ok":"badge-warn")}>STITCH: {meta?.stitch?.status || "NOT_IMPLEMENTED"}</span>
          <span className={"badge "+(meta?.equivalence?.verdict==="EQUIVALENT"?"badge-ok":"badge-warn")}>EQUIVALENCE: {meta?.equivalence?.verdict || "NOT_RUN"}</span>
        </div>
      </header>
      <div style={{display:"flex",gap:4,marginBottom:16}}>
        {["explore","dashboard","chunks","compare"].map(t =>
          <button key={t} onClick={()=>setTab(t)} className={"graph-btn"+(tab===t?" active":"")}>{t}</button>
        )}
      </div>
      {tabContent}
    </main>
  );
}