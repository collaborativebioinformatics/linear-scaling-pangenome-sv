"use client"
import { useState, useEffect } from "react"

interface Run { run_id: string; pipeline_version: string; mode: string; data_mode: string; timestamp?: string }
interface Data { data_mode: string; run?: Run; target?: Record<string,unknown>; samples?: string[]; metrics?: Record<string,Record<string,number>>; boundaries?: unknown[]; bubbles?: unknown[]; graphWindow?: Record<string,unknown> }

function Metric({ label, value }: { label: string; value?: number|null }) {
  return <div className="metric-row"><span>{label}</span><span style={{fontWeight:600}}>{value ?? "N/A"}</span></div>
}

function Card({ title, rows }: { title: string; rows: [string,number|undefined|null][] }) {
  return <div className="card"><h3 style={{margin:"0 0 12px",fontSize:16}}>{title}</h3>{rows.map(([l,v]) => <Metric key={l} label={l} value={v}/>)}</div>
}

export default function Home() {
  const [data, setData] = useState<Data|null>(null)
  const [loading, setLoading] = useState(true)
  const [err, setErr] = useState<string|null>(null)

  useEffect(() => {
    const load = (u: string) => fetch(u).then(r => { if(!r.ok) throw Error(); return r.json() })
    load("/api/data").then(setData).catch(() => load("/data/latest.json").then(setData)).catch(e => setErr(e?.message ?? "fetch failed")).finally(() => setLoading(false))
  }, [])

  if (loading) return <div className="container" style={{padding:40}}>Loading...</div>
  if (err) return <div className="container" style={{padding:40}}>Error: {err}</div>
  if (!data) return <div className="container" style={{padding:40}}>No data</div>

  const isDemo = data.data_mode === "synthetic"
  const bl = data.metrics?.baseline ?? {}
  const mg = data.metrics?.merged ?? {}

  return (
    <main className="container" style={{padding:"40px 20px"}}>
      <header className="header">
        <h1>Parallel Pangenome Graph Explorer</h1>
        <p>Comparing monolithic vs parallel pangenome graph construction</p>
        <span className={"badge " + (isDemo ? "badge-demo" : "badge-real")}>{isDemo ? "DEMO DATA" : "REAL HPRC DATA"}</span>
      </header>
      <section className="section">
        <h2>Dashboard</h2>
        <div className="grid-2">
          <Card title="Monolithic (Baseline)" rows={[["Nodes",bl.nodes as number],["Edges",bl.edges as number],["Paths",bl.paths as number]]} />
          <Card title="Parallel + Merged" rows={[["Nodes",mg.nodes as number],["Edges",mg.edges as number],["Paths",mg.paths as number]]} />
        </div>
      </section>
      <section className="section">
        <h2>Target</h2>
        <div className="card">
          <p><strong>Reference:</strong> {data.target?.reference as string || "GRCh38"} &mdash; {data.target?.chromosome as string || "chr21"}</p>
          <p><strong>Samples:</strong> {data.samples?.join(", ") || "N/A"}</p>
        </div>
      </section>
    </main>
  )
}