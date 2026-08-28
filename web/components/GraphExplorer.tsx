"use client";
import { useEffect, useRef, useState, useCallback } from "react";
import cytoscape, { Core, EventObject } from "cytoscape";
import type { SampleGraph, GraphNode, GraphEdge, NodeInfo, EdgeInfo } from "../types";

interface Props {
  graph: SampleGraph | null;
  onNodeSelect?: (info: NodeInfo | null) => void;
  onEdgeSelect?: (info: EdgeInfo | null) => void;
}

const MIN_NODE = 14;
const MAX_NODE = 46;

function nodeSize(length: number) {
  // log scale so a 1 Mb node doesn't dwarf a 1 kb node, but sizes stay readable
  if (!length) return MIN_NODE;
  const s = MIN_NODE + Math.log10(length + 1) * 4.2;
  return Math.min(MAX_NODE, Math.max(MIN_NODE, s));
}

export default function GraphExplorer({
  graph, onNodeSelect, onEdgeSelect,
}: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const cyRef = useRef<Core | null>(null);
  const [searchVal, setSearchVal] = useState("");
  const [tooltip, setTooltip] = useState<{ x: number; y: number; text: string } | null>(null);

  const buildElements = useCallback((g: SampleGraph) => {
    const els: any[] = g.nodes.map((n: GraphNode) => ({
      group: "nodes",
      data: { id: n.id, label: n.id, length: n.length,
        onPath: n.on_selected_path, onRef: n.on_reference, degree: n.degree },
      classes: n.on_reference ? "onref" : n.on_selected_path ? "onpath" : "",
    }));
    els.push(...g.edges.map((e: GraphEdge) => ({
      group: "edges",
      data: { id: e.source + "-" + e.target, source: e.source, target: e.target,
        sourceOrientation: e.source_orientation, targetOrientation: e.target_orientation,
        onPath: e.on_selected_path || false },
      classes: e.on_selected_path ? "onpath" : "",
    })));
    return els;
  }, []);

  useEffect(() => {
    if (!containerRef.current) return;
    if (cyRef.current) { cyRef.current.destroy(); cyRef.current = null; }
    if (!graph || graph.nodes.length === 0) return;

    const cy = cytoscape({
      container: containerRef.current,
      elements: buildElements(graph),
      style: [
        { selector: "node", style: {
            width: (el: any) => nodeSize(el.data("length")),
            height: (el: any) => nodeSize(el.data("length")),
            "background-color": "#cbd5e1",
            "border-width": 1.5, "border-color": "#64748b",
            // no permanent text labels — hover tooltip instead (keeps canvas readable)
            "text-opacity": 0,
        }},
        { selector: "node.onpath", style: {
            "background-color": "#3b82f6", "border-color": "#1d4ed8",
        }},
        { selector: "node.onref", style: {
            "background-color": "#22c55e", "border-color": "#15803d",
        }},
        { selector: "node:selected", style: {
            "border-color": "#f59e0b", "border-width": 4,
        }},
        { selector: "node.highlighted", style: {
            "border-color": "#f59e0b", "border-width": 4,
        }},
        { selector: "edge", style: {
            width: 1.5, "line-color": "#cbd5e1", "target-arrow-color": "#cbd5e1",
            "target-arrow-shape": "triangle", "arrow-scale": 0.8,
            "curve-style": "bezier",
        }},
        { selector: "edge.onpath", style: {
            width: 3, "line-color": "#3b82f6", "target-arrow-color": "#3b82f6",
        }},
      ],
      // "subway map" layout: reference backbone flows one direction,
      // variants split off and rejoin — far more intuitive than a force-directed ball.
      layout: {
        name: "breadthfirst",
        directed: true,
        spacingFactor: 1.15,
        avoidOverlap: true,
        animate: false,
      },
      wheelSensitivity: 0.3, minZoom: 0.05, maxZoom: 10,
    });
    cyRef.current = cy;

    cy.on("tap", "node", (evt: EventObject) => {
      const n = evt.target; const d = n.data();
      if (onNodeSelect) onNodeSelect({
        id: d.id, length: d.length || 0, degree: d.degree || 0,
        on_selected_path: !!d.onPath, on_reference: !!d.onRef, neighbors: [],
      });
      cy.elements().removeClass("highlighted");
      n.addClass("highlighted");
      n.neighborhood().addClass("highlighted");
    });

    cy.on("tap", "edge", (evt: EventObject) => {
      const e = evt.target; const d = e.data();
      if (onEdgeSelect) onEdgeSelect({
        source: d.source, target: d.target,
        source_orientation: d.sourceOrientation || "+",
        target_orientation: d.targetOrientation || "+",
        on_selected_path: !!d.onPath,
      });
    });

    cy.on("tap", (evt: EventObject) => {
      if (evt.target === cy) {
        cy.elements().removeClass("highlighted");
        if (onNodeSelect) onNodeSelect(null);
        if (onEdgeSelect) onEdgeSelect(null);
      }
    });

    cy.on("mouseover", "node", (evt: EventObject) => {
      const n = evt.target; const d = n.data();
      const pos = n.renderedPosition();
      const role = d.onRef && d.onPath ? "shared backbone + this person" :
                   d.onRef ? "shared with others" :
                   d.onPath ? "unique to this person" : "other";
      const kb = ((d.length || 0) / 1000).toFixed(1);
      setTooltip({
        x: pos.x, y: pos.y,
        text: `${d.id}\n${kb} Kb · ${role}`,
      });
    });
    cy.on("mouseout", "node", () => setTooltip(null));

    return () => { cy.destroy(); cyRef.current = null; };
  }, [graph, onNodeSelect, onEdgeSelect, buildElements]);

  const fitGraph = () => cyRef.current?.fit(undefined, 50);
  const resetGraph = () => {
    if (cyRef.current) { cyRef.current.zoom(1); cyRef.current.center(); }
  };

  const searchNode = () => {
    if (!cyRef.current || !searchVal.trim()) return;
    cyRef.current.elements().removeClass("highlighted");
    const node = cyRef.current.getElementById(searchVal.trim());
    if (node && node.length > 0) {
      node.addClass("highlighted");
      cyRef.current.animate({ fit: { eles: node, padding: 100 }, duration: 300 });
    }
  };

  if (!graph || graph.nodes.length === 0) {
    return <div className="empty-graph">Select a sample to view graph</div>;
  }

  return (
    <div>
      <div className="graph-controls">
        <button onClick={fitGraph} className="graph-btn">Fit</button>
        <button onClick={resetGraph} className="graph-btn">Reset</button>
        <input type="text" placeholder="Search node ID..."
          value={searchVal} onChange={e => setSearchVal(e.target.value)}
          onKeyDown={e => e.key === "Enter" && searchNode()}
          className="graph-search-input" />
        <button onClick={searchNode} className="graph-btn">Go</button>
        {graph.truncated && (
          <span className="graph-truncated">
            Truncated: {graph.original_counts?.nodes || "?"} nodes
          </span>
        )}
      </div>
      <div style={{ position: "relative" }}>
        <div ref={containerRef} className="graph-canvas" />
        {tooltip && (
          <div className="graph-tooltip" style={{ left: tooltip.x + 12, top: tooltip.y + 12 }}>
            {tooltip.text.split("\n").map((l, i) => <div key={i}>{l}</div>)}
          </div>
        )}
      </div>
    </div>
  );
}
