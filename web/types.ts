// Shared types for pangenome graph explorer

export interface Manifest {
  schema_version: string;
  run_id: string;
  data_mode: "synthetic" | "real" | "fixture";
  target: { reference: string; chromosome: string; start: number; end: number };
  pipeline_status: Record<string, string>;
  samples: { sample: string; haplotypes: string[]; hap_labels?: Record<string, string> }[];
  graphs: Record<string, { nodes: number; edges: number; paths: number; walks: number }>;
}

export interface GraphNode {
  id: string; length: number;
  on_selected_path: boolean; on_reference: boolean;
  degree: number; neighbors: string[];
}

export interface GraphEdge {
  source: string; target: string;
  source_orientation: string; target_orientation: string;
  on_selected_path?: boolean; on_reference?: boolean;
}

export interface SampleGraph {
  schema_version: string;
  sample: string; haplotype: string; path_name: string;
  nodes: GraphNode[]; edges: GraphEdge[];
  path: { steps: { node: string; orientation: string }[]; length_bp: number };
  metrics: Record<string, number>;
  truncated: boolean;
  original_counts?: { nodes: number; edges: number; path_steps: number };
}

export interface NodeInfo {
  id: string; length: number; degree: number;
  on_selected_path: boolean; on_reference: boolean;
  neighbors: { id: string; orientation: string }[];
}

export interface EdgeInfo {
  source: string; target: string;
  source_orientation: string; target_orientation: string;
  on_selected_path?: boolean;
}

export interface ChunkTiming {
  chunk_id: string; job_id: string; instance_type: string;
  started_running_ms: number; stopped_running_ms: number;
  wall_seconds: number; state: string;
}
