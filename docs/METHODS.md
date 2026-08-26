# Methods

## Graph Construction

Both monolithic and parallel graphs use PGGB (PanGenome Graph Builder)
with identical parameters.

## Parallel Strategy

The reference chromosome is divided into overlapping windows.
Each window is constructed independently.

## Merge Algorithm

Adjacent graphs are stitched using the GRCh38 reference path
as the primary anchor to identify and deduplicate overlapping regions.

## Validation

- Path integrity: extract sequence spelled by each path and compare
- Graph topology: nodes, edges, components
- Variant calling: vg deconstruct + Truvari comparison
