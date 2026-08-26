# Linear / Assembly-Based SV Pipeline

This is a **secondary validation workflow**. It calls structural variants
by aligning de novo assemblies against GRCh38 using `dipcall`.

## Purpose

Validates that variants discovered by the pangenome graph method
are biologically real (i.e., also detectable through independent
assembly alignment).

## Tools

- **dipcall** — primary assembly-based SV caller
- **SVIM-asm** — alternative caller (future)

## Important

This workflow is **separate** from the PGGB graph pipeline.
It is NOT the monolithic PGGB baseline.

Our core comparison is:
```
Monolithic PGGB (baseline)  vs.  Parallel PGGB + Reassembly (merged)
```

The linear assembly workflow provides orthogonal validation.