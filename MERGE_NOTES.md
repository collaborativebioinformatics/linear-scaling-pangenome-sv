# Merge Notes — read before integrating branches

## 1. `pipeline/merge/gfa.py` — take THIS version, not your branch's

Four of the five branches shipped a `Walk` (W-line) parser that does not
match the GFA 1.1 spec. It expects an eighth `step_count` column and
splits the walk field on commas:

```
W  Sample  Hap  Seq  Start  End  step_count  n1+,n2+     <-- WRONG
W  Sample  Hap  Seq  Start  End  >n1<n2                  <-- SPEC
```

Real PGGB, Minigraph-Cactus, and `vg` output uses the spec form. Any
branch that keeps its own `gfa.py` will crash with
`ValueError: invalid literal for int()` the moment it parses genuine
tool output.

Fingerprints of what was shipped (md5, first 12 chars):

| branch | gfa.py | W-line parser |
|---|---|---|
| ali-dnanexus-integration | `d77e60f5e031` | broken |
| lex-pggb-benchmarking | `d77e60f5e031` | broken |
| michael-graph-construction | `d77e60f5e031` | broken |
| quang-overlap-aware-merge | `1eddb404b78f` | broken, + 11 helper fns |
| main | `fc11937e549f` | correct, missing the 11 helpers |

Neither main nor Quang's version alone is safe: main parses W-lines
correctly but lacks `haplotype_key`, `parse_pansn`, `get_path_sequence`,
`get_walk_sequence`, `dangling_links`, `orphan_segments`, `path_length`,
`path_steps`, `segment_sequence`, `used_segments`, and `walk_steps` —
all of which `pipeline/merge/merge_graphs.py` on Quang's branch imports.
Dropping main's file onto Quang's branch fails with
`ImportError: cannot import name 'haplotype_key'`.

The `gfa.py` in this bundle is the union: main's spec-correct `Walk` plus
all 11 of Quang's helpers.

### Backward compatibility

`Walk.__init__` accepts both call shapes, so no existing code needs
editing:

```python
Walk(sample, hap, contig, start, end, path)              # spec
Walk(sample, hap, contig, start, end, step_count, path)  # legacy
```

The legacy `step_count` argument is accepted and ignored. `step_count` is
now a derived property equal to `len(path)`, so it can never drift out of
sync with the actual walk. `_walk_to_segments` also tolerates a legacy
comma-separated walk field, so old GFA files still load.

Verified: Quang's full 61-test suite passes unchanged against this file.

## 2. Two bugs worth fixing in the repo (not in this bundle)

**`scripts/setup_demo.py` ignores `config.merge.strategy`.** It imports
and calls `diagnostic_disjoint_union` directly (line ~117), so `make demo`
never exercises the overlap-aware stitch no matter what the config says.
That is why the README still reports the merge as NOT_IMPLEMENTED while
the algorithm is present and working.

**`_load_chunks` in `merge_graphs.py` hardcodes `work/chunks/{id}.gfa`**
while the demo writes to `work/demo/chunks/`. In demo mode it silently
finds zero chunks. Separately, `overlap_aware_stitch` raises
`ValueError: path ... has no subrange and no chunk manifest row` unless
`chunk_rows` is passed — it cannot be called with graphs alone.

## 3. Validated behaviour

With the merged `gfa.py` in place, running Quang's two merge strategies on
the demo data and grading them with `graph_stats.py`:

| strategy | nodes | paths | components | verdict |
|---|---|---|---|---|
| `diagnostic_disjoint_union` | 599 | 15 | 5 → 15 | DIVERGENT (2/7) |
| `overlap_aware_stitch` | 499 | 5 | 5 → 5 | EQUIVALENT (7/7) |

`component_count` is the metric that separates them, and it is the direct
quantitative answer to the research question.
