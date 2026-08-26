# lex_testing — PGGB cost benchmarking

Standalone benchmarking of PGGB graph construction. Independent of
`pipeline/` — these scripts call PGGB directly so the actual command and
its parameters stay visible.

## Ground rules

- **DNAnexus is read-only.** The only `dx` call anywhere here is
  `dx download` in `fetch_inputs.sh`. Nothing uploads, creates folders,
  or deletes.
- **All output stays under `lex_testing/`.** Nothing is written to
  `results/`, `work/`, or anywhere else in the repo.

## What gets measured

| Metric | How |
|---|---|
| Wall-clock seconds | `date` around the `docker run` |
| Peak memory (MiB) | `docker stats` polled every 5s against the running container |
| Graph size | `S` and `L` line counts plus total segment bp in the output GFA |

Peak memory is sampled, not exact — a spike between polls is missed. Treat
it as a good estimate, not a hard bound.

## Usage

```bash
export DXPROJ=project-JB6zQBj0ZQv2Bk79ggBBv76Z

bash lex_testing/fetch_inputs.sh both        # read-only download

# One graph, all sequences together (the monolithic baseline)
bash lex_testing/benchmark.sh lex_testing/inputs/smoke_1mb.fa all 16 smoke

# Reference + one haplotype at a time — N small graphs, linear in samples
bash lex_testing/benchmark.sh lex_testing/inputs/smoke_1mb.fa pairwise 16 smoke

# Reference + 1, +2, +3, +4 haplotypes — the growth curve
bash lex_testing/benchmark.sh lex_testing/inputs/smoke_1mb.fa cumulative 16 smoke

bash lex_testing/summary.sh                  # cost table so far
```

Run a single input directly if you prefer:

```bash
bash lex_testing/pggb_run.sh <input.fa> <label> [threads]
```

Tune PGGB via environment variables: `P_IDENT` (`-p`, default 90),
`S_SEG` (`-s`, 5000), `K_MIN` (`-k`, 29), `PGGB_IMAGE`.

## What lands where

```
lex_testing/
├── inputs/          downloaded FASTAs + generated subsets   (gitignored: *.fa)
├── runs/<label>/    PGGB output, logs, the GFA              (gitignored: *.gfa)
└── metrics/runs.tsv one row per run — COMMITTED
```

The repo `.gitignore` excludes `*.fa` and `*.gfa` globally, so the large
files stay local while `metrics/runs.tsv` is version-controlled. That TSV
is the actual result.

## Reading the numbers

PGGB's `wfmash` step aligns every sequence against every other, so cost
grows roughly with the **square** of the sequence count, while `pairwise`
runs grow linearly. Comparing total `pairwise` seconds against the single
`all` run is the core measurement: if 4 pairwise runs cost meaningfully
less than 1 combined run, the parallel-construction premise holds at the
sample axis.

Note this measures the **sample** axis. The project's central experiment is
the **region** axis — splitting a chromosome into chunks. Same instrument,
different inputs.
