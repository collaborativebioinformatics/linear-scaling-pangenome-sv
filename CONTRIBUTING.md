# Contributing — BCM SV Hackathon

## Branch Naming

Use feature branches with your name prefix:

```
michael/graph-construction
khoi/linear-pipeline
quang/graph-merge
ali/integration-web
alexander/docs
```

Branches describe **primary focus**, not ownership silos. Anyone can help with anything.

## Workflow

```bash
git checkout main
git pull
git checkout -b <name>/<feature>
# make changes
git add .
git commit -m "Clear description of what changed"
git push -u origin <branch>
# Open a pull request on GitHub
```

## What Not To Commit

- FASTA, GFA, VCF, or BAM files
- Large binary outputs
- DNAnexus tokens or AWS credentials
- `.env` files
- `node_modules/`

## Code Review

- Keep commits small and focused
- Write clear commit messages
- Run tests before pushing: `make test`
- Run the demo: `make demo`