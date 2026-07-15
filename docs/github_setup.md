# GitHub workflow

## Recommended repository visibility

Use a private repository while animal-level metadata and internal methods are still being organized.

## Branches

- `main`: working reproducible pipeline
- `feature/<name>`: one contained addition, such as `feature/temperature-qc`
- `fix/<name>`: bug corrections

## Commit rule

Keep organizational changes separate from scientific-method changes. For example:

```text
refactor: move scoring helpers into src/scoring
fix: enforce strict <15 s microarousal rule
feat: add temperature plateau QC
```

## Data rule

Never force-add ignored raw files. If Git reports a MAT-file as staged, remove it before pushing:

```bash
git restore --staged path/to/file.mat
```
