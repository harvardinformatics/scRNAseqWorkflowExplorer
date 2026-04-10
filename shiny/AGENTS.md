# Shiny App Notes

## Scope
This directory contains the Shiny app for comparing clustering outputs derived from the same scRNA-seq input library across multiple analysis methods.

## Local Structure
- `app.R`: app wiring, file discovery, data loading helpers, UI definitions, server reactives, and download handlers.
- `functions.R`: reusable plotting helpers and lower-level data/plot utilities.
- `data/`: local input directory scanned by the app for `.rds` and bootstrap TSV files.

## Data Contract
- Input Seurat objects are `.rds` files under `data/`.
- Method-pair comparisons require a matching bootstrap TSV beside each `.rds`.
- Bootstrap TSV files are normalized by `standardize_bootstrap_cols()` and must resolve to:
  - `clusterid`
  - `max_jaccard`
  - `bootstrap_number`
- Seurat objects are expected to contain `seurat_clusters` in `obj@meta.data`.
- Expression-oriented views assume a count layer is available in the chosen assay.

## Contribution Rules
- Put reusable plotting or data-transformation code in `functions.R`.
- Keep Shiny-specific UI glue, reactive state, and input/output synchronization in `app.R`.
- Prefer extending existing discovery/loading helpers over adding new ad hoc file scans.
- Preserve current data assumptions unless the user explicitly asks to broaden them.
- When adding a new tab or export, make its status text describe:
  - what data subset is loaded
  - the key filter/threshold settings that affect interpretation

## Labeling Rule
- If the app exposes configurable custom method labels in a Configuration tab or equivalent global state, every new method-facing feature in `shiny/` must use that mapping.
- “Method-facing” includes:
  - plot axis labels
  - legends
  - selector choices
  - table headers
  - hover text
  - plot subtitles
  - download filenames
- Raw filenames should only be used as the fallback default label when no custom label has been applied.
- Do not introduce a new tab or output that bypasses the applied custom-label mapping for one view while using it elsewhere.

## Reactive-State Guidance
- Before changing UI state synchronization, inspect how the current server code updates `selectInput()`, `textInput()`, checkbox inputs, and reactive file lists. This app is sensitive to feedback loops between reactives and server-driven input updates.
- If you add selection/configuration state, keep a clear distinction between:
  - pending UI state
  - applied state used by plots and downloads
- Keep expensive reads in reactives that validate inputs early.
- Use `validate()` / `need()` for user-facing data errors instead of allowing raw exceptions to propagate.

## Maintenance Targets
- Keep `README.md` aligned with the actual tabs and workflows implemented in `app.R`.
- If major new Shiny functionality is added, update this file with:
  - the owning helper functions
  - any new required data fields
  - whether the feature must honor global custom method labels
