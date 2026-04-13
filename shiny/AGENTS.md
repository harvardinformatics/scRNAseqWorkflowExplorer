# Agent Notes For `shiny/`

## Scope
- These notes apply only to `/Users/adamfreedman/workspace/scRNAseqWorkflowExplorer/shiny`.
- This directory currently contains a single Shiny app for comparing scRNA-seq workflow outputs stored as Seurat `.rds` files plus matching bootstrap `.tsv` files under `data/`.

## Code Layout
- `app.R` is the Shiny entrypoint and owns app wiring: package loads, file discovery, input widgets, reactive state, tab layout, `renderPlot`/`renderText`/`downloadHandler`, and any logic that binds UI inputs to plotting or data-loading functions.
- `functions.R` is the home for reusable analysis and plotting helpers. Keep plotting functions here, including future functions that may be called from plain R scripts outside the Shiny app.
- `README.md` documents input expectations and launch instructions.
- `data/` is runtime input only. Do not hardcode sample-file names into app logic.

## Current App Structure
- `app.R` contains:
  - Input-file discovery and pairing helpers for `.rds` and bootstrap `.tsv` files.
  - Label-management helpers for per-method custom names from the Configuration tab.
  - Reactive state for selected files, method choices, loaded Seurat/bootstrap data, and heatmap hover state.
  - UI tabs for Configuration, Cluster stability, Gene heatmap, Cluster Jaccard heatmap, and shared-barcode UpSet plots.
- `functions.R` contains:
  - Inter-method/intra-method stability plotting.
  - Cluster Jaccard heatmap plotting.
  - Shared-barcode UpSet plotting.
  - Gene-expression heatmap and hover-preview plotting.
  - Small shared helpers such as Jaccard computation and cluster-palette generation.

## Placement Rules
- Put new reusable plotting functions in `functions.R`.
- Put new reusable analysis/data-transformation helpers in `functions.R` unless they are tightly coupled to Shiny widget state only.
- Put Shiny-only wiring in `app.R`: new tabs, controls, observers, reactives, outputs, downloads, and session updates.
- If a helper is useful outside the app, move it to `functions.R` even if it is first introduced for a Shiny feature.

## Method Label Rule
- Custom method labels are defined in the Configuration tab and stored in `method_labels()`, keyed by the source `.rds` path.
- The canonical way to turn a method path into a user-facing label is `resolve_method_label(label_map, rds_path)`.
- Any new functionality that presents a method to the user must use the resolved custom label, not the raw filename, anywhere user-facing.

## Required Label Mapping Coverage
- When adding a new feature that references methods, map custom labels through `resolve_method_label()` for:
  - Input choices and selectors.
  - Plot titles, subtitles, axes, legends, facet labels, annotations, and hover text.
  - Status text, notifications, table headers/cells, and download filenames.
  - Any exported artifact metadata or captions shown to the user.
- Internal data structures may still use file paths as stable keys, but the UI boundary must convert them to custom labels.
- If a new feature compares methods and also supports non-Shiny use, pass display labels in as explicit function arguments rather than deriving filenames inside the plotting function.

## Data And Input Assumptions
- Treat `.rds` paths as the stable identifiers for methods inside the app.
- Expect bootstrap tables to normalize to `clusterid`, `max_jaccard`, and `bootstrap_number`.
- Expect Seurat objects to contain `seurat_clusters` in `meta.data`.
- Preserve support for recursive file discovery under `data/`.

## Contribution Guidance
- Prefer small helpers over duplicating file-discovery or label-resolution logic across tabs.
- Preserve validation messages when adding new reactives or outputs; failed input assumptions should be explicit in the UI.
- Keep user-facing labels consistent across the live plot and matching PDF download.
- If a feature adds a new method-oriented visualization, wire the same selected-file subset and custom-label behavior into it rather than creating a parallel selection system.
