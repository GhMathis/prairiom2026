# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

`prairiom2026` is an R-based PhD analysis project on multi-kingdom biodiversity in the Camargue (Arles, France) grasslands. It processes amplicon-sequencing data from four kingdoms — plants, bacteria (16S), fungi (ITS2) and viruses — together with environmental data (soil chemistry, habitat, distance to railway, land cover) to study community structure and cross-kingdom interactions.

Code is in R only (no build system, no package). Scripts are run interactively from the RStudio project (`prairiom2026.Rproj`) with the project root as the working directory. Comments and prose are mixed French/English. Translate everything in English when you detect French

## How the pipeline is wired

Scripts in `scripts/` are numbered to reflect a dependency order, not folders. The convention is:

-   `0_*` — raw → clean: read `data/<kingdom>/...` Excel/TSV, filter to Camargue (`Host_code` matches `"CAM"` and excludes `"21_CAM_1[45]"`), write tidy tables to `outputs/clean_data_to_analise/`. Bacteria and fungi additionally build `phyloseq` objects saved as `ps_bacteria.Rdata` / `ps_fungi.Rdata` (decontamination via `microDecon`, taxonomy harmonization via SILVA / UNITE reference tables in `data/bacteria/data_ref/` and `data/fungi/`).
-   `1_*` — rarefaction (`iNEXT`) and spatial buffer-size selection (`siland`).
-   `2_environmental_structure.R` — joins soil / habitat / pasture / railway-distance tables into the grid-level metadata used downstream.
-   `3_*`, `4_*`, `5_*` — per-kingdom community structure analyses (plants, viruses, fungi).
-   `bact_sbm.R`, `simple_sbm_mulit_regne.R`, `multipartite_SBM.R`, `plot_multipartite_sbm.R` — Stochastic Block Model fits (package `sbm`) on bacterial OTUs alone, on each kingdom independently at grid level, and on the multipartite cross-kingdom network. **In progress** (per recent commits): SBM is being lifted from grid level to quadrat level.
-   `functions_modified.R` — shared plotting helpers (matrix/community plots) used by the SBM scripts; source it rather than copying.

Downstream scripts read the cleaned files written by `0_*`/`2_*`, so re-running an early script invalidates everything that follows. There is no orchestrator — run scripts manually in numeric order (and `0_environnemental_data.R` / `0_cleaning_fungi_data.R` before `2_environmental_structure.R`).

## Sample identifier conventions

Sample codes follow the pattern `YY_CAM_GG_QQ` — year / locality / grid / quadrat. Many transformations rely on `str_extract(Host_code, ".._CAM_..")` to roll quadrat data up to grid level (the unit at which environmental and SBM analyses operate). Aggregation from quadrat to grid frequently multiplies covers by 2 before summing — preserve that when adding new kingdoms.

The filter pair `str_detect(..., "CAM")` + `!str_detect(..., "21_CAM_1[45]")` recurs throughout. The excluded grids 14/15 from 2021 are intentionally dropped; keep this filter when adding new kingdoms or scripts so sample sets stay aligned across kingdoms.

## Data and outputs layout

-   `data/` — raw inputs, gitignored by extension (`.xlsx`, `.txt`, `.table2`, `.shp`, …). Treat as read-only.
-   `outputs/clean_data_to_analise/` — canonical cleaned inputs that downstream scripts consume (`cover_plant_quad.txt`, `abund_plant_grid.txt`, `OTU_virus_CAM.txt`, `OTU_fungi*.txt`, `Metadata_grid_CAM.txt`, `ps_bacteria.Rdata`, `ps_fungi.Rdata`, …). Anything written here is a contract with later scripts.
-   `outputs/heavy_computation_save/` — for long-running model fits (SBM in particular). Cache fits here rather than recomputing.
-   `data/Shapefile/crop_shapefile.*` — site polygons used by the README leaflet map and by `siland` buffer analyses.

Note the historical path drift: scripts 2/3 read `data/data_clean/Metadata_*` while `0_*` writes to `outputs/clean_data_to_analise/`. When editing those scripts, expect to update the path.

## Heavy dependencies

These are non-CRAN or Bioconductor and worth knowing about before running anything:

-   Bioconductor: `phyloseq`, `microbiome`, `ComplexHeatmap` (install via `BiocManager::install(...)`).
-   R-universe: `microViz` (`repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))`).
-   GitHub: `microDecon` (`remotes::install_github("donaldtmcknight/microDecon")`).
-   CRAN, but central to the analyses: `sbm`, `iNEXT`, `iNEXT.3D`, `siland`, `vegan`, `FactoMineR`, `factoextra`, `sf`, `terra`, `tidyverse`.

There is no `renv.lock` / `DESCRIPTION` — package state is not pinned. If you add a dependency, mention it in commit messages so it can be reproduced.

## Conventions

-   Always run from the project root; every path is relative (`data/...`, `outputs/...`).
-   Keep using `read_tsv` / `read_xlsx` / `read.table` consistent with the surrounding script — mixed quoting/encoding has bitten this codebase before.
-   `main_theme` (a custom `ggplot2` theme) is redefined at the top of most plotting scripts. Match the existing definition rather than introducing a new one.
-   Commit messages in this repo use gitmoji (`:construction:`, `:mag:`, `:memo:`, `:wastebasket:`) — follow that style.
-   NEVER overwrite `.Rdata` files or saved files from a previous script except when it is explicit said. If you come to a step that need overwriting, notify clearly the user.
