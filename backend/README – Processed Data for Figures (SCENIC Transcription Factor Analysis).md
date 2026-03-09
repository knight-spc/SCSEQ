# SCENIC

This repository contains the **processed data used to generate the transcription factor analysis figures** shown in the manuscript. These results were produced using the **SCENIC-based regulon activity analysis module** implemented in the SCSEQ platform.

The files correspond to four visualization panels illustrating different aspects of transcription factor regulatory activity across cell types.

------

## Figure: UMAP (Regulon Activity)

**File:**
 `figure_umap_regulon_activity.csv`

**Description**

This file contains the **UMAP embedding coordinates and regulon activity scores** for individual cells.

Each row corresponds to a single cell and includes:

- UMAP_1
- UMAP_2
- cell_type
- selected_regulon_activity (e.g., IRF7)

**Purpose**

This dataset was used to generate the **UMAP visualization of regulon activity**, which shows the spatial distribution of transcription factor activity across the cellular landscape.

The visualization highlights how specific regulons are enriched in particular immune cell populations.

------

## Figure: DotPlot (Top Regulons Across Cell Types)

**File:**
 `figure_dotplot_top_regulons.csv`

**Description**

This file summarizes the **activity of top regulons across different cell types**.

Columns include:

- regulon_name
- cell_type
- mean_activity (average regulon activity score)
- pct_active_cells (percentage of cells with active regulon)

**Purpose**

This dataset was used to generate the **DotPlot overview figure**.

In the visualization:

- Dot size represents the **fraction of cells with active regulon activity**
- Dot color represents the **average regulon activity level**

This figure provides a **global overview of candidate transcription factors regulating specific cell populations**.

------

## Figure: Heatmap (Global Regulatory Pattern)

**File:**
 `figure_heatmap_regulon_matrix.csv`

**Description**

This file contains the **cell-type-level averaged regulon activity matrix**.

Rows represent:

- transcription factor regulons

Columns represent:

- cell types

Values represent:

- **average regulon activity scores (AUC)** across cells of the same type.

**Purpose**

This dataset was used to generate the **heatmap of global regulatory patterns**.

The heatmap reveals similarities and differences in transcriptional regulatory programs between cell populations and highlights clusters of cell types sharing similar regulatory features.

------

## Figure: Violin Plot (Key Transcription Factors)

**File:**
 `figure_violin_key_regulons.csv`

**Description**

This dataset contains **single-cell regulon activity values for selected transcription factors**.

Columns include:

- cell_id
- cell_type
- regulon_name
- regulon_activity_score

**Purpose**

This dataset was used to generate **violin plots showing the distribution of transcription factor activity across cell types**.

These plots highlight cell-type-specific regulatory activity and allow detailed inspection of regulatory heterogeneity at the single-cell level.

------

## Notes

- All regulon activity scores were derived using the **SCENIC workflow implemented in the SCSEQ platform**.
- The PBMC dataset from 10x Genomics was used as the demonstration dataset in the manuscript.
- The processed files provided here allow readers to **reproduce the figures shown in the manuscript**.