# Halstock Wood SSSI Habitat Survey Analysis

**Author:** Dan A. Holloway
**Contact:** Dan2.holloway@live.uwe.ac.uk
**Date:** 10/08/2026
**License:** MIT License

## Overview
This R script explores understorey plant community composition along a riparian distance gradient in Halstock Wood SSSI. It uses 3D Non-metric Multidimensional Scaling (NMDS), environmental vector fitting (`envfit`), PERMANOVA, and multivariate dispersion analysis (`betadisper`) to test whether, and how, distance from the river shapes vegetation community structure and its relationship to soil/canopy conditions.

## Required Data
The script expects the following two files in the working directory:
- `Species_Data.csv` – quadrat-level species abundance data (`plot_id`, `transect`, `quadrat`, plus 26 species columns)
- `Env_Factors.csv` – quadrat-level environmental data (`plot_id`, `transect`, `distance_from_river_m`, `soil_moisture_mg_g`, `soil_ph`, `canopy_cover_pct`)

See the dataset README for full column definitions.

## Required R Packages
```r
install.packages(c("vegan", "ggplot2", "ggrepel", "vegan3d", "patchwork", "plotly", "htmlwidgets"))
```
- **vegan** – NMDS, envfit, PERMANOVA (adonis2), betadisper, ordisurf
- **ggplot2** / **ggrepel** – static plots and non-overlapping vector labels
- **vegan3d** – static 3D ordination plot (`ordiplot3d`)
- **patchwork** – combining multiple ggplot figures into one panel
- **plotly** / **htmlwidgets** – interactive, rotatable 3D NMDS plot saved as HTML
