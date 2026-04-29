# Rugulopteryx okamurae colonisation dynamics

This repository contains the R script used for the statistical analyses and visualisations associated with the manuscript:

**Monitoring the colonisation dynamics of *Rugulopteryx okamurae* on algal communities in the Northwestern Mediterranean Sea**

Authors: Marie Borriglione, Thierry Thibaut, Bastien Thouroude, Aurélie Blanfuné, Frédéric Zuberer, Dorian Guillemain, Delphine Thibault, Charles-François Boudouresque, Sandrine Ruitton  
Journal: *Aquatic Conservation*  
Year: 2026

## Overview

The analyses are based on two datasets:

1. `Biomasse_rugulopteryx.xlsx`  
   Biomass sampling of *Rugulopteryx okamurae* across two depths, 5 m and 15 m, at two sites in Marseille, France.

2. `data_quadrat-colonisation.xlsx`  
   Photoquadrat data used to analyse benthic community colonisation after experimental scraping.

## Analyses included

### `Biomasse_rugulopteryx.xlsx`

The R script performs:
- Gamma GLMM analysis of *R. okamurae* dry biomass;
- model selection using AIC;
- residual diagnostics using DHARMa;
- likelihood ratio tests for fixed effects;
- post-hoc pairwise comparisons using `emmeans`;
- biomass and seawater temperature visualisation;

### `data_quadrat-colonisation.xlsx`

The R script performs:
- stacked barplots of benthic community composition;
- Bray-Curtis NMDS ordination;
- species vector fitting using `envfit`;
- PERMANOVA with constrained permutations within quadrats;
- PERMDISP analysis of multivariate dispersion;

## Input data

The script expects the following files in the working directory, available on zenodo:

```text
Biomasse_rugulopteryx.xlsx
data_quadrat-colonisation.xlsx
