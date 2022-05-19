# chem_shift_imaging_NMR

This repository contains the point-wise binned NMR spectra as well as the data analysis R-code for the study "Chemical shift imaging by NMR metabolomics using in tube extraction and slice selection." by Silvio Waschina (Kiel University), Karsten Seeger (University of Lübeck).

## Dependencies

- ***R*** (version >4.1.2)
- R-packages
  - `data.table` (v. 1.14.2)
  - `speaq` (v. 2.6.1)
  - `ggplot2` (v. 3.3.5)
  - `ggrepel` (v. 0.9.1)
  - `bit64` (v. 4.0.5)
  - `egg` (v. 0.4.5)
  - `missForest` (v. 1.4)

## Repository structure

- directory `analysis/va/scripts/` - contains all R-scripts
- directory `analysis/va/plots/` - should be created if not already there. Plots in PDF format will be saved here.
- directory `data/clean/` contains the sample meta information (e.g. coordinates of samples)
- directory `data/raw/` contains the actual data (point-wise binned spectra). Files are gzipped plain text files. 

## Data analysis

The data analysis is performed using R. To re-run the analysis, simply run the following commands from a bash terminal within the project directory.

##### Bear-shaped ham

```bash
Rscript analysis/v1/baerchi_01.R
```

##### Bierschinken

```bash
Rscript analysis/v1/beerschi_01.R
```

