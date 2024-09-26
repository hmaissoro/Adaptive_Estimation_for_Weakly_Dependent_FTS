
# Adaptive Estimation for Weakly Dependent Functional Time Series

This repository contains real data and R code for reproducing the numerical study in the paper [Adaptive Estimation for Weakly Dependent Functional Time Series](https://arxiv.org/abs/2403.13706) available on arXiv.

The real data supporting the findings of the study are publicly available \citep[see][]{misc_individual_household_electric_power_consumption_235} and are also copied in the folder `real_data_application`. The data for the Monte Carlo study can be generated using the code provided in the folder `simulations`.

## Contents

The numerical study is divided into two main parts:

- **Simulation:** The scripts to reproduce the simulation results are in the `simulations` folder.
- **Real Data Application:** The scripts to reproduce the real data analysis are in the `real_data_application` folder.

## Important Information

Each script in both folders is numbered from 01 to 07 and must be run in sequence.

At the beginning of each script, the functions from the `R` folder need to be loaded. You will see the line:

```r
all_func <- list.files("../R/")
```

You may need to adjust this line to match your working directory.

## Results

The results produced by the scripts are stored in the following folders:

- `paper_graphs`: Contains the graphs used in the paper.
- `paper_tables`: Contains the tables used in the paper.
