
# Adaptive Estimation for Weakly Dependent Functional Time Series

This repository contains real data and R code for reproducing the numerical study in the paper [Adaptive Estimation for Weakly Dependent Functional Time Series](https://arxiv.org/abs/2403.13706) available on arXiv.

The real data supporting the findings of the study are publicly available \citep[see][]{misc_individual_household_electric_power_consumption_235} and are also copied in the folder `real_data_application`. The data for the Monte Carlo study can be generated using the code provided in the folder `simulations`.

## Contents

The numerical study is divided into two main parts:

- **Simulation:** The scripts to reproduce the simulation results are in the `simulations` folder.
- **Real Data Application:** The scripts to reproduce the real data analysis are in the `real_data_application` folder.

## Installation
### install using devtools

#. Make sure you have devtools package installed :
```R
install.packages("devtools")
```

then using devtools install the package :

```R
devtools::install_github(
    repo = "allemand-instable/Adaptive_Estimation_for_Weakly_Dependent_FTS", # the repo
    ref = "HEAD", # use the latest commit
    subdir = "adaptiveFTS", # the package is not at the main root but inside the directory adaptiveFTS
    build = TRUE,
    build_manual = FALSE,
    build_vignettes = FALSE
)
```


### build manually
if you are using UNIX system, using z-shell :

```shellscript
zsh build.zsh
```

the compiled package should be available in 
```
...
┣ out
┃ ┗ adaptiveFTS_0.1.0_arm64.tar.gz
                        △
                os and architecure
```

else RUN
```zsh
cd adaptiveFTS
R CMD check .
R CMD build .
R CMD INSTALL .
Rscript -e "devtools::build_manual()"
#                             put the actual name of the compiled file
#                                             ▽
Rscript -e 'install.packages( "./adaptiveFTS_[...].tar.gz"  ,repos = NULL, type = "source")'
```

```R
#               put the actual name of the compiled file
#                               ▽
install.packages( "./adaptiveFTS_[...].tar.gz"  ,repos = NULL, type = "source")
```

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
