# Quantifying bias from dependent left truncation in survival analyses of real world data

## 1) System requirements
- The simulation scripts located in the sims directory require the R software environment, and use the survival and tidyverse packages; these are all available freely on all standard operating systems.
- For this paper, R version 4.02 was used, along with version 1.3.0 of tidyverse and 3.2.3 of survival. In addition, the package doParallel (version 1.0.15) was used in order to run simulations in parallel on a computing node.

## 2) Installation guide
- R and the required pacakges can be installed freely from the [R project](https://www.r-project.org/) and [CRAN](https://cran.r-project.org/).
- Installation time should take less than an hour.

## 3) Demo and instructions for use
- Each simulations script outputs results in an .RData file.
- Using a computing node with 36 cores, these simulations should run in a few hours; on a normal desktop, they may need to be run overnight.
- The RData results used in the paper are stored in the results directory.
- The Rmd file paper_figures.Rmd reads in the RData file and generates the figures used in the paper.

