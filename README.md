# Assessing New Predictors with Biased Data Augmented with Summary Statistics


This repository contains the simulation code for the main manuscript and the additional simulation studies reported in Web Appendix B.

## Folder structure

- `main/`: Code for the simulation studies reported in the main text.
  - `main.R`: Main script for running the simulations.
  - `main-fun-parallel.R`: Functions used for parallel simulation.
  - `parallel_other_fun.R`: Additional helper functions required for running the simulations.
  - `cpp/`: Rcpp source files used in the simulation.
  - `res/`: Folder for saving simulation results.
    - `evalua.R`: Functions for evaluating the simulation results.
    - `summary.Rmd`: R Markdown file used to summarize the final results.


- `density-ratio-violation/`: Code for the additional simulation study on density-ratio model misspecification, reported in Web Appendix B.1.
  - `main.R`: Main script for running the simulations.
  - `main-fun-parallel.R`: Functions used for parallel simulation.
  - `parallel_other_fun.R`: Additional helper functions required for running the simulations.
  - `cpp/`: Rcpp source files used in the simulations.
  - `res/`: Folder for saving simulation results.
    - `evalua.R`: Functions for evaluating the simulation results.
    - `summary.Rmd`: R Markdown file used to summarize the final results.


- `Y-model-violation/`: Code for the additional simulation study on outcome-model misspecification, reported in Web Appendix B.2.
  - `main.R`: Main script for running the simulations.
  - `main-fun-parallel.R`: Functions used for parallel simulation.
  - `parallel_other_fun.R`: Additional helper functions required for running the simulations.
  - `tune-param.R`: Functions for calibrating the intercept terms and computing the pseudo-true parameters.
  - `cpp/`: Rcpp source files used in the simulations.
  - `res/`: Folder for saving simulation results.
    - `evalua.R`: Functions for evaluating the simulation results.
    - `summary.Rmd`: R Markdown file used to summarize the final results.
