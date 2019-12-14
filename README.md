# data-processing 

MScR Project - Stellar rotation formalisms of γ Doradus stars from gravity-mode period spacings

(Updated 14/12/2019)


This repo contains the python codes used for processing the frequency outputs from GYRE and ADIPLS. 

The output files must be named and organized according to the following convention:
(if you're using the run_gyre and run_adipls tools in astero work directories, the files will be 
automatically ready for python processing)

GYRE work directory:
Outputs stored in inidividual directories named gyre_output_M_mm_O_oo_H_hh_R_rr;
inside each directory, it (only) contains a profile_puls_summary.txt for one model.
Frequencies in units of uHz.

ADIPLS work directory:
Outputs for each model named as save_M_mm_O_oo_H_hh_R_rr_mode.data

[Each model has unique set of parameters - M:mass, O:overshoot, H:central hydrogen abundance, R:rotation rate]

# Plot_mode_data.py

A simple tool for plotting period spacing of two user-specified models on a single graph,
one for GYRE, one for ADIPLS

# Plot_mode_data_bench_grid.py

The main code used in our project for comparing Period Spacings of ADIPLS/GYRE grid and benchmark.

Reads both ADIPLS/GYRE benchmark and ADIPLS/GYRE grid. It filters best-matching models with gradient 
and intercept uncertainties of the benchmark.
Extracts parameter set of best-matching models, then plots all best-matching models on a graph 
along with the benchmark.

Evaluates chi2 for all grid models, and creates 2D contour plots showing the chi2 variation 
with each parameter combination in the set.

# Plot_1st_min_rotation_adipls_gyre.py 

Reads GYRE & ADIPLS work directories, each containing models with varying rotation rates. 
Detects the first local minimum in the period spacing curve for each model,
then plots the dip locations against rotation rates of the models. 




















