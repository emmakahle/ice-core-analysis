# ice-core-analysis

## Code Description
Purpose: Produce climate data sets from measurements of the South Pole ice core

Code Requirements: MATLAB license, curve fitting toolbox, signal processing toolbox, statistics and machine learning toolbox, symbolic math toolbox

Input measurements: high-resolution water isotopes, water isotope diffusion length, annual-layer thickness, depth-age scale for ice core, delta-age (difference between ice and gas age at given depth)

Output data sets: temperature history, accumulation rate history, thinning function for ice core record

## How to run
Before running, select run settings in run_settings.m.
Run the main.m script to run the code.

## Results
Results are automatically plotted in two versions: shaded histograms for each parameter and a summary plot, which shows the mean and standard deviation for each parameter.
Output data sets are saved as .txt files for temperature, accumulation length, and thinning function. Columns in these files are age, mean, and standard deviation for each variable.

## Reference
This code is citeable using DOI: 10.5281/zenodo.4579416

This code accompanies Kahle, E. C. et al. (2020) "Temperature, accumulation rate, and layer thinning from the South Pole ice core (SPC14)" U.S. Antarctic Program (USAP) Data Center. doi: https://doi.org/10.15784/601396.
