# Mixing Length interpolator

This project is centered around creating a large Python interpolator for 
predicting the effective temperature of red giant stars, given their mass,
metallicity, and surface gravity. Our predictions are done by interpolating
in a vast, mulit-dimensional grid of stellar evolutionary models. For more 
details, see Tayar, Somers, et al. (2017) -- The Astrophysical Journal, 
Volume 840, Issue 1, article id. 17, 12.

## File descriptions

CreateCoarseGrid.py -- This program was used to produced each of the files in
ProcessedMods/NewBigGrid. These are down-sampled evolutionary models, coving
the red giant branch, from the full grid of models. They have been downsampled
to improve the speed of interpolation.

Interpolator.py -- This is the program which performs the actual interpolation.
It takes in a mass, [Fe/H], [Al/Fe], surface gravity, and effective temperature,
returning the predicted mixing length.

Runner.py -- A wrapper which repeated runs Interpolator over a series of stars 
with various stellar parameters.

ProcessedMods/NewBigGrid/aaa_EmptyTracks.txt -- This file records original YREC
files which were found to be empty. Sometimes the evolutionary calculator chokes
for some reason or another, so we save those models in order to make sure that
we skip them during interpolation.

ProcessedMods/NewBigGrid/* -- The down-sampled models from the original YREC
grid. The name says its mass, [Fe/H], helium, mixing length, and alpha element
enhancement.


