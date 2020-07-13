# FORC-Analysis
[![DOI](https://zenodo.org/badge/279344258.svg)](https://zenodo.org/badge/latestdoi/279344258)

First Order Reversal Curve (FORC) data and code

This is an archive of the code and data used to reproduce the analysis of FORC data and corresponding figures presented in "Quality Metrics for First Order Reversal Curve (FORC) Measurements and Analysis" by C.L. Dennis, S.P. Lund, and R.D. Shull.

The file "FORC02--14.dat" contains measurements collected using DC Head with a free center point and with a fixed center point.
The file "FORC011714.dat" contains measurements collected using VSM Head with autoranging. 
The file "FORC012514.dat" contains measurements collected using VSM Head with fixed range.

The file "FORC analysis and figures.R"  is the main file to load the FORC data, analyze it, and produce pdf files that contain the figures in the paper, along with many other plots.  This file sources both "LoadFORCdata2.R" and "Spline plots.R".
The file "LoadFORCdata2.R" loads and formats FORC data when sourced.
The file 'Spline plots.R' produces a set of plots specific to the sequential 1-d cubic spline analysis methods when sourced.
