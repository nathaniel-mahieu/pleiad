# Pleiad

## Introduction

Pleiad is a tool for targeted analyte analysis including both absolute comparisons and isotopologue analysis. The tool performs a few core steps:
Given a metabolite library pleiad exposes helper functions for mass, mz, and isotopologue masses.
Given a (rt, i) EIC and a (rtmin, rtmax) region within that trace pleid will align additional samples to that EIC and return the corresponding regions.
Given a set of target masses pleid will return the nearest masses from supplied scans to the target masses.
Pleiad also performs some simple normalization of spectra and consolidation of multiple scans into a single spectrum and other miscellaneous processing steps.

## Use

The core functions are simple and can be leveraged to extract data from files rapidly. This is the function of pleiad.

Taking the extracted data and doing statistics, plotting, or otherwise interpreting the data is a highly variable proess and can’t be simply codified.  Each experiment will have different experimental factors, desired comparisons, etc. and as such it is expected that later processing and plotting be customized for each experiment.

The package is supplied with an example workflow in which I have used pleiad to extract data and then plot data.  The example R workflow can be adapted to future experiments.

## Reference

To learn what a function does type ?function_name at the R prompt.

#### Processing Functions

 - dtw_getspecs()
 - subspec_nearestmasses()
 - subspec_masscorrect()
 - subspec_spectra()

#### Metabolite Helper Functions

 - m()
 - mz()
 - mass()
 - isotopologues()
 - isotopologues.combn()


#### Install

 - setRepositories(ind = 1:2)
 - devtools::install_github("nathaniel-mahieu/pleiad")

