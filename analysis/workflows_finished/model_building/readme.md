# Model Building

This directory contains the work that's gone into building and implementing models for population-specific substitution.  

## Contents

In this directory are two Microsoft exel spreadsheets that illustrate small examples of the models we've made, using dummy data to show the null hypothesis expectations.  There is also a group of report files, saved as .Rmd, .pdf, and .html.  This is R code detailing how I've worked through building the model and estimating the parameters

### Results

The results subdirectory contains two types of files:

 - LRT test output files, which contain counts of each type of substitution and likelihood calculations for all models
 - parameter files, which contain the MLE estimates of the parameters for the specified kmer and population

### Prediction

This directory is still a work in progress.

Here, you can find some of the output files for jobs I've been running on the cluster.  The results of these exploits are saved in the form of .params files, but sadly they are very large, so they're .gitignored for now.
