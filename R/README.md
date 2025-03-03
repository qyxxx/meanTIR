# Calculating mean TIR based on CGM data

An R package implementing the methods proposed in **Time-In-Range Analyses of Functional Data Subject to Missing with Applications to Inpatient Continuous Glucose Monitoring**

Please contact qi.yu2@emory.edu with any questions or comments.

## Installation

You can install the development version of meanTIR like so:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("qyxxx/meanTIR/R", build_vignettes = TRUE)
```

## Usage directions 

Please check vignette for details on usage by:

``` r
vignette("meanTIR", package = "meanTIR")
```

## A Sample dataset
Please check **dataSample.csv** in **data** folder for a simulated dataset perturbated from a real dataset.

Please note that our tool needs the input with format similar to **dataSample**, which take the long format. The Naïve estimator and the proposed estimator with non-informative follow-up duration require a minimum dataset including **patient_id**, **glucose**, and **time**. For the proposed estimator with informative follow-up duration, at least one additional covariate is necessary.
