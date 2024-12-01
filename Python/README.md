# Calculating mean TIR based on CGM data

A Python tool implementing the methods proposed in **Time-In-Range Analyses of Functional Data Subject to Missing with Applications to Inpatient Continuous Glucose Monitoring**

Please contact qi.yu2@emory.edu with any questions or comments.

## Usage directions 

Please check **tirExample.ipynb** for details on usage.

## A Sample dataset
Please check **dataSample.csv** in **data** folder for a simulated dataset perturbated from a real dataset. Descriptions of the dataset can be found in **tirExample.ipynb**

Please note that our tool needs the input with format similar to **dataSample**, which take the long format. The Naïve estimator and the proposed estimator with non-informative follow-up duration require a minimum dataset including **patient_id**, **glucose**, and **time**. For the proposed estimator with informative follow-up duration, at least one additional covariate is necessary. See more details in **tirExample.ipynb**.

## Source code
Please check **tir.py** in **src** folder for source code that implements the methods proposed in our paper.
