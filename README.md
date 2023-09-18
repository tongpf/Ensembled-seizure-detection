# Ensembled-seizure-detection

This repo contains the R source code for the unpublished paper "Ensembled Seizure Detection based on Small Training Samples" (under peer review in IEEE TSP, all right reserved). 

Fully automated seizure detection algorithm can relieve doctors from the routine work of inspecting the long-term electroencephalography EEG data. However, due to the large inter-individual variability of seizure types, the majority of researches still relies on patient specific models to achieve higher detection accuracy.

In this repo, our goal is to propose a highly efficient patient specific seizure detection algorithm using small training samples. For example, in an 24 hours continuous EEG monitoring, an experienced doctor only need to manually label one hour of the data and the algorithm can be trained and take care of the rest of the EEG data.

## Operating systems

This repo has been tested at the following operating system:

* Ubuntu 18.04
* CPU: AMD EPYC 7T83 64-Core Processor (2545 MHz)
* Memory: > 32GB
* Disk space: >30GB
* R version: 4.1.1

## Example data sets



Link of the example data. Download and unzip the data into /data/chb01/ directory. https://drive.google.com/file/d/1d5YaM1AQw_06nJbtKVS4eqfLA7QCB-6w/view?usp=sharing
