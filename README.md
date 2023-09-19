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

The well known CHB-MIT scalp EEG database (https://physionet.org/content/chbmit/1.0.0/) can be used for evaluating the model performance, which cantains 24 subjects with a total of 198 seizures. However, the original .edf files should be coverted to .csv files for R to read. The example .csv header is listed in example_head.csv, the column arangement should be the same with this file. We note that among the columns of example_head.csv, label1 = 0 means non-seizure period, and label1 = 1 means seizure period. The label2 column is not used for the CHB-MIT dataset.

We have already converted the .edf data of the first subject chb01 to .csv files. One can download and unzip the data from google drive (https://drive.google.com/file/d/1d5YaM1AQw_06nJbtKVS4eqfLA7QCB-6w/view?usp=sharing) into /data/chb01/ directory to run the example codes.

## Getting start

1. Download or clone this repo to your local environment;
2. Install all the required R packages as listed in main.R;
3. Download and unzip the example data from google drive to /data/chb01/;
4. Run main.R.

## Time complexity

For such a high frequency multi-channel EEG data, our procedure required 20 minutes to obtain the final results for one hour EEG data in an AMD Milan EPYC 7T83 (2.45 GHz) platform using 10 cores and 20 processes.
