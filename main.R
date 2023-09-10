library(doParallel)
library(foreach) 
library(data.table)
library(signal)
library(fastICA)
library(eegkit)
library(hdbinseg)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(Rfast)
library(Matrix)
library(wbs)
library(wavethresh)
library(data.table)
library(tidyverse)

library(RcppArmadillo)
library(Rcpp)
library(RcppEigen)

library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(mlr3filters)
library(mlr3fselect)
library(mlr3pipelines)
library(e1071)
library(WeightedROC)

rm(list=ls())
set.seed(2023)
datadir = './data/'
plotdir = './plot/'
codedir = './code/'
tempdir = './temp/'

allfiles = dir(paste0(datadir,'/chb01/'))

filehours = 1  # each file contains one hours EEG data
training_sample_prop = 0.05
sample_freq = 256  # sampling frequency in Hz
non_seizure_segment_length = 256  # 256 seconds one piece
seizure_segment_expand_length = 128 # expand 128 seconds for the pre-ictal and post-ictal periods
ICnum = 5 # how many independent components to keep for seizure signals
CUSUM_q = 0.2 # significance level for CUSUM
parallel_core_num = 10
diag = FALSE # if diag = TRUE, only use the variance term for CUSUM test

source(paste0(codedir,'preprocessing.R'))
source(paste0(codedir,'change_point.R'))
source(paste0(codedir,'SCICA.R'))
source(paste0(codedir,'basis_fun.R'))
source(paste0(codedir,'feature_utils.R'))
source(paste0(codedir,'classification.R'))

################################################################################
# preprocessing                                                                #
################################################################################

#We first extract the training dataset from the entire EEG records, 
#and then train and derive the corresponding spatially constrained 
#independent components for further usage

EEGdata = train_test_split(training_sample_prop,sample_freq,non_seizure_segment_length,
                           seizure_segment_expand_length,allfiles)

refrow = SC_coeff_cal(EEGdata,ICnum,verbose = TRUE)

#training data for threshold selection
#here the thresholds can be distinguished into the CUSUM threshold (q% critical value)
#and the threshold for post-processing (selection after all cps have been estimated).
#the CUSUM threshold for post-processing is more conservative.
cp_stat_list = CUSUM_stat(EEGdata,refrow,sample_freq,ICnum,CUSUM_q,diag)

save.image(".RData")

################################################################################
# feature extraction                                                           #
################################################################################

# extract and save features from the entire raw EEG series. The training and test
# sets are marked according to the indexes listed in EEGdata

source(paste0(codedir,'feature_extraction.R'))

save.image(".RData")

################################################################################
# classification                                                               #
################################################################################

features_all = lapply(purrr::transpose(features_all), function(x) do.call(rbind, x))

combined_feature = features_all[[1]]
combined_feature_fixed_segment = features_all[[2]]
fulllabels = features_all[[3]]

#centering features by hour
combined_feature = center_adj(combined_feature)
combined_feature_fixed_segment = center_adj(combined_feature_fixed_segment)

#classifier and performance

positive_sample_thr = 0.9 #seizure segment must contain the proportion of 
                          # true seizure period greater than positive_sample_thr
positive_sample_weight = 10
costs = matrix(c(-10, 20, 0.1, -0.1), nrow = 2)
max_unit_length = 200
min_event_length = 20
max_event_length = 500

EEGRF(combined_feature,positive_sample_thr,positive_sample_weight,
      costs,max_unit_length,min_event_length,max_event_length)

EEGRF(combined_feature_fixed_segment,positive_sample_thr,positive_sample_weight,
      costs,max_unit_length,min_event_length,max_event_length)
