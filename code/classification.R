EEGRF = function(EEG_feature,positive_sample_thr,positive_sample_weight,
                 costs,max_unit_length,min_event_length,max_event_length){
  training_idx = which(EEG_feature$usedfortraining == 1)
  
  EEG_feature_temp = EEG_feature[,c(1,5,6,8,10:ncol(EEG_feature))]
  EEG_feature_temp$label[which(EEG_feature_temp$label>=positive_sample_thr)] = 1
  EEG_feature_temp$label[which(EEG_feature_temp$label<positive_sample_thr)] = 0
  EEG_feature_temp$label = as.factor(EEG_feature_temp$label)

  raw_length = EEG_feature_temp$length
  EEG_feature_temp$weight = 1
  EEG_feature_temp$weight[which(EEG_feature_temp$label==1)] = positive_sample_weight
  task = as_task_classif(EEG_feature_temp,
                         target = 'label', positive = '1')
  
  # assign weight to the minority group
  task$set_col_roles("weight", roles = "weight")
  dimnames(costs) = list(response = c("1", "0"), truth = c("1", "0"))
  cost_measure = msr("classif.costs", costs = costs)
  
  learner = lrn("classif.ranger", importance = "impurity", predict_type = "prob")
  
  train_set = training_idx
  test_set = setdiff(seq_len(task$nrow), train_set)
  
  
  learner$train(task, row_ids = train_set)
  sort(learner$state$model$variable.importance,decreasing = T)

  # found the optimal threshold
  prediction = learner$predict(task, row_ids = test_set)
  tp.fp <- WeightedROC(prediction$prob[,1], 
                       as.numeric(as.character(prediction$truth)), 
                       raw_length[test_set])
  tp.fp$loss = costs[2,1]*tp.fp$FN + costs[1,2]*tp.fp$FP
  prediction$set_threshold(tp.fp$threshold[which.min(tp.fp$loss)-1])
  
  tp.fp[which.min(tp.fp$loss),]
  
  ################################################################################
  # postprocessing                                                               #
  ################################################################################
  
  testdata = EEG_feature[test_set,]
  testdata$label[which(testdata$label>=positive_sample_thr)] = 1
  testdata$label[which(testdata$label<positive_sample_thr)] = 0
  testdata$response = prediction$response
  
  #testdata$label = as.numeric(as.character(testdata$label))
  testdata$response = as.numeric(as.character(testdata$response))
  testdata$response[which(testdata$length>(max_unit_length*sample_freq))] = 0
  
  seizure_event = c(which(testdata$response==1),nrow(testdata))
  
  seizure_event_preend = -Inf
  seizure_event_gap = c()
  seizure_event_length = c()
  event_idx = c()
  event_st = 1
  event_en = 1
  for (i in seizure_event){
    gap_length = testdata$start[i]-seizure_event_preend
    seizure_event_preend = testdata$end[i]
    seizure_event_length = c(seizure_event_length,testdata$end[i]-testdata$start[i])
    seizure_event_gap = c(seizure_event_gap,gap_length)
    #merge adjacent seizure segments
    if (gap_length<sample_freq){#consecutive event
      event_idx = c(event_idx,i)
      event_en = testdata$end[i]
    }else{#new event
      if(i!=seizure_event[1]){
        if(((event_en-event_st)>(min_event_length*sample_freq))&
           ((event_en-event_st)<(max_event_length*sample_freq))){
          st = min(event_idx)
          en = max(event_idx)
          st = min(which(testdata$end>(testdata$start[st]-2*sample_freq)))
          en = max(which(testdata$start<(testdata$end[en]+2*sample_freq)))
          testdata$response[st:en] = 1
          #print(c((min(event_idx)):(i)))
        }else{
          testdata$response[min(event_idx):max(event_idx)] = 0
        }
      }
      event_st = testdata$start[i]
      event_idx = c(i)
    }
  }
  TP = sum(testdata$length[which(testdata$label==1&testdata$response==1)])
  FP = sum(testdata$length[which(testdata$label==0&testdata$response==1)])
  TN = sum(testdata$length[which(testdata$label==0&testdata$response==0)])
  FN = sum(testdata$length[which(testdata$label==1&testdata$response==0)])
  SENS = TP / (TP+FN)
  SPEC = TN / (TN+FP)
  ACC = (TP+TN) / (TP+FN+TN+FP)
  print(paste0('Sensitivity: ', SENS))
  print(paste0('Specificity: ', SPEC))
  print(paste0('Accuracy: ', ACC))
}