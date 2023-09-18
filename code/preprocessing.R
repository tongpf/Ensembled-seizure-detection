train_test_split = function(training_sample_prop,sample_freq,non_seizure_segment_length,
                            seizure_segment_expand_length,allfiles){
  # read the seizure files, prepare to split data into training and test sets
  labelinfo = list()
  data_segment_length = 0
  for (i in c(1:length(allfiles))){
    csvfile = allfiles[i]
    EEGdata_label = fread(paste0(datadir,csvfile),
                          header=TRUE,data.table=FALSE,
                          select = c('index','hour','minute','sec','hertz','label1'))
    EEGdata_label$file_idx = i
    EEGdata_label$file_idx_start = data_segment_length+1  # save the file start index
    data_segment_length = data_segment_length+nrow(EEGdata_label)
    EEGdata_label$file_idx_end = data_segment_length  # save the file end index
    labelinfo[[csvfile]] = EEGdata_label
    rm(EEGdata_label)
  }
  labelinfo = do.call(rbind, labelinfo)
  total_sample_size = nrow(labelinfo)
  row.names(labelinfo)  = NULL
  # seq_idx = seq(1,nrow(labelinfo),1000)
  # test = labelinfo[seq_idx,]
  seizure_period = which(labelinfo$label1 == 1)
  
  seizure_period_train_length = as.integer(length(seizure_period)*training_sample_prop)
  seizure_period_diff = diff(seizure_period)
  seizure_period_start = c(1,which(seizure_period_diff>1)+1)
  seizure_period_end = c(which(seizure_period_diff>1),length(seizure_period))
  seizure_period_start = seizure_period[seizure_period_start]
  seizure_period_end = seizure_period[seizure_period_end]
  
  seizure_period_train = c()
  effect_seizure_length = 0
  while(effect_seizure_length<seizure_period_train_length){
    seizure_idx = sample(c(1:length(seizure_period_start)),1)
    
    seizure_period_start_i = seizure_period_start[seizure_idx]
    seizure_period_end_i = seizure_period_end[seizure_idx]
    effect_seizure_length = 
      effect_seizure_length + seizure_period_end_i - seizure_period_start_i + 1
    
    # the expanded period should not cross files
    seizure_period_start_i = max(seizure_period_start_i - seizure_segment_expand_length*sample_freq,
                                 labelinfo$file_idx_start[seizure_period_start_i])
    seizure_period_end_i = min(seizure_period_end_i + seizure_segment_expand_length*sample_freq,
                               labelinfo$file_idx_end[seizure_period_start_i])
    
    seizure_period_train = c(seizure_period_train,c(seizure_period_start_i:seizure_period_end_i))
    
    seizure_period_start = seizure_period_start[-seizure_idx]
    seizure_period_end = seizure_period_end[-seizure_idx]
  }
  
  #evenly select files to extract non-seizure periods
  non_seizure_train_length = nrow(labelinfo)*training_sample_prop - length(seizure_period_train)
  non_seizure_train_num = ceiling(non_seizure_train_length/non_seizure_segment_length/sample_freq)
  non_seizure_train_jump = as.integer(nrow(labelinfo)/non_seizure_train_num)
  non_seizure_start_idx = sample(c(1:non_seizure_train_jump),1)
  non_seizure_train_files = seq(non_seizure_start_idx,nrow(labelinfo),non_seizure_train_jump)
  non_seizure_train_files = labelinfo$file_idx[non_seizure_train_files]
  
  labelinfo = labelinfo[seizure_period_train,]
  seizure_train_files_unique = unique(labelinfo$file_idx)
  non_seizure_train_files_unique = unique(non_seizure_train_files)
  train_files_unique = sort(unique(c(seizure_train_files_unique, non_seizure_train_files_unique)))
  
  #extract the training data
  EEGdata = list()
  for (i in train_files_unique){
    csvfile = allfiles[i]
    EEGdata_temp = fread(paste0(datadir,csvfile),
                         header=TRUE,data.table=FALSE,drop = 1)
    EEGdata_temp = EEGdata_temp[,c(-26,-30)] #remove the duplicate columns P7-T7 and T8-P8
    
    #Butterworth filter
    EEGdata_temp[8:28] = EEGdata_temp[8:28]*1e6
    for (channel_idx in c(8:ncol(EEGdata_temp))){
      series = EEGdata_temp[,channel_idx]
      #plot(series,type = 'l')
      bf <- butter(4, c(0.53*2/sample_freq,70*2/sample_freq), type="pass")
      series <- filtfilt(bf, series)
      #plot(series,type = 'l')
      bf <- butter(4, c(49*2/sample_freq,51*2/sample_freq), type="stop")
      series <- filtfilt(bf, series)
      #plot(series,type = 'l')
      EEGdata_temp[,channel_idx] = series
    }
    #remove the first and last 10 seconds after filtering
    EEGdata_temp = EEGdata_temp[c((sample_freq*10):(nrow(EEGdata_temp)-sample_freq*10)),]
    
    #select the training periods
    file_index = c()
    if (i %in% seizure_train_files_unique){
      labelinfo_i = labelinfo[which(labelinfo$file_idx==i),]
      file_index = labelinfo_i$index - labelinfo_i$file_idx_start
    }
    
    if (i %in% non_seizure_train_files_unique){
      non_seizure_num = length(which(non_seizure_train_files_unique==i))
      while (non_seizure_num > 0){
        non_seizure_start_idx = sample(c(1:(nrow(EEGdata_temp)-non_seizure_segment_length*sample_freq)),1)
        file_index_temp = c(non_seizure_start_idx:(non_seizure_start_idx+non_seizure_segment_length*sample_freq-1))
        # non_seizure_period should not contain any selected period or seizure period
        if ((!any(EEGdata_temp$label1[file_index_temp]==1))&
            (length(intersect(file_index_temp,file_index))==0)){
          file_index = sort(unique(c(file_index,file_index_temp)))
          non_seizure_num = non_seizure_num - 1
        }
      }
    }
    file_index = sort(unique(file_index))
    
    EEGdata[[csvfile]] = EEGdata_temp[file_index,]
    rm(EEGdata_temp)
  }
  EEGdata = do.call(rbind, EEGdata)
  training_sample_size = nrow(EEGdata)
  row.names(EEGdata)  = NULL
  
  training_sample_ratio = training_sample_size/total_sample_size
  print(paste0('Exact proportion of the training dataset: ', training_sample_ratio))
  return(EEGdata)
}

SC_coeff_cal = function(EEGdata, ICnum, verbose = FALSE, ifplot = FALSE){
  # derive the spatially constrained independent components and 
  # their corresponding mixing matrix coefficients
  datamat = EEGdata[,c(8:28)]
  ICAresult = fastICA(datamat,ncol(datamat), alg.typ = "parallel",fun = c('logcosh'), 
                      alpha = 1, method = "C", row.norm = FALSE, maxit = 1000, 
                      tol = 0.0001, verbose = verbose)
  Smat = ICAresult$S
  df =as.data.frame(cbind(EEGdata[,c(1)],Smat))
  colnames(df)[1] = 'X'
  
  if(ifplot){
    # plot the seizure onset zone
    myelectrodes <- c("AF7","FT7","TP7","PO7",
                      "AF3","FC3","CP3","PO3",
                      "AF4","FC4","CP4","PO4",
                      "AF8","FT8","TP8","PO8",
                      "FCz","CPz",
                      "FT9","FT10","T10")
    #Create a function to generate a continuous color palette
    rbPal <- colorRampPalette(c('blue','red'))
    Aus = ICAresult$A
    Aus = abs(Aus)  #only absolute value matters
    # This adds a column of color values based on the mixing matrix
    for (i in 1:nrow(Aus)){
      Col <- rbPal(10)[as.numeric(cut(Aus[i,],breaks = 10))]
      jpeg(paste0(plotdir,'the ',i,'-th IC topology.jpg'), width = 1000, height = 1000)
      eegcap2d(electrodes = myelectrodes,col.point = Col,cex.point = 15,
               cex.border = 15, cex.label = 3,main = paste0('the ',i,'-th IC topology'))
      dev.off()
    }
  }
  
  #using labels to select seizure IC
  #correlation
  comparing_raw_signal = EEGdata[which(EEGdata$label1==1),8:28]
  comparing_ICs = Smat[which(EEGdata$label1==1),]
  
  corrdim = c()
  corrmax = c()
  for (i in c(1:ncol(comparing_ICs))){
    corri = 0
    for (j in c(1:ncol(comparing_raw_signal))){
      corri = corri + cor(comparing_ICs[,i],comparing_raw_signal[,j],method = 'spearman')
    }
    corri = abs(corri/ncol(comparing_raw_signal))
    corrdim = c(corrdim,i)
    corrmax = c(corrmax,corri)
    print(i)
  }
  corrdim = corrdim[order(corrmax,decreasing = TRUE)]
  corrmax = corrmax[order(corrmax,decreasing = TRUE)]
  print(corrdim)
  print(corrmax)
  
  #energy
  noseizure_ICs = Smat[-which(EEGdata$label1==1),]
  energydim = c()
  energydiff = c()
  for (i in c(1:ncol(comparing_ICs))){
    meanenergy1 = mean((noseizure_ICs[,i])^2)
    meanenergy2 = mean((comparing_ICs[,i])^2)
    energydim = c(energydim,i)
    energydiff = c(energydiff,meanenergy2/meanenergy1)
  }
  energydim = energydim[order(energydiff,decreasing = TRUE)]
  energydiff = energydiff[order(energydiff,decreasing = TRUE)]
  print(energydim)
  print(energydiff)
  
  #energydim = energydim[energydiff>1]
  
  #select appropriate ICs
  i = 1
  seizure_dim = c()
  while (length(seizure_dim)<ICnum){
    seizure_dim = intersect(energydim[1:i],corrdim[1:i])
    i = i+1
  }
  seizure_dim = seizure_dim[1:ICnum]
  
  #try to build a spatial constrained ICA algorithm
  #If we want to fix the effect of column k in S, we should fix the Kth row in A
  
  refrow = matrix(ICAresult$A[seizure_dim,],nrow = length(seizure_dim))
  return(refrow)
}

CUSUM_stat = function(EEGdata,refrow,sample_freq,ICnum,CUSUM_q,diag = FALSE,ifplot = TRUE){
  final_variance_df = data.frame('label' = NA, 'variance' = NA)
  final_ratio_variance_df = data.frame('label' = NA, 'variance' = NA)
  
  final_thrdf = list()#for post processing
  final_cp_thrdf = list()#for post processing
  if (diag){
    thrdf = as.data.frame(matrix(rep(0,times = ICnum),nrow = 1))
    thrdf = cbind(thrdf,'critical value (q%)')
    thrdf[,ncol(thrdf)] = as.character(thrdf[,ncol(thrdf)])
    colnames(thrdf)[ncol(thrdf)] = 'cptype'
  }else{
    thrdf = as.data.frame(matrix(rep(0,times = ICnum*(ICnum+1)/2),nrow = 1))
    thrdf = cbind(thrdf,'critical value (q%)')
    thrdf[,ncol(thrdf)] = as.character(thrdf[,ncol(thrdf)])
    colnames(thrdf)[ncol(thrdf)] = 'cptype'
  }
  
  #determine the scales (frequencies) used for CUSUM
  #target frequency range are 2Hz-60Hz
  max_scale = ceiling(log2(sample_freq/2))
  min_scale = ceiling(log2(sample_freq/60))
  scale_range = c(max_scale:min_scale)
  
  for(i in c(1:max_scale)){
    final_thrdf[[i]] = thrdf[1,]
  }
  
  thrdf[1,ncol(thrdf)] = 0
  colnames(thrdf)[ncol(thrdf)] = 'cp_location'
  thrdf[,ncol(thrdf)] = as.numeric(thrdf[,ncol(thrdf)])
  for(i in c(1:max_scale)){
    final_cp_thrdf[[i]] = thrdf[1,]
  }
  
  cp_stat_list = list()
  cp_stat_list[['final_variance_df']] = final_variance_df
  cp_stat_list[['final_ratio_variance_df']] = final_ratio_variance_df
  cp_stat_list[['final_thrdf']] = final_thrdf
  cp_stat_list[['final_cp_thrdf']] = final_cp_thrdf
  
  trainindex_diff = diff(EEGdata$index)
  trainindex_start = c(1,which(trainindex_diff>1)+1)
  trainindex_end = c(which(trainindex_diff>1),nrow(EEGdata))
  
  trainindex_non_seizure = c()
  trainindex_seizure = c()
  
  for (i in c(1:length(trainindex_start))){
    if (any(EEGdata$label1[trainindex_start[i]:trainindex_end[i]]==1)){
      trainindex_seizure = c(trainindex_seizure,i)
    }else{
      trainindex_non_seizure = c(trainindex_non_seizure,i)
    }
  }
  
  if(length(trainindex_non_seizure)>10){
    trainindex_non_seizure_cp_thr = sort(sample(trainindex_non_seizure,10))
    trainindex_non_seizure_cp = setdiff(trainindex_non_seizure,trainindex_non_seizure_cp_thr)
  }else{
    trainindex_non_seizure_cp_thr = trainindex_non_seizure
    trainindex_non_seizure_cp = c()
  }
  
  
  #training in no seizure period
  for (training_idx in trainindex_non_seizure_cp_thr){
    print('change point threshold selection from non-seizure period')
    start_time = trainindex_start[training_idx]
    end_time = trainindex_end[training_idx]
    
    cp_stat_list = cp_threshold_cal(EEGdata,start_time,end_time,refrow,ICnum,cp_stat_list,
                                    CUSUM_q,scale_range,sample_freq,diag)
  }
  
  #obtain the critical value for CUSUM
  final_thrdf = cp_stat_list[['final_thrdf']]
  final_CUSUM_thrdf = list()
  for (i in c(3:length(final_thrdf))){
    thrdf = final_thrdf[[i]]
    #save the median of the bootstrapped critical values
    thrdf = thrdf[which(thrdf$cptype=='critical value (q%)'),-ncol(thrdf)]
    thrdf = thrdf[-1,]
    #final_CUSUM_thrdf[[i]] = as.numeric(colMeans(thrdf,na.rm = T))
    final_CUSUM_thrdf[[i]] = apply(thrdf,2, function(x) median(x, na.rm = TRUE))
  }
  
  if (length(trainindex_non_seizure_cp)>0){
    for (training_idx in trainindex_non_seizure_cp){
      start_time = trainindex_start[training_idx]
      end_time = trainindex_end[training_idx]
      
      cp_stat_list = cp_threshold_cal(EEGdata,start_time,end_time,refrow,ICnum,cp_stat_list,
                                      CUSUM_q,scale_range,sample_freq,diag,
                                      thr = final_CUSUM_thrdf)
    }
  }
  
  for (training_idx in trainindex_seizure){
    print('plot for the seizure period')
    start_time = trainindex_start[training_idx]
    end_time = trainindex_end[training_idx]
    
    cp_stat_list = cp_threshold_cal(EEGdata,start_time,end_time,refrow,ICnum,cp_stat_list,
                                    CUSUM_q,scale_range,sample_freq,diag,
                                    thr = final_CUSUM_thrdf,ifplot = ifplot)
  }
  
  final_variance_df = cp_stat_list[['final_variance_df']]
  final_ratio_variance_df = cp_stat_list[['final_ratio_variance_df']]
  final_thrdf = cp_stat_list[['final_thrdf']]
  final_cp_thrdf = cp_stat_list[['final_cp_thrdf']]
  final_variance_df = final_variance_df[-1,]
  final_ratio_variance_df = final_ratio_variance_df[-1,]
  for(scale_i in scale_range){
    final_thrdf[[scale_i]] = final_thrdf[[scale_i]][-1,]
    final_cp_thrdf[[scale_i]] = final_cp_thrdf[[scale_i]][-1,]
  }
  if(ifplot){
    plotdf = final_variance_df
    plotdf$label = as.character(plotdf$label)
    p = ggplot(plotdf)+geom_boxplot(aes(y=variance,color = label))
    p = p + theme_bw() + ylab('variance') + scale_color_npg()+
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14))
    ggsave(paste0(plotdir,'training_data_local_variance.jpg'), 
           width=15, height=12,units = "in")
    
    plotdf = final_ratio_variance_df
    plotdf$label = as.character(plotdf$label)
    p = ggplot(plotdf)+geom_boxplot(aes(y=variance,color = label))
    p = p + theme_bw() + ylab('variance') + scale_color_npg()+
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14))
    ggsave(paste0(plotdir,'training_data_local_variance_ratio.jpg'), 
           width=15, height=12,units = "in") 
    for(scale_i in scale_range){
      plotdf = reshape2::melt(final_thrdf[[scale_i]], id.vars = "cptype")
      p = ggplot(plotdf,aes(x=variable))+geom_boxplot(aes(y=value,color = cptype))
      p = p + theme_bw() + ylab('threshold') + scale_color_npg()+
        theme(axis.text=element_text(size=14),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.title=element_text(size=14,face="bold"),
              legend.text=element_text(size=14),
              legend.title=element_text(size=14))
      ggsave(paste0(plotdir,'training_data_threshold_at_level_',scale_i,'.jpg'), 
             width=15, height=12,units = "in")
    }
  }
  cp_stat_list = list()
  cp_stat_list[['final_thrdf']] = final_thrdf
  cp_stat_list[['final_cp_thrdf']] = final_cp_thrdf
  cp_stat_list[['final_CUSUM_thrdf']] = final_CUSUM_thrdf
  cp_stat_list[['scale_range']] = scale_range
  return(cp_stat_list)
}