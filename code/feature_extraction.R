keylevel = sort(cp_stat_list$scale_range)
thr = unlist(cp_stat_list$final_CUSUM_thrdf)

#build the basis function dictionary
cand_spike_freq = c(15:24)
cand_sharp_freq = c(5:14)
cand_wave_freq = c(2,3,4,6,8)
cand_spislo_sharp_freq = c(5,7,10,15,20,25)
#cand_haar_freq = c(2,4,8,16,32,64,128)
maxwaveformlength = 750
jump = 4
halfwaveformlength = as.integer(maxwaveformlength/2)
active_length = 2000
total_length = active_length+maxwaveformlength
{
  xstd_spike<- matrnorm(total_length, 1)
  xstd_spike<- data.table(xstd_spike)
  
  for (freq in cand_sharp_freq){
    cha_waveform = sharp_spike(sample_freq,freq)
    #plot(cha_waveform, type = 'l')
    cha_waveform_length = length(cha_waveform)
    xzero = rep(0,total_length)
    xnew = xzero
    laglength = as.integer((maxwaveformlength-cha_waveform_length)/2)
    xnew[(1:cha_waveform_length)+laglength] = cha_waveform
    xnew_save = xnew
    
    xstd_spike[,c(paste("V",ncol(xstd_spike)+1,sep="")) := xnew]
    
    for (loc in seq(total_length-jump+1,maxwaveformlength+2,by = -jump)){
      xnew = c(xnew_save[loc:(loc+jump-1)],xnew)[1:total_length]
      xstd_spike[,c(paste("V",ncol(xstd_spike)+1,sep="")) := xnew]
    }
  }
  for (freq in cand_spike_freq){
    cha_waveform = sharp_spike(sample_freq,freq)
    cha_waveform_length = length(cha_waveform)
    xzero = rep(0,total_length)
    xnew = xzero
    laglength = as.integer((maxwaveformlength-cha_waveform_length)/2)
    xnew[(1:cha_waveform_length)+laglength] = cha_waveform
    xnew_save = xnew
    xstd_spike[,c(paste("V",ncol(xstd_spike)+1,sep="")) := xnew]
    
    for (loc in seq(total_length-jump+1,maxwaveformlength+2,by = -jump)){
      xnew = c(xnew_save[loc:(loc+jump-1)],xnew)[1:total_length]
      xstd_spike[,c(paste("V",ncol(xstd_spike)+1,sep="")) := xnew]
    }
  }
  xstd_spike=as.matrix(xstd_spike)
  xstd_spike = xstd_spike[,-1]
  xstd_spikemask = matrix(NA,nrow = nrow(xstd_spike),ncol = ncol(xstd_spike))
  xstd_spikemask[xstd_spike!=0] = 1L
  xstd_spike = as(xstd_spike,'dgCMatrix')
  wave_num_spike = 20
  total_num_spike = 20
}

{
  xstd_spislo <- matrnorm(total_length, 1)
  xstd_spislo <- data.table(xstd_spislo)
  for (freq1 in cand_spislo_sharp_freq){
    for (freq2 in cand_wave_freq){
      cha_waveform = spike_and_wave(sample_freq,freq1,freq2)
      #plot(cha_waveform, type = 'l')
      cha_waveform_length = length(cha_waveform)
      xzero = rep(0,total_length)
      xnew = xzero
      laglength = as.integer((maxwaveformlength-cha_waveform_length)/2)
      xnew[(1:cha_waveform_length)+laglength] = cha_waveform
      xnew_save = xnew
      
      xstd_spislo[,c(paste("V",ncol(xstd_spislo)+1,sep="")) := xnew]
      
      for (loc in seq(total_length-jump+1,maxwaveformlength+2,by = -jump)){
        xnew = c(xnew_save[loc:(loc+jump-1)],xnew)[1:total_length]
        xstd_spislo[,c(paste("V",ncol(xstd_spislo)+1,sep="")) := xnew]
      }
    }
  }
  xstd_spislo=as.matrix(xstd_spislo)
  xstd_spislo = xstd_spislo[,-1]
  xstd_spislomask = matrix(NA,nrow = nrow(xstd_spislo),ncol = ncol(xstd_spislo))
  xstd_spislomask[xstd_spislo!=0] = 1L
  xstd_spislo = as(xstd_spislo,'dgCMatrix')
  wave_num_spislo = 30
  total_num_spislo = 30
}

read_file_num = parallel_core_num%/%5
`%mydo%` <- ifelse(read_file_num > 1, `%dopar%`, `%do%`)
#derive the features
cl=makeCluster(read_file_num)
registerDoParallel(cl)
features_all <- foreach (csvfile = allfiles,.inorder=FALSE,
                         .packages=c("ggplot2","reshape2","signal","hdbinseg","wbs",
                                     "wavethresh","doParallel", "Rfast",
                                     "RcppArmadillo","Matrix", 'Rcpp','RcppEigen',
                                     'data.table')) %mydo% 
  {
    EEG_feature = data.frame(
      'label' = NA,
      'usedfortraining' = NA,
      'start' = NA,
      'end' = NA,
      'length' = NA,
      
      'variance' = NA,
      'max_variance' = NA,
      
      'variance_ratio' = NA,
      'max_variance_ratio' = NA,
      
      'front_l_spike' = NA,
      'front_r_spike' = NA,
      'tempo_l_spike' = NA,
      'tempo_r_spike' = NA,
      'parietal_spike' = NA,
      'paramed_l_spike' = NA,
      'paramed_r_spike' = NA,
      'occipital_spike' = NA,
      
      'front_l_spislo' = NA,
      'front_r_spislo' = NA,
      'tempo_l_spislo' = NA,
      'tempo_r_spislo' = NA,
      'parietal_spislo' = NA,
      'paramed_l_spislo' = NA,
      'paramed_r_spislo' = NA,
      'occipital_spislo' = NA,
      
      'front_l_spike_bipolar' = NA,
      'front_r_spike_bipolar' = NA,
      'tempo_l_spike_bipolar' = NA,
      'tempo_r_spike_bipolar' = NA,
      'parietal_spike_bipolar' = NA,
      'paramed_l_spike_bipolar' = NA,
      'paramed_r_spike_bipolar' = NA,
      
      
      'front_l_spike_cor' = NA,
      'front_r_spike_cor' = NA,
      'tempo_l_spike_cor' = NA,
      'tempo_r_spike_cor' = NA,
      'parietal_spike_cor' = NA,
      'paramed_l_spike_cor' = NA,
      'paramed_r_spike_cor' = NA,
      'occipital_spike_cor' = NA,
      
      'front_l_spislo_cor' = NA,
      'front_r_spislo_cor' = NA,
      'tempo_l_spislo_cor' = NA,
      'tempo_r_spislo_cor' = NA,
      'parietal_spislo_cor' = NA,
      'paramed_l_spislo_cor' = NA,
      'paramed_r_spislo_cor' = NA,
      'occipital_spislo_cor' = NA,
      
      'front_l_spike_bipolar_cor' = NA,
      'front_r_spike_bipolar_cor' = NA,
      'tempo_l_spike_bipolar_cor' = NA,
      'tempo_r_spike_bipolar_cor' = NA,
      'parietal_spike_bipolar_cor' = NA,
      'paramed_l_spike_bipolar_cor' = NA,
      'paramed_r_spike_bipolar_cor' = NA
    )
    EEG_feature_fixed_segment = EEG_feature
    
    channel_name = c("X",'FP1-F7','F7-T7','T7-P7','P7-O1',
                     'FP1-F3','F3-C3','C3-P3','P3-O1',
                     'FP2-F4','F4-C4','C4-P4','P4-O2',
                     'FP2-F8','F8-T8','T8-P8','P8-O2',
                     'FZ-CZ','CZ-PZ','T7-FT9','FT9-FT10','FT10-T8')
    
    EEGdata = fread(paste0(datadir,csvfile),
                    header=TRUE,data.table=FALSE,drop = 1)
    EEGdata = EEGdata[,c(-26,-30)] #remove the duplicate columns P7-T7 and T8-P8
    EEGdata[8:28] = EEGdata[8:28]*1e6
    #Butterworth filter
    for (i in c(8:ncol(EEGdata))){
      series = EEGdata[,i]
      #plot(series,type = 'l')
      bf <- butter(4, c(0.53*2/sample_freq,70*2/sample_freq), type="pass")
      series <- filtfilt(bf, series)
      #plot(series,type = 'l')
      bf <- butter(4, c(49*2/sample_freq,51*2/sample_freq), type="stop")
      series <- filtfilt(bf, series)
      #plot(series,type = 'l')
      EEGdata[,i] = series
    }
    EEGdata = EEGdata[c((sample_freq*10):(nrow(EEGdata)-sample_freq*10)),]
    
    EEGdata_SCICA_1d = rowMeans(EEGdata[8:28])
    EEGdata_SCICA_sigmahat = rep(0,times = length(EEGdata_SCICA_1d))
    EEGdata_SCICA_sigmahat_ratio = rep(0,times = length(EEGdata_SCICA_1d))
    
    EEGdata_waveform_rawcoeff = matrix(0,ncol = 8*2, nrow = nrow(EEGdata))
    EEGdata_waveform_rawcoeff = as.data.frame(EEGdata_waveform_rawcoeff)
    colnames(EEGdata_waveform_rawcoeff) = 
      c('front_l_spike','front_r_spike', 'tempo_l_spike','tempo_r_spike','parietal_spike',
        'paramed_l_spike','paramed_r_spike','occipital_spike','front_l_spislo',
        'front_r_spislo','tempo_l_spislo','tempo_r_spislo','parietal_spislo',
        'paramed_l_spislo','paramed_r_spislo','occipital_spislo')
    EEGdata_waveform_corcoeff = matrix(0,ncol = 8*2, nrow = nrow(EEGdata))
    EEGdata_waveform_corcoeff = as.data.frame(EEGdata_waveform_corcoeff)
    colnames(EEGdata_waveform_corcoeff) = 
      c('front_l_spike_cor','front_r_spike_cor','tempo_l_spike_cor','tempo_r_spike_cor',
        'parietal_spike_cor','paramed_l_spike_cor','paramed_r_spike_cor','occipital_spike_cor',
        'front_l_spislo_cor','front_r_spislo_cor','tempo_l_spislo_cor','tempo_r_spislo_cor',
        'parietal_spislo_cor','paramed_l_spislo_cor','paramed_r_spislo_cor','occipital_spislo_cor')
    EEGdata_waveform_bipolar_rawcoeff = matrix(0,ncol = 7, nrow = nrow(EEGdata))
    EEGdata_waveform_bipolar_rawcoeff = as.data.frame(EEGdata_waveform_bipolar_rawcoeff)
    colnames(EEGdata_waveform_bipolar_rawcoeff) = 
      c('front_l_spike_bipolar','front_r_spike_bipolar','tempo_l_spike_bipolar',
        'tempo_r_spike_bipolar','parietal_spike_bipolar','paramed_l_spike_bipolar',
        'paramed_r_spike_bipolar')
    EEGdata_waveform_bipolar_corcoeff = matrix(0,ncol = 7, nrow = nrow(EEGdata))
    EEGdata_waveform_bipolar_corcoeff = as.data.frame(EEGdata_waveform_bipolar_corcoeff)
    colnames(EEGdata_waveform_bipolar_corcoeff) = 
      c('front_l_spike_bipolar_cor','front_r_spike_bipolar_cor','tempo_l_spike_bipolar_cor',
        'tempo_r_spike_bipolar_cor','parietal_spike_bipolar_cor','paramed_l_spike_bipolar_cor',
        'paramed_r_spike_bipolar_cor')
    
    
    #calculate the waveform coefficients first
    if (!file.exists(paste0(tempdir,csvfile,"_waveform_allchannel",".RData"))){
      waveform_allchannel = list()
      
      cl=makeCluster(5)
      registerDoParallel(cl)
      waveform_allchannel = 
        foreach(EEG_sig_full = EEGdata[,c(8:26,28)],.inorder = TRUE,
                .packages = c("signal", "Rfast","RcppArmadillo","Matrix", 'Rcpp',
                              'RcppEigen','data.table'),
                .final = function(x)setNames(x, channel_name[c(-1,-21)])) %dopar% 
        {
          sourceCpp(paste0(codedir,'sparse_eachcol_apply.cpp'))
          stopping_rule_var = 0
          artstart = 1
          
          coefficient_mat_list = list()
          corcoefficient_mat_list = list()
          xfit_spike = as.matrix(cbind(xstd_spike, 1))
          xfit_spislo = as.matrix(cbind(xstd_spislo, 1))
          
          while (artstart + total_length - 1 <= length(EEG_sig_full)) {
            #print(artstart)
            if (artstart + total_length - 1 > length(EEG_sig_full) -active_length) {
              #tail procedure
              # print('tail procedure')
              EEG_sig = EEG_sig_full[(artstart):length(EEG_sig_full)]
              total_length_tail = length(EEG_sig)
              
              xstd_spike_tail<- matrnorm(total_length_tail, 1)
              xstd_spike_tail<- data.table(xstd_spike_tail)
              for (freq in cand_sharp_freq){
                cha_waveform = sharp_spike(sample_freq,freq)
                cha_waveform_length = length(cha_waveform)
                xzero = rep(0,total_length_tail)
                xnew = xzero
                laglength = as.integer((maxwaveformlength-cha_waveform_length)/2)
                xnew[(1:cha_waveform_length)+laglength] = cha_waveform
                xnew_save = xnew
                
                xstd_spike_tail[,c(paste("V",ncol(xstd_spike_tail)+1,sep="")) := xnew]
                
                for (loc in seq(total_length_tail-jump+1,maxwaveformlength+2,by = -jump)){
                  xnew = c(xnew_save[loc:(loc+jump-1)],xnew)[1:total_length_tail]
                  xstd_spike_tail[,c(paste("V",ncol(xstd_spike_tail)+1,sep="")) := xnew]
                }
              }
              for (freq in cand_spike_freq){
                cha_waveform = sharp_spike(sample_freq,freq)
                cha_waveform_length = length(cha_waveform)
                xzero = rep(0,total_length_tail)
                xnew = xzero
                laglength = as.integer((maxwaveformlength-cha_waveform_length)/2)
                xnew[(1:cha_waveform_length)+laglength] = cha_waveform
                xnew_save = xnew
                xstd_spike_tail[,c(paste("V",ncol(xstd_spike_tail)+1,sep="")) := xnew]
                
                for (loc in seq(total_length_tail-jump+1,maxwaveformlength+2,by = -jump)){
                  xnew = c(xnew_save[loc:(loc+jump-1)],xnew)[1:total_length_tail]
                  xstd_spike_tail[,c(paste("V",ncol(xstd_spike_tail)+1,sep="")) := xnew]
                }
              }
              xstd_spike_tail=as.matrix(xstd_spike_tail)
              xstd_spike_tail = xstd_spike_tail[,-1]
              xstd_spike_tailmask = matrix(NA,nrow = nrow(xstd_spike_tail),ncol = ncol(xstd_spike_tail))
              xstd_spike_tailmask[xstd_spike_tail!=0] = 1L
              xstd_spike_tail = as(xstd_spike_tail,'dgCMatrix')
              
              
              xstd_spislo_tail <- matrnorm(total_length_tail, 1)
              xstd_spislo_tail <- data.table(xstd_spislo_tail)
              
              for (freq1 in cand_spislo_sharp_freq){
                for (freq2 in cand_wave_freq){
                  cha_waveform = spike_and_wave(sample_freq,freq1,freq2)
                  cha_waveform_length = length(cha_waveform)
                  xzero = rep(0,total_length_tail)
                  xnew = xzero
                  laglength = as.integer((maxwaveformlength-cha_waveform_length)/2)
                  xnew[(1:cha_waveform_length)+laglength] = cha_waveform
                  xnew_save = xnew
                  
                  xstd_spislo_tail[,c(paste("V",ncol(xstd_spislo_tail)+1,sep="")) := xnew]
                  
                  for (loc in seq(total_length_tail-jump+1,maxwaveformlength+2,by = -jump)){
                    xnew = c(xnew_save[loc:(loc+jump-1)],xnew)[1:total_length_tail]
                    xstd_spislo_tail[,c(paste("V",ncol(xstd_spislo_tail)+1,sep="")) := xnew]
                  }
                }
              }
              xstd_spislo_tail=as.matrix(xstd_spislo_tail)
              xstd_spislo_tail = xstd_spislo_tail[,-1]
              xstd_spislo_tailmask = matrix(NA,nrow = nrow(xstd_spislo_tail),ncol = ncol(xstd_spislo_tail))
              xstd_spislo_tailmask[xstd_spislo_tail!=0] = 1L
              xstd_spislo_tail = as(xstd_spislo_tail,'dgCMatrix')
              
              xfit_spike_tail = as.matrix(cbind(xstd_spike_tail, 1))
              xfit_spislo_tail = as.matrix(cbind(xstd_spislo_tail, 1))
              
              
              OMPresult_spike <-ompr2( y = EEG_sig, x = xstd_spike_tail, xfit = xfit_spike_tail, 
                                       xstdmask = xstd_spike_tailmask,
                                       stopping_rule_var = stopping_rule_var, xstand = FALSE,
                                       ystand = FALSE,method = 'l2_bounded_noise',
                                       basisnum = total_num_spike,specificnum = wave_num_spike,
                                       sample_freq = sample_freq,
                                       freq_prop = 1 /10,jump = jump, halfwaveformlength = halfwaveformlength)
              OMPresult_spislo <-ompr2( y = EEG_sig, x = xstd_spislo_tail, xfit = xfit_spislo_tail, 
                                        xstdmask = xstd_spislo_tailmask,
                                        stopping_rule_var = stopping_rule_var, xstand = FALSE,
                                        ystand = FALSE,method = 'l2_bounded_noise',
                                        basisnum = total_num_spislo,specificnum = wave_num_spislo,
                                        sample_freq = sample_freq,
                                        freq_prop = 1 /10,jump = jump, halfwaveformlength = halfwaveformlength)
              
              
              selaloc_spike = OMPresult_spike$info[, 1]
              selaloc_spislo = OMPresult_spislo$info[, 1]
              
              coefficient_mat = rep(0, dim(xstd_spike_tail)[2]+dim(xstd_spislo_tail)[2])
              coefficient_mat = 
                Matrix(coefficient_mat,ncol = (total_num_spike+total_num_spislo),
                       nrow = (dim(xstd_spike_tail)[2]+dim(xstd_spislo_tail)[2])/
                         (total_num_spike+total_num_spislo),
                       byrow = FALSE,sparse = T)
              corcoefficient_mat = coefficient_mat
              
              
              if (length(selaloc_spike) > 0) {
                selaregcoeff_spike = OMPresult_spike$info[, 3]
                selacorcoeff_spike = OMPresult_spike$info[, 4]
                
                coefficient_mat_spike = rep(0, dim(xstd_spike_tail)[2])
                coefficient_mat_spike[selaloc_spike] = selaregcoeff_spike
                coefficient_mat[,c(1:total_num_spike)] = 
                  Matrix(coefficient_mat_spike,ncol = (total_num_spike),
                         nrow = dim(xstd_spike_tail)[2]/total_num_spike,byrow = FALSE,sparse = T)
                corcoefficient_mat_spike = rep(0, dim(xstd_spike_tail)[2])
                corcoefficient_mat_spike[selaloc_spike] = selacorcoeff_spike
                corcoefficient_mat[,c(1:total_num_spike)] = 
                  Matrix(corcoefficient_mat_spike,ncol = (total_num_spike),
                         nrow = dim(xstd_spike_tail)[2] /total_num_spike,byrow = FALSE,sparse = T)
              }
              if (length(selaloc_spislo) > 0) {
                selaregcoeff_spislo = OMPresult_spislo$info[, 3]
                selacorcoeff_spislo = OMPresult_spislo$info[, 4]
                
                coefficient_mat_spislo = rep(0, dim(xstd_spislo_tail)[2])
                coefficient_mat_spislo[selaloc_spislo] = selaregcoeff_spislo
                coefficient_mat[,c((total_num_spike+1):ncol(coefficient_mat))] = 
                  Matrix(coefficient_mat_spislo,ncol = (total_num_spislo),
                         nrow = dim(xstd_spislo_tail)[2]/total_num_spislo,byrow = FALSE,sparse = T)
                corcoefficient_mat_spislo = rep(0, dim(xstd_spislo_tail)[2])
                corcoefficient_mat_spislo[selaloc_spislo] = selacorcoeff_spislo
                corcoefficient_mat[,c((total_num_spike+1):ncol(coefficient_mat))] = 
                  Matrix(corcoefficient_mat_spislo,ncol = (total_num_spislo),
                         nrow = dim(xstd_spislo_tail)[2] /total_num_spislo,byrow = FALSE,sparse = T)
              }
              rm(xfit_spike_tail,xfit_spislo_tail)
              
            }else{
              # print('standard procedure')
              EEG_sig = EEG_sig_full[(artstart):(artstart + total_length - 1)]
              OMPresult_spike <-ompr2( y = EEG_sig, x = xstd_spike, xfit = xfit_spike, 
                                       xstdmask = xstd_spikemask,
                                       stopping_rule_var = stopping_rule_var, xstand = FALSE,
                                       ystand = FALSE,method = 'l2_bounded_noise',
                                       basisnum = total_num_spike,specificnum = wave_num_spike,
                                       sample_freq = sample_freq,
                                       freq_prop = 1 /10,jump = jump, halfwaveformlength = halfwaveformlength)
              OMPresult_spislo <-ompr2( y = EEG_sig, x = xstd_spislo, xfit = xfit_spislo, 
                                        xstdmask = xstd_spislomask,
                                        stopping_rule_var = stopping_rule_var, xstand = FALSE,
                                        ystand = FALSE,method = 'l2_bounded_noise',
                                        basisnum = total_num_spislo,specificnum = wave_num_spislo,
                                        sample_freq = sample_freq,
                                        freq_prop = 1 /10,jump = jump, halfwaveformlength = halfwaveformlength)
              
              selaloc_spike = OMPresult_spike$info[, 1]
              selaloc_spislo = OMPresult_spislo$info[, 1]
              
              coefficient_mat = rep(0, (total_num_spike+total_num_spislo) * active_length / jump)
              coefficient_mat = 
                Matrix(coefficient_mat,ncol = (total_num_spike+total_num_spislo),
                       nrow = active_length/jump,byrow = FALSE,sparse = T)
              corcoefficient_mat = coefficient_mat
              
              if (length(selaloc_spike) > 0) {
                selaregcoeff_spike = OMPresult_spike$info[, 3]
                selacorcoeff_spike = OMPresult_spike$info[, 4]
                
                coefficient_mat_spike = rep(0, (total_num_spike) *active_length / jump)
                coefficient_mat_spike[selaloc_spike] = selaregcoeff_spike
                coefficient_mat[,c(1:total_num_spike)] = 
                  Matrix(coefficient_mat_spike,ncol = (total_num_spike),
                         nrow = active_length/jump,byrow = FALSE,sparse = T)
                corcoefficient_mat_spike = rep(0, (total_num_spike) *active_length / jump)
                corcoefficient_mat_spike[selaloc_spike] = selacorcoeff_spike
                corcoefficient_mat[,c(1:total_num_spike)] = 
                  Matrix(corcoefficient_mat_spike,ncol = (total_num_spike),
                         nrow = active_length /jump,byrow = FALSE,sparse = T)
              }
              if (length(selaloc_spislo) > 0) {
                selaregcoeff_spislo = OMPresult_spislo$info[, 3]
                selacorcoeff_spislo = OMPresult_spislo$info[, 4]
                
                coefficient_mat_spislo = rep(0, (total_num_spislo) *active_length / jump)
                coefficient_mat_spislo[selaloc_spislo] = selaregcoeff_spislo
                coefficient_mat[,c((total_num_spike+1):ncol(coefficient_mat))] = 
                  Matrix(coefficient_mat_spislo,ncol = (total_num_spislo),
                         nrow = active_length/jump,byrow = FALSE,sparse = T)
                corcoefficient_mat_spislo = rep(0, (total_num_spislo) *active_length / jump)
                corcoefficient_mat_spislo[selaloc_spislo] = selacorcoeff_spislo
                corcoefficient_mat[,c((total_num_spike+1):ncol(coefficient_mat))] = 
                  Matrix(corcoefficient_mat_spislo,ncol = (total_num_spislo),
                         nrow = active_length /jump,byrow = FALSE,sparse = T)
              }
            } 
            coefficient_mat_list[[artstart]] = coefficient_mat
            corcoefficient_mat_list[[artstart]] = corcoefficient_mat
            artstart = artstart + active_length
          }
          coefficient_mat_full_temp = do.call('rbind', coefficient_mat_list)
          corcoefficient_mat_full_temp = do.call('rbind', corcoefficient_mat_list)
          
          gc()
          
          list(coefficient_mat_full_temp, corcoefficient_mat_full_temp)
          
        }
      stopCluster(cl)
      stopImplicitCluster()
      save(waveform_allchannel, 
           file=paste0(tempdir,csvfile,"_waveform_allchannel",".RData"))
    }else{
      waveform_allchannel = list()
      load(paste0(tempdir,csvfile,"_waveform_allchannel",".RData"))
    }
    
    if (!file.exists(paste0(tempdir,csvfile,"_waveform_bipolar",".RData"))){
      waveform_bipolar = list()
      EEGdata_bipolar = EEGdata[,c(8,9,10,12,13,14,16,17,18,20,21,22,24)]
      channel_name_bipolar = c("F7","T7","P7","F3","C3",  
                               "P3","F4","C4","P4","F8", 
                               "T8","P8","CZ")
      colnames(EEGdata_bipolar) = channel_name_bipolar
      EEGdata_bipolar2 = EEGdata[,c(9,10,11,13,14,15,17,18,19,21,22,23,25)]
      colnames(EEGdata_bipolar2) = channel_name_bipolar
      EEGdata_bipolar = rbind(EEGdata_bipolar, EEGdata_bipolar2)
      rm(EEGdata_bipolar2)
      
      cl=makeCluster(3)
      registerDoParallel(cl)
      waveform_bipolar = 
        foreach(EEG_sig_full = EEGdata_bipolar,.inorder = TRUE,
                .packages = c("signal", "Rfast","RcppArmadillo","Matrix", 'Rcpp',
                              'RcppEigen','data.table'),
                .final = function(x)setNames(x, channel_name_bipolar)) %dopar% 
        {
          sourceCpp(paste0(codedir,'sparse_eachcol_apply.cpp'))
          halflength = length(EEG_sig_full)/2
          EEG_sig_full_1 = EEG_sig_full[1:halflength]
          EEG_sig_full_2 = EEG_sig_full[(halflength+1):length(EEG_sig_full)]
          stopping_rule_var = 0
          artstart = 1
          
          coefficient_mat_list = list()
          corcoefficient_mat_list = list()
          xstd_spike_bipolar = rbind(-xstd_spike,xstd_spike)
          xstd_spikemask_bipolar = rbind(xstd_spikemask,xstd_spikemask)
          xfit_spike_bipolar = as.matrix(cbind(xstd_spike_bipolar, 1))
          
          while (artstart + total_length - 1 <= halflength) {
            print(artstart)
            if (artstart + total_length - 1 > halflength -active_length) {
              #tail procedure
              # print('tail procedure')
              EEG_sig1 = EEG_sig_full_1[(artstart):halflength]
              total_length_tail = length(EEG_sig1)
              EEG_sig2 = EEG_sig_full_2[(artstart):halflength]
              EEG_sig = c(EEG_sig1,EEG_sig2)
              
              xstd_spike_tail<- matrnorm(total_length_tail, 1)
              xstd_spike_tail<- data.table(xstd_spike_tail)
              for (freq in cand_sharp_freq){
                cha_waveform = sharp_spike(sample_freq,freq)
                cha_waveform_length = length(cha_waveform)
                xzero = rep(0,total_length_tail)
                xnew = xzero
                laglength = as.integer((maxwaveformlength-cha_waveform_length)/2)
                xnew[(1:cha_waveform_length)+laglength] = cha_waveform
                xnew_save = xnew
                
                xstd_spike_tail[,c(paste("V",ncol(xstd_spike_tail)+1,sep="")) := xnew]
                
                for (loc in seq(total_length_tail-jump+1,maxwaveformlength+2,by = -jump)){
                  xnew = c(xnew_save[loc:(loc+jump-1)],xnew)[1:total_length_tail]
                  xstd_spike_tail[,c(paste("V",ncol(xstd_spike_tail)+1,sep="")) := xnew]
                }
              }
              for (freq in cand_spike_freq){
                #total_length_tail = 2000
                cha_waveform = sharp_spike(sample_freq,freq)
                cha_waveform_length = length(cha_waveform)
                xzero = rep(0,total_length_tail)
                xnew = xzero
                laglength = as.integer((maxwaveformlength-cha_waveform_length)/2)
                xnew[(1:cha_waveform_length)+laglength] = cha_waveform
                xnew_save = xnew
                xstd_spike_tail[,c(paste("V",ncol(xstd_spike_tail)+1,sep="")) := xnew]
                
                for (loc in seq(total_length_tail-jump+1,maxwaveformlength+2,by = -jump)){
                  xnew = c(xnew_save[loc:(loc+jump-1)],xnew)[1:total_length_tail]
                  xstd_spike_tail[,c(paste("V",ncol(xstd_spike_tail)+1,sep="")) := xnew]
                }
              }
              xstd_spike_tail=as.matrix(xstd_spike_tail)
              xstd_spike_tail = xstd_spike_tail[,-1]
              xstd_spike_tailmask = matrix(NA,nrow = nrow(xstd_spike_tail),ncol = ncol(xstd_spike_tail))
              xstd_spike_tailmask[xstd_spike_tail!=0] = 1L
              xstd_spike_tail = as(xstd_spike_tail,'dgCMatrix')
              
              xstd_spike_tail_bipolar = rbind(-xstd_spike_tail,xstd_spike_tail)
              xstd_spike_tailmask_bipolar = rbind(xstd_spike_tailmask,xstd_spike_tailmask)
              xfit_spike_tail_bipolar = as.matrix(cbind(xstd_spike_tail_bipolar, 1))
              
              
              OMPresult_spike <-ompr2( y = EEG_sig, x = xstd_spike_tail_bipolar, 
                                       xfit = xfit_spike_tail_bipolar, 
                                       xstdmask = xstd_spike_tailmask_bipolar,
                                       stopping_rule_var = stopping_rule_var, xstand = FALSE,
                                       ystand = FALSE,method = 'l2_bounded_noise',
                                       basisnum = total_num_spike,specificnum = wave_num_spike,
                                       sample_freq = sample_freq,
                                       freq_prop = 1 /10,jump = jump, 
                                       halfwaveformlength = halfwaveformlength)
              
              selaloc_spike = OMPresult_spike$info[, 1]
              
              coefficient_mat = rep(0, dim(xstd_spike_tail)[2])
              coefficient_mat = 
                Matrix(coefficient_mat,ncol = (total_num_spike),
                       nrow = (dim(xstd_spike_tail)[2])/(total_num_spike),
                       byrow = FALSE,sparse = T)
              corcoefficient_mat = coefficient_mat
              
              
              if (length(selaloc_spike) > 0) {
                selaregcoeff_spike = OMPresult_spike$info[, 3]
                selacorcoeff_spike = OMPresult_spike$info[, 4]
                
                coefficient_mat_spike = rep(0, dim(xstd_spike_tail)[2])
                coefficient_mat_spike[selaloc_spike] = selaregcoeff_spike
                coefficient_mat[,c(1:total_num_spike)] = 
                  Matrix(coefficient_mat_spike,ncol = (total_num_spike),
                         nrow = dim(xstd_spike_tail)[2]/total_num_spike,byrow = FALSE,sparse = T)
                corcoefficient_mat_spike = rep(0, dim(xstd_spike_tail)[2])
                corcoefficient_mat_spike[selaloc_spike] = selacorcoeff_spike
                corcoefficient_mat[,c(1:total_num_spike)] = 
                  Matrix(corcoefficient_mat_spike,ncol = (total_num_spike),
                         nrow = dim(xstd_spike_tail)[2] /total_num_spike,byrow = FALSE,sparse = T)
              }
              rm(xfit_spike_tail_bipolar)
              
            }else{
              # print('standard procedure')
              EEG_sig1 = EEG_sig_full_1[(artstart):(artstart + total_length - 1)]
              EEG_sig2 = EEG_sig_full_2[(artstart):(artstart + total_length - 1)]
              EEG_sig = c(EEG_sig1, EEG_sig2)
              OMPresult_spike <-ompr2( y = EEG_sig, x = xstd_spike_bipolar, xfit = xfit_spike_bipolar, 
                                       xstdmask = xstd_spikemask_bipolar,
                                       stopping_rule_var = stopping_rule_var, xstand = FALSE,
                                       ystand = FALSE,method = 'l2_bounded_noise',
                                       basisnum = total_num_spike,specificnum = wave_num_spike,
                                       sample_freq = sample_freq,
                                       freq_prop = 1 /10,jump = jump, 
                                       halfwaveformlength = halfwaveformlength)
              
              selaloc_spike = OMPresult_spike$info[, 1]
              
              coefficient_mat = rep(0, (total_num_spike) * active_length / jump)
              coefficient_mat = 
                Matrix(coefficient_mat,ncol = (total_num_spike),
                       nrow = active_length/jump,byrow = FALSE,sparse = T)
              corcoefficient_mat = coefficient_mat
              
              if (length(selaloc_spike) > 0) {
                selaregcoeff_spike = OMPresult_spike$info[, 3]
                selacorcoeff_spike = OMPresult_spike$info[, 4]
                
                coefficient_mat_spike = rep(0, (total_num_spike) *active_length / jump)
                coefficient_mat_spike[selaloc_spike] = selaregcoeff_spike
                coefficient_mat[,c(1:total_num_spike)] = 
                  Matrix(coefficient_mat_spike,ncol = (total_num_spike),
                         nrow = active_length/jump,byrow = FALSE,sparse = T)
                corcoefficient_mat_spike = rep(0, (total_num_spike) *active_length / jump)
                corcoefficient_mat_spike[selaloc_spike] = selacorcoeff_spike
                corcoefficient_mat[,c(1:total_num_spike)] = 
                  Matrix(corcoefficient_mat_spike,ncol = (total_num_spike),
                         nrow = active_length /jump,byrow = FALSE,sparse = T)
              }
            } 
            coefficient_mat_list[[artstart]] = coefficient_mat
            corcoefficient_mat_list[[artstart]] = corcoefficient_mat
            artstart = artstart + active_length
          }
          coefficient_mat_full_temp = do.call('rbind', coefficient_mat_list)
          corcoefficient_mat_full_temp = do.call('rbind', corcoefficient_mat_list)
          
          gc()
          
          list(coefficient_mat_full_temp, corcoefficient_mat_full_temp)
          
        }
      stopCluster(cl)
      stopImplicitCluster()
      save(waveform_bipolar, 
           file=paste0(tempdir,csvfile,"_waveform_bipolar",".RData"))
      rm(EEGdata_bipolar)
    }else{
      waveform_bipolar = list()
      load(paste0(tempdir,csvfile,"_waveform_bipolar",".RData"))
    }
    
    for (j in c(1:length(waveform_allchannel))){
      #print(j)
      strength_mat = waveform_allchannel[[j]][[1]]
      strength_mat_new = as.matrix(strength_mat)
      cor_mat = waveform_allchannel[[j]][[2]]
      cor_mat_new = as.matrix(cor_mat)
      i = 0
      for (freq in cand_sharp_freq){
        cha_waveform = sharp_spike(sample_freq,freq)
        cha_waveform_length = as.integer(length(cha_waveform)/jump)
        cha_waveform_length_half = as.integer(cha_waveform_length/2)
        i = i+1
        
        strength_mat_new = wave_feature_combine(strength_mat,strength_mat_new,i)
        cor_mat_new = wave_feature_combine(cor_mat,cor_mat_new,i)
        
      }
      for (freq in cand_spike_freq){
        cha_waveform = sharp_spike(sample_freq,freq)
        cha_waveform_length = as.integer(length(cha_waveform)/jump)
        cha_waveform_length_half = as.integer(cha_waveform_length/2)
        i = i+1
        strength_mat_new = wave_feature_combine(strength_mat,strength_mat_new,i)
        cor_mat_new = wave_feature_combine(cor_mat,cor_mat_new,i)
      }
      for (freq1 in cand_spislo_sharp_freq){
        for (freq2 in cand_wave_freq){
          cha_waveform = spike_and_wave(sample_freq,freq1,freq2)
          cha_waveform_length = as.integer(length(cha_waveform)/jump)
          cha_waveform_length_half = as.integer(cha_waveform_length/2)
          i = i+1
          strength_mat_new = wave_feature_combine(strength_mat,strength_mat_new,i)
          cor_mat_new = wave_feature_combine(cor_mat,cor_mat_new,i)
        }
      }
      strength_mat_new = as(strength_mat_new, 'dgCMatrix')
      cor_mat_new = as(cor_mat_new, 'dgCMatrix')
      waveform_allchannel[[j]][[1]] = strength_mat_new
      waveform_allchannel[[j]][[2]] = cor_mat_new
    }
    
    
    for (j in c(1:length(waveform_bipolar))){
      #print(j)
      strength_mat = waveform_bipolar[[j]][[1]]
      strength_mat_new = as.matrix(strength_mat)
      cor_mat = waveform_bipolar[[j]][[2]]
      cor_mat_new = as.matrix(cor_mat)
      i = 0
      for (freq in cand_sharp_freq){
        cha_waveform = sharp_spike(sample_freq,freq)
        cha_waveform_length = as.integer(length(cha_waveform)/jump)
        cha_waveform_length_half = as.integer(cha_waveform_length/2)
        i = i+1
        strength_mat_new = wave_feature_combine(strength_mat,strength_mat_new,i)
        cor_mat_new = wave_feature_combine(cor_mat,cor_mat_new,i)
      }
      for (freq in cand_spike_freq){
        cha_waveform = sharp_spike(sample_freq,freq)
        cha_waveform_length = as.integer(length(cha_waveform)/jump)
        cha_waveform_length_half = as.integer(cha_waveform_length/2)
        i = i+1
        strength_mat_new = wave_feature_combine(strength_mat,strength_mat_new,i)
        cor_mat_new = wave_feature_combine(cor_mat,cor_mat_new,i)
      }
      strength_mat_new = as(strength_mat_new, 'dgCMatrix')
      cor_mat_new = as(cor_mat_new, 'dgCMatrix')
      waveform_bipolar[[j]][[1]] = strength_mat_new
      waveform_bipolar[[j]][[2]] = cor_mat_new
    }
    
    #from channel level to brain regions
    for (i in c("FP1-F7","FP1-F3","F3-C3")){
      EEGdata_waveform_rawcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_rawcoeff,halfwaveformlength,1,i,1)
      EEGdata_waveform_corcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_corcoeff,halfwaveformlength,2,i,1)
    }
    for (i in c("FP2-F4","F4-C4","FP2-F8")){
      EEGdata_waveform_rawcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_rawcoeff,halfwaveformlength,1,i,2)
      EEGdata_waveform_corcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_corcoeff,halfwaveformlength,2,i,2)
    }
    for (i in c("FP1-F7","F7-T7","T7-P7","P7-O1","T7-FT9")){
      EEGdata_waveform_rawcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_rawcoeff,halfwaveformlength,1,i,3)
      EEGdata_waveform_corcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_corcoeff,halfwaveformlength,2,i,3)
    }
    for (i in c("FP2-F8","F8-T8","T8-P8","P8-O2","FT10-T8")){
      EEGdata_waveform_rawcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_rawcoeff,halfwaveformlength,1,i,4)
      EEGdata_waveform_corcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_corcoeff,halfwaveformlength,2,i,4)
    }
    for (i in c("FZ-CZ","CZ-PZ")){
      EEGdata_waveform_rawcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_rawcoeff,halfwaveformlength,1,i,5)
      EEGdata_waveform_corcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_corcoeff,halfwaveformlength,2,i,5)
    }
    for (i in c("F3-C3","C3-P3","P3-O1")){
      EEGdata_waveform_rawcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_rawcoeff,halfwaveformlength,1,i,6)
      EEGdata_waveform_corcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_corcoeff,halfwaveformlength,2,i,6)
    }
    for (i in c("F4-C4","C4-P4","P4-O2")){
      EEGdata_waveform_rawcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_rawcoeff,halfwaveformlength,1,i,7)
      EEGdata_waveform_corcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_corcoeff,halfwaveformlength,2,i,7)
    }
    for (i in c("P7-O1","P3-O1","P4-O2","P8-O2")){
      EEGdata_waveform_rawcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_rawcoeff,halfwaveformlength,1,i,8)
      EEGdata_waveform_corcoeff = 
        wave_feature_rep(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                         EEGdata_waveform_corcoeff,halfwaveformlength,2,i,8)
    }
    
    for (i in c("F3")){
      spikesig = rowSums(abs(waveform_bipolar[[i]][[1]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),1] = 
        EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),1]+spikesig
      
      spikesig = rowSums(abs(waveform_bipolar[[i]][[2]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),1] = 
        EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),1]+spikesig
    }
    for (i in c("F4")){
      spikesig = rowSums(abs(waveform_bipolar[[i]][[1]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),2] = 
        EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),2]+spikesig
      
      spikesig = rowSums(abs(waveform_bipolar[[i]][[2]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),2] = 
        EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),2]+spikesig
    }
    for (i in c("F7","T7","P7")){
      spikesig = rowSums(abs(waveform_bipolar[[i]][[1]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),3] = 
        EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),3]+spikesig
      
      spikesig = rowSums(abs(waveform_bipolar[[i]][[2]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),3] = 
        EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),3]+spikesig
    }
    for (i in c("F8","T8","P8")){
      spikesig = rowSums(abs(waveform_bipolar[[i]][[1]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),4] = 
        EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),4]+spikesig
      
      spikesig = rowSums(abs(waveform_bipolar[[i]][[2]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),4] = 
        EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),4]+spikesig
    }
    for (i in c("CZ")){
      spikesig = rowSums(abs(waveform_bipolar[[i]][[1]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),5] = 
        EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),5]+spikesig
      
      spikesig = rowSums(abs(waveform_bipolar[[i]][[2]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),5] = 
        EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),5]+spikesig
    }
    for (i in c("C3","P3")){
      spikesig = rowSums(abs(waveform_bipolar[[i]][[1]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),6] = 
        EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),6]+spikesig
      
      spikesig = rowSums(abs(waveform_bipolar[[i]][[2]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),6] = 
        EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),6]+spikesig
    }
    for (i in c("C4","P4")){
      spikesig = rowSums(abs(waveform_bipolar[[i]][[1]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),7] = 
        EEGdata_waveform_bipolar_rawcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),7]+spikesig
      
      spikesig = rowSums(abs(waveform_bipolar[[i]][[2]][,1:wave_num_spike]))
      spikesig = rep(spikesig,each = jump)
      EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),7] = 
        EEGdata_waveform_bipolar_corcoeff[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),7]+spikesig
    }
    
    #derive the local variance and variance ratio
    piecelength = non_seizure_segment_length*sample_freq
    artstart = 1
    
    final_cp_no_cpthr = c(0)
    
    #split EEGdata
    while ((artstart + piecelength - 1)<=nrow(EEGdata)){
      
      if ((artstart + 1.5*piecelength - 1)<=nrow(EEGdata)){
        art_sig = EEGdata[c(artstart:(artstart+piecelength-1)),1:28]
        EEGlabel = EEGdata[c(artstart:(artstart+piecelength-1)),c(2:4,6)]
        EEGlabel = unique(na.omit(EEGlabel))
        if ((artstart-1)%%piecelength==0){
          cal_wave = TRUE
        }else{
          cal_wave = FALSE
        }
      }else{
        art_sig = EEGdata[c(artstart:nrow(EEGdata)),1:28]
        EEGlabel = EEGdata[c(artstart:nrow(EEGdata)),c(2:4,6)]
        EEGlabel = unique(na.omit(EEGlabel))
        cal_wave = TRUE
      }
      
      datamat = art_sig[,c(8:28)]
      
      SCICAresult = SCICA(datamat,ncol(datamat), alg.typ = "parallel",fun = 'logcosh', alpha = 1, 
                          row.norm = FALSE, maxit = 200, sc.typ = 'soft',
                          tol = 0.0001, verbose = FALSE, SC = refrow)
      Smat = SCICAresult$S
      
      seizure_sig_1d = rowMeans(Smat[,1:ICnum,drop = FALSE]%*%
                                  SCICAresult$A[1:ICnum,,drop = FALSE])
      seizure_sig_1d = scale(seizure_sig_1d, scale = FALSE)
      sigmahat_kernel_DPI = local_var(seizure_sig_1d)
      
      seizure_raw_1d = rowMeans(art_sig[,8:28])
      rawsigmahat_kernel_DPI = local_var(seizure_raw_1d)
      
      ratiosigma = sigmahat_kernel_DPI/rawsigmahat_kernel_DPI
      
      #21 dimensional seizure signal
      newW = SCICAresult$K%*%SCICAresult$W
      newW = newW[,1:ICnum,drop =FALSE]
      for (i in 1:ncol(newW)){
        newW[,i] = newW[,i]/sqrt(sum(newW[,i]^2))
      }
      seizure_sig = as.matrix(datamat)%*%newW
      
      #save the SCICA results
      EEGdata_SCICA_1d[(piecelength/4+artstart-1):(nrow(art_sig)-piecelength/4+artstart-1)
      ] = seizure_sig_1d[(piecelength/4):(nrow(art_sig)-piecelength/4)]
      EEGdata_SCICA_sigmahat[(piecelength/4+artstart-1):(nrow(art_sig)-piecelength/4+artstart-1)
      ] = sigmahat_kernel_DPI[(piecelength/4):(nrow(art_sig)-piecelength/4)]
      EEGdata_SCICA_sigmahat_ratio[(piecelength/4+artstart-1):(nrow(art_sig)-piecelength/4+artstart-1)
      ] = ratiosigma[(piecelength/4):(nrow(art_sig)-piecelength/4)]
      
      print(artstart)
      #change point detection
      sbsresult_5 = sbs.alg(t(seizure_sig), cp.type=2, scales=-keylevel, trim = 64, thr = thr,
                            diag=diag, do.parallel=10, .verbose = TRUE, sq = FALSE)
      ecp_5 = sbsresult_5$ecp
      ecp_5 = ecp_5[ecp_5>=(piecelength/4)&ecp_5<(nrow(art_sig)-piecelength/4)]
      
      final_cp_no_cpthr = c(final_cp_no_cpthr,ecp_5+artstart-1)
      
      artstart = artstart + 0.5*piecelength
    }
    
    #change point detection is unstable near the segment boundary, thus we drop
    #out the start and end piecelength/4 parts
    final_cp_fixed_segment = seq(piecelength/4,nrow(EEGdata)-piecelength/4,sample_freq*2) #2 seconds length segments
    training_time = intersect(EEGdata$index,EEGdata_train_index$index)
    if (length(training_time)>0){
      training_time = training_time - EEGdata$index[1] + 1
      training_time = sort(unique(training_time))
    }else{
      training_time = c()
    }
    
    #extract the features
    cp = as.numeric(final_cp_no_cpthr)
    cp = cp[cp>(piecelength/4)&cp<(nrow(EEGdata)-piecelength/4)]
    cp = c((piecelength/4),cp,(nrow(EEGdata)-piecelength/4))
    
    #remove those cps too close of each other
    j = 1
    while (j < length(cp)){
      if ((cp[j+1]-cp[j])<(sample_freq/2-1)){
        print(paste0('remove ',cp[j+1]))
        cp = cp[-(j+1)]
      }else{
        j = j+1
      }
    }
    
    for (j in c(1:(length(cp)-1))){
      variancedata = EEGdata_SCICA_sigmahat[c(cp[j]:cp[j+1])]
      
      variance = mean(variancedata)
      max_variance = max(variancedata)
      
      variance_ratio_data = EEGdata_SCICA_sigmahat_ratio[c(cp[j]:cp[j+1])]
      
      variance_ratio = mean(variance_ratio_data)
      max_variance_ratio = max(variance_ratio_data)
      
      waveformfeature = c(
        colMeans(EEGdata_waveform_rawcoeff[c(cp[j]:cp[j+1]),]),
        colMeans(EEGdata_waveform_bipolar_rawcoeff[c(cp[j]:cp[j+1]),]),
        colMeans(EEGdata_waveform_corcoeff[c(cp[j]:cp[j+1]),]),
        colMeans(EEGdata_waveform_bipolar_corcoeff[c(cp[j]:cp[j+1]),])
      )
      
      testlabel = EEGdata$label1[c(cp[j]:cp[j+1])]
      onenum = length(which(testlabel==1))
      allnum = length(testlabel)
      startj = EEGdata$index[cp[j]]
      endj = EEGdata$index[cp[j+1]]
      len = cp[j+1]-cp[j]+1
      usedfortraining = 0
      
      if (length(training_time)>0){
        if (length(intersect(c(cp[j]:cp[j+1]),training_time))>0){
          usedfortraining = 1
        }
      }
      
      EEG_feature = 
        rbind(EEG_feature,
              as.numeric(c(onenum/allnum,usedfortraining,startj,endj,len,variance,
                           max_variance,variance_ratio,max_variance_ratio,
                           waveformfeature)))
      
    }
    
    #extract the features
    cp = as.numeric(final_cp_fixed_segment)
    cp = cp[cp>(piecelength/4)&cp<(nrow(EEGdata)-piecelength/4)]
    cp = c((piecelength/4),cp,(nrow(EEGdata)-piecelength/4))
    
    #remove those cps too close of each other
    j = 1
    while (j < length(cp)){
      if ((cp[j+1]-cp[j])<(sample_freq/2-1)){
        print(paste0('remove ',cp[j+1]))
        cp = cp[-(j+1)]
      }else{
        j = j+1
      }
    }
    
    for (j in c(1:(length(cp)-1))){
      variancedata = EEGdata_SCICA_sigmahat[c(cp[j]:cp[j+1])]
      
      variance = mean(variancedata)
      max_variance = max(variancedata)
      
      variance_ratio_data = EEGdata_SCICA_sigmahat_ratio[c(cp[j]:cp[j+1])]
      
      variance_ratio = mean(variance_ratio_data)
      max_variance_ratio = max(variance_ratio_data)
      
      waveformfeature = c(
        colMeans(EEGdata_waveform_rawcoeff[c(cp[j]:cp[j+1]),]),
        colMeans(EEGdata_waveform_bipolar_rawcoeff[c(cp[j]:cp[j+1]),]),
        colMeans(EEGdata_waveform_corcoeff[c(cp[j]:cp[j+1]),]),
        colMeans(EEGdata_waveform_bipolar_corcoeff[c(cp[j]:cp[j+1]),])
      )
      
      testlabel = EEGdata$label1[c(cp[j]:cp[j+1])]
      onenum = length(which(testlabel==1))
      allnum = length(testlabel)
      startj = EEGdata$index[cp[j]]
      endj = EEGdata$index[cp[j+1]]
      len = cp[j+1]-cp[j]+1
      usedfortraining = 0
      
      if (length(training_time)>0){
        if (length(intersect(c(cp[j]:cp[j+1]),training_time))>0){
          usedfortraining = 1
        }
      }
      EEG_feature_fixed_segment = 
        rbind(EEG_feature_fixed_segment,
              as.numeric(c(onenum/allnum,usedfortraining,startj,endj,len,variance,
                           max_variance,variance_ratio,max_variance_ratio,
                           waveformfeature)))
    }
    
    print(csvfile)
    gc()
    
    list(EEG_feature[-1,],EEG_feature_fixed_segment[-1,],
         EEGdata[(piecelength/4):(nrow(EEGdata)-piecelength/4),c('index','label1')])
  }
stopCluster(cl)
stopImplicitCluster()