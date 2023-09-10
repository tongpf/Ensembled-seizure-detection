wave_feature_combine = function(raw_wave_mat,new_wave_mat,idx){
  wave1 = raw_wave_mat[,idx]
  keyparaloc = which(wave1!=0)
  wave2 = Matrix(0,nrow = length(wave1),ncol = length(keyparaloc))
  keyparaloc_2 = rep(c(-cha_waveform_length_half:(cha_waveform_length-cha_waveform_length_half-1)),
                     each = length(keyparaloc))
  keyparaloc_3 = rep(keyparaloc,times = cha_waveform_length)
  keyparaloc_2 = keyparaloc_2 + keyparaloc_3
  keyparaloc_3 = rep(c(1:length(keyparaloc)),times = cha_waveform_length)
  keyparaloc_4 = which(keyparaloc_2>=1 & keyparaloc_2<= length(wave1))
  keyparaloc_final = cbind(keyparaloc_2[keyparaloc_4],keyparaloc_3[keyparaloc_4])
  wave2[keyparaloc_final] = 1
  wave2 = t(wave2) * (wave1[keyparaloc]/cha_waveform_length)
  wave2 = colSums(wave2)
  new_wave_mat[,idx] = wave2
  return(new_wave_mat)
}

wave_feature_rep = function(waveform_allchannel,wave_num_spike,wave_num_spislo,jump,
                            EEGdata_waveform_feature,halfwaveformlength,feature_type,
                            channel_idx,region_index){
  #feature_type = 1, signal strength
  #feature_type = 2, signal correlation
  spikesig = rowSums(abs(waveform_allchannel[[channel_idx]][[feature_type]][,1:wave_num_spike]))
  spikesig = rep(spikesig,each = jump)
  spislosig = 
    rowSums(abs(waveform_allchannel[[channel_idx]][[feature_type]][,(wave_num_spike+1):
                                                                     (wave_num_spike+wave_num_spislo)]))
  spislosig = rep(spislosig,each = jump)
  EEGdata_waveform_feature[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),region_index] = 
    EEGdata_waveform_feature[c(halfwaveformlength:(length(spikesig)+halfwaveformlength-1)),region_index]+spikesig
  EEGdata_waveform_feature[c(halfwaveformlength:(length(spislosig)+halfwaveformlength-1)),region_index+8] = 
    EEGdata_waveform_feature[c(halfwaveformlength:(length(spislosig)+halfwaveformlength-1)),region_index+8]+spislosig
  
  return(EEGdata_waveform_feature)
}

local_var = function(key_signal,ifplot = FALSE,plotindex = NULL){
  xGrid = c(1:length(key_signal))
  hDPI <- KernSmooth::dpill(x = xGrid, y = key_signal)
  # Fits
  lp1DPI <- KernSmooth::locpoly(x = xGrid, y = key_signal, bandwidth = hDPI, 
                                degree = 0, range.x = c(min(xGrid),max(xGrid)), 
                                gridsize = length(xGrid))
  eps = (key_signal - lp1DPI$y)^2
  hDPI <- KernSmooth::dpill(x = xGrid, y = eps)
  lp1DPI <- KernSmooth::locpoly(x = xGrid, y = eps, bandwidth = hDPI, degree = 0,
                                range.x = c(min(xGrid),max(xGrid)), 
                                gridsize = length(xGrid))
  key_signal_local_var = lp1DPI$y
  if(ifplot){
    jpeg(paste0(plotdir,'start index ',plotindex,
                ' epilepsy signal.jpg'), width = 2000, height = 500)
    plot(key_signal,type = 'l')
    dev.off()
    jpeg(paste0(plotdir,'start index ',plotindex,
                ' epilepsy variance.jpg'), width = 2000, height = 500)
    plot(lp1DPI$x, lp1DPI$y,ylim = c(0, max(eps)))
    points(xGrid, eps)
    rug(xGrid, side = 1); rug(eps, side = 2)
    lines(xGrid, eps, col = 1)
    lines(lp1DPI$x, lp1DPI$y, col = 3, type = "o", pch = 16, cex = 0.5)
    dev.off()
  }
  return(key_signal_local_var)
}

center_adj = function(EEG_feature){
  indexdiff = EEG_feature$start[2:nrow(EEG_feature)]-
    EEG_feature$end[1:(nrow(EEG_feature)-1)]
  hour_st = c(1,which(indexdiff>0)+1)
  hour_en = c(which(indexdiff>0),nrow(EEG_feature))
  
  for (i in c(1:length(hour_st))){
    stidx = hour_st[i]
    enidx = hour_en[i]
    EEG_feature[stidx:enidx,c(6,8,10:ncol(EEG_feature))] = 
      scale(EEG_feature[stidx:enidx,c(6,8,10:ncol(EEG_feature))],
            center = apply(EEG_feature[stidx:enidx,c(6,8,10:ncol(EEG_feature))],2,median), 
            scale = F)
  }
  return(EEG_feature)
}