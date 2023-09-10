ompr2 = function (y, x, xfit, xstdmask, stopping_rule_var, ystand = TRUE, 
                  xstand = TRUE, method = "BIC", 
                  tol = 2, basisnum = 19,specificnum = 12,freq_prop = 1/10, 
                  jump = 4, sample_freq = 256, halfwaveformlength) {
  dm <- dim(x)
  d <- dm[2]
  nrowx = dm[1]
  #n <- dm[1]
  ind <- 1:d
  runtime <- proc.time()
  basiswidth = d/basisnum #number of candidate locations
  n = basiswidth*jump
  y <- (y - mean(y))
  #xfit = as.matrix(cbind(x,1))
  if (ystand) {
    m <- sum(y)/n
    y <- (y - m)/Rfast::Var(y, std = TRUE)
  }
  #if (xstand) {
  #  x <- Rfast::standardise(x)}
  if (method == "l2_bounded_noise") {#user specific
    con <- n * log(2 * pi) + n
    vary = Rfast::Var(y[(halfwaveformlength+1):(halfwaveformlength+n)]) * (n - 1)/n
    rho <- n * log(vary) + 2 * log(n)
    r <- sparse_eachcol_apply(x, y)
    save_first_result = r
    epe <- which(is.na(r))
    ind[epe] <- 0
    sel <- which.max(abs(r))
    sela <- sel
    #res <- .lm.fit(xfit[, c(sela,d+1), drop = FALSE], y)$residuals
    res = as.numeric(Rfast::lmfit(xfit[, c(sela,d+1), drop = FALSE], y)$residuals)

    varres = Rfast::Var(res[(halfwaveformlength+1):(halfwaveformlength+n)]) * (n - 1)/n
    rho[2] <- n * log(varres) + 3 * log(n)
    ind[sel] <- 0
    r[sel] <- 0
    i <- 2
    
    while (#rho[i - 1] - rho[i] > tol & 
      i <= 2 & 0.01*vary<varres) {
      i <- i + 1
      r[sela] <- 0
      r[ind] <- sparse_eachcol_apply(x, res)[ind>0]
      sel <- which.max(abs(r))
      sela <- c(sela, sel)
      #res <- .lm.fit(xfit[, c(sela,d+1)], y)$residuals
      res = as.numeric(Rfast::lmfit(xfit[, c(sela,d+1), drop = FALSE], y)$residuals)
      varres = Rfast::Var(res[(halfwaveformlength+1):(halfwaveformlength+n)]) * (n - 1)/n
      rho[i] <- n * log(varres) + (i + 1) * log(n)
      ind[sela] <- 0
    }
    len <- length(sela)
    info <- cbind(c(sela), rho[1:len] + con)
    colnames(info) <- c("Vars", "l2_bounded_noise")
  }else{
    stop('not supported methods')
  }

  regcoef = as.numeric(Rfast::lmfit(xfit[, c(sela,d+1), drop = FALSE], y)$be)[1:length(sela)]
  corcoef = colSums(x[,sela]*scale(xstdmask[,sela]*y),na.rm = TRUE)/sqrt(colSums(xstdmask[,sela],na.rm = TRUE))
  runtime <- proc.time() - runtime
  #runtime
  list(runtime = runtime, info = cbind(info,regcoef,corcoef),first_innerprod = save_first_result)
}

#waveform detection
plot.frequency.spectrum <- function(X.k, acq.freq, xlimits=c(0,acq.freq/2)) {
  N = length(X.k)
  plot.data  <- cbind(c(0:(length(X.k)-1))/length(X.k)*acq.freq, Mod(X.k))
  
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
  colnames(plot.data) = c('freq','spectrum')
  return(plot.data[1:as.integer(N/2),])
}

#sharp or spike
sharp_spike = function(acq.freq,target.freq){
  n = as.integer(acq.freq/target.freq)
  x = seq(0,1,length.out = n)
  x1 = x[x<=0.4]
  x2 = x[x>0.4]
  y1 = x1*15/4
  y2 = x2*(-5)/2+5/2
  y = c(y1,y2)
  y = y/target.freq
  y = scale(y)/sqrt(n)
  return(y)
}

sharp_spike_triple = function(acq.freq,target.freq){
  n = as.integer(acq.freq/target.freq/3)
  x = seq(0,1,length.out = n)
  x1 = x[x<=0.4]
  x2 = x[x>0.4]
  y1 = x1*15/4
  y2 = x2*(-5)/2+5/2
  y = c(y1,y2)
  y = y/target.freq
  y = rep(y,3)
  y = scale(y)/sqrt(n)
  return(y)
}

#spike-slow wave
spike_and_wave = function(acq.freq,target.freq1,target.freq2){
  n1 = as.integer(acq.freq/target.freq1)
  n2 = as.integer(acq.freq/target.freq2)
  xa = seq(0,1,length.out = n1)
  xb = seq(0,1,length.out = n2)
  y1 = -2*abs(xa - 1/2)+1
  y2 = -4*(xb-1/2)^2+1
  # y1 = y1/target.freq1
  # y2 = y2/target.freq2
  y=c(y1,y2)
  y = scale(y)/sqrt(n1+n2)
  return(y)
}

#triple spike slow wave
spike_and_wave_triple = function(acq.freq,target.freq1,target.freq2){
  n1 = as.integer(acq.freq/target.freq1/3)
  n2 = as.integer(acq.freq/target.freq2/3)
  xa = seq(0,1,length.out = n1)
  xb = seq(0,1,length.out = n2)
  y1 = -2*abs(xa - 1/2)+1
  y2 = -4*(xb-1/2)^2+1
  y=c(y1,y2)
  y = rep(y,3)
  y = scale(y)/sqrt(n1+n2)
  return(y)
}

#haar wave
haar_wave = function(acq.freq, target.freq){
  n = as.integer(acq.freq/target.freq)
  x = seq(0,1,length.out = n)
  x1 = x[x<=0.5]
  x2 = x[x>0.5]
  y1 = x1 - x1 +1
  y2 = x2 - x2 -1
  y=c(y1,y2)
  y=scale(y)/sqrt(n)
  return(y) 
}
