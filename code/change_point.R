# modified from the hdbinseg R package
# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @keywords internal
stl_sort <- function(x) {
  .Call('_hdbinseg_stl_sort', x)
}

#' @keywords internal
rcpp_rev <- function(x) {
  .Call('_hdbinseg_rcpp_rev', x)
}

#' @keywords internal
func_coef <- function(z, scale) {
  .Call('_hdbinseg_func_coef', z, scale)
}

#' @keywords internal
func_input <- function(coef, sgn, sq, diag) {
  .Call('_hdbinseg_func_input', coef, sgn, sq, diag)
}

#' @keywords internal
func_dc <- function(z, gamma) {
  .Call('_hdbinseg_func_dc', z, gamma)
}

#' @keywords internal
func_cusum <- function(z) {
  .Call('_hdbinseg_func_cusum', z)
}

#' @keywords internal
func_cusum_vec <- function(z) {
  .Call('_hdbinseg_func_cusum_vec', z)
}

#' @keywords internal
func_dc_by <- function(z, gamma, dmby, dtby) {
  .Call('_hdbinseg_func_dc_by', z, gamma, dmby, dtby)
}

#' @keywords internal
func_mvt_ar <- function(ar, res) {
  .Call('_hdbinseg_func_mvt_ar', ar, res)
}

#test example
#x <- matrix(rnorm(20*600), nrow=20)
#cov(x[1,],x[11,])
#x[1:10, 201:400] <- x[1:10, 201:400]*4
#cov(x[1,],x[11,])
#test = sbs.alg(x, cp.type=2, scales=-1, diag=FALSE, do.parallel=0, sq = FALSE)
#test$ecp
#test$change_dim_list

sbs.alg = function (x, cp.type = c(1, 2)[1], thr = NULL, trim = NULL, height = Inf, 
          temporal = TRUE, scales = NULL, diag = FALSE, B = 1000, q = 0.01, 
          do.parallel = 4, sq = FALSE, .verbose = FALSE) {
  n <- dim(x)[1] #number of variate
  T <- dim(x)[2] #number of observations
  stopifnot(n >= 1 && T >= 1)
  stopifnot(cp.type == 1 || cp.type == 2)
  stopifnot(is.null(trim) || (trim > 0 && 2 * trim + 1 < T))
  #stopifnot(is.null(height) || 2^height < T)
  if (cp.type == 2) {
    if (!is.null(scales)) { #can add multiple scales once time
      stopifnot(sum(scales < 0) == length(scales))
      scales <- sort(scales, decreasing = TRUE)
    }
    else{
      scales <- -1
    }
  }
  if (is.null(thr)) { #bootstrap procedure is employed for the threshold selection
    stopifnot(B > 0)
    stopifnot(q >= 0 && q <= 1)
  }else {
    if (cp.type == 1) 
      stopifnot(length(thr) == n)
    if (cp.type == 2) 
      stopifnot(length(thr) == (length(scales) * n * (n + 
                1)/2 * (!diag) + length(scales) * n * diag))
  }
  if (is.null(trim)) {#length of the intervals trimmed off around the change-point candidates
    trim <- round(log(T)^1.5)}
  if (is.null(height)) { #maximum height of the binary tree
    height <- floor(log(T, 2)/2)}
  if (cp.type == 1) {
    ccp <- clean.cp(x, type = "sbs", trim = trim, height = height)
    cx <- ccp$x
    if (!temporal) {
      tau <- apply(cx, 1, mad)
    }else {tau <- apply(cx, 1, mad)/apply(cx, 1, function(z) {
      1 - acf(z, lag.max = 1, type = "correlation", 
              plot = FALSE)$acf[2, , 1]
      })
    }
    if (is.null(thr)) {
      thr <- sbs.thr(cx, interval = c(1, T), cp.type = 1, 
                     B = B, q = q/n, do.parallel = do.parallel)}
    mt <- sbs.make.tree(x, tau, thr, trim, height)
  }
  if (cp.type == 2) {
    sgn <- sign(cor(t(x)))
    input <- gen.input(x, scales, sq, diag, sgn)
    if (is.null(thr)) {
      #each dimension of the periodograms and cross-periodograms has a unique threshold
      thr <- sbs.thr(x, interval = c(1, T), cp.type = 2, 
                     scales = scales, diag = diag, sgn = sgn, B = B, 
                     q = q, do.parallel = do.parallel, sq = sq, .verbose = .verbose)
    }
    mt <- sbs.make.tree(input, rep(1, dim(input)[1]), thr, 
                        trim, height, diag = diag, n = n)
  }
  return(mt)
}

sbs.thr = function (z, interval = c(1, dim(z)[2]), cp.type = 1, do.clean.cp = TRUE, 
                    scales = NULL, diag = FALSE, sgn = NULL, B = 1000, q = 0.01, 
                    do.parallel = 4,sq = FALSE,.verbose = FALSE) {
  
  #' @keywords internal
  stl_sort <- function(x) {
    .Call('_hdbinseg_stl_sort', x)
  }
  
  #' @keywords internal
  rcpp_rev <- function(x) {
    .Call('_hdbinseg_rcpp_rev', x)
  }
  
  #' @keywords internal
  func_coef <- function(z, scale) {
    .Call('_hdbinseg_func_coef', z, scale)
  }
  
  #' @keywords internal
  func_input <- function(coef, sgn, sq, diag) {
    .Call('_hdbinseg_func_input', coef, sgn, sq, diag)
  }
  
  #' @keywords internal
  func_dc <- function(z, gamma) {
    .Call('_hdbinseg_func_dc', z, gamma)
  }
  
  #' @keywords internal
  func_cusum <- function(z) {
    .Call('_hdbinseg_func_cusum', z)
  }
  
  #' @keywords internal
  func_cusum_vec <- function(z) {
    .Call('_hdbinseg_func_cusum_vec', z)
  }
  
  #' @keywords internal
  func_dc_by <- function(z, gamma, dmby, dtby) {
    .Call('_hdbinseg_func_dc_by', z, gamma, dmby, dtby)
  }
  
  #' @keywords internal
  func_mvt_ar <- function(ar, res) {
    .Call('_hdbinseg_func_mvt_ar', ar, res)
  }
  
  if (do.parallel > 0) {
    cl <- parallel::makeCluster(do.parallel)
    doParallel::registerDoParallel(cl)
  }
  #switch between dopar and do
  `%mydo%` <- ifelse(do.parallel > 0, `%dopar%`, `%do%`)
  n <- dim(z)[1]
  len <- dim(z)[2]
  s <- interval[1]
  e <- interval[2]
  if (e - s + 1 < ceiling(log(len)^1.5)) {
    stop("Error: the interval too short!")
  }
  if (cp.type == 1) {
    d <- n
  }
  if (cp.type == 2) {
    if (is.null(sgn)) {
      sgn <- sign(cor(t(z[, s:e])))
    }
    d <- n * diag + n * (n + 1)/2 * (!diag) #number of thresholds
  }
  if (cp.type == 1 && do.clean.cp) {
    z <- clean.cp(z, type = "sbs", trim = round(log(len)^1.5), 
                  height = round(log(len, 2)/2))$x
  }
  burnin <- 100
  #the maximum AR order of periodogram
  #In Cho and Fryzlewicz 2015, the AR is forced to be AR(1), but this code
  #seems to use fitted AR(p) to simulate threshold
  p <- floor(min(sqrt(e - s + 1), 10 * log(e - s + 1, 10)))
  ar.coef <- matrix(0, nrow = d, ncol = p)
  residuals <- matrix(0, nrow = d, ncol = e - s + 1 - p)
  k <- 1
  for (i in 1:n) {
    ar.fit <- stats::ar(z[i, s:e], order.max = p, method = "yw")
    if (length(ar.fit$ar) > 0) {
      ar.coef[k, 1:length(ar.fit$ar)] <- ar.fit$ar
    }
    residuals[k, ] <- ar.fit$resid[-(1:p)]
    k <- k + 1
    if (cp.type == 2 && i < n && !diag) {#for cross-periodogram
      for (j in (i + 1):n) {
        ar.fit <- stats::ar(z[i, s:e] - sgn[i, j] * z[j, 
                                                      s:e], order.max = p, method = "yw")
        if (length(ar.fit$ar) > 0) 
          ar.coef[k, 1:length(ar.fit$ar)] <- ar.fit$ar
        residuals[k, ] <- ar.fit$resid[-(1:p)]
        k <- k + 1
      }
    }
  }
  #why greater than 0? it should be not equal 0.
  if (sum(apply(abs(ar.coef), 2, max) > 0) > 0) {#there are some non-zero AR coefs
    ar.coef <- ar.coef[, 1:max(which(apply(abs(ar.coef), 
                                           2, max) > 0)), drop = FALSE]
    p <- ncol(ar.coef)
  }else {
    p <- 0
  }
  ns <- foreach::foreach(l = iterators::iter(1:B), .combine = cbind, 
                         .packages = c("Rcpp", "RcppArmadillo", "hdbinseg"),
                         .verbose = .verbose) %mydo% {
         #sample residuals with replacement. bootstrap for x?
         bx <- residuals[, sample(dim(residuals)[2], len + 
                                    p + burnin, replace = TRUE)]
         if (p > 0) {
           bx <- func_mvt_ar(ar.coef, bx)
         }
         #burnin may be used for stabilize the result
         bx <- bx[, -(1:(p + burnin))]
         if (cp.type == 1) {
           tmp <- apply(func_cusum(bx)$acs, 1, max)
         }
         if (cp.type == 2) {
           input <- gen.input(bx, scales, sq, TRUE, sgn)
           tmp <- apply(func_cusum(input)$acs, 1, max)
         }
         tmp
       }
  if (do.parallel > 0) {
    parallel::stopCluster(cl)}
  apply(ns, 1, function(z) {
    quantile(z, 1 - q)
  })
}

#core function
sbs.make.tree = function (input, tau = rep(1, nrow(input)), thr, trim, height, diag, n) {
  N <- dim(input)[1]
  len <- dim(input)[2]
  if (is.null(height)) {
    height <- floor(log(len, 2)/2)}
  if (is.null(trim)) {
    trim <- round(log(len)^1.5)}
  tree <- list(matrix(0, 5, 1))
  mat <- matrix(NA, ncol = 0, nrow = 6)
  acs <- func_cusum(input)$acs
  stat <- apply(acs * (acs > thr)/tau, 2, sum)
  stat[c((1:trim), (len - trim):len)] <- 0 #close to the start and end point
  sb <- search.b(stat, trim)
  change_dim_list = list()
  change_thr_list = list()
  change_thr_postprocessing_list = list()
  if (sb$FLAG) {
    # find the significant dimensions
    if (diag){
      sigdim = which(acs[,sb$b]>thr)
      sigthr = acs[,sb$b]
      sigthr[which(acs[,sb$b]<=thr)] = NA
    }else{
      sigthr = acs[,sb$b]
      sigthr[which(acs[,sb$b]<=thr)] = NA
      sigdim = c()
      sigdimtemp = which(acs[,sb$b]>thr)
      ijindex = 0
      for (ii in c(1:n)){
        for (jj in c(ii:n)){
          ijindex = ijindex +1
          if (ijindex%in%sigdimtemp){
            sigdim = c(sigdim,ii,jj)
          }
        }
      }
      sigdim = unique(sigdim)
    }
    tree[[1]][1, 1] <- 1
    tree[[1]][2, 1] <- 1
    tree[[1]][3, 1] <- sb$b
    tree[[1]][4, 1] <- len
    tree[[1]][5, 1] <- sb$test.stat
    change_dim_list[[as.character(sb$b)]] = sigdim #save the change dimension
    change_thr_list[[as.character(sb$b)]] = sigthr
    mat <- cbind(mat, c(1, tree[[1]][,]))
    j <- 1
    while (length(tree) == j & j < height) {
      npc <- dim(tree[[j]])[2]
      ncc <- 0
      i <- 1
      while (i <= npc) {
        if (tree[[j]][3, i] - tree[[j]][2, i] + 1 > 2 * #first part
            trim + 1) {
          s <- tree[[j]][2, i]
          e <- tree[[j]][3, i]
          acs <- func_cusum(input[, s:e])$acs
          stat <- apply(acs * (acs > thr)/tau, 2, sum)
          stat[c((1:trim), (e - s - trim + 1):(e - s))] <- 0
          sb <- search.b(stat, trim)
          if (sb$FLAG) {
            if (length(tree) == j) {
              tree <- c(tree, list(matrix(0, 5, 0)))}
            if (diag){
              sigdim = which(acs[,sb$b]>thr)
              sigthr = acs[,sb$b]
              sigthr[which(acs[,sb$b]<=thr)] = NA
            }else{
              sigthr = acs[,sb$b]
              sigthr[which(acs[,sb$b]<=thr)] = NA
              sigdim = c()
              sigdimtemp = which(acs[,sb$b]>thr)
              ijindex = 0
              for (ii in c(1:n)){
                for (jj in c(ii:n)){
                  ijindex = ijindex +1
                  if (ijindex%in%sigdimtemp){
                    sigdim = c(sigdim,ii,jj)
                  }
                }
              }
              sigdim = unique(sigdim)
            }
            ncc <- ncc + 1
            tree[[j + 1]] <- matrix(c(tree[[j + 1]], 
                                      matrix(0, 5, 1)), 5, ncc)
            tree[[j + 1]][1, ncc] <- 2 * tree[[j]][1,i] - 1
            tree[[j + 1]][2, ncc] <- s
            tree[[j + 1]][3, ncc] <- s + sb$b - 1
            tree[[j + 1]][4, ncc] <- e
            tree[[j + 1]][5, ncc] <- sb$test.stat
            change_dim_list[[as.character(s + sb$b - 1)]] = sigdim #save the change dimension
            change_thr_list[[as.character(s + sb$b - 1)]] = sigthr
            mat <- cbind(mat, c(j + 1, tree[[j + 1]][,ncc]))
          }
        }
        if (tree[[j]][4, i] - tree[[j]][3, i] > 2 * trim + 
            1) {#last part
          s <- tree[[j]][3, i] + 1
          e <- tree[[j]][4, i]
          acs <- func_cusum(input[, s:e])$acs
          stat <- apply(acs * (acs > thr)/tau, 2, sum)
          stat[c((1:trim), (e - s - trim + 1):(e - s))] <- 0
          sb <- search.b(stat, trim)
          if (sb$FLAG) {
            if (length(tree) == j) {
              tree <- c(tree, list(matrix(0, 5, 0)))}
            if (diag){
              sigdim = which(acs[,sb$b]>thr)
              sigthr = acs[,sb$b]
              sigthr[which(acs[,sb$b]<=thr)] = NA
            }else{
              sigthr = acs[,sb$b]
              sigthr[which(acs[,sb$b]<=thr)] = NA
              sigdim = c()
              sigdimtemp = which(acs[,sb$b]>thr)
              ijindex = 0
              for (ii in c(1:n)){
                for (jj in c(ii:n)){
                  ijindex = ijindex +1
                  if (ijindex%in%sigdimtemp){
                    sigdim = c(sigdim,ii,jj)
                  }
                }
              }
              sigdim = unique(sigdim)
            }
            ncc <- ncc + 1
            tree[[j + 1]] <- matrix(c(tree[[j + 1]], 
                                      matrix(0, 5, 1)), 5, ncc)
            tree[[j + 1]][1, ncc] <- 2 * tree[[j]][1, 
                                                   i]
            tree[[j + 1]][2, ncc] <- s
            tree[[j + 1]][3, ncc] <- s + sb$b - 1
            tree[[j + 1]][4, ncc] <- e
            tree[[j + 1]][5, ncc] <- sb$test.stat
            change_dim_list[[as.character(s + sb$b - 1)]] = sigdim #save the change dimension
            change_thr_list[[as.character(s + sb$b - 1)]] = sigthr
            mat <- cbind(mat, c(j + 1, tree[[j + 1]][,ncc]))
          }
        }
        i <- i + 1
      }
      j <- j + 1
    }
  }
  rownames(mat) <- c("level_index", "child_index", 
                     "s", "b", "e", "test_stat")
  #post processing, recalculate the CUSUM statistics for each change point
  ecp = sort(mat[4, ], decreasing = FALSE)
  if(length(ecp)>0){
    ecpnew = c(0,ecp,len)
    for (i in c(2:(length(ecpnew)-1))){
      s <- ecpnew[i-1] + 1
      e <- ecpnew[i+1]
      cp = ecpnew[i]
      acs <- func_cusum(input[, s:e])$acs
      sigthr = acs[,cp+1-s]
      sigthr[which(acs[,cp+1-s]<thr)] = NA
      change_thr_postprocessing_list[[as.character(cp)]] = sigthr
    }
  }
  
  return(structure(list(tree = tree, mat = mat, ecp = sort(mat[4, 
  ], decreasing = FALSE), thr = thr, change_dim_list = change_dim_list, change_thr_list = change_thr_list,
  change_thr_postprocessing_list = change_thr_postprocessing_list
  ), class = "bin.tree"))
}

clean.cp = function (z, type = c("dcbs", "sbs"), phi = 0.5, #for mean test
          trim = NULL, height = NULL) 
{
  len <- dim(z)[2]
  n <- dim(z)[1]
  if (is.null(trim)) 
    trim <- round(log(len))
  if (is.null(height)) 
    height <- floor(log(len, 2)/2)
  if (type == "dcbs") 
    mat <- dcbs.make.tree(z, phi, trim, height)$mat
  if (type == "sbs") 
    mat <- sbs.make.tree(z, rep(1, n), apply(func_cusum(z)$acs, 
                                             1, median), trim, 2)$mat
  brks <- sort(unique(c(0, mat[4, ], len)), decreasing = FALSE)
  for (l in 1:(length(brks) - 1)) {
    int <- (brks[l] + 1):brks[l + 1]
    z[, int] <- z[, int] - apply(z[, int, drop = FALSE], 
                                 1, mean)
  }
  return(list(x = z, mat = mat))
}

gen.input = function (x, scales, sq, diag, sgn = NULL) { #can set multiple scales
  if (is.null(sgn)) {
    sgn <- sign(cor(t(x)))
  }
  input <- NULL
  for (sc in scales) {
    cc <- func_coef(x, sc) #time ordered NDWT, Haar basis
    #testwd = wd(c(x[1,],rep(0,212)), filter.number=1, family="DaubExPhase",type = 'station')
    #testwd = accessD(testwd, lev = 8)
    if (sc > min(scales)) {#only use the whole covered d
      cc <- cc[, -(1:(2^(-min(scales)) - 2^(-sc))), drop = FALSE]#drop = FALSE keeps the dataframe
    }
    #calculate the periodogram
    input <- rbind(input, t(func_input(cc, sgn, sq, diag)))
  }
  input
}

search.b = function (stat, trim) {
  FLAG <- FALSE
  b <- test.stat <- NA
  while (!FLAG && sum(stat) > 0) {
    test.stat <- max(stat)
    b <- min(which(stat == test.stat))
    int <- max(b - trim + 1, 1):min(b + trim, length(stat))
    if (sum(stat[int] == 0) > 0) {
      stat[int] <- 0
    }else {
      FLAG <- TRUE
    }
  }
  return(list(b = b, test.stat = test.stat, FLAG = FLAG))
}


cp_detection_classification = function(seizure_sig, keyscale, thr = NULL, diag, q = 0.01,
                                       final_cp_thrdf, final_thrdf, art_sig, needplot = TRUE,
                                       start_index){
  
  sbsresult = sbs.alg(t(seizure_sig), cp.type=2, scales=-keyscale, trim = 64, thr = thr,
                      diag=diag, do.parallel=10, .verbose = FALSE, sq = FALSE, q = q)
  ecp = sbsresult$ecp
  
  #save the change point locations and their thresholds
  if (length(ecp)>0) {
    for (i in c(1:length(sbsresult$tree))){
      treelevel_i = sbsresult$tree[[i]]
      if (i == 1){
        if (ncol(treelevel_i)!=1){
          stop('error in i, j interation')
        }
        treeson_ij = treelevel_i[,1]
        cp_thr = c(sbsresult$change_thr_postprocessing_list[[as.character(treeson_ij[3])]],treeson_ij[3])
        cp_thr = matrix(cp_thr,nrow = 1)
      } else if (i>1){
        for (j in c(1:ncol(treelevel_i))){
          treeson_ij = treelevel_i[,j]
          cp_pre = treeson_ij[c(2,4)]
          cp_pre = c(cp_pre,cp_pre-1)
          cp_pre = cp_thr[which(cp_thr[,ncol(cp_thr)]%in%cp_pre),-ncol(cp_thr),drop = FALSE]
          cp_now = rbind(cp_pre,sbsresult$change_thr_postprocessing_list[[as.character(treeson_ij[3])]])
          cp_now = apply(cp_now,2,min, na.rm = TRUE)
          cp_now[which(cp_now==Inf)] = NA
          cp_thr = rbind(cp_thr,c(cp_now,treeson_ij[3]))
        }
      }
    }
    cp_thr[,ncol(cp_thr)] = cp_thr[,ncol(cp_thr)]+start_index
    colnames(cp_thr) = colnames(final_cp_thrdf[[keyscale]])
    final_cp_thrdf[[keyscale]] = rbind(final_cp_thrdf[[keyscale]],cp_thr)
  }

  #plot the threshold from the margin, inner seizure and non-seizure period groups
  innerseizure = which(art_sig$label1==1)
  thrdf = as.data.frame(matrix(sbsresult$thr,nrow = 1))
  thrdf = cbind(thrdf,'critical value (q%)')
  thrdf[,ncol(thrdf)] = as.character(thrdf[,ncol(thrdf)])
  colnames(thrdf)[ncol(thrdf)] = 'cptype'
  if (diag){
    colnames(thrdf) = colnames(final_thrdf[[keyscale]])
  }
  if (length(innerseizure)>0) { #exist seizure
    if (length(ecp)>0) { #detect seizure
      first_cp = ecp[which.min(abs(ecp-min(innerseizure)))]
      last_cp = ecp[which.min(abs(ecp-max(innerseizure)))]
      inner_cp = ecp[which(ecp>first_cp&ecp<last_cp)]
      outer_cp = ecp[which(ecp<first_cp|ecp>last_cp)]
      for (i in unique(c(first_cp,last_cp))){
        thrdf = rbind(thrdf,c(sbsresult$change_thr_postprocessing_list[[as.character(i)]]))
        thrdf[nrow(thrdf),ncol(thrdf)] = 'margin'
      }
      for (i in inner_cp){
        thrdf = rbind(thrdf,c(sbsresult$change_thr_postprocessing_list[[as.character(i)]]))
        thrdf[nrow(thrdf),ncol(thrdf)] = 'inner'
      }
      for (i in outer_cp){
        thrdf = rbind(thrdf,c(sbsresult$change_thr_postprocessing_list[[as.character(i)]]))
        thrdf[nrow(thrdf),ncol(thrdf)] = 'outer'
      }
    }
  }else{
    if (length(ecp)>0) {
      outer_cp = ecp
      for (i in outer_cp){
        thrdf = rbind(thrdf,c(sbsresult$change_thr_postprocessing_list[[as.character(i)]]))
        thrdf[nrow(thrdf),ncol(thrdf)] = 'outer'
      }
    }
  }
  final_thrdf[[keyscale]] = rbind(final_thrdf[[keyscale]],thrdf)
  if (needplot){
    plotdf = reshape2::melt(thrdf, id.vars = "cptype")
    p = ggplot(plotdf,aes(x=variable))+geom_boxplot(aes(y=value,color = cptype))
    p = p + theme_bw() + ylab('threshold') + scale_color_npg()+
      theme(axis.text=element_text(size=14),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title=element_text(size=14,face="bold"),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14))
    ggsave(paste0(plotdir,'start_index_',start_index,
                  '_ecp_thr_',keyscale,'.jpg'), width=15, height=12,units = "in")
  }
  
  return(list('final_thrdf'=final_thrdf,'final_cp_thrdf'=final_cp_thrdf,'ecp' = ecp))
}

cp_threshold_cal = function(EEGdata,start_time,end_time,refrow,ICnum,cp_stat_list,
                            CUSUM_q,scale_range,sample_freq,diag,thr = NULL,ifplot = FALSE){
  
  final_variance_df = cp_stat_list[['final_variance_df']]
  final_ratio_variance_df = cp_stat_list[['final_ratio_variance_df']]
  final_thrdf = cp_stat_list[['final_thrdf']]
  final_cp_thrdf = cp_stat_list[['final_cp_thrdf']]

  hour = EEGdata$hour[start_time]
  minute = EEGdata$minute[start_time]
  sec = EEGdata$sec[start_time]
  hertz = EEGdata$hertz[start_time]
  
  data_seg = EEGdata[c(start_time:end_time),1:28]
  
  datamat = data_seg[,c(8:28)]
  
  SCICAresult = SCICA(datamat,ncol(datamat), alg.typ = "parallel",fun = 'logcosh', alpha = 1, 
                      row.norm = FALSE, maxit = 200, sc.typ = 'soft',
                      tol = 0.0001, verbose = FALSE, SC = refrow)
  
  Smat = SCICAresult$S
  
  seizure_sig_1d = rowMeans(Smat[,1:ICnum,drop = FALSE]%*%SCICAresult$A[1:ICnum,,drop = FALSE])
  seizure_sig_1d = scale(seizure_sig_1d, scale = FALSE)
  seizure_raw_1d = rowMeans(data_seg[,8:28])
  
  #local variance estimator
  sigmahat_kernel_DPI = local_var(seizure_sig_1d,ifplot=ifplot,plotindex=start_time)
  rawsigmahat_kernel_DPI = local_var(seizure_raw_1d)
  
  ratiosigma = sigmahat_kernel_DPI/rawsigmahat_kernel_DPI
  
  tempvardf = data.frame('label' = data_seg$label1, 'variance' = sigmahat_kernel_DPI)
  final_variance_df = rbind(final_variance_df,tempvardf)
  tempvardf = data.frame('label' = data_seg$label1, 'variance' = ratiosigma)
  final_ratio_variance_df = rbind(final_ratio_variance_df,tempvardf)
  
  # inverse ICA for 21 dimensional seizure signal
  newW = SCICAresult$K%*%SCICAresult$W
  newW = newW[,1:ICnum,drop =FALSE]
  for (i in 1:ncol(newW)){
    newW[,i] = newW[,i]/sqrt(sum(newW[,i]^2))
  }
  seizure_sig = as.matrix(datamat)%*%newW
  
  df =as.data.frame(cbind(data_seg[,c(1)],seizure_sig,data_seg$label1))
  colnames(df)[ncol(df)] = 'label'
  colnames(df)[1] = 'X'
  df$label = as.character(df$label)
  
  plotdf = reshape2::melt(df, id.vars = c("X","label"))
  
  for(scale_i in scale_range){
    if (is.null(thr)){
      thr_i = NULL
    }else{
      thr_i = thr[[scale_i]]
    }
    combined_output = 
      cp_detection_classification(seizure_sig = seizure_sig, keyscale = scale_i, 
                                  thr = thr_i, diag = diag, 
                                  final_cp_thrdf = final_cp_thrdf, final_thrdf=final_thrdf,
                                  art_sig = data_seg, needplot = ifplot,q = CUSUM_q,
                                  start_index = EEGdata$index[start_time]-1)
    final_thrdf = combined_output$final_thrdf
    final_cp_thrdf = combined_output$final_cp_thrdf
    ecp_i = combined_output$ecp
    
    if(ifplot){
      # plot SCICs
      mt <- ggplot(plotdf, aes(X, value)) +
        geom_line(aes(color = label), show.legend = FALSE) + facet_grid(vars(variable), scales = "free") + 
        geom_vline(xintercept = (ecp_i+df$X[1]-1),color="red", linetype="dashed") + 
        ggtitle(paste0(hour,':',minute,':',sec,'(',hertz,
                       ')_sbs_seizure_detection'))+
        scale_x_continuous(breaks = seq(from = df$X[1], to = df$X[nrow(df)], by = sample_freq),
                           minor_breaks = NULL)+ guides(x = guide_axis(angle = 90))+
        theme_minimal()+ scale_colour_brewer(palette = "Set1",direction = -1)
      #mt
      ggsave(paste0(plotdir,'start_from_',hour,'_',minute,'_',sec,'(',hertz,
                    ')_sbs_seizure_detection_',scale_i,'.jpg'), 
             width=15, height=12,units = "in")
    }
  }
  
  cp_stat_list[['final_variance_df']] = final_variance_df
  cp_stat_list[['final_ratio_variance_df']] = final_ratio_variance_df
  cp_stat_list[['final_thrdf']] = final_thrdf
  cp_stat_list[['final_cp_thrdf']] = final_cp_thrdf
  
  if(ifplot){
    #local variance boxplot
    plotdf = data.frame('label' = data_seg$label1, 'variance' = sigmahat_kernel_DPI)
    final_variance_df = rbind(final_variance_df,plotdf)
    plotdf$label = as.character(plotdf$label)
    p = ggplot(plotdf)+geom_boxplot(aes(y=variance,color = label))
    p = p + theme_bw() + ylab('variance') + scale_color_npg()+
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14))
    ggsave(paste0(plotdir,'start_from_',hour,'_',minute,'_',sec,'(',hertz,
                  ')_local_var_comp.jpg'), width=15, height=12,units = "in")
    
    ratiosigma = sigmahat_kernel_DPI/rawsigmahat_kernel_DPI
    
    #local variance boxplot
    plotdf = data.frame('label' = data_seg$label1, 'variance' = ratiosigma)
    final_ratio_variance_df = rbind(final_ratio_variance_df,plotdf)
    plotdf$label = as.character(plotdf$label)
    p = ggplot(plotdf)+geom_boxplot(aes(y=variance,color = label))
    p = p + theme_bw() + ylab('variance') + scale_color_npg()+
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14))
    ggsave(paste0(plotdir,'start_from_',hour,'_',minute,'_',sec,'(',hertz,
                  ')_local_var_ratio_comp.jpg'), width=15, height=12,units = "in")
  }
  
  return(cp_stat_list)
}