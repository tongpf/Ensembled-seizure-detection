# by tongpf 2021/04/27, tongpf@pku.edu.cn
# modified from fastICA algorithm, adaptive to spatial constraint

SCICA = function (X, n.comp, alg.typ = "parallel",#c("parallel", "deflation"),
                  sc.typ = c('hard','soft','weak'), fun = "logcosh",
                  alpha = 1, row.norm = FALSE, maxit = 200, tol = 1e-04, 
                  verbose = FALSE, w.init = NULL, SC = NULL) 
{
  # the output X is the demeaned input X
  # X = SA
  # XKW = S
  dd <- dim(X)
  d <- dd[dd != 1L]
  if (length(d) != 2L) {
    stop("data must be matrix-conformal")}
  if (length(d) != length(dd)) {
    X = matrix(X, d[1L], d[2L])
  }else {
    X = as.matrix(X)
  }
  if (alpha < 1 || alpha > 2) {
    stop("alpha must be in range [1,2]")}
  #alg.typ <- match.arg(alg.typ)
  #fun <- match.arg(fun)
  sc.typ = match.arg(sc.typ)
  n <- nrow(X)
  p <- ncol(X)
  if (n.comp > min(n, p)) {
    message("'n.comp' is too large: reset to ", min(n,p))
    n.comp <- min(n, p)
  }
  if (verbose) {
    message("Centering")}
  X <- scale(X, scale = FALSE)
  X <- if (row.norm) {
    t(scale(X, scale = row.norm))
  }else{ 
    t(X)
  }
  if (verbose) {
    message("Whitening")}
  V <- X %*% t(X)/n
  s <- La.svd(V)
  D <- diag(c(1/sqrt(s$d)))
  K <- D %*% t(s$u)
  K <- matrix(K[1:n.comp, ], n.comp, p)
  X1 <- K %*% X #whitening use Sigma^-1/2
  if (is.null(w.init)) {
    if (is.null(SC)){
      SCnum = 0
      w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
    }else{
      SCnum = nrow(SC)
      wc.qr = qr(K%*%t(SC))
      wc = qr.Q(wc.qr)
      w.init <- matrix(rnorm(n.comp*(n.comp-ncol(wc))), n.comp, (n.comp-ncol(wc)))
      w.init = cbind(wc,w.init)
      w.qr = qr(w.init)
      w.init = qr.Q(w.qr)
      w.init = t(w.init)
    }
  }else {
    SCnum = 0
    if (!is.matrix(w.init) || length(w.init) != (n.comp^2)) 
      stop("w.init is not a matrix or is the wrong size")
  }
  if (alg.typ == "parallel") {
    a <- scica.R.par(X1, n.comp, tol = tol, fun = fun, alpha = alpha, 
                     maxit = maxit, verbose = verbose, w.init = w.init, SCnum = SCnum, sc.typ = sc.typ)
  }else if (alg.typ == "deflation") {
    # not complete yet!
    stop("method deflation dose not complete yet!")
    a <- scica.R.def(X1, n.comp, tol = tol, fun = fun, alpha = alpha, 
                     maxit = maxit, verbose = verbose, w.init = w.init, SCnum = SCnum, sc.typ = sc.typ)
  }
  w <- a %*% K
  S <- w %*% X #aKX=S
  #w = w+diag(rep(10^(-5),ncol(w)))
  #print(w %*% t(w))
  #print(qr(w)$rank)
  #A <- t(w) %*% solve(w %*% t(w))
  #A <- solve(t(w) %*% w)%*%t(w)
  A <- solve(w)
  return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S)))
  #K%*%t(K) = Sigma, t(K)%*%K = Lambda
}


scica.R.par = function (X, n.comp, tol, fun, alpha, maxit, verbose, w.init, SCnum, sc.typ) 
{
  Diag <- function(d) {
    if (length(d) > 1L) {
      diag(d)
    }else {
      as.matrix(d)
    }
  }
  if (SCnum == 0|sc.typ == 'weak'){
    p <- ncol(X)
    W <- t(w.init)
    sW <- La.svd(W)
    W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W #use svd for orthogonalization
    W1 <- W
    lim <- rep(1000, maxit)
    it <- 1
    if (fun == "logcosh") {
      if (verbose) {
        message("Symmetric FastICA using logcosh approx. to neg-entropy function")}
      while (lim[it] > tol && it < maxit) {#here X has been standardized
        wx <- W %*% X #the independent component S
        #based on Kuhn-Tucker conditions
        gwx <- tanh(alpha * wx)
        v1 <- gwx %*% t(X)/p
        g.wx <- alpha * (1 - (gwx)^2)
        v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
        W1 <- v1 - v2 #Newton's method
        sW1 <- La.svd(W1)
        W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
        lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
        W <- W1
        if (verbose) {
          message("Iteration ", it, " tol = ", format(lim[it + 1]))}
        it <- it + 1
      }
    }
  }else{
    p <- ncol(X)
    W <- t(w.init)
    hcref = matrix(W[,1:SCnum],ncol = SCnum)
    #sW <- La.svd(W)
    #W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W #controlling the l2 norm of W to be 1？
    W1 <- W
    lim <- rep(1000, maxit)
    it <- 1
    if (verbose) {
      message("Symmetric FastICA using logcosh approx. to neg-entropy function")}
    while (lim[it] > tol && it < maxit) {#here X has been standardized
      hc = matrix(W[,1:SCnum],ncol = SCnum)
      hu = matrix(W[,(SCnum+1):n.comp],ncol = n.comp-SCnum)
      wx <- t(hu) %*% X #the independent component S
      #based on Kuhn-Tucker conditions
      gwx <- tanh(alpha * wx)
      v1 <- t(gwx %*% t(X)/p)
      g.wx <- alpha * (1 - (gwx)^2)
      v2 <- t(Diag(apply(g.wx, 1, FUN = mean)) %*% t(hu))
      hu1 <- v1 - v2 #Newton's method (for hu only)
      #sW1 <- La.svd(hu1)
      #hu1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% hu1
      if (sc.typ == "soft") {
        W1 = cbind(hu1,hc)
        W1.qr = qr(W1)
        W1 = qr.Q(W1.qr)
        hc = matrix(W1[,(n.comp-SCnum+1):n.comp],ncol = SCnum)
        hu1 = matrix(W1[,1:(n.comp-SCnum)],ncol = (n.comp-SCnum))
        hc = scica_sc_update(hc, hcref)
      }
      W1 = cbind(hc,hu1)
      W1.qr = qr(W1)
      W1 = qr.Q(W1.qr)
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      W <- W1
      if (verbose) {
        message("Iteration ", it, " tol = ", format(lim[it + 1]))}
      it <- it + 1
    }
  }
  t(W)
}

scica.R.def = function (X, n.comp, tol, fun, alpha, maxit, verbose, w.init, SCnum, sc.typ) 
{
  #not complete yet
  if (verbose && fun == "logcosh") 
    message("Deflation FastICA using logcosh approx. to neg-entropy function")
  p <- ncol(X)
  W <- matrix(0, n.comp, n.comp)
  for (i in 1:n.comp) {
    if (SCnum == 0|sc.typ == 'weak'){
      if (verbose) {
        message("Component ", i)}
      w <- matrix(w.init[i, ], n.comp, 1) #transpose
      if (i > 1) {
        t <- w
        t[1:length(t)] <- 0
        for (u in 1:(i - 1)) {
          k <- sum(w * W[u, ])
          t <- t + k * W[u, ]
        }
        w <- w - t #minus inner product
      }
      w <- w/sqrt(sum(w^2))
      lim <- rep(1000, maxit)
      it <- 1
      while (lim[it] > tol && it < maxit) {
        wx <- t(w) %*% X
        gwx <- tanh(alpha * wx)
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X * gwx
        v1 <- apply(xgwx, 1, FUN = mean)
        g.wx <- alpha * (1 - (tanh(alpha * wx))^2)
        v2 <- mean(g.wx) * w
        w1 <- v1 - v2
        w1 <- matrix(w1, n.comp, 1)
        it <- it + 1
        if (i > 1) {
          t <- w1
          t[1:length(t)] <- 0
          for (u in 1:(i - 1)) {
            k <- sum(w1 * W[u, ])
            t <- t + k * W[u, ]
          }
          w1 <- w1 - t
        }
        w1 <- w1/sqrt(sum(w1^2))
        lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
        if (verbose) 
          message("Iteration ", it - 1, " tol = ", 
                  format(lim[it]))
        w <- matrix(w1, n.comp, 1)
      }
      W[i, ] <- w
    }else{
      p <- ncol(X)
      hcref = matrix(w.init[,1:SCnum],ncol = SCnum)
      W <- w.init
      #sW <- La.svd(W)
      #W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W 
      W1 <- W
      lim <- rep(1000, maxit)
      it <- 1
      
      
      if (verbose) {
        message("Symmetric FastICA using logcosh approx. to neg-entropy function")}
      while (lim[it] > tol && it < maxit) {#
        hc = matrix(W[,1:SCnum],ncol = SCnum)
        hu = matrix(W[,(SCnum+1):n.comp],ncol = n.comp-SCnum)
        wx <- t(hu) %*% X #
        #
        gwx <- tanh(alpha * wx)
        v1 <- t(gwx %*% t(X)/p)#
        g.wx <- alpha * (1 - (gwx)^2)
        v2 <- t(Diag(apply(g.wx, 1, FUN = mean)) %*% t(hu))
        hu1 <- v1 - v2 #
        #sW1 <- La.svd(hu1)
        #hu1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% hu1
        if (sc.typ == "soft") {
          W1 = cbind(hu1,hc)
          W1.qr = qr(W1)
          W1 = qr.Q(W1.qr)
          hc = matrix(W1[,(n.comp-SCnum+1):n.comp],ncol = SCnum)
          hu1 = matrix(W1[,1:(n.comp-SCnum)],ncol = (n.comp-SCnum))
          hc = scica_sc_update(hc, hcref)
        }
        W1 = cbind(hc,hu1)
        W1.qr = qr(W1)
        W1 = qr.Q(W1.qr)
        lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
        W <- W1
        if (verbose) {
          message("Iteration ", it, " tol = ", format(lim[it + 1]))}
        it <- it + 1
      }
    }
  }
  W
}

scica_sc_update = function(hc,hcref,maxangle = 1/180*pi){
  hcreturn = hc
  for (i in 1:ncol(hcref)){
    a = hc[,i]
    b = hcref[,i]
    theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
    if (abs(theta)<maxangle){
      hcreturn[,i] = a
    }else{
      aa = t(a)%*%a
      bb = as.numeric(t(b)%*%b)
      ab = as.numeric(t(b)%*%a)
      c1 = (1-2*ab+(ab)^2)-((2-2*ab)*(cos(maxangle))^2)
      c2 = (2*ab-2*(ab)^2)-((2*ab-2)*(cos(maxangle))^2)
      c3 = (ab)^2-(cos(maxangle))^2
      
      temp = sqrt(c2^2-4*c1*c3)
      x1 = (-c2+temp)/2/c1
      x2 = (-c2-temp)/2/c1
      if (x1<1&x1>0){
        y = (1-x1)*a+x1*b
        y = y/sqrt(sum(y*y))
      }else{
        y = (1-x2)*a+x2*b
        y = y/sqrt(sum(y*y))
      }
      hcreturn[,i] = y
    }
  }
  return(hcreturn)
}

