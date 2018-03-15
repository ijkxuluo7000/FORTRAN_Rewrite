#Sample Data input
###############################
inspo2 = c(21,21,21,21)
po2 = c(50,53,54,55)
pco2 = c(25.6,26.3,27.2,24.6)
ph = c(7.323,7.309,7.256,7.185,7.356)
sat = c(86.4,94.4,93.0,85.0)
hb = c(19.2,20.0,20.6,20.5)
hct = c(71,71,72,73)
temp = c(36.8,36.7,36.7,36.7)
dat = as.data.frame(cbind(inspo2,po2,pco2,ph,sat,hb,hct,temp))
###############################
#Function Saturation
saturation = function(PO2,TEMP,PH,PCO2,P50){
  a1 = -8532.229;a2 = 2121.401
  a3 = -67.07399;a4 = 935960.9
  a5 = -31346.26;a6 = 2396.167
  a7 = -67.10441
  b = 0.43429 * log(40.0/PCO2,base = 10)
  xx = PO2 * 10^(0.024 * (37 - TEMP) + 0.4 * (PH - 7.4) + 0.06 * b)
  x = 26.8 * (xx/P50)
  sat = (x*(x*(x*(x+a3)+a2)+a1))/(x*(x*(x*(x+a7)+a6)+a5)+a4)
  return(100 * sat)
}

Apply = function (X, MARGIN, FUN, ...) {
    FUN <- match.fun(FUN)
    dl <- length(dim(X))
    if (!dl) 
      stop("dim(X) must have a positive length")
    if (is.object(X)) 
      X <- if (dl == 2L) 
        as.matrix(X)
    else as.array(X)
    d <- dim(X)
    dn <- dimnames(X)
    ds <- seq_len(dl)
    if (is.character(MARGIN)) {
      if (is.null(dnn <- names(dn))) 
        stop("'X' must have named dimnames")
      MARGIN <- match(MARGIN, dnn)
      if (anyNA(MARGIN)) 
        stop("not all elements of 'MARGIN' are names of dimensions")
    }
    s.call <- ds[-MARGIN]
    s.ans <- ds[MARGIN]
    d.call <- d[-MARGIN]
    d.ans <- d[MARGIN]
    dn.call <- dn[-MARGIN]
    dn.ans <- dn[MARGIN]
    d2 <- prod(d.ans)
    if (d2 == 0L) {
      newX <- array(vector(typeof(X), 1L), dim = c(prod(d.call), 
                                                   1L))
      ans <- forceAndCall(1, FUN, if (length(d.call) < 2L) newX[, 
                                                                1] else array(newX[, 1L], d.call, dn.call), ...)
      return(if (is.null(ans)) ans else if (length(d.ans) < 
                                            2L) ans[1L][-1L] else array(ans, d.ans, dn.ans))
    }
    newX <- aperm(X, c(s.call, s.ans))
    dim(newX) <- c(prod(d.call), d2)
    ans <- vector("list", d2)
    if (length(d.call) < 2L) {
      if (length(dn.call)) 
        dimnames(newX) <- c(dn.call, list(NULL))
      for (i in 1L:d2) {
        tmp <- forceAndCall(1, FUN, newX[, i], ...)
        if (!is.null(tmp)) 
          ans[[i]] <- tmp
      }
    }
    else for (i in 1L:d2) {
      tmp <- forceAndCall(1, FUN, array(newX[, i], d.call, 
                                        dn.call), ...)
      if (!is.null(tmp)) 
        ans[[i]] <- tmp
    }
    ans.list <- is.recursive(ans[[1L]])
    l.ans <- length(ans[[1L]])
    ans.names <- names(ans[[1L]])
    if (!ans.list) 
      ans.list <- any(lengths(ans) != l.ans)
    if (!ans.list && length(ans.names)) {
      all.same <- vapply(ans, function(x) identical(names(x), 
                                                    ans.names), NA)
      if (!all(all.same)) 
        ans.names <- NULL
    }
    len.a <- if (ans.list) 
      d2
    else length(ans <- unlist(ans, recursive = FALSE))
    if (length(MARGIN) == 1L && len.a == d2) {
      names(ans) <- if (length(dn.ans[[1L]])) 
        dn.ans[[1L]]
      ans
    }
    else if (len.a == d2) 
      array(ans, d.ans, dn.ans)
    else if (len.a && len.a%%d2 == 0L) {
      if (is.null(dn.ans)) 
        dn.ans <- vector(mode = "list", length(d.ans))
      dn1 <- list(ans.names)
      if (length(dn.call) && !is.null(n1 <- names(dn <- dn.call[1])) && 
          nzchar(n1) && length(ans.names) == length(dn[[1]])) 
        names(dn1) <- n1
      dn.ans <- c(dn1, dn.ans)
      array(ans, c(len.a%/%d2, d.ans), if (!is.null(names(dn.ans)) || 
                                           !all(vapply(dn.ans, is.null, NA))) 
        dn.ans)
    }
    else ans
  }

lower_bond = 10
upper_bond = 60
digit = 1
vp50 = seq(lower_bond,upper_bond,10^(-digit))

size = length(dat[,1])
satura_mat = matrix(ncol = size,nrow = length(vp50))

for(i in 1:size){
  for(j in 1:length(vp50)){
    satura_mat[j,i] = saturation(PO2 = dat$po2[i],
                               TEMP = dat$temp[i],
                               PH = dat$ph[i],
                               PCO2 = dat$pco2[i],
                               P50 = vp50[j])
  }
}

ssqvec = vector(length = length(vp50))
for(j in 1:length(vp50)){
  ssq = vector(length = size)
  for(i in 1:size){
    ssq[i] = (satura_mat[,i] - dat$sat[i])^2
  }
  ssq = sum(ssq)
  ssqvec[j] = ssq
}
            
print(ssq)
print(ssqvec)
