cnmlcd = function(x, lcd, maxit=100, tol=1e-6,
                  plot=c("null","density","logdensity","gradient")) {
  plot = match.arg(plot)
  lambda = 1e-15             # for psuedo-rank
  n = length(x)              # sample size
  xw =  x.weight(x)
  x = xw$x                   # unique values in x
  w = xw$w                   # frequencies/weights
  nx = length(x)             # number of unique x-values
  lower = x[1]               # lower boundary
  upper = x[nx]              # upper boundary
  if(nx <= 2) {
    lcd = new.lcd(alpha=0, lower=lower, upper=upper)
    return(list(lcd=lcd, ll=logLik.lcd(lcd, x, w),
                num.iterations=0, max.gradient=0, convergence=0))
  }
  attr(x, "xx") =                 # - colSums(pmax(outer(x,x,"-"),0) * w)
    rev(cumsum(rev(w))) * x - rev(cumsum(rev(x * w)))
  if(missing(lcd)) lcd = new.lcd(alpha=0, lower=lower, upper=upper)
  ll = logLik.lcd(lcd, x, w)
  convergence = 1
  ll.old = -Inf
  for(i in 1:maxit) {
    if(ll <= ll.old  + tol) {convergence = 0; break}
    lcd.old = lcd
    ll.old = ll
    switch(plot,
           density = plot(lcd, x, w),
           logdensity = plot(lcd, x, w, log=TRUE),
           gradient = plotgradient(lcd, x, w),
           null =, )
    g = maxima.gradient(lcd, x, w=w)      # finding new knots
    if(length(g$theta) != 0) {
      nsp = g$theta
      nsl = length(nsp)
      if(nsl >= 1) {
        if(plot=="gradient") points(g$theta, g$gradient, col="red", pch=20)
        lcd = new.lcd(lcd$alpha, c(lcd$theta, nsp),
                      c(lcd$pi, double(nsl)), lcd$lower, lcd$upper)
      }
    }
    knots = c(lcd$lower, lcd$theta)
    nk = length(knots)
    cpkr = lcd$cpk[,nk] - cbind(0, lcd$cpk[,-nk,drop=FALSE])  # cpk reversed
    mu = cpkr[2,] - cpkr[1,] * knots            # E{(X-theta)_+}
    grad = n * mu + attr(x, "xx")[indx(knots, x)]    # gradient vector
    mm = cpkr[3,] -                             # E{(X-theta_j)_+ (X-theta_k)_+}
      (knots + rep(knots, rep.int(nk, nk))) * cpkr[2,] +
      tcrossprod(knots) * cpkr[1,]
    mm[upper.tri(mm)] = 0
    mm = mm + t(mm)
    diag(mm) = diag(mm) / 2
    H = mm - tcrossprod(mu)                     # negative Hessian matrix
    e = eigen(H)
    v2 = sqrt(e$values[e$values >= e$values[1] * lambda])
    kr = length(v2)                             # psuedo-rank
    R = t(e$vectors[,1:kr,drop=FALSE]) * v2
    p = grad / n + drop(c(-lcd$alpha, lcd$pi) %*% H)
    b = drop(crossprod(e$vectors[,1:kr,drop=FALSE], p)) / v2
    r1 = pnnls(R,b,1)
    lcd1 = lcd
    lcd1$alpha = -r1$x[1]       # lcd1$C will be updated inside line.lcd
    lcd1$pi = r1$x[-1]
    r = line.lcd(lcd, lcd1, x, w=w, ll0=ll.old)
    lcd = r$lcd
    ll = r$ll
    ## if(any(j0 <- lcd$pi == 0)) {
    ##   j = ! j0
    ##   lcd = new.lcd(lcd$alpha, lcd$theta[j], lcd$pi[j], lcd$lower, lcd$upper)
    ## }
    if(any(lcd$pi == 0)) lcd = simplify.lcd(lcd)
  }
  list(lcd=lcd, ll=ll, num.iterations=i, max.gradient=g$gmax,
       convergence=convergence)
}

## Line search

line.lcd = function(lcd0, lcd1, x, w, ll0, tol=1e-10) {
  llt = function(alpha) {
    m = new.lcd((1 - alpha) * lcd0$alpha + alpha * lcd1$alpha, lcd1$theta,
        (1 - alpha) * lcd0$pi + alpha * lcd1$pi, lcd0$lower, lcd0$upper)
    ll = logLik.lcd(m, x, w)
    list(lcd=m, ll=ll)
  }
  grad = gradient.lcd(lcd0, x, w, "knots")
  grad[1] = - grad[1]
  delta = sum(grad * (c(lcd1$alpha, lcd1$pi) - c(lcd0$alpha, lcd0$pi))) * .333
  convergence = 0
  alpha = 1
  repeat{
    new = llt(alpha)
    if(new$ll >= ll0 + alpha * delta) break 
    if(alpha < tol) {convergence=1; new=list(lcd=lcd0, ll=ll0); break}
    alpha = alpha * .5
  }
  list(lcd=new$lcd, ll=new$ll, alpha=alpha, convergence=convergence)
}

## Create a new 'lcd' object

## alpha    first slope of the log-density
## theta    interior knots
## pi       changes of slope at the interior knots
## lower    lower boundary
## upper    upper boundary

##### Output: A list consisting of
## alpha      Initial slope on [x_1, x_2]
## C          Normalizing constant
## theta      Interior knots (knots other than x_1 and x_n)
## pi         Changes of slope at theta
## lower      Lower boundary (= x_1)
## upper      Upper boundary (= x_n)
## coef       Intercepts (in row 1) and slopes (in row 2) over all segments
## fk         Density at knots = c(lower, theta)
## dpk        Int x^o f(x) dx over each segment between knots for o = 0, 1, 2
## cpk        Cumulative sum of dpk for all segments

new.lcd = function(alpha, theta=NULL, pi=NULL, lower, upper) {
  if(length(theta) >  0) {
    k = length(theta)
    if(length(pi) < k) pi = rep(pi, len=k)
    o = order(theta)
    theta = theta[o]
    pi = pi[o]
  }
  else theta = pi = numeric(0)
  c0 = -alpha * lower + c(0, cumsum(pi * theta))  # intercepts 
  c1 = alpha - c(0, cumsum(pi))                   # slopes
  knots1 = c(lower, theta)
  knots2 = c(theta, upper)
  dk = knots2 - knots1
  nk = length(dk)
  fk1 = exp(c0 + c1 * knots1)
  fk2 = c(fk1[-1], exp(c0[nk] + c1[nk] * upper))

  dpk = matrix(0, nrow=3, ncol=nk)
  rownames(dpk) = paste0("x", 0:2)
  if(any(j <- c1 == 0)) {                       # horizontal segments
    dpk[1,j] = fk1[j] * dk[j]
    dpk[2,j] = dpk[1,j] * (knots2[j] + knots1[j]) * .5
    dpk[3,j] = dpk[1,j] * (knots2[j]^2 + knots2[j] * knots1[j] + knots1[j]^2) / 3
  }
  if(any(j2 <- ! j)) {                          # non-horizontal segments
    x1 = fk1[j2] / c1[j2]
    y1 = fk2[j2] / c1[j2]
    dpk[1,j2] = y1 - x1
    dpk[2,j2] = knots2[j2] * y1 - knots1[j2] * x1 - dpk[1,j2] / c1[j2]
    dpk[3,j2] = knots2[j2]^2 * y1 - knots1[j2]^2 * x1 - 2 * dpk[2,j2] / c1[j2]
  }
  C = sum(dpk[1,])                              # normalizing constant
  fk = fk1 / C
  cpk = dpk = dpk / C                           # normalization
  for(i in 1:nrow(dpk)) cpk[i,] = cumsum(cpk[i,])
  structure(list(alpha=alpha, C=C, theta=theta, pi=pi, lower=lower, upper=upper,
                 coef=rbind(c0,c1), fk=fk, dpk=dpk, cpk=cpk),
            class = "lcd")
}

simplify.lcd = function(lcd) {
  if(any(j0 <- lcd$pi == 0)) {
    nk = length(lcd$theta) + 1
    j = which(!j0)
    pi = lcd$pi[j]
    theta = lcd$theta[j]
    j1 = c(1,j+1)
    dpk = cpk = lcd$cpk[, c(j,nk), drop=FALSE]
    for(i in 1:nrow(dpk)) dpk[i,] = c(cpk[i,1], diff(cpk[i,]))
    lcd = structure(list(alpha=lcd$alpha, C=lcd$C, theta=theta, pi=pi,
                         lower=lcd$lower, upper=lcd$upper,
                         coef=lcd$coef[,j1,drop=FALSE],
                         fk=lcd$fk[j1], dpk=dpk, cpk=cpk),
                    class = "lcd")
  }
  lcd
}

## Log-likelihood function

logLik.lcd = function(object, x, w=NULL, ...) {
  if(is.null(w)) { xw = x.weight(x); x = xw$x; w = xw$w }
  sum(dlcd(x, object, log=TRUE) * w)
}

## Density function f(x,G,alpha)

dlcd = function(x, lcd, log=FALSE) {
  logd = rep(-Inf, length(x))
  j = x >= lcd$lower & x <= lcd$upper
  xj = x[j]
  knots = c(lcd$lower, lcd$theta)
  jk = indx(xj, knots)
  logd[j] = lcd$coef[1,jk] + lcd$coef[2,jk] * xj - log(lcd$C)
  if(log) logd else exp(logd)
}

## Probability function

plcd = function(q, lcd, lower.tail=TRUE, log.p=FALSE) {
  p = double(length(q))
  j = q >= lcd$lower & q <= lcd$upper
  p[j] = drop(cpx(lcd, q[j], order=0))
  p[q>lcd$upper] = 1
  if(!lower.tail) p = 1 - p
  if(log.p) log(p) else p
}

# knots in lcd must be a subset of x

maxima.gradient = function(lcd, x, w, tol = -Inf) {
  knots = c(lcd$lower, lcd$theta, lcd$upper)
  nk = length(knots)
  index = match(knots, x)
  grad = gradient.lcd(lcd, x, w, "x")
  ii = integer(nk - 1)
  for(i in 1:(nk-1)) {
    ima = which.max(grad[index[i]:index[i+1]])
    ii[i] = index[i] + ima - 1
  }
  theta = x[ii]
  g = grad[ii]
  gmax = max(g)         # maximum gradient over {x_1, ..., x_n}
  j = g > tol & ! theta %in% knots
  list(theta=theta[j], gradient=g[j], gmax=gmax)
}

## Gradient function

## Exact only for theta being a subset of x, for speedy computation.

## If theta = "x", then theta = x is used internally.

## If theta = "knots", then theta = c(lcd$lower,lcd$theta) is used internally.

gradient.lcd = function(lcd, x, w=NULL, theta) {
  if(length(theta) < 1) return(numeric(0))
  if(is.null(w)) { xw = x.weight(x); x = xw$x; w = xw$w }
  knots = c(lcd$lower, lcd$theta)
  nk = length(knots)
  xxt = attr(x, "xx")
  if(is.null(xxt)) xxt = rev(cumsum(rev(w))) * x -  rev(cumsum(rev(x * w)))
  if(!is.numeric(theta) && theta == "knots") {
    xxt = xxt[indx(knots, x)]
    px = cbind(0, lcd$cpk[1:2,-nk,drop=FALSE])
    xxt + sum(w) * (lcd$cpk[2,nk] - px[2,] - (1 - px[1,]) * knots)
  }
  else {
    if(is.numeric(theta)) xxt = xxt[indx(theta, x)]
    else theta = x
    cpx = cpx(lcd, theta, order=1)
    xxt + sum(w) * (lcd$cpk[2,nk] - cpx[2,] - (1 - cpx[1,]) * theta)
  }
}

## integral_lower^x t^o f(t; G, alpha) dt for o = 0:order
## Condition: x in [lcd$lower, lcd$upper]

##### Input: 
##   lcd      an object of class 'lcd'
##   x        values of x
##   order    to be evaluated for 0:order

##### Output: A (order+1) x m matrix for all m intervals
##   Row 1                    # results for o = 0
##   Row 2                    # results for o = 1
##   Row 3                    # results for o = 2

cpx = function(lcd, x, order=0) {
  knots = c(lcd$lower, lcd$theta)
  nk = length(knots)
  k = indx(x, knots)
  a = knots[k]
  dpx = matrix(0, nrow=order+1, ncol=length(x))
  rownames(dpx) = paste0("x", 0:order)
  coef = lcd$coef[,k,drop=FALSE]
  fkk = lcd$fk[k]
  if(any(j <- a != x & coef[2,] == 0)) {      # horizontal segments
    dpx[1,j] = fkk[j] * (x[j] - a[j])
    if(order >= 1) dpx[2,j] = dpx[1,j] * (x[j] + a[j]) * .5
    if(order >= 2) dpx[3,j] = dpx[1,j] * (x[j]^2 + x[j] * a[j] + a[j]^2) / 3
  }
  if(any(j2 <- a != x & coef[2,] != 0)) {     # non-horizontal segments
    x1 = fkk[j2] / coef[2,j2]
    y1 = exp(coef[1,j2] + coef[2,j2] * x[j2] - log(lcd$C)) / coef[2,j2]
    dpx[1,j2] = y1 - x1
    if(order >= 1) dpx[2,j2] = x[j2] * y1 - a[j2] * x1 - dpx[1,j2] / coef[2,j2]
    if(order >= 2)
      dpx[3,j2] = x[j2]^2 * y1 - a[j2]^2 * x1 - 2 * dpx[2,j2] / coef[2,j2]
  }
  if(nk > 1) cbind(0, lcd$cpk[1:(order+1),,drop=FALSE])[,k,drop=FALSE] + dpx
  else dpx 
}

## Plot the gradient function

plotgradient = function(lcd, x, w=NULL, col="blue", knotcol=col, pch=1,
                        lwd=1, lty=1, ...) {
  if(is.null(w)) { xw = x.weight(x); x = xw$x; w = xw$w }
  y = gradient.lcd(lcd, x, w, "x")
  plot(0, 0, xlim=range(x), ylim=range(y, 0), type="n",
       ylab="Gradient", xlab=expression(theta), ...)
  abline(h=0, col='lightgrey', lty=1)
  lines(x, y, col=col, lwd=lwd)
  knots = c(lcd$lower, lcd$theta, lcd$upper)
  points(knots, y[indx(knots, x)], col=knotcol, pch=pch)
}

## plot density or log-density function

plot.lcd = function(x, data, w=NULL, log=FALSE, col="blue", knotcol=col,
                    border="grey", lwd=1, pch=1, main=NULL, breaks=50, ylim=NULL,...) {
  if(!is.null(w)) data = rep(data, w)
  p = sort( c(seq(x$lower, x$upper, len=200), x$theta) )
  sp = c(x$lower, x$theta, x$upper)
  np = length(p)
  d = dlcd(p, lcd=x)
  sd = dlcd(sp, lcd=x)

  if(!missing(data)) h = hist(data, breaks=breaks, plot=FALSE)
  else h = hist(p, breaks=breaks, plot=FALSE)

  if(missing(main)) {
    #if(log) main = "Log density"
    #else main = "Density"
    main=NULL
  }
  if(log) {
    d = log(d)
    
    if (!is.null(ylim)) ylim=ylim
    else ylim = c(min(d), max(d))
    
    if(missing(data)) {
      plot(p, d, xlab="Data", ylab="Log-density", ylim=ylim,
           main=main, type="l", lwd=lwd, col=col, ...)
      base = min(d)
    }
    else{
      logDensity = log(h$density)
      y0 = range(logDensity, finite=TRUE)
      if (!is.null(ylim)) ylim=ylim
      else{
      ylim[1] = min(y0[1], ylim[1])
      ylim[2] = max(y0[2], ylim[2])}
      base = ylim[1] - 0.25 * diff(ylim)
      if(missing(main)) main =  "Log density estimation"
      loghist(data, h=h, main=main, ylim=ylim,
              xlab="Data", base=base, border=border, ...)
      lines(p, d, col=col, lwd=lwd)
    }
    yrange=ylim
    segments(p[1], base, p[1], d[1], col=col, lwd=lwd)
    segments(p[np], base, p[np], d[np], col=col, lwd=lwd)
    points(sp, log(sd), col=knotcol, pch=pch)
  }
  else{
    ###
    if (is.null(ylim)) ylim = c(0, max(d, h$density))
    yrange=(ylim)
    if(!missing(data)) {
      hist(data, breaks=breaks, freq=F, ylim=ylim, xlab="Data",
           main=main, border=border, ...)
      lines(p, d, col=col, lty=1, lwd=lwd)
    }
    else {
      plot(p, d, col=col, type="l", lwd=lwd,
           xlab="Data", ylab="Density", 
           main=main, ylim=ylim,...)
    }
    segments(p[1], 0, p[1], d[1], col=col, lwd=lwd)
    segments(p[np], 0, p[np], d[np], col=col, lwd=lwd)
    points(sp, sd, col=knotcol, pch=pch)
  }
  return(yrange)
}

loghist = function(x, h=NULL, breaks=30, main=NULL, ylim=NULL,
    xlim=NULL, base=NULL, xlab="x", ylab="Log-density", lty=1, border=1) {
  if(is.null(h)) h = hist(x, breaks=breaks, plot=FALSE)
  logd = log(h$density)
  if(is.null(ylim)) ylim = range(logd, finite=TRUE)
  if(is.null(base)) base = ylim[1] - 0.5 * diff(ylim)
  breaks = h$breaks
  if (is.null(xlim)) xlim = range(breaks)
  nd = length(logd)
  plot(h$mids, logd, xlim=xlim, ylim=ylim, type="n", xlab=xlab, ylab=ylab,
       main=main)
  segments(breaks[-(nd+1)], logd, breaks[-1],
           logd, lty=lty, col=border)             # horizontal
  heights = pmax(c(-Inf,logd[-nd],logd[nd]), c(logd[1],logd[-1],-Inf))
  i = heights > -Inf
  segments(breaks[i], heights[i], breaks[i], base,
           lty=lty, col=border)                   # vertical
}

## Print an 'lcd' object

print.lcd = function(x, ...) {
  a = c(alpha=x$alpha, C=x$C)
  print(a, ...)
  if(length(x$theta) > 0) {
    b = cbind(theta=x$theta, pi=x$pi)
    print(b, ...)
  }
  c = c(lower=x$lower, upper=x$upper)
  print(c, ...)
}

x.weight = function(x) {
  n = length(x)
  if(n == 1) return(list(x=x, w=1))
  x = sort(x)
  i = which(diff(x) > 0)
  i.n = c(i, n)
  list(x=x[i.n], w=i.n-c(0, i))
}

lcdmax = function(lcd) {
  # maximal value of lcd
  if (lcd$alpha<=0){
    max_pt=lcd$lower
  }
  else{
    slopes = lcd$alpha-cumsum(lcd$pi)
    if (slopes[length(slopes)]>=0){
      max_pt=lcd$upper
    }
    else{
      max_pt=lcd$theta[length(slopes[slopes>0])+1]
    }
  }
  return(dlcd(max_pt,lcd))
}


rlcd = function(lcd) {
  max_lcd = lcdmax(lcd)
  flag = 0
  while (flag == 0){
    cand = runif(1,lcd$lower,lcd$upper)
    u=runif(1)
    if (u<dlcd(cand,lcd)/max_lcd){
      flag = 1
      y = cand
    }
  }
  return(y)
  
}









