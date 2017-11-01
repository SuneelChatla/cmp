

gam.cmp <- function(formula.list,data=data,family=cmp(),weights=NULL,subset=NULL,na.action,offset=NULL, method="GCV.Cp",optimizer=c("perf","magic"),control=gam.control(), scale=0,select=FALSE,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=TRUE, paraPen=NULL,G=NULL,in.out=NULL,drop.unused.levels=TRUE,drop.intercept=NULL,...)
{
##

if(is.list(formula.list))
{
  t1 <- as.formula(formula.list[[1]])
  t2 <- as.formula(formula.list[[2]])
}else
{
  t1 <- formula.list
  dv <- all.vars(t1)[1]
  t2 <- as.formula(paste(dv,"~1"))
}
  # fixing family
  if (is.character(family)) family <- eval(parse(text=family))
  if (is.function(family)) family <- family()
  if (is.null(family$family)) stop("family not recognized")
  #
# preparing bases gam set-up
G <- gam.perf1(t1,data = data,family = cmp(),weights=NULL,subset=NULL,na.action,offset=NULL, method="GCV.Cp",optimizer=c("perf","magic"),control=gam.control(), scale=0,select=FALSE,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=FALSE, paraPen=NULL,G=NULL,in.out=NULL,drop.unused.levels=TRUE,drop.intercept=NULL,...)
G1 <- gam.perf1(t2,data = data,family = cmp(),weights=NULL,subset=NULL,na.action,offset=NULL, method="GCV.Cp",optimizer=c("perf","magic"),control=gam.control(), scale=0,select=FALSE,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=FALSE, paraPen=NULL,G=NULL,in.out=NULL,drop.unused.levels=TRUE,drop.intercept=NULL,...)
#
G$conv.tol <- control$mgcv.tol      # tolerence for mgcv
G$max.half <- control$mgcv.half
G1$conv.tol <- control$mgcv.tol      # tolerence for mgcv
G1$max.half <- control$mgcv.half
#
scale=-1
G$sig2<-scale; G1$sig2<-scale
#
G$rS <- mini.roots(G$S,G$off,ncol(G$X),G$rank)
G1$rS <- mini.roots(G1$S,G1$off,ncol(G1$X),G1$rank)

Ssp <- totalPenaltySpace(G$S,G$H,G$off,ncol(G$X))
G$Eb <- Ssp$E       ## balanced penalty square root for rank determination purposes
G$U1 <- cbind(Ssp$Y,Ssp$Z) ## eigen space basis
G$Mp <- ncol(Ssp$Z) ## null space dimension
G$UrS <- list()     ## need penalty matrices in overall penalty range space...
if (length(G$S)>0) for (i in 1:length(G$S)) G$UrS[[i]] <- t(Ssp$Y)%*%G$rS[[i]] else i <- 0
if (!is.null(G$H)) { ## then the sqrt fixed penalty matrix H is needed for (RE)ML
  G$UrS[[i+1]] <- t(Ssp$Y)%*%mroot(G$H)
}
#
Ssp1 <- totalPenaltySpace(G1$S,G1$H,G1$off,ncol(G1$X))
G1$Eb <- Ssp1$E       ## balanced penalty square root for rank determination purposes
G1$U1 <- cbind(Ssp1$Y,Ssp1$Z) ## eigen space basis
G1$Mp <- ncol(Ssp1$Z) ## null space dimension
G1$UrS <- list()     ## need penalty matrices in overall penalty range space...
if (length(G1$S)>0) for (i in 1:length(G1$S)) G1$UrS[[i]] <- t(Ssp1$Y)%*%G1$rS[[i]] else i <- 0
if (!is.null(G1$H)) { ## then the sqrt fixed penalty matrix H is needed for (RE)ML
  G1$UrS[[i+1]] <- t(Ssp1$Y)%*%mroot(G1$H)
}

### calling cmp fit
object <- gam.fit.cmp(G,G1,family=family,control=control,gamma=gamma)
###
##
object$smooth<-G$smooth
object$smooth.nu<-G1$smooth

names(object$edf) <- G$term.names
names(object$edf1) <- G$term.names
names(object$edf.nu) <- G1$term.names
names(object$edf1.nu) <- G1$term.names

if (!is.null(G$P)) { ## matrix transforming from fit to prediction parameterization
  object$coefficients <- as.numeric(G$P %*% object$coefficients)
  object$Vp <- G$P %*% object$Vp %*% t(G$P)
  object$Ve <- G$P %*% object$Ve %*% t(G$P)
  rownames(object$Vp) <- colnames(object$Vp) <- G$term.names
  rownames(object$Ve) <- colnames(object$Ve) <- G$term.names
}
#
if (!is.null(G1$P)) { ## matrix transforming from fit to prediction parameterization
  object$coefficients.nu <- as.numeric(G1$P %*% object$coefficients.nu)
  object$Vp.nu <- G1$P %*% object$Vp.nu %*% t(G1$P)
  object$Ve.nu <- G1$P %*% object$Ve.nu %*% t(G1$P)
  rownames(object$Vp.nu) <- colnames(object$Vp.nu) <- G1$term.names
  rownames(object$Ve.nu) <- colnames(object$Ve.nu) <- G1$term.names
}
#
names(object$coefficients) <- G$term.names
names(object$coefficients.nu) <- G1$term.names

##
if (!is.null(G$L)) {
  object$full.sp <- as.numeric(exp(G$L%*%log(object$sp)+G$lsp0))
  names(object$full.sp) <- names(G$lsp0)
}
if (!is.null(G1$L)) {
  object$full.sp.nu <- as.numeric(exp(G1$L%*%log(object$sp.nu)+G1$lsp0))
  names(object$full.sp.nu) <- names(G1$lsp0)
}
#
names(object$sp) <- names(G$sp)
names(object$sp.nu) <- names(G1$sp)

object$paraPen <- G$pP
object$formula <- G$formula

object$paraPen.nu <- G1$pP
object$formula.nu <- G1$formula
## store any lpi attribute of G$X for use in predict.gam...

object$var.summary <- G$var.summary
object$var.summary.nu <- G1$var.summary

object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
object$cmX.nu <- G1$cmX

object$model<-G$mf # store the model frame
object$model.nu<-G1$mf

object$na.action <- attr(G$mf,"na.action") # how to deal with NA's
object$control <- control
object$terms <- G$terms
object$terms.nu <- G1$terms
#
object$pred.formula <- G$pred.formula
object$pred.formula.nu <- G1$pred.formula

#
attr(object$pred.formula,"full") <- reformulate(all.vars(object$terms))
attr(object$pred.formula.nu,"full") <- reformulate(all.vars(object$terms))

object$pterms <- G$pterms
object$pterms.nu <- G1$pterms

object$assign <- G$assign # applies only to pterms
object$assign.nu <- G1$assign

object$contrasts <- G$contrasts
object$contrasts.nu <- G1$contrasts

object$xlevels <- G$xlevels
object$xlevels.nu <- G1$xlevels

object$offset <- G$offset
object$offset.nu <- G1$offset

if (!is.null(G$Xcentre)) object$Xcentre <- G$Xcentre
if (!is.null(G1$Xcentre)) object$Xcentre.nu <- G1$Xcentre

if (control$keepData) object$data <- data
#
object$df.residual <- nrow(G$X) - sum(object$edf)-sum(object$edf.nu)
#
object$min.edf <- G$min.edf
object$min.edf.nu <- G1$min.edf
#
object$optimizer <- "magic"

class(object) <- c("gam","glm","lm")

names(object$gcv.ubre) <- method
names(object$gcv.ubre.nu) <- method

environment(object$pred.formula) <-
environment(object$terms) <- environment(object$pterms) <- .GlobalEnv

environment(object$pred.formula.nu) <-
  environment(object$terms.nu) <- environment(object$pterms.nu) <- .GlobalEnv

if (!is.null(object$model))  environment(attr(object$model,"terms"))  <- .GlobalEnv
if (!is.null(object$model.nu))  environment(attr(object$model.nu,"terms"))  <- .GlobalEnv

if (!is.null(attr(object$pred.formula,"full"))) environment(attr(object$pred.formula,"full")) <- .GlobalEnv
if (!is.null(attr(object$pred.formula.nu,"full"))) environment(attr(object$pred.formula.nu,"full")) <- .GlobalEnv

### return
object
## end
}

####################################################
##### Function to fit extended gam for cmp
####################################################
gam.fit.cmp <- function (G, G1, start = NULL, etastart = NULL,
                          mustart = NULL,nustart=NULL, family = gaussian(),
                          control = gam.control(),gamma=1,...)
  # fitting function for a gam, modified from glm.fit.
  # note that smoothing parameter estimates from one irls iterate are carried over to the next irls iterate
  # unless the range of s.p.s is large enough that numerical problems might be encountered (want to avoid
  # completely flat parts of gcv/ubre score). In the latter case autoinitialization is requested.
  # oneStep == TRUE causes only a single IRLS step to be taken
{
  intercept<-G$intercept
  conv <- FALSE
  nobs <- NROW(G$y)
  nvars <- NCOL(G$X) # check this needed
  nvars1 <- NCOL(G1$X) # check this needed
  y<-G$y # original data
  X<-G$X # original design matrix
  Xnu <- G1$X
  if (nvars == 0||nvars1 == 0) stop("Model seems to contain no terms")
  olm <- G$am & G1$am   # only need 1 iteration as it's a pure additive model.


  # obtain average element sizes for the penalties
  n.free <-length(G$S)
  n.free1 <-length(G1$S)
  if (n.free>0)
  { S.size<-0
  for (i in 1:n.free) S.size[i]<-mean(abs(G$S[[i]]))
  }
  if (n.free1>0)
  { S.size1<-0
  for (i in 1:n.free1) S.size1[i]<-mean(abs(G1$S[[i]]))
  }
  #
  weights<-G$w # original weights
  n.score <- sum(weights!=0) ## n to use in GCV score (i.e. don't count points with no influence)
  weights1 <- G1$w
  n.score1 <- sum(weights1!=0)
  #
  offset <- G$offset
  offset1 <- G1$offset

  # cmp internal function
  cumulants <- family$cumulants;devf <- family$dev
  linkinv <- family$linkinv;linkfun <- family$linkfun
  if (!is.function(cumulants) || !is.function(devf))
    stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta))
    valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu))
    validmu <- function(mu) TRUE

  ## initialize nu
  if(is.null(nustart))
  {nustart <- rep(0.2,nobs)
  nueta <- linkfun(nustart)}
  else
  {
    if(length(nu)==1)
    {
      nustart <- rep(nu,nobs)
      nueta <- linkfun(nustart)
    }
    else nueta <- linkfun(nustart)
  }
  ##
  if (is.null(mustart))   # new from version 1.5.0
  {mustart <- (y+0.1)^nustart
  eta <- linkfun(mustart)}
  else
  {
    if(length(mustart)==1) mustart <- rep(mustart,nobs)
    eta <- linkfun(mustart)
  }
  #
  if (NCOL(y) > 1)
    stop("y must be univariate unless binomial")

  coefold <- coefold.nu <- NULL                 # 1.5.0
  mu <- linkinv(eta);nu <- linkinv(nueta)

  if (!(validmu(mu) && valideta(eta)))
    stop("Can't find valid starting values: please specify some")
  ##
  y.logfact <- sapply(y,LogFactorial)
  cum <- cumulants(y,eta,nueta,flag=1)
  devold <- devf(y, y.logfact, eta, nueta)
  #
  boundary <- FALSE
  scale<-G$sig2; scale1 <- G1$sig2
  #
  msp<-rep(-1,n.free) # free smoothing parameter vector for magic
  msp1<-rep(-1,n.free1)
  magic.control<-list(tol=G$conv.tol,step.half=G$max.half,maxit=control$maxit, rank.tol=control$rank.tol)

  ##############
  ## Main iteration of P-IRLS starts here
  ##############
  scflag=0
  for (iter in 1:(20*(control$maxit)))
  {
    good <- weights > 0
    var.y <- cum$var[good]
    mean.y <- cum$mean[good]

    if (any(is.na(var.y)))
      stop("NAs in V(y)")
    if (any(var.y == 0))
      stop("0s in V(y)")

    good <- (weights > 0) & (var.y > 0)  # note good modified here => must re-calc each iter
    if (all(!good)) {
      conv <- FALSE
      warning(paste("No observations informative at iteration",
                    iter))
      break
    }

    #####
    ## for lambda
    #####
    z<-G$y <- (eta - offset)[good] + (y - mean.y)[good]/var.y[good]
    w<- sqrt(weights[good] * var.y[good])

    G$w<-w
    G$X<-X[good,]  # truncated design matrix
    if (dim(X)[2]==1) dim(G$X)<-c(length(X[good,]),1) # otherwise dim(G$X)==NULL !!
    ngoodobs <- as.integer(nobs - sum(!good))
    ncols <- as.integer(1)
    # must set G$sig2 to scale parameter or -1 here....
    G$sig2<-scale
    mr <- magic(G$y,G$X,msp,G$S,G$off,L=G$L,lsp0=G$lsp0,G$rank,G$H,matrix(0,0,ncol(G$X)),  #G$C,
                G$w,gamma=gamma,G$sig2,G$sig2<0,
                ridge.parameter=control$irls.reg,control=magic.control,n.score=n.score,nthreads=control$nthreads)
    G$p<-mr$b;msp<-mr$sp;G$sig2<-mr$scale;G$gcv.ubre<-mr$score;

    if (any(!is.finite(G$p))) {
      conv <- FALSE
      warning(paste("Non-finite coefficients at iteration",iter))
      break
    }

    mustart <- G$p
    eta <- drop(X %*% mustart) # 1.5.0
    mu <- linkinv(eta <- eta + offset)
    eta <- linkfun(mu) # force eta/mu consistency even if linkinv truncates

    ############
    ### For nu
    ############
    cum.logy <- cumulants(y,eta,nueta,flag=2)
    mean.logy <- cum.logy$mean
    var.logy <- cum.logy$var

    if (any(is.na(var.logy)))
      stop("NAs in V(logy)")
    if (any(var.logy == 0))
      stop("0s in V(logy)")

    good1 <- (weights1 > 0) & (var.logy > 0)  # note good modified here => must re-calc each iter
    if (all(!good1)) {
      conv <- FALSE
      warning(paste("No observations informative at iteration",
                    iter))
      break
    }
    #
    Xnu1 <- as.matrix(Xnu[good1,]*nu[good1])
    z1 <- G1$y <- nueta[good1]*nu[good1]+(mean.logy-y.logfact)[good1]/var.logy[good1]
    w1 <- sqrt(weights[good1]*var.logy[good1])
    G1$w <- w1
    G1$X <- Xnu1
    # fit2 <- lm.fit(Xnu1*w1,z1*w1)
    #
    # must set G$sig2 to scale parameter or -1 here....
    G1$sig2<-scale1
    mr1 <- magic(G1$y,G1$X,msp1,G1$S,G1$off,L=G1$L,lsp0=G1$lsp0,G1$rank,G1$H,matrix(0,0,ncol(G1$X)),  #G$C,
                G1$w,gamma=gamma,G1$sig2,G1$sig2<0,
                ridge.parameter=control$irls.reg,control=magic.control,n.score=n.score1,nthreads=control$nthreads)
    G1$p<-mr1$b;msp1<-mr1$sp;G1$sig2<-mr1$scale;G1$gcv.ubre<-mr1$score;

    if (any(!is.finite(G1$p))) {
      conv <- FALSE
      warning(paste("Non-finite coefficients at iteration",iter))
      break
    }
    #
    nustart <- G1$p
    nueta <- drop(Xnu %*% nustart) # 1.5.0
    nu <- linkinv(nueta <- nueta + offset1)
    nueta <- linkfun(nu) # force eta/mu consistency even if linkinv truncates
    ##
    good2 <- good & good1
    dev <- devf(y[good2], y.logfact[good2],eta[good2], nueta[good2])
    cat("deviance",dev,"iter",iter,"\n")

    # termination of onemore loop
    if(scflag==1)
    {
      conv = TRUE
      break
    }
    #
    if (control$trace)
      cat("Deviance =", dev, "Iterations -", iter, "\n")
    boundary <- FALSE
    if (!is.finite(dev) || ((dev-devold)/(0.1+devold) > control$epsilon && iter>5)) {
      if (is.null(coefold))
        stop("no valid set of coefficients has been found:please supply starting values",
             call. = FALSE)
      warning("Step size truncated due to divergence",call.=FALSE)
      ii <- 1
      mco <- coefold; mcoo <- coefoold
      nco <- coefold.nu; ncoo <- coefoold.nu
      while (!is.finite(dev)||((dev-devold)/(0.1+devold)> control$epsilon)) {
        if (ii > control$maxit)
          stop("inner loop 1; can't correct step size")
        ii <- ii + 1
        #
        mustart <- (mustart + coefold)/2
        eta.t <- drop(X %*% mustart)
        #
        mco <- (mcoo+mco)/2
        eta <- drop(X %*% mco)

        ##
        nustart <- (nustart+coefold.nu)/2
        nueta.t <- drop(Xnu %*% nustart)
        #
        nco <- (ncoo+nco)/2
        nueta <- drop(Xnu %*% nco)
        #
        dev <- devf(y[good2], y.logfact[good2],eta.t[good2], nueta.t[good2])
      }
      boundary <- TRUE
      if(abs(dev - devold)/(0.1+abs(dev)) < control$epsilon) scflag <- 1
      coef <- mco; coef.nu  <- nco
      if (control$trace)
        cat("Step halved: new deviance =", dev, "\n")
    }
    if (!(valideta(eta) && validmu(mu) && valideta(nueta) && validmu(nu))) {
      warning("Step size truncated: out of bounds.",call.=FALSE)
      ii <- 1
      mco <- coefold; mcoo <- coefoold
      nco <- coefold.nu; ncoo <- coefoold.nu
      while (!(valideta(eta) && validmu(mu)&& valideta(nueta) && validmu(nu))) {
        if (ii > control$maxit)
          stop("inner loop 2; can't correct step size")
        ii <- ii + 1
        ##
        mustart <- (mustart + coefold)/2
        eta.t <- drop(X %*% mustart)
        #
        mco <- (mcoo+mco)/2
        eta <- drop(X %*% mco)

        ##
        nustart <- (nustart+coefold.nu)/2
        nueta.t <- drop(Xnu %*% nustart)
        #
        nco <- (ncoo+nco)/2
        nueta <- drop(Xnu %*% nco)
      }
      boundary <- TRUE
      dev <- devf(y[good2], y.logfact[good2],eta[good2], nueta[good2])
      if(abs(dev - devold)/(0.1+abs(dev)) < control$epsilon) scflag <- 1
      coef <- mco; coef.nu  <- nco
      if (control$trace)
        cat("Step halved: new deviance =", dev, "\n")
    }

    ## Test for convergence here ...
    ccrit <- abs(dev-devold)/(0.1+abs(devold))
    #
    if ( (ccrit < control$epsilon && scflag==0) || olm ) {
      conv <- TRUE
      coef <- mustart #1.5.0
      coef.nu <- nustart
      break
    }else{
      devold <- dev
      coefoold <- coefold ; coefold <- coef<-mustart
      coefoold.nu <- coefold.nu; coefold.nu <- coef.nu <- nustart
      cum <- cumulants(y,eta,nueta,flag=1)
    }

    ### end of for loop
  }

  ##
  if (!conv)
  {
    warning("Algorithm did not converge")
  }
  if (boundary) warning("Algorithm stopped at boundary value")
  eps <- 10 * .Machine$double.eps
  if (any(mu < control$eps)) warning("fitted rates numerically 0 occurred")

  residuals <- rep(NA, nobs)
  residuals[good] <- z - (eta - offset)[good]

  nr <- min(sum(good), nvars)
  nr1 <- min(sum(good1), nvars1)

  wt <- rep(0, nobs)
  wt[good] <- G$w^2

  wt1 <- rep(0, nobs)
  wt1[good1] <- G1$w^2


  # wtdmu <- if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  # nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)-sum(weights1==0)
  # nulldf <- n.ok - as.integer(intercept)
  #####################################
  #### post processing for lambda
  mv<-magic.post.proc(G$X,mr,w=G$w^2)
  G$Vp<-mv$Vb;
  G$hat<-mv$hat;
  G$Ve <- mv$Ve # frequentist cov. matrix
  G$edf<-mv$edf    #c(mv$edf,length(coef.nu))
  G$conv<-mr$gcv.info
  G$sp<-msp
  rank<-G$conv$rank

  #####################################
  #### post processing for nu
  mv1<-magic.post.proc(G1$X,mr1,w=G1$w^2)
  G1$Vp<-mv1$Vb;
  G1$hat<-mv1$hat;
  G1$Ve <- mv1$Ve # frequentist cov. matrix
  G1$edf<-mv1$edf    #c(mv$edf,length(coef.nu))
  G1$conv<-mr1$gcv.info
  G1$sp<-msp1
  rank1<-G1$conv$rank

  # nu.cmatu <- solve(G1$Vp)*G1$sig2  #t(G1$X)%*%diag(w1^2)%*%(G1$X)
  # cov.lnu <- cumulants(y,eta,nueta,flag=3)$mean-(mean.y*mean.logy)
  # cpd1 <- t(G$X)%*%diag(-1*cov.lnu)%*%(Xnu*nu)
  # #
  # ### For nu
  # nu.cmatp <- solve(nu.cmatu -(t(cpd1)%*%mv$Vb%*%cpd1))
  # nu.cmate <- solve(nu.cmatu -(t(cpd1)%*%mv$Ve%*%cpd1))
  #
  # diag(nu.cmatp)[diag(nu.cmatp)<0] <- control$rank.tol
  # diag(nu.cmate)[diag(nu.cmate)<0] <- control$rank.tol
  #
  # ### Conditional cov. for lambda
  #
  # G$Vp <- mv$Vb + mv$Vb%*%cpd1%*%nu.cmatp%*%t(cpd1)%*%mv$Vb/mr$scale
  # G$Ve <- mv$Ve + mv$Ve%*%cpd1%*%nu.cmate%*%t(cpd1)%*%mv$Vb/mr$scale
  # diag(G$Vp)[diag(G$Vp)<0] <- control$rank.tol
  # diag(G$Ve)[diag(G$Ve)<0] <- control$rank.tol
  #
  aic.model <-  dev + 2 * sum(G$edf) +2*sum(G1$edf)
  if (scale < 0) { ## deviance based GCV
    gcv.ubre.dev <- length(y)*dev/(length(y)-gamma*(sum(G$edf)+sum(G1$edf)))^2
  } else { # deviance based UBRE, which is just AIC
    gcv.ubre.dev <- dev/length(y) + 2 * gamma * sum(G$edf)/length(y) - G$sig2
  }
  #
  list(coefficients = as.vector(coef), residuals = residuals, fitted.values = mu, fitted.values.y=mean.y,
       family = family,linear.predictors = eta,linear.predictors.nu = nueta, deviance = dev,coefficients.nu=as.vector(coef.nu),
       nu=nu, iter = iter, weights = wt, weights.nu=wt1, prior.weights = weights,
       y = y, converged = conv,sig2=G$sig2,edf=G$edf,hat=G$hat,edf1=mv$edf1,
       R=mr$R,boundary = boundary,sp = G$sp,nsdf=G$nsdf,Ve=G$Ve,Vp=G$Vp,mgcv.conv=G$conv,
       gcv.ubre=G$gcv.ubre,aic=aic.model,rank=rank,gcv.ubre.dev=gcv.ubre.dev,scale.estimated = (scale < 0),sig2.nu=G1$sig2,edf.nu=G1$edf,hat.nu=G1$hat,edf1.nu=mv$edf1,
       R.nu=mr1$R,boundary = boundary,sp.nu = G1$sp,nsdf.nu=G1$nsdf,Ve.nu=G1$Ve,Vp.nu=G1$Vp,mgcv.conv.nu=G1$conv,
       gcv.ubre.nu=G1$gcv.ubre,rank.nu=rank1)
}
#end of cmpfit ####################################

####################################
######### Summary
####################################
summary.gam.cmp <- function (object, dispersion = NULL,dispersion.nu = NULL, freq = FALSE, ...) {
  ## summary method for gam object - provides approximate p values
  ## for terms + other diagnostics
  ## Improved by Henric Nilsson
  ## * freq determines whether a frequentist or Bayesian cov matrix is
  ##   used for parametric terms. Usually the default TRUE will result
  ##   in reasonable results with paraPen.
  ## If a smooth has a field 'random' and it is set to TRUE then
  ## it is treated as a random effect for some p-value dist calcs


  pinv<-function(V,M,rank.tol=1e-6) {
    ## a local pseudoinverse function
    D <- eigen(V,symmetric=TRUE)
    M1<-length(D$values[D$values>rank.tol*D$values[1]])
    if (M>M1) M<-M1 # avoid problems with zero eigen-values

    if (M+1<=length(D$values)) D$values[(M+1):length(D$values)]<-1
    D$values<- 1/D$values
    if (M+1<=length(D$values)) D$values[(M+1):length(D$values)]<-0
    res <- D$vectors%*%(D$values*t(D$vectors))  ##D$u%*%diag(D$d)%*%D$v
    attr(res,"rank") <- M
    res
  } ## end of pinv

  ########################
  ### for lambda model
  #########################
  if (is.null(object$R)) { ## Factor from QR decomp of sqrt(W)X
    warning("p-values for any terms that can be penalized to zero will be unreliable: refit model to fix this.")
    useR <- FALSE
  } else useR <- TRUE

  p.table <- pTerms.table <- s.table <- NULL

  if (freq) covmat <- object$Ve else covmat <- object$Vp
  name <- names(object$edf)
  dimnames(covmat) <- list(name, name)
  covmat.unscaled <- covmat/object$sig2
  est.disp <- object$scale.estimated
  if (!is.null(dispersion)) {
    covmat <- dispersion * covmat.unscaled
    object$Ve <- object$Ve*dispersion/object$sig2 ## freq
    object$Vp <- object$Vp*dispersion/object$sig2 ## Bayes
    est.disp <- FALSE
  } else dispersion <- object$sig2


  ## Now the individual parameteric coefficient p-values...

  se <- diag(covmat)^0.5
  residual.df <- length(object$y) - sum(object$edf)-sum(object$edf.nu)

  if (sum(object$nsdf) > 0) { # individual parameters
    if (length(object$nsdf)>1) { ## several linear predictors
      pstart <- attr(object$nsdf,"pstart")
      ind <- rep(0,0)
      for (i in 1:length(object$nsdf)) if (object$nsdf[i]>0) ind <-
        c(ind,pstart[i]:(pstart[i]+object$nsdf[i]-1))
    } else { pstart <- 1;ind <- 1:object$nsdf} ## only one lp
    p.coeff <- object$coefficients[ind]
    p.se <- se[ind]
    p.t<-p.coeff/p.se
    if (!est.disp) {
      p.pv <- 2*pnorm(abs(p.t),lower.tail=FALSE)
      p.table <- cbind(p.coeff, p.se, p.t, p.pv)
      dimnames(p.table) <- list(names(p.coeff), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    } else {
      p.pv <- 2*pt(abs(p.t),df=residual.df,lower.tail=FALSE)
      p.table <- cbind(p.coeff, p.se, p.t, p.pv)
      dimnames(p.table) <- list(names(p.coeff), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    }
  } else {p.coeff <- p.t <- p.pv <- array(0,0)}

  ## Next the p-values for parametric terms, so that factors are treated whole...

  pterms <- if (is.list(object$pterms)) object$pterms else list(object$pterms)
  if (!is.list(object$assign)) object$assign <- list(object$assign)
  npt <- length(unlist(lapply(pterms,attr,"term.labels")))
  if (npt>0)  pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0,npt)
  term.labels <- rep("",0)
  k <- 0 ## total term counter
  for (j in 1:length(pterms)) {
    tlj <- attr(pterms[[j]],"term.labels")
    nt <- length(tlj)
    if (j>1 && nt>0) tlj <- paste(tlj,j-1,sep=".")
    term.labels <- c(term.labels,tlj)
    if (nt>0) { # individual parametric terms
      np <- length(object$assign[[j]])
      ind <- pstart[j] - 1 + 1:np
      Vb <- covmat[ind,ind,drop=FALSE]
      bp <- array(object$coefficients[ind],np)

      for (i in 1:nt) {
        k <- k + 1
        ind <- object$assign[[j]]==i
        b <- bp[ind];V <- Vb[ind,ind]
        ## pseudo-inverse needed in case of truncation of parametric space
        if (length(b)==1) {
          V <- 1/V
          pTerms.df[k] <- nb <- 1
          pTerms.chi.sq[k] <- V*b*b
        } else {
          V <- pinv(V,length(b),rank.tol=.Machine$double.eps^.5)
          pTerms.df[k] <- nb <- attr(V,"rank")
          pTerms.chi.sq[k] <- t(b)%*%V%*%b
        }
        if (!est.disp)
          pTerms.pv[k] <- pchisq(pTerms.chi.sq[k],df=nb,lower.tail=FALSE)
        else
          pTerms.pv[k] <- pf(pTerms.chi.sq[k]/nb,df1=nb,df2=residual.df,lower.tail=FALSE)
      } ## for (i in 1:nt)
    } ## if (nt>0)
  }

  if (npt) {
    attr(pTerms.pv,"names") <- term.labels
    if (!est.disp) {
      pTerms.table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)
      dimnames(pTerms.table) <- list(term.labels, c("df", "Chi.sq", "p-value"))
    } else {
      pTerms.table <- cbind(pTerms.df, pTerms.chi.sq/pTerms.df, pTerms.pv)
      dimnames(pTerms.table) <- list(term.labels, c("df", "F", "p-value"))
    }
  } else { pTerms.df<-pTerms.chi.sq<-pTerms.pv<-array(0,0)}

  ## Now deal with the smooth terms....

  m <- length(object$smooth) # number of smooth terms

  df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
  if (m>0) { # form test statistics for each smooth
    ## Bayesian p-values required
    if (useR)  X <- object$R else {
      sub.samp <- max(1000,2*length(object$coefficients))
      if (nrow(object$model)>sub.samp) { ## subsample to get X for p-values calc.
        seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
        if (inherits(seed,"try-error")) {
          runif(1)
          seed <- get(".Random.seed",envir=.GlobalEnv)
        }
        kind <- RNGkind(NULL)
        RNGkind("default","default")
        set.seed(11) ## ensure repeatability
        ind <- sample(1:nrow(object$model),sub.samp,replace=FALSE)  ## sample these rows from X
        X <- predict(object,object$model[ind,],type="lpmatrix")
        RNGkind(kind[1],kind[2])
        assign(".Random.seed",seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
      } else { ## don't need to subsample
        X <- model.matrix(object)
      }
      X <- X[!is.na(rowSums(X)),] ## exclude NA's (possible under na.exclude)
    } ## end if (m>0)

    for (i in 1:m) { ## loop through smooths

      start <- object$smooth[[i]]$first.para;stop <- object$smooth[[i]]$last.para

      V <- object$Vp[start:stop,start:stop,drop=FALSE] ## Bayesian

      p <- object$coefficients[start:stop]  # params for smooth

      edf1[i] <- edf[i] <- sum(object$edf[start:stop]) # edf for this smooth
      ## extract alternative edf estimate for this smooth, if possible...
      if (!is.null(object$edf1)) edf1[i] <-  sum(object$edf1[start:stop])

      Xt <- X[,start:stop,drop=FALSE]
      fx <- if (inherits(object$smooth[[i]],"tensor.smooth")&&
                !is.null(object$smooth[[i]]$fx)) all(object$smooth[[i]]$fx) else object$smooth[[i]]$fixed
      if (!fx&&object$smooth[[i]]$null.space.dim==0&&!is.null(object$R)) { ## random effect or fully penalized term
        res <- reTest(object,i)
      } else { ## Inverted Nychka interval statistics
        df[i] <- min(ncol(Xt),edf1[i])
        if (est.disp) rdf <- residual.df else rdf <- -1
        res <- testStat(p,Xt,V,df[i],type=0,res.df = rdf)
      }
      df[i] <- res$rank
      chi.sq[i] <- res$stat
      s.pv[i] <- res$pval

      names(chi.sq)[i]<- object$smooth[[i]]$label

    }
    if (!est.disp) {
      s.table <- cbind(edf, df, chi.sq, s.pv)
      dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "Chi.sq", "p-value"))
    } else {
      s.table <- cbind(edf, df, chi.sq/df, s.pv)
      dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "F", "p-value"))
    }
  }
  w <- as.numeric(object$prior.weights)
  mean.y <- sum(w*object$y)/sum(w)
  w <- sqrt(w)
  nobs <- nrow(object$model)
  r.sq <- if (inherits(object$family,"general.family")||!is.null(object$family$no.r.sq)) NULL else
    1 - var(w*(as.numeric(object$y)-object$fitted.values))*(nobs-1)/(var(w*(as.numeric(object$y)-mean.y))*residual.df)

  ###########################################
  ######### add nu for cmp
  #############################################
  if (is.null(object$R.nu)) { ## Factor from QR decomp of sqrt(W)X
    warning("p-values for any terms that can be penalized to zero will be unreliable: refit model to fix this.")
    useR <- FALSE
  } else useR <- TRUE

  p.table.nu <- pTerms.table.nu <- s.table.nu <- NULL

  if (freq) covmat.nu <- object$Ve.nu else covmat.nu <- object$Vp.nu
  name.nu <- names(object$edf.nu)
  dimnames(covmat.nu) <- list(name.nu, name.nu)
  covmat.unscaled.nu <- covmat.nu/object$sig2.nu
  est.disp.nu <- object$scale.estimated
  if (!is.null(dispersion.nu)) {
    covmat.nu <- dispersion.nu * covmat.unscaled.nu
    object$Ve.nu <- object$Ve.nu*dispersion.nu/object$sig2.nu ## freq
    object$Vp.nu <- object$Vp.nu*dispersion.nu/object$sig2.nu ## Bayes
    est.disp.nu <- FALSE
  } else dispersion.nu <- object$sig2.nu


  ## Now the individual parameteric coefficient p-values...

  se.nu <- diag(covmat.nu)^0.5
  residual.df <- length(object$y) - sum(object$edf)-sum(object$edf.nu)

  if (sum(object$nsdf.nu) > 0) { # individual parameters
    if (length(object$nsdf.nu)>1) { ## several linear predictors
      pstart.nu <- attr(object$nsdf.nu,"pstart")
      ind.nu <- rep(0,0)
      for (i in 1:length(object$nsdf.nu)) if (object$nsdf.nu[i]>0) ind.nu <-
        c(ind.nu,pstart.nu[i]:(pstart.nu[i]+object$nsdf.nu[i]-1))
    } else { pstart.nu <- 1;ind.nu <- 1:object$nsdf.nu} ## only one lp
    p.coeff.nu <- object$coefficients.nu[ind.nu]
    p.se.nu <- se.nu[ind.nu]
    p.t.nu<-p.coeff.nu/p.se.nu
    if (!est.disp.nu) {
      p.pv.nu <- 2*pnorm(abs(p.t.nu),lower.tail=FALSE)
      p.table.nu <- cbind(p.coeff.nu, p.se.nu, p.t.nu, p.pv.nu)
      dimnames(p.table.nu) <- list(names(p.coeff.nu), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    } else {
      p.pv.nu <- 2*pt(abs(p.t.nu),df=residual.df,lower.tail=FALSE)
      p.table.nu <- cbind(p.coeff.nu, p.se.nu, p.t.nu, p.pv.nu)
      dimnames(p.table.nu) <- list(names(p.coeff.nu), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    }
  } else {p.coeff.nu <- p.t.nu <- p.pv.nu <- array(0,0)}

  ## Next the p-values for parametric terms, so that factors are treated whole...

  pterms.nu <- if (is.list(object$pterms.nu)) object$pterms.nu else list(object$pterms.nu)
  if (!is.list(object$assign.nu)) object$assign.nu <- list(object$assign.nu)
  npt.nu <- length(unlist(lapply(pterms.nu,attr,"term.labels")))
  if (npt.nu>0)  pTerms.df.nu <- pTerms.chi.sq.nu <- pTerms.pv.nu <- array(0,npt)
  term.labels.nu <- rep("",0)
  k <- 0 ## total term counter
  for (j in 1:length(pterms.nu)) {
    tlj <- attr(pterms.nu[[j]],"term.labels")
    nt <- length(tlj)
    if (j>1 && nt>0) tlj <- paste(tlj,j-1,sep=".")
    term.labels.nu <- c(term.labels.nu,tlj)
    if (nt>0) { # individual parametric terms
      np <- length(object$assign.nu[[j]])
      ind <- pstart.nu[j] - 1 + 1:np
      Vb <- covmat.nu[ind,ind,drop=FALSE]
      bp <- array(object$coefficients.nu[ind],np)

      for (i in 1:nt) {
        k <- k + 1
        ind <- object$assign.nu[[j]]==i
        b <- bp[ind];V <- Vb[ind,ind]
        ## pseudo-inverse needed in case of truncation of parametric space
        if (length(b)==1) {
          V <- 1/V
          pTerms.df.nu[k] <- nb <- 1
          pTerms.chi.sq.nu[k] <- V*b*b
        } else {
          V <- pinv(V,length(b),rank.tol=.Machine$double.eps^.5)
          pTerms.df.nu[k] <- nb <- attr(V,"rank")
          pTerms.chi.sq.nu[k] <- t(b)%*%V%*%b
        }
        if (!est.disp)
          pTerms.pv.nu[k] <- pchisq(pTerms.chi.sq.nu[k],df=nb,lower.tail=FALSE)
        else
          pTerms.pv.nu[k] <- pf(pTerms.chi.sq.nu[k]/nb,df1=nb,df2=residual.df,lower.tail=FALSE)
      } ## for (i in 1:nt)
    } ## if (nt>0)
  }

  if (npt.nu) {
    attr(pTerms.pv.nu,"names") <- term.labels.nu
    if (!est.disp) {
      pTerms.table.nu <- cbind(pTerms.df.nu, pTerms.chi.sq.nu, pTerms.pv.nu)
      dimnames(pTerms.table.nu) <- list(term.labels.nu, c("df", "Chi.sq", "p-value"))
    } else {
      pTerms.table.nu <- cbind(pTerms.df.nu, pTerms.chi.sq.nu/pTerms.df.nu, pTerms.pv.nu)
      dimnames(pTerms.table.nu) <- list(term.labels.nu, c("df", "F", "p-value"))
    }
  } else { pTerms.df.nu<-pTerms.chi.sq.nu<-pTerms.pv.nu<-array(0,0)}

  ## Now deal with the smooth terms....

  m.nu <- length(object$smooth.nu) # number of smooth terms

  df.nu <- edf1.nu <- edf.nu <- s.pv.nu <- chi.sq.nu <- array(0, m.nu)
  if (m.nu>0) { # form test statistics for each smooth
    ## Bayesian p-values required
  Xnu <- object$R.nu


    for (i in 1:m.nu) { ## loop through smooths

      start.nu <- object$smooth.nu[[i]]$first.para;stop.nu <- object$smooth.nu[[i]]$last.para

      V.nu <- object$Vp.nu[start.nu:stop.nu,start.nu:stop.nu,drop=FALSE] ## Bayesian

      p.nu <- object$coefficients.nu[start.nu:stop.nu]  # params for smooth

      edf1.nu[i] <- edf.nu[i] <- sum(object$edf.nu[start.nu:stop.nu]) # edf for this smooth
      ## extract alternative edf estimate for this smooth, if possible...
      if (!is.null(object$edf1.nu)) edf1.nu[i] <-  sum(object$edf1.nu[start.nu:stop.nu])

      Xt.nu <- Xnu[,start.nu:stop.nu,drop=FALSE]
      fx.nu <- if (inherits(object$smooth.nu[[i]],"tensor.smooth")&&
                !is.null(object$smooth.nu[[i]]$fx)) all(object$smooth.nu[[i]]$fx) else object$smooth.nu[[i]]$fixed
      ## Inverted Nychka interval statistics
        df.nu[i] <- min(ncol(Xt.nu),edf1.nu[i])
        if (est.disp.nu) rdf.nu <- residual.df else rdf.nu <- -1
        res.df <- testStat(p.nu,Xt.nu,V.nu,df.nu[i],type=0,res.df = rdf.nu)

      df.nu[i] <- res.df$rank
      chi.sq.nu[i] <- res.df$stat
      s.pv.nu[i] <- res.df$pval

      names(chi.sq.nu)[i]<- object$smooth.nu[[i]]$label

    }
    if (!est.disp) {
      s.table.nu <- cbind(edf.nu, df.nu, chi.sq.nu, s.pv.nu)
      dimnames(s.table.nu) <- list(names(chi.sq.nu), c("edf", "Ref.df", "Chi.sq", "p-value"))
    } else {
      s.table.nu <- cbind(edf.nu, df.nu, chi.sq.nu/df.nu, s.pv.nu)
      dimnames(s.table.nu) <- list(names(chi.sq.nu), c("edf", "Ref.df", "F", "p-value"))
    }

  } ## end if (m.nu>0)
  ##########################################
    # return
  ###########################################
    r.sq <- NULL
    ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,s.pv=s.pv,scale=dispersion,r.sq=NULL,family=object$family,formula=object$formula,n=nobs,edf=edf,dispersion=dispersion,pTerms.pv=pTerms.pv,pTerms.chi.sq=pTerms.chi.sq,pTerms.df = pTerms.df, cov.unscaled = covmat.unscaled, cov.scaled = covmat, p.table = p.table,pTerms.table = pTerms.table, s.table = s.table,method=object$method,sp.criterion=object$gcv.ubre,
              rank=object$rank,np=length(object$coefficients), scale=object$sig2,
  #
  p.coeff.nu=p.coeff.nu,se.nu=se.nu,p.t.nu=p.t.nu,p.pv.nu=p.pv.nu,m.nu=m.nu,chi.sq.nu=chi.sq.nu,s.pv.nu=s.pv.nu,scale.nu=dispersion.nu, edf.nu=edf.nu,dispersion.nu=dispersion.nu,pTerms.pv.nu=pTerms.pv.nu,pTerms.chi.sq.nu=pTerms.chi.sq.nu,pTerms.df.nu = pTerms.df.nu, cov.unscaled.nu = covmat.unscaled.nu, cov.scaled.nu = covmat.nu, p.table.nu = p.table.nu,pTerms.table.nu = pTerms.table.nu, s.table.nu = s.table.nu,sp.criterion.nu=object$gcv.ubre.nu,
  rank.nu=object$rank.nu,np.nu=length(object$coefficients.nu),aic=object$aic, scale.nu=object$sig2.nu )


  class(ret)<-"summary.gam.cmp"
  ret
} ## end summary.gam.cmp


print.summary.gam.cmp <- function(x, digits = max(3, getOption("digits") - 3),
                               signif.stars = getOption("show.signif.stars"), ...)
  # print method for gam summary method. Improved by Henric Nilsson
{ print(x$family)
  cat("Formula:\n")

  if (is.list(x$formula)) for (i in 1:length(x$formula)) print(x$formula[[i]]) else
    print(x$formula)
  #
  cat("\n------------------------------------------------------\n")
  cat("For Lambda model:\n")
  cat("------------------------------------------------------\n")
  #
  if (length(x$p.coeff)>0)
  { cat("\nParametric coefficients:\n")
    printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n")

  if(x$m>0)
  { cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
  }

  #
  cat("\n------------------------------------------------------\n")
  cat("For Nu model:\n")
  cat("------------------------------------------------------\n")

  if (length(x$p.coeff.nu)>0)
  { cat("\nParametric coefficients:\n")
    printCoefmat(x$p.table.nu, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n")

  if(x$m.nu>0)
  { cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$s.table.nu, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
  }


  cat("\n------------------------------------------------------\n")
  #

  if (length(x$aic)>0) cat("AIC = ",formatC(x$aic,digits=5,width=8))
  cat("  Scale est.for lambda = ",formatC(x$scale,digits=5,width=8,flag="-"),"\n")
   cat("  Scale est.for nu = ",formatC(x$scale.nu,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
  invisible(x)
} ## print.summary.gam



######
# Addiing mgcv namespace
environment(gam.cmp) <- asNamespace("mgcv")
environment(summary.gam.cmp) <- asNamespace("mgcv")
environment(print.summary.gam.cmp) <- asNamespace("mgcv")
