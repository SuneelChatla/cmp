
###

gam.perf1 <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,offset=NULL,
                method="GCV.Cp",optimizer=c("perf","magic"),control=list(),#gam.control(),
                scale=0,select=FALSE,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=TRUE,
                paraPen=NULL,G=NULL,in.out=NULL,drop.unused.levels=TRUE,drop.intercept=NULL,...) {
  ## Routine to fit a GAM to some data. The model is stated in the formula, which is then
  ## interpreted to figure out which bits relate to smooth terms and which to parametric terms.
  ## Basic steps:
  ## 1. Formula is split up into parametric and non-parametric parts,
  ##    and a fake formula constructed to be used to pick up data for
  ##    model frame. pterms "terms" object(s) created for parametric
  ##    components, model frame created along with terms object.
  ## 2. 'gam.setup' called to do most of basis construction and other
  ##    elements of model setup.
  ## 3. 'estimate.gam' is called to estimate the model. This performs further
  ##    pre- and post- fitting steps and calls either 'gam.fit' (performance
  ##    iteration) or 'gam.outer' (default method). 'gam.outer' calls the actual
  ##    smoothing parameter optimizer ('newton' by default) and then any post
  ##    processing. The optimizer calls 'gam.fit3/4/5' to estimate the model
  ##    coefficients and obtain derivatives w.r.t. the smoothing parameters.
  ## 4. Finished 'gam' object assembled.
  control <- do.call("gam.control",control)
  if (is.null(G)) {
    ## create model frame.....
    gp <- interpret.gam(formula) # interpret the formula
    cl <- match.call() # call needed in gam object for update to work
    mf <- match.call(expand.dots=FALSE)
    mf$formula <- gp$fake.formula
    mf$family <- mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-mf$select <- mf$drop.intercept <-
      mf$gamma<-mf$method<-mf$fit<-mf$paraPen<-mf$G<-mf$optimizer <- mf$in.out <- mf$...<-NULL
    mf$drop.unused.levels <- drop.unused.levels
    mf[[1]] <- quote(stats::model.frame) ## as.name("model.frame")
    pmf <- mf
    mf <- eval(mf, parent.frame()) # the model frame now contains all the data
    if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf,"terms")

    ## summarize the *raw* input variables
    ## note can't use get_all_vars here -- buggy with matrices
    vars <- all.vars(gp$fake.formula[-2]) ## drop response here
    inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

    ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
    if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data)

    dl <- eval(inp, data, parent.frame())
    names(dl) <- vars ## list of all variables needed
    var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
    rm(dl) ## save space

    ## pterms are terms objects for the parametric model components used in
    ## model setup - don't try obtaining by evaluating pf in mf - doesn't
    ## work in general (e.g. with offset)...

    if (is.list(formula)) { ## then there are several linear predictors
      environment(formula) <- environment(formula[[1]]) ## e.g. termplots needs this
      pterms <- list()
      tlab <- rep("",0)
      for (i in 1:length(formula)) {
        pmf$formula <- gp[[i]]$pf
        pterms[[i]] <- attr(eval(pmf, parent.frame()),"terms")
        tlabi <- attr(pterms[[i]],"term.labels")
        if (i>1&&length(tlabi)>0) tlabi <- paste(tlabi,i-1,sep=".")
        tlab <- c(tlab,tlabi)
      }
      attr(pterms,"term.labels") <- tlab ## labels for all parametric terms, distinguished by predictor
    } else { ## single linear predictor case
      pmf$formula <- gp$pf
      pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
      pterms <- attr(pmf,"terms") ## pmf only used for this
    }

    if (is.character(family)) family <- eval(parse(text=family))
    if (is.function(family)) family <- family()
    if (is.null(family$family)) stop("family not recognized")

    if (family$family[1]=="gaussian" && family$link=="identity") am <- TRUE
    else am <- FALSE

    if (!control$keepData) rm(data) ## save space

    ## check whether family requires intercept to be dropped...
    # drop.intercept <- if (is.null(family$drop.intercept) || !family$drop.intercept) FALSE else TRUE
    # drop.intercept <- as.logical(family$drop.intercept)
    if (is.null(family$drop.intercept)) { ## family does not provide information
      lengthf <- if (is.list(formula)) length(formula) else 1
      if (is.null(drop.intercept)) drop.intercept <- rep(FALSE, lengthf) else {
        drop.intercept <- rep(drop.intercept,length=lengthf) ## force drop.intercept to correct length
        if (sum(drop.intercept)) family$drop.intercept <- drop.intercept ## ensure prediction works
      }
    } else drop.intercept <- as.logical(family$drop.intercept) ## family overrides argument

    if (inherits(family,"general.family")&&!is.null(family$presetup)) eval(family$presetup)

    gsname <- if (is.list(formula)) "gam.setup.list" else "gam.setup"

    G <- do.call(gsname,list(formula=gp,pterms=pterms,
                             data=mf,knots=knots,sp=sp,min.sp=min.sp,
                             H=H,absorb.cons=TRUE,sparse.cons=0,select=select,
                             idLinksBases=control$idLinksBases,scale.penalty=control$scalePenalty,
                             paraPen=paraPen,drop.intercept=drop.intercept))

    G$var.summary <- var.summary
    G$family <- family

    if (ncol(G$X)>nrow(G$X)) stop("Model has more coefficients than data")

    G$terms<-terms;
    G$mf<-mf;G$cl<-cl;
    G$am <- am

    if (is.null(G$offset)) G$offset<-rep(0,G$n)

    G$min.edf <- G$nsdf ## -dim(G$C)[1]
    if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim

    G$formula <- formula
    G$pred.formula <- gp$pred.formula
    environment(G$formula)<-environment(formula)
  }

  if (!fit) return(G)

  G$conv.tol <- control$mgcv.tol      # tolerence for mgcv
  G$max.half <- control$mgcv.half # max step halving in Newton update mgcv
  if (scale==0)
  { if (family$family=="binomial"||family$family=="poisson") scale<-1 #ubre
  else scale <- -1 #gcv
  }

  G$sig2<-scale

  if(family$family[1]=="cmp")
  {
    G$rS <- mini.roots(G$S,G$off,ncol(G$X),G$rank)

    Ssp <- totalPenaltySpace(G$S,G$H,G$off,ncol(G$X))
    G$Eb <- Ssp$E       ## balanced penalty square root for rank determination purposes
    G$U1 <- cbind(Ssp$Y,Ssp$Z) ## eigen space basis
    G$Mp <- ncol(Ssp$Z) ## null space dimension
    G$UrS <- list()     ## need penalty matrices in overall penalty range space...
    if (length(G$S)>0) for (i in 1:length(G$S)) G$UrS[[i]] <- t(Ssp$Y)%*%G$rS[[i]] else i <- 0
    if (!is.null(G$H)) { ## then the sqrt fixed penalty matrix H is needed for (RE)ML
      G$UrS[[i+1]] <- t(Ssp$Y)%*%mroot(G$H)
    }
    formula.nu <- as.formula(~1)
    Xnu <- model.matrix(formula.nu,data.frame(G$X))
    object<-gam.fit.cmp1(G,Xnu,family=family,control=control,gamma=gamma,...)

    object$smooth<-G$smooth

    names(object$edf) <- G$term.names
    names(object$edf1) <- G$term.names
    if (!is.null(G$P)) { ## matrix transforming from fit to prediction parameterization
      object$coefficients <- as.numeric(G$P %*% object$coefficients)
      object$Vp <- G$P %*% object$Vp %*% t(G$P)
      object$Ve <- G$P %*% object$Ve %*% t(G$P)
      rownames(object$Vp) <- colnames(object$Vp) <- G$term.names
      rownames(object$Ve) <- colnames(object$Ve) <- G$term.names
    }
    names(object$coefficients) <- G$term.names

  }else object <- estimate.gam(G,method,optimizer,control,in.out,scale,gamma,...)


  if (!is.null(G$L)) {
    object$full.sp <- as.numeric(exp(G$L%*%log(object$sp)+G$lsp0))
    names(object$full.sp) <- names(G$lsp0)
  }
  names(object$sp) <- names(G$sp)
  object$paraPen <- G$pP
  object$formula <- G$formula
  ## store any lpi attribute of G$X for use in predict.gam...
  if (is.list(object$formula)) attr(object$formula,"lpi") <- attr(G$X,"lpi")
  object$var.summary <- G$var.summary
  object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
  object$model<-G$mf # store the model frame
  object$na.action <- attr(G$mf,"na.action") # how to deal with NA's
  object$control <- control
  object$terms <- G$terms
  object$pred.formula <- G$pred.formula
  #
  if(family$family[1]=="cmp") names(object$coefficients.nu) <- dimnames(Xnu)[[2]]
  #
  attr(object$pred.formula,"full") <- reformulate(all.vars(object$terms))

  object$pterms <- G$pterms
  object$assign <- G$assign # applies only to pterms
  object$contrasts <- G$contrasts
  object$xlevels <- G$xlevels
  object$offset <- G$offset
  if (!is.null(G$Xcentre)) object$Xcentre <- G$Xcentre
  if (control$keepData) object$data <- data
  #
  if(family$family[1]=="cmp") object$df.residual <- nrow(G$X) - sum(object$edf)-length(object$coefficients.nu)
  else object$df.residual <- nrow(G$X) - sum(object$edf)
  #
  object$min.edf <- G$min.edf
  if (G$am&&!(method%in%c("REML","ML","P-ML","P-REML"))) object$optimizer <- "magic" else object$optimizer <- optimizer
  object$call <- G$cl # needed for update() to work
  class(object) <- c("gam","glm","lm")
  if (is.null(object$deviance)) object$deviance <- sum(residuals(object,"deviance")^2)
  names(object$gcv.ubre) <- method
  environment(object$formula) <- environment(object$pred.formula) <-
    environment(object$terms) <- environment(object$pterms) <- .GlobalEnv
  if (!is.null(object$model))  environment(attr(object$model,"terms"))  <- .GlobalEnv
  if (!is.null(attr(object$pred.formula,"full"))) environment(attr(object$pred.formula,"full")) <- .GlobalEnv
  object
} ## gam


################

gam.fit.cmp1 <- function (G, Xnu, start = NULL, etastart = NULL,
                          mustart = NULL,nustart=NULL, family = gaussian(),
                          control = gam.control(),gamma=1)
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
  y<-G$y # original data
  X<-G$X # original design matrix
  if(is.null(Xnu))
  {Xnu <- matrix(1,nobs,1)
  colnames(Xnu) <- "Intercept"
  }else Xnu <- as.matrix(Xnu)
  if (nvars == 0) stop("Model seems to contain no terms")
  olm <- G$am   # only need 1 iteration as it's a pure additive model.


  # obtain average element sizes for the penalties
  n.free<-length(G$S)
  if (n.free>0)
  { S.size<-0
  for (i in 1:n.free) S.size[i]<-mean(abs(G$S[[i]]))
  }
  weights<-G$w # original weights
  n.score <- sum(weights!=0) ## n to use in GCV score (i.e. don't count points with no influence)

  offset<-G$offset

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
  {nustart <- rep(0.25,nobs)
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
  scale<-G$sig2
  #
  msp<-rep(-1,n.free) # free smoothing parameter vector for magic
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

    good1 <- (weights > 0) & (var.logy > 0) & (nu>1e-4)  # note good modified here => must re-calc each iter
    if (all(!good1)) {
      conv <- FALSE
      warning(paste("No observations informative at iteration",
                    iter))
      break
    }
    #
    Xnu1 <- as.matrix(Xnu[good1,]*nu[good1])
    z1 <- nueta[good1]*nu[good1]+(mean.logy-y.logfact)[good1]/var.logy[good1]
    w1 <- sqrt(weights[good1]*var.logy[good1])
    fit2 <- lm.fit(Xnu1*w1,z1*w1)
    #
    nustart <- fit2$coefficients
    nueta <- drop(Xnu %*% nustart) # 1.5.0
    nu <- linkinv(nueta <- nueta + offset)
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
  if (any(mu < eps)) warning("fitted rates numerically 0 occurred")

  residuals <- rep(NA, nobs)
  residuals[good] <- z - (eta - offset)[good]

  nr <- min(sum(good), nvars)

  wt <- rep(0, nobs)
  wt[good] <- G$w^2


  # wtdmu <- if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  # nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  # nulldf <- n.ok - as.integer(intercept)
  mv<-magic.post.proc(G$X,mr,w=G$w^2)
  G$Vp<-mv$Vb;
  G$hat<-mv$hat;
  G$Ve <- mv$Ve # frequentist cov. matrix
  G$edf<-mv$edf    #c(mv$edf,length(coef.nu))
  G$conv<-mr$gcv.info
  G$sp<-msp
  rank<-G$conv$rank

  ######

  nu.cmatu <- t(nu*Xnu)%*%diag(w1^2)%*%(Xnu*nu)
  cov.lnu <- cumulants(y,eta,nueta,flag=3)$mean-(mean.y*mean.logy)
  cpd1 <- t(G$X)%*%diag(-1*cov.lnu)%*%(Xnu*nu)
  #
  ### For nu
  nu.cmatp <- solve(nu.cmatu) # -(t(cpd1)%*%mv$Vb%*%cpd1))
  nu.cmate <- solve(nu.cmatu) # -(t(cpd1)%*%mv$Ve%*%cpd1))

  # diag(nu.cmatp)[diag(nu.cmatp)<0] <- control$rank.tol
  # diag(nu.cmate)[diag(nu.cmate)<0] <- control$rank.tol

  ### Conditional cov. for lambda

  G$Vp <- mv$Vb #+ mv$Vb%*%cpd1%*%nu.cmatp%*%t(cpd1)%*%mv$Vb/mr$scale
  G$Ve <- mv$Ve #+ mv$Ve%*%cpd1%*%nu.cmate%*%t(cpd1)%*%mv$Vb/mr$scale
  # diag(G$Vp)[diag(G$Vp)<0] <- control$rank.tol
  # diag(G$Ve)[diag(G$Ve)<0] <- control$rank.tol
  #
  aic.model <-  dev + 2 * sum(G$edf) +2*length(coef.nu)
  if (scale < 0) { ## deviance based GCV
    gcv.ubre.dev <- length(y)*dev/(length(y)-gamma*sum(G$edf))^2
  } else { # deviance based UBRE, which is just AIC
    gcv.ubre.dev <- dev/length(y) + 2 * gamma * sum(G$edf)/length(y) - G$sig2
  }

  list(coefficients = as.vector(coef), residuals = residuals, fitted.values = mu, fitted.values.y=mean.y,
       family = family,linear.predictors = eta, deviance = dev,coefficients.nu=as.vector(coef.nu), nu.cmatp=nu.cmatp,nu.cmate=nu.cmate,
       nu=nu, iter = iter, weights = wt, prior.weights = weights,
       y = y, converged = conv,sig2=G$sig2,edf=G$edf,hat=G$hat,edf1=mv$edf1,
       R=mr$R,boundary = boundary,sp = G$sp,nsdf=G$nsdf,Ve=G$Ve,Vp=G$Vp,mgcv.conv=G$conv,
       gcv.ubre=G$gcv.ubre,aic=aic.model,rank=rank,gcv.ubre.dev=gcv.ubre.dev,scale.estimated = (scale < 0))
}
##################


summary.gam1 <- function (object, dispersion = NULL, freq = FALSE, ...) {
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
  if(object$family$family=="cmp")
  {
    residual.df <- length(object$y) - sum(object$edf)-length(object$coefficients.nu)
  }else residual.df<-length(object$y)-sum(object$edf)
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

  # add nu for cmp
  if(object$family$family=="cmp")
  {
    if(freq) nc.mat <- object$nu.cmate else nc.mat <- object$nu.cmatp
    n.se <- sqrt(diag(nc.mat))
    n.coeff <- object$coefficients.nu
    n.t <- n.coeff/n.se
    n.pv <- 2 * pt(abs(n.t), df = residual.df, lower.tail = FALSE)
    n.table <- cbind(n.coeff, n.se, n.t, n.pv)
    dimnames(n.table) <- list(names(n.coeff), c("Estimate",  "Std. Error", "t value", "Pr(>|t|)"))

    # return
    r.sq <- NULL
    ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,
              s.pv=s.pv,scale=dispersion,r.sq=NULL,family=object$family,formula=object$formula,n=nobs,
              edf=edf,dispersion=dispersion,pTerms.pv=pTerms.pv,pTerms.chi.sq=pTerms.chi.sq,
              pTerms.df = pTerms.df, cov.unscaled = covmat.unscaled, cov.scaled = covmat, p.table = p.table,
              pTerms.table = pTerms.table, s.table = s.table,method=object$method,sp.criterion=object$gcv.ubre,
              rank=object$rank,np=length(object$coefficients),n.table=n.table)
  }else
  {
    dev.expl<-(object$null.deviance-object$deviance)/object$null.deviance
    if (object$method%in%c("REML","ML")) object$method <- paste("-",object$method,sep="")

    ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,chi.sq=chi.sq,
              s.pv=s.pv,scale=dispersion,r.sq=r.sq,family=object$family,formula=object$formula,n=nobs,
              dev.expl=dev.expl,edf=edf,dispersion=dispersion,pTerms.pv=pTerms.pv,pTerms.chi.sq=pTerms.chi.sq,
              pTerms.df = pTerms.df, cov.unscaled = covmat.unscaled, cov.scaled = covmat, p.table = p.table,
              pTerms.table = pTerms.table, s.table = s.table,method=object$method,sp.criterion=object$gcv.ubre,
              rank=object$rank,np=length(object$coefficients))
  }


  class(ret)<-"summary.gam1"
  ret
} ## end summary.gam


print.summary.gam1 <- function(x, digits = max(3, getOption("digits") - 3),
                              signif.stars = getOption("show.signif.stars"), ...)
  # print method for gam summary method. Improved by Henric Nilsson
{ print(x$family)
  cat("Formula:\n")

  if (is.list(x$formula)) for (i in 1:length(x$formula)) print(x$formula[[i]]) else
    print(x$formula)

  if (length(x$p.coeff)>0)
  { cat("\nParametric coefficients:\n")
    printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n")
  ## adding nu table
  if(x$family$family=="cmp")
  {
    cat("\nParametric coefficients for nu:\n")
    printCoefmat(x$n.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    cat("\n")
  }

  ## end of cmp nu

  if(x$m>0)
  { cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
  }
  cat("\n")
  if (!is.null(x$rank) && x$rank< x$np) cat("Rank: ",x$rank,"/",x$np,"\n",sep="")
  if (!is.null(x$r.sq)) cat("R-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5),"  ")
  if (length(x$dev.expl)>0) cat("Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%",sep="")
  cat("\n")
  if (!is.null(x$method)&&!(x$method%in%c("PQL","lme.ML","lme.REML")))
    cat(x$method," = ",formatC(x$sp.criterion,digits=5),sep="")

  cat("  Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
  invisible(x)
} ## print.summary.gam

#####
## Adding mgcv name space

environment(gam.perf1) <- asNamespace("mgcv")
environment(summary.gam1) <- asNamespace("mgcv")
environment(print.summary.gam1) <- asNamespace("mgcv")
