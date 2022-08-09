##################################################
####### predict.gam.cmp
predict.gam.cmp <- function (object, newdata, type = "link", na.action=na.pass, newdata.guaranteed=FALSE,...)
{
  
  if (type != "link" && type != "response" && type !="newdata" ) {
    warning("Unknown type, reset to response.")
    type <- "response"
  }
  
  if (missing(newdata))
    na.act <- object$na.action
  else {
    if (is.null(na.action))
      na.act <- NULL
    else {
      na.txt <- if (is.character(na.action) || is.function(na.action))
        "na.omit" #get.na.action(na.action)
      else "na.pass"
      if (na.txt == "na.pass")
        na.act <- "na.exclude"
      else if (na.txt == "na.exclude")
        na.act <- "na.omit"
      else na.act <- na.action
    }
  }
  nd.is.mf <- FALSE
  yname <- attr(attr(object$terms, "dataClasses"), "names")[attr(object$terms,
                                                                 "response")]
  if (newdata.guaranteed == FALSE) {
    if (missing(newdata)) {
      newdata <- object$model
      new.data.ok <- FALSE
      nd.is.mf <- TRUE
      response <- newdata[[yname]]
    }
    else {
      new.data.ok <- TRUE
      if (is.data.frame(newdata) && !is.null(attr(newdata,
                                                  "terms"))) {
        if (sum(!(names(object$model) %in% names(newdata))))
          stop("newdata is a model.frame: it should contain all required variables\n")
        nd.is.mf <- TRUE
      }
      else {
        resp <- get.var(yname, newdata, FALSE)
        naresp <- FALSE
        if (!is.null(object$family$predict) && !is.null(resp)) {
          if (!is.null(object$pred.formula))
            object$pred.formula <- attr(object$pred.formula,
                                        "full")
          response <- TRUE
          Terms <- terms(object)
          if (is.matrix(resp)) {
            if (sum(is.na(rowSums(resp))) > 0)
              stop("no NAs allowed in response data for this model")
          }
          else {
            if (sum(is.na(resp)) > 0) {
              naresp <- TRUE
              rar <- range(resp, na.rm = TRUE)
              thresh <- rar[1] * 1.01 - rar[2] * 0.01
              resp[is.na(resp)] <- thresh
              newdata[[yname]] <- thresh
            }
          }
        }
        else {
          response <- FALSE
          Terms <- delete.response(terms(object))
          
        }
        allNames.lam <-  all.vars(object$pred.formula)
        allNames.nu <- all.vars(object$pred.formula.nu)
        allNames <- union(allNames.lam,allNames.nu) 
        if (length(allNames) > 0) {
          ff <- reformulate(allNames)
          
          if (sum(!(allNames %in% names(newdata)))) {
            warning("not all required variables have been supplied in  newdata!\n")
          }
          newdata <- eval(model.frame(ff, data = newdata,
                                      na.action = na.act), parent.frame())
          if (naresp)
            newdata[[yname]][newdata[[yname]] <= thresh] <- NA
        }
        na.act <- attr(newdata, "na.action")
        response <- if (response)
          get.var(yname, newdata, FALSE)
        else NULL
      }
    }
  }
  else {
    na.act <- NULL
    new.data.ok = TRUE
    if (!is.null(attr(newdata, "terms")))
      nd.is.mf <- TRUE
    response <- get.var(yname, newdata, FALSE)
  }
  if (new.data.ok) {
    nn <- names(newdata)
    mn <- union(colnames(object$model),colnames(object$model.nu))
    model <- cbind(object$model,object$model.nu)[,nn]
    for (i in 1:length(newdata)) if (nn[i] %in% mn && is.factor(model[,
                                                                      nn[i]])) {
      levm <- levels(model[, nn[i]])
      levn <- if (any(is.na(levm)))
        levels(factor(newdata[[i]], exclude = NULL))
      else levels(factor(newdata[[i]]))
      if (sum(!levn %in% levm) > 0) {
        msg <- paste("factor levels", paste(levn[!levn %in%
                                                   levm], collapse = ", "), "not in original fit",
                     collapse = "")
        warning(msg)
      }
      if (is.matrix(newdata[[i]])) {
        dum <- factor(newdata[[i]], levels = levm, exclude = NULL)
        dim(dum) <- dim(newdata[[i]])
        newdata[[i]] <- dum
      }
      else newdata[[i]] <- factor(newdata[[i]], levels = levm,
                                  exclude = NULL)
    }
    if (type == "newdata")
      return(newdata)
    if (length(newdata) == 1)
      newdata[[2]] <- newdata[[1]]
    if (is.null(dim(newdata[[1]])))
      np <- length(newdata[[1]])
    else np <- dim(newdata[[1]])[1]
    nb <- length(object$coefficients)
  }
  else {
    np <- nrow(object$model)
    nb <- length(object$coefficients)
  }
  
  n.smooth <- length(object$smooth)
  fit <-  array(0, np)
  fitnu <- array(0,np)
  nnu <- length(object$coefficients.nu)
  n.smooth.nu <- length(object$smooth.nu)
  
  stop <- 0
  
  Terms <- list(delete.response(object$pterms))
  Terms.nu <- list(delete.response(object$pterms.nu))
  pterms <- list(object$pterms)
  pstart <- 1
  pind <- 1:object$nsdf
  pind.nu <- 1:object$nsdf.nu
  
  drop.intercept <- object$family$drop.intercept
  if (is.null(drop.intercept)) {
    drop.intercept <- rep(FALSE, length(Terms))
  }
  else {
    for (i in 1:length(Terms)) {
      if (drop.intercept[i] == TRUE)
        attr(Terms[[i]], "intercept") <- 1
    }
  }
  drop.ind <- attr(object$nsdf, "drop.ind")
  s.offset <- NULL
  any.soff <- FALSE
  n.blocks <- 1
  b.size <- array(np, 1)
  if (n.blocks > 0)
    for (b in 1:n.blocks) {
      start <- stop + 1
      stop <- start + b.size[b] - 1
      if (n.blocks == 1)
        data <- newdata
      else data <- newdata[start:stop, ]
      X <- matrix(0, b.size[b], nb + length(drop.ind))
      Xoff <- matrix(0, b.size[b], n.smooth)
      Xnu <- matrix(0, b.size[b], nnu + length(drop.ind))
      Xnuoff <- matrix(0, b.size[b], n.smooth.nu)
      offs <- list()
      for (i in 1:length(Terms)) {
        if (new.data.ok) {
          if (nd.is.mf)
            mf <- model.frame(data, xlev = object$xlevels)
          else {
            mf <- model.frame(Terms[[i]], data, xlev = object$xlevels)
            mf.nu <- model.frame(Terms.nu[[i]], data, xlev = object$xlevels.nu)
            if (!is.null(cl <- attr(pterms[[i]], "dataClasses")))
              .checkMFClasses(cl, mf)
          }
          oc <- if (length(object$contrasts) == 0)
            object$contrasts
          else object$contrasts[names(object$contrasts) %in%
                                  attr(Terms[[i]], "term.labels")]
          Xp <- model.matrix(Terms[[i]], mf, contrasts = oc)
          Xp.nu <- model.matrix(Terms.nu[[i]],mf.nu)
        }
        else {
          Xp <- model.matrix(Terms[[i]], object$model)
          Xp.nu <- model.matrix(Terms.nu[[i]], object$model.nu)
          mf <- newdata
        }
        offi <- attr(Terms[[i]], "offset")
        offs <- list()
        
        if (is.null(offi))
          offs[[i]] <- 0
        else {
          offs[[i]] <- mf[[names(attr(Terms[[i]], "dataClasses"))[offi +
                                                                    1]]]
        }
        #
        offi.nu <- attr(Terms.nu[[i]], "offset")
        offs.nu <- list()
        
        if (is.null(offi.nu))
          offs.nu[[i]] <- 0
        else {
          offs.nu[[i]] <- mf.nu[[names(attr(Terms[[i]], "dataClasses"))[offi.nu +
                                                                          1]]]
        }
        
        #
        if (drop.intercept[i]) {
          xat <- attributes(Xp)
          ind <- xat$assign > 0
          Xp <- Xp[, xat$assign > 0, drop = FALSE]
          xat$assign <- xat$assign[ind]
          xat$dimnames[[2]] <- xat$dimnames[[2]][ind]
          xat$dim[2] <- xat$dim[2] - 1
          attributes(Xp) <- xat
        }
        if (object$nsdf[i] > 0)
          X[, pstart[i] - 1 + 1:object$nsdf[i]] <- Xp
        if(object$nsdf.nu[i] > 0)
          Xnu[,pstart[i]-1+1:object$nsdf.nu[i]] <- Xp.nu
      }
      if (!is.null(drop.ind))
        X <- X[, -drop.ind]
      if (n.smooth)
        for (k in 1:n.smooth) {
          klab <- object$smooth[[k]]$label
          Xfrag <- PredictMat(object$smooth[[k]], data)
          X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
          Xfrag.off <- attr(Xfrag, "offset")
          if (!is.null(Xfrag.off)) {
            Xoff[, k] <- Xfrag.off
            any.soff <- TRUE
          }
          
        }
      #
      if (n.smooth.nu)
        for (k in 1:n.smooth.nu) {
          klab <- object$smooth.nu[[k]]$label
          Xfragnu <- PredictMat(object$smooth.nu[[k]], data)
          Xnu[, object$smooth.nu[[k]]$first.para:object$smooth.nu[[k]]$last.para] <- Xfragnu
          Xfragnu.off <- attr(Xfragnu, "offset")
          if (!is.null(Xfragnu.off)) {
            Xnuoff[, k] <- Xfragnu.off
            any.soff <- TRUE
          }
          
        }
      
      
      #
      lpi <- if (is.list(object$formula)) 
        attr(object$formula, "lpi")
      else NULL
      nlp <- if (is.null(lpi)) 
        1
      else length(lpi)
      fam <- object$family
      k <- attr(attr(object$model, "terms"), "offset")
      knu <-   attr(attr(object$model.nu, "terms"), "offset")
      
      offs <- if (is.null(k))
        rowSums(Xoff)
      else rowSums(Xoff) + model.offset(mf)
      
      offs.nu <- if (is.null(knu))
        rowSums(Xnuoff)
      else rowSums(Xnuoff) + model.offset(mf.nu)
      
      fit[start:stop] <- X %*% object$coefficients +
        offs
      #
      fitnu[start:stop] <- Xnu %*% object$coefficients.nu +
        offs.nu
      #
      fit1 <- NULL
      
      
      if (type == "response") {
        linkinv <- fam$linkinv
        if (is.null(fam$predict)) {
          cumulants <- CMPCumulants(start:stop,llambda = fit,lnu=fitnu,flag = 1)
          fitmean <- cumulants$mean
          
        }
        
        
      }
    }
  if(type == "link") {
    fitmean <- list(fit=fit,fitnu=fitnu, Xmat=X, Xnumat=Xnu)
    #rm(X)
    #rm(Xnu)
    H <- fitmean
    return(H)
  }else
  {
    H <- fitmean
    rn <- rownames(newdata)
    if (is.null(nrow(H))) names(H) <- rn else rownames(H) <- rn
    H <- napredict(na.act, H)
    return(list(fit=H, Xmat=X, Xnumat=Xnu))
  }
}


