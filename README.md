# CMP Distribution
Conway-Maxwell-Poisson GLM and GAM

The package contains two routines to fit CMP generalized additive models. They only differ interms of the dispersion parameter. 

1.  gam.cmp(): allows models for both lambda and nu parameters
    main inputs: 
    - formula: could be a list of two formulas for both lambda and nu or just a single formula for lambda and in this case
    the parameter nu is treated as constant.
    - dataset
    - family=cmp
    
2.  gam.perf1(): allows only a constant dispersion parameter.
    - formula: allows only a single formula for lambda model and nu will be treated as constant.
    - dataset
    - family=cmp or negbin or poisson
    
    For the remaining parameters, the default values are good enough to get the desired results.
    
If the formulas do not contain smooth terms (s()), then the model fitted is simply a generalized linear model.

For more details about the model formulations and estimation procedures, please refer to 

**Suneel Babu Chatla, Galit Shmueli**, Efficient estimation of COM–Poisson regression and a generalized additive model,
*Computational Statistics & Data Analysis*, Volume 121, 2018, Pages 71-88, ISSN 0167-9473, https://doi.org/10.1016/j.csda.2017.11.011.
(http://www.sciencedirect.com/science/article/pii/S0167947317302608)

## Example:

###  Data simulation example
  `set.seed(123)`\
  `n=200`\
  `sdata=data.frame(x0 = runif(n, 0, 1))`\
  `sdata$x1 <- runif(n, 0, 1)`\
  `sdata$x2 <- runif(n, 0, 1)`\
  `sdata$x3 <- runif(n, 0, 1)`
 
  `f0 <- function(x) sin(pi * x)`\
  `f1 <- function(x) exp( x)`\
  `f2 <- function(x) 0.02 * x^2 * ( (1 - x)) + (0.5 * x)^2 * (1 - x)^3`\
  `f3 <- function(x) 2*x-(x^2)`

  `sdata$f <- 2*f3(sdata$x3) +1*f1(sdata$x1) +1* f2(sdata$x2)`\
  `lambda=exp(sdata$f)`\
  `nu= exp(f0(sdata$x0))`\
  `s=rep(0,n)`\
  `y=.C("cmpsim_all",lambda=as.double(lambda),nu=as.double(nu),n=as.integer(n),y=as.double(s))$y`\
  `sdata$y=y`

## Model formulations
`m1 <- as.formula(y~s(x3)+s(x1)+s(x2))`\
`m21 <- as.formula(y~s(x0))`\
`m22 <- as.formula(y~x0)`

  ### CMP GAM with additive models for both lambda and nu
  `cmpgam <- gam.cmp(list(m1,m21),data = sdata,family = cmp)`\
  `summary.gam.cmp(cmpgam)`
  
  
  ### CMP GAM with additive model for lambda and linear model for nu
  `cmpglm <- gam.cmp(list(m1,m22),data = sdata,family = cmp)`\
  `summary.gam.cmp(cmpglm)`
  
  
  ### CMP GAM with additive model for lambda and constant model for nu
  `cmpglm1 <- gam.perf1(m1,data = sdata,family = cmp)`\
  `summary.gam1(cmpglm1)`
  
  ###  Negative Binomial GAM
  `nbgam <- gam.perf1(m1,data=sdata,family = negbin(c(1,10000)))`\
  `summary.gam1(nbgam)`
