## Estimating HIV Transmission rates 
rm(list = ls())
# have a genealogy reconcstructed from history of transmission
# at a single point in time a random sample of 75 individual
# taken
library(rcolgem)
tree <- read.tree(system.file( 'extdata/sirModel0.nwk', package='rcolgem'))


plot(tree)
library(rjson)
epidata <- fromJSON(file=system.file('extdata/sirModel0.json', package='rcolgem')) 

file.show( system.file('extdata/sirModel0.xml', package='rcolgem'))


# creating a list to store the true parameter values
# going to focus on estimating beta and the nuisance paramter I0
parms_truth <- list( beta = .00020002, gamma = 1, S0 = 9999, t0 = 0 )

sampleTimes <- rep(12, 75)
names(sampleTimes) <- tree$tip.label
bdt <- binaryDatedTree( tree, sampleTimes=sampleTimes)

births <- c(I = 'parms$beta*S*I') # total rate that new infectiosn are generated
deaths <- c(I = 'parms$gamma*I') # rate that lineages are terminated
nonDemeDynamics <- c(S = '-parms$beta*S*I') # state variables not direclty involved in the genealogy 

x0 <- c(I = 1, S = unname(parms_truth$S0))
# this should before the root of the tree 
t0 <- bdt$maxSampleTime-max(bdt$heights) - 1

# calculate the likelihood of the tree 
coalescent.log.likelihood(
  bdt
  , births, deaths, nonDemeDynamics
  , t0, x0
  , parms=parms_truth
  , fgyResolution = 1000
  , integrationMethod = 'rk4')


# Fit the model using the maximum likelihood 
library(bbmle)
obj_fun <- function(lnbeta, lnI0)
{
  beta <- exp(lnbeta)
  I0 <- exp(lnI0)
  parms <- parms_truth
  parms$beta <- beta
  x0 <- c(I=unname(I0), S = unname(parms$S0) )
  mll <- -coalescent.log.likelihood(
    bdt
    , births, deaths, nonDemeDynamics
    ,  t0, x0
    , parms=parms
    , fgyResolution = 1000
    , integrationMethod = 'rk4')
  print(paste(mll, beta, I0))
  mll
}

# Fitting the model
fit <- mle2(
  obj_fun
  , start=list(lnbeta=log(parms_truth$beta*.75), lnI0=log(1))
  , method='Nelder-Mead'
  , optimizer='optim' 
  , control=list(trace=6, reltol=1e-8)
)


# Comparing the fitted model to the true number of infected through time
AIC(fit)
logLik(fit)
exp(coef(fit))


exp(coef(fit)['lnbeta']) - parms_truth$beta

## Testing the fit model to the truth 
beta <- exp(coef(fit)['lnbeta'])
I0 <- exp(coef(fit)['lnI0'])
parms <- parms_truth
parms$beta <- beta


x0 <- c(I=unname(I0), S = unname(parms$S0))

# Solving model with estimated parameters
o <- solve.model.unstructured(t0,bdt$maxSampleTime, x0
                              , births
                              ,  deaths
                              , nonDemeDynamics, parms)

# What the trajectory actually would look like
otruth <- solve.model.unstructured(t0, bdt$maxSampleTime, x0
                                   , births
                                   ,  deaths
                                   , nonDemeDynamics, parms_truth)

plot(epidata$t, epidata$I, type='line'
     , ylim=c(0, 100+max(max(o[,2]),max(epidata$I))))
lines(o[,1], o[,2], col='red' )
lines(otruth[,1], otruth[,2], col='blue' )











# keeping track of birth and migration events between demes

INFECTEDNAMES <- c('I0', 'I1', 'I2')

# the birth events between demes expressed in a 3x3 matrix F
# there are zero rates in the 2nd and third columns since all new 
# hosts start  out in the first infection


# Births, migrations, and deaths are all part of the 
# coalescent model 
births <- rbind(
    c('beta0*S*I0 / (S+I0+I1+I2)', '0','0'),
    c('beta1*S*I1/(S+I0+I1+I2)', '0', '0'),
    c('beta2*S*I2/(S+I0+I1+I2)', '0', '0')
)
rownames(births) = colnames(births) <- INFECTEDNAMES
births

# This matrix tells the progression from the different
# stages of infection
migrations <- rbind(
  c('0', 'gamma0 * I0', '0'),
  c('0', '0', 'gamma1 *I1'),
  c('0','0','0')
)
rownames(migrations) = colnames(migrations) <- INFECTEDNAMES


# What also terminates a lineage

deaths <- c(
  'mu*I0',
  'mu*I1',
  'mu*I2 + gamma2*I2'
)
names(deaths) <- INFECTEDNAMES


nonDemeDynamics <- c(S = '-mu*S + mu*(S+I0+I1+I2)-
S*(beta0*I0+beta1*I1+beta2*I2)/(S+I0+I1+I2)'
                                     )

# This creates the model and we need to provide the parameters
# The SDE treats the equations as rates within a stochastic 
# differentaiil equation
demo.model <- build.demographic.process (
  births
  , nonDemeDynamics
  , migrations = migrations
  , deaths = deaths
  , parameterNames = c(
    'beta0'
    , 'beta1'
    , 'beta2'
    , 'gamma0'
    , 'gamma1'
    , 'gamma2'
    , 'mu')
  , rcpp = TRUE
  , sde = TRUE
)

# not you an simulate the demogrpahic process

# theta name vectors of the parameters
theta <- c(gamma0 = 1
           , gamma1 = 1/7
           , gamma2 = 1/2
           , mu = 1/30
           , beta0 = 12./10
           , beta1 = 3./100
           , beta2 = 9./100
           )

# Time of inital and ending time conditions
t0 <- 0
t1 <- 50

#x0 named vector of the initial conditions

x0 <- c(S = 999, I0 = 1, I1 = .1, I2 = .1)

# can be modified not to be matplot, change to ggplot
show.demographic.process(demo.model, theta, x0, t0, t1)


# Specifying the time that each lineage is sampled and the state of the lineage
# at the time of sampling

# Vector of uniformly spaced sample times

n <- 100
sampleTimes <-seq(15, 25, length.out = n)


# constructing the state: each element corresponds to the probability
# that lineage i is in deme j 

sampleStates <- t(rmultinom(n, size = 1, prob = c(0.25, .9, .075)))
head(sampleStates)

tree <- sim.co.tree(theta, demo.model, x0, t0, sampleTimes, sampleStates, res = 1e3)

tree



## Vignette 2! 
### Estimating HIV transmission rates given a structured coalescent model
# and knowing the stage of infection when patients are sampled

# Creating list of the parameter values

parms_truth <- list(gamma0 = 1
                     , gamma1 = 1/7
                     , gamma2 = 1/2
                     , mu = 1/30
                     , b = .036
                     , beta0 = 12./10
                     , beta1 = 3./100
                     , beta2 = 9./100
                     , S_0 = 3000
                     , I0_0= 1, I1_0 = 0.01, I2_0=0.01
                     , m = 3, mm = 1)


INFECTEDNAMES <- c('I0', 'I1', 'I2')

births <- rbind(
  c('parms$beta0 *S*I0/(S+I0+I1+I2)','0','0'),
  c('parms$beta1*S*I1/(S+I0+I1+I2)', '0','0'),
  c('parms$beta2*S*I2/(S+I0+I1+I2)', '0', '0')
)
rownames(births) = colnames(births) <- INFECTEDNAMES

migrations <- rbind(
   c('0', 'parms$gamma0 * I0', '0'),
   c('0', '0', 'parms$gamma1 * I1'),
   c('0', '0', '0')
   )
rownames(migrations)=colnames(migrations) <- INFECTEDNAMES

deaths <- c(
   'parms$mu*I0'
   , 'parms$mu*I1'
   , 'parms$mu*I2 + parms$gamma2 * I2'
   )

names(deaths) <- INFECTEDNAMES

nonDemeDynamics <- paste(sep='',
                            '-parms$mu*S + parms$mu*(S + I0 + I1 + I2)'
                            ,'- S * (parms$beta0*I0+parms$beta1*I1+parms$beta2*I2) / (S + I0 + I1 + I2)'
                            )
names(nonDemeDynamics) <- 'S'

# This model will be fitted to a binary tree with dated tips
tree <- read.tree(system.file('extdata/hivSimulation.nwk', package = 'rcolgem'))

# All tips were sampled at the same time 50 years ago
# Have info on the tip label (individual) and the time it was sampled
sampleTimes <- rep(50, length(tree$tip.label))
names(sampleTimes) <- tree$tip.label


#creating a tree with date tips and internal nodes 
# inferring sampe states from tip labels
bdt <- binaryDatedTree(tree
                       , sampleTimes
                       , sampleStatesAnnotations=INFECTEDNAMES)


coalescent.log.likelihood(bdt
                          , births, deaths, nonDemeDynamics
                          , t0 = 0
                          , x0=c(I0=1, I1=.01, I2=.01, S = parms_truth$S_0)
                          , migrations = migrations
                          , parms=parms_truth
                          , fgyResolution=1000
                          , integrationMethod='euler'
                          )

library(bbmle)

# Creating objective function to be minimized
# Focus on estimating three transmission rates of the system
# nuisance parameter that controls initial conditions, t0
# time of the origin of the epidemic

# The objective function uses log-transformation for variables 
# that must be positive 
obj_fun <- function(lnbeta0, lnbeta1, lnbeta2, t0)
{
  parms <- parms_truth
  parms$beta0 <- exp(lnbeta0)
  parms$beta1 <- exp(lnbeta1)
  parms$beta2 <- exp(lnbeta2)
  mll <- -coalescent.log.likelihood(bdt
                                    , births, deaths, nonDemeDynamics
                                    , t0 = 0
                                    , x0=c(I0=1, I1=.01, I2=.01, S = parms_truth$S_0)
                                    , migrations = migrations
                                    , parms=parms_truth
                                    , fgyResolution=1000
                                    , integrationMethod='rk4'
  )
  
  # track progress
  print(c(mll, exp(c(lnbeta0, lnbeta1, lnbeta2) ), t0))
  mll
}

# Fitting the model 
fit <- mle2(obj_fun
               , start=list(lnbeta0=log(.6), lnbeta1=log(.2), lnbeta2=log(.05), t0=0)
               , method='Nelder-Mead', optimizer='optim'
               , control=list(trace=6, reltol=1e-8))

## How well did the optimizer work?
AIC(fit)
logLik(fit)
coef(fit)
exp(coef(fit))


#Compare the fitted model to the true number
profbeta <- profile(fit, which = 1, alpha  = .05, std.err = .5,
                    trace = TRUE, tol.newmin=1)

c(exp(confint(profbeta)), TrueVal = parms_truth$beta0)
plot(profbeta)
