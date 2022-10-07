## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(MASS)
library(PanelCount)
set.seed(1)
N = 200
periods = 10
rho = 0.25
tau = 0.5

id = rep(1:N, each=periods)
time = rep(1:periods, N)
x = rnorm(N*periods)
w = rnorm(N*periods)

# correlated random effects at the individual level
r = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), nrow=2))
r1 = rep(r[,1], each=periods)
r2 = rep(r[,2], each=periods)

# correlated error terms at the individual-time level
e = mvrnorm(N*periods, mu=c(0,0), Sigma=matrix(c(1,tau,tau,1), nrow=2))
e1 = e[,1]
e2 = e[,2]

# selection
z = as.numeric(1+x+w+r1+e1>0)
# outcome
y = rpois(N*periods, exp(-1+x+r2+e2))
y[z==0] = NA
sim = data.frame(id,time,x,w,z,y)
head(sim)

## -----------------------------------------------------------------------------
m1 = PoissonRE(y~x, data=sim[!is.na(sim$y), ], id.name='id', verbose=-1)
round(m1$estimates, digits=3)

## -----------------------------------------------------------------------------
m2 = PLN_RE(y~x, data=sim[!is.na(sim$y), ], id.name='id', verbose=-1)
round(m2$estimates, digits=3)

## -----------------------------------------------------------------------------
m3 = ProbitRE(z~x+w, data=sim, id.name='id', verbose=-1)
round(m3$estimates, digits=3)

## -----------------------------------------------------------------------------
m4 = ProbitRE_PoissonRE(z~x+w, y~x, data=sim, id.name='id', verbose=-1)
round(m4$estimates, digits=3)

## -----------------------------------------------------------------------------
# The estimation may take up to 1 minute
m5 = ProbitRE_PLNRE(z~x+w, y~x, data=sim, id.name='id', verbose=-1)
round(m5$estimates, digits=3)

