LL_PoissonRE = function(par,y,x,group,H,verbose=1){
    if(length(par) != ncol(x)+1)  stop("Number of parameters incorrect")
    beta = par[1:ncol(x)]
    sigma = par[length(par)]

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    xb = as.vector(x %*% beta)
    Li = rep(0, N)
    for(h in 1:H){
        lamh = exp(xb + sqrt(2)*sigma*v[h])
        ix = !is.na(y)
        Pz = rep(1, obs)
        Pz[ix] = dpois(y[ix], lamh[ix])
        Pz_prod = as.vector(groupProd(Pz, group))
        Li = Li + Omega[h] * Pz_prod
    }
    LL = sum(log(Li/sqrt(pi)))
    if(verbose>=1){
        writeLines(paste("==== Iteration ", tmp.env$iter, ": LL=",round(LL,digits=5)," =====", sep=""))
        print(round(par,digits=3))
    }
    tmp.env$iter = tmp.env$iter + 1
    if(is.na(LL) || !is.finite(LL)){
        if(verbose>=2) writeLines("NA or infinite likelihood, will try others")
        LL = -1e300
    }
    return (LL)
}


# 2. Gradient function
Gradient_PoissonRE = function(par,y,x,group,H,variance=F,verbose=1){
    if(length(par) != ncol(x)+1)  stop("Number of parameters incorrect")
    beta = par[1:ncol(x)]
    sigma = par[length(par)]

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    xb = as.vector(x %*% beta)
    x_ext = cbind(x,sigma=0)

    Li = rep(0, N)
    sumH_beta = matrix(0, N, ncol(x_ext))
    for(h in 1:H){
        lamh = exp(xb + sqrt(2)*sigma*v[h])
        ix = !is.na(y)
        Pz = rep(1, obs)
        Pz[ix] = dpois(y[ix], lamh[ix])
        Pz_prod = as.vector(groupProd(Pz, group))
        Li = Li + Omega[h] * Pz_prod

        zyl = rep(0, obs)
        zyl[ix] = y[ix] - lamh[ix]
        sumT_beta = matVecProdSum(x_ext, sqrt(2)*v[h], zyl, group)
        sumH_beta = sumH_beta + matVecProd(sumT_beta, Omega[h]*Pz_prod)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    # whether rho is included depends on length of par, pi cancells out
    dLogLi = matVecProd(sumH_beta, 1/Li)
    colnames(dLogLi) = names(par)
    gradient = colSums(dLogLi)
    if(variance){
        if(verbose>=1){
            writeLines("----Converged, the gradient at optimum:")
            print(round(gradient,digits=3))
        }
        var = solve(crossprod(dLogLi))
        return (list(g=gradient, var=var, I=crossprod(dLogLi)))
    }
    if(verbose>=2){
        writeLines("----Gradient:")
        print(round(gradient,digits=3))
    }
    if(any(is.na(gradient) | !is.finite(gradient))){
        if(verbose>=2) writeLines("NA or infinite gradient, reset to all -1")
        gradient = rep(-1, length(gradient))
    }
    return (gradient)
}

#' A Poisson Model with Random Effects
#' @description Estimate a Poisson model with random effects in panel counting data. Note this model is different with the Poisson Lognormal model for counting data.
#' @param formula Formula of the model
#' @param data Input data, a data frame
#' @param id A vector that represents the identity of individuals, numeric or character
#' @param method Searching algorithm, don't change default unless you know what you are doing
#' @param par Starting values for estimates
#' @param sigma Variance of random effects on the individual level
#' @param lower Lower bound for estiamtes
#' @param upper Upper bound for estimates
#' @param H A vector of length 2, specifying the number of points for inner and outer Quadratures
#' @param accu 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim
#' @param verbose Level of output during estimation. Lowest is 0.
#' @return A list containing the results of the estimated model
#' @examples
#' data(rt)
#' est = PoissonRE(num.words~fans+tweets+as.factor(tweet.id),
#'                      id=rt$user.id[rt$isRetweet==1],
#'                      data=rt[rt$isRetweet==1,])
#' @export
#' @family PanelCount
PoissonRE = function(formula, id, data=NULL, par=NULL, sigma=1, method='BFGS',lower=NULL,upper=NULL,H=20,accu=1e10,verbose=0){
    # 1.1 Sort data based on id
    ord = order(id)
    data = data[ord,]
    id = id[ord]
    group = c(0,cumsum(table(as.integer(factor(id)))))
    # 1.2 Parse x and y
    mf = model.frame(formula, data=data, na.action=NULL, drop.unused.levels=T)
    y = model.response(mf, "numeric")
    x = model.matrix(attr(mf, "terms"), data=mf)
    if(length(id)!=nrow(x)) stop("Error: id and x length don't match, potentially due to removal of data points with missing values!")
    # 1.3 use Poisson estimates as initial values for parameters
    if(is.null(par)){
        writeLines("========Initializing estimates with Possion model=======")
        psn = summary(glm(formula, family=poisson(link="log"), data=data))
        par = c(psn$coefficients[,1], sigma=sigma)
        if(verbose>=2){
            writeLines('Initial parameters')
            print(par)
        }
    }
    # 2. Estimation
    # Typical values for factr are: 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy.
    tmp.env <<- new.env()
    tmp.env$iter = 1
    begin = Sys.time()
    if(method=="L-BFGS-B"){
        res = optim(par=par, fn=LL_PoissonRE, gr=Gradient_PoissonRE, method="L-BFGS-B", control=list(factr=accu,fnscale=-1), lower=lower, upper=upper, y=y, x=x, group=group, H=H, verbose=verbose)
    } else {
        res = optim(par=par, fn=LL_PoissonRE, gr=Gradient_PoissonRE, method="BFGS", control=list(factr=accu,fnscale=-1), y=y, x=x, group=group, H=H, verbose=verbose)
    }
    # 3. Likelihood, standard error, and p values
    res$LL = LL_PoissonRE(res$par,y,x,group,H,verbose=verbose-1)
    gvar = Gradient_PoissonRE(res$par,y,x,group,H,variance=T,verbose=verbose-1)
    res$var = gvar$var
    res$g = gvar$g
    res$gtHg = matrix(res$g, nrow=1) %*% res$var %*% matrix(res$g, ncol=1)
    res$se = sqrt(diag(res$var))
    res$z = res$par/res$se
    res$p = 1 - pchisq(res$z^2, 1)
    res$estimates = cbind(estimates=res$par,se=res$se,z=res$z,p=res$p)
    writeLines(paste("\n *** Estimation of Poisson RE model finished, LL=",res$LL," ***", sep=""))
    print(res$estimates)
    print(Sys.time()-begin)
    # 4. Convergence criterions
    res$scgrad=tryCatch(solve(chol(gvar$I),gvar$g), error=function(e)e)
    if('numeric' %in% class(res$scgrad)[1]) {
        max_grad = round(max(abs(res$scgrad)), digits=3)
        max_grad_var = names(which.max(abs(res$scgrad)))
        if(verbose>=1)
            writeLines(paste0('Max absolute scaled gradient: ', max_grad, ' on ', paste(max_grad_var, collapse=', ')))
    } else {
        writeLines('Error occured while computing scaled gradient, details below:')
        print(res$scgrad)
    }
    if(verbose>=1) writeLines(paste0('Convergence criterion gtHg: ', round(res$gtHg, digits=3)))
    rm(tmp.env)
    return (res)
}
