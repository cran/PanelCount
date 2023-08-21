LL_ProbitRE = function(par,z,w,group,H,verbose=1){
    par = transformToBounded(par, c(delta='exp'))
    alpha = par[1:ncol(w)]
    delta = par['delta']

    N = length(group)-1
    obs = length(z)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    wa = as.vector(w %*% alpha)
    d = 2*z - 1

    Li = rep(0, N)
    for(k in 1:H){
        q = wa + sqrt(2)*delta*v[k]
        Phi = pnorm(d*q)
        Phi_prod = as.vector(groupProd(Phi, group)) # @@@
        Li = Li + Omega[k]*Phi_prod
    }
    LL = sum(log(Li/sqrt(pi)))

    if(verbose>=1){
        cat(sprintf('==== Likelihood function call %d: LL=%.5f ====\n', panel.count.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}


# 2. Gradient function
Gradient_ProbitRE = function(par,z,w,group,H,variance=FALSE,verbose=1){
    par = transformToBounded(par, c(delta='exp'))
    alpha = par[1:ncol(w)]
    delta = par['delta']

    N = length(group)-1
    obs = length(z)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    wa = as.vector(w %*% alpha)
    d = 2*z - 1
    w_ext = cbind(w,delta=0)

    Li = rep(0, N)
    sumK_alp = matrix(0, N, ncol(w_ext))
    for(k in 1:H){
        dq = d * (wa + sqrt(2)*delta*v[k])
        phi = dnorm(dq)
        Phi = pnorm(dq)
        Phi_prod = as.vector(groupProd(Phi, group)) # @@@
        Li = Li + Omega[k]*Phi_prod

        sumT_alp = matVecProdSum(w_ext, sqrt(2)*v[k], d*phi/Phi, group) # @@@
        sumK_alp = sumK_alp + matVecProd(sumT_alp, Omega[k]*Phi_prod)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    dLogLi = matVecProd(sumK_alp, 1/Li) # sqrt(pi) cancels out
    colnames(dLogLi) = names(par)
    # transformation
    dLogLi = transformDerivative(dLogLi, tail(par, 1), c(delta='exp'))
    gradient = colSums(dLogLi)

    if(verbose>=2){
        cat("----Gradient:\n")
        print(gradient,digits=3)
    }
    if(any(is.na(gradient) | !is.finite(gradient))) gradient = rep(NA, length(gradient))
    if(variance){
        var = tryCatch( solve(crossprod(dLogLi)), error = function(e){
            cat('BHHH cross-product not invertible: ', e$message, '\n')
            diag(length(par)) * NA
        } )
        return (list(g=gradient, var=var, I=crossprod(dLogLi)))
    }
    return(gradient)
}

#' A Probit Model with Random Effects
#' @description Estimate a Probit model with random effects at the individual level.
#' \deqn{z_{it}=1(\boldsymbol{\alpha}\mathbf{w_{it}}'+\delta u_i+\xi_{it} > 0)}{z_it=1(\alpha*w_it'+\delta*u_i+\xi_it > 0)}
#' Notations:
#' * \eqn{w_{it}}{w_it}: variables influencing the selection decision \eqn{z_{it}}{z_it}, which could be a mixture of time-variant variables, time-invariant variables, and time dummies
#' * \eqn{u_i}: individual level random effect
#' * \eqn{\xi_{it}}{\xi_it}: error term
#' @md
#' @param formula Formula of the model
#' @param data Input data, a data.frame object
#' @param id.name The name of the column representing id. Data will be sorted by id to improve estimation speed.
#' @param par Starting values for estimates. Default to estimates of Probit model.
#' @param method Optimization method used by optim. Defaults to 'BFGS'.
#' @param se_type Report Hessian or BHHH standard errors. Defaults to Hessian.
#' @param delta Starting value for delta. Defaults to 1 and will be ignored if par is provided.
#' @param H Number of Quadrature points used for numerical integration using the Gaussian-Hermite Quadrature method. Defaults to 20.
#' @param reltol Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
#' @param verbose A integer indicating how much output to display during the estimation process.
#' * <0 - No ouput
#' * 0 - Basic output (model estimates)
#' * 1 - Moderate output, basic ouput + parameter and likelihood in each iteration
#' * 2 - Extensive output, moderate output + gradient values on each call
#' @return A list containing the results of the estimated model, some of which are inherited from the return of optim
#' * estimates: Model estimates with 95% confidence intervals
#' * par: Point estimates
#' * var_bhhh: BHHH covariance matrix, inverse of the outer product of gradient at the maximum
#' * var_hessian: Inverse of negative Hessian matrix (the second order derivative of likelihood at the maximum)
#' * se_bhhh: BHHH standard errors
#' * g: Gradient function at maximum
#' * gtHg: \eqn{g'H^-1g}, where H^-1 is approximated by var_bhhh. A value close to zero (e.g., <1e-3 or 1e-6) indicates good convergence.
#' * LL: Likelihood
#' * AIC: AIC
#' * BIC: BIC
#' * n_obs: Number of observations
#' * time: Time takes to estimate the model
#' * partial: Average partial effect at the population level
#' * paritalAvgObs: Partial effect for an individual with average characteristics
#' * predict: A list with predicted participation probability (prob), predicted potential outcome (outcome), and predicted actual outcome (actual_outcome).
#' * counts: From optim. A two-element integer vector giving the number of calls to fn and gr respectively. This excludes those calls needed to compute the Hessian, if requested, and any calls to fn to compute a finite-difference approximation to the gradient.
#' * message: From optim. A character string giving any additional information returned by the optimizer, or NULL.
#' * convergence: From optim. An integer code. 0 indicates successful completion.
#' Note that the list inherits all the complements in the output of optim. See the documentation of optim for more details.
#' * estimates model estimates with 95% confidence intervals
#' * par point estimates
#' * var_bhhh BHHH covariance matrix, inverse of the outer product of gradient at the maximum
#' * var_hessian Inverse of negative Hessian matrix (the second order derivative of likelihood at the maximum)
#' * se_bhhh BHHH standard errors
#' * g graident function at maximum
#' * LL likelihood
#' * AIC AIC
#' * BIC BIC
#' * n_obs Number of observations
#' * counts A two-element integer vector giving the number of calls to fn and gr respectively. This excludes those calls needed to compute the Hessian, if requested, and any calls to fn to compute a finite-difference approximation to the gradient.
#' * time Time takes to estimate the model
#' * message A character string giving any additional information returned by the optimizer, or NULL.
#' * convergence An integer code. 0 indicates successful completion.
#' Note that the list inherits all the complements in the output of optim. See the documentation of optim for more details.
#' @examples
#' # Use the simulated dataset, in which the true coefficients of x and w are 1.
#' data(sim)
#' res = ProbitRE(z~x+w, data=sim, id.name='id', verbose=-1)
#' res$estimates
#' @export
#' @family PanelCount
#' @references
#' 1. Peng, J., & Van den Bulte, C. (2023). Participation vs. Effectiveness in Sponsored Tweet Campaigns: A Quality-Quantity Conundrum. Management Science (forthcoming). Available at SSRN: <https://www.ssrn.com/abstract=2702053>
#'
#' 2. Peng, J., & Van den Bulte, C. (2015). How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. 2015 International Conference on Information Systems. <https://aisel.aisnet.org/icis2015/proceedings/SocialMedia/24/>
ProbitRE = function(formula, data, id.name, par=NULL, delta=NULL, method='BFGS', se_type=c('Hessian', 'BHHH')[1], H=20, reltol=sqrt(.Machine$double.eps), verbose=0){
    # 1.1 Sort data based on id
    data = data[order(data[, id.name]), ]
    group = c(0,cumsum(table(as.integer(factor(data[, id.name])))))

    # 1.2 Parse w and z
    mf = model.frame(formula, data=data, na.action=NULL, drop.unused.levels=TRUE)
    z = model.response(mf, "numeric")
    w = model.matrix(attr(mf, "terms"), data=mf)

    # 1.3 use probit estimates as initial values for parameters
    if(is.null(par)){
        probit = summary(glm(formula, family=binomial(link='probit'), data=data))
        par = c(probit$coefficients[,1], delta=1)
        if(!is.null(delta)) par['delta'] = delta
    }
    # check parameters
    if(length(par) != ncol(w)+1)  stop(sprintf('Number of parameters incorrect: desired %d, provided %d', ncol(w)+1, length(par)))
    if(names(par)[length(par)] != 'delta') stop('Last parameter must be named delta')
    # transform to unbounded parameters
    if(method != 'L-BFGS-B')
        par = transformToUnbounded(par, c(delta='exp'))

    # 2. Estimation
    resetIter()
    updateMethod(method)
    res = optim(par=par, fn=LL_ProbitRE, gr=Gradient_ProbitRE, method=method, hessian=(se_type=='Hessian'), control=list(reltol=reltol,fnscale=-1), z=z, w=w, group=group, H=H, verbose=verbose)

    # 3. Likelihood, standard error, and p values
    res$n_obs = length(z)
    gvar = Gradient_ProbitRE(res$par,z,w,group,H,variance=TRUE,verbose=verbose-1)
    res = compileResults(res, gvar, se_type, trans_vars=c(delta='delta'), trans_types=c('exp'))

    # 4. Print results
    if(verbose>=0){
        cat(sprintf('==== Probit RE Model converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$counts[1], res$LL, res$gtHg))
        print(res$estimates, digits=3)
        print(res$time <- Sys.time() - panel.count.env$begin)
    }
    return (res)
}
