LL_PLN_RE = function(par,y,x,group,H,verbose=1){
    par = transformToBounded(par, c(sigma='exp', gamma='exp'))
    beta = par[1:ncol(x)]
    sigma = par['sigma']
    gamma = par['gamma']

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    xb = as.vector(x %*% beta)

    Li = rep(0, N)
    for(h in 1:H){
        sumK = rep(0, obs)
        for(k in 1:H){
            lamh = exp(xb + sqrt(2)*sigma*v[h] + sqrt(2)*gamma*v[k])
            PhiP = Omega[k] / sqrt(pi) * dpois(y, lamh)
            sumK = sumK + PhiP
        }
        prodT = groupProd(sumK, group)
        Li = Li + Omega[h] * prodT
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    LL = sum(log(Li/sqrt(pi)))

    if(verbose>=1){
        cat(sprintf('==== Likelihood function call %d: LL=%.5f ====\n', panel.count.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}

integrate_PLN_RE_LL = function(y, xb, group, sigma, gamma, mu, scale, H){
    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    N = length(group)-1
    obs = length(y)

    Li = rep(0, N)
    for(h in 1:H){
        tau = sqrt(2)*v[h]*scale + mu
        tau_rep = rep(tau, times=diff(group)) # i,t level

        # product of probabilities by individual
        sumK = rep(0, obs)
        for(k in 1:H){
            lam = exp(xb + sigma*tau_rep + sqrt(2)*gamma*v[k])
            PhiP = Omega[k] / sqrt(pi) * dpois(y, lam)
            sumK = sumK + PhiP
        }
        prodT = as.vector(groupProd(sumK, group))
        Li = Li + Omega[h] * exp(v[h]^2) * dnorm(tau) * prodT
    }
    Li = pmax(sqrt(2)*scale*Li, 1e-100)
    sum(log(Li))
}

LL_PLN_RE_AGQ = function(par,y,x,group,H,verbose=1){
    par = transformToBounded(par, c(sigma='exp', gamma='exp'))
    beta = par[1:ncol(x)]
    sigma = par['sigma']
    gamma = par['gamma']

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    xb = as.vector(x %*% beta)

    # update mode and scale
    mu = panel.count.env$mu
    scale = panel.count.env$scale
    n_count = 1
    # Usually converge in a few iterations
    while(panel.count.env$stopUpdate==FALSE){
        Li = rep(0, N)
        mu_new = rep(0, N)
        scale_new = rep(0, N)
        for(h in 1:H){
            tau = sqrt(2)*v[h]*scale + mu
            tau_rep = rep(tau, times=diff(group)) # i,t level

            # product of probabilities by individual
            sumK = rep(0, obs)
            for(k in 1:H){
                lam = exp(xb + sigma*tau_rep + sqrt(2)*gamma*v[k])
                PhiP = Omega[k] / sqrt(pi) * dpois(y, lam)
                sumK = sumK + PhiP
            }
            prodT = as.vector(groupProd(sumK, group))

            # scale is a vector, calculate products of scale first
            Lih = sqrt(2)*Omega[h]*exp(v[h]^2) * scale * dnorm(tau) * prodT
            Li = Li + Lih
            mu_new = mu_new + tau*Lih
            scale_new = scale_new + tau^2*Lih
        }
        Li = pmax(Li, 1e-100)
        mu_new = mu_new / Li
        scale_new = sqrt(pmax(scale_new / Li - mu_new^2, 1e-6))

        # Update only when LL improves
        LL = sum(log(Li))
        if(is.na(LL) || LL < panel.count.env$LL) break

        # Stop on error (takes longer than Poisson RE to converge)
        # Stop on small n_count threshold can dramatically improve speed and has no impact on converged mu
        if(any(is.na(scale_new)) || any(is.na(mu_new)) || any(scale_new<0) || n_count >= 20){
            if(verbose>=2){
                print(sprintf('--- Failing to get reliable adaptive parameters after %d inner iteration ----', n_count))
                print(par)
                print(max(abs(scale_new - scale)/abs(scale)))
                print(max(abs(mu_new - mu)/abs(mu)))
                print(min(scale_new))
                # print(mu_new)
                # print(scale_new)
            }
            break
        }

        # Check if mu and scale converged
        if(all(abs(mu_new - mu) < pmax(1e-2*abs(mu), 1e-3))
           && all(abs(scale_new - scale) < pmax(1e-2*abs(scale), 1e-3))){
            # Without gradient, LL has to be numerically differentiated
            # update only when LL improves (big impact on results)
            updateMu(mu_new)
            updateScale(scale_new)
            if(verbose>=2) print(sprintf('Adaptive parameters converged after %d iterations', n_count))
            break
        }
        n_count = n_count + 1
        mu = mu_new
        scale = scale_new
    }

    LL = integrate_PLN_RE_LL(y, xb, group, sigma, gamma, panel.count.env$mu, panel.count.env$scale, H)

    if(verbose>=1){
        cat(sprintf('==== Likelihood function call %d: LL=%.5f ====\n', panel.count.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    if(!is.na(LL) && LL > panel.count.env$LL) {
        # update LL only when improves
        if(panel.count.env$stopUpdate==FALSE && (LL - panel.count.env$LL < 1e-6 * abs(panel.count.env$LL)) ) {
            updateStopUpdate(TRUE)
            if(verbose>=2) print('~~ Stop updating mu and scale now')
        }
        updateLL(LL)
    }
    return (LL)
}


# 2. Gradient function
Gradient_PLN_RE_AGQ = function(par,y,x,group,H,variance=FALSE,verbose=1){
    par = transformToBounded(par, c(sigma='exp', gamma='exp'))
    beta = par[1:ncol(x)]
    sigma = par['sigma']
    gamma = par['gamma']

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes
    mu = panel.count.env$mu
    scale = panel.count.env$scale

    xb = as.vector(x %*% beta)
    x_ext = cbind(x, gamma=0, sigma=0)

    Li = rep(0, N)
    dLogLi = matrix(0, N, ncol(x_ext))
    for(h in 1:H){
        tau = sqrt(2)*v[h]*scale + mu
        tau_rep = rep(tau, times=diff(group)) # i,t level

        # product of probabilities by individual
        P_int = rep(0, obs) # P_int, integral over e_it (same as sumK in LL function)
        dP_int = matrix(0, obs, ncol(x_ext)) # dP_int/dtheta
        for(k in 1:H){
            lam = exp(xb + sigma*tau_rep + sqrt(2)*gamma*v[k])
            wP = Omega[k] / sqrt(pi) * dpois(y, lam)
            P_int = P_int + wP
            # switch gamma and sigma positions to accommodate matVecProdSumExt
            dP_int = dP_int + matVecProdSumExt(x_ext, sqrt(2)*v[k]*gamma, tau_rep*sigma, wP*(y-lam), numeric(0))
        }
        prodT = as.vector(groupProd(P_int, group))
        Lih = sqrt(2)*Omega[h]*exp(v[h]^2) * scale * dnorm(tau) * prodT
        Li = Li + Lih

        # if P_int=0, the corresponding rows in dP_int are likely to be zero as well
        sumT = matVecProdSum(dP_int, numeric(0), 1/pmax(P_int,1e-100), group)
        dLogLi = dLogLi + matVecProd(sumT, Lih)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    dLogLi = matVecProd(dLogLi, 1/Li)
    # switch gamma and sigma back
    dLogLi = dLogLi[, c(1:(length(par)-2), length(par), length(par)-1)]
    colnames(dLogLi) = names(par)
    # transformation
    dLogLi = transformDerivative(dLogLi, tail(par, 2), c(sigma='exp', gamma='exp'))
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


Gradient_PLN_RE = function(par,y,x,group,H,variance=FALSE,verbose=1){
    par = transformToBounded(par, c(sigma='exp', gamma='exp'))
    beta = par[1:ncol(x)]
    sigma = par['sigma']
    gamma = par['gamma']

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    xb = as.vector(x %*% beta)
    x_ext = cbind(x,sigma=0,gamma=0)

    Li = rep(0, N)
    sumH = matrix(0, N, ncol(x_ext))
    for(h in 1:H){
        sumK = rep(0, obs)
        sumK_ext = matrix(0, obs, ncol(x_ext))
        for(k in 1:H){
            lamh = exp(xb + sqrt(2)*sigma*v[h] + sqrt(2)*gamma*v[k])
            PhiP = Omega[k] / sqrt(pi) * dpois(y, lamh)
            sumK = sumK + PhiP
            sumK_ext = sumK_ext + matVecProdSum(x_ext, c(sqrt(2)*v[h]*sigma,sqrt(2)*v[k]*gamma), PhiP*(y-lamh), numeric(0))
        }
        prodT = groupProd(sumK, group)
        Li = Li + Omega[h] * prodT
        # if sumK=0, the corresponding rows in sumK_ext are likely to be zero as well
        sumT = matVecProdSum(sumK_ext, numeric(0), 1/pmax(sumK,1e-100), group)
        sumH = sumH + matVecProd(sumT, Omega[h] * prodT)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    dLogLi = matVecProd(sumH, 1/Li)
    colnames(dLogLi) = names(par)
    # transformation
    dLogLi = transformDerivative(dLogLi, tail(par, 2), c(sigma='exp', gamma='exp'))
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

#' A Poisson Lognormal Model with Random Effects
#' @description Estimate a Poisson model with random effects at the individual and individual-time levels.
#' \deqn{E[y_{it}|x_{it},v_i,\epsilon_{it}] = exp(\boldsymbol{\beta}\mathbf{x_{it}}' + \sigma v_i + \gamma \epsilon_{it})}{E[y_it | x_it,v_i,\epsilon_it] = exp(\beta*x_it' + \sigma*v_i + \gamma*\epsilon_it)}
#' Notations:
#' * \eqn{x_{it}}{x_it}: variables influencing the selection decision \eqn{y_{it}}{y_it}, which could be a mixture of time-variant variables, time-invariant variables, and time dummies
#' * \eqn{v_i}: individual level random effect
#' * \eqn{\epsilon_{it}}{\epsilon_it}: individual-time level random effect
#'
#' \eqn{v_i} and \eqn{\epsilon_{it}}{\epsilon_it} can both account for overdispersion.
#' @md
#' @param formula Formula of the model
#' @param data Input data, a data.frame object
#' @param id.name The name of the column representing id. Data will be sorted by id to improve estimation speed.
#' @param method Optimization method used by optim. Defaults to 'BFGS'.
#' @param se_type Report Hessian or BHHH standard errors. Defaults to BHHH.
#' @param par Starting values for estimates. Default to estimates of Poisson RE model.
#' @param adaptiveLL Whether to use Adaptive Gaussian Quadrature. Defaults to TRUE because it is more reliable (though slower) for long panels.
#' @param stopUpdate Whether to disable update of Adaptive Gaussian Quadrature parameters. Defaults to FALSE.
#' @param sigma Starting value for sigma. Defaults to 1 and will be ignored if par is provided.
#' @param gamma Starting value for gamma. Defaults to 1 and will be ignored if par is provided.
#' @param H Number of Quadrature points used for numerical integration using the Gaussian-Hermite Quadrature method. Defaults to 20.
#' @param psnH Number of Quadrature points for Poisson RE model
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
#' @examples
#' \donttest{
#' # Use the simulated dataset, in which the true coefficient of x is 1.
#' # Estimated coefficient is biased due to omission of self-selection
#' data(sim)
#' res = PLN_RE(y~x, data=sim[!is.na(sim$y), ], id.name='id', verbose=-1)
#' res$estimates
#' }
#' @export
#' @family PanelCount
#' @references
#' 1. Peng, J., & Van den Bulte, C. (2023). Participation vs. Effectiveness in Sponsored Tweet Campaigns: A Quality-Quantity Conundrum. Management Science (forthcoming). Available at SSRN: <https://www.ssrn.com/abstract=2702053>
#'
#' 2. Peng, J., & Van den Bulte, C. (2015). How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. 2015 International Conference on Information Systems. <https://aisel.aisnet.org/icis2015/proceedings/SocialMedia/24/>
PLN_RE = function(formula, data, id.name, par=NULL, sigma=NULL, gamma=NULL, method='BFGS', adaptiveLL=TRUE, stopUpdate=FALSE, se_type=c('BHHH', 'Hessian')[1], H=12, psnH=12, reltol=sqrt(.Machine$double.eps), verbose=0){
    # 1.1 Sort data based on id
    data = data[order(data[, id.name]), ]
    group = c(0,cumsum(table(as.integer(factor(data[, id.name])))))

    # 1.2 Parse x and y
    mf = model.frame(formula, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y = model.response(mf, "numeric")
    x = model.matrix(attr(mf, "terms"), data=mf)

    # 1.3 Initial values
    if(is.null(par)){
        psn_re = PoissonRE(formula, id.name=id.name, data=data, sigma=sigma, method=method, H=psnH, reltol=reltol, verbose=verbose-1)
        par = c(psn_re$estimates[,1], gamma=1)
        if(!is.null(sigma)) par['sigma'] = sigma
        if(!is.null(gamma)) par['gamma'] = gamma
    }
    # check parameters
    if(length(par) != ncol(x)+2)  stop(sprintf('Number of parameters incorrect: desired %d, provided %d', ncol(x)+2, length(par)))
    if(names(par)[length(par)-1] != 'sigma') stop('Second last parameter must be named sigma')
    if(names(par)[length(par)] != 'gamma') stop('Last parameter must be named gamma')
    # transform to unbounded parameters
    if(method != 'L-BFGS-B')
        par = transformToUnbounded(par, c(sigma='exp', gamma='exp'))

    # 2. Estimation
    resetIter()
    updateMethod(method)
    updateMu(rep(0, length(group)-1))
    updateScale(rep(1, length(group)-1))
    updateStopUpdate(stopUpdate)
    if(exists('psn_re')){
        # psn re mu and scale are very close to pln re
        updateMu(psn_re$mu)
        updateScale(psn_re$scale)
    }

    gr = ifelse(adaptiveLL==TRUE, Gradient_PLN_RE_AGQ, Gradient_PLN_RE)
    res = optim(par=par, fn=ifelse(adaptiveLL==TRUE, LL_PLN_RE_AGQ, LL_PLN_RE), gr=gr, method=method, hessian=(se_type=='Hessian'), control=list(reltol=reltol,fnscale=-1), y=y, x=x, group=group, H=H, verbose=verbose)

    # 3. Likelihood, standard error, and p values
    res$n_obs = length(y)
    gvar = gr(res$par,y,x,group,H,variance=TRUE,verbose=verbose-1)
    res = compileResults(res, gvar, se_type, trans_vars=c(sigma='sigma', gamma='gamma'), trans_types=c('exp', 'exp'))

    res$mu = panel.count.env$mu
    res$scale = panel.count.env$scale
    res$psn_re = psn_re

    if(verbose>=0){
        cat(sprintf('==== PLN RE Model converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$counts[1], res$LL, res$gtHg))
        print(res$estimates, digits=3)
        print(res$time <- Sys.time() - panel.count.env$begin)
    }
    return (res)
}
