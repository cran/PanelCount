# mu and scale are defined as the user level
# xb and y are defined at the user-time level
integrate_LL = function(y, xb, group, sigma, mu, scale, H){
    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    N = length(group)-1
    obs = length(y)

    Li = rep(0, N)
    for(h in 1:H){
        tau = sqrt(2) * scale * v[h] + mu
        lam = exp(xb + sigma * rep(tau, times=diff(group)))
        ix = !is.na(y)
        Pz = rep(1, obs)
        Pz[ix] = dpois(y[ix], lam[ix])
        Pz_prod = as.vector(groupProd(Pz, group))

        Li = Li + Omega[h] * exp(v[h]^2) * dnorm(tau) * Pz_prod
    }
    Li = pmax(sqrt(2)*scale*Li, 1e-100)
    sum(log(Li))
}


# Adaptive Gaussian Quadrature
LL_PoissonRE_AGQ = function(par, y, x, group, H, verbose=1){
    par = transformToBounded(par, c(sigma='exp'))
    beta = par[1:ncol(x)]
    sigma = par['sigma']

    N = length(group)-1
    obs = length(y)

    xb = as.vector(x %*% beta)

    # update mode and scale
    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes
    mu = panel.count.env$mu
    scale = panel.count.env$scale

    n_count = 1
    # Usually converge in a few iterations
    while(panel.count.env$stopUpdate==FALSE){
        Li = rep(0, N)
        mu_new = rep(0, N)
        scale_new = rep(0, N)
        for(h in 1:H){
            tau = sqrt(2) * scale * v[h] + mu
            lam = exp(xb + sigma * rep(tau, times=diff(group)))
            ix = !is.na(y)
            Pz = rep(1, obs)
            Pz[ix] = dpois(y[ix], lam[ix])
            Pz_prod = as.vector(groupProd(Pz, group))

            # scale is a vector, calculate products of scale first
            tmp = sqrt(2)*Omega[h]*exp(v[h]^2) * scale * dnorm(tau) * Pz_prod
            Li = Li + tmp
            mu_new = mu_new + tau*tmp
            scale_new = scale_new + tau^2*tmp
        }
        Li = pmax(Li, 1e-100)
        mu_new = mu_new / Li
        scale_new = sqrt(pmax(scale_new / Li - mu_new^2, 1e-6))

        # Update only when LL improves
        LL = sum(log(Li))
        if(is.na(LL) || LL < panel.count.env$LL) break

        # Stop on error
        if(any(is.na(scale_new)) || any(is.na(mu_new)) || any(scale_new<0) || n_count >= 100){
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
        if(all(abs(mu_new - mu) < pmax(1e-4*abs(mu), 1e-6))
           && all(abs(scale_new - scale) < pmax(1e-4*abs(scale), 1e-6))){
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

    LL = integrate_LL(y, xb, group, sigma, panel.count.env$mu, panel.count.env$scale, H)

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

Gradient_PoissonRE_AGQ = function(par,y,x,group,H,variance=FALSE,verbose=1){
    par = transformToBounded(par, c(sigma='exp'))
    beta = par[1:ncol(x)]
    sigma = par['sigma']

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes
    mu = panel.count.env$mu
    scale = panel.count.env$scale

    xb = as.vector(x %*% beta)
    x_ext = cbind(x, sigma=0)

    Li = rep(0, N)
    dLogLi = matrix(0, N, ncol(x_ext))
    for(h in 1:H){
        tau = sqrt(2) * scale * v[h] + mu
        tau_rep = rep(tau, times=diff(group)) # i,t level
        lam = exp(xb + sigma * tau_rep)
        ix = !is.na(y)
        Pz = rep(1, obs)
        Pz[ix] = dpois(y[ix], lam[ix])
        Pz_prod = as.vector(groupProd(Pz, group))

        # The scalor sqrt(2) * scale is ignored as it will be canceled out
        hx = Omega[h] * exp(v[h]^2) * dnorm(tau) * Pz_prod
        Li = Li + hx

        # Gradient
        yl = rep(0, obs)
        yl[ix] = y[ix] - lam[ix]
        # matVecProdSumExt: replace last column by sigma * tau_rep and multiple by y-lam, then sum by group. The second parameter is set to empty
        sumT = matVecProdSumExt(x_ext, numeric(0), sigma * tau_rep, yl, group)
        dLogLi = dLogLi + matVecProd(sumT, hx)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    dLogLi = matVecProd(dLogLi, 1/Li)
    colnames(dLogLi) = names(par)
    # transformation
    dLogLi = transformDerivative(dLogLi, tail(par, 1), c(sigma='exp'))
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


#' A Poisson Model with Random Effects
#' @description Estimate a Poisson model with random effects at the individual level.
#' \deqn{E[y_{it}|x_{it},v_i] = exp(\boldsymbol{\beta}\mathbf{x_{it}}' + \sigma v_i)}{E[y_it | x_it,v_i] = exp(\beta*x_it' + \sigma*v_i)}
#' Notations:
#' * \eqn{x_{it}}{x_it}: variables influencing the outcome \eqn{y_{it}}{y_it}, which could be a mixture of time-variant variables, time-invariant variables, and time dummies
#' * \eqn{v_i}: individual level random effect
#' @md
#' @param formula Formula of the model
#' @param data Input data, a data.frame object
#' @param id.name The name of the column representing id. Data will be sorted by id to improve estimation speed.
#' @param method Optimization method used by optim. Defaults to 'BFGS'.
#' @param se_type Report Hessian or BHHH standard errors. Defaults to Hessian.
#' @param par Starting values for estimates. Default to estimates of Poisson Model
#' @param stopUpdate Whether to disable update of Adaptive Gaussian Quadrature parameters. Defaults to FALSE.
#' @param sigma Starting value for sigma. Defaults to 1 and will be ignored if par is provided.
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
#' @examples
#' # Use the simulated dataset, in which the true coefficient of x is 1.
#' # Estimated coefficient is biased primarily due to omission of self-selection
#' data(sim)
#' res = PoissonRE(y~x, data=sim[!is.na(sim$y), ], id.name='id', verbose=-1)
#' res$estimates
#' @export
#' @family PanelCount
#' @references
#' 1. Peng, J., & Van den Bulte, C. (2023). Participation vs. Effectiveness in Sponsored Tweet Campaigns: A Quality-Quantity Conundrum. Management Science (forthcoming). Available at SSRN: <https://www.ssrn.com/abstract=2702053>
#'
#' 2. Peng, J., & Van den Bulte, C. (2015). How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. 2015 International Conference on Information Systems. <https://aisel.aisnet.org/icis2015/proceedings/SocialMedia/24/>
PoissonRE = function(formula, data, id.name, par=NULL, sigma=NULL, method='BFGS', stopUpdate=FALSE, se_type=c('Hessian', 'BHHH')[1], H=20, reltol=sqrt(.Machine$double.eps), verbose=0){
    # 1.1 Sort data based on id
    data = data[order(data[, id.name]), ]
    group = c(0,cumsum(table(as.integer(factor(data[, id.name])))))

    # 1.2 Parse x and y
    mf = model.frame(formula, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y = model.response(mf, "numeric")
    x = model.matrix(attr(mf, "terms"), data=mf)

    # 1.3 use Poisson estimates as initial values for parameters
    if(is.null(par)){
        psn = summary(glm(formula, family=poisson(link="log"), data=data))
        par = c(psn$coefficients[,1], sigma=1)
        if(!is.null(sigma)) par['sigma'] = sigma
    }
    if(length(par) != ncol(x)+1)  {
        stop(sprintf('Number of parameters incorrect: desired %d, provided %d', ncol(x)+1, length(par)))
    }
    if(names(par)[length(par)] != 'sigma') stop('Last parameter must be named sigma')
    # transform to unbounded parameters
    if(method != 'L-BFGS-B')
        par = transformToUnbounded(par, c(delta='sigma'))

    # 2. Estimation
    resetIter()
    updateMethod(method)
    updateMu(rep(0, length(group)-1))
    updateScale(rep(1, length(group)-1))
    updateStopUpdate(stopUpdate)
    res = optim(par=par, fn=LL_PoissonRE_AGQ, gr=Gradient_PoissonRE_AGQ, method=method, hessian=(se_type=='Hessian'), control=list(reltol=reltol,fnscale=-1), y=y, x=x, group=group, H=H, verbose=verbose)

    # 3. Likelihood, standard error, and p values
    res$n_obs = length(y)
    gvar = Gradient_PoissonRE_AGQ(res$par,y,x,group,H,variance=TRUE,verbose=verbose-1)
    res = compileResults(res, gvar, se_type, trans_vars=c(sigma='sigma'), trans_types=c('exp'))

    res$mu = panel.count.env$mu
    res$scale = panel.count.env$scale

    if(verbose>=0){
        cat(sprintf('==== Poisson RE Model converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$counts[1], res$LL, res$gtHg))
        print(res$estimates, digits=3)
        print(res$time <- Sys.time() - panel.count.env$begin)
    }
    return (res)
}
