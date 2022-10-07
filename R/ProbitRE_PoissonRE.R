LL_CRE = function(par,y,z,x,w,group,rule,offset_w=NULL,offset_x=NULL,verbose=1){
    par = transformToBounded(par, c(delta='exp', sigma='exp', rho='correlation'))
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    delta = par['delta']
    sigma = par['sigma']
    rho = par['rho']

    # Omega & v (inner), Weight & r (outer)
    Omega = rule$Omega
    v = rule$v
    Weight = rule$Weight
    r = rule$r

    N = length(group)-1
    obs = length(y)

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    if(!is.null(offset_w)) wa = wa + offset_w
    if(!is.null(offset_x)) xb = xb + offset_x
    d = 2*z - 1

    Li = rep(0, N)
    for(h in 1:length(r)){
        lamh = exp(xb + sqrt(2)*sigma*r[h])
        ix = (z==1)
        Pz = rep(1, obs)
        Pz[ix] = dpois(y[ix], lamh[ix])
        Pz_prod = as.vector(groupProd(Pz, group))
        sumk = rep(0, N)
        for(k in 1:length(v)){
            q = wa + sqrt(2)*delta*rho*r[h] + sqrt(2*(1-rho^2))*delta*v[k]
            Phi = pnorm(d*q)
            Phi_prod = as.vector(groupProd(Phi, group)) # @@@
            sumk = sumk + Omega[k]*Phi_prod
        }
        Li = Li + Weight[h] * Pz_prod * sumk
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    LL = sum(log(Li/pi))

    if(verbose>=1){
        cat(sprintf('==== Likelihood function call %d: LL=%.5f ====\n', panel.count.env$iter, LL))
        print(par,digits=3)
        print(Sys.time() - panel.count.env$begin)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}


# 2. Gradient function
Gradient_CRE = function(par,y,z,x,w,group,rule,offset_w=NULL,offset_x=NULL,variance=FALSE,verbose=1){
    par = transformToBounded(par, c(delta='exp', sigma='exp', rho='correlation'))
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    delta = par['delta']
    sigma = par['sigma']
    rho = par['rho']

    N = length(group)-1
    obs = length(y)

    # Omega & v (inner), Weight & r (outer)
    Omega = rule$Omega
    v = rule$v
    Weight = rule$Weight
    r = rule$r

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    if(!is.null(offset_w)) wa = wa + offset_w
    if(!is.null(offset_x)) xb = xb + offset_x
    d = 2*z - 1

    w_ext = cbind(w,delta=0,rho=0) # preallocate memory for delta and rho
    x_ext = cbind(x,sigma=0)
    # renames variables for uniqueness
    colnames(w_ext) = c(paste0("w",1:ncol(w)),"delta","rho")
    colnames(x_ext) = c(paste0("x",1:ncol(x)),"sigma")

    Li = rep(0, N)
    sumH_alp = matrix(0, N, ncol(w_ext))
    sumH_beta = matrix(0, N, ncol(x_ext))
    colnames(sumH_alp) = colnames(w_ext)
    colnames(sumH_beta) = colnames(x_ext)

    # saved coef
    a1 = sqrt(2)*rho
    a2 = sqrt(2*(1-rho^2))
    a3 = sqrt(2)*delta
    a4 = sqrt(2/(1-rho^2))*rho*delta

    b1 = sqrt(2)*delta*rho
    b2 = sqrt(2*(1-rho^2))*delta

    for(h in 1:length(r)){
        lamh = exp(xb + sqrt(2)*sigma*r[h])
        ix = (z==1)
        Pz = rep(1, obs)
        Pz[ix] = dpois(y[ix], lamh[ix])
        Pz_prod = as.vector(groupProd(Pz, group))

        zyl = rep(0, obs)
        zyl[ix] = y[ix] - lamh[ix]

        # meta data for beta and sigma
        sumK_beta = rep(0, N)
        sumT_beta = matVecProdSum(x_ext, sqrt(2)*r[h], zyl, group)

        # meta data for alpha, delta, and rho
        sumK_alp = matrix(0, N, ncol(w_ext))
        for(k in 1:length(v)){
            dq = d * (b1*r[h] + b2*v[k] + wa)
            phi = dnorm(dq)
            Phi = pnorm(dq)
            Phi_prod = as.vector(groupProd(Phi, group)) # @@@
            sumK_beta = sumK_beta + Omega[k]*Phi_prod

            # delta and rho
            sumT_alp = matVecProdSum(w_ext,c(a1*r[h] + a2*v[k], a3*r[h] - a4*v[k]), d*phi/Phi, group) # @@@
            sumK_alp = sumK_alp + matVecProd(sumT_alp, Omega[k]*Phi_prod)
        }
        Li = Li + Weight[h] * Pz_prod * sumK_beta
        sumH_alp = sumH_alp + matVecProd(sumK_alp, Weight[h] * Pz_prod)
        sumH_beta = sumH_beta + matVecProd(sumT_beta, Weight[h]*Pz_prod*sumK_beta)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    # whether rho is included depends on length of par, pi cancells out
    dLogLi = matVecProd(cbind(sumH_alp, sumH_beta), 1/Li)
    colnames(dLogLi) = c(colnames(sumH_alp), colnames(sumH_beta))
    dLogLi = dLogLi[, c(paste0("w",1:ncol(w)),paste0("x",1:ncol(x)),names(par[(ncol(x)+ncol(w)+1):length(par)]))]
    colnames(dLogLi) = names(par)

    # transformation
    dLogLi = transformDerivative(dLogLi, tail(par, length(par)-ncol(x)-ncol(w)), c(delta='exp', sigma='exp', rho='correlation'))

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

# 3. Partial Effects on the fit sample
Partial_CRE = function(res,w,xnames,offset_w=NULL,intercept=FALSE){
    wnames = colnames(w)
    par = res$estimates[, 1] # delta, sigma, etc. already transformed

    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:length(xnames)]
    delta = par['delta']
    sigma = par['sigma']
    rho = par["rho"]

    com = xnames[xnames %in% wnames & xnames!="(Intercept)"]
    if(intercept){ # no effect on other variables
        xunq = c("(Intercept)", xnames[!xnames %in% wnames])
        wunq = c("(Intercept)", wnames[!wnames %in% xnames])
    } else {
        xunq = xnames[!xnames %in% wnames]
        wunq = wnames[!wnames %in% xnames]
    }

    wa = as.vector(w %*% alpha)
    if(!is.null(offset_w)) wa = wa + offset_w

    # Relative Partial Effects (G)
    s = (wa + delta*rho*sigma) / sqrt(1 + delta^2)
    ratio = dnorm(s) / (sqrt(1 + delta^2) * pnorm(s)) # c_it
    Gw = mean(ratio) * alpha
    Gx = beta
    G = Gw[com]+Gx[com]
    if(length(wunq)>0) G = c(G, Gw[wunq])
    if(length(xunq)>0) G = c(G, Gx[xunq])
    names(G) = c(com, wunq, xunq)

    # Average Gradient of Partial Effects (J)
    # The latter part of Jw_1 does not vary across obs and hence mean of Jw_1 is the same as mean of ratio
    Jw_1 = mean(ratio) * cbind(diag(length(alpha)), -delta/(1+delta^2)*alpha, matrix(0,length(alpha),2))
    # The first part of Jw_2 is the constant vector alpha and hence mean of Jw_2 is the same as mean of the latter part.
    w_ext = cbind(w,delta=rho*sigma-s*delta/sqrt(1+delta^2),sigma=delta*rho,rho=delta*sigma)
    Jw_2 = outer(alpha, colMeans(matVecProd(w_ext, -ratio * (s/sqrt(1+delta^2) + ratio))))
    Jw = Jw_1 + Jw_2
    # Now insert zeros for beta
    Jw = cbind(Jw[, 1:(ncol(Jw)-3), drop=FALSE], matrix(0,length(alpha),length(beta)), Jw[, tail(1:ncol(Jw), 3)])
    rownames(Jw) = wnames
    Jx = cbind(matrix(0,length(beta),length(alpha)), diag(length(beta)), matrix(0,length(beta),3))
    rownames(Jx) = xnames

    J = Jx[com, , drop=FALSE] + Jw[com, , drop=FALSE]
    if(length(wunq)>0) J = rbind(J, Jw[wunq, , drop=FALSE])
    if(length(xunq)>0) J = rbind(J, Jx[xunq, , drop=FALSE])
    rownames(J) = c(com, wunq, xunq)
    colnames(J) = c(wnames,xnames,"delta","sigma","rho")

    # rearrange to align with covariance matrix, remove missing structural parameters
    plain_par = length(wnames) + length(xnames)
    J = cbind(J[,1:plain_par,drop=FALSE], J[, tail(names(par), length(par)-plain_par),drop=FALSE])
    se = sqrt(diag(J %*% res$var_bhhh %*% t(J)))
    z = G/se
    p = 1 - pchisq(z^2, 1)
    pe = cbind(estimates=G,se=se,z=z,p=p)
    return (pe)
}


#' Predictions of CRE model on new sample
#' @description Predictions of CRE model on new sample. Please make sure the factor variables in the test data do not have levels not shown in the training data.
#' @param par Model estimates
#' @param sel_form Formula for selection equation, a Probit model with random effects
#' @param out_form Formula for outcome equation, a Poisson Lognormal model with random effects
#' @param data Input data, a data.frame object
#' @param offset_w_name Offset variables in selection equation, if any.
#' @param offset_x_name Offset variables in outcome equation, if any.
#' @return A list with three sets of predictions
#' * prob: Predicted probability to participate
#' * outcome: Predicted potential outcome
#' * actual_outcome: Predicted actual outcome
#' @md
#' @export
predict_ProbitRE_PoissonRE = function(par,sel_form,out_form,data,offset_w_name=NULL,offset_x_name=NULL){
    # keep unused levels in case test data have less levels, but can lead to error if test data have more levels than training data
    mf = model.frame(sel_form, data=data, na.action=NULL, drop.unused.levels=FALSE)
    w = model.matrix(attr(mf, "terms"), data=mf)
    mf2 = model.frame(out_form, data=data, na.action=NULL, drop.unused.levels=FALSE)
    x = model.matrix(attr(mf2, "terms"), data=mf2)

    if(length(par) != ncol(x)+ncol(w)+3){
        cat('Variable names in par\n')
        print(names(par))
        cat('Variable names in w and x\n')
        print(c(colnames(w), colnames(x)))
        stop(sprintf('Number of parameters incorrect: desired %d, provided %d. Check if the test data have more or less levels than the training data in any factor variables', ncol(x)+ncol(w)+3, length(par)))
    }

    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    delta = par['delta']
    sigma = par['sigma']
    rho = par['rho']

    if(!identical(names(alpha), colnames(w)) || !identical(names(beta), colnames(x))){
        cat('Variable names in par\n')
        print(names(par))
        cat('Variable names in w and x\n')
        print(c(colnames(w), colnames(x)))
        stop('Parameter variable names with supplied data mismatch. Check if the test data have more or less levels than the training data in any factor variables')
    }

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    if(!is.null(offset_w_name)) wa = wa + data[, offset_w_name]
    if(!is.null(offset_x_name)) xb = xb + data[, offset_x_name]

    prob = pnorm(wa/sqrt(1+delta^2))
    outcome = exp(xb+sigma^2/2)
    actual_outcome = outcome * pnorm((wa+rho*sigma*delta)/sqrt(1+delta^2))
    return (list(prob=prob, outcome=outcome, actual_outcome=actual_outcome))
}

#' Poisson RE model with Sample Selection
#' @description Estimates the following two-stage model \cr \cr
#' Selection equation (ProbitRE - Probit model with individual level random effects):
#' \deqn{z_{it}=1(\boldsymbol{\alpha}\mathbf{w_{it}}'+\delta u_i+\xi_{it} > 0)}{z_it=1(\alpha*w_it'+\delta*u_i+\xi_it > 0)}
#' Outcome Equation (PoissonRE - Poisson with individual level random effects):
#' \deqn{E[y_{it}|x_{it},v_i] = exp(\boldsymbol{\beta}\mathbf{x_{it}}' + \sigma v_i)}{E[y_it | x_it,v_i] = exp(\beta*x_it' + \sigma*v_i)}
#' Correlation (self-selection at individual level):
#' * \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}.
#'
#' Notations:
#' * \eqn{w_{it}}{w_it}: variables influencing the selection decision \eqn{z_{it}}{z_it}, which could be a mixture of time-variant variables, time-invariant variables, and time dummies
#' * \eqn{x_{it}}{x_it}: variables influencing the outcome \eqn{y_{it}}{y_it}, which could be a mixture of time-variant variables, time-invariant variables, and time dummies
#' * \eqn{u_i}: individual level random effect in the selection equation
#' * \eqn{v_i}: individual level random effect in the outcome equation
#' * \eqn{\xi_{it}}{\xi_it}: error term in the selection equation
#' @param sel_form Formula for selection equation, a Probit model with random effects
#' @param out_form Formula for outcome equation, a Poisson model with random effects
#' @param data Input data, a data.frame object
#' @param testData Test data for prediction, a data.frame object
#' @param id.name The name of the column representing id. Data will be sorted by id to improve estimation speed.
#' @param method Optimization method used by optim. Defaults to 'BFGS'.
#' @param se_type Report Hessian or BHHH standard errors. Defaults to BHHH.
#' @param par Starting values for estimates. Default to estimates of standalone selection and outcome models.
#' @param delta Starting value for delta. Will be ignored if par is provided.
#' @param sigma Starting value for sigma. Will be ignored if par is provided.
#' @param rho Starting value for rho. Defaults to 0 and will be ignored if par is provided.
#' @param H A integer vector of length 2, specifying the number of points for inner and outer Quadratures
#' @param psnH Number of Quadrature points for Poisson RE model
#' @param prbH Number of Quddrature points for Probit RE model
#' @param reltol Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
#' @param offset_w_name An offset variable whose coefficient is assumed to be 1 in the selection equation
#' @param offset_x_name An offset variable whose coefficient is assumed to be 1 in the outcome equation
#' @param verbose A integer indicating how much output to display during the estimation process.
#' * <0 - No ouput
#' * 0 - Basic output (model estimates)
#' * 1 - Moderate output, basic ouput + parameter and likelihood in each iteration
#' * 2 - Extensive output, moderate output + gradient values on each call
#' @return A list containing the results of the estimated model, some of which are inherited from the return of optim
#' * estimates: Model estimates with 95% confidence intervals
#' * par: Point estimates
#' * var_bhhh: BHHH covariance matrix, inverse of the outer product of gradient at the maximum
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
#' @md
#' @examples
#' \donttest{
#' # Use the simulated dataset, in which the true coefficients of x and w are 1 in both stages.
#' # The simulated dataset includes self-selection at both individual and individual-time level,
#' # but this model only considers self-selection at the individual level.
#' data(sim)
#' res = ProbitRE_PoissonRE(z~x+w, y~x, data=sim, id.name='id')
#' res$estimates
#' }
#' @export
#' @family PanelCount
#' @references
#' 1. Peng, J., & Van den Bulte, C. (2022). Participation vs. Effectiveness in Sponsored Tweet Campaigns: A Quality-Quantity Conundrum. Available at SSRN: https://ssrn.com/abstract=2702053
#'
#' 2. Peng, J., & Van Den Bulte, C. (2015). How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. 2015 International Conference on Information Systems. https://aisel.aisnet.org/icis2015/proceedings/SocialMedia/24
ProbitRE_PoissonRE = function(sel_form, out_form, data, id.name, testData=NULL, par=NULL, delta=NULL, sigma=NULL, rho=NULL, method='BFGS', se_type=c('BHHH', 'Hessian')[1], H=c(10,10), psnH=20, prbH=20, reltol=sqrt(.Machine$double.eps), verbose=1, offset_w_name=NULL, offset_x_name=NULL){
    # 1.1 Sort data based on id
    data = data[order(data[, id.name]), ]
    group = c(0,cumsum(table(as.integer(factor(data[, id.name])))))
    offset_w = offset_x = NULL
    if(!is.null(offset_w_name)) offset_w = data[, offset_w_name]
    if(!is.null(offset_x_name)) offset_x = data[, offset_x_name]

    # 1.2 Quadrature rules
    rule1 = gauss.quad(H[1], "hermite")
    rule2 = gauss.quad(H[2], "hermite")
    rule = list(Weight=rule1$weights, r=rule1$nodes, Omega=rule2$weights, v=rule2$nodes)

    # 1.3 parse w,z: need to make sure if data is reordered
    mf = model.frame(sel_form, data=data, na.action=NULL, drop.unused.levels=TRUE)
    z = model.response(mf, "numeric")
    w = model.matrix(attr(mf, "terms"), data=mf)

    # 1.4 parse x,y
    mf2 = model.frame(out_form, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y = model.response(mf2, "numeric")
    x = model.matrix(attr(mf2, "terms"), data=mf2)

    # 1.5 Initialize parameters
    if(is.null(par)){
        if(verbose>=1) cat('---Estimating outcome equation alone...\n')
        psn_re = PoissonRE(out_form, data=data[!is.na(y),], id.name=id.name, sigma=sigma, H=psnH, method='BFGS', se_type='BHHH', reltol=reltol, verbose=verbose-1)
        if(verbose>=1) cat('---Estimating selection equation alone...\n')
        probit_re = ProbitRE(sel_form, data=data, id.name=id.name, delta=delta, method='BFGS', se_type='BHHH', H=prbH, reltol=reltol, verbose=verbose-1)
        sel_est = probit_re$estimates[,1]
        out_est = psn_re$estimates[,1]
        par = c(sel_est[-length(sel_est)], out_est[-length(out_est)], sel_est[length(sel_est)], out_est[length(out_est)], rho=0)

        # update structural parameters to given initial values if provided
        if(!is.null(delta)) par['delta'] = delta
        if(!is.null(sigma)) par['sigma'] = sigma
        if(!is.null(rho)) par['rho'] = rho
    }
    # check parameter length and names
    if(length(par) != ncol(x)+ncol(w)+3){
        stop(sprintf('Number of parameters incorrect: desired %d, provided %d', ncol(x)+ncol(w)+3, length(par)))
    }
    if(names(par)[length(par)-2] != 'delta') stop('Third last parameter must be named delta')
    if(names(par)[length(par)-1] != 'sigma') stop('Second last parameter must be named sigma')
    if(names(par)[length(par)] != 'rho') stop('Last parameter must be named rho')
    # transform to unbounded parameters
    if(method != 'L-BFGS-B')
        par = transformToUnbounded(par, c(delta='exp', sigma='exp', rho='correlation'))

    # 2. Estimation
    resetIter()
    updateMethod(method)
    res = optim(par=par, fn=LL_CRE, gr=Gradient_CRE, method=method, hessian=(se_type=='Hessian'), control=list(reltol=reltol,fnscale=-1), y=y, z=z, x=x, w=w, group=group, rule=rule, offset_w=offset_w, offset_x=offset_x, verbose=verbose)

    # 3. Likelihood, standard error, and p values
    res$n_obs = length(z)
    gvar = Gradient_CRE(res$par,y,z,x,w,group,rule,offset_w=offset_w,offset_x=offset_x,variance=TRUE,verbose=verbose-1)
    res = summary.panel.count(res, gvar, se_type=se_type, trans_vars=c(sigma='sigma', delta='delta', rho='rho'), trans_types=c('exp', 'exp', 'correlation'))

    # 4.1 making predictions on new data
    if(!is.null(testData)) data = testData
    res$predict = predict_ProbitRE_PoissonRE(res$estimates[, 1],sel_form,out_form,data,offset_w_name=offset_w_name,offset_x_name=offset_x_name)
    par_no_cor = c(sel_est[-length(sel_est)], out_est, sel_est[length(sel_est)], rho=0)
    res$predict_no_cor = predict_ProbitRE_PoissonRE(par_no_cor,sel_form,out_form,data,offset_w_name=offset_w_name,offset_x_name=offset_x_name)

    # 4.2 Partial Effects on the fit sample
    res$partial = Partial_CRE(res,w,colnames(x),offset_w=offset_w)
    wavg = t(colMeans(w)) # convert vector to matrix with a single row (names kept)
    avg_offset_w = NULL
    if(!is.null(offset_w)) avg_offset_w = mean(offset_w)
    res$partialAvgObs = Partial_CRE(res,wavg,colnames(x),offset_w=avg_offset_w)


    # 5. Meta data
    res$input = list(sel_form=sel_form, out_form=out_form, par=par, delta=delta, sigma=sigma, rho=rho, method=method, se_type=se_type, H=H,psnH=psnH,prbH=prbH,reltol=reltol,verbose=verbose,offset_w_name=offset_w_name,offset_x_name=offset_x_name)
    res$psn_re = psn_re
    res$probit_re = probit_re

    # LR test statistics can be negative due to numerical integration
    res$LR_stat = 2 * ( res$LL - probit_re$LL - psn_re$LL )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)

    if(verbose>=0){
        cat(sprintf('==== CRE Model converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$counts[1], res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$estimates, digits=3)
        print(res$time <- Sys.time() - panel.count.env$begin)
    }
    return (res)
}
