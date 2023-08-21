# a few utility functions
panel.count.env = new.env(parent = emptyenv())
panel.count.env$iter = 1

resetIter = function(){
    panel.count.env$begin = Sys.time()
    panel.count.env$iter = 1
    panel.count.env$LL = -.Machine$double.xmax
    panel.count.env$mu = NULL
    panel.count.env$scale = NULL
    panel.count.env$stopUpdate = NULL
    panel.count.env$method = NULL
}

addIter = function(){
    panel.count.env$iter = panel.count.env$iter + 1
}

updateMethod = function(method){
    panel.count.env$method = method
}

updateMu = function(mu){
    panel.count.env$mu = mu
}

updateScale = function(scale){
    panel.count.env$scale = scale
}

updateLL = function(LL){
    panel.count.env$LL = LL
}

updateStopUpdate = function(stopUpdate){
    panel.count.env$stopUpdate = stopUpdate
}

.onUnload <- function (libpath) {
    library.dynam.unload("PanelCount", libpath)
}

getVarSE = function(res, gvar, se_type=c('BHHH', 'Hessian')[1]){
    # Greene: For small and moderate sized samples, Hessian is preferable
    # https://en.wikipedia.org/wiki/Information_matrix_test
    # https://stats.stackexchange.com/questions/122149/name-for-outer-product-of-gradient-approximation-of-hessian
    # BHHH SE
    res$var_bhhh = gvar$var
    res$g = gvar$g
    rownames(res$var_bhhh) = rownames(res$hessian)
    colnames(res$var_bhhh) = colnames(res$hessian)
    res$se_bhhh = sqrt(diag(res$var_bhhh))

    # Hessian SE
    if(!is.null(res$hessian)){
        res$var_hessian = tryCatch( solve(-res$hessian), error = function(e){
            cat('Hessian not invertible: ', e$message, '\n')
            m = NA * diag(length(res$par))
            rownames(m) = colnames(m) = names(res$par)
            m
        } )
        res$se_hessian = sqrt(diag(res$var_hessian))
    }

    if(se_type=='BHHH') {
        res$se = res$se_bhhh
        res$gtHg = matrix(res$g, nrow=1) %*% res$var_bhhh %*% matrix(res$g, ncol=1)
    }
    if(se_type=='Hessian') {
        if(is.null(res$hessian)) stop("Error: cannot report Hessian matrix because Hessian matrix is not calculated")
        res$se = res$se_hessian
        res$gtHg = matrix(res$g, nrow=1) %*% res$var_hessian %*% matrix(res$g, ncol=1)
    }
    names(res$se) = names(res$par)
    res
}

# Transform to derivative w.r.t. original parameter: trans_vars stores <var_name, value> pairs
# Trans_types store <var_name, type> pairs
transformDerivative = function(dr, trans_vars, trans_types){
    trans_types = trans_types[names(trans_vars)] # in case some parameter is disabled
    if(length(trans_vars)!=length(trans_types))
        stop('Numbers of transformation variables and types are inconsistent are matching')
    if(panel.count.env$method == 'L-BFGS-B') return(dr) # no transformation needed

    ix = names(trans_vars) %in% colnames(dr)
    trans_vars = trans_vars[ix]
    trans_types = trans_types[ix]
    if(length(trans_vars)==0) return(dr)

    for(i in 1:length(trans_vars)){
        type = trans_types[i]
        var_name = names(trans_vars)[i]
        value = trans_vars[i]
        if(type == 'exp'){
            dr[, var_name] = value * dr[, var_name]
        }else if(type == 'correlation'){
            dr[, var_name] = (2 * exp(value) / (exp(value)+1)^2) * dr[, var_name]
        }else{
            stop(sprintf('Derivative transformation for type %s not implemented yet', type))
        }
    }
    dr
}

# Transform to unbounded parameter: trans_dict stores <var_name, type> pairs
transformToUnbounded = function(par, trans_dict){
    trans_dict = trans_dict[names(trans_dict) %in% names(par)]
    if(length(trans_dict)==0) return(par)
    for(i in 1:length(trans_dict)){
        type = trans_dict[i]
        var_name = names(trans_dict)[i]
        value = par[var_name]
        if(type == 'exp'){
            if(value <= 0){
                warning(sprintf('The parameter %s with a value of %.3f is non-positive and will be set to 1e-6', var_name, value))
                value = 1e-6
            }
            par[var_name] = log(value)
        }else if(type == 'correlation'){
            if(value >= 1){
                warning(sprintf('The parameter %s with a value of %.3f is greater than or equal to 1 and will be set to 1-1e-6', var_name, value))
                value = 1 - 1e-6
            }
            if(value <= -1){
                warning(sprintf('The parameter %s with a value of %.3f is less than or equal to -1 and will be set to -1+1e-6', var_name, value))
                value = -1 + 1e-6
            }
            par[var_name] = 0.5*log((1+value)/(1-value))
        }else if(type == 'logit'){
            if(value <= 0){
                warning(sprintf('The parameter %s with a value of %.3f is non-positive and will be set to 1e-6', var_name, value))
                value = 1e-6
            }
            if(value >= 1){
                warning(sprintf('The parameter %s with a value of %.3f is greater than or equal to 1 and will be set to 1-1e-6', var_name, value))
                value = 1 - 1e-6
            }
            par[var_name] = log(value/(1-value))
        }
    }
    par
}

# Transform to original bounded parameter: trans_dict stores <var_name, type> pairs
transformToBounded = function(par, trans_dict){
    if(panel.count.env$method == 'L-BFGS-B') return(par) # no transformation needed
    trans_dict = trans_dict[names(trans_dict) %in% names(par)]
    if(length(trans_dict)==0) return(par)
    for(i in 1:length(trans_dict)){
        type = trans_dict[i]
        var_name = names(trans_dict)[i]
        if(type == 'exp'){
            par[var_name] = exp(par[var_name])
        }else if(type == 'correlation'){
            par[var_name] = 1 - 2/(exp(par[var_name])+1)
        }else if(type == 'logit'){
            par[var_name] = plogis(par[var_name])
        }
    }
    par
}

# summarize results
compileResults = function(res, gvar, se_type=c('BHHH', 'Hessian')[1], trans_vars=NULL, trans_types=NULL){
    if(length(trans_vars)!=length(trans_types))
        stop('Numbers of transformation variables and types are inconsistent')

    # Miscellaneous
    res$LL = res$value
    res$AIC = -2*res$LL + 2 * length(res$par)
    res$BIC = -2*res$LL + log(res$n_obs) * length(res$par)
    res$n_par = length(res$par)
    res$LL_calls = panel.count.env$iter # res$counts store counts by optim

    # get SE
    res = getVarSE(res, gvar, se_type=se_type)

    # Retain the raw estimate and se
    par = res$par
    se = res$se

    lci = par - qnorm(0.975)*se
    uci = par + qnorm(0.975)*se

    if(panel.count.env$method != 'L-BFGS-B' && length(trans_vars)>=1){
        for(i in 1:length(trans_vars)){
            var_name = trans_vars[i]
            type = trans_types[i]

            # Note: delta method is not reliable when the original parameter is not normal
            # https://www.statisticshowto.com/delta-method-definition/#:~:text=Disadvantages,LePage%20%26%20Billard%2C%201992).
            if(type=='exp'){ # >0
                par[var_name] = exp(par[var_name])
                se[var_name] = par[var_name] * se[var_name]

                lci[var_name] = exp(lci[var_name])
                uci[var_name] = exp(uci[var_name])
            }else if(type=='correlation'){ # [-1, 1]
                tau = par[var_name]
                par[var_name] = 1 - 2/(exp(tau)+1)
                se[var_name] = 2*exp(tau)/(exp(tau)+1)^2 * se[var_name]

                lci[var_name] = 1 - 2/(exp(lci[var_name])+1)
                uci[var_name] = 1 - 2/(exp(uci[var_name])+1)
            }else if(type=='logit'){ # [0,1]
                q_eta = par[var_name]
                par[var_name] = plogis(q_eta)
                se[var_name] = dlogis(q_eta) * se[var_name]

                lci[var_name] = plogis(lci[var_name])
                uci[var_name] = plogis(uci[var_name])
            }else{
                cat(sprintf('No transformation is done for parameter %s with type %s\n', var_name, type))
            }
            ix = which(names(par)==var_name)
            names(se)[ix] = names(par)[ix] = names(lci)[ix] = names(uci)[ix] = names(trans_vars)[i]
        }
    }
    z = par/se
    p = 1 - pchisq(z^2, 1)
    res$estimates = cbind(estimate=par,se=se,z=z,p=p,lci=lci,uci=uci)
    res
}
