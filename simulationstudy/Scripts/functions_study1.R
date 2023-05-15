# all functions for study 1

library(tidyr)
library(dplyr)
library(purrr)
library(lavaan)
library(ggplot2)

#### create_pop_lavaan #######################################################

# Create_pop_lavaan() produces model syntax according to population values 
# for data generation.
create_pop_lavaan <- function(condition,
                              wSigma,
                              Phi,
                              Psi) {
  
  # Generate default variable names
  name_var <- LETTERS[1:2]
  
  # Create matrix of names for observed variable, within, and between components
  name_obs <- sapply(name_var, paste0, 1:condition[["waves"]])
  name_within <- sapply(name_var, function(x) {
    paste0("w", x, 1:condition[["waves"]])
  })
  name_RI <- paste0("RI_", name_var)
  
  # Create population parameter table
  pop_tab <- rbind(
    lav_RI(condition = condition, name_RI, name_obs),
    pop_RI_var(condition = condition, name_RI),
    pop_RI_cor(condition = condition, name_RI),
    pop_within(condition = condition, name_within, name_obs),
    pop_lagged(condition = condition, name_within, Phi),
    pop_within_var1(condition = condition, name_within),
    pop_within_cov1(condition = condition, name_within, wSigma),
    pop_within_var2(condition = condition, name_within, Psi),
    pop_within_cov2(condition = condition, name_within, Psi),
    pop_MEx(condition = condition, name_obs),
    pop_MEy(condition = condition, name_obs)
  )
  rownames(pop_tab) <- NULL
  
  
  
  # Create lavaan syntax
  pop_synt <- paste0( # Paste over parameters
    paste0( # Paste over columns
      pop_tab[, 1],
      pop_tab[, 2],
      pop_tab[, 3],
      pop_tab[, 4],
      pop_tab[, 5]
    ),
    collapse = "\n"
  )
  
  # Create condition list with extra element space
  list(
    waves = condition[["waves"]],
    sample_size = condition[["sample_size"]],
    ME_proportion= condition[["ME_prop"]],
    RI_var = condition[["RI_var"]],
    within_var = condition[["within_var"]],
    crosslagged = condition[["crosslagged"]],
    pop_synt = pop_synt,
    pop_tab = pop_tab
  )
}

lav_RI <- function(condition, name_RI, name_obs) {
  lhs <- rep(name_RI, each = condition[["waves"]])
  op <- rep("=~", times = 2 * condition[["waves"]])
  pv <- rep("1", times = 2 * condition[["waves"]]) # HERE MAYBE ??
  con <- rep("*", times = 2 * condition[["waves"]])
  rhs <- c(unlist(name_obs))
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}

pop_RI_var <- function(condition, name_RI) {
  lhs <- rhs <- name_RI
  op <- "~~"
  con <- "*"
  pv <- condition[["RI_var"]]
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


pop_RI_cor <- function(condition, name_RI) {
  lhs <- name_RI[1]
  rhs <- name_RI[2]
  op <- "~~"
  con <- "*"
  pv <- condition[["RI_cov"]]
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


pop_within <- function(condition, name_within, name_obs) {
  lhs <- c(name_within)
  op <- "=~"
  pv <- condition[["within_var"]]
  con <- "*"
  rhs <- c(name_obs)
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}



pop_lagged <- function(condition, name_within, Phi) {
  lhs <- rep(c(t(name_within))[-(1:2)], each = 2)
  op <- "~"
  con <- "*"
  pv <- c(t(Phi))
  free <- FALSE
  rhs <- c(apply(name_within[-condition[["waves"]], ], 1, rep, times = 2))
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}

pop_within_var1 <- function(condition, name_within) {
  lhs <- rhs <- c(t(name_within[1, ]))
  op <- "~~"
  con <- "*"
  pv <- condition[["within_var"]]
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}



pop_within_cov1 <- function(condition, name_within, wSigma) {
  lhs <- name_within[1, "A"]
  rhs <- name_within[1, "B"]
  op <- "~~"
  con <- "*"
  pv <- c(wSigma[lower.tri(wSigma)]) # Get covariances
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


pop_within_var2 <- function(condition, name_within, Psi) {
  lhs <- rhs <- c(name_within[-1, ])
  op <- "~~"
  con <- "*"
  pv <- rep(diag(Psi), each = (condition[["waves"]] - 1))
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}



pop_within_cov2 <- function(condition, name_within, Psi) {
  lhs <- name_within[-1, 1]
  rhs <- name_within[-1, 2]
  op <- "~~"
  con <- "*"
  pv <- c(Psi[lower.tri(Psi)])
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}



pop_MEx <- function(condition, name_obs) {
  lhs <- rhs <-name_obs[,"A"]
  op <- "~~"
  pv <- unlist(condition[["MEx_pop"]])
  free <- FALSE
  con <- "*"
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}

pop_MEy <- function(condition, name_obs) {
  lhs <- rhs <- name_obs[,"B"]
  op <- "~~"
  pv <- unlist(condition[["MEy_pop"]])
  free <- FALSE
  con <- "*"
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}



##### create_lavaan ###########################################################


# create_lavaan() produces model syntax for RI-CLPM estimation with
# ME variances fixed to different values
create_lavaan <- function(conditions,
                          constraints)
{
  
  # Generate default variable names
  name_var <- LETTERS[24:25] # X and Y
  
  
  # Create matrix of names for observed variable, within, measurement error, and between components
  name_obs <- sapply(name_var, paste0, 1:conditions[["waves"]])
  name_within <- sapply(name_var, function(x) {
    paste0("w", x, 1:conditions[["waves"]])
  })
  
  
  name_RI <- paste0("RI_", name_var)
  
  
  # Create estimation parameter table
  est_tab <- rbind(
    lav_RI(conditions = conditions, name_RI, name_obs),
    est_RI_var(name_RI),
    est_RI_cov(name_RI), 
    est_within(name_within, name_obs, constraints),
    est_lagged(conditions = conditions, name_within, constraints),
    est_within_var1(name_within, constraints),
    est_within_cov1(name_within, constraints),
    est_within_var2(conditions = conditions, name_within, constraints),
    est_within_cov2(conditions = conditions, name_within, constraints),
    fix_MEx(conditions = conditions, name_obs),
    fix_MEy(conditions = conditions, name_obs)
  )
  rownames(est_tab) <- NULL
  
  # Create lavaan syntax
  
  est_synt <- paste0( # Paste over parameters
    paste0( # Paste over columns
      est_tab[, 1],
      est_tab[, 2],
      est_tab[, 3],
      est_tab[, 4],
      est_tab[, 5]
    ),
    collapse = "\n"
  )
  
  
  # Create condition list with extra element space
  list(
    waves = conditions[["waves"]],
    ME_proportion_in_x = conditions[["MEpropx"]],
    ME_proportion_in_y = conditions[["MEpropy"]],
    est_synt = est_synt,
    est_tab = est_tab,
    estimates = NA,
    uncertainty = NA,
    errors = NA,
    not_converged = NA,
    inadmissible = NA
  )
}


lav_RI <- function(conditions, name_RI, name_obs) {
  lhs <- rep(name_RI, each = conditions[["waves"]])
  op <- rep("=~", times = 2 * conditions[["waves"]])
  pv <- rep("1", times = 2 * conditions[["waves"]])
  con <- rep("*", times = 2 * conditions[["waves"]])
  rhs <- c(unlist(name_obs))
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


est_RI_var <- function(name_RI) {
  lhs <- rhs <- name_RI
  op <- "~~"
  con <- "*"
  pv <- "NA"
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}

est_RI_cov <- function(name_RI) {
  lhs <- name_RI[[1]]
  rhs <- name_RI[[2]]
  op <- "~~"
  con <- "*"
  pv <- "NA"
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


est_within <- function(name_within, name_obs, constraints) {
  lhs <- c(name_within)
  op <- "=~"
  if (constraints == "stationarity") {
    pv <- "NA"
    free <- TRUE
  } else {
    pv <- "1"
    free <- FALSE
  }
  con <- "*"
  rhs <- c(name_obs)
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


est_lagged <- function(conditions, name_within, constraints) {
  lhs <- rep(c(t(name_within))[-(1:2)], each = 2)
  op <- "~"
  con <- "*"
  if (constraints == "lagged" || constraints == "within" ||
      constraints == "stationarity"){
    pv <- c("a", "b", "c", "d") # Labels for constraints
    free <- TRUE
    rhs <- c(apply(name_within[-conditions[["waves"]], ], 1, rep, times = 2))
  } else {
    pv <- "NA"
    free <- TRUE
    rhs <- c(apply(name_within[-conditions[["waves"]], ], 1, rep, times = 2))
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE))
}



est_within_var1 <- function(name_within, constraints) {
  lhs <- rhs <- c(t(name_within[1, ]))
  op <- "~~"
  con <- "*"
  if (constraints == "stationarity") {
    pv <-     pv <- c(
      paste0("varx"),
      paste0("vary"))
    free <- TRUE
    
  } else {
    pv <- "NA"
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


est_within_cov1 <- function(name_within, constraints) {
  lhs <- name_within[1, 1]
  rhs <- name_within[1, 2]
  op <- "~~"
  con <- "*"
  if (constraints == "stationarity") { # Label
    pv <- "cov1"
    free <- TRUE
  } else { # Freely estimate
    pv <- "NA"
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE))
}


est_within_var2 <- function(conditions, name_within, constraints) {
  lhs <- rhs <- c(name_within[-1, ])
  op <- "~~"
  con <- "*"
  if (constraints == "residuals" || constraints == "within") { # Constrain over time
    pv <- rep(c("rvarx", "rvary"), each = (conditions[["waves"]] - 1))
    free <- TRUE
  } #else if (constraints == "stationarity") {
  # pv <- c(
  #  paste0("rvarx", 2:conditions[["waves"]]),
  # paste0("rvary", 2:conditions[["waves"]]))
  #  free <- TRUE
  #  }
  else { # Freely estimate
    pv <- "NA"
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE))
}


est_within_cov2 <- function(conditions, name_within, constraints) {
  lhs <- name_within[-1, 1]
  rhs <- name_within[-1, 2]
  op <- "~~"
  con <- "*"
  if (constraints == "residuals" || constraints == "within") { # Constrain over time
    pv <- "rcov"
    free <- TRUE
  }# else if (constraints == "stationarity") { # Label
  #   pv <- paste0("rcov", 2:conditions[["waves"]])
  #    free <- TRUE
  else { # Freely estimate
    pv <- "NA"
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE))
}


fix_MEx <- function(conditions, name_obs){
  lhs <- rhs <- name_obs[, "X"]
  op <- "~~"
  con <- "*"
  pv <- unlist(conditions[["MEx"]])
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE))
}


fix_MEy <- function(conditions, name_obs){
  lhs <- rhs <- name_obs[, "Y"]
  op <- "~~"
  con <- "*"
  pv <- unlist(conditions[["MEy"]])
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE))
}



# sensRICLPM input needs to be a data frame with all observations of one variable X in the first few columns and the other variable Y in the second half of columns
# data: (example:) if there are three waves of data, column 1-3 should contain the observed scores on X, and column 4-6 should contain the observed scores on Y)
# ME_prop_x: the steps in the proportions of measurement error in the observed scores of X that will be included in the sensitivity analysis
# ME_prop_y: the steps in the proportions of measurement error in the observed scores of Y that will be included in the sensitivity analysis
# constraints: can be set to "none", "stationarity", "lagged", "within", "residual", and "equal ME variances"

sensRICLPM.v2 <- function(data, 
                          ME_prop_x = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9), 
                          ME_prop_y = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9), 
                          # var_name_x = "X",
                          # var_name_y = "Y",
                          p,
                          constraints = "none"){
  
  p()
  # add input checks
  
  conditions <- data.frame(MEpropx = ME_prop_x, 
                           MEpropy = ME_prop_y) # only equal ME proportions
  
  # add number of waves and sample size to conditions
  conditions$waves <- timepoints(data)
  conditions$samplesize <- nrow(data)
  
  # consistent names for observed variables
  var_names <- sapply(c("X", "Y"), paste0, 1:conditions[["waves"]][[1]])
  names(data) <- var_names
  
  conditions$constraints <- constraints
  
  # Compute ME variances from proportions
  conditions$MEx <- purrr::map(conditions$MEpropx, .f = compute_ME_var_x, data, constraints)
  conditions$MEy <- purrr::map(conditions$MEpropy, .f = compute_ME_var_y, data, constraints)
  
  
  # generate lavaan syntax for all combinations of ME proportions
  syntaxes <- lapply(X = asplit(conditions, MARGIN = 1), 
                     FUN = create_lavaan, 
                     constraints = constraints
  )
  # extract all model syntaxes
  conditions$syntax <- purrr::map(syntaxes, function(i){
    pluck(i, "est_synt")
  })
  
  
  message(rlang::format_error_bullets(c(
    i = "Fitting lavaan models..."
  )))
  
  # fit all models
  fits <- suppressWarnings(purrr::map(conditions$syntax, lavaan::lavaan, data))
  
  
  message(rlang::format_error_bullets(c(
    i ="Saving estimates...")))
  
  # check convergence
  conditions$converged <- suppressWarnings(
    purrr::map_lgl(fits, function(i){
      lavInspect(i, what ="converged")}))
  
  # check admissibility
  conditions$admissible <- suppressWarnings(
    purrr::map_lgl(fits, function(i){
      lavInspect(i, what ="post.check")}))
  
  
  # save parameter estimates
  conditions$results <- purrr::map(fits, function(i) parameterEstimates(i, remove.nonfree = TRUE))
  
  # create dataframe for estimates of parameters
  
  parameterlist <- list(
    
    alphas <- purrr::map(conditions$results, function(i) # extract all alphas (autoregressive)
      dplyr::filter(i, startsWith(lhs, prefix = "wX"), op == "~", startsWith(rhs, prefix = "wX"))
    ),
    betas <- purrr::map(conditions$results, function(i)  # extract all betas (cross-lagged)
      dplyr::filter(i, startsWith(lhs, prefix = "wY"), op == "~", startsWith(rhs, prefix = "wX"))
    ),
    gammas <- purrr::map(conditions$results, function(i)  # extract all gammas (cross-lagged)
      dplyr::filter(i, startsWith(lhs, prefix = "wX"), op == "~", startsWith(rhs, prefix = "wY"))
    ),
    deltas <- purrr::map(conditions$results, function(i)  # extract all deltas (autoregressive)
      dplyr::filter(i, startsWith(lhs, prefix = "wY"), op == "~", startsWith(rhs, prefix = "wY"))
    ),
    RIX <- purrr::map(conditions$results, function(i) # extract all random intercept variance estimates for X
      dplyr::filter(i, startsWith(lhs, prefix = "RI_X"), op == "~~", startsWith(rhs, prefix = "RI_X"))
    ),
    RIY <- purrr::map(conditions$results, function(i) # extract all random intercept variance estimates for Y
      dplyr::filter(i, startsWith(lhs, prefix = "RI_Y"), op == "~~", startsWith(rhs, prefix = "RI_Y"))
    ),
    wX <- purrr::map(conditions$results, function(i) # extract all within variance estimates for X
      dplyr::filter(i, startsWith(lhs, prefix = "wX"), op == "~~", startsWith(rhs, prefix = "wX"))
    ),
    wY <- purrr::map(conditions$results, function(i) # extract all within variance estimates for Y
      dplyr::filter(i, startsWith(lhs, prefix = "wY"), op == "~~", startsWith(rhs, prefix = "wY"))
    ),
    RIcov <- purrr::map(conditions$results, function(i) # extract all alphas (autoregressive)
      dplyr::filter(i, startsWith(lhs, prefix = "RI_X"), op == "~~", startsWith(rhs, prefix = "RI_Y"))
    ),
    wcov <- purrr::map(conditions$results, function(i) # extract all alphas (autoregressive)
      dplyr::filter(i, startsWith(lhs, prefix = "wX"), op == "~~", startsWith(rhs, prefix = "wY"))
    ))
  
  # create list containing estimates for all lagged effects, SEs & CIs (grouped by wave)

  outputtable <- create_outputtable(parameterlist, conditions)
  
  result <- list(outputtable = outputtable)
 
} 

### measurement error computations & others #############################################

## ME_var_pop2() computes the population values for measurement error variances
#  from the target ME proportion, the randon intercept variance, the within variance,
# and the number of measurement timepoints.
ME_var_pop2 <- function(ME.p, RI.var, W.var, waves){
  matrix(
    (ME.p%*%RI.var + ME.p%*%W.var)/(1 - ME.p), nrow=waves, byrow = T)
}



# compute_ME_var_x() computes values that ME variances of x need to be fixed to from different proportions of ME & total #variance in the data
compute_ME_var_x <- function(ME_prop_x, data, constraints){
  waves <- timepoints(data)
  
  if (constraints == "equal_ME_variances"){ # ME variances constrained to be equal over waves
    ME_var_x_mat <-matrix(
      (sum(diag(cov(data[1:waves]))) / waves) %*% ME_prop_x, # computes ME variance from average total variance in x
      nrow = length(diag(cov(data[1:waves]))), byrow=T)
  }else{
    ME_var_x_mat <-matrix(
      diag(cov(data[1:waves])) %x% ME_prop_x, # computes ME variance for x at each wave
      nrow=length(diag(cov(data[1:waves]))), byrow=T)
  }
}


# compute_ME_var_y() computes values that ME variances of y need to be fixed to from different proportions of ME & total #variance in the data
compute_ME_var_y <- function(ME_prop_y, data, constraints){
  waves <- timepoints(data)
  if (constraints == "equal_ME_variances"){ # ME variances constrained to be equal over waves
    matrix(
      (sum(diag(cov(data[-(1:waves)]))) / waves) %*% ME_prop_y, # computes ME variance from average total variance in y
      nrow = length(diag(cov(data[-(1:waves)]))), byrow=T)
  }else{
    matrix(
      diag(cov(data[-(1:waves)])) %x% ME_prop_y,  # computes ME variance for y at each wave
      nrow=length(diag(cov(data[-(1:waves)]))), byrow=T)
  }
}


# timepoints() saves the number of measurement waves
timepoints<- function(d){
  nwaves <- length(diag(cov(d)))/2
  return(nwaves)
}



# create_outputtable() creates the outputtable from a lavaan.data.frame (list) object 
create_outputtable <- function(parameter, conditions){
  
  parframe <- as.data.frame(purrr::map(parameter, create_plotdata)) # grab relevant estimates
  
  # consitent column names: 
  parnames <- c("alpha", "beta", "gamma", "delta", "RIX", "RIY", "wX", "wY", "RIcov", "wcov")
  estnames <- c("est", "se", "pvalue", "cilower", "ciupper")
  colnames(parframe) <-  as.vector(outer(estnames, parnames, paste, sep="_")) # add column names
  
  plotframe <- cbind(conditions, parframe) # join condition table and data frame containing estimates
  
  plotframe <- tidyr:: unnest_wider(plotframe,!c(MEpropx, MEpropy, waves, samplesize, constraints, MEx, MEy, converged, 
                                                 admissible, results, syntax), names_sep = ".") # unnest wave estimates
  
  
  outputtable<- pivot_longer(plotframe, names_to = c(".value", "parameter", "wave"), !c(MEpropx, MEpropy, waves, samplesize, constraints, MEx, MEy, converged, 
                                                                                          admissible, results, syntax),
                           names_pattern = "(.+)_(.+).(\\d)")  # make longer format
  
}

# create_plotdata() creates a long format table from a lavaan.data.frame object, including estimates for lagged effects, SEs, pvalues, and confidence intervals grouped by wave, 
create_plotdata <- function(parameter){
  
  plotdata <- cbind(
    est <- purrr::map(parameter, function(i){
      pluck(i, "est") # extract estimates
    }),
    se <- purrr::map(parameter, function(i){
      pluck(i, "se")  # extract SEs
    }),
    pvalue <- purrr::map(parameter, function(i){
      pluck(i, "pvalue")  # extract pvalues
    }),
    ci_lower <- purrr::map(parameter, function(i){
      pluck(i, "ci.lower") # extract lower bounds of confidence intervals
    }),
    ci_higher <- purrr::map(parameter, function(i){
      pluck(i, "ci.upper") # extract upper bounds of confidence intervals
    })
  )
  plotdata <- as.data.frame(plotdata)

}


# extra function
ME_var_pop2 <- function(ME.p, RI.var, W.var, waves){
  matrix(
    (ME.p%*%RI.var + ME.p%*%W.var)/(1 - ME.p), nrow=waves, byrow = T)
}

