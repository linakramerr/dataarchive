# data generation study 1 script

# creates 100 data sets for each of 144 conditions (14.000 data sets in total)

library(tidyr)
library(purrr)
library(furrr)
library(lavaan)
library(ggplot2)
library(future)
library(tictoc)
library(progressr)

set.seed(15)


future::plan(multisession, workers = 12)
Sys.setenv("OMP_THREAD_LIMIT" = 46)
### make conditions ############################################################

simconds <- list("ME_prop" = c(0, 0.2, 0.4, 0.6),
                 "sample_size" = c(300, 1000), 
                 "waves" = c(3, 8), 
                 "crosslagged" = c("beta =.2, gamma = 0", "beta =-.2, gamma = .2", "beta =-.2, gamma = -.2"), 
                 "between_var_prop" = c(70, 50, 30))

simconds <- expand.grid(simconds)


###  Create these variables for each condition:  
### RI_var,
### RI_cov,  
### (co-) variances of within components at wave 1, 
### lagged parameters (alpha, beta, gamma, delta),
### Population measurement error variances required to reach each proportion

# RI_var (the random intercept variances) 
# (proportion of between variance not in the total variance VAR(within) + VAR(between) + VAR(error), 
# but in the reliable variance VAR(within) + VAR(between))
simconds$RI_var[simconds$between_var_prop == "70"] <- 1 *1.3
simconds$RI_var[simconds$between_var_prop == "50"] <- 1*1
simconds$RI_var[simconds$between_var_prop == "30"] <- 1*0.7

# RI_cov (covariance (here also the correlation) between the random intercepts
simconds$RI_cov <- 0.5

# within variances and covariacnes at timepoint 1
simconds$within_var <- 1
simconds$within_cov <- 0.3

# alpha, beta, gamma, delta
#simconds$alpha <- 0.4
#simconds$beta[simconds$crosslagged == "beta =.2, gamma = 0"] <- 0.2
#simconds$beta[simconds$crosslagged == "beta =-.2, gamma = .2" | simconds$crosslagged =="beta =-.2, gamma = -.2"] <- -0.2
#simconds$gamma[simconds$crosslagged == "beta =.2, gamma = 0"] <- 0
#simconds$gamma[simconds$crosslagged == "beta =-.2, gamma = .2"] <- 0.2
#simconds$gamma[simconds$crosslagged == "beta =-.2, gamma = -.2"] <- -0.2
#simconds$delta <- 0.3

# Measurement error variances for each proportion 
simconds$MEx_pop <- purrr::pmap(list(simconds$ME_prop, simconds$RI_var, simconds$within_var, simconds$waves),
                               ME_var_pop2) # function ME_var_pop2() computes the necessary ME variance to reach a certain proportion
simconds$MEy_pop <- purrr::pmap(list(simconds$ME_prop, simconds$RI_var,simconds$within_var, simconds$waves),
                                ME_var_pop2) 

wSigma <- matrix(c(1, 0.3, 0.3, 1), ncol=2, byrow=T)

# the conditions differ in their phi (& subsequent psi) matrices, so we are splitting up the conditions into subsets that share the same phi
simconds_1 <- simconds %>% subset(crosslagged == "beta =.2, gamma = 0") 
simconds_2 <- simconds %>% subset(crosslagged == "beta =-.2, gamma = .2")
simconds_3 <- simconds %>% subset(crosslagged == "beta =-.2, gamma = -.2")

Phi_1 <- matrix(c(.4, .2, 0, .3), ncol=2, byrow=F) # alpha, beta, gamma, delta
Phi_2 <- matrix(c(.4, -.2, .2, .3), ncol=2, byrow=F)
Phi_3 <- matrix(c(.4, -.2, -.2, .3), ncol=2, byrow=F)

Psi_1 <- matrix((diag(length(wSigma)) - Phi_1 %x% Phi_1) %*% c(wSigma), nrow = nrow(Phi_1), byrow=T)
Psi_2 <- matrix((diag(length(wSigma)) - Phi_2 %x% Phi_2) %*% c(wSigma), nrow = nrow(Phi_2), byrow=T)
Psi_3 <- matrix((diag(length(wSigma)) - Phi_3 %x% Phi_3) %*% c(wSigma), nrow = nrow(Phi_3), byrow=T)


### Obtain lavaan syntaxes for data generation ################################


# Create lavaan syntax for all conditions:
simconds_1$popsynt <- lapply(X = asplit(simconds_1, MARGIN = 1), # apply create_pop_lavaan() function to all rows
                      FUN = create_pop_lavaan, 
                      wSigma,
                      Phi_1,
                      Psi_1) %>% purrr::map(function(i){ # grab only the lavaan syntax string for each condition and add to condition table
                        pluck(i, "pop_synt")
                      })

simconds_2$popsynt <- lapply(X = asplit(simconds_2, MARGIN = 1), 
                    FUN = create_pop_lavaan, 
                    wSigma,
                    Phi_2,
                    Psi_2) %>% purrr::map(function(i){
                      pluck(i, "pop_synt")
                      })

simconds_3$popsynt <- lapply(X = asplit(simconds_3, MARGIN = 1), 
                    FUN = create_pop_lavaan, 
                    wSigma,
                    Phi_3,
                    Psi_3) %>% purrr::map(function(i){
                      pluck(i, "pop_synt")
                      })

# recombine the previously separated simulation condition table, now inlcuding lavaan syntax in "popsynt" for data generation
simconds <- rbind(simconds_1, simconds_2, simconds_3) 

saveRDS(simconds, file = "simulationconditions.rds") # save condition file

### Generate 100 data sets for each of the 144 conditions #####################

# splitting up data generation for N = 500 and N = 10.000:

simconds_n300 <-  simconds %>% subset(sample_size == 300) %>% as_tibble() 
# frame contains the 72 conditions where  N=500

dat_n300 <- replicate(100, purrr::map(simconds_n300$popsynt, 
                                      lavaan::simulateData, 
                                      meanstructure = F,
                                      sample.nobs = 300), simplify = F)
saveRDS(dat_n300, "dat_n300.rds") 
# this is a list of 100 lists that each contain 72 data sets (in the order of table simcons_n500)
# -> 7200 data sets

simconds_n1000 <- simconds %>% subset(sample_size == 1000) %>% as_tibble()
# frame contains the 72 conditions where N = 10.000


dat_n1000 <- replicate(100, furrr::future_map(simconds_n1000$popsynt,
                                       lavaan::simulateData,
                                       meanstructure = F,
                                       sample.nobs = 1000,
                                       .options= furrr_options(seed = T, scheduling = 2L)), simplify = F)


saveRDS(dat_n1000, "dat_n1000.rds")
# list of 100 lists that each contain 72 data sets (in the order of the table simconds_n10000)
# -> 7200 data sets

listdat_n300 <- unlist(dat_n300, recursive = F)

# function that runs the simulation for a list of data sets
runsims_study1 <- function(x) {
  p <- progressor(steps = length(x))
  
furrr::future_map(x,
                  sensRICLPM.v2, 
                        ME_prop_x = c(0, .2, .4, .6, .8), 
                        ME_prop_y = c(0, .2, .4, .6, .8),
                        p = p,
                        .options = furrr::furrr_options(seed = TRUE)) ### fits 72 data sets to 10 conditions 

}

# run simulation for n = 300
tic()
with_progress({
  simresults_study1_n300 <- runsims_study1(listdat_n300)
})

toc()
saveRDS(simresults_study1_n300, "study1fits_n300.rds")


# run simulation for n = 10.000

listdat_n1000 <- unlist(dat_n1000, recursive = F)

tic()
with_progress({
  simresults_study1_n1000 <- runsims_study1(listdat_n1000)
})

toc()
saveRDS(simresults_study1_n1000, "study1fits_n1000.rds")

outputtables <-  c(study1fits.v2_n300, study1fits_n1000) %>%  
  purrr::map(function(i){
  pluck(i, "outputtable") # extract estimates
})# join results for both sample sizes

saveRDS(outputtables, "study1_fits.rds")


res_study1 <- data.frame(matrix(nrow=14400))
res_study1$id <- c(rep(c(1:72, 73:144), each = 100))


outputtables<- purrr::map(study1_fits, function(i){
  pluck(i, "outputtable") # extract estimates
})

res_study1$admissible <- purrr::map(outputtables, function(i){
  pluck(i, "admissible") # extract estimates
})

res_study1$MEpropx <- purrr::map(outputtables, function(i){
  pluck(i, "MEpropx") # extract estimates
})

res_study1$MEpropy <- purrr::map(outputtables, function(i){
  pluck(i, "MEpropy") # extract estimates
})
res_study1$parameters <- purrr::map(outputtables, function(i){
  pluck(i, "parameter") # extract estimates
})
res_study1$est <- purrr::map(outputtables, function(i){
  pluck(i, "est") # extract estimates
})
res_study1$se <- purrr::map(outputtables, function(i){
  pluck(i, "se") # extract estimates
})
res_study1$pvalue <- purrr::map(outputtables, function(i){
  pluck(i, "pvalue") # extract estimates
})
res_study1$cilower <- purrr::map(outputtables, function(i){
  pluck(i, "cilower") # extract estimates
})
res_study1$ciupper <- purrr::map(outputtables, function(i){
  pluck(i, "ciupper") # extract estimates
})

res_study1$waves <- purrr::map(outputtables, function(i){
  pluck(i, "waves") # extract estimates
})




res_study1$MEprop_sample <- c(rep(simconds_n300$ME_prop, 100), rep(simconds_n1000$ME_prop, 100))
res_study1$crosslagged <- c(rep(simconds_n300$crosslagged, 100), rep(simconds_n1000$crosslagged, 100))
res_study1$between_var_prop <- c(rep(simconds_n300$between_var_prop, 100), rep(simconds_n1000$between_var_prop, 100))
res_study1$RI_var <- c(rep(simconds_n300$RI_var, 100), rep(simconds_n1000$RI_var, 100))
res_study1$within_var <- c(rep(simconds_n300$within_var, 100), rep(simconds_n1000$within_var, 100))
res_study1$RI_cov <- c(rep(simconds_n300$RI_cov, 100), rep(simconds_n1000$RI_cov,  100))
res_study1$sample_size <- c(rep(simconds_n300$sample_size, 100), rep(simconds_n1000$sample_size,100))
res_study1$within_cov <- c(rep(simconds_n300$within_cov, 100), rep(simconds_n1000$within_cov, 100))

l_res_study1 <- res_study1 %>% unnest_longer(col = c(admissible, MEpropx, MEpropy, parameters, est, cilower, ciupper, pvalue, se, waves))
