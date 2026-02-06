## This script will be used to run all of the 
## Bayesian cp models for the pfactor-nonlinear manuscript

## load library(s)
#library(doParallel)
library(rstan)

# > names(combined)
# [1] "src_subject_id"      "eventname"
# [3] "site_id_l"           "rel_family_id"
# [5] "interview_age"       "asr_scr_thought_r"
# [7] "asr_scr_attention_r" "asr_scr_internal_r"
# [9] "asr_scr_external_r"

# Now run through all of these model options in a for parallel loop
mod.dv <- c("asr_scr_internal_r", "asr_scr_external_r", "asr_scr_attention_r", "asr_scr_thought_r")
mod.iv <- mod.dv
iter <- 1
all.mods <- expand.grid(mod.dv, mod.iv, iter)
all.mods$Var1 <- as.character(all.mods$Var1)
all.mods$Var2 <- as.character(all.mods$Var2)
all.mods <- all.mods[-which(as.character(all.mods$Var1)==as.character(all.mods$Var2)),]
all.mods <- unique(all.mods)
rowID <- as.integer(commandArgs(1))
i <- all.mods[rowID,3]
## Now go through the wave value
## create the data
if(i == 1){
  in.dat <- read.csv("./data/ABCD_CBCL_BL.csv") ## This file path needs to be updated with the path to the adult's csv file
}
# if(i == 2){
#   in.dat <- read.csv("./data/ABCD_CBCL_Y1.csv")
# }
# if(i == 3){
#   in.dat <- read.csv("./data/ABCD_CBCL_Y2.csv")
# }
# if(i == 4){
#   in.dat <- read.csv("./data/ABCD_CBCL_Y3.csv")
# }
# if(i > 1){
#   in.dat.bp <- read.csv("./data/ABCD_CBCL_BL.csv")
#   in.dat.bp.family <- in.dat.bp[,c("src_subject_id", "site_id_l", "rel_family_id")]
#   in.dat <- merge(in.dat, in.dat.bp.family, by=c("src_subject_id", "site_id_l"), suffixes = c(".x",""))
# }
tmp.dat <- in.dat
tmp.dat <- in.dat[complete.cases(tmp.dat[,c(all.mods$Var1[rowID], all.mods$Var2[rowID], "interview_age")]),]
data_jags <- list(
  N = nrow(tmp.dat),
  x = tmp.dat[,all.mods$Var2[rowID]],
  y = tmp.dat[,all.mods$Var1[rowID]],
  x2 = tmp.dat$interview_age,
  N_group = as.numeric(factor(tmp.dat$site_id_l)),
  count_group = length(unique(as.numeric(factor(tmp.dat$site_id_l)))),
  N_group2 = as.numeric(factor(tmp.dat$rel_family_id)),
  count_group2 = length(unique(as.numeric(factor(tmp.dat$rel_family_id))))
)

## Now make the scaled values here
all.dat <- data_jags
file.out <- paste("./data/brmsModsOut/model_rawX_NB_NOcp_allmods_", rowID, ".RDS", sep='')
stanmonitor = c("alpha", "beta", "phi", "sigma_p", "sigma_p2", "log_lik")
if(!file.exists(file.out)){
  result_case = stan(file="./scripts/stan_models/quick_c0_test.stan", 
                     data = all.dat, cores=2,chains=2, refresh = 10, 
                     pars = stanmonitor, 
                     iter=30000, warmup = 20000, thin = 2,control = list(max_treedepth=9))
  saveRDS(result_case, file.out)
}else{
  print("Done")
  result_case <- readRDS(file.out)
}

summary(do.call(rbind, 
                args = get_sampler_params(result_case, inc_warmup = FALSE)),
        digits = 2)
## Now do the logLik loo call here
loo::elpd(result_case)
