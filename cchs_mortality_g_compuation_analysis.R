##################
## Below are example codes for g-computation analysis
##################

## Set up the stage--please EDIT ME before running the codes
##################
library(gfoRmula)
library(data.table)
library(Hmisc)
library(RColorBrewer)
##################

## Read in the data
##################
## This is a mock dataset of 1000 individuals without meaningful association
## It is used to test the codes and estimate the amount of running time
## The codes to generate this dataset is in "cchs_mortality_g_compuation_data_generating.R"
dt <- read.csv("1000individual simulated dataset.csv", as.is = TRUE)
##################

## Intervention
##################
## function for relative reduction with lower limit at 1.8 for change
relative_reduction <- function(newdf, pool, intvar, intvals, time_name, t){ 
  temp <- newdf[[intvar]] # save an original copy
  newdf[, (intvar) := eval(as.name(intvar)) + log(intvals[[1]]) ] # reduce exposure 
  
  loc1 <- which( newdf[[intvar]] < intvals[[2]] & newdf[[intvar]] >= (intvals[[2]] + log(intvals[[1]])) )
  newdf[loc1, (intvar) := intvals[[2]]] # change back to 1.8 if reduced below
  
  loc2 <- which(newdf[[intvar]] < (intvals[[2]] + log(intvals[[1]])))
  newdf[loc2, (intvar) := temp[loc2]]   # change back to the original value
}

APvar.threshold.1 <- log(12)      # US national standards
APvar.threshold.2 <- log(8.8)     # Canadian national standards
APvar.threshold.3 <- log(8.8*0.8) # 20% below Canadian national standards
APvar.threshold.4 <- log(5)       # revised WHO guideline value
APvar.threshold.5 <- log(4)       # 1 ug/m3 below revised WHO guideline
APvar.threshold.6 <- log(1.8)     # presumed background conc ug/m3

interventions <- list(list(c(threshold, -Inf, APvar.threshold.1)), # all time-points with levels above CA standard set to 12
                      list(c(threshold, -Inf, APvar.threshold.2)), # all time-points with levels above US standard set to 8.8
                      list(c(threshold, -Inf, APvar.threshold.3)), # all time-points with levels above theoretical minimum risk level of GDB set to 80% below Canadian standards
                      list(c(threshold, -Inf, APvar.threshold.4)), # all time-points with levels above theoretical minimum risk level of GDB set to 5
                      list(c(threshold, -Inf, APvar.threshold.5)), # all time-points with levels above theoretical minimum risk level of GDB set to 4
                      list(c(relative_reduction, 0.95, APvar.threshold.6)), # all time-points set to 95% of its original value
                      list(c(relative_reduction, 0.9, APvar.threshold.6)) # all time-points set to 90% of its original value
)

int_descript <- c("Lower or equal to 12", "Lower or equal to 8.8", 
                  "Lower or equal to 80% of 8.8",                
                  "Lower or equal to 5", "Lower or equal to 4",
                  "Reduced to 95% of original",
                  "Reduced to 90% of original"
)
##################


## Main analysis without bootstrap
##################
## set the exposure
dt$APvar <- log(dt$PM3yr)
expsr <- "APvar" # log transformed exposure

## the length of simulated follow-up time 
time_points <- max(dt$time) + 1 # length of follow-up till the longest observation--"time" starts at 0

## Create restricted cubic splines for age--3 internal knots for 25-85 age cohort
rcs <- data.frame(rcspline.eval(dt$age, nk=5, inclx = TRUE)) 
names(rcs) <- paste0("AGE", 1:4)
dt <- cbind(dt, rcs)

## transform categorical variables into factors
ft.vb <- c("bmi5c", "smk6c", "alc6c", "fvg3c", "exc3c", "lbs3c", "edu4c", 
           "abor", "ms3c", "airsh", "urban")
for (i in ft.vb) {
  dt[, i] <- as.factor(dt[, i])
}

## gformula package requires data.table
setDT(dt)

## list of time varying covarites
covnames <- c("csize", "incqu_new_imp", "ethcon", "depriv", "instab", "depend", expsr)
covtypes <- c("bounded normal", "bounded normal", "bounded normal", "bounded normal",
              "bounded normal", "bounded normal",
              "normal")

covlabels <- covnames
names(covlabels) <- c("Community size", "Income quintile", "Ethnic concentration", "Material deprivation", 
               "Residential instability", "Population dependency", "PM2.5 concentration")

## list of baseline covariates for covariate model and outcome model
basecovs.all <- c("yrsincan", "AGE1", "AGE2", "AGE3", "AGE4", 
                  "abor", "ms3c", "airsh", "urban",
                  "bmi5c", "smk6c", "alc6c", "fvg3c", "exc3c",
                  "lbs3c", "edu4c",
                  "DHH_SEX", "immigrant", "vismin")


# Outcome model
f0 <- reformulate(c(expsr,
                    "DHH_SEX", "AGE1", "AGE2", "AGE3", "AGE4", 
                    "abor", "ms3c", "airsh", "urban", 
                    "immigrant", "I(yrsincan * immigrant)", 
                    "bmi5c", "smk6c", "alc6c", "fvg3c", "exc3c",
                    "vismin",   
                    "lbs3c", "edu4c",
                    "incqu_new_imp", "csize", "depend", "depriv", "ethcon", "instab",
                    paste0(expsr, " : as.factor(time)"), "time", "I(time*time)"
), response = "nonacc1")

## time-varying covariates could depend on time, baseline covariates excluing behavior ones, 
## previous measurement of the time-varying covariate being simulated
f.cov <- sapply(1:length(covnames), function(a) {
  reformulate(c(
    basecovs.all[-c(1)], "I(yrsincan * immigrant)", # yrsincan need to be modeled differently
    paste0("lag1_", covnames[a]), # include previous value of the time-varying covariate of interest
    if (a > 1) {covnames[1:(a-1)]},
    "time", "I(time * time)"),
    response = covnames[a])
})

## main analysis: gformula till the end of max follow-up time without bootstrapping
gform_1 <- gformula_survival(
  obs_data = dt, 
  # outcome_type = "survival",
  id = "uniqid",
  time_points = time_points,
  time_name = "time",
  covnames = covnames,
  outcome_name = "nonacc1",
  ymodel = f0, # this is where I specified the outcome model to use--saturated one
  covtypes = covtypes,
  covparams = list(covmodels = f.cov), # this is where I specified the covariate model to use
  intvars = list(expsr, expsr, expsr, expsr, expsr, expsr, expsr),
  ref_int = 0, # so that the reference is the natural course
  interventions = interventions,
  int_descript = int_descript,
  histories = c(lagged),
  histvars = list(covnames),
  basecovs = basecovs.all,
  sim_data_b = TRUE, # this needs to be false when running bootstrapping
  nsimul = 10000,
  seed = 1234
)
gform_1
write.csv(gform_1$result, "gformula_results.csv", row.names = FALSE)

## below are plots for adjusted survival curve, risk ratio and risk difference without confidence interval
d <- as.data.frame(gform_1$result)
names(d) <- c("k", "Intervention", "NP_risk", "g_form_risk", "Risk_ratio", "Risk_dif")
d$surv <- 1 - d$g_form_risk
sc <- reshape(d[, c("k", "Intervention", "surv", "Risk_ratio", "Risk_dif")], v.names = c("surv", "Risk_ratio", "Risk_dif"), timevar = "Intervention", idvar = "k", direction = "wide")
sc$time <- sc$k + 1
sc <- rbind(sc, c(0, rep(1, length(unique(d$Intervention))*3), 0)) # add the original value and move start time to 1
sc <- sc[order(sc$time), ]

## interventions to plot
plot.ints <- 1:length(int_descript)
names(plot.ints) <- int_descript

## Survival curves
colr.line <- brewer.pal(n = length(plot.ints), name = "Set2")
sl <- names(plot.ints)
plot(sc$time, sc$surv.0, type = "n",
     family = "sans", main = "Adjusted survival curve", 
     ylab = "Survival probability", xlab = "Year since enrollment", cex.lab=1.5)
for (k in length(plot.ints):1) {
  lines(sc$time, sc[, paste0("surv.", plot.ints[k])], type="o", col = colr.line[k], lwd=2, 
        pch=k+15)
}
legend("topright", legend=sl, pch = 15 + (1:length(sl)), col=colr.line, xpd=T, lwd = 2)

## risk ratio
sc <- sc[sc$time>0, ]
colr.line2 <- colr.line[-1] 
sl1 <- paste0("Risk_ratio.", plot.ints[-1])
sl2 <- paste0("Risk_dif.", plot.ints[-1])
par(mar = c(5, 4, 4, 4) + 0.3)
plot(sc$time, sc[, sl1[1]], type = "n", 
     ylim = c(min(sc[, sl1])*0.95, max(sc[, sl1])*1.05),
     family = "sans",
     ylab = "", xlab = "Follow-up year", cex.lab=1.5)
for (k in length(sl1):1) {
  lines(sc$time, sc[, sl1[k]], type="o", col = colr.line2[k], lwd=2, pch= k+16)
}
abline(h=1, lwd=1.5)
mtext("Cumulative incidence ratio", side = 2, line = 3, cex=1.5)
legend("topright", legend=names(plot.ints)[-1], pch = 16 + (1:length(sl1)), col=colr.line2, xpd=T, lwd = 2)

## risk difference
plot(sc$time, sc[, sl2[1]], type = "n", 
     ylim = c(min(sc[, sl2])*1.05, max(sc[, sl2])*0.95),
     family = "sans", 
     ylab = "", xlab = "Follow-up year", cex.lab=1.5)
for (k in length(sl2):1) {
  lines(sc$time, sc[, sl2[k]], type="o", col = colr.line2[k], lwd=2, pch=k+16)
}
abline(h=0, lwd=1.5)
mtext("Risk Difference", side = 2, line = 3, cex=1.5) 
legend("bottomleft", legend=names(plot.ints)[-1], pch = 16 + (1:length(sl2)), col=colr.line2, xpd=T, lwd = 2)
##################

## Bootstrap for main analysis
library(doParallel)
library(dplyr)
##################
## generate a folder for bootstrapping results
outdir_bt <- "bootstrap"
if (!dir.exists(outdir_bt)) dir.create(outdir_bt)

## the length of simulated follow-up time 
time_points <- max(dt$time) + 1 # length of follow-up till the longest observation--"time" starts at 0

## list of time varying covarites
covnames <- c("csize", "incqu_new_imp", "ethcon", "depriv", "instab", "depend", expsr)
covtypes <- c("bounded normal", "bounded normal", "bounded normal", "bounded normal",
              "bounded normal", "bounded normal",
              "normal")  

## list of baseline covariates for covariate model and outcome model
basecovs.all <- c("yrsincan", "AGE1", "AGE2", "AGE3", "AGE4", 
                  "abor", "ms3c", "airsh", "urban",
                  "bmi5c", "smk6c", "alc6c", "fvg3c", "exc3c",
                  "lbs3c", "edu4c",
                  "DHH_SEX", "immigrant", "vismin")


# Outcome model
f0 <- reformulate(c(expsr,
                    "DHH_SEX", "AGE1", "AGE2", "AGE3", "AGE4", 
                    "abor", "ms3c", "airsh", "urban", 
                    "immigrant", "I(yrsincan * immigrant)", 
                    "bmi5c", "smk6c", "alc6c", "fvg3c", "exc3c",
                    "vismin",   
                    "lbs3c", "edu4c",
                    "incqu_new_imp", "csize", "depend", "depriv", "ethcon", "instab",
                    paste0(expsr, " : as.factor(time)"), "time", "I(time*time)"
), response = "nonacc1")

## time-varying covariates could depend on time, baseline covariates excluing behavior ones, 
## previous measurement of the time-varying covariate being simulated
f.cov <- sapply(1:length(covnames), function(a) {
  reformulate(c(
    basecovs.all[-c(1)], "I(yrsincan * immigrant)", # yrsincan need to be modeled differently
    paste0("lag1_", covnames[a]), # include previous value of the time-varying covariate of interest
    if (a > 1) {covnames[1:(a-1)]},
    "time", "I(time * time)"),
    response = covnames[a])
})

sub.id <- unique(dt$uniqid)
range <- 1:2 # e.g., 1:100 and 101:200, please use different range for each job so that the seed is different for each bootstrap
job.id <- "all"  # e.g., 1, 2 or all

## IMPORTANT: codes below are for parallel running
ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
foo <- NULL
foo <- foreach(i = range,
               .combine = "rbind", .inorder = FALSE,
               .packages=c("gfoRmula", "data.table", "Hmisc")) %dopar% {
                 # for (i in range){
                 st_ <- Sys.time()
                 set.seed(1234 + i)
                 ids <- as.data.table(sample(sub.id, length(sub.id), replace = TRUE))
                 ids[, bid := 1:nrow(ids)]
                 colnames(ids) <- c("uniqid", "bid")
                 new <- copy(dt)
                 setkey(new, "uniqid")
                 new <- new[J(ids), allow.cartesian = TRUE]
                 new[, uniqid := new$bid]
                 new[, bid := NULL]
                 
                 # Recode rcs(age)
                 new[, AGE1 := NULL]
                 new[, AGE2 := NULL]
                 new[, AGE3 := NULL]
                 new[, AGE4 := NULL]
                 rcs <- data.frame(rcspline.eval(new$age, nk=5, inclx = TRUE))
                 names(rcs) <- paste0("AGE", 1:4)
                 new <- cbind(new, rcs)
                 
                 gform_2 <- gformula_survival(
                   obs_data = new, 
                   id = "uniqid",
                   time_points = time_points,
                   time_name = "time",
                   covnames = covnames,
                   outcome_name = "nonacc1",
                   ymodel = f0, # this is where I specified the outcome model to use--saturated one
                   covtypes = covtypes,
                   ##below are specifications related to glm
                   covparams = list(covmodels = f.cov), # this is where I specified the covariate model to use
                   ##below are for continuous exposure
                   intvars = list(expsr, expsr, expsr, expsr, expsr, expsr, expsr),
                   ref_int = 0, # so that the reference is the natural course
                   interventions = interventions,
                   int_descript = int_descript,
                   histories = c(lagged),
                   histvars = list(covnames), #exposure does not need lag
                   basecovs = basecovs.all,
                   nsimul = 10000,
                   seed = 1234 + i
                 )
                 saveRDS(gform_2$result, file.path(outdir_bt, paste0("gformula_results_bootstrap.", i, ".rds")))
               }
stopCluster(cl)
rm(foo)

## Codes below combined individual bootstrapping results from each cycle into one file
## also checking for duplicates in the code
risks_bs <- NULL
for (i in range) {
  risk_ <- readRDS(file.path(outdir_bt, paste0("gformula_results_bootstrap.", i, ".rds")))
  risks_bs <- rbind(risks_bs, risk_)
}
names(risks_bs) <- c("k", "Intervention", "NP_risk", "g_form_risk", "Risk_ratio", "Risk_dif")

## Check if bootstrapping lead to duplicated results
temp <- risks_bs[c(which(duplicated(risks_bs$g_form_risk)), 
                   which(duplicated(risks_bs$g_form_risk, fromLast = TRUE))), ]
if (nrow(temp) > 0) {
  cat("Error: duplicated results found in bootstrapping", "\n")
  if (sum(unique(temp$k) != 0) > 0 ) {
    cat("Error requires checking: duplicated results not just found in the beginning", "\n")
    write.csv(temp, file.path(outdir_bt, paste0("Error_duplicated_results.csv")), row.names = FALSE) ## Edit me if want to save error message somewhere else
  } else {
    cat("Error fixed: duplicated results only in the beginning", "\n")
  }
}

## Transform bootstrapping results to get CI and SE 
l_risk_bs <- risks_bs %>% group_by(Intervention, k) %>% 
  summarize(NP_risk.bsmean=mean(NP_risk, na.rm=T),
            
            g_form_risk.bsmean=mean(g_form_risk, na.rm=T),
            g_form_risk.se=sd(g_form_risk, na.rm = TRUE), 
            g_form_risk.lb=quantile(g_form_risk, probs = c(0.025), na.rm = TRUE),
            g_form_risk.ub=quantile(g_form_risk, probs = c(0.975), na.rm = TRUE),
            
            Risk_ratio.bsmean=mean(Risk_ratio, na.rm=T),
            Risk_ratio.se=sd(Risk_ratio, na.rm = TRUE), 
            Risk_ratio.lb=quantile(Risk_ratio, probs = c(0.025), na.rm = TRUE),
            Risk_ratio.ub=quantile(Risk_ratio, probs = c(0.975), na.rm = TRUE),
            
            Risk_dif.bsmean=mean(Risk_dif, na.rm=T),
            Risk_dif.se=sd(Risk_dif, na.rm = TRUE), 
            Risk_dif.lb=quantile(Risk_dif, probs = c(0.025), na.rm = TRUE),
            Risk_dif.ub=quantile(Risk_dif, probs = c(0.975), na.rm = TRUE)
  )
l_risk_bs <- as.data.frame(l_risk_bs)

write.csv(l_risk_bs, "gformula_results_bootstrap.csv", row.names = FALSE)
##################