# main_summer.R


# Main script for the mortality burden during summer months.


# Project : paper_heat_cold_mortality_covid
# Author  : Jeremie Boudreault
# Email   : Jeremie.Boudreault@inrs.ca
# Depends : R [v4.4.2], rjutils [v0.1]
# Imports : See below.
# License : CC BY-NC-SA 4.0 DEED


# Packages ---------------------------------------------------------------------


library(data.table)
library(dlnm)
library(splines)
library(mixmeta)
library(rjutils) 


# Functions --------------------------------------------------------------------


source("R/funs/create_lagged_mat.R")
source("R/funs/fit_meta_stepwise_aic.R")


# Imports ----------------------------------------------------------------------


# Daily health data by HR.
health_data <- fread("PATH TO MORTALITY DATA")

# Daily lagged Daymet data by HR.
daymet_data <- fread("PATH TO LAGGED DAYMET WEATHER DATA")

# Daily lagged ECCC data by HR.
eccc_data <- fread("PATH TO LAGGED ECCC WEATHER DATA")

# Meta-predictors for each HR.
meta_data <- fread("PATH TO META-DATA")


# Table of analyses ------------------------------------------------------------


# Number of entries for the table.
n_tbl <- 20L

# Generate initial table.
tbl <- data.table(
    ID      = 1:n_tbl,
    YEARS   = rep("2000:2024", n_tbl),
    YVAR    = rep("COUNT_TOT", n_tbl),
    TVAR    = rep("T_MEAN", n_tbl),
    RELH    = rep(TRUE, n_tbl),
    MONTHS  = rep("5:9", n_tbl),
    LAG     = rep(7L, n_tbl),
    KNOTS   = rep("c(50, 90)", n_tbl),
    DF      = rep(4L, n_tbl)
)

# Updated the analysis.
tbl[1L, NAME := "Main"]
tbl[2L, `:=`(NAME = "2000–2024 (w/o Covid)", YVAR = "COUNT_WO_COVID") ]
tbl[3L, `:=`(NAME = "2000–2019", YEARS = "2000:2019") ]
tbl[4L, `:=`(NAME = "2005–2024 (w Covid)", YEARS = "2005:2024") ]
tbl[5L, `:=`(NAME = "2005–2024 (w/o Covid)", YEARS = "2005:2024", YVAR = "COUNT_WO_COVID") ]
tbl[6L, `:=`(NAME = "Tmean (ECCC)", TVAR = "T_MEAN_ECCC") ]
tbl[7L, `:=`(NAME = "Tmin (Daymet)", TVAR = "T_MIN") ]
tbl[8L, `:=`(NAME = "Tmin (ECCC)", TVAR = "T_MIN_ECCC") ]
tbl[9L, `:=`(NAME = "Tmax (Daymet)", TVAR = "T_MAX") ]
tbl[10L, `:=`(NAME = "Tmax (ECCC)", TVAR = "T_MAX_ECCC") ]
tbl[11L, `:=`(NAME = "Without RH adj.", RELH = FALSE) ]
tbl[12L, `:=`(NAME = "June–August", MONTHS = "6:8") ]
tbl[13L, `:=`(NAME = "April–October", MONTHS = "4:10") ]
tbl[14L, `:=`(NAME = "Lag 3 days", LAG = 3L) ]
tbl[15L, `:=`(NAME = "Lag 10 days", LAG = 10L) ]
tbl[16L, `:=`(NAME = "Knots (25, 75)", KNOTS = "c(25, 75)")]
tbl[17L, `:=`(NAME = "Knots (25, 90)", KNOTS = "c(25, 90)")]
tbl[18L, `:=`(NAME = "Knots (50, 75)", KNOTS = "c(50, 75)")]
tbl[19L, `:=`(NAME = "Seas. adj. 3 df", DF = 3L)]
tbl[20L, `:=`(NAME = "Seas. adj. 5 df", DF = 5L)]

# Look up the results.
tbl


# Parameters -------------------------------------------------------------------


# Loop on all parametrs.
for (tbl_i in c(1:20)) {

# Message.
message("Running ", tbl_i, ".")

# Maximum lag (7 days for summer, 21 days for winter)
max_lag <- tbl[tbl_i, LAG]

# Number of degrees of freedom for date variable (4 df for 5 months)
df_date <- tbl[tbl_i, DF]

# Number of knots for lag dimension.
nk_xlag <- 2

# Position of knots for main exposure.
knots_xvar <- eval(parse(text = tbl[tbl_i, KNOTS]))/100

# Quantiles for the MMT search.
q_mmt <- seq(5L, 95L, by = 1L)/100

# Months for the analysis.
months <- eval(parse(text = tbl[tbl_i, MONTHS]))

# Years for the analysis.
years <- eval(parse(text = tbl[tbl_i, YEARS]))

# Health variables.
hvar <- tbl[tbl_i, YVAR]

# Temperature exposure variable.
tvar <- tbl[tbl_i, TVAR]

# Adjustment for relative humidity.
relh <- tbl[tbl_i, RELH]

# Number of simulations for AN/AF calculation.
nsim <- 1000L
cols_sim <- paste0("SIM", 1:nsim)

# Filename.
fname <- paste0("heat_", tbl_i, "_")

# Verbose.
verbose <- FALSE


# ------------------------------------------------------------------------------
# Step 0 : Data preparation ----------------------------------------------------
# ------------------------------------------------------------------------------


# Merge data health and daymet data.
merged_data <- merge(
    x     = health_data[, c(hvar, "DATE", "HR", "WEEKDAY", "HOL"), with = FALSE], 
    y     = daymet_data, 
    by    = c("HR", "DATE"), 
    all.x = TRUE
)[order(HR, DATE), ]

# Merge data with ECCC.
merged_data <- merge(merged_data, eccc_data, all.x = TRUE, by = c("HR", "DATE"))[order(HR, DATE), ]

# Rename the health variable to "COUNT".
setnames(merged_data, hvar, "COUNT")

# Extract all HR to be treated from 1:9 and 11:16.
hr_sel <- c(1:9, 11:16)
merged_data <- merged_data[HR %in% hr_sel, ]

# Convert weekday and year variables to factor.
merged_data[, WEEKDAY_F := factor(WEEKDAY, levels = 1:7)]
merged_data[, YEAR_F := factor(YEAR)]

# Fix <DOY> for leap years (only for summer analysis).
merged_data[YEAR %% 4 == 0L, DOY := DOY - 1]

# Create lagged matrix of tvar and relh.
tvar_mat_l <- create_lagged_mat(tvar)
if (relh) relh_mat_l <- create_lagged_mat("RELH_MEAN")

# Create final dataset using selected months and years only.
data_final <- merged_data[MONTH %in% months & YEAR %in% years, ]
set(data_final, j = tvar, value = data_final[[paste0(tvar, "0")]])

# Compute heat threshold. 
heat_thresh <- merged_data[HR %in% hr_sel & YEAR %in% 2000:2024, .(Q = quantile(get(paste0(tvar, "0")), 0.975)), by = "HR"]

# Check that matrices fits with results data.
sum(sapply(tvar_mat_l, nrow)) == nrow(data_final)
sum(sapply(relh_mat_l, nrow)) == nrow(data_final)


# ------------------------------------------------------------------------------
# Step 1 DLNM model fitting by region (HR) -------------------------------------
# ------------------------------------------------------------------------------


# Model specification for temperature variable and control (relh).
arglag_tvar <- list(fun = "ns", knots = logknots(max_lag, nk = nk_xlag))
argvar_ctrl <- list(fun = "lin")
arglag_ctrl <- arglag_tvar  

# Empty elements to store results.
coef_list  <- as.list(rep(NA, length.out = length(hr_sel)))
vcov_list  <- as.list(rep(NA, length.out = length(hr_sel)))
cr_list    <- as.list(rep(NA, length.out = length(hr_sel)))

# Fit DLNM by region (HR).
for (hr in hr_sel){ 

    # Message.
    if (verbose) message("Fitting DLNM on HR ", hr, ".")
    
    # Filter data based on <HR>.
    data_hr <- data_final[HR == hr, ]
    nyears <- length(unique(data_hr$YEAR))
    
    # Extract indicator <hr_i>.
    hr_i <- which(hr_sel == hr)
    
    # Extract temperature and control variables.
    tvar_mat <- tvar_mat_l[[hr_i]]
    relh_mat <- relh_mat_l[[hr_i]]
    
    # Argument for the tvar cross-basis.
    argvar_tvar <- list(
        fun   = "ns", 
        knots = quantile(data_hr[[tvar]], probs = knots_xvar), 
        Bound = range(data_hr[[tvar]])  
    ) 
    
    # Create cross-basis.
    cb_tvar <- crossbasis(tvar_mat, lag = max_lag, argvar = argvar_tvar, arglag = arglag_tvar)
    if (relh) cb_relh <- crossbasis(relh_mat, lag = max_lag, argvar = argvar_ctrl, arglag = arglag_ctrl) 
    
    # Create formula.
    if (relh) form <- paste0("COUNT ~ cb_tvar + cb_relh + WEEKDAY_F + HOL + YEAR_F * ns(DOY, df = df_date)")
    if (!relh) form <- paste0("COUNT ~ cb_tvar + WEEKDAY_F + HOL + YEAR_F * ns(DOY, df = df_date)")
    
    # Fitting DLNM with a quasiPoisson dsitribution.
    dlnm_HR <- glm(
        formula = as.formula(form), 
        family  = quasipoisson(), 
        data    = data_hr
    ) 
    
    # Reduce DLNM effect using mean values as MMT.
    cb_tvar_reduce <- crossreduce(
        basis = cb_tvar, 
        model = dlnm_HR, 
        cen   = mean(data_hr[[tvar]]),
        at    = quantile(data_hr[[tvar]], probs = q_mmt)
    )
    
    # Extract MMT and extreme temperature threshold.
    mmt <- quantile(data_hr[[tvar]], probs = q_mmt[which.min(cb_tvar_reduce$RRfit)])
    ext <- max(mmt, heat_thresh[HR == hr, Q])
    
    # Reduced function at MMT for plotting.
    cb_tvar_reduce <- crossreduce(
        basis = cb_tvar, 
        model = dlnm_HR, 
        cen   = mmt,
        at    = seq(min(data_hr[[tvar]]), max(data_hr[[tvar]]), by = 0.1)
    )
    
    # Save information for future steps.
    cr_list[[hr_i]]    <- cb_tvar_reduce
    coef_list[[hr_i]]  <- coef(cb_tvar_reduce) 
    vcov_list[[hr_i]]  <- vcov(cb_tvar_reduce)

}

# Note : 'cr_list' countains the reduced function for each HR (i.e., regional DLNM).
# For example : 
# plot(cr_list[[6L]], xlab = tvar, exp = TRUE, y = "RR")

# Convert saved coefficient into a matrix.
coef_matrix <- matrix(
    data     = unlist(coef_list),
    nrow     = length(coef_list), 
    ncol     = length(coef_list[[1L]]),
    byrow    = TRUE, 
    dimnames = list(hr_sel)
)


# ------------------------------------------------------------------------------
# Step 2 : Regional BLUPs from meta-regression ---------------------------------
# ------------------------------------------------------------------------------


# Fit meta-regression with stepwise variable selection.
meta_reg <- fit_stepwise_meta_aic(
    coef       = coef_matrix, 
    vcov       = vcov_list,
    meta_preds = names(meta_data)[-1]
)
meta_reg_summ <- summary(meta_reg)

# I2 statistic.
meta_reg_summ$i2stat[1L]

# P-value for the Cochran Q-test.
meta_reg_summ$qstat$pvalue[1L]

# Formula of the meta-regression
meta_reg$formula

# Extract BLUPs.
blup_meta <- blup(meta_reg, vcov = TRUE) 

# Appending results to the following vectors.
qmin_hr <- rep(NA, length(hr_sel))
mmt_hr <- rep(NA, length(hr_sel))

# Extract MMT based on BLUP for each region.
for (hr in hr_sel){ 
    
    # Filter data based on <HR>.
    data_hr <- data_final[HR == hr, ]
    hr_i <- which(hr_sel == hr)
    
    # Extract tvar values to extract MMT based on q_mmt.
    predvar <- quantile(data_hr[[tvar]], probs = q_mmt) 
    
    # Argument for the tvar cross-basis.
    argvar_tvar <- list(
        x     = predvar,            
        fun   = "ns", 
        knots = quantile(data_hr[[tvar]], probs = knots_xvar), 
        Bound = range(data_hr[[tvar]]) 
    ) 
    
    # Apply one basis to the tvar.
    bvar <- do.call(onebasis, argvar_tvar) 
    
    # Extract quantile of MMT and MMT. 
    qmin_hr[hr_i] <- q_mmt[which.min((bvar %*% blup_meta[[hr_i]]$blup))]
    mmt_hr[hr_i] <- quantile(data_hr[[tvar]], qmin_hr[hr_i]) 
    
}

# Extract BLUPs of reduced functions at MMT for each HR.
red_blups_mmt <- lapply(hr_sel, function(hr) {
    
    # Filter data based on <HR>.
    data_hr <- data_final[HR == hr, ]
    hr_i <- which(hr_sel == hr)
    
    # Argument for the tvar cross-basis.
    argvar_tvar <- list(
        fun   = "ns", 
        knots = quantile(data_hr[[tvar]], probs = knots_xvar), 
        Bound = range(data_hr[[tvar]])
    ) 
    
    # Extract coef and vcov from the BLUPs of meta-regression.
    coef <- blup_meta[[hr_i]]$blup
    vcov <- blup_meta[[hr_i]]$vcov
    
    # Extract MMT.
    mmt <- mmt_hr[hr_i]
    
    # Centered basis functions.
    bvar <- do.call(onebasis, c(list(x = data_hr[[tvar]]), argvar_tvar))
    cenvec <- do.call(onebasis, c(list(x = mmt), argvar_tvar))
    bvarcen <- scale(bvar, center = cenvec, scale = F)
    
    # Predict function values.
    red_blup <- crosspred(bvarcen, coef = coef, vcov = vcov, cen = mmt)
    
})

# Set names.
names(red_blups_mmt) <- hr_sel

# Note : 'red_blups_mmt' contains the BLUPs estimates for each HR at MMT.
# plot(red_blups_mmt[[6L]], xlab = tvar, exp = TRUE, ylab = "RR")

# Exports red_blups_mmt.
fwrite(do.call(rbind, lapply(1:length(hr_sel), function(hr_i) {
    data.table(
        HR     = hr_sel[[hr_i]],
        MODEL  = "Summer",
        ID     = tbl_i, 
        NAME   = tbl[tbl_i, NAME],
        X      = red_blups_mmt[[hr_i]]$predvar,
        Y      = exp(red_blups_mmt[[hr_i]]$allfit),
        Y_LOW  = exp(red_blups_mmt[[hr_i]]$alllow),
        Y_HIGH = exp(red_blups_mmt[[hr_i]]$allhigh)
    )
})), paste0("out/", fname, "blups.csv"))


# ------------------------------------------------------------------------------
# Step 3 : Burden quantification with reduced regional BLUPs -------------------
# ------------------------------------------------------------------------------


# Get years.
years <- unique(data_final$YEAR)
nyears <- length(years)

# Compute simulated AN for each HR based on BLUPs.
data_sim <- do.call(rbind, lapply(hr_sel, function(hr) {
    
    # Message.
    if (verbose) message("Computing AN for region ", hr, ".")
    
    # Filter data based on <HR>.
    data_hr <- data_final[HR == hr, ]
    hr_i <- which(hr_sel == hr)
    
    # Argument for the tvar cross-basis.
    argvar_tvar <- list(
        fun   = "ns", 
        knots = quantile(data_hr[[tvar]], probs = knots_xvar), 
        Bound = range(data_hr[[tvar]]) 
    ) 
    
    # Extract MMT and threshold for extreme temperature
    mmt <- mmt_hr[hr_i]
    ext <- heat_thresh[HR == hr, Q]
    
    # Extract coef and vcov from the meta-regression.
    coef <- blup_meta[[hr_i]]$blup
    vcov <- blup_meta[[hr_i]]$vcov
    
    # Centered basis functions.
    bvar <- do.call(onebasis, c(list(x = data_hr[[tvar]]), argvar_tvar))
    cenvec <- do.call(onebasis, c(list(x = mmt), argvar_tvar))
    bvarcen <- scale(bvar, center = cenvec, scale = F)
    
    # Indicators for heat and extreme heat days.
    data_hr[, heat_all     := get(tvar) > mmt]
    data_hr[, heat_mod     := get(tvar) > mmt & get(tvar) <= max(mmt, ext)]
    data_hr[, heat_extreme := get(tvar) > max(mmt, ext)]
    
    # Compute the contribution of daily values (average value).
    an_mean <- (1 - exp(-bvarcen %*% coef)) * data_hr[["COUNT"]]
    
    # Mean AN values for heat and extreme heat.
    an_heatall <- sum(an_mean[data_hr$heat_all])/nyears
    an_heatmod <- sum(an_mean[data_hr$heat_mod])/nyears
    an_heatext <- sum(an_mean[data_hr$heat_extreme])/nyears
    
    # Sampling coefficient from BLUP assuming mvnorm.
    set.seed(2912L)
    coefsim <- MASS::mvrnorm(nsim, coef, vcov)
    
    # Compute simulations.
    an_sim <- sapply(1:nsim, function(s) {
        (1 - exp(-bvarcen %*% coefsim[s, ])) * data_hr[["COUNT"]]
    })
    
    # Export data.
    colnames(an_sim) <- cols_sim
    data_sim <- cbind(data_hr[, .(HR, DATE, YEAR, heat_mod, heat_extreme, heat_all)], an_sim)
    
    # Return data_sim.
    return(data_sim)
    
}))

# Aggregate simulations by years.
x_year <- do.call(rbind, lapply(c("heat_mod", "heat_extreme", "heat_all"), function(trange) {
    
    # Compute mean, lower and higher value.
    do.call(rbind, lapply(years, function(year) {
        
        # Extract values of interest.
        val <- ul(data_sim[YEAR == year & get(trange) == TRUE, lapply(.SD, sum), .SDcols = cols_sim])
        
        # Extract mean values and quantiles.
        return(data.table(
            HR     = 99L,
            YEAR    = year,
            TRANGE  = trange,
            AN      = mean(val),
            AN_LOW  = quantile(val, probs = 0.025),
            AN_HIGH = quantile(val, probs = 0.975)
        ))
        
    }))
    
}))

# Aggregate simulations by region (HR).
x_hr <- do.call(rbind, lapply(c("heat_mod", "heat_extreme", "heat_all"), function(trange) {
    
    # Compute mean, lower and higher value.
    do.call(rbind, lapply(hr_sel, function(hr) {
        
        # Extract values of interest.
        val <- ul(data_sim[HR == hr & get(trange) == TRUE, lapply(.SD, sum), .SDcols = cols_sim])
        
        # Extract mean values and quantiles.
        return(data.table(
            HR      = hr,
            YEAR    = 9999L,
            TRANGE  = trange,
            AN      = mean(val)/nyears,
            AN_LOW  = quantile(val, probs = 0.025)/nyears,
            AN_HIGH = quantile(val, probs = 0.975)/nyears
        ))
        
    }))
    
}))

# Aggregate simulations for the whole province.
x_qc <- do.call(rbind, lapply(c("heat_mod", "heat_extreme", "heat_all"), function(trange) {
    
    # Extract values of interest.
    val <- ul(data_sim[get(trange) == TRUE, lapply(.SD, sum), .SDcols = cols_sim])
    
    # Extract mean values and quantiles.
    return(data.table(
        HR     = 99L,
        YEAR    = 9999L,
        TRANGE  = trange,
        AN      = mean(val)/nyears,
        AN_LOW  = quantile(val, probs = 0.025)/nyears,
        AN_HIGH = quantile(val, probs = 0.975)/nyears
    ))
    
}))

# Merge all tables.
an_tbl <- rbind(x_qc, x_hr, x_year)

# Count number per year and region (HR).
count_year_hr <- rbind(
    data_final[, .(COUNT = sum(COUNT)/length(unique(YEAR)), YEAR = 9999L), by = c("HR")],
    data_final[, .(COUNT = sum(COUNT), HR = 99L), by = c("YEAR")],
    data_final[, .(COUNT = sum(COUNT)/length(unique(YEAR)), YEAR = 9999L, HR = 99L)]
)

# Compute AF results.
an_tbl <- merge(an_tbl, count_year_hr, by = c("HR", "YEAR"), all.x = TRUE)
an_tbl[, `:=`(AF = AN / COUNT, AF_LOW = AN_LOW / COUNT, AF_HIGH = AN_HIGH / COUNT)]
an_tbl[, COUNT := NULL]

# Note : 'an_tbl' contains AN/AF by HR/Year based on regional BLUPS.
an_tbl[HR == 99 & YEAR == 9999, ]

# Export an_tbl.
fwrite(an_tbl[, `:=`(ID = tbl_i, NAME = tbl[tbl_i, NAME])][], paste0("out/", fname, "an_tbl.csv"))


# ------------------------------------------------------------------------------
# Step 4 : Pooled cumulative functions at MMT ----------------------------------
# ------------------------------------------------------------------------------


# Extract full range of temperature.
range_temp <- range(data_final[[tvar]])

# Extract tested values of temperature within Q25-Q98.
predvar <- quantile(data_final[[tvar]], probs = q_mmt) 
argvar_tvar <- list(
    fun    = "ns", 
    knots  = quantile(data_final[[tvar]], probs = knots_xvar),
    Bound  = range_temp
) 
bvar <- do.call("onebasis", c(list(x = predvar), argvar_tvar))

# Extract coef and vcov from meta-regression.
coef_meta <- coef(meta_reg)[1:ncol(coef_matrix)]
vcov_meta <- vcov(meta_reg)[1:ncol(coef_matrix), 1:ncol(coef_matrix)]

# Get pooled effect (without MMT).
cpall_prelim <- crosspred(
    basis      = bvar,
    coef       = coef_meta,
    vcov       = vcov_meta,
    model.link = "log",
    cen        = mean(mmt_hr),
    by         = 0.1
) 

# Get MMT.
mmt_pooled <- cpall_prelim$predvar[which.min(cpall_prelim$matRRfit)]

# Refit pooled effect at all temperature values.
predvar <- seq(range_temp[1L], range_temp[2L], by = 0.1)
argvar_tvar <- list(
    fun    = "ns", 
    knots  = quantile(data_final[[tvar]], probs = knots_xvar),
    Bound  = range_temp
) 
bvar <- do.call("onebasis", c(list(x = predvar), argvar_tvar))

# Get pooled effect (with MMT).
cpall_final <- crosspred(
    basis      = bvar,
    coef       = coef(meta_reg)[1:ncol(coef_matrix)],
    vcov       = vcov(meta_reg)[1:ncol(coef_matrix), 1:ncol(coef_matrix)],
    model.link = "log",
    by         = 0.1,
    from       = range_temp[1],
    to         = range_temp[2], 
    cen        = mmt_pooled
) 

# Note : 'cpall_final' is the pooled (across Quebec) cumulative function at MMT.
# plot(cpall_final, xlab = tvar, exp = TRUE, ylab = "RR")

# Export cpall_final
fwrite(data.table(
    ID     = tbl_i, 
    NAME   = tbl[tbl_i, NAME],
    X      = cpall_final$predvar,
    Y      = cpall_final$allRRfit,
    Y_LOW  = cpall_final$allRRlow,
    Y_HIGH = cpall_final$allRRhigh
), paste0("out/", fname, "cpall_final.csv"))

# Export MMT values.
fwrite(data.table(
    ID     = tbl_i, 
    NAME   = tbl[tbl_i, NAME],
    HR     = c(hr_sel, 99L),
    MMT    = c(mmt_hr, mmt_pooled)
), paste0("out/", fname, "mmt.csv"))

# Export results of meta-regression.
fwrite(data.table(
    ID     = tbl_i, 
    NAME   = tbl[tbl_i, NAME],
    I2     = meta_reg_summ$i2stat[1L],
    Q      = meta_reg_summ$qstat$pvalue[1L],
    FORM   = paste0(as.character(meta_reg$formula), collapse = "_")
), paste0("out/", fname, "meta.csv"))


}
