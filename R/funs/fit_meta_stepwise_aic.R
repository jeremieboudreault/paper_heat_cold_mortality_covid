# fit_stepwise_meta_aic.R

fit_stepwise_meta_aic <- function(coef, vcov, meta_preds) {
    
    # Function to AIC from meta-regression.
    .get_meta_aic <- function(meta) summary(meta)$AIC
    
    # Initial parameters.
    meta_pred_avail <- meta_preds
    form_init <- "coef ~ 1L"
    meta0 <- mixmeta(coef ~ 1L, S = vcov, method = "reml")
    aic0 <- .get_meta_aic(meta0)
    loop <- TRUE
    
    # Loop to get best model.
    while(loop) {
        
        # Fit meta-reg with all available predictors.
        aic_preds <- sapply(meta_pred_avail, function(pred) {
            form <- as.formula(paste0(form_init, " + meta_data$", pred))
            meta <- mixmeta(form, S = vcov, method = "reml")
            return(.get_meta_aic(meta))
        })
        
        # Which meta-predictors is select.
        delta_aic <-  aic0 - aic_preds
        pred_to_add <- delta_aic[max(delta_aic) == delta_aic]
        
        # Stop criteria.
        if (length(pred_to_add) == 0 | pred_to_add < 0)
            loop = FALSE
        else {
            form_init <- paste0(form_init, " + meta_data$", names(pred_to_add))
            aic0 <- aic_preds[names(aic_preds) == names(pred_to_add)]
            meta_pred_avail <- meta_pred_avail[-which(names(pred_to_add) == meta_pred_avail)]
        }
        
    }
    
    # Refit final model.
    meta_final <-  mixmeta(as.formula(form_init), S = vcov, method = "reml")
    
    # Return.
    return(meta_final)
    
}