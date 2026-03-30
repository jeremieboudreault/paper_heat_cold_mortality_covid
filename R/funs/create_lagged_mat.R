# create_lagged_mat.R

create_lagged_mat <- function(var) {
    
    # Loop on all RSS.
    mat_l <- lapply(hr_sel, function(hr_i) {
        
        return(as.matrix(merged_data[
            i    = (HR == hr_i & MONTH %in% months & YEAR %in% years),
            j    = paste0(var, 0:max_lag), 
            with = FALSE
        ]))
        
    })
    
    # Add names.
    names(mat_l) <- hr_sel
    
    # Return.
    return(mat_l)
    
}
