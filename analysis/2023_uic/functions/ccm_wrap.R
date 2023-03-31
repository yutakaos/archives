#--------------------------------------------------------------------------------#
# Wrapper function for CCM
#--------------------------------------------------------------------------------#

ccm_wrap = function (
    block, lib = c(1, NROW(block)), pred = lib, lib_var = 1, tar_var = 2,
    E = 1, tau = 1, tp = 0)
{
    require(rEDM) # v0.7.5
    op_ccm <- rEDM::ccm(
        block, lib=lib, pred=pred, norm=2, E=E, tau=tau, tp=tp,
        num_neighbors="e+1", lib_sizes=c(E+1, nrow(block)), num_samples=100,
        lib_column=lib_var, target_column=tar_var, silent=TRUE) 
    op <- rEDM::ccm_means(op_ccm, FUN=mean, na.rm=TRUE)
    te <- with(op, log(rmse[1])-log(rmse[2]))
    op <- data.frame(op[2,c("E","tau","tp","nn","num_pred","rmse")], te=te)
    rownames(op) <- NULL
    op
}

# End