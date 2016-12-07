#' RV Table function
#'
#' Function to calculate the RV coefficient between sets of variables
#' @param mfa the mfa object
#' @param sets list of vectors containing the variables between which the RV
#'             coefficient
#' will be calculated
#' @keywords RV similarity measure table set
#' @examples
#' data(wine)                                    
#' i <- grep("V", colnames(wine))                
#' data <- wine[,i]                              
#' sets <- list(1:6, 7:12, 13:18, 19:23, 24:29,  
#'             30:34, 35:38, 39:44, 45:49, 50:53)
#' MFA <- mfa(data, sets)
#' RV_table(MFA, sets = list(1:6, 19:23, 24:29))
#' @export
RV_table <- function(mfa, sets = NULL){
    check_mfa(mfa)
    if(is.null(sets)){
        sets <- attributes(mfa)$sets
    }
    check_sets(sets)
    dataset <- mfa$commonFactorScores %*% t(mfa$Q)
    rvMat <- matrix(rep(0, length(sets)^2), length(sets))
    for(i in 1:length(sets)){
        for(j in i:length(sets)){
            rvMat[i,j] <- RV(dataset[, sets[[i]]], dataset[, sets[[j]]])
            if(i != j){
                rvMat[j, i] <- rvMat[i, j]
            }
        }
    }
    rvMat
}
