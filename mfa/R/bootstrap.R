#' bootstrap
#'
#' Perform a boostrap approximation for the compromise of the
#' given mfa object
#' @param mfa an mfa object
#' @param B the number of bootstrap samples to run
#' @param seed a seed for the random number generator
#' @return a copy of the mfa object with the bootstrap as an
#'         attribute.
#' @export
bootstrap <- function(mfa, B = 100, seed = 0){

    ## run the necessary checks
    check_mfa(mfa)
    check_B(B)
    check_seed(seed)
    
    ## set the seed
    set.seed(seed)
    
    ## extract all the necessary info from the mfa object
    data <- mfa$commonFactorScores %*% t(mfa$Q)
    sets <- attributes(mfa)$sets
    nSets <- length(sets)
    ids <- attributes(mfa)$ids
    ncomps <- attributes(mfa)$ncomps
    M <- diag(attributes(mfa)$rowWeights)
    A <- diag(attributes(mfa)$colWeights)

    ## these are important for keeping track of
    ## each bootstrap sample, so make sure they're useful
    ## (not null)
    if(is.null(ids)){
        ids <- 1:nrow(data)
    }
    
    ## make a temporary boot strap function (just a segment
    ## of the mfa function)
    boot_fun <- function(bootSets){
        Xk <- lapply(bootSets, function(g){
            as.matrix(scale(data[, g]))
        })

        X <- do.call(cbind, Xk)
        S <- X %*% A %*% t(X)

        ## perform GSVD
        SVD <- svd(S)
        sqrt(solve(M)) %*% SVD$u
    }

    ## run the bootstrap
    pBoot <- lapply(1:B, function(x) NA)
    b <- 0
    while(b < B){
        i <- sample(1:nSets, nSets, replace = TRUE)
        ## some samples will lead to singular decompositions
        ## so wrap each attempt in a 'try' function and only
        ## keep the successful ones.
        boot <- try(boot_fun(sets[i]), silent = TRUE)
        if(!identical(class(boot), "try-error")){
            b <- b + 1
            boot <- as.data.frame(boot[,1:ncomps])
            boot$id <- ids
            pBoot[[b]] <- boot
        }
    }
    
    attr(mfa, "boot") <- do.call(rbind, pBoot)
    mfa
}
