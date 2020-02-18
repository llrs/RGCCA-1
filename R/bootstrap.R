#' Compute bootstrap
#'
#' Bootstrapping individuals for obtaining confidence intervals for weights of RGCCCA (draw with replacement).
#' Requires then the use of get_bootstrap, to extract the component and block of interest.
#' @inheritParams rgcca
#' @inheritParams plot_var_2D
#' @param n_boot A integer for the number of boostrap
#' @param n_cores An integer for the number of cores used in parallelization 
#' @param ... other RGCCA parameters # TODO
#' @return A list of RGCCA bootstrap weights
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' bootstrap(rgcca_out, n_boot = 2, n_cores = 1)
#' bootstrap(rgcca_out, n_boot = 2, n_cores = 1, blocks = lapply(blocks, scale),
#'  superblock = FALSE)
#' @export
#' @seealso get_bootstrap, plot_bootstrap_1D
bootstrap <- function(
    rgcca_res,
    n_boot = 5,
    n_cores = parallel::detectCores() - 1,
    ...) {

    # TODO : nboot > 1
    stopifnot(!missing(rgcca_res))

    if (n_cores == 0)
        n_cores <- 1

    # if (any(unlist(lapply(rgcca_res$call$blocks, NCOL) > 1000)))
    #     verbose <- TRUE

    cat("Bootstrap in progress...")

    W <- parallel::mclapply(
        seq(n_boot), 
        function(x) bootstrap_k(rgcca_res, ...), 
        mc.cores = n_cores)

    cat("OK.\n", append = TRUE)

    return(list(bootstrap = W, rgcca = rgcca_res))
}
