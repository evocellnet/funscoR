#' Function to log scale resulting functional scores
#'
#' @param d vector of scores
#'
#' @return log-scaled scores
#' @export
log_scaling <- function(d) {
    l <- log10(d)
    return((l - min(l)) / (max(l) - min(l)))
}
