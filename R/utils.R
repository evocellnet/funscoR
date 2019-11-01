log_scaling <- function(d) {
    l <- log10(d)
    return((l - min(l)) / (max(l) - min(l)))
}
