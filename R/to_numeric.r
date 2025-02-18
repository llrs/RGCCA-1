to_numeric <- function(df) {
    df <- as.matrix(df) # deals with data.frame case
    matrix(
        sapply(
            seq(NROW(df) * NCOL(df)),
            function(i)
                tryCatch(
                    as.numeric(df[i]),
                    warning = function(e) NA
                )),
        NROW(df),
        NCOL(df),
        dimnames = list(row.names(df), colnames(df))
    )
}
