
longDF2ratetable <-
  function(DF,
           value.var = "haz",
           by.vars = setdiff(names(DF), value.var)) {
    univals <- lapply(DF[, by.vars], unique)
    names(univals) <- NULL
    dimvec <- sapply(DF[, by.vars], function(x) {
      length(unique(x))
    },
    simplify = TRUE)
    ar <- array(DF[, value.var], dim = dimvec)
    dimnames(ar) <- univals
    attr(ar, "class") <- "ratetable"
    attr(ar, "dimid") <- colnames(DF)
    ar
    invisible()
  }
