Csurvsplit2 <- function(){
  newfunc <- eval(parse(text = paste0(
    "Csurvsplit", '<-survival:::', "Csurvsplit"
  )))
   newfunc
}
