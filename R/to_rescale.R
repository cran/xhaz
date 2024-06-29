to_rescale2 <- function(){
  newfunc <- eval(parse(text = paste0(
    "mexhazAlpha", '<-mexhaz:::', "mexhazAlpha"
  )))
  newfunc
}

