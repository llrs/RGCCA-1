shorten.label = function(x){
  if (stringr::str_length(x) >25){
    return(stringr::str_c(str_sub(x,1,25),"..."))
  } else {
    return(x)
  }
}