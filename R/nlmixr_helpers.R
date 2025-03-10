#' Identify residual error parameters
#' 
#' Identify the error parameters in prop(X) and add(X) in the nlmixr2 model.
#' 
#' NULL
#' 
#' @param fit_obj Fitted object with nlmixr2.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nlmixr_err_identifier <- function(fit_obj){
  model_name <- fit_obj$modelName
  
  if(!exists(model_name)){
    stop("The model function in the nlmixr2 fit must be R global environment.")
  }
  
  model_text <- deparse(get(model_name))
  
  prop_err_parms <- regexpr("prop\\([^\\)]+\\)",model_text) %>%
    regmatches(model_text,.) %>%
    {gsub("prop\\(|\\(|\\)","",.)} %>%
    unlist()
  add_err_parms <- regexpr("add\\([^\\)]+\\)",model_text) %>%
    regmatches(model_text,.) %>%
    {gsub("add\\(|\\(|\\)","",.)} %>%
    unlist()
  
  err_parms <- c(prop_err_parms,
                 add_err_parms)
  
  return(err_parms)
}
