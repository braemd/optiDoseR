#' Set up the nlmixr2 optiProject
#' 
#' Function to generate data and model from an optiProject-object
#' 
#' Data and model will be saved in the global environment, \emph{_data} will be added to the name 
#' for the data and \emph{_model} to the name for the model.
#' 
#' @param optiProject A optiProject-object with a model, penalizations, dose levels, and parameters.
#' For requirements regarding penalizations, dose levels, and parameters.
#' @param name Name under which optiProject will be saved.
#' @param pop Whether multiple subjects should be used (TRUE), either from added individual parameters or
#' simulated from population parameters, or if only one subject with population parameters should be used
#' (FALSE).
#' @param iiv iiv Whether doses should include inter-individual variability (TRUE) or not (FALSE).
#' @param sim If \code{pop=TRUE}, whether individuals should be simulated (TRUE) or added individual
#' parameters should be used (FALSE).
#' @param n_id If \code{pop=TRUE} and \code{sim=TRUE}, number of individuals to be simulated.
#' @param const_err If doses are estimated with inter-individual variability for a population of individual
#' subjects, \code{const_err} defines the balance between fulfilling the objectives (constraints, target values etc.)
#' and having the doses log-normal distributed.
#' @param run_model Whether the generated nlmixr2 optiProject should directly be run.
#' @param control If \code{run_model=TRUE}, list of control arguments given to nlmxir2. If \code{iiv=FALSE},
#' control for "bobyqa" and if \code{iiv=TRUE} control for "saem" must be given.
#' @importFrom magrittr %>%
#' @author Dominic Bräm
nlmixrSetAll <- function(optiProject,
                     name,
                     pop,
                     iiv,
                     sim,
                     n_id,
                     const_err,
                     run_model,
                     control){
  data <- DataGen(optiProject,pop = pop, sim = sim, n_id = n_id) %>%
    mutate(AMT = as.numeric(ifelse(AMT==".",0,AMT)))
  data_name <- paste0(name,"_data")
  assign(data_name,data,pos=".GlobalEnv")
  print(paste0("Data has been saved in the global environment under: ",data_name))
  
  optiProject <- optiProject %>%
    nlmixrUpdateIni(pop,iiv,const_err) %>%
    nlmixrUpdateModelInis() %>% 
    nlmixrUpdateModelParms(iiv) %>% 
    nlmixrUpdateModelBioav() %>%
    nlmixrUpdateModelPens() %>%
    nlmixrUpdateModelOut()
  
  model_name <- paste0(name,"_model")
  model_fun <- paste(optiProject$Model,collapse="\n") %>%
    {eval(str2lang(.))}
  assign(model_name,model_fun,pos=".GlobalEnv")
  print(paste0("Model has been saved in the global environment under: ",model_name))
  
  if(run_model){
    if(!pop){
      fit <- nlmixr2::nlmixr2(object = get(model_name,pos=".GlobalEnv"),
                              data = get(data_name,pos=".GlobalEnv"),
                              est = "bobyqa",
                              control = control)
    } else{
      fit <- nlmixr2::nlmixr2(object = get(model_name,pos=".GlobalEnv"),
                              data = get(data_name,pos=".GlobalEnv"),
                              est = "saem",
                              control = control)
    }
    fit_name <- paste0(name,"_fit")
    assign(fit_name,fit,pos=".GlobalEnv")
    print(paste0("Fit has been saved in the global environment under: ",fit_name))
  }
}
