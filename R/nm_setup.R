#' Set up the NONMEM optiProject
#' 
#' Function to generate data and model from an optiProject-object
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object with a model, penalizations, dose levels, and parameters.
#' For requirements regarding penalizations, dose levels, and parameters.
#' @param name Name under which optiProject will be saved. .ctl will be added to name.
#' @param save_path Path, where optiProject should be saved.
#' @param pop Whether multiple subjects should be used (TRUE), either from added individual parameters or
#' simulated from population parameters, or if only one subject with population parameters should be used
#' (FALSE).
#' @param iiv iiv Whether doses should include inter-individual variability (TRUE) or not (FALSE).
#' @param sim If \code{pop=TRUE}, whether individuals should be simulated (TRUE) or added individual
#' parameters should be used (FALSE).
#' @param n_id If \code{pop=TRUE} and \code{sim=TRUE}, number of individuals to be simulated.
#' @param seperate If \code{pop=TRUE} and \code{iiv=FALSE}, whether individuals should be fitted separately (TRUE) or
#' all individuals should be fitted together.
#' @importFrom magrittr %>%
#' @author Dominic Bräm
nmSetAll <- function(optiProject,
                     name,
                     save_path,
                     pop,
                     iiv,
                     sim,
                     n_id,
                     separate){
  
  data <- DataGen(optiProject,pop = pop, sim = sim, n_id = n_id)
  data_name <- paste0(save_path,"/",name,"_data.csv")
  write.csv(data,data_name,row.names = F,quote = F)
  print(paste0("Data has been saved under: ",data_name))
  
  optiProject <- optiProject %>%
    nmUpdateThetaEta() %>%
    nmUpdateInputs() %>%
    nmUpdateData(name) %>%
    nmUpdateSubrout() %>%
    nmUpdateNcomp() %>%
    nmUpdatePK(iiv) %>%
    nmUpdateDES() %>%
    nmUpdateError(iiv) %>%
    nmUpdateTheta() %>%
    nmUpdateOmega(iiv) %>%
    nmUpdateSigma(iiv) %>%
    nmUpdateEst(iiv) %>%
    nmUpdateTable(name,iiv)
  
  if(pop & !iiv & separate){
    optiProject <- optiProject %>%
      nmUpdateMulti(name,sim,n_id)
    nm_multi_sep_code(optiProject,name,save_path)
  }
  
  model_name <- paste0(save_path,"/",name,"_model.ctl")
  writeLines(optiProject$Model,con=model_name)
  print(paste0("Model has been saved under: ",model_name))
}