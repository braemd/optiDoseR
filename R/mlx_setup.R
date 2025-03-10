#' Set up the Monolix optiProject
#' 
#' Function to generate data and model from an optiProject-object
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object with a model, penalizations, dose levels, and parameters.
#' For requirements regarding penalizations, dose levels, and parameters.
#' @param name Name under which optiProject will be saved. .mlxtran will be added to name.
#' @param save_path Path, where optiProject should be saved.
#' @param pop Whether multiple subjects should be used (TRUE), either from added individual parameters or
#' simulated from population parameters, or if only one subject with population parameters should be used
#' (FALSE).
#' @param iiv iiv Whether doses should include inter-individual variability (TRUE) or not (FALSE).
#' @param sim If \code{pop=TRUE}, whether individuals should be simulated (TRUE) or added individual
#' parameters should be used (FALSE).
#' @param n_id If \code{pop=TRUE} and \code{sim=TRUE}, number of individuals to be simulated.
#' @param parametrisation If a Monolix-optiProject was created and a population of individuals is simulated
#' (\code{pop=TRUE} and \code{sim=TRUE}), the parametrisation of the parameters in the Monolix model
#' must be given. Must be in "log","normal",or "logit". If all parameters have same parametrisation,
#' a single string can be given. If different parametrisations for different parameters are required, a vector
#' with parametrisation for each parameter is required. Default is "log" for all parameters.
#' @param lower_parm_limits If one or multiple parameters in the model are "logit" distributed and individuals
#' are simulated, the lower limit(s) for the "logit"-distributed parameter(s).
#' @param upper_parm_limits If one or multiple parameters in the model are "logit" distributed and individuals
#' are simulated, the upper limit(s) for the "logit"-distributed parameter(s).
#' @param const_err If doses are estimated with inter-individual variability for a population of individual
#' subjects, \code{const_err} defines the balance between fulfilling the objectives (constraints, target values etc.)
#' and having the doses log-normal distributed.
#' @param create_mlx If optiProject-software is "Monolix", whether directly a Monolix project (.mlxtran file) should
#' be created.
#' @importFrom magrittr %>%
#' @author Dominic Bräm
mlxSetAll <- function(optiProject,
                      name,
                      save_path,
                      pop,
                      iiv,
                      sim,
                      n_id,
                      parametrisation,
                      lower_parm_limits,
                      upper_parm_limits,
                      const_err,
                      create_mlx){
  data <- DataGen(optiProject,pop = pop, sim = sim, n_id = n_id, parametrisation = parametrisation,
                  upper_parm_limits = upper_parm_limits, lower_parm_limits = lower_parm_limits)
  data_name <- paste0(save_path,"/",name,"_data.csv")
  write.csv(data,data_name,row.names = F)
  print(paste0("Data has been saved under: ",data_name))
  
  optiProject <- optiProject %>%
    mlxUpdateParmNames() %>%
    mlxUpdateInputs() %>%
    mlxUpdatePK() %>%
    mlxUpdateBioav() %>%
    mlxUpdateEquation() %>%
    mlxUpdateOutput()
  
  model_name <- paste0(save_path,"/",name,"_model.txt")
  writeLines(optiProject$Model,con=model_name)
  print(paste0("Model has been saved under: ",model_name))
  
  if(create_mlx){
    header_types <- mlx_header_types(data)
    mlxtran_name <- paste0(save_path,"/",name,"_mlx.mlxtran")
    mlxSetMlxtran(optiProject,
                  mlxtran_name,
                  model_name,
                  data_name,
                  header_types,
                  iiv,
                  const_err)
  }
}

#' Set up the .mlxtran file
#' 
#' Set up the .mlxtran file from an optiProject-object
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object with a model, penalizations, dose levels, and parameters.
#' For requirements regarding penalizations, dose levels, and parameters.
#' @param mlxtran_name Name under which the Monolix file should be saved. Must include .mlxtran
#' @param model_name Name of Monolix structural model file to be used.
#' @param data_name Name of data file to be used.
#' @param header_types Header types of columns in the data file.
#' @param iiv Whether doses should include inter-individual variability (TRUE) or not (FALSE).
#' @param const_err If doses are estimated with inter-individual variability for a population of individual
#' subjects, \code{const_err} defines the balance between fulfilling the objectives (constraints, target values etc.)
#' and having the doses log-normal distributed.
#' @author Dominic Bräm
mlxSetMlxtran <- function(optiProject,
                          mlxtran_name,
                          model_name,
                          data_name,
                          header_types,
                          iiv=FALSE,
                          const_err=0.0001){
  
  newProject(modelFile = model_name, data = list(dataFile=data_name, headerTypes=header_types))
  
  parms <- paste0("F",1:length(optiProject$Doses))
  
  #Variability
  if(!iiv){
    parms_var <- as.list(rep(FALSE,length(parms))) %>%
      `names<-`(parms)
  } else{
    parms_var <- as.list(rep(TRUE,length(parms))) %>%
      `names<-`(parms)
  }
  setIndividualParameterVariability(parms_var)
  
  #Population initial values
  parm_inis <- lapply(optiProject$Doses,function(x) x$Ini) %>%
    unlist()
  if(is.null(parm_inis)){
    parm_inis <- 1
  }
  
  mlx_inis <- data.frame(name = paste0(parms,"_pop"),
                               initialValue = parm_inis,
                               method = "MLE")
  
  #Omega initial values
  if(iiv){
    omega_inis <- data.frame(name = paste0("omega_",parms),
                             initialValue = 1,
                             method = "MLE")
    
    mlx_inis <- rbind(mlx_inis,
                      omega_inis)
  }
  
  #Error model
  error_inis <- data.frame(name = c("a","b","c"),
                           initialValue = c(const_err,0,0),
                           method = c("FIXED","FIXED","FIXED"))
  
  mlx_inis <- rbind(mlx_inis,
                    error_inis)
  
  setPopulationParameterInformation(mlx_inis)
  
  setErrorModel(DV = "constant")
  
  saveProject(mlxtran_name)
  print(paste0("Monolix file saved under: ",mlxtran_name))
  
}

