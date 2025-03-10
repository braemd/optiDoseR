#' Generate data for an optiProject
#' 
#' Transforming targets, dose levels, and parameters of an optiProject to a data file
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object with a model, penalizations, dose levels, and parameters.
#' For requirements regarding penalizations, dose levels, and parameters.
#' @param pop Whether multiple subjects should be used (TRUE), either from added individual parameters or
#' simulated from population parameters, or if only one subject with population parameters should be used
#' (FALSE).
#' @param iiv Whether doses should include inter-individual variability (TRUE) or not (FALSE).
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
#' @importFrom magrittr %>%
#' @import checkmate
#' @import dplyr
#' @author Dominic Bräm
DataGen <- function(optiProject,
                    pop=FALSE,
                    sim=FALSE,
                    n_id=1,
                    parametrisation="log",
                    lower_parm_limits=NULL,
                    upper_parm_limits=NULL){
  
  checkmate::assertClass(optiProject,"optiProject")
  checkmate::assert(checkmate::checkTRUE(pop),checkmate::checkFALSE(pop))
  checkmate::assert(checkmate::checkTRUE(sim),checkmate::checkFALSE(sim))
  checkmate::assertNumeric(n_id,len=1,lower=1)
  checkmate::assertCharacter(parametrisation,pattern="log|normal|logit")
  checkmate::assert(checkmate::checkNumeric(upper_parm_limits,lower=1e-3),checkmate::checkNull(upper_parm_limits))
  checkmate::assert(checkmate::checkNumeric(lower_parm_limits,lower=0),checkmate::checkNull(lower_parm_limits))
  if(length(upper_parm_limits) != sum(parametrisation=="logit") | length(upper_parm_limits) != sum(parametrisation=="logit")){
    stop("For each parameter with logit-parametrisation, a lower and upper limit must be given")
  }
  if(!is.null(lower_parm_limits)){
    if(any(lower_parm_limits > upper_parm_limits)){
      stop("Lower parameter limits must be smaler than upper parameter limits")
    }
  }
  if(is.null(lower_parm_limits)){
    lower_parm_limits <- NA
  }
  if(is.null(upper_parm_limits)){
    upper_parm_limits <- NA
  }
  
  software <- optiProject$Software
  
  dose_data_single <- DoseDataGen(optiProject$Doses)
  
  const_data_single <- ConstraintDataGen(optiProject$Constraints)
  if(!("ADM" %in% colnames(dose_data_single))){
    const_data_single <- const_data_single %>%
      dplyr::select(-ADM)
  }
  
  if(!pop){
    if(sim | n_id>1){
      warning("No population will be generated, sim and n_id is ignored")
    }
    parm_data <- PParmDataGen(optiProject$Parms,software)
  } else{
    if(sim){
      parm_data <- IParmDataSim(optiProject$Parms,n_id,software,parametrisation,
                                lower_parm_limits,upper_parm_limits)
    } else{
      parm_data <- optiProject$IParms %>%
        dplyr::rename("ID" = "id")
      n_id <- nrow(parm_data)
    }
  }
  
  n_rows_per_id <- sum(nrow(dose_data_single),
                       nrow(const_data_single))
  
  if(is.factor(parm_data$ID)){
    parm_data$ID <- as.numeric(as.character(parm_data$ID))
  }
  
  if(pop & sim){
    data <- rbind(dose_data_single,
                  const_data_single) %>%
      `[`(rep(1:n_rows_per_id,n_id),) %>%
      dplyr::mutate(ID=rep(1:n_id,each=n_rows_per_id)) %>%
      dplyr::arrange(ID,TIME) %>%
      dplyr::left_join(parm_data,by="ID")
  } else if(pop & !sim){
    data <- rbind(dose_data_single,
                  const_data_single) %>%
      `[`(rep(1:n_rows_per_id,n_id),) %>%
      dplyr::mutate(ID=rep(unique(parm_data$ID),each=n_rows_per_id)) %>%
      dplyr::arrange(ID,TIME) %>%
      dplyr::left_join(parm_data,by="ID")
  } else{
    data <- rbind(dose_data_single,
                  const_data_single) %>%
      dplyr::arrange(ID,TIME) %>%
      dplyr::left_join(parm_data,by="ID")
  }
  
  if(software == "NONMEM"){
    data <- data %>%
      dplyr::mutate(C = ".") %>%
      dplyr::relocate("C")
  }
  
  if(software %in% c("Monolix","NONMEM")){
    data_cols <- colnames(data)
    cols_to_paste <- !(data_cols %in% c("C","ID","TIME","DV","AMT","EVID","ADM","DOSE_NR"))
    data_cols[cols_to_paste] <- paste0(data_cols[cols_to_paste],"_cov")
    colnames(data) <- data_cols
  }
  
  if(software == "NONMEM" & "ADM" %in% colnames(data)){
    data <- data %>%
      dplyr::rename("CMT" = "ADM")
  }
  
  if(software == "nlmixr2"){
    data <- data %>%
      dplyr::mutate(MDV = ifelse(DV==".",1,0),
                    DV = as.numeric(ifelse(DV==".",0,DV)),
                    DOSE_NR = as.numeric(ifelse(DOSE_NR==".",0,DOSE_NR))) %>%
      dplyr::relocate(MDV,.after=EVID)
  }
  
  if(software == "nlmixr2" & "ADM" %in% colnames(data)){
    data <- data %>%
      dplyr::rename("CMT" = "ADM")
  }
  
  return(data)
}

#' Generate dose data
#' 
#' Generate dose data from the dose level information of an optiProject-object
#' 
#' NULL
#' 
#' @param doses Dose levels in an optiProject-object added through \code{addDoseLevel}.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
DoseDataGen <- function(doses){
  dose_times <- lapply(doses,
                       function(x) x$Times) %>%
    unlist()
  
  dose_adms <- lapply(doses,
                      function(x) rep(x$ADM,length(x$Times))) %>%
    unlist()
  
  if(!(all(is.null(dose_adms)) | all(!is.null(dose_adms)))){
    stop("ADM must be given either for all or none dose level")
  }
  
  dose_numbers <- lapply(1:length(doses),
                         function(x) rep(x,length(doses[[x]]$Times))) %>%
    unlist()
  
  if(all(is.null(dose_adms))){
    dose_data <- data.frame(ID = 1,
                            TIME = dose_times,
                            DV = ".",
                            AMT = 1,
                            EVID = 1,
                            DOSE_NR = dose_numbers)
  } else{
    dose_data <- data.frame(ID = 1,
                            TIME = dose_times,
                            DV = ".",
                            AMT = 1,
                            EVID = 1,
                            ADM = dose_adms,
                            DOSE_NR = dose_numbers)
  }
  
  dose_data <- dose_data %>%
    dplyr::arrange(TIME)

  return(dose_data)
}

#' Generate penalization data
#' 
#' Generate penalization data from the constraints information of an optiProject-object
#' 
#' NULL
#' 
#' @param constraints Penalization constraints/targets in an optiProject-object.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
ConstraintDataGen <- function(constraints){
  constraint_times <- lapply(constraints,
                             function(x) x$Time) %>%
    unlist() %>%
    unique()
  
  const_data <- data.frame(ID = 1,
                           TIME = constraint_times,
                           DV = 0,
                           AMT = ".",
                           EVID = 0,
                           ADM = ".",
                           DOSE_NR = ".") %>%
    dplyr::arrange(TIME)
  
  return(const_data)
}

#' Generate parameter data for typical subject
#' 
#' Generate  parameter data for the typcial subject from the population information of
#' an optiProject-object
#' 
#' NULL
#' 
#' @param parms Population parameter data frame with columns for parameter name and value,
#' including at least \emph{_pop} parameters for Monolix or \emph{THETA} parameters for NONMEM
#' @param software Software of the optiProject-object
#' @author Dominic Bräm
#' @importFrom magrittr %>%
#' @importFrom tibble column_to_rownames
#' @import dplyr
PParmDataGen <- function(parms,software){
  if(software == "NONMEM"){
    pparms <- parms %>%
      dplyr::mutate(value = ifelse(grepl("THETA",parameter),value,0)) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      as.data.frame()
  } else if(software == "Monolix"){
    pparms <- parms %>%
      dplyr::filter(grepl("_pop",parameter)) %>%
      dplyr::mutate(parameter = gsub("_pop","",parameter)) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      as.data.frame()
  } else if(software == "nlmixr2"){
    pparms <- parms %>%
      dplyr::mutate(value = ifelse(grepl("eta\\.",parameter),0,value)) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      as.data.frame()
  }
  
  pparm_data <- data.frame(ID=1) %>%
    cbind(pparms)
  
  return(pparm_data)
}

#' Generate parameter data for individual subjects
#' 
#' Generate parameter data for individual subjects by simulating from the population 
#' information of an optiProject-object
#' 
#' NULL
#' 
#' @param parms Population parameters to simulate from.
#' @param n_id Number of subjects to be simulated.
#' @param software Software of the optiProject-object.
#' @param parametrisation If \code{software="Monolix"}, the parametrisation of the parameters in the Monolix model
#' must be given. Must be in "log","normal",or "logit". If all parameters have same parametrisation,
#' a single string can be given. If different parametrisations for different parameters are required, a vector
#' with parametrisation for each parameter is required. Default is "log" for all parameters.
#' @param lower_parm_limits If one or multiple parameters in the model are "logit" distributed and individuals
#' are simulated, the lower limit(s) for the "logit"-distributed parameter(s).
#' @param upper_parm_limits If one or multiple parameters in the model are "logit" distributed and individuals
#' are simulated, the upper limit(s) for the "logit"-distributed parameter(s).
#' @author Dominic Bräm
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
IParmDataSim <- function(parms,n_id,software,parametrisation,
                         lower_parm_limit,upper_parm_limit){
  if(software == "Monolix"){
    if(length(parametrisation) != 1 & length(parametrisation) != sum(grepl("_pop",parms$parameter))){
      stop("parametrisation must be of length 1 or n_id")
    }
    
    sim_parm_data <- parms %>%
      dplyr::mutate(parm_type = dplyr::case_when(grepl("_pop",parameter) ~ "pop",
                                                 grepl("omega_",parameter) ~ "omega",
                                                 TRUE ~ "0")) %>%
      dplyr::mutate(parameter = gsub("_pop|omega_","",parameter)) %>%
      dplyr::filter(parm_type != "0") %>%
      tidyr::pivot_wider(values_from = value,
                         names_from = parm_type,
                         values_fill = 0) %>%
      dplyr::mutate(pop = as.numeric(pop),
             omega = as.numeric(omega)) %>%
      dplyr::mutate(parmet = parametrisation,
                    upper = NA,
                    lower = NA) %>%
      `[<-`(.$parmet=="logit","upper",value=c(upper_parm_limit)) %>%
      `[<-`(.$parmet=="logit","lower",value=c(lower_parm_limit)) %>%
      `rownames<-`(NULL) %>%
      {suppressWarnings({
        apply(.,1,function(x) {
          varname <- x["parameter"]
          pop <- as.numeric(x["pop"])
          omega <- as.numeric(x["omega"])
          if(x["parmet"]=="log"){
            out <- data.frame(id = 1:n_id) %>%
                          dplyr::mutate(!!varname := pop * exp(rnorm(n_id,0,omega))) %>%
              dplyr::select(-id)
          } else if(x["parmet"]=="normal"){
            out <- data.frame(id = 1:n_id) %>%
                          dplyr::mutate(!!varname := pop + rnorm(n_id,0,omega)) %>%
              dplyr::select(-id)
          } else if(x["parmet"]=="logit"){
            up <- as.numeric(x["upper"])
            lo <- as.numeric(x["lower"])
            eta <- rnorm(n_id,0,omega)
            out <- data.frame(id = 1:n_id) %>%
              dplyr::mutate(!!varname := (up*(pop-lo)/(up-pop)*exp(eta)+lo)/(1+(pop-lo)/(up-pop)*exp(eta))) %>%
              dplyr::select(-id)
          }
        })
      })} %>%
      do.call(cbind,.) %>%
      dplyr::mutate(ID = 1:n_id) %>%
      dplyr::relocate("ID")
      
    return(sim_parm_data)
  } else if(software == "NONMEM"){
    thetas <- parms %>%
      dplyr::filter(grepl("THETA",parameter)) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      as.data.frame() %>%
      `[`(rep(1,n_id),)
    
    etas <- parms %>%
      dplyr::filter(!grepl("THETA",parameter)) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      as.data.frame() %>%
      apply(2,function(x) {
        out <- rnorm(n_id,0,as.numeric(x))
        return(out)
      })
    
    sim_parm_data <- data.frame(ID = 1:n_id) %>%
      cbind(thetas,etas)
    return(sim_parm_data)
  } else if(software == "nlmixr2"){
    thetas <- parms %>%
      dplyr::filter(!grepl("eta\\.",parameter)) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      as.data.frame() %>%
      `[`(rep(1,n_id),)
    
    etas <- parms %>%
      dplyr::filter(grepl("eta\\.",parameter)) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      as.data.frame() %>%
      apply(2,function(x) {
        out <- rnorm(n_id,0,as.numeric(x))
        return(out)
      })
    
    sim_parm_data <- data.frame(ID = 1:n_id) %>%
      cbind(thetas,etas)
    return(sim_parm_data)
  }
}

