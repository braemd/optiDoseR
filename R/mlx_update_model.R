#' Update structural model parameters
#' 
#' Define structural model parameters through regressors.Values for structural model 
#' parameters are given through regressors in the data, thus, in the PK section, 
#' structural model parameters are defined as parmX = parmX_cov where
#' parmX_cov is a regressor in the data. 
#' 
#' Should be called first.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
mlxUpdateParmNames <- function(optiProject){
  out_model <- optiProject$Model
  
  if(!any(grepl("input *= *\\{[^\\}]+\\}",optiProject$Model))){
    input_start <- grep("input *= *\\{",optiProject$Model)
    input_end <- min(grep("\\}",optiProject$Model))
    input_line <- paste(optiProject$Model[input_start:input_end],collapse = "")
    out_model <- out_model[-input_end]
  } else{
    input_start <- grep("input *= *\\{",optiProject$Model)
    input_line <- optiProject$Model[input_start]
  }
  
  pmx_parms <- gsub("input *= *\\{([^\\}]*)\\}","\\1",input_line) %>%
    strsplit(",") %>%
    `[[`(1)
  
  parms_def <- lapply(pmx_parms,function(x) {
    paste0(x," = ",x,"_cov")
  }) %>%
    unlist()
  
  pk_line <- grep("PK:",out_model)
  out_model <- append(out_model,parms_def,after=pk_line)
  pk_line <- grep("PK:",out_model)
  out_model <- append(out_model,"",after=pk_line)
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update input
#' 
#' Update input to the Monolix model with F1-FX for dose levels, DOSE_NR as regressor to define 
#' when which dose level should be applied, and parmX_cov as regressor for parmX in the 
#' structural model.
#' 
#' Should be called after \code{mlxUpdateParmNames}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
mlxUpdateInputs <- function(optiProject){
  out_model <- optiProject$Model
  
  if(!any(grepl("input *= *\\{[^\\}]+\\}",optiProject$Model))){
    input_start <- grep("input *= *\\{",optiProject$Model)
    input_end <- min(grep("\\}",optiProject$Model))
    input_line <- paste(optiProject$Model[input_start:input_end],collapse = "")
    out_model <- out_model[-input_end]
  } else{
    input_start <- grep("input *= *\\{",optiProject$Model)
    input_line <- optiProject$Model[input_start]
  }
  
  pmx_parms <- gsub("input *= *\\{([^\\}]*)\\}","\\1",input_line) %>%
    strsplit(",") %>%
    `[[`(1) %>%
    paste0("_cov")
  
  n_dose_levels <- length(optiProject$Doses)
  dose_parms <- paste0("F",1:n_dose_levels)
  
  updated_input_line <- paste0("input = {",
                               paste(dose_parms,collapse = ","),
                               ",",
                               paste(pmx_parms,collapse = ","),
                               ",",
                               "DOSE_NR}")
  
  dose_regressor <- "DOSE_NR = {use=regressor}"
  parm_regressor <- paste0(pmx_parms," = {use=regressor}")
  
  regressor <- c(dose_regressor,
                 parm_regressor)
  
  out_model[input_start] <- updated_input_line
  out_model <- append(out_model,regressor,after=input_start)
  
  optiProject$Model <- out_model
  
  return(optiProject)
}

#' Update Equation section
#' 
#' Update EQUATION: section of the Monolix model by adding the penalization code of
#' the optiProject-object constraints
#' 
#' Should be called after \code{mlxUpdateBioav}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
mlxUpdateEquation <- function(optiProject){
  out_model <- optiProject$Model
  
  out_line <- grep("OUTPUT:",out_model)
  eq_line <- grep("EQUATION:",out_model)
  if(length(eq_line)==0){
    out_model <- append(out_model,"EQUATION:",after=out_line-1)
  }
  
  constr_code <- lapply(optiProject$Constraints,function(x) x$Code) %>%
    unlist()
  
  pens <- lapply(optiProject$Constraints,function(x) x$Name) %>%
    unlist()
  
  auc_constr <- lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)) %>%
    unlist()
  if(any(auc_constr)){
    inis <- lapply(optiProject$Constraints[auc_constr],function(x) x$Ini) %>%
      unlist()
    
    auc_code <- lapply(optiProject$Constraints[auc_constr],function(x) x$AUC_Code) %>%
      unlist()
    
    out_line <- grep("OUTPUT",out_model)
    out_model <- append(out_model,inis,after=out_line-1)
    out_line <- grep("OUTPUT",out_model)
    out_model <- append(out_model,"",after=out_line-1)
    out_line <- grep("OUTPUT",out_model)
    out_model <- append(out_model,auc_code,after=out_line-1)
    out_line <- grep("OUTPUT",out_model)
    out_model <- append(out_model,"",after=out_line-1)
    out_line <- grep("OUTPUT",out_model)
  }
  
  out_model <- append(out_model,constr_code,after=out_line-1)
  out_line <- grep("OUTPUT",out_model)
  out_model <- append(out_model,"",after=out_line-1)
  out_line <- grep("OUTPUT",out_model)
  
  final_pen <- paste0("pen = ",paste(pens,collapse = " + "))
  out_model <- append(out_model,final_pen,after=out_line-1)
  out_line <- grep("OUTPUT",out_model)
  out_model <- append(out_model,"",after=out_line-1)
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update dosing macros
#' 
#' Update dosing macros, such as depot, absorption, iv, and pkmodel, by adding 
#' bioavailability for dose levels.
#' 
#' Should be called after \code{mlxUpdatePK}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
mlxUpdateBioav <- function(optiProject){
  out_model <- optiProject$Model
  
  admin_macros <- c("depot\\([^\\)]+\\)","absorption\\([^\\)]+\\)","iv\\([^\\)]+\\)","pkmodel\\([^\\)]+\\)")
  adm_line_nr <- grep(paste(admin_macros,collapse = "|"),out_model)
  adm_lines <- out_model[adm_line_nr]
  
  model_macros <- regexpr("([^\\)]+)",adm_lines) %>%
    {regmatches(adm_lines,.)}
  adjusted_macros <- paste0(model_macros,",p=F)")
  
  out_model[adm_line_nr] <- adjusted_macros
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update PK section
#' 
#' Update PK section to define estimated dose levels FX as bioavailability factor
#' F.
#' 
#' Should be called after \code{mlxUpdateInputs}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
mlxUpdatePK <- function(optiProject){
  out_model <- optiProject$Model
  
  dose_levels <- length(optiProject$Doses)
  dose_code <- lapply(1:dose_levels,function(x) {
    out <- list(paste0("if DOSE_NR==",x),
                paste0("   F = F",x),
                "end")
  }) %>% unlist()
  
  pk_line <- grep("PK:",out_model)
  out_model <- append(out_model,dose_code,after=pk_line)
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update output
#' 
#' Update output of the structural model by removing previous outputs and adding
#' penalization.
#' 
#' Should be called after \code{mlxUpdateEquation}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
mlxUpdateOutput <- function(optiProject){
  out_model <- optiProject$Model
  
  out_line_nr <- grep("output *= *",out_model)
  og_outputs <- gsub("output *= *","",out_model[out_line_nr])
  out_model[out_line_nr] <- "output = pen"
  out_model <- append(out_model,paste0("table = ",og_outputs),after=out_line_nr)
  
  optiProject$Model <- out_model
  return(optiProject)
}
