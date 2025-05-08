#' Update ini section of model
#' 
#' Remove model parameters and add dosing parameters.
#' 
#' Should be called first.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param iiv Whether doses should contain IIV.
#' @param const_err Constant error for error model.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nlmixrUpdateIni <- function(optiProject,pop,iiv,const_err){
  out_model <- optiProject$Model
  
  ini_start <- grep("ini\\(\\{",out_model)
  ini_end <- min(grep("}",out_model)[grep("\\}",out_model) > ini_start])
  
  out_model <- out_model[-c((ini_start+1):(ini_end-1))]
  
  ini_start <- grep("ini\\(\\{",out_model)
  out_model <- append(out_model,"",after=ini_start)
  if(pop){
    out_model <- append(out_model,paste0("        add.err <- fix(",const_err,")"),after=ini_start)
  } else{
    out_model <- append(out_model,paste0("        add.err <- c(0,",const_err,",",const_err*1.1,")"),after=ini_start)
  }
  
  dose_levels <- length(optiProject$Doses)
  inis <- lapply(optiProject$Doses,function(x) x$Ini) %>%
    unlist()
  
  if(iiv){
    eta_doses <- paste0("        eta.F",1:dose_levels," ~ 0.1")
    out_model <- append(out_model,eta_doses,after=ini_start)
    out_model <- append(out_model,"",after=ini_start)
  }
  
  theta_doses <- paste0("        lF",1:dose_levels," <- log(",inis,")")
  out_model <- append(out_model,theta_doses,after=ini_start)
  out_model <- append(out_model,"",after=ini_start)
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update model section initial states
#' 
#' Adding initial states for AUC constraints to the model section.
#' 
#' Should be called after \code{nlmixrUpdateModelInis}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
nlmixrUpdateModelInis <- function(optiProject){
  out_model <- optiProject$Model
  
  model_start <- grep("model\\(\\{",out_model)
  
  auc_constr <- lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)) %>%
    unlist()
  
  if(any(auc_constr)){
    inis <- lapply(optiProject$Constraints[auc_constr],function(x) x$Ini) %>%
      unlist()
    
    out_model <- append(out_model,paste0("        ",inis),after=model_start)
    optiProject$Model <- out_model
  }
  
  return(optiProject)
}

#' Update model section dose definition
#' 
#' Adding dose definitions to the model section.
#' 
#' Should be called after \code{nlmixrUpdateModelInis}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param iiv Whether doses should contain IIV.
#' @author Dominic Bräm
nlmixrUpdateModelParms <- function(optiProject,iiv){
  out_model <- optiProject$Model
  
  model_start <- grep("model\\(\\{",out_model)
  dose_levels <- length(optiProject$Doses)
  
  if(!iiv){
    doses <- paste0("        F",1:dose_levels," = exp(lF",1:dose_levels,")")
  } else{
    doses <- paste0("        F",1:dose_levels," = exp(lF",1:dose_levels," + eta.F",1:dose_levels,")")
  }
  
  out_model <- append(out_model,c(doses,""),after=model_start)
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update model section bioavailability definition
#' 
#' Adding bioavailability to the model section.
#' 
#' Should be called after \code{nlmixrUpdateModelParms}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nlmixrUpdateModelBioav <- function(optiProject){
  out_model <- optiProject$Model
  
  model_end <- grep("\\}",out_model) %>% `[`(length(.)-1)
  
  dose_levels <- length(optiProject$Doses)
  
  adms <- lapply(optiProject$Doses,function(x) x$ADM) %>%
    unlist()
  
  comps <- regexpr("d/dt\\([^\\)]+\\)",out_model) %>%
    {regmatches(out_model,.)} %>%
    {gsub("d/dt|\\(|\\)","",.)}
  
  if(is.null(adms)){
    bioav_lines <- lapply(1:dose_levels,function(x){
      out <- list(paste0("        if(DOSE_NR == ",x,"){"),
                  paste0("           F0 = F",x),
                  paste0("        }"))
      return(out)
    }) %>% unlist()
  } else{
    bioav_lines <- lapply(1:dose_levels,function(x){
      out <- list(paste0("        if(DOSE_NR == ",x,"){"),
                  paste0("           F0_",adms[x]," = F",x),
                  paste0("        }"))
      return(out)
    }) %>% unlist()
  }
  
  if(is.null(adms)){
    bioav_lines <- c(bioav_lines,paste0("f(",comps[1],") = F0"))
  } else{
    boav_lines <- c(bioav_lines,paste0("f(",comps[adms],") = F0_",adms))
  }
  
  
  out_model <- append(out_model,c("",bioav_lines),after=model_end-1)
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update model section penalisation code
#' 
#' Adding penalisation code to the model section.
#' 
#' Should be called after \code{nlmixrUpdateModelBioav}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nlmixrUpdateModelPens <- function(optiProject){
  out_model <- optiProject$Model
  
  model_end <- grep("\\}",out_model) %>% `[`(length(.)-1)
  
  auc_constr <- lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)) %>%
    unlist()
  if(any(auc_constr)){
    auc_code <- lapply(optiProject$Constraints[auc_constr],function(x) x$AUC_Code) %>%
      unlist() %>%
      {paste0("        ",.)}
    
    out_model <- append(out_model,c(auc_code,""),after=model_end-1)
  }
  
  model_end <- grep("\\}",out_model) %>% `[`(length(.)-1)
  
  constr_code <- lapply(optiProject$Constraints,function(x) x$Code) %>%
    unlist() %>%
    {paste0("        ",.)}
  
  out_model <- append(out_model,c(constr_code,""),after=model_end-1)
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update model section output code
#' 
#' Removing previous model output and adding penalization output to the model section.
#' 
#' Should be called after \code{nlmixrUpdateModelPens}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nlmixrUpdateModelOut <- function(optiProject){
  out_model <- optiProject$Model
  
  og_out_lines <- grep("prop \\(|add\\(",out_model)
  out_model <- out_model[-og_out_lines]
  
  model_end <- grep("\\}",out_model) %>% `[`(length(.)-1)
  
  pens <- lapply(optiProject$Constraints,function(x) x$Name) %>%
    unlist()
  
  out_lines <- list(paste0("        PEN = ",paste(pens,collapse = " + ")),
                    "        PEN ~ add(add.err)") %>%
    unlist()
  out_model <- append(out_model,out_lines,after=model_end-1)
  
  optiProject$Model <- out_model
  return(optiProject)
  
}

#' Update model if population without IIV
#' 
#' If a single dose level is estimated with a population but without IIV, the estimation
#' algorithm "bobyqa" seems to struggle. Therfore, a dummy parameter "tmp" is introduced.
#' 
#' Should be called only if pop=TRUE and IIV=FALSE, after \code{nlmixrUpdateModelOut}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
nlmixrUpdatePopNoIIV <- function(optiProject){
  out_model <- optiProject$Model
  
  ini_start <- grep("ini\\(\\{",out_model)
  out_model <- append(out_model,"        tmp <- c(0,1,2)", after = ini_start)
  
  model_start <- grep("model\\(\\{",out_model)
  out_model <- append(out_model,"        tmptmp <- exp(tmp)", after = model_start)
  
  optiProject$Model <- out_model
  return(optiProject)
}