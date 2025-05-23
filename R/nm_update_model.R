#' Update Thetas and Etas in Model
#' 
#' Update THETA(X) and ETA(Y) to THETAX_cov and ETAY_cov in the NONMEM model.
#' 
#' Should be called first to avoid updating dose THETAs and ETAs.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateThetaEta <- function(optiProject){
  out_model <- optiProject$Model
  
  out_model <- gsub("THETA\\((\\d+)\\)","THETA\\1_cov",out_model) %>%
    {gsub("ETA\\((\\d+)\\)","ETA\\1_cov",.)}
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $INPUT in a NONMEM model
#' 
#' Update the $INPUT line in a NONMEM model
#' 
#' Should be called after \code{nmUpdateThetaEta}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateInputs <- function(optiProject){
  out_model <- optiProject$Model
  
  input_line_nr <- grep("\\$INPUT",out_model)
  
  if(!is.null(optiProject$Parms)){
    parms <- optiProject$Parms$parameter
  } else if(!is.null(optiProject$IParms)){
    parms <- colnames(optiProject$IParms)[-1]
  } else{
    stop("No parameters provided in the optiProject")
  }
  
  input_parms <- paste0(parms,"_cov")
  adms <- lapply(optiProject$Doses,function(x) x$ADM) %>%
    unlist()
  
  if(is.null(adms)){
    new_input <- paste0("$INPUT C   ID   TIME   DV   AMT   EVID   DOSE_NR   ",
                        paste(input_parms,collapse = "   "))
  } else{
    new_input <- paste0("$INPUT C   ID   TIME   DV   AMT   EVID   DOSE_NR   CMT   ",
                        paste(input_parms,collapse = "   "))
  }
  
  out_model[input_line_nr] <- new_input
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $DATA in NONMEM model
#' 
#' Update the $DATA line in a NONMEM model to new optiProject data
#' 
#' Should be called after \code{nmUpdateInputs}
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param new_data_name Name of the new OptiProject-data generated through
#' \code{DataGen}.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateData <- function(optiProject,name){
  new_data_name <- data_name <- paste0(name,"_data.csv")
  out_model <- optiProject$Model
  
  data_line_nr <- grep("\\$DATA",out_model)
  
  out_model[data_line_nr] <- paste0("$DATA ",new_data_name)
  if(!any(grepl("IGNORE=C",out_model))){
    out_model[data_line_nr] <- paste0("$DATA ",new_data_name," IGNORE=C")
  }
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $SUBROUTINES in a NONMEM model
#' 
#' Update the $SUBROUTINGES line in a NONMEM model.
#' 
#' Should be called after \code{nmUpdateData}. Tolerances will be set to 15.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateSubrout <- function(optiProject){
  out_model <- optiProject$Model
  
  subrout_line_nr <- grep("\\$SUBROUTINES",out_model)
  
  if(grepl("TRANS1",out_model[subrout_line_nr])){
    out_model[subrout_line_nr] <- "$SUBROUTINES ADVAN13 TRANS1 TOL=12 ATOL=10"
  } else{
    out_model[subrout_line_nr] <- "$SUBROUTINES ADVAN13 TOL=12 ATOL=10"
  }
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update NCOMPARTMENTS in a NONMEM model
#' 
#' Update the NCOMPARTMENTS argument in a NONMEM model.
#' 
#' Should be called after \code{nmUpdateSubrout}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateNcomp <- function(optiProject){
  out_model <- optiProject$Model
  
  ncomp_line_nr <- grep("NCOMPARTMENTS",out_model)
  if(length(ncomp_line_nr)!=0){
    ncomp_line <- out_model[ncomp_line_nr]
    new_comps <- sum(unlist(lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type))))
    new_ncomp <- regexpr("NCOMPARTMENTS *= *\\d+",ncomp_line) %>%
      {regmatches(ncomp_line,.)} %>%
      {gsub("NCOMPARTMENTS *= *(\\d+)","\\1",.)} %>%
      as.numeric() %>%
      `+`(new_comps)
    new_ncomp_line <- gsub("NCOMPARTMENTS *= *\\d+",paste0("NCOMPARTMENTS=",new_ncomp),ncomp_line)
    out_model[ncomp_line_nr] <- new_ncomp_line
  } else{
    model_line_nr <- grep("$MODEL",out_model,ignore.case = T)
    new_comps <- sum(unlist(lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type))))
    new_ncomp <- gregexpr("A\\(\\d\\)",out_model) %>%
      {regmatches(out_model,.)} %>%
      unlist() %>%
      unique() %>%
      length() %>%
      `+`(new_comps)
    out_model[model_line_nr] <- paste0("$MODEL NCOMPARTMENTS=",new_ncomp)
  }
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $PK in a NONMEM model
#' 
#' Update the $PK line in a NONMEM model, defining the dose levels as THETAs and ETAs.
#' 
#' If model is generated with IIV, MU-referencing will be applied.
#' Should be called after \code{nmUpdateNcomp}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param iiv Whether inter-individual variability should be added to the estimated doses.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdatePK <- function(optiProject,iiv=FALSE){
  out_model <- optiProject$Model
  
  dose_levels <- length(optiProject$Doses)
  
  adms <- lapply(optiProject$Doses,function(x) x$ADM) %>%
    unlist()
  
  if(!iiv){
    if(is.null(adms)){
      dose_code <- lapply(1:dose_levels,function(x) {
        out <- paste0("IF(DOSE_NR.EQ.",x,") F1=THETA(",x,")")
      }) %>% unlist()
    } else{
      dose_code <- lapply(1:dose_levels,function(x) {
        out <- paste0("IF(DOSE_NR.EQ.",x,") F",adms[x],"=THETA(",x,")")
      }) %>% unlist()
    }
  } else{
    if(is.null(adms)){
      dose_code <- lapply(1:dose_levels,function(x) {
        out <- list(paste0("LF",x," = LOG(THETA(",x,"))"),
                    paste0("MU_",x," = LF",x),
                    paste0("TF",x," = EXP(MU_",x," + ETA(",x,"))"),
                    paste0("IF(DOSE_NR.EQ.",x,") F1=TF",x)) %>%
          unlist()
      }) %>% unlist()
    } else{
      dose_code <- lapply(1:dose_levels,function(x) {
        out <- list(paste0("LF",x," = LOG(THETA(",x,"))"),
                    paste0("MU_",x," = LF",x),
                    paste0("TF",x," = EXP(MU_",x," + ETA(",x,"))"),
                    paste0("IF(DOSE_NR.EQ.",x,") F",adms[x],"=TF",x)) %>%
          unlist()
      }) %>% unlist()
    }
  }
  
  pk_line <- grep("\\$PK",out_model)
  out_model <- append(out_model,"",after=pk_line)
  pk_line <- grep("\\$PK",out_model)
  out_model <- append(out_model,dose_code,after=pk_line)
  
  auc_constr <- lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)) %>%
    unlist()
  if(any(auc_constr)){
    inis <- lapply(optiProject$Constraints[auc_constr],function(x) x$Ini) %>%
      unlist()
    
    des_line <- grep("\\$DES",out_model)
    out_model <- append(out_model,inis,after=des_line-1)
    out_line <- grep("\\$DES",out_model)
    out_model <- append(out_model,"",after=out_line-1)
  }
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $DES in a NONMEM model
#' 
#' Update the $DES section in a NONMEM model, adding the constraint/target-codes
#' to the $DES section
#' 
#' Should be called after \code{nmUpdatePK}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateDES <- function(optiProject){
  out_model <- optiProject$Model
  
  out_line <- grep("\\$ERROR",out_model)
  
  auc_constr <- lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)) %>%
    unlist()
  if(any(auc_constr)){
    auc_code <- lapply(optiProject$Constraints[auc_constr],function(x) x$AUC_Code) %>%
      unlist()
    
    out_model <- append(out_model,auc_code,after=out_line-1)
    out_line <- grep("\\$ERROR",out_model)
    out_model <- append(out_model,"",after=out_line-1)
    out_line <- grep("\\$ERROR",out_model)
  }
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $ERROR in a NONMEM model
#' 
#' Update the $ERROR section in a NONMEM model, removing previous $ERROR code and 
#' adding $ERROR code for penalization.
#' 
#' Should be called after \code{nmUpdateDES}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param iiv Whether inter-individual variability should be added to the estimated doses.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateError <- function(optiProject,iiv=FALSE){
  out_model <- optiProject$Model
  
  constr_code <- lapply(optiProject$Constraints,function(x) x$Code) %>%
    unlist()
  
  error_start <- grep("\\$ERROR",out_model)
  error_end <- grep("\\$THETA",out_model)
  
  out_model <- out_model[-c((error_start+1):(error_end-1))]
  
  pens <- lapply(optiProject$Constraints,function(x) x$Name) %>%
    unlist()
  
  out_line <- grep("\\$ERROR",out_model)
  out_model <- append(out_model,"",after=out_line)
  out_line <- grep("\\$ERROR",out_model)
  
  if(!iiv){
    out_model <- append(out_model,"Y = PEN",after=out_line)
  } else{
    out_model <- append(out_model,"Y = PEN + EPS(1)",after=out_line)
  }
  
  final_pen <- paste0("PEN = ",paste(pens,collapse = " + "))
  out_model <- append(out_model,final_pen,after=out_line)
  out_line <- grep("\\$ERROR",out_model)
  out_model <- append(out_model,"",after=out_line)
  out_line <- grep("\\$ERROR",out_model)
  out_model <- append(out_model,constr_code,after=out_line)
  out_line <- grep("\\$ERROR",out_model)
  out_model <- append(out_model,"",after=out_line)
  out_line <- grep("\\$ERROR",out_model)
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $THETA in a NONMEM model
#' 
#' Update the $THETA section in a NONMEM model, removing previous $THETA code and
#' adding initial values for dose levels.
#' 
#' Should be called after \code{nmUpdateError}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateTheta <- function(optiProject){
  out_model <- optiProject$Model
  
  theta_start <- grep("\\$THETA",out_model)
  theta_end <- min(grep("\\$",out_model)[grep("\\$",out_model)>theta_start])
  
  out_model <- out_model[-c((theta_start+1):(theta_end-1))]
  
  dose_levels <- length(optiProject$Doses)
  inis <- lapply(optiProject$Doses,function(x) x$Ini) %>%
    unlist()
  
  new_thetas <- lapply(1:dose_levels,function(x){
    out <- paste0("(0.0001, ",inis[x],") ; [F",x,"]")
  }) %>% unlist()
  
  out_model <- append(out_model,"",after=theta_start)
  out_model <- append(out_model,new_thetas,after=theta_start)
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $OMEGA in a NONMEM model
#' 
#' Update the $OMEGA section in a NONMEM model, either removing $OMEGA section
#' entirely if \code{iiv=FALSE}, or removing previous $OMEGA code and adding
#' initial values for dose levels.
#' 
#' Should be called after \code{nmUpdateTheta}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param iiv Whether inter-individual variability should be added to the estimated doses.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateOmega <- function(optiProject,iiv=FALSE){
  out_model <- optiProject$Model
  
  dose_levels <- length(optiProject$Doses)
  new_omegas <- lapply(1:dose_levels,function(x){
    out <- paste0("0.1 ; [etaF",x,"]")
  }) %>% unlist()
  
  omega_start <- grep("\\$OMEGA",out_model)
  
  if(length(omega_start) != 0 & !iiv){
    omega_end <- min(grep("\\$",out_model)[grep("\\$",out_model)>omega_start])
    
    out_model <- out_model[-c((omega_start):(omega_end-1))]
    
    optiProject$Model <- out_model
  } else if(length(omega_start) == 0 & !iiv){
    optiProject$Model <- out_model
  } else if(length(omega_start) != 0 & iiv){
    omega_end <- min(grep("\\$",out_model)[grep("\\$",out_model)>omega_start])
    
    out_model <- out_model[-c((omega_start+1):(omega_end-1))]
    
    out_model <- append(out_model,"",after=omega_start)
    out_model <- append(out_model,new_omegas,after=omega_start)
    
    optiProject$Model <- out_model
  } else if(length(omega_start) == 0 & iiv){
    est_start <- grep("\\$EST",out_model)
    
    out_model <- append(out_model,"",after=omega_start)
    out_model <- append(out_model,c("$OMEGA",new_omegas),after=est_start)
    
    optiProject$Model <- out_model
  }
  return(optiProject)
}

#' Update $SIGMA in a NONMEM model
#' 
#' Update the $SIGMA section in a NONMEM model, either removing $SIGMA section
#' entirely if \code{iiv=FALSE}, or removing previous $SIGMA code and adding
#' initial values for residual error parameter (not yet implemented).
#' 
#' Should be called after \code{nmUpdateOmega}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param iiv Whether inter-individual variability should be added to the estimated doses.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateSigma <- function(optiProject,iiv=FALSE){
  out_model <- optiProject$Model
  
  sigma_start <- grep("\\$SIGMA",out_model)
  sigma_end <- min(grep("\\$",out_model)[grep("\\$",out_model)>sigma_start])
  
  if(!iiv){
    out_model <- out_model[-c((sigma_start):(sigma_end-1))]
    
    optiProject$Model <- out_model
  } else{
    out_model <- out_model[-c((sigma_start+1):(sigma_end-1))] %>%
      append("0.001 FIX",after=sigma_start) %>%
      append("",after=sigma_start+1)
    
    optiProject$Model <- out_model
  }
  
  return(optiProject)
}

#' Update $EST in a NONMEM model
#' 
#' Update the $EST line in a NONMEM model, removing previous $EST lines and adding
#' $EST line for optiProject (with iiv not yet implemented).
#' 
#' If model is with IIV, SAEM will be used. Should be called after \code{nmUpdateSigma}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param iiv Whether inter-individual variability should be added to the estimated doses.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateEst <- function(optiProject,iiv=FALSE){
  out_model <- optiProject$Model
  
  est_line_nr <- grep("\\$EST",out_model)
  
  out_model <- out_model[-est_line_nr]
  
  if(!iiv){
    new_est_line <- "$EST -2LL PRINT=50 MAXEVAL=99999 NSIG=6"
    out_model <- append(out_model,new_est_line,after=min(est_line_nr)-1)
  } else{
    new_est_line <- "$EST METHOD=SAEM AUTO=1 NITER=500 PRINT=20"
    out_model <- append(out_model,new_est_line,after=min(est_line_nr)-1)
  }
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update $TABLE in a NONMEM model
#' 
#' Update the $TABLE line in a NONMEM model, removing previous $TABLE lines and adding
#' outputs TIME PEN
#' 
#' Should be called after \code{nmUpdateEst}.
#' 
#' @param optiProject A optiProject-object with a model added.
#' @param iiv Whether inter-individual variability should be added to the estimated doses.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateTable <- function(optiProject,name,iiv=FALSE){
  out_model <- optiProject$Model
  
  table_line_nr <- grep("\\$TABLE",out_model)
  
  out_model <- out_model[-table_line_nr]
  
  if(!iiv){
    new_table_line <- paste0("$TABLE TIME PEN NOAPPEND FILE=",name,".tab")
    out_model <- append(out_model,new_table_line,after=min(table_line_nr)-1)
  } else{
    new_table_line <- paste0("$TABLE TIME PEN NOAPPEND FILE=",name,".tab")
    out_model <- append(out_model,new_table_line,after=min(table_line_nr)-1)
  }
  
  optiProject$Model <- out_model
  return(optiProject)
}

#' Update NONMEM model for multiple seperate individuals
#' 
#' RECS=ID is added to the model $DATA line and INCLUDE of the multiple fitting
#' text file at the end of the model.
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object
#' @param name Name under which the optiProject is saved.
#' @param sim Whether individuals were simulated (TRUE) or IParms was used (FALSE).
#' @param n_id Number of simulated individuals.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nmUpdateMulti <- function(optiProject,name,sim,n_id){
  out_model <- optiProject$Model
  
  if(sim & n_id>1){
    reps <- n_id - 1
  } else if(sim & n_id <= 1){
    stop("If sim=TRUE, n_id must be greater than 1.")
  } else if(!sim & !is.null(optiProject$IParms)){
    if(nrow(optiProject$IParms) > 1){
      reps <- nrow(optiProject$IParms) - 1
    } else{
      stop("IParm must include more than 1 subject.")
    }
  }
  
  data_line_nr <- grep("\\$DATA",out_model)
  out_model[data_line_nr] <- paste0(out_model[data_line_nr]," RECS=ID")
  
  multi_name <- paste0(name,"_multi.txt")
  
  
  out_model <- out_model %>%
    append("") %>%
    append(paste0("INCLUDE ",multi_name," ",reps))
  
  optiProject$Model <- out_model
  return(optiProject)
}