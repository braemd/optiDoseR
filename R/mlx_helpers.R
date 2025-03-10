#' Check for administration macros
#' 
#' The optiDose package only works if administration macros, i.e., depot, absorption,
#' iv, or pkmodel were used in the structural model. This function checks whether this is 
#' the case and if multiple administration IDs are present in the dose levels, whether 
#' these were already defined in the structural model.
#' 
#' NULL
#' 
#' @param model Monolix structural model.
#' @author Dominic Bräm
mlx_admin_checker <- function(model){
  admin_macros <- c("depot","absorption","iv","pkmodel")
  model <- model[grep("\\[LONGITUDINAL\\]",model):length(model)]
  
  admin_line_nrs <- grep(paste(admin_macros,collapse = "|"),model)
  if(length(admin_line_nrs)==1){
  } else if(length(admin_line_nrs)==0){
    stop("No administration macro identified in the model")
  } else if(length(admin_line_nrs)>1){
    adms <- grepl("adm",model[admin_line_nrs])
    if(!all(adms)){
      stop("If multiple administration macros are present, all need a adm argument")
    }
  }
  
  ps <- grepl("p *= *",model[admin_line_nrs])
  if(any(ps)){
    warning("Bioavailability parameters in the model will be removed for dose optimization")
  }
}

#' Identify header types
#' 
#' Defines header types of optiProject-data.
#' 
#' NULL
#' 
#' @param data Dataframe of the optiProject, generated through DataGen.
#' @author Dominic Bräm
mlx_header_types <- function(data){
  col_names <- colnames(data)
  col_names[col_names == "ID"] <- "id"
  col_names[col_names == "TIME"] <- "time"
  col_names[col_names == "DV"] <- "observation"
  col_names[col_names == "AMT"] <- "amount"
  col_names[col_names == "EVID"] <- "evid"
  col_names[col_names == "ADM"] <- "admid"
  col_names[col_names == "DOSE_NR"] <- "regressor"
  col_names[grepl("_cov",col_names)] <- "regressor"
  
  return(col_names)
}
