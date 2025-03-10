#' Get maximal number of compartments
#' 
#' Get the maximal number of compartments in the structural model. Needed to define
#' the compartment number for penalizations with AUC.
#' 
#' Looking for highest DADT(X) in the model.
#' 
#' @param model_text Structural NONMEM model.
#' @author Dominic Bräm
nm_nstate_extractor <- function(model_text){
  spoter <- regexpr("DADT\\(\\d+\\)",model_text)
  matches <- regmatches(model_text,spoter)
  numbers <- unlist(gsub("DADT|\\(|\\)","",matches))
  max_number <- max(as.numeric(numbers))
  return(max_number)
}

#' Update Thetas and Etas in Model
#' 
#' Update THETA(X) and ETA(Y) to THETAX_cov and ETAY_cov in the NONMEM model.
#' 
#' Should be called first to avoid updating dose THETAs and ETAs.
#' 
#' @param model_text A structural NONMEM model.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
update_nm_model_parms <- function(model_text){
  updated_model <- gsub("THETA\\((\\d+)\\)","THETA\\1",model_text) %>%
    {gsub("ETA\\((\\d+)\\)","ETA\\1",.)}
  
  return(updated_model)
}

#' Generate code for multiple seperate subjects in NONMEM
#' 
#' To estimate multiple subject seperately in NONMEM $PROB, $INPUT, $DATA, 
#' $THETA, $EST, and $TABLE can be saved in a seperate file and included in the
#' control file
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object
#' @param name Name under which the optiProject is saved.
#' @param save_path Path where optiProject is saved.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
nm_multi_sep_code <- function(optiProject,name,save_path){
  out_model <- optiProject$Model
  
  prob <- out_model[grep("\\$PROB",out_model)]
  input <- out_model[grep("\\$INPUT",out_model)]
  data <- out_model[grep("\\$DATA",out_model)] %>%
    paste0(" NOREWIND")
  
  theta_start <- grep("\\$THETA",out_model)
  theta_end <- min(grep("\\$",out_model[(theta_start+1):length(out_model)]))
  theta <- out_model[theta_start:(theta_start+theta_end)]
  
  est <- out_model[grep("\\$EST",out_model)]
  table <- out_model[grep("\\$TABLE",out_model)]
  
  if(length(prob)>0){
    multi_sep_text <- list(prob)
  } else{
    multi_sep_text <- list()
  }
  
  multi_sep_text <- multi_sep_text %>%
    append(input) %>%
    append("") %>%
    append(data) %>%
    append("") %>%
    append(theta) %>%
    append("") %>%
    append(est) %>%
    append("") %>%
    append(table) %>%
    unlist()
  
  multi_name <- paste0(save_path,"/",name,"_multi.txt")
  writeLines(multi_sep_text,con = multi_name)
  
  print(paste0("Multi-fitting file has been saved under: ",multi_name))
}
