#' General estimates extractor
#' 
#' Gets estimated population and individual parameters from a Monolix or NONMEM run.
#' 
#' NULL
#' 
#' @param pmx_file Name of either a Monolix .mlxtran file or a NONMEM results file.
#' @param software Initialized software of the optiProject-object.
#' @param phi_file If \code{software="NONMEM"}, a .phi file of the NONMEM run to extract
#' individual etas.
#' @author Dominic Bräm
#' @importFrom dplyr select
#' @importFrom dplyr relocate
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom tools file_path_sans_ext
#' @importFrom tibble column_to_rownames
parm_extractor <- function(pmx_file,software,phi_file=NULL){
  if(software == "Monolix"){
    parm_file <- tools::file_path_sans_ext(pmx_file) %>%
      paste0("/populationParameters.txt")
    pparms <- parm_extractor_mlx(parm_file) %>%
      dplyr::select(parameter,value)
    
    iparm_file <- tools::file_path_sans_ext(pmx_file) %>%
      paste0("/IndividualParameters/estimatedIndividualParameters.txt")
    iparms <- iparm_extractor_mlx(iparm_file)
  } else if(software == "NONMEM"){
    pparms <- parm_extractor_nm(pmx_file)
    
    etas <- eta_extractor_nm(phi_file)
    thetas <- theta_extractor_nm(pmx_file) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      `[`(rep(1,length(unique(etas$id))),) %>%
      as.data.frame() %>%
      dplyr::mutate(id = etas$id) %>%
      dplyr::relocate(id)
    iparms <- dplyr::left_join(thetas,etas,by="id")
  } else if(software == "nlmixr2"){
    err_parms <- nlmixr_err_identifier(pmx_file)
    
    pparms <- parm_extractor_nlmixr(pmx_file,err_parms)
    
    etas <- eta_extractor_nlmixr(pmx_file)
    thetas <- theta_extractor_nlmixr(pmx_file,err_parms) %>%
      tibble::column_to_rownames("parameter") %>%
      `colnames<-`(NULL) %>%
      t() %>%
      `[`(rep(1,length(unique(etas$id))),) %>%
      as.data.frame() %>%
      dplyr::mutate(id = etas$id) %>%
      dplyr::relocate(id)
    iparms <- dplyr::left_join(thetas,etas,by="id")
  }
  
  out <- list("PParms" = pparms,
              "IParms" = iparms)
  return(out)
}

#' Population estimates extractor for Monolix
#' 
#' Extracts population estimates from a Monolix run.
#' 
#' NULL
#' 
#' @param parm_file File of estimated population parameters.
#' @author Dominic Bräm
parm_extractor_mlx <- function(parm_file){
  parms <- read.table(parm_file,sep=",",header=T)
  return(parms)
}

#' Individual estimates extractor for Monolix
#' 
#' Extracts individual estimates (EBEs) from a Monolix run.
#' 
#' NULL
#' 
#' @param iparm_file File of estimated individual parameters.
#' @author Dominic Bräm
iparm_extractor_mlx <- function(iparm_file){
  ind_file_path <- iparm_file
  
  if(!file.exists(ind_file_path)){
    stop("No individual estimates available for provided Monolix file")
  }
  
  parm_table <- read.table(ind_file_path,header=T,sep=",")
  
  ind_parm_table <- parm_table[,grepl("id|_mode",colnames(parm_table))]
  colnames(ind_parm_table) <- gsub("_mode","",names(ind_parm_table))
  
  return(ind_parm_table)
}

#' Population estimates extractor for NONMEM
#' 
#' Extracts population estimates from a NONMEM results file.
#' 
#' NULL
#' 
#' @param res_file NONMEM results file.
#' @author Dominic Bräm
parm_extractor_nm <- function(res_file){
  parms <- rbind(theta_extractor_nm(res_file),
                 omega_extractor_nm(res_file))
  return(parms)
}

#' Theta estimates extractor for NONMEM
#' 
#' Extracts Theta estimates from a NONMEM results file.
#' 
#' NULL
#' 
#' @param res_path NONMEM results file.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
theta_extractor_nm <- function(res_path){
  res_file <- readLines(res_path)
  
  theta_est_start <- grep("THETA - VECTOR OF FIXED EFFECTS PARAMETERS",res_file)
  omega_est_start <- grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS",res_file)
  
  theta_est_lines <- res_file[(theta_est_start+1):(omega_est_start-1)]
  
  theta_est_names_match <- gregexpr("TH\\s*\\d+",theta_est_lines)
  theta_est_names <- unlist(regmatches(theta_est_lines,theta_est_names_match)) %>%
    {gsub("TH","THETA",.)} %>% 
    {gsub(" *","",.)}
  theta_est_values_match <- gregexpr("-?\\d+\\.?\\d+E[+-]?\\d+",theta_est_lines)
  theta_est_values <- unlist(regmatches(theta_est_lines,theta_est_values_match))
  
  df <- data.frame(parameter = theta_est_names,
                   value = as.numeric(theta_est_values))
  rownames(df) <- NULL
  return(df)
  
}

#' Omega estimates extractor for NONMEM
#' 
#' Extracts Omega estimates from a NONMEM results file.
#' 
#' NULL
#' 
#' @param res_path NONMEM results file.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
omega_extractor_nm <- function(res_path){
  res_file <- readLines(res_path)
  
  omega_est_start <- grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS",res_file)
  sigma_est_start <- grep("SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS",res_file)
  
  omega_est_lines <- res_file[(omega_est_start+1):(sigma_est_start-1)] %>% 
    {gsub(" ","",.)} %>% 
    `[`(!grepl("ET",.))
  
  omega_names_lines <- res_file[(omega_est_start+1):(sigma_est_start-1)] %>%
    `[`(1:(min(grep("\\+",.))-2))
  
  omega_est_values <- which(substr(omega_est_lines,1,1)=="+") %>% 
    lapply(function(x) omega_est_lines[x:length(omega_est_lines)]) %>% 
    lapply(function(x) {
      empties <- which(lapply(x,function(y) nchar(y)) %>% unlist() == 0)
      min_empty <- min(empties)
      return(x[1:(min_empty-1)])
    }) %>%
    lapply(function(x) paste(x,collapse = "")) %>%
    lapply(function(x) substr(x,(nchar(x)-7),nchar(x))) %>%
    unlist() %>%
    as.numeric()
  
  omega_est_names_match <- gregexpr("ET[A]*\\s*\\d+",omega_names_lines)
  omega_est_names <- unlist(regmatches(omega_names_lines,omega_est_names_match)) %>%
    {gsub("A","",.)} %>%
    {gsub("ET","ETA",.)}
  
  df <- data.frame(parameter = omega_est_names,
                   value = omega_est_values)
  rownames(df) <- NULL
  return(df)
  
}

#' Individual estimates extractor for NONMEM
#' 
#' Extracts individual Etas from a NONMEM .phi file.
#' 
#' NULL
#' 
#' @param phi_file NONMEM .phi file
#' @author Dominic Bräm
#' @importFrom magrittr %>%
eta_extractor_nm <- function(phi_file){
  if(!file.exists(phi_file)){
    stop("Provided .phi file from NONMEM run doesn't exist")
  }
  
  if(tools::file_ext(phi_file) != "phi"){
    stop("Please provide a .phi file from NONMEM run")
  }
  
  phi_text <- readLines(phi_file)[-1]
  phi_standardised <- gsub("\\s+"," ",phi_text)
  phi_split <- lapply(phi_standardised,function(x) unlist(strsplit(x,split=" ")))
  
  id_place <- grep("ID",phi_split[[1]])
  etas_place <- grep("ETA",phi_split[[1]])
  etas_names <- phi_split[[1]][etas_place] %>%
    {gsub("\\(|\\)","",.)}
  
  ids <- unlist(lapply(phi_split[-1],function(x) x[id_place]))
  etas <- lapply(phi_split[-1],function(x) as.numeric(x[etas_place])) %>%
    do.call(rbind,.) %>%
    `colnames<-`(etas_names)
  
  
  
  id_df <- data.frame(id = ids)
  out <- cbind(id_df,etas)
  
  return(out)
}

#' Population estimates extractor for nlmixr2
#' 
#' Extracts population estimates from a nlmixr2 run.
#' 
#' NULL
#' 
#' @param fit_obj Fitted object with nlmixr2.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
parm_extractor_nlmixr <- function(fit_obj,err_parms){
  thetas <- fit_obj$fixef %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column() %>% 
    `colnames<-`(c("parameter","value"))
  etas <- fit_obj$omega %>%
    diag()%>% 
    as.data.frame() %>% 
    tibble::rownames_to_column() %>% 
    `colnames<-`(c("parameter","value"))
  
  parms <- rbind(thetas,etas) %>%
    dplyr::filter(!(parameter %in% err_parms))
    
  return(parms)
}

#' Fixed effects extractor for nlmixr2
#' 
#' Extract estimated fixed effects from nlmixr2 fit.
#' 
#' NULL
#' 
#' @param fit_obj Fitted object with nlmixr2.
#' @param err_parms Names of residual error parameters in the model, extracted
#' through \link[=func]{nlmixr_err_identifier} to be excluded from fixed effects.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
theta_extractor_nlmixr <- function(fit_obj,err_parms){
  thetas <- fit_obj$fixef %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column() %>% 
    `colnames<-`(c("parameter","value")) %>%
    dplyr::filter(!(parameter %in% err_parms))
  
  return(thetas)
}

#' Etas extractor for nlmixr2
#' 
#' Extract estimated individual Etas from nlmixr2 fit.
#' 
#' NULL
#' 
#' @param fit_obj Fitted object with nlmixr2.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
eta_extractor_nlmixr <- function(fit_obj){
  ranefs <- coef(fit_obj)$random %>%
    dplyr::rename("id" = "ID")
  
  return(ranefs)
}
