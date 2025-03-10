#' Create new optiProject
#' 
#' Initializes a new optiProject-object for a software to which a model, parameters,
#' dose levels, and constraints/targets can be added.
#' 
#' NULL
#' 
#' @param software Software for which the optiProject should be initialized (currently
#' "Monolix" and "NONMEM" available)
#' @author Dominic Bräm
#' @returns A optiProject-object, to which additional elements of the optiDose approach
#' can be added.
#' @importFrom checkmate assertChoice
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix")
#' }
#' @export
newOptiProject <- function(software){
  checkmate::assertChoice(software,c("Monolix","NONMEM","nlmixr2"))
  
  out <- list("Model" = NULL,
              "Parms" = NULL,
              "IParms" = NULL,
              "Covs" = NULL,
              "Software" = software,
              "Doses" = list(),
              "Constraints" = list())
  
  attr(out, "class") <- "optiProject"
  return(out)
}

#' Adding an already fitted pharmacometric model
#' 
#' Adding an already fitted pharmacometric model to the optiProject-object, including
#' structural model, population and individual parameter estimates.
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object
#' @param pmx_file An already fitted pharmacometric model, i.e., a \emph{.mlxtran} file
#' if optiProject-object is initialized with "Monolix", or a NONMEM control file if initialized
#' with "NONMEM".
#' @param lixoft_path If initialized with "Monolix" and the structural model used in the \emph{.mlxtran} file is
#' a model from the library, \code{lixoft_path} is needed. Standard is "C:/ProgramData/Lixoft/MonolixSuiteXXXX".
#' @param res_file A NONMEM results file if optiProject is initialized with "NONMEM" (in order
#' to extract population estimates).
#' @param phi_file A NONMEM .phi file if optiProject is initialized with "NONMEM" (in order
#' to extract individual estimates).
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' 
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2")
#' }
#' @importFrom checkmate assertClass
#' @export
addPmxRun <- function(optiProject,
                      pmx_file,
                      lixoft_path = NULL,
                      res_file = NULL,
                      phi_file = NULL){
  
  checkmate::assertClass(optiProject,"optiProject")
  
  software <- optiProject$Software
  if(software == "Monolix"){
    model_info <- newOptiMlx(pmx_file,lixoft_path)
    ext_parms <- parm_extractor(pmx_file,software)
  } else if(software == "NONMEM"){
    if(is.null(res_file) | is.null(phi_file)){
      stop("Please provide a results file for NONMEM")
    }
    model_info <- list(model_text=readLines(pmx_file))
    ext_parms <- parm_extractor(pmx_file=res_file,software=software,phi_file=phi_file)
  } else if(software == "nlmixr2"){
    model_info <- newOptiNlmixr(pmx_file)
    ext_parms <- parm_extractor(pmx_file,software)
  }
  
  out <- optiProject
  out$Model <- model_info$model_text
  out$Parms <- ext_parms$PParms
  out$IParms <- ext_parms$IParms
  
  return(out)
}

#' Import structural model
#' 
#' Import the structural model used in a \emph{.mlxtran} run.
#' 
#' NULL
#' 
#' @param pmx_file A Monolix \emph{.mlxtran} file.
#' @param lixoft_path If initialized with "Monolix" and the structural model used in the \emph{.mlxtran} file is
#' a model from the library, \code{lixoft_path} is needed. Standard is "C:/ProgramData/Lixoft/MonolixSuiteXXXX".
#' @author Dominic Bräm
#' @importFrom magrittr %>%
#' @importFrom tools file_path_as_absolute
newOptiMlx <- function(pmx_file,lixoft_path){
  model_file <- readLines(pmx_file) %>%
    `[`(grep("\\[LONGITUDINAL\\]",.):length(.)) %>%
    `[`(1:grep("DEFINITION",.)) %>%
    `[`(grep("file",.)) %>%
    unlist() %>%
    {gsub("file *= *","",.)} %>%
    {gsub("\\'","",.)}
  
  if(grepl("lib:",model_file)){
    if(is.null(lixoft_path)){
      stop("Library model was used. Please provide path to Lixoft installation")
    }
    
    all_lib_path <- lixoft_path %>%
      paste0("/factory/library")
    for(i in dir(all_lib_path)){
      if(gsub("lib:","",model_file) %in% dir(paste0(all_lib_path,"/",i))){
        lib_folder <- i
      }
    }
    if(!exists("lib_folder")){
      stop("Model file could not be found in library")
    } else{
      model_file <- model_file %>%
        {gsub("lib:",paste0(all_lib_path,"/",lib_folder,"/"),.)}
    }
  } else{
    pmx_file_abs <- tools::file_path_as_absolute(dirname(pmx_file))
    model_file <- paste0(pmx_file_abs,"/",model_file)
  }
  
  model_text <- readLines(model_file)
  
  out <- list("model_file" = model_file,
              "model_text" = model_text)
  return(out)
}

#' Import structural model from nlmixr2
#' 
#' Import the structural model used in a nlmixr2-model run.
#' 
#' NULL
#' 
#' @param pmx_file A fitted nlmixr2 object.
#' @author Dominic Bräm
#' @importFrom magrittr %>%
newOptiNlmixr <- function(pmx_file){
  model_name <- pmx_file$modelName
  
  if(!exists(model_name)){
    stop("The model function in the nlmixr2 fit must be R global environment.")
  }
  
  model_info <- list(model_text=deparse(get(model_name)))
  
  return(model_info)
}

#' Add structural model to optiProject-object
#' 
#' If not an already fitted pharmacometric model is added through addPmxRun, an explicit structural
#' model can be added to the optiProject-object with this function.
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object
#' @param pmx_model (Path/)File name of structural model or rxode2 style model function 
#' to be added to the optiProject-object.If software is "Monolix", a \emph{.txt} file is 
#' required (as in a Monolix project), if software is "NONMEM", the NONMEM control file 
#' name can be added, if software is "nlmixr2", a rxode2 style model function that would be
#' used in the nlmixr2 fit can be given.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxModel("../inst/oral1_1cpt_IndirectModelInhibitionKin_TlagkaVClR0koutImaxIC50.txt")
#' }
#' @importFrom checkmate assertClass
#' @importFrom tools file_ext
#' @export
addPmxModel <- function(optiProject,
                        pmx_model){
  
  checkmate::assertClass(optiProject,"optiProject")
  
  software <- optiProject$Software
  
  if(software %in% c("Monolix","NONMEM")){
    if(software == "Monolix" & tools::file_ext(pmx_model) != "txt"){
      stop("The PMX model for Monolix must be a .txt file")
    }
    
    model_text <- readLines(pmx_model)
  } else if(software == "nlmixr2"){
    model_text <- deparse(pmx_model)
  }
  
  out <- optiProject
  out$Model <- model_text
  
  return(out)
}

#' Add population estimates to optiProject-object
#' 
#' If not an already fitted pharmacometric model is added through addPmxRun, explicit population
#' parameter estimates can be added to the optiProject-object with this function.
#' 
#' If software ist "NONMEM" and a data frame is given, the parameter names must be of type THETAX and ETAX with X the number 
#' in THETA(X) and ETA(X) in the control file. This function is not applicable for nlmixr2 as
#' a fitted nlmixr2 object can directly added with \link[=func]{addPmxRun}.
#' 
#' @param optiProject A optiProject-object
#' @param pmx_parms Either a file name (populationEstimates.txt for "Monolix" or results file for
#' "NONMEM") containing population estimates, or a data frame with a column "parameter" for parameter
#' names and "value" for the population estimates. 
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxModel("../inst/oral1_1cpt_IndirectModelInhibitionKin_TlagkaVClR0koutImaxIC50.txt") %>%
#'                   addPopParms("../inst/warfarinPKPD_project/populationParameters.txt")
#' }
#' @importFrom checkmate assertClass
#' @export
addPopParms <- function(optiProject,
                        pmx_parms){
  
  checkmate::assertClass(optiProject,"optiProject")
  if(is.null(optiProject$Model)){
    stop("Please add first a pharmacometric model before adding parameters")
  }
  
  out <- optiProject
  software <- optiProject$Software
  if(class(pmx_parms) == "data.frame"){
    df_cols <- colnames(pmx_parms)
    if(!all(df_cols == c("parameter","value"))){
      stop("The data frame must have only a 'parameter' and a 'value' column")
    }
    out$Parms <- pmx_parms
    return(out)
  } else if(class(pmx_parms) == "character"){
    if(!file.exists(pmx_parms)){
      stop("The file you provided does not exist")
    }
    if(software == "Monolix"){
      parms <- read.table(pmx_parms,header=T,sep=",")
      df_cols <- colnames(parms)
      if(!all(df_cols == c("parameter","value"))){
        stop("The data frame must have only a 'parameter' and a 'value' column")
      }
      out$Pars <- parms
      return(out)
    } else if(software == "NONMEM"){
      parms <- parm_extractor_nm(pmx_parms)
      out$Parms <- parms
      return(out)
    }
  }
}

#' Add individual estimates to optiProject-object
#' 
#' If not an already fitted pharmacometric model is added through addPmxRun, explicit individual
#' parameter estimates can be added to the optiProject-object with this function.
#' 
#' If software ist "NONMEM", the parameter names must be of type THETAX and ETAX with X the number 
#' in THETA(X) and ETA(X) in the control file. This function is not applicable for nlmixr2 as
#' a fitted nlmixr2 object can directly added with \link[=func]{addPmxRun}.
#' 
#' @param optiProject A optiProject-object
#' @param pmx_parms Either a file name (estimatedIndividualParameters.txt for "Monolix" or results file for
#' "NONMEM") containing population estimates, or a data frame with a column "id" for subject IDs and columns 
#' for each parameter (if "Monolix", the explicit individual parameter, if "NONMEM", THETAs and ETAs).
#' @param phi_file If software is "NONMEM" and pmx_parms not a data frame, a .phi file of the fitted model
#' must be given to extract individual ETA estimates.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxModel("../inst/oral1_1cpt_IndirectModelInhibitionKin_TlagkaVClR0koutImaxIC50.txt") %>%
#'                   addIndParms("../inst/warfarinPKPD_project/IndividualParameters/estimatedIndividualParameters.txt")
#' }
#' @importFrom checkmate assertClass
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr left_join
#' @importFrom tibble column_to_rownames
#' @export
addIndParms <- function(optiProject,
                        pmx_parms,
                        phi_file=NULL){
  checkmate::assertClass(optiProject,"optiProject")
  if(is.null(optiProject$Model)){
    stop("Please add first a pharmacometric model before adding parameters")
  }

  out <- optiProject
  software <- optiProject$Software
  
  if(class(pmx_parms) == "data.frame"){
    df_cols <- colnames(pmx_parms)
    if(software == "Monolix"){
      if(!any(grepl("input *= *\\{[^\\}]+\\}",optiProject$Model))){
        input_start <- grep("input *= *\\{",optiProject$Model)
        input_end <- min(grep("\\}",optiProject$Model))
        input_line <- paste(optiProject$Model[input_start:input_end],collapse = "")
      } else{
        input_line <- optiProject$Model[grep("input *= *\\{",optiProject$Model)]
      }
      parm_names <- input_line %>%
        {gsub("input *= *\\{|\\}","",.)} %>%
        strsplit(",") %>%
        unlist()
    } else if(software == "NONMEM"){
      parm_names <- regexpr("THETA\\(\\d+\\)|ETA\\(\\d+\\)",optiProject$Model) %>%
        {regmatches(optiProject$Model,.)} %>%
        {gsub("\\(|\\)","",.)}
    }
    if(!(all(c("id",parm_names) %in% df_cols) & all(df_cols %in% c("id",parm_names)))){
      stop("If a dataframe is provided, only columns for id and parameters are allowed and required")
    }
    out$IParms <- pmx_parms
    return(out)
  } else if(class(pmx_parms) == "character"){
    if(software == "Monolix"){
      iparms <- iparm_extractor_mlx(pmx_parms)
    } else if(software == "NONMEM"){
      if(is.null(phi_file)){
        stop("For adding individual parameters to NONMEM from a file, a .phi file is required")
      } else if(!file.exists(phi_file)){
        stop("Provided .phi file doesn't exist")
      }
      etas <- eta_extractor_nm(phi_file)
      thetas <- theta_extractor_nm(pmx_parms) %>%
        tibble::column_to_rownames("parameter") %>%
        `colnames<-`(NULL) %>%
        t() %>%
        `[`(rep(1,length(unique(etas$id))),) %>%
        as.data.frame() %>%
        mutate(id = etas$id) %>%
        relocate(id)
      iparms <- dplyr::left_join(thetas,etas,by="id")
    }
    
    out$IParms <- iparms
    return(out)
  }
   
}

#' Add dose level to optiProject-object
#' 
#' Add a dose level to be estimated to a optiProject-object
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object
#' @param time Time of first dose for this dose level.
#' @param n_dose Number of doses for this dose level.
#' @param ini_est Initial estimate for this dose level.
#' @param dose_interval If \code{n_dose>1}, a dosing interval is required.
#' @param adm If different administration IDs are present in the structural model,
#' the administration ID of this dose level must be defined.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxModel("../inst/oral1_1cpt_IndirectModelInhibitionKin_TlagkaVClR0koutImaxIC50.txt") %>%
#'                   addPopParms("../inst/warfarinPKPD_project/populationParameters.txt") %>%
#'                   addDoseLevel(time=0,ini_est = 100)
#' }
#' @importFrom checkmate assertClass
#' @importFrom checkmate assertNumeric
#' @export
addDoseLevel <- function(optiProject,
                         time,
                         n_dose = 1,
                         ini_est = 1,
                         dose_interval = NULL,
                         adm = NULL){
  
  checkmate::assertClass(optiProject,"optiProject")
  checkmate::assertNumeric(ini_est,lower = 0,len = 1)
  
  out <- optiProject
  
  if(n_dose > 1 & is.null(dose_interval)){
    stop("If more than one dose is given, please provide the dosing interval")
  }
  
  if(n_dose == 1){
    out$Doses[[length(out$Doses)+1]] <- list("Times" = time,
                                             "ADM" = adm,
                                             "Ini" = ini_est)
  } else if(n_dose > 1){
    out$Doses[[length(out$Doses)+1]] <- list("Times" = seq(time,time + (n_dose-1)*dose_interval, by=dose_interval),
                                             "ADM" = adm,
                                             "Ini" = ini_est)
  }
  
  return(out)
}
