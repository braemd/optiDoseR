#' Realize and save your OptiProject
#' 
#' This function generates the required data, the model file from the OptiProject object,
#' and if it is a Monolix-OptiProject, the corresponding .mlxtran file
#' 
#' The optiProject needs to have added a \emph{Model}, at least one penalization (i.e., constraint, 
#' secondary target, target value or function), and at least one dose level. Additionally, either
#' population parameters or individual parameters need to be added.
#' If doses are estimated with inter-individual variability, the objective function consists of term
#' that tries to minimize the penalization and a term that aims at having normal distributed random effects.
#' The term that minimizes the penalization is scaled by the residual error term, thus by fixing a small
#' constant error, emphasize is given to minimizing the penalization. So if simulations with estimated
#' individual doses show that targets are not met, reducing \code{const_err} will result in better meeting
#' the targets.
#' 
#' @param optiProject A optiProject-object with a model, penalizations, dose levels, and parameters.
#' For requirements regarding penalizations, dose levels, and parameters, see details.
#' @param name Name under which optiProject will be saved (.mlxtran will be added for Monolix, and 
#' .ctl for NONMEM).
#' @param save_path Path, where optiProject should be saved.
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
#' @param const_err If doses are estimated with inter-individual variability for a population of individual
#' subjects, \code{const_err} defines the balance between fulfilling the objectives (constraints, target values etc.)
#' and having the doses log-normal distributed. For more detail, see \emph{Details}.
#' @param create_mlx If optiProject-software is "Monolix", whether directly a Monolix project (.mlxtran file) should
#' be created. Default is TRUE.
#' @param seperate If \code{pop=TRUE} and \code{iiv=FALSE}, whether individuals should be fitted separately (TRUE) or
#' all individuals should be fitted together.
#' @param run_nlmixr2 If optiProject-object is for nlmixr2, whether the generated nlmixr2 optiProject should directly be run.
#' @param control If optiProject-object if for nlmixr2 and \code{run_nlmixr2=TRUE}, list of control arguments given to nlmxir2. If \code{iiv=FALSE},
#' control for "bobyqa" and if \code{iiv=TRUE} control for "saem" must be given.
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addConstraintAUC("R",limit=30,pen_time = c(0,200),
#'                                   eval_time=200,pen_values="below",gamma=1000) %>%
#'                   addSecondaryAUC("Cc",pen_values = "low",pen_time = c(0,500),eval_time=500) %>%
#'                   addDoseLevel(time=0,ini_est = 100) %>%
#'                   addDoseLevel(time=100,ini_est = 100)
#'                   
#' saveOptiProject(warf_proj,"warf_optiDose",save_path = "~",
#'                 parametrisation = c("log","log","log","log","log","log","logit","log"),
#'                 lower_parm_limits = 0,
#'                 upper_parm_limits = 1)
#' }
#' @author Dominic Bräm
#' @export
saveOptiProject <- function(optiProject,
                            name,
                            save_path = NULL,
                            pop = FALSE,
                            iiv = FALSE,
                            sim = FALSE,
                            n_id = 1,
                            parametrisation = "log",
                            lower_parm_limits = NULL,
                            upper_parm_limits = NULL,
                            const_err = 0.0001,
                            create_mlx = TRUE,
                            separate = FALSE,
                            run_nlmixr2 = FALSE,
                            control = list()){
  
  checkmate::assertClass(optiProject,"optiProject")
  checkmate::assert(checkmate::checkTRUE(pop),checkmate::checkFALSE(pop))
  checkmate::assert(checkmate::checkTRUE(sim),checkmate::checkFALSE(sim))
  checkmate::assertNumeric(n_id,len=1,lower=1)
  checkmate::assertCharacter(parametrisation,pattern="log|normal|logit")
  checkmate::assert(checkmate::checkNull(save_path),checkmate::checkString(save_path))
  checkmate::assertString(name)
  checkmate::assert(checkmate::checkTRUE(create_mlx),checkmate::checkFALSE(create_mlx))
  checkmate::assert(checkmate::checkTRUE(separate),checkmate::checkFALSE(separate))
  checkmate::assert(checkmate::checkTRUE(run_nlmixr2),checkmate::checkFALSE(run_nlmixr2))
  
  if(is.null(optiProject$Model)){
    stop("Please add PMX model")
  }
  
  if(pop & sim & is.null(optiProject$Parms)){
    stop("For simulating individuals, please add population parameters")
  }
  
  if(pop & !sim & is.null(optiProject$IParms)){
    stop("If multiple subjects should be included, please add individual parameters or set sim=TRUE")
  } else if(pop & !sim & !is.null(optiProject$IParms)){
    if(nrow(optiProject$IParms)<=1){
      stop("IParm must include more than 1 subject.")
    }
  }
  
  if(iiv & !pop){
    stop("If IIV should be included, pop must be TRUE")
  }
  
  if(pop & sim & n_id<=1){
    stop("If sim=TRUE, n_id must be greater than 1.")
  }
  
  if(is.null(save_path)){
    save_path <- getwd()
  }
  
  software <- optiProject$Software
  if(software == "Monolix" & create_mlx){
    if(separate){
      warning("Option separate not available for Monolix, will be ignored")
    }
    
    if(!("lixoftConnectors" %in% row.names(installed.packages()))){
      stop("lixoftConnectors is needed for pmxOptiD in Monolix, please install this package first")
    }
    if(!("lixoftConnectors" %in% .packages())){
      stop("Please first load lixoftConnectors and initialize your MonolixSuite")
    }
    if(is.null(suppressMessages(lixoftConnectors::isProjectLoaded()))){
      stop("Please initialize lixoftConnectors first")
    }
    
    mlxSetAll(optiProject,
              name,
              save_path,
              pop,
              iiv,
              sim,
              n_id,
              parametrisation,
              lower_parm_limits = lower_parm_limits,
              upper_parm_limits = upper_parm_limits,
              const_err = const_err,
              create_mlx = create_mlx)
  } else if(software == "NONMEM"){
    if(iiv & separate){
      warning("iiv and separate was set to TRUE, model will be with iiv and not separate.")
    }
    
    nmSetAll(optiProject,
             name,
             save_path,
             pop,
             iiv,
             sim,
             n_id,
             separate)
  } else if(software == "nlmixr2"){
    if(run_nlmixr2){
      if(!("nlmixr2" %in% row.names(installed.packages()))){
        stop("Package nlmixr2 is needed to directly run the model, please install this package first")
      }
      if(!("nlmixr2" %in% .packages())){
        stop("Please first load the package nlmixr2")
      }
    }
    nlmixrSetAll(optiProject,
                 name,
                 pop,
                 iiv,
                 sim,
                 n_id,
                 const_err,
                 run_nlmixr2,
                 control)
  }
}
