#' Add constraint AUC to optiProject
#' 
#' Add a constraint, i.e., penalizing state-values below or above a certain limit, 
#' over a time-frame.
#' 
#' The area between limit and state is calculated if state-values are above or below 
#' the limit, within the defined penalization time-frame.
#' 
#' @param optiProject A optiProject-object
#' @param state State in the structural model to be penalized.
#' @param limit The limit, under/above which the state gets penalized.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param pen_values What values should be penalized; "below" if limit is the lower and "above" if
#' limit is the upper limit.
#' @param gamma Penalization strength.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addConstraint("R",limit=30,pen_time = c(0,200),
#'                                 eval_time=200,pen_values="below",gamma=1000)
#' }
#' @import checkmate
#' @export
addConstraintAUC <- function(optiProject,
                             state,
                             limit,
                             pen_time,
                             eval_time,
                             pen_values,
                             gamma = 100){
  
  checkmate::assertClass(optiProject,"optiProject")
  
  checkmate::assertNumeric(pen_time,len=2)
  checkmate::assertNumeric(eval_time,len=1)
  checkmate::assertChoice(pen_values,c("below","above"))
  
  software <- optiProject$Software
  n_constraint <- length(optiProject$Constraints)
  out <- optiProject
  
  if(software == "Monolix"){
    out$Constraints[[n_constraint+1]] <- mlxConstraintAUC(state,limit,pen_time,eval_time,
                                                          pen_values,n_constraint,gamma)
  } else if(software == "NONMEM"){
    model_text <- optiProject$Model
    if(is.null(model_text)){
      stop("To add AUC constraint to NONMEM, first a model code must be added to the optiProject")
    }
    n_state <- nm_nstate_extractor(model_text) + sum(unlist(lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)))) + 1
    out$Constraints[[n_constraint+1]] <- nmConstraintAUC(state,limit,pen_time,eval_time,
                                                         pen_values,n_state,n_constraint,gamma)
  } else if(software == "nlmixr2"){
    out$Constraints[[n_constraint+1]] <- nlmixrConstraintAUC(state,limit,pen_time,eval_time,
                                                          pen_values,n_constraint,gamma)
  }
  
  return(out)
}

#' Add AUC of divergence from target value to optiProject
#' 
#' Add a target value, i.e., penalizing AUC from state-values different from this value, 
#' over a time-frame.
#' 
#' The area between state-curve and target value is calculated and penalized.
#' 
#' @param optiProject A optiProject-object
#' @param state State in the structural model to be penalized.
#' @param target_value Target value of the state.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param gamma Penalization strength.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addTargetValueAUC("R",target_value = 20, pen_time = c(30,50),
#'                                     eval_time = 50)
#' }
#' @import checkmate
#' @export
addTargetValueAUC <- function(optiProject,
                              state,
                              target_value,
                              pen_time,
                              eval_time,
                              gamma = 100){
  
  checkmate::assertClass(optiProject,"optiProject")
  
  checkmate::assertNumeric(pen_time,len=2)
  checkmate::assertNumeric(eval_time,len=1)
  
  software <- optiProject$Software
  n_constraint <- length(optiProject$Constraints)
  out <- optiProject
  
  if(software == "Monolix"){
    out$Constraints[[n_constraint+1]] <- mlxTargetValueAUC(state,target_value,pen_time,eval_time,
                                                           n_constraint,gamma)
  } else if(software == "NONMEM"){
    model_text <- optiProject$Model
    if(is.null(model_text)){
      stop("To add AUC constraint to NONMEM, first a model code must be added to the optiProject")
    }
    n_state <- nm_nstate_extractor(model_text) + sum(unlist(lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)))) + 1
    out$Constraints[[n_constraint+1]] <- nmTargetValueAUC(state,target_value,pen_time,eval_time,
                                                          n_state,n_constraint,gamma)
  } else if(software == "nlmixr2"){
    out$Constraints[[n_constraint+1]] <- nlmixrTargetValueAUC(state,target_value,pen_time,eval_time,
                                                           n_constraint,gamma)
  }
  
  return(out)
}

#' Add AUC of divergence from target function to optiProject
#' 
#' Add a target function, i.e., penalizing AUC from state-values different from this function, 
#' over a time-frame.
#' 
#' The area between state-curve and target function is calculated and penalized.
#' 
#' @param optiProject A optiProject-object
#' @param state State in the structural model to be penalized.
#' @param target_fun Target function for the state.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param gamma Penalization strength.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addTargetFunAUC("R",target_fun = "100 * exp(-0.2 * t)", pen_time = c(0,30),
#'                                     eval_time = 50)
#' }
#' @import checkmate
#' @export
addTargetFunAUC <- function(optiProject,
                            state,
                            target_fun,
                            pen_time,
                            eval_time,
                            gamma=1){
  checkmate::assertClass(optiProject,"optiProject")
  
  checkmate::assertString(target_fun)
  
  software <- optiProject$Software
  n_constraint <- length(optiProject$Constraints)
  out <- optiProject
  
  if(software == "Monolix"){
    out$Constraints[[n_constraint+1]] <- mlxTargetFunAUC(state,target_fun,pen_time,eval_time,
                                                         n_constraint,gamma)
  } else if(software == "NONMEM"){
    model_text <- optiProject$Model
    if(is.null(model_text)){
      stop("To add AUC constraint to NONMEM, first a model code must be added to the optiProject")
    }
    n_state <- nm_nstate_extractor(model_text) + sum(unlist(lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)))) + 1
    out$Constraints[[n_constraint+1]] <- nmTargetFunAUC(state,target_fun,pen_time,eval_time,
                                                        n_state,n_constraint,gamma)
  } else if(software == "nlmixr2"){
    out$Constraints[[n_constraint+1]] <- nlmixrTargetFunAUC(state,target_fun,pen_time,eval_time,
                                                         n_constraint,gamma)
  }
  
  return(out)
  
}

#' Add secondary target AUC to optiProject
#' 
#' Add a secondary target, i.e., penalizing low or high state-AUC, 
#' over a time-frame.
#' 
#' The area under the state-curve is calculated and either high or low AUCs are penalized.
#' 
#' @param optiProject A optiProject-object
#' @param state State in the structural model to be penalized.
#' @param pen_values What values should be penalized; "high" if high state.values and "low" if
#' low state-values should be penalized.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param gamma Penalization strength.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addConstraintAUC("R",limit = 30,pen_time = c(0,50), eval_time = 50, pen_values="below") %>%
#'                   addSecondaryAUC("Cc",pen_values = "low",pen_time = c(0,50), eval_time = 50)
#' }
#' @import checkmate
#' @export
addSecondaryAUC <- function(optiProject,
                            state,
                            pen_values,
                            pen_time,
                            eval_time,
                            gamma = 1){
  
  checkmate::assertClass(optiProject,"optiProject")
  
  checkmate::assertNumeric(pen_time,len=2)
  checkmate::assertNumeric(eval_time,len=1)
  checkmate::assertChoice(pen_values,c("high","low"))
  
  software <- optiProject$Software
  n_constraint <- length(optiProject$Constraints)
  out <- optiProject
  
  if(software == "Monolix"){
    out$Constraints[[n_constraint+1]] <- mlxSecondaryAUC(state,pen_values,pen_time,eval_time,
                                                          n_constraint,gamma)
  } else if(software == "NONMEM"){
    model_text <- optiProject$Model
    if(is.null(model_text)){
      stop("To add AUC constraint to NONMEM, first a model code must be added to the optiProject")
    }
    n_state <- nm_nstate_extractor(model_text) + sum(unlist(lapply(optiProject$Constraints,function(x) grepl("AUC",x$Type)))) + 1
    out$Constraints[[n_constraint+1]] <- nmSecondaryAUC(state,pen_values,pen_time,eval_time,
                                                         n_state,n_constraint,gamma)
  } else if(software == "nlmixr2"){
    out$Constraints[[n_constraint+1]] <- nlmixrSecondaryAUC(state,pen_values,pen_time,eval_time,
                                                         n_constraint,gamma)
  }
  
  return(out)
}
