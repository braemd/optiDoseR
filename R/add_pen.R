#' Add constraint to optiProject
#' 
#' Add a constraint, i.e., penalizing state-values below or above a certain limit, 
#' at a specific time-point.
#' 
#' Constraint at a specific time-points may be useful if only a specific time-point is of
#' interest, e.g., trough-concentrations. For broader application, see \link[=func]{addConstraintAUC}.
#' 
#' @param optiProject A optiProject-object
#' @param state State in the structural model to be penalized.
#' @param limit The limit, under/above which the state gets penalized.
#' @param time Time-point of penalization (one numeric value).
#' @param pen_values What values should be penalized; "below" if limit is the lower and "above" if
#' limit is the upper limit.
#' @param gamma Penalization strength.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addConstraint("R",limit=30,pen_time = 50,pen_values="below")
#' }
#' @import checkmate
#' @export
addConstraint <- function(optiProject,
                          state,
                          limit,
                          time,
                          pen_values,
                          gamma = 100){
  
  checkmate::assertClass(optiProject,"optiProject")
  
  checkmate::assertNumeric(time,len=1)
  checkmate::assertChoice(pen_values,c("below","above"))
  
  software <- optiProject$Software
  n_constraint <- length(optiProject$Constraints)
  out <- optiProject
  
  if(software == "Monolix"){
    out$Constraints[[n_constraint+1]] <- mlxConstraint(state,limit,time,pen_values,n_constraint,gamma)
  } else if(software == "NONMEM"){
    out$Constraints[[n_constraint+1]] <- nmConstraint(state,limit,time,pen_values,n_constraint,gamma)
  } else if(software == "nlmixr2"){
    out$Constraints[[n_constraint+1]] <- nlmixrConstraint(state,limit,time,pen_values,n_constraint,gamma)
  }
  
  return(out)
}

#' Add secondary target to optiProject
#' 
#' Add a secondary target, i.e., penalizing low or high state-values, 
#' at a specific time-point.
#' 
#' Usually, secondary targets are required when a constraint is applied. 
#' Secondary targets at a specific time-points may be useful if only a specific time-point is of
#' interest, e.g., trough-concentrations. For broader application, see \link[=func]{addSecondaryAUC}.
#' 
#' @param optiProject A optiProject-object
#' @param state State in the structural model to be penalized.
#' @param pen_values What values should be penalized; "high" if high state.values and "low" if
#' low state-values should be penalized.
#' @param time Time-point of penalization (one numeric value).
#' @param gamma Penalization strength.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addConstraint("R",limit = 30,time = 50,pen_values="below") %>%
#'                   addSecondary("Cc",pen_values = "low",time = 30)
#' }
#' @import checkmate
#' @export
addSecondary <- function(optiProject,
                         state,
                         pen_values,
                         time,
                         gamma=1){
  
  checkmate::assertClass(optiProject,"optiProject")
  
  checkmate::assertNumeric(time,len=1)
  checkmate::assertChoice(pen_values,c("high","low"))
  
  software <- optiProject$Software
  n_constraint <- length(optiProject$Constraints)
  out <- optiProject
  
  if(software == "Monolix"){
    out$Constraints[[n_constraint+1]] <- mlxSecondary(state,time,pen_values,n_constraint,gamma)
  } else if(software == "NONMEM"){
    out$Constraints[[n_constraint+1]] <- nmSecondary(state,time,pen_values,n_constraint,gamma)
  } else if(software == "nlmixr2"){
    out$Constraints[[n_constraint+1]] <- nlmixrSecondary(state,time,pen_values,n_constraint,gamma)
  }
  
  return(out)
}

#' Add a target value to optiProject
#' 
#' Add a target value, i.e., penalizing state-values different from this value, 
#' at a specific time-point.
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object
#' @param state State in the structural model to be penalized.
#' @param target_value Target value of the state.
#' @param time Time-point of penalization (one numeric value).
#' @param gamma Penalization strength.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addTargetValue("R",target_value = 20, time = 50)
#' }
#' @import checkmate
#' @export
addTargetValue <- function(optiProject,
                           state,
                           target_value,
                           time,
                           gamma=1){
  
  checkmate::assertClass(optiProject,"optiProject")
  
  checkmate::assertNumeric(target_value,len=1)
  checkmate::assertNumeric(time,len=1)
  
  software <- optiProject$Software
  n_constraint <- length(optiProject$Constraints)
  out <- optiProject
  
  if(software == "Monolix"){
    out$Constraints[[n_constraint+1]] <- mlxTargetValue(state,target_value,time,n_constraint,gamma)
  } else if(software == "NONMEM"){
    out$Constraints[[n_constraint+1]] <- nmTargetValue(state,target_value,time,n_constraint,gamma)
  } else if(software == "nlmixr2"){
    out$Constraints[[n_constraint+1]] <- nlmixrTargetValue(state,target_value,time,n_constraint,gamma)
  }
  
  return(out)
}

#' Add a target function to optiProject
#' 
#' Add a target function, i.e., penalizing state-values different from the target function values, 
#' at a specific time-point.
#' 
#' NULL
#' 
#' @param optiProject A optiProject-object
#' @param state State in the structural model to be penalized.
#' @param target_fun Target function of the state.
#' @param time Time-point of penalization (one numeric value).
#' @param gamma Penalization strength.
#' @author Dominic Bräm
#' @examples
#' \dontrun{
#' warf_proj <- newOptiProject("Monolix") %>%
#'                   addPmxRun("../inst/warfarinPKPD_project.mlxtran",
#'                             lixoft_path = "C:/ProgramData/Lixoft/MonolixSuite2021R2") %>%
#'                   addTargetFun("R",target_fun = "100 \* exp(-0.2 \* t)", time = 50)
#' }
#' @import checkmate
#' @export
addTargetFun <- function(optiProject,
                         state,
                         target_fun,
                         time,
                         gamma=1){
  checkmate::assertClass(optiProject,"optiProject")
  
  checkmate::assertString(target_fun)
  
  software <- optiProject$Software
  n_constraint <- length(optiProject$Constraints)
  out <- optiProject
  
  if(software == "Monolix"){
    out$Constraints[[n_constraint+1]] <- mlxTargetFun(state,target_fun,time,n_constraint,gamma)
  } else if(software == "NONMEM"){
    out$Constraints[[n_constraint+1]] <- nmTargetFun(state,target_fun,time,n_constraint,gamma)
  } else if(software == "nlmixr2"){
    out$Constraints[[n_constraint+1]] <- nlmixrTargetFun(state,target_fun,time,n_constraint,gamma)
  }
  
  return(out)
  
}