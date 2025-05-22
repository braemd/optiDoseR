#' Generate Monolix code for target function penalization
#' 
#' Generate Monolix code for target function penalization, called by \link[=func]{addTargetFun} if
#' software is "Monolix".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param target_fun Target function for the state.
#' @param time Time-point of penalization (one numeric value).
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
mlxTargetFun <- function(state,
                           target_fun,
                           time,
                           n_constraint,
                           gamma){
  pen_name <- paste0("pen",n_constraint)
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * (",
                 state,
                 " - ",
                 target_fun,
                 ")^2")
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "TargetFun")
  
  return(out)
}

#' Generate Monolix code for target function AUC penalization
#' 
#' Generate Monolix code for target function AUC penalization, called by \link[=func]{addTargetFunAUC} if
#' software is "Monolix".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param target_fun Target function for the state.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
mlxTargetFunAUC <- function(state,
                         target_fun,
                         pen_time,
                         eval_time,
                         n_constraint,
                         gamma){
  pen_name <- paste0("pen",n_constraint)
  
  ini_state <- paste0("AUC",
                      n_constraint,
                      "_0 = 0")
  
  auc_code <- list(paste0("d",
                          pen_name,
                          " = 0"),
                   paste0("if (t>=",
                          pen_time[1],
                          " & t<=",
                          pen_time[2],
                          ")"),
                   paste0("   d",
                          pen_name,
                          " = (",
                          state,
                          " - ",
                          target_fun,
                          ")^2"),
                   "end",
                   paste0("ddt_AUC",
                          n_constraint,
                          " = d",
                          pen_name)
  )
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * AUC",
                 n_constraint)
  
  out <- list("Name" = pen_name,
              "Time" = eval_time,
              "AUC_Code" = auc_code,
              "Code" = code,
              "Ini" = ini_state,
              "Type" = "TargetFunAUC")
  
  return(out)
}

#' Generate NONMEM code for target function penalization
#' 
#' Generate NONMEM code for target function penalization, called by \link[=func]{addTargetFun} if
#' software is "NONMEM".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param target_fun Target function for the state.
#' @param time Time-point of penalization (one numeric value).
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
nmTargetFun <- function(state,
                          target_fun,
                          time,
                          n_constraint,
                          gamma){
  pen_name <- paste0("PEN",n_constraint)
  
  code <- list(paste0("REF_FUN",
                      n_constraint,
                      " = ",
                      target_fun),
               paste0(pen_name,
                      " = ",
                      gamma,
                      " * (",
                      state,
                      " - REF_FUN)**2"))
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "TargetFun")
  
  return(out)
}

#' Generate NONMEM code for target function AUC penalization
#' 
#' Generate NONMEM code for target function AUC penalization, called by \link[=func]{addTargetFunAUC} if
#' software is "NONMEM".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param target_fun Target function for the state.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param n_state The index of the state for AUC (number of states in the structural model plus 1).
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
nmTargetFunAUC <- function(state,
                        target_fun,
                        pen_time,
                        eval_time,
                        n_state,
                        n_constraint,
                        gamma){
  pen_name <- paste0("PEN",n_constraint)
  
  ini_state <- paste0("A_0(",n_state,") = 0")
  
  auc_code <- list(paste0("REF_FUN",
                          n_constraint,
                          " = ",
                          target_fun),
                   paste0("DAUC",
                          n_state,
                          " = 0"),
                   paste0("IF(TIME.GE.",
                          pen_time[1],
                          ".AND.TIME.LE.",
                          pen_time[2],
                          ") THEN"
                   ),
                   paste0("   DAUC",
                          n_state,
                          " = (",
                          state,
                          " - REF_FUN)**2"
                   ),
                   "ENDIF",
                   paste0("DADT(",
                          n_state,
                          ") = DAUC",
                          n_state)
  )
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * A(",
                 n_state,
                 ")")
  
  out <- list("Name" = pen_name,
              "Time" = eval_time,
              "AUC_Code" = auc_code,
              "Code" = code,
              "Ini" = ini_state,
              "Type" = "TargetFunAUC")
  
  return(out)
}

#' Generate nlmixr2 code for target function penalization
#' 
#' Generate nlmixr2 code for target function penalization, called by \link[=func]{addTargetFun} if
#' software is "nlmixr2".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param target_fun Target function for the state.
#' @param time Time-point of penalization (one numeric value).
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
nlmixrTargetFun <- function(state,
                        target_fun,
                        time,
                        n_constraint,
                        gamma){
  pen_name <- paste0("PEN",n_constraint)
  
  code <- list(paste0("REF_FUN",
                      n_constraint,
                      " = ",
                      target_fun),
               paste0(pen_name,
                      " = ",
                      gamma,
                      " * (",
                      state,
                      " - REF_FUN)^2"))
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "TargetFun")
  
  return(out)
}

#' Generate nlmixr2 code for target function AUC penalization
#' 
#' Generate nlmixr2 code for target function AUC penalization, called by \link[=func]{addTargetFunAUC} if
#' software is "nlmixr2".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param target_fun Target function for the state.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
nlmixrTargetFunAUC <- function(state,
                           target_fun,
                           pen_time,
                           eval_time,
                           n_constraint,
                           gamma){
  pen_name <- paste0("PEN",n_constraint)
  
  ini_state <- paste0("AUC",n_constraint,"(0) = 0")
  
  auc_code <- list(paste0("REF_FUN",
                          n_constraint,
                          " = ",
                          target_fun),
                   paste0("dAUC",
                          n_constraint,
                          " = 0"),
                   paste0("if(t >= ",
                          pen_time[1],
                          " & t <= ",
                          pen_time[2],
                          "){"
                   ),
                   paste0("   dAUC",
                          n_constraint,
                          " = (",
                          state,
                          " - REF_FUN)^2"
                   ),
                   "}",
                   paste0("d/dt(AUC",
                          n_constraint,
                          ") = dAUC",
                          n_constraint)
  )
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * AUC",
                 n_constraint)
  
  out <- list("Name" = pen_name,
              "Time" = eval_time,
              "AUC_Code" = auc_code,
              "Code" = code,
              "Ini" = ini_state,
              "Type" = "TargetFunAUC")
  
  return(out)
}
