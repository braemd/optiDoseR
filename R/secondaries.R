#' Generate Monolix code for secondary target
#' 
#' Generate Monolix code for secondary target, called by \link[=func]{addSecondary} if
#' software is "Monolix".
#' 
#' Penalization function is just \emph{state} if high state-values are penalized, and
#' exp(-0.01 * state) if low state-values are penalized.
#' 
#' @param state State in the structural model to be penalized.
#' @param time Time-point of penalization (one numeric value).
#' @param pen_values What values should be penalized; "high" if high state.values and "low" if
#' low state-values should be penalized.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @param decay_high Rate of penalization-decay with increasing state-values.
#' @author Dominic Bräm
mlxSecondary <- function(state,
                         time,
                         pen_values,
                         n_constraint,
                         gamma,
                         decay_high=0.01){
  
  pen_name <- paste0("pen",n_constraint)
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * ",
                 ifelse(pen_values == "high",
                        state,
                        paste0(
                          "exp(-",
                          decay_high,
                          " * ",
                          state,
                          ")")))
  
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "Secondary")
  
  return(out)
}

#' Generate Monolix code for secondary target AUC
#' 
#' Generate Monolix code for secondary target, called by \link[=func]{addSecondaryAUC} if
#' software is "Monolix".
#' 
#' Penalization function is just \emph{AUC_state} if high AUCs are penalized, and
#' exp(-0.01 * AUC_state) if low AUCs are penalized.
#' 
#' @param state State in the structural model to be penalized.
#' @param pen_values What values should be penalized; "high" if high AUCs and "low" if
#' low AUCs should be penalized.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @param decay_high Rate of penalization-decay with increasing state-values.
#' @author Dominic Bräm
mlxSecondaryAUC <- function(state,
                         pen_values,
                         pen_time,
                         eval_time,
                         n_constraint,
                         gamma,
                         decay_high=0.01){
  
  pen_name <- paste0("pen",n_constraint)
  
  ini_state <- paste0("AUC",
                      n_constraint,
                      "_0 = 0")
  
  auc_code <- list(paste0("d",
                          pen_name,
                          " = 0"),
                   paste0("if (t>",
                          pen_time[1],
                          " & t<",
                          pen_time[2],
                          ")"),
                   paste0("   d",
                          pen_name,
                          " = ",
                          state),
                   "end",
                   paste0("ddt_AUC",
                          n_constraint,
                          " = d",
                          pen_name)
  )
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * ",
                 ifelse(pen_values == "high",
                        paste0("AUC",
                               n_constraint),
                        paste0(
                          "exp(-",
                          decay_high,
                          " * ",
                          "AUC",
                          n_constraint,
                          ")")))
  
  
  out <- list("Name" = pen_name,
              "Time" = eval_time,
              "AUC_Code" = auc_code,
              "Code" = code,
              "Ini" = ini_state,
              "Type" = "SecondaryAUC")
  
  return(out)
}

#' Generate NONMEM code for secondary target
#' 
#' Generate NONMEM code for secondary target, called by \link[=func]{addSecondary} if
#' software is "NONMEM".
#' 
#' Penalization function is just \emph{state} if high state-values are penalized, and
#' exp(-0.01 * state) if low state-values are penalized.
#' 
#' @param state State in the structural model to be penalized.
#' @param time Time-point of penalization (one numeric value).
#' @param pen_values What values should be penalized; "high" if high state.values and "low" if
#' low state-values should be penalized.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @param decay_high Rate of penalization-decay with increasing state-values.
#' @author Dominic Bräm
nmSecondary <- function(state,
                         time,
                         pen_values,
                         n_constraint,
                         gamma,
                         decay_high=0.01){
  
  pen_name <- paste0("PEN",n_constraint)
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * ",
                 ifelse(pen_values == "high",
                        state,
                        paste0(
                          "2.7182**(-",
                          decay_high,
                          " * ",
                          state,
                          ")")))
  
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "Secondary")
  
  return(out)
}

#' Generate NONMEM code for secondary target AUC
#' 
#' Generate NONMEM code for secondary target, called by \link[=func]{addSecondaryAUC} if
#' software is "NONMEM".
#' 
#' Penalization function is just \emph{AUC_state} if high AUCs are penalized, and
#' exp(-0.01 * AUC_state) if low AUCs are penalized.
#' 
#' @param state State in the structural model for which AUC is calculated to be penalized.
#' @param pen_values What values should be penalized; "high" if high AUCs and "low" if
#' low AUCs should be penalized.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param n_state The index of the state for AUC (number of states in the structural model plus 1).
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @param decay_high Rate of penalization-decay with increasing state-values.
#' @author Dominic Bräm
nmSecondaryAUC <- function(state,
                           pen_values,
                           pen_time,
                           eval_time,
                           n_state,
                           n_constraint,
                           gamma,
                           decay_high=0.01){
  
  pen_name <- paste0("PEN",n_constraint)
  
  ini_state <- paste0("A_0(",n_state,") = 0")
  
  auc_code <- list(paste0("DAUC",
                          n_state,
                          " = 0"),
                   paste0("IF(TIME.GT.",
                          pen_time[1],
                          ".AND.TIME.LT.",
                          pen_time[2],
                          ") THEN"
                   ),
                   paste0("   DAUC",
                          n_state,
                          " = ",
                          state
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
                 " * ",
                 ifelse(pen_values == "high",
                        paste0("A(",
                               n_state,
                               ")"),
                        paste0(
                          "2.7182**(-",
                          decay_high,
                          " * ",
                          paste0("A(",
                                 n_state,
                                 ")"),
                          ")")))
  
  
  out <- list("Name" = pen_name,
              "Time" = eval_time,
              "AUC_Code" = auc_code,
              "Code" = code,
              "Ini" = ini_state,
              "Type" = "SecondaryAUC")
  
  return(out)
}

#' Generate nlmixr2 code for secondary target
#' 
#' Generate nlmixr2 code for secondary target, called by \link[=func]{addSecondary} if
#' software is "nlmixr2".
#' 
#' Penalization function is just \emph{state} if high state-values are penalized, and
#' exp(-0.01 * state) if low state-values are penalized.
#' 
#' @param state State in the structural model to be penalized.
#' @param time Time-point of penalization (one numeric value).
#' @param pen_values What values should be penalized; "high" if high state.values and "low" if
#' low state-values should be penalized.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @param decay_high Rate of penalization-decay with increasing state-values.
#' @author Dominic Bräm
nlmixrSecondary <- function(state,
                        time,
                        pen_values,
                        n_constraint,
                        gamma,
                        decay_high=0.01){
  
  pen_name <- paste0("PEN",n_constraint)
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * ",
                 ifelse(pen_values == "high",
                        state,
                        paste0(
                          "2.7182^(-",
                          decay_high,
                          " * ",
                          state,
                          ")")))
  
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "Secondary")
  
  return(out)
}

#' Generate nlmixr2 code for secondary target AUC
#' 
#' Generate nlmixr2 code for secondary target, called by \link[=func]{addSecondaryAUC} if
#' software is "nlmixr2".
#' 
#' Penalization function is just \emph{AUC_state} if high AUCs are penalized, and
#' exp(-0.01 * AUC_state) if low AUCs are penalized.
#' 
#' @param state State in the structural model for which AUC is calculated to be penalized.
#' @param pen_values What values should be penalized; "high" if high AUCs and "low" if
#' low AUCs should be penalized.
#' @param pen_time Time-range when penalization should be applied (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @param decay_high Rate of penalization-decay with increasing state-values.
#' @author Dominic Bräm
nlmixrSecondaryAUC <- function(state,
                           pen_values,
                           pen_time,
                           eval_time,
                           n_constraint,
                           gamma,
                           decay_high=0.01){
  
  pen_name <- paste0("PEN",n_constraint)
  
  ini_state <- paste0("AUC",n_constraint,"(0) = 0")
  
  auc_code <- list(paste0("dAUC",
                          n_constraint,
                          " = 0"),
                   paste0("if(t > ",
                          pen_time[1],
                          " & t < ",
                          pen_time[2],
                          "){"
                   ),
                   paste0("   dAUC",
                          n_constraint,
                          " = ",
                          state
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
                 " * ",
                 ifelse(pen_values == "high",
                        paste0("AUC",
                               n_constraint),
                        paste0(
                          "2.7182^(-",
                          decay_high,
                          " * ",
                          paste0("AUC",
                                 n_constraint),
                          ")")))
  
  
  out <- list("Name" = pen_name,
              "Time" = eval_time,
              "AUC_Code" = auc_code,
              "Code" = code,
              "Ini" = ini_state,
              "Type" = "SecondaryAUC")
  
  return(out)
}
