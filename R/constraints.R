#' Generate Monolix code for constraint
#' 
#' Generate Monolix code for constraint, called by \link[=func]{addConstraint} if
#' software is "Monolix".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param limit The limit, under/above which the state gets penalized.
#' @param time Time-point of penalization (one numeric value).
#' @param pen_values What values should be penalized; "below" if limit is the lower and "above" if
#' limit is the upper limit.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
mlxConstraint <- function(state,
                          limit,
                          time,
                          pen_values,
                          n_constraint,
                          gamma){
  pen_name <- paste0("pen",n_constraint)
  
  code <- paste0(pen_name,
                 " = ",
                 gamma,
                 " * max(0,",
                 ifelse(pen_values=="above",
                        paste0(state," - ",limit),
                        paste0(limit," - ",state)),
                 ")")
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "Constraint")
  
  return(out)
}

#' Generate Monolix code for constraint AUC
#' 
#' Generate Monolix code for constraint, called by \link[=func]{addConstraintAUC} if
#' software is "Monolix".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param limit The limit, under/above which the AUC gets calculated for penalization.
#' @param pen_time Time-range of penalization (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param eval_time Time-point at which the constraint AUC should be evaluated.
#' @param pen_values What values should be penalized; "below" if limit is the lower and "above" if
#' limit is the upper limit.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
mlxConstraintAUC <- function(state,
                          limit,
                          pen_time,
                          eval_time,
                          pen_values,
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
                      " = max(0,",
                      ifelse(pen_values=="above",
                             paste0(state," - ",limit),
                             paste0(limit," - ",state)),
                      ")"),
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
              "Type" = "ConstraintAUC")
  
  return(out)
}

#' Generate NONMEM code for constraint
#' 
#' Generate NONMEM code for constraint, called by \link[=func]{addConstraint} if
#' software is "NONMEM".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param limit The limit, under/above which the state gets penalized.
#' @param time Time-point of penalization (one numeric value).
#' @param pen_values What values should be penalized; "below" if limit is the lower and "above" if
#' limit is the upper limit.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
nmConstraint <- function(state,
                         limit,
                         time,
                         pen_values,
                         n_constraint,
                         gamma){
  pen_name <- paste0("PEN",n_constraint)
  
  code <- list(
    paste0(pen_name,
           " = 0"),
    paste0("IF (",
           ifelse(pen_values=="above",
                  paste0(state," - ",limit),
                  paste0(limit," - ",state)),
           " > 0) ",
           pen_name,
           " = ",
           gamma,
           " * (",
           ifelse(pen_values=="above",
                  paste0(state," - ",limit),
                  paste0(limit," - ",state)),
           ")**2")
  )
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "Constraint")
  
  return(out)
  
}

#' Generate NONMEM code for constraint AUC
#' 
#' Generate NONMEM code for constraint, called by \link[=func]{addConstraintAUC} if
#' software is "NONMEM".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param limit The limit, under/above which the AUC gets calculated for penalization.
#' @param pen_time Time-range of penalization (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param pen_values What values should be penalized; "below" if limit is the lower and "above" if
#' limit is the upper limit.
#' @param n_state The index of the state for AUC (number of states in the structural model plus 1)
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
nmConstraintAUC <- function(state,
                         limit,
                         pen_time,
                         eval_time,
                         pen_values,
                         n_state,
                         n_constraint,
                         gamma){
  pen_name <- paste0("PEN",n_constraint)
  
  ini_state <- paste0("A_0(",n_state,") = 0")
  
  auc_code <- list(paste0("DAUC",
                      n_state,
                      " = 0"),
               paste0("D",
                      pen_name,
                      " = 0"),
               paste0("IF (",
                      ifelse(pen_values=="above",
                             paste0(state," - ",limit),
                             paste0(limit," - ",state)),
                      " > 0) D",
                      pen_name,
                      " = (",
                      ifelse(pen_values=="above",
                             paste0(state," - ",limit),
                             paste0(limit," - ",state)),
                      ")**2"),
               paste0("IF(TIME.GE.",
                      pen_time[1],
                      ".AND.TIME.LE.",
                      pen_time[2],
                      ") THEN"
                      ),
               paste0("   DAUC",
                      n_state,
                      " = D",
                      pen_name
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
              "Type" = "ConstraintAUC")
  
  return(out)
  
}

#' Generate nlmixr2 code for constraint
#' 
#' Generate nlmixr2 code for constraint, called by \link[=func]{addConstraint} if
#' software is "nlmixr2".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param limit The limit, under/above which the state gets penalized.
#' @param time Time-point of penalization (one numeric value).
#' @param pen_values What values should be penalized; "below" if limit is the lower and "above" if
#' limit is the upper limit.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
nlmixrConstraint <- function(state,
                         limit,
                         time,
                         pen_values,
                         n_constraint,
                         gamma){
  pen_name <- paste0("PEN",n_constraint)
  
  code <- list(
    paste0(pen_name,
           " = 0"),
    paste0("if(",
           ifelse(pen_values=="above",
                  paste0(state," - ",limit),
                  paste0(limit," - ",state)),
           " > 0) ",
           pen_name,
           " = ",
           gamma,
           " * (",
           ifelse(pen_values=="above",
                  paste0(state," - ",limit),
                  paste0(limit," - ",state)),
           ")^2")
  )
  
  out <- list("Name" = pen_name,
              "Time" = time,
              "Code" = code,
              "Type" = "Constraint")
  
  return(out)
  
}

#' Generate nlmixr2 code for constraint AUC
#' 
#' Generate nlmixr2 code for constraint, called by \link[=func]{addConstraintAUC} if
#' software is "nlmixr2".
#' 
#' NULL
#' 
#' @param state State in the structural model to be penalized.
#' @param limit The limit, under/above which the AUC gets calculated for penalization.
#' @param pen_time Time-range of penalization (numeric vector of length 2, i.e., c(min_time, max_time)).
#' @param pen_values What values should be penalized; "below" if limit is the lower and "above" if
#' limit is the upper limit.
#' @param n_constraint The index of the penalization (number of penalizations in the optiProject-object 
#' plus 1)
#' @param gamma Penalization strength.
#' @author Dominic Bräm
nlmixrConstraintAUC <- function(state,
                            limit,
                            pen_time,
                            eval_time,
                            pen_values,
                            n_constraint,
                            gamma){
  pen_name <- paste0("PEN",n_constraint)
  
  ini_state <- paste0("AUC",n_constraint,"(0) = 0")
  
  auc_code <- list(paste0("dAUC",
                          n_constraint,
                          " = 0"),
                   paste0("d",
                          pen_name,
                          " = 0"),
                   paste0("if(",
                          ifelse(pen_values=="above",
                                 paste0(state," - ",limit),
                                 paste0(limit," - ",state)),
                          " > 0) d",
                          pen_name,
                          " = (",
                          ifelse(pen_values=="above",
                                 paste0(state," - ",limit),
                                 paste0(limit," - ",state)),
                          ")^2"),
                   paste0("if(t >= ",
                          pen_time[1],
                          " & t <= ",
                          pen_time[2],
                          "){"
                   ),
                   paste0("   dAUC",
                          n_constraint,
                          " = d",
                          pen_name
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
              "Type" = "ConstraintAUC")
  
  return(out)
  
}
