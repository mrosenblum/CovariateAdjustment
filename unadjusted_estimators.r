# dim_estimate: estimate the difference in means
# outcome: ordinal outcome variable (Y)
# treatment: treatment assignment variable (A)
# outcome_levels: a vector of possible values of the outcome, indicating their\
#   order. If `outcome` is not numeric, `outcome_levels` must be specified.
# treatment_levels: a vector of possible values of the treatment assignment. If
#   this is not specified, the unique values are sorted - the first is taken to
#   be reference treatment.
# utilities: a vector of utilities for each level of outcome_levels
#
# outcome_levels allows users to specify the entire range and ordering of values
#   that the outcome variable can take, even if not observed in the sample.
#
# If an outcome does not have numeric levels (e.g. 'Excellent', 'Good', 'Fair',
# 'Poor'), utilities must be provided.

dim_estimate <-
  function(
    outcome,
    treatment,
    outcome_levels = NULL,
    treatment_levels = NULL,
    utilities = NULL
  ){
    
    if (is.null(outcome_levels) & !is.numeric(outcome)) {
      stop("'outcome_levels' should be specified for non-numeric outcomes.")
    }
    
    if(is.null(outcome_levels))
      outcome_levels <- sort(unique(outcome))
    
    if(is.null(treatment_levels)) 
      treatment_levels <- sort(unique(treatment))
    
    if(length(treatment_levels) < 1 | length(treatment_levels) > 2)
      stop("Only binary treatments implemented.")
    
    if(is.null(utilities)){
      if(is.numeric(outcome)) {
        utilities = outcome_levels
      } else{
        stop("Non-numeric outcome supplied without supplying a utility.")
      }
    }
    
    if(length(utilities) != length(outcome_levels))
      stop("Length of utilities (", length(utilities), ") does not match",
           "number of outcome levels (", length(outcome_levels), ")")
    
    pmf <-
      prop.table(
        table(
          factor(outcome, levels = outcome_levels),
          factor(treatment, levels = treatment_levels)
        ),
        margin = 2
      )
    
    matrix(data = utilities, nrow = 1) %*%
      (pmf[, colnames(pmf)[2]] - pmf[, colnames(pmf)[1]])
  }




# mw_estimate: Mann-Whitney estimate
# outcome: ordinal outcome variable (Y)
# treatment: treatment assignment variable (A)
# outcome_levels: a vector of possible values of the outcome, indicating their\
#   order. If `outcome` is not numeric, `outcome_levels` must be specified.
# treatment_levels: a vector of possible values of the treatment assignment. If
#   this is not specified, the unique values are sorted - the first is taken to
#   be reference treatment.
#
# outcome_levels allows users to specify the entire range and ordering of values
#   that the outcome variable can take, even if not observed in the sample. This
#   can be used for non-numeric levels, or to indicate if higher values are
#   less preferable to lower values.
mw_estimate <-
  function(
    outcome,
    treatment,
    outcome_levels = NULL,
    treatment_levels = NULL
  ){
    
    if (is.null(outcome_levels) & !is.numeric(outcome)) {
      stop("'outcome_levels' should be specified for non-numeric outcomes.")
    }
    
    if(is.null(outcome_levels))
      outcome_levels <- sort(unique(outcome))
    
    if(is.null(treatment_levels)) 
      treatment_levels <- sort(unique(treatment))
    
    if(length(treatment_levels) < 1 | length(treatment_levels) > 2)
      stop("Only binary treatments implemented.")
    
    pmf <-
      prop.table(
        table(
          factor(outcome, levels = outcome_levels),
          factor(treatment, levels = treatment_levels)
        ),
        margin = 2
      )
    
    cdf <-
      apply(X = head(pmf, -1),
            MARGIN = 2,
            FUN = cumsum)
    
    sum((c(0, cdf[, colnames(cdf)[1]]) +
           0.5*pmf[, colnames(pmf)[1]])*pmf[, colnames(pmf)[2]])
  }

lor_estimate <-
  function(
    outcome,
    treatment,
    outcome_levels = NULL,
    treatment_levels = NULL
  ){
    if (is.null(outcome_levels) & !is.numeric(outcome)) {
      stop("'outcome_levels' should be specified for non-numeric outcomes.")
    }
    
    if(is.null(outcome_levels))
      outcome_levels <- sort(unique(outcome))
    
    if(is.null(treatment_levels)) 
      treatment_levels <- sort(unique(treatment))
    
    if(length(treatment_levels) < 1 | length(treatment_levels) > 2)
      stop("Only binary treatments implemented.")
    
    pmf <-
      prop.table(
        table(
          factor(outcome, levels = outcome_levels),
          factor(treatment, levels = treatment_levels)
        ),
        margin = 2
      )
    
    cdf <-
      apply(X = head(pmf, -1),
            MARGIN = 2,
            FUN = cumsum)
    
    mean(log((cdf[, colnames(cdf)[2]]*(1 - cdf[, colnames(cdf)[1]]))/
               (cdf[, colnames(cdf)[1]]*(1 - cdf[, colnames(cdf)[2]]))))
  }