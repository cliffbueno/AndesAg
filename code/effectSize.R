# make effect size pie chart as in Darcy et al. (2018)
eta_sq <- function(aov_in){
  aovtab <- summary(aov_in)[[1]]
  n_terms <- length(aovtab[["Sum Sq"]])
  SSt <- sum(aovtab[["Sum Sq"]])
  output <- aovtab[["Sum Sq"]] / SSt
  names(output) <- rownames(aovtab)
  return(output)
}

eta_sq_adonis <- function(aov_in){
  aovtab <- as.data.frame(aov_in)
  n_terms <- length(aovtab[["SumOfSqs"]])
  SSt <- sum(aovtab[["SumOfSqs"]])
  output <- aovtab[["SumOfSqs"]] / SSt
  names(output) <- rownames(aovtab)
  return(output)
}

eta_sq_glm <- function(glm_in){
  aovtab <- Anova(glm_in, test.statistic = "F")
  n_terms <- length(aovtab[["Sum Sq"]])
  SSt <- sum(aovtab[["Sum Sq"]])
  output <- aovtab[["Sum Sq"]] / SSt
  names(output) <- rownames(aovtab)
  return(output)
}