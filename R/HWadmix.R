#' Carries out likelihood ratio test and chi-squared testfor Hardy Weinberg equilibrium for 
#' X chromosomal markers, accounting for the possibility of sex-biased admixture
#'
#' @param X a named vector containing the genotype counts (names should be AA, AB, BB, A, B)
#'
#' @return list with elements
#'   lrt_pval: the p-value of the LRT test
#'   prev_p_f, prev_p_m: the ML estimates of the previous generation allele frequencies in females and males
#'   F_ic: the ML estimate of the inbreeding coefficient
#'   current_p_f, current_p_m: the ML estimates of the current generation allele frequencies in females and males
#'   current_p_AA, current_p_AB, current_p_BB: the ML estimates of the current generation genotype frequencies in females
#'   cs_ml_pval: the p-value of the chi-squared (maximum likelihood) test
#'   you_lrt0_unconstrained_pval: the p-value of the You et al. LRT0 test
#'        where the inbreeding coefficient is not constrained to be in [0,1]
#'
#' @examples
#' snp <- c(A = 122, B = 278, AA = 88, AB = 242, BB = 70)
#' HWadmix(snp)
#'
#' @export
HWadmix <- 
function(X) {
  
  if (!any(c("integer", "numeric") %in% class(X))) {
    stop("X must be a numeric or integer vector")
  }
  if (!all(c("AA", "AB", "BB", "A", "B") %in% names(X))) {
    stop("X must have entries AA, AB, BB, A and B")
  }
  
  grid_base <- expand.grid(p_m = seq(0, 1, by = 0.01), 
                      p_f = seq(0, 1, by = 0.01), 
                      F_ic = seq(-1, 1, by = 0.01),
                      LL = NA, 
                      index = NA)
  grid <- grid_base
  grid[, "index"] <- 1:nrow(grid)
  
  grid[, "LL"] <- LogLikelihood(grid, X)
  simple_grid <- grid[grid[, "F_ic"] == 0, ]
  simple_which_max <- which.max(simple_grid[, "LL"])
   
  simple_model_logLL <- simple_grid[simple_which_max, "LL"]

  complex_which_max <- which.max(grid[, "LL"])
  
  complex_model_logLL <- grid[complex_which_max, "LL"]
  
  estimates <- as.list(grid[complex_which_max, c("p_m", "p_f", "F_ic")])
  statistic <- -2 * (simple_model_logLL - complex_model_logLL)
  pval <- pchisq(statistic, 1, lower.tail = F)
  prev_p_f <- estimates[['p_f']]
  prev_p_m <- estimates[['p_m']]
  
  chi_squared_mm_pval <- ChiSquaredTest(dat = X)
  chi_squared_ml_pval <- ChiSquaredTest(dat = X, 
                                        p_Afprev = simple_grid[simple_which_max, "p_f"], 
                                        p_Amprev = simple_grid[simple_which_max, "p_m"])
  
  F_ic <- estimates[['F_ic']]
  current_p_f <- (prev_p_f + prev_p_m) / 2
  current_p_m <- prev_p_f
  Q <- F_ic * current_p_f * (1 - current_p_f)
  current_p_AA <- current_p_f ^ 2 + Q
  current_p_BB <- (1 - current_p_f) ^ 2 + Q
  current_p_AB <- 1 - current_p_AA - current_p_BB
  
  grid_you <- grid_base
  grid_you[, "LL"] <- LogLikelihoodLRT0(grid_you, X)
  
  complex_you_which_max <- which.max(grid_you[, "LL"])
  complex_you_logLL <- grid_you[complex_you_which_max, "LL"]
  
  simple_you_logLL <- LogLikelihoodLRT0Null(X)
  you_statistic <- -2 * (simple_you_logLL - complex_you_logLL)
  pval_LRT0 <- pchisq(you_statistic, 2, lower.tail = F)
  
  return(list(lrt_pval = pval,
              prev_p_f = prev_p_f,
              prev_p_m = prev_p_m,
              F_ic = F_ic, 
              current_p_f = current_p_f,
              current_p_m = current_p_m,
              current_p_AA = current_p_AA, 
              current_p_AB = current_p_AB, 
              current_p_BB = current_p_BB, 
              cs_ml_pval = chi_squared_ml_pval,
              you_lrt0_unconstrained_pval = pval_LRT0))
}

CalcLL <- function(n_AA, n_AB, n_BB, n_A, n_B, 
                   p_AA, p_AB, p_BB, p_A, p_B) {
  
  df <- list(n = c(n_AA = n_AA, n_AB = n_AB, n_BB = n_BB, n_A = n_A, n_B = n_B), 
             p = list(p_AA = p_AA, p_AB = p_AB, p_BB = p_BB, p_A = p_A, p_B = p_B))
  return(CalcLLDf(df))
}

SafeLog <- function(x) {
  negs <- x < 0
  non_negs <- x >= 0
  logx <- rep(NA, length(x))
  logx[non_negs] <- log(x[non_negs])
  return(logx)
}

CalcLLDf <- function(df) {
  
  terms_to_keep <- names(df$n)[which(df$n > 0)]
  result <- 0
  for (n in terms_to_keep) {
    result <- result + df$n[n] * SafeLog(df$p[[gsub("n_", "p_", n)]])
  }
  return(result)
}

LogLikelihood <- function(grid, dat) {
  # dat is a vector with entries AA, AB, BB, A and B
  # grid is a matrix with columns p_f, p_m and F_ic
  
  # returns numeric vector of log-likelihoods of the data, one for each row of the matrix grid
  
  n_AA <- unname(dat['AA'])
  n_AB <- unname(dat['AB'])
  n_BB <- unname(dat['BB'])
  n_A <- unname(dat['A'])
  n_B <- unname(dat['B'])
  
  p_f <- grid[, "p_f"]
  p_m <- grid[, "p_m"]
  F_ic <- grid[, "F_ic"]
  
  Q <- p_m * (1 - p_f) + p_f * (1 - p_m)
  
  p_AA <- p_m * p_f + F_ic / 2 * Q
  p_AB <- (1 - F_ic) * Q
  p_BB <- (1 - p_m) * (1 - p_f) + F_ic / 2 * Q
  p_A <- p_f
  p_B <- 1 - p_f
  LL <- CalcLL(n_AA = n_AA, n_AB = n_AB, n_BB = n_BB, n_A = n_A, n_B = n_B, 
               p_AA = p_AA, p_AB = p_AB, p_BB = p_BB, p_A = p_A, p_B = p_B)
  LL[p_AA < 0 | p_BB < 0] <- -Inf
  return(LL)
}

LogLikelihoodLRT0 <- function(grid, dat) {
  # dat is a vector with entries AA, AB, BB, A and B
  # grid is a matrix with columns p_f, p_m and F_ic
  
  # returns numeric vector of log-likelihoods of the data, one for each row of the matrix grid
  
  n_AA <- unname(dat['AA'])
  n_AB <- unname(dat['AB'])
  n_BB <- unname(dat['BB'])
  n_A <- unname(dat['A'])
  n_B <- unname(dat['B'])
  
  p_f <- grid[, "p_f"]
  p_m <- grid[, "p_m"]
  F_ic <- grid[, "F_ic"]
  
  Q <- p_f * (1 - p_f)
  
  p_AA <- p_f * p_f + F_ic * Q
  p_AB <- 2 * (1 - F_ic) * Q
  p_BB <- (1 - p_f) ^ 2 + F_ic * Q
  p_A <- p_m
  p_B <- 1 - p_m
  LL <- CalcLL(n_AA = n_AA, n_AB = n_AB, n_BB = n_BB, n_A = n_A, n_B = n_B, 
               p_AA = p_AA, p_AB = p_AB, p_BB = p_BB, p_A = p_A, p_B = p_B)
  LL[p_AA < 0 | p_BB < 0] <- -Inf
  return(LL)
}

CalculateChiSquared <- function(o, e) {
  # calculates chi-squared statistic
  # o: vector of observed counts
  # e: vector of expected counts
  
  exclude <- o == 0 & e == 0
  o_include <- o[!exclude]
  e_include <- e[!exclude]
  
  return(sum((o_include - e_include) ^ 2 / e_include))
}

ChiSquaredTest <- function(dat, p_Afprev = NULL, p_Amprev = NULL) {
  n_AA <- unname(dat['AA'])
  n_AB <- unname(dat['AB'])
  n_BB <- unname(dat['BB'])
  n_A <- unname(dat['A'])
  n_B <- unname(dat['B'])
  
  n_m <- n_A + n_B
  n_f <- n_AA + n_AB + n_BB
  n <- n_m + n_f
  p_Af <- (n_AA + 0.5 * n_AB) / n_f
  p_Am <- n_A / n_m
  if (is.null(p_Amprev)) {
    p_Amprev <- 2 * p_Af - p_Am     # allele frequency among males in previous generation
  }
  if (is.null(p_Afprev)) {
    p_Afprev <- p_Am     # allele frequency among females in previous generation
  }
  
  e_AA <- n_f * p_Afprev * p_Amprev
  e_AB <- n_f * ((1 - p_Afprev) * p_Amprev + p_Afprev * (1 - p_Amprev))
  e_BB <- n_f * (1 - p_Afprev) * (1 - p_Amprev)
  
  observed <- c(n_AA, n_AB, n_BB, n_A, n_B)
  expected = c(e_AA, e_AB, e_BB, n_m * p_Afprev, n_m * (1 - p_Afprev))
  
  statistic <- CalculateChiSquared(o = observed, 
                                   e = expected)
  
  pval <- pchisq(statistic, 1, lower.tail = F)
  if (any(expected == 0 & observed > 0)) {
    pval <- 0
  } 
  if (is.na(pval)) browser()
  return(pval)
}

LogLikelihoodLRT0Null <- function(dat) {
  # dat is a vector with entries AA, AB, BB, A and B
  
  # returns log-likelihood of the data under You et al. LRT0 test null hypothesis
  
  n_AA <- unname(dat['AA'])
  n_AB <- unname(dat['AB'])
  n_BB <- unname(dat['BB'])
  n_A <- unname(dat['A'])
  n_B <- unname(dat['B'])
  
  tot_A <- n_A + 2 * n_AA + n_AB
  tot_B <- n_B + 2 * n_BB + n_AB
  tot <- tot_A + tot_B
  p_A <- tot_A / tot
  p_AA <- p_A ^ 2
  p_AB <- 2 * p_A * (1 - p_A)
  p_BB <- (1 - p_A) ^ 2
  p_B <- 1 - p_A
  LL <- CalcLL(n_AA = n_AA, n_AB = n_AB, n_BB = n_BB, n_A = n_A, n_B = n_B, 
               p_AA = p_AA, p_AB = p_AB, p_BB = p_BB, p_A = p_A, p_B = p_B)
  return(LL)
}