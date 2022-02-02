# Original file name: "instabilityTest.R"
# Created: 2020.12.16
# Last modified: 2022.02.02
# License: MIT
# Written by: Astaf'ev Sergey <seryymail@mail.ru>
# This is a part of instabilityTest R package.

#' @export genWkSample
#' @export getInstabilityQuantile

#' @importFrom stats runif
#' @importFrom stats quantile

#' @title genWkSample
#' @inheritParams instabilityQuantileTopEst
#' @param N - размер выборки
genWkSample = function(fi, delta, sigma, kappa, a, b, k, N)
{
  sample = rep(0, N)
  for(i in 1:N){
    sample[i] = genWk(fi, delta, sigma, kappa, a, b, runif(k))
  }
  return(sample)
}

#' @title getInstabilityQuantile
#' @param alpha - уровень значимости
#' @inheritParams instabilityQuantileTopEst
getInstabilityQuantile = function(alpha, fi, delta, sigma, kappa, a, b, k)
{
  sample = genWkSample(fi, delta, sigma, kappa, a, b, k, 1/min(alpha, 1-alpha)*1000)
  return(quantile(sample, names = FALSE, probs = alpha))
}
