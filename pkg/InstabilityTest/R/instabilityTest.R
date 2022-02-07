# Original file name: "instabilityTest.R"
# Created: 2020.12.16
# Last modified: 2022.02.07
# License: MIT
# Written by: Astaf'ev Sergey <seryymail@mail.ru>
# This is a part of instabilityTest R package.

#' @export genWkSample
#' @export getInstabilityQuantile

#' @importFrom stats runif
#' @importFrom stats quantile

#' @title genWkSample
#' @param N - размер выборки
#' @param sample - выборка, ранее уже сгенерированная этой функцией.
genWkSample = function(fi, delta, sigma, kappa, a, b, k, N, sample = NULL)
{
  if(is.null(sample)){
    sample = rep(0, N)
  }
  for(i in 1:N){
    sample[i] = genWk(sample[i], fi, delta, sigma, kappa, a, b, runif(k))
  }
  return(sample)
}

#' @title getInstabilityQuantile
#' @param alpha - уровень значимости
getInstabilityQuantile = function(alpha, fi, delta, sigma, kappa, a, b, k)
{
  sample = genWkSample(fi, delta, sigma, kappa, a, b, k, ceiling(1/min(alpha, 1-alpha)*1000))
  return(quantile(sample, names = FALSE, probs = alpha))
}
