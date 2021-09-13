# Original file name: "instabilityTest.R"
# Created: 2020.12.16
# Last modified: 2020.12.17
# License: MIT
# Written by: Astaf'ev Sergey <seryymail@mail.ru>
# This is a part of instabilityTest R package.

#' instabilityTest
#' @useDynLib instabilityTest, .registration = TRUE
#' @importFrom Rcpp sourceCpp

#' @export instabilityTest
#' @export instabilityQuantile

#' @name instabilityTest
#' @description Test simulation result on instability.
#' @param w increment of Markov chain in simulation time.
#' @param fi upper estimate for the absolute value of Markov chain increment in one step.
#' @param kappa lower estimate for the initial state of Markov chain.
#' @param delta the downward drift that a process must exhibit in order to be stable. Significance level analogue.
#' @param tau simulation time.
#' @return TRUE, if  Markov chain is instable. FALSE doesn't mean that chain is stable.
instabilityTest = function(w, fi, kappa, delta, tau)
{
  if (w > findZ(0, fi, kappa, delta, 1, tau)){
    return (TRUE)
  }else{
    return (FALSE)
  }
}
