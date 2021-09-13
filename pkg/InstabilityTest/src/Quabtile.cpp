// Original file name: "instabilityTest.R"
// Created: 2020.12.17
// Last modified: 2021.09.13
// License: MIT
// Written by: Astaf'ev Sergey <seryymail@mail.ru>
// This is a part of instabilityTest R package.

#include <Rcpp.h>
#include <cmath>
#include <cfloat>

// Хвост распределения из статьи Detecting Markov chain instability: a Monte Carlo Approach
// @param w - значение функции в начальном состояние цепи, из которой мы начинали симуляцию.
// @param z - значение случайной величины, стохастически мажорирующей приращение.
// @param fi - оценка сверху модуля приращения функции за 1 шаг.
// @param kappa - оценка снизу модуля начального состояния цепи (т.е. X_{k-1})
// @param delta - "уровень значимости"
// @param sigma - оценка снизу числа шагов алгоритма поиска (почти всегда можно брать 1)
// @param tau - длина симуляции
double distTail(double w, double z, double fi, double kappa, double delta, double sigma, double tau)
{
  double n_w = ceil(tau/sigma);
  double alpha_1 = sigma*fi-sigma*n_w*delta;
  double alpha_2 = pow(((fi+delta)*sigma),2)*n_w;
  double alpha_3 = sigma*fi-w+kappa;
  double alpha_4 = pow((fi*sigma),2)*n_w;

  double exp1 = -pow((z-alpha_1),2)/(2*alpha_2);
  double exp2 = -pow((z-alpha_3),2)/(2*alpha_4);
  double prob = exp(exp1)+n_w*exp(exp2);
  return prob > 1.0 ? 1.0 : prob;
}

// [[Rcpp::export]]
double findZ(double alpha, double w, double fi, double kappa, double delta, double sigma, double tau)
{
  double left = 0;
  double right = DBL_MAX;
  double value = (tau/sigma - left)/2;
  double prev = value;

  // Двоичный поиск. Мы можем его использовать, поскольку
  // хвост распределения это монотонна убывающая функция
  do
  {
    double prob = distTail(w, value, fi, kappa, delta, sigma, tau);
    if(prob <= alpha){
      right = value;
    } else{
      left = value;
    }
    prev = value;
    value = left + (right - left)/2;

  }while(fabs(prev - value) > 0.0000000001);
  return value;
}


//' instabilityQuantile
//' @description Calculation of a quantile with given parameters for a random variable stochastically dominating over a Markov chain.
//' @param alpha the significance level.
//' @param w initial state of Markov chain.
//' @param fi upper estimate for the absolute value of Markov chain increment in one step.
//' @param delta the downward drift that a process must exhibit in order to be stable.
//' @param sigma the lower estimate for the number of steps.
//' @param kappa the lower estimate for initial state.
//' @param a a parameter for time calculation. Time is calculated by a formula: a*(X_{i-1})+b.
//' @param b b parameter for time calculation. Time is calculated by a formula: a*(X_{i-1})+b.
//' @param k number of steps to select a new point.
//' @return Quantile with given parameters
// [[Rcpp::export]]
double instabilityQuantile(double alpha, double w, double fi, double delta, double sigma, double kappa, double a, double b, double k)
{
  double res = 0;
  double last = w;
  for(unsigned int i = 0; i < k; i++){
    last = findZ(alpha, w, fi, kappa, delta, sigma, res*a+b);
    res = res + last;
  }
  return res;
}
