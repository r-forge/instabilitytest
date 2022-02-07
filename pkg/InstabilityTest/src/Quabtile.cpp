// Original file name: "instabilityTest.R"
// Created: 2020.12.17
// Last modified: 2022.02.07
// License: MIT
// Written by: Astaf'ev Sergey <seryymail@mail.ru>
// This is a part of instabilityTest R package.

#include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <iostream>

//' @useDynLib instabilityTest, .registration = TRUE
//' @importFrom Rcpp sourceCpp

// Хвост распределения из статьи Detecting Markov chain instability: a Monte Carlo Approach
// @param w - предыдущее значение величины
// @param z - значение случайной величины, стохастически мажорирующей приращение.
// @param fi - оценка сверху модуля приращения функции за 1 шаг.
// @param kappa - оценка снизу модуля начального состояния цепи
// @param delta - "уровень значимости"
// @param sigma - оценка снизу числа шагов алгоритма поиска (почти всегда можно брать 1)
// @param tau - длина симуляции


inline double e_1(double w, double z, double fi, double kappa, double delta, double sigma, double n_w)
{
  double alpha_1 = sigma*fi-sigma*n_w*delta;
  double alpha_2 = pow(((fi+delta)*sigma),2)*n_w;
  return exp(-pow((z-alpha_1),2)/(2*alpha_2));
}

inline double e_2(double w, double z, double fi, double kappa, double delta, double sigma, double n_w)
{
  double alpha_3 = sigma*fi-w+kappa;
  double alpha_4 = pow((fi*sigma),2)*n_w;
  return n_w*exp(-pow((z-alpha_3),2)/(2*alpha_4));
}


inline double findZ(double alpha, double w, double fi, double kappa, double delta, double sigma, double tau)
{
  double n_w = ceil(tau/sigma);
  double z1_1 = sigma*(fi-n_w*delta+sqrt(2)*sqrt(-n_w*log(alpha)*pow(fi+delta, 2)));
  double z1_2 = sigma*fi+kappa-w+sqrt(2)*sigma*fi*sqrt(log(n_w/alpha)*n_w);
  double z = std::max(z1_1, z1_2);
  
  double prev;
  do{
    prev = z;
    double exp1 = e_1(w, z, fi, kappa, delta, sigma, n_w);
    double exp2 = e_2(w, z, fi, kappa, delta, sigma, n_w);
    double num = (exp1 + exp2 - alpha)*n_w*pow(fi*sigma*(fi+delta),2);
    double den = fi*fi*(sigma*n_w*delta-sigma*fi+z)*exp1+
      pow(fi+delta,2)*(w - sigma*fi - kappa + z)*exp2;
    z = z + num/den;
  }while(abs(z-prev) > 0.0000000001);
  if((z < 0) || std::isnan(z)){
    z = 0;
  }
  return z;
}

// [[Rcpp::export]]
double genWk(double w, double fi, double delta, double sigma, double kappa, double a, double b, Rcpp::NumericVector u)
{
  double res = w;
  for(uint_fast64_t i = 0; i < u.size(); i++){
    res +=  findZ(u[i], res, fi, kappa, delta, sigma, res*a+b);
  }
  return res;
}

