// Original file name: "instabilityTest.R"
// Created: 2020.12.17
// Last modified: 2021.09.13
// License: MIT
// Written by: Astaf'ev Sergey <seryymail@mail.ru>
// This is a part of instabilityTest R package.

#include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include <algorithm>

//' @useDynLib instabilityTest, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export instabilityQuantileTopEst

// Хвост распределения из статьи Detecting Markov chain instability: a Monte Carlo Approach
// @param w - предыдущее значение величины
// @param z - значение случайной величины, стохастически мажорирующей приращение.
// @param fi - оценка сверху модуля приращения функции за 1 шаг.
// @param kappa - оценка снизу модуля начального состояния цепи
// @param delta - "уровень значимости"
// @param sigma - оценка снизу числа шагов алгоритма поиска (почти всегда можно брать 1)
// @param tau - длина симуляции
inline double distTail(double w, double z, double fi, double kappa, double delta, double sigma, double n_w)
{
  double alpha_1 = sigma*fi-sigma*n_w*delta;
  double alpha_2 = pow(((fi+delta)*sigma),2)*n_w;
  double alpha_3 = sigma*fi-w+kappa;
  double alpha_4 = pow((fi*sigma),2)*n_w;

  double exp1 = -pow((z-alpha_1),2)/(2*alpha_2);
  double exp2 = -pow((z-alpha_3),2)/(2*alpha_4);
  double prob = exp(exp1)+n_w*exp(exp2);
  return prob;
}

inline double findZ(double alpha, double w, double fi, double kappa, double delta, double sigma, double tau)
{
  double n_w = ceil(tau/sigma);
  // Сужаем границы хитрым способом, что-бы ускорить сходимость без потери точности.
  double left = std::max(0.0, std::min((fi-delta*n_w)*sigma, sigma*fi+kappa-w));
  double r1 = (-delta*n_w+fi+2*sqrt(537)*sqrt(log(2))*sqrt(n_w*pow(fi+delta, 2)))*sigma;
  double r2 = sigma*fi+kappa-w+sqrt(2)*sqrt(n_w*pow(fi*sigma,2)*(1074*log(2)+log(n_w)));
  double right = std::min(DBL_MAX, std::max(r1, r2));
  // Здесь мы немного хитрим, чтобы ускорить сходимость в большинстве случаев.
  //double value = 2*w;
  double v1 = (-delta*n_w+fi+2*sqrt(15)*sqrt(log(2))*sqrt(n_w*pow(fi+delta, 2)))*sigma;
  double v2 = sigma*fi+kappa-w+sqrt(2)*sqrt(n_w*pow(fi*sigma,2)*(30*log(2)+log(n_w)));
  double value = std::max(abs(v1), abs(v2));
  double prev = value;

  if(distTail(w, left, fi, kappa, delta, sigma, n_w) < alpha){
    return left;
  }
  // Двоичный поиск. Мы можем его использовать, поскольку
  // хвост распределения это монотонна убывающая функция
  do
  {
    double prob = distTail(w, value, fi, kappa, delta, sigma, n_w);
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


//' @title instabilityQuantileTopEst
//' @param alpha - уровень значимости
//' @param fi - оценка сверху модуля приращения функции за 1 шаг.
//' @param delta - отрицательный снос. В большинстве случаев можно задать равным 0.01 или 0.05
//' @param sigma - оценка снизу числа шагов алгоритма поиска нестабильного лямбда (почти всегда можно брать 1)
//' @param kappa - оценка снизу модуля начального состояния цепи
//' @param a - множитель, определящий время следующей симуляции
//' @param b - добавка, определяющая время следующей симуляции
//' @param k - число шагов алгоритма.
// [[Rcpp::export]]
double instabilityQuantileTopEst(double alpha, double fi, double delta, double sigma, double kappa, double a, double b, double k)
{
  // У нас хвост распределения. По этому всё наоборот.
  alpha = 1 - alpha;
  double res = 0;
  for(unsigned int i = 0; i < k; i++){
    res += findZ(alpha, res, fi, kappa, delta, sigma, res*a+b);
  }
  return res;
}

// [[Rcpp::export]]
double genWk(double fi, double delta, double sigma, double kappa, double a, double b, Rcpp::NumericVector u)
{
  double res = 0;
  for(uint_fast64_t i = 0; i < u.size(); i++){
    res +=  findZ(u[i], res, fi, kappa, delta, sigma, res*a+b);
  }
  return res;
}

