#include "downhill.hh"
#include "powell.hh"

#include <cmath>
#include <iostream>

size_t evaluation_counter;

// nD functions:

// optimal value 0 at the origin
double ackley(const std::vector<double> &x) {
  ++evaluation_counter;
  static double a = 20, b = 0.2, c = 2 * M_PI;
  size_t d = x.size();
  double sqsum = 0, cosum = 0;
  for (double xi : x) {
    sqsum += xi * xi;
    cosum += std::cos(c * xi);
  }
  return -a * std::exp(-b * std::sqrt(sqsum / d)) - std::exp(cosum / d) + a + std::exp(1);
}

// optimal value in 2D: -1.8011 at ~ [2.20, 1.57]
double michaelewicz(const std::vector<double> &x) {
  ++evaluation_counter;
  static double m = 10;
  size_t d = x.size();
  double sum = 0;
  for (size_t i = 1; i <= d; ++i) {
    double v = x[i-1];
    sum += std::sin(v) * std::pow(std::sin(i * v * v / M_PI), 2 * m);
  }
  return -sum;
}

// 2D functions:

// optimal value ~0.397887
double branin(const std::vector<double> &x) {
  ++evaluation_counter;
  static double a = 1, b = 5.1 / (4 * M_PI * M_PI), c = 5 / M_PI, r = 6, s = 10, t = 1 / (8 * M_PI);
  return a * std::pow(x[1] - b * x[0] * x[0] + c * x[0] - r, 2) + s * (1 - t) * cos(x[0]) + s;
}

// global minimum should be at the origin, but undefined
double flower(const std::vector<double> &x) {
  ++evaluation_counter;
  static double a = 1, b = 1, c = 4;
  double norm = std::sqrt(x[0] * x[0] + x[1] * x[1]);
  return a * norm + b * std::sin(c * std::atan(x[1] / x[0]));
}

int main() {
  size_t maxit = 1000;
  double tolerance = 1.0e-8;
  std::vector<double> x;

  std::cout << "[Ackley]" << std::endl;

  std::cout << "Nelder-Mead:" << std::endl;
  evaluation_counter = 0;
  x = { 0.3, 0.7, -0.4, 0.2 };
  NelderMead::optimize(ackley, x, maxit, tolerance, 1.0);
  std::cout << evaluation_counter << " evaluations, result: " << ackley(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;

  std::cout << "Powell:" << std::endl;
  evaluation_counter = 0;
  x = { 0.3, 0.7, -0.4, 0.2 };
  Powell::optimize(ackley, x, maxit, tolerance, 100, 1.0e-4);
  std::cout << evaluation_counter << " evaluations, result: " << ackley(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;

  std::cout << "[Michaelewicz]" << std::endl;

  std::cout << "Nelder-Mead:" << std::endl;
  evaluation_counter = 0;
  x = { 0.3, 0.7 };
  NelderMead::optimize(michaelewicz, x, maxit, tolerance, 1.0);
  std::cout << evaluation_counter << " evaluations, result: " << michaelewicz(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;

  std::cout << "Powell:" << std::endl;
  evaluation_counter = 0;
  x = { 0.3, 0.7 };
  Powell::optimize(michaelewicz, x, maxit, tolerance, 100, 1.0e-4);
  std::cout << evaluation_counter << " evaluations, result: " << michaelewicz(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;

  std::cout << "[Branin]" << std::endl;

  std::cout << "Nelder-Mead:" << std::endl;
  evaluation_counter = 0;
  x = { 0.3, 0.7 };
  NelderMead::optimize(branin, x, maxit, tolerance, 1.0);
  std::cout << evaluation_counter << " evaluations, result: " << branin(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;

  std::cout << "Powell:" << std::endl;
  evaluation_counter = 0;
  x = { 0.3, 0.7 };
  Powell::optimize(branin, x, maxit, tolerance, 100, 1.0e-4);
  std::cout << evaluation_counter << " evaluations, result: " << branin(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;

  std::cout << "[Flower]" << std::endl;

  std::cout << "Nelder-Mead:" << std::endl;
  evaluation_counter = 0;
  x = { 0.3, 0.7 };
  NelderMead::optimize(flower, x, maxit, tolerance, 1.0);
  std::cout << evaluation_counter << " evaluations, result: " << flower(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;

  std::cout << "Powell:" << std::endl;
  evaluation_counter = 0;
  x = { 0.3, 0.7 };
  Powell::optimize(flower, x, maxit, tolerance, 100, 1.0e-4);
  std::cout << evaluation_counter << " evaluations, result: " << flower(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;
}
