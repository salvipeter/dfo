#include "downhill.hh"
#include "mads.hh"
#include "powell.hh"

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

enum class Optimizer { MADS, NELDER_MEAD, POWELL };

using Function = std::function<double (const std::vector<double> &)>;

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

void tryMethod(Optimizer optimizer, Function f, std::vector<double> x) {
  static size_t maxit = 1000;
  static double tolerance = 1.0e-8;

  std::chrono::steady_clock::time_point start, stop;
  std::string name = "";
  evaluation_counter = 0;
  start = std::chrono::steady_clock::now();
  switch (optimizer) {
  case Optimizer::MADS:
    name = "MADS";
    MADS::optimize(f, x, maxit, tolerance, 1.0);
    break;
  case Optimizer::NELDER_MEAD:
    name = "Nelder-Mead";
    NelderMead::optimize(f, x, maxit, tolerance, 1.0);
    break;
  case Optimizer::POWELL:
    name = "Powell";
    Powell::optimize(f, x, maxit, tolerance, 100, 1.0e-4);
    break;
  }
  stop = std::chrono::steady_clock::now();
  std::cout << name << ":" << std::endl;
  std::cout << evaluation_counter << " evaluations, result: " << f(x) << std::endl;
  std::cout << "  at:";
  for (double xi : x)
    std::cout << ' ' << xi;
  std::cout << std::endl;
  std::cout << "  evaluation time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;
}

void tryFunction(std::string name, Function f, std::vector<double> init) {
  std::cout << "[" << name << "]" << std::endl;
  tryMethod(Optimizer::MADS, f, init);
  tryMethod(Optimizer::NELDER_MEAD, f, init);
  tryMethod(Optimizer::POWELL, f, init);
}

int main() {
  tryFunction("Ackley", ackley, { 3.0, 12.5, -4.5, 5.7, 10.1 });
  tryFunction("Michaelewicz", michaelewicz, { 1.3, 3.7 });
  tryFunction("Branin", branin, { 12.0, 14.5 });
  tryFunction("Flower", flower, { 1.3, 2.7 });

  std::random_device rd;
  std::default_random_engine re(rd());
  std::uniform_real_distribution<double> rgen(-15, 15);
  std::vector<double> large;
  for (size_t i = 0; i < 100; ++i)
    large.push_back(rgen(re));

  tryFunction("Ackley - large", ackley, large);
  tryFunction("Michaelewicz - large", michaelewicz, large);
}
