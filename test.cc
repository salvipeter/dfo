#include "direct.hh"
#include "mads.hh"
#include "nelder-mead.hh"
#include "powell.hh"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>

#ifndef M_PI
static constexpr auto M_PI = std::acos(-1);
#endif

enum class Optimizer { DIRECT, MADS, NELDER_MEAD, POWELL };

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

double mckinnon(const std::vector<double> &x) {
  ++evaluation_counter;
  double tau = 2, theta = 6, phi = 60;
  return (x[0] <= 0 ? phi : 1) * theta * std::pow(std::abs(x[0]), tau) + x[1] + x[1] * x[1];
}

void tryMethod(Optimizer optimizer, Function f, std::vector<double> x, bool print_values) {
  static size_t maxit = 1000;
  static double tolerance = 1.0e-8;

  std::chrono::steady_clock::time_point start, stop;
  std::string name = "";
  evaluation_counter = 0;
  start = std::chrono::steady_clock::now();
  switch (optimizer) {
  case Optimizer::DIRECT:
    name = "DIRECT";
    {
      std::vector<double> a(x.size(), -3), b(x.size(), 5);
      x = DIRECT::optimize(f, a, b, 30, tolerance);
    }
    break;
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
  if (print_values) {
    std::cout << "  at:";
    for (double xi : x)
      std::cout << ' ' << xi;
    std::cout << std::endl;
  }
  std::cout << "  evaluation time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;
}

void tryFunction(std::string name, Function f, std::vector<double> init, bool print_values) {
  std::cout << "[" << name << "]" << std::endl;
  tryMethod(Optimizer::DIRECT, f, init, print_values);
  tryMethod(Optimizer::MADS, f, init, print_values);
  tryMethod(Optimizer::NELDER_MEAD, f, init, print_values);
  tryMethod(Optimizer::POWELL, f, init, print_values);
}

void numericTests() {
  tryFunction("Ackley", ackley, { 3.0, 12.5, -4.5, 5.7, 10.1 }, true);
  tryFunction("Michaelewicz", michaelewicz, { 1.3, 3.7 }, true);
  tryFunction("Branin", branin, { 12.0, 14.5 }, true);
  tryFunction("Flower", flower, { 1.3, 2.7 }, true);
  tryFunction("McKinnon", mckinnon, { 1.3, 2.7 }, true);

  std::mt19937 re(1);
  std::uniform_real_distribution<double> rgen(-15, 15);
  std::vector<double> large;
  for (size_t i = 0; i < 100; ++i)
    large.push_back(rgen(re));

  tryFunction("Ackley - large", ackley, large, false);
  tryFunction("Michaelewicz - large", michaelewicz, large, false);
}

void imageTest() {
  auto fun = michaelewicz;
  std::vector<double> min = { -3, -8 }, max = { 4, 4 }, x = { 0, 1 };
  // auto fun = ackley;
  // std::vector<double> min = { -10, -10 }, max = { 8, 8 }, x = { -7, 3 };
  // auto fun = branin;
  // std::vector<double> min = { 0, 0 }, max = { 10, 10 }, x = { 3, 7 };
  // auto fun = flower;
  // std::vector<double> min = { -4, -4 }, max = { 3, 3 }, x = { -2, -3 };
  // auto fun = mckinnon;
  // std::vector<double> min = { -0.3, -1.5 }, max = { 1, 1 }, x = { 0.8, -1 };

  double step = 0.01;
  {
    std::ofstream f("/tmp/evaluations.txt");
    auto printer = [&](const std::vector<double> &x) {
                     f << x[0] << ' ' << x[1] << std::endl;
                     return fun(x);
                   };
    // DIRECT::optimize(printer, min, max, 30, 1.0e-8);
    MADS::optimize(printer, x, 1000, 1.0e-8, 1.0);
    // NelderMead::optimize(printer, x, 1000, 1.0e-8, 1.0);
    // Powell::optimize(printer, x, 1000, 1.0e-8, 100, 1.0e-4);
  }

  double vmin = std::numeric_limits<double>::infinity(), vmax = -vmin;
  std::vector<double> data;
  for (double x = min[0]; x <= max[0]; x += step)
    for (double y = min[1]; y <= max[1]; y += step) {
      data.push_back(fun({ x, y }));
      if (data.back() < vmin)
        vmin = data.back();
      if (data.back() > vmax)
        vmax = data.back();
    }

  std::ofstream f("/tmp/image.txt");
  size_t index = 0;
  for (double x = min[0]; x <= max[0]; x += step)
    for (double y = min[1]; y <= max[1]; y += step, ++index)
      f << x << ' ' << y << ' ' << (data[index] - vmin) * 8.0 / (vmax - vmin) << std::endl;
}

int main() {
  // numericTests();
  imageTest();
}
