#pragma once

#include <functional>
#include <vector>

namespace Powell {

  using Point = std::vector<double>;
  using Function = std::function<double (const Point &)>;

  bool optimize(const Function &f, Point &x, size_t max_iteration, double tolerance,
                size_t max_iteration_1d, double tolerance_1d);

  using Function1D = std::function<double (double)>;

  double optimize(const Function1D &f, double x, double step,
                  size_t max_iteration, double tolerance);

  double optimizeBracketed(const Function1D &f, double a, double b,
                           size_t max_iteration, double tolerance);

}
