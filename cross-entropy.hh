#pragma once

#include <functional>
#include <vector>

namespace CrossEntropy {

  using Point = std::vector<double>;
  using Function = std::function<double (const Point &)>;

  void optimize(const Function &f, Point &x, const Point &max_search_radii,
                size_t max_iteration, size_t samples, size_t elite_samples, double tolerance);

  void optimize(const Function &f, Point &x, double search_radius,
                size_t max_iteration, size_t samples, size_t elite_samples, double tolerance);

  Point optimizeBox(const Function &f, const Point &a, const Point &b,
                    size_t max_iteration, size_t samples, size_t elite_samples, double tolerance);

}
