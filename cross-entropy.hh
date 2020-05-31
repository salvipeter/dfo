#pragma once

#include <functional>
#include <vector>

namespace CrossEntropy {

  using Point = std::vector<double>;
  using Function = std::function<double (const Point &)>;

  Point optimize(const Function &f, Point &mean, const Point &stddev,
                 size_t max_iteration, size_t samples, size_t elite_samples, double tolerance);

  Point optimize(const Function &f, Point &mean, double stddev,
                 size_t max_iteration, size_t samples, size_t elite_samples, double tolerance);

  // This is just a convenience function; search is not restricted to the box
  Point optimizeBox(const Function &f, const Point &a, const Point &b,
                    size_t max_iteration, size_t samples, size_t elite_samples, double tolerance);

}
