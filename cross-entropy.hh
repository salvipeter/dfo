#pragma once

#include <functional>
#include <vector>

namespace CrossEntropy {

  using Point = std::vector<double>;
  using Function = std::function<double (const Point &)>;

  Point optimize(const Function &f, Point mean, Point stddev,
                 size_t max_iteration, size_t samples, size_t elite_samples,
                 double tolerance, double relaxation = 0.9);

  Point optimize(const Function &f, Point &mean, double stddev,
                 size_t max_iteration, size_t samples, size_t elite_samples,
                 double tolerance, double relaxation = 0.9);

  // This is just a convenience function; search is not restricted to the box
  Point optimizeBox(const Function &f, const Point &a, const Point &b,
                    size_t max_iteration, size_t samples, size_t elite_samples,
                    double tolerance, double relaxation = 0.9);

}
