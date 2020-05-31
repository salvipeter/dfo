#pragma once

#include <functional>
#include <vector>

namespace CrossEntropy {

  using Point = std::vector<double>;
  using Function = std::function<double (const Point &)>;

  Point optimize(const Function &f, const Point &a, const Point &b,
                 size_t iteration, size_t samples, size_t elite_samples, double tolerance);

}
