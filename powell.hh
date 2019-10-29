#pragma once

#include <functional>
#include <vector>

namespace Powell {

  using Point = std::vector<double>;
  using Function = std::function<double (const Point &)>;

  bool optimize(const Function &f, Point &x, size_t max_iteration, double tolerance,
                size_t max_1d_iteration, double tolerance_1d);

}
