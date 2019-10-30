#pragma once

#include <functional>
#include <vector>

namespace DIRECT {

  using Point = std::vector<double>;
  using Function = std::function<double (const Point &)>;

  Point optimize(const Function &f, const Point &a, const Point &b,
                 size_t iteration, double tolerance);

}
