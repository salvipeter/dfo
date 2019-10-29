#include "mads.hh"

#include <algorithm>
#include <random>

namespace MADS {

  std::vector<Point> generatePositiveSpanningSet(std::default_random_engine &re,
                                                 double alpha, size_t n) {
    std::vector<Point> result;
    int delta = std::round(1.0 / std::sqrt(alpha));
    static std::uniform_int_distribution<int> sign(0, 1), value(-delta + 1, delta - 1);
    static double signs[] = { -1.0, 1.0 };
    for (size_t i = 0; i < n; ++i) {
      Point row(n, 0);
      row[i] = delta * signs[sign(re)];
      for (size_t j = 0; j < i; ++j)
        row[j] = value(re);
      result.push_back(row);
    }

    // Shuffle
    std::vector<size_t> indices;
    for (size_t i = 0; i < n; ++i)
      indices.push_back(i);
    std::shuffle(indices.begin(), indices.end(), re);
    for (size_t j = 0; j < n - 1; ++j)
      if (indices[j] > j)
        for (size_t i = 0; i < n; ++i)
          std::swap(result[i][j], result[i][indices[j]]);
    std::shuffle(result.begin(), result.end(), re);

    // Add an extra row
    Point row(n, 0);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        row[j] -= result[i][j];
    result.push_back(row);

    return result;
  }

  Point addmul(const Point &x, const Point &d, double alpha) {
    Point result = x;
    for (size_t i = 0; i < result.size(); ++i)
      result[i] += d[i] * alpha;
    return result;
  }

  bool optimize(const Function &f, Point &x, size_t max_iteration, double tolerance,
                double step_size) {
    std::random_device rd;
    std::default_random_engine re(rd());
    double alpha = 1.0, y = f(x);
    size_t n = x.size();
    for (size_t iter = 0; iter < max_iteration && alpha >= tolerance; ++iter) {
      bool improved = false;
      auto span = generatePositiveSpanningSet(re, alpha, n);
      for (size_t i = 0; i <= n; ++i) {
        const auto &d = span[i];
        auto x1 = addmul(x, d, alpha * step_size);
        double y1 = f(x1);
        if (y1 < y) {
          x = x1;
          y = y1;
          improved = true;
          x1 = addmul(x, d, 3 * alpha * step_size);
          y1 = f(x1);
          if (y1 < y) {
            x = x1;
            y = y1;
          }
          break;
        }
      }
      alpha = improved ? std::min(4 * alpha, 1.0) : alpha / 4;
    }
    return alpha < tolerance;
  }

} // namespace
