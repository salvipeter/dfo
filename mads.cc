#include "mads.hh"

#include <algorithm>
#include <random>
#include <tuple>

namespace MADS {

  using Permutation = std::vector<size_t>;

  std::tuple<Point, Permutation, Permutation>
  generatePositiveSpanningSet(std::default_random_engine &re, double alpha, size_t n) {
    Point result((n + 1) * n);
    int delta = std::round(1.0 / std::sqrt(alpha));
    static std::uniform_int_distribution<size_t> sign(0, 1);
    static double signs[] = { -1.0, 1.0 };
    std::uniform_int_distribution<int> value(-delta + 1, delta - 1);
    for (size_t i = 0; i < n; ++i) {
      size_t index = i * n;
      for (size_t j = 0; j < i; ++j)
        result[index+j] = value(re);
      result[index+i] = delta * signs[sign(re)];
    }

    // Shuffle
    std::vector<size_t> rows, cols;
    for (size_t i = 0; i < n; ++i)
      rows.push_back(i);
    cols = rows;
    rows.push_back(n);
    std::shuffle(rows.begin(), rows.end(), re);
    std::shuffle(cols.begin(), cols.end(), re);

    // Add an extra row
    size_t index = n * n;
    for (size_t i = 0, k = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j, ++k)
        result[index+j] -= result[k];

    return { result, rows, cols };
  }

  Point addmul(const Point &x, const double *d, const Permutation &cols, double alpha) {
    Point result = x;
    for (size_t i = 0; i < result.size(); ++i)
      result[i] += d[cols[i]] * alpha;
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
      auto [span, rows, cols] = generatePositiveSpanningSet(re, alpha, n);
      for (size_t i = 0; i <= n; ++i) {
        const double *d = &span[rows[i] * n];
        auto x1 = addmul(x, d, cols, alpha * step_size);
        double y1 = f(x1);
        if (y1 < y) {
          x = x1;
          y = y1;
          improved = true;
          x1 = addmul(x, d, cols, 3 * alpha * step_size);
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
