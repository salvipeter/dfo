#include "cross-entropy.hh"

#include <random>

namespace CrossEntropy {

  Point centroid(const std::vector<Point> &points) {
    size_t n = points[0].size(), m = points.size();
    Point result(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < m; ++j)
        result[i] += points[j][i];
      result[i] /= m;
    }
    return result;
  }

  Point fitStdDev(const std::vector<Point> &points, const Point &mean) {
    size_t n = points[0].size(), m = points.size();
    Point result(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < m; ++j)
        result[i] += std::pow(points[j][i] - mean[i], 2);
      result[i] = std::sqrt(result[i] / (m - 1));
    }
    return result;
  }

  Point deviation(const Point &p, const Point &q) {
    Point result = p;
    for (size_t i = 0; i < q.size(); ++i)
      result[i] -= q[i];
    return result;
  }

  double distance(const Point &p, const Point &q) {
    double result = 0;
    for (size_t i = 0; i < p.size(); ++i)
      result += std::pow(p[i] - q[i], 2);
    return std::sqrt(result);
  }

  Point optimize(const Function &f, Point &mean, const Point &stddev,
                 size_t max_iteration, size_t samples, size_t elite_samples, double tolerance) {
    // Set up the multivariate Gaussian distribution
    std::random_device rd;
    std::default_random_engine re(rd());
    std::vector<std::normal_distribution<>> P;
    size_t n = mean.size();
    for (size_t i = 0; i < n; ++i)
      P.push_back(std::normal_distribution<>{ mean[i], stddev[i] });

    // Main iteration
    std::vector<Point> s(samples);
    std::vector<double> fs(samples);
    auto last_mean = mean;
    for (size_t iter = 0; iter < max_iteration; ++iter) {
      // Generate samples
      for (size_t j = 0; j < samples; ++j) {
        s[j].clear();
        for (size_t i = 0; i < n; ++i)
          s[j].push_back(P[i](re));
        fs[j] = f(s[j]);
      }

      // Sort values
      for (size_t j = 1; j < samples; ++j)
        for (size_t k = j; k > 0 && fs[k-1] > fs[k]; --k) {
          std::swap(s[k], s[k-1]);
          std::swap(fs[k], fs[k-1]);
        }

      // Fit a new distribution on the best values
      s.resize(elite_samples);  // does not change the capacity
      mean = centroid(s);
      auto stddev = fitStdDev(s, mean);
      for (size_t i = 0; i < n; ++i)
        P[i] = std::normal_distribution<>{ mean[i], stddev[i] };
      s.resize(samples);

      // Check for convergence
      if (distance(mean, last_mean) < tolerance)
        return mean;
      last_mean = mean;
    }

    return mean;
  }

  Point optimize(const Function &f, Point &mean, double stddev,
                 size_t max_iteration, size_t samples, size_t elite_samples, double tolerance) {
    Point stddevs(mean.size(), stddev);
    return optimize(f, mean, stddevs, max_iteration, samples, elite_samples, tolerance);
  }

  Point optimizeBox(const Function &f, const Point &a, const Point &b,
                    size_t max_iteration, size_t samples, size_t elite_samples, double tolerance) {
    auto x = centroid({ a, b });
    auto radii = deviation(a, b);
    size_t n = x.size();
    for (size_t i = 0; i < n; ++i)
      radii[i] = std::abs(radii[i]) / 6;
    return optimize(f, x, radii, max_iteration, samples, elite_samples, tolerance);
  }

}
