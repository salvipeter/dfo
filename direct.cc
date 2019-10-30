#include "direct.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <queue>

namespace DIRECT {

  struct Interval {
    Point c;                    // center
    double y;                   // central value
    std::vector<size_t> depths; // number of divisions in each dimension

    Interval(const Point &c, double y, const std::vector<size_t> depths) 
      : c(c), y(y), depths(depths) { }
    size_t min_depth() const { return *std::min_element(depths.begin(), depths.end()); }
    bool operator<(const Interval &i) const { return y < i.y; }
  };
  using Queue = std::priority_queue<Interval>;
  using Intervals = std::map<size_t, Queue>;

  void addInterval(Intervals &intervals, const Interval &interval) {
    size_t d = interval.min_depth();
    if (intervals.find(d) == intervals.end())
      intervals[d] = Queue();
    intervals[d].push(interval);
  }

  Point deviation(const Point &p, const Point &q) {
    Point result = p;
    for (size_t i = 0; i < q.size(); ++i)
      result[i] -= q[i];
    return result;
  }

  Point addmul(const Point &x, const Point &d, double alpha) {
    Point result = x;
    for (size_t i = 0; i < d.size(); ++i)
      result[i] += d[i] * alpha;
    return result;
  }

  Point addmul(const Point &x, const Point &d, const Point &alphas) {
    Point result = x;
    for (size_t i = 0; i < d.size(); ++i)
      result[i] += d[i] * alphas[i];
    return result;
  }

  std::vector<Interval> getOptIntervals(const Intervals &intervals,
                                        double tolerance, double y_best) {
    size_t max_depth = 0;
    for (const auto &kv : intervals)
      if (kv.first > max_depth)
        max_depth = kv.first;

    std::vector<Interval> stack;
    stack.push_back(intervals.at(max_depth).top());

    size_t d = max_depth;
    while (d > 0) {
      --d;
      if (intervals.find(d) == intervals.end() || intervals.at(d).empty())
        continue;
      const auto &interval = intervals.at(d).top();
      double x = 0.5 * std::pow(3.0, -(int)interval.min_depth()), y = interval.y;
      while (!stack.empty()) {
        const auto &interval1 = stack.back();
        double x1 = 0.5 * std::pow(3.0, -(int)interval1.min_depth()), y1 = interval1.y;
        double l1 = (y - y1) / (x - x1);
        if (y1 - l1 * x1 > y_best - tolerance || y < y1)
          stack.pop_back();
        else if (stack.size() > 1) {
          const auto &interval2 = stack[stack.size()-2];
          double x2 = 0.5 * std::pow(3.0, -(int)interval2.min_depth()), y2 = interval2.y;
          double l2 = (y1 - y2) / (x1 - x2);
          if (l2 > l1)
            stack.pop_back();
          else
            break;
        } else
          break;
      }
      stack.push_back(interval);
    }

    return stack;
  }

  Point basis(size_t i, size_t n) {
    Point row(n, 0);
    row[i] = 1;
    return row;
  }

  std::vector<size_t> sortPermutation(std::vector<double> &x) {
    size_t n = x.size();
    std::vector<size_t> perm;
    for (size_t i = 0; i < n; ++i)
      perm.push_back(i);

    // Insertion sort
    for (size_t i = 1; i < n; ++i)
      for (size_t j = i; j > 0 && x[j-1] > x[j]; --j) {
        std::swap(x[j], x[j-1]);
        std::swap(perm[j], perm[j-1]);
      }
    return perm;
  }

  std::vector<Interval> divide(const Function &f, const Interval &interval) {
    const Point &c = interval.c;
    size_t d = interval.min_depth();
    size_t n = c.size();
    std::vector<size_t> dirs;
    for (size_t i = 0; i < n; ++i)
      if (interval.depths[i] == d)
        dirs.push_back(i);
    std::vector<std::pair<Point, Point>> cs;
    std::vector<std::pair<double, double>> vs;
    std::vector<double> minvals;
    for (size_t i : dirs) {
      cs.push_back({ addmul(c, basis(i, n), std::pow( 3.0, -(int)d-1)),
                     addmul(c, basis(i, n), std::pow(-3.0, -(int)d-1)) });
      vs.push_back({ f(cs.back().first), f(cs.back().second) });
      minvals.push_back(std::min(vs.back().first, vs.back().second));
    }

    std::vector<Interval> intervals;
    std::vector<size_t> depths = interval.depths;
    for (size_t j : sortPermutation(minvals)) {
      depths[dirs[j]] += 1;
      intervals.push_back({ cs[j].first, vs[j].first, depths });
      intervals.push_back({ cs[j].second, vs[j].second, depths });
    }
    intervals.push_back({ c, interval.y, depths });
    return intervals;
  }

  Point optimize(const Function &f, const Point &a, const Point &b,
                 size_t iteration, double tolerance) {
    auto Delta = deviation(b, a);
    auto g = [&](const Point &x) { return f(addmul(a, x, Delta)); };

    Intervals intervals;
    size_t n = a.size();
    Point c(n, 0.5);
    std::vector<size_t> zero(n, 0);
    Interval interval(c, g(c), zero);
    addInterval(intervals, interval);

    Point c_best = interval.c;
    double y_best = interval.y;

    for (size_t k = 0; k < iteration; ++k) {
      auto S = getOptIntervals(intervals, tolerance, y_best);
      std::vector<Interval> to_add;
      for (const auto &interval : S) {
        auto divided = divide(g, interval);
        to_add.insert(to_add.end(), divided.begin(), divided.end());
        intervals[interval.min_depth()].pop();
      }
      for (const auto &interval : to_add) {
        addInterval(intervals, interval);
        if (interval.y < y_best) {
          c_best = interval.c;
          y_best = interval.y;
        }
      }
    }

    return addmul(a, c_best, Delta);
  }

} // namespace
