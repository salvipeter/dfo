#include "powell.hh"

#include <cmath>
#include <limits>

namespace Powell {

  using Function1D = std::function<double (double)>;

  const double phi = (3.0 - std::sqrt(5.0)) / 2.0;
  
  std::vector<Point> basis(size_t n) {
    std::vector<Point> result;
    for (size_t i = 0; i < n; ++i) {
      Point row(n, 0);
      row[i] = 1;
      result.push_back(row);
    }
    return result;
  }

  Point deviation(const Point &p, const Point &q) {
    Point result = p;
    for (size_t i = 0; i < q.size(); ++i)
      result[i] -= q[i];
    return result;
  }

  double norm(const Point &p) {
    double sum = 0.0;
    for (double x : p)
      sum += x * x;
    return std::sqrt(sum);
  }

  Point addmul(const Point &x, const Point &d, double alpha) {
    Point result = x;
    for (size_t i = 0; i < d.size(); ++i)
      result[i] += d[i] * alpha;
    return result;
  }

  double sign(double x) {
    return x < 0 ? -1 : 1;
  }

  std::pair<double, double> bracket_minimum(const Function1D &f, double a = 0, double b = 1) {
    static const double tiny = std::numeric_limits<double>::epsilon() * 100.0;
    static const double glimit = 100; // maximum magnification in the parabolic-fit step
    double fa = f(a), fb = f(b);
    if (fb > fa) {
      std::swap(a, b);
      std::swap(fa, fb);
    }
    double c = b + phi * (b - a), fc = f(c);
    while (fb > fc) {
      double r = (b - a) * (fb - fa), q = (b - c) * (fb - fa);
      double u =
        b - ((b - c) * q - (b - a) * r) / (2 * std::max(std::abs(q - r), tiny) * sign(q - r));
      double fu;
      double ulim = b + glimit * (c - b);
      if ((b - u) * (u - c) > 0) {
        fu = f(u);
        if (fu < fc)
          return { b, c };
        if (fu > fb)
          return { c, b };
        u = c + phi * (c - b); fu = f(u);
      } else if ((c - u) * (u - ulim) > 0) {
        fu = f(u);
        if (fu < fc) {
          b = c; fb = fc;
          c = u; fc = fu;
          u = c + phi * (c - b); fu = f(u);
        }
      } else if ((u - ulim) * (ulim - c) >= 0) {
        u = ulim; fu = f(u);
      } else {
        u = c + phi * (c - b); fu = f(u);
      }
      a = b; fa = fb;
      b = c; fb = fc;
      c = u; fc = fu;
    }

    return { a, c };
  }

  // As in Chapter 5, Section 8 (pp. 79-80) of
  //  R.P. Brent: Algorithms for minimization without derivatives, Prentice Hall, 1973.
  double brent(const Function1D &f, double a, double b,
               size_t max_iteration, double eps, double t = 1.0e-10) {
    if (a > b)
      std::swap(a, b);
    double x = a + phi * (b - a), v = x, w = x, d, e = 0;
    double fx = f(x), fv = fx, fw = fx;
    for (size_t iter = 0; iter < max_iteration; ++iter) {
      // Main loop
      double m = 0.5 * (a + b);
      double tol = eps * std::abs(x) + t, t2 = 2 * tol;
      // Check stopping criterion
      if (std::abs(x - m) <= t2 - 0.5 * (b - a))
        break;
      double p = 0, q = 0, r = 0;
      if (std::abs(e) > tol) {
        // Fit parabola
        r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = 2 * (q - r);
        if (q > 0)
          p = -p;
        else
          q = -q;
        r = e;
        e = d;
      }
      if (std::abs(p) < std::abs(0.5 * q * r) && p > q * (a - x) && p < q * (b - x)) {
        // A "parabolic interpolation" step
        d = p / q;
        double u = x + d;
        // f must not be evaluated too close to a or b
        if (u - a < t2 || b - u < t2)
          d = x < m ? tol : -tol;
      } else {
        // A "golden section" step
        e = (x < m ? b : a) - x;
        d = phi * e;
      }
      // f must not be evaluated too close to x
      double u = x + (std::abs(d) >= tol ? d : (d > 0 ? tol : -tol)), fu = f(u);
      // Update a, b, v, w and x
      if (fu <= fx) {
        if (u < x)
          b = x;
        else
          a = x;
        v = w; fv = fw;
        w = x; fw = fx;
        x = u; fx = fu;
      } else {
        if (u < x)
          a = u;
        else
          b = u;
        if (fu <= fw || w == x) {
          v = w; fv = fw;
          w = u; fw = fu;
        } else if (fu <= fv || v == x || v == w) {
          v = u; fv = fu;
        }
      }
    }

    return x;
  }

  Point line_search(const Function &f, const Point &x, const Point &d,
                    size_t max_iteration, double tolerance) {
    Function1D objective = [&](double alpha) { return f(addmul(x, d, alpha)); };
    auto [a, b] = bracket_minimum(objective);
    double alpha = brent(objective, a, b, max_iteration, tolerance);
    return addmul(x, d, alpha);
  }

  bool optimize(const Function &f, Point &x, size_t max_iteration, double tolerance,
                size_t max_iteration_1d, double tolerance_1d) {
    size_t n = x.size();
    auto U = basis(n);
    double Delta = std::numeric_limits<double>::infinity();

    for (size_t iter = 0; iter < max_iteration && Delta >= tolerance; ++iter) {
      auto x1 = x;
      for (size_t i = 0; i < n; ++i) {
        auto d = U[i];
        x1 = line_search(f, x1, d, max_iteration_1d, tolerance_1d);
      }
      auto d = deviation(x1, x);
      for (size_t i = 1; i < n; ++i)
        U[i-1] = U[i];
      U[n-1] = d;
      x1 = line_search(f, x1, d, max_iteration_1d, tolerance_1d);
      Delta = norm(deviation(x1, x));
      x = x1;
    }

    return Delta < tolerance;
  }

} // namespace
