#include "powell.hh"

#include <cmath>
#include <limits>

namespace Powell {

  using Function1D = std::function<double (double)>;

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
    return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
  }

  // As in Chapter 10 (pp. 491-492) of
  //   W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flanner:
  //     Numerical Recipes - The Art of Scientific Computing (3rd Ed).
  //       Cambridge University Press, 2007.
  std::pair<double, double> bracket_minimum(const Function1D &f, double a = 0, double b = 1) {
    static const double gold = (1.0 + std::sqrt(5.0)) / 2.0;
    static const double tiny = std::numeric_limits<double>::epsilon() * 100.0;
    static const double glimit = 100; // maximum magnification in the parabolic-fit step
    double fa = f(a), fb = f(b);
    if (fb > fa) {              // switch roles so that we can go donwhill from a to b
      std::swap(a, b);
      std::swap(fa, fb);
    }
    double c = b + gold * (b - a), fc = f(c); // first guess for c
    while (fb > fc) {                         // keep returning here until we bracket
      double r = (b - a) * (fb - fc), q = (b - c) * (fb - fa);
      double u =                // compute u by parabolic extrapolation from a, b, c
        b - ((b - c) * q - (b - a) * r) / (2 * std::max(std::abs(q - r), tiny) * sign(q - r));
      double fu;
      double ulim = b + glimit * (c - b); // we won't go farther than this
      // Test various possibilities:
      if ((b - u) * (u - c) > 0) { // parabolic u is between b and c
        fu = f(u);
        if (fu < fc)            // got a minimum between b and c
          return b < c ? std::make_pair(b, c) : std::make_pair(c, b);
        if (fu > fb)            // got a minimum between a and u
          return a < u ? std::make_pair(a, u) : std::make_pair(u, a);
        u = c + gold * (c - b); fu = f(u); // parabolic fit was no use, use default magnification
      } else if ((c - u) * (u - ulim) > 0) { // parabolic fit between c and the allowed limit
        fu = f(u);
        if (fu < fc) {
          b = c; fb = fc;
          c = u; fc = fu;
          u = c + gold * (c - b); fu = f(u); // (2nd Ed. has this correctly)
        }
      } else if ((u - ulim) * (ulim - c) >= 0) { // limit parabolic u to maximum allowed value
        u = ulim; fu = f(u);
      } else {                  // reject parabolic u, use default magnification
        u = c + gold * (c - b); fu = f(u);
      }
      // Eliminate oldest point and continue
      a = b; fa = fb;
      b = c; fb = fc;
      c = u; fc = fu;
    }

    return a < c ? std::make_pair(a, c) : std::make_pair(c, a);
  }

  // As in Chapter 5, Section 8 (pp. 79-80) of
  //  R.P. Brent: Algorithms for minimization without derivatives, Prentice Hall, 1973.
  std::pair<double, double> brent(const Function1D &f, double a, double b,
                                  size_t max_iteration, double eps, double t = 1.0e-10) {
    static const double c = (3.0 - std::sqrt(5.0)) / 2.0;
    double x = a + c * (b - a), v = x, w = x, d, e = 0;
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
        d = c * e;
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

    return { x, fx };
  }

  std::pair<Point, double> line_search(const Function &f, const Point &x, const Point &d,
                                       size_t max_iteration, double tolerance) {
    Function1D objective = [&](double alpha) { return f(addmul(x, d, alpha)); };
    auto [a, b] = bracket_minimum(objective);
    auto [alpha, fx] = brent(objective, a, b, max_iteration, tolerance);
    return { addmul(x, d, alpha), fx };
  }

  // Direction update based on Powell's heuristic (see Numerical Recipes in C, Section 10.7.3)
  bool optimize(const Function &f, Point &x, size_t max_iteration, double tolerance,
                size_t max_iteration_1d, double tolerance_1d) {
    size_t n = x.size();
    double fx = f(x);
    auto U = basis(n);
    double Delta = std::numeric_limits<double>::infinity();

    for (size_t iter = 0; iter < max_iteration && Delta >= tolerance; ++iter) {
      // Minimize in all directions and save the best direction
      Point x1 = x;
      double fx_prev = fx, fx1 = 0, largest_decrease = 0;
      size_t best_dir = 0;
      for (size_t i = 0; i < n; ++i) {
        auto d = U[i];
        std::tie(x1, fx1) = line_search(f, x1, d, max_iteration_1d, tolerance_1d);
        if (fx_prev - fx1 > largest_decrease) {
          largest_decrease = fx_prev - fx1;
          best_dir = i;
        }
        fx_prev = fx1;
      }

      // Update directions
      auto d = deviation(x1, x), extrapolated = addmul(x1, d, 1);
      double fextrapolated = f(extrapolated);
      // Only update directions when the overall deviation direction seems good
      if (fextrapolated < fx &&
          2 * (fx - 2 * fx1 + fextrapolated) * std::pow(fx - fx1 - largest_decrease, 2) -
          largest_decrease * std::pow(fx - fextrapolated, 2) < 0) {
        // Replace the best direction with the deviation direction
        U[best_dir] = U[n-1];
        U[n-1] = d;
        // Minimize in the new direction
        std::tie(x1, fx1) = line_search(f, x1, d, max_iteration_1d, tolerance_1d);
      }

      Delta = norm(deviation(x1, x));
      x = x1; fx = fx1;
    }

    return Delta < tolerance;
  }

} // namespace
