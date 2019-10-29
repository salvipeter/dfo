// based on the algorithm in "Numerical Recipes in C"

#include "powell.hh"

#include <cmath>
#include <limits>
#include <algorithm>

namespace Powell {

  using Function1D = std::function<double (double)>;

  const double GOLD = (1.0 + sqrt(5.0)) / 2.0;
  const double CGOLD = (3.0 - sqrt(5.0)) / 2.0;
  const double VERY_TINY = std::numeric_limits<double>::epsilon() * 100.0;
  const double GLIMIT = 100.0;  // maximum magnification in the parabolic-fit step
  const double ZEPS = 1.0e-10;  // protects against trying to achieve fractional accuracy for 0

  // Global variables (yuck)
  size_t ITMAX;                 // maximum iterations for the Brent method
  double TOL;                   // tolerance for the Brent method
  Point pcom, xicom;
  Function nrfunc;
  // Warning: This class uses global variables for data transfer.
  // Do not run two instances of this class simultaneously!

  double sign(double x) {
    return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
  }

  double sqr(double x) {
    return x * x;
  }

  void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc,
              const Function1D &func) {
    // ax and bx are input parameters
    fa = func(ax);
    fb = func(bx);
    if(fb > fa) {
      std::swap(ax, bx);
      std::swap(fa, fb);
    }
    cx = bx + GOLD * (bx - ax);
    fc = func(cx);
    while(fb > fc) {
      double const r = (bx - ax) * (fb - fc);
      double const q = (bx - cx) * (fb - fa);
      double u, fu;
      u = bx - ((bx - cx) * q - (bx - ax) * r) /
        (2.0 * std::max(std::abs(q - r), VERY_TINY) * sign(q - r));
      double const ulim = bx + GLIMIT * (cx - bx);
      if((bx - u) * (u - cx) > 0.0) {
        fu = func(u);
        if(fu < fc) {
          ax = bx;
          bx = u;
          fa = fb;
          fb = fu;
          return;
        }
        else if(fu > fb) {
          cx = u;
          fc = fu;
          return;
        }
        u = cx + GOLD * (cx - bx);
        fu = func(u);
      } else if((cx - u) * (u - ulim) > 0.0) {
        fu = func(u);
        if(fu < fc) {
          bx = cx; cx = u; u = cx + GOLD * (cx - bx);
          fb = fc; fc = fu; fu = func(u);
        }
      } else if((u - ulim) * (ulim - cx) >= 0.0) {
        u = ulim;
        fu = func(u);
      } else {
        u = cx + GOLD * (cx - bx);
        fu = func(u);
      }
      ax = bx; bx = cx; cx = u;
      fa = fb; fb = fc; fc = fu;
    }
  }

  double brent(double ax, double bx, double cx, const Function1D &f, double tol, double &xmin) {
    double e = 0.0, a, b, d, w, v, u, fw, fv, fu;
    a = std::min(ax, cx);
    b = std::max(ax, cx);
    double x = w = v = bx;
    double fx = fw = fv = f(bx);
    for(unsigned iter = 1; iter <= ITMAX; ++iter) {
      double const xm = 0.5 * (a + b);
      double const tol1 = tol * std::abs(x) + ZEPS;
      double const tol2 = 2.0 * tol1;
      if(std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
        xmin = x;
        return fx;
      }
      if(std::abs(e) > tol1) {
        double r, p, q;
        r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = 2.0 * (q - r);
        if(q > 0.0)
          p = -p;
        q = std::abs(q);
        double const etemp = e;
        e = d;
        if(std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
          e = (x >= xm ? a - x : b - x);
          d = CGOLD * e;
        } else {
          d = p / q;
          u = x + d;
          if(u - a < tol2 || b - u < tol2)
            d = tol1 * sign(xm - x);
        }
      } else {
        e = (x >= xm ? a - x : b - x);
        d = CGOLD * e;
      }
      u = (std::abs(d) >= tol1 ? x + d : x + tol1 * sign(d));
      fu = f(u);
      if(fu <= fx) {
        if(u >= x)
          a = x;
        else
          b = x;
        v = w; w = x; x = u;
        fv = fw; fw = fx; fx = fu;
      } else {
        if(u < x)
          a = u;
        else
          b = u;
        if(fu <= fw || w == x) {
          v = w; w = u;
          fv = fw; fw = fu;
        } else if(fu <= fv || v == x || v == w) {
          v = u;
          fv = fu;
        }
      }
    }
    xmin = x;
    return fx;
  }

  double f1dim(double x) {
    size_t const n = pcom.size();
    std::vector<double> xt(n);
    for(size_t j = 0; j < n; ++j)
      xt[j] = pcom[j] + x * xicom[j];
    return nrfunc(xt);
  }

  void linmin(Point &p, Point &xi, double &fret, const Function &func) {
    size_t const n = p.size();
    pcom = p;
    xicom = xi;
    nrfunc = func;

    double ax = 0.0, xx = 1.0, bx, fa, fx, fb, xmin;
    mnbrak(ax, xx, bx, fa, fx, fb, f1dim);
    fret = brent(ax, xx, bx, f1dim, TOL, xmin);

    for(size_t j = 0; j < n; ++j) {
      xi[j] *= xmin;
      p[j] += xi[j];
    }
  }

  bool powell(Point &p, Point &xi, double ftol, unsigned int &iter, double &fret,
              const Function &func, size_t max_iteration) {
    size_t const n = p.size();
    std::vector<double> pt(p), ptt(n), xit(n);
    fret = func(p);
    for(iter = 1; ; ++iter) {
      double fp = fret, del = 0.0;
      size_t ibig = 0;
      for(size_t i = 0; i < n; ++i) {
        for(size_t j = 0, k = i * n; j < n; ++j, ++k)
          xit[j] = xi[k];
        double const fptt = fret;
        linmin(p, xit, fret, func);
        if(fptt - fret > del) {
          del = fptt - fret;
          ibig = i;
        }
      }
      if(2.0 * (fp - fret) <= ftol * (fabs(fp) + fabs(fret)) + VERY_TINY)
        return true;
      if(iter == max_iteration)
        return false;
      for(size_t j = 0; j < n; ++j) {
        ptt[j] = 2.0 * p[j] - pt[j];
        xit[j] = p[j] - pt[j];
        pt[j] = p[j];
      }
      double const fptt = func(ptt);
      if(fptt < fp) {
        double const t = 2.0 * (fp - 2.0 * fret + fptt) *
          sqr(fp - fret - del) - del * sqr(fp -fptt);
        if(t < 0.0) {
          linmin(p, xit, fret, func);
          for(size_t j = 0, k = ibig * n, l = (n - 1) * n; j < n;
              ++j, ++k, ++l) {
            xi[k] = xi[l];
            xi[l] = xit[j];
          }
        }
      }
    }
    return false;               // should not come here
  }

  bool optimize(const Function &f, Point &x, size_t max_iteration, double tolerance,
                size_t max_1d_iteration, double tolerance_1d) {
    ITMAX = max_1d_iteration;
    TOL = tolerance_1d;
    size_t const n = x.size();
    Point xi(n * n, 0.0);
    for(size_t j = 0; j < n; ++j)
      xi[j * n + j] = 1.0;
    unsigned int iter;
    double fret;
    return powell(x, xi, tolerance, iter, fret, f, max_iteration);
  }

} // namespace
