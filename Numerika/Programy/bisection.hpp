#ifndef BISECTION_H
#define BISECTION_H
template <int maxit=100, double ftol=1.0e-16, double xtol=1.0e-16, bool verbose = true>
double find_root(double a, double b, std::function<double(const double&)> f) {
  int it {1};
  double c {}, cprev {};
  do {
    cprev = (it==1) ? a : c; // Set it to the middle of the interval as well in the first iteration, it does not matter
    c = (a+b)/2;
    if (f(a)*f(c) > 0) {a = c;} else {b = c;}
    if (verbose) {
      std::cout << "it=" << it << " xprev=" << cprev << " x=" << c << " f(xprev)=" << f(cprev) << " f(x)=" << f(c) << "\n";
    }
    ++it;
  } while (((std::abs(f(c)) > ftol) || (std::abs(cprev-c) > xtol)) && (it <= maxit));
  std::cout << "find_root(): " << it-1 << " iteration" << ((it-1>1)?"s":"") << " \n";
  return c;
}
#endif
