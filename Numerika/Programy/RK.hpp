#ifndef RK_HPP
#define RK_HPP

#include <iostream>
#include <functional>
#include <cmath>
#include <ranges>
#include <vector>

class RK_solver {
  /*
    Solve dy/dx = f(x,y) on [a,b] divided into N intervals (N+1 grid points) given y(a)
  */
private:
  double a {}, b {}, h {};
  int N {};
  std::function<double(double,double)> f; //
  
  double RK4_step (const double& xn, const double& yn) {
    double K1 {h * f( xn, yn)};
    double K2 {h * f( xn + h/2, yn + K1/2)};
    double K3 {h * f( xn + h/2, yn + K2/2)};
    double K4 {h * f( xn + h, yn + K3)};
    double ynp1 = yn + K1/6 + K2/3 + K3/3 + K4/6; // +O(h^5)
    return ynp1;
  }
  
public:
  std::vector<double> x {}, y {};
  RK_solver (const double a, const double b, const int N, std::function<double(double,double)> f, const double y0)
    : a(a), b(b), h ((b-a)/N), f (f)
  {
    if (a>b || a==b) {
      std::cout << "a, b = " << a << ", " << b << ", abort()\n";
      abort();}
    std::cout << "RK4: local error O(" << std::pow(h,5) << "), cummulative error O(" << std::pow(h,4) << ")\n";
    y.push_back(y0);
    for (int i=0; i<=N; ++i) {
      x.push_back(a + i * h);
      if (i>0) {
	y.push_back( RK4_step( x[i-1], y[i-1]));
      }
    }
  }
};
#endif // RK_HPP
