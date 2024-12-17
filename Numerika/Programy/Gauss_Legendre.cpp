#include <numbers>
#include <cmath>
#include <vector>
#include "Gauss_Legendre.hpp"

Gauss_Legendre_mesh::Gauss_Legendre_mesh (const int N, const double x1, const double x2)
  : x(N,0), w(N,0)
{
  double z1, pp, p3, p2, p1;

  const int m{(N + 1) / 2};
  const double
    xm{0.5 * (x2 + x1)},
    xl{0.5 * (x2 - x1)},
    pi{std::numbers::pi_v<double>}, // Requires -std=c++20
    tol{3.0e-14};
  
  for (int i {1}; i<=m; ++i) {
    double z = std::cos(pi * (i - 0.25) / (N + 0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (int j{1}; j<=N; ++j) {
	p3 = p2;
	p2 = p1;
	p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      pp = N * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    } while (std::abs(z1 - z) > tol);
    // Scaling from [-1,1] to [x1,x2]
    x.at(i-1) = xm - xl * z;
    x.at(N+1-(i+1)) = xm + xl * z;
    w.at(i-1) = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w.at(N+1-(i+1)) = w.at(i-1);
    // TODO
    // Scaling to infinite interval
    // Adapt x,w -> t,u for finite:
    // t = 0.5 * (x + 1) * (b - a) + a
    // u = w * 0.5 * (b - a)
    // [-1,1] -> [0,1] -> (-oo,oo)
    // scale = 100 (?)
    // pi_over_4 = np.pi / 4
    // t = scale * np.tan(pi_over_4 * (x + 1))
    // u = scale * pi_over_4 / np.cos(pi_over_4 * (x + 1))**2 * w
  }
};
