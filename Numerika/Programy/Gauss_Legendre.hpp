#ifndef GAULEG_H
#define GAULEG_H
//#include <numbers>
//#include <cmath>
#include <vector>

class Gauss_Legendre_mesh {
  /* Adapted FORTRAN subroutine I got from Petr Navrátil
     Gauss_Legendre_mesh(N,-1,1) gives the same as numpy.polynomial.legendre.leggauss(N) */
  public:
  std::vector<double> x, w; // points and weights
  Gauss_Legendre_mesh (const int N, const double x1, const double x2);
  };
#endif
