#ifndef NUMEROV_H
#define NUMEROV_H

#include <iostream>
#include <functional>
#include <vector>
#include <ranges>
#include <algorithm>
#include <format>
#include <cmath>

class Numerov_solver {
  /*
    Solve y''(x) + k^2(x) = F(x) on [a,b] with y(a) = alpha, y(b) = beta
  */
  // private:
public:
  double a, b; // interval [a,b]
  int N; // number of intervals = numper of points - 1
  double h; // step size (b-a)/N
  std::vector<double> grid, dkk, dF; // coordinate grid and discretized k^2(), F()
  
  std::vector<double> make_grid(const double a, const double b, const int N) {
    if (a>=b) {std::cout<<"Numerov_solver: a>=b, abort()\n"; abort();}
    std::vector<double> grid{};
    for (int i=0; i<N+1; ++i) {
      grid.push_back(a+i*(b-a)/N);
    }
    return grid;
  }
  
  std::vector<double> discretize(std::function<double(double)> fun) {
    std::vector<double> dfun;
    for (const auto& x : grid) {
      dfun.push_back(fun(x));
    }
    return dfun;
  }
  
  inline double step(const double yi, const double yim1,
		     const double kkim1, const double kki, const double kkip1,
		     const double Fim1, const double Fi, const double Fip1) {
    double yip1{
      1.0/(1.0 + h*h/12.0 * kkip1)*(2 * yi * (1 - 5 * h*h/12 * kki)
				    - yim1 * (1 + h*h/12 * kkim1)
				    + h*h/12 * (Fip1 + 10*Fi + Fim1))
    };
    return yip1;
  }
  
  // public:
  Numerov_solver() : a {}, b {}, N {}, h {}, grid {}, dkk {}, dF {} {}
  
  Numerov_solver(const double a, const double b, const int N, const std::function<double(double)> kk, const std::function<double(double)> F)
    : a{a}, b{b},
      N {N},
      h {(b-a)/N},
      grid {make_grid(a,b,N)},
      dkk {discretize(kk)},
      dF {discretize(F)}
  {
  }

  std::vector<double> integrate_outwards(const double alpha, const double delta) {
    std::vector<double> y(grid.size(), 0.0);
    y.at(0) = alpha; // y(a) = alpha
    y.at(1) = delta; // arbitrary (?) delta \approx f'(a) <= find it by shooting at y(b)=beta
    for (int i {1}; i<N; ++i) {
      y.at(i+1) = step(y[i], y[i-1], dkk[i-1], dkk[i], dkk[i+1], dF[i-1], dF[i], dF[i+1]);
    }
    return y;
  }

  std::vector<double> integrate_inwards(const double beta, const double delta) {
    std::vector<double> y(grid.size(), 0.0);
    y.at(N) = beta; // y(b)=beta
    y.at(N-1) = delta; // arbitrary (?) delta \approx f'(a) <= find it by shooting at y(a)=alpha
    for (int i {N-1}; i>0; --i) {
      y.at(i-1) = step(y.at(i), y.at(i+1), dkk.at(i+1), dkk.at(i), dkk.at(i-1), dF.at(i+1), dF.at(i), dF.at(i-1));
    }
    return y;
  }
  
  void print_grid() {
    std::cout << "i " << "x(i)" <<"\n";
    for (const auto& [i,e] : std::views::enumerate(grid))
      std::cout << i << " " << e << "\n";
    std::cout << "h " << h << "\n";
  }

  void print_sol(std::vector<double> y) {
    std::cout<< "i " << "y(i)" <<"\n";
    for (const auto& [i,e] : std::views::enumerate(y))
      std::cout << std::format("{} {}\n", grid[i], e);
  }

  int grid_index_closest_to(const double& val) {
    auto absSubtValCompare = [&val] (const auto& a, const auto& b) {return std::abs(a-val) < std::abs(b-val);};
    auto iterator = std::ranges::min_element(grid, absSubtValCompare);
    auto position = std::ranges::distance(grid.begin(), iterator);
    return position;
  }
};
#endif
