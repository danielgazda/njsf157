#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <ranges>
#include <algorithm>
#include <functional>

class Matrix {
  // Simple matrix class -> use Eigen::Matrix
public:
  int rows {}, cols {};
  std::vector<double> storage {};

  Matrix () {}
  Matrix (const int r, const int c) : rows {r}, cols {c}, storage(rows*cols, 0.0) {}
  Matrix (const int d) : Matrix(d, d) {}
  // Row-wise storage:
  // double operator() (const int r, const int c) const {return storage.at(c + cols*r);}
  // double& operator() (const int r, const int c) {return storage.at(c + cols*r);}
  // Column-wise storage:
  double operator[] (const int r, const int c) const {return storage.at(r + rows*c);}
  double& operator[] (const int r, const int c) {return storage.at(r + rows*c);}
  double operator() (const int r, const int c) const {return storage.at(r + rows*c);}
  double& operator() (const int r, const int c) {return storage.at(r + rows*c);}
  // private:
};

Matrix discretize(const std::function<double(double,double)> f, const std::vector<double> grid) {
  Matrix df(grid.size());
  // Needs c++23
  for (const auto [ipa,pa] : std::views::enumerate(grid)) {
    for (const auto [ipb,pb] : std::views::enumerate(grid)) {
      df(ipa,ipb) = f(pa,pb);
    }
  }
  // for (std::vector<double>::size_type ipp {0}; ipp<grid.size(); ++ipp) {
  //   for (std::vector<double>::size_type ip {0}; ip<grid.size(); ++ip) {
  //     df(ipp,ip) = f(grid[ipp],grid[ip]);
  //   }
  // }
  return df;
}
#endif // MATRIX_HPP
