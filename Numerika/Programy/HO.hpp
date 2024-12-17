#ifndef HO_H
#define HO_H

#include <vector>

class radial_HO_wave_functions {

private:
  int nmax, lmax, maxgam;
  double fgamal(const int arg);
  double fdsq(const int n, const int l);
  
public:
  
  class u_type {
  private:
    int nmax, lmax;
    std::vector<double> storage;
  public:
    u_type ();
    u_type (const int nmax, const int lmax);
    double operator[] (const int n, const int l) const;
    double& operator[] (const int n, const int l);
    double operator() (const int n, const int l) const;
    double& operator() (const int n, const int l);
  } u;

  class ur_type {
  private:
    int nmax, lmax;
    std::vector<double> rgrid;
    std::vector<double> storage;
  public:
    ur_type ();
    ur_type (const int nmax, const int lmax, std::vector<double> rgrid);
    double operator[] (const int n, const int l, const int ir) const;
    double& operator[] (const int n, const int l, const int ir);
    double operator() (const int n, const int l, const int ir) const;
    double& operator() (const int n, const int l, const int ir);
  } ur;
  
  radial_HO_wave_functions (const int n_re, const int l_re, const double anu, const std::vector<double> rgrid);
  radial_HO_wave_functions (const int n_re, const int l_re, const double anu, const double r);
};
#endif // HO_H
