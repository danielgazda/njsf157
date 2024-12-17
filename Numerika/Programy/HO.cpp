#include <numbers>
#include <cmath>
//#include <algorithm>
#include <vector>
#include <iostream>
#include <ranges>
#include "HO.hpp"

double radial_HO_wave_functions::fgamal(const int arg) {
  switch (arg) {
  case 2:
    return 0;
  case 1:
    return 0.5*std::log(std::numbers::pi_v<double>);
  default:
    return std::log(static_cast<double>(arg)/2 - 1) + fgamal(arg-2);
  }
}

double radial_HO_wave_functions::fdsq(const int n, const int l) {
  return std::sqrt(n*(l+n+0.5));
}

// void gamasub(const int nmax, const int lmax) {
//   for (int arg{1}; arg<=maxgam; arg++) {
//     gamal.at(arg) = fgamal(arg);
//   }
//   for (int l{0}; l<=lmax; ++l) {
//     for (int n{1}; n<=nmax; ++n) {
//    dsq[n,l] = fdsq(n,l);
//     };
//   };
// }
 
radial_HO_wave_functions::u_type::u_type () : nmax(-1), lmax(-1), storage{} {}

radial_HO_wave_functions::u_type::u_type (const int nmax, const int lmax) : nmax(nmax), lmax(lmax), storage((nmax+1)*(lmax+1), 0.0) {}

double radial_HO_wave_functions::u_type::operator[] (const int n, const int l) const {
  return storage.at(n+(nmax+1)*l);
}

double& radial_HO_wave_functions::u_type::operator[] (const int n, const int l) {
  return storage.at(n+(nmax+1)*l);
}

double radial_HO_wave_functions::u_type::operator() (const int n, const int l) const {
  return storage.at(n+(nmax+1)*l);
}

double& radial_HO_wave_functions::u_type::operator() (const int n, const int l) {
  return storage.at(n+(nmax+1)*l);
}

radial_HO_wave_functions::ur_type::ur_type () : nmax(-1), lmax(-1), rgrid{}, storage{} {}

radial_HO_wave_functions::ur_type::ur_type (const int nmax, const int lmax, std::vector<double> rgrid)
  : nmax(nmax), lmax(lmax), rgrid(rgrid), storage((nmax+1)*(lmax+1)*rgrid.size(), 0.0) {}

double radial_HO_wave_functions::ur_type::operator[] (const int n, const int l, const int ir) const {
  return storage.at(ir + n*rgrid.size() + l*rgrid.size()*(nmax+1));
}

double& radial_HO_wave_functions::ur_type::operator[] (const int n, const int l, const int ir) {
  return storage.at(ir + n*rgrid.size() + l*rgrid.size()*(nmax+1));
}

double radial_HO_wave_functions::ur_type::operator() (const int n, const int l, const int ir) const {
  return storage.at(ir + n*rgrid.size() + l*rgrid.size()*(nmax+1));
}
double& radial_HO_wave_functions::ur_type::operator() (const int n, const int l, const int ir) {
  return storage.at(ir + n*rgrid.size() + l*rgrid.size()*(nmax+1));
}

radial_HO_wave_functions::radial_HO_wave_functions (const int n_re, const int l_re, const double anu, const std::vector<double> rgrid)
  : nmax(n_re), lmax(l_re), maxgam(2*l_re+5+1),
    // gamal(maxgam),
    // dsq(n_re,l_re),
    u(),
    ur(n_re,l_re,rgrid)
{
  for (const auto& [ir,r] : std::views::enumerate(rgrid)) {
    // TODO FIX: BAD - the code is just copied from below
    std::vector<long double> laguer_sav_qp(n_re+1, 0.0l), scale(n_re+1, 1.0l);
    //gamasub(n_re,l_re);
    for (int l{0}; l<=l_re; ++l) {
      double dlanu{std::log(anu)}, zz{anu*r*r};
      double wavel = 0.25*dlanu - zz/2.0 + (l+1) * (0.5 * dlanu + std::log(r));

      long double guerpq = std::exp(0.5*(std::log(2.0)-fgamal(2*l+3)));
      // guerpq = std::exp(0.5*(std::log(2.0)-gamal[2*l+3]));
      laguer_sav_qp.at(0) = guerpq;

      long double zzq = static_cast<long double>(zz);
      if (n_re>0) { // I adedd this
	guerpq = std::exp(0.5l*(std::log(2.0l) - fgamal(2*l+5))) * (static_cast<long double>(l) + 1.5l - zzq);
	// guerpq = std::exp(0.5l*(std::log(2.0l) - gamal[2*l+5])) * (static_cast<long double>(l) + 1.5l - zzq);
	laguer_sav_qp.at(1) = guerpq;
      }
     
      long double aq = std::exp(0.5l*(std::log(2.0l)-fgamal(2*l+3)));
      // aq = std::exp(0.5l*(std::log(2.0l)-gamal[2*l+3]));
   
      for (int nnn{2}; nnn<=n_re; ++nnn) {
	if (std::abs(aq)>1.0e+290l || abs(guerpq)>1.0e+290l) {
	  for (int i{nnn}; i<=n_re; ++i) {scale.at(i) *= 1.0e+290l;}
	  aq = aq * 1.0e-290l;
	  guerpq = guerpq * 1.0e-290l;
	}
	long double bq = (static_cast<long double>(l+2*nnn) - 0.5l - zzq) / fdsq(nnn,l) * guerpq - fdsq(nnn-1,l) / fdsq(nnn,l) * aq;
	// bq = (static_cast<long double>(l+2*nnn) - 0.5l - zzq) / dsq[nnn,l] * guerpq - dsq[nnn-1,l] / dsq[nnn,l] * aq;
	aq = guerpq;
	guerpq = bq;
	laguer_sav_qp.at(nnn) = guerpq;
      };
     
      for (int nnn{0}; nnn<=n_re; ++nnn) {
	guerpq = laguer_sav_qp.at(nnn);
	long double sig{(guerpq>=0.0l) ? 1.0l : -1.0l};
	guerpq = std::log(std::abs(guerpq)) + std::log(scale.at(nnn));
	ur[nnn,l,ir] = sig*std::exp(wavel + guerpq);
      };
    } 
  }
}
  
radial_HO_wave_functions::radial_HO_wave_functions (int n_re, int l_re, double anu, double r)
  : nmax(n_re), lmax(l_re), maxgam(2*l_re+5+1),
    // gamal(maxgam),
    // dsq(n_re,l_re),
    u(n_re,l_re), ur()
{
  std::vector<long double> laguer_sav_qp(n_re+1, 0.0l), scale(n_re+1, 1.0l);
  //gamasub(n_re,l_re);
  for (int l{0}; l<=l_re; ++l) {
    double dlanu{std::log(anu)}, zz{anu*r*r};
    double wavel = 0.25*dlanu - zz/2.0 + (l+1) * (0.5 * dlanu + std::log(r));
    
    long double guerpq = std::exp(0.5*(std::log(2.0)-fgamal(2*l+3)));
    // guerpq = std::exp(0.5*(std::log(2.0)-gamal[2*l+3]));
    laguer_sav_qp.at(0) = guerpq;

    long double zzq = static_cast<long double>(zz);
    if (n_re>0) { // I adedd this
      guerpq = std::exp(0.5l*(std::log(2.0l) - fgamal(2*l+5))) * (static_cast<long double>(l) + 1.5l - zzq);
      // guerpq = std::exp(0.5l*(std::log(2.0l) - gamal[2*l+5])) * (static_cast<long double>(l) + 1.5l - zzq);
      laguer_sav_qp.at(1) = guerpq;
    }
     
    long double aq = std::exp(0.5l*(std::log(2.0l)-fgamal(2*l+3)));
    // aq = std::exp(0.5l*(std::log(2.0l)-gamal[2*l+3]));
   
    for (int nnn{2}; nnn<=n_re; ++nnn) {
      if (std::abs(aq)>1.0e+290l || abs(guerpq)>1.0e+290l) {
	for (int i{nnn}; i<=n_re; ++i) {scale.at(i) *= 1.0e+290l;}
	aq = aq * 1.0e-290l;
	guerpq = guerpq * 1.0e-290l;
      }
      long double bq = (static_cast<long double>(l+2*nnn) - 0.5l - zzq) / fdsq(nnn,l) * guerpq - fdsq(nnn-1,l) / fdsq(nnn,l) * aq;
      // bq = (static_cast<long double>(l+2*nnn) - 0.5l - zzq) / dsq[nnn,l] * guerpq - dsq[nnn-1,l] / dsq[nnn,l] * aq;
      aq = guerpq;
      guerpq = bq;
      laguer_sav_qp.at(nnn) = guerpq;
    };
     
    for (int nnn{0}; nnn<=n_re; ++nnn) {
      guerpq = laguer_sav_qp.at(nnn);
      long double sig{(guerpq>=0.0l) ? 1.0l : -1.0l};
      guerpq = std::log(std::abs(guerpq)) + std::log(scale.at(nnn));
      u[nnn,l] = sig*std::exp(wavel + guerpq);
    };
  }
}
