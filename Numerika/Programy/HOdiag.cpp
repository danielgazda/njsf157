#include <cmath>
#include <vector>
#include <ranges>
#include <algorithm>
#include <functional>
#include <iostream>

#include "Gauss_Legendre.hpp"
#include "constants.hpp"
#include "potentials.hpp"
#include "HO.hpp"
#include "matrix.hpp"

#include <Eigen/Dense>

int delta(int i, int j) {return (i==j)? 1 : 0;}

class HOBasis {
public:
  int Nmax {};
  int dim {};
  int parity {};
  // Different storage schemes
  std::vector<int> n {}, l {};
  std::vector<std::pair<int,int>> nl{};
  std::vector<std::tuple<int,int,int>> inl{};

  HOBasis (const int Nmax, const int pi) : Nmax(Nmax), parity(pi) {
    int ind {0};
    const int NNmin {(parity==-1) ? 1 : 0};
    // Check if we have odd Nmax for negative parity
    if (parity == -1 && Nmax % 2 == 0) {
      std::cout << "parity = " << parity << ", abort()";
      abort();
    }
    const int dNN {(parity==-1 || parity==1)? 2 : 1}; // Allows to build a basis with mixed parity
    for (int NN {NNmin}; NN <= Nmax; NN += dNN) {
      for (int _n {0}; _n<=Nmax/2; ++_n) {
	for (int _l {0}; _l<=Nmax; ++_l) {
	  if (_l != 0) continue; // Deuteron without tensor force, L=0 s-wave, no s-d coupling
	  if (2*_n + _l == NN) {
	    n.push_back(_n);
	    l.push_back(_l);
	    // std::vector<int> _nl {_n,_l};
	    std::pair<int,int> _nl {_n, _l};
	    nl.push_back(_nl);
	    std::tuple<int, int, int> _inl {ind, _n, _l};
	    inl.push_back(_inl);
	    ++ind;
	  }
	}
      }
    }
    dim = nl.size();
  }
  
  void print() {
    std::cout << "HO states (idex,n,l):\n";
    for (auto [i,n,l] : inl) {
      std::cout << i << ": " << n << " " << l << "\n";
    }
  }
};

void findSingleChannelBoundStateHO(const int& Ntotmax, const double& hbo) {

  // Kinetic energy matrix elements
  auto T = [&hbo] (const int na, const int nb, const int l)->double{
    if (na==nb) {
      return hbo/2 * (2*na + l + 1.5);
    } else if (na==nb-1) {
      return hbo/2 * std::sqrt((na + 1)*(na + l + 1.5));
    } else if (na==nb+1) {
      return hbo/2 * std::sqrt((nb + 1)*(nb + l + 1.5));
    } else {
      return 0.0;
    }
  };

  // Construct 2-body HO basis
  HOBasis sts2b( Ntotmax, +1);
  std::cout << "Basis dimension = " << sts2b.dim << "\n";
  // sts2b.print();

  using constants::Mp, constants::Mn, constants::hbc;
  double M {Mp*Mn/(Mp+Mn)};

#ifdef MOMENTUMSPACE
  /*
    For momentum-space potential
  */

  // Momentum grid
  double
    kCore {6.0}, // 1/fm
    kInf {0.0}; // 1/fm
  int kN {100};
  Gauss_Legendre_mesh kmesh(kN, kInf, kCore);
  auto GLikw = std::views::zip(std::views::iota(0,kN), kmesh.x, kmesh.w);
  // Oscillator length
  double bb {M * hbo / (hbc*hbc)}; // In momentum space, b is inverted and there is additional (-1)^n * u_{n,l}(r)
  std::cout << "HO length scale b = sqrt(M*omega/hbar) = " << std::sqrt(bb) << "\n";
  radial_HO_wave_functions HOk(Ntotmax/2, Ntotmax, 1/bb, kmesh.x);
#else
  /*
    For coordinate-space potential (default)
  */

  // Radial coordinate grid
  const double
    rCore {1.0e-9}, // fm
    rInf {22.0}; // fm
  const int rN {4*200}; // 200 was not enough for square radial well
  // Use finite grid
  std::cout << "Integrate, " << rN << " points, [" << rCore << ", " << rInf << "] fm\n";
  Gauss_Legendre_mesh rmesh(rN, rCore, rInf);
  auto GLirw = std::views::zip(std::views::iota(0,rN), rmesh.x, rmesh.w);

  // nu = (Ma + Mb) / (Ma * Mb) / hbo * hbarc^2
  // = b^2, with b the HO length scale
  double bb  = hbc*hbc / M / hbo; // M in MeV, bb in fm^2
  std::cout << "HO length scale b = sqrt(hbar/(M*omega)) = " << std::sqrt(bb) << " fm \n";
  radial_HO_wave_functions HOr(Ntotmax/2, Ntotmax, 1/bb, rmesh.x);
  // I need the factor (-1)^n * u(n,l,r) in momentum space
#endif

  auto v = [](const double r){return potentials::spherical_square_well(-38.5, 1.93, r);}; // pheno np 3S1
  // auto v = [](const double r){return potentials::spherical_square_well(-14.3, 2.50, r);}; // pheno np 1S0
  // Use Minnesota potential to get deuteron
  // auto v = [](const double& r){return potentials::Minnesota( 1, 0, r);}; // Minnesota in the deuteron channel

#ifdef MOMENTUMSPACE
  // Radial Fourier transform of v(r) to integrate in momentum space - testing only
  // It is very slow, use small Ntotmax
  const double rCore {1.0e-9}, rInf {22.0};
  const int rN {200};
  Gauss_Legendre_mesh rmesh(rN, rCore, rInf);
  auto FTv = [&v, &rmesh](const double& ka, const double& kb)->double{
    const int l {0}; // s-wave only
    double I {0.0};
    for (auto [r,w] : std::views::zip( rmesh.x, rmesh.w)) {
      I += w * r*r
	* 2/constants::PI
	* ka * std::sph_bessel( l, ka*r)
	* v(r)
	* kb * std::sph_bessel( l, kb*r);
    }
    return I;
  };
#endif

  // Potential matrix elements in HO basis

  Matrix V(sts2b.dim);
  for (auto& [ia,na,la] : sts2b.inl) {
    for (auto& [ib,nb,lb] : sts2b.inl) {
      if (ib>ia) continue; // Lower triangle only, V must be symmetric

      double I {0.0};
#ifdef MOMENTUMSPACE
      // Integrate over the momentum grid
      for (auto [ikwa,ka,wa] : GLikw) {
	for (auto [ikwb,kb,wb] : GLikw) {
	  // u(ipa,na,la)*w(ipa)*dV(ipa,ipb)*w(ipb)*u(ipb,nb,lb) // With discretized V
	  I += HOk.ur[na,la,ikwa] * wa
	    * FTv(ka,kb)
	    * wb * HOk.ur(nb,lb,ikwb);
	}
      }
      I *= std::pow(-1,na)*std::pow(-1,nb); // Phase factors from the radial wave functions
#else
      // Integrate over the coordinate grid
      for (auto [i,r,w] : GLirw) {
	I += w
	  * HOr.ur(na,la,i)
	  * v(r)
	  * HOr.ur(nb,lb,i);
      }
#endif
      V(ia,ib) = I;
      V(ib,ia) = I;
    }
  }

  // Hamiltonian matrix
  Matrix H(sts2b.dim);
  for (auto& [ia,na,la] : sts2b.inl) {
    for (auto& [ib,nb,lb] : sts2b.inl) {
      if (ib>ia) continue; // Lower triangle only
      double t {(la==lb) ? T(na,nb,la) : 0.0};
      H(ia,ib) = t + V(ia,ib);
      H(ib,ia) = H(ia,ib); // H is real and symmetric
    }
  }

  // Diagonalize H
  auto X = Eigen::Map<Eigen::MatrixXd>(H.storage.data(), sts2b.dim, sts2b.dim);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(X);
  const int nev {3}; // Number of eigenvalues to print
  std::cout << "The smallest eigenvalues (" << std::min(sts2b.dim, nev) << ")\n";
  for (int i : std::views::iota(0,std::min(sts2b.dim,nev))) {
    std::cout << eigensolver.eigenvalues()(i) << "\n";
  }
}

int main () {
  {
    double hbo {20.0}; // MeV, hbar omega
    for (int Ntotmax {20}; Ntotmax<=100; Ntotmax+=4) {
      std::cout << "Ntotmax=" << Ntotmax << "\n";
      findSingleChannelBoundStateHO(Ntotmax, hbo);
      std::cout << std::endl;
    }
  }
  return 0;
}
