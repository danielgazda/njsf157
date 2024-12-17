#include <iostream>
#include <cmath>
#include <ranges>
#include <vector>

#include "constants.hpp"
#include "potentials.hpp"
#include "RK.hpp"

double phase_shift_vpa (const double& E, const bool& verbose) {
  /*
    Solve
    
    d/dr \delta_l(k,r) = - V(r) / (k * \hbar^2 / (2 \mu)) * [\cos(\delta_l) \hat{j}_l(kr) - \sin(\delta_l) \hat{n}_l(kr)]^2
    Taken from https://arxiv.org/pdf/2403.19173, quite a random reference

    Takes energy E>0 in MeV,
    verbose = true turns on more printing
  */

  if (E<0) {
    std::cout << "E=" << E << ", E<0 does not make sense here, abort()";
    abort();
  }

  // Define potential
  auto V = [](const double& r){return potentials::spherical_square_well(-38.5, 1.93, r);}; // Pheno np in 3S1 channel
  // auto V = [](const double r){return potentials::spherical_square_well(-14.3, 2.50, r);}; // Pheno np in 1S0 channel
  // auto V = [](const double r) {return potentials::Minnesota(1,0,r);}; // 2H deuteron 3S1

  // Define r.h.s.
  using constants::Mp;
  using constants::Mn;
  double mu {Mn*Mp/(Mn+Mp)}; // Reduced mass
  double k {std::sqrt(2*mu*E)/constants::hbc}; // momentum 1/fm
  int l {0}; // s-wave
  if (verbose) {
    std::cout << "E=" << E << " MeV\n"
	      << "k=" << k << " 1/fm\n"
	      << "mu=" << mu << " MeV\n"
	      << "l=" << l << "\n";
  }
  auto hatj = [](const int l, const double z){return z*std::sph_bessel(l,z);}; // \hat{j}(z)
  auto hatn = [](const int l, const double z){return z*std::sph_neumann(l,z);}; // \hat{n}(z)
  auto rhs = [&V, &mu, &k, &hatj, &hatn, &l] (const double r, const double delta) {
    return -V(r) / (k* constants::hbc*constants::hbc / (2*mu))
      * std::pow((std::cos(delta) * hatj(l,k*r)
		  - std::sin(delta) * hatn(l,k*r)), 2);
  };

  // Solve the equation
  constexpr double rCore {1.0e-16}, rMax {3.0};
  constexpr int N {1000};
  constexpr double delta0 {0.0};
  RK_solver rk(rCore, // a, min radius
	       rMax,  // b, max radius - Must be larger than the range of the interaction
	       N,     // N, number of intervals
	       rhs,   // r.h.s
	       delta0 // initial value, delta_l(k,r) = 0 for r = 0 and all l,k
	       );
  /*
    We should look at large radius, where the potential vanishes, to
    extract delta_l(k) from delta_l(r,k). The solution should converge
    as delta_l(r,k) -> delta_l(k) for r >> range of the interaction.
  */

  if (verbose) {
    // Print convergence with r
    std::cout << "( r, delta_l(r,k)) for convergence analysis:\n";
    for (auto [r,d] : std::views::zip(rk.x,rk.y)) {
      std::cout << r << " " << d << "\n";
    }
    std::cout << std::endl;
  }
  
  // Return phase shift at the largest radius
  return rk.y.back();
}

int main () {
  constexpr double E {0.1}; // MeV
  constexpr bool verbose {true};
  double delta = phase_shift_vpa(E,verbose);
  std::cout << "Phase shift tan(delta_l(k))=" << std::tan(delta) << " for E=" << E << " MeV" << std::endl;
  return 0;
}
