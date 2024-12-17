#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <cmath>
#include "constants.hpp"

namespace potentials {

  double separableGaussian(const double& pp, const double& p) {
    // For p',p in fm^{-1}, return v(p',p) in fm^2. For Lambda -> oo, it
    // should approach zero-range contact interaction.
    constexpr double gamma {-1}, // fm^2
      Lambda {500 / constants::hbc}; // Regulator cutoff momentum, 1/fm
    auto g = [&Lambda] (const double& q)->double{return std::exp(-q*q / (Lambda*Lambda)); };
    // Sharper cutoff
    // alpha = 1
    // g = lambda q: np.exp(-q**(2*alpha) / Lambda**(2*alpha))
    // g = lambda q: 1 / ((2 * np.pi)**(3/2) * Lambda**3) * np.exp(-q**2 / (2 * Lambda**2)) # g(x;Lambda) -> delta^{3}(\vec{x}) for Lambda -> 0
    return gamma * g(pp) * g(p);
  }
  
  double separable_Yamaguchi(const double& pp, const double& p) {
    // For p',p in fm^{-1}, return v(p',p) in fm^2.
    constexpr double gamma {-1}, // fm^2
      Lambda {500 / constants::hbc}; // Regulator cutoff momentum in 1/fm
    auto g = [&Lambda] (const double& q)->double{ return Lambda*Lambda / (q*q + Lambda*Lambda); };
    return gamma * g(pp) * g(p);
  }

  double spherical_square_well(const double V0, const double R0, const double r) {
    return (r<R0) ? V0 : 0.0;
  }

  double Minnesota(const int& s, const int& t, const double& r){
    const double
      rmu1 {1.487},
      rmu2 {0.639},
      rmu3 {0.465},
      VMn1 {200.0},
      VMn2 {-178.0},
      VMn3 {-91.85};
    const double rr {r*r};
    const double
      sgns {std::pow(-1,s)},
      sgnt {std::pow(-1,t)},
      sgnspt {std::pow(-1,s+t)};
    return
      0.5 * VMn1 * std::exp(-rmu1*rr) * (1 - sgnspt)
      + 0.25 * VMn2 * std::exp(-rmu2*rr) * (1 + (-sgns + sgnt - sgnspt))
      + 0.25 * VMn3 * std::exp(-rmu3*rr) * (1 + sgns - sgnt - sgnspt);
  }
  
  double MalflietTjon(const double& r) {
    const double
      rmu1 {3.11},
      rmu2 {1.55},
      VMT1 {1458.05},
      VMT2 {-578.09};
    return VMT1 * std::exp(-rmu1*r)/r + VMT2 * std::exp(-rmu2*r)/r;
  }

  double harmonic_oscillator(const double& m, const double& hbo, const double& r) {
    using constants::hbc;
    return 1.0/2 * m * hbo*hbo / (hbc*hbc) * r*r; // in MeV for M, hbo in MeV and r in fm
  }

}
#endif
