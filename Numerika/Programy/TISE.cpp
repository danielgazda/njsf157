#include "numerov.hpp"
#include "bisection.hpp"
#include "constants.hpp"
#include "potentials.hpp"

void bound_state () {
  /*
    Bound state case, E < 0
    
    Simplified solution without matching - eigenvalue search by
    shooting at u(rCore)=rCore^{l+1} from rInf or
    u(rInf)=\exp{-\sqrt{2M(-E)}rInf} from rCore. Both seem to work
    fine.
  */
  using potentials::spherical_square_well;
  using constants::hbc;
  auto V = [](const double r){return spherical_square_well(-38.5, 1.93, r);};
  // auto V = [](const double r){return potentials::Minnesota(1,0,r);}; // Minnesota potential in the S=1, T=0 deuteron channel

  std::vector<double> rgrid, u; // To store the (reduced) radial wave function

  auto solve_by_shooting = [&V, &rgrid, &u] (const int& N,
					     const double& rCore,
					     const double& rInf,
					     const double& Emin,
					     const double& Emax) {

    auto shooting = [&V, &N, &rCore, &rInf, &rgrid, &u] (const double& E) {
      constexpr double M {constants::Mp * constants::Mn / (constants::Mp + constants::Mn)}; // = M*c^2, (reduced) mass in MeV/cc
      const int l {0}; // S-wave partial wave
// k^2(x)

      // Define k^2(x)
      auto kk = [&M, &E, &V, &l] (const double& r) -> double {
	using constants::hbc;
	return - l*(l+1)/(r*r) - 2*M*V(r)/(hbc*hbc) + 2*M*E/(hbc*hbc);
      }; // 1/fm^2

      // Define r.h.s. - F(x)
      auto zero = [](const double r){return 0.0;};

      // k - E, momentum<-energy relation
      auto k = [&M] (const double E){return std::sqrt(2 * M * E / (hbc*hbc));}; // 1/fm

      // Integrate outwards from rCore to rInf
      Numerov_solver nsOutw(rCore, rInf, N, kk, zero);
      auto uOutw = nsOutw.integrate_outwards(std::pow(nsOutw.grid[0], l+1), // It should not matter, u(0)=0, asymptotics is u = A*r^(l+1)
					     std::pow(nsOutw.grid[1], l+1)); // It should not matter
      // Integrate inwards from rInf to rCore
      // Numerov_solver nsInw(rCore, rInf, N, kk, zero);
      // auto uInw = nsInw.integrate_inwards(std::exp(-k(-E)*nsInw.grid.back()), // Asymptotics is u ~ A * e^{-k*r}, say A=1, E is negative
      // 					  std::exp(-k(-E)*nsInw.grid.rbegin()[1]));

      
      // double u_rCore {uInw.at(0)};
      double u_rInf {uOutw.rbegin()[1]};

      // Store the solution
      rgrid.clear();
      u.clear();
      // Fill by the outward solution
      for (const auto& val : uOutw) {u.push_back(val);}
      for (const auto& val : nsOutw.grid) {rgrid.push_back(val);}
      // Fill by the inward solution
      // for (const auto& val : uInw) {u.push_back(val);}
      // for (const auto& val : nsInw.grid) {rgrid.push_back(val);}
      
      // return u_rCore - std::pow(rCore,l+1);
      return u_rInf - std::exp(-k(-E)*rInf);
    };
      
    double Eb = find_root<100, 1.0e-9, 1.0e-9, true>(Emin, Emax, shooting);
    std::cout << "E=" << Eb << "\n";
  };

  int N = 4000; // Number of points for solution
  double Emin = -4.0; // Lower bound for energy eigenvalue search
  double Emax = -1.0; // Upper bound for energy eigenvalue search
  double rCore = 1.0e-9; // Lower bound for radius. Avoid the centrifugal singularity for l>0
  double rInf = 22.0; // Upper bound for radius.
  solve_by_shooting(N, rCore, rInf, Emin, Emax);

  // Print the solution
  std::cout << "Wave function (r, u(r)):\n";
  for (const auto& [i,uval] : std::views::enumerate(u)) {
    std::cout << rgrid[i] << " " << uval << "\n";
  }
  std::cout << std::endl;
}

void analytic_radial_square_well_bound_state () {
  /*
    Semi-analytical bound-state solution for radial square well
  */
  double M {constants::Mp * constants::Mn / (constants::Mp + constants::Mn)};
  double V0 {38.5}, // MeV
    a {1.93}; // fm
  auto f = [&M, &V0, &a] (const double& E) {
    using constants::hbc;
    double k {std::sqrt(2*M * (V0 - std::abs(E)) )/ hbc}; // 1/fm
    double kappa {std::sqrt((2*M*V0)/(hbc*hbc) - k*k)}; // 1/fm
    return kappa + k * 1/std::tan(k*a);
  };
  double Eanalytical = find_root<100,1.0e-9,1.0e-9,true>(-10.0, -0.1, f);
  std::cout << "analytic_radial_square_well_bound_state():\n";
  std::cout << "E analytical = " << Eanalytical << "\n\n";
}

double tan_phase_shift(const double E, double rCore, double rInf, int N, int l, double approx_r1, double approx_r2) {
  /*
    Scattering case, E > 0
  */
  if (E<0) {std::cout << "E<0, abort()\n"; abort();}

  const double M {constants::Mp * constants::Mn / (constants::Mp + constants::Mn)};
  using potentials::spherical_square_well;
  using constants::hbc;
  auto V = [](const double r){return spherical_square_well(-38.5, 1.93, r);}; // pheno np 3S1
  // auto V = [](const double r){return spherical_square_well(-14.3, 2.50, r);}; // pheno np 1S0
  // auto V = [](const double r){return potentials::Minnesota(1,0, r);}; // Minnesota potential in the S=1, T=0 deuteron channel

  auto kk = [&M, &E, &V, &l] (const double& r) -> double {
    return - l*(l+1)/(r*r) - 2*M*V(r)/(hbc*hbc) + 2*M*E/(hbc*hbc);
  }; // 1/fm^2
  auto k = [&M] (const double E){return std::sqrt(2 * M * E / (hbc*hbc));}; // 1/fm
  auto zero = [](const double& x){return 0.0;};

  Numerov_solver ns(rCore, rInf, N, kk, zero);
  // Log-derivative of the asymptotic, r->oo, solution
  // auto uInf = [](const double k, const double r, const double delta) {
  //   return std::cos(delta) * k*r * std::sph_bessel(k*r) - std::sin(delta) * k*r * std::sph_neumann(k*r);
  // }
  auto u = ns.integrate_outwards(std::pow(ns.grid[0], l+1), std::pow(ns.grid[1], l+1));

  // Find grid values and indices closest to requested r1, r2
  int i1 {ns.grid_index_closest_to(approx_r1)},
    i2 {ns.grid_index_closest_to(approx_r2)};
  double r1 {ns.grid.at(i1)}, r2 {ns.grid.at(i2)};
  // std::cout << "Matching at r1, r2 = " << r1 << ", " << r2 << "\n";

  double beta {r1*u[i2] / (r2*u[i1])};
  double tandelta {
    (beta * std::sph_bessel(l, k(E)*r1) - std::sph_bessel(l, k(E)*r2))
    / (beta * std::sph_neumann(l, k(E)*r1) - std::sph_neumann(l, k(E)*r2))
  };
  return tandelta; // Return tan(phase shift)
}

double analytic_radial_square_well_scatterring(const double E) {
  /*
    Analytical scattering (phase shift) solution for radial square well
   */
  const double M {constants::Mp * constants::Mn / (constants::Mp + constants::Mn)};
  const double V0 {38.5}, a {1.93}; // Pheno np 3S1
  // const double V0 {14.3}, a {2.50}; // Pheno np 3S1
  using constants::hbc;
  double k {
    std::sqrt(2 * M * E / (hbc*hbc))
  };
  double alpha {
    std::sqrt(2 * M * (E+V0) / (hbc*hbc))
  };
  using std::sph_bessel;
  using std::sph_neumann;
  const int l {0};
  double tandelta {
    (k*sph_bessel(l+1,k*a)*sph_bessel(l,alpha*a) - alpha*sph_bessel(l,k*a)*sph_bessel(l+1,alpha*a))
    / (k*sph_neumann(l+1,k*a)*sph_bessel(l,alpha*a) - alpha*sph_neumann(l,k*a)*sph_bessel(l+1,alpha*a))
  };
  return tandelta;
}
  

int main () {

  analytic_radial_square_well_bound_state();
  bound_state();

  for (double E {0.1}; E<10.0; E+=+0.2){
    double tandelta = tan_phase_shift(E,      // C-o-m energy in MeV
				      1.0e-9, // rCore in fm
				      200.0,  // rInf in fm
				      10'000,   // N, number of intervals
				      0,      // l,  S-wave
				      180,    // approximate r1 radius
				      190     // approximate r2 radius
				      );
    double analytic_tandelta = analytic_radial_square_well_scatterring(E);
    // std::cout << E << " " << std::atan(tandelta) << "\n";
    std::cout << E << " " << std::atan(tandelta) << " " << std::atan(analytic_tandelta) << " diff=" << std::atan(tandelta)-std::atan(analytic_tandelta) << "\n";
    // std::cout << E << " " << std::atan(tandelta) << "\n";
  }

  return 0;
}
