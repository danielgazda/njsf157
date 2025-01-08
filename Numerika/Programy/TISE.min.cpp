#define OPTIM_ENABLE_EIGEN_WRAPPERS

#include "optim.hpp"
#include "numerov.hpp"
#include "bisection.hpp"
#include "constants.hpp"
#include "potentials.hpp"

void bound_state () {
  /*
    Bound state case, E < 0
    
    Simplified solution without matching - eigenvalue search by
    shooting at u(xCore)=xCore^{l+1} from xInf or
    u(xInf)=\exp{-\sqrt{2M(-E)}xInf} from xCore. Both seem to work
    fine.
  */
  using potentials::spherical_square_well;
  using constants::hbc;
  auto V = [](const double x){return spherical_square_well(-38.5, 1.93, r);};
  // auto V = [](const double x){return potentials::Minnesota(1,0,r);}; // Minnesota potential in the S=1, T=0 deuteron channel

  std::vector<double> xgrid, u; // To store the (reduced) radial wave function

  auto solve_by_minimization = [&V, &xgrid, &u] (const int& N,
						 const double& xCore,
						 const double& xInf,
						 const double& Emin,
						 const double& Emax) {
    
    auto tail = [&V, &xCore, &xInf, &N, &xgrid, &u] (const Eigen::VectorXd& Einp, Eigen::VectorXd* grad_out, void* opt_data) {
      constexpr double M {constants::Mp * constants::Mn / (constants::Mp + constants::Mn)};
      const int l {0};
      auto kk = [&M, &Einp, &V, &l] (const double& r) -> double {
	using constants::hbc;
	double E {Einp(0)};
	return - l*(l+1)/(r*r) - 2*M*V(r)/(hbc*hbc) + 2*M*E/(hbc*hbc);
      };
      auto zero = [](const double x){return 0.0;};
      // auto k = [&M] (const double E){return std::sqrt(2 * M * E / (hbc*hbc));};
      Numerov_solver nsOutw(xCore, xInf, N, kk, zero);
      auto uOutw = nsOutw.integrate_outwards(std::pow(nsOutw.grid[0], l+1),
					     std::pow(nsOutw.grid[1], l+1));

      // Integrate up to a very large radius to see the divergence
      Numerov_solver nsOutwLarge(xCore, 10*xInf, 10*N, kk, zero);
      auto uOutwLarge = nsOutwLarge.integrate_outwards(std::pow(nsOutw.grid[0], l+1),
					     std::pow(nsOutw.grid[1], l+1));
      
      xgrid.clear();
      u.clear();
      // for (const auto& val : uOutw) {u.push_back(val);}
      // for (const auto& val : nsOutw.grid) {xgrid.push_back(val);}
      for (const auto& val : uOutwLarge) {u.push_back(val);}
      for (const auto& val : nsOutwLarge.grid) {xgrid.push_back(val);}

      return uOutw.back()*uOutw.back(); // Minimize u(xInf)^2
    };

    Eigen::VectorXd E0(1);
    E0 << -2.1;// Initial energy
    optim::algo_settings_t settings;
    Eigen::VectorXd lb(1), ub(1);
    lb << Emin;
    ub << Emax;
    settings.vals_bound = true;
    settings.lower_bounds = lb;
    settings.upper_bounds = ub;
    bool success = optim::nm(E0, tail, nullptr, settings);
    if (success) {
      std::cout << E0 << "\n";
    } else {
      std::cout << "Optimization failed, abort()\n";
      abort();
    }
  };
    
  int N = 20'000; // Number of points for solution
  double Emin = -4.0; // Lower bound for energy eigenvalue search
  double Emax = -1.0; // Upper bound for energy eigenvalue search
  double xCore = 1.0e-9; // Lower bound for radius. Avoid the centrifugal singularity for l>0
  double xInf = 30.0; // Upper bound for radius.
  // std::cin >> N >> xInf;
  solve_by_minimization(N, xCore, xInf, Emin, Emax);

  // Print the solution
  std::cout << "Wave function (r, u(r)):\n";
  for (const auto& [i,uval] : std::views::enumerate(u)) {
    std::cout << xgrid[i] << " " << uval << "\n";
  }
  std::cout << std::endl;
}

int main () {
  bound_state();
  return 0;
}
