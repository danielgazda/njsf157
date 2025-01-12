#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

#include "numerov.hpp"
#include "bisection.hpp"
#include "constants.hpp"
#include "potentials.hpp"

#include <iostream>
#include <vector>
#include <fstream>

void bound_state () {
  /*
    Bound state case, E < 0
    
    Eigenvalue search by integrating from xCore, u(xCore)=xCore^{l+1},
    and from xInf, u(xInf)=\exp{-\sqrt{2M(-E)}xInf}, and matching the
    derivative at xMatch.
  */
  using potentials::spherical_square_well;
  using constants::hbc;
  auto V = [](const double x){return spherical_square_well(-38.5, 1.93, x);};
  // auto V = [](const double x){return potentials::Minnesota(1,0,r);}; // Minnesota potential in the S=1, T=0 deuteron channel

  std::vector<double> xgrid, u; // To store the (reduced) radial wave function

  auto solve_by_matching = [&V, &xgrid, &u] (const int& N,
						 const double& xCore,
						 const double& xInf,
						 const double& Emin,
						 const double& Emax) {
    
    auto derivatives_diff = [&V, &xCore, &xInf, &N, &xgrid, &u] (const Eigen::VectorXd& Einp, void* opt_data) {
      constexpr double M {constants::Mp * constants::Mn / (constants::Mp + constants::Mn)};
      const int l {0};
      auto kk = [&M, &Einp, &V, &l] (const double& r) -> double {
	using constants::hbc;
	double E {Einp(0)};
	return - l*(l+1)/(r*r) - 2*M*V(r)/(hbc*hbc) + 2*M*E/(hbc*hbc);
      };
      auto zero = [](const double x){return 0.0;};
      auto k = [&M] (const double E){return std::sqrt(2 * M * (-E) / (hbc*hbc));};

      Numerov_solver nsInw(xCore, xInf, N, kk, zero);
      auto uInw = nsInw.integrate_inwards(std::exp( -k(Einp[0]) * nsInw.grid.at(N)),
					  std::exp( -k(Einp[0]) * nsInw.grid.at(N-1)));
      
      // Find turning point - index
      int iMatch {-1};
      for (int i {N}; i>=1; --i) {
	if (uInw.at(i-1) < uInw.at(i)) {
	  iMatch = i;
	  std::cout << "u(x) turns down at x=" << nsInw.grid.at(i) << "=xMatch (i=" <<  iMatch << ")\n";
	  break;
	}
      }
      if (iMatch == -1) {
	std::cout << "Matching radius not found, abort()\n";
	abort();
      }

      Numerov_solver nsOutw(xCore, xInf, N, kk, zero);
      auto uOutw = nsOutw.integrate_outwards(std::pow(nsOutw.grid[0], l+1),
					     std::pow(nsOutw.grid[1], l+1));

      // Stich the solutions at xMatch (to be continous)
      {
	double C {uOutw[iMatch]/uInw[iMatch]};
	for (auto& u: uInw) {u *= C;}
      }

      // Integrate up to a very large radius to see the divergence
      // Numerov_solver nsOutwLarge(xCore, 10*xInf, 10*N, kk, zero);
      // auto uOutwLarge = nsOutwLarge.integrate_outwards(std::pow(nsOutw.grid[0], l+1),
      // 					     std::pow(nsOutw.grid[1], l+1));
      
      xgrid.clear();
      u.clear();
      for (int i {0}; i<iMatch; ++i) {
	xgrid.push_back(nsOutw.grid[i]);
	u.push_back(uOutw[i]);
      }
      for (int i {iMatch}; i<N; ++i) {
	xgrid.push_back(nsInw.grid[i]);
	u.push_back(uInw[i]);
      }
      // for (const auto& val : uOutw) {u.push_back(val);}
      // for (const auto& val : nsOutw.grid) {xgrid.push_back(val);}
      // for (const auto& val : uInw) {u.push_back(val);}
      // for (const auto& val : nsInw.grid) {xgrid.push_back(val);}
      // for (const auto& val : uOutwLarge) {u.push_back(val);}
      // for (const auto& val : nsOutwLarge.grid) {xgrid.push_back(val);}

      Eigen::VectorXd res(1);
      res <<  uInw[iMatch-1] + uOutw[iMatch+1] -2*uOutw[iMatch];
      return res;
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
    settings.rel_objfn_change_tol = 1e-6;
    settings.rel_sol_change_tol = 1e-6;
    settings.print_level = 1;
    bool success = optim::broyden_df(E0, derivatives_diff, nullptr, settings);
    if (success) {
      std::cout << "E=" << E0 << " MeV\n";
    } else {
      std::cout << "Optimization failed, abort()\n";
      abort();
    }
  };
    
  int N = 5'000; // Number of points for solution
  double Emin = -4; // Lower bound for energy eigenvalue search
  double Emax = -1; // Upper bound for energy eigenvalue search
  double xCore = 1e-9; // Lower bound for radius. Avoid the centrifugal singularity for l>0
  double xInf = 40; // Upper bound for radius.
  std::cin >> N >> xCore >> xInf;
  solve_by_matching(N, xCore, xInf, Emin, Emax);

  // Print the solution
  std::ofstream wf("TISE.match.wf.txt");
  wf << "# Wave function (r, u(r)):\n";
  for (const auto& [xval,uval] : std::views::zip(xgrid,u)) {
    wf << xval << " " << uval << "\n";
  }
  wf << std::endl;
}

int main () {
  bound_state();
  return 0;
}
