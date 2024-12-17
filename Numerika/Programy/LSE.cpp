#include <cmath>
#include <vector>
#include <ranges>
#include <algorithm>
#include <functional>
// #include <numeric> // provides std::reduce()
#include <iostream>

#include "Gauss_Legendre.hpp"
#include "constants.hpp"
#include "potentials.hpp"
#include "HO.hpp"
#include "bisection.hpp"
#include "matrix.hpp"

#include <Eigen/Dense>

int delta(int i, int j) {return (i==j)? 1 : 0;}

void findSingleChannelBoundState () {
  // Momentum grid
  constexpr double
    kCore {6.0}, // 1/fm
    kInf {1.0e-9}; // 1/fm
  constexpr int N {2*200};
  // Use finite grid
  Gauss_Legendre_mesh mesh(N, kInf, kCore);
  std::vector<double> k(mesh.x);
  std::vector<double> omega(mesh.w);

  // Test with Fourier-transformed Minnesota or square well
  using constants::Mp, constants::Mn, constants::hbc;
  double M {2 * Mp * Mn / (Mp + Mn)}; // = 2 * reduced mass! (To be consistent with Haftel-Tabakhin)

  auto Vr = [](const double r){return potentials::spherical_square_well(-38.5, 1.93, r);}; // pheno np 3S1
  // auto Vr = [](const double r){return potentials::spherical_square_well(-14.3, 2.50, r);}; // pheno np 1S0
  // auto Vr = [](const double& r)->double{return potentials::Minnesota(1,0,r);};

  // Fourier transform Vr (this is very slow!)
  constexpr double
    rCore {1.0e-9}, // fm
    rInf {22.0}; // fm
  constexpr int rN {4*200};
  Gauss_Legendre_mesh rmesh(rN, rCore, rInf);
  auto Vk = [&Vr, &rmesh, &M](const double& ka, const double& kb)->double{
    constexpr int l {0}; // s-wave only
    double I {0.0};
    for (auto [r,w] : std::views::zip( rmesh.x, rmesh.w)) {
      I += w * r*r
	* std::sph_bessel(l,ka*r)
	* Vr(r) / constants::hbc // Vr in MeV -> 1/fm
	* std::sph_bessel(l,kb*r) // MeV fm
	* M * constants::hbc / (constants::hbc*constants::hbc); // Because of the H.-T. potential partial waves definition
    }
    return I;
  };

  std::cout << "Doing ourier transform of V ...\n";
  Matrix V = discretize(Vk,k);
  
  auto FD = [&V, &k, &omega] (const double& kD) -> Matrix {
    Matrix X(N);
    for (int i {0}; i<N; ++i) {
      for (int j {0}; j<N; ++j) {
	X(i,j) = delta(i,j) + 2/constants::PI * k[j]*k[j] * omega[j] / (k[j]*k[j] + kD*kD) * V(i,j);
      }
    }
    return X;
  };

  auto detFD = [&FD] (const double& kD)->double{
    // Returns det FD(kD)
    return Eigen::Map<Eigen::MatrixXd>( FD(kD).storage.data(), N, N).determinant();
  };

  auto kFromE = [&M](const double& E)->double{
    // Converts energy to momentum
    return std::sqrt(-E*M/(constants::hbc*constants::hbc));
  };

  {
    // Search for bound states by det FD = 0
    double Emin {-3.0}, Emax {-1.0};
    double kD = find_root<>(kFromE(Emin), kFromE(Emax), detFD);
    auto EFromk = [&M](const double& k)->double{
      return -k*k * (constants::hbc*constants::hbc) / M;
    };
    std::cout << "Bound-state E=" << EFromk(kD) << " MeV\n";
  }
}

double computeRMatrix(const double& k0) {
  using constants::PI;
  constexpr double
    kCore {6.0}, // 1/fm
    kInf {0.0}; // 1/fm
  constexpr int N {150}; // Number of grid points for momentum discretization / integration
  // Integration grid points and weights
  // TODO Switch to (0,oo) interval by some transformation
  Gauss_Legendre_mesh mesh(N, kInf, kCore);
  auto k = Eigen::Map<Eigen::Array<double,N,1>>(mesh.x.data());
  auto omega = Eigen::Map<Eigen::Array<double,N,1>>(mesh.w.data());
  // Append k0 as (N+1)st point
  Eigen::ArrayXd kp(N+1,1), omegap(N+1,1);
  kp << k, k0;
  omegap(Eigen::seq(0,N-1)) = 2/PI * k*k * omega / (k*k - k0*k0);
  omegap(N) = -2/PI * k0*k0 * (omega / (k*k - k0*k0)).sum();
  // Potential
  Eigen::MatrixXd V(N+1,N+1);

  // for (int i : std::views::iota(0,N+1)) {
  //   for (int j : std::views::iota(0,N+1)) {
  //     V(i,j) = potentials::separableGaussian(kp(i),kp(j));
  //   }
  // }
  
  {
    // Test with Fourier-transformed Minnesota
    using constants::Mp, constants::Mn, constants::hbc;
    double M {2*Mp*Mn/(Mp+Mn)}; // = 2 * reduced mass
    auto Vr = [](const double& r)->double{return potentials::Minnesota(1,0,r);};
    // Fourier transform Vr
    // This is very slow
    const double
      rCore {1.0e-9}, // fm
      rInf {22.0}; // fm
    const int rN {150};
    Gauss_Legendre_mesh rmesh(rN, rCore, rInf);
    auto Vk = [&Vr, &rmesh, &M](const double& ka, const double& kb)->double{
      const int l {0}; // s-wave only
      double I {0.0};
      for (auto [r,w] : std::views::zip( rmesh.x, rmesh.w)) {
	I += w * r*r
	  * std::sph_bessel(l,ka*r)
	  * Vr(r) / constants::hbc // Vr in MeV -> 1/fm
	  * std::sph_bessel(l,kb*r) // MeV fm
	  * M * constants::hbc / (constants::hbc*constants::hbc); // Because of the H.-T. potential partial waves definition
      }
      return I;
    };
    // Matrix Vktmp = discretize(Vk,kp.data())
    // V = Eigen::Map<Eigen::MatrixXd>(Vtmp);
    // or:
    std::cout << "Doing ourier transform of V ...\n";
    for (int i : std::views::iota(0,N+1)) {
      for (int j : std::views::iota(0,N+1)) {
	V(i,j) = Vk(kp(i),kp(j));
      }
    }
  }

  // Construct the F matrix
  Eigen::MatrixXd F(N+1,N+1);
  for (int i : std::views::iota(0,N+1)) {
    for (int j : std::views::iota(0,N+1)) {
      F(i,j) = delta(i,j) + omegap(j) * V(i,j);
    }
  }
  // Invert the F matrix
  Eigen::MatrixXd Finv = F.inverse();
  // Test the inversion, F * F^{-1} = 1
  // std::cout << "Tr abs F * Finv = " << (F * Finv).array().abs().matrix().trace() << " (should be " << N+1 << ")\n";
  std::cout << "Tr (F * Finv) = " << (F * Finv).trace() << " (should be " << N+1 << ")\n";
  // R-matrix R(k_i, k0)
  Eigen::ArrayXd R = Finv * V.col(N);

  // On-shell R-matrix element R(k0,k0) is the last element of R
  return R(R.size()-1);
}

int main () {

  // Find bound states
  findSingleChannelBoundState();

  // Compute phase shifts
  {
    constexpr double E {0.1}; // C-o-m energy in MeV
    using constants::Mp;
    using constants::Mn;
    using constants::hbc;
    constexpr double m {Mn*Mp/(Mn+Mp)}; // Reduced mass in MeV
    constexpr double k0 {std::sqrt(2*m*E)/hbc}; // in 1/fm
      double R_k0 = computeRMatrix(k0);
      double tan_delta_k0 = - k0 * R_k0;
      std::cout << "tan delta(k)=" << tan_delta_k0 << " for E=" << E << " MeV (k=" << k0 << " 1/fm)\n";
  }
  return 0;
}
