// Solution of FHE

#include <frlap.hpp>
#include <frlap/gdm/strategies/centered_3_point_periodized.hpp>
#include "linspace.hpp"
#include "toeplitz/fsp.hpp"
#include "toeplitz/solvers/levinson.hpp"
#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <functional>

namespace fhe {
  std::vector<double> advance(std::vector<double>& Un,
               std::vector<double> const & Hm,
               std::vector<double> const & Hp)
  //               std::vector<double> const & Rn1)
  {
    auto Qn1 = fast_symm_prod(Hm, Un);
    return toeplitz::solvers::levinson<double>(Hp, Qn1);
  }
}

namespace frlap { namespace gdm { namespace strategies {
  using c3point =
    typename frlap::gdm::strategies::centered_3_point_periodized;
    }}} // frlap::gdm::strategies

struct frheateq_implicit_euler {
  double mealpha;
  double mdeltax;
  double mdeltat;
  double mdiffus;
  std::pair<double, double> mdomain;
  std::function<double(double)> mic;
  std::vector<double> Hp;
  
  void set_ealpha(double const& ealpha) {
    if (ealpha != mealpha) {
      mealpha = ealpha;
      //perform operations
    }
  }
  void set_domain(double const& a, double const& b) {
    // check arguments
    mdomain.first  = a;
    mdomain.second = b;
    // perform_operations
  }
  void set_deltax(double const& deltax) {
    // check arguments
    mdeltax = deltax;
    // perform_operations
  }
  void set_deltat(double const& deltat) {
    // check arguments
    mdeltat = deltat;
    // perform_operations
  }
  void set_diffus(double const& diffus) {
    // check arguments
    mdiffus = diffus;
    // perform_operations
  }
  void set_u0(std::function<double(double)> const& u0) {
    // check arguments
    mic = u0;
    // perform operations
  }

  std::size_t prepare() {
    using std::cout;
    using std::endl;
    
    cout << "\n\tIn function prepare..." << endl;
    
    std::size_t const J =
      std::ceil((mdomain.second - mdomain.first) / mdeltax) + 1;
    mdeltax = (mdomain.second - mdomain.first) / (J - 1);
    std::vector<double> r(J);
    frlap::gdm::strategies::c3point()
      .generate_coefficients(mealpha, mdeltax, r);
    cout << "MU = {"; for (auto mu : r) cout << mu << ", "; cout << "\b\b}" << endl;
    { // Multiply the coefficients by diffus * deltat / 2:
      double mult = mdiffus * mdeltat;
      std::for_each(r.begin(), r.end(), [&mult](double & c) { c *= mult; });
      cout << "r = {"; for (auto c : r) cout << c << ", "; cout << "\b\b}" << endl;
    }
    
    // Generate matrix H+:
    Hp = r; Hp.at(0) += 1;
    cout << "Hp = {"; for (auto c : Hp) cout << c << ", "; cout << "\b\b}" << endl;  

    return J;
  }
  
  
  
  std::vector<double> advance(std::vector<double>& Un)
  //               std::vector<double> const & Rn1)
  {
    return toeplitz::solvers::levinson<double>(Hp, Un);
  }
  
};




int main(int argc, char * argv[]) {
  using std::cout;
  using std::endl;

  frheateq_implicit_euler solver;
  double const ealpha = 0.4;
  
  auto u0 = [ealpha](double const & x) -> double {
    double const sgn = ( std::isgreater(ealpha, 1.0) ? -1.0 : 1.0 );
    return sgn * std::pow(1.0 + x*x, -0.5+ealpha/2.0);
  };
  
  solver.set_ealpha(1.4);
  solver.set_domain(-4.0, 4.0);
  solver.set_deltax(0.8);
  solver.set_deltat(0.1);
  solver.set_diffus(1);
  solver.set_u0(u0);
  
  std::size_t J = solver.prepare();

  cout << "numer of nodes = " << J << endl;
  cout << "deltax = " << solver.mdeltax << endl;

  auto Un = linspace(-4.0, 4.0, J);

  cout << "X = {"; for (auto x : Un) cout << x << ", "; cout << "\b\b}" << endl;

  std::transform(Un.cbegin(), Un.cend(), Un.begin(), u0);

  cout << "U0 = {"; for (auto x : Un) cout << x << ", "; cout << "\b\b}" << endl;

  auto Un1 = solver.advance(Un);
  cout << "U1 = {"; for (auto x : Un1) cout << x << ", "; cout << "\b\b}" << endl;

  

  
  
  return 0;
}

