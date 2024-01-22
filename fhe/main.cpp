// Solution of FHE

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "frheateq.hpp"
#include "closest.hpp"

int main(int argc, char * argv[]) {
  using std::cout;
  using std::endl;
  double const ealpha = 1.4;
  double const a = -1000.1;
  double const b =  1000.1;
  double const a0 = -1.1;
  double const b0 =  1.1;
  double const deltat = 0.1;
  double const deltax = 0.1;
  double const max_solution_time = 1.0;
  std::string const filename{"bump.dat"};
  
//  frheateq::solvers::explicit_euler solver;
//  frheateq::solvers::implicit_euler solver;
  frheateq::solvers::crank_nicolson solver;

  
  auto u0 = [ealpha](double const & x) -> double {
    double const sgn = ( std::isgreater(ealpha, 1.0) ? -1.0 : 1.0 );
    return sgn * std::pow(1.0 + x*x, -0.5+ealpha/2.0);
  };

  auto phi = [](double const & x) -> double {
    double const norm = std::fabs(x);
    if (!(norm < 1))
      return 0;
    return std::exp(-1.0/(1-norm*norm));    
  };

  
  solver.set_ealpha(ealpha);
  solver.set_domain(a,b);
  solver.set_deltax(deltax);
  solver.set_deltat(deltat);
  solver.set_diffus(1);
  
  std::size_t J = solver.prepare();

  cout << "numer of nodes = " << J << endl;
  cout << "deltax = " << solver.get_deltax() << endl;

  std::vector<double> const X{linspace(a, b, J)};
  std::vector<double> U0(X.size());
  std::transform(X.cbegin(), X.cend(), U0.begin(), phi);


  
  std::ofstream ofs(filename);
  if (!ofs.is_open()) {
    std::cerr << "cannot open file " << filename << std::endl;
    std::abort();
  }
  ofs << std::scientific;
  ofs << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  std::size_t ja = closest<double>(a0, X);
  std::size_t jb = closest<double>(b0, X);

  /*
  std::cout << "ja = " << ja << std::endl;
  std::cout << "jb = " << jb << std::endl;
  std::cout << "X[ja] = " << X[ja] << std::endl;
  std::cout << "X[jb] = " << X[jb] << std::endl;
  return 0;
  */
  
  for (std::size_t j = ja; j <= jb; ++j) {
    ofs << std::setw(25) << 0.0       << ' ';
    ofs << std::setw(25) << X.at(j)   << ' ';
    ofs << std::setw(25) << U0.at(j) << '\n';
  }
  ofs << std::endl;

  
  // PRECONDITIONS
  std::vector<double> Un{U0};
  std::vector<double> Un1;
  bool   program_continue = false;
  double solution_time    = solver.get_deltat();
  while (!(solution_time > max_solution_time)) {

    std::cout << "solution time = " << solution_time << std::endl;
        
    Un1 = solver.advance(Un); // TODO: pass Un1 by reference
    
    for (std::size_t j = ja; j <= jb; ++j) {
      ofs << std::setw(25) << solution_time << ' ';
      ofs << std::setw(25) << X.at(j)       << ' ';
      ofs << std::setw(25) << Un1.at(j)     << '\n';
    }
    ofs << std::endl;
    
    /*
    // TODO: program continuation when I/O operation fails
    if (!ofs.good() && !program_continue) {
      std::cout << "Some error occured while writing to file `" << filename << "`.";
      // TODO: implement program continuation....
      if (!program_continue) std::abort();
    }
    */
    
    Un  = Un1; // TODO: swap
    solution_time += solver.get_deltat();
  }
  
  
  ofs.close();
  return EXIT_SUCCESS;
}

