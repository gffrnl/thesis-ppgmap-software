#include "fsp-v2.hpp"
#include <iostream>
#include <numeric>

std::vector<double> const yexact = {
  1.3306224490559279693969,
  0.9145886970054426523546,
  1.1300942935645541531642,
  0.8439319542836019039100, 
  1.7043853495260707919812
};

std::vector<double> const c = {
    0.7560439,
    0.0002211,
    0.3303271,
    0.6653811,
    0.6283918
};

std::vector<double> const x = {
    0.8497452,
    0.6857310,
    0.8782165,
    0.0683740,
    0.5608486
};



int main(int argc, char * argv[])
{
  double const tol = 10e-5;
  double sumabserr;

  std::cout << "tol = " << tol << std::endl;
  
  auto y = fast_symm_prod(c, x);
  
  std::cout << "y = \n";
  for (auto elem : y)
    std::cout << " " << elem << "\n";
  std::cout << std::endl;

  sumabserr = 
    std::inner_product(yexact.cbegin(), yexact.cend(),
                       y.cbegin(),
                       static_cast<double>(0),
                       std::plus<double>(),
                       [](double const & a,
                          double const & b) -> double
                       {
                         double const c = a - b;
                         return (c < 0)? -c : c;
                       });
  std::cout << "sumabserr = " << sumabserr << std::endl;
  
  if (sumabserr > tol)
    std::cout << "fsp-test NOT passed!" << std::endl;
  else
    std::cout << "fsp-test PASSED!"     << std::endl;

      
  return 0;
}
