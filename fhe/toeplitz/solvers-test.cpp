#include <iostream>

#include <vector>
#include <algorithm>  // transform
#include <numeric>    // inner_product, transform_reduce
#include <cmath>      // fabs

namespace toeplitz {
  namespace solvers {

    template<typename Real>
    std::vector<Real> durbin(std::vector<Real> const & r) {
      // Matrix computations, Golub & Van Loan (2012)
      // Algorithm 4.7.1 (Durbin)

      Real alpha, beta;
      std::vector<Real> y(r.size());

      y[0]  = -r[0];
      alpha =  y[0];
      beta  =   1.0;

      for (std::size_t k = 0; k < r.size()-1; ++k) {
        beta *= - (alpha*alpha - 1);
        alpha = - (r[k+1] +
                   std::transform_reduce(//execution::par_unseq,
                     std::next(r.crbegin(), r.size()-1-k), r.crend(),
                     y.cbegin(),
                     static_cast<Real>(0))
          ) / beta;

        std::vector<Real> z(r.size());
        std::transform(//execution::par_unseq,
          y.cbegin(), std::next(y.cbegin(), k+1),
          std::next(y.crbegin(), r.size()-1-k),
          z.begin(),
          [&alpha](Real const & a, Real const & b) -> double {
            return a + alpha * b;
          });
        std::move(z.begin(), std::prev(z.end(), 1), y.begin());
        y[k+1] = alpha;
      }
      return y;
    }

    
    template<typename Real>
    std::vector<Real> levinson(std::vector<Real> const & r,
                               std::vector<Real> const & b) {
      // Matrix computations, Golub & Van Loan (2012)
      // Algorithm 4.7.2 (Durbin)
      // Modified version

      // Check arguments
      if (r.size() != b.size())
        throw std::invalid_argument("in r and b have different sizes");
      if (std::fabs(r[0]) < std::numeric_limits<Real>::epsilon())
        throw std::invalid_argument("null diagonal");

      std::vector<Real> r2(r.size());
      std::vector<Real> b2(r.size());

      { // Rescale Toeplitz and rhs vectors:
        Real const r_0 = r[0];
        auto f = [&r_0](Real const& a) -> Real { return a/r_0; };
        std::transform(//execution::par_unseq,
          r.cbegin(), r.cend(), r2.begin(), f);
        r2[0] = 1.0;
        std::transform(//execution::par_unseq,
          b.cbegin(), b.cend(), b2.begin(), f);
      }
	
      Real alpha, beta, mu;
      std::vector<Real> y(r2.size());
      std::vector<Real> x(r2.size());
      y[0]  = -r2[1];
      x[0]  =  b2[0];
      alpha =  y[0];
      beta  =   1.0;
      for (std::size_t k = 0; k < r2.size()-1; ++k) {
        beta *= - (alpha*alpha - 1);
        mu    =   (b2[k+1] -
                   transform_reduce(//execution::par_unseq,
                     std::next(x.crbegin(), r2.size()-1-k), x.crend(),
                     std::next(r2.cbegin(), 1),
                     static_cast<Real>(0))
          ) / beta;

        std::vector<Real> v(r2.size());
        std::transform(//execution::par_unseq,
          x.cbegin(), std::next(x.cbegin(), k+1),
          std::next(y.crbegin(), r2.size()-1-k),
          v.begin(),
          [&mu](Real const & a, Real const & b) -> Real {
            return a + mu * b;
          });
        std::move(v.begin(), std::prev(v.end(), 1), x.begin());
        x[k+1] = mu;

        if (k < r2.size()-2) {
          alpha = - (r2[k+2] +
                     std::transform_reduce(//execution::par_unseq,
                       std::next(r2.crbegin(), r2.size()-(k+2)), r2.crend(),
                       y.cbegin(),
                       static_cast<Real>(0))
            ) / beta;

          std::vector<Real> z(r2.size());
          std::transform(//execution::par_unseq,
            y.cbegin(), std::next(y.cbegin(), k+1),
            std::next(y.crbegin(), r2.size()-1-k),
            z.begin(),
            [&alpha](Real const &a, Real const &b) -> Real {
              return a + alpha * b;
            });
          std::move(z.begin(), prev(z.end(), 1), y.begin());
          y[k+1] = alpha;
        }
      }

      return x;
    }
    
  }
}

#include <cmath>
#include <functional>
namespace toeplitz { namespace tests {

    template<typename Real = double>
    bool durbin_solver() {
      using namespace toeplitz::solvers;
      Real const tol =
        1000 * std::numeric_limits<Real>::epsilon();      
      Real sumabserr;
      std::vector<Real> y;
      {
        std::vector<Real> const y_exact = {
          - static_cast<Real>(75) / 140,
          + static_cast<Real>(12) / 140,
          - static_cast<Real>( 5) / 140
        };
        std::vector<Real> r = {
          static_cast<Real>(1) / 2 ,
          static_cast<Real>(1) / 5 ,
          static_cast<Real>(1) / 10
        };
        y = durbin(r);
        sumabserr = 
          std::inner_product(y_exact.cbegin(), y_exact.cend(),
                             y.cbegin(),
                             static_cast<Real>(0),
                             std::plus<Real>(),
                             [](Real const & a,
                                Real const & b) -> Real 
                             {
                               Real const c = a - b;
                               return (c < 0)? -c : c;
                             });
                              
      }
      if (sumabserr > tol) return false;
      return true;
    }


    template<typename Real = double>
    bool levinson_solver() {
      using namespace toeplitz::solvers;
      Real const tol =
        1000 * std::numeric_limits<Real>::epsilon();      
      Real sumabserr;
      std::vector<Real> y;
      {
        std::vector<Real> const y_exact = {
          + static_cast<Real>(355) / 56,
          - static_cast<Real>(376) / 56,
          + static_cast<Real>(285) / 56
        };        
        std::vector<Real> r = {
          static_cast<Real>(1)    ,
          static_cast<Real>(1) / 2,
          static_cast<Real>(1) / 5
        };
        std::vector<Real> b = { 4, -1, 3 };
        y = levinson(r, b); 
        sumabserr = 
          std::inner_product(y_exact.cbegin(), y_exact.cend(),
                             y.cbegin(),
                             static_cast<Real>(0),
                             std::plus<Real>(),
                             [](Real const & a,
                                Real const & b) -> Real 
                             {
                               Real const c = a - b;
                               return (c < 0)? -c : c;
                             });
                              
      }
      if (sumabserr > tol) return false;
      return true;
    }
    
  }}






int main(int argc, char * argv[])
{
  using namespace std;
  using namespace toeplitz;
  
  bool const durbin_passed = tests::durbin_solver<float>();
  bool const levins_passed = tests::levinson_solver<float>();

  if (durbin_passed)
    std::cout << "Durbin's solver test PASSED!"     << std::endl;
  else
    std::cout << "Durbin's solver test NOT passed!" << std::endl;

  if (levins_passed)
    std::cout << "Levinson's solver test PASSED!"     << std::endl;
  else
    std::cout << "Levinson's solver test NOT passed!" << std::endl;

  
  return 0;
}

