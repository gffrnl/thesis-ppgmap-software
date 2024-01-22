#pragma once

#include <frlap.hpp>
#include <frlap/gdm/strategies/centered_3_point_periodized.hpp>
#include "linspace.hpp"
#include "toeplitz/fsp.hpp"
#include "toeplitz/solvers/levinson.hpp"
#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>
//#include <functional>

namespace frlap { namespace gdm { namespace strategies {
  using c3point =
    typename frlap::gdm::strategies::centered_3_point_periodized;
    }}} // frlap::gdm::strategies

namespace frheateq {

  class solver {
  public:
    virtual ~solver() {}
    virtual std::size_t prepare() = 0;
    virtual std::vector<double> advance(std::vector<double>& Un) = 0;
  public:
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
    double get_deltax() 
    {
      return mdeltax;
    }
    double get_deltat()
    {
      return mdeltat;
      
    }
    
    
  protected:
    double mealpha;
    double mdeltax;
    double mdeltat;
    double mdiffus;
    std::pair<double, double> mdomain;
  };

  
  namespace solvers {

    class explicit_euler : public frheateq::solver {
    public:
      std::size_t prepare() override {
        using std::cout;
        using std::endl;
    
        //  cout << "\n\tIn function prepare..." << endl;
    
        std::size_t const J =
          std::ceil((mdomain.second - mdomain.first) / mdeltax) + 1;
        mdeltax = (mdomain.second - mdomain.first) / (J - 1);
        std::vector<double> r(J);
        frlap::gdm::strategies::c3point()
          .generate_coefficients(mealpha, mdeltax, r);
        //cout << "MU = {"; for (auto mu : r) cout << mu << ", "; cout << "\b\b}" << endl;
        { // Multiply the coefficients by - diffus * deltat:
          double mult = - mdiffus * mdeltat;
          std::for_each(r.begin(), r.end(), [&mult](double & c) { c *= mult; });
          //cout << "r = {"; for (auto c : r) cout << c << ", "; cout << "\b\b}" << endl;
        }
    
        // Generate matrix H-:
        Hm = r;
        Hm.at(0) += 1;
        //cout << "Hm = {"; for (auto c : Hm) cout << c << ", "; cout << "\b\b}" << endl;
    
        return J;
      }
  
  
  
      std::vector<double> advance(std::vector<double>& Un) override
      //               std::vector<double> const & Rn1)
      {
        return fast_symm_prod(Hm, Un);
      }
    private:
      std::vector<double> Hm;
    };

    class implicit_euler : public frheateq::solver {
    public:
      std::size_t prepare() override {
        using std::cout;
        using std::endl;
    
        //cout << "\n\tIn function prepare..." << endl;
    
        std::size_t const J =
          std::ceil((mdomain.second - mdomain.first) / mdeltax) + 1;
        mdeltax = (mdomain.second - mdomain.first) / (J - 1);
        std::vector<double> r(J);
        frlap::gdm::strategies::c3point()
          .generate_coefficients(mealpha, mdeltax, r);
        //cout << "MU = {"; for (auto mu : r) cout << mu << ", "; cout << "\b\b}" << endl;
        { // Multiply the coefficients by diffus * deltat / 2:
          double mult = mdiffus * mdeltat;
          std::for_each(r.begin(), r.end(), [&mult](double & c) { c *= mult; });
          //cout << "r = {"; for (auto c : r) cout << c << ", "; cout << "\b\b}" << endl;
        }
    
        // Generate matrix H+:
        Hp = r; Hp.at(0) += 1;
        //cout << "Hp = {"; for (auto c : Hp) cout << c << ", "; cout << "\b\b}" << endl;  

        return J;
      }

      std::vector<double> advance(std::vector<double>& Un) override
        //               std::vector<double> const & Rn1)
      {
        return toeplitz::solvers::levinson<double>(Hp, Un);
      }
    private:
      std::vector<double> Hp;
    };

    class crank_nicolson : public frheateq::solver {
    public:
        std::size_t prepare() override  {
          using std::cout;
          using std::endl;
    
          //cout << "\n\tIn function prepare..." << endl;
    
          std::size_t const J =
            std::ceil((mdomain.second - mdomain.first) / mdeltax) + 1;
          mdeltax = (mdomain.second - mdomain.first) / (J - 1);
          std::vector<double> r(J);
          frlap::gdm::strategies::c3point()
            .generate_coefficients(mealpha, mdeltax, r);
          //cout << "MU = {"; for (auto mu : r) cout << mu << ", "; cout << "\b\b}" << endl;

          { // Multiply the coefficients by diffus * deltat / 2:
            double mult = mdiffus * mdeltat / 2;
            std::for_each(r.begin(), r.end(), [&mult](double & c) { c *= mult; });
            //cout << "r = {"; for (auto c : r) cout << c << ", "; cout << "\b\b}" << endl;
          }
          // Generate matrix H+:
          Hp = r; Hp.at(0) += 1;
          //cout << "Hp = {"; for (auto c : Hp) cout << c << ", "; cout << "\b\b}" << endl;  

          { // Multiply the coefficients by diffus * deltat / 2:
            std::for_each(r.begin(), r.end(), [](double & c) { c *= -1; });
            //cout << "r = {"; for (auto c : r) cout << c << ", "; cout << "\b\b}" << endl;
          }
          // Generate matrix H-:
          Hm = r;
          Hm.at(0) += 1;
          //cout << "Hm = {"; for (auto c : Hm) cout << c << ", "; cout << "\b\b}" << endl;
          return J;
        }

      std::vector<double> advance(std::vector<double>& Un) override
        //               std::vector<double> const & Rn1)
      {
        auto Qn1 = fast_symm_prod(Hm, Un);
        return toeplitz::solvers::levinson<double>(Hp, Qn1);
      }
    private:
      std::vector<double> Hm, Hp;
    };
    
  }
  
}
