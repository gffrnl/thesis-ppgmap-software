#include <iostream>
#include <vector>

struct toeplitz_matrix
{
  std::vector<double> rv;
  std::vector<double> cv;
  
  toeplitz_matrix(std::initializer_list<double> r) : rv(r)
  {
    
  }
  
  bool is_symmetric() const
  {
    if (cv.size() == 0)
      return true;
    return false;
  }

  std::vector<double> operator=(toeplitz_matrix& m)
  {
    if (!m.is_symmetric())
      throw "toeplitz matrix not symmetric";
    return rv;
  }
  std::vector<double> operator=(toeplitz_matrix&& m)
  {
    if (!m.is_symmetric())
      throw "toeplitz matrix not symmetric";
    return rv;
  }
  
};


  

int main() 
{
  toeplitz_matrix A = {
    0.7560439,
    0.0002211,
    0.3303271,
    0.6653811,
    0.6283918
  };
  std::cout << "A.is_symmetric() = "
            << std::boolalpha << A.is_symmetric()
            << std::endl;

  std::vector<double> v1 = A;
  std::vector<double> v2;

  v2 = A;

  
  return 0;
}
