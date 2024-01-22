#pragma once
#include <cstddef>
#include <cmath>
#include <vector>

template<typename T>
std::size_t closest(T val, std::vector<T> const & x)
{
  std::size_t n = x.size();
  std::size_t idx = 0;

  if ( (val < x[0]) || (val > x[n-1]) )
    return n;

  n--;
  for (std::size_t i = 1; i < n; i++)
    if ( std::fabs(x[idx]-val) > std::fabs(x[i]-val) )
      {
        idx = i;
        if ( std::fabs(x[i+1]-val) > std::fabs(x[idx]-val) )
          return idx;
      }

  if ( std::fabs(x[n]-val) < std::fabs(x[idx]-val) )
    return n;

  return idx;
}
