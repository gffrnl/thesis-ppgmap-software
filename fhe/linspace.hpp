#pragma once
#include <cstddef>
#include <vector>

template<typename T>
std::vector<T> linspace(T start, T end, std::size_t num)
{
  std::vector<double> linspaced;

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  --num;
  for (std::size_t i = 0; i < num; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}
