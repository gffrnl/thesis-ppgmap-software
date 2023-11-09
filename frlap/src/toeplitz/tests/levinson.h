#pragma once
#include <stddef.h>

int
levinson (const size_t n,
    const double A1[restrict const static n],
          double x [restrict const static n],
    const double b [restrict const static n]);
    
int
levinson_in_place (const size_t n,
          double A1[restrict const static n],
          double x [restrict const static n],
          double b [restrict const static n]);
