#pragma once
#include <stddef.h>

void
durbin (const size_t n,
    const double r[restrict const static n],
          double y[restrict const static n]);
