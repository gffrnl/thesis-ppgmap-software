/* fltcmp.h */
#pragma once

#include <stdbool.h>
#include <math.h>

bool
fltcmp_dig (float dbl1, float dbl2, unsigned dig)
{
  if (isless( fabsf((dbl1-dbl2)/fminf(dbl1, dbl2)),
              5 / powf(10.0, dig))
      == 0)
    return false;
  return true;
}

bool
dblcmp_dig (double dbl1, double dbl2, unsigned dig)
{
  if (isless( fabs((dbl1-dbl2)/fmin(dbl1, dbl2)),
              5 / pow(10.0, dig))
      == 0)
    return false;
  return true;
}

