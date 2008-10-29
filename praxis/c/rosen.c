#include <math.h>
#include "machine.h"

double rosen(x, n)
double x[];
int n;
{
   double f1, f2;

   f1 = x[0] - 1.0;
   f2 = 10.0*(x[1] - x[0]*x[0]);
   return f1*f1 + f2*f2;
}
