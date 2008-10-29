#include <math.h>
#include "machine.h"

extern double rosen();
extern double praxis();
extern double tol;
extern int prin, maxfun;

main()
{
   double x[N];
   double fx;
   int n;
   
   x[0] = -1.2;
   x[1] = 1.0;  
   n = 2;
   tol = 1.0e-20;
   printf("print parameter: ");
   scanf("%d", &prin);
   fx = praxis(rosen, x, n);
   printf("Minimum = %e at point %e,%e\n", fx, x[0], x[1]);
}
