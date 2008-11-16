#include <lopti/hook-jeevs.h>
using lopti::hook_jeevs_minimizer;
using lopti::make_loft;

typedef   array<double,2>   A;		// see boost::array
double   rosenberg_fn (A& X)    {  return  100 * (X[1]-X[0]*X[0])*(X[1]-X[0]*X[0]) + (1-X[0])*(1-X[0]);  };

int main() {
	hook_jeevs_minimizer<A>	  m;
	m.loft	( make_loft<A> (&rosenberg_fn,"rosenberg") );  // converts plain function to loft and pass it as minimizer param
	A	X = {{ -1.2,  1   }};		// X0 - starting point
	A	S = {{  0.2,  0.2 }};		// step0 - initial step
	m.x0	(X);
	m.step0	(S);

	A	X_min	=  m.argmin();
	double	y_min	=  m.ymin();
	cout  <<  "Minimum: "  <<  y_min  <<  "    at point: "  <<  X_min  <<  endl;
	m.print(); // almost the same as 3 lines above
 }
