
// select minimizer with include. 
#include <lopti/hook-jeevs.h>
using lopti::hook_jeevs_minimizer;
using lopti::make_loft;

// Array-type (see boost::array) 
typedef   array<double,2>   A;

// Objective Function
double   rosenberg_fn (A& X)    {  return  100 * (X[1]-X[0]*X[0])*(X[1]-X[0]*X[0]) + (1-X[0])*(1-X[0]);  };

int main() {

	hook_jeevs_minimizer<A>	  m;

	// configure minimizer
	m.loft	( make_loft<A> (&rosenberg_fn,"rosenberg") );  // converts plain function to loft and pass it to minimizer
	A	X = {{ -1.2,  1   }};		// X0 - starting point. 
	A	S = {{  0.2,  0.2 }};		// step0 - initial step
	m.x0	(X);
	m.step0	(S);
	m.max_iter(500);

	// Result
	A	X_min	=  m.argmin();
	double	y_min	=  m.ymin();
	cout  <<  "Minimum: "  <<  y_min  <<  "    at point: "  <<  X_min  <<  endl;
 }
