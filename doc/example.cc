
#include <lopti/hook-jeevs.h>
using lopti::hook_jeevs_minimizer;			// minimizer class
using lopti::make_objective;			

// Array-type 
typedef   lvv::array<double,2>   V;			// == double	X[2]

// Objective Function
double   rosenberg_fn (V& X)    {  return  100 * (X[1]-X[0]*X[0])*(X[1]-X[0]*X[0]) + (1-X[0])*(1-X[0]);  };

int main() {

	V	X = {{ -1.2,  1   }};                             // X0 - starting point. 
	V	S = {{  0.2,  0.2 }};                             // step0 - initial step

	hook_jeevs_minimizer<V>	  m;

	// minimizer config 
	m.objective     ( make_objective<V> (&rosenberg_fn,"rosenberg") );
	m.x0       	(X);
	m.step0    	(S);
	m.max_iter 	(500);

	// Result
	V	X_min	=  m.argmin();
	double	y_min	=  m.ymin();
	cout  <<  "Minimum: "  <<  y_min  <<  "    at point: "  <<  X_min  <<  endl;
 }
