	#include <lvv/math.h>
	using lvv::pow2;
	#include	<lopti/hook-jeevs.h>
	using lopti::hook_jeevs_minimizer;
	#include	<lopti/object_function.h>
	using lopti::make_loft;


typedef 	array<double,2> 	A;

double   	rosenberg_fn (A& X)   {  return  100 * pow2(X[1]-pow2(X[0])) + pow2(1-X[0]); };

int main() {

	A		X = {{ -1.2,  1  }};	// X0 - starting point
	A		S = {{  0.2,  0.2 }};	// step0 - initial step
	
	hook_jeevs_minimizer<A>	  m;
		m.loft		( make_loft<A> (&rosenberg_fn,"rosenberg") );
		m.x0		(X);
		m.step0		(S);

	A	X_min	=  m.argmin();
	double	y_min	=  m.ymin();

	cout  <<  "minimu at: "  <<  X_min  <<  "   with value: "  <<  y_min  <<  endl;
 }
