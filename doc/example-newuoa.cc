	#include <lvv/math.h>
		using lvv::pow2;
	#include <lopti/newuoa-wrap.h>
		using lopti::MINIMIZER;			// current minimizer #define-ed in above include 
		using lopti::make_objective;

	typedef   array<double,2,1>   V;		// see boost::array
double   rosenbrock_fn (V& X)    {  return  100 * pow2(X[2]-pow2(X[1])) + pow2(1-X[1]);  };

int main() {
	MINIMIZER<V>	  m;
	m.objective	( make_objective<V> (&rosenbrock_fn,"rosenbrock") );  // converts plain function to objective class and pass it as minimizer param
	V	X = {{ -1.2,  1  }};		// X0 - starting point
	m.x0	(X);
	m.rho_begin	(1);
	m.rho_end	(0.001);
	m.verbose	(true);
	m.max_iter	(200);

	V	X_min	=  m.argmin();
	double	y_min	=  m.ymin();
	cout  <<  "Minimum: "  <<  y_min  <<  "    at point: "  <<  X_min  <<  endl;
	m.print(); // almost the same as 3 lines above
 }
