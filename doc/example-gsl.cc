	#include <lvvlib/math.h>
		using lvv::pow2;
	#include <lopti/gsl-nelder-mead-wrap.h>
		using lopti::MINIMIZER;
		using lopti::make_loft;
	typedef   array<double,2>   A;		// see boost::array
double   rosenbrock_fn (A& X)    {  return  100 * pow2(X[1]-pow2(X[0])) + pow2(1-X[0]);  };

int main() {
	MINIMIZER<A>	  m;
	m.loft	( make_loft<A> (&rosenbrock_fn,"rosenbrock") );  // converts plain function to loft and pass it as minimizer param
	A	X = {{ -1.2,  1  }};		// X0 - starting point
	A	S = {{  0.2,  0.2 }};		// step0 - initial step
	m.x0	(X);
	m.step0	(S);

	A	X_min	=  m.argmin();
	double	y_min	=  m.ymin();
	cout  <<  "Minimum: "  <<  y_min  <<  "    at point: "  <<  X_min  <<  endl;
	m.print(); // almost the same as 3 lines above
 }
