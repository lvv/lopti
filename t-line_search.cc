

#include <iostream>
	using namespace std;

//#include <lvv/lvv.h>
//	using namespace lvv;

	#include	<functional>
		using std::unary_function;
	//#include	<boost/function.hpp>
	//	using boost::function;
		#include   	<lopti/hook-jeevs.h>

#include "lopti.h"
	using namespace lopti;

#include "line_search.h"


int main() { 
	typedef 	array<double,_N,0>		V;	
	objective_v_t	objective (trace<V>(rosenbrock<V>()));
	line_search_backtracking_t<V>	ls(trace<V>(rosenbrock<V>()));
	V		X0 = {{ 1, 0 }};
	V		DX = {{ 0, 1 }};
	cout << "starting point:    " << X0 << endl;
	cout << "direction:         " << DX << endl;
	cout << "line search found: " << ls.find(X0, DX) << endl;
	
	cout << endl;
}
