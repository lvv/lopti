

#include  <iostream>
	using namespace std;

//#include <lvv/lvv.h>
//	using namespace lvv;

	#include  <functional>
		using std::unary_function;

#include  "lopti.h"
	using namespace lopti;

#include  "line_search.h"


int main() { 

		typedef 	array<float,2,1>		V;	
	//	objective1<V>*	objective_v  = objective1<V>*(trace<V>(rosenbrock<V>()));
	//	line_search_backtracking_t<V>	lls(trace<V>(rosenbrock<V>()));
		line_search_backtracking_t<V>	ls = line_search_backtracking_t<V>(rosenbrock<V>());
	//	line_search_backtracking_t<V>	ls(*new rosenbrock<V>);
	//	line_search_backtracking_t<V>	lls(*objective_v);
		const V		X0 = {{ 1, 0 }};
		const V		DX = {{ 0, 1 }};

	cout << "starting point:    " << X0 << endl;
	cout << "direction:         " << DX << endl;
	cout << "line search found: " << ls.find (X0, DX) << endl;
	
	cout << endl;
}
