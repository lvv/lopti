
	#ifndef LOPTI_LINE_SEARCH
	#define LOPTI_LINE_SEARCH

	#include <lopti/lopti.h>

	#include <cassert>
		//using namespace std;

	#include <lvv/array.h>
		using lvv::array;

	//  ALGORITHM FROM [[[6]]] Stefen Boyd, "Backtraking Line Search, Armijo-Goldstein condition" (ee364a lecture 15,  18min)
	//  	α  ∈  (0, 0.5);        Zarowsky: (0.1, 0.3)
	//  	β  ∈  (0, 1.0)  ≈0.5;  Zarowsky: (0.1, 0.5)
	//  	t0 ≈ 1;
	//
	//  	t=t0;		
	//  	Δx = -∇f(x)	
	//  	do  {
	//  		t=t/β;  
	//  	} until ( x + t * Δx  <  f(x) + α * t * ∇f(x) * Δx) 	// lvv:  relace ∇f with normalized ∇f/|∇f| 
	//
	//  IDIA: t0 for next step:
	//  
	//  	if (steps in last iter > 1)
	//  		decrease t0
	//  	else
	//  		if  (too close to f_mod1_low()  in last iter)   increase t0
	//  		if  (too close to f_mod1_high() in last iter)   decrease t0

	// gnuplot:
	// 	set grid; set yrange[-1:4]; f(x)=(x-1)**2 +1; df(x)=2*(x-1); x0=0; a=0.3; ;dx=-df(x0)/abs(df(x0));   plot [-1:2]  f(x0+dx*x),  f(x0)+dx*df(x0)*x, f(x0) +a*dx*df(x0)*x
	// 	a=.5; x0=3; set grid; e=2.71828;  f(x)=e**(1-x)+x; df(x)=1-e**(1-x);  f_mod_low(x)=f(x0)+df(x0)*(x-x0);   f_mod_high(x)=f(x0)+a*df(x0)*(x-x0);  plot [-1.8:12] [-2:12]  f(x) w p , df(x),  f_mod_low(x), f_mod_high(x)
	// 
	// 	for quadratic F(),  optimal step is F'/2
	// 		plot [0:4] x**2, 2*x
	// 	
namespace lopti  {

		template<typename V>
struct 	line_search_backtracking_t   {	
						OBJECTIVE_TYPES;
				objective_p_t	objective_v;	
				const T		alpha;
				const T		beta;
				const T		t0;
	line_search_backtracking_t (
		//virtual void		objective	(objective_cref_t ref)	{ wrapped_objective_v = objective_p_t(&ref.clone());	name_ = ref.name(); };
		objective_cref_t ref,
		const T		alpha = 0.5,
		const T		beta  = 0.5,
		const T		t0    = 1. 
	) :
		objective_v	(&ref.clone()),
		alpha		(alpha),
		beta		(beta),
		t0		(t0)
	{};


	V&&	find( const V&	X0,	const V& DX) {
				T  t       = t0;
				T  f0      ( objective_v->eval0 (X0));
				V  G0      ( objective_v->eval1 (X0));
				V  X;

		for (int i = 1;  i< 50;  i++) {
			X = X0 + t * DX;
			T f = objective_v->eval0(X);
			if ( f  <  f0 + alpha * t * dot(G0,DX)) 	 {
				return  std::move(X);
			}
			t = beta * t;
		}

		assert(false);
		return  std::move(X);
	};
 };

	} // namespace lopti
	#endif  // LOPTI_H
