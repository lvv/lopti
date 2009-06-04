
	// NOTE: there is 2nd implementation:  http://www.netlib.org/opt/hooke.c
	#ifndef LOPTI_LINE_SEARCH
	#define LOPTI_LINE_SEARCH

	//#include <lopti/lopti.h>


	#include <cassert>
		//using namespace std;

	//#include <lvv/array.h>
		using lvv::array;

	//  ALGORITHM FROM [[[6]]] Stefen Boyd, "Backtraking Line Search, Armijo-Goldstein condition" (ee364a lecture 15,  18min)
	//  	α  ∈  (0, 0.5)
	//  	β  ∈  (0, 1.0)  ≈0.5
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
	// 	a=0.5; x0=2.9; set grid; e=2.71828;  f(x)=e**(1-x)+x; df(x)=1-e**(1-x);  f_mod_low(x)=f(x0)+df(x0)*(x-x0);   f_mod_high(x)=f(x0)+a*df(x0)*(x-x0);  plot [-1:4] [-5:5]  f(x) w p , df(x),  f_mod_low(x), f_mod_high(x)
namespace lopti  {

struct 	line_search_backtracking {
		int eval_cnt;
		objective
	line_search_backtracking (){};
	find


 }

}


	} // namespace lopti
	#endif  // LOPTI_H
