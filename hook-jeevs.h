
	// NOTE: there is 2nd implementation:  http://www.netlib.org/opt/hooke.c
	#ifndef LOPTI_HOOK_JEEVS_H
	#define LOPTI_HOOK_JEEVS_H

	#include <lopti/lopti.h>

	#ifndef		MINIMIZER
		#define	MINIMIZER	hook_jeevs_minimizer
	#endif

	#include <cassert>
	#include <algorithm>
	#include <iterator>
		//using namespace std;

	#include <lvv/array.h>
		//using namespace lvv;
		using lvv::array;

namespace lopti  {
					// Examine the function f in the vicinity of the BASE point
					// by making tentative steps fro/back along the every coordinate.
					// When the function is found to decrease, the point BASE is correspondingly
					// updated. Examination() returns the minimal function value found in
					// the investigated region.

			template<typename V> 
struct   hook_jeevs_minimizer   : minimizer<V>  {
						MINIMIZER_MEMBERS;  OBJECTIVE_TYPES;
					#define  EPSILON 0.000001
						T	tau_;	   // Termination criterium        

	explicit			hook_jeevs_minimizer	()		:  tau_(10*EPSILON) { S.assign(0.2);};
	const string		name			() 	const	{  return minimizer<V>::mk_name("hook-jeevs"); };
	minimizer<V>&	 	tau			(T _tau)	{ tau_ = _tau; return *this; };
	minimizer<V>&		step0			(V& _S)		{ S=_S; return *this; };

	V&		 	argmin			()		{
					const T	threshold = 1e-12;	   // Threshold for the function   
					const T	step_reduce_factor = 10;	// decay to be treated as significant                  
					//T  ymin_;				   // Min function value found     

		BASE = X ;
		f_base = ymin_ = (*objective_v) (X);	
		iter_++;


		for (;;) {										   // Main iteration loop        // X is a next approximation to min             
			if (examination() < ymin_ - threshold) {					   // Function falled markedly     
				do	{									   // Perform search on patter     // BASE - pattern, X - new min appr 
						V t;
						((t = BASE)  *= 2.) -= X;

						X = BASE;
						BASE = t;


					ymin_ = f_base;
					f_base = (*objective_v) (BASE);	
					iter_++;

				}	while (examination() < ymin_ - threshold);	// Continue search until f doesn't  decrease         

				BASE = X;
				f_base = ymin_;						// Save the best approx found   

			} else {							// Function didn't fall markedly 
				// upon examination near BASE    
				int success = 1;					// Try to reduce the step       
				typename V::iterator	s_it   =  S.begin();
				typename V::iterator	base_it =  BASE.begin();

				for (int i = S.ibegin();  i < S.iend();  i++) {
					*s_it /= step_reduce_factor;
					success &= (*s_it++ / (1 + fabs(*base_it++)) < tau_);
				}

				if (success) {
					found_ = true;
					ymin_ = f_base;
					Xmin_ = X;
					return  Xmin_;
				}

				if (iter_  >  max_iter_) {
					return  Xmin_;
				}
			}

		}
 }

  private:

			V	S;				//  Steps along axes        
			V	BASE;				//  Base point              
			T	f_base;				//  Function value at it    

	T   examination() {
		
		for (int  i = BASE.ibegin();  i < BASE.iend(); i++) {	// Perform step along a coordinate              
			 T basi_old = BASE[i];		   // Old coordinate value         
			 T f_new;

			BASE[i] = basi_old + S[i];
			f_new = (*objective_v) (BASE); 
			iter_++;

			if	(f_new < f_base)  {
				f_base = f_new;								// Step caused f to decrease, OK 
			} else {
				BASE[i] = basi_old - S[i];
				f_new = (*objective_v) (BASE); 
				iter_++;
				if (f_new < f_base) f_base = f_new;
				else               BASE[i] = basi_old;					 // No fall was found along this coord 
			}
		}
		return f_base;
	 }

 }; // class 

	} // namespace lopti
	#endif  // LOPTI_H
