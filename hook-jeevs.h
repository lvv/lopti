
	#ifndef LOPTI_HOOK_JEEVS_H
	#define LOPTI_HOOK_JEEVS_H

	#include <lopti/lopti.h>

	#ifndef		MINIMIZER
		#define	MINIMIZER	hook_jeevs_minimizer
	#endif

	#include <cassert>
	#include <algorithm>
	#include <iterator>
		using namespace std;

	#include <lvv/array.h>
		using namespace lvv;
		using lvv::array;

namespace lopti  {
					// Examine the function f in the vicinity of the BASE point
					// by making tentative steps fro/back along the every coordinate.
					// When the function is found to decrease, the point BASE is correspondingly
					// updated. Examination() returns the minimal function value found in
					// the investigated region.

			template<typename V> 
struct   hook_jeevs_minimizer   : minimizer<V>  {
						MINIMIZER_MEMBERS;  LOFT_TYPES;
					#define  EPSILON 0.000001
						fp_t	tau_;	   // Termination criterium        

	hook_jeevs_minimizer	()			: minimizer<V>("hook-jeevs"), tau_(10*EPSILON) {};
	virtual	minimizer<V>&	 	tau		(fp_t _tau)	{ tau_ = _tau; return *this; };

	virtual minimizer<V>&		step0(V& _S)	 { S=_S; return *this; };

	V&		 	argmin			()		{
					const fp_t	threshold = 1e-12;	   // Threshold for the function   
					const fp_t	step_reduce_factor = 10;	// decay to be treated as significant                  
					//fp_t  ymin_;				   // Min function value found     

		BASE = X ;
		f_base = ymin_ = (*loft_v) (X);


		for (;;) {										   // Main iteration loop        // X is a next approximation to min             
			if (examination() < ymin_ - threshold) {					   // Function falled markedly     
				do	{									   // Perform search on patter     // BASE - pattern, X - new min appr 
						V t;
						((t = BASE)  *= 2.) -= X;

						X = BASE;
						BASE = t;


					ymin_ = f_base;
					f_base = (*loft_v) (BASE);
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
					return  X;
				}
			}

		}
 }

  private:

			V	S;				//  Steps along axes        
			V	BASE;				//  Base point              
			fp_t	f_base;				//  Function value at it    

	fp_t   examination() {
		
		for (int  i = BASE.ibegin();  i < BASE.iend(); i++) {	// Perform step along a coordinate              
			 fp_t basi_old = BASE[i];		   // Old coordinate value         
			 fp_t f_new;

			BASE[i] = basi_old + S[i];
			f_new = (*loft_v) (BASE); 

			if	(f_new < f_base)  {
				f_base = f_new;								// Step caused f to decrease, OK 
			} else {
				BASE[i] = basi_old - S[i];
				f_new = (*loft_v) (BASE); 
				if (f_new < f_base) f_base = f_new;
				else               BASE[i] = basi_old;					 // No fall was found along this coord 
			}
		}
		return f_base;
	 }

 }; // class 

	} // namespace lopti
	#endif  // LOPTI_H
