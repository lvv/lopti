 /************************************************************************
 *      function HJMIN - find local minimum of a given function
 *	   	        by the Hook-Jeevse method
 *
 * Input
 *	double hjmin(X,h0,funct)
 *	VECTOR X 			On input contains the initial
 *					guess to the minimum location
 *					On return has the vector
 *					specifying for the minimum location
 *	const VECTOR h0			Initial values for the steps along
 *					each directions
 *	double f(const VECTOR x)	Procedure to compute a function
 *					value at the specified point
 *
 * Output
 *	Hjmin returns the function value at the point of minimum
 *
 * Algorithm
 *	Hook-Jeevse method of direct search for a function minimum
 *	The method is of the 0. order (i.e. requiring no gradient computation)
 *	See
 *	B.Bondi. Methods of optimization. An Introduction - M.,
 *	"Radio i sviaz", 1988 - 127 p. (in Russian)
 *	
 *
 * Notes
 *	static VECTORs S and BASE are used as work arrays. They are not
 *	destroyed (freed) on exit so that next call to Hjmin could use
 *	them if they are still appropriate for that call.
 *
 ************************************************************************/
 


	#include <cassert>
	#include <algorithm>
	#include <iterator>
		using namespace std;

	#include <lvv/array.h>
		using namespace lvv;
		using lvv::array;
	//#include <lvv/math.S>
	//using lvv::powi;

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

	virtual minimizer<V>&		step(V& _S)	 { S=_S; return *this; };

	V&		 	argmin			()		{
					const fp_t	threshold = 1e-12;	   // Threshold for the function   
					const fp_t	step_reduce_factor = 10;	// decay to be treated as significant                  
					fp_t  f_min;				   // Min function value found     

		BASE = X ;
		f_base = f_min = (*loft_v) (X);


		for (;;) {										   // Main iteration loop        // X is a next approximation to min             
			if (examination() < f_min - threshold) {					   // Function falled markedly     
				do {									   // Perform search on patter     // BASE - pattern, X - new min appr 
					typename V::iterator  base_it =  BASE.begin();
					typename V::iterator  min_it =  X.begin();

					for (int i = BASE.ibegin();  i < BASE.iend();  i++) {
						 fp_t t = (*base_it - *min_it) + *base_it;

						*min_it++ = *base_it;
						*base_it++ = t;
					}
					f_min = f_base;
					f_base = (*loft_v) (BASE);
				}	while (examination() < f_min - threshold);	// Continue search until f doesn't  decrease         

				BASE = X;
				f_base = f_min;						// Save the best approx found   

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
					ymin_  =  f_min;
					return  X;
				}
			}

		}
 }

  private:

	V	S;				//  Steps along axes        
	V	BASE;				//  Base point              
	fp_t	f_base;				//  Function value at it    

	fp_t examination() {
		
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
