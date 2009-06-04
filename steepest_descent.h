
	// NOTE: there is 2nd implementation:  http://www.netlib.org/opt/hooke.c
	#ifndef LOPTI_STEEPEST_DESCENT
	#define LOPTI_STEEPEST_DESCENT

	#include <lopti/lopti.h>

	#ifndef		MINIMIZER
		#define	MINIMIZER	steepest_descent_minimizer
	#endif

	#include <cassert>
		//using namespace std;

	#include <lvv/array.h>
		//using namespace lvv;
		using lvv::array;

namespace lopti  {

			template<typename V> 
struct   MINIMIZER   : minimizer<V>  {
						MINIMIZER_MEMBERS;  OBJECTIVE_TYPES;

	explicit			MINIMIZER	()		: minimizer<V>("steepest_descent") {};
	//virtual	minimizer<V>&	 	tau			(T _tau)	{ tau_ = _tau; return *this; };

	V&		 	argmin			()		{

										/*
										def steepest(x0):

										    i = 0
										    iMax = 10
										    x = x0
										    Delta = 1
										    alpha = 1

										    while i<iMax and Delta>10**(-5):
											p = -Jacobian(x)
											xOld = x
											x = x + alpha*p
											Delta = sum((x-xOld)**2)
											print x
											i += 1
										*/
			    int	i = 0;
			    int	iMax = 10;
			    X = X0;
			    Delta = 1
			    alpha = 1

			    while i<iMax and Delta>10**(-5):
				p = -Jacobian(x)
				xOld = x
				x = x + alpha*p
				Delta = sum((x-xOld)**2)
				print x
				i += 1



		}

  private:


 }; // class 

	} // namespace lopti
	#endif  // LOPTI_H
