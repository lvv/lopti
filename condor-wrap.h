	
	#ifndef LOPTI_CONDOR_H
	#define LOPTI_CONDOR_H

	// there is a namespace CONDOR, to avoid name conflict:
	#ifdef		CONDOR
		#undef	CONDOR
	#endif

	#include <lopti/lopti.h>
	#include <lopti/object_function.h>

	#ifndef		MINIMIZER
		#define	MINIMIZER	condor_minimizer
	#endif

	#ifndef		V_IBEGIN
		#define	V_IBEGIN 0
	#endif
	
	#include <lopti/convert-condor.h>
		using lvv::array;

	#include <stdio.h>
	#include <stdlib.h>
	#include <cmath>
	#include <string>

	#include <condor/Solver.h>
	#include <condor/tools.h>
	#include <condor/ObjectiveFunction.h>
		using CONDOR::sqr; using CONDOR::Vector; using CONDOR::CONDORSolver;
	#include <lvv/lvv.h>

	extern int globalPrintLevel;

namespace lopti {
			template<typename V>
class	c_of_t  : public CONDOR::ObjectiveFunction { public:
						OBJECTIVE_TYPES;
				int eval_cnt;
				bool verbose;
				objective_p_t objective_v;

	void	init (objective_p_t objective_p, V& X0) {  
			objective_v = objective_p;
			strcpy(name,"condor_of_wrap");
			xOptimal.setSize(X0.size());
			xStart.  setSize(X0.size());
			xStart << X0;
			//this->var =  var;
			eval_cnt = 0;
			verbose = false;
	};

	double  eval (CONDOR::Vector cX, int *nerror=NULL) {  		// condor use this to eval
			V X;	
			X << cX;
			T y = (*objective_v)(X);
			updateCounter(y,cX);
			eval_cnt++;

			if (verbose) {
				printf("%5d 	 %19.15g  ",  eval_cnt, y);
				for  (typename V::iterator x=X.begin();   x != X.end();   x++)    printf("\t%19.15g",  *x);
				cout << endl;
			}

			return	y;
	};
 };

/////////////////////////////////////////////////////////////// 
                 template<typename V>
struct	condor_minimizer: trust_region_minimizer<V>   {
						MINIMIZER_MEMBERS;  TR_MINIMIZER_MEMBERS;  OBJECTIVE_TYPES;
			c_of_t<V>			c_of;	// condor object func wrap	
									//CONDOR::ObjectiveFunction*	c_of;		
									//CONDOR::ObjectiveFunction*	c_rof_wrap;	// condor object func rescaled wrap	
			CONDOR::Vector			cX;		// param in condor format
			CONDOR::Vector			cR;		// rescale divider
									//	void		rescale			(V& R) 		{ cR << R; c_rof_wrap = new CONDOR::CorrectScaleOF(2, c_of, cR); };
	//explicit 			condor_minimizer	()			: trust_region_minimizer<V>("condor") {};
	minimizer<V>&	 	verbose			(bool flag)	{ verbose_     = flag; c_of.verbose = flag; return *this; };
	const string		name			() 	const	{  return minimizer<V>::mk_name("condor"); };
	V&		 	argmin			()		{
							assert(objective_v != 0 );  assert(!isnan(rho_begin_) && !isnan(rho_end_));
		c_of.init(objective_v, X);
		globalPrintLevel = 10;		// off
		CONDOR::CONDORSolver(rho_begin_, rho_end_, max_iter_,  &c_of);
		Xmin_ << (c_of.xBest);
		ymin_ = c_of.valueBest; 				//c_of.printStats();
		iter_ = c_of.getNFE(); 
		return Xmin_;
	};
};
	} // namespace lopti
	#endif // LOPTI_CONDOR_H
