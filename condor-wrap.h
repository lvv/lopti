	
	#ifndef LVV_LOPTI_CONDOR_H
	#define LVV_LOPTI_CONDOR_H

	#include <lopti.h>

	#undef CONDOR
	#define  LOPTI_CONDOR

	#ifndef		MINIMIZER
		#define	MINIMIZER	condor_minimizer
	#endif
	
	#include <lvv/convert-condor.h>
		using lvv::array;

	#include <stdio.h>
	#include <stdlib.h>
	#include <cmath>
	#include <string>

	#include <condor/Solver.h>
	#include <condor/tools.h>
	#include <condor/ObjectiveFunction.h>
		using CONDOR::sqr;
		using CONDOR::Vector;
		using CONDOR::CONDORSolver;

	#include <lvv/lvv.h>
		using std::cout;
		using std::cerr;
		using std::endl;

	extern int globalPrintLevel;


			template<typename V>
class	c_of_t  : public CONDOR::ObjectiveFunction { public:

					typedef  typename V::value_type fp_t;
					//typedef		loft_base<V>*	loft_v_t;	// virt_ptr to  loft_base
					typedef 	const loft_base<V>&		loft_v_t;		// virt ref

				int eval_cnt;
				bool verbose;
				loft_v_t loft_v;
				//loft_base<V>&  loft_v;

/*	void	init (loft_v_t _loft_v, V& X0) {  			// we use this CTOR to construct
			strcpy(name,"condor_of");
			xOptimal.setSize(X0.size());
			xStart.  setSize(X0.size());
			xStart << X0;
			//loft_v =  _loft_v;
			loft_vp =  &_loft_v;
			//this->var =  var;
			eval_cnt = 0;
			verbose = false;
	};*/

	void	init (V& X0) {  			// we use this CTOR to construct
			xOptimal.setSize(X0.size());
			xStart.  setSize(X0.size());
			xStart << X0;
	};

	explicit	c_of_t	 (loft_v_t _loft_v)    : loft_v(_loft_v)  {  			// we use this CTOR to construct
			strcpy(name,"condor_loft_wrapper");
			//loft_v =  _loft_v;
			//loft_vp =  &_loft_v;
			//this->var =  var;
			eval_cnt = 0;
			verbose = false;
	};

	double  eval (CONDOR::Vector cX, int *nerror=NULL) {  		// condor use this to eval
			V X;	
			X << cX;
			fp_t y = (const_cast< loft_base<V>& > (loft_v))(X);
			updateCounter(y,cX);
			eval_cnt++;

			if (verbose) {
				FMT("%5d \t %19.15g  ")  %eval_cnt  %y;
				for  (typename V::iterator x=X.begin();   x != X.end();   x++)    FMT("\t%19.15g")   %(*x);
				cout << endl;
			}

			return	y;
	};
 };

/////////////////////////////////////////////////////////////// 
                 template<typename V>
class	condor_minimizer: public trust_region_minimizer<V> { public:
				typedef  typename minimizer<V>::fp_t		fp_t;
				typedef 	const loft_base<V>&		loft_v_t;

				using minimizer<V>::X;  		// without this we woun't see minimizer members
				using minimizer<V>::X0;
				using minimizer<V>::max_iter_;
				using minimizer<V>::iter_;
				using minimizer<V>::verbose_;
				using minimizer<V>::verbose_;
				using minimizer<V>::ymin_;
				using minimizer<V>::Xmin_;
				using minimizer<V>::name;
				using minimizer<V>::loft_v;
				using trust_region_minimizer<V>::rho_begin_;
				using trust_region_minimizer<V>::rho_end_;

			// m-vars
			//V	Xmin;
			c_of_t<V>			c_of;	// condor object func wrap	
			//CONDOR::ObjectiveFunction*	c_of;		
			//CONDOR::ObjectiveFunction*	c_rof_wrap;	// condor object func rescaled wrap	

			CONDOR::Vector			cX;		// param in condor format
			CONDOR::Vector			cR;		// rescale divider


			//	void		rescale			(V& R) 		{ cR << R; c_rof_wrap = new CONDOR::CorrectScaleOF(2, c_of, cR); };
	//explicit 			condor_minimizer	(V& _X)			: trust_region_minimizer<V>( _X, "condor") {};
	explicit 			condor_minimizer	(loft_v_t _loft_v)
		:
			trust_region_minimizer<V>( _loft_v, "condor"),
			c_of(loft_v)	// this is cloned loft from base
		{};



	virtual	minimizer<V>&	 	verbose			(bool flag)	{
		verbose_     = flag;
		c_of.verbose = flag;
		return *this;
	};

	virtual V&		 	argmin			()		{
						// TODO assert(!loft_v.empty());  
						assert(!isnan(rho_begin_) && " rho_begin not initialised ");
						assert(!isnan(rho_end_)   && " rho_end   not initialised ");
		c_of.init(X);
		globalPrintLevel = 10;		// off
		CONDOR::CONDORSolver(rho_begin_, rho_end_, max_iter_,  &c_of);
		Xmin_ << (c_of.xBest);
		ymin_ = c_of.valueBest; 				//c_of.printStats();
		iter_ = c_of.getNFE(); 
		return Xmin_;
	};

};
	#endif // LVV_LOPTI_CONDOR_H
