	
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
			function<fp_t(V&)>		oco;

			int eval_cnt;
			bool verbose;

	//c_of_t (fp_t of(V&, void*), V X0, void* var) {  			// we use this CTOR to construct
	void	init (function<fp_t(V&)> of, V& X0) {  			// we use this CTOR to construct
			strcpy(name,"condor_of");
			xOptimal.setSize(X0.size());
			xStart.  setSize(X0.size());
			xStart << X0;
			oco =  of;
			//this->var =  var;
			eval_cnt = 0;
			verbose = false;
	};

	double  eval (CONDOR::Vector cX, int *nerror=NULL) {  		// condor use this to eval
			V X;	
			X << cX;
			fp_t y = oco(X);
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
				//typedef  typename minimizer<V>::function<fp_t(V&)>	function<fp_t(V&)>;


				using minimizer<V>::X;  		// without this we woun't see minimizer members
				using minimizer<V>::max_iter_;
				using minimizer<V>::iter_;
				using minimizer<V>::verbose_;
				using minimizer<V>::oco;
				using minimizer<V>::verbose_;
				using minimizer<V>::ymin_;
				using minimizer<V>::Xmin_;
				using minimizer<V>::name;
				using trust_region_minimizer<V>::rho_begin_;
				using trust_region_minimizer<V>::rho_end_;

	explicit 		condor_minimizer(V& _X) : trust_region_minimizer<V>( _X, "condor") {};


	// m-vars
	//V	Xmin;
	c_of_t<V>			c_of;	// condor object func wrap	
	//CONDOR::ObjectiveFunction*	c_of;		
	//CONDOR::ObjectiveFunction*	c_rof_wrap;	// condor object func rescaled wrap	

	CONDOR::Vector			cX;		// param in condor format
	CONDOR::Vector			cR;		// rescale divider


	//	void		rescale			(V& R) 		{ cR << R; c_rof_wrap = new CONDOR::CorrectScaleOF(2, c_of, cR); };

	virtual minimizer<V>&	 verbose		(bool flag)	{
		#undef GP_F
		#define  GP_F "splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		cout << "# :gnuplot: set view 0,0,1.7;   set font \"arial,6\"; set dgrid3d;  set key off;"
			"  set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 3:4:2:1 with labels; \n";

		verbose_     = flag;
		c_of.verbose = flag;
		return *this;
	};

	virtual V&		 argmin			()		{
						assert(!oco.empty());
						assert(!isnan(rho_begin_) && " rho_begin not initialised ");
						assert(!isnan(rho_end_)   && " rho_end   not initialised ");
		c_of.init(oco, X);
		globalPrintLevel = 10;		// off
		CONDOR::CONDORSolver(rho_begin_, rho_end_, max_iter_,  &c_of);
		Xmin_ << (c_of.xBest);
		ymin_ = c_of.valueBest; 				//c_of.printStats();
		iter_ = c_of.getNFE(); 
		return Xmin_;
	};

};
	#endif // LVV_LOPTI_CONDOR_H
