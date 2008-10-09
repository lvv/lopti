	
	#define  OPTI_CONDOR
	#undef CONDOR
	
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
	#include <lvv/array.h>
		using lvv::array;

	extern int globalPrintLevel;


	template<typename V>
class	of_wrap  : public CONDOR::ObjectiveFunction { public:

			typedef  typename V::value_type v;
			v (*of)(V&, void*);
			void* var;
			int eval_cnt;
			bool verbose;

	of_wrap (v of(V&, void*), V X0, void* var) {  			// we use this CTOR to construct
		strcpy(name,"condor_of_wrap");
		xOptimal.setSize(X0.size());
		xStart.  setSize(X0.size());
		xStart << X0;
		this->of =  of;
		this->var =  var;
		eval_cnt = 0;
		verbose = false;
	};

	double  eval (CONDOR::Vector cX, int *nerror=NULL) {  		// condor use this to eval
		V X;	
		X << cX;
		v y = of(X, var);
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
 
                 template<typename V>
class	minimizer { public:
		typedef  typename V::value_type v;
		int				max_iter_;
		//bool				verbose_;
		v 				c_rho_start;	// r(rho) start
		v 				c_rho_end;	// r end
		of_wrap<V>*			c_of_wrap;	// condor object func wrap	
		//CONDOR::ObjectiveFunction*	c_of_wrap;		
		CONDOR::ObjectiveFunction*	c_rof_wrap;	// condor object func rescaled wrap	
		CONDOR::Vector			cX;		// param in condor format
		CONDOR::Vector			cR;		// rescale divider
		V				Xmin;

			minimizer		(v (*of)(V&, void*), V& X, void* _var=NULL)        :max_iter_(500), c_rof_wrap(NULL)  { c_of_wrap = new of_wrap<V>(of, X, _var); };
			~minimizer		()		{ delete   c_of_wrap; };

	void		rescale			(V& R) 		{ cR << R; c_rof_wrap = new CONDOR::CorrectScaleOF(2, c_of_wrap, cR); };

	void		condor_rho_start	(v rho)		{ c_rho_start = rho; };
	void		condor_rho_end		(v rho)		{ c_rho_end   = rho; };
	void		max_iter		(int mx)	{ max_iter_   = mx;  };
	v 	 	ymin			()		{  return c_of_wrap->valueBest; };
	v 	 	iter			()		{  return c_of_wrap->getNFE(); };

	void 	 verbose		(bool flag)	{
		#define  GP_F "splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		cout << "# :gnuplot: set view 0,0,1.7;   set font \"arial,6\"; set dgrid3d;  set key off;"
			"  set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 3:4:2:1 with labels; \n";

		c_of_wrap->verbose =flag;
	};

	V&		 argmin			()		{
						assert(c_of_wrap);
		globalPrintLevel = 10;		// off
		CONDOR::CONDORSolver(c_rho_start, c_rho_end, max_iter_,  c_rof_wrap == NULL ? c_of_wrap : c_rof_wrap);
		Xmin << (c_of_wrap->xBest);
						//c_of_wrap->printStats();
		return Xmin;
	};
};
