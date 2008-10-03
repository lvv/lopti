	extern int globalPrintLevel;
	
	#include <condor/Solver.h>
	#include <condor/tools.h>
	#include <condor/ObjectiveFunction.h>
		using CONDOR::sqr;
		using CONDOR::Vector;
		using CONDOR::CONDORSolver;



	template<typename V>
class	of_wrap  : public CONDOR::ObjectiveFunction { public:

			typedef  typename V::value_type v;
			v (*F)(V);
			int eval_cnt;
			bool verbose;

	of_wrap (v F(V), V X0) {
		strcpy(name,"condor_of_wrap");
		xOptimal.setSize(X0.size());
		xStart.  setSize(X0.size());
		for (int i=X0.ibegin();  i < X0.iend();  i++)      xStart[i] = X0[i];
		this->F =  F;
		eval_cnt = 0;
		verbose = false;
	};

			//typeof(*X.begin())	eval (V X) {
	double  eval (CONDOR::Vector cX, int *nerror=NULL) {
		V X;	
		for (int i=X.ibegin();  i < X.iend();  i++)      X[i] = cX[i];
		v y = F(X);
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
		v 				c_rho_start;
		v 				c_rho_end;
		of_wrap<V>*			c_of_wrap;		
		//CONDOR::ObjectiveFunction*	c_of_wrap;		
		CONDOR::ObjectiveFunction*	c_rof_wrap;		
		CONDOR::Vector			cX;
		CONDOR::Vector			cR;
		V				Xmin;

			minimizer		(v (*F)(V), V X)      :max_iter_(500), c_rof_wrap(NULL)  { c_of_wrap = new of_wrap<array_t>(F, X); };
			~minimizer		()		{ delete   c_of_wrap; };

	minimizer&	rescale			(V R)  {
		cR.setSize(R.size());
		for (int i = R.ibegin();  i<R.iend();  i++)     cR[i] = R[i];
		c_rof_wrap = new CONDOR::CorrectScaleOF(2, c_of_wrap, cR);
		return *this;
	};


	minimizer& 	condor_rho_start	(v rho)		{ c_rho_start = rho;	return *this; };
	minimizer& 	condor_rho_end		(v rho)		{ c_rho_end   = rho;	return *this; };
	minimizer& 	max_iter		(int mx)	{ max_iter_   = mx;	return *this; };
	v 	 	ymin			()		{  return c_of_wrap->valueBest; };
	v 	 	iter			()		{  return c_of_wrap->getNFE(); };

	minimizer& 	 verbose		(bool flag)	{
		#define  GP_F "splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		cout << "# :gnuplot: set view 0,0,1.7;   set font \"arial,6\"; set dgrid3d;  set key off;"
			"  set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 3:4:2:1 with labels; \n";

		c_of_wrap->verbose =flag;	return *this;
	};

	V&		 argmin			()		{
						assert(c_of_wrap);
		globalPrintLevel = 10;		// off
		CONDOR::CONDORSolver(c_rho_start, c_rho_end, max_iter_,  c_rof_wrap == NULL ? c_of_wrap : c_rof_wrap);
		for (int i = Xmin.ibegin();  i<Xmin.iend();  i++)     Xmin[i] = c_of_wrap->xBest[i];
						//c_of_wrap->printStats();
		return Xmin;

	};
};

///////////////////////////////////////////////////////////////////////////////////////////////  GENERIC
