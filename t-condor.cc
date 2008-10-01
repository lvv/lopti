
	// for sleep:
	#include <unistd.h>
	    
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include <string.h>

	#include <lvv/lvv.h>
	#include <lvv/math.h>
	using lvv::pow2;
	#include <lvv/array.h>
	using lvv::array;

typedef array<double,2>		array_t;	
//template<typename V> typeof(*V::begin())  of_rb(V X)   {  double y=100*pow2(X[1]-pow2(X[0]))+pow2(1-X[0]);  FMT("y=%g %15t x%f \n") %y %X; return y; };
template<typename V>  typename V::value_type  of_rb(V X)   {  double y=100*pow2(X[1]-pow2(X[0]))+pow2(1-X[0]);  FMT("y=%g %15t x%f \n") %y %X; return y; };


///////////////////////////////////////////////////////////////////////////////////  CONDOR specific
	#include <condor/Solver.h>
	#include <condor/tools.h>
	#include <condor/ObjectiveFunction.h>
	using CONDOR::sqr;
	using CONDOR::Vector;
	using CONDOR::CONDORSolver;


							/*
							class Rosenbrock : public CONDOR::ObjectiveFunction { public:

								Rosenbrock() {
									strcpy(name,"ROSEN");
									
									xOptimal.setSize(2);
									xStart.setSize(2);
									
									xOptimal[0]=1.0;
									xOptimal[1]=1.0;
									
									valueOptimal=0.0;
									
									xStart[0]=-1.2;
									xStart[1]=1.0;
								}

								double eval(CONDOR::Vector X, int *nerror=NULL) {
									double *x=X, r=100*sqr(x[1]-sqr(x[0]))+sqr(1-x[0]);
									updateCounter(r,X);
									return r;
								}
							}; */

	template<typename V>
class	of_wrap  : public CONDOR::ObjectiveFunction { public:

		typedef  typename V::value_type v;
		v (*F)(V);

	of_wrap (v F(V), V X0) {
		strcpy(name,"condor_of_wrap");
		xOptimal.setSize(X0.size());
		xStart.  setSize(X0.size());
		for (int i=X0.ibegin();  i < X0.iend();  i++)      xStart[i] = X0[i];
		this->F =  F;
	};


			//typeof(*X.begin())	eval (V X) {
	double eval(CONDOR::Vector cX, int *nerror=NULL) {
		V X;	
		for (int i=X.ibegin();  i < X.iend();  i++)      X[i] = cX[i];
		v y = F(X);
		updateCounter(y,cX);
		return	y;
	};
 };

                 template<typename V>
class	minimizer { public:
		typedef  typename V::value_type v;
		int				max_iter_;
		//V				S;
		v 				c_rho_start;
		v 				c_rho_end;
		CONDOR::ObjectiveFunction*	c_of_wrap;		
		CONDOR::Vector			cX;
		V				Xmin;
		//typedef  typeof(*Xmin.begin()) v;

	minimizer	(v (*F)(V), V X)
		:max_iter_(500)
	{
		//const int			n =   sizeof(X)/sizeof(*X.begin());
		c_of_wrap = new of_wrap<array_t>(F, X);
									//strcpy(name,"ROSEN");
									//xOptimal[0]=1.0;
									//xOptimal[1]=1.0;
									//valueOptimal=0.0;
	};

	~minimizer	()		{ delete   c_of_wrap; };

	minimizer& 	 condor_rho_start	(v rho)	{ c_rho_start = rho;	return *this; };
	minimizer& 	 condor_rho_end		(v rho)	{ c_rho_end   = rho;	return *this; };

	minimizer& 	 max_iter		(int mx) { max_iter_ = mx;	return *this; };
	V&		 argmin			()	{
						assert(c_of_wrap);
						// rescale
						//CorrectScaleOF(int _t, ObjectiveFunction *_of, Vector _rescaling);
						//CorrectScaleOF(int _t, ObjectiveFunction *_of);
		CONDOR::CONDORSolver(c_rho_start, c_rho_end, max_iter_, c_of_wrap);
		for (int i = Xmin.ibegin();  i<Xmin.iend();  i++)  Xmin[i] = c_of_wrap->xBest[i];
		c_of_wrap->printStats();
		return Xmin;

	};
};

///////////////////////////////////////////////////////////////////////////////////////////////  GENERIC
int main(int argc, char **argv) {
	
	//double     rhoStart = 1e-0;
	//double     rhoEnd   = 1e-7;
	//int        max_iter    = 1000;
	//Rosenbrock of;
				// of->setSaveFile();
	array_t		X0 = {{ -1.2, 1}};

	minimizer<array_t>	mzr(of_rb, X0);
	mzr.condor_rho_start	(1);
	mzr.condor_rho_end	(1e-7);
	array_t	Xmin = mzr.argmin();
	cout << "Result: " << Xmin << endl;

				//CONDORSolver(rhoStart, rhoEnd, niter, &of);
				//of_wrap.printStats();
				// FILE * ff = fopen("resultsOne.txt", "w");
				// fprintf(ff, "%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
				// fclose(ff);
	return 0;
 }


