	// for sleep:
	#include <unistd.h>
	    
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include <string.h>
	#include <condor/Solver.h>
	#include <condor/tools.h>
	#include <condor/ObjectiveFunction.h>
	    
	    

//using namespace CONDOR;
using CONDOR::sqr;
using CONDOR::Vector;
using CONDOR::CONDORSolver;

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
};

		 21 class   condor_of_wrap { public:
		 24         condor_of_wrap(v *F(V), V X0, V S0) {};
		 26 }

                 template<typename V>
class	minimizer { public:
	minimizer(v *F(V), V X0, V S0) {};
	argmin(v *F(V), V X0, V S0) {};
	typedef typeof(*V.begin())      v;
	const int n =   sizeof(V)/sizeof(*V.begin());
};


int main(int argc, char **argv) {
	
	double     rhoStart = 1e-0;
	double     rhoEnd   = 1e-7;
	int        niter    = 1000;
	Rosenbrock of;
				// of->setSaveFile();
	CONDORSolver(rhoStart, rhoEnd, niter, &of);
	of.printStats();
				// FILE * ff = fopen("resultsOne.txt", "w");
				// fprintf(ff, "%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
				// fclose(ff);
	return 0;
 }


