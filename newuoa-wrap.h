
#include <lopti.h>
#define LOPTI_NEWUOA
#undef  MINIMIZER                                                                                                                                  
#define	MINIMIZER	newuoa_minimizer
#undef        NEWUOA

// reverse indexes for fortran compat
#define		BMAT(i,j)	BMAT[j][i]
#define		ZMAT(i,j)	ZMAT[j][i]
#define		XPT(i,j)	XPT[j][i]

#include        <cmath>
		using   std::sqrt;  
#include        <algorithm>
		using   std::max;  
		using   std::min;  

#include	<lvv/math.h>
		using lvv::pow2;
		using lvv::pow3;
		using lvv::abs;
#include	<lvv/array.h>
		using lvv::array;
		using lvv::vector;
		using lvv::matrix;

extern "C" void  trsapp_  (int* N, int* NPT, double* XOPT, double* XPT, double* GQ, double* HQ, double* PQ, double* DELTA, double* D, double* W, double* /*W[NP]*/, double* /*W[NP+N]*/, double* /*W[NP+2*N]*/, double* CRVMIN);
extern "C" void  biglag_  (int* N, int* NPT, double* XOPT, double* XPT, double* BMAT, double* ZMAT, int* IDZ, int* NDIM, int* KNEW, double* DSTEP, double* D, double* ALPHA, double* VLAG, double* /*VLAG[NPT+1]*/, double* W, double* /*W[NP]*/, double* /*W[NP+N]*/);
extern "C" void  bigden_  (int* N, int* NPT, double* XOPT, double* XPT, double* BMAT, double* ZMAT, int* IDZ, int* NDIM, int* KOPT, int*  KNEW, double* D, double* W, double* VLAG, double* BETA, double* XNEW, double* /*W[NDIM+1]*/, double* /*W[6*NDIM+1]*/);
extern "C" void  update_  (int* N, int* NPT, double* BMAT, double* ZMAT, int* IDZ, int* NDIM, double* VLAG, double* BETA, int* KNEW, double* W);



		template<typename V, int NPT=2*V::sz+1>
class newuoa_minimizer:  public trust_region_minimizer<V> { public:
	typedef  typename minimizer<V>::fp_t	fp_t;
	typedef  typename minimizer<V>::of_ptr_t	of_ptr_t;

	using minimizer<V>::X;  		// without this we woun't see minimizer members
	using minimizer<V>::max_iter_;
	using minimizer<V>::iter_;
	using minimizer<V>::verbose_;
	using minimizer<V>::of_;
	using minimizer<V>::verbose_;
	using minimizer<V>::ymin_;
	using minimizer<V>::Xmin_;
	using trust_region_minimizer<V>::rho_begin_;
	using trust_region_minimizer<V>::rho_end_;

	explicit 		newuoa_minimizer	(of_ptr_t of, V& _X):   trust_region_minimizer<V>(of, _X)  {};
	virtual const char*	name			()	const 		{ return "newuoa"; };
	virtual V&		argmin			();
};

		template<typename V, int NPT> 	V&
newuoa_minimizer<V,NPT>::argmin () {
						assert(!isnan(rho_begin_) && " rho_begin definition ");
						assert(!isnan(rho_end_)   && " rho_end   definition ");
						assert(V::ibg == 1        && " 1st vector index == 1 "); //  newuoa have index:  [1:N]
	const int N = V::sz;
	const int NP = N+1;
	const int NPTM = NPT-NP;
	const int NDIM = NPT+N;

	if ((NPT < N+2) || ( NPT > ((N+2)*NP)/2))  { cout << "error: NPT is not in the required interval\n"; exit(33); }
	if (verbose_) FMT("olduoa:  N =%d and NPT =%d   ----------------------------------------------------------\n")  % N  % NPT;

	//array<double,N,1>		XBASE;
	//V&	XBASE = *new V;
	V		XBASE;
	V		XOPT; 
	V		XNEW;
	V		GQ;
	V		D;
	vector<double,NPT>		FVAL;
	vector<double,NPT>		PQ;
	vector<double,NDIM>		VLAG;
	vector<double,(N*NP)/2>		HQ;

	matrix<double,NDIM,N>		BMAT;
	matrix<double,NPT,N>		XPT;  
	matrix<double,NPT,NPTM>		ZMAT;   
		
	array<double,1000,1>		W; //ws

	// lvv: fortran auto declared
	double CRVMIN;
	double DSQ;
	double DSTEP;
	double DISTSQ;
	double ALPHA;
	double BETA;
	double BSUM;
	double DELTA;
	double DETRAT;
	double DIFF;
	double DIFFA;
	double DIFFB;
	double DIFFC;
	double DNORM  = 0;
	double DX;
	double FSAVE;
	double F;
	double FOPT;
	double FBEG   = 0;
	double GQSQ;
	double RATIO  = 0;
	double RHO;
	double SUM;
	double SUMA;
	double SUMB;
	double SUMZ;
	double TEMP;
	double TEMPQ;
	double VQUAD;
	double XIPT   = 0;
	double XJPT   = 0;
	double XOPTSQ = 0;
	int    IDZ;
	int    IH;
	int    IP;
	int    IPT    = 0;
	int    ITEMP;
	int    ITEST  = 0;
	int    JPT    = 0;
	int    KNEW;
	int    KSAVE;
	int    KTEMP;
	int    NFSAV  = 0;
	int    KOPT   = 0;

	int    _ndim  = NDIM;                 //const err
	int    _n     = N;
	int    _npt   = NPT;

     	double	HALF   = 0.5          ; 		// newuob start
     	double	ONE    = 1.0          ; 
     	double	TENTH  = 0.1          ; 
     	double	ZERO   = 0.0          ; 
     	int	NH     = (N*NP)/2     ; 


	//  Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.

 	for (int J=1; J<=N; J++)  {
		XBASE[J] = X[J];
	 	for (int K=1; K<=NPT; K++)  	XPT(K,J) = ZERO;
	 	for (int I=1; I<=NDIM; I++)  	BMAT(I,J) = ZERO;
	}

 	for (int IH=1; IH<=NH; IH++)  	HQ[IH] = ZERO;

 	for (int K=1; K<=NPT; K++)  {
		PQ[K] = ZERO;
	 	for (int J=1; J<=NPTM; J++)  	ZMAT(K,J) = ZERO;
	}

     	double	RHOSQ = rho_begin_*rho_begin_;
     	double	RECIP = ONE/RHOSQ;
     	double	RECIQ = sqrt(HALF)/RHOSQ;

     	iter_ = 0;

fill_xpt_50:

	//  Begin the initialization procedure. iter_ becomes one more than the number
	//  of function values so far. The coordinates of the displacement of the
	//  next initial interpolation point from XBASE are set in XPT(iter_,.).

	int  NFM = iter_;
     	int  NFMM = iter_-N;
     	iter_++;
     	if (NFM <= 2*N)  {
     	    if (NFM >= 1  &&  NFM <= N)	XPT(iter_,NFM) = rho_begin_;
     	    else if (NFM > N) 			XPT(iter_,NFMM) = -rho_begin_;
	} else {
     	    ITEMP = (NFMM-1)/N;
     	    JPT = NFM-ITEMP*N-N;
     	    IPT = JPT+ITEMP;
     	    if (IPT > N)  {
     	        ITEMP = JPT;
     	        JPT = IPT-N;
     	        IPT = ITEMP;
     	    }
     	    XIPT = rho_begin_;
     	    if (FVAL[IPT+NP] < FVAL[IPT+1]) XIPT = -XIPT;
     	    XJPT = rho_begin_;
     	    if (FVAL[JPT+NP] < FVAL[JPT+1]) XJPT = -XJPT;
     	    XPT(iter_,IPT) = XIPT;
     	    XPT(iter_,JPT) = XJPT;
     	}

	//  Calculate the next value of F, label 70 being reached immediately
	//  after this calculation. The least function value so far and its index
	//  are required.

 	for (int J=1; J<=N; J++)  	X[J] = XPT(iter_,J)+XBASE[J];

     	goto eval_f_310;

return_to_init_from_eval_70:

   	FVAL[iter_] = F;

     	if (iter_ == 1)  {
     	    FBEG = F;
     	    FOPT = F;
     	    KOPT = 1;
     	} else if (F < FOPT)  {
     	    FOPT = F;
     	    KOPT = iter_;
     	}

	//  Set the nonzero initial elements of BMAT and the quadratic model in
	//  the cases when iter_ is at most 2*N+1.

     	if (NFM <= 2*N)  {

     		if (NFM >= 1  &&  NFM <= N)  {
     		    GQ[NFM] = (F-FBEG)/rho_begin_;
     		    if (NPT < iter_+N)  {
     		        BMAT(1,NFM) = -ONE/rho_begin_;
     		        BMAT(iter_,NFM) = ONE/rho_begin_;
     		        BMAT(NPT+NFM,NFM) = -HALF*RHOSQ;
     		    }

     		} else if (NFM > N)  {
     		    BMAT(iter_-N,NFMM) = HALF/rho_begin_;
     		    BMAT(iter_,NFMM) = -HALF/rho_begin_;
     		    ZMAT(1,NFMM) = -RECIQ-RECIQ;
     		    ZMAT(iter_-N,NFMM) = RECIQ;
     		    ZMAT(iter_,NFMM) = RECIQ;
     		    IH = (NFMM*(NFMM+1))/2;
     		    TEMP = (FBEG-F)/rho_begin_;
     		    HQ[IH] = (GQ[NFMM]-TEMP)/rho_begin_;
     		    GQ[NFMM] = HALF*(GQ[NFMM]+TEMP);
     		}

		//  Set the off-diagonal second derivatives of the Lagrange functions and
		//  the initial quadratic model.

     	} else {
     		IH = (IPT*(IPT-1))/2+JPT;
     		if (XIPT < ZERO) IPT = IPT+N;
     		if (XJPT < ZERO) JPT = JPT+N;
     		ZMAT(1,NFMM) = RECIP;
     		ZMAT(iter_,NFMM) = RECIP;
     		ZMAT(IPT+1,NFMM) = -RECIP;
     		ZMAT(JPT+1,NFMM) = -RECIP;
     		HQ[IH] = (FBEG-FVAL[IPT+1]-FVAL[JPT+1]+F)/(XIPT*XJPT);
     	}

     	if (iter_ < NPT) goto fill_xpt_50;

	//  Begin the iterative procedure, because the initial model is complete.

     	RHO = rho_begin_;
     	DELTA = RHO;
     	IDZ = 1;
     	DIFFA = ZERO;
     	DIFFB = ZERO;
     	ITEST = 0;
     	XOPTSQ = ZERO;

 	for (int I=1; I<=N; I++)  {
		XOPT[I] = XPT(KOPT,I);
	   	XOPTSQ = XOPTSQ+pow2(XOPT[I]);
	}

begin_iter_90:

	NFSAV = iter_;

	//  Generate the next trust region step and test its length. Set KNEW
	//  to -1 if the purpose of the next F will be to improve the model.

gen_tr_100:

  	KNEW = 0;
	trsapp_ (&_n, &_npt, (double*)&XOPT, (double*)&XPT, (double*)&GQ, (double*)&HQ, (double*)&PQ, &DELTA,
		 	(double*)&D, (double*)&W, &W[NP], &W[NP+N], &W[NP+2*N], &CRVMIN);  
     	DSQ = ZERO;
 	for (int I=1; I<=N; I++)  DSQ = DSQ+pow2(D[I]);
     	DNORM = min(DELTA,sqrt(DSQ));
     	if (DNORM < HALF*RHO)  {
     	    KNEW = -1;
     	    DELTA = TENTH*DELTA;
     	    RATIO = -1.0;
     	    if (DELTA <= 1.5*RHO) DELTA = RHO;
     	    if (iter_ <= NFSAV+2) goto L460;
     	    TEMP = 0.125*CRVMIN*RHO*RHO;
     	    if (TEMP <= max(DIFFA, max(DIFFB,DIFFC))) goto L460;
     	    goto new_rho_490;
     	}

shift_xbase_120:
	//  Shift XBASE if XOPT may be too far from XBASE. First make the changes
	//  to BMAT that do not depend on ZMAT.

  	if (DSQ <= 1.0 - 3*XOPTSQ)  {
		TEMPQ = 0.25*XOPTSQ;

 		for (int K=1; K<=NPT; K++)  {
			SUM = ZERO;
	 		for (int I=1; I<=N; I++)  	SUM = SUM+XPT(K,I)*XOPT[I];
			TEMP = PQ[K]*SUM;
			SUM = SUM-HALF*XOPTSQ;
			W[NPT+K] = SUM;
	 		for (int I=1; I<=N; I++)  {
				GQ[I] = GQ[I]+TEMP*XPT(K,I);
				XPT(K,I) = XPT(K,I)-HALF*XOPT[I];
				VLAG[I] = BMAT(K,I);
				W[I] = SUM*XPT(K,I)+TEMPQ*XOPT[I];
				IP = NPT+I;
		 		for (int J=1; J<=I; J++)  	BMAT(IP,J) = BMAT(IP,J)+VLAG[I]*W[J]+W[I]*VLAG[J];
		 	}
		}

		//  Then the revisions of BMAT that depend on ZMAT are calculated.

 		for (int K=1; K<=NPTM; K++)  {
			SUMZ = ZERO;

	 		for (int I=1; I<=NPT; I++)  {
				SUMZ = SUMZ+ZMAT(I,K);
		  		W[I] = W[NPT+I]*ZMAT(I,K);
		  	}

	 		for (int J=1; J<=N; J++)  {
				SUM = TEMPQ*SUMZ*XOPT[J];
		 		for (int I=1; I<=NPT; I++)	SUM += W[I]*XPT(I,J);
				VLAG[J] = SUM;
				if (K < IDZ) SUM = -SUM;
		 		for (int I=1; I<=NPT; I++)  	BMAT(I,J) += SUM*ZMAT(I,K);
		 	}

	 		for (int I=1; I<=N; I++)  {
				IP = I+NPT;
				TEMP = VLAG[I];
				if (K < IDZ) TEMP = -TEMP;
		 		for (int J=1; J<=I; J++)  	BMAT(IP,J) = BMAT(IP,J)+TEMP*VLAG[J];
		 	}
	  	}

		//  The following instructions complete the shift of XBASE, including
		//  the changes to the parameters of the quadratic model.

		IH = 0;

 		for (int J=1; J<=N; J++)  {
			W[J] = ZERO;

	 		for (int K=1; K<=NPT; K++)  {
				W[J] = W[J]+PQ[K]*XPT(K,J);
		  		XPT(K,J) -= HALF*XOPT[J];
		  	}

	 		for (int I=1; I<=J; I++)  {
				IH++;
				if (I < J) GQ[J] += HQ[IH]*XOPT[I];
				GQ[I] += HQ[IH]*XOPT[J];
				HQ[IH] += W[I]*XOPT[J]+XOPT[I]*W[J];
		  		BMAT(NPT+I,J) = BMAT(NPT+J,I);
		  	}
		}

 		for (int J=1; J<=N; J++)  {
			XBASE[J] = XBASE[J]+XOPT[J];
	  		XOPT[J] = ZERO;
	  	}

		XOPTSQ = ZERO;
     	}

	//  Pick the model step if KNEW is positive. A different choice of D
	//  may be made later, if the choice of D by BIGLAG causes substantial
	//  cancellation in DENOM.

     	if (KNEW > 0) 
     	    //CALL BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,DSTEP,D,ALPHA,VLAG,VLAG[NPT+1],W,W[NP],W[NP+N])
	    biglag_ (&_n, &_npt, (double*)&XOPT, (double*)&XPT, (double*)&BMAT, (double*)&ZMAT, &IDZ, &_ndim,
	     &KNEW, &DSTEP, (double*)&D, &ALPHA, (double*)&VLAG, &VLAG[NPT+1], (double*)&W, &W[NP], &W[NP+N]);

	//  Calculate VLAG and BETA for the current choice of D. The first NPT
	//  components of W_check will be held in W.

 	for (int K=1; K<=NPT; K++)  {
		SUMA = ZERO;
		SUMB = ZERO;
		SUM = ZERO;
	 	for (int J=1; J<=N; J++)  {
			SUMA += XPT(K,J)*D[J];
			SUMB += XPT(K,J)*XOPT[J];
		  	SUM  += BMAT(K,J)*D[J];
		}
		W[K] = SUMA*(HALF*SUMA+SUMB);
	  	VLAG[K] = SUM;
	}

     	BETA = ZERO;

 	for (int K=1; K<=NPTM; K++)  {
		SUM = ZERO;
	 	for (int I=1; I<=NPT; I++)  	SUM = SUM+ZMAT(I,K)*W[I];

		if (K < IDZ)  {
		    BETA = BETA+SUM*SUM;
		    SUM = -SUM;
		} else {
		    BETA = BETA-SUM*SUM;
		}
	 	for (int I=1; I<=NPT; I++)  	VLAG[I] += SUM*ZMAT(I,K);
	}

     	BSUM = ZERO;
     	DX = ZERO;

 	for (int J=1; J<=N; J++)  {
		SUM = ZERO;
	 	for (int I=1; I<=NPT; I++)  	SUM = SUM+W[I]*BMAT(I,J);
		BSUM = BSUM+SUM*D[J];
		int JP = NPT+J;
	 	for (int K=1; K<=N; K++)  	SUM = SUM+BMAT(JP,K)*D[K];
		VLAG[JP] = SUM;
		BSUM = BSUM+SUM*D[J];
	  	DX = DX+D[J]*XOPT[J];
	}

     	BETA = DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM;
     	VLAG[KOPT] += ONE;

	//  If KNEW is positive and if the cancellation in DENOM is unacceptable,
	//  then BIGDEN calculates an alternative model step, XNEW being used for
	//  working space.

     	if (KNEW > 0)  {
     	    TEMP = ONE+ALPHA*BETA/pow2(VLAG[KNEW]);
     	    if (abs(TEMP) <= 0.8)  {
     	        //CALL BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,KNEW,D,W,VLAG,BETA,XNEW,W[NDIM+1],W[6*NDIM+1])
		bigden_( &_n, &_npt, (double*)&XOPT, (double*)&XPT, (double*)&BMAT, (double*)&ZMAT, &IDZ, &_ndim, &KOPT,  &KNEW,
		(double*)&D, (double*)&W, (double*)&VLAG, &BETA, (double*)&XNEW, &W[NDIM+1], &W[6*NDIM+1]);

     	    }
	}
	//  Calculate the next value of the objective function.

 xnew_290:
  	for (int I=1; I<=N; I++) {
		XNEW[I] = XOPT[I]+D[I];
	  	X[I] = XBASE[I]+XNEW[I];
	}

     	iter_ = iter_+1;

eval_f_310:

  	if (iter_ > max_iter_)  {
		iter_ = iter_-1;
		if (verbose_) FMT("\n    error:  CALFUN has been called max_iter_ times.");
     	    goto exit_530;
     	}

	// CALL CALFUN (N,X,F)
	// calfun_ (&_n, (double*)&X, &F);
	//F = (*of_)(X, NULL);
	F = of_(X);

	if (verbose_) FMT ("%d \t %18.10g  \t  %18.10g \n")  %iter_  %F  %X;

	if (iter_ <= NPT) goto return_to_init_from_eval_70;
	if (KNEW == -1) goto exit_530;

	//  Use the quadratic model to predict the change in F due to the step D,
	//  and set DIFF to the error of this prediction.

     	VQUAD = ZERO;
     	IH = 0;

 	for (int J=1; J<=N; J++)  {
		VQUAD = VQUAD+D[J]*GQ[J];
	 	for (int I=1; I<=J; I++)  {
			IH++;
			TEMP = D[I]*XNEW[J]+D[J]*XOPT[I];
			if (I == J)   TEMP *= HALF;
		  	VQUAD += TEMP*HQ[IH];
		}
	}

 	for (int K=1; K<=NPT; K++)  	VQUAD += PQ[K]*W[K];

     	DIFF = F-FOPT-VQUAD;
     	DIFFC = DIFFB;
     	DIFFB = DIFFA;
     	DIFFA = abs(DIFF);

     	if (DNORM > RHO) NFSAV = iter_;

	//  Update FOPT and XOPT if the new F is the least value of the objective
	//  function so far. The branch when KNEW is positive occurs if D is not
	//  a trust region step.

     	FSAVE = FOPT;

     	if (F < FOPT)  {
     		FOPT = F;
     		XOPTSQ = ZERO;
 		for (int I=1; I<=N; I++)  {
			XOPT[I] = XNEW[I];
		  	XOPTSQ = XOPTSQ+pow2(XOPT[I]);
		}
     	}

     	KSAVE = KNEW;

     	if (KNEW > 0) goto update_410;

	//  Pick the next value of DELTA after a trust region step.

     	if (VQUAD >= ZERO)  {
		if (verbose_)     FMT ("\n        error: trust region step has failed to reduce Q.");
		goto exit_530;
     	}

     	RATIO = (F-FSAVE)/VQUAD;

     	if 	(RATIO <= TENTH)	DELTA = HALF*DNORM;
     	else if (RATIO <= 0.7)		DELTA = max(HALF*DELTA,DNORM);
     	else				DELTA = max(HALF*DELTA,DNORM+DNORM);

     	if (DELTA <= 1.5*RHO)		DELTA = RHO;

	//  Set KNEW to the index of the next interpolation point to be deleted.

     	RHOSQ = pow2(max(TENTH*DELTA,RHO));
     	KTEMP = 0;
     	DETRAT = ZERO;

     	if (F >= FSAVE) { 
		KTEMP = KOPT;
		DETRAT = ONE;
     	}

 	for (int K=1; K<=NPT; K++)  {

		double HDIAG = ZERO;

	 	for (int J=1; J<=NPTM; J++)  {
			TEMP = ONE;
			if (J < IDZ) TEMP = -ONE;
		  	HDIAG = HDIAG+TEMP*pow2(ZMAT(K,J));
		}

		TEMP = abs(BETA*HDIAG+pow2(VLAG[K]));
		DISTSQ = ZERO;

	 	for (int J=1; J<=N; J++)  	DISTSQ = DISTSQ+pow2((XPT(K,J)-XOPT[J]));

		if (DISTSQ > RHOSQ)		TEMP = TEMP*pow3(DISTSQ/RHOSQ);

		if (TEMP > DETRAT  &&  K != KTEMP)  {
		    DETRAT = TEMP;
		    KNEW = K;
		}
	}

     	if (KNEW == 0) goto L460;

	//  Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
	//  can be moved. Begin the updating of the quadratic model, starting
	//  with the explicit second derivative term.

update_410:

	// CALL UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
	update_  (&_n, &_npt, (double*)&BMAT, (double*)&ZMAT, &IDZ, &_ndim, (double*)&VLAG, &BETA, &KNEW, (double*)&W);

     	FVAL[KNEW] = F;
     	IH = 0;

 	for (int I=1; I<=N; I++)  {
		TEMP = PQ[KNEW]*XPT(KNEW,I);
	 	for (int J=1; J<=I; J++)  {
			IH = IH+1;
		  	HQ[IH] = HQ[IH]+TEMP*XPT(KNEW,J);
		}
	}

     	PQ[KNEW] = ZERO;

	//  Update the other second derivative parameters, and then the gradient
	//  vector of the model. Also include the new interpolation point.

 	for (int J=1; J<=NPTM; J++)  {
		TEMP = DIFF*ZMAT(KNEW,J);
		if (J < IDZ) 			TEMP = -TEMP;
	 	for (int K=1; K<=NPT; K++)  	PQ[K] = PQ[K]+TEMP*ZMAT(K,J);
	}

     	GQSQ = ZERO;
 	for (int I=1; I<=N; I++)  {
		GQ[I] = GQ[I]+DIFF*BMAT(KNEW,I);
		GQSQ = GQSQ+pow2(GQ[I]);
	  	XPT(KNEW,I) = XNEW[I];
	}

	//  If a trust region step makes a small change to the objective function,
	//  then calculate the gradient of the least Frobenius norm interpolant at
	//  XBASE, and store it in W, using VLAG for a vector of right hand sides.

     	if (KSAVE == 0  &&  DELTA == RHO)  {
     		if (abs(RATIO) > 0.01)  {
     			ITEST = 0;
     		} else {
			for (int K=1; K<=NPT; K++)  	VLAG[K] = FVAL[K]-FVAL[KOPT];
     			double GISQ = ZERO;

 			for (int I=1; I<=N; I++)  {
				SUM = ZERO;
	 			for (int K=1; K<=NPT; K++)  	SUM += BMAT(K,I)*VLAG[K];
				GISQ = GISQ+SUM*SUM;
	  			W[I] = SUM;
	  		}

			//  Test whether to replace the new quadratic model by the least Frobenius
			//  norm interpolant, making the replacement if the test is satisfied.

     			ITEST = ITEST+1;
     			if (GQSQ < 100*GISQ) ITEST = 0;

     			if (ITEST >= 3)  {
 				for (int I=1; I<=N; I++)  	GQ[I] = W[I];
 				for (int IH=1; IH<=NH; IH++)  	HQ[IH] = ZERO;

 				for (int J=1; J<=NPTM; J++)  {
					W[J] = ZERO;
	 				for (int K=1; K<=NPT; K++)  	W[J] = W[J]+VLAG[K]*ZMAT(K,J);
	  				if (J < IDZ)			W[J] = -W[J];
	  			}

 				for (int K=1; K<=NPT; K++)  {
					PQ[K] = ZERO;
	 				for (int J=1; J<=NPTM; J++)  	PQ[K] = PQ[K]+ZMAT(K,J)*W[J];
	  			}

				ITEST = 0;
     			}
     		}
     	}

     	if (F < FSAVE)		KOPT = KNEW;

	//  If a trust region step has provided a sufficient decrease in F, then
	//  branch for another trust region calculation. The case KSAVE>0 occurs
	//  when the new function value was calculated by a model step.

     	if (F <= FSAVE+TENTH*VQUAD)	goto gen_tr_100;
     	if (KSAVE > 0)			goto gen_tr_100;

	//  Alternatively, find out if the interpolation points are close enough
	//  to the best point so far.

     	KNEW = 0;

  L460:

  	DISTSQ = 4.0*DELTA*DELTA;

 	for (int K=1; K<=NPT; K++)  {
		SUM = ZERO;
	 	for (int J=1; J<=N; J++)  	SUM += pow2(XPT(K,J)-XOPT[J]);

		if (SUM > DISTSQ)  { 
		    KNEW = K;
		    DISTSQ = SUM;
		}
	}

	//  If KNEW is positive, then set DSTEP, and branch back for the next
	//  iteration, which will generate a "model step".

     	if (KNEW > 0)  {
     	    DSTEP = max(min(TENTH*sqrt(DISTSQ),HALF*DELTA),RHO);
     	    DSQ = DSTEP*DSTEP;
     	    goto shift_xbase_120;
     	}

     	if (RATIO > ZERO)		goto gen_tr_100;
     	if (max(DELTA,DNORM) > RHO)	goto gen_tr_100;


new_rho_490: 

	//  The calculations with the current value of RHO are complete. Pick the
	//  next values of RHO and DELTA.

  	if (RHO > rho_end_)  { 

		DELTA = HALF*RHO;
		RATIO = RHO/rho_end_;
		
		if	(RATIO <= 16.0)	RHO = rho_end_;
		else if	(RATIO <= 250.0)	RHO = sqrt(RATIO)*rho_end_;
		else				RHO = TENTH*RHO;
		
		DELTA = max(DELTA,RHO);

		if (verbose_) FMT("-- (%d) RHO =%9.6g \t F =%9.6g   X%.8g \n")   %iter_  %RHO  %FOPT  %X;
		goto  begin_iter_90; 
     	}

	//  Return from the calculation, after another Newton-Raphson step, if
	//  it is too short to have been tried before.

     	if (KNEW == -1)		goto xnew_290;

  exit_530:

	if (FOPT <= F)  {
	 	for (int I=1; I<=N; I++)  	X[I] = XBASE[I]+XOPT[I];
		F = FOPT;
	}

	if (verbose_)  FMT("-- (%d) RETURNED: \t F =%.15g    X is: %.15g\n\n")  %iter_ %F  %X;
	ymin_ = F;
	Xmin_ = X;
	return X;
}; // newuoa

