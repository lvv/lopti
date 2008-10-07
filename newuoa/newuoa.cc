// g/\]\[/s/\vBMAT *\[(.{-})\]\[(.{-})\]/BMAT(\2,\1)/g

#define BMAT(i,j)  BMAT[j][i]
#define ZMAT(i,j)  ZMAT[j][i]
#define XPT(i,j)  XPT[j][i]

#include        <cmath>
	using   std::sqrt;  
#include        <algorithm>
	using   std::max;  
	using   std::min;  
	using   std::swap;  

#include <lvv/lvv.h>
#include <lvv/math.h>
	using lvv::pow2;
	using lvv::pow3;
	using lvv::abs;
#include <lvv/array.h>
	using lvv::array;

extern "C" void  trsapp_  (int* N, int* NPT, double* XOPT, double* XPT, double* GQ, double* HQ, double* PQ, double* DELTA, double* D, double* W, double* /*W[NP]*/, double* /*W[NP+N]*/, double* /*W[NP+2*N]*/, double* CRVMIN);
extern "C" void  biglag_  (int* N, int* NPT, double* XOPT, double* XPT, double* BMAT, double* ZMAT, int* IDZ, int* NDIM, int* KNEW, double* DSTEP, double* D, double* ALPHA, double* VLAG, double* /*VLAG[NPT+1]*/, double* W, double* /*W[NP]*/, double* /*W[NP+N]*/);
extern "C" void  bigden_  (int* N, int* NPT, double* XOPT, double* XPT, double* BMAT, double* ZMAT, int* IDZ, int* NDIM, int* KOPT, int*  KNEW, double* D, double* W, double* VLAG, double* BETA, double* XNEW, double* /*W[NDIM+1]*/, double* /*W[6*NDIM+1]*/);
extern "C" void  update_  (int* N, int* NPT, double* BMAT, double* ZMAT, int* IDZ, int* NDIM, double* VLAG, double* BETA, int* KNEW, double* W);
extern "C" void  calfun_  (int* N, double* X, double* F);



		template<int N, int NPT>
void newuoa (array<double,N, 1>& X,  double RHOBEG,  double RHOEND,  int IPRINT,  int MAXFUN) {

	const int NP=N+1;
	const int NPTM=NPT-NP;

	if ((NPT < N+2) || ( NPT > ((N+2)*NP)/2))  { cout << "Return from NEWUOA because NPT is not in the required interval\n"; exit(33); }

	// work space (former W)
	const int NDIM=NPT+N;
	//W[IXB], W(IXO), W(IXN), W(IXP), W(IFV), W(IGQ), W(IHQ), W(IPQ), W(IBMAT), W(IZMAT), NDIM, W(ID), W(IVL), W(IW)
	//XBASE , XOPT  , XNEW  , XPT   , FVAL  , GQ    , HQ    , PQ    , BMAT    , ZMAT    , NDIM, D    , VLAG  , W)
	array<double,N,1>		XBASE;
	array<double,N,1>		XOPT; 
	array<double,N,1>		XNEW;
	array<double,NPT,1>		FVAL;
	array<double,N,1>		GQ;
	array<double,(N*NP)/2,1>	HQ;
	array<double,NPT,1>		PQ;
	array<double,N,1>		D;
	array<double,NDIM,1>		VLAG;

	//  DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),FVAL(*),  GQ(*),HQ(*),PQ(*),D(*),VLAG(*),W(*)
	//  DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),XPT(NPT,*)
	array<array<double,N,1>,NPT,1>	XPT;  
	array<array<double,N,1>,NDIM,1>	BMAT;
	array<array<double,NPTM,1>,NPT,1>ZMAT;   
		
	array<double,1000,1>		W; //ws


	//  Set some constants.

     	const double 	HALF=0.5;
     	const double 	ONE=1.0;
     	const double 	TENTH=0.1;
     	const double 	ZERO=0.0;
     	const int	NH=(N*NP)/2;
     	const int 	NFTEST=max(MAXFUN,1);

	//  Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.

	for(int J=1; J<=N;  J++)  {
		XBASE[J]=X[J];
		for (int K=1;  K<=NPT;  K++) 	XPT(J,K)  = ZERO;
     		for (int I=1; I<=NDIM; I++) 	BMAT(J,I) = ZERO;
   	}

	for(int IH=1; IH<=NH;  IH++) 	HQ[IH]=ZERO;

	for(int K=1; K<=NPT; K++)  { 
		PQ[K]=ZERO;
		for(int J=1; J<=NPTM; J++)  ZMAT(J,K)=ZERO;
   	}

	double  RHOSQ=RHOBEG*RHOBEG;
	double  RECIP=ONE/RHOSQ;
	double  RECIQ=sqrt(HALF)/RHOSQ;
	int NF=0;
	
	// lvv: comp error
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
	double TEMP;
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

	//  Begin the initialization procedure. NF becomes one more than the number
	//  of function values so far. The coordinates of the displacement of the
	//  next initial interpolation point from XBASE are set in XPT(NF,.).

fill_xpt_50: //do {
		int NFM=NF;
		int NFMM=NF-N;
		NF++;

		if (NFM <= 2*N) {
			if (NFM >= 1 && NFM <= N) {
				XPT(NFM,NF) = RHOBEG;
			} else if (NFM > N) {
				XPT(NFMM,NF) = -RHOBEG;
			}
		} else {
			ITEMP=(NFMM-1)/N;
			JPT=NFM-ITEMP*N-N;
			IPT=JPT+ITEMP;

			if (IPT > N) {
			    ITEMP=JPT;
			    JPT=IPT-N;
			    IPT=ITEMP;
			}

			XIPT=RHOBEG;
			if (FVAL[IPT+NP] < FVAL[IPT+1]) XIPT=-XIPT;
			XJPT=RHOBEG;
			if (FVAL[JPT+NP] < FVAL[JPT+1]) XJPT=-XJPT;
			XPT(IPT,NF)=XIPT;
			XPT(JPT,NF)=XJPT;
		}
		
		//     Calculate the next value of F, label 70 being reached immediately
		//     after this calculation. The least function value so far and its index  are required.
		
		for (int J=1; J<=N; J++)  X[J]=XPT(J,NF)+XBASE[J];

		goto eval_f_310;

return_to_init_from_eval_70:

		FVAL[NF]=F;

		if (NF == 1) {
			FBEG=F;
			FOPT=F;
			KOPT=1;
		} else if (F < FOPT) {
			FOPT=F;
			KOPT=NF;
		}

		
		//     Set the nonzero initial elements of BMAT and the quadratic model in
		//     the cases when NF is at most 2*N+1.
		
		if (NFM <= 2*N) {
			if (NFM >= 1  &&  NFM <= N)  {
				GQ[NFM] = (F-FBEG)/RHOBEG;
				if  (NPT < NF+N)  {
				  	BMAT(NFM,1)       = -ONE/RHOBEG;
				  	BMAT(NFM,NF)      = ONE/RHOBEG ;
				  	BMAT(NFM,NPT+NFM) = -HALF*RHOSQ;
				}
			} else if (NFM > N) {
				BMAT(NFMM,NF-N) = HALF/RHOBEG           ; 
				BMAT(NFMM,NF)   = -HALF/RHOBEG          ; 
				ZMAT(NFMM,1)    = -RECIQ-RECIQ          ; 
				ZMAT(NFMM,NF-N) = RECIQ                 ; 
				ZMAT(NFMM,NF)   = RECIQ                 ; 
				IH           = (NFMM*(NFMM+1))/2     ; 
				double TEMP      = (FBEG-F)/RHOBEG       ; 
				HQ[IH]           = (GQ[NFMM]-TEMP)/RHOBEG; 
				GQ[NFMM]         = HALF*(GQ[NFMM]+TEMP)  ; 
			}
			
		//     Set the off-diagonal second derivatives of the Lagrange functions and
		//     the initial quadratic model.
		
		} else {
			IH=(IPT*(IPT-1))/2+JPT;
			if (XIPT < ZERO) IPT = IPT+N;
			if (XJPT < ZERO) JPT = JPT+N;
			ZMAT(NFMM,1)     = RECIP;
			ZMAT(NFMM,NF)    = RECIP;
			ZMAT(NFMM,IPT+1) = -RECIP;
			ZMAT(NFMM,JPT+1) = -RECIP;
			HQ[IH]            = (FBEG-FVAL[IPT+1]-FVAL[JPT+1]+F)/(XIPT*XJPT);
		}

      if (NF < NPT) goto fill_xpt_50;



//     Begin the iterative procedure, because the initial model is complete.


	RHO    = RHOBEG;
	DELTA  = RHO;
	IDZ    = 1;
	DIFFA  = ZERO;
	DIFFB  = ZERO;
	ITEST  = 0;
	XOPTSQ = ZERO;
	
	for (int I=1; I<=N; I++)  {
		XOPT[I] = XPT(I,KOPT); 
		XOPTSQ  = XOPTSQ+pow2(XOPT[I]); 
	}

begin_iter_90:
	NFSAV=NF;

	//  Generate the next trust region step and test its length. Set KNEW
	//  to -1 if the purpose of the next F will be to improve the model.

gen_tr_100:

	KNEW=0;
	double CRVMIN;
	double DSQ;
	trsapp_ (&_n, &_npt, (double*)&XOPT, (double*)&XPT, (double*)&GQ, (double*)&HQ, (double*)&PQ, &DELTA, (double*)&D, (double*)&W, &W[NP], &W[NP+N], &W[NP+2*N], &CRVMIN);
	DSQ=ZERO;
	for (int I=1; I<=N; I++)   DSQ += pow2(D[I]);

	DNORM = min(DELTA,sqrt(DSQ));

	if (DNORM < HALF*RHO)  {
	    KNEW=-1;
	    DELTA=TENTH*DELTA;
	    RATIO=-1.0;
	    if (DELTA <= 1.5e0*RHO) DELTA=RHO;
	    if (NF <= NFSAV+2) goto L460;
	    TEMP = 0.125 * CRVMIN * RHO * RHO;
	    if (TEMP <= max(DIFFA,max(DIFFB,DIFFC))) goto L460;
	    goto new_rho_490;
	}

//	Shift XBASE if XOPT may be too far from XBASE. First make the changes
//	to BMAT that do not depend on ZMAT.

shift_xbase_120:

	if (DSQ <= 1.0e-3 * XOPTSQ)  {

        	double TEMPQ = 0.25 * XOPTSQ;

        	for (int K=1;  K<=NPT;  K++) {
			SUM=ZERO;
			for (int  I=1;  I<=N;  I++)    	 SUM  +=  XPT(I,K) * XOPT[I];
			TEMP=PQ[K]*SUM;
			SUM=SUM-HALF*XOPTSQ;
			W[NPT+K]=SUM;
			for (int I=1;  I<=N;  I++)  {
				GQ[I]=GQ[I]+TEMP*XPT(I,K);
				XPT(I,K)=XPT(I,K)-HALF*XOPT[I];
				VLAG[I]=BMAT(I,K);
				W[I]=SUM*XPT(I,K)+TEMPQ*XOPT[I];
				IP=NPT+I;
				for (int J=1; J<=I; J++)   BMAT(J,IP)  +=  VLAG[I]*W[J]+W[I]*VLAG[J];
			}
		}

		//  Then the revisions of BMAT that depend on ZMAT are calculated.

        	for (int K=1;  K<=NPTM; K++)  {

			double SUMZ=ZERO;

			for (int I=1; I<=NPT; I++)  {
				SUMZ=SUMZ+ZMAT(K,I);
		  		W[I]=W[NPT+I]*ZMAT(K,I);
		  	}

			for (int J=1; J<=N; J++)  {
				SUM=TEMPQ*SUMZ*XOPT[J];
				for (int I=1; I<=NPT;  I++)    SUM=SUM+W[I]*XPT(J,I);
				VLAG[J]=SUM;
				if (K < IDZ)  SUM=-SUM;
				for (int I=1; I<=NPT; I++) 	BMAT(J,I) +=  SUM*ZMAT(K,I);
		  	}

			for (int I=1; I<=N; I++)  {
				IP=I+NPT;
				TEMP=VLAG[I];
				if (K < IDZ) TEMP=-TEMP;
				for (int J=1; J<=I; J++) 	BMAT(J,IP) += TEMP*VLAG[J];
			}
		}

		//  The following instructions complete the shift of XBASE, including
		//  the changes to the parameters of the quadratic model.

     		int IH=0;

     		for (int J=1; J<=N; J++) {
			W[J]=ZERO;

			for (int K=1; K<=NPT; K++) {
				W[J]=W[J]+PQ[K]*XPT(J,K);
		  		XPT(J,K)=XPT(J,K)-HALF*XOPT[J];
		  	}

			for (int I=1; I<=J; I++)  {
				IH=IH+1;
				if (I < J) GQ[J]=GQ[J]+HQ[IH]*XOPT[I];
				GQ[I] += HQ[IH] * XOPT[J];
				HQ[IH] += W[I]*XOPT[J]+XOPT[I]*W[J];
		  		BMAT(J,NPT+I)=BMAT(I,NPT+J);
		  	}
	  	}

     		for (int J=1; J<=N; J++)  { 
			XBASE[J] = XBASE[J]+XOPT[J]; 
			XOPT[J]  = ZERO;
	  	}

     		XOPTSQ=ZERO;
	}


	//  Pick the model step if KNEW is positive. A different choice of D
	//  may be made later, if the choice of D by BIGLAG causes substantial
	//  cancellation in DENOM.

	double DSTEP;
	double ALPHA;
	if (KNEW > 0) {
          	biglag_ (&_n, &_npt, (double*)&XOPT, (double*)&XPT, (double*)&BMAT, (double*)&ZMAT, &IDZ, &_ndim, &KNEW, &DSTEP, (double*)&D, &ALPHA, (double*)&VLAG, &VLAG[NPT+1], (double*)&W, &W[NP], &W[NP+N]);
	}

	//  Calculate VLAG and BETA for the current choice of D. The first NPT
	//  components of W_check will be held in W.

	for (int K=1; K<=NPT; K++) {
		double SUMA=ZERO;
		double SUMB=ZERO;
		SUM=ZERO;
		for (int J=1; J<=N; J++) {
			SUMA=SUMA+XPT(J,K)*D[J];
			SUMB=SUMB+XPT(J,K)*XOPT[J];
			SUM=SUM+BMAT(J,K)*D[J];
		}
		W[K]=SUMA*(HALF*SUMA+SUMB);
	  	VLAG[K]=SUM;
	}

     	BETA=ZERO;
     	for (int K=1; K<=NPTM; K++)  {
		SUM=ZERO;
		for (int I=1; I<=NPT; I++) SUM += ZMAT(K,I)*W[I];

		if (K < IDZ)  {
		    BETA +=  SUM*SUM;
		    SUM   = -SUM;
		} else {
		    BETA -= SUM*SUM;
		}

		for (int I=1; I<=NPT; I++) 	VLAG[I] += SUM*ZMAT(K,I);
	}

     	BSUM = ZERO;
     	DX   = ZERO;
     	int    JP;

     	for (int J=1; J<=N; J++) {
		SUM=ZERO;
		for (int I=1; I<=NPT; I++)	SUM += W[I]*BMAT(J,I);

		BSUM += SUM*D[J];
		JP=NPT+J;
		for (int K=1; K<=N; K++) 	SUM += BMAT(K,JP)*D[K];
		VLAG[JP]=SUM;
		BSUM += SUM*D[J];
	  	DX += D[J]*XOPT[J];
	}

     	BETA = DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM;
     	VLAG[KOPT] = VLAG[KOPT]+ONE;

	//  If KNEW is positive and if the cancellation in DENOM is unacceptable,
	//  then BIGDEN calculates an alternative model step, XNEW being used for  working space.

     	if (KNEW > 0)  {
     	    TEMP = ONE+ALPHA*BETA/pow2(VLAG[KNEW]);
     	    if (abs(TEMP) <= 0.8)
     	        bigden_ (
			&_n, &_npt, (double*)&XOPT, (double*)&XPT, (double*)&BMAT, (double*)&ZMAT, &IDZ, &_ndim, &KOPT,  &KNEW,
			(double*)&D, (double*)&W, (double*)&VLAG, &BETA, (double*)&XNEW, &W[NDIM+1], &W[6*NDIM+1]
		);
     	}

	//  Calculate the next value of the objective function.

xnew_290:

	for (int I=1; I<=N; I++) {
		XNEW[I] = XOPT[I]+D[I]    ; 
		X[I]    = XBASE[I]+XNEW[I]; 
	}

     	NF++;

eval_f_310:

 	if (NF > NFTEST)  {
     		NF--;
     		if (IPRINT > 0) FMT("\n    Return from NEWUOA because CALFUN has been called MAXFUN times.");
     		goto exit_530;
     	}


     	//CALL CALFUN (N,X,F)
     	calfun_ (&_n, (double*)&X, &F);

     	if (IPRINT == 3) FMT ("\n	Function number %d    F =%18.10g    The corresponding X is:  %18.10g \n")  %NF  %F  %X;

     	if (NF <= NPT) goto return_to_init_from_eval_70;
     	if (KNEW == -1) goto exit_530;
//
//    	Use the quadratic model to predict the change in F due to the step D,
//    	and set DIFF to the error of this prediction.
//
     	VQUAD=ZERO;
     	IH=0;
     	for (int J=1; J<=N; J++) {
		VQUAD += D[J]*GQ[J];
		for (int I=1; I<=J; I++) {
			IH=IH+1;
			TEMP=D[I]*XNEW[J]+D[J]*XOPT[I];
			if (I == J) TEMP=HALF*TEMP;
			VQUAD += TEMP*HQ[IH];
		}
	}

     	for (int K=1; K<=NPT; K++)	VQUAD += PQ[K]*W[K];

     	DIFF  = F-FOPT-VQUAD; 
     	DIFFC = DIFFB       ; 
     	DIFFB = DIFFA       ; 
     	DIFFA = abs(DIFF)   ; 
     	if (DNORM > RHO) NFSAV=NF;

	//  Update FOPT and XOPT if the new F is the least value of the objective
	//  function so far. The branch when KNEW is positive occurs if D is not
	//  a trust region step.

     	//KTEMP = 0;
     	//DETRAT=ZERO;

     	FSAVE=FOPT;
     	if (F < FOPT)  {
     		FOPT=F;
     		XOPTSQ=ZERO;
     		for (int I=1; I<=N; I++)  {
			XOPT[I] = XNEW[I];
			XOPTSQ += pow2(XOPT[I]);
		}
     	}

     	KSAVE=KNEW;
     	if (KNEW > 0) goto update_410;

	//  Pick the next value of DELTA after a trust region step.

     	if (VQUAD >= ZERO) {
     	    if (IPRINT > 0) 	FMT ("\n	Return from NEWUOA because a trust region step has failed to reduce Q.");
     	    goto exit_530;
     	}

     	RATIO=(F-FSAVE)/VQUAD;
     	if (RATIO <= TENTH) 
     	    DELTA = HALF*DNORM;
     	else if (RATIO <= 0.7)
     	    DELTA = max (HALF*DELTA,DNORM);
     	else
     	    DELTA = max (HALF*DELTA,DNORM+DNORM);
     	
     	if (DELTA <= 1.5*RHO) DELTA=RHO;

	//    	Set KNEW to the index of the next interpolation point to be deleted.

     	RHOSQ = pow2(max(TENTH*DELTA,RHO));
     	KTEMP = 0;
     	DETRAT=ZERO;
	double  HDIAG;
	double  DISTSQ;

     	if (F >= FSAVE) {
     	    KTEMP  = KOPT; 
     	    DETRAT = ONE ; 
     	}

     	for (int K=1; K<=NPT; K++) {
		HDIAG=ZERO;

		for (int J=1; J<=NPTM; J++) {
			TEMP=ONE;
			if (J < IDZ) TEMP=-ONE;
		  	HDIAG += TEMP*pow2(ZMAT(J,K));
		}

		TEMP=abs(BETA*HDIAG+pow2(VLAG[K]));
		DISTSQ=ZERO;

		for (int J=1; J<=N; J++) 	DISTSQ += pow2(XPT(J,K)-XOPT[J]);

		if (DISTSQ > RHOSQ) TEMP *= pow3(DISTSQ/RHOSQ);
		if (TEMP > DETRAT &&  K != KTEMP)  {
		    DETRAT=TEMP;
		    KNEW=K;
		}
	}

     	if (KNEW == 0) goto L460;

//    	Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
//    	can be moved. Begin the updating of the quadratic model, starting
//    	with the explicit second derivative term.

update_410:

	//CALL UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
	update_  (&_n, &_npt, (double*)&BMAT, (double*)&ZMAT, &IDZ, &_ndim, (double*)&VLAG, &BETA, &KNEW, (double*)&W);

     	FVAL[KNEW]=F;
     	IH=0;
     	for (int I=1; I<=N; I++)  {
		TEMP = PQ[KNEW]*XPT(I,KNEW);
		for (int J=1; J<=I; J++)  {
			IH++;
		  	HQ[IH] += TEMP*XPT(J,KNEW);
		}
	}

     	PQ[KNEW] = ZERO;

	//  Update the other second derivative parameters, and then the gradient
	//  vector of the model. Also include the new interpolation point.

     	for (int J=1; J<=NPTM; J++)  {
		TEMP = DIFF*ZMAT(J,KNEW);
		if (J < IDZ) TEMP=-TEMP;
		for (int K=1; K<=NPT; K++) 	PQ[K] += TEMP*ZMAT(J,K);
	}

     	GQSQ=ZERO; 

     	for (int I=1; I<=N; I++) {
		GQ[I] += DIFF*BMAT(I,KNEW);
		GQSQ  += pow2(GQ[I]);
	  	XPT(I,KNEW) = XNEW[I];
	}

	//  If a trust region step makes a small change to the objective function,
	//  then calculate the gradient of the least Frobenius norm interpolant at
	//  XBASE, and store it in W, using VLAG for a vector of right hand sides.

     	if (KSAVE == 0   &&  DELTA == RHO)   {

		if (abs(RATIO) > 0.01) {
			ITEST=0;

     	    	} else {

     			for (int K=1; K<=NPT; K++)   VLAG[K] = FVAL[K]-FVAL[KOPT];
     			double GISQ=ZERO;

     			for (int I=1; I<=N; I++)  {
				SUM=ZERO;
				for (int K=1; K<=NPT; K++) 	SUM += BMAT(I,K)*VLAG[K];
				GISQ += SUM*SUM;
	  			W[I] = SUM;
	  		}
     
     			//    	Test whether to replace the new quadratic model by the least Frobenius
     			//    	norm interpolant, making the replacement if the test is satisfied.
     		
     			ITEST++;

     			if (GQSQ < 100*GISQ) ITEST=0;

     			if (ITEST >= 3) {

     				for (int I=1; I<=N; I++)	GQ[I] = W[I];
     				for (int IH=1; IH<=NH; IH++)	HQ[IH] = ZERO;

     				for (int J=1; J<=NPTM; J++) {
					W[J]=ZERO;
					for (int K=1; K<=NPT; K++) 	W[J] += VLAG[K]*ZMAT(J,K);
	  				if (J < IDZ) W[J] =- W[J];
	  			}

     				for (int K=1; K<=NPT; K++) {
					PQ[K] = ZERO;
					for (int J=1; J<=NPTM; J++)	PQ[K] += ZMAT(J,K)*W[J];
				}

     				ITEST=0;
     			}
     	    	}
     	}


     	if (F < FSAVE)  KOPT=KNEW;

	//  If a trust region step has provided a sufficient decrease in F, then
	//  branch for another trust region calculation. The case KSAVE>0 occurs
	//  when the new function value was calculated by a model step.

	if (F <= FSAVE+TENTH*VQUAD) goto gen_tr_100;
	if (KSAVE > 0) goto gen_tr_100;

	//  Alternatively, find out if the interpolation points are close enough
	//  to the best point so far.

     	KNEW=0;
L460:	DISTSQ = 4.0*DELTA*DELTA;

     	for (int K=1; K<=NPT; K++)  {
		SUM=ZERO;
		for (int J=1; J<=N; J++)	SUM += pow2(XPT(J,K)-XOPT[J]);
		if (SUM > DISTSQ)  {
		    KNEW=K;
		    DISTSQ=SUM;
		}
	 }

	//  If KNEW is positive, then set DSTEP, and branch back for the next
	//  iteration, which will generate a "model step".

     	if (KNEW > 0)  {
     	    DSTEP = max (min(TENTH*sqrt(DISTSQ),HALF*DELTA),RHO);
     	    DSQ=DSTEP*DSTEP;
     	    goto shift_xbase_120;
     	}

     	if (RATIO > ZERO)		goto gen_tr_100;
     	if (max (DELTA,DNORM) > RHO)	goto gen_tr_100;

	//  The calculations with the current value of RHO are complete. Pick the
	//  next values of RHO and DELTA.

new_rho_490:	if (RHO >  RHOEND)  {
     		DELTA=HALF*RHO;
     		RATIO=RHO/RHOEND;

     		if (RATIO <= 16.0)
     		    RHO = RHOEND;
     		else if (RATIO <= 250.0) 
     		    RHO = sqrt(RATIO)*RHOEND;
     		else
     		    RHO = TENTH*RHO;

     		DELTA = max (DELTA,RHO);

     		if (IPRINT >= 2) {
     			if (IPRINT >= 3) cout << "	";
			FMT("\n	New RHO =%g   Number of function values =%d")   %RHO  %NF;
  			FMT("\n	Least value of F =%.10g  The corresponding X is:")  %FOPT;
  			for (int i=XBASE.ibegin(); i<XBASE.iend(); i++)  FMT("%.10g")   %(XBASE[i]+XOPT[i]);
			cout << endl;
     		}

     		goto begin_iter_90;
     	}

	//  Return from the calculation, after another Newton-Raphson step, if
	//  it is too short to have been tried before.

     	if (KNEW == -1) goto xnew_290;

 exit_530:
 	if (FOPT <= F) {
     	    for (int I=1; I<=N; I++)	X[I]=XBASE[I]+XOPT[I];
     	    F=FOPT;
     	}

     	if (IPRINT >= 1)  {
  		FMT("\n	At the return from NEWUOA:  Number of function values =%d")  %NF;
		FMT("\n	Least value of F =%.15g  \n	The corresponding X is: %.15g")   %F  %X;
		cout << endl;
     	}

}; // newuoa


int main() {

	
	/*
	for (int N=2;  N <= 8;  N+=2)  {
                int NPT=2*N+1;

		for (int I=1; I<=N; I++) X[I]=I/double(N+1);

                double RHOBEG=0.2e0 * X[1];

                FMT("Results with N =%d and NPT =%d")  % N  % NPT;

              //CALL NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
              newuoa (N, NPT, (double*)&X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}*/
	int IPRINT=2;
	int MAXFUN=5000;
	double RHOEND=1.0e-6;
	{
		const int N=2;
		const int NPT=2*N+1;
		array<double,N, 1> X;
		for (int I=1; I<=N; I++) X[I]=I/double(N+1);
		double RHOBEG=0.2 * X[1];
		FMT("Results with N =%d and NPT =%d")  % N  % NPT;
		newuoa<N, NPT>(X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}

	{
		const int N=4;
		const int NPT=2*N+1;
		array<double,N, 1> X;
		for (int I=1; I<=N; I++) X[I]=I/double(N+1);
		double RHOBEG=0.2 * X[1];
		FMT("Results with N =%d and NPT =%d")  % N  % NPT;
		newuoa<N, NPT>(X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}

	{
		const int N=6;
		const int NPT=2*N+1;
		array<double,N, 1> X;
		for (int I=1; I<=N; I++) X[I]=I/double(N+1);
		double RHOBEG=0.2 * X[1];
		FMT("Results with N =%d and NPT =%d")  % N  % NPT;
		newuoa<N, NPT>(X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}

	{
		const int N=8;
		const int NPT=2*N+1;
		array<double,N, 1> X;
		for (int I=1; I<=N; I++) X[I]=I/double(N+1);
		double RHOBEG=0.2 * X[1];
		FMT("Results with N =%d and NPT =%d")  % N  % NPT;
		newuoa<N, NPT>(X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}
return 0 ;}

