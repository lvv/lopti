/*
    NLEQ_RES - RESidual-based Damped Newton algorithm
    
 *  Written by        L. Weimann 
 *  Purpose           Solution of systems of nonlinear equations
 *  Method            RESidual-based Damped Newton algorithm
                      (see reference below)
 *  Category          F2a. - Systems of nonlinear equations
 *  Keywords          Nonlinear equations, Newton methods
 *  Version           1.1.1
 *  Revision          May 2006
 *  Latest Change     June 2006
 *  Library           NewtonLib
 *  Code              C, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 *  Copyright     (c) Konrad-Zuse-Zentrum fuer
                      Informationstechnik Berlin (ZIB)
                      Takustrasse 7, D-14195 Berlin-Dahlem
                      phone : + 49/30/84185-0
                      fax   : + 49/30/84185-125
 *  Contact           Lutz Weimann
                      ZIB, Division Scientific Computing, 
                           Department Numerical Analysis and Modelling
                      phone : + 49/30/84185-185
                      fax   : + 49/30/84185-107
                      e-mail: weimann@zib.de
 
 *    References:
 
      /1/ P. Deuflhard:
          Newton Methods for Nonlinear Problems. -
          Affine Invariance and Adaptive Algorithms.
          Series Computational Mathematics 35, Springer (2004)
 
   ---------------------------------------------------------------
 
 * Licence
     You may use or modify this code for your own non commercial
     purposes for an unlimited time. 
     In any case you should not deliver this code without a special 
     permission of ZIB.
     In case you intend to use the code commercially, we oblige you
     to sign an according licence agreement with ZIB.
 
 * Warranty 
     This code has been tested up to a certain level. Defects and
     weaknesses, which may be included in the code, do not establish
     any warranties by ZIB. ZIB does not take over any liabilities
     which may follow from acquisition or application of this code.
 
 * Software status 
     This code is under care of ZIB and belongs to ZIB software class 2.
 
      ------------------------------------------------------------
 
 *    Parameters description
      ======================
 
      The calling interface looks as follows:

      extern void nleq_res(struct NLEQ_FUN funs, int n, double *x,
                           struct NLEQ_OPT *opt, 
                           struct NLEQ_INFO *info)

      The structures used within the parameter list are defined
      as follows:
      ---
      struct NLEQ_FUN
      {
        NLEQ_FFUN *fun;
        NLEQ_JFUN *jac;
      };
      
      where the types used within this structure are defined by
      typedef void NLEQ_FFUN(int*,double*,double*,int*);
      and
      typedef void NLEQ_JFUN(int*,int*,double*,double*,int*);
      ---
      struct NLEQ_OPT
      {
         double rtol;
         int maxiter;
         LOGICAL restricted, nleqcalled; 
         PRINT_LEVEL errorlevel, monitorlevel, datalevel;
         FILE *errorfile, *monitorfile, *datafile;
         PROBLEM_TYPE nonlin;
         double *scale;
      };
      
      where the types used within this structure are defined by
      typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL;
      typedef enum {False=0, True=1} LOGICAL ;
      typedef enum {Mildly_Nonlinear=2, Highly_Nonlinear=3,
                    Extremely_Nonlinear=4} PROBLEM_TYPE ;
      ---
      struct NLEQ_INFO
      {
         double precision, normdx;
         double *fx;
         int iter, rcode, subcode, nofunevals, nojacevals;
      };
      ---
      
      A detailed description of the parameters follows: 
      
      struct NLEQ_FUN funs :
      
      The field funs.fun must contain a pointer to the problem function fun -
      The required parameters interface of fun is described in detail below
      
      The field funs.jac must either contain a pointer to the Jacobian function
      jac or a NULL pointer. If a NULL pointer is supplied, then the Jacobian
      will be approximately computed by an internal function of the nleq_res
      package.
      
      int n :
      The number of equations and unknown variables of the nonlinear system.
      
      double *x :
      A pointer to an array of double values of size n.
      The array must contain on input an initial guess of the problems solution,
      which is used as the start-vector of the damped Newton iteration.
      On output, the pointed array contains an approximate solution vector x*,
      which fits the small residual condition
      || funs.fun(x*) || <= opt->tol,
      where ||...|| is a scaled Euclidian norm.
      
      struct NLEQ_OPT *opt:
      A pointer to an options structure. The pointed fields of the structure
      contain input options to nleq_res.
      
      opt->tol is of type double and must contain the residuum threshold
      which funs.fun(x*) must fit for the solution vector x*.
      
      opt->maxiter is of type int and must contain the maximum number of allowed
      iterations. if a zero or negative value is supplied, then the maximum
      iteration count will be set to 50.
      
      opt->nonlin is of type PROBLEM_TYPE and must classify the problem to
      be solved. The following classifications may be used:
      Mildly_Nonlinear: The problem is considered to be mildly nonlinear and
                        nleq_res starts up with dampingfactor=1.
      Highly_Nonlinear: The problem is considered to be highly nonlinear and
                        nleq_res starts up with dampingfactor=1.0e-4.
      Extremely_Nonlinear: The problem is considered to be extremely nonlinear
                        and nleq_res starts up with dampingfactor=1.0e-6.
                        Moreover, opt->restricted is set automatically to True.

      opt->restricted is of type LOGICAL.
      If set to True, then the restricted monotonicity test will be applied for
      determination whether the next iterate (and the associate damping factor
      lambda) will be accepted. This means, with
      theta = ||F(x(k+1))|| / ||F(x(k))||, the condition
      theta <= 1.0 - lambda/4 must be fit
      If set to False, then the standard monotonicity test will be applied, i.e.
      the following condition must be fit:
      theta < 1.0.
      
      opt->nleqcalled is of type LOGICAL and only used internally. This field
      should always be set to False.
      
      opt->errorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no error message will be printed if an error condition occurs.
      If it is set to level Minimum or any higher level, then error messages
      will be printed, if appropriate.
      
      opt->monitorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no monitor output will be printed.
      If it is set to level Minimum, a few infomation will be printed.
      If set to level Verbose, then some infomation about each Newton iteration
      step, fitting into a single line, will be printed. The higher level Debug
      is reserved for future additional information output.
      
      opt->datalevel is of type PRINT_LEVEL. If it is set to level None,
      then no data output will be printed.
      If it is set to level Minimum, then the values of the initial iteration
      vector x and the final vector x will be printed.
      If set to level Verbose, then the iteration vector x will be printed for
      each Newton step. The higher level Debug is reserved for future additional
      information output.
      
      opt->errorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->errorfile will be set to stdout. The error 
      messages will be printed to opt->errorfile.
      
      opt->monitorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->monitorfile will be set to stdout. The monitor 
      output will be printed to opt->monitorfile.
      
      opt->datafile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, a file named "nleq_res.data" will be opened by a
      fopen call and opt->datafile will be set to the filepointer which the
      fopen returns. The data output will be printed to opt->datafile.
      
      opt->iterfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the iteration vector will be written out to the associated file, for
      each Newton iteration step. If opt->iterfile is set to NULL, no such 
      data will be written out.
      
      opt->resfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the residuum vector will be written out to the associated file, for
      each Newton iteration step. If opt->resfile is set to NULL, no such 
      data will be written out.
      
      opt->miscfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number, an
      identification number of the calling code (1 for NLEQ_RES), the norm
      of the residuum, the norm of the Newton correction, a zero value 
      as a dummy placeholder value, the accepted damping factor, and another
      zero value as a dummy placeholder value will be written out, for
      each Newton iteration step. If opt->miscfile is set to NULL, no such 
      data will be written out. For additional information on the file output,
      refer to the description of this option in the QNRES documentation.
     
      Note: The output to the files opt->iterfile, opt->resfile and
            opt->miscfile is written as a single for each iteration step.
            Such, the data in these files are suitable as input to the
            graphics utility GNUPLOT.
      
      opt->scale is of type pointer to a double array of size n. 
      This array must, if present, contain positive scaling values, which are
      used in computations of scaled norms and Jacobian scaling, as follows:
      || f || = squareroot(sum(1 to n) ( f_i/scale_i )^2)
      The pointer may be initialized with a NULL pointer. In this case, all
      scaling values are internally set to 1.
      
      opt->scaleopt is of type enum{...}. This option is only meaningful, if
      the user has not supplied a scaling vector via opt->scale.
      In this case, if opt->scaleopt is set to StandardScale, all scaling
      vector components are set to 1. If it is set to StartValueScale,
      then the setting is fscale_i = max(1,abs(f0_i)) for i=0,...,n-1,
      where f0=(f0_0,...,f0_(n-1))=problem_function(x0), x0=start vector.
      
      struct NLEQ_INFO *info:
      A pointer to an info structure. The pointed fields of the structure
      are set output info of nleq_res.
      
      info->precision is of type double and is set to the achieved scaled norm
      of the residuum at the final iterate.
      
      info->normdx is of type double and is set to the unscaled norm of the
      last Newton correction.
      
      info->fx is a pointer to a double array of size n, which contains the
      final residuum vector.
      
      info->iter is set to number of Newton iteration steps done.
      
      info->nofunevals is set to the number of done calls to the problem
      function funs.fun.
      
      info->nojacevals is set to the number of done calls to the Jacobian
      function funs.jac.
      
      info->rcode is set to the return-code of nleq_res. A return-code 0
      means that nleq_res has terminated sucessfully. For the meaning of a
      nonzero return-code, see the error messages list below.
      
      info->subcode is set for certain failure conditions to the error code
      which has been returned by a routine called from nleq_res.


      Parameter definitions of the required problem routine funs.fun
      and the optional Jacobian routine funs.jac
      --------------------------------------------------------------
      
      void fun(int *n, double *x, double *f, int *fail);
        int    *n     input  Number of vector components.
        double *x     input  Vector of unknowns, of size *n .
        double *f     output Vector of function values.
        int    *fail  output fun evaluation-failure indicator.
                      On input:  undefined.
                      On output: Indicates failure of fun evaluation, 
                                 if having a value <= 2.
                      If <0 or >2: nleq_res will be terminated with
                                   error code = 82, and *fail will be stored to
                                   info->subcode.
                      If =1: A new trial Newton iterate will
                             computed, with the damping factor
                             reduced to it's half.
                      If =2: A new trial Newton iterate will computed, with the
                             damping factor reduced by a reduction factor, which
                             must be output through f[0] by fun, and it's value
                             must be >0 and <1.
                      Note, that if IFAIL = 1 or 2, additional conditions 
                      concerning the damping factor, e.g. the minimum damping
                      factor may also influence the value of the reduced 
                      damping factor.
 
      void jac(int *n, int *ldjac, double *x, double *dfdx, int *fail);
                        Ext    Jacobian matrix subroutine
        int    *n     input  Number of vector components.
        int    *ldjac input  Leading dimension of Jacobian array, i.e.
                             the total row length for C-style two-dimensional
                             arrays, or the total column length for 
                             Fortran-style two-dimensional arrays.
                             See Note below!
        double *x     input  Vector of unknowns, of size *n .
        double *dfdx  output dfdx[i][k]: partial derivative of i-th component
                             of output parameter *f from fun with respect 
                             to x[k].
        int    *fail  output fun evaluation-failure indicator.
                      On input:  undefined.
                      On output: Indicates failure of jac evaluation
                      and causes termination of nleq_res, f set to a nonzero
                      value on output.
                                 
      Note: The calling interfaces of the user routines fun and jac has
            been designed to be compatible with routines programmed for
            use with the Fortran codes NLEQ1 and NLEQ2. However, note
            that the Fortran matrix storage mode is columnwise while
            the C matrix storage mode is rowwise. If you intend to link 
            a Jacobian routine, which has been programmed in Fortran for
            use with NLEQ1 or NLEQ2, you must either transpose the Jacobian,
            or you must compile the nleq_res package for use with Fortran
            matrix storage mode, by setting the C preprocessor flag FMAT,
            i.e. setting the gcc compiler option -DFMAT, when compiling the
            file jacobian_and_linalg.c .


      The following error conditions may occur: (returned via info->rcode)
      --------------------------------------------------------------------
      
      -999 routine nleq_fwalloc failed to allocate double memory via malloc.
      -998 routine nleq_iwalloc failed to allocate int memory via malloc.
      -997 routine nleq_pfwalloc failed to allocate double pointer memory
           via malloc.
      -995 Internal i/o control block could not be allocated via malloc.
      -994 Internally used data structure could not be allocated via malloc.
      -989 Default data-output file could not be opened via fopen call.
       -99 NULL pointer obtained from funs.fun field - the problem function
           must be defined!
         1 Singular Jacobian matrix (detected by routine nleq_linfact),
           nleq_res cannot proceed the iteration.
         2 Maximum number of Newton iteration (as set by opt->maxiter) exceeded.
         3 No convergence of Newton iteration, damping factor became too small.
        20 Nonpositive input for dimensional parameter n.
        21 Nonpositive value for opt->tol supplied.
        22 Negative scaling value for some component of vector opt->scale
           supplied.
        80 nleq_linfact returned with an error other than singular Jacobian.
           Check info->subcode for the nleq_linfact failure code.
        81 nleq_linsol returned with an error.
           Check info->subcode for the nleq_linsol failure code.
        82 The user defined problem function funs.fun returned a nonzero code
           other than 1 or 2. 
           Check info->subcode for the user-function failure code.
        83 The user defined Jacobian function funs.jac returned a nonzero code.
           Check info->subcode for the Jacobian-function failure code.
         
      Summary of changes:
      -------------------
      
      Version  Date        Changes
      1.1.1    2006/06/06  Missing int return code in function 
                           nleqres_initscale, bug fixed.
      1.1      2006/06/02  Added the output data files iterfile, resfile and
                           miscfile, where optional data output is provided,
                           a single line, starting with the iteration number,
                           for each iteration step. 
      1.0      2006/05/30  Initial release.
      
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nleq.h"

int nleqres_initscale(int n, double **scale, struct NLEQ_OPT *opt);
void nleqres_monitor(int k, int n, double normdx, double normf,
                     double lambda);

#define THETA_MAX 0.25
#define MAX_ITER_DEFAULT 50
#define LAMBDA_START_DEFAULT 1.0e-2
#define LAMBDA_START_EXTREMELY_DEFAULT 1.0e-4
#define LAMBDA_MIN_DEFAULT 1.0e-4
#define LAMBDA_MIN_EXTREMELY_DEFAULT 1.0e-8
#define NONLIN opt->nonlin

extern struct NLEQ_IO *nleq_ioctl;

extern void nleq_res(struct NLEQ_FUN fun, int n, double *x,
                  struct NLEQ_OPT *opt, struct NLEQ_INFO *info)
{  int inc_one=1;
   double one=1.0, minus_one=-1.0;
   double ftol=opt->tol;
   double lambda, lambda_new, mue, normfk, normfkm1, normfkp1, 
          reduction_factor, normdx, theta, s;
   double lambda_min;
   int i, j, k=0, fail=0, nrhs=1, ldjac,
            max_iter=opt->maxiter;
   LOGICAL     qnres_iter=False,
               saved_nleqcalled=opt->nleqcalled,
               restricted=opt->restricted,
               io_allocated=False,
               reducted;
   PRINT_LEVEL error_level;
   double *dx, *fxk, *fxkp1, *fscale, *xkp1, *w;
   JACOBIAN *jac;
   int scale_allocated=0, nfcn=0, njac=0;
   NLEQ_FFUN *f = fun.fun;
   NLEQ_JFUN *jc = fun.jac;
   struct NLEQ_DATA *data=malloc(sizeof(struct NLEQ_DATA));
   if (!nleq_ioctl) nleq_ioctl=malloc(sizeof(struct NLEQ_IO));
   if (!nleq_ioctl) 
     { fprintf(stderr,"\n could not allocate output controlblock\n");
       RCODE=-995; return; }
   else
     io_allocated = True;
   if (!data)
     { fprintf(stderr,"\n could not allocate struct data\n");
       RCODE=-994; return; };
   data->codeid    = NLEQ_RES;
   data->normdxbar = 0.0;
   data->theta     = 0.0;
   data->mode      = Initial;
   ERRORLEVEL   = opt->errorlevel;
   MONITORLEVEL = opt->monitorlevel;
   DATALEVEL    = opt->datalevel;
   error_level  = opt->errorlevel;
   ERROR    = opt->errorfile;
   MONITOR  = opt->monitorfile;
   DATA     = opt->datafile;
   FITER    = opt->iterfile;
   FRES     = opt->resfile;
   FMISC    = opt->miscfile;
   if ( !ERROR && ERRORLEVEL>0 )     ERROR   = stdout;
   if ( !MONITOR && MONITORLEVEL>0 ) MONITOR = stdout;
   if ( !DATA && DATALEVEL>0 )
     { DATA=fopen("nleq_res.data","w");
       if (!DATA && ERRORLEVEL>0)
         { fprintf(ERROR,"\n fopen of file nleq_res.data failed\n");
           RCODE=-989; return;
         };
     };
   opt->errorfile   = ERROR;
   opt->monitorfile = MONITOR;
   opt->datafile    = DATA;
   if ( MONITORLEVEL > 0 ) fprintf(MONITOR,"\n NLEQ_RES - Version 1.1\n");
   RCODE = nleq_parcheck_and_print(n,opt,fun,1);
   if ( RCODE !=0 ) 
     { if (io_allocated) {free(nleq_ioctl); nleq_ioctl=NULL;};
       if (data) free(data);
       return;
     };
   opt->nleqcalled = True;
   if ( max_iter <= 0 ) max_iter = MAX_ITER_DEFAULT;
   if      ( NONLIN==Mildly_Nonlinear ) 
     { lambda = 1.0; lambda_min = LAMBDA_MIN_DEFAULT; }
   else if ( NONLIN==Highly_Nonlinear ) 
     { lambda = LAMBDA_START_DEFAULT; lambda_min = LAMBDA_MIN_DEFAULT; }
   else if ( NONLIN==Extremely_Nonlinear ) 
     { lambda = LAMBDA_START_EXTREMELY_DEFAULT;
       lambda_min = LAMBDA_MIN_EXTREMELY_DEFAULT;
       restricted = True;
     } ;
   RCODE = nleq_jacalloc(n,&jac,&ldjac,opt);  if ( RCODE !=0 ) return;
   RCODE = nleq_fwalloc(n,&dx,"dx");          if ( RCODE !=0 ) return;
   RCODE = nleq_fwalloc(n,&xkp1,"xkp1");      if ( RCODE !=0 ) return;
   RCODE = nleq_fwalloc(n,&fxk,"fxk");        if ( RCODE !=0 ) return;
   RCODE = nleq_fwalloc(n,&fxkp1,"fxkp1");    if ( RCODE !=0 ) return;
   RCODE = nleq_fwalloc(n,&w,"w");            if ( RCODE !=0 ) return;
   data->fx = fxk;
   data->dx = dx;
   f(&n,x,fxk,&fail);  nfcn++;
   if (fail != 0) { RCODE=82; goto errorexit ;};
   fscale = opt->scale;
   if ( fscale == NULL ) 
     {RCODE=nleqres_initscale(n,&fscale,opt);
      if (RCODE !=0) goto errorexit;
      scale_allocated = 1; opt->scale = fscale; };
   if ( opt->scaleopt == StartValueScale && scale_allocated == 1 ) 
     for (i=0;i<n;i++) fscale[i] = MAX(fabs(fxk[i]),1.0);
   if ( MONITORLEVEL > 1 ) 
      fprintf(MONITOR,"\n iter      norm(dx)  norm_scl(fk)    lambda \n\n");
   normfk = nleq_scaled_norm2(n,fxk,fscale);
   normdx = 0.0;
   RCODE = 2;

   do 
     { 
       if ( normfk <= ftol )
          {RCODE=0; break;}; /* stop, x contains the solution! */
       if (*jc != NULL) 
         { jc(&n,&ldjac,x,jac,&fail);  njac++;
           if (fail !=0) {RCODE = 83; goto errorexit;};
         }
       else 
         { fail=nleq_numjac(f,n,x,fxk,NULL,jac,&nfcn,opt);
           if (fail !=0) { RCODE = 82; goto errorexit; };
         };
       /* scale jacobian and change sign of it */
       nleq_jacrow_scale(n,jac,fscale,opt);
       fail = nleq_linfact(n,jac,opt);  /* compute LU-factorization of Jacobian */
       if ( fail < 0 )      { RCODE=80; goto errorexit; }
       else if ( fail > 0 ) { RCODE=1;  goto errorexit; };
       /* compute newton correction */
       nleq_scale(n,fxk,dx,fscale);
       fail = nleq_linsol(n,dx,opt);
       if ( fail != 0 ) { RCODE=81; goto errorexit; };
       normdx = nleq_norm2(n,dx);
       if ( k>0 )
         { mue = (normfkm1/normfk)*mue; lambda = MIN(1.0,mue); };
       if ( MONITORLEVEL > 1 ) 
          normfkp1 = normfk; /* set normfkp1 for print only */
       reducted = False;
       checkregularity:
       if ( MONITORLEVEL > 1 && fail==0 ) 
         nleqres_monitor(k,n,normdx,normfkp1,lambda);
       if ( lambda < lambda_min ) 
         { RCODE=3; break; };  /* stop, convergence failure! */
       for (i=0;i<n;i++) xkp1[i]=x[i]+lambda*dx[i]; /* new trial iterate */
       f(&n,xkp1,fxkp1,&fail);  nfcn++;
       if ( fail<0 || fail>2 ) { RCODE=82; goto errorexit; }
       else if ( fail==1 || fail== 2 )
         { if ( fail==1 ) reduction_factor = 0.5;
           else           reduction_factor = fxkp1[0];
           if ( reduction_factor <= 0.0 || reduction_factor >= 1.0 )
             { RCODE=82; goto errorexit; };
           if ( MONITORLEVEL>1 )
             fprintf(MONITOR," %4i  FUN could not be evaluated  %7f\n",
                             k,lambda);
           if ( lambda > lambda_min )
             lambda = MAX(lambda*reduction_factor,lambda_min);
           else 
             lambda = lambda*reduction_factor;
           reducted = True;
           goto checkregularity;  
         };
       normfkp1 = nleq_scaled_norm2(n,fxkp1,fscale);
       theta = normfkp1/normfk;
       s = 1.0-lambda;
       for (i=0;i<n;i++) w[i] = fxkp1[i]-s*fxk[i];
       mue = (0.5*normfk*lambda*lambda) / nleq_scaled_norm2(n,w,fscale);
       if ( ( !restricted && theta >= 1.0 ) || 
            (  restricted && theta > 1.0-lambda/4.0) )
         { lambda_new = MIN(mue,0.5*lambda);
           if ( lambda <= lambda_min ) lambda = lambda_new;
           else                        lambda = MAX(lambda_new,lambda_min);
           reducted = True;
           goto checkregularity; };
       lambda_new = MIN(1.0,mue);
       if ( lambda==1.0 && lambda_new==1.0 && theta < THETA_MAX )
         qnres_iter = True;
       else
         { if( lambda_new >= 4.0*lambda && !reducted ) 
             { lambda=lambda_new; goto checkregularity; };
         };
       data->normf = normfk;         data->normdx = normdx;
       data->lambda = lambda;
       nleq_dataout(k,n,x,data);
       data->mode = Intermediate;
       for (i=0;i<n;i++) x[i]=xkp1[i];  /* accept new iterate */
       /* next step */
       k++;  normfkm1 = normfk;  normfk   = normfkp1;
       /* start qnres only, if residuum not yet small enough  */ 
       qnres_iter = qnres_iter && normfk > ftol;
       /* perform qnres steps, if choosen */
       if (qnres_iter)
         { info->fx=fxk;
           info->iter=k;
           info->normdx=normdx;
           opt->maxiter=max_iter-k+1;
           opt->errorlevel = 0;
           qnres(fun,n,x,opt,info);
           nfcn += info->nofunevals;
           k = info->iter;
           normfkp1=info->precision;
           opt->errorlevel = error_level;
           ERRORLEVEL      = error_level;
           /* if QNRES failed, try to continue NLEQ_RES */
           if ( RCODE != 0 ) {RCODE = 2; qnres_iter=False; };
         }
       else
         for (i=0;i<n;i++) fxk[i]=fxkp1[i];
     }
   while ( k <= max_iter && RCODE == 2 );

   if ( !qnres_iter )
     { if ( MONITORLEVEL > 1 ) 
         nleqres_monitor(k,n,normdx,normfk,lambda);
       data->normf = normfk;         data->normdx = normdx;
       data->lambda = lambda;
       data->mode = ( RCODE==0 ? Solution : Final );
       nleq_dataout(k,n,x,data);
     };
     
   errorexit:
   if ( ERRORLEVEL > 0 && RCODE != 0 )
     {
       switch ( RCODE )
        {
         case     1:
           fprintf(ERROR,"\n Error return from nleq_linfact: Singular Jacobian\n");
           break;
         case     2:
           fprintf(ERROR,"\n Error - Maximum allowed number of iterations exceeded\n");
           break;
         case     3:
           fprintf(ERROR,"\n Error - no convergence, damping factor became too small\n");
           break;
         case    80:
           fprintf(ERROR,"\n Error return from nleq_linfact: fail=%i\n",fail);
           break;
         case    81:
           fprintf(ERROR,"\n Error return from nleq_linsol: fail=%i\n",fail);
           break;
         case    82:
           fprintf(ERROR,"\n Error return from problem function: fail=%i\n",
                         fail);
           break;
         case    83:
           fprintf(ERROR,"\n Error return from Jacobian function: fail=%i\n",
                         fail);
           break;
         default   :
           fprintf(ERROR,"\n Error, code=%i,  subcode=%i\n",RCODE,fail);
        };
     };
   info->subcode = fail;
   if (io_allocated) {free(nleq_ioctl); nleq_ioctl=NULL;};
   free(data);
   free(dx); free(xkp1); if (qnres_iter) free(fxk);
   free(fxkp1); free(w);
   if (scale_allocated) { free(fscale);  opt->scale = NULL; };
   nleq_linfree();
   if (!qnres_iter) 
     { info->precision = normfk;
       info->normdx    = normdx;
     };
   /* restore original values */
   opt->maxiter    = max_iter;
   opt->nleqcalled = saved_nleqcalled;
   info->iter       = k; 
   info->nofunevals = nfcn;
   info->nojacevals = njac;
   if (!qnres_iter) info->fx         = fxk;
}

int nleqres_initscale(int n, double **scale, struct NLEQ_OPT *opt)
{  int i, rcode;
   rcode = nleq_fwalloc(n,scale,"scale");  if ( rcode != 0 ) return rcode;
   for (i=0;i<n;i++) (*scale)[i]=1.0;  return 0;
}

void nleqres_monitor(int k, int n, double normdx, double normf,
                     double lambda)
{  fprintf(MONITOR," %4i  %12e  %12e  %7f\n",k,normdx,normf,lambda);
}
