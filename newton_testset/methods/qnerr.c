/*
    QNERR - ERRor-based Quasi-Newton algorithm
    
 *  Written by        L. Weimann 
 *  Purpose           Solution of systems of nonlinear equations
 *  Method            ERRor-based Quasi-Newton algorithm
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

      extern void qnerr(struct NLEQ_FUN funs, int n, double *x,
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
      will be approximately computed by an internal function of the qnerr
      package.
      
      int n :
      The number of equations and unknown variables of the nonlinear system.
      
      double *x :
      A pointer to an array of double values of size n.
      The array must contain on input an initial guess of the problems solution,
      which is used as the start-vector of the Quasi-Newton iteration.
      On output, the pointed array contains an approximate solution vector x^k,
      which fits the error condition
      || x^k - x* || <= opt->tol,
      where x* denotes the exact solution, and ||...|| is a scaled 
      Euclidian norm.
      
      struct NLEQ_OPT *opt:
      A pointer to an options structure. The pointed fields of the structure
      contain input options to qnerr.
      
      opt->tol is of type double and must contain the error tolerance which
      the final approximate solution x^k must fit.
      
      opt->maxiter is of type int and must contain the maximum number of allowed
      iterations. if a zero or negative value is supplied, then the maximum
      iteration count will be set to 50.
      
      opt->nleqcalled is of type LOGICAL and only used internally. This field
      must always be set to False.
      
      opt->errorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no error message will be printed if an error condition occurs.
      If it is set to level Minimum or any higher level, then error messages
      will be printed, if appropriate.
      
      opt->monitorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no monitor output will be printed.
      If it is set to level Minimum, a few infomation will be printed.
      If set to level Verbose, then some infomation about each Quasi-Newton
      iteration step, fitting into a single line, will be printed. The higher
      level Debug is reserved for future additional information output.
      
      opt->datalevel is of type PRINT_LEVEL. If it is set to level None,
      then no data output will be printed.
      If it is set to level Minimum, then the values of the initial iteration
      vector x and the final vector x will be printed.
      If set to level Verbose, then the iteration vector x will be printed for
      each Quasi-Newton step. The higher level Debug is reserved for future
      additional information output.
      
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
      If it is set to NULL, a file named "qnerr.data" will be opened by a
      fopen call and opt->datafile will be set to the filepointer which the
      fopen returns. The data output will be printed to opt->datafile.
      
      opt->iterfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the iteration vector will be written out to the associated file, for
      each Quasi-Newton iteration step. If opt->iterfile is set to NULL,
      no such data will be written out.
      
      opt->resfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the residuum vector will be written out to the associated file, for
      each Quasi-Newton iteration step. If opt->resfile is set to NULL,
      no such data will be written out.
      
      opt->miscfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number, an
      identification number of the calling code (2 for QNERR), the norm
      of the residuum, the norm of the Quasi-Newton correction, a zero value 
      as a dummy placeholder value, a 1.0 as damping factor value, and
      the control parameter theta will be written out, for each Quasi-Newton
      iteration step. If opt->miscfile is set to NULL, no such data will be
      written out.
     
      Note: The output to the files opt->iterfile, opt->resfile and
            opt->miscfile is written as a single for each iteration step.
            Such, the data in these files are suitable as input to the
            graphics utility GNUPLOT.
      
      opt->scale is of type pointer to a double array of size n. 
      This array must, if present, contain positive scaling values, which are
      used in computations of scaled norms and Jacobian scaling, as follows:
      || x || = squareroot(sum(1 to n) ( x_i/scale_i )^2)
      The pointer may be initialized with a NULL pointer. In this case, all
      scaling values are internally set to 1.
      
      struct NLEQ_INFO *info:
      A pointer to an info structure. The pointed fields of the structure
      are set output info of qnerr.
      
      info->precision is of type double and is set to the achieved scaled norm
      of the error of the final iterate.
      
      info->normdx is of type double and is set to the unscaled norm of the
      last Quasi-Newton correction.
      
      info->fx is a pointer to a double array of size n, which contains the
      final residuum vector.
      
      info->iter is set to number of Quasi-Newton iteration steps done.
      
      info->nofunevals is set to the number of done calls to the problem
      function funs.fun.
      
      info->nojacevals is set to the number of done calls to the Jacobian
      function funs.jac.
      
      info->rcode is set to the return-code of qnerr. A return-code 0
      means that qnerr has terminated sucessfully. For the meaning of a
      nonzero return-code, see the error messages list below.
      
      info->subcode is set for certain failure conditions to the error code
      which has been returned by a routine called from qnerr.


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
                                 if set to a nonzero value. In this case,
                                 qnerr will be terminated with error code=82,
                                 and *fail will be stored to info->subcode.
 
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
                      On output: Indicates failure of jac evaluation,
                                 if set to a nonzero value. In this case,
                                 qnerr will be terminated with error code=83,
                                 and *fail will be stored to info->subcode.
                                 
      Note: The calling interfaces of the user routines fun and jac has
            been designed to be compatible with routines programmed for
            use with the Fortran codes NLEQ1 and NLEQ2. However, note
            that the Fortran matrix storage mode is columnwise while
            the C matrix storage mode is rowwise. If you intend to link 
            a Jacobian routine, which has been programmed in Fortran for
            use with NLEQ1 or NLEQ2, you must either transpose the Jacobian,
            or you must compile the qnerr package for use with Fortran
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
           qnerr cannot proceed the iteration.
         2 Maximum number of Quasi-Newton iteration (as set by opt->maxiter)
           exceeded.
         3 No convergence of Quasi-Newton iteration, theta became too large.
         6 Ill conditioned update (kappa became too large)
        80 nleq_linfact returned with an error other than singular Jacobian.
           Check info->subcode for the nleq_linfact failure code.
        81 nleq_linsol returned with an error.
           Check info->subcode for the nleq_linsol failure code.
        82 The user defined problem function funs.fun returned a nonzero code.
           Check info->subcode for the user-function failure code.
        83 The user defined Jacobian function funs.jac returned a nonzero code.
           Check info->subcode for the Jacobian-function failure code.
         
      Summary of changes:
      -------------------
      
      Version  Date        Changes
      1.1.1    2006/06/06  Missing int return code in function 
                           qne_initscale, bug fixed.
      1.1      2006/06/02  Added the output data files iterfile, resfile and
                           miscfile, where optional data output is provided,
                           a single line, starting with the iteration number,
                           for each iteration step. 
      1.0      2006/05/30  Initial release.
      
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nleq.h"

int qne_initscale(int n, double **scale, double *x, double tol);
void qne_monitor(int k, int n, double normdx, double normf,
                 char qnerr_id[], char fail_reason[]);

#define THETA_MAX 0.5
#define KAPPA_MAX 1.0e5
#define MAX_ITER_DEFAULT 50

extern struct NLEQ_IO *nleq_ioctl;

extern void qnerr(struct NLEQ_FUN fun, int n, double *x,
                  struct NLEQ_OPT *opt, struct NLEQ_INFO *info)
{  int inc_one=1;
   double one=1.0, minus_one=-1.0;
   double thetak, xtol=opt->tol;
   double s, normfk, normdx, alpha_bar, alpha_kp1;
   int i, k=0, k1, kprint, fail=0, ldjac, 
       max_iter=opt->maxiter;
   LOGICAL nleqcalled = opt->nleqcalled, skipstep, 
           io_allocated = False;
   double *dxk, *fxk, *sigma, *v;
   double *xscale;
   double **dx;
   JACOBIAN *jac;
   int scale_allocated=0, nfcn=0, njac=0;
   char  qnerr_id[9]="        \0", fail_reason[7]="      \0";
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
   data->codeid    = QNERR;
   data->normdxbar = 0.0;
   data->lambda    = 1.0;
   ERRORLEVEL   = opt->errorlevel;
   MONITORLEVEL = opt->monitorlevel;
   DATALEVEL    = opt->datalevel;
   ERROR    = opt->errorfile;
   MONITOR  = opt->monitorfile;
   DATA     = opt->datafile;
   FITER    = opt->iterfile;
   FRES     = opt->resfile;
   FMISC    = opt->miscfile;
   if ( !ERROR && ERRORLEVEL>0 )     ERROR   = stdout;
   if ( !MONITOR && MONITORLEVEL>0 ) MONITOR = stdout;
   if ( !DATA && DATALEVEL>0 )
     { DATA=fopen("qnerr.data","w");
       if (!DATA && ERRORLEVEL>0)
         { fprintf(ERROR,"\n fopen of file qnerr.data failed\n");
           RCODE=-989; return;
         };
     };
   opt->errorfile   = ERROR;
   opt->monitorfile = MONITOR;
   opt->datafile    = DATA;
   if ( !nleqcalled )
     { if ( MONITORLEVEL > 0 )
         fprintf(MONITOR,"\n QNERR - Version 1.1\n");
       RCODE = nleq_parcheck_and_print(n,opt,fun,2);
       if ( RCODE !=0 ) 
         { if (io_allocated) {free(nleq_ioctl); nleq_ioctl=NULL;};
           if (data) free(data);
           return;
         };
     };
   if ( max_iter <= 0 ) max_iter = MAX_ITER_DEFAULT;
   RCODE = nleq_pfwalloc(max_iter+2,&dx,"dx");        if ( RCODE !=0 ) return;
   RCODE = nleq_fwalloc(max_iter+2,&sigma,"sigma");   if ( RCODE !=0 ) return;
   RCODE = nleq_fwalloc(n,&fxk,"fxk");                if ( RCODE !=0 ) return;
   RCODE = nleq_fwalloc(n,&v,"v");                    if ( RCODE !=0 ) return;
   xscale = opt->scale;
   if ( xscale == NULL ) 
     {RCODE=qne_initscale(n,&xscale,x,xtol); if (RCODE !=0) goto errorexit;
      scale_allocated = 1;};
   skipstep=nleqcalled;
   if(!skipstep)
     { RCODE = nleq_jacalloc(n,&jac,&ldjac,opt);      if ( RCODE !=0 ) return;
       f(&n,x,fxk,&fail);  nfcn++;
       if (fail != 0) {RCODE=82; goto errorexit ;};
       normfk = nleq_norm2(n,fxk);
       if (*jc != NULL) 
         { jc(&n,&ldjac,x,jac,&fail);  njac++;
           if (fail !=0) {RCODE = 83; goto errorexit;};
         }
       else 
         { fail=nleq_numjac(f,n,x,fxk,xscale,jac,&nfcn,opt);
           if (fail !=0) {RCODE = 82; goto errorexit;};
         };
       /* column scale jacobian and change sign of it */
       nleq_jaccolumn_scale(n,jac,xscale,opt);
       fail = nleq_linfact(n,jac,opt);  /* compute LU-factorization of Jacobian */
       if ( fail < 0 ) { RCODE=80; goto errorexit; }
       else if ( fail > 0 ) { RCODE=1; goto errorexit; };
       RCODE = nleq_fwalloc(n,&dx[0],"dx[0]");      if ( RCODE !=0 ) return;
       dxk = dx[0];
       for (i=0;i<n;i++) dxk[i]=fxk[i];
       fail = nleq_linsol(n,dx[0],opt);  /* compute first quasi newton correction */
       if ( fail != 0 ) { RCODE=81; goto errorexit; };
       nleq_descale(n,dx[0],dx[0],xscale);
       sigma[0] = nleq_scaled_sprod(n,dx[0],dx[0],xscale);
       normdx = sqrt(sigma[0]);
       data->normf = normfk;  data->normdx = normdx;  data->theta = 0.0; 
       data->dx = dxk;
       data->mode = Initial;
       nleq_dataout(k,n,x,data);
       data->mode = Intermediate;
       if ( MONITORLEVEL > 1 ) 
         { fprintf(MONITOR,"\n iter     norm_scl(dx)      norm(fk)\n\n");
           qne_monitor(k,n,normdx,normfk,qnerr_id,fail_reason);
         };
     }
   else
     { strcpy(qnerr_id,"   QNERR\0");  data->mode = Intermediate;
       dx[0] = info->dx;   s = info->normdx;   sigma[0] = s*s;
     };

   do 
     { 
       if (!skipstep) 
         { daxpy_(&n,&one,dx[k],&inc_one,x,&inc_one); /* x(k+1)=x(k)+dx(k) */
           if ( sigma[k] <= xtol*xtol ) { k++; RCODE=0; break; };
         };
       skipstep = False;
       RCODE = nleq_fwalloc(n,&dx[k+1],"dx[k+1]");    if ( RCODE !=0 ) return;
       RCODE = 2;
       f(&n,x,fxk,&fail);  nfcn++;
       if (fail != 0) {RCODE=82; break;};
       data->fx = fxk;
       normfk = nleq_norm2(n,fxk);
       for (i=0;i<n;i++) v[i]=fxk[i];
       fail = nleq_linsol(n,v,opt);
       if ( fail != 0 ) { RCODE=81; break; };
       nleq_descale(n,v,v,xscale);
       for (i=1;i<=k;i++)
         { alpha_bar = nleq_scaled_sprod(n,v,dx[i-1],xscale) / sigma[i-1];
           /* v=v+alpha_bar*dx(i) */
           daxpy_(&n,&alpha_bar,dx[i],&inc_one,v,&inc_one); 
         };
       alpha_kp1 = nleq_scaled_sprod(n,v,dx[k],xscale) / sigma[k];  
       thetak = sqrt( nleq_scaled_sprod(n,v,v,xscale) / sigma[k] );
       if ( thetak > THETA_MAX ) { k++; RCODE=3; break; };
       s = 1.0-alpha_kp1;
       dxk = dx[k+1];
       for (i=0;i<n;i++) dxk[i] = v[i]/s;
       sigma[k+1] = nleq_scaled_sprod(n,dx[k+1],dx[k+1],xscale);
       normdx = sqrt(sigma[k+1]);
       k++;
       kprint = ( nleqcalled ? (info->iter)+k-1 : k );
       data->normf = normfk;  data->normdx = normdx;  data->theta = thetak;  
       data->dx = dxk;
       nleq_dataout(kprint,n,x,data);
       if ( MONITORLEVEL > 1 ) 
          qne_monitor(kprint,n,normdx,normfk,qnerr_id,fail_reason);
     }
   while ( k <= max_iter && RCODE == 2 );

   kprint = ( nleqcalled ? (info->iter)+k-1 : k );
   if ( RCODE != 2 && RCODE != 0 )
     { if ( MONITORLEVEL > 1 ) 
         { if (nleqcalled && RCODE==3 ) strcpy(fail_reason,"THETA!\0");
           qne_monitor(kprint,n,normdx,normfk,qnerr_id,fail_reason);
         };
     };
   data->normf = normfk;  data->normdx = normdx;  data->theta = thetak;
   data->dx = dxk;
   data->mode = ( RCODE==0 ? Solution : Final );
   if ( RCODE==0 || !nleqcalled ) nleq_dataout(kprint,n,x,data);
   
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
           fprintf(ERROR,"\n Error - QNERR failed - theta = %e\n",thetak);
           break;
         case    80:
           fprintf(ERROR,"\n Error return from nleq_linfact: fail=%i\n",fail);
           break;
         case    81:
           fprintf(ERROR,"\n Error return from nleq_linsol: fail=%i\n",fail);
           break;
         case    82:
           fprintf(ERROR,"\n Error return from problem function\n");
           break;
         case    83:
           fprintf(ERROR,"\n Error return from Jacobian evaluation function\n");
           break;
         default   :
           fprintf(ERROR,"\n Error, code=%i,  subcode=%i\n",RCODE,fail);
        };
     };
   info->subcode = fail;
   if ( !nleqcalled && io_allocated ) { free(nleq_ioctl); nleq_ioctl=NULL; };
   free(data);
   if (scale_allocated) free(xscale);
   k1 = ( RCODE == 2 ? k : k+1 );
   info->fx         = fxk;
   for (i=1;i<=k1;i++) free(dx[i]);
   if (!nleqcalled) {free(dx[0]); nleq_linfree();};
   free(dx);  free(v); free(sigma);
   info->precision  = normdx;
   info->normdx     = normdx;
   info->iter       = (nleqcalled ? (info->iter)+k-1 : k );
   info->nofunevals = nfcn;
   info->nojacevals = njac;
}

int qne_initscale(int n, double **scale, double *x, double tol)
{  int i, rcode;
   rcode = nleq_fwalloc(n,scale,"scale");     if ( rcode !=0 ) return rcode;
   for (i=0;i<n;i++) (*scale)[i] = 1.0;
   return 0;
}

void qne_monitor(int k, int n, double normdx, double normf,
                 char qnerr_id[9], char fail_reason[7])
{  fprintf(MONITOR," %4i     %12e  %12e  %s  %s\n"
                 ,k,normdx,normf,qnerr_id,fail_reason);
}
