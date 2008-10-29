/*
    CNLSOL -   Call Residuum based and Error based Damped Newton
               and Quasi-Newton methods
               Interface routine for a calling Fortran program.
    
 *  Written by        L. Weimann 
 *  Purpose           Solution of systems of nonlinear equations
 *  Method            Various Newton-type methods
                      (see reference below)
 *  Category          F2a. - Systems of nonlinear equations
 *  Keywords          Nonlinear equations, Newton methods
 *  Version           1.1
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
  
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nleq.h"

typedef void NLEQ_SOLVER(struct NLEQ_FUN,int,double*,struct NLEQ_OPT*,
                         struct NLEQ_INFO*);

static struct NLEQ_FUN  *problem = NULL;
static struct NLEQ_OPT  *opt     = NULL;
static struct NLEQ_INFO *ninfo   = NULL;

extern void cinit_(char filename[6][80])
{ int i,j;
  char null='\0';
  LOGICAL err_eq_mon, err_eq_dat, mon_eq_dat, blank_found;
  if (!problem) problem = malloc(sizeof(struct NLEQ_FUN));
  if (!opt)     opt     = malloc(sizeof(struct NLEQ_OPT));
  if (!ninfo)   ninfo   = malloc(sizeof(struct NLEQ_INFO));
  if (!problem || !opt || !ninfo) 
    { fprintf(stdout,"\n CINIT: Memory allocation error!\n"); exit(-10);};
  opt->errorfile    = NULL;
  opt->monitorfile  = NULL;
  opt->datafile     = NULL;
  opt->iterfile     = NULL;
  opt->resfile      = NULL;
  opt->miscfile     = NULL;
  for (i=0;i<6;i++) 
    { filename[i][79]=' ';
      blank_found = False;
      for(j=0; j<80 && !blank_found ;j++)
        { if ( filename[i][j]==' ' )
            { filename[i][j] = '\0'; blank_found = True; };
        };
     };
  err_eq_mon = !strncmp(filename[0],filename[1],80);
  err_eq_dat = !strncmp(filename[0],filename[2],80);
  mon_eq_dat = !strncmp(filename[1],filename[2],80);
  if ( filename[0][0] != null ) 
    opt->errorfile   = fopen(filename[0],"w");
  if ( filename[1][0] != null && !err_eq_mon )
    opt->monitorfile = fopen(filename[1],"w");
  else if ( err_eq_mon ) 
    opt->monitorfile = opt->errorfile;
  if ( filename[2][0] != null && !mon_eq_dat && !err_eq_dat )
    opt->datafile = fopen(filename[2],"w");
  else if ( mon_eq_dat ) 
    opt->datafile = opt->monitorfile;
  else if ( err_eq_dat ) 
    opt->datafile = opt->errorfile;
  if ( filename[3][0] != null ) 
    opt->iterfile   = fopen(filename[3],"w");
  if ( filename[4][0] != null ) 
    opt->resfile   = fopen(filename[4],"w");
  if ( filename[5][0] != null ) 
    opt->miscfile   = fopen(filename[5],"w");
  if ( !opt->errorfile )   opt->errorfile   = stderr;
  if ( !opt->monitorfile ) opt->monitorfile = stdout;
  if ( !opt->datafile )    opt->datafile    = stdout;
}

extern void cnlsol_(int *n, NLEQ_FFUN *fun, NLEQ_JFUN *jac,
                      double *x, double *scale, double *tol,
                      int *iopt, char pronam[48], int *info)
{ int fail;
  NLEQ_SOLVER *solver=NULL;
  char proname[49];
  double fnorm;
  double *f = (double*) malloc( (*n)*sizeof(double) );
  if (!f) 
    { fprintf(stdout,"\n CNLSOL: Memory allocation error!\n"); exit(-10);};
  strncpy(proname,pronam,48);   proname[48]='\0';
  problem->fun = fun;
  if ( iopt[0] ) problem->jac = jac;
  else           problem->jac = NULL;
  opt->tol = *tol;
  opt->maxiter      = iopt[1];
  opt->nonlin       = iopt[2];
  opt->restricted   = iopt[3];
  opt->errorlevel   = iopt[4];
  opt->monitorlevel = iopt[5];
  opt->datalevel    = iopt[6];
  opt->scaleopt     = iopt[7]; 
  if ( iopt[7] == 2 ) opt->scale = scale;
  if      ( iopt[8] == 1 )    solver = &nleq_err;
  else if ( iopt[8] == 2 )    solver = &nleq_res;
  else if ( iopt[8] == 3 )    solver = &qnerr;
  else if ( iopt[8] == 4 )    solver = &qnres;
  fprintf(opt->monitorfile," problem: %s\n\n",proname);
  fprintf(opt->datafile,   " problem: %s\n\n",proname);
  if (!solver) { info[0] = -2000; return; };
  solver(*problem,*n,x,opt,ninfo);
  fprintf(opt->monitorfile," return code:   % 16i\n",ninfo->rcode);
  fprintf(opt->monitorfile," subcode:       % 16i\n",ninfo->subcode);
  info[0] = ninfo->rcode;
  info[1] = ninfo->subcode;
  info[2] = ninfo->iter;
  info[3] = ninfo->nofunevals;
  info[4] = ninfo->nojacevals;
  *tol    = ninfo->precision;
  fun(n,x,f,&fail);
  fnorm = nleq_norm2(*n,f);
  fprintf(opt->monitorfile," norm of residuum = %12.4e\n\n",fnorm);
  free(f);
  fflush(opt->errorfile);
  fflush(opt->monitorfile);
  fflush(opt->datafile);
  if (opt->iterfile) fflush(opt->iterfile);
  if (opt->resfile)  fflush(opt->resfile);
  if (opt->miscfile) fflush(opt->miscfile);
}


