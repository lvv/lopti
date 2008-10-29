/*
    Common declarations for programs of the NewtonLib package.
    
 *  Written by        L. Weimann 
 *  Version           1.1
 *  Revision          May 2006
 *  Latest Change     May 2006
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
#define RCODE info->rcode
#define MIN(A,B)  ( A < B ? A : B )
#define MAX(A,B)  ( A > B ? A : B )
#define SIGN(A)   ( A > 0 ? 1 : -1 )

#define SMALL  1.0e-150
#define EPMACH 1.0e-17

#include "jacobian.h"

typedef void NLEQ_FFUN(int*,double*,double*,int*);
typedef void NLEQ_JFUN(int*,int*,double*,JACOBIAN*,int*);
typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL ;
typedef enum {False=0, True=1} LOGICAL ;
typedef enum {Mildly_Nonlinear=2, Highly_Nonlinear=3, Extremely_Nonlinear=4}
             PROBLEM_TYPE ;

struct NLEQ_FUN
{
   NLEQ_FFUN *fun;
   NLEQ_JFUN *jac;
};

struct NLEQ_OPT
{
   double tol;
   int maxiter;
   LOGICAL restricted, nleqcalled; 
   enum {StandardScale=0,StartValueScale=1} scaleopt;
   PRINT_LEVEL errorlevel, monitorlevel, datalevel;
   FILE *errorfile, *monitorfile, *datafile,
        *iterfile, *resfile, *miscfile;
   PROBLEM_TYPE nonlin;
   double *scale;
};

struct NLEQ_INFO
{
   double precision, normdx;
   double *fx, *dx;
   int iter, rcode, subcode, nofunevals, nojacevals;
};

struct NLEQ_DATA
{
  double *fx, *dx;
  double normf, normdx, normdxbar, lambda, theta;
  enum { QNRES=0, NLEQ_RES=1, QNERR=2,NLEQ_ERR=3 } codeid;
  enum {Initial=1,Intermediate=2,Solution=3,Final=4} mode;
};

struct NLEQ_IO
{
   FILE *errfile, *monfile, *datfile,
        *iterfile, *resfile, *miscfile;
   PRINT_LEVEL errlevel, monlevel, datlevel;
};

#define ERRORLEVEL   nleq_ioctl->errlevel
#define MONITORLEVEL nleq_ioctl->monlevel
#define DATALEVEL    nleq_ioctl->datlevel
#define ERROR        nleq_ioctl->errfile
#define MONITOR      nleq_ioctl->monfile
#define DATA         nleq_ioctl->datfile
#define FITER        nleq_ioctl->iterfile
#define FRES         nleq_ioctl->resfile
#define FMISC        nleq_ioctl->miscfile

extern void daxpy_(int *n, double *alpha, double *x, int *incx,
                   double *y, int *incy);

/* routines defined in utils.c */
int    nleq_fwalloc(int size, double **ptr, char vname[]);
int    nleq_iwalloc(int size, int **ptr, char vname[]);
int    nleq_pfwalloc(int size, double ***ptr, char vname[]);
double nleq_scaled_norm2(int n, double *v, double *scale);
double nleq_scaled_sprod(int n, double *v1, double *v2, double *scale);
double nleq_norm2(int n, double *v);
void   nleq_scale(int n, double *v1, double *v2, double *scale);
void   nleq_descale(int n, double *v1, double *v2, double *scale);
void   nleq_dataout(int k, int n, double *x, struct NLEQ_DATA *data);
int    nleq_parcheck_and_print(int n, struct NLEQ_OPT *opt,
                               struct NLEQ_FUN fun,
                               int nleq_code);
                               
/* routines defined in jacobian_and_linalg.c */
int    nleq_jacalloc(int n, JACOBIAN **jac, int *ldjac,
                  struct NLEQ_OPT *opt);
int    nleq_linfact(int n, JACOBIAN *jac, struct NLEQ_OPT *opt);
int    nleq_linsol(int n, double *b, struct NLEQ_OPT *opt);
int    nleq_numjac(NLEQ_FFUN *f, int n, double *x, double *fx,
                   double *scale, JACOBIAN *jac, int *nfcn,
                   struct NLEQ_OPT *opt);
void    nleq_linfree();
void    nleq_jacrow_scale(int n, JACOBIAN *jac, double *scale, 
                          struct NLEQ_OPT *opt);
void    nleq_jaccolumn_scale(int n, JACOBIAN *jac, double *scale, 
                             struct NLEQ_OPT *opt);

/* main routines ( see file <routine-name>.c ) */
void qnres(struct NLEQ_FUN fun, int n, double *x,
           struct NLEQ_OPT *opt, struct NLEQ_INFO *info);
void nleq_res(struct NLEQ_FUN fun, int n, double *x,
              struct NLEQ_OPT *opt, struct NLEQ_INFO *info);
void qnerr(struct NLEQ_FUN fun, int n, double *x,
           struct NLEQ_OPT *opt, struct NLEQ_INFO *info);
void nleq_err(struct NLEQ_FUN fun, int n, double *x,
              struct NLEQ_OPT *opt, struct NLEQ_INFO *info);
