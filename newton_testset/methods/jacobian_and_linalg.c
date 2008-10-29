/*
    Jacobian related routines for the NewtonLib package.
    
 *  Written by        L. Weimann 
 *  Purpose           Jacobian and linear algebra related routines
 *  Category          ???. - Utilities
 *  Keywords          Linear system solution, Numerical differentiation,
                      Jacobian scaling
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
 
   The implementation of the routines within this file depend on the
   chosen Jacobian format and the selected linear system solver. 
   The JACOBIAN type used within this file is defined in the header 
   file "jacobian.h". By changing the contents of this header file,
   you may be able to implement any type of Jacobians, for example,
   sparse Jacobians. The current implementation supports Jacobians 
   stored in full storage mode, and uses the Fortran LAPACK solver 
   routines DGETRF for the LU-decomposition and DGETRS for the 
   solution (i.e. backsubstitution).
   
   The routines and functions are the following:
   ---------------------------------------------
   
   int nleq_jacalloc(int n, JACOBIAN **jac, int *ldjac,
                     struct NLEQ_OPT *opt)
   This function must allocate memory for the Jacobian, using either
   directly the malloc function or the functions nleq_fwalloc,
   nleq_iwalloc, nleq_pfwalloc (see file "utils.c"), as suitable. 
   The parameters are:
   int      n        (input)   The number of equations and unknowns
   JACOBIAN **jac    (output)  A pointer of type (JACOBIAN *) returned
                               by the memory allocation function, and
                               pointing to the reserved memory.
   int     *ldjac    (output)  A positive int number specifying the
                               leading dimension size of the Jacobian
                               (asuming the Jacobian to be a two dimensional
                                array of (double)). 
   struct NLEQ_OPT *opt (input) The input options of the main routine.
   
   The int return code of nleq_jacalloc must be 0, if everything worked fine
   with the memory allocation, and otherwise a number from -999 to -996.
   ---
   Note: You may expand the struct NLEQ_OPT, which is defined in the "nleq.h"
         header file, by adding fields which may describe characteristics
         of your Jacobian, for example, such as lower and upper bandwidth
         of a Jacobian which will be stored in band mode.
   ---
   int nleq_linfact(int n, JACOBIAN *jac, struct NLEQ_OPT *opt)
   This function must call the linear algebra routine for computation of
   the (LU-)decomposition.
   The parameters are:
   int      n        (input)   The number of equations and unknowns
   JACOBIAN *jac     (in/out)  On input:  The Jacobian
                               On output: The decomposition of the Jacobian.
   struct NLEQ_OPT *opt (input) The input options of the main routine.
   
   The int return code must be 0, if the matrix decomposition was successful.
   A singular matrix must be indicated by a positive return code, any other
   error by a negative return code.
   ---
   Note: To keep all needed information of the matrix decomposition, it
         may be necessary to allocate additional memory on the first call
         of nleq_linfact, i.e. to store pivot indices or other information.
   ---
   int nleq_linsol(int n, double *b, struct NLEQ_OPT *opt)
   This function must call the linear algebra routine for the solution
   (i.e. backsubstitution) of the linear system, after the Jacobian has
   been factorized through nleq_linfact.
   The parameters are:
   int      n        (input)   The number of equations and unknowns
   double   *b       (in/out)  On input: The right hand side vector of the
                                         linear equations system.
                               On output: The solution of the linear system.
   struct NLEQ_OPT *opt (input) The input options of the main routine.
   
   A successful completion must be indicated by a return code 0, a failure
   condition by any nonzero return code.
   ---
   void nleq_linfree()
   This routine must free any memory, which has been allocated either by 
   the routine nleq_jacalloc or the routine nleq_linfact.
   ---
   int nleq_numjac(NLEQ_FFUN *f, int n, double *x, double *fx,
                   double *scale, JACOBIAN *jac, int *nfcn,
                   struct NLEQ_OPT *opt)
   This function should compute the Jacobian by numerical difference
   approximation. The main codes will also work if this function just
   does nothing more than returning a nonzero value, but for this case 
   the Jacobian routine must be always supplied by the user.
   The parameters are:
   NLEQ_FFUN *f       (input)   A pointer to the user problem function f.
   int       n        (input)   The number of equations and unknowns
                                and size of any vectors.
   double    *x       (input)   The input vector where to compute the Jacobian.
                                The content of x must be preserved on output!
   double    *fx      (input)   The vector f(x).
                                The content of fx must be preserved on output!
   double    *scale   (input)   The scaling vector.
                                The content of scale must be preserved on
                                output!
   JACOBIAN *jac      (output)  A pointer to the Jacobian storage.
   int      *nfcn     (in/out)  The value of *nfcn must be incremented by one
                                for each call of the user problem function.
   struct NLEQ_OPT *opt (input) The input options of the main routine.

   A successful completion must be indicated by a return code 0, a failure
   condition by any nonzero return code.
   ---
   void nleq_jacrow_scale(int n, JACOBIAN *jac, double *scale, 
                          struct NLEQ_OPT *opt)
   This routine must scale the rows of the Jacobian and change the sign
   of all Jacobian elements. The function is used by the main codes
   nleq_res and qnres.
   The scaling must accomplish the following transformation:
   
   jac[i][j] = - jac[i][j]/scale[i] for i=0,...,n-1  and j=0,...,n-1.
   
   The parameters are:
   int       n        (input)   The number of equations and unknowns
                                and size of any vectors.
   JACOBIAN *jac      (in/out)  A pointer to the Jacobian storage.
   double    *scale   (input)   The scaling vector.
                                The content of scale must be preserved on
                                output!
   struct NLEQ_OPT *opt (input) The input options of the main routine.
   ---
   void nleq_jaccolumn_scale(int n, JACOBIAN *jac, double *scale, 
                             struct NLEQ_OPT *opt)
   This routine must scale the columns of the Jacobian and change the sign
   of all Jacobian elements. The funtion is used by the main codes
   nleq_err and qnerr.
   The scaling must accomplish the following transformation:
   
   jac[i][j] = - jac[i][j]*scale[j] for j=0,...,n-1  and i=0,...,n-1.
   
   The parameters are:
   int       n        (input)   The number of equations and unknowns
                                and size of any vectors.
   JACOBIAN *jac      (in/out)  A pointer to the Jacobian storage.
   double    *scale   (input)   The scaling vector.
                                The content of scale must be preserved on
                                output!
   struct NLEQ_OPT *opt (input) The input options of the main routine.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nleq.h"

#ifdef FMAT
#define JMODE 'N'
#else
#define JMODE 'T'
#endif

#define AJDEL 1.0e-8
#define AJMIN 1.0e-4


extern void dgetrf_(int *m, int *n, double *a, int *lda,
                    int *ipiv, int *info);
extern void dgetrs_(char *trans, int *n, int *nrhs, double *a,
                    int *lda, int *ipiv, double *b,
                    int *ldb, int *info);

static int *pivot=NULL;
static JACOBIAN *mat=NULL;
static double *v=NULL;

int nleq_jacalloc(int n, JACOBIAN **jac, int *ldjac,
                  struct NLEQ_OPT *opt)
{  *ldjac = n;
   return nleq_fwalloc(n*n,jac,"jac");
}

int nleq_linfact(int n, JACOBIAN *jac, struct NLEQ_OPT *opt)
{ int fail=0; 
  if (!pivot) 
     { fail=nleq_iwalloc(n,&pivot,"pivot");
       if ( fail !=0 ) return fail;
     };
  mat=jac;
  dgetrf_(&n,&n,mat,&n,pivot,&fail);
  return fail;
}

int nleq_linsol(int n, double *b, struct NLEQ_OPT *opt)
{ int fail=0, nrhs=1;
  char mode=JMODE;
  if (mat) dgetrs_(&mode,&n,&nrhs,mat,&n,pivot,b,&n,&fail);
  else     fail = -980;
  return fail;
}

void nleq_linfree()
{ if(pivot) free(pivot); pivot = NULL;
  if (v)    free(v);     v     = NULL;
  if (mat)  free(mat);   mat   = NULL;
}


int nleq_numjac(NLEQ_FFUN *f, int n, double *x, double *fx,
                double *scale, JACOBIAN *jac, int *nfcn, struct NLEQ_OPT *opt)
{ int i, k, kn;
  int fail=0;
  double w, u;
  if (!v) fail=nleq_fwalloc(n,&v,"v");
  if (fail!=0) return -992;
  for (k=0;k<n;k++)
    { w = x[k];
      if (scale) u = MAX( fabs(x[k]) , scale[k] );
      else       u = fabs(x[k]);
      u = MAX(u,AJMIN)*AJDEL*SIGN(x[k]);
      x[k] = w+u;
      f(&n,x,v,&fail);  (*nfcn)++;
      if ( fail != 0 ) return fail;
      x[k] = w;
#ifndef FMAT 
      for (i=0;i<n;i++) jac[i*n+k] = (v[i]-fx[i]) / u ;  
#else    
      kn = k*n;
      for (i=0;i<n;i++) jac[kn+i] = (v[i]-fx[i]) / u ;
#endif
    };
  return fail;
}

void nleq_jacrow_scale(int n, JACOBIAN *jac, double *scale, 
                       struct NLEQ_OPT *opt)
{ int i,j;
  double s;
       /* scale jacobian and change sign of it */
  for (i=0;i<n;i++)
    { s = -scale[i];
#ifndef FMAT
      for (j=i*n;j<(i+1)*n;j++) jac[j] /= s;
#else 
      for (j=i;j<n*n;j+=n) jac[j] /= s;
#endif
    };
}

void nleq_jaccolumn_scale(int n, JACOBIAN *jac, double *scale, 
                          struct NLEQ_OPT *opt)
{ int i,j;
  double s;
       /* scale jacobian and change sign of it */
  for (i=0;i<n;i++)
    { s = -scale[i];
#ifdef FMAT
      for (j=i*n;j<(i+1)*n;j++) jac[j] *= s;
#else 
      for (j=i;j<n*n;j+=n) jac[j] *= s;
#endif
    };
}
