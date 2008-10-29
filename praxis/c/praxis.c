/*********************************************************************/
/* 	f u n c t i o n     p r a x i s                              */
/*                                                                   */
/* praxis is a general purpose routine for the minimization of a     */
/* function in several variables. the algorithm used is a modifi-    */
/* cation of conjugate gradient search method by powell. the changes */
/* are due to r.p. brent, who gives an algol-w program, which served */
/* as a basis for this function.                                     */
/*                                                                   */
/* references:                                                       */
/*     - powell, m.j.d., 1964. an efficient method for finding       */
/*       the minimum of a function in several variables without      */
/*       calculating derivatives, computer journal, 7, 155-162       */
/*     - brent, r.p., 1973. algorithms for minimization without      */
/*       derivatives, prentice hall, englewood cliffs.               */
/*                                                                   */
/*     problems, suggestions or improvements are always wellcome     */
/*                       karl gegenfurtner   07/08/87                */
/*                                           c - version             */
/*********************************************************************/
/*                                                                   */
/* usage: min = praxis(fun, x, n);                                   */
/*                                                                   */
/*  fun        the function to be minimized. fun is called from      */
/*             praxis with x and n as arguments                      */
/*  x          a double array containing the initial guesses for     */
/*             the minimum, which will contain the solution on       */
/*             return                                                */
/*  n          an integer specifying the number of unknown           */
/*             parameters                                            */
/*  min        praxis returns the least calculated value of fun      */
/*                                                                   */
/* some additional global variables control some more aspects of     */
/* the inner workings of praxis. setting them is optional, they      */
/* are all set to some reasonable default values given below.        */
/*                                                                   */
/*   prin      controls the printed output from the routine.         */
/*             0 -> no output                                        */
/*             1 -> print only starting and final values             */
/*             2 -> detailed map of the minimization process         */
/*             3 -> print also eigenvalues and vectors of the        */
/*                  search directions                                */
/*             the default value is 1                                */
/*  tol        is the tolerance allowed for the precision of the     */
/*             solution. praxis returns if the criterion             */
/*             2 * ||x[k]-x[k-1]|| <= sqrt(macheps) * ||x[k]|| + tol */
/*             is fulfilled more than ktm times.                     */
/*             the default value depends on the machine precision    */
/*  ktm        see just above. default is 1, and a value of 4 leads  */
/*             to a very(!) cautious stopping criterion.             */
/*  step       is a steplength parameter and should be set equal     */
/*             to the expected distance from the solution.           */
/*             exceptionally small or large values of step lead to   */
/*             slower convergence on the first few iterations        */
/*             the default value for step is 1.0                     */
/*  scbd       is a scaling parameter. 1.0 is the default and        */
/*             indicates no scaling. if the scales for the different */
/*             parameters are very different, scbd should be set to  */
/*             a value of about 10.0.                                */
/*  illc       should be set to true (1) if the problem is known to  */
/*             be ill-conditioned. the default is false (0). this    */
/*             variable is automatically set, when praxis finds      */
/*             the problem to be ill-conditioned during iterations.  */
/*  maxfun     is the maximum number of calls to fun allowed. praxis */
/*             will return after maxfun calls to fun even when the   */
/*             minimum is not yet found. the default value of 0      */
/*             indicates no limit on the number of calls.            */
/*             this return condition is only checked every n         */
/*             iterations.                                           */
/*                                                                   */
/*********************************************************************/

#include <math.h>
#include <stdio.h>
#include "machine.h"


/* control parameters */
double tol = SQREPSILON,
       scbd = 1.0,
       step = 1.0;
int    ktm = 1,
       prin = 2,
       maxfun = 0,
       illc = 0;
       
/* some global variables */
static int i, j, k, k2, nl, nf, kl, kt;
static double s, sl, dn, dmin,
       fx, f1, lds, ldt, sf, df,
       qf1, qd0, qd1, qa, qb, qc,
       m2, m4, small, vsmall, large, 
       vlarge, ldfac, t2;
static double d[N], y[N], z[N],
       q0[N], q1[N], v[N][N];

/* these will be set by praxis to point to it's arguments */
static int n;
double *x;
double (*fun)();

/* these will be set by praxis to the global control parameters */
static double h, macheps, t;

double 
random()	/* return random no between 0 and 1 */
{
   return (double)(rand()%(8192*2))/(double)(8192*2);
}

sort()		/* d and v in descending order */
{
   int k, i, j;
   double s;

   for (i=0; i<n-1; i++) {
       k = i; s = d[i];
       for (j=i+1; j<n; j++) {
           if (d[j] > s) {
	      k = j;
	      s = d[j];
	   }
       }
       if (k > i) {
	  d[k] = d[i];
	  d[i] = s;
	  for (j=0; j<n; j++) {
	      s = v[j][i];
	      v[j][i] = v[j][k];
	      v[j][k] = s;
	  }
       }
   }
}

print()		/* print a line of traces */
{
   int i;

   printf("\n");
   printf("... chi square reduced to ... %20.10e\n", fx);
   printf("... after %u function calls ...\n", nf);
   printf("... including %u linear searches ...\n", nl);
   vecprint("... current values of x ...", x, n);
}

matprint(s, v, n)
char *s;
double v[N][N];
{
   int k, i;
   
   printf("%s\n", s);
   for (k=0; k<n; k++) {
       for (i=0; i<n; i++) {
           printf("%20.10e ", v[k][i]);
       }
       printf("\n");
   }
}

vecprint(s, x, n)
char *s;
double x[N];
{
   int i;
   
   printf("%s\n", s);
   for (i=0; i<n; i++)
       printf("%20.10e ", x[i]);
   printf("\n");
}

#ifdef MSDOS
static double tflin[N];
#endif

double
flin(l, j)
double l;
{
   int i;
#ifndef MSDOS
   double tflin[N];
#endif   

   if (j != -1) {		/* linear search */
      for (i=0; i<n; i++)
          tflin[i] = x[i] + l *v[i][j];
   }
   else {			/* search along parabolic space curve */
      qa = l*(l-qd1)/(qd0*(qd0+qd1));
      qb = (l+qd0)*(qd1-l)/(qd0*qd1);
      qc = l*(l+qd0)/(qd1*(qd0+qd1));
      for (i=0; i<n; i++)
          tflin[i] = qa*q0[i]+qb*x[i]+qc*q1[i];
   }
   nf++;
   return (*fun)(tflin, n);
}

min(j, nits, d2, x1, f1, fk)
double *d2, *x1, f1;
{
   int k, i, dz;
   double x2, xm, f0, f2, fm, d1, t2,
          s, sf1, sx1;

   sf1 = f1; sx1 = *x1;
   k = 0; xm = 0.0; fm = f0 = fx; dz = *d2 < macheps;
   /* find step size */
   s = 0;
   for (i=0; i<n; i++) s += x[i]*x[i];
   s = sqrt(s);
   if (dz)
      t2 = m4*sqrt(fabs(fx)/dmin + s*ldt) + m2*ldt;
   else
      t2 = m4*sqrt(fabs(fx)/(*d2) + s*ldt) + m2*ldt;
   s = s*m4 + t;
   if (dz && t2 > s) t2 = s;
   if (t2 < small) t2 = small;
   if (t2 > 0.01*h) t2 = 0.01 * h;
   if (fk && f1 <= fm) {
      xm = *x1;
      fm = f1;
   }
   if (!fk || fabs(*x1) < t2) {
      *x1 = (*x1 > 0 ? t2 : -t2);
      f1 = flin(*x1, j);
   }
   if (f1 <= fm) {
      xm = *x1;
      fm = f1;
   }
next:
   if (dz) {
      x2 = (f0 < f1 ? -(*x1) : 2*(*x1));
      f2 = flin(x2, j);
      if (f2 <= fm) {
         xm = x2;
	 fm = f2;
      }
      *d2 = (x2*(f1-f0) - (*x1)*(f2-f0))/((*x1)*x2*((*x1)-x2));
   }
   d1 = (f1-f0)/(*x1) - *x1**d2; dz = 1;
   if (*d2 <= small) {
      x2 = (d1 < 0 ? h : -h);
   }
   else {
      x2 = - 0.5*d1/(*d2);
   }
   if (fabs(x2) > h)
      x2 = (x2 > 0 ? h : -h);
try:
   f2 = flin(x2, j);
   if ((k < nits) && (f2 > f0)) {
      k++;
      if ((f0 < f1) && (*x1*x2 > 0.0))
         goto next;
      x2 *= 0.5;
      goto try;
   }
   nl++;
   if (f2 > fm) x2 = xm; else fm = f2;
   if (fabs(x2*(x2-*x1)) > small) {
      *d2 = (x2*(f1-f0) - *x1*(fm-f0))/(*x1*x2*(*x1-x2));
   }
   else {
      if (k > 0) *d2 = 0;
   }
   if (*d2 <= small) *d2 = small;
   *x1 = x2; fx = fm;
   if (sf1 < fx) {
      fx = sf1;
      *x1 = sx1;
   }
   if (j != -1)
      for (i=0; i<n; i++)
          x[i] += (*x1)*v[i][j];
}

quad()	/* look for a minimum along the curve q0, q1, q2	*/
{
   int i;
   double l, s;

   s = fx; fx = qf1; qf1 = s; qd1 = 0.0;
   for (i=0; i<n; i++) {
       s = x[i]; l = q1[i]; x[i] = l; q1[i] = s;
       qd1 = qd1 + (s-l)*(s-l);
   }
   s = 0.0; qd1 = sqrt(qd1); l = qd1;
   if (qd0>0.0 && qd1>0.0 &&nl>=3*n*n) {
      min(-1, 2, &s, &l, qf1, 1);
      qa = l*(l-qd1)/(qd0*(qd0+qd1));
      qb = (l+qd0)*(qd1-l)/(qd0*qd1);
      qc = l*(l+qd0)/(qd1*(qd0+qd1));
   }
   else {
      fx = qf1; qa = qb = 0.0; qc = 1.0;
   }
   qd0 = qd1;
   for (i=0; i<n; i++) {
       s = q0[i]; q0[i] = x[i];
       x[i] = qa*s + qb*x[i] + qc*q1[i];
   }
}

double
praxis(_fun, _x, _n)
double (*_fun)();
double _x[N];
{
   /* init global extern variables and parameters */
   macheps = EPSILON; h = step; t = tol;
   n = _n; x = _x; fun = _fun;
   small = macheps*macheps; vsmall = small*small;
   large = 1.0/small; vlarge = 1.0/vsmall;
   m2 = sqrt(macheps); m4 = sqrt(m2);
   ldfac = (illc ? 0.1 : 0.01);
   nl = kt = 0; nf = 1; fx = (*fun)(x, n); qf1 = fx;
   t2 = small + fabs(t); t = t2; dmin = small;

   if (h < 100.0*t) h = 100.0*t;
   ldt = h;
   for (i=0; i<n; i++) for (j=0; j<n; j++)
       v[i][j] = (i == j ? 1.0 : 0.0);
   d[0] = 0.0; qd0 = 0.0;
   for (i=0; i<n; i++) q1[i] = x[i];
   if (prin > 1) {
      printf("\n------------- enter function praxis -----------\n");
      printf("... current parameter settings ...\n");
      printf("... scaling ... %20.10e\n", scbd);
      printf("...   tol   ... %20.10e\n", t);
      printf("... maxstep ... %20.10e\n", h);
      printf("...   illc  ... %20u\n", illc);
      printf("...   ktm   ... %20u\n", ktm);
      printf("... maxfun  ... %20u\n", maxfun);
   }
   if (prin) print();

mloop:
   sf = d[0];
   s = d[0] = 0.0;

   /* minimize along first direction */
   min(0, 2, &d[0], &s, fx, 0);
   if (s <= 0.0)
      for (i=0; i < n; i++)
          v[i][0] = -v[i][0];
   if ((sf <= (0.9 * d[0])) || ((0.9 * sf) >= d[0]))
      for (i=1; i<n; i++)
          d[i] = 0.0;
   for (k=1; k<n; k++) {
       for (i=0; i<n; i++)
           y[i] = x[i];
       sf = fx;
       illc = illc || (kt > 0);
next:
       kl = k;
       df = 0.0;
       if (illc) {        /* random step to get off resolution valley */
          for (i=0; i<n; i++) {
              z[i] = (0.1 * ldt + t2 * pow(10.0,(double)kt)) * (random() - 0.5);
              s = z[i];
              for (j=0; j < n; j++)
                  x[j] += s * v[j][i];
  	  }
          fx = (*fun)(x, n);
          nf++;
       }
       /* minimize along non-conjugate directions */ 
       for (k2=k; k2<n; k2++) {  
           sl = fx;
           s = 0.0;
           min(k2, 2, &d[k2], &s, fx, 0);
           if (illc) {
	      double szk = s + z[k2];
              s = d[k2] * szk*szk;
	   }
           else 
	      s = sl - fx;
           if (df < s) {
              df = s;
              kl = k2;
           }
       }
       if (!illc && (df < fabs(100.0 * macheps * fx))) {
          illc = 1;
          goto next;
       }
       if ((k == 1) && (prin > 1))
          vecprint("\n... New Direction ...",d,n);
       /* minimize along conjugate directions */ 
       for (k2=0; k2<=k-1; k2++) {
           s = 0.0;
           min(k2, 2, &d[k2], &s, fx, 0);
       }
       f1 = fx;
       fx = sf;
       lds = 0.0;
       for (i=0; i<n; i++) {
           sl = x[i];
           x[i] = y[i];
           y[i] = sl - y[i];
           sl = y[i];
           lds = lds + sl*sl;
       }
       lds = sqrt(lds);
       if (lds > small) {
          for (i=kl-1; i>=k; i--) {
              for (j=0; j < n; j++)
                  v[j][i+1] = v[j][i];
                  d[i+1] = d[i];
              }
              d[k] = 0.0;
              for (i=0; i < n; i++)
                  v[i][k] = y[i] / lds;
              min(k, 4, &d[k], &lds, f1, 1);
              if (lds <= 0.0) {
                 lds = -lds;
                 for (i=0; i<n; i++)
                     v[i][k] = -v[i][k];
              }
       }
       ldt = ldfac * ldt;
       if (ldt < lds)
          ldt = lds;
       if (prin > 1)
          print();
       t2 = 0.0;
       for (i=0; i<n; i++)
           t2 += x[i]*x[i];
       t2 = m2 * sqrt(t2) + t;
       if (ldt > (0.5 * t2))
          kt = 0;
       else 
	  kt++;
       if (kt > ktm)
          goto fret; 
   }
   /*  try quadratic extrapolation in case    */
   /*  we are stuck in a curved valley        */
   quad();
   dn = 0.0;
   for (i=0; i<n; i++) {
       d[i] = 1.0 / sqrt(d[i]);
       if (dn < d[i])
          dn = d[i];
   }
   if (prin > 2)
      matprint("\n... New Matrix of Directions ...",v,n);
   for (j=0; j<n; j++) {
       s = d[j] / dn;
       for (i=0; i < n; i++)
           v[i][j] *= s;
   }
   if (scbd > 1.0) {       /* scale axis to reduce condition number */
      s = vlarge;
      for (i=0; i<n; i++) {
          sl = 0.0;
          for (j=0; j < n; j++)
              sl += v[i][j]*v[i][j];
          z[i] = sqrt(sl);
          if (z[i] < m4)
             z[i] = m4;
          if (s > z[i])
             s = z[i];
      }
      for (i=0; i<n; i++) {
          sl = s / z[i];
          z[i] = 1.0 / sl;
          if (z[i] > scbd) {
             sl = 1.0 / scbd;
             z[i] = scbd;
          }
      }
   }
   for (i=1; i<n; i++)
       for (j=0; j<=i-1; j++) {
           s = v[i][j];
           v[i][j] = v[j][i];
           v[j][i] = s;
       }
   minfit(n, macheps, vsmall, v, d);
   if (scbd > 1.0) {
      for (i=0; i<n; i++) {
          s = z[i];
          for (j=0; j<n; j++)
              v[i][j] *= s;
      }
      for (i=0; i<n; i++) {
          s = 0.0;
          for (j=0; j<n; j++)
              s += v[j][i]*v[j][i];
          s = sqrt(s);
          d[i] *= s;
          s = 1.0 / s;
          for (j=0; j<n; j++)
              v[j][i] *= s;
      }
   }
   for (i=0; i<n; i++) {
       if ((dn * d[i]) > large)
          d[i] = vsmall;
       else if ((dn * d[i]) < small)
          d[i] = vlarge;
       else 
          d[i] = pow(dn * d[i],-2.0);
   }
   sort();               /* the new eigenvalues and eigenvectors */
   dmin = d[n-1];
   if (dmin < small)
      dmin = small;
   illc = (m2 * d[0]) > dmin;
   if ((prin > 2) && (scbd > 1.0))
      vecprint("\n... Scale Factors ...",z,n);
   if (prin > 2)
      vecprint("\n... Eigenvalues of A ...",d,n);
   if (prin > 2)
      matprint("\n... Eigenvectors of A ...",v,n);

   if ((maxfun > 0) && (nl > maxfun)) {
      if (prin)
	 printf("\n... maximum number of function calls reached ...\n");
      goto fret;
   }
   goto mloop; 	 /* back to main loop */

fret:
   if (prin > 0) {
         vecprint("\n... Final solution is ...", x, n);
         printf("\n... ChiSq reduced to %20.10e ...\n", fx);
	 printf("... after %20u function calls.\n", nf);
   }
   
   return(fx);
}

