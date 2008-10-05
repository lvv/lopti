#include <string.h>

#include <R.h>
#include <Rdefines.h>

SEXP R_UObyQA(SEXP x, SEXP fn, SEXP control, SEXP rho);

int is_valid(SEXP x);

void UObyQA_a(unsigned int *n,double *x,double *rhobeg,double *rhoend,
              unsigned int *iprint,unsigned int *maxfun,double *w,
              double *fnscale,double *parscale);
void UObyQA_b(unsigned int *n,double *x,double *rhobeg,
              double *rhoend,unsigned int *iprint,
              unsigned int *maxfun,unsigned int *npt,
              double *xbase,double *xopt,double *xnew,double *xpt,
              double *pq,double *pl,double *h__,double *g,
              double *d__,double *vlag,double *w,
              double *fnscale,double *parscale);
void lagmax(unsigned int *n,double *g,double *h__,double *rho,
            double *d__,double *v, double *vmax);
void trstep(unsigned int *n,double *g,double *h__,
            double *delta,double *tol,double *d__,
            double *gg,double *td,double *tn,
            double *w,double *piv,double *z__,double *evalue);
double d_sign(double *a,double *b);
void getHessian(double *pq);

unsigned int verbose;
SEXP par, expr, env, fval, H;
unsigned int *nf;
double f, *Ha;
unsigned int npar;
