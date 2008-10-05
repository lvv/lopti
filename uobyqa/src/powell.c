#include "powell.h"

#ifndef _MAX
#define _MAX(a, b) (a > b ? a : b)
#endif

#ifndef _MIN
#define _MIN(a, b) (a < b ? a : b)
#endif

SEXP R_UObyQA(SEXP x, SEXP fn, SEXP control, SEXP rho)
{
  SEXP ret, nfun;
  double *p_par, *w, rhobeg, rhoend, fnscale, *parscale;
  unsigned int nw, maxfun;
  par = x;
  env = rho;
  expr = fn;
  verbose = INTEGER(VECTOR_ELT(control, 0))[0];
  rhobeg = REAL(VECTOR_ELT(control, 1))[0];
  rhoend = REAL(VECTOR_ELT(control, 2))[0];
  maxfun = (unsigned int)REAL(VECTOR_ELT(control, 3))[0];
  fnscale = REAL(AS_NUMERIC(VECTOR_ELT(control, 4)))[0];
  parscale = REAL(AS_NUMERIC(VECTOR_ELT(control, 5)));
  p_par = NUMERIC_POINTER(par);
  PROTECT(nfun = allocVector(INTSXP, 1));
  nf = (unsigned int*)INTEGER_POINTER(nfun);
  npar = (unsigned int)LENGTH(x);
  PROTECT(H = allocVector(REALSXP, npar * npar));
  Ha = NUMERIC_POINTER(H);
  nw = (unsigned int)(pow((double)npar, 4) +
                      pow((double)8 * npar, 3) +
                      23 * npar * npar + 42 * npar +
                      _MAX(2 * npar * npar + 4, 18 * npar))/4 + 1;
  w = malloc(nw * sizeof(double));
  UObyQA_a(&npar, p_par, &rhobeg, &rhoend, &verbose, &maxfun, w, &fnscale, parscale);
  free(w);
  PROTECT(ret = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ret, 0, par);
  SET_VECTOR_ELT(ret, 1, fval);
  SET_VECTOR_ELT(ret, 2, nfun);
  SET_VECTOR_ELT(ret, 3, H);
  UNPROTECT(3);

  return ret;
}

int is_valid(SEXP x)
{
  long i, n;
  double *y = NUMERIC_POINTER(x);
  n = LENGTH(x);
  for(i = 0; i < n; i++)
    if(!R_FINITE(y[i])) return 0;
  return 1;
}

// Unconstrained Optimization by Quadratic Approximation
void UObyQA_a(unsigned int *n,double *x,double *rhobeg,double *rhoend,
              unsigned int *iprint,unsigned int *maxfun,double *w,
              double *fnscale,double *parscale)
{
  // http://plato.asu.edu/topics/problems/nlounres.html

  // This subroutine seeks the least value of a function of many variables,
  // by a trust region method that forms quadratic models by interpolation.
  // The algorithm is described in "UOBYQA: unconstrained optimization by
  // quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA14,
  // University of Cambridge. The arguments of the subroutine are as follows.

  // N must be set to the number of variables and must be at least two.
  // Initial values of the variables must be set in X(1),X(2),...,X(N). They
  //   will be changed to the values that give the least calculated F.
  // RHOBEG and RHOEND must be set to the initial and final values of a trust
  //   region radius, so both must be positive with RHOEND<=RHOBEG. Typically
  //   RHOBEG should be about one tenth of the greatest expected change to a
  //   variable, and RHOEND should indicate the accuracy that is required in
  //   the final values of the variables.
  // The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
  //   amount of printing. Specifically, there is no output if IPRINT=0 and
  //   there is output only at the return if IPRINT=1. Otherwise, each new
  //   value of RHO is printed, with the best vector of variables so far and
  //   the corresponding value of the objective function. Further, each new
  //   value of F with its variables are output if IPRINT=3.
  // MAXFUN must be set to an upper bound on the number of calls of CALFUN.
  // The array W will be used for working space. Its length must be at least
  //   ( N**4 + 8*N**3 + 23*N**2 + 42*N + max [ 2*N**2 + 4, 18*N ] ) / 4.
  // H_2der will provide the second derivatives that TRSTEP and LAGMAX require.
  // G_1der will provide the first derivatives that TRSTEP and LAGMAX require.

  // Function CALFUN (N,X,F) must be provided by the user. It must set F to
  // the value of the objective function for the variables X(1),X(2),...,X(N).

  // Partition the working space array, so that different parts of it can be
  // treated separately by the subroutine that performs the main calculation.

  static unsigned int id, ig, ih, iw;
  static unsigned int ixb, ipl, ipq, ivl, ixn, ixo, ixp, npt;

  // Parameter adjustments
  --w;
  --x;
  --parscale;

  // Function Body
  npt = (*n * *n + *n * 3 + 2) / 2;
  ixb = 1;
  ixo = ixb + *n;
  ixn = ixo + *n;
  ixp = ixn + *n;
  ipq = ixp + *n * npt;
  ipl = ipq + npt - 1;
  ih = ipl + (npt - 1) * npt;
  ig = ih + *n * *n;
  id = ig + *n;
  ivl = ih;
  iw = id + *n;
  UObyQA_b(n,&x[1],rhobeg,rhoend,iprint,maxfun,
           &npt,&w[ixb],&w[ixo],&w[ixn],&w[ixp],&w[ipq],&w[ipl],
           &w[ih],&w[ig],&w[id],&w[ivl],&w[iw],
           fnscale,&parscale[1]);

  // get the Hessian for confidence interval
  getHessian(&w[ipq+*n]);

  return;
}

void UObyQA_b(unsigned int *n,double *x,double *rhobeg,
              double *rhoend,unsigned int *iprint,
              unsigned int *maxfun,unsigned int *npt,
              double *xbase,double *xopt,double *xnew,double *xpt,
              double *pq,double *pl,double *h__,double *g,
              double *d__,double *vlag,double *w,
              double *fnscale,double *parscale)
{
  // The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to
  //   the corresponding arguments in SUBROUTINE UOBYQA.
  // NPT is set by UOBYQA to (N*N+3*N+2)/2 for the above dimension statement.
  // XBASE will contain a shift of origin that reduces the contributions from
  //   rounding errors to values of the model and Lagrange functions.
  // XOPT will be set to the displacement from XBASE of the vector of
  //   variables that provides the least calculated F so far.
  // XNEW will be set to the displacement from XBASE of the vector of
  //   variables for the current calculation of F.
  // XPT will contain the interpolation point coordinates relative to XBASE.
  // PQ will contain the parameters of the quadratic model.
  // PL will contain the parameters of the Lagrange functions.
  // H will provide the second derivatives that TRSTEP and LAGMAX require.
  // G will provide the first derivatives that TRSTEP and LAGMAX require.
  // D is reserved for trial steps from XOPT, except that it will contain
  //   diagonal second derivatives during the initialization procedure.
  // VLAG will contain the values of the Lagrange functions at a new point X.
  // The array W will be used for working space. Its length must be at least
  // max [ 6*N, ( N**2 + 3*N + 2 ) / 2 ].

  // Table of constant values
  static double c_b35 = 1.5;

  // System generated locals
  int xpt_dim1,xpt_offset,pl_dim1,pl_offset,h_dim1,h_offset;
  double d1,d2;

  // Local variables
  static double diff, half;
  static int knew;
  static double temp, fopt, sumg, sumh;
  static int kopt;
  static double zero, vmax;
  static unsigned int i, j, k, nptm;
  static double fbase, delta, fsave, tempa;
  static int ksave;
  static double ratio, dnorm, vquad;
  static unsigned int ktemp;
  static double estim, rhosq, wmult;
  static unsigned int ih, ip, iq, iw;
  static double ddknew, evalue, detrat;
  static unsigned int nftest;
  static double errtol, sixthm;
  static double tworsq, one, rho;
  static unsigned int nnp;
  static double tol, sum, two;
  static int jswitch;
  static double distest;

  // Parameter adjustments
  h_dim1 = *n;
  h_offset = 1 + h_dim1 * 1;
  h__ -= h_offset;
  --x;
  pl_dim1 = *npt;
  pl_offset = 1 + pl_dim1 * 1;
  pl -= pl_offset;
  xpt_dim1 = *npt;
  xpt_offset = 1 + xpt_dim1 * 1;
  xpt -= xpt_offset;
  --xbase;
  --xopt;
  --xnew;
  --pq;
  --g;
  --d__;
  --vlag;
  --w;
  --parscale;

  // Function Body
  one = 1.;
  two = 2.;
  zero = 0.;
  half = .5;
  tol = .01;
  nnp = *n + *n + 1;
  nptm = *npt - 1;
  nftest = _MAX(*maxfun,(unsigned int)1);

  // Initialization. NF is the number of function calculations so far.
  rho = *rhobeg;
  rhosq = rho * rho;
  *nf = 0;
  for (i = 1; i <= *n; ++i) {
    xbase[i] = x[i];
    for (k = 1; k <= *npt; ++k) {
      xpt[k + i * xpt_dim1] = zero;
    }
  }
  for (k = 1; k <= *npt; ++k) {
    for (j = 1; j <= nptm; ++j) {
      pl[k + j * pl_dim1] = zero;
    }
  }

  // The branch to label 120 obtains a new value of the objective function
  // and then there is a branch back to label 50, because the new function
  // value is needed to form the initial quadratic model. The least function
  // value so far and its index are noted below.
 L30:
  for (i = 1; i <= *n; ++i) {
    x[i] = xbase[i] + xpt[*nf + 1 + i * xpt_dim1];
  }
  goto L120;
 L50:
  if (*nf == 1) {
    fopt = f;
    kopt = *nf;
    fbase = f;
    j = 0;
    jswitch = -1;
    ih = *n;
  } else {
    if (f < fopt) {
      fopt = f;
      kopt = *nf;
    }
  }

  // Form the gradient and diagonal second derivatives of the initial
  // quadratic model and Lagrange functions.
  if (*nf <= nnp) {
    jswitch = -jswitch;
    if (jswitch > 0) {
      if (j >= 1) {
        ih += j;
        if (w[j] < zero) {
          d__[j] = (fsave + f - two * fbase) / rhosq;
          pq[j] = (fsave - f) / (two * rho);
          pl[ih * pl_dim1 + 1] = -two / rhosq;
          pl[*nf - 1 + j * pl_dim1] = half / rho;
          pl[*nf - 1 + ih * pl_dim1] = one / rhosq;
        } else {
          pq[j] = (fsave * 4. - fbase * 3. - f) / (two * rho);
          d__[j] = (fbase + f - two * fsave) / rhosq;
          pl[j * pl_dim1 + 1] = -1.5 / rho;
          pl[ih * pl_dim1 + 1] = one / rhosq;
          pl[*nf - 1 + j * pl_dim1] = two / rho;
          pl[*nf - 1 + ih * pl_dim1] = -two / rhosq;
        }
        pq[ih] = d__[j];
        pl[*nf + j * pl_dim1] = -half / rho;
        pl[*nf + ih * pl_dim1] = one / rhosq;
      }

      // Pick the shift from XBASE to the next initial interpolation point
      // that provides diagonal second derivatives.
      if (j < *n) {
        ++j;
        xpt[*nf + 1 + j * xpt_dim1] = rho;
      }
    } else {
      fsave = f;
      if (f < fbase) {
        w[j] = rho;
        xpt[*nf + 1 + j * xpt_dim1] = two * rho;
      } else {
        w[j] = -rho;
        xpt[*nf + 1 + j * xpt_dim1] = -rho;
      }
    }
    if (*nf < nnp) {
      goto L30;
    }
    // Form the off-diagonal second derivatives of the initial quadratic model.
    ih = *n;
    ip = 1;
    iq = 2;
  }
  ++ih;
  if (*nf > nnp) {
    temp = one / (w[ip] * w[iq]);
    tempa = f - fbase - w[ip] * pq[ip] - w[iq] * pq[iq];
    pq[ih] = (tempa - half * rhosq * (d__[ip] + d__[iq])) * temp;
    pl[ih * pl_dim1 + 1] = temp;
    iw = ip + ip;
    if (w[ip] < zero) {
      ++iw;
    }
    pl[iw + ih * pl_dim1] = -temp;
    iw = iq + iq;
    if (w[iq] < zero) {
      ++iw;
    }
    pl[iw + ih * pl_dim1] = -temp;
    pl[*nf + ih * pl_dim1] = temp;

    // Pick the shift from XBASE to the next initial interpolation point
    // that provides off-diagonal second derivatives.
    ++ip;
  }
  if (ip == iq) {
    ++ih;
    ip = 1;
    ++iq;
  }
  if (*nf < *npt) {
    xpt[*nf + 1 + ip * xpt_dim1] = w[ip];
    xpt[*nf + 1 + iq * xpt_dim1] = w[iq];
    goto L30;
  }

  // Set parameters to begin the iterations for the current RHO.
  sixthm = zero;
  delta = rho;
 L60:
  // Computing 2nd power
  d1 = two * rho;
  tworsq = d1 * d1;
  rhosq = rho * rho;

  // Form the gradient of the quadratic model at the trust region centre.
 L70:
  knew = 0;
  ih = *n;
  for (j = 1; j <= *n; ++j) {
    xopt[j] = xpt[kopt + j * xpt_dim1];
    g[j] = pq[j];
    for (i = 1; i <= j; ++i) {
      ++ih;
      g[i] += pq[ih] * xopt[j];
      if (i < j) {
        g[j] += pq[ih] * xopt[i];
      }
      h__[i + j * h_dim1] = pq[ih];
    }
  }

  // Generate the next trust region step and test its length. Set KNEW
  // to -1 if the purpose of the next F will be to improve conditioning,
  // and also calculate a lower bound on the Hessian term of the model Q.
  trstep(n, &g[1], &h__[h_offset], &delta, &tol, &d__[1], &w[1], &w[*n + 1],
         &w[(*n << 1) + 1], &w[*n * 3 + 1], &w[(*n << 2) + 1], &w[*n * 5 + 1], &evalue);
  temp = zero;
  // Computing 2nd power
  for (i = 1; i <= *n; ++i) {
    d1 = d__[i];
    temp += d1 * d1;
  }
  // Computing MIN
  d1 = delta, d2 = sqrt(temp);
  dnorm = _MIN(d1,d2);
  errtol = -one;
  if (dnorm < half * rho) {
    knew = -1;
    errtol = half * evalue * rho * rho;
    if (*nf <= *npt + 9) {
      errtol = zero;
    }
    goto L290;
  }

  // Calculate the next value of the objective function.
 L100:
  for (i = 1; i <= *n; ++i) {
    xnew[i] = xopt[i] + d__[i];
    x[i] = xbase[i] + xnew[i];
  }
 L120:
  if (*nf >= nftest) {
    if (*iprint >= 0) {
      char message[100];
      sprintf(message, "Return from optimizer because model function has been called MAX iterations = %d times", *maxfun);
      warning(message);
    }
    goto L420;
  }
  ++(*nf);

  // function evaluation and validation for new par
  for(i = 1; i <= *n; ++i) {
    x[i] /= parscale[i];
  }
  defineVar(install("..par.."), par, env);
  fval = eval(expr, env);
  if(!is_valid(fval)) {
    char message[100];
    sprintf(message, "Function returned NA or Inf at step %d.", *nf);
    error(message);
  }
  if(LENGTH(fval) != 1) {
    if(*nf == 1) warning("Function does not return a scalar");
    f = REAL(fval)[0];
  } else {
    f = REAL(fval)[0];
  }
  f /= *fnscale;
  // end function evaluation and validation

  if (*iprint >= 2 && (*nf % 10 == 0 || *nf == 1)) {
    char message[100];
    sprintf(message, "iteration %5d;  F = %10.3f;  rho = %10.3f\n", *nf, f, rho);
    printf("%s", message);
  }
  if (*nf <= *npt) {
    goto L50;
  }
  if (knew == -1) {
    goto L420;
  }

  // Use the quadratic model to predict the change in F due to the step D,
  // and find the values of the Lagrange functions at the new point.
  vquad = zero;
  ih = *n;
  for (j = 1; j <= *n; ++j) {
    w[j] = d__[j];
    vquad += w[j] * pq[j];
    for (i = 1; i <= j; ++i) {
      ++ih;
      w[ih] = d__[i] * xnew[j] + d__[j] * xopt[i];
      if (i == j) {
        w[ih] = half * w[ih];
      }
      vquad += w[ih] * pq[ih];
    }
  }
  for (k = 1; k <= *npt; ++k) {
    temp = zero;
    for (j = 1; j <= nptm; ++j) {
      temp += w[j] * pl[k + j * pl_dim1];
    }
    vlag[k] = temp;
  }
  vlag[kopt] += one;

  // Update SIXTHM, which is a lower bound on one sixth of the greatest
  // third derivative of F.
  diff = f - fopt - vquad;
  sum = zero;
  for (k = 1; k <= *npt; ++k) {
    temp = zero;
    // Computing 2nd power
    for (i = 1; i <= *n; ++i) {
      d1 = xpt[k + i * xpt_dim1] - xnew[i];
      temp += d1 * d1;
    }
    temp = sqrt(temp);
    sum += (d1 = temp * temp * temp * vlag[k], fabs(d1));
  }
  // Computing MAX
  d1 = sixthm, d2 = fabs(diff) / sum;
  sixthm = _MAX(d1,d2);

  // Update FOPT and XOPT if the new F is the least value of the objective
  // function so far. Then branch if D is not a trust region step.
  fsave = fopt;
  if (f < fopt) {
    fopt = f;
    for (i = 1; i <= *n; ++i) {
      // value x[i]=xbase[i]+xopt[i]
      xopt[i] = xnew[i];
    }
  }
  ksave = knew;
  if (knew > 0) {
    goto L240;
  }

  // Pick the next value of DELTA after a trust region step.
  if (vquad >= zero) {
    if (*iprint >= 0) {
      warning("Return from optimizer because a trust region step has failed to reduce");
    }
    goto L420;
  }
  ratio = (f - fsave) / vquad;
  if (ratio <= .1) {
    delta = half * dnorm;
  } else if (ratio <= .7) {
    // Computing MAX
    d1 = half * delta;
    delta = _MAX(d1,dnorm);
  } else {
    // Computing MAX
    d1 = delta, d2 = dnorm * 1.25, d1 = _MAX(d1,d2), d2 = dnorm + rho;
    delta = _MAX(d1,d2);
  }
  if (delta <= rho * 1.5) {
    delta = rho;
  }

  // Set KNEW to the index of the next interpolation point to be deleted.
  ktemp = 0;
  detrat = zero;
  if (f >= fsave) {
    ktemp = kopt;
    detrat = one;
  }
  for (k = 1; k <= *npt; ++k) {
    sum = zero;
    // Computing 2nd power
    for (i = 1; i <= *n; ++i) {
      d1 = xpt[k + i * xpt_dim1] - xopt[i];
      sum += d1 * d1;
    }
    temp = (d1 = vlag[k], fabs(d1));
    if (sum > rhosq) {
      d1 = sum / rhosq;
      temp *= pow(d1, c_b35);
    }
    if (temp > detrat && k != ktemp) {
      detrat = temp;
      ddknew = sum;
      knew = k;
    }
  }
  if (knew == 0) {
    goto L290;
  }

  // Replace the interpolation point that has index KNEW by the point XNEW,
  // and also update the Lagrange functions and the quadratic model.
 L240:
  for (i = 1; i <= *n; ++i) {
    xpt[knew + i * xpt_dim1] = xnew[i];
  }
  temp = one / vlag[knew];
  for (j = 1; j <= nptm; ++j) {
    pl[knew + j * pl_dim1] = temp * pl[knew + j * pl_dim1];
    pq[j] += diff * pl[knew + j * pl_dim1];
  }
  for (k = 1; k <= *npt; ++k) {
    if (knew>=0 && k != (unsigned int)knew) {
      temp = vlag[k];
      for (j = 1; j <= nptm; ++j) {
        pl[k + j * pl_dim1] -= temp * pl[knew + j * pl_dim1];
      }
    }
  }

  // Update KOPT if F is the least calculated value of the objective
  // function. Then branch for another trust region calculation. The
  // case KSAVE>0 indicates that a model step has just been taken.
  if (f < fsave) {
    kopt = knew;
    goto L70;
  }
  if (ksave > 0) {
    goto L70;
  }
  if (dnorm > two * rho) {
    goto L70;
  }
  if (ddknew > tworsq) {
    goto L70;
  }

  // Alternatively, find out if the interpolation points are close
  // enough to the best point so far.
 L290:
  for (k = 1; k <= *npt; ++k) {
    w[k] = zero;
    // Computing 2nd power
    for (i = 1; i <= *n; ++i) {
      d1 = xpt[k + i * xpt_dim1] - xopt[i];
      w[k] += d1 * d1;
    }
  }
 L310:
  knew = -1;
  distest = tworsq;
  for (k = 1; k <= *npt; ++k) {
    if (w[k] > distest) {
      knew = k;
      distest = w[k];
    }
  }

  // If a point is sufficiently far away, then set the gradient and Hessian
  // of its Lagrange function at the centre of the trust region, and find
  // half the sum of squares of components of the Hessian.
  if (knew > 0) {
    ih = *n;
    sumh = zero;
    for (j = 1; j <= *n; ++j) {
      g[j] = pl[knew + j * pl_dim1];
      for (i = 1; i <= j; ++i) {
        ++ih;
        temp = pl[knew + ih * pl_dim1];
        g[j] += temp * xopt[i];
        if (i < j) {
          g[i] += temp * xopt[j];
          sumh += temp * temp;
        }
        h__[i + j * h_dim1] = temp;
      }
      sumh += half * temp * temp;
    }

    // If ERRTOL is positive, test whether to replace the interpolation point
    // with index KNEW, using a bound on the maximum modulus of its Lagrange
    // function in the trust region.
    if (errtol > zero) {
      w[knew] = zero;
      sumg = zero;
      // Computing 2nd power
      for (i = 1; i <= *n; ++i) {
        d1 = g[i];
        sumg += d1 * d1;
      }
      estim = rho * (sqrt(sumg) + rho * sqrt(half * sumh));
      wmult = sixthm * pow(distest, c_b35);
      if (wmult * estim <= errtol) {
        goto L310;
      }
    }

    // If the KNEW-th point may be replaced, then pick a D that gives a large
    // value of the modulus of its Lagrange function within the trust region.
    // Here the vector XNEW is used as temporary working space.
    lagmax(n, &g[1], &h__[h_offset], &rho, &d__[1], &xnew[1], &vmax);
    if (errtol > zero) {
      if (wmult * vmax <= errtol) {
        goto L310;
      }
    }
    goto L100;
  }
  if (dnorm > rho) {
    goto L70;
  }

  // Prepare to reduce RHO by shifting XBASE to the best point so far,
  // and make the corresponding changes to the gradients of the Lagrange
  // functions and the quadratic model.
  if (rho > *rhoend) {
    ih = *n;
    for (j = 1; j <= *n; ++j) {
      xbase[j] += xopt[j];
      for (k = 1; k <= *npt; ++k) {
        xpt[k + j * xpt_dim1] -= xopt[j];
      }
      for (i = 1; i <= j; ++i) {
        ++ih;
        pq[i] += pq[ih] * xopt[j];
        if (i < j) {
          pq[j] += pq[ih] * xopt[i];
          for (k = 1; k <= *npt; ++k) {
            pl[k + j * pl_dim1] += pl[k + ih * pl_dim1] * xopt[i];
          }
        }
        for (k = 1; k <= *npt; ++k) {
          pl[k + i * pl_dim1] += pl[k + ih * pl_dim1] * xopt[j];
        }
      }
    }

    // Pick the next values of RHO and DELTA.
    delta = half * rho;
    ratio = rho / *rhoend;
    if (ratio <= 16.) {
      rho = *rhoend;
    } else if (ratio <= 250.) {
      rho = sqrt(ratio) * *rhoend;
    } else {
      rho *= .1;
    }
    delta = _MAX(delta,rho);
    if(*iprint >= 3) {
      char message[100];
      sprintf(message, "iteration %5d;  F = %10.3f;  rho = %10.3f***\n", *nf, f, rho);
      printf(message);
    }
    goto L60;
  }

  // Return from the calculation, after another Newton-Raphson step, if
  // it is too short to have been tried before.
  if (errtol >= zero) {
    goto L100;
  }
 L420:
  if (fopt <= f) {
    for (i = 1; i <= *n; ++i) {
      x[i] = xbase[i] + xopt[i];
    }
    f = fopt;
  }
  if(*iprint >= 3) {
    printf("*** rho reduction step\n");
  }
  if(*iprint >= 1) {
    char message[100];
    printf("\nFinished:\n");
    sprintf(message, "iteration %5d;  F = %10.3f;  rho = %10.3f\n", *nf, f, rho);
    printf(message);
  }

  return;
}

void lagmax(unsigned int *n,double *g,double *h__,double *rho,
            double *d__,double *v, double *vmax)
{
  // N is the number of variables of a quadratic objective function, Q say.
  // G is the gradient of Q at the origin.
  // H is the symmetric Hessian matrix of Q. Only the upper triangular and
  //   diagonal parts need be set.
  // RHO is the trust region radius, and has to be positive.
  // D will be set to the calculated vector of variables.
  // The array V will be used for working space.
  // VMAX will be set to |Q(0)-Q(D)|.

  // Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| .LEQ. RHO
  // requires of order N**3 operations, but sometimes it is adequate if
  // |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
  // subroutine provides such a solution in only of order N**2 operations,
  // where the claim of accuracy has been tested by numerical experiments.

  // System generated locals
  int h_dim1, h_offset;
  double d1, d2, d3, d4;

  // Local variables
  static double half, dlin, hmax, temp, vlin, wcos, zero, wsin, sumv;
  static unsigned int i, j, k;
  static double scale, tempa, tempb, tempc, tempd, ratio, gnorm, tempv, vnorm, dd, gd, gg, vv;
  static double halfrt, dhd, ghg, one, vhg, dsq, vhv, sum, whw, vhw, vmu, vsq, wsq;

  // Parameter adjustments
  h_dim1 = *n;
  h_offset = 1 + h_dim1 * 1;
  h__ -= h_offset;
  --g;
  --d__;
  --v;

  // Function Body
  half = .5;
  halfrt = sqrt(half);
  one = 1.;
  zero = 0.;

  // Pick V such that ||HV|| / ||V|| is large
  hmax = zero;
  for (i = 1; i <= *n; ++i) {
    sum = zero;
    for (j = 1; j <= *n; ++j) {
      h__[j + i * h_dim1] = h__[i + j * h_dim1];
      // Computing 2nd power
      d1 = h__[i + j * h_dim1];
      sum += d1 * d1;
    }
    if (sum > hmax) {
      hmax = sum;
      k = i;
    }
  }
  for (j = 1; j <= *n; ++j) {
    v[j] = h__[k + j * h_dim1];
  }

  // Set D to a vector in the subspace spanned by V and HV that maximizes
  // |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel
  // The vector that has the name D at label 60 used to be the vector W
  vsq = zero;
  vhv = zero;
  dsq = zero;
  for (i = 1; i <= *n; ++i) {
    // Computing 2nd power
    d1 = v[i];
    vsq += d1 * d1;
    d__[i] = zero;
    for (j = 1; j <= *n; ++j) {
      d__[i] += h__[i + j * h_dim1] * v[j];
    }
    vhv += v[i] * d__[i];
    // Computing 2nd power
    d1 = d__[i];
    dsq += d1 * d1;
  }
  if (vhv * vhv <= dsq * .9999 * vsq) {
    temp = vhv / vsq;
    wsq = zero;
    for (i = 1; i <= *n; ++i) {
      d__[i] -= temp * v[i];
      // Computing 2nd power
      d1 = d__[i];
      wsq += d1 * d1;
    }
    whw = zero;
    ratio = sqrt(wsq / vsq);
    for (i = 1; i <= *n; ++i) {
      temp = zero;
      for (j = 1; j <= *n; ++j) {
        temp += h__[i + j * h_dim1] * d__[j];
      }
      whw += temp * d__[i];
      v[i] = ratio * v[i];
    }
    vhv = ratio * ratio * vhv;
    vhw = ratio * wsq;
    temp = half * (whw - vhv);
    // Computing 2nd power
    d2 = temp;
    // Computing 2nd power
    d3 = vhw;
    d1 = sqrt(d2 * d2 + d3 * d3);
    d4 = whw + vhv;
    temp += d_sign(&d1, &d4);
    for (i = 1; i <= *n; ++i) {
      d__[i] = vhw * v[i] + temp * d__[i];
    }
  }

  // We now turn our attention to the subspace spanned by G and D. A multiple
  // of the current D is returned if that choice seems to be adequate.
  gg = zero;
  gd = zero;
  dd = zero;
  dhd = zero;
  for (i = 1; i <= *n; ++i) {
    // Computing 2nd power
    d1 = g[i];
    gg += d1 * d1;
    gd += g[i] * d__[i];
    // Computing 2nd power
    d1 = d__[i];
    dd += d1 * d1;
    sum = zero;
    for (j = 1; j <= *n; ++j) {
      sum += h__[i + j * h_dim1] * d__[j];
    }
    dhd += sum * d__[i];
  }
  temp = gd / gg;
  vv = zero;
  d1 = *rho / sqrt(dd);
  d2 = gd * dhd;
  scale = d_sign(&d1, &d2);
  for (i = 1; i <= *n; ++i) {
    v[i] = d__[i] - temp * g[i];
    // Computing 2nd power
    d1 = v[i];
    vv += d1 * d1;
    d__[i] = scale * d__[i];
  }
  gnorm = sqrt(gg);
  if (gnorm * dd <= *rho * .005 * fabs(dhd) || vv / dd <= 1e-4) {
    *vmax = (d1 = scale * (gd + half * scale * dhd), fabs(d1));
  } else {
    // G and V are now orthogonal in the subspace spanned by G and D. Hence
    // we generate an orthonormal basis of this subspace such that (D,HV) is
    // negligible or zero, where D and V will be the basis vectors.
    ghg = zero;
    vhg = zero;
    vhv = zero;
    for (i = 1; i <= *n; ++i) {
      sum = zero;
      sumv = zero;
      for (j = 1; j <= *n; ++j) {
        sum += h__[i + j * h_dim1] * g[j];
        sumv += h__[i + j * h_dim1] * v[j];
      }
      ghg += sum * g[i];
      vhg += sumv * g[i];
      vhv += sumv * v[i];
    }
    vnorm = sqrt(vv);
    ghg /= gg;
    vhg /= vnorm * gnorm;
    vhv /= vv;
    // Computing MAX
    d1 = fabs(ghg), d2 = fabs(vhv);
    if (fabs(vhg) <= _MAX(d1,d2) * .01) {
      vmu = ghg - vhv;
      wcos = one;
      wsin = zero;
    } else {
      temp = half * (ghg - vhv);
      // Computing 2nd power
      d2 = temp;
      // Computing 2nd power
      d3 = vhg;
      d1 = sqrt(d2 * d2 + d3 * d3);
      vmu = temp + d_sign(&d1, &temp);
      // Computing 2nd power
      d1 = vmu;
      // Computing 2nd power
      d2 = vhg;
      temp = sqrt(d1 * d1 + d2 * d2);
      wcos = vmu / temp;
      wsin = vhg / temp;
    }
    tempa = wcos / gnorm;
    tempb = wsin / vnorm;
    tempc = wcos / vnorm;
    tempd = wsin / gnorm;
    for (i = 1; i <= *n; ++i) {
      d__[i] = tempa * g[i] + tempb * v[i];
      v[i] = tempc * v[i] - tempd * g[i];
    }

    // The final D is a multiple of the current D, V, D+V or D-V. We make the
    // choice from these possibilities that is optimal.
    dlin = wcos * gnorm / *rho;
    vlin = -wsin * gnorm / *rho;
    tempa = fabs(dlin) + half * (d1 = vmu + vhv, fabs(d1));
    tempb = fabs(vlin) + half * (d1 = ghg - vmu, fabs(d1));
    tempc = halfrt * (fabs(dlin) + fabs(vlin)) + (d1 = ghg + vhv, fabs(d1)) * .25;
    if (tempa >= tempb && tempa >= tempc) {
      d1 = dlin * (vmu + vhv);
      tempd = d_sign(rho, &d1);
      tempv = zero;
    } else if (tempb >= tempc) {
      tempd = zero;
      d1 = vlin * (ghg - vmu);
      tempv = d_sign(rho, &d1);
    } else {
      d1 = halfrt * *rho;
      d2 = dlin * (ghg + vhv);
      tempd = d_sign(&d1, &d2);
      d1 = halfrt * *rho;
      d2 = vlin * (ghg + vhv);
      tempv = d_sign(&d1, &d2);
    }
    for (i = 1; i <= *n; ++i) {
      d__[i] = tempd * d__[i] + tempv * v[i];
    }
    // Computing MAX
    d1 = _MAX(tempa,tempb);
    *vmax = *rho * *rho * _MAX(d1,tempc);
  }

  return;
}

void trstep(unsigned int *n,double *g,double *h__,
            double *delta,double *tol,double *d__,
            double *gg,double *td,double *tn,
            double *w,double *piv,double *z__,double *evalue)
{
  // N is the number of variables of a quadratic objective function, Q say.
  // G is the gradient of Q at the origin.
  // H is the Hessian matrix of Q. Only the upper triangular and diagonal
  //   parts need be set. The lower triangular part is used to store the
  //   elements of a Householder similarity transformation.
  // DELTA is the trust region radius, and has to be positive.
  // TOL is the value of a tolerance from the open interval (0,1).
  // D will be set to the calculated vector of variables.
  // The arrays GG, TD, TN, W, PIV and Z will be used for working space.
  // EVALUE will be set to the least eigenvalue of H if and only if D is a
  // Newton-Raphson step. Then EVALUE will be positive, but otherwise it
  // will be set to zero.

  // Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| .LEQ. DELTA,
  // and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
  // We take the view that any D is acceptable if it has the properties
  //         ||D|| .LEQ. DELTA  and  ACTRED .LEQ. (1-TOL)*MAXRED.

  // The calculation of D is done by the method of Section 2 of the paper
  // by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
  // after transforming H to tridiagonal form.

  // System generated locals
  int h_dim1, h_offset;
  double d1, d2, d3;

  // Local variables
  static double phil, parl;
  static int ksav;
  static double temp, phiu, paru, zero, wwsq;
  static unsigned int i, j, k;
  static double scale;
  static int iterc;
  static double tempa, delsq, tempb;
  static int ksave;
  static double tdmin, shift, dnorm, gnorm, hnorm, slope, pivot;
  static unsigned int jp, nm, kp;
  static double posdef, wz, shfmin, shfmax, pivksv, dhd, gam, dtg, phi, one, par, dsq;
  static int kpp;
  static double gsq, dtz, sum, two, wsq, zsq, parlest, paruest;

  // Parameter adjustments
  h_dim1 = *n;
  h_offset = 1 + h_dim1 * 1;
  h__ -= h_offset;
  --g;
  --d__;
  --gg;
  --td;
  --tn;
  --w;
  --piv;
  --z__;

  // Function Body
  one = 1.;
  two = 2.;
  zero = 0.;
  delsq = *delta * *delta;
  *evalue = zero;
  for (i = 1; i <= *n; ++i) {
    d__[i] = zero;
    td[i] = h__[i + i * h_dim1];
    for (j = 1; j <= i; ++j) {
      h__[i + j * h_dim1] = h__[j + i * h_dim1];
    }
  }

  // Apply Householder transformations to obtain a tridiagonal matrix that
  // is similar to H, and put the elements of the Householder vectors in
  // the lower triangular part of H. Further, TD and TN will contain the
  // diagonal and other nonzero elements of the tridiagonal matrix.
  nm = *n - 1;
  for (k = 1; k <= nm; ++k) {
    kp = k + 1;
    sum = zero;
    if (kp < *n) {
      kpp = kp + 1;
      for (i = kpp; i <= *n; ++i) {
        // Computing 2nd power
        d1 = h__[i + k * h_dim1];
        sum += d1 * d1;
      }
    }
    if (sum == zero) {
      tn[k] = h__[kp + k * h_dim1];
      h__[kp + k * h_dim1] = zero;
    } else {
      temp = h__[kp + k * h_dim1];
      d1 = sqrt(sum + temp * temp);
      tn[k] = d_sign(&d1, &temp);
      h__[kp + k * h_dim1] = -sum / (temp + tn[k]);
      // Computing 2nd power
      d1 = h__[kp + k * h_dim1];
      temp = sqrt(two / (sum + d1 * d1));
      for (i = kp; i <= *n; ++i) {
        w[i] = temp * h__[i + k * h_dim1];
        h__[i + k * h_dim1] = w[i];
        z__[i] = td[i] * w[i];
      }
      wz = zero;
      for (j = kp; j <= nm; ++j) {
        jp = j + 1;
        for (i = jp; i <= *n; ++i) {
          z__[i] += h__[i + j * h_dim1] * w[j];
          z__[j] += h__[i + j * h_dim1] * w[i];
        }
        wz += w[j] * z__[j];
      }
      wz += w[*n] * z__[*n];
      for (j = kp; j <= *n; ++j) {
        td[j] += w[j] * (wz * w[j] - two * z__[j]);
        if (j < *n) {
          jp = j + 1;
          for (i = jp; i <= *n; ++i) {
            h__[i + j * h_dim1] = h__[i + j * h_dim1] - w[i] * z__[j] - w[j] * (z__[i] - wz * w[i]);
          }
        }
      }
    }
  }

  // Form GG by applying the similarity transformation to G.
  gsq = zero;
  for (i = 1; i <= *n; ++i) {
    gg[i] = g[i];
    // Computing 2nd power
    d1 = g[i];
    gsq += d1 * d1;
  }
  gnorm = sqrt(gsq);
  for (k = 1; k <= nm; ++k) {
    kp = k + 1;
    sum = zero;
    for (i = kp; i <= *n; ++i) {
      sum += gg[i] * h__[i + k * h_dim1];
    }
    for (i = kp; i <= *n; ++i) {
      gg[i] -= sum * h__[i + k * h_dim1];
    }
  }

  // Begin the trust region calculation with a tridiagonal matrix by
  // calculating the norm of H. Then treat the case when H is zero.
  hnorm = fabs(td[1]) + fabs(tn[1]);
  tdmin = td[1];
  tn[*n] = zero;
  for (i = 2; i <= *n; ++i) {
    temp = (d1 = tn[i - 1], fabs(d1)) + (d2 = td[i], fabs(d2)) + (d3 = tn[i], fabs(d3));
    hnorm = _MAX(hnorm,temp);
    // Computing MIN
    d1 = tdmin, d2 = td[i];
    tdmin = _MIN(d1,d2);
  }
  if (hnorm == zero) {
    if (gnorm == zero) {
      goto L400;
    }
    scale = *delta / gnorm;
    for (i = 1; i <= *n; ++i) {
      d__[i] = -scale * gg[i];
    }
    goto L370;
  }

  // Set the initial values of PAR and its bounds.
  // Computing MAX
  d1 = zero, d2 = -tdmin, d1 = _MAX(d1,d2), d2 = gnorm / *delta  - hnorm;
  parl = _MAX(d1,d2);
  parlest = parl;
  par = parl;
  paru = zero;
  paruest = zero;
  posdef = zero;
  iterc = 0;

  // Calculate the pivots of the Cholesky factorization of (H+PAR*I).
 L140:
  ++iterc;
  ksav = 0;
  piv[1] = td[1] + par;
  k = 1;
 L150:
  if (piv[k] > zero) {
    // Computing 2nd power
    d1 = tn[k];
    piv[k + 1] = td[k + 1] + par - d1 * d1 / piv[k];
  } else {
    if (piv[k] < zero || tn[k] != zero) {
      goto L160;
    }
    ksav = k;
    piv[k + 1] = td[k + 1] + par;
  }
  ++k;
  if (k < *n) {
    goto L150;
  }
  if (piv[k] < zero) {
    goto L160;
  }
  if (piv[k] == zero) {
    ksav = k;
  }

  // Branch if all the pivots are positive, allowing for the case when
  // G is zero.
  if (ksav == 0 && gsq > zero) {
    goto L230;
  }
  if (gsq == zero) {
    if (par == zero) {
      goto L370;
    }
    paru = par;
    paruest = par;
    if (ksav == 0) {
      goto L190;
    }
  }
  k = ksav;

  // Set D to a direction of nonpositive curvature of the given tridiagonal
  // matrix, and thus revise PARLEST.
 L160:
  d__[k] = one;
  if ((d1 = tn[k], fabs(d1)) <= (d2 = piv[k], fabs(d2))) {
    dsq = one;
    dhd = piv[k];
  } else {
    temp = td[k + 1] + par;
    if (temp <= (d1 = piv[k], fabs(d1))) {
      d1 = -tn[k];
      d__[k + 1] = d_sign(&one, &d1);
      dhd = piv[k] + temp - two * (d1 = tn[k], fabs(d1));
    } else {
      d__[k + 1] = -tn[k] / temp;
      dhd = piv[k] + tn[k] * d__[k + 1];
    }
    // Computing 2nd power
    d1 = d__[k + 1];
    dsq = one + d1 * d1;
  }
 L170:
  if (k > 1) {
    --k;
    if (tn[k] != zero) {
      d__[k] = -tn[k] * d__[k + 1] / piv[k];
      // Computing 2nd power
      d1 = d__[k];
      dsq += d1 * d1;
      goto L170;
    }
    for (i = 1; i <= k; ++i) {
      d__[i] = zero;
    }
  }
  parl = par;
  parlest = par - dhd / dsq;

  // Terminate with D set to a multiple of the current D if the following
  // test suggests that it suitable to do so.
 L190:
  temp = paruest;
  if (gsq == zero) {
    temp *= one - *tol;
  }
  if (paruest > zero && parlest >= temp) {
    dtg = zero;
    for (i = 1; i <= *n; ++i) {
      dtg += d__[i] * gg[i];
    }
    d1 = *delta / sqrt(dsq);
    scale = -d_sign(&d1, &dtg);
    for (i = 1; i <= *n; ++i) {
      d__[i] = scale * d__[i];
    }
    goto L370;
  }

  // Pick the value of PAR for the next iteration.
 L220:
  if (paru == zero) {
    par = two * parlest + gnorm / *delta;
  } else {
    par = (parl + paru) * .5;
    par = _MAX(par,parlest);
  }
  if (paruest > zero) {
    par = _MIN(par,paruest);
  }
  goto L140;

  // Calculate D for the current PAR in the positive definite case.
 L230:
  w[1] = -gg[1] / piv[1];
  for (i = 2; i <= *n; ++i) {
    w[i] = (-gg[i] - tn[i - 1] * w[i - 1]) / piv[i];
  }
  d__[*n] = w[*n];
  for (i = nm; i >= 1; --i) {
    d__[i] = w[i] - tn[i] * d__[i + 1] / piv[i];
  }

  // Branch if a Newton-Raphson step is acceptable.
  dsq = zero;
  wsq = zero;
  for (i = 1; i <= *n; ++i) {
    // Computing 2nd power
    d1 = d__[i];
    dsq += d1 * d1;
    // Computing 2nd power
    d1 = w[i];
    wsq += piv[i] * (d1 * d1);
  }
  if (par == zero && dsq <= delsq) {
    goto L320;
  }

  // Make the usual test for acceptability of a full trust region step.
  dnorm = sqrt(dsq);
  phi = one / dnorm - one / *delta;
  temp = *tol * (one + par * dsq / wsq) - dsq * phi * phi;
  if (temp >= zero) {
    scale = *delta / dnorm;
    for (i = 1; i <= *n; ++i) {
      d__[i] = scale * d__[i];
    }
    goto L370;
  }
  if (iterc >= 2 && par <= parl) {
    goto L370;
  }
  if (paru > zero && par >= paru) {
    goto L370;
  }

  // Complete the iteration when PHI is negative.
  if (phi < zero) {
    parlest = par;
    if (posdef == one) {
      if (phi <= phil) {
        goto L370;
      }
      slope = (phi - phil) / (par - parl);
      parlest = par - phi / slope;
    }
    slope = one / gnorm;
    if (paru > zero) {
      slope = (phiu - phi) / (paru - par);
    }
    temp = par - phi / slope;
    if (paruest > zero) {
      temp = _MIN(temp,paruest);
    }
    paruest = temp;
    posdef = one;
    parl = par;
    phil = phi;
    goto L220;
  }

  // If required, calculate Z for the alternative test for convergence.
  if (posdef == zero) {
    w[1] = one / piv[1];
    for (i = 2; i <= *n; ++i) {
      temp = -tn[i - 1] * w[i - 1];
      w[i] = (d_sign(&one, &temp) + temp) / piv[i];
    }
    z__[*n] = w[*n];
    for (i = nm; i >= 1; --i) {
      z__[i] = w[i] - tn[i] * z__[i + 1] / piv[i];
    }
    wwsq = zero;
    zsq = zero;
    dtz = zero;
    for (i = 1; i <= *n; ++i) {
      // Computing 2nd power
      d1 = w[i];
      wwsq += piv[i] * (d1 * d1);
      // Computing 2nd power
      d1 = z__[i];
      zsq += d1 * d1;
      dtz += d__[i] * z__[i];
    }

    // Apply the alternative test for convergence.
    tempa = (d1 = delsq - dsq, fabs(d1));
    tempb = sqrt(dtz * dtz + tempa * zsq);
    gam = tempa / (d_sign(&tempb, &dtz) + dtz);
    temp = *tol * (wsq + par * delsq) - gam * gam * wwsq;
    if (temp >= zero) {
      for (i = 1; i <= *n; ++i) {
        d__[i] += gam * z__[i];
      }
      goto L370;
    }
    // Computing MAX
    d1 = parlest, d2 = par - wwsq / zsq;
    parlest = _MAX(d1,d2);
  }

  // Complete the iteration when PHI is positive.
  slope = one / gnorm;
  if (paru > zero) {
    if (phi >= phiu) {
      goto L370;
    }
    slope = (phiu - phi) / (paru - par);
  }
  // Computing MAX
  d1 = parlest, d2 = par - phi / slope;
  parlest = _MAX(d1,d2);
  paruest = par;
  if (posdef == one) {
    slope = (phi - phil) / (par - parl);
    paruest = par - phi / slope;
  }
  paru = par;
  phiu = phi;
  goto L220;

  // Set EVALUE to the least eigenvalue of the second derivative matrix if
  // D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.
 L320:
  shfmin = zero;
  pivot = td[1];
  shfmax = pivot;
  for (k = 2; k <= *n; ++k) {
    // Computing 2nd power
    d1 = tn[k - 1];
    pivot = td[k] - d1 * d1 / pivot;
    shfmax = _MIN(shfmax,pivot);
  }

  // Find EVALUE by a bisection method, but occasionally SHFMAX may be
  // adjusted by the rule of false position.
  ksave = 0;
 L340:
  shift = (shfmin + shfmax) * .5;
  k = 1;
  temp = td[1] - shift;
 L350:
  if (temp > zero) {
    piv[k] = temp;
    if (k < *n) {
      // Computing 2nd power
      d1 = tn[k];
      temp = td[k + 1] - shift - d1 * d1 / temp;
      ++k;
      goto L350;
    }
    shfmin = shift;
  } else {
    if (ksave>=0 && k < (unsigned int)ksave) {
      goto L360;
    }
    if (ksave>=0 && k == (unsigned int)ksave) {
      if (pivksv == zero) {
        goto L360;
      }
      if (piv[k] - temp < temp - pivksv) {
        pivksv = temp;
        shfmax = shift;
      } else {
        pivksv = zero;
        shfmax = (shift * piv[k] - shfmin * temp) / (piv[k] - temp);
      }
    } else {
      ksave = k;
      pivksv = temp;
      shfmax = shift;
    }
  }
  if (shfmin <= shfmax * .99) {
    goto L340;
  }
 L360:
  *evalue = shfmin;

  // Apply the inverse Householder transformations to D.
 L370:
  nm = *n - 1;
  for (k = nm; k >= 1; --k) {
    kp = k + 1;
    sum = zero;
    for (i = kp; i <= *n; ++i) {
      sum += d__[i] * h__[i + k * h_dim1];
    }
    for (i = kp; i <= *n; ++i) {
      d__[i] -= sum * h__[i + k * h_dim1];
    }
  }

  // Return from the subroutine.
 L400:
  return;
}

double d_sign(double *a,double *b)
{
  double x;
  x = (*a >= 0 ? *a : - *a);

  return(*b >= 0 ? x : -x);
}

void getHessian(double *pq)
{
  // Hessian matrix [(a,b),(c,d)] upper triangle (column major in Powell, row major in IMSL)
  // Hessian is in the quadratic model and triangel symmetric (n*(n+1)/2 data points)
  //   unsigned int npar((*m_model).m_npar);
  //   double *Ha=new double[npar*npar];
  //   (*m_model).m_invH=new double[npar*npar];

  unsigned int i, j, k;
  // upper triangle + diagonal
  for (j=0,k=0; j<npar; ++j) {
    for (i=0; i<=j; ++i) {
      Ha[j*npar+i]=pq[k++];
    }
  }
  // lower triangle
  for (j=1; j<npar; ++j) {
    for (i=0; i<j; ++i) {
      Ha[i*npar+j]=Ha[j*npar+i];
    }
  }

  return;
}
// End of Unconstrained Optimization by Quadratic Approximation
