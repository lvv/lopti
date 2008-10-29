#include <stdio.h>
#include <math.h>
#include "machine.h"

#ifdef MSDOS
static double e[N];	/* to save stack space */
#endif

/* singular value decomposition */

minfit(n, eps, tol, ab, q)
int n;
double eps, tol, ab[N][N], q[N];
{
   int l, kt, l2, i, j, k;
   double c, f, g, h, s, x, y, z;
#ifndef MSDOS
   double e[N];		/* plenty of stack on a vax */
#endif

   /* householder's reduction to bidiagonal form */
   x = g = 0.0;
   for (i=0; i<n; i++) {
       e[i] = g; s = 0.0; l = i+1;
       for (j=i; j<n; j++)
	   s += ab[j][i] * ab[j][i];
       if (s < tol) {
	  g = 0.0;
       }
       else {
	  f = ab[i][i];
          if (f < 0.0) 
	     g = sqrt(s);
	  else
	     g = -sqrt(s);
	  h = f*g - s; ab[i][i] = f - g;
	  for (j=l; j<n; j++) {
	      f = 0.0;
	      for (k=i; k<n; k++)
		  f += ab[k][i] * ab[k][j];
	      f /= h;
	      for (k=i; k<n; k++)
		  ab[k][j] += f * ab[k][i];
	  }
       }
       q[i] = g; s = 0.0;
       if (i < n)
	  for (j=l; j<n; j++)
	      s += ab[i][j] * ab[i][j];
       if (s < tol) {
	  g = 0.0;
       }
       else {
	  f = ab[i][i+1];
	  if (f < 0.0)
	     g = sqrt(s);
	  else 
	     g = - sqrt(s);
	  h = f*g - s; ab[i][i+1] = f - g;
	  for (j=l; j<n; j++)
	      e[j] = ab[i][j]/h;
	  for (j=l; j<n; j++) {
	      s = 0;
	      for (k=l; k<n; k++) s += ab[j][k]*ab[i][k];
	      for (k=l; k<n; k++) ab[j][k] += s * e[k];
	  }
       }
       y = fabs(q[i]) + fabs(e[i]);
       if (y > x) x = y;
   }
   /* accumulation of right hand transformations */
   for (i=n-1; i >= 0; i--) {
       if (g != 0.0) {
          h = ab[i][i+1]*g;
	  for (j=l; j<n; j++) ab[j][i] = ab[i][j] / h;
	  for (j=l; j<n; j++) {
              s = 0.0;
	      for (k=l; k<n; k++) s += ab[i][k] * ab[k][j];
	      for (k=l; k<n; k++) ab[k][j] += s * ab[k][i];
	  }
       }
       for (j=l; j<n; j++)
           ab[i][j] = ab[j][i] = 0.0;
       ab[i][i] = 1.0; g = e[i]; l = i;
   }
   /* diagonalization to bidiagonal form */
   eps *= x;
   for (k=n-1; k>= 0; k--) {
       kt = 0;
TestFsplitting:
       if (++kt > 30) {
          e[k] = 0.0;
	  fprintf(stderr, "\n+++ qr failed\n");
       }
       for (l2=k; l2>=0; l2--) {
           l = l2;
	   if (fabs(e[l]) <= eps)
	      goto TestFconvergence;
	   if (fabs(q[l-1]) <= eps)
   	      break;	/* goto Cancellation; */
       }
Cancellation:
       c = 0.0; s = 1.0;
       for (i=l; i<=k; i++) {
           f = s * e[i]; e[i] *= c;
	   if (fabs(f) <= eps)
	      goto TestFconvergence;
	   g = q[i];
   	   if (fabs(f) < fabs(g)) {
	      double fg = f/g;
	      h = fabs(g)*sqrt(1.0+fg*fg);
	   }
	   else {
	      double gf = g/f;
	      h = (f!=0.0 ? fabs(f)*sqrt(1.0+gf*gf) : 0.0);
	   }
	   q[i] = h;
	   if (h == 0.0) { h = 1.0; g = 1.0; }
	   c = g/h; s = -f/h;
       }
TestFconvergence:
       z = q[k];
       if (l == k)
          goto Convergence;
       /* shift from bottom 2x2 minor */
       x = q[l]; y = q[k-l]; g = e[k-1]; h = e[k];
       f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.0*h*y);
       g = sqrt(f*f+1.0);
       if (f <= 0.0)
          f = ((x-z)*(x+z) + h*(y/(f-g)-h))/x;
       else
          f = ((x-z)*(x+z) + h*(y/(f+g)-h))/x;
       /* next qr transformation */
       s = c = 1.0;
       for (i=l+1; i<=k; i++) {
           g = e[i]; y = q[i]; h = s*g; g *= c;
	   if (fabs(f) < fabs(h)) {
	      double fh = f/h;
	      z = fabs(h) * sqrt(1.0 + fh*fh);
	   }
	   else {
	      double hf = h/f;
	      z = (f!=0.0 ? fabs(f)*sqrt(1.0+hf*hf) : 0.0);
	   }
	   e[i-1] = z;
	   if (z == 0.0) 
 	      f = z = 1.0;
	   c = f/z; s = h/z;
	   f = x*c + g*s; g = - x*s + g*c; h = y*s;
	   y *= c;
	   for (j=0; j<n; j++) {
	       x = ab[j][i-1]; z = ab[j][i];
	       ab[j][i-1] = x*c + z*s;
	       ab[j][i] = - x*s + z*c;
	   }
	   if (fabs(f) < fabs(h)) {
	      double fh = f/h;
	      z = fabs(h) * sqrt(1.0 + fh*fh);
	   }
	   else {
	      double hf = h/f;
	      z = (f!=0.0 ? fabs(f)*sqrt(1.0+hf*hf) : 0.0);
	   }
           q[i-1] = z;
	   if (z == 0.0) z = f = 1.0;
	   c = f/z; s = h/z;
	   f = c*g + s*y; x = - s*g + c*y;
       }
       e[l] = 0.0; e[k] = f; q[k] = x;
       goto TestFsplitting;
Convergence:
       if (z < 0.0) {
          q[k] = - z;
	  for (j=0; j<n; j++) ab[j][k] = - ab[j][k];
       }
   }
}
