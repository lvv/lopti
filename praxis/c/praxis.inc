(****************************************************************************)
(*              P R O C E D U R E      P R A X I S                          *)
(*                                                                          *)
(*   P R A X I S  is for the minimization of a function in several          *)
(*   variables. The algorithm used is a modification of a conjugate         *)
(*   gradient method developed by Powell. Changes are due to Brent,         *)
(*   who gives an ALGOL W program.                                          *)
(*   Users who are interested in more of the details should read            *)
(*          - Powell, M. J. D., 1962. An efficient method for finding       *)
(*            the minimum of a function of several variables without        *)
(*            calculating derivatives, Computer Journal 7, 155-162.         *)
(*          - Brent, R. P., 1973. Algorithms for minimization without       *)
(*            derivatives. Prentice Hall, Englewood Cliffs.                 *)
(*   If you have any further comments, questions, etc. please feel free     *)
(*   and contact me.                                                        *)
(*                                                                          *)
(*                           Karl Gegenfurtner     02/17/86                 *)
(*                                                 TURBO - P Version 1.0    *)
(****************************************************************************)
(*  The use of PRAXIS is fairly simple. There are only two parameters:      *)
(*      - X is of type vector and contains on input an initial guess        *)
(*        of the solution. on output X holds the final solution of the      *)
(*        system of equations.                                              *)
(*      - N is of type integer and gives the number of unknown parameters   *)
(*  The result of PRAXIS is the least calculated value of the function F    *)
(*  F is one of the global parameters used by PRAXIS:                       *)
(*    - F(X,N) is declared forward and is the function to be minimized      *)
(*  All other globals are optional, i.e. you can change them or else        *)
(*  PRAXIS will use some default values, which are adequate for most        *)
(*  problems under consideration.                                           *)
(*    - Prin controls the printout from PRAXIS                              *)
(*           0:  no printout at all                                         *)
(*           1:  only initial and final values                              *)
(*           2:  detailed map of the minimization process                   *)
(*           3:  also prints eigenvalues and eigenvectors of the            *)
(*               direction matices in use (for insiders only).              *)
(*    - T is the tolerance for the precision of the solution                *)
(*      PRAXIS returns if the criterion                                     *)
(*             2 * ||x(k)-x(k-1)|| <= Sqrt(MachEps) * ||x(k)|| + T          *)
(*             is fulfilled more than KTM times, where                      *)
(*    - KTM is another global parameter. It's default value is 1 and        *)
(*      a value of 4 leads to a very(!) cautious stopping criterion         *)
(*    - MachEps is the relative machine precision and is                    *)
(*      1.0 E-15 using an 8087 and TURBO-87 and                             *)
(*      1.0 E-06 without 8087.                                              *)
(*    - H is a steplength parameter and shoul be set to the expected        *)
(*      distance to the solution. An exceptional large or small value       *)
(*      of H leads to slower convergence on the first few iterations.       *)
(*    - Scbd is a scaling parameter and should be set to about 10.          *)
(*      The default is 1 and with that value no scaling is done at all      *)
(*      It's only necessary when the parameters are scaled very different   *)
(*    - IllC is a Boolean variable and should be set to TRUE if the         *)
(*      problem is known to be illconditioned.                              *)
(****************************************************************************)


CONST MaxPar = 20;

TYPE  Vector = ARRAY[1..MaxPar] OF REAL;
      Matrix = ARRAY[1..MaxPar] OF Vector;
      PrString = STRING[80];

CONST Macheps: REAL = 1.0E-08;  (* set to 1.0E-13 for TURBO-87 *)
      H:       REAL = 1.0E+00;
      T:       REAL = 1.0E-05;
      Prin:    INTEGER = 2;
      IllC:    BOOLEAN = FALSE;
      Scbd:    REAL = 1;
      Ktm:     INTEGER = 2;


FUNCTION F(VAR x:Vector;n:INTEGER):REAL; FORWARD;

FUNCTION Power(a,b: REAL):REAL;
BEGIN    Power := Exp(b*Ln(a));
END;     (* Power *)


PROCEDURE MinFit(N:INTEGER;Eps,Tol:REAL;VAR ab:Matrix;VAR q:Vector);
LABEL     TestFsplitting,
          TestFconvergence,
          Convergence,
          Cancellation;
VAR       l, kt,
          l2,
          i, j, k: INTEGER;
          c, f, g,
          h, s, x,
          y, z:    REAL;
          e:       Vector;
BEGIN     (* Householders reduction to bidiagonal form *)
          x:= 0; g:= 0;
          FOR i:= 1 TO N DO
          BEGIN e[i]:= g; s:= 0; l:= i+1;
                FOR j:= i TO N DO
                    s:= s + Sqr(ab[j,i]);
                IF s<Tol THEN g:= 0 ELSE
                BEGIN f:= ab[i,i];
                      IF f<0
                         THEN g:= Sqrt(s)
                         ELSE g:= - Sqrt(s);
                      h:= f*g-s;  ab[i,i]:= f - g;
                      FOR j:= l TO N DO
                      BEGIN f:= 0;
                            FOR k:= i TO N DO
                                f:= f + ab[k,i]*ab[k,j];
                            f:= f/h;
                            FOR k:= i TO N DO
                                ab[k,j]:= ab[k,j] + f*ab[k,i];
                      END; (* j *)
                END; (* IF *)
                q[i]:= g; s:= 0;
                IF i<=N
                   THEN FOR j:= l TO N DO
                            s:= s + Sqr(ab[i,j]);
                IF s<Tol THEN g:= 0 ELSE
                BEGIN f:= ab[i,i+1];
                      IF f<0
                         THEN g:= Sqrt(s)
                         ELSE g:= - Sqrt(s);
                      h:= f*g-s;  ab[i,i+1]:= f-g;
                      FOR j:= l TO N DO e[j]:= ab[i,j]/h;
                      FOR j:= l TO N DO
                      BEGIN s:= 0;
                            FOR k:= l TO N DO s:= s + ab[j,k]*ab[i,k];
                            FOR k:= l TO N DO ab[j,k]:= ab[j,k] + s*e[k];
                      END; (* J *)
                END; (* IF *)
                y:= Abs(q[i])+Abs(e[i]);
                IF y > x THEN x:= y;
          END; (* I *)
          (* accumulation of right hand transformations *)
          FOR i:= N DOWNTO 1 DO
          BEGIN IF g<>0.0 THEN
                BEGIN h:= ab[i,i+1]*g;
                      FOR j:= l TO N DO ab[j,i]:= ab[i,j]/h;
                      FOR j:= l TO N DO
                      BEGIN s:= 0;
                            FOR k:= l TO N DO s:= s + ab[i,k]*ab[k,j];
                            FOR k:= l TO N DO ab[k,j]:= ab[k,j] + s*ab[k,i];
                      END; (* J *)
                END; (* IF *)
                FOR j:= l TO N DO
                BEGIN ab[j,i]:= 0;
                      ab[i,j]:= 0;
                END;
                ab[i,i]:= 1; g:= e[i]; l:= i;
          END; (* I *)
          (* diagonalization to bidiagonal form *)
          eps:= eps*x;
          FOR k:= N DOWNTO 1 DO
          BEGIN kt:= 0;
TestFsplitting: kt:= kt + 1;
                IF kt > 30 THEN
                BEGIN e[k]:= 0;
                      WriteLn('+++ QR failed');
                END;
                FOR l2:= k DOWNTO 1 DO
                BEGIN l:= l2;
                      IF Abs(e[l])<=Eps THEN Goto TestFconvergence;
                      IF Abs(q[l-1])<=Eps THEN Goto Cancellation;
                END; (* L2 *)
Cancellation:   c:= 0; s:= 1;
                FOR i:= l TO k DO
                BEGIN f:= s*e[i]; e[i]:= c*e[i];
                      IF Abs(f)<=Eps THEN Goto TestFconvergence;
                      g:= q[i];
                      IF Abs(f) < Abs(g)
                         THEN h:= Abs(g)*Sqrt(1+sqr(f/g))
                         ELSE IF f <> 0
                                 THEN h:= Abs(f)*Sqrt(1+Sqr(g/f))
                                 ELSE h:= 0;
                      q[i]:= h;
                      IF h = 0 THEN
                      BEGIN h:= 1;
                            g:= 1;
                      END;
                      c:= g/h; s:= -f/h;
                END; (* I *)
TestFconvergence:
                z:= q[k];
                IF l=k THEN Goto convergence;
                (* shift from bottom 2*2 minor *)
                x:= q[l]; y:= q[k-1]; g:= e[k-1];  h:= e[k];
                f:= ((y-z)*(y+z) + (g-h)*(g+h))/(2*h*y);
                g:= Sqrt(Sqr(f)+1);
                IF f<=0
                   THEN f:= ((x-z)*(x+z)+h*(y/(f-g)-h))/x
                   ELSE f:= ((x-z)*(x+z)+h*(y/(f+g)-h))/x;
                (* next QR transformation *)
                s:= 1; c:= 1;
                FOR i:= l+1 TO k DO
                BEGIN g:= e[i]; y:= q[i]; h:= s*g; g:= g*c;
                      IF Abs(f)<Abs(h)
                         THEN z:= Abs(h)*Sqrt(1+Sqr(f/h))
                         ELSE IF f<>0
                                 THEN z:= Abs(f)*Sqrt(1+Sqr(h/f))
                                 ELSE z:= 0;
                      e[i-1]:= z;
                      IF z=0 THEN
                      BEGIN f:= 1;
                            z:= 1;
                      END;
                      c:= f/z;  s:= h/z;
                      f:= x*c + g*s;  g:= - x*s + g*c; h:= y*s;
                      y:= y*c;
                      FOR j:= 1 TO N DO
                      BEGIN x:= ab[j,i-1]; z:= ab[j,i];
                            ab[j,i-1]:= x*c + z*s;
                            ab[j,i]:= - x*s + z*c;
                      END;
                      IF Abs(f)<Abs(h)
                         THEN z:= Abs(h)*Sqrt(1+Sqr(f/h))
                         ELSE IF f<>0
                                 THEN z:= Abs(f)*Sqrt(1+Sqr(h/f))
                                 ELSE z:= 0;
                      q[i-1]:= z;
                      IF z=0 THEN
                      BEGIN z:= 1;
                            f:= 1;
                      END;
                      c:= f/z; s:= h/z;
                      f:= c*g + s*y; x:= -s*g + c*y;
                END; (* I *)
                e[l]:= 0; e[k]:= f; q[k]:= x;
                Goto TestFsplitting;
Convergence:    IF z<0 THEN
                BEGIN q[k]:= -z;
                      FOR j:= 1 TO N DO ab[j,k]:= - ab[j,k];
                END;
          END; (* K *)
END;      (* Minfit *)

FUNCTION Praxis(VAR X:Vector; N:INTEGER):REAL;
LABEL    l0, l1, l2;
VAR      i, j,
         k, k2,
         nl, nf, kl, kt: INTEGER;
         s, sl, dn, dmin,
         fx, f1,
         lds, ldt, sf, df,
         qf1, qd0, qd1,
         qa, qb, qc,
         m2, m4,
         small, vsmall,
         large, vlarge,
         ldfac, t2:      REAL;
         d, y, z,
         q0, q1:         Vector;
         v:              Matrix;

PROCEDURE SORT;                         (* sort d and v in descending order *)
VAR       k, i, j:       INTEGER;
          s:             REAL;
BEGIN     FOR i:= 1 TO N-1 DO
          BEGIN k:= i; s:= d[i];
                FOR j:= i+1 TO N DO
                    IF d[j] > s THEN
                    BEGIN k:= j; s:= d[j];
                    END;
                IF k>i THEN
                BEGIN d[k]:= d[i]; d[i]:= s;
                      FOR j:= 1 TO N DO
                      BEGIN s:= v[j,i];
                            v[j,i]:= v[j,k];
                            v[j,k]:= s;
                      END;
                END; (* IF *)
          END; (* FOR *)
END;      (* sort *)

PROCEDURE Print;                                   (* print a line of traces *)
VAR       i:             INTEGER;
BEGIN     WriteLn('------------------------------------------------------');
          WriteLn('Chi Square reduced to ... ', fx);
          WriteLn(' ... after ',nf,' function calls ...');
          WriteLn(' ... including ',nl,' linear searches.');
          WriteLn('Current Values of X ...');
          FOR i:= 1 TO N DO
              Write(X[i]);
          WriteLn;
END;      (* print *)

PROCEDURE MatPrint(s:PrString;v:Matrix;n,m:INTEGER);
VAR       k, i:          INTEGER;
BEGIN     WriteLn;
          WriteLn(s);
          FOR k:= 1 TO N DO
          BEGIN FOR i:= 1 TO m DO
                    Write(v[k,i]);
                WriteLn;
          END;
          WriteLn;
END;      (* MatPrint *)

PROCEDURE VecPrint(s:PrString;v:Vector;n:INTEGER);
VAR       i:             INTEGER;
BEGIN     WriteLn;
          WriteLn(s);
          FOR i:= 1 TO n DO
              Write(v[i]);
          WriteLn;
END;      (* VecPrint *)

PROCEDURE Min(j, Nits:INTEGER; VAR d2, x1:REAL; f1:REAL;fk:BOOLEAN);
LABEL     l0, l1;
VAR       k, i:          INTEGER;
          dz:            BOOLEAN;
          x2, xm, f0,
          f2, fm, d1, t2,
          s, sf1, sx1:   REAL;
 FUNCTION Flin(l:REAL):REAL;
 VAR      i:             INTEGER;
          t:             Vector;
 BEGIN    IF j>0 THEN                     (* linear search *)
             FOR i:= 1 TO N DO
                 t[i]:= x[i]+l*v[i,j]
             ELSE BEGIN              (* search along parabolic space curve *)
                  qa:= l*(l-qd1)/(qd0*(qd0+qd1));
                  qb:= (l+qd0)*(qd1-l)/(qd0*qd1);
                  qc:= l*(l+qd0)/(qd1*(qd0+qd1));
                  FOR i:= 1 TO N DO
                      t[i]:= qa*q0[i]+qb*x[i]+qc*q1[i];
             END; (* ELSE *)
          nf:= nf+1;
          Flin:= F(t, N);
 END;     (* Flin *)
BEGIN     (* Min *)
          sf1:= f1; sx1:= x1;
          k:= 0; xm:= 0; fm:= fx; f0:= fx; dz:= (d2<MachEps);
          (* Find step size *)
          s:= 0; FOR i:= 1 TO N DO s:= s + Sqr(X[i]);
          s:= Sqrt(s);
          IF dz
             THEN t2:= m4*Sqrt(Abs(fx)/Dmin + s*ldt) + m2*ldt
             ELSE t2:= m4*Sqrt(Abs(fx)/D2 + s*ldt) + m2*ldt;
          s:= m4*s + t;
          IF dz AND (t2>s) THEN t2:= s;
          IF (t2<small) THEN t2:= small;
          IF (t2>(0.01*h)) THEN t2:= 0.01*h;
          IF fk AND (f1<=fm) THEN BEGIN xm:= x1; fm:= f1; END;
          IF NOT fk OR (Abs(x1)<t2) THEN
          BEGIN IF x1>0 THEN x1:= t2 ELSE x1:= - t2;
                f1:= Flin(x1);
          END;
          IF f1<=fm THEN BEGIN xm:= x1; fm:= f1; END;
l0:       IF dz THEN
          BEGIN IF f0<f1 THEN x2:= - x1 ELSE x2:= 2*x1;
                f2:= Flin(x2);
                IF f2 <= fm THEN BEGIN xm:= x2; fm:= f2; END;
                d2:= (x2*(f1-f0) - x1*(f2-f0))/(x1*x2*(x1-x2));
          END;
          d1:= (f1-f0)/x1 - x1*d2; dz:= TRUE;
          IF d2 <= small
             THEN IF d1<0
                     THEN x2:= h
                     ELSE x2:= -h
             ELSE x2:= -0.5*d1/d2;
          IF Abs(x2)>h
             THEN IF x2>0
                     THEN x2:= h
                     ELSE x2:= -h;
l1:       f2:= Flin(x2);
          IF (k<Nits) AND (f2>f0) THEN
          BEGIN k:= k + 1;
                IF (f0<f1) AND ((x1*x2)>0) THEN GOTO l0;
                x2:= 0.5*x2; GOTO l1;
          END;
          nl:= nl + 1;
          IF f2>fm THEN x2:= xm ELSE fm:= f2;
          IF Abs(x2*(x2-x1))>small
             THEN d2:= (x2*(f1-f0) - x1*(fm-f0))/(x1*x2*(x1-x2))
             ELSE IF k>0
                     THEN d2:= 0;
          IF d2<=small THEN d2:= small;
          x1:= x2; fx:= fm;
          IF sf1<fx THEN BEGIN fx:= sf1; x1:= sx1; END;
          IF j>0
             THEN FOR i:= 1 TO N DO
                      x[i]:= x[i] + x1*v[i,j];
END;      (* Min *)

PROCEDURE Quad;        (* look for a minimum along the curve q0, q1, q2 *)
VAR       i:             INTEGER;
          l, s:          REAL;
BEGIN     s:= fx; fx:= qf1; qf1:= s; qd1:= 0;
          FOR i:= 1 TO N DO
          BEGIN s:= x[i];  l:= q1[i];  x[i]:= l;  q1[i]:= s;
                qd1:= qd1 + Sqr(s - l);
          END;
          s:= 0; qd1:= Sqrt(qd1); l:= qd1;
          IF (qd0>0) AND (qd1>0) AND (nl >= (3*N*N)) THEN
          BEGIN Min(0, 2, s, l, qf1, TRUE);
                qa:= l*(l-qd1)/(qd0*(qd0+qd1));
                qb:= (l+qd0)*(qd1-l)/(qd0*qd1);
                qc:= l*(l+qd0)/(qd1*(qd0+qd1));
          END ELSE
          BEGIN fx:= qf1; qa:= 0; qb:= 0; qc:= 1;
          END;
          qd0:= qd1;
          FOR i:= 1 TO N DO
          BEGIN s:= q0[i];  q0[i]:= x[i];
                x[i]:= qa*s + qb*x[i] + qc*q1[i];
          END;
END;      (* quad *)

BEGIN     (****  P R A X I S  ****)
          (* initialization *)
          small:= Sqr(MachEps); vsmall:= Sqr(small);
          large:= 1.0/small;    vlarge:= 1.0/vsmall;
          m2:= Sqrt(MachEps);   m4:= Sqrt(m2);
          IF IllC THEN ldfac:= 0.1 ELSE ldfac:= 0.01;
          nl:= 0; kt:= 0; nf:= 1; fx:= f(X, N); qf1:= fx;
          t2:= small + Abs(t); t:= t2; dmin:= small;
          IF h<(100*t) THEN h:= 100*t; ldt:= h;
          FOR i:= 1 TO N DO FOR j:= 1 TO N DO
              IF i=j THEN v[i,j]:= 1 ELSE v[i,j]:= 0;
          d[1]:= 0; qd0:= 0;
          FOR i:= 1 TO N DO q1[i]:= x[i];
          IF Prin > 1 THEN
          BEGIN WriteLn;WriteLn('---------- enter function praxis ------------');
                WriteLn('Current paramter settings are:');
                Writeln('... scaling ... ',scbd);
                WriteLn('... MachEps ... ',MachEps);
                WriteLn('...   Tol   ... ',t);
                WriteLn('... MaxStep ... ',h);
                WriteLn('...   IllC  ... ',IllC);
                WriteLn('...   Ktm   ... ',ktm);
          END;
          IF Prin > 0 THEN Print;

l0:       (* main loop *)
          sf:= d[1]; s:= 0; d[1]:= 0;
          (* minimize along first direction *)
          Min(1, 2, d[1], s, fx, FALSE);
          IF s<= 0
             THEN FOR i:= 1 TO N DO
                      v[i,1]:= - v[i,1];
          IF (sf<= (0.9*d[1])) OR ((0.9*sf)>=d[1])
             THEN FOR i:= 2 TO N DO
                      d[i]:= 0;
          FOR k:= 2 TO N DO
          BEGIN FOR i:= 1 TO N DO y[i]:= x[i];
                sf:= fx;
                IllC:= IllC OR (kt>0);
l1:             kl:= k; df:= 0;
                IF IllC THEN   (* random step to get off resolution valley *)
                BEGIN FOR i:= 1 TO N DO
                      BEGIN z[i]:= (0.1*ldt + t2*Power(10,kt))*(Random-0.5);
                            s:= z[i];
                            FOR j:= 1 TO N DO x[j]:= x[j]+s*v[j,i];
                      END; (* I *)
                      fx:= F(x, N);   nf:= nf + 1;
                END; (* IF *)
                FOR k2:= k TO N DO  (* minimize along non-conjugate directions *)
                BEGIN sl:= fx; s:= 0;
                      Min(k2, 2, d[k2], s, fx, FALSE);
                      IF IllC
                         THEN s:= d[k2]*Sqr(s+z[k2])
                         ELSE s:= sl - fx;
                      IF df<s THEN BEGIN df:= s;  kl:= k2; END;
                END; (* K2 *)
                IF NOT IllC AND (df < Abs(100*MachEps*fx)) THEN
                BEGIN IllC:= TRUE; GOTO l1;
                END;
                IF (k=2) AND (Prin>1) THEN VecPrint('New Direction ...', d, N);
                FOR k2:= 1 TO k-1 DO (* minimize along conjugate directions *)
                BEGIN s:= 0;
                      Min(k2, 2, d[k2], s, fx, FALSE);
                END; (* K2 *)
                f1:= fx; fx:= sf; lds:= 0;
                FOR i:= 1 TO N DO
                BEGIN sl:= x[i]; x[i]:= y[i]; y[i]:= sl - y[i]; sl:= y[i];
                      lds:= lds + Sqr(sl);
                END;
                lds:= Sqrt(lds);
                IF lds>small THEN
                BEGIN FOR i:= kl-1 DOWNTO k DO
                      BEGIN FOR j:= 1 TO N DO v[j,i+1]:= v[j,i];
                            d[i+1]:= d[i];
                      END;
                      d[k]:= 0;
                      FOR i:= 1 TO N DO v[i,k]:= y[i]/lds;
                      Min(k, 4, d[k], lds, f1, TRUE);
                      IF lds<=0 THEN
                      BEGIN lds:= -lds;
                            FOR i:= 1 TO N DO v[i,k]:= -v[i,k];
                      END;
                END; (* IF *)
                ldt:= ldfac*ldt; IF ldt<lds THEN ldt:= lds;
                IF Prin > 1 THEN Print;
                t2:= 0; FOR i:= 1 TO N DO t2:= t2 + Sqr(x[i]);
                t2:= m2*Sqrt(t2) + t;
                IF ldt > (0.5*t2) THEN kt:= 0 ELSE kt:= kt+1;
                IF kt>ktm THEN GOTO l2;
          END; (* K *)
          (*  try quadratic extrapolation in case    *)
          (*  we are stuck in a curved valley        *)
          Quad;
          dn:= 0;
          FOR i:= 1 TO N DO
          BEGIN d[i]:= 1.0/Sqrt(d[i]);
                IF dn<d[i] THEN dn:= d[i];
          END;
          IF Prin>2 THEN MatPrint('New Matrix of Directions ...', v, n, n);
          FOR j:= 1 TO N DO
          BEGIN s:= d[j]/dn;
                FOR i:= 1 TO N DO v[i,j]:= s*v[i,j];
          END;
          IF scbd > 1 THEN  (* scale axis to reduce condition number *)
          BEGIN s:= vlarge;
                FOR i:= 1 TO N DO
                BEGIN sl:= 0;
                      FOR j:= 1 TO N DO sl:= sl + Sqr(v[i,j]);
                      z[i]:= Sqrt(sl);
                      IF z[i]<m4 THEN z[i]:= m4;
                      IF s>z[i] THEN s:= z[i];
                END; (* I *)
                FOR i:= 1 TO N DO
                BEGIN sl:= s/z[i];  z[i]:= 1.0/sl;
                      IF z[i] > scbd THEN
                      BEGIN sl:= 1/scbd;
                            z[i]:= scbd;
                      END;
                END;
          END; (* IF *)
          FOR i:= 2 TO N DO  FOR j:= 1 TO i-1 DO
          BEGIN s:= v[i,j]; v[i,j]:= v[j,i]; v[j,i]:= s;
          END;
          MinFit(N, MachEps, Vsmall, v, d);
          IF scbd>1 THEN
          BEGIN FOR i:= 1 TO N DO
                BEGIN s:= z[i];
                      FOR j:= 1 TO n DO v[i,j]:= v[i,j]*s;
                END;
                FOR i:= 1 TO N DO
                BEGIN s:= 0;
                      FOR j:= 1 TO N DO s:= s + Sqr(v[j,i]);
                      s:= Sqrt(s);  d[i]:= s*d[i];  s:= 1.0/s;
                      FOR j:= 1 TO N DO v[j,i]:= s*v[j,i];
                END;
          END;
          FOR i:= 1 TO N DO
          BEGIN IF (dn*d[i])>Large
                   THEN d[i]:= vsmall
                   ELSE IF (dn*d[i])<small
                           THEN d[i]:= Vlarge
                           ELSE d[i]:= Power(dn*d[i],-2);
          END;
          Sort;   (* the new eigenvalues and eigenvectors *)
          Dmin:= d[n];
          IF Dmin < small THEN Dmin:= small;
          IllC:= (m2*d[1]) > Dmin;
          IF (Prin > 2) AND (scbd > 1)
             THEN VecPrint('Scale Factors ...', z, n);
          IF Prin > 2 THEN VecPrint('Eigenvalues of A ...', d, N);
          IF Prin > 2 THEN MatPrint('Eigenvectors of A ...', v, n, n);
          GOTO l0;  (* back to main loop *)
l2:       IF Prin > 0 THEN
          BEGIN VecPrint('Final solution is ...', x, n);
                WriteLn; WriteLn('ChiSq reduced to ',fx,' after ',nf, ' function calls.');
          END;
          Praxis:= fx;
END;      (****  Praxis  ****)
