(****************************************************************************)
(*              p r o c e d u r e      p r a x i s                          *)
(*                                                                          *)
(*   p r a x i s  is for the minimization of a function in several          *)
(*   variables. the algorithm used is a modification of a conjugate         *)
(*   gradient method developed by powell. changes are due to brent,         *)
(*   who gives an algol w program.                                          *)
(*   users who are interested in more of the details should read            *)
(*          - powell, m. j. d., 1962. an efficient method for finding       *)
(*            the minimum of a function of several variables without        *)
(*            calculating derivatives, computer journal 7, 155-162.         *)
(*          - brent, r. p., 1973. algorithms for minimization without       *)
(*            derivatives. prentice hall, englewood cliffs.                 *)
(*   if you have any further comments, questions, etc. please feel free     *)
(*   and contact me.                                                        *)
(*                                                                          *)
(*                           karl gegenfurtner     02/17/86                 *)
(*                                                 turbo - p version 1.0    *)
(****************************************************************************)
(*  the use of praxis is fairly simple. there are only two parameters:      *)
(*      - x is of type vector and contains on input an initial guess        *)
(*        of the solution. on output x holds the final solution of the      *)
(*        system of equations.                                              *)
(*      - n is of type integer and gives the number of unknown parameters   *)
(*  the result of praxis is the least calculated value of the function f    *)
(*  f is one of the global parameters used by praxis:                       *)
(*    - f(x,n) is declared forward and is the function to be minimized      *)
(*  all other globals are optional, i.e. you can change them or else        *)
(*  praxis will use some default values, which are adequate for most        *)
(*  problems under consideration.                                           *)
(*    - prin controls the printout from praxis                              *)
(*           0:  no printout at all                                         *)
(*           1:  only initial and final values                              *)
(*           2:  detailed map of the minimization process                   *)
(*           3:  also prints eigenvalues and eigenvectors of the            *)
(*               direction matices in use (for insiders only).              *)
(*    - t is the tolerance for the precision of the solution                *)
(*      praxis returns if the criterion                                     *)
(*             2 * ||x(k)-x(k-1)|| <= sqrt(macheps) * ||x(k)|| + t          *)
(*             is fulfilled more than ktm times, where                      *)
(*    - ktm is another global parameter. it's default value is 1 and        *)
(*      a value of 4 leads to a very(!) cautious stopping criterion         *)
(*    - macheps is the relative machine precision and is                    *)
(*      1.0 e-15 using an 8087 and turbo-87 and                             *)
(*      1.0 e-06 without 8087.                                              *)
(*    - h is a steplength parameter and shoul be set to the expected        *)
(*      distance to the solution. an exceptional large or small value       *)
(*      of h leads to slower 3 on the first few iterations.       *)
(*    - scbd is a scaling parameter and should be set to about 10.          *)
(*      the default is 1 and with that value no scaling is done at all      *)
(*      it's only necessary when the parameters are scaled very different   *)
(*    - illc is a boolean variable and should be set to true if the         *)
(*      problem is known to be illconditioned.                              *)
(****************************************************************************)


function power(a,b: real):real;
begin    power := exp(b*ln(a));
end;     (* power *)


procedure minfit(n:integer;eps,tol:real;var ab:matrix;var q:vector);
label     1,
          2,
          3,
          4;
var       l, kt,
          l2,
          i, j, k: integer;
          c, f, g,
          h, s, x,
          y, z:    real;
          e:       vector;
begin     (* householders reduction to bidiagonal form *)
          x:= 0; g:= 0;
          for i:= 1 to n do
          begin e[i]:= g; s:= 0; l:= i+1;
                for j:= i to n do
                    s:= s + sqr(ab[j,i]);
                if s<tol then g:= 0 else
                begin f:= ab[i,i];
                      if f<0
                         then g:= sqrt(s)
                         else g:= - sqrt(s);
                      h:= f*g-s;  ab[i,i]:= f - g;
                      for j:= l to n do
                      begin f:= 0;
                            for k:= i to n do
                                f:= f + ab[k,i]*ab[k,j];
                            f:= f/h;
                            for k:= i to n do
                                ab[k,j]:= ab[k,j] + f*ab[k,i];
                      end; (* j *)
                end; (* if *)
                q[i]:= g; s:= 0;
                if i<=n
                   then for j:= l to n do
                            s:= s + sqr(ab[i,j]);
                if s<tol then g:= 0 else
                begin f:= ab[i,i+1];
                      if f<0
                         then g:= sqrt(s)
                         else g:= - sqrt(s);
                      h:= f*g-s;  ab[i,i+1]:= f-g;
                      for j:= l to n do e[j]:= ab[i,j]/h;
                      for j:= l to n do
                      begin s:= 0;
                            for k:= l to n do s:= s + ab[j,k]*ab[i,k];
                            for k:= l to n do ab[j,k]:= ab[j,k] + s*e[k];
                      end; (* j *)
                end; (* if *)
                y:= abs(q[i])+abs(e[i]);
                if y > x then x:= y;
          end; (* i *)
          (* accumulation of right hand transformations *)
          for i:= n downto 1 do
          begin if g<>0.0 then
                begin h:= ab[i,i+1]*g;
                      for j:= l to n do ab[j,i]:= ab[i,j]/h;
                      for j:= l to n do
                      begin s:= 0;
                            for k:= l to n do s:= s + ab[i,k]*ab[k,j];
                            for k:= l to n do ab[k,j]:= ab[k,j] + s*ab[k,i];
                      end; (* j *)
                end; (* if *)
                for j:= l to n do
                begin ab[j,i]:= 0;
                      ab[i,j]:= 0;
                end;
                ab[i,i]:= 1; g:= e[i]; l:= i;
          end; (* i *)
          (* diagonalization to bidiagonal form *)
          eps:= eps*x;
          for k:= n downto 1 do
          begin kt:= 0;
1: kt:= kt + 1;
                if kt > 30 then
                begin e[k]:= 0;
                      writeln('+++ qr failed');
                end;
                for l2:= k downto 1 do
                begin l:= l2;
                      if abs(e[l])<=eps then goto 2;
                      if abs(q[l-1])<=eps then goto 4;
                end; (* l2 *)
4:   c:= 0; s:= 1;
                for i:= l to k do
                begin f:= s*e[i]; e[i]:= c*e[i];
                      if abs(f)<=eps then goto 2;
                      g:= q[i];
                      if abs(f) < abs(g)
                         then h:= abs(g)*sqrt(1+sqr(f/g))
                         else if f <> 0
                                 then h:= abs(f)*sqrt(1+sqr(g/f))
                                 else h:= 0;
                      q[i]:= h;
                      if h = 0 then
                      begin h:= 1;
                            g:= 1;
                      end;
                      c:= g/h; s:= -f/h;
                end; (* i *)
2:
                z:= q[k];
                if l=k then goto 3;
                (* shift from bottom 2*2 minor *)
                x:= q[l]; y:= q[k-1]; g:= e[k-1];  h:= e[k];
                f:= ((y-z)*(y+z) + (g-h)*(g+h))/(2*h*y);
                g:= sqrt(sqr(f)+1);
                if f<=0
                   then f:= ((x-z)*(x+z)+h*(y/(f-g)-h))/x
                   else f:= ((x-z)*(x+z)+h*(y/(f+g)-h))/x;
                (* next qr transformation *)
                s:= 1; c:= 1;
                for i:= l+1 to k do
                begin g:= e[i]; y:= q[i]; h:= s*g; g:= g*c;
                      if abs(f)<abs(h)
                         then z:= abs(h)*sqrt(1+sqr(f/h))
                         else if f<>0
                                 then z:= abs(f)*sqrt(1+sqr(h/f))
                                 else z:= 0;
                      e[i-1]:= z;
                      if z=0 then
                      begin f:= 1;
                            z:= 1;
                      end;
                      c:= f/z;  s:= h/z;
                      f:= x*c + g*s;  g:= - x*s + g*c; h:= y*s;
                      y:= y*c;
                      for j:= 1 to n do
                      begin x:= ab[j,i-1]; z:= ab[j,i];
                            ab[j,i-1]:= x*c + z*s;
                            ab[j,i]:= - x*s + z*c;
                      end;
                      if abs(f)<abs(h)
                         then z:= abs(h)*sqrt(1+sqr(f/h))
                         else if f<>0
                                 then z:= abs(f)*sqrt(1+sqr(h/f))
                                 else z:= 0;
                      q[i-1]:= z;
                      if z=0 then
                      begin z:= 1;
                            f:= 1;
                      end;
                      c:= f/z; s:= h/z;
                      f:= c*g + s*y; x:= -s*g + c*y;
                end; (* i *)
                e[l]:= 0; e[k]:= f; q[k]:= x;
                goto 1;
3:    if z<0 then
                begin q[k]:= -z;
                      for j:= 1 to n do ab[j,k]:= - ab[j,k];
                end;
          end; (* k *)
end;      (* minfit *)

function praxis(var x:vector; n:integer):real;
label    5, 6, 7;
var      i, j,
         k, k2,
         nl, nf, kl, kt: integer;
         s, sl, dn, dmin,
         fx, f1,
         lds, ldt, sf, df,
         qf1, qd0, qd1,
         qa, qb, qc,
         m2, m4,
         small, vsmall,
         large, vlarge,
         ldfac, t2:      real;
         d, y, z,
         q0, q1:         vector;
         v:              matrix;

procedure sort;                         (* sort d and v in descending order *)
var       k, i, j:       integer;
          s:             real;
begin     for i:= 1 to n-1 do
          begin k:= i; s:= d[i];
                for j:= i+1 to n do
                    if d[j] > s then
                    begin k:= j; s:= d[j];
                    end;
                if k>i then
                begin d[k]:= d[i]; d[i]:= s;
                      for j:= 1 to n do
                      begin s:= v[j,i];
                            v[j,i]:= v[j,k];
                            v[j,k]:= s;
                      end;
                end; (* if *)
          end; (* for *)
end;      (* sort *)

procedure print;                                   (* print a line of traces *)
var       i:             integer;
begin     writeln('------------------------------------------------------');
          writeln('chi square reduced to ... ', fx);
          writeln(' ... after ',nf,' function calls ...');
          writeln(' ... including ',nl,' linear searches.');
          writeln('current values of x ...');
          for i:= 1 to n do
              write(x[i]);
          writeln;
end;      (* print *)

procedure matprint(s:prstring;v:matrix;n,m:integer);
var       k, i:          integer;
begin     writeln;
          writeln(s);
          for k:= 1 to n do
          begin for i:= 1 to m do
                    write(v[k,i]);
                writeln;
          end;
          writeln;
end;      (* matprint *)

procedure vecprint(s:prstring;v:vector;n:integer);
var       i:             integer;
begin     writeln;
          writeln(s);
          for i:= 1 to n do
              write(v[i]);
          writeln;
end;      (* vecprint *)

procedure min(j, nits:integer; var d2, x1:real; f1:real;fk:boolean);
label     5, 6;
var       k, i:          integer;
          dz:            boolean;
          x2, xm, f0,
          f2, fm, d1, t2,
          s, sf1, sx1:   real;
 function flin(l:real):real;
 var      i:             integer;
          t:             vector;
 begin    if j>0 then                     (* linear search *)
             for i:= 1 to n do
                 t[i]:= x[i]+l*v[i,j]
             else begin              (* search along parabolic space curve *)
                  qa:= l*(l-qd1)/(qd0*(qd0+qd1));
                  qb:= (l+qd0)*(qd1-l)/(qd0*qd1);
                  qc:= l*(l+qd0)/(qd1*(qd0+qd1));
                  for i:= 1 to n do
                      t[i]:= qa*q0[i]+qb*x[i]+qc*q1[i];
             end; (* else *)
          nf:= nf+1;
          flin:= f(t, n);
 end;     (* flin *)
begin     (* min *)
          sf1:= f1; sx1:= x1;
          k:= 0; xm:= 0; fm:= fx; f0:= fx; dz:= (d2<macheps);
          (* find step size *)
          s:= 0; for i:= 1 to n do s:= s + sqr(x[i]);
          s:= sqrt(s);
          if dz
             then t2:= m4*sqrt(abs(fx)/dmin + s*ldt) + m2*ldt
             else t2:= m4*sqrt(abs(fx)/d2 + s*ldt) + m2*ldt;
          s:= m4*s + t;
          if dz and (t2>s) then t2:= s;
          if (t2<small) then t2:= small;
          if (t2>(0.01*h)) then t2:= 0.01*h;
          if fk and (f1<=fm) then begin xm:= x1; fm:= f1; end;
          if not fk or (abs(x1)<t2) then
          begin if x1>0 then x1:= t2 else x1:= - t2;
                f1:= flin(x1);
          end;
          if f1<=fm then begin xm:= x1; fm:= f1; end;
5:       if dz then
          begin if f0<f1 then x2:= - x1 else x2:= 2*x1;
                f2:= flin(x2);
                if f2 <= fm then begin xm:= x2; fm:= f2; end;
                d2:= (x2*(f1-f0) - x1*(f2-f0))/(x1*x2*(x1-x2));
          end;
          d1:= (f1-f0)/x1 - x1*d2; dz:= true;
          if d2 <= small
             then if d1<0
                     then x2:= h
                     else x2:= -h
             else x2:= -0.5*d1/d2;
          if abs(x2)>h
             then if x2>0
                     then x2:= h
                     else x2:= -h;
6:       f2:= flin(x2);
          if (k<nits) and (f2>f0) then
          begin k:= k + 1;
                if (f0<f1) and ((x1*x2)>0) then goto 5;
                x2:= 0.5*x2; goto 6;
          end;
          nl:= nl + 1;
          if f2>fm then x2:= xm else fm:= f2;
          if abs(x2*(x2-x1))>small
             then d2:= (x2*(f1-f0) - x1*(fm-f0))/(x1*x2*(x1-x2))
             else if k>0
                     then d2:= 0;
          if d2<=small then d2:= small;
          x1:= x2; fx:= fm;
          if sf1<fx then begin fx:= sf1; x1:= sx1; end;
          if j>0
             then for i:= 1 to n do
                      x[i]:= x[i] + x1*v[i,j];
end;      (* min *)

procedure quad;        (* look for a minimum along the curve q0, q1, q2 *)
var       i:             integer;
          l, s:          real;
begin     s:= fx; fx:= qf1; qf1:= s; qd1:= 0;
          for i:= 1 to n do
          begin s:= x[i];  l:= q1[i];  x[i]:= l;  q1[i]:= s;
                qd1:= qd1 + sqr(s - l);
          end;
          s:= 0; qd1:= sqrt(qd1); l:= qd1;
          if (qd0>0) and (qd1>0) and (nl >= (3*n*n)) then
          begin min(0, 2, s, l, qf1, true);
                qa:= l*(l-qd1)/(qd0*(qd0+qd1));
                qb:= (l+qd0)*(qd1-l)/(qd0*qd1);
                qc:= l*(l+qd0)/(qd1*(qd0+qd1));
          end else
          begin fx:= qf1; qa:= 0; qb:= 0; qc:= 1;
          end;
          qd0:= qd1;
          for i:= 1 to n do
          begin s:= q0[i];  q0[i]:= x[i];
                x[i]:= qa*s + qb*x[i] + qc*q1[i];
          end;
end;      (* quad *)

begin     (****  p r a x i s  ****)
          (* initialization *)
          small:= sqr(macheps); vsmall:= sqr(small);
          large:= 1.0/small;    vlarge:= 1.0/vsmall;
          m2:= sqrt(macheps);   m4:= sqrt(m2);
          if illc then ldfac:= 0.1 else ldfac:= 0.01;
          nl:= 0; kt:= 0; nf:= 1; fx:= f(x, n); qf1:= fx;
          t2:= small + abs(t); t:= t2; dmin:= small;
          if h<(100*t) then h:= 100*t; ldt:= h;
          for i:= 1 to n do for j:= 1 to n do
              if i=j then v[i,j]:= 1 else v[i,j]:= 0;
          d[1]:= 0; qd0:= 0;
          for i:= 1 to n do q1[i]:= x[i];
          if prin > 1 then
          begin writeln;writeln('---------- enter function praxis ------------');
                writeln('current paramter settings are:');
                writeln('... scaling ... ',scbd);
                writeln('... macheps ... ',macheps);
                writeln('...   tol   ... ',t);
                writeln('... maxstep ... ',h);
                writeln('...   illc  ... ',illc);
                writeln('...   ktm   ... ',ktm);
          end;
          if prin > 0 then print;

5:       (* main loop *)
          sf:= d[1]; s:= 0; d[1]:= 0;
          (* minimize along first direction *)
          min(1, 2, d[1], s, fx, false);
          if s<= 0
             then for i:= 1 to n do
                      v[i,1]:= - v[i,1];
          if (sf<= (0.9*d[1])) or ((0.9*sf)>=d[1])
             then for i:= 2 to n do
                      d[i]:= 0;
          for k:= 2 to n do
          begin for i:= 1 to n do y[i]:= x[i];
                sf:= fx;
                illc:= illc or (kt>0);
6:             kl:= k; df:= 0;
                if illc then   (* random step to get off resolution valley *)
                begin for i:= 1 to n do
                      begin z[i]:= (0.1*ldt + t2*power(10,kt))*(random(1)-0.5);
                            s:= z[i];
                            for j:= 1 to n do x[j]:= x[j]+s*v[j,i];
                      end; (* i *)
                      fx:= f(x, n);   nf:= nf + 1;
                end; (* if *)
                for k2:= k to n do  (* minimize along non-conjugate directions *)
                begin sl:= fx; s:= 0;
                      min(k2, 2, d[k2], s, fx, false);
                      if illc
                         then s:= d[k2]*sqr(s+z[k2])
                         else s:= sl - fx;
                      if df<s then begin df:= s;  kl:= k2; end;
                end; (* k2 *)
                if not illc and (df < abs(100*macheps*fx)) then
                begin illc:= true; goto 6;
                end;
                if (k=2) and (prin>1) then vecprint('new direction ...', d, n);
                for k2:= 1 to k-1 do (* minimize along conjugate directions *)
                begin s:= 0;
                      min(k2, 2, d[k2], s, fx, false);
                end; (* k2 *)
                f1:= fx; fx:= sf; lds:= 0;
                for i:= 1 to n do
                begin sl:= x[i]; x[i]:= y[i]; y[i]:= sl - y[i]; sl:= y[i];
                      lds:= lds + sqr(sl);
                end;
                lds:= sqrt(lds);
                if lds>small then
                begin for i:= kl-1 downto k do
                      begin for j:= 1 to n do v[j,i+1]:= v[j,i];
                            d[i+1]:= d[i];
                      end;
                      d[k]:= 0;
                      for i:= 1 to n do v[i,k]:= y[i]/lds;
                      min(k, 4, d[k], lds, f1, true);
                      if lds<=0 then
                      begin lds:= -lds;
                            for i:= 1 to n do v[i,k]:= -v[i,k];
                      end;
                end; (* if *)
                ldt:= ldfac*ldt; if ldt<lds then ldt:= lds;
                if prin > 1 then print;
                t2:= 0; for i:= 1 to n do t2:= t2 + sqr(x[i]);
                t2:= m2*sqrt(t2) + t;
                if ldt > (0.5*t2) then kt:= 0 else kt:= kt+1;
                if kt>ktm then goto 7;
          end; (* k *)
          (*  try quadratic extrapolation in case    *)
          (*  we are stuck in a curved valley        *)
          quad;
          dn:= 0;
          for i:= 1 to n do
          begin d[i]:= 1.0/sqrt(d[i]);
                if dn<d[i] then dn:= d[i];
          end;
          if prin>2 then matprint('new matrix of directions ...', v, n, n);
          for j:= 1 to n do
          begin s:= d[j]/dn;
                for i:= 1 to n do v[i,j]:= s*v[i,j];
          end;
          if scbd > 1 then  (* scale axis to reduce condition number *)
          begin s:= vlarge;
                for i:= 1 to n do
                begin sl:= 0;
                      for j:= 1 to n do sl:= sl + sqr(v[i,j]);
                      z[i]:= sqrt(sl);
                      if z[i]<m4 then z[i]:= m4;
                      if s>z[i] then s:= z[i];
                end; (* i *)
                for i:= 1 to n do
                begin sl:= s/z[i];  z[i]:= 1.0/sl;
                      if z[i] > scbd then
                      begin sl:= 1/scbd;
                            z[i]:= scbd;
                      end;
                end;
          end; (* if *)
          for i:= 2 to n do  for j:= 1 to i-1 do
          begin s:= v[i,j]; v[i,j]:= v[j,i]; v[j,i]:= s;
          end;
          minfit(n, macheps, vsmall, v, d);
          if scbd>1 then
          begin for i:= 1 to n do
                begin s:= z[i];
                      for j:= 1 to n do v[i,j]:= v[i,j]*s;
                end;
                for i:= 1 to n do
                begin s:= 0;
                      for j:= 1 to n do s:= s + sqr(v[j,i]);
                      s:= sqrt(s);  d[i]:= s*d[i];  s:= 1.0/s;
                      for j:= 1 to n do v[j,i]:= s*v[j,i];
                end;
          end;
          for i:= 1 to n do
          begin if (dn*d[i])>large
                   then d[i]:= vsmall
                   else if (dn*d[i])<small
                           then d[i]:= vlarge
                           else d[i]:= power(dn*d[i],-2);
          end;
          sort;   (* the new eigenvalues and eigenvectors *)
          dmin:= d[n];
          if dmin < small then dmin:= small;
          illc:= (m2*d[1]) > dmin;
          if (prin > 2) and (scbd > 1)
             then vecprint('scale factors ...', z, n);
          if prin > 2 then vecprint('eigenvalues of a ...', d, n);
          if prin > 2 then matprint('eigenvectors of a ...', v, n, n);
          goto 5;  (* back to main loop *)
7:        if prin > 0 then
          begin vecprint('final solution is ...', x, n);
                writeln; writeln('chisq reduced to ',fx,' after ',nf, ' function calls.');
          end;
          praxis:= fx;
end;      (****  praxis  ****)
