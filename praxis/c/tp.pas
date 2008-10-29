PROGRAM TestPraxis;

{$IPraxis.inc}

FUNCTION F;
VAR f1, f2, f3, f4: REAL;
BEGIN f1:= x[1]-1.0;
      f2:= 10*(x[2]-Sqr(x[1]));
      F:= (Sqr(f1)+Sqr(f2));
END;


VAR x: Vector;
    n: INTEGER;
    fin:REAL;

BEGIN
      X[1]:= -1.2;
      x[2]:= 1;
      N:= 2;
      Fin:= Praxis(X,N);
END.
