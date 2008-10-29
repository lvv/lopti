program testpraxis(input, output);


const maxpar = 20;



type  vector = array[1..maxpar] of real;
      matrix = array[1..maxpar, 1..maxpar] of real;
      prstring = array[1..80] of char;


var x: vector;
    n: integer;
    macheps:real;
    h:real;
    t:real;
    fin:real;
    prin:integer;
    illc:boolean;
    scbd:integer;
    ktm:integer;

function f(x:vector; n:integer):real;
var f1, f2: real;
begin f1:= x[1]-1.0;
      f2:= (10*(x[2]-sqr(x[1])));
      f:= sqr(f1)+sqr(f2);
end;

#include "praxis.i"


begin
      x[1]:= -1.2;
      x[2]:= 1;
      n:= 2;

      macheps  := 1.0e-8;
      h        := 1.0e+00;
      t        := 1.0e-08;
      prin     := 2;
      illc     := false;
      scbd     := 1;
      ktm      := 2;

      write('enter print paramter: '); readln(prin);
      fin:= praxis(x, n);
end.

