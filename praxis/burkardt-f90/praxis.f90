function flin ( n, j, l, f, x, nf, v, q0, q1 )
!
!*******************************************************************************
!
!! FLIN is the function of one variable to be minimized by MINNY.
!
!
!  Discussion:
!
!    In fact, what is happening is that the scalar function F(X),
!    where X is an N dimensional vector, is being minimized along a 
!    fixed line.
!
!  Modified:
!
!    05 March 2002
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer N, the number of variables.
!
!    Input, integer J, indicates the kind of search.
!    If J is nonzero, then the search is linear in the direction of V(*,J).
!    If J is zero, then the search is parabolic, based on X, Q0 and Q1.
!
!    Input, real L, is the parameter determining the particular point
!    at which F is to be evaluated.  
!    For a linear search ( J is nonzero ), L is the size of the step
!    along the direction V(*,J).
!    For a quadratic search ( J is zero ), L is a parameter which specifies
!    a point in the plane of X, Q0 and Q1.
!
!    Input, external F, is the name of the function to be minimized.
!    The function should have the form 
!      function f(x,n)
!      integer n
!      real f
!      real x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    Input, real X(N), ?
!
!    Input/output, integer NF, the function evaluation counter.
!
!    Input, real V(N,N), a matrix whose columns constitute search directions.
!    If J is nonzero, then a linear search is being carried out in the
!    direction of V(*,J).
!
!    Input, real Q0(N), Q1(N), two auxiliary points used to determine
!    the plane when a quadratic search is performed.
!
!    Output, real FLIN, the value of the function at the given point.
!
  implicit none
!
  integer n
!
  real, external :: f
  real flin
  integer i
  integer j
  real l
  integer nf
  real q0(n)
  real q1(n)
  real qa
  real qb
  real qc
  real qd0
  real qd1
  real qf1
  real t(n)
  real v(n,n)
  real x(n)
!
  common /q/ qa,qb,qc,qd0,qd1,qf1
!
  if ( j /= 0 ) then
!
!  The search is linear.
!
    t(1:n) = x(1:n) + l * v(1:n,j)

  else
!
!  The search is along a parabolic space curve.
!
    qa = ( l * ( l - qd1 ) ) / ( qd0 * ( qd0 + qd1 ) )
    qb = ( ( l + qd0 ) * ( qd1 - l ) ) / ( qd0 * qd1 )
    qc = ( l * ( l + qd0 ) ) / ( qd1 * ( qd0 + qd1 ) )

    t(1:n) = qa * q0(1:n) + qb * x(1:n) + qc * q1(1:n)

  end if
!
!  The function evaluation counter NF is incremented.
!
  nf = nf + 1
!
!  Evaluate the function.
!
  flin = f(t,n)

  return
end
subroutine minfit ( m, n, tol, ab, q )
!
!*******************************************************************************
!
!! MINFIT computes the singular value decomposition of an N by N array.
!
!
!  Discussion:
!
!    This is an improved version of the EISPACK routine MINFIT
!    restricted to the case M = N and P = 0.
!
!    The singular values of the array AB are returned in Q.  AB is
!    overwritten with the orthogonal matrix V such that u.diag(q) = ab.v,
!    where U is another orthogonal matrix.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer M, the leading dimension of AB, which must be
!    at least N.
!
!    Input, integer N, the order of the matrix AB.
!
!    Input, real TOL, a tolerance which determines when a vector
!    (a column or part of a column of the matrix) may be considered
!    "essentially" equal to zero.
!
!    Input/output, real AB(M,N).  On input, an N by N array whose
!    singular value decomposition is desired.  On output, ?
!
!    Input, real Q(N), ?
!
  implicit none
!
  integer m
  integer n
!
  real ab(m,n)
  real c
  real e(n)
  real eps
  real f
  real g
  real h
  integer i
  integer ii
  integer j
  integer k
  integer kk
  integer kt
  integer, parameter :: kt_max = 30
  integer l
  integer l2
  integer ll2
  real q(n)
  real s
  real temp
  real tol
  real x
  real y
  real z
!
!  Householder's reduction to bidiagonal form.
!
  if ( n == 1 ) then
    q(1) = ab(1,1)
    ab(1,1) = 1.0E+00
    return
  end if

  eps = epsilon ( eps )
  g = 0.0E+00
  x = 0.0E+00

  do i = 1, n

    e(i) = g
    l = i + 1

    s = sum ( ab(1:n,i)**2 )

    g = 0.0E+00

    if ( s >= tol ) then

      f = ab(i,i)

      g = sqrt ( s )
      if ( f >= 0.0E+00 ) then
        g = -g
      end if

      h = f * g - s
      ab(i,i) = f - g

      do j = l, n

        f = dot_product ( ab(i:n,i), ab(i:n,j) ) / h

        ab(i:n,j) = ab(i:n,j) + f * ab(i:n,i)

      end do 

    end if

    q(i) = g

    s = sum ( ab(i,l:n)**2 )

    g = 0.0E+00

    if ( s >= tol ) then

      if ( i /= n ) then
        f = ab(i,i+1)
      end if

      g = sqrt ( s )
      if ( f >= 0.0E+00 ) then
        g = - g
      end if

      h = f * g - s

      if ( i /= n ) then

        ab(i,i+1) = f - g
        e(l:n) = ab(i,l:n) / h

        do j = l, n

          s = dot_product ( ab(j,l:n), ab(i,l:n) )

          ab(j,l:n) = ab(j,l:n) + s * e(l:n)

        end do

      end if

    end if

    y = abs ( q(i) ) + abs ( e(i) )

  end do

  x = max ( x, y )
!
!  Accumulation of right-hand transformations.
!
  ab(n,n) = 1.0E+00
  g = e(n)
  l = n

  do ii = 2, n

    i = n - ii + 1

    if ( g /= 0.0E+00 ) then

      h = ab(i,i+1) * g

      do j = l, n
        ab(j,i) = ab(i,j) / h
      end do

      do j = l, n

        s = dot_product ( ab(i,l:n), ab(l:n,j) )

        ab(l:n,j) = ab(l:n,j) + s * ab(l:n,i)

      end do

    end if

    ab(i,l:n) = 0.0E+00
    ab(l:n,i) = 0.0E+00
    ab(i,i) = 1.0E+00

    g = e(i)

  end do

  l = i
!
!  Diagonalization of the bidiagonal form.
!
  eps = eps * x

  do kk = 1, n

    k = n - kk + 1
    kt = 0

101 continue

    kt = kt + 1

    if ( kt > kt_max ) then
      e(k) = 0.0E+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MINFIT - Warning!'
      write ( *, '(a)' ) '  The QR algorithm failed to converge.'
    end if

    do ll2 = 1, k

      l2 = k - ll2 + 1
      l = l2

      if ( abs ( e(l) ) <= eps ) then
        go to 120
      end if

      if ( l /= 1 ) then
        if ( abs ( q(l-1) ) <= eps ) then
          exit
        end if
      end if

    end do
!
!  Cancellation of E(L) if L>1.
!
    c = 0.0E+00
    s = 1.0E+00

    do i = l, k

      f = s * e(i)
      e(i) = c * e(i)
      if ( abs ( f ) <= eps ) then
        go to 120
      end if
      g = q(i)
!
!  q(i) = h = sqrt(g*g + f*f).
!
      if ( abs ( f ) < abs ( g ) ) then
        h = abs ( g ) * sqrt ( 1.0E+00 + ( f / g )**2 )
      else if ( f == 0.0E+00 ) then
        h = 0.0E+00
      else
        h = abs ( f ) * sqrt ( 1.0E+00 + ( g / f )**2 )
      end if

      q(i) = h

      if ( h == 0.0E+00 ) then
        g = 1.0E+00
        h = 1.0E+00
      end if

      c = g / h
      s = - f / h

    end do
!
!  Test for convergence.
!
120 continue

    z = q(k)

    if ( l == k ) then
      if ( z < 0.0E+00 ) then
        q(k) = - z
        ab(1:n,k) = - ab(1:n,k)
      end if
      cycle
    end if
!
!  Shift from bottom 2*2 minor.
!
    x = q(l)
    y = q(k-1)
    g = e(k-1)
    h = e(k)
    f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0E+00 * h * y )

    g = sqrt ( f * f + 1.0E+00 )

    if ( f < 0.0E+00 ) then
      temp = f - g
    else
      temp = f + g
    end if

    f = ( ( x - z ) * ( x + z ) + h * ( y / temp - h ) ) / x
!
!  Next QR transformation.
!
    c = 1.0E+00
    s = 1.0E+00

    do i = l+1, k

      g = e(i)
      y = q(i)
      h = s * g
      g = g * c

      if ( abs ( f ) < abs ( h ) ) then
        z = abs ( h ) * sqrt ( 1.0E+00 + ( f / h )**2 )
      else if ( f == 0.0E+00 ) then
        z = 0.0E+00
      else
        z = abs ( f ) * sqrt ( 1.0E+00 + ( h / f )**2 )
      end if

      e(i-1) = z

      if ( z == 0.0E+00 ) then
        f = 1.0E+00
        z = 1.0E+00
      end if

      c = f / z
      s = h / z
      f =   x * c + g * s
      g = - x * s + g * c
      h = y * s
      y = y * c

      do j = 1, n
        x = ab(j,i-1)
        z = ab(j,i)
        ab(j,i-1) = x * c + z * s
        ab(j,i) = - x * s + z * c
      end do

      if ( abs ( f ) < abs ( h ) ) then
        z = abs ( h ) * sqrt ( 1.0E+00 + ( f / h ) **2 )
      else if ( f == 0.0E+00 ) then
        z = 0.0E+00
      else
        z = abs ( f ) * sqrt ( 1.0E+00 + ( h / f )**2 )
      end if

      q(i-1) = z

      if ( z == 0.0E+00 ) then
        f = 1.0E+00
        z = 1.0E+00
      end if

      c = f / z
      s = h / z
      f = c * g + s * y
      x = - s * g + c * y

    end do

    e(l) = 0.0E+00
    e(k) = f
    q(k) = x
    go to 101

  end do

  return
end
subroutine minny ( n, j, nits, d2, x1, f1, fk, f, x, t, h, v, q0, q1 )
!
!*******************************************************************************
!
!! MINNY minimizes a scalar function of N variables along a line.
!
!
!  Discussion:
!
!    MINNY minimizes F along the line from X in the direction V(*,j) unless
!    j is less than 1, when a quadratic search is made in the plane
!    defined by q0,q1,x.
!
!    If fk = .true., then f1 is flin(x1).  Otherwise x1 and f1 are ignored
!    on entry unless final fx is greater than f1.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer N, the number of variables.
!
!    Input, integer J, indicates the kind of search.
!    If J is nonzero, then the search is linear in the direction of V(*,J).
!    If J is zero, then the search is parabolic, based on X, Q0 and Q1.
!
!    Input, integer NITS, the maximum number of times the interval may be 
!    halved to retry the calculation.
!
!    Input, real D2, is either zero, or an approximation to the value
!    of (1/2) times the second derivative of F.
!
!    Input/output, real X1, on entry, an estimate of the distance from X 
!    to the minimum along V(*,J), or, if J = 0, a curve.  On output, 
!    the distance between X and the minimizer that was found.
!
!    Input/output, real F1, ?
!
!    Input, logical FK, ?
!
!    Input, external F, is the name of the function to be minimized.
!    The function should have the form 
!      function f(x,n)
!      integer n
!      real f
!      real x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    ?, real X(N), ?
!
!    ?, real T, ?
!
!    ?, real H, ?
!
!    Input, real V(N,N), a matrix whose columns are direction vectors
!    along which the function may be minimized.
!
!    ?, real Q0(N), ?
!
!    ?, real Q1(N), ?
!
  implicit none
!
  real d1
  real d2
  real dmin
  logical dz
  real, external :: f
  real f0
  real f1
  real f2
  logical fk
  real flin
  real fm
  real fx
  real h
  integer i
  integer j
  integer k
  real ldt
  real m2
  real m4
  real machep
  integer n
  integer nf
  integer nits
  integer nl
  real q0(n)
  real q1(n)
  real qa
  real qb
  real qc
  real qd0
  real qd1
  real qf1
  real s
  real sf1
  real small
  real sx1
  real t
  real t2
  real temp
  real v(n,n)
  real x(n)
  real x1
  real x2
  real xm
!
  common /global/ fx,ldt,dmin,nf,nl
  common /q/ qa,qb,qc,qd0,qd1,qf1
!
  machep = epsilon ( machep )
  small = machep**2
  m2 = sqrt ( machep )
  m4 = sqrt ( m2 )
  sf1 = f1
  sx1 = x1
  k = 0
  xm = 0.0E+00
  fm = fx
  f0 = fx
  dz = ( d2 < machep )
!
!  Find the step size.
! 
  s = sqrt ( sum ( x(1:n)**2 ) )

  if ( dz ) then
    temp = dmin
  else
    temp = d2
  end if

  t2 = m4 * sqrt ( abs ( fx ) / temp + s * ldt ) + m2 * ldt
  s = m4 * s + t
  if ( dz .and. t2 > s ) then
    t2 = s
  end if

  t2 = max ( t2, small )
  t2 = min ( t2, 0.01E+00 * h )

  if ( fk .and. f1 <= fm ) then
    xm = x1
    fm = f1
  end if

  if ( .not. fk .or. abs ( x1 ) < t2 ) then

    if ( x1 >= 0.0E+00 ) then
      temp = 1.0E+00
    else
      temp = - 1.0E+00
    end if

    x1 = temp * t2
    f1 = flin ( n, j, x1, f, x, nf, v, q0, q1 )

  end if

  if ( f1 <= fm ) then
    xm = x1
    fm = f1
  end if
!
!  Evaluate FLIN at another point and estimate the second derivative.
!
4 continue

  if ( dz ) then

    if ( f0 >= f1 ) then
      x2 = 2.0E+00 * x1
    else
      x2 = - x1
    end if

    f2 = flin ( n, j, x2, f, x, nf, v, q0, q1 )

    if ( f2 <= fm ) then
      xm = x2
      fm = f2
    end if

    d2 = ( x2 * ( f1 - f0 ) - x1 * ( f2 - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )

  end if
!
!  Estimate the first derivative at 0.
!
  d1 = ( f1 - f0 ) / x1 - x1 * d2
  dz = .true.
!
!  Predict the minimum.
!
  if ( d2 <= small ) then

    if ( d1 >= 0.0E+00 ) then
      x2 = - h
    else
      x2 = h
    end if

  else

    x2 = ( - 0.5E+00 * d1 ) / d2

  end if

  if ( abs ( x2 ) > h ) then

    if ( x2 <= 0.0E+00 ) then
      x2 = - h
    else
      x2 = h
    end if

  end if
!
!  Evaluate F at the predicted minimum.
!
  do

    f2 = flin ( n, j, x2, f, x, nf, v, q0, q1 )

    if ( k >= nits .or. f2 <= f0 ) then
      exit
    end if

    k = k + 1

    if ( f0 < f1 .and. x1 * x2 > 0.0E+00 ) then
      go to 4
    end if

    x2 = 0.5E+00 * x2

  end do
!
!  Increment the one-dimensional search counter.
!
  nl = nl + 1

  if ( f2 > fm ) then
    x2 = xm
  else
    fm = f2
  end if
!
!  Get a new estimate of the second derivative.
!
  if ( abs ( x2 * ( x2 - x1 ) ) > small ) then
    d2 = ( x2 * ( f1 - f0 ) - x1 * ( fm - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )
  else
    if ( k > 0 ) then
      d2 = 0.0E+00
    end if
  end if

  d2 = max ( d2, small )

  x1 = x2
  fx = fm

  if ( sf1 < fx ) then
    fx = sf1
    x1 = sx1
  end if
!
!  Update X for linear but not parabolic search.
!
  if ( j /= 0 ) then

    x(1:n) = x(1:n) + x1 * v(1:n,j)

  end if

  return
end
function praxis ( t0, h0, n, prin, x, f )
!
!*******************************************************************************
!
!! PRAXIS seeks an N-dimensional minimizer X of a scalar function F(X).
!
!
!  Discussion:
!
!    PRAXIS returns the minimum of the function F(X,N) of N variables
!    using the principal axis method.  The gradient of the function is
!    not required.
!
!    The approximating quadratic form is
!
!      Q(x') = F(x,n) + (1/2) * (x'-x)-transpose * A * (x'-x)
!
!    where X is the best estimate of the minimum and 
!
!      A = inverse(V-transpose) * D * inverse(V)
!
!   V(*,*) is the matrix of search directions; d(*) is the array
!   of second differences.  
!
!   If F(X) has continuous second derivatives near X0, then A will tend 
!   to the hessian of F at X0 as X approaches X0.
!
!  Modified:
!
!    25 February 2002
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Finding Zeros and Extrema of Functions Without
!      Calculating Derivatives,
!    Stanford University Technical Report STAN-CS-71-198.
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, real T0, is a tolerance.  PRAXIS attempts to return 
!    praxis = f(x) such that if X0 is the true local minimum near X, then
!    norm ( x - x0 ) < T0 + sqrt ( EPSILON ( X ) ) * norm ( X ),
!    where EPSILON ( X ) is the machine precision for X.
!
!    Input, real H0, is the maximum step size.  H0 should be set to about the
!    maximum distance from the initial guess to the minimum.
!    If H0 is set too large or too small, the initial rate of
!    convergence may be slow.
!
!    Input, integer N, the number of variables.
!
!    Input, integer PRIN, controls the printing of intermediate results.
!    0, nothing is printed.
!    1, F is printed after every n+1 or n+2 linear minimizations.  
!       final X is printed, but intermediate X is printed only 
!       if N is at most 4.
!    2, the scale factors and the principal values of the approximating 
!       quadratic form are also printed.
!    3, X is also printed after every few linear minimizations.
!    4, the principal vectors of the approximating quadratic form are 
!       also printed.
!
!    Input/output, real X(N), is an array containing on entry a guess 
!    of the point of minimum, on return the estimated point of minimum.
!
!    Input, external F, is the name of the function to be minimized.
!    The function should have the form 
!      function f(x,n)
!      integer n
!      real f
!      real x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    Output, real PRAXIS, the function value at the minimizer.
!
  implicit none
!
  integer n
!
  real d(n)
  real df
  real dmin
  real dn
  real dni
  real, external :: f
  real f1
  real fx
  real h
  real h0
  integer i
  integer ii
  logical illc
  integer, save :: iseed = 1234567
  integer j
  integer k
  integer k2
  integer kl
  integer kt
  integer ktm
  real large
  real ldfac
  real lds
  real ldt
  real m2
  real m4
  real machep
  integer nits
  integer nl
  integer nf
  real praxis
  integer prin
  real q0(n)
  real q1(n)
  real qa
  real qb
  real qc
  real qd0
  real qd1
  real qf1
  real r
  real s
  real scbd
  real sf
  real sl
  real small
  real t
  real t0
  real t2
  real v(n,n)
  real value
  real vlarge
  real vsmall
  real x(n)
  real y(n)
  real z(n)
!
  common /global/ fx,ldt,dmin,nf,nl
  common /q/ qa,qb,qc,qd0,qd1,qf1
!
!  Initialization.
!
  machep = epsilon ( machep )
  small = machep * machep
  vsmall = small * small
  large = 1.0E+00 / small
  vlarge = 1.0E+00 / vsmall
  m2 = sqrt ( machep )
  m4 = sqrt ( m2 )
!
!  heuristic numbers:
!
!  If the axes may be badly scaled (which is to be avoided if
!  possible), then set scbd = 10.  otherwise set scbd=1.
!
!  If the problem is known to be ill-conditioned, set ILLC = true.
!
!  KTM is the number of iterations without improvement before the
!  algorithm terminates.  KTM = 4 is very cautious; usually KTM = 1
!  is satisfactory.
!
  scbd = 1.0E+00
  illc = .false.
  ktm = 1

  if ( illc ) then
    ldfac = 0.1E+00
  else
    ldfac = 0.01E+00
  end if

  kt = 0
  nl = 0
  nf = 1
  fx = f(x,n)
  qf1 = fx
  t = small + abs ( t0 )
  t2 = t
  dmin = small
  h = h0
  h = max ( h, 100.0E+00 * t )
  ldt = h
!
!  The initial set of search directions V is the identity matrix.
!
  v(1:n,1:n) = 0.0E+00
  do i = 1, n
    v(i,i) = 1.0E+00
  end do

  d(1) = 0.0E+00
  qd0 = 0.0E+00
  q0(1:n) = x(1:n)
  q1(1:n) = x(1:n)

  if ( prin > 0 ) then
    call print2 ( n, x, prin, fx, nf, nl )
  end if
!
!  The main loop starts here.
!
  do

    sf = d(1)
    d(1) = 0.0E+00
    s = 0.0E+00
!
!  Minimize along the first direction v(*,1).
!
    nits = 2
    value = fx

    call minny ( n, 1, nits, d(1), s, value, .false., f, x, t, h, &
      v, q0, q1 )

    if ( s <= 0.0E+00 ) then
      v(1:n,1) = - v(1:n,1)
    end if

    if ( sf <= 0.9E+00 * d(1) .or. 0.9E+00 * sf >= d(1) ) then
      d(2:n) = 0.0E+00
    end if
!
!  The inner loop starts here.
!
    do k = 2, n

      y(1:n) = x(1:n)

      sf = fx

      if ( kt > 0 ) then
        illc = .true.
      end if

80    continue

      kl = k
      df = 0.0E+00
!
!  A random step follows (to avoid resolution valleys).
!  PRAXIS assumes that the random number generator returns a random 
!  number uniformly distributed in (0,1).
!
      if ( illc ) then

        do i = 1, n
          call r_random ( 0.0E+00, 1.0E+00, r )
          s = ( 0.1E+00 * ldt + t2 * 10.0E+00**kt ) * ( r - 0.5E+00 )
          z(i) = s
          x(1:n) = x(1:n) + s * v(1:n,i)
        end do

        fx = f(x,n)
        nf = nf + 1

      end if
!
!  Minimize along the "non-conjugate" directions V(*,K),...,V(*,N).
!
      do k2 = k, n

        sl = fx
        s = 0.0E+00
        nits = 2
        value = fx

        call minny ( n, k2, nits, d(k2), s, value, .false., f, x, t, &
          h, v, q0, q1 )

        if ( illc ) then
          s = d(k2) * ( ( s + z(k2) )**2 )
        else
          s = sl - fx
        end if

        if ( df <= s ) then
          df = s
          kl = k2
        end if

      end do
!
!  If there was not much improvement on the first try, set
!  ILLC = true and start the inner loop again.
!
      if ( .not. illc ) then
        if ( df < abs ( 100.0E+00 * machep * fx ) ) then
          illc = .true.
          go to 80
        end if
      end if

      if ( k == 2 .and. prin > 1 ) then
        call rvec_print ( n, d, '  The second difference array' )
      end if
!
!  Minimize along the "conjugate" directions V(*,1),...,V(*,K-1).
!
      do k2 = 1, k-1

        s = 0.0E+00
        nits = 2
        value = fx

        call minny ( n, k2, nits, d(k2), s, value, .false., f, x, t, &
          h, v, q0, q1 )

      end do

      f1 = fx
      fx = sf
      lds = 0

      do i = 1, n
        sl = x(i)
        x(i) = y(i)
        sl = sl - y(i)
        y(i) = sl
        lds = lds + sl**2
      end do

      lds = sqrt ( lds )
!
!  Discard direction V(*,kl).
!
!  If no random step was taken, V(*,KL) is the "non-conjugate"
!  direction along which the greatest improvement was made.
!
      if ( lds > small ) then

        do ii = 1, kl-k
          i = kl - ii
          do j = 1, n
            v(j,i+1) = v(j,i)
          end do
          d(i+1) = d(i)
        end do

        d(k) = 0
        v(1:n,k) = y(1:n) / lds
!
!  Minimize along the new "conjugate" direction V(*,k), which is
!  the normalized vector:  (new x) - (old x).
!
        nits = 4
        value = f1

        call minny ( n, k, nits, d(k), lds, value, .true., f, x, t, &
          h, v, q0, q1 )

        if ( lds <= 0.0E+00 ) then
          lds = - lds
          v(1:n,k) = - v(1:n,k)
        end if

      end if

      ldt = ldfac * ldt
      ldt = max ( ldt, lds )

      if ( prin > 0 ) then
        call print2 ( n, x, prin, fx, nf, nl )
      end if

      t2 = m2 * sqrt ( sum ( x(1:n)**2 ) ) + t
!
!  See whether the length of the step taken since starting the
!  inner loop exceeds half the tolerance.
!
      if ( ldt > 0.5E+00 * t2 ) then
        kt = -1
      end if

      kt = kt + 1

      if ( kt > ktm ) then
        go to 400
      end if

    end do
!
!  The inner loop ends here.
!
!  Try quadratic extrapolation in case we are in a curved valley.
!
171 continue

    call quad ( n, f, x, t, h, v, q0, q1 )

    d(1:n) = 1.0E+00 / sqrt ( d(1:n) )

    dn = maxval ( d(1:n) )

    if ( prin > 3 ) then
      call rmat_print ( v, n, n, n, '  The new direction vectors' )
    end if

    do j = 1, n
      v(1:n,j) = ( d(j) / dn ) * v(1:n,j)
    end do
!
!  Scale the axes to try to reduce the condition number.
!
    if ( scbd > 1.0E+00 ) then

      do i = 1, n
        z(i) = sqrt ( sum ( v(i,1:n)**2 ) )
        z(i) = max ( z(i), m4 )
      end do

      s = minval ( z(1:n) )

      do i = 1, n

        sl = s / z(i)
        z(i) = 1.0E+00 / sl

        if ( z(i) > scbd ) then
          sl = 1.0E+00 / scbd
          z(i) = scbd
        end if

        v(i,1:n) = sl * v(i,1:n)

      end do

    end if
!
!  Calculate a new set of orthogonal directions before repeating
!  the main loop.
!
!  Transpose V for MINFIT:
!
    v(1:n,1:n) = transpose ( v(1:n,1:n) )
!
!  Call MINFIT to find the singular value decomposition of V.
!
!  This gives the principal values and principal directions of the
!  approximating quadratic form without squaring the condition number.
!
    call minfit ( n, n, vsmall, v, d )
!
!  Unscale the axes.
!
    if ( scbd > 1.0E+00 ) then

      do i = 1, n
        v(i,1:n) = z(i) * v(i,1:n)
      end do

      do i = 1, n

        s = sqrt ( sum ( v(1:n,i)**2 ) )

        d(i) = s * d(i)
        v(1:n,i) = v(1:n,i) / s

      end do

    end if

    do i = 1, n

      dni = dn * d(i)

      if ( dni > large ) then
        d(i) = vsmall
      else if ( dni < small ) then
        d(i) = vlarge
      else
        d(i) = 1.0E+00 / dni**2
      end if

    end do
!
!  Sort the eigenvalues and eigenvectors.
!
    call sort ( n, n, d, v )
!
!  Determine the smallest eigenvalue.
!
    dmin = max ( d(n), small )
!
!  The ratio of the smallest to largest eigenvalue determines whether
!  the system is ill conditioned.
!
    if ( m2 * d(1) > dmin ) then
      illc = .true.
    else
      illc = .false.
    end if

    if ( prin > 1 ) then

      if ( scbd > 1.0E+00 ) then
        call rvec_print ( n, z, '  The scale factors' )
      end if 

      call rvec_print ( n, d, '  Principal values of the quadratic form' )

    end if

    if ( prin > 3 ) then
      call rmat_print ( v, n, n, n, '  The principal axes:' )
    end if
!
!  The main loop ends here.
!
  end do

400   continue

  if ( prin > 0 ) then
    call rvec_print ( n, x, '  X:' )
  end if

  praxis = fx

  return
end
subroutine print2 ( n, x, prin, fx, nf, nl )
!
!*******************************************************************************
!
!! PRINT2 prints certain data about the progress of the iteration.
!
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer N, the number of variables.
!
!    Input, real X(N), the current estimate of the minimizer.
!
!    Input, integer PRIN, ?
!
!    Input, real FX, the smallest value of F(X) found so far.
!
!    Input, integer NF, the number of function evaluations.
!
!    Input, integer NL, the number of linear searches.
!
  implicit none
!
  integer n
!
  real fx
  integer nf
  integer nl
  integer prin
  real x(n)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Linear searches      ', nl
  write ( *, '(a,i6)' ) '  Function evaluations ', nf 

  if ( n <= 4 .or. prin > 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'X:'
    write ( *, '(5g14.6)' ) x(1:n)
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) 'FX: ', fx

  return
end
subroutine quad ( n, f, x, t, h, v, q0, q1 )
!
!*******************************************************************************
!
!! QUAD seeks to minimize the scalar function F along a particular curve.
!
!
!  Discussion:
!
!    The minimizer to be sought is required to lie on a curve defined
!    by Q0, Q1 and X.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer N, the number of variables.
!
!    Input, external F, is the name of the function to be minimized.
!    The function should have the form 
!      function f(x,n)
!      integer n
!      real f
!      real x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    Input, real X(N), ?
!
!    Input, real T, ?
!
!    Input, real H, ?
!
!    Input, real V(N,N), ?
!
!    Input, real Q0(N), Q1(N), ?
!
  implicit none
!
  integer n
!
  real dmin
  real, external :: f
  real fx
  real h
  integer i
  real l
  real ldt
  real machep
  integer nf
  integer nits
  integer nl
  real q0(n)
  real q1(n)
  real qa
  real qb
  real qc
  real qd0
  real qd1
  real qf1
  real s
  real t
  real v(n,n)
  real value
  real x(n)
!
  common /global/ fx,ldt,dmin,nf,nl
  common /q/ qa,qb,qc,qd0,qd1,qf1
!
  machep = epsilon ( machep )

  call r_swap ( fx, qf1 )

  call rvec_swap ( n, x, q1 )

  qd1 = sqrt ( sum ( ( x(1:n) - q1(1:n) )**2 ) )

  l = qd1
  s = 0.0E+00

  if ( qd0 <= 0.0E+00 .or. qd1 <= 0.0E+00 .or. nl < 3 * n**2 ) then

    fx = qf1
    qa = 0.0E+00
    qb = 0.0E+00
    qc = 1.0E+00

  else

    nits = 2
    value = qf1

    call minny ( n, 0, nits, s, l, value, .true., f, x, t, &
      h, v, q0, q1 )

    qa = ( l * ( l - qd1 ) ) / ( qd0 * ( qd0 + qd1 ) )
    qb = ( ( l + qd0 ) * ( qd1 - l ) ) / ( qd0 * qd1 )
    qc = ( l * ( l + qd0 ) ) / ( qd1 * ( qd0 + qd1 ) )

  end if

  qd0 = qd1

  do i = 1, n
    s = q0(i)
    q0(i) = x(i)
    x(i) = ( qa * s + qb * x(i) ) + qc * q1(i)
  end do

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Output, real R, the randomly chosen value.
!
  implicit none
!
  real r
  real rhi
  real rlo
  logical, save :: seed = .false.
  real t
!
  if ( .not. seed ) then
    call random_seed ( )
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine rmat_print ( v, lda, m, n, label )
!
!*******************************************************************************
!
!! RMAT_PRINT prints out a matrix.
!
!
!  Discussion:
!
!    The matrix is printed out five columns at a time.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real V(LDA,N), an M by N matrix.
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer M, N, the number of rows and columns in the matrix.
!
!    Input, character ( len = * ) LABEL, a label for the matrix.
!
  implicit none
!
  integer lda
  integer n
!
  integer i
  integer j
  integer jhi
  integer jlo
  character ( len = * ) label
  integer m
  real v(lda,n)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) label

  do jlo = 1, n, 5

    jhi = min ( jlo + 4, n )

    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(5g14.6)' ) v(i,jlo:jhi)
    end do

  end do

  return
end
subroutine rvec_print ( n, a, title )
!
!*******************************************************************************
!
!! RVEC_PRINT prints a real vector, with an optional title.
!
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) title
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine rvec_swap ( n, x, y )
!
!*******************************************************************************
!
!! RVEC_SWAP switches two real vectors.
!
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  integer n
!
  real x(n)
  real y(n)
  real z(n)
!
  z(1:n) = x(1:n)
  x(1:n) = y(1:n)
  y(1:n) = z(1:n)

  return
end
subroutine sort ( m, n, d, v ) 
!
!*******************************************************************************
!
!! SORT sorts a vector D and adjusts the corresponding columns of a matrix V.
!
!
!  Discussion:
!
!    A simple bubble sort is used on D.
!
!    In our application, D contains eigenvalues, and the columns of V are
!    the corresponding eigenvectors.
!
!  Modified:
!
!    25 February 2002
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer M, the row dimension of V, which must be at least N.
!
!    Input, integer N, the length of D, and the order of V.
!
!    Input/output, real D(N), the vector to be sorted.  On output, the
!    entries of D are in descending order.
!
!    Input/output, real V(M,N), an N by N array to be adjusted as D is
!    sorted.  In particular, if the value that was in D(I) on input is
!    moved to D(J) on output, then the input column V(*,I) is moved to
!    the output column V(*,J).
!
  implicit none
!
  integer m
  integer n
!
  real d(n)
  integer i
  integer k(1)
  real v(m,n)
!
  do i = 1, n-1
!
!  Find K, the index of the largest entry in D(I:N).
!  MAXLOC apparently requires its output to be an array.
!
    k = maxloc ( d(i:n) )
!
!  If K > I, swap D(K) and D(I), and columns K and I of V.
!
    if ( k(1) > i ) then

      call r_swap ( d(i), d(k(1)) )

      call rvec_swap ( n, v(1:n,i), v(1:n,k(1)) )

    end if

  end do

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
