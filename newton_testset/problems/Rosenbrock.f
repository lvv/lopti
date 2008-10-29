      subroutine Rosenbrock(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      implicit none
      character*(*) task
      integer n, ldfjac, np, nx0, ldx, nbound, ncon, nsol, ifail
      double precision x(*), f(*), dfdx(ldfjac,*), p(*), x0(ldx,*),
     $                 xlow(ldx,*), xup(ldx,*), conxl(ldx,*),
     $                 conxu(ldx,*), xsol(ldx,*), xthrsh(*), fthrsh(*),
     $                 rtol
c-----------------------------------------------------------------------
c
c   Rosenbrock function.
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Rosenbrock.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     J.J. More, B.S. Garbow, K.E. Hillstrom:
c     Testing unconstrained optimization software.
c     ACM Trans. Math. Soft. 7, 17-41, 1981
c
c-----------------------------------------------------------------------
      ifail = 0
      if (task .eq. 'XS') then
        n = 2
        nx0 = 1
        x0(1,1) = -1.2d0
        x0(2,1) = 1.0d0
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1) = 1.0d0 - x(1)
        f(2) = 10.0d0*(x(2) - x(1)**2)
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        dfdx(1,1) = -1.0d0
        dfdx(1,2) = 0.0d0
        dfdx(2,1) = -20.0d0*X(1)
        dfdx(2,2) = 10.0d0
      endif
      return
      end
