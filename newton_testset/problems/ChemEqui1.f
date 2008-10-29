      subroutine ChemEqui1(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Chemical equilibrium 1
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             ChemEqui1.f
c   Version          1.0
c   Latest Change    7th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     K.L. Hiebert:
c     An evaluation of mathematical software that
c     solves systems of nonlinear equations.
c     ACM Trans. Math. Soft. 8, 5-20, 1982
c
c-----------------------------------------------------------------------
      ifail = 0
      if (task .eq. 'XS') then
        n = 2
        x0(1,1)=1.d0
        x0(2,1)=1.d0
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1) = x(2) - 10.d0  
        f(2) = x(1)* x(2) - 5.0d+04
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		dfdx(1,1) = 0.d0  
		dfdx(2,1) = x(2)
		dfdx(1,2) = 1.d0  
		dfdx(2,2) = x(1)
      endif
      return
      end
