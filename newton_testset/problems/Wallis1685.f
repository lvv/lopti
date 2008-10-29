      subroutine Wallis1685(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Wallis function, used to illustate Newtons method in 1685
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Wallis1685.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     N. Kollerstrom:
c     Newton’s method of approximation, an enduring myth.
c     British Journal for History of Science,1992, 25, 347-354
c
c-----------------------------------------------------------------------
      if (task .eq. 'XS') then
        n = 1
		x0(1,1) = 0.82 
		x0(1,2) = -2.5 
		x0(1,3) = 0.0
		x0(1,4) = 5.0
		x0(1,5) = 1.0
		nx0 = 5
		xthrsh(1) = 1.0d0
		fthrsh(1) = 1.0d0
		rtol = 1.0d-12
c range of interest/solution
		xlow(1,1) =  -5.0
		xup (1,1) =   5.0
		nbound = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1) = x(1)**3 - 2.0d0*x(1) - 5.0d0
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        dfdx(1,1) = 3*x(1)**2 - 2.0d0
      endif
      ifail = 0
      return
      end
