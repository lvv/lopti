      subroutine HelicalValley(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Helical valley function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             HelicalValley.f
c   Version          1.0
c   Latest Change    9th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     J.J. More, B.S. Garbow, K.E. Hillstrom:
c     Testing unconstrained optimization software.
c     ACM Trans. Math. Soft. 7, 17-41, 1981
c
c-----------------------------------------------------------------------
      double precision temp, temp1, temp2, tpi
      double precision c7, c8
      data c7, c8 /2.5d-1,5.0d-1/
c
      ifail = 0
      if (task .eq. 'XS') then
        n = 3
        x0(1,1) = -1.0d0
        x0(2,1) = 0.0d0
        x0(3,1) = 0.0d0
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		tpi = 8.0d0*datan(1.0d0)
		temp1 = dsign(c7,x(2))
		if (x(1) .gt. 0.0d0) temp1 = datan(x(2)/x(1))/tpi
		if (x(1) .lt. 0.0d0) temp1 = datan(x(2)/x(1))/tpi + c8
		temp2 = dsqrt(x(1)**2+x(2)**2)
		f(1) = 10.0d0*(x(3) - 10.0d0*temp1)
		f(2) = 10.0d0*(temp2 - 1.0d0)
		f(3) = x(3)
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		tpi = 8.0d0*datan(1.0d0)
		temp = x(1)**2 + x(2)**2
		temp1 = tpi*temp
		temp2 = dsqrt(temp)
		dfdx(1,1) = 100.0d0*x(2)/temp1
		dfdx(1,2) = -100.0d0*x(1)/temp1
		dfdx(1,3) = 10.0d0
		dfdx(2,1) = 10.0d0*x(1)/temp2
		dfdx(2,2) = 10.0d0*x(2)/temp2
		dfdx(2,3) = 0.0d0
		dfdx(3,1) = 0.0d0
		dfdx(3,2) = 0.0d0
		dfdx(3,3) = 1.0d0
      endif
      return
      end
