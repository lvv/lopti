      subroutine BroydenBanded(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Broyden banded function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             BroydenBanded.f
c   Version          1.0
c   Latest Change    7th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     J.J. More, B.S. Garbow, K.E. Hillstrom:
c     Testing unconstrained optimization software.
c     ACM Trans. Math. Soft. 7, 17-41, 1981
c
c-----------------------------------------------------------------------
c
      integer j, k, ml, mu, k1, k2
      double precision temp
c
      ifail = 0
      if (task .eq. 'XS') then
        if (n.eq.0) n = 10
		do j = 1, n
		   x0(j,1) = -1.0d0
		enddo
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		ml = 5
		mu = 1
		do k = 1, n
		   k1 = max0(1,k-ml)
		   k2 = min0(k+mu,n)
		   temp = 0.0d0
		   do j = k1, k2
			  if (j .ne. k) temp = temp + x(j)*(1.0d0 + x(j))
		   enddo
		   f(k) = x(k)*(2.0d0 + 5.0d0*x(k)**2) + 1.0d0 - temp
		enddo
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		ml = 5
		mu = 1
		do k = 1, n
		   do j = 1, n
			  dfdx(k,j) = 0.0d0
		   enddo
		   k1 = max0(1,k-ml)
		   k2 = min0(k+mu,n)
		   do j = k1, k2
			  if (j .ne. k) dfdx(k,j) = -(1.0d0 + 2.0d0*x(j))
		   enddo
		   dfdx(k,k) = 2.0d0 + 15.0d0*x(k)**2
		enddo
      endif
      return
      end
