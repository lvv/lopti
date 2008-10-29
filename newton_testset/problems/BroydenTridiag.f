      subroutine BroydenTridiag(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Broyden tridiagonal function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             BroydenTridiag.f
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
      integer j, k
      double precision temp, temp1, temp2
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
		do k = 1, n
		   temp = (3.0d0 - 2.0d0*x(k))*x(k)
		   temp1 = 0.0d0
		   if (k .ne. 1) temp1 = x(k-1)
		   temp2 = 0.0d0
		   if (k .ne. n) temp2 = x(k+1)
		   f(k) = temp - temp1 - 2.0d0*temp2 + 1.0d0
		enddo
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		do k = 1, n
		   do j = 1, n
			  dfdx(k,j) = 0.0d0
		   enddo
		   dfdx(k,k) = 3.0d0 - 4.0d0*x(k)
		   if (k .ne. 1) dfdx(k,k-1) = -1.0d0
		   if (k .ne. n) dfdx(k,k+1) = -2.0d0
		enddo
      endif
      return
      end
