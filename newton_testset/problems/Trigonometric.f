      subroutine Trigonometric(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Trigonometric function.
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Trigonometric.f
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
      integer j, k
      double precision h, sum, temp
c
      ifail = 0
      if (task .eq. 'XS') then
        if (n.eq.0) n = 10
		h = 1.0d0/dble(n)
		do j = 1, n
		   x0(j,1) = h
		enddo
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		sum = 0.0d0
		do j = 1, n
		   f(j) = dcos(x(j))
		   sum = sum + f(j)
		enddo
		do k = 1, n
		   f(k) = dble(n+k) - dsin(x(k)) - sum - dble(k)*f(k)
		enddo
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		do j = 1, n
		   temp = dsin(x(j))
		   do k = 1, n
			  dfdx(k,j) = temp
	       enddo
		   dfdx(j,j) = dble(j+1)*temp - dcos(x(j))
	    enddo
      endif
      return
      end
