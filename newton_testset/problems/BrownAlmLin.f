      subroutine BrownAlmLin(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Browns almost-linear function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             BrownAlmLin.f
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
      double precision sum, prod, temp
c
      ifail = 0
      if (task .eq. 'XS') then
        if (n.eq.0) n = 10
		do j = 1, n
		   x0(j,1) = 0.5d0
		enddo
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		sum = -dble(n+1)
		prod = 1.0d0
		do j = 1, n
		   sum = sum + x(j)
		   prod = x(j)*prod
		enddo
		do k = 1, n
		   f(k) = x(k) + sum
		enddo
		f(n) = prod - 1.0d0
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		prod = 1.0d0
		do j = 1, n
		   prod = x(j)*prod
		   do k = 1, n
			  dfdx(k,j) = 1.0d0
		   enddo
		   dfdx(j,j) = 2.0d0
		enddo
		do j = 1, n
		   temp = x(j)
		   if (temp .eq. 0.0d0) then
			 temp = 1.0d0
			 prod = 1.0d0
			 do k = 1, n
				if (k .ne. j) prod = x(k)*prod
			 enddo
		   endif
		   dfdx(n,j) = prod/temp
		enddo
      endif
      return
      end
