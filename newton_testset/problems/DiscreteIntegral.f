      subroutine DiscreteIntegral(n, x, f, dfdx, ldfjac, task, p, np,
     $           ldx, x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
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
c   Discrete integral equation function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             DiscreteIntegral.f
c   Version          1.0
c   Latest Change    8th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     J.J. More, B.S. Garbow, K.E. Hillstrom:
c     Testing unconstrained optimization software.
c     ACM Trans. Math. Soft. 7, 17-41, 1981
c
c-----------------------------------------------------------------------
      integer j, k, kp1
      double precision h, temp, tj, tk, sum1, sum2
c
      ifail = 0
      if (task .eq. 'XS') then
        if (n.eq.0) n = 10
		h = 1.0d0/dble(n+1)
		do j = 1, n
		   tj = dble(j)*h
		   x0(j,1) = tj*(tj - 1.0d0)
	    enddo
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		h = 1.0d0/dble(n+1)
		do k = 1, n
		   tk = dble(k)*h
		   sum1 = 0.0d0
		   do j = 1, k
			  tj = dble(j)*h
			  temp = (x(j) + tj + 1.0d0)**3
			  sum1 = sum1 + tj*temp
		   enddo
		   sum2 = 0.0d0
		   kp1 = k + 1
		   if (n .lt. kp1) go to 250
		   do j = kp1, n
			  tj = dble(j)*h
			  temp = (x(j) + tj + 1.0d0)**3
			  sum2 = sum2 + (1.0d0 - tj)*temp
		   enddo
  250      continue
		   f(k) = x(k) + h*((1.0d0 - tk)*sum1 + tk*sum2)/2.0d0
		enddo
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		h = 1.0d0/dble(n+1)
		do k = 1, n
		   tk = dble(k)*h
		   do j = 1, n
			  tj = dble(j)*h
			  temp = 3.0d0*(x(j) + tj + 1.0d0)**2
			  dfdx(k,j)=h*dmin1(tj*(1.0d0-tk),tk*(1.0d0-tj))*temp/2.0d0
		   enddo
		   dfdx(k,k) = dfdx(k,k) + 1.0d0
		enddo
      endif
      return
      end
