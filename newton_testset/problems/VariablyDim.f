      subroutine VariablyDim(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Variably dimensioned function.
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             VariablyDim.f
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
		   x0(j,1) = 1.0d0 - dble(j)*h
	    enddo
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		sum = 0.0d0
		do j = 1, n
		   sum = sum + dble(j)*(x(j) - 1.0d0)
		enddo
		temp = sum*(1.0d0 + 2.0d0*sum**2)
		do k = 1, n
		   f(k) = x(k) - 1.0d0 + dble(k)*temp
		enddo
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		sum = 0.0d0
		do j = 1, n
		   sum = sum + dble(j)*(x(j) - 1.0d0)
		enddo
		temp = 1.0d0 + 6.0d0*sum**2
		do k = 1, n
		   do j = k, n
			  dfdx(k,j) = dble(k*j)*temp
			  dfdx(j,k) = dfdx(k,j)
		   enddo
		   dfdx(k,k) = dfdx(k,k) + 1.0d0
		enddo
      endif
      return
      end
