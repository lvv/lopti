      subroutine Chebyquad(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Chebyquad function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Chebyquad.f
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
      integer i, j, k, iev
      double precision h, temp, temp1, temp2, temp3, temp4, ti, tk
c
      ifail = 0
      if (task .eq. 'XS') then
        if (n.eq.0) n = 9
        h = 1.0d0/dble(n+1)
        do j = 1, n
           x0(j,1) = dble(j)*h
        enddo
        nx0 = 1
        do j = 1, n
          xlow(j,1) = 0.0d0
          xup (j,1) = 1.0d0
        enddo
        nbound = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		do k = 1, n
		   f(k) = 0.0d0
		enddo
		do j = 1, n
		   temp1 = 1.0d0
		   temp2 = 2.0d0*x(j) - 1.0d0
		   temp = 2.0d0*temp2
		   do i = 1, n
			  f(i) = f(i) + temp2
			  ti = temp*temp2 - temp1
			  temp1 = temp2
			  temp2 = ti
		   enddo
		enddo
		tk = 1.0d0/dble(n)
		iev = -1
		do k = 1, n
		   f(k) = tk*f(k)
		   if (iev .gt. 0) f(k) = f(k) + 1.0d0/(dble(k)**2 - 1.0d0)
		   iev = -iev
		enddo
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		tk = 1.0d0/dble(n)
		do j = 1, n
		   temp1 = 1.0d0
		   temp2 = 2.0d0*x(j) - 1.0d0
		   temp = 2.0d0*temp2
		   temp3 = 0.0d0
		   temp4 = 2.0d0
		   do k = 1, n
			  dfdx(k,j) = tk*temp4
			  ti = 4.0d0*temp2 + temp*temp4 - temp3
			  temp3 = temp4
			  temp4 = ti
			  ti = temp*temp2 - temp1
			  temp1 = temp2
			  temp2 = ti
		   enddo
		enddo
      endif
      return
      end
