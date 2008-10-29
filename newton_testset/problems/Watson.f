      subroutine Watson(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Watson function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Watson.f
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
      integer i, j, k
      double precision temp, temp1, temp2, ti, tj, tk, sum1, sum2
      double precision c9
      data c9 /2.9d1/
c
      ifail = 0
      if (task .eq. 'XS') then
        if (n.eq.0) n = 10
        do j = 1, n
          x0(j,1) = 0.0d0
        enddo
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		do k = 1, n
		   f(k) = 0.0d0
		enddo
		do i = 1, 29
		   ti = dfloat(i)/c9
		   sum1 = 0.0d0
		   temp = 1.0d0
		   do j = 2, n
			  sum1 = sum1 + dfloat(j-1)*temp*x(j)
			  temp = ti*temp
		   enddo
		   sum2 = 0.0d0
		   temp = 1.0d0
		   do j = 1, n
			  sum2 = sum2 + temp*x(j)
			  temp = ti*temp
		   enddo
		   temp1 = sum1 - sum2**2 - 1.0d0
		   temp2 = 2.0d0*ti*sum2
		   temp = 1.0d0/ti
		   do k = 1, n
			  f(k) = f(k) + temp*(dfloat(k-1) - temp2)*temp1
			  temp = ti*temp
		   enddo
		enddo
		temp = x(2) - x(1)**2 - 1.0d0
		f(1) = f(1) + x(1)*(1.0d0 - 2.0d0*temp)
		f(2) = f(2) + temp
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        do k = 1, n
           do j = k, n
              dfdx(k,j) = 0.0d0
           enddo
        enddo
        do i = 1, 29
           ti = dfloat(i)/c9
           sum1 = 0.0d0
           temp = 1.0d0
           do j = 2, n
              sum1 = sum1 + dfloat(j-1)*temp*x(j)
              temp = ti*temp
           enddo
           sum2 = 0.0d0
           temp = 1.0d0
           do j = 1, n
              sum2 = sum2 + temp*x(j)
              temp = ti*temp
           enddo
           temp1 = 2.0d0*(sum1 - sum2**2 - 1.0d0)
           temp2 = 2.0d0*sum2
           temp = ti**2
           tk = 1.0d0
           do k = 1, n
              tj = tk
              do j = k, n
                 dfdx(k,j) = dfdx(k,j)
     *                       + tj*((dfloat(k-1)/ti - temp2)
     *                           *(dfloat(j-1)/ti - temp2) - temp1)
                 tj = ti*tj
              enddo
              tk = temp*tk
           enddo
        enddo
        dfdx(1,1) = dfdx(1,1) + 6.0d0*x(1)**2 - 2.0d0*x(2) + 3.0d0
        dfdx(1,2) = dfdx(1,2) - 2.0d0*x(1)
        dfdx(2,2) = dfdx(2,2) + 1.0d0
        do k = 1, n
           do j = k, n
              dfdx(j,k) = dfdx(k,j)
           enddo
        enddo
      endif
      return
      end
