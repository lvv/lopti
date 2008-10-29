      subroutine Wood(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Wood function.
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Wood.f
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
      integer j,k
      double precision temp1, temp2
      double precision c3, c4, c5, c6
      data c3, c4, c5, c6 /2.0d2,2.02d1,1.98d1,1.8d2/
c
      ifail = 0
      if (task .eq. 'XS') then
        n = 4
		x0(1,1) = -3.0d0
		x0(2,1) = -1.0d0
		x0(3,1) = -3.0d0
		x0(4,1) = -1.0d0
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        temp1 = x(2) - x(1)**2
        temp2 = x(4) - x(3)**2
        f(1) = -c3*x(1)*temp1 - (1.0d0 - x(1))
        f(2) = c3*temp1 + c4*(x(2) - 1.0d0) + c5*(x(4) - 1.0d0)
        f(3) = -c6*x(3)*temp2 - (1.0d0 - x(3))
        f(4) = c6*temp2 + c4*(x(4) - 1.0d0) + c5*(x(2) - 1.0d0)
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        do k = 1, 4
          do j = 1, 4
            dfdx(k,j) = 0.0d0
          enddo
        enddo
        temp1 = x(2) - 3.0d0*x(1)**2
        temp2 = x(4) - 3.0d0*x(3)**2
        dfdx(1,1) = -c3*temp1 + 1.0d0
        dfdx(1,2) = -c3*x(1)
        dfdx(2,1) = -2.0d0*c3*x(1)
        dfdx(2,2) = c3 + c4
        dfdx(2,4) = c5
        dfdx(3,3) = -c6*temp2 + 1.0d0
        dfdx(3,4) = -c6*x(3)
        dfdx(4,2) = c5
        dfdx(4,3) = -2.0d0*c6*x(3)
        dfdx(4,4) = c6 + c4
      endif
      return
      end
