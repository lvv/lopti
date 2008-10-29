      subroutine ChemEqui3(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Chemical equilibrium 3
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             ChemEqui3.f
c   Version          1.0
c   Latest Change    7th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     K.L. Hiebert:
c     An evaluation of mathematical software that
c     solves systems of nonlinear equations.
c     ACM Trans. Math. Soft. 8, 5-20, 1982
c
c-----------------------------------------------------------------------
      integer i,j, ivar, ntemp
      double precision rrr
      save rrr
      double precision xhh, xhh2, tot, tmp, temp, t5, t6, t7, t8, t9,
     $                 t10, xh(10), dmaxex
c
      ifail = 0
      if (np.eq.0) then
        ivar = 0
        rrr = 10.0d0
        xhh = 1.0d0
        p(1) = dble(ivar)
        p(2) = rrr
        p(3) = xhh
        np = 3
      else
        ivar = int(p(1))
        rrr  = p(2)
        xhh  = p(3)
      endif
      if (task .eq. 'XS') then
        n = 10
        if (ivar.eq.0 .or. ivar.eq.2) then
          do i=1,n
            x0(i,1)=xhh 
          enddo
          do i=1,n
            xlow(i,1) =  0.0d0
            xup (i,1) =  1.0d0
          enddo
        else if (ivar.eq.1) then
          xhh2 = dlog(xhh)
          do i=1,n
            x0(i,1)=xhh2 
          enddo
          do i=1,n
            xlow(i,1) = -1.0d0
            xup (i,1) =  1.0d0
          enddo
        else if (ivar.eq.3) then
          xhh2=dsqrt(xhh)
          do i=1,n
            x0(i,1)=xhh2 
          enddo
          do i=1,n
            xlow(i,1) = -1.0d0
            xup (i,1) =  1.0d0
          enddo
        endif
        nx0 = 1
		do i=1,n
		  xlow(i,1) =  1.0d-3
		  xup (i,1) =  1.0d1
		enddo
        nbound = 1
        return
      endif

      T5  = 0.193D+0
      T6  = 2.597*1.0D-03
      T7  = 3.448*1.0D-03
      T8  = 1.799*1.0D-05
      T9  = 2.155*1.0D-04
      T10 = 3.846*1.0D-05
	  TOT = 0.0d0
	  ntemp = 0
	  DO I= 1,n
		 xh(i)=x(i)
		 if (ivar.eq.0 .or. ivar.eq.1) then
		   dmaxex = 350.0d0
		   if (x(i).gt.dmaxex) then
			 ifail = 1
			 goto 1000
		   endif
		 endif
		 ntemp = i
		 if (ivar.eq.0 .or. ivar.eq.1) x(i)=dexp(x(i))
		 if (ivar.eq.3) x(i)=x(i)*x(i)
		 TOT = TOT + X(I)
	  enddo
c
c       if (ivar.eq.2) then
          if ( X(2)*X(4)*TOT.le.0.0d0 .or. X(2).le.0.0d0 .or.
     $         X(1)*X(4)*TOT.le.0.0d0 .or. X(4).le.0.0d0 .or.
     $         X(3)*TOT.le.0.0d0                              ) then
             ifail = 1
             goto 1000
          endif
c       endif
C
c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
c
        f(1) = -3.0d+0 + x(1) + x(4)
        f(2) =-rrr+2.0d0*x(1)+x(2)+x(4)+x(7)+x(8)+x(9)+2.0d0*x(10)
        f(3) = -8.0d+0 + 2.0d0*x(2) + 2.0d0*x(5) + x(6) + x(7)
        f(4) = -4.0d+0*rrr + 2.0d0*x(3) + x(5)
        f(5) = x(1)*x(5) - t5* x(2)*x(4)
        f(6) = dsqrt( x(2) )* x(6) - t6*dsqrt( x(2)*x(4) * tot)
        f(7) = dsqrt( x(4) )* x(7) - t7*dsqrt( x(1)*x(4) * tot)
        f(8) = x(4) * x(8) - t8 * x(2) * tot
        f(9) = x(4) * x(9) - t9 * x(1) * dsqrt( x(3) * tot)
        f(10)= (x(4)**2) * x(10) - t10 *( x(4)**2) * tot
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        do i = 1, n
           do j = 1, n
              dfdx(j,i) = 0.0d0
           enddo
        enddo
c
c  f(1) = -3.0d+0 + x(1) + x(4)
        dfdx(1,1) = 1.0d0
        dfdx(1,4) = 1.0d0
c
c  f(2) = -r + 2.0d0*x(1) + x(2) + x(4) + x(7) + x(8)+ x(9) + 2.0d0*x(10)
        dfdx(2,1) = 2.0d0
        dfdx(2,2) = 1.0d0
        dfdx(2,10) = 2.0d0
        dfdx(2,4) = 1.0d0
        dfdx(2,7) = 1.0d0
        dfdx(2,8) = 1.0d0
        dfdx(2,9) = 1.0d0
c
c  f(3) = -8.0d+0 + 2.0d0*x(2) + 2.0d0*x(5) + x(6) + x(7)
        dfdx(3,2) = 2.0d0
        dfdx(3,5) = 2.0d0
        dfdx(3,6) = 1.0d0
        dfdx(3,7) = 1.0d0
c
c  f(4) = -4.0d+0*r + 2.0d0*x(3) + x(5)
        dfdx(4,3) = 2.0d0
        dfdx(4,5) = 1.0d0
c
c  f(5) = x(1)*x(5) - t5* x(2)*x(4)
        dfdx(5,1) = x(5)
        dfdx(5,5) = x(1)
        dfdx(5,2) = -t5*x(4)
        dfdx(5,4) = -t5*x(2)
c
c  f(6) = dsqrt( x(2) )* x(6) - t6*dsqrt( x(2)*x(4) * tot)
		tmp  = 0.5d0*t6 / dsqrt( x(2)*x(4)*tot )
		temp = tmp *x(2)*x(4)
		dfdx(6,1) = -temp
		dfdx(6,2) = 0.5d0*x(6) / dsqrt(x(2)) -tmp*x(4)*(tot + x(2) )
		dfdx(6,3) = -temp
		dfdx(6,4) = -tmp * x(2) * (tot + x(4) )
		dfdx(6,5) = -temp
		dfdx(6,6) = -temp + dsqrt( x(2) )
		dfdx(6,7) = -temp
		dfdx(6,8) = -temp
		dfdx(6,9) = -temp
		dfdx(6,10) =-temp
c
c  f(7) = dsqrt( x(4) )* x(7) - t7*dsqrt( x(1)*x(4) * tot)
		tmp  = 0.5d0*t7 / dsqrt( x(1)*x(4)*tot )
		temp = tmp*x(1)*x(4)
		dfdx(7,1) = -tmp*x(4)* (tot + x(1) )
		dfdx(7,2) = -temp
		dfdx(7,3) = -temp
		dfdx(7,4) = 0.5d0*x(7) / dsqrt(x(4)) -tmp*x(1)*(tot+x(4) )
		dfdx(7,5) = -temp
		dfdx(7,6) = -temp
		dfdx(7,7) = -temp + dsqrt( x(4) )
		dfdx(7,8) = -temp
		dfdx(7,9) = -temp
		dfdx(7,10) = -temp
c
c  f(8) = x(4) * x(8) - t8 * x(2) * tot
		temp = t8 * x(2)
		dfdx(8,1) = -temp
		dfdx(8,2) = -(temp + t8 * tot)
		dfdx(8,3) = -temp
		dfdx(8,4) = x(8) - temp
		dfdx(8,5) = -temp
		dfdx(8,6) = -temp
		dfdx(8,7) = -temp
		dfdx(8,8) = x(4) - temp
		dfdx(8,9) = -temp
		dfdx(8,10) = -temp
c
c  f(9) = x(4) * x(9) - t9 * x(1) * dsqrt( x(3) * tot)
		tmp = 0.5d0* t9 *x(1) / dsqrt( x(3)*tot )
	    temp= tmp*x(3)
		dfdx(9,1) = -t9*dsqrt( x(3)*tot ) - temp
		dfdx(9,2) = -temp
		dfdx(9,3) = -tmp*( tot + x(3) )
		dfdx(9,4) = x(9) - temp
		dfdx(9,5) = -temp
		dfdx(9,6) = -temp
		dfdx(9,7) = -temp
		dfdx(9,8) = -temp
		dfdx(9,9) = x(4) - temp
		dfdx(9,10) = -temp
c
c  f(10)= (x(4)**2) * x(10) - t10 *( x(4)**2) * tot
		temp = t10*( x(4)**2 )
		dfdx(10,1) = -temp
		dfdx(10,2) = -temp
		dfdx(10,3) = -temp
		dfdx(10,4) = 2.0d0*x(4)*x(10) - temp -2.0d0*x(4)*t10*tot
		dfdx(10,5) = -temp
		dfdx(10,6) = -temp
		dfdx(10,7) = -temp
		dfdx(10,8) = -temp
		dfdx(10,9) = -temp
		dfdx(10,10) = x(4)**2 - temp
        if (ivar.eq.0 .or. ivar.eq.1) then
          do i = 1, 10
             do j = 1, 10
                dfdx(j,i) = dfdx(j,i)*x(i)
             enddo
          enddo
        else if (ivar.eq.3) then
          do i = 1, 10
             do j = 1, 10
                dfdx(j,i) = dfdx(j,i)*2.0d0*xh(i)
             enddo
          enddo
        endif
c
      endif
1000  continue
      do i= 1,ntemp
         x(i)=xh(i)
      enddo
      return
      end
