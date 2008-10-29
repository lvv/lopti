      subroutine DewPoint(n, xa, fvec, fjac, ldfjac, task, p, np, ldx, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      implicit none
      character*(*) task
      integer n, ldfjac, np, nx0, ldx, nbound, ncon, nsol, ifail
      double precision xa(*), fvec(*), fjac(ldfjac,*), p(*), x0(ldx,*),
     $                 xlow(ldx,*), xup(ldx,*), conxl(ldx,*),
     $                 conxu(ldx,*), xsol(ldx,*), xthrsh(*), fthrsh(*),
     $                 rtol
c-----------------------------------------------------------------------
c
c   Subroutine DewPoint
c
c   This subroutine computes the function and Jacobian matrix of
c   the - dew point temperature for a 20% isobutanol and 80%
c   water mixture - problem. (reduced formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             DewPoint.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     R. Luus:
c     Iterative Dynamic Programming.
c     Chapman & Hall/CRC, 2000
c     page 28ff (example 6)
c
c-----------------------------------------------------------------------
c     Compute the standard starting point if task = 'XS'
      double precision x(16)
      integer i, j
      logical qprint
      common /dewpri/ qprint
c
      ifail = 0
      if (task .eq. 'XS') then
        x(1) = 0.1d0
        x(2) = 0.9d0
        x(3) = 85.0d0
        n = 3
        nx0 = 1
        do i=1,n
          x0(i,1) = x(i)
        enddo
        xlow(1,1) =  0.01d0
        xup (1,1) =  0.99d0
        xlow(2,1) =  0.01d0
        xup (2,1) =  0.99d0
        xlow(3,1) =  80.0d0
        xup (3,1) =  90.0d0
        nbound = 1
        return
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
        x(1) = xa(1)
        x(2) = xa(2)
        x(3) = xa(3)
        x(4) = 1.0d0-x(1)
        x(5) = 1.0d0-x(2)
        x(6) = 10.0d0**(7.62231d0-1417.90d0/(191.15d0+x(3)))
        x(7) = 10.0d0**(8.10765d0-1750.29d0/(235.00d0+x(3)))
        x(8) = 10.0d0**(1.7d0*x(4)**2/(1.7d0*x(1)/0.7d0+x(4))**2)
        x(9) = 10.0d0**(0.7d0*x(1)**2/(x(1)+0.7d0*x(4)/1.7d0)**2)
        x(10) = 10.0d0**(1.7d0*x(5)**2/(1.7d0*x(2)/0.7d0+x(5))**2)
        x(11) = 10.0d0**(0.7d0*x(2)**2/(x(2)+0.7d0*x(5)/1.7d0)**2)
        x(12) = x(6)*x(8)/760.0d0
        x(13) = x(7)*x(9)/760.0d0
        x(14) = x(10)*x(6)/760.0d0
        x(15) = x(11)*x(7)/760.0d0
        x(16) = (0.2d0/x(1)-x(12)/x(14))/(1.0d0-x(12)/x(14))
        if ( x(12) .eq. x(14) ) then
          ifail = 1
          return
        endif
        if (qprint) then
          write(20,'(a,25(4(1x,d18.10),/))') 'x00 =',(x(i),i=1,16)
        endif
        fvec(1) = x(12)*x(1)+x(13)*x(4)-1.0d0
        fvec(2) = x(2)-(x(1)*x(12))/x(14)
        fvec(3) = x(4)-0.8d0/(x(16)+(1.0d0-x(16))*x(13)/x(15))
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ifail = -999
      end if

      end
      
      block data dewcom
      logical qprint
      common /dewpri/ qprint
      data qprint /.false./
      end
            
