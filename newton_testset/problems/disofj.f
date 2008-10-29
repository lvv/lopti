      subroutine disofj(n, x, fvec, fjac, ldfjac, task, p, np, ldx, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      implicit none
      character*(*) task
      integer n, ldfjac, np, nx0, ldx, nbound, ncon, nsol, ifail
      double precision x(*), fvec(*), fjac(ldfjac,*), p(*), x0(ldx,*),
     $                 xlow(ldx,*), xup(ldx,*), conxl(ldx,*),
     $                 conxu(ldx,*), xsol(ldx,*), xthrsh(*), fthrsh(*),
     $                 rtol
c-----------------------------------------------------------------------
c
c   Subroutine disofj
c
c   This subroutine computes the function and Jacobian matrix of
c   the - dew point temperature for a 20% isobutanol and 80%
c   water mixture - problem. (full formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             disofj.f
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
      double precision ln10in
      integer i, j
c
      ifail = 0
      if (task .eq. 'XS') then
        x(1) = 0.1d0
        x(2) = 0.9d0
        x(3) = 85.0d0
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
        n = 16
        nx0 = 2
        do i=1,n
          x0(i,1) = x(i)
        enddo
        x(1) = 0.0226
        x(2) = 0.686
        x(3) = 88.5
        x(4) = 0.977
        x(5) = 0.313
        x(6) = 357.05
        x(7) = 498.65
        x(8) = 33.36
        x(9) = 1.004
        x(10) = 1.102
        x(11) = 3.134
        x(12) = 15.675
        x(13) = 0.659
        x(14) = 0.518
        x(15) = 2.056
        x(16) = 0.733
        do i=1,n
          x0(i,2) = x(i)
        enddo
        do i=1,n
           xlow(i,1) =  0.0d0
           xup (i,1) =  5.0d1
        enddo
        nbound = 1
        return
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      do i=6,11
        if (x(i).le.0.0) then
          ifail = 1
          return
        endif
      enddo
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        fvec(1) = 1.0d0-x(1)-x(4)
        fvec(2) = 1.0d0-x(2)-x(5)
        fvec(3) = 7.62231d0-1417.90d0/(191.15d0+x(3))-dlog10(x(6))
        fvec(4) = 8.10765d0-1750.29d0/(235.00d0+x(3))-dlog10(x(7))
        fvec(5) = 1.7d0*x(4)**2/(1.7d0*x(1)/0.7d0+x(4))**2-dlog10(x(8))
        fvec(6) = 0.7d0*x(1)**2/(x(1)+0.7d0*x(4)/1.7d0)**2-dlog10(x(9))
        fvec(7) = 1.7d0*x(5)**2/(1.7d0*x(2)/0.7d0+x(5))**2-dlog10(x(10))
        fvec(8) = 0.7d0*x(2)**2/(x(2)+0.7d0*x(5)/1.7d0)**2-dlog10(x(11))
        fvec(9)  = x(6)*x(8)/760.0d0-x(12)
        fvec(10) = x(7)*x(9)/760.0d0-x(13)
        fvec(11) = x(10)*x(6)/760.0d0-x(14)
        fvec(12) = x(11)*x(7)/760.0d0-x(15)
        fvec(13) = (0.2d0/x(1)-x(12)/x(14))/(1.0d0-x(12)/x(14))-x(16)
        fvec(14) = x(12)*x(1)+x(13)*x(4)-1.0d0
        fvec(15) = x(2)-(x(1)*x(12))/x(14)
        fvec(16) = x(4)-0.8d0/(x(16)+(1.0d0-x(16))*x(13)/x(15))
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ln10in = 1.0d0/dlog(10.0d0)
        do 20 i=1,16
          do 20 j=1,16
            fjac(i,j) = 0.0d0
20      continue
        fjac(1,1) = -1.0d0
        fjac(1,4) = -1.0d0
        fjac(2,2) = -1.0d0
        fjac(2,5) = -1.0d0
        fjac(3,3) = 1417.90d0/(191.15d0+x(3))**2
        fjac(3,6) = -ln10in/x(6)
        fjac(4,3) = 1750.29d0/(235.00d0+x(3))**2
        fjac(4,7) = -ln10in/x(7)
        fjac(5,1) = -2.0d0*1.7d0/0.7d0*1.7d0*x(4)**2 /
     $              (1.7d0*x(1)/0.7d0+x(4))**3
        fjac(5,4) = 3.4d0*x(4)/(1.7d0*x(1)/0.7d0+x(4))**2 - 
     $              2.0d0*1.7d0*x(4)**2/(1.7d0*x(1)/0.7d0+x(4))**3
        fjac(5,8) = -ln10in/x(8)
        fjac(6,1) = 1.4d0*x(1)/(x(1)+0.7d0*x(4)/1.7d0)**2-
     $              0.7d0*x(1)**2*2.0d0/(x(1)+0.7d0*x(4)/1.7d0)**3
        fjac(6,4) = -2.0d0*0.7d0/1.7d0*0.7d0*x(1)**2 /
     $               (x(1)+0.7d0*x(4)/1.7d0)**3
        fjac(6,9) = -ln10in/x(9)
        fjac(7,5) = 3.4d0*x(5)/(1.7d0*x(2)/0.7d0+x(5))**2 - 
     $              3.4d0*x(5)**2/(1.7d0*x(2)/0.7d0+x(5))**3
        fjac(7,2) = -3.4d0/0.7d0*1.7d0*x(5)**2/
     $              (1.7d0*x(2)/0.7d0+x(5))**3
        fjac(7,10)  = -ln10in/x(10)
        fjac(8,2)   = 1.4d0*x(2)/(x(2)+0.7d0*x(5)/1.7d0)**2 -
     $                1.4d0*x(2)**2/(x(2)+0.7d0*x(5)/1.7d0)**3
        fjac(8,5)   = -2.0d0*0.7d0/1.7d0*0.7d0*x(2)**2/
     $                (x(2)+0.7d0*x(5)/1.7d0)**3
        fjac(8,11)  = -ln10in/x(11)
        fjac(9,6)   = x(8)/760.0d0
        fjac(9,8)   = x(6)/760.0d0
        fjac(9,12)  = -1.0d0
        fjac(10,7)  = x(9)/760.0d0
        fjac(10,9)  = x(7)/760.0d0
        fjac(10,13) = -1.0d0
        fjac(11,10) = x(6)/760.0d0
        fjac(11,6)  = x(10)/760.0d0
        fjac(11,14) = -1.0d0
        fjac(12,11) = x(7)/760.0d0
        fjac(12,7)  = x(11)/760.0d0
        fjac(12,15) = -1.0d0
        fjac(13,1)  = -0.2d0/x(1)**2/(1.0d0-x(12)/x(14))
        fjac(13,12) = 1.0d0/x(14)*(-(1.0d0-x(12)/x(14))+
     $                (0.2d0/x(1)-x(12)/x(14)))/(1.0d0-x(12)/x(14))**2
        fjac(13,14) = -x(12)/x(14)**2*(-(1.0d0-x(12)/x(14))+
     $                (0.2d0/x(1)-x(12)/x(14)))/(1.0d0-x(12)/x(14))**2
        fjac(13,16) = -1.0d0
        fjac(14,1)  = x(12)
        fjac(14,12) = x(1)
        fjac(14,13) = x(4)
        fjac(14,4)  = x(13)
        fjac(15,2)  = 1.0d0
        fjac(15,1)  = -x(12)/x(14)
        fjac(15,12) = -x(1)/x(14)
        fjac(15,14) = (x(1)*x(12))/x(14)**2
        fjac(16,4)  = 1.0d0
        fjac(16,16) = 0.8d0*(1.0d0-x(13)/x(15))/
     $                (x(16)+(1.0d0-x(16))*x(13)/x(15))**2
        fjac(16,13) = 0.8d0*(1.0d0-x(16))/x(15)/
     $                (x(16)+(1.0d0-x(16))*x(13)/x(15))**2
        fjac(16,15) = -0.8d0*(1.0d0-x(16))*x(13)/x(15)**2/
     $                (x(16)+(1.0d0-x(16))*x(13)/x(15))**2
      end if

      end
