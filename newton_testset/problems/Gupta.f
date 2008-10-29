      subroutine Gupta(n, xa, fvec, fjac, ldfjac, task, p, np, ldx, 
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
c   Subroutine Gupta
c
c   This subroutine computes the function and Jacobian matrix of
c   the Gupta problem. (reduced formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Gupta.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     R. Luus:
c     Iterative Dynamic Programming.
c     Chapman & Hall/CRC, 2000
c     page 21ff (example 2)
c
c-----------------------------------------------------------------------
      double precision s, x(12)
      integer i,j

      ifail = 0
c     Compute the standard starting point if task = 'XS'

      if (task .eq. 'XS') then
         n = 2
         nx0 = 2
         x0(1,1) = 0.1d0
         x0(2,1) = 0.1d0
         x0(1,2) = 0.5d0
         x0(2,2) = 0.5d0
         do i=1,n
           xlow(i,1) =  1.0d-4
           xup (i,1) =  1.0d1
         enddo
         nbound = 1
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
        x(1) = xa(1)
        x(2) = xa(2)
        x(11) = x(1)/3.846d-5
        x(10) = 3.0d0-x(2)
        x(6) = 2.597d-3*dsqrt(x(2)*x(11))
        x(7) = 3.448d-3*dsqrt(x(10)*x(11))
        x(4) = (4.0d0-0.5d0*(x(6)+x(7)))/(1.0d0+0.193d0*x(2)/x(10))
        x(5) = 0.193d0*x(2)*x(4)/x(10)
        x(12) = 40.0d0
        x(3) = 2.0d0*x(12)-0.5d0*x(5)
        if (x(2)*x(11).le.0.0d0 .or. x(10)*x(11).le.0.0d0 .or.
     $      x(3)*x(11).le.0.0d0) then
          ifail = 1
          return
        endif
        x(8) = (1.799d0*x(4)*x(1))/(3.846d0*x(2))
        x(9) = 2.155d-4*x(10)*dsqrt(x(3)*x(11))/x(2)
        fvec(1) = 2.0d0*(x(1)+x(10))+x(2)+x(4)+x(7)+x(8)+x(9)-x(12)
        s = 0.0d0
        do 10 i=1,10
          s = s + x(i)
10      continue
        fvec(2) = s-x(11)
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ifail = -999
      end if

      end
