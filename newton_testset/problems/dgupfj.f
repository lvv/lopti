      subroutine dgupfj(n, x, fvec, fjac, ldfjac, task, p, np, ldx, 
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
c   Subroutine dgupfj
c
c   This subroutine computes the function and Jacobian matrix of
c   the Gupta problem. (full formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             dgupfj.f
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
      double precision s
      integer i,j

      ifail = 0
c     Compute the standard starting point if task = 'XS'

      if (task .eq. 'XS') then
        x(1) = 0.1d0
        x(2) = 0.1d0
        x(11) = x(1)/3.846d-5
        x(10) = 3.0d0-x(2)
        x(6) = 2.597d-3*dsqrt(x(2)*x(11))
        x(7) = 3.448d-3*dsqrt(x(10)*x(11))
        x(4) = (4.0d0-0.5d0*(x(6)+x(7)))/(1.0d0+0.193d0*x(2)/x(10))
        x(5) = 0.193d0*x(2)*x(4)/x(10)
        x(12) = 40.0d0
        x(3) = 2.0d0*x(12)-0.5d0*x(5)
        x(8) = (1.799d0*x(4)*x(1))/(3.846d0*x(2))
        x(9) = 2.155d-4*x(10)*dsqrt(x(3)*x(11))/x(2)
        n = 12
        nx0 = 1
        do i=1,n
          x0(i,1) = x(i)
        enddo
         do i=1,n
           xlow(i,1) =  0.0d0
           xup (i,1) =  1.0d1
         enddo
         nbound = 1
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
        if (x(2)*x(11).le.0.0d0 .or. x(10)*x(11).le.0.0d0 .or.
     $      x(3)*x(11).le.0.0d0) then
          ifail = 1
          return
        endif
        fvec(1) = x(1)/3.846d-5-x(11)
        fvec(2) = 3.0d0-x(2)-x(10)
        fvec(3) = 2.597d-3*dsqrt(x(2)*x(11))-x(6)
        fvec(4) = 3.448d-3*dsqrt(x(10)*x(11))-x(7)
        fvec(5) = (4.0d0-0.5d0*(x(6)+x(7)))/(1.0d0+0.193d0*x(2)/x(10))-
     $            x(4)
        fvec(6) = 0.193d0*x(2)*x(4)/x(10)-x(5)
        fvec(7) = 40.0d0-x(12)
        fvec(8) = 2.0d0*x(12)-0.5d0*x(5)-x(3)
        fvec(9) = (1.799d0*x(4)*x(1))/(3.846d0*x(2))-x(8)
        fvec(10) = 2.155d-4*x(10)*dsqrt(x(3)*x(11))/x(2)-x(9)
        fvec(11) = 2.0d0*(x(1)+x(10))+x(2)+x(4)+x(7)+x(8)+x(9)-x(12)
        s = 0.0d0
        do 10 i=1,10
          s = s + x(i)
10      continue
        fvec(12) = s-x(11)
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        do 20 i=1,12
          do 20 j=1,12
            fjac(i,j) = 0.0d0
20      continue
        fjac(1,1)  = 1.0d0/3.846d-5
        fjac(1,11) = -1.0d0
        fjac(2,2)  = -1.0d0
        fjac(2,10) = -1.0d0
        fjac(3,2)  = 2.597d-3*0.5d0*x(11)/dsqrt(x(2)*x(11))
        fjac(3,6)  = -1.0d0
        fjac(3,11) = 2.597d-3*0.5d0*x(2)/dsqrt(x(2)*x(11))
        fjac(4,7)  = -1.0d0
        fjac(4,10) = 3.448d-3*0.5d0*x(11)/dsqrt(x(10)*x(11))
        fjac(4,11) = 3.448d-3*0.5d0*x(10)/dsqrt(x(10)*x(11))
        fjac(5,2)  = - (4.0d0-0.5d0*(x(6)+x(7)))/
     $                (1.0d0+0.193d0*x(2)/x(10))**2 * 0.193d0/x(10)
        fjac(5,6)  = -0.5d0/(1.0d0+0.193d0*x(2)/x(10))
        fjac(5,7)  = -0.5d0/(1.0d0+0.193d0*x(2)/x(10))
        fjac(5,10) =  (4.0d0-0.5d0*(x(6)+x(7)))/
     $             (1.0d0+0.193d0*x(2)/x(10))**2 * 0.193d0*x(2)/x(10)**2
        fjac(5,4)  = -1.0d0
        fjac(6,2)  = 0.193d0*x(4)/x(10)
        fjac(6,4)  = 0.193d0*x(2)/x(10)
        fjac(6,10) = -0.193d0*x(2)*x(4)/x(10)**2
        fjac(6,5)  = -1.0d0
        fjac(7,12) = -1.0d0
        fjac(8,12) = 2.0d0
        fjac(8,5)  = -0.5d0
        fjac(8,3)  = -1.0d0
        fjac(9,1)  = (1.799d0*x(4))/(3.846d0*x(2))
        fjac(9,4)  = (1.799d0*x(1))/(3.846d0*x(2))
        fjac(9,2)  = -3.846d0*(1.799d0*x(4)*x(1))/(3.846d0*x(2))**2
        fjac(9,8)  = -1.0d0
        fjac(10,10)= 2.155d-4*dsqrt(x(3)*x(11))/x(2)
        fjac(10,3) = 2.155d-4*x(10)*0.5d0*x(11)/dsqrt(x(3)*x(11))/x(2)
        fjac(10,11)= 2.155d-4*x(10)*0.5d0*x(3)/dsqrt(x(3)*x(11))/x(2)
        fjac(10,2) = -2.155d-4*x(10)*dsqrt(x(3)*x(11))/x(2)**2
        fjac(10,9) = -1.0d0
        fjac(11,1) = 2.0d0
        fjac(11,10)= 2.0d0
        fjac(11,2) = 1.0d0
        fjac(11,4) = 1.0d0
        fjac(11,7) = 1.0d0
        fjac(11,8) = 1.0d0
        fjac(11,9) = 1.0d0
        fjac(11,12) = -1.0d0
        do 30 i=1,10
          fjac(12,i) = 1.0d0
30      continue
        fjac(12,11) = -1.0d0
      end if

      end
