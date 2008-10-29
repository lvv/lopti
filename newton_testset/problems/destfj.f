      subroutine desti(n,x)
      integer n
      double precision x(n)
c
c     x(1) and x(2) must have been set prior calling this 
c     starting point (initialization) routine
        x(9) = x(1)+x(2)
        x(3) = x(9)*(6.875d0*x(9)+x(2))/(5.797d0*x(9)+x(2))
        x(4) = 1.25d0*x(1)+2.25d0*x(2)-1.054d0*x(3)
        x(7) = 0.025714*x(3)-0.007346189*x(2)-0.2720539*x(4)
        x(5) = 12.128256*x(3)-x(2)-14.0d0*x(4)-1.75d0*x(7)-75.85949d0
        x(6) = 61.177d0-8.7614*x(3)-x(2)+91.0*x(4)-x(5)+x(7)
        x(8) = 1.0986*x(3)-x(2)-9.0d0*x(4)-x(5)-x(6)+x(7)
      return
      end
c
      subroutine destfj(n, x, fvec, fjac, ldfjac, task, p, np, ldx, 
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
c   Subroutine destfj
c
c   This subroutine computes the function and Jacobian matrix of
c   the esterification reaction problem. (full formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             destfj.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     R. Luus:
c     Iterative Dynamic Programming.
c     Chapman & Hall/CRC, 2000
c     page 25ff (example 4)
c
c-----------------------------------------------------------------------
c     Compute the standard starting point if task = 'XS'
      integer i,j
      ifail = 0
c
      if (task .eq. 'XS') then
        x(1) = 20.0d0
        x(2) = 20.0d0
        x(9) = x(1)+x(2)
        x(3) = x(9)*(6.875d0*x(9)+x(2))/(5.797d0*x(9)+x(2))
        x(4) = 1.25d0*x(1)+2.25d0*x(2)-1.054d0*x(3)
        x(7) = 0.025714*x(3)-0.007346189*x(2)-0.2720539*x(4)
        x(5) = 12.128256*x(3)-x(2)-124.0d0*x(4)-1.75d0*x(7)-75.85949d0
        x(6) = 61.177d0-8.7614*x(3)-x(2)+91.0*x(4)-x(5)+x(7)
        x(8) = 1.0986*x(3)-x(2)-9.0d0*x(4)-x(5)-x(6)+x(7)
        n = 9
        nx0 = 1
        do i=1,n
          x0(i,1) = x(i)
        enddo
        nbound = 1
        do i=1,n
          xlow(i,1) = -40.0d0
          xup (i,1) =  40.0d0
        enddo
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
        fvec(1) = x(1)+x(2)-x(9)
        fvec(2) = x(9)*(6.875d0*x(9)+x(2))/(5.797d0*x(9)+x(2))-x(3)
        fvec(3) = 1.25d0*x(1)+2.25d0*x(2)-1.054d0*x(3)-x(4)
        fvec(4) = 0.025714*x(3)-0.007346189*x(2)-0.2720539*x(4)-x(7)
        fvec(5) = 12.128256*x(3)-x(2)-124.0d0*x(4)-1.75d0*x(7)
     $            -75.85949d0-x(5)
        fvec(6) = 61.177d0-8.7614*x(3)-x(2)+91.0*x(4)-x(5)+x(7)-x(6)
        fvec(7) = 1.0986*x(3)-x(2)-9.0d0*x(4)-x(5)-x(6)+x(7)-x(8)
        fvec(8) = x(4)*(x(3)-x(2)-x(1))**2
     $            -2056.4*(0.0986d0*x(3)-x(4))**2
        fvec(9) = 5.5d0*(0.0986d0*x(3)-x(2)-x(5)
     $            -0.01d0*x(8))*(x(5)+x(6))-x(5)*x(8)
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        do 20 i=1,9
          do 20 j=1,9
            fjac(i,j) = 0.0d0
20      continue
        fjac(1,1) = 1.0d0
        fjac(1,2) = 1.0d0
        fjac(1,9) = -1.0d0
        fjac(2,2) = (x(9)*(5.797d0*x(9)+x(2))-x(9)*(6.875d0*x(9)+x(2)))/
     $              (5.797d0*x(9)+x(2))**2
        fjac(2,9) = (2.0d0*6.875d0*x(9)*(5.797d0*x(9)+x(2))-
     $              x(9)*(6.875d0*x(9)+x(2))*5.797d0)/
     $              (5.797d0*x(9)+x(2))**2
        fjac(2,3) = -1.0d0
        fjac(3,1) = 1.25d0
        fjac(3,2) = 2.25d0
        fjac(3,3) = -1.054d0
        fjac(3,4) = -1.0d0
        fjac(4,3) = 0.025714
        fjac(4,2) = -0.007346189
        fjac(4,4) = -0.2720539
        fjac(4,7) = -1.0d0
        fjac(5,3) = 12.128256
        fjac(5,2) = -1.0d0
        fjac(5,4) = -124.0d0
        fjac(5,7) = -1.75d0
        fjac(5,5) = -1.0d0
        fjac(6,3) = -8.7614
        fjac(6,2) = -1.0d0
        fjac(6,4) = 91.0
        fjac(6,5) = -1.0d0
        fjac(6,7) = 1.0d0
        fjac(6,6) = -1.0d0
        fjac(7,3) = 1.0986
        fjac(7,2) = -1.0d0
        fjac(7,4) = -9.0d0
        fjac(7,5) = -1.0d0
        fjac(7,6) = -1.0d0
        fjac(7,7) = 1.0d0
        fjac(7,8) = -1.0d0
        fjac(8,1) = -2.0d0*x(4)*(x(3)-x(2)-x(1))
        fjac(8,2) = -2.0d0*x(4)*(x(3)-x(2)-x(1))
        fjac(8,3) = 2.0d0*x(4)*(x(3)-x(2)-x(1))-
     $              0.0986d0*2.0d0*2056.4*(0.0986d0*x(3)-x(4))
        fjac(8,4) = (x(3)-x(2)-x(1))**2+
     $              2.0d0*2056.4*(0.0986d0*x(3)-x(4))
        fjac(9,3) = 5.5d0*0.0986d0*(x(5)+x(6))
        fjac(9,2) = -5.5d0*(x(5)+x(6))
        fjac(9,5) = -5.5d0*(x(5)+x(6))+5.5d0*(0.0986d0*x(3)-x(2)-x(5)
     $            -0.01d0*x(8))-x(8)
        fjac(9,8) = -5.5d0*0.01d0*(x(5)+x(6))-x(5)
        fjac(9,6) = 5.5d0*(0.0986d0*x(3)-x(2)-x(5)-0.01d0*x(8))
      end if

      end
