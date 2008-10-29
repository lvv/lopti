      subroutine EsterificRea(n, xa, fvec, fjac, ldfjac, task, p,np,ldx,
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
c   Subroutine EsterificRea
c
c   This subroutine computes the function and Jacobian matrix of
c   the esterification reaction problem. (reduced formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             EsterificRea.f
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
      double precision x(9)
      integer i,j
      ifail = 0
c
      if (task .eq. 'XS') then
        x(1) = 20.0d0
        x(2) = 20.0d0
        n = 2
        nx0 = 1
        do i=1,n
          x0(i,1) = x(i)
        enddo
        nbound = 2
        do i=1,n
          xlow(i,1) = 1.0d-4
          xup (i,1) = 1.0d2
        enddo
        do i=1,n
          xlow(i,2) = 1.0d-4
          xup (i,2) = 1.0d3
        enddo
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
        x(1) = xa(1)
        x(2) = xa(2)
        x(9) = x(1)+x(2)
        x(3) = x(9)*(6.875d0*x(9)+x(2))/(5.797d0*x(9)+x(2))
        x(4) = 1.25d0*x(1)+2.25d0*x(2)-1.054d0*x(3)
        x(7) = 0.025714*x(3)-0.007346189*x(2)-0.2720539*x(4)
        x(5) = 12.128256*x(3)-x(2)-124.0d0*x(4)-1.75d0*x(7)-75.85949d0
        x(6) = 61.177d0-8.7614*x(3)-x(2)+91.0*x(4)-x(5)+x(7)
        x(8) = 1.0986*x(3)-x(2)-9.0d0*x(4)-x(5)-x(6)+x(7)
        fvec(1) = x(4)*(x(3)-x(2)-x(1))**2
     $            -2056.4*(0.0986d0*x(3)-x(4))**2
        fvec(2) = 5.5d0*(0.0986d0*x(3)-x(2)-x(5)
     $            -0.01d0*x(8))*(x(5)+x(6))-x(5)*x(8)
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ifail = -999
      end if

      end
