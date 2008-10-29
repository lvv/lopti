      subroutine InterSEllHy(n, xa, fvec, fjac, ldfjac, task, p, np,ldx,
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
c     Subroutine InterSEllHy
c
c     This subroutine computes the function and Jacobian matrix of
c     the intersection of an ellipsoid with a hyperboloid problem.
c     (reduced formulation)
c
c     Coded by         L. Weimann, 
c                      Zuse Institute Berlin (ZIB),
c                      Dep. Numerical Analysis and Modelling
c     File             InterSEllHy.f
c     Version          1.0
c     Latest Change    13th December 2005
c     Code             Fortran 77 with common extensions
c                      Double precision
c
c     Reference:
c       R. Luus:
c       Iterative Dynamic Programming.
c       Chapman & Hall/CRC, 2000
c       page 27f (example 5)
c
c-----------------------------------------------------------------------
c     Compute the standard starting point if task = 'XS'
      double precision x(7)
      integer i, j

      ifail = 0
      if (task .eq. 'XS') then
        x(1) = 1.0d0
        x(2) = 1.0d0
        x(3) = 1.0d0
        n = 3
        nx0 = 1
        do i=1,n
          x0(i,1) = x(i)
        enddo
        nbound = 1
        do i=1,n
          xlow(i,1) = -2.0d0
          xup (i,1) =  2.0d0
        enddo
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
        x(1) = xa(1)
        x(2) = xa(2)
        x(3) = xa(3)
        x(7) = x(3)/(2.0d0*x(3)+0.2*x(2))
        x(6) = 4.0d0*(x(2)-0.2d0)+0.1d0*x(1)+0.2d0*x(3)
        x(5) = (x(6)*x(7)-x(2))/(x(2)+2.0d0*x(6)*x(7))
        x(4) = (4.0d0*x(5)-2.0d0)*x(7)
        fvec(1) = 2.0d0*x(1)+x(4)*(8.0d0*(x(1)-0.5d0)+0.1*x(2))
     $            +4.0d0*x(1)*x(5)
        fvec(2) = 2.0d0*(x(1)**2-x(3)**2-1.0d0)+x(2)**2
        fvec(3) = 4.0d0*(x(1)-0.5d0)**2+2.0d0*(x(2)-0.2d0)**2+x(3)**2
     $            +0.1d0*x(1)*x(2)+0.2d0*x(2)*x(3)-16.0d0
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ifail = -999
      end if

      end
