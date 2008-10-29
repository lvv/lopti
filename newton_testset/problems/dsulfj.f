      subroutine dsulfj(n, x, fvec, fjac, ldfjac, task, p, np, ldx, 
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
c   Subroutine dsulfj
c
c   This subroutine computes the function and Jacobian matrix of
c   the sulphur dioxide to sulphur trioxide problem.
c   (full formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             dsulfj.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     R. Luus:
c     Iterative Dynamic Programming.
c     Chapman & Hall/CRC, 2000
c     page 19ff (example 1)
c
c-----------------------------------------------------------------------
      integer i
      double precision arg

      ifail = 0
c     Compute the standard starting point if task = 'XS'

      if (task .eq. 'XS') then
        x(1) = 0.4d0
        x(2) = (43260.0d0*x(1)+105128.0d0)/(1.84d0*x(1)+77.3d0)
        x(3) = dexp(42300.0d0/x(2)-24.2d0+0.17*dlog(x(2)))
        n = 3
        nx0 = 1
        do i=1,n
          x0(i,1) = x(i)
        enddo
        do i=1,n
          xlow(i,1) =  0.0d0
          xup (i,1) =  1.0d3
        enddo
        xlow(1,2) = 0.0d0
        xup (1,2) = 1.0d0
        xlow(2,2) = 1000.0d0
        xup (2,2) = 2000.0d0
        xlow(3,2) = 0.0d0
        xup (3,2) = 100.0d0
        nbound = 2
        return
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (x(1).eq.1.0d0) then
        ifail = 1
        return
      endif
      if (x(2).le.0.0d0) then
        ifail = 1
        return
      endif
      if (x(3).eq.0.0d0) then
        ifail = 1
        return
      endif
      arg = 42300.0d0/x(2)-24.2d0+0.17*dlog(x(2))
      if (arg.gt.350.0d0) then
        ifail = 1
        return
      endif

      if (task .eq. 'F' .or. task .eq. 'FJ') then
        fvec(1) = (0.91d0-0.5d0*x(1))/(9.1d0-0.5d0*x(1)) -
     $            x(1)**2/((1.0d0-x(1))**2*x(3))
        fvec(2) = x(2)*(1.84d0*x(1)+77.3d0)-43260.0d0*x(1)-105128.0d0
        fvec(3) = dexp(arg)-x(3)
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        fjac(1,1) = -4.095/(9.1-0.5*x(1))**2 
     $             - 2.0d0*x(1)/((1.0d0-x(1))**3*x(3))
        fjac(1,2) = 0.0d0
        fjac(1,3) = x(1)**2 / ((1.0d0-x(1))**2*x(3)**2)
        fjac(2,1) = 1.84d0*x(2)-43260.0d0
        fjac(2,2) = 1.84d0*x(1)+77.3d0
        fjac(2,3) = 0.0d0
        fjac(3,1) = 0.0d0
        fjac(3,2) = dexp(arg)*(-42300.0d0/x(2)**2+0.17d0/x(2))
        fjac(3,3) = -1.0d0
      end if

      end
