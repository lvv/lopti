      subroutine SulpDiox(n, xa, fvec, fjac, ldfjac, task, p, np, ldx, 
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
c   Subroutine SulpDiox
c
c   This subroutine computes the function and Jacobian matrix of
c   the sulphur dioxide to sulphur trioxide problem.
c   (reduced formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             SulpDiox.f
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
      double precision arg,x(3)

      ifail = 0
c     Compute the standard starting point if task = 'XS'

      if (task .eq. 'XS') then
        n = 1
        nx0 = 2
        x0(1,1) = 0.4d0
        x0(1,2) = 0.0d0
        xlow(1,1) =  0.0d0
        xup (1,1) =  1.0d3
        xlow(1,2) =  0.0d0
        xup (1,2) =  1.0d0
        nbound = 2
        return
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (xa(1).eq.1.0d0) then
        ifail = 1
        return
      endif

      if (task .eq. 'F' .or. task .eq. 'FJ') then
        x(1) = xa(1)
        x(2) = (43260.0d0*x(1)+105128.0d0)/(1.84d0*x(1)+77.3d0)
		if (x(2).le.0.0d0) then
		  ifail = 1
		  return
		endif
		arg = 42300.0d0/x(2)-24.2d0+0.17*dlog(x(2))
		if (arg.gt.350.0d0) then
		  ifail = 1
		  return
		endif
        x(3) = dexp(arg)
		if (x(3).eq.0.0d0) then
		  ifail = 1
		  return
		endif
        fvec(1) = (0.91d0-0.5d0*x(1))/(9.1d0-0.5d0*x(1)) -
     $            x(1)**2/((1.0d0-x(1))**2*x(3))
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ifail = -999
      end if

      end
