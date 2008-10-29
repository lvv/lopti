      subroutine CRSteadyState5(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Chemical Reactor Steady State
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             CRSteadyState5.f
c   Version          1.0
c   Latest Change    8th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     M. Shacham:
c     Recent developments in solution techniques for systems
c     of nonlinear equations.
c     In A.W. Westerberg, H.H. Chien (eds.):
c     Proceedings of the second international conference
c     on foundations of computer-aided
c     process design. CACHE Publications, 1983
c
c-----------------------------------------------------------------------
      double precision y, T, k, dkdT, arg
c
      ifail = 0
      if (task .eq. 'XS') then
        n = 2
c       x0a
        x0(1,1) = 0.5
        x0(2,1) =320.0
c       x0b
        x0(1,2) = 0.0
        x0(2,2) = 300.0
c       x0c
        x0(1,3) = 0.0
        x0(2,3) = 350.0
c       x0d
        x0(1,4) = 1.0
        x0(2,4) = 400.0
c       very close to the solution
        x0(1,5) = 0.964
        x0(2,5) = 338.0
        x0(1,6) = 0.9
        x0(2,6) = 350.0
        x0(1,7) = 0.9
        x0(2,7) = 310.0
        x0(1,8) = 0.9
        x0(2,8) = 390.0
        x0(1,9) = 0.906948356846D+00
        x0(2,9) = 0.308103350833D+03
        nx0 = 9
        xthrsh(1) = 1.0d0
        xthrsh(2) = 1.0d0
        fthrsh(1) = 1.0d0
        fthrsh(2) = 1.0d0
        rtol = 1.0d-10
c       range of interest/solution
        xlow(1,1) = 0.0
        xlow(2,1) = 300.0
        xup (1,1) = 1.0
        xup (2,1) = 400.0
        nbound = 1
        return
      endif
      y = x(1)
      T = x(2)
      arg = 12581.0*(T-298.0)/(298.0*T)
      if (arg.gt.395.0d0) then
        ifail = 1
        return
      endif
      k = 0.12 * exp( arg )

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1) = 120.0*y - 75.0*k*(1.0-y)
        f(2) = -y*(873.0-T) + 11.0*(T-300.0)
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        dkdT = k * (12581.0 *298.0*T - 12581.0*(T-298.0)*298.0) /
     $        (298.0*T)**2
        dfdx(1,1) = 120.0 +75.0*k
        dfdx(1,2) = -75.0*dkdT*(1.0-y)
        dfdx(2,1) = -(873.0-T)
        dfdx(2,2) = y + 11.0
      endif
      return
      end

