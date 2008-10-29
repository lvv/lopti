      subroutine CREquiConv4(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Chemical Reactor Equilibrium Conversion
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             CREquiConv4.f
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
      double precision y, T, k1, k2
c
      ifail = 0
      if (task .eq. 'XS') then
        n = 2
c       close to solution, but large initial residual
        x0(1,1) = 0.5
        x0(2,1) = 1700.0
c       hard to converge due to singular points
c       x0a
        x0(1,2) = 0.0
        x0(2,2) = 1600.0
c       x0c
        x0(1,3) = 0.9
        x0(2,3) = 1600.0
c       x0b
        x0(1,4) = 0.0
        x0(2,4) = 1650.0
c       x0d
        x0(1,5) = 0.9
        x0(2,5) = 1700.0
c
c       Luus I.V.
        x0(1,6) = 0.0
        x0(2,6) = 1360.0
        x0(1,7) = 0.0
        x0(2,7) = 0.0
        x0(1,8) = 0.0
        x0(2,8) = 1650.0
        nx0 = 8
        xthrsh(1) = 1.0d0
        xthrsh(2) = 1.0d0
        fthrsh(1) = 1.0d0
        fthrsh(2) = 1.0d0
        rtol = 1.0d-8
        xlow(1,1) = 0.0d0
        xlow(2,1) = 1500.0d0
        xup (1,1) = 1.0d0
        xup (2,1) = 2000.0d0
        xlow(1,2) = 0.0d0
        xlow(2,2) = 1600.0d0
        xup (1,2) = 1.0d0
        xup (2,2) = 1700.0d0
        nbound = 2
        return
      endif

      y = x(1)
      if (y.gt.1.0d0) then
        ifail=1
        return
      endif
      T = x(2)
      if (T.le.0.0d0) then
        ifail=1
        return
      endif
      k1 =  -149750.0/T + 92.5
      k2 = 42300.0/T - 24.2 + 0.17*dlog(T)
      if ( dmax1(k1,k2).gt.350.0d0 ) then
        ifail=1
        return
      endif
      k1 = dexp( k1 )
      k2 = dexp( k2 )
      
c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1) = k1*sqrt(1.0-y) * ( (0.91-0.5*y)/(9.1 - 0.5*y)
     $                               -  y**2/((1-y)**2*k2) )     
        f(2) = T*(1.84*y+77.3) - 43260.0*y - 105128.0
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        dfdx( 1, 1) = -k1*0.5/sqrt(1.0-y) *( (0.91-0.5*y)/(9.1 - 0.5*y)
     $                               -  y**2/((1-y)**2*k2) ) +
     $                 k1*sqrt(1.0-y) * ( -4.095/(9.1-0.5*y)**2 -
     $                                     2.0*y/((1-y)**3*k2) )
        dfdx( 1, 2) = k1*149750.0/T**2 * sqrt(1.0-y) * ( 
     $               (0.91-0.5*y)/(9.1 - 0.5*y)  - y**2/((1-y)**2*k2) ) 
     $                 + k1*sqrt(1-y) *  y**2/((1-y)**2*k2)*
     $                 (-42300.0/T**2 +   0.17/T) 
        dfdx( 2,1 ) =  1.84*T - 43260.0
        dfdx( 2,2 ) =  1.84*y + 77.3
      endif
      return
      end
