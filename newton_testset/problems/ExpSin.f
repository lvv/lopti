      subroutine ExpSin(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Exponential sine
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             ExpSin.f
c   Version          1.0
c   Latest Change    9th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     U. Nowak, L. Weimann:
c     A family of Newton codes for systems of highly nonlinear
c     equations.
c     Technical Report TR91-10, ZIB, 1991
c    
c-----------------------------------------------------------------------
      double precision arg, dmaxex
c
      ifail = 0
      if (np.eq.0) then
        p(1) = 3.0d0
        p(2) = 3.0d0
        np = 2
      endif
      if (task .eq. 'XS') then
        n = 2
        x0(1,1)=0.81d0  
        x0(2,1)=0.82d0  
        x0(1,2)= 2.99714682530400d0
        x0(2,2)=-2.95330710241260d0    
        x0(1,3)= 2.44169014107600d0
        x0(2,3)=-2.51379895001340d0    
        x0(1,4)= 2.51913768310200d0
        x0(2,4)=-2.85296012826840d0    
        x0(1,5)=-0.226027327648D+01  
        x0(2,5)=-0.241295786268D+01
        nx0 = 5
        xlow(1,1) = -1.5d0
        xlow(2,1) = -1.5d0
        xup (1,1) =  1.5d0
        xup (2,1) =  1.5d0
        xlow(1,2) = -3.0d0
        xlow(2,2) = -3.0d0
        xup (1,2) =  3.0d0
        xup (2,2) =  3.0d0
        nbound = 2
        return
      endif
      arg=x(1)**2+x(2)**2
      dmaxex = 350.0d0
      if (arg.gt.dmaxex) then
        ifail = 1
        return
      endif
c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1)=dexp(arg)-p(1)
        f(2)=x(1)+x(2)-dsin(p(2)*(x(1)+x(2)))
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		dfdx(1,1)=2.d0*x(1)*dexp(arg)
		dfdx(1,2)=2.d0*x(2)*dexp(arg)
		dfdx(2,1)=1.d0-p(2)*dcos(p(2)*(x(1)+x(2)))
		dfdx(2,2)=dfdx(2,1)
      endif
      return
      end
