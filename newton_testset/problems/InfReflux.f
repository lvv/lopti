      subroutine InfReflux(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   InfReflux function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             InfReflux.f
c   Version          1.0
c   Latest Change    9th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     Paterson, W. R.:
c     Chem. Engng. Sci., 41, 1986, p.1935 ;
c     http://www.polymath-software.com/library/nle/Oneeq10.htm
c
c-----------------------------------------------------------------------
c
      double precision arg1,arg2
c
      ifail = 0
      if (task .eq. 'XS') then
        n = 1
c       default initial guess (critical as close to J=0 for x=0.229)
        x0(1,1) = 0.23
c       three alternative initial guesses
c       an even harder initial point, and two relatively easy points
        x0(1,2) = 0.228
        x0(1,3) = 0.6
        x0(1,4) = 0.01
        nx0 = 4
        xthrsh(1) = 1.0d0
        fthrsh(1) = 1.0d0
        rtol = 1.0d-10
c       range of interest/solution
        xlow(1,1) = 0.01
        xup (1,1) = 0.7
        nbound = 1
c       constrained range 
        conxl(1,1) = 0.0
        conxu(1,1) = 1.0
        ncon = 1
        return
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      arg1 = 1./(1.-x(1))
      arg2 = 0.95-x(1)
      if (x(1).le.0.0d0 .or. arg1.le.0.0d0 .or. arg2 .le. 0.0d0) then
        ifail = 1
        return
      endif
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1) = (1./63.)*dlog(x(1))+(64./63.)*dlog(arg1)+
     $         dlog(arg2)-dlog(0.9d0)
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        dfdx(1,1) = (1./63.)/x(1)+(64./63.)/(1.-x(1))-1./(0.95-x(1))
      endif
      ifail = 0
      return
      end

