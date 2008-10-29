      subroutine Order10to11(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Order10to11 function
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Order10to11.f
c   Version          1.0
c   Latest Change    9th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c   Reference:
c     Shacham, M. and Kehat, E.:
c     Chem. Engng. Sci, 27, 1972, p.2099
c     http://www.polymath-software.com/library/nle/Oneeq3.htm
c
c-----------------------------------------------------------------------
      double precision T
c
      if (task .eq. 'XS') then
        n = 1
c default initial guess 
        x0(1,1) = 555.0
        nx0 = 1
        xthrsh(1) = 1.0d0
        fthrsh(1) = 1.0d10
        rtol = 1.0d-10

c range of interest/solution
        xlow(1,1) = 550.0
        xup (1,1) = 560.0
        nbound = 1
c constrained range 
        conxl(1,1) = 1.0d-10
        conxu(1,1) = 1.0d+100
        ncon = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        T = x(1)
        f(1) = dexp(21000./T)/T**2 - 1.11d11
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        T = x(1)
        dfdx(1,1) = (-dexp(21000./T)*21000. - exp(21000./T)*2.*T)/T**4
      endif
      ifail = 0
      return
      end
