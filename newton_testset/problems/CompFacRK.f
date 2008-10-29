      subroutine CompFacRK(n, x, f, dfdx, ldfjac, task, par, np, ldx, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      implicit none
      character*(*) task
      integer n, ldfjac, np, nx0, ldx, nbound, ncon, nsol, ifail
      double precision x(*), f(*), dfdx(ldfjac,*), par(*), x0(ldx,*),
     $                 xlow(ldx,*), xup(ldx,*), conxl(ldx,*),
     $                 conxu(ldx,*), xsol(ldx,*), xthrsh(*), fthrsh(*),
     $                 rtol
c-----------------------------------------------------------------------
c
c  Compressibility factor from the RK equation - three solutions
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             CompFacRK.f
c   Version          1.0
c   Latest Change    8th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     Cutlip, M. B. and Shacham, M:
c     Problem Solving in Chemical Engineering with Numerical
c     Methods (2nd Ed).
c     Prentice Hall  Inc., 1999
c
c-----------------------------------------------------------------------
      double precision P, Pc, T, Tc, Pr, Tr, Asqr, B, Q, r, z
c
      ifail = 0
      if (task .eq. 'XS') then
        n = 1
c       default initial guess (critical as not far from J=0, should find x*(2))
        x0(1,1) = 0.65
c       three alternative initial guesses
c       from first one should run to x*(1), from second easily to x*(3)     
c       from last to x*(2) ??
        x0(1,2) = -0.5
        x0(1,3) = 1.0
        x0(1,4) = -0.02
        nx0 = 4
        xthrsh(1) = 1.0d0
        fthrsh(1) = 1.0d0
        rtol = 1.0d-10
c       range of interest/solution
        xlow(1,1) = -0.50
        xup (1,1) =  1.20
        nbound = 1
        return
      endif
c
	  P = 200 
	  Pc=33.5 
	  T= 631*2 
	  Tc=126.2
	  Pr=P/Pc
	  Tr=T/Tc  
	  Asqr=0.4278*Pr/(Tr**2.5)  
	  B=0.0867*Pr/Tr   
	  Q=B**2+B-Asqr  
	  z = x(1) 

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        r=Asqr*B   
        f(1) = z**3 - z**2 - Q*z - r
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        dfdx(1,1) = 3*z**2 - 2*z - Q
      endif
      ifail = 0
      return
      end
