      subroutine sst(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   SST nonlinearity term.
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             sst.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     R.F. Sincovec, N. K. Madsen:
c     Software for Nonlinear Partial Differential Equations.
c     ACM Transactions on Mathematical Software Vol. 1, No. 3, 
c     p.232-260, 1975

c
c-----------------------------------------------------------------------
      integer iprob
      double precision sst1
c
      ifail = 0
      if (np.eq.0) then
        iprob = 1
        p(1) = dble(iprob)
        np = 1
      else
        iprob = int(p(1))
      endif
      if (task .eq. 'XS') then
        n = 4
		x0(1,1)=1.0d9
		x0(2,1)=1.0d9
		x0(3,1)=1.0d13
		x0(4,1)=1.0d7
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		if (iprob.eq.1) then
		  sst1=360.0d0
		else if (iprob.eq.2) then
		  sst1=3250.0d0 
		else
		  ifail = -1000
		  return
		endif
        f(1) = 4.0d5 - 272.443800016*x(1) + 1.0d-4*x(2)
     $        + 0.007d0*x(4) - 3.67d-16*x(1)*x(2) - 4.13d-12*x(1)*x(4)
        f(2) = 272.4438d0*x(1) - 1.00016d-4*x(2)
     $        + 3.67d-16*x(1)*x(2) - 3.57d-15*x(2)*x(3)
        f(3) = -1.6d-8*x(3) + 0.007d0*x(4)
     $        + 4.1283d-12*x(1)*x(4) - 3.57d-15*x(2)*x(3) + 800.0d0+sst1
        f(4) = -7.000016d-3*x(4) + 3.57d-15*x(2)*x(3)
     $        - 4.1283d-12*x(1)*x(4) + 800.0d0
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		dfdx(1,1) =  - 272.443800016 - 3.67d-16*x(2) - 4.13d-12*x(4)
		dfdx(2,1) = 272.4438d0  + 3.67d-16*x(2) 
		dfdx(3,1) =  4.1283d-12*x(4) 
		dfdx(4,1) = - 4.1283d-12*x(4) 
		dfdx(1,2) =  1.0d-4 - 3.67d-16*x(1)
		dfdx(2,2) = - 1.00016d-4 + 3.67d-16*x(1) - 3.57d-15*x(3)
		dfdx(3,2) = - 3.57d-15*x(3) 
		dfdx(4,2) =  + 3.57d-15*x(3)
		dfdx(1,3) = 0.0d0
		dfdx(2,3) = - 3.57d-15*x(2)
		dfdx(3,3) = -1.6d-8 - 3.57d-15*x(2) 
		dfdx(4,3) = 3.57d-15*x(2)
		dfdx(1,4) =  0.007d0  - 4.13d-12*x(1)
		dfdx(2,4) = 0.0d0
		dfdx(3,4) =  0.007d0 + 4.1283d-12*x(1)
		dfdx(4,4) = -7.000016d-3 - 4.1283d-12*x(1)
      endif
      return
      end
