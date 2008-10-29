      subroutine ChemEqui2(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Chemical equilibrium 2
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             ChemEqui2.f
c   Version          1.0
c   Latest Change    7th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     K.L. Hiebert:
c     An evaluation of mathematical software that
c     solves systems of nonlinear equations.
c     ACM Trans. Math. Soft. 8, 5-20, 1982
c
c-----------------------------------------------------------------------
      integer i
c
      ifail = 0
      if (task .eq. 'XS') then
		 n = 6
		 x0(1,1)=3.0d-4
		 x0(2,1)=3.0d-4
		 x0(3,1)=27.5d0
		 x0(4,1)=3.0d-4
		 x0(5,1)=27.5d0
		 x0(6,1)=27.5d0
		 do i=1,n
			x0(i,2)=1.0d0
		 enddo
         nx0 = 2
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
		f(1)=x(1)+x(2)+x(4)-1.0d-3 
		f(2)=x(5)+x(6)-55.0d0 
		f(3)=x(1)+x(2)+x(3)+2.0d0*x(5)+x(6)-110.001d0 
		f(4)=x(1)-1.0d-1*x(2) 
		f(5)=x(1)-1.0d4*x(3)*x(4) 
		f(6)=x(5)-55.0d14*x(3)*x(6) 
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
		dfdx(1,1)=1.0d0
		dfdx(1,2)=1.0d0
		dfdx(1,3)=0.0d0
		dfdx(1,4)=1.0d0
		dfdx(1,5)=0.0d0
		dfdx(1,6)=0.0d0
c
		dfdx(2,1)=0.0d0
		dfdx(2,2)=0.0d0
		dfdx(2,3)=0.0d0
		dfdx(2,4)=0.0d0
		dfdx(2,5)=1.0d0
		dfdx(2,6)=1.0d0
c
		dfdx(3,1)=1.0d0
		dfdx(3,2)=1.0d0
		dfdx(3,3)=1.0d0
		dfdx(3,4)=0.0d0
		dfdx(3,5)=2.0d0
		dfdx(3,6)=1.0d0 
c
		dfdx(4,1)=1.0d0
		dfdx(4,2)=-0.1d0
		dfdx(4,3)=0.0d0
		dfdx(4,4)=0.0d0
		dfdx(4,5)=0.0d0
		dfdx(4,6)=0.0d0
c
		dfdx(5,1)=1.0d0
		dfdx(5,2)=0.0d0
		dfdx(5,3)=-1.0d4*x(4)
		dfdx(5,4)=-1.0d4*x(3)
		dfdx(5,5)=0.0d0
		dfdx(5,6)=0.0d0
c
		dfdx(6,1)=0.0d0
		dfdx(6,2)=0.0d0
		dfdx(6,3)=-55.0d14*x(6)
		dfdx(6,4)=0.0d0
		dfdx(6,5)=1.0d0
		dfdx(6,6)=-55.0d14*x(3)
      endif
      return
      end
