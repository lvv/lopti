      subroutine Cutlip(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Cutlips steady state for reaction rate equations,
c   parameter set 1
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             Cutlip.f
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
      double precision p_k1, p_k2, p_k3, p_kr1, p_kr2
      double precision x00(7,6)
      integer npar, i, j
c
      data (x00(1,i),i=1,6) /0.99, 0.05, 0.05, 0.99, 0.05, 0.0/
      data (x00(2,i),i=1,6) /0.05, 0.99, 0.05, 0.05, 0.99, 0.0/
      data (x00(3,i),i=1,6) /0.97, 0.98, 0.06, 0.99, 0.0, 0.0/ 
      data (x00(4,i),i=1,6) /3.56e-2, 3.57e-1, 1.92, 3.60e-2, 8.84e-2,
     $                       8.76e-1/
      data (x00(5,i),i=1,6) /3.62e-2, 3.57e-1, 1.92, 9.36e-2, 3.40e-2,
     $                       8.72e-1/
      data (x00(6,i),i=1,6) /1.03, 1.02, -6.65e-2, 1.10e-4, 1.00,
     $                       -1.03e-3/
      data (x00(7,i),i=1,6) / 150.994, 1066.746, 1466.816, 0.438855,
     $                        0.54453, 0.016612/
c
      npar = np
      if (np.ge.0 .and. np.le.2) then
        if (np.eq.0) then
		  p_k1 = 31.24
		  p_k2 = 2.062
		  p_k3 = 303.03
		  p_kr1 = 0.272
		  p_kr2 =  0.02
		else if (np.eq.1) then
          p_k1 = 17.721
          p_k2  = 3.483
          p_kr1= 0.118
          p_kr2 =  0.033
          p_k3 = 505.051
        else
          p_k1 = 17.721
          p_k2  = 6.966
          p_kr1= 0.118
          p_kr2 =  333.333
          p_k3 = 505.051
		endif
		p(1) = p_k1
		p(2) = p_k2
		p(3) = p_k3
		p(4) = p_kr1
		p(5) = p_kr2
		np = 5
      else
		p_k1 = p(1)
		p_k2 = p(2)
		p_k3 = p(3)
		p_kr1 = p(4)
		p_kr2 = p(5)
      endif  
c
      if (task .eq. 'XS') then
        n = 6
        if (npar.eq.0 .or. npar.eq.1) then
		  nx0 = 6
		  do j=1,nx0
			do i=1,n
			  x0(i,j) = x00(j,i)
			enddo
		  enddo
        else
		  nx0 = 3
		  do j=1,nx0-1
			do i=1,n
			  x0(i,j) = x00(j,i)
			enddo
		  enddo
	      do i=1,n
			x0(i,3) = x00(7,i)
		  enddo
        endif
c
        do i=1,n
          xthrsh(i) = 1.0d-4
          fthrsh(i) = 1.0d0
        enddo
        rtol = 1.0d-10
c
c range of interest/solution
        do i=1,6
          xlow(i,1) = 0.0d0
          xup (i,1) = 10.0d0
        enddo
        nbound = 1
        nsol = 0
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1) = 1 - x(1) - p_k1*x(1)*x(6) + p_kr1*x(4)
        f(2) = 1 - x(2) - p_k2*x(2)*x(6) + p_kr2*x(5)
        f(3) = - x(3) + 2*p_k3*x(4)*x(5)
        f(4) = p_k1*x(1)*x(6) - p_kr1*x(4) - p_k3*x(4)*x(5)
        f(5) = 1.5*(p_k2*x(2)*x(6) - p_kr2*x(5)) - p_k3*x(4)*x(5)
        f(6) = 1 - x(4) - x(5) - x(6)
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        do j=1,6
          do i=1,6
            dfdx(i,j) = 0.0d0
          enddo
        enddo
        dfdx( 1, 1) = -1.0 - p_k1*x(6)
        dfdx( 1, 4) = p_kr1
        dfdx( 1, 6) = -p_k1*x(1)
        dfdx( 2, 2) = -1.0 - p_k2*x(6)
        dfdx( 2, 5) = p_kr2
        dfdx( 2, 6) = -p_k2*x(2)
        dfdx( 3, 3) = -1.0
        dfdx( 3, 4) = 2*p_k3*x(5)
        dfdx( 3, 5) = 2*p_k3*x(4)
        dfdx( 4, 1) = p_k1*x(6)
        dfdx( 4, 4) = -p_kr1 - p_k3*x(5)
        dfdx( 4, 5) = -p_k3*x(4)
        dfdx( 4, 6) = p_k1*x(1)
        dfdx( 5, 2) = 1.5*p_k2*x(6)
        dfdx( 5, 4) = -p_k3*x(5)
        dfdx( 5, 5) = -1.5*p_kr2 - p_k3*x(4)
        dfdx( 5, 6) = 1.5*p_k2*x(2)
        dfdx( 6, 4) = -1.0
        dfdx( 6, 5) = -1.0
        dfdx( 6, 6) = -1.0
      endif
      ifail = 0
      return
      end
      