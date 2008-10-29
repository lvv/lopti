      subroutine LargeExp(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Semiconductor
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             LargeExp.f
c   Version          1.0
c   Latest Change    9th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     J. Molenaar, P.W. Hemker:
c     A multigrid approach for the solution of the 
c     2D semiconductor equations.
c     IMPACT Comput. Sci. Eng. 2, No. 3, 219-243, 1990.
c
c-----------------------------------------------------------------------
      integer i, j
      double precision epsilo, q, x00, dot, dotr, dotl, xnh, xph,
     $                 h1, h2, h3, h4, axnh1, axph1, axnh2, axph2,
     $                 dmaxex
      double precision p1,p2l,p2r,alfa,xni,volt
      save p1,p2l,p2r,alfa,xni,volt
c
      ifail = 0
      if (np.eq.0) then
        volt = 100.0d0
        p(1) = volt
        np = 1
      else
        volt = p(1)
      endif
      if (task .eq. 'XS') then
         n = 6
         alfa=38.683d0
         epsilo=1.036d-12
         xni=1.22d10
         q=1.6d-19
         dot=1.d17
         dotr=dot
         dotl=-dot
         p1=xni*q/epsilo
         p2r=dotr/xni
         p2l=dotl/xni
         x00=1.d0
         do i=1,n
            x0(i,1)=x00 
         enddo 
         nx0 = 1
         return
      endif
      
	  h1 = alfa*(x(1)-x(2))
	  h2 = alfa*(x(3)-x(1))
	  h3 = alfa*(x(4)-x(5))
	  h4 = alfa*(x(6)-x(4))
	  dmaxex = 350.0d0
	  if (dmax1(h1,h2,h3,h4).gt.dmaxex) then
		ifail = 1
		return
	  endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
c
c   left boundary
c
        f(2)=x(2)-0.d0
        f(3)=x(3)-0.d0
        xnh=dexp(h1)
        xph=dexp(h2)
        f(1)=xph-xnh+p2l
c
c   right boundary
c
        f(5)=x(5)-volt
        f(6)=x(6)-volt
        xnh=dexp(h3)
        xph=dexp(h4)
        f(4)=xph-xnh+p2r
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        axnh1=alfa*dexp(h1)
        axph1=alfa*dexp(h2)
        axnh2=alfa*dexp(h3)
        axph2=alfa*dexp(h4)
        do i=1,n
          do j=1,n
            dfdx(j,i) = 0.0d0
          enddo
        enddo
        dfdx(1,1) = -axph1-axnh1
        dfdx(1,2) = axnh1
        dfdx(2,2) = 1.0d0
        dfdx(1,3) = axph1
        dfdx(3,3) = 1.0d0
        dfdx(4,4) = -axph2-axnh2
        dfdx(4,5) = axnh2
        dfdx(5,5) = 1.0d0
        dfdx(4,6) = axph2
        dfdx(6,6) = 1.0d0
      endif
      return
      end
