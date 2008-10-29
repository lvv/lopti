      subroutine PartMethaneOx(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Chemical equilibrium resulting from a Partial Methane Oxidation
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             PartMethaneOx.f
c   Version          1.0
c   Latest Change    9th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     M. Shacham:
c     Comparing Software for the Solution of Systems
c     of Nonlinear Algebraic Equations Arising in
c     Chemical Engineering.
c     Computers & Chemical Engineering 9, 103-112, 1985 
c
c-----------------------------------------------------------------------
      integer i,j
      double precision x00(7,5)
c     the first hard initial guess
      data (x00(i,1),i=1,7) / 0.5,0.0,0.0,0.5,0.0,0.5,2.0 /
c     the second hard initial guess
      data (x00(i,2),i=1,7) / 0.22,0.075,0.001,0.58,0.125,0.436,2.35 /
c     an initial guess near the feasible solution
      data (x00(i,3),i=1,7) / 0.3,0.01,0.05,0.6,0.004,0.6,3.0/
c     an initial guess near Luus infeasible solution
      data (x00(i,4),i=1,7) / 1.5,-1.13,1.33,-0.66,-0.0007,0.8,3.0/
c     an initial guess hard for global residual oriented (2xx iterates)
      data (x00(i,5),i=1,7) / 1.5,0.001,1.33,1.e-3,1.e-4,0.8,3.0/
c
      ifail = 0
      if (task .eq. 'XS') then
        n = 7
        nx0 = 5
        do j=1,nx0
          do i=1,n
            x0(i,j) = x00(i,j)
          enddo
        enddo
        do i=1,n
          xthrsh(i) = 1.0d0
          fthrsh(i) = 1.0d0
        enddo
        rtol = 1.0d-9
c       range of interest/solution
        do i=1,n
          xlow(i,1) = 1.0d-4
        enddo
        do i=1,5
          xup(i,1) = 1.0d0
        enddo
        xup(6,1) = 5.0d0
        xup(7,1) = 5.0d0
        do i=1,n
          xup (i,2) = 5.0d0
          xlow(i,2) = -xup(i,2)
        enddo
        nbound = 2
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        f(1) = 0.5*x(1) + x(2) + 0.5*x(3) - x(6)/x(7)
        f(2) = x(3) + x(4) + 2*x(5) - 2/x(7) 
        f(3) = x(1) + x(2) + x(5) - 1/x(7) 
        f(4) = -28837*x(1) - 139009*x(2) - 78213*x(3) + 18927*x(4)
     $            + 8427*x(5) + 13492/x(7) - 10690*x(6)/x(7)
        f(5) = x(1) + x(2) + x(3) + x(4) + x(5) - 1
        f(6) = 400*x(1)*x(4)**3 - 1.7837e5*x(3)*x(5)
        f(7) = x(1)*x(3) - 2.6058*x(2)*x(4)
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        dfdx(1,1) = 0.5
        dfdx(1,2) = 1.0
        dfdx(1,3) = 0.5
        dfdx(1,6) = -1/x(7)
        dfdx(1,7) = x(6)/x(7)**2
        dfdx(2,3) = 1.0
        dfdx(2,4) = 1.0
        dfdx(2,5) = 2.0
        dfdx(2,7) = 2.0/x(7)**2
        dfdx(3,1) = 1.0
        dfdx(3,2) = 1.0
        dfdx(3,5) = 1.0
        dfdx(3,7) = 1.0/x(7)**2
        dfdx(4,1) = -28837
        dfdx(4,2) = -139009
        dfdx(4,3) = -78213
        dfdx(4,4) =  18927
        dfdx(4,5) =  8427
        dfdx(4,6) = -10690/x(7)
        dfdx(4,7) = -(13492 - 10690*x(6))/x(7)**2
        dfdx(5,1) = 1.0
        dfdx(5,2) = 1.0
        dfdx(5,3) = 1.0
        dfdx(5,4) = 1.0
        dfdx(5,5) = 1.0
        dfdx(6,1) = 400*x(4)**3
        dfdx(6,3) = -1.7837e5*x(5)
        dfdx(6,4) = 1200*x(1)*x(4)**2
        dfdx(6,5) = -1.7837e5*x(3)
        dfdx(7,1) = x(3)
        dfdx(7,2) = -2.6058*x(4)
        dfdx(7,3) = x(1)
        dfdx(7,4) = -2.6058*x(2)
      endif
      return
      end
