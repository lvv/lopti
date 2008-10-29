      subroutine PartOxiOfMet(n, xa, fvec, fjac, ldfjac, task, p,np,ldx, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      implicit none
      character*(*) task
      integer n, ldfjac, np, nx0, ldx, nbound, ncon, nsol, ifail
      double precision xa(*), fvec(*), fjac(ldfjac,*), p(*), x0(ldx,*),
     $                 xlow(ldx,*), xup(ldx,*), conxl(ldx,*),
     $                 conxu(ldx,*), xsol(ldx,*), xthrsh(*), fthrsh(*),
     $                 rtol
c-----------------------------------------------------------------------
c
c   Subroutine PartOxiOfMet
c
c   This subroutine computes the function and Jacobian matrix of
c   the - chemical equilibrium resulting from partial oxidation
c   of methane with oxygen - problem. (reduced formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             PartOxiOfMet.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     R. Luus:
c     Iterative Dynamic Programming.
c     Chapman & Hall/CRC, 2000
c     page 23ff. (example 3)
c
c-----------------------------------------------------------------------
      integer i,j
      double precision x(7)

      ifail = 0
c     Compute the standard starting point if task = 'XS'

      if (task .eq. 'XS') then
        x(1) = 0.22d0
        x(2) = 0.075d0
        x(3) = 0.001d0
        n = 3
        nx0 = 1
        do i=1,n
          x0(i,1) = x(i)
        enddo
        do i=1,n
           xlow(i,1) =  0.0d0
           xup (i,1) =  0.95d0
        enddo
        nbound = 1
        return
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      do i=2,3
        if (xa(i).eq.0.0d0) then
          ifail = 1
          return
        endif
      enddo
      
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        x(1) = xa(1)
        x(2) = xa(2)
        x(3) = xa(3)
        x(4) = x(1)*x(3)/(2.6058d0*x(2))
        x(5) = 400.0d0*x(1)*x(4)**3/(1.7837d5*x(3))
        x(7) = 2.0d0/(x(3)+x(4)+2.0d0*x(5))
        x(6) = x(7)*(0.5d0*(x(1)+x(3))+x(2))
        if (x(7).eq.0.0d0) then
          ifail = 1
          return
        endif
        fvec(1) = x(1)+x(2)+x(5)-1.0d0/x(7)
        fvec(2) = -28837.0d0*x(1)-139009.0d0*x(2)-78213.0d0*x(3)
     $            +18927.0d0*x(4)+8427.0d0*x(5) 
     $            +(13492.0d0-10690.0d0*x(6))/x(7)
        fvec(3) = x(1)+x(2)+x(3)+x(4)+x(5)-1.0d0
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ifail = -999
      end if

      end
