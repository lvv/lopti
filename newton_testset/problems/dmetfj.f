      subroutine dmetfj(n, x, fvec, fjac, ldfjac, task, p, np, ldx, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      implicit none
      character*(*) task
      integer n, ldfjac, np, nx0, ldx, nbound, ncon, nsol, ifail
      double precision x(*), fvec(*), fjac(ldfjac,*), p(*), x0(ldx,*),
     $                 xlow(ldx,*), xup(ldx,*), conxl(ldx,*),
     $                 conxu(ldx,*), xsol(ldx,*), xthrsh(*), fthrsh(*),
     $                 rtol
c-----------------------------------------------------------------------
c
c   Subroutine dmetfj
c
c   This subroutine computes the function and Jacobian matrix of
c   the - chemical equilibrium resulting from partial oxidation
c   of methane with oxygen - problem. (full formulation)
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             dmetfj.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     R. Luus:
c     Iterative Dynamic Programming.
c     Chapman & Hall/CRC, 2000
c     page 23ff (example 3)
c
c-----------------------------------------------------------------------
      integer i,j

      ifail = 0
c     Compute the standard starting point if task = 'XS'

      if (task .eq. 'XS') then
        x(1) = 0.22d0
        x(2) = 0.075d0
        x(3) = 0.001d0
        x(4) = x(1)*x(3)/(2.6058d0*x(2))
        x(5) = 400.0d0*x(1)*x(4)**3/(1.7837d5*x(3))
        x(7) = 2.0d0/(x(3)+x(4)+2.0d0*x(5))
        x(6) = x(7)*(0.5d0*(x(1)+x(3))+x(2))
        n = 7
        nx0 = 3
        do i=1,n
          x0(i,1) = x(i)
        enddo
        do i=1,n
          x0(i,2) = 0.5d0
        enddo
        x0(1,3) = 0.322
        x0(2,3) = 0.0092
        x0(3,3) = 0.046
        x0(4,3) = 0.6181
        x0(5,3) = 0.0037
        x0(6,3) = 0.5767
        x0(7,3) = 2.978
        do i=1,n
           xlow(i,1) =  0.0d0
           xup (i,1) =  1.0d1
        enddo
        nbound = 1
        return
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      do i=2,3
        if (x(i).eq.0.0d0) then
          ifail = 1
          return
        endif
      enddo
      if (x(7).eq.0.0d0) then
        ifail = 1
        return
      endif
      
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        fvec(1) = x(1)*x(3)/(2.6058d0*x(2))-x(4)
        fvec(2) = 400.0d0*x(1)*x(4)**3/(1.7837d5*x(3))-x(5)
        fvec(3) = 2.0d0/(x(3)+x(4)+2.0d0*x(5))-x(7)
        fvec(4) = x(7)*(0.5d0*(x(1)+x(3))+x(2))-x(6)
        fvec(5) = x(1)+x(2)+x(5)-1.0d0/x(7)
        fvec(6) = -28837.0d0*x(1)-139009.0d0*x(2)-78213.0d0*x(3)
     $            +18927.0d0*x(4)+8427.0d0*x(5) 
     $            +(13492.0d0-10690.0d0*x(6))/x(7)
        fvec(7) = x(1)+x(2)+x(3)+x(4)+x(5)-1.0d0
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        do 20 i=1,7
          do 20 j=1,7
            fjac(i,j) = 0.0d0
20      continue
        fjac(1,1) = x(3)/(2.6058d0*x(2))
        fjac(1,3) = x(1)/(2.6058d0*x(2))
        fjac(1,2) = x(1)*x(3)/(2.6058d0*x(2))**2*2.6058d0
        fjac(1,4) = -1.0d0
        fjac(2,1) = 400.0d0*x(4)**3/(1.7837d5*x(3))
        fjac(2,4) = 400.0d0*x(1)*3.0d0*x(4)**2/(1.7837d5*x(3))
        fjac(2,3) = -400.0d0*x(1)*x(4)**3/(1.7837d5*x(3))**2*1.7837d5
        fjac(2,5) = -1.0d0
        fjac(3,3) = -2.0d0/(x(3)+x(4)+2.0d0*x(5))**2
        fjac(3,4) = -2.0d0/(x(3)+x(4)+2.0d0*x(5))**2
        fjac(3,5) = -4.0d0/(x(3)+x(4)+2.0d0*x(5))**2
        fjac(3,7) = -1.0d0
        fjac(4,7) = 0.5d0*(x(1)+x(3))+x(2)
        fjac(4,1) = 0.5d0*x(7)
        fjac(4,3) = 0.5d0*x(7)
        fjac(4,2) = x(7)
        fjac(4,6) = -1.0d0
        fjac(5,1) = 1.0d0
        fjac(5,2) = 1.0d0
        fjac(5,5) = 1.0d0
        fjac(5,7) = 1.0d0/x(7)**2
        fjac(6,1) = -28837.0d0
        fjac(6,2) = -139009.0d0
        fjac(6,3) = -78213.0d0
        fjac(6,4) = 18927.0d0
        fjac(6,5) = 8427.0d0
        fjac(6,6) = -10690.0d0/x(7)
        fjac(6,7) = -(13492.0d0-10690.0d0*x(6))/x(7)**2
        do 30 i=1,5
          fjac(7,i) = 1.0d0
30      continue
      end if

      end
