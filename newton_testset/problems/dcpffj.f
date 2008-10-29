      subroutine dcpffj(n, x, fvec, fjac, ldfjac, task, par, np, ldx, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      implicit none
      character*(*) task
      integer n, ldfjac, np, nx0, ldx, nbound, ncon, nsol, ifail
      double precision x(*), fvec(*), fjac(ldfjac,*), par(*), x0(ldx,*),
     $                 xlow(ldx,*), xup(ldx,*), conxl(ldx,*),
     $                 conxu(ldx,*), xsol(ldx,*), xthrsh(*), fthrsh(*),
     $                 rtol
c-----------------------------------------------------------------------
c
c     Subroutine dcpffj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the combustion of propane (full formulation) problem.
c
c     Adapted by       L. Weimann, 
c                      Zuse Institute Berlin (ZIB),
c                      Dep. Numerical Analysis and Modelling
c     File             dcpffj.f
c     Version          1.0
c     Latest Change    13th December 2005
c     Code             Fortran 77 with common extensions
c                      Double precision
c
c     Reference:
c       M.B. Averick, R.G. Carter, J.J. More, G.-L. Xue:   
c       The MINPACK-2 test problem collection.
c       Argonne National Laboratory, Preprint MCS-P153-0692
c
c-----------------------------------------------------------------------
      double precision one, p, p5, rr, ten, two, zero
      parameter (zero=0.0d0,p5=0.5d0,one=1.0d0,two=2.0d0,ten=1.0d1)
      parameter (p=4.0d1,rr=ten)

      integer i, j
      double precision pdx, sqpdx, sqrtp, xfrac, xtau
      double precision k(10)

      data k/zero, zero, zero, zero, 1.930d-1, 2.597d-3, 3.448d-3,
     +     1.799d-5, 2.155d-4, 3.846d-5/

      ifail = 0

      sqrtp = sqrt(p)

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         n = 11
         x0(1,1) = 5.0d0
         x0(2,1) = 2.5d0
         x0(3,1) = 5.0d0
         x0(4,1) = 1.0d-1
         x0(5,1) = 5.0d-2*k(5)
         x0(6,1) = k(6)/sqrtp
         x0(7,1) = 5.0d1*k(7)/sqrtp
         x0(8,1) = 1.0d3*k(8)/p
         x0(9,1) = 5.0d2*k(9)/sqrtp
         x0(10,1) = 5.0d4*k(10)/p
         x0(11,1) = 2.0d1
c
         x0( 1,2) = 1.0614059029506D-22
         x0( 2,2) = 9.7801821638273D-18    
         x0( 3,2) = 5.7688208649105D-20
         x0( 4,2) = 3.8913971394121D-02    
         x0( 5,2) = 8.0180330879019D-20    
         x0( 6,2) = 4.6586053320405D-14
         x0( 7,2) = 2.3836251700297D-15    
         x0( 8,2) = 1.7315685102394D-13    
         x0( 9,2) = 2.8555467179170D-22
         x0(10,2) = 5.8824393637661D-04    
         x0(11,2) = 8.9581084045863D-28
c
         nx0 = 2
c
         do 10 i = 1, n
            xlow(i,1) = zero
            xup (i,1) = 1.0d3
c           xthrsh(i) = 1.0d-19
   10    continue
         nbound = 1
         return

      end if

c     Check input arguments for errors.

      if (n .ne. 11) then
         task = 'ERROR: N .NE. 11 IN DCPFFJ'

         return

      end if
      do i=1, 4
        if (x(i).le.0.0d0) then
          ifail = 1
          return
        endif
      enddo
	  if (x(11).le.0.0d0) then
		ifail = 1
		return
	  endif
      

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      pdx = p/x(11)
      sqpdx = sqrt(pdx)

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         xtau = zero
         do 20 i = 1, n - 1
            xtau = xtau + x(i)
   20    continue
         fvec(1) = x(1) + x(4) - 3.0d0
         fvec(2) = two*x(1) + x(2) + x(4) + x(7) + x(8) + x(9) +
     +             two*x(10) - rr
         fvec(3) = two*x(2) + two*x(5) + x(6) + x(7) - 8.0d0
         fvec(4) = two*x(3) + x(9) - 4.0d0*rr
         fvec(5) = k(5)*x(2)*x(4) - x(1)*x(5)
         fvec(6) = k(6)*sqrt(x(2)*x(4)) - sqrt(x(1))*x(6)*sqpdx
         fvec(7) = k(7)*sqrt(x(1)*x(2)) - sqrt(x(4))*x(7)*sqpdx
         fvec(8) = k(8)*x(1) - x(4)*x(8)*pdx
         fvec(9) = k(9)*x(1)*sqrt(x(3)) - x(4)*x(9)*sqpdx
         fvec(10) = k(10)*x(1)**2 - (x(4)**2)*x(10)*pdx
         fvec(11) = x(11) - xtau
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         do 40 j = 1, n
            do 30 i = 1, n - 1
               fjac(i,j) = zero
   30       continue
            fjac(n,j) = -one
   40    continue
         fjac(n,n) = one

         xfrac = one/(sqrt(x(11))**3)

         fjac(1,1) = one
         fjac(2,1) = two
         fjac(5,1) = -x(5)
         fjac(6,1) = -p5*x(6)*sqpdx/sqrt(x(1))
         fjac(7,1) = p5*k(7)*sqrt(x(2))/sqrt(x(1))
         fjac(8,1) = k(8)
         fjac(9,1) = k(9)*sqrt(x(3))
         fjac(10,1) = two*k(10)*x(1)

         fjac(2,2) = one
         fjac(3,2) = two
         fjac(5,2) = k(5)*x(4)
         fjac(6,2) = p5*k(6)*sqrt(x(4))/sqrt(x(2))
         fjac(7,2) = p5*k(7)*sqrt(x(1))/sqrt(x(2))

         fjac(4,3) = two
         fjac(9,3) = p5*k(9)*x(1)/sqrt(x(3))

         fjac(1,4) = one
         fjac(2,4) = one
         fjac(5,4) = k(5)*x(2)
         fjac(6,4) = p5*k(6)*sqrt(x(2))/sqrt(x(4))
         fjac(7,4) = -p5*x(7)*sqpdx/sqrt(x(4))
         fjac(8,4) = -x(8)*pdx
         fjac(9,4) = -x(9)*sqpdx
         fjac(10,4) = -two*x(4)*x(10)*pdx

         fjac(3,5) = two
         fjac(5,5) = -x(1)

         fjac(3,6) = one
         fjac(6,6) = -sqrt(x(1))*sqpdx

         fjac(2,7) = one
         fjac(3,7) = one
         fjac(7,7) = -sqrt(x(4))*sqpdx

         fjac(2,8) = one
         fjac(8,8) = -x(4)*pdx

         fjac(2,9) = one
         fjac(4,9) = one
         fjac(9,9) = -x(4)*sqpdx

         fjac(2,10) = two
         fjac(10,10) = -(x(4)**2)*pdx

         fjac(6,11) = p5*sqrt(x(1))*x(6)*sqrtp*xfrac
         fjac(7,11) = p5*sqrt(x(4))*x(7)*sqrtp*xfrac
         fjac(8,11) = x(4)*x(8)*p/(x(11)**2)
         fjac(9,11) = p5*x(4)*x(9)*sqrtp*xfrac
         fjac(10,11) = x(4)**2*x(10)*p/(x(11)**2)
      end if

      end
