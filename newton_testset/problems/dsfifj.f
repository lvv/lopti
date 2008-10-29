      subroutine dsfifj(n, x, fvec, fjac, ldfjac, task, p, np, ldx, 
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
c     Subroutine dsfifj
c
c     This subroutine computes the function and Jacobian matrix of
c     the solid fuel ignition problem.
c
c     Adapted by       L. Weimann, 
c                      Zuse Institute Berlin (ZIB),
c                      Dep. Numerical Analysis and Modelling
c     File             dsfifj.f
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
      integer nx, ny
      double precision lambda
      double precision four, one, two, zero
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)

      integer i, j, k
      double precision hx, hxdhy, hxhy, hy, hydhx, temp, temp1, u, ub,
     +                 ul, ur, ut, uxx, uyy
      ifail = 0
      if (np.eq.0) then
        nx = 10
        ny = 10
        lambda = 5.0d0
        p(1) = dble(nx)
        p(2) = dble(ny)
        p(3) = lambda
        np = 3
      else
        nx = int(p(1))
        ny = int(p(2))
        lambda = p(3)
      endif
      n = nx*ny
      hx = one/dble(nx+1)
      hy = one/dble(ny+1)
      hydhx = hy/hx
      hxdhy = hx/hy
      hxhy = hx*hy

c     Compute the standard starting point if task = 'XS'

      if (task .eq. 'XS') then
         temp1 = lambda/(lambda+one)
         do 20 j = 1, ny
            temp = dble(min(j,ny-j+1))*hy
            do 10 i = 1, nx
               k = nx*(j-1) + i
               x(k) = temp1*sqrt(min(dble(min(i,nx-i+1))*hx,temp))
   10       continue
   20    continue
         nx0 = 1
         do i=1,n
           x0(i,1) = x(i)
         enddo
         do i=1,n
           xlow(i,1) =  0.0d0
           xup (i,1) =  1.0d0
         enddo
         nbound = 1

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         do 40 j = 1, ny
            do 30 i = 1, nx
               k = (j-1)*nx + i
               ut = zero
               ub = zero
               ul = zero
               ur = zero
               u = x(k)
               if (i .ne. 1) ul = x(k-1)
               if (i .ne. nx) ur = x(k+1)
               if (j .ne. 1) ub = x(k-nx)
               if (j .ne. ny) ut = x(k+nx)
               uxx = (-ur+two*u-ul)*hydhx
               uyy = (-ut+two*u-ub)*hxdhy
               fvec(k) = uxx + uyy - hxhy*lambda*exp(x(k))
   30       continue
   40    continue
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         do 60 j = 1, n
            do 50 i = 1, n
               fjac(i,j) = zero
   50       continue
   60    continue

c        Evaluate the Jacobian at x.

         do 80 j = 1, ny
            do 70 i = 1, nx
               k = (j-1)*nx + i
               if (i .ne. 1) fjac(k,k-1) = -hydhx
               if (i .ne. nx) fjac(k,k+1) = -hydhx
               if (j .ne. 1) fjac(k,k-nx) = -hxdhy
               if (j .ne. ny) fjac(k,k+nx) = -hxdhy
               fjac(k,k) = two*(hydhx+hxdhy) - hxhy*lambda*exp(x(k))
   70       continue
   80    continue
      end if

      end
