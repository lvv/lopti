      subroutine dficfj(n, x, fvec, fjac, ldfjac, task, p, np, ldx, 
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
c     Subroutine dficfj
c
c     This subroutine computes the function and Jacobian matrix of the
c     flow in a channel problem.
c
c
c     Adapted by       L. Weimann, 
c                      Zuse Institute Berlin (ZIB),
c                      Dep. Numerical Analysis and Modelling
c     File             dficfj.f
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
      integer nint
      double precision r
      integer bc, cpts, deg, dim, npi
      parameter (bc=2,cpts=4,deg=4,dim=deg+cpts-1,npi=cpts+deg)
      double precision one, six, three, twelve, two, zero
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0,
     +          twelve=12.0d0)

      logical feval, jeval
      integer i, j, k, m, eqn, var
      double precision h, hm, nf, rhoijh, xt
      double precision dw(deg+1,cpts+deg),
     +                 rhnfhk(cpts,0:dim,0:dim,0:deg), rho(cpts),
     +                 w(deg+1)

      data (rho(i),i=1,cpts)/0.694318413734436035d-1,
     +     0.330009490251541138d0, 0.669990539550781250d0,
     +     0.930568158626556396d0/

      ifail = 0
      if (np.eq.0) then
         nint = 10
         r = 100.0d0
         p(1) = dble(nint)
         p(2) = r
         np = 2
      else
         nint = int(p(1))
         r = p(2)
      endif

c     Initialization.

      h = one/dble(nint)

c     Compute the standard starting point if task = 'XS'

      if (task .eq. 'XS') then
         n = 8*nint
c        The standard starting point corresponds to the solution of the
c        flow in a channel problem with R = 0.

         xt = zero
         do 20 i = 1, nint
            var = (i-1)*npi
            x(var+1) = xt*xt*(three-two*xt)
            x(var+2) = six*xt*(one-xt)
            x(var+3) = six*(one-two*xt)
            x(var+4) = -twelve
            do 10 j = 1, cpts
               x(var+deg+j) = zero
   10       continue
            xt = xt + h
   20    continue
         nx0 = 1
         do i=1,n
           x0(i,1) = x(i)
         enddo
         do i=1,n
           xlow(i,1) = -2.0d1
           xup (i,1) =  2.0d1
         enddo
         nbound = 1
         return

      end if

c     Check input arguments for errors.

      if (n .ne. 8*nint) then
         task = 'ERROR: N .NE. 8*NINT IN DFICFJ'

         return

      end if
      
c     Store all possible combinations of rho, h, and n factorial.

      hm = one
      do 60 m = 0, deg
         do 50 i = 1, cpts
            rhoijh = hm
            do 40 j = 0, dim
               nf = one
               do 30 k = 0, dim
                  rhnfhk(i,j,k,m) = rhoijh/nf
                  nf = nf*dble(k+1)
   30          continue
               rhoijh = rhoijh*rho(i)
   40       continue
   50    continue
         hm = hm*h
   60 continue

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         feval = .true.
      else
         feval = .false.
      end if
      if (task .eq. 'J' .or. task .eq. 'FJ') then
         jeval = .true.
      else
         jeval = .false.
      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

c     Initialize arrays.

      do 80 j = 1, n
         if (feval) fvec(j) = zero
         if (jeval) then
            do 70 i = 1, n
               fjac(i,j) = zero
   70       continue
         end if
   80 continue
      do 100 k = 1, npi
         do 90 j = 1, deg + 1
            dw(j,k) = zero
   90    continue
  100 continue

c     Set up the boundary equations at t = 0.  u(0) = 0, u'(0) = 0.

      if (feval) then
         fvec(1) = x(1)
         fvec(2) = x(2)
      end if
      if (jeval) then
         fjac(1,1) = one
         fjac(2,2) = one
      end if

      do 220 i = 1, nint
         var = (i-1)*npi

c        Set up the collocation equations.

         eqn = var + bc
         do 150 k = 1, cpts
            do 130 m = 1, deg + 1
               w(m) = zero
               do 110 j = m, deg
                  w(m) = w(m) + rhnfhk(k,j-m,j-m,j-m)*x(var+j)
                  dw(m,j) = rhnfhk(k,j-m,j-m,j-m)
  110          continue
               do 120 j = 1, cpts
                  w(m) = w(m) + rhnfhk(k,deg+j-m,deg+j-m,deg-m+1)*
     +                   x(var+deg+j)
                  dw(m,deg+j) = rhnfhk(k,deg+j-m,deg+j-m,deg-m+1)
  120          continue
  130       continue
            if (feval) fvec(eqn+k) = w(5) - r*(w(2)*w(3)-w(1)*w(4))
            if (jeval) then
               do 140 j = 1, npi
                  fjac(eqn+k,var+j) = dw(5,j) -
     +                                r*(dw(2,j)*w(3)+w(2)*dw(3,j)-
     +                                dw(1,j)*w(4)-w(1)*dw(4,j))
  140          continue
            end if
  150    continue

c        Set up the continuity equations.

         eqn = var + bc + cpts
         do 180 m = 1, deg
            w(m) = zero
            do 160 j = m, deg
               w(m) = w(m) + rhnfhk(1,0,j-m,j-m)*x(var+j)
               dw(m,j) = rhnfhk(1,0,j-m,j-m)
  160       continue
            do 170 j = 1, cpts
               w(m) = w(m) + rhnfhk(1,0,deg+j-m,deg-m+1)*x(var+deg+j)
               dw(m,deg+j) = rhnfhk(1,0,deg+j-m,deg-m+1)
  170       continue
  180    continue
         if (i .eq. nint) go to 230
         if (feval) then
            do 190 m = 1, deg
               fvec(eqn+m) = x(var+cpts+deg+m) - w(m)
  190       continue
         end if
         if (jeval) then
            do 210 m = 1, deg
               fjac(eqn+m,var+cpts+deg+m) = one
               do 200 j = 1, npi
                  fjac(eqn+m,var+j) = -dw(m,j)
  200          continue
  210       continue
         end if
  220 continue

c     Set up the boundary equations at t = 1.  u(1) = 1, u'(1) = 0.

  230 continue
      if (feval) then
         fvec(n-1) = w(1) - one
         fvec(n) = w(2)
      end if
      if (jeval) then
         var = n - npi
         do 240 j = 1, npi
            fjac(n-1,var+j) = dw(1,j)
            fjac(n,var+j) = dw(2,j)
  240    continue
      end if

      end
