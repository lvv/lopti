      subroutine dhhdfj(n, x, fvec, fjac, ldfjac, task, p, np, ldx, 
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
c     Subroutine dhhdfj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the human heart dipole problem.
c
c     Adapted by       L. Weimann, 
c                      Zuse Institute Berlin (ZIB),
c                      Dep. Numerical Analysis and Modelling
c     File             dhhdfj.f
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
      integer prob
      double precision one, three, twenty, two, zero
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,
     +          twenty=2.0d1)

      integer i
      double precision a, b, c, d, suma, sumb, sumc, sumd, sume, sumf,
     +                 summx, summy, t, temp, ts3vs, tsvs, tt, tv, u,
     +                 us3ws, usws, uu, uw, v, vs3ts, vv, w, ws3us, ww

c     Initialization.
      ifail = 0
      if (np.eq.0) then
         prob = 1
         p(1) = dble(prob)
         np =1
      else
         prob = int(p(1))
      endif

      if (prob .eq. 1) then
         summx = 0.485d0
         summy = -0.0019d0
         suma = -0.0581d0
         sumb = 0.015d0
         sumc = 0.105d0
         sumd = 0.0406d0
         sume = 0.167d0
         sumf = -0.399d0
      else if (prob .eq. 2) then
         summx = -0.69d0
         summy = -0.044d0
         suma = -1.57d0
         sumb = -1.31d0
         sumc = -2.65d0
         sumd = 2.0d0
         sume = -12.6d0
         sumf = 9.48d0
      else if (prob .eq. 3) then
         summx = -0.816d0
         summy = -0.017d0
         suma = -1.826d0
         sumb = -0.754d0
         sumc = -4.839d0
         sumd = -3.259d0
         sume = -14.023d0
         sumf = 15.467d0
      else if (prob .eq. 4) then
         summx = -0.809d0
         summy = -0.021d0
         suma = -2.04d0
         sumb = -0.614d0
         sumc = -6.903d0
         sumd = -2.934d0
         sume = -26.328d0
         sumf = 18.639d0
      else if (prob .eq. 5) then
         summx = -0.807d0
         summy = -0.021d0
         suma = -2.379d0
         sumb = -0.364d0
         sumc = -10.541d0
         sumd = -1.961d0
         sume = -51.551d0
         sumf = 21.053d0
      else
         ifail = -1000
         return
      end if

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         n = 8
         if (prob .eq. 1) then
            x(1) = 0.299d0
            x(2) = 0.186d0
            x(3) = -0.0273d0
            x(4) = 0.0254d0
            x(5) = -0.474d0
            x(6) = 0.474d0
            x(7) = -0.0892d0
            x(8) = 0.0892d0
         else if (prob .eq. 2) then
            x(1) = -0.3d0
            x(2) = -0.39d0
            x(3) = 0.3d0
            x(4) = -0.344d0
            x(5) = -1.2d0
            x(6) = 2.69d0
            x(7) = 1.59d0
            x(8) = -1.5d0
         else if (prob .eq. 3) then
            x(1) = -0.041d0
            x(2) = -0.775d0
            x(3) = 0.03d0
            x(4) = -.047d0
            x(5) = -2.565d0
            x(6) = 2.565d0
            x(7) = -0.754d0
            x(8) = 0.754d0
         else if (prob .eq. 4) then
            x(1) = -0.056d0
            x(2) = -0.753d0
            x(3) = 0.026d0
            x(4) = -0.047d0
            x(5) = -2.991d0
            x(6) = 2.991d0
            x(7) = -0.568d0
            x(8) = 0.568d0
         else if (prob .eq. 5) then
            x(1) = -0.074d0
            x(2) = -0.733d0
            x(3) = 0.013d0
            x(4) = -0.034d0
            x(5) = -3.632d0
            x(6) = 3.632d0
            x(7) = -0.289d0
            x(8) = 0.289d0
         end if
         nx0 = 2
         do i=1,n
           x0(i,1) = x(i)
         enddo
         x0(1,2) = -0.5298801800D-02
         x0(2,2) =  0.4902988018D+00
         x0(3,2) = -0.2920994772D-02
         x0(4,2) =  0.1020994772D-02  
         x0(5,2) = -0.8543213045D-01
         x0(6,2) = -0.9844531015D-01
         x0(7,2) = -0.4190180792D+01
         x0(8,2) = -0.1805682305D-01  
         do i=1,n
           xlow(i,1) = -twenty
           xup (i,1) =  twenty
         enddo
         xup(1,1) = zero
         nbound = 1
         return

      end if

c     Check input arguments for errors.

      if (n .ne. 8) then
         task = 'ERROR: N .NE. 8 IN DHHDFJ'

         return

      end if

      a = x(1)
      b = x(2)
      c = x(3)
      d = x(4)
      t = x(5)
      u = x(6)
      v = x(7)
      w = x(8)
      tv = t*v
      tt = t*t
      vv = v*v
      tsvs = tt - vv
      ts3vs = tt - three*vv
      vs3ts = vv - three*tt
      uw = u*w
      uu = u*u
      ww = w*w
      usws = uu - ww
      us3ws = uu - three*ww
      ws3us = ww - three*uu

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         fvec(1) = a + b - summx
         fvec(2) = c + d - summy
         fvec(3) = t*a + u*b - v*c - w*d - suma
         fvec(4) = v*a + w*b + t*c + u*d - sumb
         fvec(5) = a*tsvs - two*c*t*v + b*usws - two*d*u*w - sumc
         fvec(6) = c*tsvs + two*a*t*v + d*usws + two*b*u*w - sumd
         fvec(7) = a*t*ts3vs + c*v*vs3ts + b*u*us3ws + d*w*ws3us - sume
         fvec(8) = c*t*ts3vs - a*v*vs3ts + d*u*us3ws - b*w*ws3us - sumf
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         fjac(1,1) = one
         fjac(1,2) = one
         fjac(1,3) = zero
         fjac(1,4) = zero
         fjac(1,5) = zero
         fjac(1,6) = zero
         fjac(1,7) = zero
         fjac(1,8) = zero

         fjac(2,1) = zero
         fjac(2,2) = zero
         fjac(2,3) = one
         fjac(2,4) = one
         fjac(2,5) = zero
         fjac(2,6) = zero
         fjac(2,7) = zero
         fjac(2,8) = zero

         fjac(3,1) = t
         fjac(3,2) = u
         fjac(3,3) = -v
         fjac(3,4) = -w
         fjac(3,5) = a
         fjac(3,6) = b
         fjac(3,7) = -c
         fjac(3,8) = -d

         fjac(4,1) = v
         fjac(4,2) = w
         fjac(4,3) = t
         fjac(4,4) = u
         fjac(4,5) = c
         fjac(4,6) = d
         fjac(4,7) = a
         fjac(4,8) = b

         fjac(5,1) = tsvs
         fjac(5,2) = usws
         fjac(5,3) = -two*tv
         fjac(5,4) = -two*uw
         fjac(5,5) = two*(a*t-c*v)
         fjac(5,6) = two*(b*u-d*w)
         fjac(5,7) = -two*(a*v+c*t)
         fjac(5,8) = -two*(b*w+d*u)

         fjac(6,1) = two*tv
         fjac(6,2) = two*uw
         fjac(6,3) = tsvs
         fjac(6,4) = usws
         fjac(6,5) = two*(c*t+a*v)
         fjac(6,6) = two*(d*u+b*w)
         fjac(6,7) = two*(a*t-c*v)
         fjac(6,8) = two*(b*u-d*w)

         fjac(7,1) = t*ts3vs
         fjac(7,2) = u*us3ws
         fjac(7,3) = v*vs3ts
         fjac(7,4) = w*ws3us
         fjac(7,5) = three*(a*tsvs-two*c*tv)
         fjac(7,6) = three*(b*usws-two*d*uw)
         fjac(7,7) = -three*(c*tsvs+two*a*tv)
         fjac(7,8) = -three*(d*usws+two*b*uw)

         fjac(8,1) = -v*vs3ts
         fjac(8,2) = -w*ws3us
         fjac(8,3) = t*ts3vs
         fjac(8,4) = u*us3ws
         fjac(8,5) = three*(c*tsvs+two*a*tv)
         fjac(8,6) = three*(d*usws+two*b*uw)
         fjac(8,7) = three*(a*tsvs-two*c*tv)
         fjac(8,8) = three*(b*usws-two*d*uw)
      end if

      end
