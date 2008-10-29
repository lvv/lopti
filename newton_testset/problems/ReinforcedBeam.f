      subroutine ReinforcedBeam(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Model of a typical rectangular singly reinforced concrete beam  
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             ReinforcedBeam.f
c   Version          1.0
c   Latest Change    13th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     C.J. Egelhoff, D.M. Blackketter and J.L. Benson:
c     Algorithms for Solving Nonlinear Equation Systems Assist
c      Students to Become Better Problem Solvers.
c     29th ASEE/IEEE Frontiers in Education Conference, p.12a4-18ff. 
c
c-----------------------------------------------------------------------
      integer i
      double precision p_b, p_h, p_d, p_Es, p_Ec, p_fy, p_j, p_ph, p_FS,
     $                 p_Mcrack
      double precision  Ac, Iconc, m1, FR, ytop, Ibeam, ytens, a, Mnom,
     $                  rho, k, Mservice, pMnom, fc, As, fyc
c
      ifail = 0
      if (np.eq.0) then
        p_b = 24.  
        p_h = 18.
        p_d = 15.5
        p_Es = 29e6
        p_Ec = 3.6e6
        p_fy = 60e3
        p_j = 0.855
        p_ph = 0.9
        p_FS = 2.75
        p_Mcrack = 787732.
        p(1) = p_b
        p(2) = p_h
        p(3) = p_d
        p(4) = p_Es
        p(5) = p_Ec
        p(6) = p_fy
        p(7) = p_j
        p(8) = p_ph
        p(9) = p_FS
        p(10) = p_Mcrack
      else
         p_b =  p(1) 
         p_h =  p(2) 
         p_d =  p(3) 
         p_Es = p(4) 
         p_Ec = p(5) 
         p_fy = p(6) 
         p_j =  p(7) 
         p_ph = p(8) 
         p_FS = p(9) 
         p_Mcrack = p(10) 
      endif
c
      if (task .eq. 'XS') then
        n = 13
c       initial guess
        do i=1,n
          x0(i,1) = 1.0d0
        enddo
c       an initial guess close to the solution
        x0( 1,2) = 474.
        x0( 2,2) = 9.73
        x0( 3,2) = 13726.
        x0( 4,2) = 8.26
        x0( 5,2) = 5.73
        x0( 6,2) = 5.91e6
        x0( 7,2) = 0.02097
        x0( 8,2) = 0.436
        x0( 9,2) = 1.93e6
        x0(10,2) = 5.32e6
        x0(11,2) = 1800.
        x0(12,2) = 7.8
        x0(13,2) = 4000.
c
        x0( 1,3) = 100.
        x0( 2,3) = 10.
        x0( 3,3) = 10000.
        x0( 4,3) = 10.
        x0( 5,3) = 10.
        x0( 6,3) = 1.e6
        x0( 7,3) = 1.e-2
        x0( 8,3) = 1.
        x0( 9,3) = 1.e6
        x0(10,3) = 1.e6
        x0(11,3) = 1000.
        x0(12,3) = 10.
        x0(13,3) = 1000.
        nx0 = 3
        do i=1,n
          xthrsh(i) = 1.0d0
          fthrsh(i) = 1.0d0
        enddo
        rtol = 1.0d-9
c       range of interest/solution
c       that is o.k.
        xlow( 1,1) = 100.
        xlow( 2,1) = 1.
        xlow( 3,1) = 5000.
        xlow( 4,1) = 1.
        xlow( 5,1) = 1.
        xlow( 6,1) = 2.e6
        xlow( 7,1) = 1.e-2
        xlow( 8,1) = 0.1
        xlow( 9,1) = 1.e6
        xlow(10,1) = 1.e6
        xlow(11,1) = 500.
        xlow(12,1) = 1.
        xlow(13,1) = 1000.
c
        xup ( 1,1) = 1000.
        xup ( 2,1) = 50.
        xup ( 3,1) = 20000.
        xup ( 4,1) = 10.
        xup ( 5,1) = 50.
        xup ( 6,1) = 9.e6
        xup ( 7,1) = 1.e-1
        xup ( 8,1) = 1.
        xup ( 9,1) = 8.e6
        xup (10,1) = 9.e6
        xup (11,1) = 5000.
        xup (12,1) = 50.
        xup (13,1) = 6000.
        do i=1,n
          x0(i,4) = xlow(i,1)
        enddo
c       that is quite hard (local method diverges)
		xlow( 1,2) = 10.
		xlow( 2,2) = 1.
		xlow( 3,2) = 1000.
		xlow( 4,2) = 1.
		xlow( 5,2) = 1.
		xlow( 6,2) = 1.e6
		xlow( 7,2) = 1.e-3
		xlow( 8,2) = 0.01
		xlow( 9,2) = 1.e6
		xlow(10,2) = 1.e6
		xlow(11,2) = 100.
		xlow(12,2) = 1.
		xlow(13,2) = 1000.
c
        xup ( 1,1) = 1000
        xup ( 2,1) = 100
        xup ( 3,1) = 100000
        xup ( 4,1) = 10
        xup ( 5,1) = 100
        xup ( 6,1) = 1.e7
        xup ( 7,1) = 1.e-1
        xup ( 8,1) = 1.
        xup ( 9,1) = 1.e7
        xup (10,1) = 1.e7
        xup (11,1) = 10000
        xup (12,1) = 100
        xup (13,1) = 10000
        nbound = 2
        do i=1,n
          x0(i,5) = xlow(i,2)
        enddo
        nx0 = 5
      endif
c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        Ac = p_b*p_h 
        Iconc = p_b*p_h**3/12.0
c       Iconc = p_b*p_h*p_h*p_h/12.0
        m1 = p_Es/p_Ec
        FR = x(1)
        ytop = x(2)
        Ibeam = x(3)
        ytens = x(4)
        a = x(5)
        Mnom = x(6)
        rho = x(7)
        k = x(8)
        Mservice = x(9)
        pMnom = x(10)
        fc = x(11)
        As = x(12)
        fyc = x(13)
        if (fyc.lt.0.0d0 .or. 2*rho*m1 + (rho*m1)**2.lt.0.0d0) then
          ifail=1
          return
        endif
        f(1) = ytop -(Ac*p_h/2 + (m1-1)*As*p_d)/(Ac+(m1-1)*As)
        f(2) = FR - 7.5*dsqrt(fyc)
        f(3) = rho - As/(p_b*p_d)
        f(4) = k - (-rho*m1 + dsqrt(2*rho*m1 + (rho*m1)**2))
        f(5) = ytens - (p_h - ytop)
        f(6) = fc - 0.45*fyc
        f(7) = p_FS - pMnom/Mservice
        f(8) = Mservice - ((fc/2)*k*p_b*p_j*p_d**2)
        f(9) = Mnom - (As*p_fy*(p_d-a/2))
        f(10) = a - (As*p_fy)/(0.85*fyc*p_b)
        f(11) = pMnom - p_ph*Mnom
        f(12) = p_Mcrack - FR*Ibeam/ytens
        f(13) = Ibeam - (Iconc + Ac*(p_h/2-ytop)**2 + 
     $            (m1-1)*As*(p_d-ytop)**2)
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ifail = -999
      endif
      return
      end
