      subroutine HodgkinHuxley1(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Hodgkin-Huxley neuronal membrane model
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             HodgkinHuxley1.f
c   Version          1.0
c   Latest Change    14th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     Hodgkin, A. L. and Huxley, A. F.:
c     A Quantitative Description of Membrane Current and 
c     its Application to Conduction and Excitation in Nerve
c     Journal of Physiology 117: 500-544, 1952
c
c-----------------------------------------------------------------------
      double precision p_Cinv, p_g_Na, p_g_K, p_g_L, p_v_Na,
     $                 p_v_K, p_v_L, Iext
      double precision v, m1, m2, h, mdot, ndot, hdot
      double precision alpham, betam, alphan, betan, alphah, betah
      integer i
c
      alpham(v) = ( 2.5 - 0.1*v) / ( exp( 2.5 - 0.1*v) - 1.)
      betam(v) = 4. * exp( - v / 18.)
      alphan(v) = ( 0.1 - 0.01*v) / ( exp( 1. - 0.1*v) - 1.)
      betan(v) = 0.125 * exp( - v / 80.)
      alphah(v) = 0.07 * exp( - v / 20.)
      betah(v) = 1. / ( exp( 3. - 0.1*v) + 1.)
c
      ifail = 0
      if (np.eq.0) then
		p_Cinv = 1.
		p_g_Na = 120.
		p_g_K  =  36.
		p_g_L  =   0.3
		p_v_Na = 115.
		p_v_K  = -12.
		p_v_L  =  10.6
		p(1) = p_Cinv
		p(2) = p_g_Na
		p(3) = p_g_K
		p(4) = p_g_L
		p(5) = p_v_Na
		p(6) = p_v_K
		p(7) = p_v_L
		np = 7
      else
		p_Cinv = p(1)
		p_g_Na = p(2)
		p_g_K  = p(3)
		p_g_L  = p(4)
		p_v_Na = p(5)
		p_v_K  = p(6)
		p_v_L  = p(7)
      endif
c
      if (task .eq. 'XS') then
        n = 4
        x0( 1, 1) = 0.1
        x0( 2, 1) = 0.1
        x0( 3, 1) = 0.1
        x0( 4, 1) = 0.1
        x0( 1, 2) = 0.e0
        x0( 2, 2) = alpham(0.0d0)/(alpham(0.0d0)+betam(0.0d0))
        x0( 3, 2) = alphan(0.0d0)/(alphan(0.0d0)+betan(0.0d0))
        x0( 4, 2) = alphah(0.0d0)/(alphah(0.0d0)+betah(0.0d0))
        nx0 = 2
        do i=1,n
          xthrsh(i) = 1.0d0
          fthrsh(i) = 1.0d0
        enddo
        rtol = 1.0d-10
        do i=1,n
          xlow(i,1) =  0.0d0
          xup (i,1) =  1.0d0
        enddo
        nbound = 1
        nsol = 0
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
        Iext = 0.
        v  = x(1)
        m1 = x(2)
        m2 = x(3)
        h  = x(4)
        mdot = alpham(v) * ( 1. - m1) - betam( v) * m1
        ndot = alphan(v) * ( 1. - m2) - betan( v) * m2
        hdot = alphah(v) * ( 1. - h) - betah( v) * h
        f(1) = p_Cinv * (Iext -  ( p_g_Na*m1**3*h*(v-p_v_Na) + 
     $                p_g_K*m2**4*(v-p_v_K) + p_g_L*(v-p_v_L) ) )
        f(2) = mdot
        f(3) = ndot
        f(4) = hdot
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        ifail = -999
      endif
      return
      end
