      subroutine DistillColumn(n, x, f, dfdx, ldfjac, task, p, np, ldx, 
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
c   Distillation column problems
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             DistillColumn.f
c   Version          1.0
c   Latest Change    8th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c   Reference:
c     More, J.J.:
c     A Collection of Nonlinear Model Problems.
c     Preprint MCS-P60-0289, 
c     Mathematics and Computer Science Division, 
c     Argonne National Laboratory (1989)
c
c-----------------------------------------------------------------------
      integer nmax, mmax, ndat, nzmax
      parameter ( nmax=40, mmax=3, ndat=2 )
      parameter( nzmax = (11*mmax+6)*(nmax-2)+11*mmax+5 )
      integer irowt(nzmax), icolt(nzmax)
      double precision dfzt(nzmax) 
      integer i, j, iprob, iptyp, nfillt
      integer nd, md
      save nd, md
c
      integer k, ntot
      double precision a(mmax), b(mmax), c(mmax),
     $                 alph(mmax), alphi(mmax), alphii(mmax),
     $                 beta(mmax), betai(mmax), betaii(mmax),
     $                 fl(mmax), fv(mmax),
     $                 tf, bf, df, qf,
     $                 pi(0:nmax)
      common /cdis/ a, b, c, alph, alphi, alphii,
     $                beta, betai, betaii, fl, fv, tf, bf, df, qf, pi,
     $                k, ntot
      save /cdis/
c
      double precision ah(mmax,ndat),bh(mmax,ndat),ch(mmax,ndat),
     $             alh(mmax,ndat), alih(mmax,ndat), aliih(mmax,ndat),
     $             beth(mmax,ndat), betih(mmax,ndat), betiih(mmax,ndat),
     $             flh(mmax,ndat), fvh(mmax,ndat),
     $             tfh(ndat),bfh(ndat),dfh(ndat),qfh(ndat),
     $             pih(0:nmax,ndat)
      data ah(1,1)/9.647/,bh(1,1)/-2998.00/,ch(1,1)/230.66/,
     $     ah(2,1)/9.953/,bh(2,1)/-3448.10/,ch(2,1)/235.88/,
     $     ah(3,1)/9.466/,bh(3,1)/-3347.25/,ch(3,1)/215.31/
      data alh(1,1)/0./,alih(1,1)/37.6/,aliih(1,1)/0./,
     $     alh(2,1)/0./,alih(2,1)/48.2/,aliih(2,1)/0./,
     $     alh(3,1)/0./,alih(3,1)/45.4/,aliih(3,1)/0./
      data beth(1,1)/8425. /,betih(1,1)/24.2/,betiih(1,1)/0./,
     $     beth(2,1)/9395. /,betih(2,1)/35.6/,betiih(2,1)/0./,  
     $     beth(3,1)/10466./,betih(3,1)/31.9/,betiih(3,1)/0./
      data flh(1,1)/30./,flh(2,1)/30./,flh(3,1)/40./,
     $     fvh(1,1)/0./ ,fvh(2,1)/0./ ,fvh(3,1)/0./ ,
     $     tfh(1), bfh(1), dfh(1), qfh(1) /100.,40.,60.,2500000./,
     $     (pih(i,1),i=0,nmax) /41*1./
      data ah(1,2)/18.5751/,bh(1,2)/-3632.649 /,ch(1,2)/239.2/,
     $     ah(2,2)/18.3443/,bh(2,2)/-3841.2203/,ch(2,2)/228. /
      data alh(1,2)/0./,alih(1,2)/15.97/,aliih(1,2)/.0422/,
     $     alh(2,2)/0./,alih(2,2)/18.1 /,aliih(2,2)/0./
      data beth(1,2)/9566.67 /,betih(1,2)/-1.59/,betiih(1,2)/.0422/,
     $     beth(2,2)/10834.67/,betih(2,2)/8.74 /,betiih(2,2)/0.   /
      data flh(1,2)/451.25/,flh(2,2)/684.25/,fvh(1,2)/0./,fvh(2,2)/0./,
     $     tfh(2), bfh(2), dfh(2), qfh(2) /89.,693.37,442.13,8386200./,
     $     (pih(i,2),i=0,7) 
     $     /1210.,1200.,1190.,1180.,1170.,1160.,1150.,1140./
c
      ifail = 0
      if (np.eq.0) then
        iprob = 1
        p(1) = dble(iprob)
        np = 1
      else
        iprob = int(p(1))
      endif
      
      if (task .eq. 'XS') then
c       write(6,'(a,i1)') ' iprob=',iprob
        if (iprob.eq.1) then
c         Hydrocarbon-6
          nd =  6
          md =  3
          k  =  2
          iptyp = 1
        else if (iprob.eq.2) then
c         Hydrocarbon-20
          nd = 20
          md =  3
          k  =  9
          iptyp = 1
        else if (iprob.eq.3) then
c         Methanol-8
          nd =  8
          md =  2
          k  =  2
          iptyp = 2
        else if (iprob.eq.4) then
c         Hydrocarbon-40
          nd = 40
          md =  3
          k  =  19
          iptyp = 1
        else
          ifail = -1000
          return
        endif
        do j=1,md
          a(j)      = ah(j,iptyp)
          b(j)      = bh(j,iptyp)
          c(j)      = ch(j,iptyp)
          alph(j)   = alh(j,iptyp)
          alphi(j)  = alih(j,iptyp)
          alphii(j) = aliih(j,iptyp)
          beta(j)   = beth(j,iptyp)
          betai(j)  = betih(j,iptyp)
          betaii(j) = betiih(j,iptyp)
          fl(j)     = flh(j,iptyp)
          fv(j)     = fvh(j,iptyp)
        enddo
        tf = tfh(iptyp)
        bf = bfh(iptyp)
        df = dfh(iptyp)
        qf = qfh(iptyp)
        do i=0,nd-1
          pi(i) = pih(i,iptyp)
        enddo 
        n = (md+2)*nd-1
c       write(6,'(3(a,i2))') ' inidis: md=',md,'  nd=',nd,'  n=',n
        call inidis(nd-1, md, x0(1,1), x0(nd*md+1,1), x0(nd*(md+1)+1,1),
     $              iprob)
        nx0 = 1
      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.
      if (task .eq. 'F' .or. task .eq. 'FJ') then
c       write(6,'(3(a,i2))') ' fcndis: md=',md,'  nd=',nd
        call fcndis( nd, nd-1, md, x(1), x(nd*md+1), x(nd*(md+1)+1),
     $               f(1), f(md+1), f((nd-1)*md+1), f(nd*md+1),
     $               f((md+1)*nd+1), f((md+1)*nd+2) ) 
      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
        nfillt = nzmax
c       write(6,'(3(a,i2))') ' jacdis: md=',md,'  nd=',nd
        call jacdis( nd-1, md, x(1), x(nd*md+1), x(nd*(md+1)+1),
     $               dfzt, irowt, icolt, nfillt, ifail )
        if (ifail.ne.0) return
        call matcop(nfillt, dfzt, irowt, icolt, n, ldfjac, dfdx, ifail)
      endif
      return
      end
c
      subroutine inidis( nm1, m, x, t, v, nprobh ) 
      integer nm1, m, nprobh
      double precision x(0:nm1, m), t(0:nm1), v(0:*)
      parameter ( nmax=40, mmax=3, ndata=4 )
      double precision x0(0:nmax,mmax,ndata), t0(0:nmax,ndata),
     $                 v0(0:nmax,ndata)
      data ((x0(i,j,1),j=1,3),i=0,5)  /0.,.2,.9,0.,.2,.8,.05,.3,.8,
     $     .1,.3,.6,.3,.5,.3,.6,.6,0. /,
     $     (t0(i,1),i=0,5)/6*100./,
     $     (v0(i,1),i=0,4)/5*300./
      data ((x0(i,j,2),j=1,3),i=0,19) /0.,.3,1.,0.,.3,.9,.01,.3,.9,
     $     .02,.4,.8,.05,.4,.8,.07,.45,.8,.09,.5,.7,.1,.5,.7,.15,.5,.6,
     $     .2,.5,.6,.25,.6,.5,.3,.6,.5,.35,.6,.5,.4,.6,.4,.4,.7,.4,
     $     .42,.7,.3,.45,.75,.3,.45,.75,.2,.5,.8,.1,.5,.8,.0/,
     $     (t0(i,2),i=0,19)/20*100./,
     $     (v0(i,2),i=0,18)/19*300./
      data ((x0(i,j,3),j=1,2),i=0,7)  /
     $     .09203,.908,.1819,.8181,.284,.716,.3051,.6949,
     $     .3566,.6434,.468,.532,.6579,.3421,.8763,.1237 /,
     $     (t0(i,3),i=0,7) /120.,110.,100.,88.,86.,84.,80.,76./,
     $     (v0(i,3),i=0,6) 
     $       /886.37,910.01,922.52,926.45,935.56,952.83,975.73/
      data ((x0(i,j,4),j=1,3),i=0,39) /0., .3, 1., 0., .3, 1., 0., .3,
     $      .9, 0., .3, .9, .01, .3, .9, .01, .3, .9,
     $      .02, .4, .8, .02, .4, .8, .05, .4, .8, .05, .4, .8, .07, 
     $      .45, .8, .07, .45, .8,
     $      .09, .5, .7, .09, .5, .7, .1, .5, .7, .1, .5, .7, .15, .5,
     $      .6, .15, .5, .6,
     $      .2, .5, .6, .2, .5, .6, .25, .6, .5, .25, .6, .5, .3, .6,
     $      .5, .3, .6, .5,
     $      .35, .6, .5, .35, .6, .5, .4, .6, .4, .4, .6, .4, .4, .7,
     $      .4, .4, .7, .4,
     $      .42, .7, .3, .42, .7, .3, .45, .75, .3, .45, .75, .3, .45,
     $      .75, .2, .45, .75, .2,
     $      .5, .5, .8, .8, .1, .1, .5, .5, .8, .8, .0, .0/,
     $     (t0(i,4),i=0,39)/40*100./,
     $     (v0(i,4),i=0,38)/39*300./
      do 10 i=0,nm1
        do 11 j=1,m
          x(i,j) = x0(i,j,nprobh)
11      continue
10    continue
      do 20 i=0,nm1
        t(i) = t0(i,nprobh)
20    continue
      do 30 i=0,nm1-1
        v(i) = v0(i,nprobh)
30    continue
      return
      end
c
      subroutine fcndis(n, nm1, m, x, t, v, fz1, fz2, fz3, fz7, fz8,
     $                  fz9)
      implicit double precision (a-h,o-z)
      integer n, nm1, m
      double precision x(0:nm1, m), t(0:nm1), v(0:*), fz1(m), fz2(m,*),
     $                 fz3(m), fz7(n), fz8, fz9(*)
c
c     problem data common:
c
      parameter ( nmax=40, mmax=3 )
      integer k, ntot
      double precision a(mmax), b(mmax), c(mmax),
     $                 alph(mmax), alphi(mmax), alphii(mmax),
     $                 beta(mmax), betai(mmax), betaii(mmax),
     $                 fl(mmax), fv(mmax),
     $                 tf, bf, df, qf,
     $                 pi(0:nmax)
      common /cdis/ a, b, c, alph, alphi, alphii,
     $                beta, betai, betaii, fl, fv, tf, bf, df, qf, pi,
     $                k, ntot
c
c     internal storage of subroutine fcnint
c
      double precision l(nmax), hh(0:nmax), h(0:nmax), hhf, hf, kf, y
c
c     (2.10)
      kf(i,j) = dexp( a(j)+b(j)/(c(j)+t(i)) ) / pi(i)
c     (2.6)
      y(i,j) = kf(i,j)*x(i,j)
c     (2.4)
      do 10 i=1,k
        l(i) = v(i-1)+bf
10    continue
c     (2.5)
      do 20 i=k+1,n-1
        l(i) = v(i-1)-df
20    continue
c
      do 30 j=1,m
c       (2.1)
        fz1(j) = v(0)*y(0,j)+bf*x(0,j)-l(1)*x(1,j)
c       (2.2)
        do 40 i=1,n-2
          fz2(j,i) = v(i-1)*y(i-1,j)+l(i+1)*x(i+1,j)
     $              -v(i)*y(i,j)-l(i)*x(i,j)
40      continue
        fz2(j,k) = fz2(j,k) + fl(j)
        fz2(j,k+1) = fz2(j,k+1) + fv(j)
c       (2.3)
        fz3(j) = y(n-2,j)-x(n-1,j)
30    continue
c     (2.7)
      do 100 i=0,n-1      
        sh = 0.0d0
        do 110 j=1,m
          sh = sh + y(i,j) 
110     continue
        fz7(i+1) = sh - 1.0d0
100   continue
c     (2.11)
      do 120 i=0,n-1
        sh = 0.0d0
        do 121 j=1,m
          sh = sh + x(i,j)*(alph(j)+alphi(j)*t(i)+alphii(j)*t(i)*t(i))
121     continue
        h(i) = sh
120   continue
c     (2.12)
      do 130 i=0,n-2
        sh = 0.0d0
        do 131 j=1,m
          sh = sh + y(i,j)*(beta(j)+betai(j)*t(i)+betaii(j)*t(i)*t(i))
131     continue
        hh(i) = sh
130   continue
c     (2.13)
      hf = 0.0d0
      do 140 j=1,m
        hf = hf + fl(j)*(alph(j)+alphi(j)*tf+alphii(j)*tf*tf)
140   continue
c     (2.14)
      hhf = 0.0d0
      do 150 j=1,m
        hhf = hhf + fv(j)*(beta(j)+betai(j)*tf+betaii(j)*tf*tf)
150   continue
c     (2.8)
      fz8 = v(0)*hh(0)+bf*h(0)-l(1)*h(1)-qf
c     (2.9)
      do 170 i=1,n-2
        ih=i
        fz9(i) = v(i-1)*hh(i-1)+l(i+1)*h(i+1)-v(ih)*hh(i)-l(i)*h(i)
170   continue
      fz9(k) = fz9(k)+hf
      fz9(k+1) = fz9(k+1)+hhf
      return
      end
c
      subroutine jacdis(nm1, m, x, t, v, dfz, irow, icol, nfill, ifail )
      implicit double precision (a-h,o-z)
      integer nm1, m, nfill
      double precision x(0:nm1, m), t(0:nm1), v(0:*), dfz(nfill)
      integer irow(nfill), icol(nfill)
c
c     problem data common:
c
      parameter ( nmax=40, mmax=3 )
      integer k, ntot
      double precision a(mmax), b(mmax), c(mmax),
     $                 alph(mmax), alphi(mmax), alphii(mmax),
     $                 beta(mmax), betai(mmax), betaii(mmax),
     $                 fl(mmax), fv(mmax),
     $                 tf, bf, df, qf,
     $                 pi(0:nmax)
      common /cdis/ a, b, c, alph, alphi, alphii,
     $                beta, betai, betaii, fl, fv, tf, bf, df, qf, pi,
     $                k, ntot
c
c     internal storage of subroutine jacint
c
      double precision l(nmax), hh(0:nmax), h(0:nmax), hhf, hf, kf, kfi,
     $                 y, hhi(0:nmax), hi(0:nmax)
c
      kf(i,j) = dexp( a(j)+b(j)/(c(j)+t(i)) ) / pi(i)
      kfi(i,j) = -b(j)/(c(j)+t(i))**2 * kf(i,j)
      y(i,j) = kf(i,j)*x(i,j)
c
c     column index mappings:
c
      ix(i,j) = (j-1)*n+i+1
      it(i) = m*n+i+1
      iv(i) = (m+1)*n+i+1
c
      n=nm1+1
c     (2.4)
      do 10 i=1,k
        l(i) = v(i-1)+bf
10    continue
c     (2.5)
      do 20 i=k+1,n-1
        l(i) = v(i-1)-df
20    continue
c
      nz = 0
c
      do 50 j=1,m
c       (2.1)
        ir = j
c
        nz = nz+1
        dfz(nz)  = v(0)*kf(0,j) + bf
        irow(nz) = ir
        icol(nz) = ix(0,j)
c
        nz = nz+1
        dfz(nz)  = -v(0)-bf
        irow(nz) = ir
        icol(nz) = ix(1,j)
c
        nz = nz+1
        dfz(nz)  = v(0)*kfi(0,j)*x(0,j)
        irow(nz) = ir
        icol(nz) = it(0)
c
        nz = nz+1
        dfz(nz)  = y(0,j)-x(1,j)
        irow(nz) = ir
        icol(nz) = iv(0)
c
c       (2.2)
        do 55 i=1,n-2
          ir = i*m+j
c
          nz = nz+1
          dfz(nz)  = v(i-1)*kf(i-1,j)
          irow(nz) = ir
          icol(nz) = ix(i-1,j)
c
          nz = nz+1
          dfz(nz)  = l(i+1)
          irow(nz) = ir
          icol(nz) = ix(i+1,j)
c
          nz = nz+1
          dfz(nz)  = -v(i)*kf(i,j)-l(i)
          irow(nz) = ir
          icol(nz) = ix(i,j)
c
          nz = nz+1
          dfz(nz)  = v(i-1)*kfi(i-1,j)*x(i-1,j)
          irow(nz) = ir
          icol(nz) = it(i-1)
c
          nz = nz+1
          dfz(nz)  = -v(i)*kfi(i,j)*x(i,j)
          irow(nz) = ir
          icol(nz) = it(i)
c
          nz = nz+1
          dfz(nz)  = y(i-1,j)-x(i,j)
          irow(nz) = ir
          icol(nz) = iv(i-1)
c
          nz = nz+1
          dfz(nz)  = x(i+1,j)-y(i,j)
          irow(nz) = ir
          icol(nz) = iv(i)
c
55      continue
c     (2.3)
        ir = (n-1)*m+j
c
        nz = nz+1
        dfz(nz)  = kf(n-2,j)
        irow(nz) = ir
        icol(nz) = ix(n-2,j)
c
        nz = nz+1
        dfz(nz)  = -1.0d0
        irow(nz) = ir
        icol(nz) = ix(n-1,j)
c
        nz = nz+1
        dfz(nz)  = kfi(n-2,j)*x(n-2,j)
        irow(nz) = ir
        icol(nz) = it(n-2)
c     
50    continue
c
c     (2.7)
      do 60 i=0,n-1
        ir = n*m+i+1
c
        do 62 j=1,m
c
          nz = nz+1
          dfz(nz)  = kf(i,j)
          irow(nz) = ir
          icol(nz) = ix(i,j)
c
62      continue
c
        sh = 0.0d0
        do 65 j=1,m
          sh = sh + kfi(i,j)*x(i,j)
65      continue
c
        nz = nz+1
        dfz(nz)  = sh
        irow(nz) = ir
        icol(nz) = it(i)
c
60    continue
c
c     (2.11)
      do 120 i=0,n-1
        sh = 0.0d0
        do 121 j=1,m
          sh = sh + x(i,j)*(alph(j)+alphi(j)*t(i)+alphii(j)*t(i)*t(i))
121     continue
        h(i) = sh
120   continue
c     (2.12)
      do 130 i=0,n-2
        sh = 0.0d0
        do 131 j=1,m
          sh = sh + y(i,j)*(beta(j)+betai(j)*t(i)+betaii(j)*t(i)*t(i))
131     continue
        hh(i) = sh
130   continue
c     (2.13)
      hf = 0.0d0
      do 140 j=1,m
        hf = hf + fl(j)*(alph(j)+alphi(j)*tf+alphii(j)*tf*tf)
140   continue
c     (2.14)
      hhf = 0.0d0
      do 150 j=1,m
        hhf = hhf + fv(j)*(beta(j)+betai(j)*tf+betaii(j)*tf*tf)
150   continue
c     derivatives of (2.11)
      do 160 i=0,n-1
        sh = 0.0d0
        do 161 j=1,m
          sh = sh + x(i,j)*(alphi(j)+2.0d0*alphii(j)*t(i))
161     continue
        hi(i) = sh
160   continue
c     derivatives of (2.12)
      do 170 i=0,n-2
        sh = 0.0d0
        do 171 j=1,m
          sh = sh + y(i,j)*(betai(j)+2.0d0*betaii(j)*t(i)-b(j)*
     $       (beta(j)+betai(j)*t(i)+betaii(j)*t(i)*t(i))/(c(j)+t(i))**2)
171     continue
        hhi(i) = sh
170   continue
c
c     (2.8)
      ir = (m+1)*n+1 
c 
      do 200 j=1,m
c
        nz = nz+1
        dfz(nz)  =   v(0)*kf(0,j)*
     $               (beta(j)+betai(j)*t(0)+betaii(j)*t(0)*t(0))
     $             + bf*(alph(j)+alphi(j)*t(0)+alphii(j)*t(0)*t(0))
        irow(nz) = ir
        icol(nz) = ix(0,j)
c
        nz = nz+1
        dfz(nz)  = -(v(0)+bf)*
     $              (alph(j)+alphi(j)*t(1)+alphii(j)*t(1)*t(1))
        irow(nz) = ir
        icol(nz) = ix(1,j)
c
200   continue 
c
      nz = nz+1
      dfz(nz)  = v(0)*hhi(0)+bf*hi(0)
      irow(nz) = ir
      icol(nz) = it(0)
c
      nz = nz+1
      dfz(nz)  = -(v(0)+bf)*hi(1)
      irow(nz) = ir
      icol(nz) = it(1)
c
      nz = nz+1
      dfz(nz)  = hh(0)-h(1)
      irow(nz) = ir
      icol(nz) = iv(0)
c
c     (2.9)
      do 250 i=1,n-2
        ir = (m+1)*n+1+i
        do 255 j=1,m
c
          nz = nz+1
          dfz(nz)  = v(i-1)*kf(i-1,j)*
     $               (beta(j)+betai(j)*t(i-1)+betaii(j)*t(i-1)*t(i-1))
          irow(nz) = ir
          icol(nz) = ix(i-1,j)
c
          nz = nz+1
          dfz(nz)  = l(i+1)*
     $               (alph(j)+alphi(j)*t(i+1)+alphii(j)*t(i+1)*t(i+1))
          irow(nz) = ir
          icol(nz) = ix(i+1,j)
c
          nz = nz+1
          dfz(nz)  = -v(i)*kf(i,j)*
     $                (beta(j)+betai(j)*t(i)+betaii(j)*t(i)*t(i))
     $               -l(i)*(alph(j)+alphi(j)*t(i)+alphii(j)*t(i)*t(i))
          irow(nz) = ir
          icol(nz) = ix(i,j)
c
255   continue
c
        nz = nz+1
        dfz(nz)  = v(i-1)*hhi(i-1)
        irow(nz) = ir
        icol(nz) = it(i-1)
c
        nz = nz+1
        dfz(nz)  = l(i+1)*hi(i+1)
        irow(nz) = ir
        icol(nz) = it(i+1)
c
        nz = nz+1
        dfz(nz)  = -v(i)*hhi(i)-l(i)*hi(i)
        irow(nz) = ir
        icol(nz) = it(i)
c
        nz = nz+1
        dfz(nz)  = hh(i-1)-h(i)
        irow(nz) = ir
        icol(nz) = iv(i-1)
c
        nz = nz+1
        dfz(nz)  = h(i+1)-hh(i)
        irow(nz) = ir
        icol(nz) = iv(i)
c
250   continue
c
      ifail = 0
      if (nz.gt.nfill) ifail=nfill-nz
      nfill = nz
      return
      end
c
      subroutine matcop(nz, dfzt, irowt, icolt, nn, ldjac, dfzf, ifail)
      integer nz, nn, ldjac, ifail
      integer irowt(nz), icolt(nz)
      double precision dfzt(nz), dfzf(ldjac,nn)
c
      ifail = 0
      do 5 i=1,ldjac
        do 6 j=1,nn
          dfzf(i,j)=0.0d0
6       continue
5     continue
      do 10 l=1,nz
        ir = irowt(l)
        ic = icolt(l) 
        if ( ir.ge.1 .and. ir.le.nn .and. ic.ge.1 .and. ic.le.nn ) then 
          dfzf(ir,ic) = dfzt(l)
        else
          ifail = -1000
          return
        endif
10    continue
      return
      end

