      program testset
      implicit none
c --------------------------------------------------------------------
c
c*  Title
c
c     testset: a testset program to test solver routines for
c              systems of nonlinear equations with various test 
c              problems.
c 
c*  Written by        L. Weimann 
c*  Category          F2a. - Systems of nonlinear equations
c*  Keywords          Nonlinear equations, Newton methods
c*  Version           1.1
c*  Revision          June 2006
c*  Latest Change     June 2006
c*  Code              Fortran 77, Double Precision
c*  Environment       Standard Fortran 77 environment on PC's,
c                     workstations and hosts.
c*  Copyright     (c) Konrad-Zuse-Zentrum fuer
c                     Informationstechnik Berlin (ZIB)
c                     Takustrasse 7, D-14195 Berlin-Dahlem
c                     phone : + 49/30/84185-0
c                     fax   : + 49/30/84185-125
c*  Contact           Lutz Weimann
c                     ZIB, Division Scientific Computing, 
c                          Department Numerical Analysis and Modelling
c                     phone : + 49/30/84185-185
c                     fax   : + 49/30/84185-107
c                     e-mail: weimann@zib.de
c
c  ---------------------------------------------------------------
      integer max_n, max_con, max_par, max_x0, max_bound, max_sol,
     $        max_iter, n_solver
      parameter ( max_n     = 200,
     $            max_con   =  10, 
     $            max_par   = 100,
     $            max_x0    =  10,
     $            max_bound =  10,
     $            max_sol   =  20,
     $            max_iter  = 100,
     $            n_solver  =   6)
      character problem*48, directive*80, solver*16, control*16
      integer n, nx0, nbound, nsol, ifail, ifail2
      integer ix0, ibound, icon, isolver, i, ii, ie, ijac, lsolver,
     $        gsolver, nonlin, l1, derfun, gnonlin
      double precision x0(max_n,max_x0),xlow(max_n,max_bound), 
     $                 xup(max_n,max_bound), xsol(max_n,max_sol),
     $                 xthrsh(max_n), fthrsh(max_n), rtol, grtol,
     $                 gxthr, gfthr
      double precision f(max_n)
c     i/o units
      integer lucon, lumon, lusol, lusum
      parameter(lucon=10, lumon=20, lusol=21, lusum=30)
c     options and workspace for nonlinear solver
      integer liopt, lrwk, liwk
      parameter (liopt=50, lrwk=70000, liwk=1000)
      integer iopt(liopt), iwk(liwk), info(5)
      double precision rwk(lrwk), xscale(max_n), fscale(max_n), x(max_n)
      character*80 fname(6)
      double precision fnorm 
      logical qinit
      external f_nleq, j_nleq
c     statistics
      integer niter, nfcn, njac, nsolvd, nfaild, nfctot, njctot
      real stime, etime
c
      character*32 modnam
      common /modul/ modnam
c
      integer np, ncon, nsave
      double precision p(max_par), conxl(max_n,max_con),
     $                 conxu(max_n,max_con)
      common /newlab/ np, p, ncon, conxl, conxu, nsave
c
      integer nfev, njev
      common /refnum/ nfev, njev
c
      qinit = .false.
      nsolvd = 0
      nfaild = 0
      nfctot = 0
      njctot = 0
      write(6,*) ' Enter control file id (name without ext. ".control")'
      read(5,'(a16)') control
      if (control.eq.' ') control = 'testset'
      l1 = index(control,' ')-1
      if (l1.eq.-1) l1=16
      open(lucon,file=control(:l1)//'.control',status='old',err=9100)
      open(lusum,file=control(:l1)//'.sumout',status='unknown',err=9120)
      write(lusum,'(a7,41x,4(1x,a5),6x,a7,2x,a8)') 
     $          'problem','rcode','#iter',
     $          ' #fcn',' #jac','resnorm','cpu-time'
      grtol = 0.0d0
      gsolver = 0
      gxthr = -1.0d0
      gfthr = -1.0d0
      gnonlin = 3
      read(lucon,'(a32,a48)',end=9000) modnam, problem
      if (modnam.ne.'*global ') goto 1050
600   continue
      read(lucon,'(a80)',end=9000) directive
      if (directive(1:1).eq.'#') goto 600
      if (directive(1:5).eq.' end ') then
        goto 1000
      else if (directive(1:8).eq.' solver ') then
        read(directive(9:80),*) gsolver
        if (gsolver.gt.n_solver) then
          write(6,*) ' error in definition: problem = ',problem
          write(6,'(a,i2,a,i2)') 'gsolver=',gsolver,' .gt. n_solver=',
     $                           n_solver
          gsolver = 0
        endif
      else if (directive(1:6).eq.' rtol ') then
        read(directive(7:80),*) grtol
      else if (directive(1:9).eq.' xthresh ') then
        read(directive(10:80),*) gxthr
      else if (directive(1:9).eq.' fthresh ') then
        read(directive(10:80),*) gfthr
      else if (directive(1:8).eq.' nonlin ') then
        read(directive(9:80),*) gnonlin
      else if (directive(1:1).ne.'#') then
        write(6,*) ' error in *global definition'
        write(6,*) ' invalid directive: ',directive(1:60)
      endif
      goto 600
1000  continue
      read(lucon,'(a32,a48)',end=9000) modnam, problem
      if (modnam(1:1).eq.' ') goto 1000
1050  continue
      if ( gsolver.ge.3 .and. .not.qinit ) then
        fname(1)=control(:l1)//'.monout'
        fname(2)=control(:l1)//'.monout'
        fname(3)=control(:l1)//'.solout'
        fname(4)=' '
        fname(5)=' '
        fname(6)=' '
        call cinit(fname)
      else if ( .not.qinit ) then
        open(lumon,file=control(:l1)//'.monout',status='unknown',
     1       err=9110)
        open(lusol,file=control(:l1)//'.solout',status='unknown',
     1       err=9130)
      endif
      qinit = .true.
      do i=1,max_n
        xthrsh(i) = 1.0d0
        fthrsh(i) = 1.0d0
      enddo
      n = 0
      np = 0
      nonlin = gnonlin
      ibound = 1
      call init(n, x0, xlow, xup, xsol, xthrsh, fthrsh, 
     $          nx0, nbound, nsol, rtol, ifail)
      if (ifail.ne.0) then
          write(6,*) ' error in init: problem = ',problem
          write(6,'(a,i10)') ' ifail = ',ifail
          write(6,*) ' this problem will be skipped'
          goto 1000
      endif
      isolver = 1
      ix0 = 1
      if (gsolver.gt.0) isolver = gsolver
      if (grtol.gt.0.0d0) rtol = grtol
1100  continue
      read(lucon,'(a80)',end=2000) directive
      if (directive(1:1).eq.'#') goto 1100
      if (directive(1:5).eq.' end ') then
        goto 2000
      else if (directive(1:11).eq.' dimension ') then
        read(directive(12:80),*) n
        call init(n, x0, xlow, xup, xsol, xthrsh, fthrsh, 
     $            nx0, nbound, nsol, rtol, ifail)
        if (ifail.ne.0) then
          write(6,*) ' error in init: problem = ',problem
          write(6,'(a,i10)') ' ifail = ',ifail
          write(6,*) ' this problem will be skipped'
          goto 1000
        endif
      else if (directive(1:8).eq.' xstart ') then
        read(directive(9:80),*) ix0
        if (ix0.gt.nx0) then
          write(6,*) ' error in definition: problem = ',problem
          write(6,'(a,i2,a,i2)') ' ix0=',ix0,' .gt. nx0=',nx0
          ix0 = 1
        endif
      else if (directive(1:8).eq.' solver ') then
        read(directive(9:80),*) isolver
        if (isolver.gt.n_solver) then
          write(6,*) ' error in definition: problem = ',problem
          write(6,'(a,i2,a,i2)') ' isolver=',isolver,' .gt. n_solver=',
     $                           n_solver
          isolver = 1
        endif
      else if (directive(1:6).eq.' rtol ') then
        read(directive(7:80),*) rtol
      else if (directive(1:8).eq.' nonlin ') then
        read(directive(9:80),*) nonlin
      else if (directive(1:5).eq.' par(') then
        ie = 6
        do i=6,80
          if (directive(i:i).eq.')') then
            ie = i-1
            goto 1150
          endif
        enddo
1150    continue
        read(directive(6:ie),*) ii
        if (ii.gt.0 .and. ii.le.max_par)
     $     read(directive(ie+3:80),*) p(ii)
      else
          write(6,*) ' error in definition: problem = ',problem
          write(6,*) ' invalid directive: ',directive(1:60)
      endif
      goto 1100
2000  continue
      if (n.lt.1 .or. n.gt.max_n) then
          write(6,*) ' error in init: problem = ',problem
          write(6,'(a,i10)') ' invalid n = ',n
          write(6,*) ' this problem will be skipped'
          goto 1000
      endif
      if (nbound.eq.0) ibound = 0
      if ( isolver.le.2 ) then
         write(lumon,'(a,a,//)') ' problem: ',problem
         write(lusol,'(a,a,//)') ' problem: ',problem
      endif
      do i=1,liopt
        iopt(i) = 0
      enddo
      do i=1,liwk
        iwk(i) = 0
      enddo
      ijac = 1
      call jac(n, n, x, rwk, p, ifail)
      if (ifail.eq.-999) ijac = 0
      do i=1,lrwk
        rwk(i) = 0.0d0
      enddo
      do i=1,n
        x(i) = x0(i,ix0)
      enddo
      if ( gxthr .ge. 0.0d0 ) then
        do i=1,n
          xscale(i) = gxthr
        enddo
      else
        do i=1,n
          xscale(i) = xthrsh(i)
        enddo
      endif
      if ( gfthr .ge. 0.0d0 ) then
        do i=1,n
          fscale(i) = gfthr
        enddo
      else
        do i=1,n
          fscale(i) = fthrsh(i)
        enddo
      endif
      if (isolver.eq.1) then
c       isolver = 1 : call nleq1
        solver = 'nleq1'
        lsolver = 5
        if (ijac.eq.1) then
          iopt(3) = 1
        else
          iopt(3) = 2
        endif
        iopt(11) = 2
        iopt(12) = lumon
        iopt(13) = 3
        iopt(14) = lumon
        iopt(15) = 1
        iopt(16) = lusol
c       iopt(19) = 1
c       iopt(20) = lumon
c       assumed nonlinearity level of problem
        iopt(31) = nonlin
c       maximum number of iterations
        iwk(31) = max_iter
        ifail = 0
        call zibsec(stime,ifail2)
        call nleq1(n,f_nleq,j_nleq,x,xscale,rtol,iopt,ifail,
     $             liwk,iwk,lrwk,rwk)
        call zibsec(etime,ifail2)
        niter = iwk(1)
        nfcn  = iwk(4)
        njac  = iwk(5)
      else if (isolver.eq.2) then
c       isolver = 2 : call nleq2
        solver = 'nleq2'
        lsolver = 5
        if (ijac.eq.1) then
          iopt(3) = 1
        else
          iopt(3) = 2
        endif
        iopt(11) = 2
        iopt(12) = lumon
        iopt(13) = 3
        iopt(14) = lumon
        iopt(15) = 1
        iopt(16) = lusol
c       iopt(19) = 1
c       iopt(20) = lumon
c       assumed nonlinearity level of problem
        iopt(31) = nonlin
c       maximum number of iterations
        iwk(31) = max_iter
c       Maximum permitted sub-condition number
c       rwk(25) = 1.0d10
        ifail = 0
        call zibsec(stime,ifail2)
        call nleq2(n,f_nleq,j_nleq,x,xscale,rtol,iopt,ifail,
     $             liwk,iwk,lrwk,rwk)
        call zibsec(etime,ifail2)
        niter = iwk(1)
        nfcn  = iwk(4)
        njac  = iwk(5)
      else if  (isolver.ge.3.and.isolver.le.6) then
c       isolver = 3,4,5,6 : call cnlsol
        if      (isolver.eq.3) then
          solver = 'nleq_err'
          lsolver = 8
        else if (isolver.eq.4) then
          solver = 'nleq_res'
          lsolver = 8
        else if (isolver.eq.5) then
          solver = 'qnerr'
          lsolver = 5
        else
          solver = 'qnres'
          lsolver = 5
        endif
        do i=1,8
          iopt(i)=0
        enddo
c       jacobian: 0=computed by numerical differentation
c                 1=supplied by user routine jac
        iopt(1) = ijac
c       maximum number of iterations
        iopt(2)=max_iter
c       2 = mildly nonlinear , 3 = highly nonlinear
c       4 = extremely nonlinear
        iopt(3)=nonlin
c       monotonicity check: 0=normal, 1=restricted
        iopt(4)=0
c       set mprerr, mprmon, mprsol
        iopt(5)=2
        iopt(6)=3
        iopt(7)=1
c       The solver to call:
c         nleq_err=1, nleq_res=2, qnerr=3, qnres=4
        iopt(9)=isolver-2
c       call nonlinear solver
        if ( isolver.eq.3 .or. isolver.eq.5 ) then
c         default scaling=0 or 1, use xscale as lower threshold=2
          iopt(8)=2
          call cnlsol(n,f_nleq,j_nleq,x,xscale,rtol,iopt,problem,info)
        else
c         no scaling=0, f(x^0) based scaling=1, use fscale=2
          iopt(8)=1
          call cnlsol(n,f_nleq,j_nleq,x,fscale,rtol,iopt,problem,info)
        endif
        ifail =   info(1)
        niter =   info(3)
        nfcn  =   info(4)
        njac  =   info(5)
      endif
      call f_nleq(n,x,f,ifail2)
      fnorm = 0.0d0
      do i=1,n
        fnorm = fnorm+f(i)*f(i)
      enddo
      fnorm = dsqrt(fnorm/dble(n))
      if (isolver.le.2) 
     1   write(lumon,'(" norm of residuum = ",d12.4,/)') fnorm
      if (ifail.eq.0) then
        nsolvd = nsolvd + 1
      else
        nfaild = nfaild + 1
      endif
      nfctot = nfctot + nfcn
      njctot = njctot + njac
c     write summary statistics
      etime = etime - stime
      if (niter.ge.0) then
        write(lusum,'(a48,4(3x,i3),3x,e10.2,3x,f7.3)') problem, ifail,
     $                                 niter, nfcn, njac, fnorm, etime
      else
        write(lusum,8900) problem, ifail, 'n/a', nfcn, njac, fnorm,
     $                    etime
8900    format(a48,3x,i3,3x,a3,2(3x,i3),3x,d10.2,3x,f7.3)
      endif
      goto 1000
9000  continue
      write(lusum,10100) nsolvd, nfaild, nfctot, njctot
      close(lucon)
      close(lusum)
      if (gsolver.le.2) then
        close(lumon)
        close(lusol)
      endif
      goto 9999
9100  write(6,*) ' +++ fatal error - could not open ',
     $           control(1:l1)//'.control'
      stop
9110  write(6,*) ' +++ fatal error - could not open ',
     $           control(1:l1)//'.monout'
      stop
9120  write(6,*) ' +++ fatal error - could not open ',
     $           control(1:l1)//'.sumout'
      stop
9130  write(6,*) ' +++ fatal error - could not open ',
     $           control(1:l1)//'.solout'
      stop
9999  continue
10100 format(//,' Number of problem solved:             ',i10,/
     1          ' Numer of problems failed:             ',i10,/,
     2          ' Total number of Function evaluations: ',i10,/,
     3          ' Total number of Jacobian evaluations: ',i10)
      end
c
      subroutine f_nleq(n,x,f,ifail)
      implicit none
c
c     function evaluation interface routine for NLEQ1
c
      integer n, ifail
      double precision x(n), f(n)
c
      integer max_n, max_con, max_par
c      parameter ( max_n   = 100,
      parameter ( max_n   = 200,
     $            max_con =  10, 
     $            max_par = 100)
c
      integer np, ncon, nsave
      double precision par(max_par), conxl(max_n,max_con),
     $                 conxu(max_n,max_con)
      common /newlab/ np, par, ncon, conxl, conxu, nsave
c
      integer nfev, njev
      common /refnum/ nfev, njev
c
      call fcn(n, x, f, par, ifail)
      nfev = nfev + 1
      return
      end
c
      subroutine j_nleq(n,ldfjac,x,dfdx,ifail) 
      implicit none
c
c     Jacobian evaluation interface routine for NLEQ1
c
      integer ldfjac, n, ifail
      double precision x(n), dfdx(ldfjac,n)
c
      integer max_n, max_con, max_par
c      parameter ( max_n   = 100,
      parameter ( max_n   = 200,
     $            max_con =  10, 
     $            max_par = 100)
c
      integer np, ncon, nsave
      double precision par(max_par), conxl(max_n,max_con),
     $                 conxu(max_n,max_con)
      common /newlab/ np, par, ncon, conxl, conxu, nsave
c
      integer nfev, njev
      common /refnum/ nfev, njev
c
      call jac(ldfjac, n, x, dfdx, par, ifail)
      njev = njev + 1
      return
      end
