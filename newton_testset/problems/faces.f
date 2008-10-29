c-----------------------------------------------------------------------
c
c   Switching subroutines to problem routines.
c
c   Coded by         L. Weimann, 
c                    Zuse Institute Berlin (ZIB),
c                    Dep. Numerical Analysis and Modelling
c   File             faces.f
c   Version          1.0.1
c   Latest Change    14th December 2005
c   Code             Fortran 77 with common extensions
c                    Double precision
c
c-----------------------------------------------------------------------
      subroutine init(n, x0, xlow, xup, xsol, xthrsh, fthrsh,
     $           nx0, nbound, nsol, rtol, ifail)
      implicit none
      integer max_n, max_con, max_par
c      parameter ( max_n   = 100,
      parameter ( max_n   = 200,
     $            max_con =  10, 
     $            max_par = 100)
      integer n, nx0, nbound, nsol, ifail
      double precision x0(max_n,*),xlow(max_n,*), xup(max_n,*),
     $                 xsol(max_n,*), xthrsh(max_n), fthrsh(max_n), rtol
c     **********
      character*80 task
      integer ldfjac
      double precision x(max_n), f(1), dfdx(1,1)
c
      character*32 modnam
      common /modul/ modnam
c
      integer np, ncon, nsave
      double precision p(max_par), conxl(max_n,max_con),
     $                 conxu(max_n,max_con)
      common /newlab/ np, p, ncon, conxl, conxu, nsave
c
      ldfjac = 1
      nx0    = 1
      nbound = 0
      ncon   = 0
      nsol   = 0
c     rtol = 1.0d-10
      task = 'XS'
      if (modnam .eq. 'InfReflux') then
        call InfReflux(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Order10to11') then
        call Order10to11(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ReinforcedBeam') then
        call ReinforcedBeam(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CompFacRK') then
        call CompFacRK(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CREquiConv4') then
        call CREquiConv4(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PartMethaneOx') then
        call PartMethaneOx(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CRSteadyState5') then
        call CRSteadyState5(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'HodgkinHuxley1') then
        call HodgkinHuxley1(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Cutlip') then
        call Cutlip(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Wallis1685') then
        call Wallis1685(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsfifj') then
        call dsfifj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'delhfj') then
        call delhfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'disofj') then
        call disofj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsfdfj') then
        call dsfdfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dierfj') then
        call dierfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsulfj') then
        call dsulfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dgupfj') then
        call dgupfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dmetfj') then
        call dmetfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'destfj') then
        call destfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dhhdfj') then
        call dhhdfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dficfj') then
        call dficfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dcprfj') then
        call dcprfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dfdcfj') then
        call dfdcfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dcpffj') then
        call dcpffj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail) 
      else if (modnam .eq. 'BrownAlmLin') then
        call BrownAlmLin(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'BroydenBanded') then
        call BroydenBanded(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'BroydenTridiag') then
        call BroydenTridiag(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Chebyquad') then
        call Chebyquad(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DiscreteBoundary') then
        call DiscreteBoundary(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DiscreteIntegral') then
        call DiscreteIntegral(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'HelicalValley') then
        call HelicalValley(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PowellBadlyScaled') then
        call PowellBadlyScaled(n, x, f, dfdx, ldfjac, task, p, np,max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PowellSingular') then
        call PowellSingular(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Rosenbrock') then
        call Rosenbrock(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Trigonometric') then
        call Trigonometric(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'VariablyDim') then
        call VariablyDim(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Watson') then
        call Watson(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Wood') then
        call Wood(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui2') then
        call ChemEqui2(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'LargeExp') then
        call LargeExp(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ExpSin') then
        call ExpSin(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui3') then
        call ChemEqui3(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui1') then
        call ChemEqui1(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DistillColumn') then
        call DistillColumn(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'sst') then
        call sst(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'InterSEllHy') then
        call InterSEllHy(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'EsterificRea') then
        call EsterificRea(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Gupta') then
        call Gupta(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DewPoint') then
        call DewPoint(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PartOxiOfMet') then
        call PartOxiOfMet(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'SulpDiox') then
        call SulpDiox(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else
        write(6,*) ' +++ init: invalid module name: ', modnam
        stop ' +++ init stopped due to severe error +++'
      endif
      nsave = n
c      
      return
      end
c
      subroutine fcn(n, x, f, p, ifail)
      implicit none
      integer n, ifail
      double precision x(n), f(n), p(*)
c
      integer max_n, max_con, max_par
c      parameter ( max_n   = 100,
      parameter ( max_n   = 200,
     $            max_con =  10, 
     $            max_par = 100)
c
      character*80 task
      integer ldfjac, nx0, nbound, nsol
      double precision dfdx, x0, xlow, xup, xsol, 
     $                 xthrsh, fthrsh, rtol
c
      character*32 modnam
      common /modul/ modnam
c
      integer np, ncon, nsave
      double precision par(max_par), conxl(max_n,max_con),
     $                 conxu(max_n,max_con)
      common /newlab/ np, par, ncon, conxl, conxu, nsave
c
      task = 'F'
      ldfjac = 1
      if (modnam .eq. 'InfReflux') then
        call InfReflux(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Order10to11') then
        call Order10to11(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ReinforcedBeam') then
        call ReinforcedBeam(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CompFacRK') then
        call CompFacRK(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CREquiConv4') then
        call CREquiConv4(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PartMethaneOx') then
        call PartMethaneOx(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CRSteadyState5') then
        call CRSteadyState5(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'HodgkinHuxley1') then
        call HodgkinHuxley1(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Cutlip') then
        call Cutlip(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Wallis1685') then
        call Wallis1685(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsfifj') then
        call dsfifj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'delhfj') then
        call delhfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'disofj') then
        call disofj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsfdfj') then
        call dsfdfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dierfj') then
        call dierfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsulfj') then
        call dsulfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dgupfj') then
        call dgupfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dmetfj') then
        call dmetfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'destfj') then
        call destfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dhhdfj') then
        call dhhdfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dficfj') then
        call dficfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dcprfj') then
        call dcprfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dfdcfj') then
        call dfdcfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dcpffj') then
        call dcpffj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail) 
      else if (modnam .eq. 'BrownAlmLin') then
        call BrownAlmLin(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'BroydenBanded') then
        call BroydenBanded(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'BroydenTridiag') then
        call BroydenTridiag(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Chebyquad') then
        call Chebyquad(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DiscreteBoundary') then
        call DiscreteBoundary(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DiscreteIntegral') then
        call DiscreteIntegral(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'HelicalValley') then
        call HelicalValley(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PowellBadlyScaled') then
        call PowellBadlyScaled(n, x, f, dfdx, ldfjac, task, p, np,max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PowellSingular') then
        call PowellSingular(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Rosenbrock') then
        call Rosenbrock(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Trigonometric') then
        call Trigonometric(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'VariablyDim') then
        call VariablyDim(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Watson') then
        call Watson(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Wood') then
        call Wood(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui2') then
        call ChemEqui2(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'LargeExp') then
        call LargeExp(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ExpSin') then
        call ExpSin(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui3') then
        call ChemEqui3(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui1') then
        call ChemEqui1(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DistillColumn') then
        call DistillColumn(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'sst') then
        call sst(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'InterSEllHy') then
        call InterSEllHy(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'EsterificRea') then
        call EsterificRea(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Gupta') then
        call Gupta(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DewPoint') then
        call DewPoint(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PartOxiOfMet') then
        call PartOxiOfMet(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'SulpDiox') then
        call SulpDiox(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else
        write(6,*) ' +++ fcn: invalid module name: ', modnam
        stop ' +++ fcn stopped due to severe error +++'
      endif
      if (task(1:5).eq.'ERROR') ifail = -1000
      return
      end
c
      subroutine jac(ldfjac, n, x, dfdx, p, ifail)
      implicit none
      integer ldfjac, n, ifail
      double precision x(n), dfdx(ldfjac,n), p(*)
c
      integer max_n, max_con, max_par
c      parameter ( max_n   = 100,
      parameter ( max_n   = 200,
     $            max_con =  10, 
     $            max_par = 100)
c
      character*80 task
      integer nx0, nbound, nsol
      double precision f, x0, xlow, xup, xsol, 
     $                 xthrsh, fthrsh, rtol
c
      character*32 modnam
      common /modul/ modnam
c
      integer np, ncon, nsave
      double precision par(max_par), conxl(max_n,max_con),
     $                 conxu(max_n,max_con)
      common /newlab/ np, par, ncon, conxl, conxu, nsave
c
      task = 'J'
      if (modnam .eq. 'InfReflux') then
        call InfReflux(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Order10to11') then
        call Order10to11(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ReinforcedBeam') then
        call ReinforcedBeam(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CompFacRK') then
        call CompFacRK(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CREquiConv4') then
        call CREquiConv4(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PartMethaneOx') then
        call PartMethaneOx(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'CRSteadyState5') then
        call CRSteadyState5(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'HodgkinHuxley1') then
        call HodgkinHuxley1(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Cutlip') then
        call Cutlip(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Wallis1685') then
        call Wallis1685(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsfifj') then
        call dsfifj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'delhfj') then
        call delhfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'disofj') then
        call disofj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsfdfj') then
        call dsfdfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dierfj') then
        call dierfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dsulfj') then
        call dsulfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dgupfj') then
        call dgupfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dmetfj') then
        call dmetfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'destfj') then
        call destfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dhhdfj') then
        call dhhdfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dficfj') then
        call dficfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dcprfj') then
        call dcprfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dfdcfj') then
        call dfdcfj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'dcpffj') then
        call dcpffj(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail) 
      else if (modnam .eq. 'BrownAlmLin') then
        call BrownAlmLin(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'BroydenBanded') then
        call BroydenBanded(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'BroydenTridiag') then
        call BroydenTridiag(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Chebyquad') then
        call Chebyquad(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DiscreteBoundary') then
        call DiscreteBoundary(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DiscreteIntegral') then
        call DiscreteIntegral(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'HelicalValley') then
        call HelicalValley(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PowellBadlyScaled') then
        call PowellBadlyScaled(n, x, f, dfdx, ldfjac, task, p, np,max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PowellSingular') then
        call PowellSingular(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Rosenbrock') then
        call Rosenbrock(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Trigonometric') then
        call Trigonometric(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'VariablyDim') then
        call VariablyDim(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Watson') then
        call Watson(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Wood') then
        call Wood(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui2') then
        call ChemEqui2(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'LargeExp') then
        call LargeExp(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ExpSin') then
        call ExpSin(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui3') then
        call ChemEqui3(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'ChemEqui1') then
        call ChemEqui1(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DistillColumn') then
        call DistillColumn(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'sst') then
        call sst(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'InterSEllHy') then
        call InterSEllHy(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'EsterificRea') then
        call EsterificRea(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'Gupta') then
        call Gupta(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'DewPoint') then
        call DewPoint(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'PartOxiOfMet') then
        call PartOxiOfMet(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else if (modnam .eq. 'SulpDiox') then
        call SulpDiox(n, x, f, dfdx, ldfjac, task, p, np, max_n, 
     $           x0, xlow, xup, conxl, conxu, xsol, xthrsh, fthrsh,
     $           nx0, nbound, ncon, nsol, rtol, ifail)
      else
        write(6,*) ' +++ jac: invalid module name: ', modnam
        stop ' +++ jac stopped due to severe error +++'
      endif
      if (task(1:5).eq.'ERROR') ifail = -1000
      return
      end
