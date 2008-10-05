!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-11-09  Time: 22:01:17

PROGRAM Test_uobyqa
!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6,8.

USE Powell_Optimize
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp)  :: rhobeg, rhoend, x(10)
INTEGER    :: i, iprint, maxfun, n

iprint=2
maxfun=5000
rhoend=1.0D-8
DO  n=2,8,2
  DO  i=1,n
    x(i)=REAL(i)/REAL(n+1)
  END DO
  rhobeg=0.2D0*x(1)
  WRITE(*, 20) n
  20 FORMAT (//t6, '******************'/ t6, 'Results with N =', i2 / t6, '******************')
  CALL uobyqa (n, x, rhobeg, rhoend, iprint, maxfun)
END DO
STOP
END PROGRAM Test_uobyqa

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calfun.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
SUBROUTINE calfun(n, x, f)

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f

! COMMON a(40,20), b(40,20), e(41)
REAL (dp)  :: sum, y(10,10)
INTEGER    :: i, iw, j, np

DO  j = 1, n
  y(1,j) = 1.0D0
  y(2,j) = 2.0D0 * x(j) - 1.0D0
END DO
DO  i = 2, n
  DO  j = 1, n
    y(i+1,j) = 2.0D0 * y(2,j) * y(i,j) - y(i-1,j)
  END DO
END DO
f = 0.0D0
np = n + 1
iw = 1
DO  i = 1, np
  sum = 0.0D0
  DO  j = 1, n
    sum = sum + y(i,j)
  END DO
  sum = sum / REAL(n)
  IF (iw > 0) sum = sum + 1.0_dp / REAL(i*i - 2*i)
  iw = -iw
  f = f + sum * sum
END DO
RETURN
END SUBROUTINE calfun
