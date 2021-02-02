FUNCTION fBed(x0,y0)
  USE types
  IMPLICIT NONE

  REAL(kind=dp) :: x0, y0, fBed, a, b

! droite descendante
  a = -0.001_dp
  b = -100.0_dp

  fBed = a * x0 + b

END FUNCTION fBed
