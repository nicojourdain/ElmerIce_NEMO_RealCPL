!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    to be compiled with the associated geometry function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION BedPlusEpsilon ( Model, nodenumber, x) RESULT(yb)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  INTEGER :: nodenumber
  REAL(KIND=dp) :: x, y, yb, fbed          

  x = Model % Nodes % x ( nodenumber )
  y = Model % Nodes % y ( nodenumber )

  yb = fBed( x , y ) + 10

END FUNCTION BedPlusEpsilon


FUNCTION Bed ( Model, nodenumber, x) RESULT(yb)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  INTEGER :: nodenumber
  REAL(KIND=dp) :: x, y, yb, fbed          

  x = Model % Nodes % x ( nodenumber )
  y = Model % Nodes % y ( nodenumber )

  yb = fBed( x , y )

END FUNCTION Bed

