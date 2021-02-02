SUBROUTINE GetArea( Model,Solver,dt,TransientSimulation )

!DEC$ATTRIBUTES DLLEXPORT :: GetHydrostaticLoads
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Modification of GetHydrostaticLoads to get the area of each element
!                               (Basile 30/05/2012)
!  
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
 
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: PointerToVariable
  TYPE(Nodes_t) :: Nodes
  TYPE(GaussIntegrationPoints_t) :: IP

  LOGICAL :: AllocationsDone = .FALSE., stat

  INTEGER :: i, n, m, t, p, Nn, istat
  INTEGER, POINTER :: Permutation(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp) :: s
  REAL(KIND=dp) :: detJ                                    

  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), &
                                ddBasisddx(:,:,:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
       

  SAVE AllocationsDone, SolverName
  SAVE Nodes, Basis, dBasisdx, ddBasisddx
  
  !------------------------------------------------------------------------------

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------

  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN

    WRITE(SolverName, '(A)') 'GetArea'
    n = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
    m = Model % Mesh % NumberOfNodes
    IF (AllocationsDone) DEALLOCATE( Basis, dBasisdx, ddBasisddx)

    ALLOCATE( Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), STAT=istat )
         
    IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
    END IF

    AllocationsDone = .TRUE.
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )

  END IF  

  !--------------------------------------------
  ! Apply a unity force at each element
  !--------------------------------------------

  VariableValues = 0.0_dp

  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
    n = GetElementNOFNodes()

    !
    ! Integration of element surface
    ! 

    CALL GetElementNodes( Nodes )
    IP = GaussPoints( Element )
    DO p = 1, IP % n

      stat = ElementInfo( Element, Nodes, IP % U(p), IP % V(p), &
      IP % W(p), detJ, Basis, dBasisdx, ddBasisddx, .FALSE.)        
      s = detJ * IP % S(p)                            
 

      DO i = 1, n
         Nn = Permutation(Element % NodeIndexes(i))
         VariableValues(Nn) = VariableValues(Nn) + s * Basis(i) 
      END DO
      
    END DO

  END DO
    
  IF ( ParEnv % PEs.GT.1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues)

!------------------------------------------------------------------------------
END SUBROUTINE GetArea
!------------------------------------------------------------------------------



