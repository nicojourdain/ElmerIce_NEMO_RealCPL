SUBROUTINE GetLateralFrictionForce( Model,Solver,dt,TransientSimulation )

!DEC$ATTRIBUTES DLLEXPORT :: GetHydrostaticLoads
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Computes nodal value of the force resulting from the 
!  lateral friction law (USF_LateralFriction.f90)
!  work only in 2D (no sense in 3D)
!  Compute also the nodal area of each nodes belonging in the ice-shelf
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
  TYPE(ValueList_t), POINTER :: BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable, FlowVariable
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Nodes_t) :: Nodes
  TYPE(GaussIntegrationPoints_t) :: IP

  LOGICAL :: AllocationsDone = .FALSE., GotIt, stat

  INTEGER :: i, j, n, m, t, p, Nn, istat, DIM
  INTEGER, POINTER :: Permutation(:), FlowPerm(:), NodeIndexes(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), FlowValues(:)
  REAL(KIND=dp) :: Norm, Velo(3), Kspring, FVector(3), s
  REAL(KIND=dp) :: detJ, mm, Xg, density                                    

  REAL(KIND=dp), ALLOCATABLE :: Kl(:), Basis(:), dBasisdx(:,:), NodalDensity(:),&
                                ddBasisddx(:,:,:), LocalVelo(:,:), x(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FlowSolverName
       

  SAVE AllocationsDone, DIM, SolverName, FlowSolverName, mm 
  SAVE Nodes, Basis, dBasisdx, ddBasisddx
  SAVE LocalVelo, Kl, x, NodalDensity
  
  !------------------------------------------------------------------------------

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------

  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN

    DIM = CoordinateSystemDimension()
    WRITE(SolverName, '(A)') 'GetLateralFrictionForce'
    n = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
    m = Model % Mesh % NumberOfNodes
    IF (AllocationsDone)  &
          DEALLOCATE(Kl, Basis, dBasisdx, ddBasisddx, LocalVelo, x, &
                     NodalDensity)

    ALLOCATE(Kl(n), Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), &
             LocalVelo(3,n), x(n), NodalDensity(n),  STAT=istat )
         
    IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
    END IF

    AllocationsDone = .TRUE.
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )


  END IF  

  !--------------------------------------------
  ! Calculate the friction force for each elements
  !--------------------------------------------
      Element => GetActiveElement(1)
      BodyForce => GetBodyForce()  
      FlowSolverName = GetString( BodyForce, 'Flow Solver Name', GotIt )    
      IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'

      FlowSolverName = 'Flow Solution'
      FlowVariable => VariableGet( Model % Variables, FlowSolverName )
      IF ( ASSOCIATED( FlowVariable ) ) THEN
         FlowPerm    => FlowVariable % Perm
         FlowValues  => FlowVariable % Values
      ELSE
         CALL Info('LateralFrictionForce', &
                    & 'No variable for velocity associated.', Level=4)
      END IF



      Xg = GetConstReal ( Model % Constants, 'Position GrdL', GotIt )


  VariableValues = 0.0_dp
  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes

    BodyForce => GetBodyForce()  
    Material => GetMaterial()  

    mm = GetConstReal( BodyForce, 'Lateral Friction Exponent', GotIt )
    IF (.Not.GotIt ) CALL FATAL('LateralFrictionForce', &
                'Lateral Friction Exponent must be defined')
!    
! Get the nodale velocities
    LocalVelo = 0.0_dp
    DO i=1, DIM
      LocalVelo(i,1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + i)
    END DO
!
! Get the nodale value of K

    Kl(1:n) = GetReal( BodyForce, &
                   'Lateral Friction Coefficient', GotIt )
    IF (.Not.GotIt ) CALL FATAL('LateralFrictionForce', &
                'Lateral Friction Coefficient must be defined')
!
! Get the nodale density
    NodalDensity(1:n) = GetReal(Material , &
                   'Density', GotIt )
    IF (.Not.GotIt ) CALL FATAL('LateralFrictionForce', &
                'Density must be defined')
  
! Force computed only for elements where all x > xg                   
    x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
    If (MINVAL(x) < Xg) Kl = 0.0 
!
! Integration
! 
    CALL GetElementNodes( Nodes )
    IP = GaussPoints( Element )
    DO p = 1, IP % n

      stat = ElementInfo( Element, Nodes, IP % U(p), IP % V(p), &
      IP % W(p), detJ, Basis, dBasisdx, ddBasisddx, .FALSE.)          
      s = detJ * IP % S(p)                  

!
!  Value of K * |u|^{mm-1} at integartion point  
!
      DO i=1, DIM
         Velo(i) = SUM(LocalVelo(i,1:n)*Basis(1:n))
      END DO
      Kspring = SUM(Kl(1:n)*Basis(1:n)) 
      density = SUM(NodalDensity(1:n)*Basis(1:n))
      Kspring = density * Kspring
      Norm = 0.0
      DO i=1, DIM
        norm = norm + Velo(i)**2.0
      END DO
      norm = SQRT(norm)
      IF ((norm > 1.0e-6).AND.(mm/=1.0))  Kspring = Kspring * norm**(mm-1.0) 
!
! Compute nodal Fx, nodal Fy, and nodal area   
!
      FVector(1:DIM) = - Kspring * Velo(1:DIM)

      DO i = 1, n
      Nn = Permutation(Element % NodeIndexes(i))
        DO j = 1, DIM
           VariableValues(3*(Nn-1)+j) = VariableValues(3*(Nn-1)+j) + FVector(j) * s * Basis(i)
      END DO
          If (Kspring > 0.0) VariableValues(3*Nn) = VariableValues(3*Nn) + s * Basis(i)
      END DO

    END DO

  END DO

!------------------------------------------------------------------------------
END SUBROUTINE GetLateralFrictionForce
!------------------------------------------------------------------------------



