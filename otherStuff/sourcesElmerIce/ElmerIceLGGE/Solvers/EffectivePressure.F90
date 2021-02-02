!!! Calculate the effective pressure N = pice - pwater and export a variable
!
!!! This version consider that pwater = rhow.g.(zsea-zb)
!   OUTPUT Variables:
!     N
!   
!   INPUT Variable:
!     H
!     bedrock (optionnal)
!
! PARAMETERS:
!   Constants: 
!     gravity
!   BodyForce :
!     Water Pressure = Variable Zb
!       Real Procedure "..../USF_WaterPressure" "WaterPressure"
!   Material:
!      SSA Mean Density
!     
!------------------------------------------------------------------------------
SUBROUTINE EffectivePressure( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE MeshUtils
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Variable_t),POINTER :: ZbVar=>NULL(),ZsVar=>NULL()
  TYPE(Variable_t),POINTER :: NeffVar=>NULL()
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t),POINTER :: Material, BodyForce

  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: Density
  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: Pwater
  REAL(KIND=dp) :: rhoi,g,pw
  REAL(KIND=dp) :: zb,zs,Neff


  INTEGER, POINTER,SAVE :: NodeIndexes(:)
  INTEGER :: GroundedNode
  INTEGER :: t,i,n

  LOGICAL,SAVE :: Initialized = .FALSE.
  LOGICAL :: Found

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='EffectivePressure'

!------------------------------------------------------------------------------
  Mesh => Model % Mesh

!!! Where the variable will be stored
  NeffVar => VariableGet (Model % Mesh % Variables, 'Effective Pressure')
  IF (.NOT.ASSOCIATED(NeffVar)) THEN
     Message='Effective Pressure Variable not found'
     CALL FATAL(SolverName,Message)
  END IF 

!!! get required variables Zb,Zs
  zbVar => VariableGet( Model % Mesh % Variables, 'Zb')
  IF (.NOT.ASSOCIATED(zbVar)) THEN
     Message='Zb  not found'
     CALL FATAL(SolverName,Message)
  END IF
  zsVar => VariableGet( Model % Mesh % Variables, 'Zs')
  IF (.NOT.ASSOCIATED(zsVar)) THEN
     Message='Zs  not found'
     CALL FATAL(SolverName,Message)
  END IF

!!! Do some initialisation/allocation
  IF ((.NOT.Initialized).OR.Mesh%Changed) THEN

    IF (Initialized) deallocate(Density,Pwater)
    
    N=Model % MaxElementNodes
    allocate(Density(N),Pwater(N))

    Initialized = .TRUE.
  END IF
!!

 g = GetCReal ( Model % Constants, 'Gravity', Found )
 If (.NOT.Found) THEN
    WRITE(Message,'(A)') 'Constant >Gravity< not found. &
       &Setting to -9.746289e15'
    CALL INFO(SolverName, Message, level=3)
    g = -9.746289e15
 End if


 Do t=1,Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes
! Get Body Forces and Material properties
    Material => GetMaterial(Element)
    BodyForce => GetBodyForce(Element)

    Density(1:n) = ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,Found)
    IF (.NOT.Found) THEN
       Message='<SSA Mean Density> not found in material'
       CALL FATAL(SolverName, Message)
    END IF
  
    Pwater (1:n) = GetReal(BodyForce, 'Water Pressure', Found, Element)
    IF (.NOT.Found) THEN
       Message='<Water Pressure> not found in material'
       CALL FATAL(SolverName, Message)
    END IF
    

    Do i=1,n 

       rhoi=Density(i)
       pw= Pwater(i)
       zs=ZsVar%Values(ZsVar%Perm(NodeIndexes(i)))
       zb=ZbVar%Values(ZbVar%Perm(NodeIndexes(i)))
       
       Neff=rhoi*abs(g)*(zs-zb)- pw

       NeffVar%Values(NeffVar%Perm(NodeIndexes(i)))=Neff
    End do
 End do
!------------------------------------------------------------------------------
END SUBROUTINE EffectivePressure
!------------------------------------------------------------------------------
