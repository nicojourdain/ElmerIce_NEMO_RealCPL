!SeaDamagePressure.f90
!Last modif : 13 November 2014
! 
! Solver Variable : Sea Damage Pressure
!
!The mesh must be vertically extruded, i.e. the nodes must be vertically aligned
!
!  Serial/Parallel   and 2D/3D
!
! The solver has to be executed in the main body (1)
!
! Requires some inputs in the sif file.
!    In the Solver Section:
!     - Variable = "PSeaD"
!     - Opening Damage = Real ; Critical value of Damage beyond which crevasses opened, i.e. beyond which sea pressure has to be taken into account for damage 
!     - Active Coordinate = integer ; the extrusion direction (2 in 2d and 3 in 3d in our case because it has to be VERTICALY extruded
!
!    In the Material Section:
!     - Sea Level = Real
!
!    In the Constants Section:
!     - Water Density = Real 
!
!
! Required Variables:
!     - GroundedMask
!     - Damage
!   
! WARNING:
!   The vertical direction is necessary associated to the maximum value of CoordinateSystemDimension:
!      - 2 in 2D	 
!      - 3 in 3D
!
!******************************************************************************
SUBROUTINE SeaDamagePressure( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t),TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
! 
  TYPE(ValueList_t), POINTER :: SolverParams, Material, ParentMaterial, BC, Constants, BodyForce !pointer required to access information in the .sif

!!!!  Usefull variables for element and basic functions
  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t) :: ElementNodes
  INTEGER, POINTER :: NodeIndexes(:)

!!!!! Elmer variables
  TYPE(Variable_t), POINTER :: GMSol,DamageVariable,PseaD
  REAL(KIND=dp), POINTER ::  GM(:),DamageValues(:),PseaDValues(:)
  INTEGER, POINTER :: GMPerm(:),DamagePerm(:),PseaDPerm(:)

!!!! Variables for DetectExtrudedStructure
  TYPE(Solver_t), POINTER :: PSolver  
  TYPE(Variable_t), POINTER :: Var
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),Upointer(:),DownPointer(:)

!!!! Other variables
  real(kind=dp) :: g,rho_w,lw
  real(kind=dp) :: alti
  real(kind=dp) :: OpD

  integer, allocatable :: BotNodes(:),TMP(:)
  integer :: NMAX,k,n,i,C,l
  integer :: DIM, Up

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  Logical, allocatable :: VisitedNodes(:)  
  Logical ::  Firsttime=.true.
  Logical ::  Found


  save Firsttime,DIM
  save ElementNodes
  save SolverName
  SAVE BotNodes
  save OpD,rho_w,g
  save C
  SAVE TopPointer,BotPointer,Upointer,DownPointer

  !!! Where the variable Sea Damage Pressure will be stored
  PSeaD => Solver % Variable
  PSeaDValues => PSeaD % Values
  PSeaDPerm => PSeaD % Perm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! INITIALIZATION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  If (Firsttime) then
      WRITE(SolverName, '(A)') 'SeaDamagePressure'
      DIM = CoordinateSystemDimension()

      !!!! 
      PSolver => Solver
      CALL DetectExtrudedStructure(Solver % Mesh , PSolver, Var, &
           TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
           UpNodePointer = Upointer , DownNodePointer = DownPointer)

      NMAX=Solver % Mesh % NumberOfNodes
      
      allocate(VisitedNodes(NMAX),TMP(NMAX))

!------------------------------------------------------------------
! Get usefull constants
!------------------------------------------------------------------
      !!! Get Water Density
      Constants => GetConstants()

      rho_w =  GetConstReal( Constants,'Water Density', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >WaterDensity< not found in section >Constants<')
               CALL WARN(SolverName,'Taking default value >1.03225e-18<')
               rho_w=1.03225e-18_dp
          END IF

      !!! Get the gravity


      IF (DIM==2) THEN   
        g = -GetConstReal( Model % BodyForces(1) % Values, 'Flow BodyForce 2',found)
      ELSE
        g = -GetConstReal( Model % BodyForces(1) % Values, 'Flow BodyForce 3',found)
      END IF
      
      !!! Check that g is positive
      IF (g < 0) THEN
        WRITE (Message, '(A,A,A)') 'Wrong sign for gravity in section >Body Force<'
        CALL FATAL (SolverName, Message)
      END IF
      !!! Get Solver Variables
      SolverParams => GetSolverParams()

      OpD =  GetConstReal( SolverParams,'Opening Damage', Found)
            IF(.NOT.Found) THEN
               WRITE (Message, '(A,A,A)') 'Keyword >Opening Damage< not found in section >Solver<'
               CALL FATAL (SolverName, Message)
            END IF
!------------------------------------------------------------------
! Get the basale nodes and store them in BotNodes()
!------------------------------------------------------------------ 
 
     VisitedNodes=.false.
     C=0.0
DoBoundaryElements: DO k=1,Solver % Mesh % NumberOfBoundaryElements
              Element => GetBoundaryElement(k)
              n = GetElementNOFNodes()
              NodeIndexes => Element % NodeIndexes
  
              DoNodes: Do i=1,n
                  If ( VisitedNodes(NodeIndexes(i)) ) CYCLE
                       VisitedNodes(NodeIndexes(i))=.true. 
                  If( BotPointer(NodeIndexes(i)) /= NodeIndexes(i) ) CYCLE
                    C=C+1
                    TMP(C)=NodeIndexes(i)
                  End Do DoNodes

  END DO DoBoundaryElements 

  allocate (BotNodes(C))
      BotNodes(1:C)=TMP(1:C) 
  deallocate (TMP, VisitedNodes)
 !!! End of First visit
      Firsttime=.false.
  Endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! END OF INITIALIZATION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!------------------------------------------------------------------
! Get the Sea Level
!------------------------------------------------------------------

!  Material => GetMaterial()
  ParentMaterial => Model % Materials(1) % Values
  
  lw =  GetConstReal( ParentMaterial,'Sea level', Found )
      IF(.NOT.Found) THEN        
         CALL WARN(SolverName,'Keyword >Sea level< not found in section >Material<')
         CALL WARN(SolverName,'Taking default value >0.0<')
         lw=0.0_dp
      END IF 
      
!------------------------------------------------------------------
!Get the variables needed by the solver
!------------------------------------------------------------------

!!! Get the grounded mask

  GMSol => VariableGet( Solver % Mesh % Variables, 'GroundedMask' )
      IF (ASSOCIATED(GMSol)) THEN
         GM => GMSol % Values
         GMPerm => GMSol % Perm
      ELSE
         WRITE(Message,'(A,A,A)') 'No variable >GroundedMask< found'
         CALL FATAL(SolverName,Message)
      END IF
        
!!! Get the damage field
  DamageVariable => VariableGet( Solver % Mesh % Variables, 'Damage' )
      IF (ASSOCIATED(DamageVariable)) THEN
         DamageValues => DamageVariable % Values
         DamagePerm => DamageVariable % Perm
      ELSE
         WRITE(Message,'(A,A,A)') 'No variable >Damage< found'
         CALL FATAL(SolverName,Message)
      END IF

!-----------------------------------------------------------------
!Apply the water pressure in the submarine crevasses
!-----------------------------------------------------------------

!!!Go through the bottom nodes and check if damage is higher than Opening damage. If yes, apply water pressure and look at the above node. Repeat the procedure until damage is less than opening Damage or top node is reached.  
  PSeaDValues =0.0_dp

  DO l=1,C 
      IF (GM(GMPerm(BotNodes(l))) > 0.5) THEN
            PSeaDValues(PSeaDPerm(BotNodes(l)))=0.0_dp
      ELSE
         IF (DIM .eq. 2) THEN    !!!!! In this IF loop we apply PSeaD on every nodes located under the shelf, even if D < OpD
            PSeaDValues(PSeaDPerm(BotNodes(l)))= MAX (0.0_dp, rho_w * g * (lw - Model % Nodes % y (BotNodes (l))))
         ELSE
            PSeaDValues(PSeaDPerm(BotNodes(l)))= MAX (0.0_dp, rho_w * g * (lw - Model % Nodes % z (BotNodes (l))))
         END IF
         Up = Upointer (BotNodes(l))
         Do While (DamageValues(DamagePerm(Up)) >= OpD .AND. Up /= TopPointer(BotNodes(l))) 
            IF (DIM .eq. 2) THEN
               alti = Model % Nodes % y (Up)
            ELSE
               alti = Model % Nodes % z (Up)
            END IF
            PSeaDValues(PSeaDPerm(Up))= MAX (0.0_dp, rho_w * g * (lw - alti))
            Up = Upointer (Up)
         End Do
      ENDIF 
  End Do  
!------------------------------------------------------------------------------
END SUBROUTINE SeaDamagePressure
!------------------------------------------------------------------------------


