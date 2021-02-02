!Computation of a coupling parameter dependant on the parameters that lead the
!wall melting and ice creeping in the channels.
!-----------------------------------------------------------------------------
FUNCTION RothPhi(Model,nodenumber,x) RESULT(Coupling)

  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription

  IMPLICIT NONE
  !-----------------------------------------------
  !    External variables
  !-----------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(KIND=dp) :: Coupling, x

  !----------------------------------------------- 
  !    Local variables
  !-----------------------------------------------

  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: WaterLoadGrad1, WaterLoadGrad2, &
       CELWaterLoad, TimeVar
  TYPE(ValueList_t), POINTER :: Constant, Material, BodyForce


  REAL(KIND=dp), POINTER :: WLG1Values(:), WLG2Values(:), &
       CELLoadValues(:), Hwrk(:,:,:)
  REAL(KIND=dp), ALLOCATABLE :: GlenFactor(:), GlenExp(:), &
       g(:,:), gravity(:), OldCoupling(:,:), &
       IceDensity(:), WaterDensity(:), Hmax(:), &
       Transmitivity(:,:,:)
  REAL(KIND=dp) :: LFusion, TimeStep, EffectivePress, PressRatio, &
       GradFactor, CELLoad, tps, &
       Temperature, R, Tlimit, A1, A2, Q1, Q2


  INTEGER, POINTER :: WLG1Perm(:), WLG2Perm(:), &
       CELLoadPerm(:), NodeIndexes(:)
  INTEGER :: N, istat, bf_id, i, j, M, &
       material_id, DIM

  CHARACTER(LEN=20) :: SolverName, ViscosityFlag

  LOGICAL :: AllocationsDone=.FALSE., FirstTime=.TRUE., &
       Found

  SAVE GlenFactor,      &
       GlenExp,         &
       g,               &
       gravity,         &
       AllocationsDone, &
       OldCoupling,     &
       IceDensity,      &
       WaterDensity,    &
       Hmax,            &
       Transmitivity,   &
       FirstTime


  IF (FirstTime)THEN
     FirstTime = .FALSE.
     IF ( AllocationsDone ) THEN
        DEALLOCATE(GlenFactor, &
             GlenExp,          &
             g,                &
             gravity,          &
             OldCoupling,      &
             IceDensity,       &
             WaterDensity,     &
             Hmax,             &
             Transmitivity)

     END IF

     M = Model % Nodes % NumberOfNodes
     N = Model % MaxElementNodes
     ALLOCATE(GlenFactor( N ),  &
          GlenExp( N ),         &
          g( 3,N ),             &
          gravity( N ),         &
          OldCoupling( M, 3 ),  &
          IceDensity( N ),      &
          WaterDensity( N ),    &
          Hmax( N ),            &
          Transmitivity(3,3,N), &
          STAT=istat)

     IF ( istat /= 0 ) THEN
        CALL FATAL(  'USF_TransferCoupling', 'Memory allocation error' )
     ELSE
        CALL INFO('Memory allocation done in USF_TransferCoupling', Message, level=1 )
     END IF
     AllocationsDone = .TRUE. 
     OldCoupling(:,:) = 0.0
  END IF

  !----------------------------------------
  !Get Variables needed for the computation
  !----------------------------------------

  Constant => GetConstants()
  DIM = CoordinateSystemDimension()
  SolverName = "USF_TransferCoupling"
  Timevar => VariableGet( Model % Variables,'Time')
  tps = TimeVar % Values(1)
  !Get Latent heat of fusion per unit mass of ice
  !----------------------------------------------
  LFusion = GetConstReal( Constant, 'Massic Latent heat of Fusion for Ice', Found )
  IF (.NOT.Found) THEN
     LFusion = 3.33e20
  END IF

  !Get the CEL water load gradient
  !-------------------------------
  WaterLoadGrad1 => VariableGet( Model % Variables, 'CEL Wload Grad 1' )
  IF ( ASSOCIATED( WaterLoadGrad1 ) ) THEN
     WLG1Perm    => WaterLoadGrad1 % Perm
     WLG1Values  => WaterLoadGrad1 % Values
  ELSE
     CALL FATAL('RothPhi', 'Need  to compute gradient in FluxSolver with target variable CEL Wload !!!') 
  END IF

  WaterLoadGrad2 => VariableGet( Model % Variables, 'CEL Wload Grad 2' )
  IF ( ASSOCIATED( WaterLoadGrad2 ) ) THEN
     WLG2Perm    => WaterLoadGrad2 % Perm
     WLG2Values  => WaterLoadGrad2 % Values
  ELSE
     CALL FATAL('RothPhi', 'Need  to compute gradient in FluxSolver with target variable CEL Wload !!!') 
  END IF

  !Get CEL  water load
  !-----------------------------
  CELWaterLoad => VariableGet( Model % Variables, 'CEL Wload' )
  IF ( ASSOCIATED( CELWaterLoad ) ) THEN
     CELLoadPerm    => CELWaterLoad % Perm
     CELLoadValues  => CELWaterLoad % Values
  ELSE
     CALL FATAL('RothPhi', 'Need  to compute CEL water load in Hydrosolver !!!') 
  END IF 

  !-----------------------------------------
  !Get parameters needed for the computation
  !-----------------------------------------
  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  n = GetElementNOFNodes()

  !Get Glen parameters
  !--------------------
  
  !Computation of Arhenius factor

  ViscosityFlag = ListGetString( Material,'Viscosity Model', Found)
  IF(.NOT.Found)THEN
     GlenFactor = 7.568d1
     WRITE(Message,'(a)') 'No Viscosity model fixing Glen A to 7.568d1'
  ELSE
     SELECT CASE( ViscosityFlag )
     CASE('glen')
       !Assuming null Temperature to compute the Glen Factor
        Temperature = 0.0
        
        R = GetConstReal(Constant,'Gas Constant',Found)
        IF (.NOT.Found) R = 8.314_dp
        ! lets for the time being have this hardcoded
        Tlimit = GetConstReal(Material, 'Limit Temperature', Found)
        IF (.NOT.Found) THEN
           Tlimit = -10.0_dp
           CALL INFO('EffectiveViscosity','Limit Temperature not found. Setting to -10', Level=5)
        END IF
        A1 = GetConstReal(Material, 'Rate Factor 1', Found)
        IF (.NOT.Found) THEN
           A1 = 3.985d-13
           CALL INFO('EffectiveViscosity','Rate Factor 1 not found. Setting to 3.985e-13', Level=5)
        END IF
        A2 = GetConstReal(Material, 'Rate Factor 2', Found)
        IF (.NOT.Found) THEN
           A2 = 1.916d03
           CALL INFO('EffectiveViscosity','Rate Factor 2 not found. Setting to 1.916E03', Level=5)
        END IF
        Q1 = GetConstReal(Material, 'Activation Energy 1', Found)
        IF (.NOT.Found) THEN
           Q1 = 60.0d03
           CALL INFO('EffectiveViscosity','Activation Energy 1 not found. Setting to 60.0E03', Level=5)
        END IF
        Q2 = GetConstReal(Material, 'Activation Energy 2', Found)
        IF (.NOT.Found) THEN
           Q2 = 139.0d03
           CALL INFO('EffectiveViscosity','Activation Energy 2 not found. Setting to 139.0d03', Level=5)
        END IF
        
        IF (Temperature .LE. Tlimit) THEN
           GlenFactor = A1 * EXP( -Q1/(R * (273.15 + Temperature)))
        ELSE IF((Tlimit .LT. Temperature) .AND. (Temperature .LE. 0.0_dp)) THEN
           GlenFactor = A2 * EXP( -Q2/(R * (273.15 + Temperature)))
        ELSE
           GlenFactor = A2 * EXP( -Q2/(R * (273.15)))
           CALL INFO('EffectiveViscosity','Positive Temperature detected in Glen - limiting to zero!', Level = 5)
        END IF

     CASE DEFAULT 
        GlenFactor = 7.568d1
        CALL WARN('EffectiveViscosity','Unknown material model, Fixing GlenFactor to 7.568d1')
     END SELECT
  END IF

  GlenExp(1:n) = ListGetReal( Material, 'Glen Exponent',  n, Element % NodeIndexes, Found )
  IF (.NOT.Found) THEN
     GlenExp = 3.0
     WRITE(Message,'(a)') 'Keyword >Glen Exponent< not found Set to 3 '
     CALL INFO(SolverName,Message,Level=4)
  END IF

  !Get Ice Density
  !---------------        
  IceDensity(1:n) = ListGetReal( Material, 'Density',  n, Element % NodeIndexes, Found )
  IF (.NOT.Found) THEN
     IceDensity = 910.0D00
     WRITE(Message,'(a,a,i5)') 'Keyword >Density< not found Set to 910 kg/m3 '
     CALL INFO(SolverName,Message,Level=4)
  END IF

  !Get Water Density
  !-----------------
  WaterDensity(1:n) = ListGetReal( Material, 'Water Density',  n, Element % NodeIndexes, Found )
  IF (.NOT.Found) THEN
     WaterDensity = 1000.0D00
     WRITE(Message,'(a)') 'Keyword >Water Density< not found Set to 1000 kg/m3 '
     CALL INFO(SolverName,Message,Level=4)
  END IF

  !Get acceleration of gravity
  !---------------------------
  BodyForce => GetBodyForce()    
  IF ( ASSOCIATED( BodyForce ) ) THEN
     bf_id = GetBodyForceId()
     g = 0.0_dp  
     g(1,1:n) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 1',Found)
     g(2,1:n) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Found)
     g(3,1:n) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Found)
     Gravity(1:n) = SQRT(SUM(g(:,1:n)**2.0/n))
  END IF

  ! Get Max water load (normaly normal basal constraint)
  !-----------------------------------------------------
  Hmax(1:n) = ListGetReal( Material, 'Hwater Upper Limit',  n, Element % NodeIndexes, Found )
  IF (.NOT.Found) THEN
     CALL FATAL('RothPhi', 'Need a water load upper bound to compute Effective Pressure !!!') 
  END IF

  !Get CEL transmitivity
  !---------------------
  CALL ListGetRealArray( Material,'CEL Transmitivity',Hwrk,N, Element % NodeIndexes )
  Transmitivity = 0.0D0

  IF (SIZE(Hwrk,1) .EQ. 1) THEN
     DO i=1,3
        DO j=1,N
           Transmitivity(i,i,j) = Hwrk(1,1,j)
        END DO
     END DO
  ELSE
     WRITE(Message,'(a)') 'Keyword >CEL Transmitivity< should be isotrop '
     CALL INFO(SolverName,Message,Level=4)

  END IF
  ! Get the timeStep
  !-----------------
  Timestep = Model % Solver % dt 

  !--------------------
  !Start the real Stuff
  !--------------------

  !Computation of the effective pressure using the upper limit of water load which should be Snn
  !---------------------------------------------------------------------------------------------

  DO i=1,N

     IF (Element % Nodeindexes(i).EQ.nodenumber)THEN
        
        IF(tps.NE.OldCoupling(nodenumber,2))THEN
           OldCoupling(nodenumber,2) = tps
           OldCoupling(nodenumber,1) = OldCoupling(nodenumber,3)
        END IF
        
        IF(DIM.EQ.2)THEN 


           CELLoad = (CELLoadValues(CELLoadPerm(Element % Nodeindexes(i))) &
                -(Model % Nodes % y(Element % NodeIndexes(i)))) &
                * WaterDensity(i) * Gravity(i)

           EffectivePress = (Hmax(i)-(Model % Nodes % y(Element % NodeIndexes(i)))) &
                * WaterDensity(i) * Gravity(i) &
                - MAX(0.0,CELLoad)

        ELSEIF(DIM.EQ.3)THEN

           CELLoad = (CELLoadValues(CELLoadPerm(Element % Nodeindexes(i))) &
                -(Model % Nodes % z(Element %NodeIndexes(i)))) &
                * WaterDensity(i) * Gravity(i)

           EffectivePress = (Hmax(i)-(Model%Nodes%z(Element % NodeIndexes(i)))) &
                * WaterDensity(i) * Gravity(i) &
                - MAX(0.0,CELLoad)

        END IF

        IF (ABS(EffectivePress).GT.1e-6)THEN
           
           PressRatio = (EffectivePress/GlenExp(i))**(-GlenExp(i))
           
        ELSE
           PressRatio = (1e-6/GlenExp(i))**(-GlenExp(i))
        END IF
           

        GradFactor =  WaterDensity(i) * Gravity(i) * Transmitivity(1,1,i) & 
             *(WLG1Values(WLG1Perm(Element %Nodeindexes(i))) &
             + WLG2Values(WLG2Perm(Element %Nodeindexes(i))))**2.0


        Coupling = (Timestep / (2.0 * Timestep * GlenFactor(i) + PressRatio)) &
             * (GradFactor + OldCoupling(nodenumber,1)/Timestep) &
             * PressRatio
        

        OldCoupling(nodenumber,3) = Coupling

     END IF
  END DO

END FUNCTION RothPhi
