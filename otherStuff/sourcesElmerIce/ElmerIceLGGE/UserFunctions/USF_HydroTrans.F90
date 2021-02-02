  ! Computation of the transmitivities depending on the state of the efficient layer
  !-closed efficient system have sediment transmitivity
  !-active efficient layer have a mean(geometric) transmitivity 
  !-efficient effectiv layer have pipe transmitivity
  !Transmitivity is Isotrop
  !--------------------------------------------------------------------------------
  FUNCTION OpenCELVarTrans(Model,nodenumber,x) RESULT(Transmitivity)
  
  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  
  IMPLICIT NONE
!------------------------------------------------------------  
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: x, Transmitivity
!------------------------------------------------------------  
  TYPE(Element_t), POINTER :: Element 
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Variable_t), POINTER ::CELState
  REAL(KIND=dp) ::MaxTrans, MinTrans, SedThick, CELThick

  LOGICAL :: FirstTime = .TRUE., GotIt
  INTEGER :: n, i
  INTEGER, POINTER :: NodeIndexes(:), CELPerm(:)
  REAL(KIND=dp), POINTER :: CELValues(:)
  REAL(KIND=dp), ALLOCATABLE ::  AuxReal(:)
 
  SAVE :: FirstTime, AuxReal 
  
  IF (FirstTime) THEN
     FirstTime = .FALSE.
     n = Model % MaxElementNodes
      ALLOCATE(AuxReal(n))
   ENDIF
  
   Element => Model % CurrentElement 
   n = GetElementNOFNodes(Element)
   NodeIndexes => Element % NodeIndexes
   Material => GetMaterial(Element)

   !Read Max and min Transmitivities, max should 
   !correspond to the transmitivity of the efficient layer 
   !and min to the one of the uneficient layer
   !-----------------------------------------------------

   MaxTrans = GetConstReal( Material, 'Maximum CEL Transmitivity', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('OpenCELVarTrans', 'Need a Maximum CEL Transmitivity')
   END IF

   MinTrans = GetConstReal( Material, 'Minimum CEL Transmitivity', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('OpenCELVarTrans', 'Need a Minimum CEL Transmitivity')
   END IF

   !Getting parameters to compute the transmitivity
   !-----------------------------------------------------
   
   AuxReal(1:n) = ListGetReal( Material, &
        'Sediment Thickness', n, NodeIndexes, GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(A)') 'Sediment Thickness not found. &
           Setting to 10.0 m'
      CALL INFO('OpenCELVarTrans', Message, Level = 20)
      AuxReal(1:n) = 10.0_dp
   END IF 
   DO i=1, n
      IF (NodeNumber .EQ. NodeIndexes(i)) EXIT 
   END DO
   SedThick = AuxReal(i)

   AuxReal(1:n) = ListGetReal( Material, &
        'CEL Thickness', n, NodeIndexes, GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(A)') 'CEL Thickness not found. &
           Setting to 10.0 m'
      CALL INFO('OpenCELVarTrans', Message, Level = 20)
      AuxReal(1:n) = 1.0_dp
   END IF 
   DO i=1, n
      IF (NodeNumber .EQ. NodeIndexes(i)) EXIT 
   END DO
   CELThick = AuxReal(i)

   !Getting the variables from HydroSolver
   !-----------------------------------------------------

   CELState => VariableGet( Model % Variables, 'Open CEL' )
   IF ( ASSOCIATED( CELState ) ) THEN
      CELPerm   => CELState % Perm
      CELValues => CELState % Values
   ELSE
       CALL FATAL('OpenCELVarTrans', 'Need an Open CEL State Variable')   
   END IF


   !Computing Transmitivities
   !------------------------------------------------------
     
   !pipe is open and active, hight transmitivity (the one of the pipe)
   IF ( CELValues(CELPerm(nodenumber)).GE.2.0) THEN
      
      Transmitivity = MaxTrans
      
      !pipe is open and not active, mean transmitivity beetween pipe and sed transmitivity
   ELSEIF(CELValues(CELPerm(nodenumber)).EQ.1.0)THEN
      
      Transmitivity = (MaxTrans + MinTrans)/(CELThick + Sedthick)
            
      !CEL is closed but a node of the element is open,mean transmitivity   
   ELSEIF(SUM(CELValues(CELPerm(Element % NodeIndexes(1:n)))).NE.0.0 &
        .AND.(CELValues(CELPerm(nodenumber)).EQ.0.0)) THEN

         
      Transmitivity = (MaxTrans + MinTrans)/(CELThick + Sedthick)
      
      
      
      ! CEL is closed, low transmitivity (the one of the sediment)
   ELSEIF(SUM(CELValues(CELPerm(Element % NodeIndexes(1:n)))).EQ.0.0) THEN
      
      Transmitivity = MinTrans
      
   END IF

 END FUNCTION OpenCELVarTrans



 !--------------------------------------------------------------------------------
 FUNCTION NVarTrans(Model,nodenumber,x) RESULT(Transmitivity)
  
   USE DefUtils
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   
   IMPLICIT NONE
   !------------------------------------------------------------  
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(kind=dp) :: x, Transmitivity
   !------------------------------------------------------------  
   TYPE(Element_t), POINTER :: Element 
   TYPE(ValueList_t), POINTER :: Material
   TYPE(Variable_t), POINTER ::CELState, SedLoad, DepthVariable
   REAL(KIND=dp) ::MaxTrans, MinTrans, SedThick, CELThick, &
        Density, WaterDensity, ratio
   LOGICAL :: FirstTime = .TRUE., GotIt
   INTEGER :: n, i
   INTEGER, POINTER :: NodeIndexes(:), CELPerm(:), SedLoadPerm(:), &
        DepthPerm(:)
   REAL(KIND=dp), POINTER :: CELValues(:), SedLoadValues(:), &
        DepthValues(:)
   REAL(KIND=dp), ALLOCATABLE ::  AuxReal(:)
   
   SAVE :: FirstTime, AuxReal 
   
   IF (FirstTime) THEN
      FirstTime = .FALSE.
      n = Model % MaxElementNodes
      ALLOCATE(AuxReal(n))
   ENDIF
   
   Element => Model % CurrentElement 
   n = GetElementNOFNodes(Element)
   NodeIndexes => Element % NodeIndexes
   Material => GetMaterial(Element)

   !Read Max and min Transmitivities, max should 
   !correspond to the transmitivity of the efficient layer 
   !and min to the one of the uneficient layer
   !-----------------------------------------------------

   MaxTrans = GetConstReal( Material, 'Maximum CEL Transmitivity', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('OpenCELVarTrans', 'Need a Maximum CEL Transmitivity')
   END IF

   MinTrans = GetConstReal( Material, 'Minimum CEL Transmitivity', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('OpenCELVarTrans', 'Need a Minimum CEL Transmitivity')
   END IF

   !Getting parameters to compute the transmitivity
   !-----------------------------------------------------
   
   AuxReal(1:n) = ListGetReal( Material, &
        'Sediment Thickness', n, NodeIndexes, GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(A)') 'Sediment Thickness not found. &
           Setting to 10.0 m'
      CALL INFO('NVarTrans', Message, Level = 20)
      AuxReal(1:n) = 10.0_dp
   END IF 
   DO i=1, n
      IF (NodeNumber .EQ. Model % CurrentElement % NodeIndexes(i)) EXIT 
   END DO
   SedThick = AuxReal(i)

   AuxReal(1:n) = ListGetReal( Material, &
        'CEL Thickness', n, NodeIndexes, GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(A)') 'CEL Thickness not found. &
           Setting to 10.0 m'
      CALL INFO('NVarTrans', Message, Level = 20)
      AuxReal(1:n) = 1.0_dp
   END IF 
   DO i=1, n
      IF (NodeNumber .EQ. Model % CurrentElement % NodeIndexes(i)) EXIT 
   END DO
   CELThick = AuxReal(i)
   
   AuxReal(1:n) = GetReal( Material,'Density', GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(A)') 'Density not found. &
           Setting to 910 kg/m3'
      CALL INFO('Stress_bound', Message, Level = 20)
      AuxReal(1:n) = 910.0_dp
   END IF
   DO i=1, n
      IF (NodeNumber.EQ.Model % CurrentElement % NodeIndexes( i )) EXIT 
   END DO
   Density = AuxReal(i)

   AuxReal(1:n) = GetReal( Material,'Water Density', GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(A)') 'Water Density not found. &
           Setting to 1000 kg/m3'
      CALL INFO('Stress_bound', Message, Level = 20)
      AuxReal(1:n) = 1000.0_dp
   END IF
   DO i=1, n
      IF (NodeNumber.EQ.Model % CurrentElement % NodeIndexes( i )) EXIT 
   END DO
   WaterDensity = AuxReal(i)

   !Getting the variables from HydroSolver
   !-----------------------------------------------------

   CELState => VariableGet( Model % Variables, 'Open CEL' )
   IF ( ASSOCIATED( CELState ) ) THEN
      CELPerm   => CELState % Perm
      CELValues => CELState % Values
   ELSE
       CALL FATAL('NVarTrans', 'Need an Open CEL State Variable')   
   END IF
   SedLoad => VariableGet( Model % Variables, 'Sediment Wload' )
   IF ( ASSOCIATED( CELState ) ) THEN
      SedLoadPerm   => SedLoad % Perm
      SedLoadValues => SedLoad % Values
   ELSE
       CALL FATAL('NVarTrans', 'Need a Sedient Wload')   
   END IF  
   DepthVariable => VariableGet( Model % Variables, 'Depth' )
   IF ( ASSOCIATED( DepthVariable ) ) THEN
      DepthPerm    => DepthVariable % Perm
      DepthValues  => DepthVariable % Values
   ELSE
      CALL FATAL('NVarTrans', 'Need Flowdepth Solver , Depth not associated !!!') 
   END IF
   
   !Computing Transmitivities
   !------------------------------------------------------
     
   !pipe is open and active, hight transmitivity (the one of the pipe)
   IF ( CELValues(CELPerm(nodenumber)).GE.2.0) THEN
      
      Transmitivity = MaxTrans
      
   ELSE
      ratio = SedLoadValues(SedLoadPerm(nodenumber))/&
           ((Density * DepthValues(DepthPerm(Nodenumber))/ WaterDensity) &
           + (Model%Nodes%z(nodenumber)))
          

      Transmitivity = CELThick * ((MaxTrans * ratio + MinTrans)/(CELThick + Sedthick))
       
   END IF

 END FUNCTION NVarTrans
