!-------------------------------------------------------------------------
!Initialisation of the pore Pressure, Pressure at the base of the sediment
!-------------------------------------------------------------------------

FUNCTION InitPP(Model,nodenumber,x) RESULT(Press)
  
  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  
  IMPLICIT NONE

!------------------------------------------------------------  
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: x, Press
!------------------------------------------------------------  
  TYPE(Element_t), POINTER :: Element 
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp) :: gravity, g(3), rhow, st

  LOGICAL :: FirstTime = .TRUE., GotIt 
  INTEGER :: n, i, bf_id 
  INTEGER, POINTER :: NodeIndexes(:)
  REAL(KIND=dp), ALLOCATABLE :: WaterDensity(:), SedThick(:)

  SAVE :: FirstTime, gravity 
  SAVE :: WaterDensity, SedThick


   IF (FirstTime) THEN
      FirstTime = .FALSE.
      bf_id = ListGetInteger( Model % Bodies(1) % Values,&
           'Body Force',  minv=1, maxv=Model % NumberOFMaterials)
      g = 0.0_dp   
      g(1) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 1',Gotit)
      g(2) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Gotit)
      g(3) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Gotit)
      gravity = SQRT(SUM(g**2.0))
        
      
      n = Model % MaxElementNodes
      ALLOCATE(WaterDensity(n), SedThick(n))
   ENDIF

   Element => Model % CurrentElement 
   n = GetElementNOFNodes(Element)
   NodeIndexes => Element % NodeIndexes
   Material => GetMaterial(Element)

   WaterDensity(1:n) = ListGetReal( Material, &
                         'Water Density', n, NodeIndexes, GotIt )
   IF (.NOT.GotIt) THEN
       WRITE(Message,'(A)') 'Water Density not found. &
                            Setting to 1000 kg/m3'
       CALL INFO('InitPP', Message, Level = 20)
       WaterDensity(1:n) = 1000.0_dp
   END IF

   SedThick(1:n) = ListGetReal( Material, &
                         'Sediment Thickness', n, NodeIndexes, GotIt )
   IF (.NOT.GotIt) THEN
       WRITE(Message,'(A)') 'Sediment Thickness not found. &
                            Setting to 10.0 m'
       CALL INFO('InitPP', Message, Level = 20)
       SedThick(1:n) = 10.0_dp
   END IF
      
   DO i= 1, n
     IF (NodeIndexes(i)==nodenumber) THEN
             st = SedThick(i)
             rhow = WaterDensity(i)
             EXIT
     END IF 
   END DO
     
   Press = gravity * rhow * st 

END FUNCTION InitPP


!------------------------------------------
!Initialisation of the free surface
!------------------------------------------

FUNCTION DyIni ( Model, nodenumber, x) RESULT(Dy)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i
   REAL(KIND=dp) :: x,   Dy
   REAL(KIND=dp), ALLOCATABLE :: Dy0(:)
   LOGICAL :: FirstTime=.True.

   SAVE FirstTime
   SAVE Dy0

   IF (FirstTime) THEN
          FirstTime=.False.
          NMAX = Model % NumberOfNodes
          ALLOCATE(Dy0(NMAX))
          Do i = 1, NMAX
          Dy0(i) = Model % Nodes % y (i)
          EndDo
   END IF

   Dy = Dy0(nodenumber)

END FUNCTION DyIni

!----------------------------------------------
!Minimum value for the altitude of the surface
!TO be done
!----------------------------------------------

FUNCTION Input ( Model, nodenumber, Depth) RESULT(flux)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i
   REAL(KIND=dp) :: Depth,   flux

   if((Model%Nodes%y(nodenumber)+Depth-1070).GT.0.0)THEN
      flux = 0.02
   elseif((Model%Nodes%y(nodenumber)+Depth-1070).GT.-150.0)THEN
      flux = 3.72-((Model%Nodes%y(nodenumber)+Depth-1070)*0.0082)
   else
      flux = 0.02-((Model%Nodes%y(nodenumber)+Depth-1070)*0.0082)
   end if

 END FUNCTION Input

!---------------------------------------------------------------------------
!Fixing an upper bound to the water load in agreement with the normal stress
!---------------------------------------------------------------------------

FUNCTION Stress_bound (Model,nodenumber,Depth) RESULT(hwat_bound)


  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL (KIND=dp) :: Depth              
  INTEGER :: nodenumber
  TYPE(ValueList_t), POINTER :: Material,BC
  TYPE(Variable_t), POINTER :: StressVariable, NormalVar, FlowVariable, DepthVariable
  REAL(KIND=dp), POINTER :: StressValues(:), NormalValues(:), FlowValues(:), DepthValues(:)
  INTEGER, POINTER :: StressPerm(:), NormalPerm(:), FlowPerm(:), DepthPerm(:)
  INTEGER :: DIM, i, j, Ind(3,3), n, bf_id
  REAL (KIND=dp) :: gravity, g(3), WaterDensity, Density
  REAL (KIND=dp) :: Snn, MaxWaterPress, hwat_bound
  REAL (KIND=dp), ALLOCATABLE :: Sig(:,:), normal(:), Sn(:), AuxReal(:) 
  LOGICAL :: GotIt, FirstTime = .TRUE., Cauchy

  SAVE :: Sig, normal, DIM, Ind, Sn, gravity, FirstTime

  IF (FirstTime) THEN
     FirstTime = .FALSE.  
     DIM = CoordinateSystemDimension()
     IF ((DIM == 2).OR.(DIM == 3))  THEN
        ALLOCATE(Sig(DIM,DIM),normal(DIM), Sn(DIM))
     ELSE
        CALL FATAL('Stress_bound', 'Bad dimension of the problem')
     END IF
     DO i=1, 3
        Ind(i,i) = i
     END DO
     Ind(1,2) = 4
     Ind(2,1) = 4
     Ind(2,3) = 5
     Ind(3,2) = 5
     Ind(3,1) = 6
     Ind(1,3) = 6
     
     bf_id = ListGetInteger( Model % Bodies(1) % Values,&
          'Body Force',  minv=1, maxv=Model % NumberOFMaterials)
     g = 0.0_dp  
     g(1) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 1',Gotit)
     g(2) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Gotit)
     g(3) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Gotit)
     gravity = SQRT(SUM(g**2.0))
  END IF


     n = GetElementNOFNodes()
     Material => GetMaterial(Model % CurrentElement)
     !get the water density

     ALLOCATE(AuxReal(n))
     AuxReal(1:n) = GetReal( Material,'Water Density', GotIt )
     IF (.NOT.GotIt) THEN
        WRITE(Message,'(A)') 'Water Density not found. &
             Setting to 1000 kg/m3'
        CALL INFO('Stress_bound', Message, Level = 20)
        AuxReal(1:n) = 1000.0_dp
     END IF
     DO i=1, n
        IF (NodeNumber== Model % CurrentElement % NodeIndexes( i )) EXIT 
     END DO
     WaterDensity = AuxReal(i)

     AuxReal(1:n) = GetReal( Material,'Density', GotIt )
     IF (.NOT.GotIt) THEN
        WRITE(Message,'(A)') 'Density not found. &
             Setting to 910 kg/m3'
        CALL INFO('Stress_bound', Message, Level = 20)
        AuxReal(1:n) = 910.0_dp
     END IF
     DO i=1, n
        IF (NodeNumber== Model % CurrentElement % NodeIndexes( i )) EXIT 
     END DO
     Density = AuxReal(i)

     DEALLOCATE(AuxReal)

     ! Get the variables to compute the normal stress
     StressVariable => VariableGet( Model % Variables, 'Stress' )
     IF ( ASSOCIATED( StressVariable ) ) THEN
        StressPerm    => StressVariable % Perm
        StressValues  => StressVariable % Values
     ELSE
        CALL FATAL('Stress_bound', 'Need ComputeDevStressNS Solver, Stress not associated !!!') 
     END IF
     !And the depth
     DepthVariable => VariableGet( Model % Variables, 'Depth' )
     IF ( ASSOCIATED( DepthVariable ) ) THEN
        DepthPerm    => DepthVariable % Perm
        DepthValues  => DepthVariable % Values
     ELSE
        CALL FATAL('Stress_bound', 'Need Flowdepth Solver , Depth not associated !!!') 
     END IF
     ! Cauchy or deviatoric stresses ?
     !
     Cauchy = ListGetLogical( material , 'Cauchy', Gotit )

     ! Get the variables to compute deviatoric part of the stress

     FlowVariable => VariableGet( Model % Variables, 'Flow Solution' )
     IF ( ASSOCIATED( FlowVariable ) ) THEN
        FlowPerm    => FlowVariable % Perm
        FlowValues  => FlowVariable % Values
     ELSE
        CALL FATAL('Stress_bound', 'Need NS Solver, Flow Solution not associated !!!')
     END IF

     ! Get the variable to compute the normal
     NormalVar =>  VariableGet(Model % Variables,'Normal Vector')
     IF ( ASSOCIATED( NormalVar ) ) THEN
        NormalPerm => NormalVar % Perm
        NormalValues => NormalVar % Values

     ELSE
        CALL FATAL('Stress_bound', 'Need ComputeNormal Solver, Normal Vector not associated !!!') 
     END IF


     DO i=1, DIM
        normal(i) = -NormalValues(DIM*(NormalPerm(Nodenumber)-1) + i)
     END DO

     DO i=1, DIM
        DO j= 1, DIM
           Sig(i,j) =  &
                StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(i,j) )
        END DO


        IF (.NOT.Cauchy) THEN 
           Sig(i,i) = Sig(i,i) - FlowValues((DIM+1)*FlowPerm(Nodenumber))

        END IF
     END DO

     ! Stress vector Sn       
     DO i=1, DIM
        Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
     END DO

     Snn = SUM( Sn(1:DIM) * normal(1:DIM) )


!For the first timestep, when Sigma is not known we use the hydrostatic hypotesis
     IF ((Snn.EQ.0.0).OR.(GetTimeStep().EQ.1)) THEN
        Snn = - gravity * Density * DepthValues(DepthPerm(Nodenumber))

!Faking Snn for first and last nodes

     ELSEIF((Model%Nodes%x(nodenumber).EQ.MAXVAL(Model%Nodes%x))&
          .OR.Model%Nodes%x(nodenumber).EQ.MINVAL(Model%Nodes%x))THEN
        Snn = - gravity * Density * DepthValues(DepthPerm(Nodenumber))
     END IF


     !From the normal stress compute the maximum water load

     MaxWaterPress = -Snn


     IF(DIM.EQ.2)THEN 

        hwat_bound =(MaxWaterPress / (WaterDensity * gravity)) + (Model%Nodes%y(nodenumber))


     ELSEIF(DIM.EQ.3)THEN
        
        hwat_bound =(MaxWaterPress / (WaterDensity * gravity)) + (Model%Nodes%z(nodenumber))
     END IF



END FUNCTION Stress_bound


FUNCTION InitCEL(Model,nodenumber,x) RESULT(Snout)
  
  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  
  IMPLICIT NONE

!------------------------------------------------------------  
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: x, Snout
!------------------------------------------------------------  
 
  IF(Model%Nodes%y(nodenumber).EQ.400)THEN
     Snout = 3.0
  ELSE
     Snout = 0.0
  END IF


END FUNCTION InitCEL


FUNCTION RandBound(Model,nodenumber,x) RESULT(Bound)
  
  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  
  IMPLICIT NONE

!------------------------------------------------------------  
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: x, Bound
!------------------------------------------------------------  
  
  INTEGER :: M
  LOGICAL :: FirstTime = .TRUE.
  REAL(KIND=dp) :: res
  REAL(KIND=dp), ALLOCATABLE ::BoundArray(:)

  SAVE :: FirstTime 
  SAVE :: BoundArray

  IF (FirstTime)THEN
     FirstTime = .FALSE.
     IF(ALLOCATED(BoundArray))DEALLOCATE(BoundArray)
     M = Model % NumberOfNodes
     Allocate(BoundArray(M))
     BoundArray = -1.0
  END IF
     
  IF(BoundArray(nodenumber).LT.0.0)THEN

     res = RAND(1)
     BoundArray(nodenumber) = Model%Nodes%y(nodenumber) + res
  END IF

  Bound = BoundArray(nodenumber)

END FUNCTION RandBound


