!/*****************************************************************************/! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - Scientific Computing Ltd., Finland
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module containing solvers and routines for standard ice dynamic problems
! *
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Mikko Lyly, Juha Ruokolainen
! *  Email:   Thomas.Zwinger@csc.fi, Juha.Ruokolainen@csc.fi 
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - Scientific Computing Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 14 May 2007
! *
! *****************************************************************************/
! * 2012/09/26 : Rewriting from Scratch step by step from HydroSolver at CSC
! *****************************************************************************/

RECURSIVE SUBROUTINE HydroSolver( Model,Solver,Timestep,TransientSimulation )
  !------------------------------------------------------------------------------

  USE DiffuseConvective
  USE DiffuseConvectiveGeneral
  USE Differentials
  USE MaterialModels
  USE DefUtils
  !------------------------------------------------------------------------------
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  !  Solve the convection diffusion equation with limiters!
  !
  !  ARGUMENTS:
  !
  !  TYPE(Model_t) :: Model,  
  !     INPUT: All model information (mesh,materials,BCs,etc...)
  !
  !  TYPE(Solver_t) :: Solver
  !     INPUT: Linear equation solver options
  !
  !  REAL(KIND=dp) :: Timestep
  !     INPUT: Timestep size for time dependent simulations
  !
  !******************************************************************************

  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Matrix_t), POINTER :: SystemMatrix
  TYPE(Element_t), POINTER :: Element
  TYPE(Variable_t), POINTER :: WloadSol, VarWload, VarWloadResidual, &
       WaterPressure
  TYPE(ValueList_t), POINTER :: Constants, SolverParams, Equation, &
       Material, BodyForce,BC
  TYPE(Nodes_t) :: ElementNodes  

  REAL(KIND=dp), POINTER :: Wload(:), PointerToResidualVector(:), &
       ForceVector(:), Hwrk(:,:,:), WloadHomologous(:), &
       ResidualVector(:), Wpress(:)

  REAL(KIND=dp), ALLOCATABLE :: Transmitivity(:,:,:) 
  REAL(KIND=dp), ALLOCATABLE :: g(:,:), Nochange(:,:)
  REAL(KIND=dp), ALLOCATABLE :: UpperLimit(:), SedComp(:), &
       Porosity(:), Gravity(:), SedThick(:) , Density(:), &
       StoringCoef(:), Viscosity(:), C0(:), C1(:), Zero(:), &
       VNull(:), Work(:), Pressure(:), TransferCoeff(:), &
       WloadExt(:), OldValues(:), OldRHS(:),StiffVector(:), &
       Influx(:), CELToSed(:)
  REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), &
       FORCE(:),LOAD(:), TimeForce(:)

  REAL(KIND=dp) :: LinearTol, NonlinearTol, & !CumulativeTime,
       totat, totst, st
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

  REAL(KIND=dp) :: WatComp, Norm, PrevNorm, RelativeChange

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, VariableName

  INTEGER, POINTER :: WloadPerm(:), PTRPerm(:), HomolPerm(:),&
       WpPerm(:)
  INTEGER :: DIM, NonlinearIter, LocalNodes, &
       PenIter, k, N, M, L, t, i, j, istat
  INTEGER :: body_id, material_id, bf_id, bc_id

  LOGICAL, ALLOCATABLE :: ActiveSed(:)

  LOGICAL ::Found = .FALSE., &
       Stabilize = .TRUE. ,&
       UseBubbles = .FALSE., &
       AllocationsDone = .FALSE., &
       FluxBC = .FALSE., &
       IsPeriodicBC = .FALSE., &
       ApplyDirichlet = .FALSE.,&
       FirstTime = .TRUE.

  SAVE LinearTol,       &
       NonlinearTol,    &
       SolverName,      &
       VariableName,    &
       ElementNodes,    &
       ResidualVector,  &
       Wpress,          &
       Transmitivity,   &
       g,               &
       Nochange,        &
       UpperLimit,      &
       SedComp,         &
       Porosity,        &
       Gravity,         &
       SedThick,        &
       Density,         &
       StoringCoef,     &
       Viscosity,       &
       C0,              &
       C1,              &
       Zero,            &
       VNull,           &
       Work,            &
       Pressure,        &
       TransferCoeff,   &
       WloadExt,        &
       OldValues,       &
       OldRHS,          &
       StiffVector,     &
       Influx,          &
       CELToSed,        &
       MASS,            &
       STIFF,           &
       FORCE,           &
       LOAD,            &
       TimeForce,       &
       ActiveSed,       & 
       Hwrk,            &
       AllocationsDone
  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  SolverName = 'HydroSolver ('// TRIM(Solver % Variable % Name) // ')'
  VariableName = TRIM(Solver % Variable % Name)

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  PointerToSolver => Solver
  SystemMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % RHS

  WloadSol => Solver % Variable
  WloadPerm  => WloadSol % Perm
  Wload => WloadSol % Values

  LocalNodes = COUNT( WloadPerm .GT. 0 )
  IF ( LocalNodes .LE. 0 ) RETURN  

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     N = Solver % Mesh % MaxElementNodes
     M = Model % Mesh % NumberOfNodes
     K = SIZE( SystemMatrix % Values )
     L = SIZE( SystemMatrix % RHS )

     IF ( AllocationsDone ) THEN  

        DEALLOCATE(            &
             ResidualVector,   &
             Transmitivity,    &
             g,                &
             Nochange,         &
             UpperLimit,       &
             SedComp,          &
             Porosity,         &
             Gravity,          &
             SedThick,         &
             Density,          &
             StoringCoef,      &
             Viscosity,        &
             C0,               &
             C1,               &
             Zero,             &
             VNull,            &
             Work,             &
             Pressure,         &
             TransferCoeff,    &
             WloadExt,         &
             OldValues,        &
             OldRHS,           &
             StiffVector,      &
             Influx,           &
             CELToSed,         &
             MASS,             &
             STIFF,            &
             FORCE,            &
             LOAD,             &
             TimeForce,        &
             ActiveSed,        &
             ElementNodes % x, &
             ElementNodes % y, &
             ElementNodes % z)          
     END IF
     ALLOCATE(                   &
          ResidualVector( L ),   &
          Transmitivity( 3,3,N ),&
          g( 3,N ),              &
          Nochange( 3,N ),       &
          UpperLimit( M ),       &
          SedComp( N ),          &
          Porosity( N ),         &
          Gravity( N ),          &
          SedThick( N ),         &
          Density( N ),          &
          StoringCoef( N ),      &
          Viscosity( N ),        &
          C0( N ),               &
          C1( N ),               &
          Zero( N ),             &
          VNull( N ),            &
          Work( N ),             &
          Pressure( N ),         &
          TransferCoeff( N ),    &
          WloadExt( N ),         &
          OldValues( K ),        &
          OldRHS( L ),           &
          StiffVector( L ),      &
          Influx( N ),           &
          CELToSed( N ),         &
          MASS( 2*N,2*N ),       &
          STIFF( 2*N,2*N ),      &
          FORCE( 2*N ),          &
          LOAD( N ),             &
          TimeForce( 2*N ),      &
          ActiveSed( M ),        &
          ElementNodes % x( N ), &
          ElementNodes % y( N ), &
          ElementNodes % z( N ), &
          STAT=istat )

     ActiveSed(:)=.FALSE.

     IF ( istat .NE. 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error in Sediment Solve' )
     ELSE
        CALL INFO(SolverName, 'Memory allocation done in Sediment Solve', level=1 )
     END IF
     AllocationsDone = .TRUE.
  END IF

  !------------------------------------------------------------------------------
  !    Say hello
  !------------------------------------------------------------------------------
  WRITE(Message,'(A,A)')&
       'Limited diffusion Solver for variable ', VariableName
  CALL INFO(SolverName,Message,Level=1)
  !------------------------------------------------------------------------------
  !    Read Iteration related constants
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  Stabilize = GetLogical( SolverParams,'Stabilize', Found )
  IF (.NOT. Found) THEN
     Stabilize = .FALSE.
  END IF

  UseBubbles = GetLogical( SolverParams,'Bubbles', Found )

  LinearTol = GetConstReal( SolverParams, &
       'Linear System Convergence Tolerance', Found )
  IF ( .NOT.Found ) THEN
     CALL FATAL(SolverName, 'No >Linear System Convergence Tolerance< found')
  END IF

  NonlinearIter = GetInteger( SolverParams, &
       'Nonlinear System Max Iterations', Found )
  IF ( .NOT.Found ) THEN
     NonlinearIter = 1
  END IF

  NonlinearTol = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance', Found )
  IF ( .NOT.Found ) THEN
     NonlinearTol = 1e-6
  END IF

  IF (FirstTime)THEN
     FirstTime = .FALSE.
     CALL CalculateNodalWeights(Solver,.TRUE.)
  END IF
  !------------------------------------------------------------------------------
  !    Read Physical constants
  !------------------------------------------------------------------------------
  Constants => GetConstants()

  WatComp = GetConstReal( Constants, &
       'Water Compressibility', Found)
  IF ( .NOT.Found ) THEN
     WatComp = 5.04e-4
  END IF

  totat = 0.0d0
  totst = 0.0d0

  !Variables needed to compute the residual and water pressure
  !----------------------------------------

  VarWloadResidual => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Residual' )
  IF (ASSOCIATED(VarWloadResidual)) THEN
     PointerToResidualVector => VarWloadResidual % Values
     PTRPerm => VarWloadResidual % Perm
  ELSE
     WRITE(Message,'(A)') '>' // TRIM(Solver % Variable % Name) // ' Residual< not associated'
     CALL FATAL( SolverName, Message)
  END IF

  VarWload => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Homologous' )
  IF (ASSOCIATED(VarWload)) THEN
     WloadHomologous => VarWload % Values
     HomolPerm => VarWload % Perm
  ELSE
     WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' VarWload not associated'
     CALL FATAL( SolverName, Message)
  END IF

  WaterPressure => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Pressure' )
  IF (ASSOCIATED(WaterPressure)) THEN
     Wpress => WaterPressure % Values
     WpPerm  => WaterPressure % Perm
  ELSE
     WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' WaterPressure not associated'
     CALL FATAL( SolverName, Message)
  END IF

  !Checking wether to use upper limit
  !------------------------------------
  SolverParams => GetSolverParams()

  ApplyDirichlet = GetLogical( SolverParams, &
       'Apply Dirichlet', Found)
  IF ( .NOT.Found ) THEN
     ApplyDirichlet = .FALSE.
  END IF
  IF (.NOT.ApplyDirichlet) THEN
     ActiveSed(:) = .FALSE.
  END IF

  !------------------------------------------------------------------------------
  !       non-linear system iteration loop
  !------------------------------------------------------------------------------
  DO Peniter=1,NonlinearIter

     !------------------------------------------------------------------------------
     ! print out some information
     !------------------------------------------------------------------------------
     at  = CPUTime()
     at0 = RealTime()
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     WRITE( Message,'(A,A,I3,A,I3)') &
          TRIM(Solver % Variable % Name),  ' iteration no.', Peniter,' of ',NonlinearIter
     CALL Info( SolverName, Message, Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, 'Starting Assembly...', Level=4 )

     !------------------------------------------------------------------------------
     ! lets start
     !------------------------------------------------------------------------------
     CALL DefaultInitialize()

      
     ! Get  Upper limit:
     !-----------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(t)
        
        n = GetElementNOFNodes()
        CALL GetElementNodes( ElementNodes )
        Material => GetMaterial()
        
        ! upper limit
        !------------
        UpperLimit(Element % Nodeindexes(1:n)) = ListGetReal(Material,TRIM(VariableName) // & 
             ' Upper Limit',n,Element % NodeIndexes, Found)
        
        IF (.NOT. Found) THEN
           WRITE(Message,'(a,i10)') 'No upper limit of solution for element no. ', t
           CALL INFO(SolverName, Message, level=10)
        END IF
     END DO

     !------------------------------------------------------------------------------
     ! write some info on max/min values
     !------------------------------------------------------------------------------    
     WRITE(Message,'(a,e13.6,a,e13.6)') &
          'Max/min values of sediment Water load: ', MAXVAL(Wload(:)), &
          '/',MINVAL(Wload(:))
     CALL INFO(SolverName,Message,Level=4)

     !------------------------------------------------------------------------------
     body_id = -1
     NULLIFY(Material)

     !------------------------------------------------------------------------------
     ! Bulk elements Loop
     !------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        !------------------------------------------------------------------------------
        ! Check if this element belongs to a body where scalar equation
        ! should be calculated
        !------------------------------------------------------------------------------

        Element => GetActiveElement(t,Solver) 

        IF (.NOT.ASSOCIATED(Element)) CYCLE 

        IF ( Element % BodyId .NE. body_id ) THEN
           Equation => GetEquation()
           IF (.NOT.ASSOCIATED(Equation)) THEN
              WRITE (Message,'(A,I3)') 'No Equation  found for boundary element no. ', t
              CALL FATAL(SolverName,Message)
           END IF

           Material => GetMaterial()
           IF (.NOT.ASSOCIATED(Material)) THEN
              WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
              CALL FATAL(SolverName,Message)
           ELSE
              material_id = GetMaterialId( Element, Found)
              IF(.NOT.Found) THEN
                 WRITE (Message,'(A,I3)') 'No Material ID found for boundary element no. ', t
                 CALL FATAL(SolverName,Message)
              END IF
           END IF
        END IF

        !------------------------------------------------------------------------------
        ! Get element material parameters
        !------------------------------------------------------------------------------       
        N = GetElementNOFNodes(Element)
        CALL GetElementNodes(ElementNodes,Element)
       
        SedComp(1:N) = listGetReal( Material,'Sediment Compressibility', N, Element % NodeIndexes, Found )
        IF (.NOT.Found) THEN
           SedComp(1:N) = 1.0D-2
           WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Compressibility', &
                '< not found for element ', t, ' material ', material_id
           CALL INFO(SolverName,Message,Level=4)
        END IF

        Porosity(1:N) = listGetReal( Material,'Sediment Porosity', N, Element % NodeIndexes, Found )
        IF (.NOT.Found) THEN
           Porosity(1:N) = 0.4D00
           WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Porosity', &
                '< not found for element ', t, ' material ', material_id
           CALL INFO(SolverName,Message,Level=4)
        END IF

        BodyForce => GetBodyForce()
        IF ( ASSOCIATED( BodyForce ) ) THEN
           bf_id = GetBodyForceId()
           g = 0.0_dp  
           g(1,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 1',Found)
           g(2,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Found)
           g(3,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Found)
           Gravity(1:N) = SQRT(SUM(g**2.0/N))
        END IF

        SedThick(1:N) = listGetReal( Material,'Sediment Thickness', N, Element % NodeIndexes, Found )
        IF (.NOT.Found) THEN
           SedThick(1:N) = 10.0D00
           WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Thickness', &
                '< not found for element ', t, ' material ', material_id
           CALL INFO(SolverName,Message,Level=4)
        END IF

        Density(1:N) = ListGetReal( Material, 'Water Density',  N, Element % NodeIndexes, Found )
        IF (.NOT.Found) THEN
           Density(1:N) = 1.0055e-18
           WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Water Density< not found for element ',&
                t, ' material ', material_id
           CALL INFO(SolverName,Message,Level=4)
        END IF


        !Computing the Storing coeficient and transmitivity of the layer 
        !----------------------------------------------------------------------------------
        CALL ListGetRealArray( Material,'Sediment Transmitivity',Hwrk,N, Element % NodeIndexes )
        Transmitivity = 0.0D0
        IF ( SIZE(Hwrk,1) == 1 ) THEN
           DO i=1,3
              Transmitivity( i,i,1:N ) = Hwrk( 1,1,1:N)
           END DO
        ELSE
           WRITE(Message,'(a,a,a)') 'Keyword >Sediment Transmitivity< should be isotrop '
           CALL INFO(SolverName,Message,Level=4)

        END IF

        StoringCoef(1:N) = SedThick(1:N) * &
             Gravity(1:N) * Porosity(1:N) * Density(1:N) * &
             (WatComp + SedComp(1:N)/Porosity(1:N))

        !------------------------------------------
        ! NB.: viscosity needed for strain heating
        !      but Newtonian flow is assumed
        !------------------------------------------
        Viscosity = 0.0D00
        !------------------------------------------------------------------------------
        ! no contribution proportional to Water load by default
        ! No convection by default give C1 equal 0
        !------------------------------------------------------------------------------
        C0 = 0.0d00
        C1 = 0.0d00
        !------------------------------------------------------------------------------
        ! add water contribution from input volume and transfer
        !------------------------------------------------------------------------------
        LOAD = 0.0D00

        BodyForce => GetBodyForce()
        IF (ASSOCIATED(BodyForce)) THEN
           bf_id = GetBodyForceId()
           Influx(1:N) = ListGetReal( BodyForce, TRIM(Solver % Variable % Name) // ' Source flux', &
                N, Element % NodeIndexes(1:N),Found ) 
           IF (.NOT.Found) THEN
              WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Source Flux', &
                   '< not found for element ', t, ' Body Force ', bf_id
              CALL INFO(SolverName,Message,Level=4)
              Influx(1:N) = 0.0
           END IF
           CELToSed(1:N) = ListGetReal( BodyForce,'CELToSed Transfer', &
                N, Element % NodeIndexes(1:N),Found ) 
           IF (.NOT.Found) THEN
              WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', 'CELToSed Transfer', &
                   '< not found for element ', t, ' Body Force ', bf_id
              CALL INFO(SolverName,Message,Level=4)
              CELToSed(1:N) = 0.0 
           END IF
        END IF

        LOAD(1:N) = Influx(1:N) + CELToSed(1:N)

        !------------------------------------------------------------------------------
        ! dummy input array for faking   heat capacity, density, temperature, 
        !                                enthalpy and viscosity
        ! Also faking the velocities to compute the Water load
        !------------------------------------------------------------------------------
        Work = 1.0d00
        Zero = 0.0D00
        VNull = 0.0D00
        Nochange = 0.0D00
        !------------------------------------------------------------------------------
        ! Get element local matrices, and RHS vectors
        !------------------------------------------------------------------------------
        MASS = 0.0d00
        STIFF = 0.0d00
        FORCE = 0.0D00

        ! cartesian coords
        !----------------
        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           CALL DiffuseConvectiveCompose( MASS, STIFF, FORCE, LOAD, &
                StoringCoef(1:N), C0, C1(1:N), Transmitivity, &
                .FALSE., Zero, Zero, VNull, VNull, VNull, &
                Nochange(1,1:N),Nochange(2,1:N),Nochange(3,1:N),&
                Viscosity, Density, Pressure, Zero, Zero,&
                .FALSE., Stabilize, UseBubbles, Element, N, ElementNodes )


           ! special coords (account for metric)
           !-----------------------------------
        ELSE
           CALL DiffuseConvectiveGenCompose( &
                MASS, STIFF, FORCE, LOAD, &
                StoringCoef(1:N), C0, C1(1:N), Transmitivity, &
                .FALSE., Zero, Zero, VNull, VNull, VNull, &
                Nochange(1,1:N),Nochange(2,1:N),Nochange(3,1:N), Viscosity,&
                Density, Pressure, Zero, Zero,.FALSE.,&
                Stabilize, Element, N, ElementNodes )

        END IF
        !------------------------------------------------------------------------------
        ! If time dependent simulation add mass matrix to stiff matrix
        !------------------------------------------------------------------------------

        TimeForce = FORCE
        IF ( TransientSimulation ) THEN
           IF ( UseBubbles ) FORCE = 0.0d0
           CALL Default1stOrderTime( MASS,STIFF,FORCE )
        END IF
        !------------------------------------------------------------------------------
        !  Update global matrices from local matrices
        !-----------------------------------------------------------------------------
        IF (  UseBubbles ) THEN
           CALL Condensate( N, STIFF, FORCE, TimeForce )
           IF (TransientSimulation) CALL DefaultUpdateForce( TimeForce )
        END IF
        CALL DefaultUpdateEquations( STIFF, FORCE )


        IF (DIM.EQ.2)THEN
           Wpress(WpPerm(Element % NodeIndexes(1:N))) = Density(1:N) * Gravity(1:N) &
                *(Wload(WloadPerm(Element % NodeIndexes(1:N)))-ElementNodes % y(1:N))
        ELSEIF(DIM.EQ.3)THEN
           Wpress(WpPerm(Element % NodeIndexes(1:N))) = Density(1:N) * Gravity(1:N) &
                *(Wload(WloadPerm(Element % NodeIndexes(1:N)))-ElementNodes % z(1:N))
        END IF
        
        WHERE(Wpress(WpPerm(Element % NodeIndexes(1:N))).LT.0.0)&
             Wpress(WpPerm(Element % NodeIndexes(1:N)))=0.0

     END DO     !  Bulk elements

     ! This was introduced to enable the computation of loads in the new system 
     
     CALL DefaultFinishBulkAssembly()
     
     !------------------------------------------------------------------------------
     ! Neumann & Newton boundary conditions
     !------------------------------------------------------------------------------
     DO t=1, Solver % Mesh % NumberOfBoundaryElements

        ! get element information
        Element => GetBoundaryElement(t)

        IF ( .NOT.ActiveBoundaryElement() ) CYCLE
        n = GetElementNOFNodes()

        BC => GetBC()
        bc_id = GetBCId( Element )
        CALL GetElementNodes( ElementNodes,Element )            

        IF ( ASSOCIATED( BC ) ) THEN

           ! Check that we are on the correct boundary part!
           STIFF=0.0D00
           FORCE=0.0D00
           MASS=0.0D00
           LOAD=0.0D00
           TransferCoeff = 0.0D00
           WloadExt = 0.0D00
           FluxBC = .FALSE.
           FluxBC =  GetLogical(BC,TRIM(Solver % Variable % Name) // ' Flux BC', Found)

           IF (FluxBC) THEN

              !BC: -k@Hw/@n = a(Hw - HwExt)
              !Check it if you want to use it this would represent a flux through a leaking media
              !Checking of equations, parameters or units have not been done
              !----------------------------------------------------------------------------------

              TransferCoeff(1:n) = GetReal( BC, TRIM(Solver % Variable % Name) //  ' Transfer Coefficient', Found )
              IF ( ANY(TransferCoeff(1:n).NE.0.0d0) ) THEN
                 WloadExt(1:n) = GetReal( BC, TRIM(Solver % Variable % Name) // ' External Value', Found )   
                 DO j=1,n
                    LOAD(j) = LOAD(j) +  TransferCoeff(j) * WloadExt(j)
                 END DO

              END IF
              !---------------
              !BC: -k@T/@n = q
              !---------------
              LOAD(1:n)  = LOAD(1:n) + &
                   GetReal( BC, TRIM(Solver % Variable % Name) // ' Water Flux', Found )
              !-------------------------------------
              ! set boundary due to coordinate system
              ! -------------------------------------
              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                 CALL DiffuseConvectiveBoundary( STIFF,FORCE, &
                      LOAD,TransferCoeff,Element,n,ElementNodes )
              ELSE
                 CALL DiffuseConvectiveGenBoundary(STIFF,FORCE,&
                      LOAD,TransferCoeff,Element,n,ElementNodes ) 
              END IF
           END IF
        END IF
        !------------------------------------------------------------------------------
        ! Update global matrices from local matrices
        !------------------------------------------------------------------------------
        IF ( TransientSimulation ) THEN
           MASS = 0.d0
           CALL Default1stOrderTime( MASS, STIFF, FORCE )
        END IF

        CALL DefaultUpdateEquations( STIFF, FORCE )

     END DO ! Neumann & Newton BCs

     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()

     OldValues = SystemMatrix % Values
     OldRHS = ForceVector 

     !------------------------------------------------------------------------------
     ! Dirichlet method - matrix and force-vector manipulation 
     !------------------------------------------------------------------------------
     IF (ApplyDirichlet) THEN

        ! manipulation of the matrix
        !---------------------------   
        DO i=1,Model % Mesh % NumberofNodes
           k = WloadPerm(i)  
           IF (k > 0) THEN
              IF (ActiveSed(i)) THEN
                 CALL ZeroRow( SystemMatrix, k ) 
                 CALL SetMatrixElement( SystemMatrix, k, k, 1.0d0 )
                 SystemMatrix % RHS(k) = UpperLimit(i)
              END IF
           END IF
        END DO
     END IF

     CALL Info( TRIM(SolverName) // ' Sediment Solver', 'Assembly done', Level=4 )
     !------------------------------------------------------------------------------
     !     Solve the system
     !------------------------------------------------------------------------------

     at = CPUTime() - at
     st = CPUTime()

     PrevNorm = Solver % Variable % Norm
     Norm = DefaultSolve()

     st = CPUTime()-st
     totat = totat + at
     totst = totst + st
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',Peniter,' Assembly: (s)', at, totat
     CALL Info( SolverName, Message, Level=4 )
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',Peniter,' Solve:    (s)', st, totst
     CALL Info( SolverName, Message, Level=4 )
     
     
     IF ( PrevNorm + Norm /= 0.0d0 ) THEN
        RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
     ELSE
        RelativeChange = 0.0d0
     END IF
     
     WRITE( Message, * ) 'Result Norm   : ',Norm
     CALL Info( SolverName, Message, Level=4 )
     WRITE( Message, * ) 'Relative Change : ',RelativeChange
     CALL Info( SolverName, Message, Level=4 )

     SystemMatrix % Values = OldValues
     ForceVector = OldRHS
     !------------------------------------------------------------------------------
     ! compute residual
     !------------------------------------------------------------------------------ 
     IF ( ParEnv % PEs .GT. 1 ) THEN 

        CALL ParallelInitSolve( SystemMatrix, Wload, ForceVector, ResidualVector )
        CALL ParallelMatrixVector( SystemMatrix, Wload, StiffVector, .TRUE. )
        ResidualVector =  StiffVector - ForceVector
        CALL ParallelSumVector( SystemMatrix, ResidualVector )
     ELSE 

        CALL CRS_MatrixVectorMultiply( SystemMatrix, Wload, StiffVector)
        ResidualVector =  StiffVector - ForceVector
     END IF
     !-----------------------------
     ! determine "active" nodes set 
     !-----------------------------


     DO i=1,Model % Mesh % NumberOfNodes 
        k = HomolPerm(i)
        l = WloadPerm(i)
        IF ( (k.LE.0).OR.(l.LE.0)) CYCLE

        WloadHomologous(k) = Wload(l) - UpperLimit(i)

        IF ((ApplyDirichlet).AND.(WloadHomologous(k).GE.0.0)) THEN
           !----------------------------------------------------------
           ! if upper limit is exceeded, manipulate matrix in any case
           !----------------------------------------------------------
           ActiveSed(i) = .TRUE.

           WloadHomologous(k) = LinearTol
           !---------------------------------------------------
           ! if there is "heating", don't manipulate the matrix
           !---------------------------------------------------
           IF( (ResidualVector(l).GE.0.0) &
                .AND.( Peniter.GT.1)) ActiveSed(i) = .FALSE.
        END IF

        IF( .NOT.ActiveSed(i) ) THEN
           PointerToResidualVector(PTRPerm(i)) = 0.0D00
        ELSE
           PointerToResidualVector(PTRPerm(i)) = ResidualVector(l)
        END IF
     END DO
     !------------------------------------------
     ! special treatment for periodic boundaries
     !------------------------------------------

     k=0
     DO t=1, Solver % Mesh % NumberOfBoundaryElements

        ! get element information
        Element => GetBoundaryElement(t)
        IF ( .NOT.ActiveBoundaryElement() ) CYCLE
        n = GetElementNOFNodes()
        BC => GetBC()
        bc_id = GetBCId( Element )

        CALL GetElementNodes( ElementNodes,Element )

        IF ( ASSOCIATED( BC ) ) THEN    
           IsPeriodicBC = GetLogical(BC,'Periodic BC ' // TRIM(Solver % Variable % Name), Found)
           IF (.NOT.Found) IsPeriodicBC = .FALSE.
           IF (IsPeriodicBC) THEN 
              DO i=1,N
                 IF  (ActiveSed(Element % NodeIndexes(i))) THEN
                    k = k + 1
                    ActiveSed(Element % NodeIndexes(i)) = .FALSE.
                 END IF
              END DO
           END IF
        END IF
     END DO
     
     !Checking for convergence
     !------------------------

     IF (RelativeChange.LT.NonlinearTol) THEN
        EXIT

     END IF

     IF((PenIter.EQ.NonlinearIter).AND.(RelativeChange.GT.NonlinearTol))THEN
        Write(*,*)'NOT CONVERGED'
        STOP
     END IF
  END DO


END SUBROUTINE HydroSolver
