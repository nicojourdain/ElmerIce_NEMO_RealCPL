
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
! * Migrating to subglacial hydrology of the TemperateIceSolver.
! * Starting on the 20 of may 2008
! * Values to define are : Transmitivity of the terrain which is Hydraulic 
! *                            Conductivity multiplied by the thickness 
! *                            of the sediment
! *                        Storing Coeficient of the Aquifer which is define 
! *                            has a compression coefficient times porosity times 
! *                            thickness time gravitation acceleration
! *                        Density of water
! *                        The sink source flux
! * For the boundary conditions: External Value: value of the Water load out of the domain
! *                              Transfer coefficient at the boundary, roughly it is the 
! *                                   transmitivity over a length value that could be 
! *                                   assimilate at the thickness of the interface layer
! *                              Water Flux: boundary condition is a given flux
! *
! * Hwater =0 when the top of the aquifer is at the sea level
! *
! *
!------------------------------------------------------------------------------
! *2012/04/06 : Reprise du code a grenoble
! *2012/04/19 : Nettoyage pour prendre en compte les nouvelles fonctions d'appel
! *2012/04/20 : Nettoyage de variables, FreeEnd est squeeze pour une meilleur gestion 
! *             de la condition limite par matc directement avec OpenCEL
! *2012/04/20 : Copie de la version standard pour mode dev
! * =>  05/15 : Quelques modifs depuis la derniere version sur le calcul et reinjection 
! *             de closure flux
! *2012/05/15 : Ajout du calcul des gradient de charge (pipe et Sed)
! *2012/05/16 : Copie de dev pour nettoyage
! *2012/05/?? : Calcul des flux, calcul de la transmitivite par USF
! *2012/05/30 : Modification du calcul de la lame d'eau réinjectee dans la CEL, maintenant
! *             avec un solver externe GetArea
! *2012/06/05 : Modification du calcul du volume d'eau reinjecter en fin de CEL
! *2012/06/07 : Passage des conductivites en transmitivites
! *2012/06/07 : Calcul des flux par un solver externe
! *2012/06/12 : Changement de la prise en compte des conditions limites
! *2012/06/12 : Debut de modification de WaterTransfer
! *2012/06/29 : All Pipe stuff replaced by CEL stuff
! *2012/07/08 : Calcul externe du parametre de couplage
! *2012/08/17 : Changement des condition de calcul CEL pour parallelisation
! *2012/08/20 : Nettoyage
!------------------------------------------------------------------------------

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
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: Element, CurrentElement
  TYPE(Variable_t), POINTER :: WloadSol, SedimentLoad, CELLoad, VarWload, &
       VarWloadResidual, VarSedResidual, VarCELResidual, WaterPressure, &
       PorePressure, Piping, ResVolume, ResLayer

  TYPE(ValueList_t), POINTER :: Material, SolverParams, BodyForce, &
       Constants, BC 

  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, SolverName

  INTEGER :: i,j,k,l,m,n,t,p,iter,material_id, &
       istat, LocalNodes,bf_id, DIM, NonlinearIter, &
       TimeStepNum, time, bc_id

  INTEGER, POINTER :: NodeIndexes(:), WloadPerm(:), PpPerm(:), WpPerm(:), &
       PTRPerm(:), PTSRPerm(:), PTPRPerm(:), SedLoadPerm(:), &
       CELLoadPerm(:), PipingPerm(:), VolumePerm(:), ResLayPerm(:)

  LOGICAL :: Stabilize = .TRUE., UseBubbles = .FALSE., &
       Found,AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., &
       FirstTime=.TRUE., ApplyDirichlet, FirstTime1=.TRUE.,FirstTime2=.TRUE., &
       Draining = .FALSE., CELDrainage

  LOGICAL, ALLOCATABLE :: ActiveCEL(:), ActiveSed(:)

  REAL(KIND=dp) :: NonlinearTol, LinearTol, Relax, &
       dt, CumulativeTime, RelativeChange, &
       SedimentNorm, Norm, CELNorm, PrevNorm, PrevSedNorm, &
       PrevCELNorm,WatComp, IceBase, &
       HydroTimeStep, tps, totat, st, totst

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

  REAL(KIND=dp), POINTER :: Wload(:), &
       PointerToResidualVector(:), Ppress(:), Wpress(:), &
       PointerToSedResidual(:), PointerToCELResidual(:), &
       SedLoad(:), DrainageLoad(:), OpenCEL(:), VolumeValues(:), &
       ResLayValues(:),Hwrk(:,:,:)

  REAL(KIND=dp), ALLOCATABLE ::Density(:), Gravity(:), SedThick(:), &
       Bedrock(:), WaterDensity(:), Grav(:), g(:,:), InitPP(:), &
       PreviousSed(:,:), PreviousCEL(:,:), SedToCEL(:),&
       ClosureFlux(:), CurrentLoad(:), Transmitivity(:,:,:), &
       CELThick(:), Buf(:)


  SAVE &
       ElementNodes,    &
       Density,         &
       ActiveCEL,       &
       ActiveSed,       &
       AllocationsDone, &
       FirstTime,       &
       FirstTime1,      &
       FirstTime2,      &
       VariableName,    &
       SolverName,      &
       NonLinearTol,    & 
       M,               &
       InitPP,          &
       Gravity,         &
       SedThick,        &
       Bedrock,         &
       WaterDensity,    &
       g,               &
       Grav,            &
       PreviousSed,     &
       PreviousCEL,     &
       SedToCEL,        &     
       ClosureFlux,     &
       CurrentLoad,     &
       Hwrk,            &
       Transmitivity,   &
       CELThick,        &
       Buf

  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  SolverName = 'HydroSolver ('// TRIM(Solver % Variable % Name) // ')'
  VariableName = TRIM(Solver % Variable % Name)

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  PointerToSolver => Solver

  WloadSol => Solver % Variable
  WloadPerm  => WloadSol % Perm
  Wload => WloadSol % Values
  LocalNodes = COUNT( WloadPerm > 0 )
  IF ( LocalNodes <= 0 ) RETURN


  !------------------------------------------------------------------------------
  !    Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     N = Solver % Mesh % MaxElementNodes
     M = LocalNodes
     IF ( AllocationsDone ) THEN

        DEALLOCATE(            &
             ElementNodes % x, &
             ElementNodes % y, &
             ElementNodes % z, &
             Density,          &
             ActiveCEL,        &
             ActiveSed,        &
             Gravity,          &
             InitPP,           &
             SedThick,         &
             Bedrock,          &
             WaterDensity,     &
             Grav,             &
             g,                &
             PreviousSed,      &
             PreviousCEL,      &
             SedToCEL,         &
             ClosureFlux,      &
             CurrentLoad,      &
             Transmitivity,    &
             CELThick,         &
             Buf)              

     END IF

     ALLOCATE(                   &
          ElementNodes % x( N ), &
          ElementNodes % y( N ), &
          ElementNodes % z( N ), &
          Density( N ),          &
          ActiveCEL( M ),        &
          ActiveSed( M ),        &
          Gravity( N ),          &
          InitPP( M ),           &
          SedThick( N ),         &
          Bedrock( M ),          &
          WaterDensity( M ),     &
          Grav( M ),             &
          g( 3,N ),              &
          PreviousSed( M,1 ),    &
          PreviousCEL( M,1 ),    &
          SedToCEL( M ),         &
          ClosureFlux( M ),      &
          CurrentLoad( M ),      &
          Transmitivity(3,3,N),  &
          CELThick( N ),         &
          Buf( N ),              &
          STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error 1' )
     ELSE
        CALL INFO(SolverName, 'Memory allocation done 1', level=1 )
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
  !    Read physical and numerical constants and initialize 
  !------------------------------------------------------------------------------
  Constants => GetConstants()
  SolverParams => GetSolverParams()

  Stabilize = GetLogical( SolverParams,'Stabilize',Found )
  IF (.NOT. Found) Stabilize = .FALSE.
  UseBubbles = GetLogical( SolverParams,'Bubbles',Found )
  IF ( .NOT.Found .AND. (.NOT.Stabilize)) UseBubbles = .TRUE.

  LinearTol = GetConstReal( SolverParams, &
       'Linear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) THEN
     CALL FATAL(SolverName, 'No >Linear System Convergence Tolerance< found')
  END IF


  NonlinearIter = GetInteger(   SolverParams, &
       'Nonlinear System Max Iterations', Found )
  IF ( .NOT.Found ) NonlinearIter = 1

  NonlinearTol = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) NonlinearTol = 1e-6

!!$  Relax = GetConstReal( SolverParams, &
!!$       'Nonlinear System Relaxation Factor',Found )
!!$  IF ( .NOT.Found ) Relax = 1.0D00

  ApplyDirichlet = GetLogical( SolverParams, &
       'Apply Dirichlet', Found)
  IF ( .NOT.Found ) THEN
     ApplyDirichlet = .FALSE.
  END IF

  WatComp = GetConstReal( Constants, &
       'Water Compressibility', Found)
  IF ( .NOT.Found ) THEN
     WatComp = 5.04e-4
  END IF

  CELDrainage = GetLogical( SolverParams, &
       'Apply CEL Drainage', Found)
  IF ( .NOT.Found ) THEN
     CELDrainage = .FALSE.
  END IF

  Relax = 0.7
  CumulativeTime = 0.0d0
  ActiveCEL = .FALSE.
  ActiveSed = .FALSE.
  ClosureFlux = 0.0
  FirstTime = .TRUE.

  !------------------------------------------------------------------------------
  !       The first time around this has been done by the caller...
  !------------------------------------------------------------------------------

  IF ( TransientSimulation .AND.(.NOT.FirstTime)) THEN
     CALL InitializeTimestep( Solver )
  END IF

  FirstTime = .FALSE.

  totat = 0.0d0
  totst = 0.0d0
  !------------------------------------------------------------------------------
  !       Get externally declared DOFs
  !------------------------------------------------------------------------------

  IF (.NOT.ApplyDirichlet) THEN
     ActiveCEL = .FALSE.
     ActiveSed = .FALSE.
  END IF

  !Variables needed to compute the pressure
  !----------------------------------------

  WaterPressure => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Pressure' )
  IF (.NOT.ASSOCIATED(WaterPressure)) THEN
     WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' WaterPressure not associated'
     CALL FATAL( SolverName, Message)
  END IF
  Wpress => WaterPressure % Values
  WpPerm  => WaterPressure % Perm

  PorePressure => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' PorePressure' )
  IF (.NOT.ASSOCIATED(PorePressure)) THEN
     WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' PorePressure not associated'
     CALL FATAL( SolverName, Message)
  END IF
  Ppress => PorePressure % Values
  PpPerm  => PorePressure % Perm

  IF (FirstTime1) THEN
     InitPP = PorePressure % Values
     FirstTime1 = .FALSE.
  END IF

  VarWload => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Homologous' )
  IF (.NOT.ASSOCIATED(VarWload)) THEN
     WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' VarWload not associated'
     CALL FATAL( SolverName, Message)
  END IF

  !Dummy array for the residual
  !----------------------------

  VarWloadResidual => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Residual' )
  IF (.NOT.ASSOCIATED(VarWloadResidual)) THEN
     WRITE(Message,'(A)') '>' // TRIM(Solver % Variable % Name) // ' Residual< not associated'
     CALL FATAL( SolverName, Message)
  END IF
  PointerToResidualVector => VarWloadResidual % Values
  PTRPerm => VarWloadResidual % Perm

  !Residual for the sediment layer
  !-------------------------------

  VarSedResidual => VariableGet( Model % Mesh % Variables,'Sediment Residual' )
  IF (.NOT.ASSOCIATED(VarSedResidual)) THEN
     WRITE(Message,'(A)') '>Sediment Residual< not associated'
     CALL FATAL( SolverName, Message)
  END IF
  PointerToSedResidual => VarSedResidual % Values
  PTSRPerm => VarSedResidual % Perm

  !Residual for the Draining layer
  !-------------------------------

  VarCELResidual => VariableGet( Model % Mesh % Variables,'CEL Residual' )
  IF (.NOT.ASSOCIATED(VarCELResidual)) THEN
     WRITE(Message,'(A)') '>CEL Residual< not associated'
     CALL FATAL( SolverName, Message)
  END IF
  PointerToCELResidual => VarCELResidual % Values
  PTPRPerm => VarCELResidual % Perm

  !Residual needed to convert water load from penalisation method to volumes
  !-------------------------------------------------------------------------

  ResLayer => VariableGet ( Solver % Mesh % Variables, 'Residual Layer' )
  IF ( ASSOCIATED( ResLayer ) ) THEN
     ResLayPerm    => ResLayer % Perm
     ResLayValues  => ResLayer % Values
  ELSE
     CALL Info('VolumeToLayer', &
          & 'No variable for Layer associated.', Level=4)
  END IF

  !Water load in the Sediment
  !--------------------------

  SedimentLoad => VariableGet( Model % Mesh % Variables,'Sediment Wload' )
  IF (.NOT.ASSOCIATED(SedimentLoad)) THEN
     WRITE(Message,'(A)') '>Sediment Wload< not associated'
     CALL FATAL( SolverName, Message)
  END IF
  SedloadPerm  => SedimentLoad % Perm
  Sedload => SedimentLoad % Values

  !Water load in the CELs
  !-----------------------

  CELLoad => VariableGet( Model % Mesh % Variables,'CEL Wload' )
  IF (.NOT.ASSOCIATED(CELLoad)) THEN
     WRITE(Message,'(A)') '>CEL Load< not associated'
     CALL FATAL( SolverName, Message)
  END IF
  CELLoadPerm  => CELLoad % Perm
  DrainageLoad => CELLoad % Values

  !Is the pipe open
  !----------------

  Piping => VariableGet( Model % Mesh % Variables,'Open CEL' )
  IF (.NOT.ASSOCIATED(Piping)) THEN
     WRITE(Message,'(A)') '>Open CEL< not associated'
     CALL FATAL( SolverName, Message)
  END IF
  PipingPerm  => Piping % Perm
  OpenCEL => Piping % Values

  !Initialisation of the CEL Boundary
  !------------------------------------

  DO t=1, Solver % Mesh % NumberOfBoundaryElements

     ! get element information
     Element => GetBoundaryElement(t)
     n = GetElementNOFNodes()
     IF (ALL(PipingPerm(Element % NodeIndexes(1:N)).LE.0))CYCLE
     BC => GetBC()

     CALL GetElementNodes( ElementNodes,Element )

     !Open CEL =3 for an open channel and =0 for a closed one

     IF ( ASSOCIATED( BC ) ) THEN 
        bc_id = GetBCId( Element ) 
        DO i=1,N
           IF(PipingPerm(Element % NodeIndexes(i)).LE.0)CYCLE

           Buf(1:N) = ListGetReal( BC, 'Open CEL',n ,Element % NodeIndexes, Found )

           IF(Found)OpenCEL(PipingPerm(Element % NodeIndexes(i))) &
                = MAX(Buf(i),OpenCEL(PipingPerm(Element % NodeIndexes(i))))

        END DO
     END IF

  END DO

  !Saving values of the previous timestep
  !---------------------------------------

  PreviousSed(:,1) = SedLoad
  PreviousCEL(:,1) = DrainageLoad 

  !initiation of closure flux
  ClosureFlux = 0.0
  
  !------------------------------------------------------------------------------
  !       non-linear system iteration loop
  !------------------------------------------------------------------------------
  DO iter=1,NonlinearIter

     !------------------------------------------------------------------------------
     ! print out some information
     !------------------------------------------------------------------------------
     at  = CPUTime()
     at0 = RealTime()
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     WRITE( Message,'(A,A,I3,A,I3)') &
          TRIM(Solver % Variable % Name),  ' iteration no.', iter,' of ',NonlinearIter
     CALL Info( SolverName, Message, Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, 'Starting Assembly...', Level=4 )


     at = CPUTime() - at
     st = CPUTime()

     IF( Draining) THEN
        Solver % Variable % Norm = SedimentNorm + CELNorm
     ELSE
        Solver % Variable % Norm = SedimentNorm
     END IF

     PrevSedNorm = SedimentNorm
     PrevCELNorm = CELNorm
     PrevNorm = Solver % Variable % Norm

     !Computing the amount of water that should pass from the sediment to the CEL
     !----------------------------------------------------------------------------

     IF (Draining.AND.(iter.GT.1))THEN
        
        CALL WaterTransfer(Model, Solver, SedLoad, DrainageLoad, &
             Piping, SedToCEL, Solvername, LocalNodes, Relax)
        
     END IF

     !Solving the Hydro equation for the sediment layer
     !-------------------------------------------------

     CALL SedimentSolve(Model, Solver, VarWload, TransientSimulation, &
          Stabilize, ApplyDirichlet, UseBubbles, Draining, ActiveSed, Piping, LinearTol, &
          WloadSol, WatComp, SedimentNorm, InitPP, VarWloadResidual, VarSedResidual, &
          VarCELResidual, Ppress,  SedimentLoad, VariableName, SolverName, PreviousSed, SedToCEL, &
          ClosureFlux, iter)


     !If necessary, Solving the hydro Equation for the CEL equivalent Layer
     !-Necessary when sediment load is at floatation or any of the node of
     !the efficient drainage is activated
     !-When computation is not done, water load in the pipe is fixed at the 
     !lowest point of the bed
     !----------------------------------------------------------------------
     
     !Sequential version (safe)
!     IF ((ANY(PointerToSedResidual.LT.0.0).AND.(CELDrainage)) & 
!          .OR.((ANY((OpenCEL.GT.0).AND.(OpenCEL.LT.3)).AND.(CELDrainage))))THEN 
     
        !parallel version (to check)

     IF (((ParallelReduction(MINVAL(PointerToSedResidual)).LT.0.0).AND.(CELDrainage)) & 
          .OR.((ParallelReduction(SUM(OpenCEL)).GT.0.0).AND.(ParallelReduction(MAXVAL(OpenCEL)).LT.3.0).AND.(CELDrainage))& 
          .OR.((ParallelReduction(SUM(OpenCEL)).GT.3.0)).AND.(CELDrainage))THEN  

        !PrevSedNorm = ParallelReduction(PrevSedNorm)
        !SedimentNorm = ParallelReduction(SedimentNorm)
       
write(126,*)Parenv % mype,PrevSedNorm,ParallelReduction(PrevSedNorm)

        IF(2.0d0*ABS(PrevSedNorm-SedimentNorm)/(PrevSedNorm+SedimentNorm)&
             .LT.NonlinearTol)THEN

           Draining = .TRUE.
           CALL CELSolve(Model, Solver, VarWload, TransientSimulation, Stabilize, &
                ApplyDirichlet, UseBubbles, ActiveCEL, Piping, LinearTol, WloadSol, &
                WatComp, CELNorm, VarWloadResidual, VarSedResidual, VarCELResidual, CELLoad, &
                Ppress, SolverName, PreviousCEL, SedToCEL, ClosureFlux, CurrentLoad, &
                ResLayer, iter)
           
        END IF
     ELSE 
        IF(DIM.EQ.3)THEN
           DrainageLoad  = ParallelReduction(minval(Model % Nodes % z (:)))
        ELSEIF(DIM.EQ.2)THEN
           DrainageLoad  = ParallelReduction(minval(Model % Nodes % y (:)))
        END IF
     END IF

     !If the pipe is not open everywhere compute the flow at the end
     !--------------------------------------------------------------

     IF((ParallelReduction(SUM(OpenCEL)).GT.1.0).AND.(ParallelReduction(MINVAL(OpenCEL)).EQ.0.0)&
          .AND.(CELDrainage))THEN

        CALL EndOfCEL (Model, SedLoad, DrainageLoad, WloadPerm, Piping, &
             ClosureFlux, CurrentLoad, LocalNodes)
     ELSE

        ClosureFlux = 0.0
     END IF

     !Checking for convergence
     !------------------------

     IF(Draining) THEN
        Norm = SedimentNorm + CELNorm
     ELSE
        Norm = SedimentNorm
     END IF

     st = CPUTime()-st
     totat = totat + at
     totst = totst + st
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly tot: (s)', at, totat
     CALL Info( SolverName, Message, Level=4 )
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
     CALL Info( SolverName, Message, Level=4 )


     IF ((PrevNorm + Norm).NE.0.0d0 ) THEN
        RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
     ELSE
        RelativeChange = 0.0d0
     END IF

     WRITE( Message, * ) 'Result Norm   : ',Norm
     CALL Info( SolverName, Message, Level=4 )
     WRITE( Message, * ) 'Relative Change : ',RelativeChange
     CALL Info( SolverName, Message, Level=4 )

   WRITE(Message,'(a,e13.6,a,e13.6)') &
         'parared Max/min values of sediment Water load:', ParallelReduction(MAXVAL(Wload(:)),2), &
         '/',ParallelReduction(MAXVAL(Wload(:)),1)
    CALL INFO(SolverName,Message,Level=4)
    WRITE(Message,'(a,e13.6,a,e13.6)') &
         'Max/min values of sedimentload:', ParallelReduction(MAXVAL( SedLoad(:))), &
         '/',ParallelReduction(MINVAL( SedLoad(:)))
    CALL INFO(SolverName,Message,Level=4)

     IF (RelativeChange.LT.NonlinearTol) THEN
        EXIT
     ELSE
        IF (ApplyDirichlet) THEN
           WRITE(Message,'(a,i10)') 'Deactivated Periodic BC nodes:', k
           CALL INFO(SolverName,Message,Level=1)
           WRITE(Message,'(a,i10)') 'Number of constrained points for the sediment:', COUNT(ActiveSed)
           CALL INFO(SolverName,Message,Level=1)

        END IF

     END IF

     !Brutal stop if ever the solver does not converge
     !------------------------------------------------

     IF((iter.EQ.NonlinearIter).AND.(RelativeChange.GT.NonlinearTol))THEN
        Write(*,*)'NOT CONVERGED'
        STOP
     END IF

  END DO ! of the nonlinear iteration



  !------------------------------------------------------------------------------
  !Compute the water pressure at the base of the ice (WaterPressure)
  !and the pressure at the base of the sediment (PorePressure)
  !------------------------------------------------------------------------------

  IF (FirstTime2) THEN
     FirstTime2 = .FALSE.

     DO t=1,Solver % NumberOfActiveElements
        CurrentElement => GetActiveElement(t) 
        Material => GetMaterial(CurrentElement)
        n = GetElementNOFNodes() 

        !get the density 
        !---------------

        Density(1:n) = ListGetReal( Material, 'Water Density',  n, Element % NodeIndexes, Found )
        IF (.NOT.Found) THEN
           Density = 1000.0D00
           WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Water Density< not found for element ',&
                t, ' material,Set to 1000 kg/m3 ', material_id
           CALL INFO(SolverName,Message,Level=4)
        END IF

        SedThick(1:n) = listGetReal( Material,'Sediment Thickness', n, Element % NodeIndexes, Found )
        IF (.NOT.Found) THEN
           SedThick = 10.0D00
           WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Thickness', &
                '< not found for element ', t, ' material ', material_id
           CALL INFO(SolverName,Message,Level=4)
        END IF

        !Get gravity
        !-----------

        BodyForce => GetBodyForce()    
        IF ( ASSOCIATED( BodyForce ) ) THEN
           bf_id = GetBodyForceId()
           g = 0.0_dp  
           g(1,1:n) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 1',Found)
           g(2,1:n) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Found)
           g(3,1:n) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Found)
           Gravity(1:n) = SQRT(SUM(g(:,1:n)**2.0/n))
        END IF

        DO j=1,n

           k = CurrentElement % NodeIndexes(j)
           WaterDensity(WloadPerm(k)) = Density(j)
           Grav(WloadPerm(k)) = Gravity(j)

           !On the first time saving the position of the bedrock
           !----------------------------------------------------
           IF(DIM.EQ.3)THEN
              IceBase = Model%Nodes%z(k)
           ELSEIF (DIM.EQ.2)THEN
              IceBase = Model%Nodes%y(k)
           END IF

           !That is introduced to take into account a lower free surface not included 
           !so far and to check
           !-------------------------------------------------------------------------
           Bedrock(WloadPerm(k))=IceBase

           !Wpress is the pressure at the base of the ice, computed from the Water load
           !It's the height of the water column other the bed multiplied by G and rho wat
           !-----------------------------------------------------------------------------
           Wpress(WpPerm(k)) = (WaterDensity(WloadPerm(k)) * Grav(WloadPerm(k))) &
                *(SedLoad(SedLoadPerm(k))-IceBase) 

           !Ppress is the pressure at the base of the Sediment, computed from the Water load
           !The thickness of the sediment layer is added under the bedrock
           !--------------------------------------------------------------------------------
           Ppress(PpPerm(k)) =  (WaterDensity(WloadPerm(k)) * Grav(WloadPerm(k))) &
                *(SedLoad(SedLoadPerm(k))-Bedrock(WloadPerm(k))-SedThick(j))
        END DO
     END DO

  ELSE

     DO t=1,Solver % NumberOfActiveElements

        CurrentElement => GetActiveElement(t) 
        Material => GetMaterial(CurrentElement)
        n = GetElementNOFNodes()

        SedThick(1:n) = listGetReal( Material,'Sediment Thickness', n, Element % NodeIndexes, Found )
        IF (.NOT.Found) THEN
           SedThick = 10.0D00
           WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Thickness', &
                '< not found for element ', t, ' material ', material_id
           CALL INFO(SolverName,Message,Level=4)
        END IF

        DO j=1,n

           k = CurrentElement % NodeIndexes(j)

           IF(DIM.EQ.3)THEN
              IceBase = Model%Nodes%z(k)
           ELSEIF (DIM.EQ.2)THEN
              IceBase = Model%Nodes%y(k)
           END IF

           Wpress(WpPerm(k)) = (WaterDensity(WloadPerm(k)) * Grav(WloadPerm(k))) &
                *(SedLoad(SedLoadPerm(k))-IceBase)

           Ppress(PpPerm(k)) = (WaterDensity(WloadPerm(k)) * Grav(WloadPerm(k))) &
                *(SedLoad(SedLoadPerm(k))-Bedrock(WloadPerm(k))-SedThick(j))

        END DO
     END DO
  END IF

  WHERE (Wpress(:).LT.0.)
     Wpress(:) = 0.
  END WHERE

  WHERE (Ppress(:).LT.0.)
     Ppress(:) = 0.
  END WHERE

  SubroutineVisited = .TRUE.
  !------------------------------------------------------------------------------


CONTAINS

  SUBROUTINE SedimentSolve(Model, Solver, VarWload, TransientSimulation, &
       Stabilize, ApplyDirichlet, UseBubbles, Draining, ActiveSed, Piping, LinearTol, &
       WloadSol, WatComp, SedimentNorm, InitPP, VarWloadResidual, VarSedResidual, &
       VarCELResidual, Ppress,  SedimentLoad, VariableName, SolverName, PreviousSed, SedToCEL, &
       ClosureFlux, iter)

    !-------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Variable_t), POINTER ::  VarWload, Piping, WloadSol, VarWloadResidual, &
         VarSedResidual, VarCELResidual, SedimentLoad 
    LOGICAL :: TransientSimulation, Stabilize, ApplyDirichlet, &
         UseBubbles, Draining
    LOGICAL, ALLOCATABLE :: ActiveSed(:)
    REAL(KIND=dp) :: LinearTol, WatComp, SedimentNorm
    REAL(KIND=dp), POINTER ::Ppress(:)
    REAL(KIND=dp), ALLOCATABLE :: InitPP(:), PreviousSed(:,:), SedToCEL(:), ClosureFlux(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, SolverName
    INTEGER :: iter
    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------
    TYPE(Nodes_t) :: ElementNodes  
    TYPE(Matrix_t), POINTER :: Systemmatrix
    TYPE(Element_t),POINTER :: Element
    TYPE(ValueList_t), POINTER :: Equation,Material, BC,BodyForce

    INTEGER :: i,j,k,l,m,n,t,p, body_id, material_id, &
         bf_id, bc_id, istat

    INTEGER, POINTER :: NodeIndexes(:)

    LOGICAL ::Found, FluxBC, IsPeriodicBC=.FALSE.,&
         AllocationsDone = .FALSE.

    REAL(KIND=dp), POINTER :: WloadHomologous(:), Hwrk(:,:,:), &
         ResidualVector(:), ForceVector(:)

    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), &
         Transmitivity(:,:,:), FORCE(:),TimeForce(:), &
         TransferCoeff(:), Work(:), C1(:), C0(:), Zero(:), Viscosity(:),&
         UpperLimit(:), StoringCoef(:), Density(:), WloadExt(:), &
         StiffVector(:), OldValues(:), OldRHS(:), VNull(:), Nochange(:,:), & 
         Sv(:,:,:), Gravity(:), SedComp(:), SedThick(:) , &
         Porosity(:), Pressure(:), g(:,:), Basis(:)

    REAL(KIND=dp) :: Dist, NOpenLim

    SAVE                        &
         OldValues,             &
         OldRHS,                &
         VNull,                 &
         Nochange,              &
         Pressure,              &
         ElementNodes,          &
         Work,Zero,             &
         Viscosity,             &
         StoringCoef,           &
         Density,               &
         WloadExt,              &
         C1,                    &
         C0,                    &
         TransferCoeff,         &
         Transmitivity,         &
         Sv,                    &
         MASS,                  &
         STIFF,LOAD,            &
         FORCE,                 &
         TimeForce,             &
         StiffVector,           &
         ResidualVector,        &
         UpperLimit,            &
         AllocationsDone,       &
         Hwrk,                  &
         Porosity,              &
         Gravity,               &
         SedComp,               &
         SedThick,              &
         g,                     &
         Basis


    SystemMatrix => Solver % Matrix
    ForceVector => Solver % Matrix % RHS


    IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
       N = Solver % Mesh % MaxElementNodes
       M = LocalNodes
       K = SIZE( SystemMatrix % Values )
       L = SIZE( SystemMatrix % RHS )

       IF ( AllocationsDone ) THEN
          DEALLOCATE(                    &
               OldValues,                &
               OldRHS,                   &
               VNull,                    &
               Nochange,                 &
               Pressure,                 &
               ElementNodes % x,         &
               ElementNodes % y,         &
               ElementNodes % z,         &
               Work,Zero,                &
               Viscosity,                &
               StoringCoef,              &
               Density,                  &
               WloadExt,                 &
               C1,                       &
               C0,                       &
               TransferCoeff,            &
               Transmitivity,            &
               Sv,                       &
               MASS,                     &
               STIFF,                    &
               LOAD,                     &
               FORCE,                    &
               TimeForce,                &
               StiffVector,              &
               ResidualVector,           &
               UpperLimit,               &
               Porosity,                 &
               Gravity,                  &
               SedComp,                  &
               SedThick,                 &
               g,                        &
               Basis)              
       END IF

       ALLOCATE(                                  &
            OldValues( K ),                       &
            OldRHS( L ),                          &
            VNull( N ),                           &
            Nochange( 3,N ),                      &
            Pressure( N ),                        &
            ElementNodes % x( N ),                &
            ElementNodes % y( N ),                &
            ElementNodes % z( N ),                &
            Work( N ), Zero( N ),                 &
            Viscosity( N ),                       &
            StoringCoef( N ),                     &
            Density( N ),                         &
            WloadExt( N ),                        &
            C1( N ),                              &
            C0( N ),                              &
            TransferCoeff( N ),                   &
            Transmitivity( 3,3,N ),               &
            Sv( 3,3,N ),                          &
            MASS(  2*N,2*N ),                     &
            STIFF( 2*N,2*N ),                     &
            LOAD( N ),                            &
            FORCE( 2*N ),                         &
            TimeForce( 2*N ),                     &
            StiffVector( L ),                     &
            ResidualVector(L),                    &
            UpperLimit( M ),                      &
            Porosity( N ),                        &
            Gravity( N ),                         &
            SedComp( N ),                         &
            SedThick( N ),                        &
            g( 3,N ),                             &
            Basis( 2*N ),                         &
            STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL FATAL( SolverName, 'Memory allocation error 1bis' )
       ELSE
          CALL INFO(SolverName, 'Memory allocation done here 1bis', level=1 )
       END IF
       AllocationsDone = .TRUE.

    END IF

    !---------------------------------------
    !Coppying the reality to the dummy
    !---------------------------------------

    IF ( TransientSimulation ) THEN
       WloadSol % PrevValues = PreviousSed
    END IF

    WloadPerm = SedLoadPerm
    Wload = SedLoad

    PointerToResidualVector = PointerToSedResidual
    PTRPerm = PTSRPerm
    !------------------------------------------------------------------------------
    ! lets start
    !------------------------------------------------------------------------------
    CALL DefaultInitialize()
    !-----------------------------------------------------------------------------
    ! Get Upper limit:
    !-----------------------------------------------------------------------------
    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
       n = GetElementNOFNodes()
       CALL GetElementNodes( ElementNodes,Element )
       Material => GetMaterial()

       UpperLimit(WloadPerm(Element % Nodeindexes(1:n))) = ListGetReal(Material,TRIM(VariableName) // & 
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
         'Max/min values of sediment Water load:', ParallelReduction(MAXVAL( Wload(:))), &
         '/',ParallelReduction(MINVAL( Wload(:)))
    CALL INFO(SolverName,Message,Level=4)
    WRITE(Message,'(a,e13.6,a,e13.6)') &
         'Max/min values of sedimentload:', ParallelReduction(MAXVAL( SedLoad(:))), &
         '/',ParallelReduction(MINVAL( SedLoad(:)))
    CALL INFO(SolverName,Message,Level=4)
    !------------------------------------------------------------------------------

    body_id = -1
    NULLIFY(Material)

    !------------------------------------------------------------------------------
    ! Bulk elements
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

       k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
            minv=1, maxv=Model % NumberOFEquations )

       !------------------------------------------------------------------------------
       ! Get element material parameters
       !------------------------------------------------------------------------------              

       N = GetElementNOFNodes(Element)
       CALL GetElementNodes( ElementNodes,Element )

       SedComp(1:N) = listGetReal( Material,'Sediment Compressibility', N, Element % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          SedComp = 1.0D-2
          WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Compressibility', &
               '< not found for element ', t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       Porosity(1:N) = listGetReal( Material,'Sediment Porosity', N, Element % NodeIndexes, Found )

       IF (.NOT.Found) THEN
          Porosity = 0.4D00
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
          SedThick = 10.0D00
          WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Thickness', &
               '< not found for element ', t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       Density(1:N) = ListGetReal( Material, 'Water Density',  N, Element % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          Density = 1.0055e-18
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
       ! add water contribution from sif file, closure and transfer routines
       !------------------------------------------------------------------------------
       LOAD = 0.0D00

       IF (Draining)THEN
          DO i=1,N
             LOAD(i) = ClosureFlux(WloadPerm(Element % NodeIndexes(i)))&
                  - SedToCEL(WloadPerm(Element % NodeIndexes(i)))
          END DO
       ELSE
          LOAD = 0.0
       END IF

       BodyForce => GetBodyForce()
       IF (ASSOCIATED(BodyForce)) THEN
          bf_id = GetBodyForceId()
          LOAD(1:N) = LOAD(1:N) +   &
               ListGetReal( BodyForce, &
               TRIM(Solver % Variable % Name) // ' Source flux', &
               N, Element % NodeIndexes(1:N),Found )
       END IF


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
               StoringCoef, C0, C1(1:N), Transmitivity, &
               .FALSE., Zero, Zero, VNull, VNull, VNull, &
               Nochange(1,1:N),Nochange(2,1:N),Nochange(3,1:N),&
               Viscosity, Density, Pressure, Zero, Zero,&
               .FALSE., Stabilize, UseBubbles, Element, N, ElementNodes )


          ! special coords (account for metric)
          !-----------------------------------
       ELSE
          CALL DiffuseConvectiveGenCompose( &
               MASS, STIFF, FORCE, LOAD, &
               StoringCoef, C0, C1(1:N), Transmitivity, &
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
       !------------------------------------------------------------------------------

       IF (  UseBubbles ) THEN
          CALL Condensate( N, STIFF, FORCE, TimeForce )
          IF (TransientSimulation) CALL DefaultUpdateForce( TimeForce )
       END IF
       CALL DefaultUpdateEquations( STIFF, FORCE )

    END DO     !  Bulk elements

    !------------------------------------------------------------------------------
    ! Neumann & Newton boundary conditions
    !------------------------------------------------------------------------------
    DO t=1, Solver % Mesh % NumberOfBoundaryElements

       ! get element information
       Element => GetBoundaryElement(t)
       bc_id = GetBCId( Element )

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

             TransferCoeff(1:n) = GetReal( BC, TRIM(Solver % Variable % Name) //  ' Transfer Coefficient',Found )
             IF ( ANY(TransferCoeff(1:n).NE.0.0d0) ) THEN
                WloadExt(1:n) = GetReal( BC, TRIM(Solver % Variable % Name) // ' External Value',Found )   
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

    END DO   ! Neumann & Newton BCs

    !------------------------------------------------------------------------------

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
             IF (ActiveSed(k)) THEN
                CALL ZeroRow( SystemMatrix, k ) 
                CALL SetMatrixElement( SystemMatrix, k, k, 1.0d0 )
                SystemMatrix % RHS(k) = UpperLimit(k)
             END IF
          END IF
       END DO
    END IF

    CALL Info( TRIM(SolverName) // ' Sediment layer', 'Assembly done', Level=4 )


    !------------------------------------------------------------------------------
    !     Solve the system
    !------------------------------------------------------------------------------

    SedimentNorm = Solver % Variable % Norm
    SedimentNorm = DefaultSolve()

    SystemMatrix % Values = OldValues
    ForceVector = OldRHS
    !------------------------------------------------------------------------------
    ! compute residual
    !------------------------------------------------------------------------------ 
    IF ( ParEnv % PEs .GT. 1 ) THEN !!!!!!!!!!!!!!!!!!!!!! we have a parallel run
       CALL ParallelInitSolve( SystemMatrix, Wload, ForceVector, ResidualVector )
       CALL ParallelMatrixVector( SystemMatrix, Wload, StiffVector, .TRUE. )
       ResidualVector =  StiffVector - ForceVector
       CALL ParallelSumVector( SystemMatrix, ResidualVector )
    ELSE !!!!!!!!!!!!!!!!!!!!!! serial run 

       CALL CRS_MatrixVectorMultiply( SystemMatrix, Wload, StiffVector)
       ResidualVector =  StiffVector - ForceVector

    END IF
    !-----------------------------
    ! determine "active" nodes set 
    !-----------------------------

!!$    IF (ASSOCIATED(VarWload)) THEN
!!$       WloadHomologous => VarWload % Values
!!$       DO i=1,Model % Mesh % NumberOfNodes ! <______________IS THIS OK IN PARALLEL????????
!!$          k = VarWload % Perm(i)
!!$          l = WloadPerm(i)
!!$          IF ( (k.LE.0).OR.(l.LE.0)) CYCLE
!!$
!!$          WloadHomologous(k) = Wload(l) - UpperLimit(l)
!!$
!!$          IF (ApplyDirichlet) THEN
!!$             !----------------------------------------------------------
!!$             ! if upper limit is exceeded, manipulate matrix in any case
!!$             !----------------------------------------------------------
!!$             IF (WloadHomologous(k).GE.0.0 ) THEN
!!$                ActiveSed(l) = .TRUE.
!!$
!!$                WloadHomologous(k) = LinearTol
!!$
!!$             END IF
!!$             !---------------------------------------------------
!!$             ! if there is "heating", don't manipulate the matrix
!!$             !---------------------------------------------------
!!$             IF( (ResidualVector(l).GT.- LinearTol) &
!!$                  .AND.( iter.GT.1)) ActiveSed(l) = .FALSE.
!!$          END IF
!!$          IF( .NOT.ActiveSed(l) ) THEN
!!$             PointerToResidualVector(PTRPerm(i)) = 0.0D00
!!$          ELSE
!!$             PointerToResidualVector(PTRPerm(i)) = ResidualVector(l)
!!$
!!$             IF (OpenL(PipingPerm(i)).LT.2.0) OpenCEL(PipingPerm(i)) = 1.0
!!$          END IF
!!$          write(*,*)'insed', i, PointerToResidualVector(PTRPerm(i))
!!$       END DO
!!$    ELSE
!!$       WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' Homologous not associated'
!!$       CALL WARN( SolverName, Message)
!!$    END IF


   IF (ASSOCIATED(VarWload)) THEN
       WloadHomologous => VarWload % Values
       DO t=1,Solver % NumberOfActiveElements
          Element => GetActiveElement(t,Solver)
          
          IF (.NOT.ASSOCIATED(Element)) CYCLE
          N = GetElementNOFNodes(Element)
          CALL GetElementNodes( ElementNodes,Element )

          DO i=1,N
             k = VarWload % Perm(Element % NodeIndexes(i))
             l = WloadPerm(Element % NodeIndexes(i))
             IF ( (k.LE.0).OR.(l.LE.0)) CYCLE
             
             WloadHomologous(k) = Wload(l) - UpperLimit(l)
             
             IF (ApplyDirichlet) THEN
                !----------------------------------------------------------
                ! if upper limit is exceeded, manipulate matrix in any case
                !----------------------------------------------------------
                IF (WloadHomologous(k).GE.0.0 ) THEN
                   ActiveSed(l) = .TRUE.
                   
                   WloadHomologous(k) = LinearTol
                   
                END IF
                !---------------------------------------------------
                ! if there is "heating", don't manipulate the matrix
                !---------------------------------------------------
                IF((ResidualVector(l).GT.-LinearTol) &
                     .AND.( iter.GT.1)) ActiveSed(l) = .FALSE.
             END IF
             IF( .NOT.ActiveSed(l) ) THEN
                PointerToResidualVector(PTRPerm(Element % NodeIndexes(i))) = 0.0D00
             ELSE
                PointerToResidualVector(PTRPerm(Element % NodeIndexes(i))) = ResidualVector(l)
                
                IF (OpenCEL(PipingPerm(i)).LT.2.0) OpenCEL(PipingPerm(Element % NodeIndexes(i))) = 1.0
             END IF
          END DO
       END DO
    ELSE
       WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' Homologous not associated'
       CALL WARN( SolverName, Message)
    END IF


    !----------------------------------------------------
    ! Activating CEL if N is less than a prescribed value 120504
    !----------------------------------------------------
    !That is The N Upper limit which should be fix in an other way
    NOpenLim=0.2

    DO t=1,Solver % NumberOfActiveElements

       Element => GetActiveElement(t,Solver)
       IF (.NOT.ASSOCIATED(Element)) CYCLE
       IF ( Element % BodyId /= body_id ) THEN
          N = GetElementNOFNodes(Element)
          CALL GetElementNodes( ElementNodes,Element )

          BodyForce => GetBodyForce()
          IF ( ASSOCIATED( BodyForce ) ) THEN
             bf_id = GetBodyForceId()
             g = 0.0_dp  
             g(1,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 1',Found)
             g(2,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Found)
             g(3,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Found)
             Gravity(1:N) = SQRT(SUM(g**2.0/N))
          END IF

          Density(1:N) = ListGetReal( Material, 'Water Density',  N, Element % NodeIndexes, Found )
          IF (.NOT.Found) THEN
             Density = 1000.0D00
             WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Water Density< not found for element ',&
                  t, ' material ', material_id
             CALL INFO(SolverName,Message,Level=4)
          END IF

          DO i=1,N

             IF (((UpperLimit(WloadPerm(Element % NodeIndexes(i)))-&
                  Wload(WloadPerm(Element % NodeIndexes(i))))*Density(i)*Gravity(i)).LE.NOpenLim)THEN

                IF (OpenCEL(PipingPerm(Element % NodeIndexes(i))).LT.2.0)&
                     OpenCEL(PipingPerm(Element % NodeIndexes(i))) = 1.0   


             END IF
          END DO
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
          IsPeriodicBC = GetLogical(BC,'Periodic BC ' // TRIM(Solver % Variable % Name),Found)
          IF (.NOT.Found) IsPeriodicBC = .FALSE.
          IF (IsPeriodicBC) THEN 
             DO i=1,N
                IF  (ActiveSed(WloadPerm(Element % NodeIndexes(i)))) THEN
                   k = k + 1
                   ActiveSed(WloadPerm(Element % NodeIndexes(i))) = .FALSE.
                END IF
             END DO
          END IF
       END IF
    END DO

    !---------------------------------------
    !Coppying the dummy to the reality
    !---------------------------------------
    SedLoadPerm = WloadPerm
    SedLoad = Wload

    PointerToSedResidual = PointerToResidualVector 
    PTSRPerm = PTRPerm

   WRITE(Message,'(a,e13.6,a,e13.6)') &
         'End Max/min values of sediment Water load:', ParallelReduction(MAXVAL( Wload(:))), &
         '/',ParallelReduction(MINVAL( Wload(:)))
    CALL INFO(SolverName,Message,Level=4)
    WRITE(Message,'(a,e13.6,a,e13.6)') &
         'Max/min values of sedimentload:', ParallelReduction(MAXVAL( SedLoad(:))), &
         '/',ParallelReduction(MINVAL( SedLoad(:)))
    CALL INFO(SolverName,Message,Level=4)

  END SUBROUTINE SedimentSolve



  SUBROUTINE CELSolve(Model, Solver, VarWload, TransientSimulation, Stabilize, &
       ApplyDirichlet, UseBubbles, ActiveCEL, Piping, LinearTol, WloadSol, &
       WatComp, CELNorm, VarWloadResidual, VarSedResidual, VarCELResidual, CELLoad, &
       Ppress, SolverName, PreviousCEL, SedToCEL, ClosureFlux, CurrentLoad, &
       ResLayer, iter)
    !------------------------------------------------------------------------------------------------
    !We are solving here the Hydro equation for a pipe equivalent layer, This is done to deal with the
    !upper bound of the sediment layer. When water load get to the Upper bound, the exess of water is 
    !routed in these pipes 
    !------------------------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Variable_t), POINTER ::  VarWload, Piping, VarWloadResidual, VarSedResidual, &
         VarCELResidual, CELLoad, WloadSol, ResLayer
    LOGICAL :: TransientSimulation, Stabilize, ApplyDirichlet, UseBubbles
    LOGICAL, ALLOCATABLE :: ActiveCEL(:)
    REAL(KIND=dp) :: LinearTol, WatComp, CELNorm
    REAL(KIND=dp), POINTER ::Ppress(:)
    REAL(KIND=dp), ALLOCATABLE :: PreviousCEL(:,:), SedToCEL(:), ClosureFlux(:), &
         CurrentLoad(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    INTEGER :: iter
    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------

    TYPE(Nodes_t) :: ElementNodes
    TYPE(Matrix_t), POINTER :: Systemmatrix
    TYPE(Element_t),POINTER :: Element
    TYPE(ValueList_t), POINTER :: Equation,Material,BodyForce,BC
    INTEGER :: i,j,k,l,m,n,t,p, body_id, material_id, &
         bf_id, bc_id, istat

    INTEGER, POINTER :: NodeIndexes(:)

    CHARACTER(LEN=MAX_NAME_LEN) :: AreaName

    LOGICAL ::Found, FluxBC, IsPeriodicBC=.FALSE.,&
         AllocationsDone = .FALSE.

    REAL(KIND=dp), POINTER :: WloadHomologous(:), Hwrk(:,:,:), &
         ResidualVector(:), ForceVector(:)

    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), &
         Transmitivity(:,:,:), FORCE(:),TimeForce(:), &
         TransferCoeff(:), Work(:), C1(:), C0(:), Zero(:), Viscosity(:),&
         StoringCoef(:), Density(:), WloadExt(:), &
         StiffVector(:), OldValues(:), OldRHS(:), VNull(:), Nochange(:,:), & 
         Gravity(:), SedComp(:), SedThick(:) , &
         CELThick(:) ,Porosity(:), Pressure(:), g(:,:),&
         LocalMassMatrix(:,:), LocalStiffMatrix(:,:), LocalForce(:),&
         LowerLimit(:), LocalArea(:)

    REAL(KIND=dp) :: Unorm, Qnode

    SAVE                        &
         OldValues,             &
         OldRHS,                &
         VNull,                 &
         Nochange,              &
         Pressure,              &
         ElementNodes,          &
         Work,Zero,             &
         Viscosity,             &
         StoringCoef,           &
         Density,               &
         WloadExt,              &
         C1,                    &
         C0,                    &
         TransferCoeff,         &
         Transmitivity,         &
         MASS,                  &
         STIFF,LOAD,            &
         FORCE,                 &
         TimeForce,             &
         StiffVector,           &
         ResidualVector,        &
         AllocationsDone,       &
         Hwrk,                  &
         Porosity,              &
         Gravity,               &
         SedComp,               &
         SedThick,              &
         CELThick,              &
         g,                     &
         LocalMassMatrix,       &
         LocalStiffMatrix,      &
         LocalForce,            &
         LowerLimit,            &
         LocalArea

    SystemMatrix => Solver % Matrix
    ForceVector => Solver % Matrix % RHS


    IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
       N = Solver % Mesh % MaxElementNodes
       M = LocalNodes
       K = SIZE( SystemMatrix % Values )
       L = SIZE( SystemMatrix % RHS )

       IF ( AllocationsDone ) THEN
          DEALLOCATE(                    &
               OldValues,                &
               OldRHS,                   &
               VNull,                    &
               Nochange,                 &
               Pressure,                 &
               ElementNodes % x,         &
               ElementNodes % y,         &
               ElementNodes % z,         &
               Work,Zero,                &
               Viscosity,                &
               StoringCoef,              &
               Density,                  &
               WloadExt,                 &
               C1,                       &
               C0,                       &
               TransferCoeff,            &
               Transmitivity,            &
               MASS,                     &
               STIFF,                    &
               LOAD,                     &
               FORCE,                    &
               TimeForce,                &
               StiffVector,              &
               ResidualVector,           &
               Porosity,                 &
               Gravity,                  &
               SedComp,                  &
               SedThick,                 &
               CELThick,                 &
               g,                        &
               LocalMassMatrix,          &
               LocalStiffMatrix,         &
               LocalForce,               &
               LowerLimit,               &
               LocalArea)              
       END IF

       ALLOCATE(                         &
            OldValues( K ),              &
            OldRHS( L ),                 &
            VNull( N ),                  &
            Nochange( 3,N ),             &
            Pressure( N ),               &
            ElementNodes % x( N ),       &
            ElementNodes % y( N ),       &
            ElementNodes % z( N ),       &
            Work( N ), Zero( N ),        &
            Viscosity( N ),              &
            StoringCoef( N ),            &
            Density( N ),                &
            WloadExt( N ),               &
            C1( N ),                     &
            C0( N ),                     &
            TransferCoeff( N ),          &
            Transmitivity( 3,3,N ),      &
            MASS(  2*N,2*N ),            &
            STIFF( 2*N,2*N ),            &
            LOAD( N ),                   &
            FORCE( 2*N ),                &
            TimeForce( 2*N ),            &
            StiffVector( L ),            &
            ResidualVector(L),           &
            Porosity( N ),               &
            Gravity( N ),                &
            SedComp( N ),                &
            SedThick( N ),               &
            CELThick( N ),              &
            g( 3,N ),                    & 
            LocalMassMatrix( 2*N,2*N ),  &
            LocalStiffMatrix( 2*N,2*N ), &
            LocalForce( 2*N ),           &
            LowerLimit( M ),             &
            LocalArea( N ),              &
            STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL FATAL( SolverName, 'Memory allocation error 2' )
       ELSE
          CALL INFO(SolverName, 'Memory allocation done here 2', level=1 )
       END IF
       AllocationsDone = .TRUE.

    END IF

    !---------------------------------------
    !Coppying the reality to the dummy
    !---------------------------------------
    IF ( TransientSimulation ) THEN
       WloadSol % PrevValues = PreviousCEL
    END IF

    Solver % Variable % Name = 'CEL'

    WloadPerm = CELLoadPerm
    Wload = DrainageLoad

    PointerToResidualVector = PointerToCELResidual
    PTRPerm = PTPRPerm


    !------------------------------------------------------------------------------
    ! lets start
    !------------------------------------------------------------------------------
    CALL DefaultInitialize() 
    !-----------------------------------------------------------------------------
    ! Get lower limit:
    !-----------------------------------------------------------------------------
    DO t=1,Solver % NumberOfActiveElements
       
       Element => GetActiveElement(t)
       n= GetElementNOFNodes()
       CALL GetElementNodes( ElementNodes,Element )
       Material => GetMaterial()
       
       LowerLimit(WloadPerm(Element % Nodeindexes(1:n))) = ListGetReal(Material,TRIM(VariableName) // & 
            ' Lower Limit',n,Element % NodeIndexes, Found) 
    END DO
    
    AreaName = GetString(Constants,'Area Name', Found)
    IF (.NOT.Found) THEN
       CALL WARN('HydroSolver', 'No Keyword >Area Name< defined. Using >Area< as default.')
       WRITE(AreaName,'(A)') 'Area'
    ELSE
       WRITE(Message,'(a,a)') 'Variable Name for Element area: ', AreaName
       CALL INFO('HydroSolver',Message,Level=12)
    END IF
    !------------------------------------------------------------------------------
    ! write some info on max/min values
    !------------------------------------------------------------------------------

    WRITE(Message,'(a,e13.6,a,e13.6)') &
         'Max/min values of CELd Water load:', ParallelReduction(MAXVAL( Wload(:))), &
         '/', ParallelReduction(MINVAL( Wload(:)))
    CALL INFO(SolverName,Message,Level=4)  
    body_id = -1
    NULLIFY(Material)

    !------------------------------------------------------------------------------
    ! Bulk elements
    !------------------------------------------------------------------------------

    DO t=1,Solver % NumberOfActiveElements
    
       !------------------------------------------------------------------------------
       ! Check if this element belongs to a body where scalar equation
       ! should be calculated
       !------------------------------------------------------------------------------

       Element => GetActiveElement(t,Solver)
       IF (.NOT.ASSOCIATED(Element)) CYCLE
       IF ( Element % BodyId /= body_id ) THEN
!!$          
!!$          !Okay in serial??
!!$          IF (Element % PartIndex .NE. Parenv % mype) CYCLE 
          
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

       k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
            minv=1, maxv=Model % NumberOFEquations )

       !--------------------------------
       ! Get element material parameters
       !--------------------------------              

       N = GetElementNOFNodes(Element)
       CALL GetElementNodes( ElementNodes,Element )

       Porosity(1:N) = listGetReal( Material,'Sediment Porosity', N, Element % NodeIndexes, Found )

       IF (.NOT.Found) THEN
          Porosity = 0.4D00
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

       SedComp(1:N) = listGetReal( Material,'Sediment Compressibility', N, Element % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          SedComp = 1.0D-2
          WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Compressibility', &
               '< not found for element ', t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       SedThick(1:N) = listGetReal( Material,'Sediment Thickness', N, Element % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          SedThick = 10.0D00
          WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Thickness', &
               '< not found for element ', t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       CELThick(1:N) = listGetReal( Material,'CEL Thickness', N, Element % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          CELThick = 1.0D00
          WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' CEL Thickness', &
               '< not found for element ', t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       Density(1:N) = ListGetReal( Material, 'Water Density',  N, Element % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          Density = 1000.0D00
          WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Water Density< not found for element ',&
               t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       !Getting transmitivity of the two CEL
       !-------------------------------------------
       CALL ListGetRealArray( Material,'CEL Transmitivity',Hwrk,N, Element % NodeIndexes )
       Transmitivity = 0.0D0


       IF (SIZE(Hwrk,1) .EQ. 1) THEN
          DO i=1,3
             DO j=1,N
                Transmitivity(i,i,j) = Hwrk(1,1,j)
             END DO
          END DO
       ELSE
          WRITE(Message,'(a,a,a)') 'Keyword >CEL Transmitivity< should be isotrop '
          CALL INFO(SolverName,Message,Level=4)

       END IF

       !Computation of storing coefficient 
      !----------------------------------
       StoringCoef(1:N) = Gravity(1:N) * Porosity(1:N) * &
            Density(1:N) * CELThick(1:N) * &
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
       ! Add sources, flux at the closure of the pipe, transfer from the sediment and 
       ! flux from the overload of the sediment. They are given in water load
       !------------------------------------------------------------------------------

       LOAD = 0.0D00
       CurrentLoad(WloadPerm(Element % NodeIndexes(1:N))) = 0.0D00

       CALL GetScalarLocalSolution(LocalArea,AreaName)


       DO i=1,N
         
          ResLayValues(ResLayPerm(Element % NodeIndexes(i))) = &
               PointerToSedResidual(PTSRPerm(Element % NodeIndexes(i)))/LocalArea(i)
          
          LOAD(i) = SedToCEL(WloadPerm(Element % NodeIndexes(i)))&
               - ResLayValues(ResLayPerm(Element % NodeIndexes(i)))

          CurrentLoad(WloadPerm(Element % NodeIndexes(i))) = LOAD(i) + CurrentLoad(WloadPerm(Element % NodeIndexes(i)))


          LOAD(i) =  LOAD(i) - ClosureFlux(WloadPerm(Element % NodeIndexes(i)))
              
          
       END DO

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
          CALL DiffuseConvectiveCompose( &
               MASS, STIFF, FORCE, LOAD, &
               StoringCoef, C0, C1(1:N), Transmitivity, &
               .FALSE., Zero, Zero, VNull, VNull, VNull, &
               Nochange(1,1:N),Nochange(2,1:N),Nochange(3,1:N),&
               Viscosity, Density, Pressure, Zero, Zero,&
               .FALSE., Stabilize, UseBubbles, Element, N, ElementNodes )


          ! special coords (account for metric)
          !-----------------------------------
       ELSE
          CALL DiffuseConvectiveGenCompose( &
               MASS, STIFF, FORCE, LOAD, &
               StoringCoef, C0, C1(1:N), Transmitivity, &
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
       !------------------------------------------------------------------------------
       IF (  UseBubbles ) THEN
          CALL Condensate( N, STIFF, FORCE, TimeForce )
          IF (TransientSimulation) CALL DefaultUpdateForce( TimeForce )
       END IF

       CALL DefaultUpdateEquations( STIFF, FORCE )


    END DO     ! Bulk elements

    !------------------------------------------------------------------------------
    ! Neumann & Newton boundary conditions
    !------------------------------------------------------------------------------
    DO t=1, Solver % Mesh % NumberOfBoundaryElements

       ! get element information
       Element => GetBoundaryElement(t)
       bc_id = GetBCId( Element )

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
             !------------------------------
             !BC: -k@Hw/@n = a(Hw - HwExt)
             !Check it if you want to use it
             !------------------------------
             TransferCoeff(1:n) = GetReal( BC, TRIM(Solver % Variable % Name) //  ' Transfer Coefficient',Found )
             IF ( ANY(TransferCoeff(1:n) /= 0.0d0) ) THEN
                WloadExt(1:n) = GetReal( BC, TRIM(Solver % Variable % Name) // ' External Value',Found )   
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


    END DO   ! Neumann & Newton BCs
    !------------------------------------------------------------------------------

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
       DO i=1,Model % Mesh % NumberOfNodes

          k = WloadPerm(i)
          IF (k > 0) THEN
             IF (ActiveCEL(k)) THEN
                CALL ZeroRow( SystemMatrix, k ) 
                CALL SetMatrixElement( SystemMatrix, k, k, 1.0d0 )
                SystemMatrix % RHS(k) = LowerLimit(k)

             END IF
          END IF
       END DO
    END IF

    CALL Info( TRIM(SolverName) // ' CEL equivalent layer', ' Assembly done', Level=4 )

    !------------------------------------------------------------------------------
    !     Solve the system 
    !------------------------------------------------------------------------------

    CELNorm = DefaultSolve()

    SystemMatrix % Values = OldValues
    ForceVector = OldRHS

    !------------------------------------------------------------------------------
    ! compute residual
    !------------------------------------------------------------------------------ 
    IF ( ParEnv % PEs.GT.1 ) THEN !!!!!!!!!!!!!!!!!!!!!! we have a parallel run
       CALL ParallelInitSolve( SystemMatrix, Wload, ForceVector, ResidualVector )
       CALL ParallelMatrixVector( SystemMatrix, Wload, StiffVector, .TRUE. )
       ResidualVector =  StiffVector - ForceVector
       CALL ParallelSumVector( SystemMatrix, ResidualVector )
    ELSE !!!!!!!!!!!!!!!!!!!!!! serial run 

       CALL CRS_MatrixVectorMultiply( SystemMatrix, Wload, StiffVector)
       ResidualVector =  StiffVector - ForceVector

    END IF

    !-----------------------------
    ! determine "active" nodes set 
    !-----------------------------

!!$    IF (ASSOCIATED(VarWload)) THEN
!!$       WloadHomologous => VarWload % Values
!!$       DO i=1,Model % Mesh % NumberOfNodes ! <______________IS THIS OK IN PARALLEL????????
!!$
!!$          k = VarWload % Perm(i)
!!$          l = WloadPerm(i)
!!$          IF ( (k.LE.0).OR.(l.LE.0) ) CYCLE
!!$
!!$          WloadHomologous(k) = Wload(l) - LowerLimit(l)
!!$
!!$          IF (ApplyDirichlet) THEN
!!$             !----------------------------------------------------------
!!$             ! if upper limit is exceeded, manipulate matrix in any case
!!$             !----------------------------------------------------------
!!$             IF (WloadHomologous(k).LE. 0.0 )THEN
!!$                ActiveCEL(l) = .TRUE.
!!$                WloadHomologous(k) = LinearTol
!!$             END IF
!!$
!!$             !---------------------------------------------------
!!$             ! if there is "heating", don't manipulate the matrix
!!$             !---------------------------------------------------
!!$             IF (ResidualVector(l).LT.-LinearTol &
!!$                  .AND. iter.GT.1) ActiveCEL(l) = .FALSE.
!!$          END IF
!!$          IF(.NOT.ActiveCEL(l)) THEN
!!$             PointerToResidualVector(PTRPerm(i)) = 0.0D00
!!$          ELSE
!!$             PointerToResidualVector(PTRPerm(i)) = ResidualVector(l)
!!$          END IF
!!$       END DO
!!$    ELSE
!!$       WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' Homologous not associated'
!!$       CALL WARN( SolverName, Message)
!!$    END IF

    IF (ASSOCIATED(VarWload)) THEN
       WloadHomologous => VarWload % Values
       DO t=1,Solver % NumberOfActiveElements
          Element => GetActiveElement(t,Solver)
          IF (.NOT.ASSOCIATED(Element)) CYCLE
          N = GetElementNOFNodes(Element)
          CALL GetElementNodes( ElementNodes,Element )

          DO i=1,N
             k = VarWload % Perm(Element % NodeIndexes(i))
             l = WloadPerm(Element % NodeIndexes(i))
             IF ( (k.LE.0).OR.(l.LE.0)) CYCLE
             
             WloadHomologous(k) = Wload(l) - LowerLimit(l)
             
             IF (ApplyDirichlet) THEN
                !----------------------------------------------------------
                ! if upper limit is exceeded, manipulate matrix in any case
                !----------------------------------------------------------
                IF (WloadHomologous(k).LE.0.0 ) THEN
                   ActiveSed(l) = .TRUE.
                   WloadHomologous(k) = LinearTol
                END IF
                !---------------------------------------------------
                ! if there is "heating", don't manipulate the matrix
                !---------------------------------------------------
                IF( (ResidualVector(l).LT.-LinearTol) &
                     .AND.( iter.GT.1)) ActiveCEL(l) = .FALSE.
             END IF
             IF( .NOT.ActiveCEL(l) ) THEN
                PointerToResidualVector(PTRPerm(Element % NodeIndexes(i))) = 0.0D00
             ELSE
                PointerToResidualVector(PTRPerm(Element % NodeIndexes(i))) = ResidualVector(l)      
             END IF
          END DO
       END DO
    ELSE
       WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' Homologous not associated'
       CALL WARN( SolverName, Message)
    END IF

    !------------------------------------------
    ! special treatment for periodic boundaries
    !------------------------------------------

    k=0
    DO t=1, Solver % Mesh % NumberOfBoundaryElements

       ! get element information
       Element => GetBoundaryElement(t)
       IF ( .NOT.ActiveBoundaryElement() ) CYCLE
       n = GetElementNOFNodes()
       !IF ( GetElementFamily() == 1 ) CYCLE
       BC => GetBC()
       bc_id = GetBCId( Element )

       CALL GetElementNodes( ElementNodes,Element )


       IF ( ASSOCIATED( BC ) ) THEN    
          IsPeriodicBC = GetLogical(BC,'Periodic BC ' // TRIM(Solver % Variable % Name),Found)
          IF (.NOT.Found) IsPeriodicBC = .FALSE.
          IF (IsPeriodicBC) THEN 
             DO i=1,N
                IF  (ActiveCEL(WloadPerm(Element % NodeIndexes(i)))) THEN
                   k = k + 1
                   ActiveCEL(WloadPerm(Element % NodeIndexes(i))) = .FALSE.
                END IF
             END DO
          END IF
       END IF
    END DO


    !---------------------------------------
    !Coppying the dummy to the reality
    !---------------------------------------

    CELLoadPerm = WloadPerm
    DrainageLoad = Wload

    PointerToCELResidual = PointerToResidualVector 
    PTPRPerm = PTRPerm
    Solver % Variable % Name = 'Hwater'

   WRITE(Message,'(a,e13.6,a,e13.6)') &
         'EndMax/min values of CELd Water load:', ParallelReduction(MAXVAL( Wload(:))), &
         '/',ParallelReduction(MINVAL( Wload(:)))
    CALL INFO(SolverName,Message,Level=4)  

  END SUBROUTINE CELSolve
  !-----------------------------------------------------------------------------


  SUBROUTINE WaterTransfer(Model, Solver, SedLoad, DrainageLoad, &
       Piping, SedToCEL ,Solvername, LocalNodes, Relax)
    !-----------------------------------------------------------------------------
    !When the water load in the sediment is higher than the one in the pipes, some
    !of the water pass from the sediment to the pipe.

    TYPE(Model_t)  :: Model
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Variable_t) :: Piping
    REAL(KIND=dp),ALLOCATABLE :: SedToCEL(:)
    REAL(KIND=dp), POINTER :: SedLoad(:), DrainageLoad(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    INTEGER :: LocalNodes
    REAL(KIND=dp) :: Relax
    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------

    TYPE(Element_t),POINTER :: CurrentElement, Element
    TYPE(ValueList_t), POINTER :: Material, BodyForce
    TYPE(Variable_t), POINTER:: TimeVar
    INTEGER :: i,j,k,n,t,m,istat, bc_id, material_id
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: ReferenceLoad
    REAL(KIND=dp), ALLOCATABLE :: Taupe(:), StoringCoef(:), Density(:), &
         Gravity(:), SedComp(:), SedThick(:),CELThick(:), Porosity(:), &
         UpperLimit(:), g(:,:), Coupling(:), OldTrans(:)
    LOGICAL :: AllocationsDone = .FALSE., Changed = .TRUE.


    SAVE AllocationsDone, &
         Taupe,           & 
         Changed,         &
         StoringCoef,     &
         Density,         &
         Gravity,         &
         SedComp,         &
         SedThick,        &
         CELThick,        &
         Porosity,        &
         UpperLimit,      &
         g,               &
         Coupling,        &
         OldTrans


    IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
       N = Solver % Mesh % MaxElementNodes
       M = LocalNodes

       IF ( AllocationsDone ) THEN

          DEALLOCATE(Taupe, &
               StoringCoef, &
               Density,     &
               Gravity,     &
               SedComp,     &
               SedThick,    &
               CELThick,    &
               Porosity,    &
               UpperLimit,  &
               g,           &
               Coupling,    &
               OldTrans)              

       END IF

       ALLOCATE(Taupe( N ),   &
            StoringCoef( N ), &
            Density( N ),     &
            Gravity( N ),     &
            SedComp( N ),     &
            SedThick( N ),    &
            CELThick( N ),    &
            Porosity( N ),    &
            UpperLimit( M ),  &
            g( 3,N ),         &
            Coupling( N ),    &
            OldTrans( M ),    &
            STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL FATAL( SolverName, 'Memory allocation error 3' )
       ELSE
          CALL INFO(SolverName, 'Memory allocation done 3', level=1 )
       END IF

       AllocationsDone = .TRUE.
    END IF

    !Creating a tab to know if the pipe is linked to the "exterior". This is set by the argument 
    !"Open CEL" in BC part of the sif, this tab is piping with value of two for an open pipe 
    !linked to the "exterior"

    OldTrans = SedToCEL

    Changed = .TRUE.

    DO WHILE (Changed)

       Changed = .FALSE.
       DO t=1, Solver % Mesh % NumberOfBoundaryElements

          ! get element information
          Element => GetBoundaryElement(t)
          IF ( .NOT.ActiveBoundaryElement() ) CYCLE
          n = GetElementNOFNodes()

          CALL GetElementNodes( ElementNodes,Element )

          !some case did not need treatment, not opened (0), every node active (ALL 2)
          !or no active node(ANY 2)
          !---------------------------------------------------------------------------

          IF (SUM(OpenCEL(PipingPerm(Element % NodeIndexes(1:N)))).EQ.0.0) CYCLE
          IF (ALL(OpenCEL(PipingPerm(Element % NodeIndexes(1:N))).EQ.2.0)) CYCLE
          IF (.NOT.(ANY(OpenCEL(PipingPerm(Element % NodeIndexes(1:N))).GE.2.0))) CYCLE

          !If a node is efficient and upstream of other nodes of the element, these nodes
          !become active
          !------------------------------------------------------------------------------

          DO i=1,N 
             IF (OpenCEL(PipingPerm(Element % NodeIndexes(i))).GE.2.0) THEN

                ReferenceLoad = DrainageLoad(CELLoadPerm(Element % NodeIndexes(i)))
                DO j = 1,N
                   IF (i.NE.j)THEN ! ajout du 19/06

                      IF((OpenCEL(PipingPerm(Element % NodeIndexes(j))).EQ.1.0).AND. &
                           (DrainageLoad(CELloadPerm(Element % NodeIndexes(j))).GE.ReferenceLoad))THEN

                         Changed = .TRUE.
                         OpenCEL(PipingPerm(Element % NodeIndexes(j))) = 2.0
                      END IF
                   END IF
                END DO
             END IF
          END DO

       END DO

    END DO

    DO t=1,Solver % NumberOfActiveElements
       CurrentElement => GetActiveElement(t)
       Material => GetMaterial(CurrentElement)
       n = GetElementNOFNodes()
       material_id = GetMaterialId(CurrentElement, Found)

       UpperLimit(WloadPerm(CurrentElement % Nodeindexes(1:n))) = ListGetReal(Material,TRIM(VariableName) // & 
            ' Upper Limit',n,CurrentElement % NodeIndexes, Found)
       !Getting the variables to compute storingcoefficient
       !---------------------------------------------------

       Porosity(1:N) = listGetReal( Material,'Sediment Porosity', N, CurrentElement % NodeIndexes, Found )

       IF (.NOT.Found) THEN
          Porosity = 0.4D00
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

       SedComp(1:N) = listGetReal( Material,'Sediment Compressibility', N, CurrentElement % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          SedComp = 1.0D-2
          WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Compressibility', &
               '< not found for element ', t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       SedThick(1:N) = listGetReal( Material,'Sediment Thickness', N, CurrentElement % NodeIndexes, Found )

       IF (.NOT.Found) THEN
          SedThick = 10.0D00
          WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' Sediment Thickness', &
               '< not found for element ', t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       CELThick(1:N) = listGetReal( Material,'CEL Thickness', N, CurrentElement % NodeIndexes, Found )

       IF (.NOT.Found) THEN
          CELThick = 1.0D00
          WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', ' CEL Thickness', &
               '< not found for element ', t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF


       Density(1:N) = ListGetReal( Material, 'Water Density',  N, CurrentElement % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          Density = 1000.0D00
          WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Water Density< not found for element ',&
               t, ' material ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       !Storing coeficient is not complete as it is not multiplied by any layer thickness so far
       !----------------------------------------------------------------------------------------
       StoringCoef(1:N) = Gravity(1:N) * Porosity(1:N) * Density(1:N) * &
            (WatComp + SedComp(1:N)/Porosity(1:N))

       !Getting the time constant
       !---------------------------------------------------------------------
       Taupe(1:n) = ListGetReal( Material, 'Transfer time constant',  n, CurrentElement % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          Taupe = 1000.0D00
          WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Transfer time constant< not found for element ',&
               t, ' material,Set to 1e-3 a-1 ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       !Getting the Coupling Parameter
       !---------------------------------------------------------------------
       Coupling(1:n) = ListGetReal( Material, 'Coupling Parameter',  n, CurrentElement % NodeIndexes, Found )
       IF (.NOT.Found) THEN
          Coupling = 1
          WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Coupling Parameter< not found for element ',&
               t, ' material,Set to 1e-3 a-1 ', material_id
          CALL INFO(SolverName,Message,Level=4)
       END IF

       !Storing coeffcicient is multiplied by one or the other thickness. Sediment thickness if the water
       !pass from the sediment to the pipe and pipe thickness in the other case.
       !-------------------------------------------------------------------------------------------------

       DO j=1,n

          k = CurrentElement % NodeIndexes(j)

          IF(OpenCEL(PipingPerm(k)).LT.2.0)THEN

             IF((SedLoad(SedloadPerm(k)).LT.DrainageLoad(CELloadPerm(k))).AND.&
                  (SedLoad(SedloadPerm(k)).LT.UpperLimit(WloadPerm(k))))THEN

                SedToCEL(SedloadPerm(k)) = Coupling(j) &
                     *(SedLoad(SedloadPerm(k)) - DrainageLoad(CELloadPerm(k))) &
                     *StoringCoef(j)*CELThick(j)/Taupe(j)
             END IF

          ELSE

             IF((SedLoad(SedloadPerm(k)).GE.DrainageLoad(CELloadPerm(k))))THEN

                SedToCEL(SedloadPerm(k)) = Coupling(j) &
                     *(SedLoad(SedloadPerm(k)) - DrainageLoad(CELloadPerm(k))) &
                     *StoringCoef(j)*SedThick(j)/Taupe(j)

                SedToCEL(SedloadPerm(k)) = MIN(SedToCEL(SedloadPerm(k)),&
                     ((SedLoad(SedloadPerm(k)) - DrainageLoad(CELloadPerm(k)))*StoringCoef(j)*SedThick(j)))

             ELSE 

                SedToCEL(SedloadPerm(k)) = Coupling(j) &
                     *(SedLoad(SedloadPerm(k)) - DrainageLoad(CELloadPerm(k))) &
                     *StoringCoef(j)*CELThick(j)/Taupe(j)

             END IF
          END IF
       END DO
    END DO

    SedToCEL = (SedToCEL * Relax) + (1-Relax)*OldTrans

  END SUBROUTINE WaterTransfer
  !-----------------------------------------------------------------------------



  SUBROUTINE EndOfCEL (Model, SedLoad, DrainageLoad, WloadPerm, Piping, &
       ClosureFlux, CurrentLoad, LocalNodes)
    !----------------------------------------------------------------------------
    !If at a point, the pipe is closed, all the water that as entered the pipe 
    !should go back to the sediment (the slope that is considered for the 
    !direction of the flux is the slope of the water load in the pipe.)
    !----------------------------------------------------------------------------

    !External Variables
    !------------------
    TYPE(Model_t)  :: Model
    TYPE(Variable_t) :: Piping
    REAL(KIND=dp), POINTER ::  SedLoad(:), DrainageLoad(:)
    INTEGER, POINTER ::WloadPerm(:)
    REAL(KIND=dp), ALLOCATABLE :: ClosureFlux(:), CurrentLoad(:)
    INTEGER :: LocalNodes

    !Local Variables
    !---------------
    TYPE(Element_t),POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: InputNumber
    INTEGER :: i, N, t,k, M
    LOGICAL :: Changed=.TRUE., AllocationsDone=.FALSE.
    REAL(KIND=dp), ALLOCATABLE ::UpperLimit(:), Prefered(:)


    SAVE Changed, AllocationsDone, UpperLimit, Prefered

    !----------Allocation----------

    IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
       M = LocalNodes

       IF ( AllocationsDone ) THEN
          DEALLOCATE(UpperLimit, Prefered)              
       END IF

       ALLOCATE(UpperLimit( M ), Prefered( M ), STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL FATAL( SolverName, 'Memory allocation error 4' )
       ELSE
          CALL INFO(SolverName, 'Memory allocation done 4', level=1 )
       END IF

       AllocationsDone = .TRUE.
    END IF

    !----------Get Upper Limit----------

    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
       n = GetElementNOFNodes()
       CALL GetElementNodes( ElementNodes,Element )
       Material => GetMaterial()

       UpperLimit(WloadPerm(Element % Nodeindexes(1:n))) = ListGetReal(Material,TRIM(VariableName) // & 
            ' Upper Limit',n,Element % NodeIndexes, Found)

       IF (.NOT. Found) THEN
          WRITE(Message,'(a,i10)') 'No upper limit of solution for element no. ', t
          CALL INFO(SolverName, Message, level=10)
       END IF
    END DO

    k=0

    !In a first time we build a tab that sum the loads in the open pipes 
    !at the last opened node of the pipe
    ClosureFlux = 0.0
    Changed = .TRUE.
    Prefered = 0.0

    DO WHILE (Changed)
       k=k+1
       Changed = .FALSE.

       DO t=1, Solver % Mesh % NumberOfBoundaryElements
          Element => GetBoundaryElement(t)
          IF ( .NOT.ActiveBoundaryElement() ) CYCLE
          N = GetElementNOFNodes()

          CALL GetElementNodes( ElementNodes,Element )

          IF((SUM(OpenCEL(PipingPerm(Element % NodeIndexes(1:N)))).GT.0.0).AND.&
               (SUM(CurrentLoad(WloadPerm(Element % NodeIndexes(1:N)))).GT.0.0))THEN

             DO i=1,N
                
                IF(PipingPerm(Element % NodeIndexes(i)).LE.0)CYCLE
                
                IF((Prefered(PipingPerm(Element % NodeIndexes(i))).NE.1.0).AND.&
                     (ANY(Prefered(PipingPerm(Element % NodeIndexes(1:N))).EQ.1.0)))CYCLE

                IF(ANY(OpenCEL(PipingPerm(Element % NodeIndexes(1:N))).EQ.0.0))THEN

                   IF ((Model % Nodes % z (Element % NodeIndexes(i))).EQ. &
                        MINVAL(Model % Nodes % z (Element % NodeIndexes(1:N))).AND.&
                        OpenCEL(PipingPerm(Element % NodeIndexes(i))).EQ.0.0)THEN


                      ClosureFlux(WloadPerm(Element % NodeIndexes(i))) = SUM(CurrentLoad(WloadPerm(Element % NodeIndexes(1:N))))

                      IF(ClosureFlux(WloadPerm(Element % NodeIndexes(i))).NE.&
                           CurrentLoad(WloadPerm(Element % NodeIndexes(i)))) THEN 
                         Prefered(PipingPerm(Element % NodeIndexes(i))) = 1.0
                         Changed = .TRUE.
                      END IF
                      CurrentLoad(WloadPerm(Element % NodeIndexes(1:N))) = ClosureFlux(WloadPerm(Element % NodeIndexes(1:N)))
                      ClosureFlux(WloadPerm(Element % NodeIndexes(1:N))) = 0.0

                   END IF
                ELSEIF ((Model % Nodes % z (Element % NodeIndexes(i))).EQ. &
                     MINVAL(Model % Nodes % z (Element % NodeIndexes(1:N))))THEN
                   
                   ClosureFlux(WloadPerm(Element % NodeIndexes(i))) = SUM(CurrentLoad(WloadPerm(Element % NodeIndexes(1:N))))
                   
                   IF(ClosureFlux(WloadPerm(Element % NodeIndexes(i))).NE.&
                        CurrentLoad(WloadPerm(Element % NodeIndexes(i)))) THEN 
                      Changed = .TRUE.
                   END IF
                   CurrentLoad(WloadPerm(Element % NodeIndexes(1:N))) = ClosureFlux(WloadPerm(Element % NodeIndexes(1:N)))
                   ClosureFlux(WloadPerm(Element % NodeIndexes(1:N))) = 0.0
                   
                END IF
             END DO
          END IF
       END DO

       IF (k.GT.100) THEN
          write(*,*)"Too many iterations in End of CEL, Stoping"
          STOP
       END IF
    END DO

    !In a second time, we input the preceding load at the first closed point

    DO t=1, Solver % Mesh % NumberOfBoundaryElements

       ! get element information
       Element => GetBoundaryElement(t)
       IF ( .NOT.ActiveBoundaryElement() ) CYCLE
       N = GetElementNOFNodes()

       CALL GetElementNodes( ElementNodes,Element )

       IF (SUM(CurrentLoad(WloadPerm(Element % NodeIndexes(1:N)))).GT.0.0)THEN

          InputNumber = REAL(N)
          DO i=1,N
             IF(((OpenCEL(PipingPerm(Element % NodeIndexes(i)))).GT.0.0).OR.& 
                  (SedLoad(SedloadPerm(Element % NodeIndexes(i))).GE.UpperLimit(WloadPerm(Element % NodeIndexes(i)))))THEN 
                InputNumber = InputNumber - 1
                Prefered(PipingPerm(Element % NodeIndexes(i))) = 0.0
             END IF
          END DO

          DO i = 1,N
             IF(OpenCEL(PipingPerm(Element % NodeIndexes(i))).EQ.0.0)THEN

                IF(SedLoad(SedloadPerm(Element % NodeIndexes(i))).LT. &
                     UpperLimit(WloadPerm(Element % NodeIndexes(i))))THEN

                   IF(InputNumber.EQ.0)THEN
                      OpenCEL(PipingPerm(Element % NodeIndexes(i))) = 1
                      Prefered(PipingPerm(Element % NodeIndexes(i))) = 0.0

                   ELSE
                      ClosureFlux(WloadPerm(Element % NodeIndexes(i))) = &
                           SUM(CurrentLoad(WloadPerm(Element % NodeIndexes(1:N))))/InputNumber

                   END IF
                END IF
             END IF
          END DO
       END IF
    END DO


  END SUBROUTINE EndOfCEL


  !------------------------------------------------------------------------------

END SUBROUTINE HydroSolver
!------------------------------------------------------------------------------
