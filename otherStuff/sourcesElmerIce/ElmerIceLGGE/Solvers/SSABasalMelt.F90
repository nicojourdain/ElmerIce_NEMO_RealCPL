!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini             
! *  Email:   gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 30. April 2010
! * 
! *****************************************************************************
!> This solver compute either the basal melt and/or the heat induced by friction
!> melt = (G + ub.taub)/(rhoi L)
!> heat = ub.taub
!> work only if the basal velocity are from the SSA solution. 
SUBROUTINE SSABasalMelt( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
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
  TYPE(Solver_t),TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, Material
  TYPE(Variable_t), POINTER :: VeloSol, MeltSol, HeatSol, NSol

  LOGICAL, SAVE :: AllocationsDone = .FALSE., Heat, Melt, GeoFlux = .FALSE.
  LOGICAL :: Found, GotIt, FirstTime = .True.

  INTEGER :: i, j, n, M, t, istat, DIM, VeloSTDOFs
  INTEGER, POINTER :: VeloPerm(:), MeltPerm(:), HeatPerm(:), NPerm(:)

  REAL(KIND=dp), POINTER :: Velo(:), MeltValue(:), HeatValue(:), Nval(:)
  REAL(KIND=dp), SAVE :: LatentHeat 
  REAL(KIND=dp) :: fm, NVelo, alpha, Xi, PostPeak, tbub, CN
  REAL(KIND=dp), ALLOCATABLE :: LocalVelo(:,:), LocalBeta(:), LocalC(:), &
     LocalG(:), LocalN(:), LocalDens(:)
  INTEGER :: iFriction
  INTEGER, POINTER :: NodeIndexes(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: Friction

  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: SolverName = 'SSABasalMelt'

  SAVE :: LocalVelo, LocalBeta, LocalC, LocalG, LocalN, LocalDens

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF (FirstTime) THEN
     FirstTime = .False. 
     Heat = GetLogical(Solver % Values, 'Compute SSA Friction Heat',Found)
     Melt = GetLogical(Solver % Values, 'Compute SSA Friction Melt',Found)
     IF ((.Not.Heat).AND.(.Not.Melt)) CALL FATAL(SolverName,'Nothing to do!')
     IF (Heat) CALL INFO( SolverName, 'Solve for basal friction heat',Level=5 )
     IF (Melt) THEN
        CALL INFO( SolverName, 'Solve for basal melt induced by friction',Level=5 )
        GeoFlux = GetLogical(Solver % Values, 'Add Geothermal Heat Flux',Found)
        IF (GeoFlux) CALL INFO( SolverName, 'Geothermal HeatFlux will be added &
               & to compute the melt',Level=5 )
     END IF
  END IF

 ! Get Pointer to auxiliary variables
  VeloSol => VariableGet( Solver % Mesh % Variables, 'SSAVelocity' )
  IF (ASSOCIATED(VeloSol)) THEN
      Velo=> VeloSol % Values
      VeloPerm => VeloSol % Perm
      VeloSTDOFs = VeloSol % DOFs
     
  ELSE
      CALL FATAL(SolverName,'Could not find variable >SSAVelocity<')
  END IF

  IF (Heat) THEN
     HeatSol => VariableGet( Solver % Mesh % Variables, 'SSAFrictionHeat' )
     IF (ASSOCIATED(HeatSol)) THEN
        HeatValue => HeatSol % Values
        HeatPerm => HeatSol % Perm
     ELSE
        CALL FATAL(SolverName,'Could not find variable >SSAFrictionHeat<')
     END IF
  END IF

  IF (Melt) THEN
     MeltSol => VariableGet( Solver % Mesh % Variables, 'SSAFrictionMelt' )
     IF (ASSOCIATED(MeltSol)) THEN
        MeltValue => MeltSol % Values
        MeltPerm => MeltSol % Perm
     ELSE
        CALL FATAL(SolverName,'Could not find variable >SSAFrictionMelt<')
     END IF
  END IF

  NSol => VariableGet( Solver % Mesh % Variables, 'Effective Pressure' )
  IF (ASSOCIATED(NSol)) THEN
     Nval => NSol % Values
     NPerm => NSol % Perm
  END IF

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     ! Allocate
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(LocalVelo, LocalBeta, LocalC, &
                        LocalG, LocalN, LocalDens)

     ALLOCATE(LocalVelo(2,M), LocalBeta(M), LocalC(M),  &
           LocalG(M), LocalN(M), LocalDens(M), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

  IF (Melt) THEN 
     ! Get Latent heat  and water density from the constant section
     LatentHeat = GetCReal( Model % Constants, 'Latent Heat', Found )
     If (.NOT.Found) THEN
       WRITE(Message,'(A)') 'Constant >Latent Heat< not found. &
         &Setting to 3.35e5_dp'
       CALL INFO(SolverName, Message, level=3)
       LatentHeat = 3.35e5_dp
     End if
  END IF   

  ! Loop over elements (but only nodal evaluation needed)
  DO t=1,Solver % NumberOFActiveElements
     CurrentElement => GetActiveElement(t)
     n = GetElementNOFNodes()
     NodeIndexes => CurrentElement % NodeIndexes
     Material => GetMaterial()

     ! Read Geothermal Heat Flux 
     LocalG(1:n) = 0.0_dp
     IF (GeoFlux) THEN
        LocalG(1:n) = ListGetReal( Material, &
                         'Geothermal Heat Flux', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) CALL FATAL(SolverName,'Could not find >Geothermal &
             & Heat Flux< in Material')
     END IF
     LocalDens(1:n) = 0.0_dp
     IF (Melt) THEN
        LocalDens(1:n) = ListGetReal( Material, &
                         'Density', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) CALL FATAL(SolverName,'Could not find >Density< in Material')
     END IF

     !
     LocalVelo = 0.0_dp
     DO i=1, VeloSTDOFs
        LocalVelo(i,1:n) = Velo((VeloSTDOFs)*(VeloPerm(NodeIndexes(1:n))-1) + i)
     END DO

     ! Read the friction law information
     Friction = GetString(Material, 'SSA Friction Law', Found)
     IF (.NOT.Found) &
        CALL FATAL(SolverName,'Could not find Material keyword >SSA Friction Law<')
     SELECT CASE(Friction)
        CASE('linear')
           iFriction = 1
           fm = 1.0_dp
        CASE('weertman')
           iFriction = 2
        CASE('coulomb')
           iFriction = 3
        CASE DEFAULT
        CALL FATAL(SolverName,'Friction should be linear, Weertman or Coulomb')
     END SELECT

! for all friction law
     LocalBeta = 0.0_dp
     LocalBeta(1:n) = GetReal( Material, 'SSA Friction Parameter', Found, CurrentElement)
     IF (.NOT.Found) &
        CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Parameter<')
     IF (iFriction > 1) THEN
        fm = ListGetConstReal( Material, 'SSA Friction Exponent', Found )
        IF (.NOT.Found) &
           CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Exponent<')
     END IF

! only for Coulomb friction
     IF (iFriction > 2) THEN
        PostPeak = ListGetConstReal( Material, 'SSA Friction Post-Peak', Found )
        IF (.NOT.Found) &
           CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Post-Peak<')
        LocalC = 0.0_dp
        LocalC(1:n) = GetReal(Material, 'SSA Friction Maximum Value', Found, CurrentElement)
        IF (.NOT.Found) &
           CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Maximum Value<')
   ! Get the effective pressure
        IF (ASSOCIATED(NSol)) THEN
           LocalN(1:n) = Nval(NPerm(NodeIndexes(1:n)))
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Effective Pressure<')
        END IF
     END IF

     DO j=1,n
        NVelo = LocalVelo(1,j)**2.0
        IF (VeloSTDOFs > 1) NVelo = NVelo + LocalVelo(2,j)**2.0
        NVelo = SQRT(NVelo)
        IF (iFriction < 3) THEN
           tbub = LocalBeta(j) * NVelo**(fm+1.0) 
        ELSE
           IF (PostPeak.NE.1.0_dp) THEN
              alpha = (PostPeak-1.0_dp)**(PostPeak-1.0_dp) / PostPeak**PostPeak
           ELSE
              alpha = 1.0_dp
           END IF
           CN = LocalN(j)*LocalC(j)
           IF (CN > 1.0e-10) THEN
              Xi = NVelo*(LocalBeta(j)/CN)**(1.0/fm)
              tbub = NVelo * CN * (Xi / (1.0 + alpha*Xi**PostPeak))**fm 
           ELSE
              tbub = 0.0_dp
           END IF
        END IF

        ! If Needed, compute the Heat 
        IF (Heat) HeatValue(HeatPerm(NodeIndexes(j))) = tbub
        ! If Needed, compute the melt
        IF (Melt) THEN 
           MeltValue(MeltPerm(NodeIndexes(j))) = &
                           (LocalG(j) + tbub)/(LocalDens(j)*LatentHeat)
        END IF
     END DO
  END DO


!------------------------------------------------------------------------------
END SUBROUTINE SSABasalMelt
!------------------------------------------------------------------------------

