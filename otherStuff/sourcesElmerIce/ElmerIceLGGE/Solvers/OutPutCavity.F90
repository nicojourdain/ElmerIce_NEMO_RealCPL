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
! *  Original Date: 12 July 2010
! * 
! *****************************************************************************
!>   Solver to produce scalar outputs for the cavity problem                   
SUBROUTINE OutPutCavity( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!   Replace the node on the bedrock 
!   This is only valid for very small adjustments 
!   Only the nodes on the bed are moved
!   Mesh Update and Mesh Velocity are updated accordingly
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

  TYPE(Variable_t), POINTER :: residuSol, MaskSol, BedSol
  TYPE(Variable_t), POINTER :: TimeVar
  REAL(KIND=dp), POINTER :: residu(:), Mask(:), Bed(:)
  INTEGER, POINTER :: ResiduPerm(:), MaskPerm(:), BedPerm(:)

  LOGICAL :: Found, AllocationsDone = .FALSE., FirstTime = .TRUE.
  LOGICAL, ALLOCATABLE :: Cavity(:)

  REAL(KIND=dp) :: x, y, b, hmean, CavityLength, x1, x2, &
       TotalForce(2)

  INTEGER :: i, j, k, t, n, M, DIM, STDOFs, NodeOnBed, NodeInCavity, istat
  INTEGER, POINTER :: NodeIndexes(:)
       
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'OutPutCavity', FileName

  SAVE :: DIM, Cavity, FirstTime, FileName
!------------------------------------------------------------------------------

  IF (FirstTime) THEN 
     FirstTime = .FALSE.
     FileName = GetString( Solver % Values, 'FileName', Found )
     IF (.Not.Found) CALL FATAL(SolverName, 'FileName not given')
  END IF

  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(cavity)
     ALLOCATE(cavity(M), STAT=istat )
     IF ( istat /= 0 ) &
        CALL FATAL( SolverName, 'Memory allocation error.' )
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

  ResiduSol => VariableGet( Model % Mesh % Variables, 'Flow Solution Loads' )
  IF (ASSOCIATED(ResiduSol)) THEN
     ResiduPerm => ResiduSol % Perm
     Residu => ResiduSol % Values
  ELSE
     CALL FATAL(SolverName, 'No > Flow Solution Loads < associated')
  END IF

  MaskSol => VariableGet( Model % Mesh % Variables, 'GroundedMask' )
  IF (ASSOCIATED(MaskSol)) THEN
     MaskPerm => MaskSol % Perm
     Mask => MaskSol % Values
  ELSE
     CALL FATAL(SolverName, 'No > GroundedMask < associated')
  END IF

  BedSol => VariableGet( Model % Mesh % Variables, 'bedrock' )
  IF (ASSOCIATED(BedSol)) THEN
     BedPerm => BedSol % Perm
     Bed => BedSol % Values
  ELSE
     CALL FATAL(SolverName, 'No > Bedrock < associated')
  END IF
  
  Cavity = .TRUE.   
  TotalForce = 0.0_dp
  NodeOnBed = 0.0_dp
  CavityLength = 0
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t) 
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes
      
     Do i = 1, n
        k = NodeIndexes(i)
        ! we cycle nodes in the cavity
        IF (Mask(MaskPerm(k)) < 0.0_dp) CYCLE
        ! First time we meet that node
        IF (Cavity(k)) THEN
           Cavity(k) = .FALSE. 
           j = (DIM+1)*ResiduPerm(k)
           TotalForce(1) = TotalForce(1) + Residu(j+1)
           TotalForce(2) = TotalForce(2) + Residu(j+2)
           NodeOnBed = NodeOnBed + 1
        END IF
     END DO
     
     IF (ANY(Cavity(NodeIndexes(1:n)))) THEN
        x1 = Model % Nodes % x(NodeIndexes(1))
        x2 = Model % Nodes % x(NodeIndexes(2))
        CavityLength = CavityLength + ABS(x2 - x1) 
     END IF


  END DO 

  ! Compute the mean elevation of the cavity
  hmean = 0.0_dp
  NodeInCavity = 0
  DO i=1, Solver % Mesh % NumberOfNodes     
     IF (MaskPerm(i) == 0) CYCLE
     IF (.Not.Cavity(i)) CYCLE
     y = Model % Nodes % y(i)
     b = bed(bedPerm(i)) 
     hmean = hmean + (y-b)
     NodeInCavity = NodeInCavity + 1
  END DO 
  hmean = hmean / NodeInCavity


  TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
  OPEN (10, FILE=TRIM(FileName), POSITION='APPEND')
  WRITE(10, '(5(e14.8,2x),i4,2x,i4)')TimeVar % values(1), TotalForce(1), TotalForce(2), &
      CavityLength, hmean, NodeOnBed, NodeInCavity 
  CLOSE(10)    


!------------------------------------------------------------------------------
END SUBROUTINE OutPutCavity
!------------------------------------------------------------------------------


