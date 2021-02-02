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
!>   Solver to replace the node exactly on the bedrock                        
SUBROUTINE InterpBedrock( Model, Solver, dt, TransientSimulation )
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
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement, &
                             BoundaryElement

  TYPE(Variable_t), POINTER :: MeshSol, ZbSol, MVSol, GMSol, BedSol

  REAL(KIND=dp), POINTER :: MeshUpdate(:), Zb(:), MV(:), GM(:), Bed(:)

  INTEGER, POINTER :: MeshPerm(:), ZbPerm(:), MVPerm(:), GMPerm(:), BedPerm(:)

  LOGICAL :: Found, CalvingOccurs, NoGroundedMask = .False.

  REAL(KIND=dp) :: x, y, ynew, GroundingLine

  INTEGER :: i, j, k, DIM, STDOFs, LocalNodes   
       
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FlowSolverName

!------------------------------------------------------------------------------

  DIM = CoordinateSystemDimension()

  ! Get variables
  MeshSol => VariableGet( Model % Mesh % Variables, 'Mesh Update' )
  MeshPerm => MeshSol % Perm
  STDOFs =  MeshSol % DOFs
  MeshUpdate => MeshSol % Values

  MVSol => VariableGet( Model % Mesh % Variables, 'Mesh Velocity' )
  MVPerm => MVSol % Perm
  MV => MVSol % Values

  ZbSol => VariableGet( Model % Mesh % Variables, 'Zb' )
  ZbPerm => ZbSol % Perm
  Zb => ZbSol % Values

  BedSol => VariableGet( Model % Mesh % Variables, 'bedrock' )
  BedPerm => BedSol % Perm
  Bed => BedSol % Values

  GMSol => VariableGet( Model % Mesh % Variables, 'GroundedMask' )
    IF ( ASSOCIATED( GMSol ) ) THEN
      GMPerm => GMSol % Perm
      GM => GMSol % Values
    ELSE
      CALL INFO( "InterpBedrock", 'No GroundedMask, considering fully grounded ice')
      NoGroundedMask = .True.
    END IF

  CalvingOccurs = ListGetLogical( Model % Simulation, 'CalvingOccurs', Found)
  IF(.NOT. Found) THEN
     CALL INFO("InterpBedrock","Can't find CalvingOccurs Logical, assuming false")
     CalvingOccurs = .FALSE.
  END IF

  GroundingLine = ListGetConstReal( Model % Simulation, 'GroundingLine', Found)
  IF(.NOT. Found) THEN
     CALL INFO("InterpBedrock","Can't find GroundingLine real")
  END IF

  !Correct vertical nodes positions
  LocalNodes = COUNT( MeshPerm > 0 )
  IF ( LocalNodes <= 0 ) RETURN
  
  DO i=1, Solver % Mesh % NumberOfNodes     
     j = ZbPerm(i)
     IF (j == 0) CYCLE
     k = DIM*MeshPerm(i)
     x = Model % Nodes % x(i)
     y = Model % Nodes % y(i)
     ynew = Bed(BedPerm(i))
!   PRINT*,'y =', y
!   PRINT*,'ynew =', ynew
!   PRINT*,'============='
     ! If there is a calving event, the solver TwoMeshes does not interpolate
     ! the groundedmask and so it is wrong. The right pre-calving value of GM is
     ! stored in GroundingLine ConstReal. Thus, we have to make a double check
     ! in the case there is a calving event :
     IF (NoGroundedMask) THEN
       MeshUpdate(k) = MeshUpdate(k) + (ynew - y) 
       IF (ynew > y) THEN
         MV(k) = MV(k) + (ynew - y)/dt
       ENDIF
       Model % Nodes % y(i) = ynew
     ELSE
       IF ( CalvingOccurs ) THEN
         IF (x .LE. GroundingLine) THEN
           MeshUpdate(k) = MeshUpdate(k) + (ynew - y) 
           IF (ynew > y) THEN
             MV(k) = MV(k) + (ynew - y)/dt
           ENDIF
           Model % Nodes % y(i) = ynew
         ENDIF
       ELSE
       IF (GM(GMPerm(i)) < -0.5_dp) CYCLE ! Not interested in floating nodes
         MeshUpdate(k) = MeshUpdate(k) + (ynew - y) 
         IF (ynew > y) THEN
           MV(k) = MV(k) + (ynew - y)/dt
         ENDIF
         Model % Nodes % y(i) = ynew
       ENDIF
     ENDIF  
  END DO 


!------------------------------------------------------------------------------
END SUBROUTINE InterpBedrock
!------------------------------------------------------------------------------
