!/*****************************************************************************/
! *
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
!/*****************************************************************************/
! *
! *  SmoothFS Solver : A solver to smooth the free surface oscillations given
! *                   a chararteristic period
! *  No Variable needed
! *
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini
! *  Email:   gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://www.csc.fi/elmer
! *
! *  Original Date: 28 August 2009 
! * 
! *****************************************************************************
SUBROUTINE SmoothFSSolver( Model,Solver,dt,TransientSimulation )
!******************************************************************************
!
!  Smooth the free surface solution ...
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
  TYPE(Element_t),POINTER :: CurrentElement, Element
  TYPE(ValueList_t), POINTER :: SolverParams, Material
  TYPE(Variable_t), POINTER :: FreeSurfSol 
  TYPE(Nodes_t)   :: ElementNodes

  LOGICAL :: AllocationsDone = .FALSE., Found

  INTEGER :: i, j, k, l, n, m, t, istat, DIM, NbPt 
  INTEGER, POINTER :: FreeSurfPerm(:)

  REAL(KIND=dp), POINTER :: FreeSurfValue(:), PrevFreeSurf(:,:) 

  REAL(KIND=dp) :: ca, cb, cc, x0, x, y, Deno, dx 
  REAL(KIND=dp) :: SX, SX2, SX3, SX4, SY, SYX, SYX2

  CHARACTER(LEN=MAX_NAME_LEN) :: FreeSurfVarName

  SAVE ElementNodes
       


!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
        DIM = CoordinateSystemDimension()
        SolverParams => GetSolverParams()

        FreeSurfVarName =  GetString( SolverParams,'Free Surface Name', Found)
        IF(.NOT.Found) THEN        
           CALL WARN('SmoothFS','Keyword >Free Surface Name< not found in section >Equation<')
           CALL WARN('Smooth','Taking default value >Zs<')
           WRITE(FreeSurfVarName,'(A)') 'Zs'
        END IF
        FreeSurfSol => VariableGet( Solver % Mesh % Variables, freeSurfVarName )
        IF ( ASSOCIATED( FreeSurfSol ) ) THEN
           FreeSurfPerm     => FreeSurfSol % Perm
           FreeSurfValue => FreeSurfSol % Values
           PrevFreeSurf => FreeSurfSol % PrevValues
        ELSE
            WRITE(Message,'(A,A,A)') &
                 'SmoothFS >',FreeSurfVarName,'< not found'
            CALL FATAL('SmoothFS',Message)              
        END IF

        dx = GetConstReal( SolverParams, 'Smooth FS Length', Found)
        IF(.NOT. Found) dx = 1.0e4_dp 
        WRITE(Message,'(a,F8.2)') 'Smooth FS Length=', dx
        CALL INFO('SmoothFS', Message, Level=4)


  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
! IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
!    WRITE(SolverName, '(A)') 'SmoothFS'
!    N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
!    M = Model % Mesh % NumberOfNodes
!    IF (AllocationsDone) DEALLOCATE()
!
!    ALLOCATE( FORCE(N) )
!    IF ( istat /= 0 ) THEN
!       CALL Fatal( SolverName, 'Memory allocation error.' )
!    END IF
!
!
!    AllocationsDone = .TRUE.
!    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
! END IF


! No non-linear iteration, no time dependency  
! No assembly, no system to be solved

     DO i=1,Model % Mesh % NumberOfNodes
        k = FreeSurfPerm(i)           
        IF (k == 0) CYCLE
        x0 = Solver % Mesh % Nodes % x (i)
        SX = 0.0
        SX2 = 0.0
        SX3 = 0.0
        SX4 = 0.0
        SY = 0.0
        SYX = 0.0
        SYX2 = 0.0
        NbPt = 0

        DO j=1,Model % Mesh % NumberOfNodes
           l = FreeSurfPerm(j)           
           IF (l == 0) CYCLE
           x = Solver % Mesh % Nodes % x (j)
           y = FreeSurfValue(l)
           IF (ABS(x0-x)>dx) CYCLE
           NbPt = NbPt+1
           SX = SX + x
           SX2 = SX2 + x**2.0         
           SX3 = SX3 + x**3.0
           SX4 = SX4 + x**4.0
           SY = SY + y        
           SYX = SYX + y * x  
           SYX2 = SYX2 + y * x**2.0
        END DO 
        Deno = SX4 * (NbPt * SX2 - SX**2.0) + SX3 * (SX * SX2 - NbPt * SX3) +&
             SX2 * (SX * SX3 - SX2**2.0)
        cc = (SX4 * (SX2 * SY - SX * SYX) + SX3 * (SX * SYX2 - SX3 * SY) +&
               SX2 * (SX3 * SYX - SX2 * SYX2)) / Deno
        cb = (SX4 * (NbPt * SYX - SY * SX) + SX3 * (SY * SX2 - NbPt * SYX2) +&
                 SX2 * (SYX2 * SX - SYX * SX2)) / Deno
        ca = (SYX2 * (NbPt * SX2 - SX**2.0) + SYX * (SX * SX2 - NbPt * SX3) +&
                SY * (SX * SX3 - SX2**2.0)) / Deno
! New FreeSurfValue for point i
!       write(*,*)j,x0,FreeSurfValue(k), ca * x0**2.0 + cb * x0 + cc
        FreeSurfValue(k) = ca * x0**2.0 + cb * x0 + cc
     END DO
      



!------------------------------------------------------------------------------
END SUBROUTINE SmoothFSSolver
!------------------------------------------------------------------------------



