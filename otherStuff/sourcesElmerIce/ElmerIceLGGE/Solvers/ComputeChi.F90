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
! ******************************************************************************
! *
! *  Compute the Chi value for the Damage model                          
! *  Only where Chi > 0 ice is damaged               
! *
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini
! *  Email:   gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://www.csc.fi/elmer
! *
! *  Original Date: 14 Dec 2012 
! * 
! *****************************************************************************


! *****************************************************************************
SUBROUTINE ComputeChi( Model,Solver,dt,TransientSimulation )
!DEC$ATTRIBUTES DLLEXPORT :: ComputeChi
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
  USE GeneralUtils

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
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable

  TYPE(Variable_t), POINTER :: StressVariable, FlowVariable, DamageVariable
  REAL(KIND=dp), POINTER :: StressValues(:), FlowValues(:), DamageValues(:)
  INTEGER, POINTER :: StressPerm(:), FlowPerm(:), DamagePerm(:)
  INTEGER :: Ind(3,3), DIM 
  REAL(KIND=dp) :: Sig(3,3), SigDev(3,3), EigVect(3,3), EigValues(3), tmp
  REAL(KIND=dp) :: SigmaI, SigmaII, D, Chi, B, sigmath, lambdah, r 
  LOGICAL :: GotIt, FirstTime = .TRUE., Cauchy

  REAL(KIND=DP) :: EI(3),Dumy(1),Work(24)
  INTEGER :: infor 

  INTEGER :: i, j, ni, n, t, NodeNumber, ord(3)
  INTEGER, POINTER :: Permutation(:) 

  REAL(KIND=dp), POINTER :: VariableValues(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
       

  SAVE SolverName
  SAVE :: Ind, DIM, sigmath, lambdah, B
  SAVE :: FirstTime, Cauchy
!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  WRITE(SolverName, '(A)') 'ComputeChi'


  IF (FirstTime) THEN
      FirstTime = .FALSE.  
      DIM = CoordinateSystemDimension()

      DO i=1, 3
         Ind(i,i) = i
      END DO
      Ind(1,2) = 4
      Ind(2,1) = 4
      Ind(2,3) = 5
      Ind(3,2) = 5
      Ind(3,1) = 6
      Ind(1,3) = 6

   Material => GetMaterial()
   !Read the coefficients B, sigmath, lambdah, r  

      B = GetConstReal( Material, 'Damage Enhancement Factor', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Damage Enhancement Factor not found. &
              &Setting to 0.1'
         CALL INFO('Damage Source', Message, level=2)
         B = 0.1_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Enhancement Factor = ', B
         CALL INFO('Damage Source', Message, level=2)
      ENDIF

      sigmath = GetConstReal( Material, 'Damage Parameter sigmath', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Damage Parameter sigmath not found. &
              &Setting to 0.1'
         CALL INFO('Damage Source', Message, level=2)
         sigmath = 0.044_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Parameter sigmath = ', sigmath 
         CALL INFO('Damage Source', Message, level=2)
      ENDIF

      lambdah = GetConstReal( Material, 'Damage Parameter lambdah', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Damage Parameter lambdah not found. &
              &Setting to 0.4'
         CALL INFO('Damage Source', Message, level=2)
         lambdah = 0.4_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Parameter lambdah = ', lambdah 
         CALL INFO('Damage Source', Message, level=2)
      ENDIF
   
      r = GetConstReal( Material, 'Damage Criterion Exponent', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Damage Criterion Exponent not found. &
              &Setting to 0.43'
         CALL INFO('Damage Source', Message, level=2)
         r = 0.43_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Criterion Exponent = ', r
         CALL INFO('Damage Source', Message, level=2)
      ENDIF

! Cauchy or deviatoric stresses ?
      Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )
      write(*,*)'Cauchy',Cauchy
   END IF


   ! Get the Stress                     
   StressVariable => VariableGet( Model % Variables, 'Stress' )
   IF ( ASSOCIATED( StressVariable ) ) THEN
      StressPerm    => StressVariable % Perm
      StressValues  => StressVariable % Values
   ELSE
      CALL FATAL(SolverName, 'Need ComputeDevStress Solver, Stress not associated !!!')
   END IF

   ! Get the variables to compute the hydrostatic pressure  
   FlowVariable => VariableGet( Model % Variables, 'Flow Solution' )
   IF ( ASSOCIATED( FlowVariable ) ) THEN
      FlowPerm    => FlowVariable % Perm
      FlowValues  => FlowVariable % Values
   ELSE
      CALL FATAL(SolverName, 'Need NS Solver, Flow Solution not associated !!!')
   END IF

   ! Get the Damage variable
   DamageVariable => VariableGet( Model % Variables, 'Damage' )
   IF ( ASSOCIATED( DamageVariable ) ) THEN
      DamagePerm    => DamageVariable % Perm
      DamageValues  => DamageVariable % Values
   ELSE
      CALL FATAL(SolverName, 'Need Damage Solver, Damage Variable not associated !!!')
   END IF


! No non-linear iteration, no time dependency  
  VariableValues = -9999.0_dp

  ! Just a loop
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()
     
     ! Loop over all nodes of the element t
     DO ni  = 1, n
        ! Test if we already compute Chi for this node
        NodeNumber = Element % NodeIndexes(ni)
        IF (VariableValues(Permutation(Nodenumber)).NE.-9999.0_dp) CYCLE

        Sig = 0.0
        DO i=1, DIM
           DO j= 1, DIM
              Sig(i,j) =  &
                      StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(i,j) )
           END DO
        END DO
        IF (DIM==2) Sig(3,3) = StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(3,3))

  ! S = sigma - p
  ! Need Cauchy Stress and Deviatoric Stress 
         IF (.NOT.Cauchy) THEN ! If Deviatoric Stress is computed, then, get the
                         ! Cauchy Stress
            SigDev = Sig
            DO i=1,3  
               Sig(i,i) = SigDev(i,i) - FlowValues((DIM+1)*FlowPerm(Nodenumber))
            END DO
         ELSE ! If the Cauchy Stress is computed, then get the Deviatoric Stress 
            DO i=1,3  
               SigDev(i,i) = Sig(i,i) + FlowValues((DIM+1)*FlowPerm(Nodenumber))
            END DO
         END IF

! Compute the maximum principal stress and the second invariant of deviatoric stress

   ! Get the principal stresses
     CALL DGEEV('N','N',3,Sig,3,EigValues,EI,Dumy,1,Dumy,1,Work,24,infor )
     IF (infor.ne.0) &
         CALL FATAL(SolverName, 'Failed to compute EigenValues') 

   ! Get ordered the EigenValues      
     ord = (/(i,i=1,3)/)
     CALL SortD(3,EigValues,ord)
     
     SigmaI = EigValues(3)
     SigmaII = EigValues(2)

     !SigmaI = MAXVAL(EigValues)

     D = DamageValues(DamagePerm(NodeNumber))

     Chi = 1.0_dp / (1.0_dp - D) * (SigmaI) - sigmath
   
     VariableValues(Permutation(NodeNumber)) = Chi

    END DO ! nodes ni
  END DO   ! element t
  
!------------------------------------------------------------------------------
END SUBROUTINE ComputeChi
!------------------------------------------------------------------------------


