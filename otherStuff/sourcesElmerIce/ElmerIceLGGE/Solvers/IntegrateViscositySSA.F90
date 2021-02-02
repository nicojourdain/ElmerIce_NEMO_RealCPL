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
! *  Authors: Julien Brondex
! *
! *  Original Date: 10 Jun 2015
! *
! *****************************************************************************
!> This solver enables to integrate the viscosity with Glen's flow law or Power law   
!> It is also possible to use this solver in order to prescribe directly the Glen's 
!flow law in the sif (as it is done in FullStokes) even if integration is not required 
!(constant T, no damage). In that case (integration not needded) "integration required"
!must be set to False to save computationnal time.                            
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
!******************************************************************************
SUBROUTINE IntegrateViscositySSA( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t),TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
! 
  TYPE(ValueList_t), POINTER :: SolverParams, Material !pointer required to access information in the .sif

!!!!  Usefull variables for element and basic functions
  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t) :: ElementNodes
  INTEGER, POINTER :: NodeIndexes(:)

!!!!! Elmer variables
  TYPE(Variable_t), POINTER :: IntegratedViscosity
  REAL(KIND=dp), POINTER ::  Values(:)
  INTEGER, POINTER :: Perm(:)

!!!! Variables for DetectExtrudedStructure
  TYPE(Solver_t), POINTER :: PSolver  
  TYPE(Variable_t), POINTER :: Var
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),Upointer(:),DownPointer(:)

!!!! Other variables
  REAL(KIND=dp) :: Viscosity
  REAL(KIND=dp) :: alti, altiUp, altiTop, altiBot
  REAL(KIND=dp), ALLOCATABLE :: IntValues(:)
  INTEGER, ALLOCATABLE :: BotNodes(:),TMP(:)
  INTEGER :: NMAX,k,n,i,C,l,j
  INTEGER :: mat
  INTEGER :: DIM, Up

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  LOGICAL, ALLOCATABLE :: VisitedNodes(:)  
  LOGICAL ::  Firsttime=.true.
  LOGICAL ::  GotIt
  LOGICAL ::  Integrate, Top


  SAVE Firsttime,DIM
  SAVE ElementNodes
  SAVE SolverName
  SAVE Material
  SAVE BotNodes 
  SAVE C, mat
  SAVE Integrate, Top
  SAVE TopPointer,BotPointer,Upointer,DownPointer
  SAVE IntValues

  !!! Where the variable Integrated viscosity will be stored
  IntegratedViscosity => Solver % Variable
  Values => IntegratedViscosity % Values
  Perm => IntegratedViscosity % Perm

!-------------------------------------------------------------
! Initialization
!-------------------------------------------------------------

  IF (Firsttime) THEN
    WRITE(SolverName, '(A)') 'IntegrateViscositySSA'
    DIM = CoordinateSystemDimension()

! Get the number of the body material

    Element => GetActiveElement(1)
    Material => GetMaterial (Element)
    !!!! 
    PSolver => Solver
    CALL DetectExtrudedStructure(Solver % Mesh , PSolver, Var, &
         TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
         UpNodePointer = Upointer , DownNodePointer = DownPointer)

    NMAX=Solver % Mesh % NumberOfNodes
      
  
    ALLOCATE(VisitedNodes(NMAX),TMP(NMAX),IntValues(NMAX))
 
  
    VisitedNodes=.FALSE.
    C=0.0
    DoBoundaryElements: DO k=1,Solver % Mesh % NumberOfBoundaryElements
              Element => GetBoundaryElement(k)
              n = GetElementNOFNodes()
              NodeIndexes => Element % NodeIndexes
  
              DoNodes: DO i=1,n
                  IF ( VisitedNodes(NodeIndexes(i)) ) CYCLE
                       VisitedNodes(NodeIndexes(i))=.true. 
                  IF ( BotPointer(NodeIndexes(i)) /= NodeIndexes(i) ) CYCLE
                    C=C+1
                    TMP(C)=NodeIndexes(i)
                  END DO DoNodes
    END DO DoBoundaryElements

    ALLOCATE (BotNodes(C))
      BotNodes(1:C)=TMP(1:C)
    DEALLOCATE (TMP, VisitedNodes)

!Check if integration is required (not the case if viscosity uniform) 
    SolverParams => GetSolverParams()

    Integrate = ListGetLogical ( SolverParams,'Integration Required', GotIt)
    IF (.NOT. GotIt) THEN
      Integrate = .TRUE.
      WRITE(Message, '(A)') 'Integration Required set to True'
      CALL INFO(SolverName, Message, Level = 20 )
    END IF 

! Check if we want the integrated viscosity at the top nodes (at the bottom nodes by default)
    Top = ListGetLogical ( SolverParams,'Integrated Value Top', GotIt)
    IF (.NOT. GotIt) THEN
      Top = .FALSE.
      WRITE(Message, '(A)') 'Integrated Value Top set to False. Mean viscosity given at bottom nodes'
      CALL INFO(SolverName, Message, Level = 20 )
    END IF 
!!!End of First time visit
     Firsttime=.FALSE.
  END IF

!-------------------------------------------------------------
!End of initialization
!-------------------------------------------------------------
      

  Values = 0.0_dp



  DO l=1,C !!Loop on the bottom nodes
    IF (.NOT. Integrate) THEN    
      IF (.NOT. Top) THEN
        Values(Perm(BotNodes(l)))=Visco(BotNodes(l))
      ELSE 
        Values(Perm(TopPointer(BotNodes(l))))=Visco(BotNodes(l))
      END IF
    ELSE
      IntValues(BotNodes(l))= 0.0_dp
      Up = BotNodes(l)
      DO WHILE (Up /= TopPointer(BotNodes(l))) !Loop on the nodes column
        IF (DIM .eq. 2) THEN
          alti = Model % Nodes % y (Up)
          altiUp = Model % Nodes % y (Upointer(Up))
        ELSE
          alti = Model % Nodes % z (Up)
          altiUp = Model % Nodes % z (Upointer(Up))
        END IF
        IntValues (Upointer(Up))= 0.5 * (Visco(Up)+Visco(Upointer(Up)))*(altiUp-alti) + IntValues(Up)
        Up = Upointer(Up)
      END DO !End of the loop on the nodes column
        IF (DIM .eq. 2) THEN
          altiBot = Model % Nodes % y (BotNodes(l))
          altiTop = Model % Nodes % y (TopPointer(BotNodes(l)))
        ELSE
          altiBot = Model % Nodes % z (BotNodes(l))
          altiTop = Model % Nodes % z (TopPointer(BotNodes(l)))
        END IF
      IF (.NOT. Top) THEN
          Values(Perm(BotNodes(l)))= IntValues (TopPointer(BotNodes(l))) / (altiTop - altiBot)
      ELSE !If we want the mean viscosity at the top nodes 
          Values(Perm(TopPointer(BotNodes(l))))= IntValues (TopPointer(BotNodes(l))) / (altiTop - altiBot)
      END IF
    END IF
  END DO
!  PRINT*,'fin loop nodes'
!!*****************************************************************************
!END SUBROUTINE IntegrateViscositySSA
!*****************************************************************************
  CONTAINS
  FUNCTION Visco(nodenumber) RESULT(eta0)
    !------------------------------------------
  !   USE DefUtils
     IMPLICIT NONE
    !------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: nodenumber
     REAL(KIND=dp)  :: eta0
    !------------------------------------------

     TYPE(Element_t),POINTER :: Element
     TYPE(ValueList_t), POINTER :: Material

    !------------------------------------------------------------------------------

     LOGICAL :: GotIt

     CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityFlag, TemperatureName
     REAL(KIND=dp) :: c1, c4, Temp, Tlimit, &
          A1, A2, Q1, Q2, R, EhF, ArrheniusFactor

     ! Temperature is needed for thermal models
     TYPE(Variable_t), POINTER :: TempSol 
     REAL(KIND=dp), POINTER :: Temperature(:)
     INTEGER, POINTER :: TempPerm(:)

     LOGICAL :: SetArrheniusFactor=.FALSE.
     
     eta0 = 0.0_dp

    Element => GetActiveElement(1)
    Material => GetMaterial (Element)

     ViscosityFlag = ListGetString( Material,'Viscosity Model', GotIt)

     IF(.NOT. GotIt) RETURN


     SELECT CASE( ViscosityFlag )

     CASE('glen')
        c1 = ListGetRealAtNode( Material, 'Glen Exponent', nodenumber, GotIt ) ! this is the real exponent, n, not 1/n
        IF (.NOT.GotIt) c1 = 3.0_dp
        
        SetArrheniusFactor = GetLogical(Material, 'Set Arrhenius Factor', GotIt)
        IF ( (.NOT.GotIt) .OR. .NOT.(SetArrheniusFactor)) THEN
           Temp = ListGetRealAtNode(Material, 'Constant Temperature', nodenumber, GotIt) !we are happy as is
           IF(.NOT.GotIt) THEN !we have to find a temperature field

              TemperatureName = GetString(Material, 'Temperature Field Variable', GotIt)
              IF (.NOT.GotIt) WRITE(TemperatureName,'(A)') 'Temperature'
              TempSol => VariableGet( CurrentModel % Variables,TRIM(TemperatureName))
              IF ( ASSOCIATED( TempSol) ) THEN
                 TempPerm    => TempSol % Perm
                 Temperature => TempSol % Values   
                 Temp =  Temperature(TempPerm(nodenumber))
              ELSE
                 WRITE(Message, '(A,A,A)') 'Could not find variable ',&
                      TRIM(TemperatureName),' to inquire temperatur field for Glen'
                 CALL FATAL('IntegrateViscositySSA',Message)
              END IF
           END IF
        
           R = GetConstReal( CurrentModel % Constants,'Gas Constant',GotIt)
           IF (.NOT.GotIt) R = 8.314_dp
           ! lets for the time being have this hardcoded
           Tlimit = GetConstReal(Material, 'Limit Temperature', GotIt)
           IF (.NOT.GotIt) THEN
              Tlimit = -10.0_dp
              CALL INFO('IntegrateViscositySSA','Limit Temperature not found. Setting to -10', Level=5)
           END IF
           A1 = GetConstReal(Material, 'Rate Factor 1', GotIt)
           IF (.NOT.GotIt) THEN
              A1 = 3.985d-13
              CALL INFO('IntegrateViscositySSA','Rate Factor 1 not found. Setting to 3.985e-13', Level=5)
           END IF
           A2 = GetConstReal(Material, 'Rate Factor 2', GotIt)
           IF (.NOT.GotIt) THEN
              A2 = 1.916d03
              CALL INFO('IntegrateViscositySSA','Rate Factor 2 not found. Setting to 1.916E03', Level=5)
           END IF
           Q1 = GetConstReal(Material, 'Activation Energy 1', GotIt)
           IF (.NOT.GotIt) THEN
              Q1 = 60.0d03
              CALL INFO('IntegrateViscositySSA','Activation Energy 1 not found. Setting to 60.0E03', Level=5)
           END IF
           Q2 = GetConstReal(Material, 'Activation Energy 2', GotIt)
           IF (.NOT.GotIt) THEN
              Q2 = 139.0d03
              CALL INFO('IntegrateViscositySSA','Activation Energy 2 not found. Setting to 139.0d03', Level=5)
           END IF
        
           IF (Temp.LE. Tlimit) THEN
              ArrheniusFactor = A1 * EXP( -Q1/(R * (273.15 + Temp)))
           ELSE IF((Tlimit<Temp) .AND. (Temp .LE. 0.0_dp)) THEN
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15 + Temp)))
           ELSE
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15)))
              CALL INFO('IntegrateViscositySSA','Positive Temperature detected in Glen - limiting to zero!', Level = 5)
           END IF
        ELSE
           ArrheniusFactor = GetConstReal(Material,'Arrhenius Factor', GotIt)
           IF (.NOT.(GotIt)) THEN 
              CALL FATAL('IntegrateViscositySSA','<Set Arrhenius Factor> is TRUE, but no value <Arrhenius Factor> found')
           END IF
        END IF

        EhF =  ListGetRealAtNode( Material, 'Glen Enhancement Factor', nodenumber, GotIt )
        IF (.NOT.GotIt) EhF = 1.0_dp
       
        ! compute the viscosity eta0
        eta0 = (2 * EhF * ArrheniusFactor)**(-1.0_dp/c1);

     CASE('power law')
           
        eta0 = ListGetRealAtNode(Material, 'Viscosity', nodenumber, GotIt) 
        IF(.NOT.GotIt) THEN 
           WRITE(Message,'(A)')'Variable Viscosity not found. Setting to 1.0'
           CALL INFO('IntegrateViscositySSA',Message, Level = 20)
           eta0 = 1.0_dp
        END IF
        
     CASE DEFAULT
        CALL WARN('IntegrateViscositySSA','Unknown material model')

     END SELECT

  END FUNCTION Visco
!******************************************************************************
END SUBROUTINE IntegrateViscositySSA
!******************************************************************************
