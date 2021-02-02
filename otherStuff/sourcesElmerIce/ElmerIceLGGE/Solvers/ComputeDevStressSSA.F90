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
! *  Authors: Juha Ruokolainen, Fabien Gillet-Chaulet, Olivier Gagliardini
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 08 Jun 1997
! *  Date of modification: 13/10/05 from version 1.5
! *
! *****************************************************************************
!> Module containing a solver for computing deviatoric or Cauchy   
!>          stress from flow solution (Restricted to NS solver)    
!> 2D SDOFs = 4 (S11, S22, S33, S12)                               
!> 3D SDOFs = 6 (S11, S22, S33, S12, S23, S31)                     
!> Keywords : Cauchy (Logical),                                    
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE ComputeDevStressSSA( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  
  USE DefUtils
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve stress equations for one timestep
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations (NOTE: Not used
!            currently)
!
!******************************************************************************
  
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: PSolver
  
  TYPE(Matrix_t),POINTER :: StiffMatrix
  
  INTEGER :: i, j, k, l, n, t, NDeg, M
  INTEGER :: dim, STDOFs, StressDOFs, LocalNodes, istat
  
  TYPE(ValueList_t),POINTER :: Material, BC, BodyForce
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  
  REAL(KIND=dp) :: UNorm
  
  REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
  REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:)
  REAL(KIND=dp) :: u, v, w, detJ
  REAL(KIND=dp), ALLOCATABLE :: Density(:)
  REAL(KIND=dp) :: rho, rhoM, Gravity, alti, altiM, altiTop, altiTopM, g, Pcryo
  
  LOGICAL :: stat, CSymmetry 
  
  TYPE(Variable_t), POINTER :: StressSol,Var, FlowVariable 
  TYPE(Mesh_t), POINTER :: Mesh 
  
  REAL(KIND=dp), POINTER ::  Stress(:), Solution(:), &
       ForceVector(:), FlowValues(:) 
  
  INTEGER, POINTER :: StressPerm(:), NodeIndexes(:), &
       FlowPerm(:)
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),UpPointer(:)
  
  INTEGER :: body_id
  INTEGER :: old_body = -1
  
  LOGICAL :: Isotropic, AllocationsDone = .FALSE.,  &
       Requal0, StiffMatrixAssemblyDone = .FALSE.
  LOGICAL :: GotIt,  Cauchy = .FALSE.
  
  REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LocalForce(:), &
       LocalP(:), LocalVelo(:,:)
  
  INTEGER :: NumberOfBoundaryNodes, COMP
  INTEGER, POINTER :: BoundaryReorder(:)
  
  REAL(KIND=dp), POINTER :: BoundaryNormals(:,:), &
       BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName, StressSolverName
  
  
!------------------------------------------------------------------------------
  SAVE NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals, &
       BoundaryTangent1, BoundaryTangent2
  
  SAVE Basis, dBasisdx, ddBasisddx
  SAVE TopPointer, BotPointer, UpPointer
  SAVE Gravity, Density, g
  SAVE LocalMassMatrix, LocalStiffMatrix, LocalForce, &
       ElementNodes,  &
       AllocationsDone,  &
       old_body, &
       StiffMatrixAssemblyDone, &
       Cauchy
  
  SAVE LocalVelo, LocalP, dim
  
!  NULLIFY(StressSol, FlowVariable)

  IF ( CurrentCoordinateSystem() .NE. Cartesian) THEN
    CALL FATAL('ComputeDevStressSSA', &
         & 'Cartesian coordinates required')
  END IF

  Mesh => Model % Mesh 

!------------------------------------------------------------------------------
!  Read the name of the SSA Bulk Velocities Solver
!------------------------------------------------------------------------------
  FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', GotIt )    
  IF (.NOT.Gotit) FlowSolverName = 'SSABulkVelocity'
  FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName )
  IF ( ASSOCIATED( FlowVariable ) ) THEN
     FlowPerm    => FlowVariable % Perm
     FlowValues  => FlowVariable % Values
  ELSE
     CALL FATAL('ComputeDevStressSSA', &
          & 'No variable for velocity associated. SSABulkVelocities solver required')
  END IF

!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  
  Solution => Solver % Variable % Values
  STDOFs   =  Solver % Variable % DOFs
  
  IF ( STDOFs /=1 ) THEN
     CALL Fatal( 'ComputeDevStressSSA', 'DOF must be equal to 1' )
  END IF
  
  StressSolverName = GetString( Solver % Values, 'Stress Variable Name', GotIt )    
  IF (.NOT.Gotit) CALL FATAL('ComputeDevStressSSA', & 
       'Stress Variable Name not defined')
  
  StressSol => VariableGet( Solver % Mesh % Variables, TRIM(StressSolverName) )
  IF ( .NOT. ASSOCIATED(StressSol) ) THEN
     WRITE(Message,'(A,A)') 'Could not find Stress Variable: ', StressSolverName
     CALL FATAL('ComputeDevStressSSA',Message)
  END IF
  StressPerm => StressSol % Perm
  StressDOFs = StressSol % DOFs
  Stress => StressSol % Values
  
  dim = CoordinateSystemDimension()
  IF (StressDOfs /= 2*dim) THEN
     CALL Fatal( 'ComputeDevStressSSA', 'Bad dimension of Stress Variable (4 in 2D, 6 in 3D)' )
  ENDIF
  
  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS
!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
     N = Model % MaxElementNodes
     M = Model % Mesh % NumberOfNodes
   
     PSolver => Solver
     CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, BotNodePointer = BotPointer, &
                                 TopNodePointer = TopPointer, UpNodePointer = UpPointer)
   
!!! Get the gravity
     CurrentElement => GetActiveElement(1)
     BodyForce => GetBodyForce(CurrentElement)

     Gravity = 0._dp
     IF (dim ==2) THEN
        Gravity = GetConstReal(BodyForce, 'Flow BodyForce 2', GotIt)
        IF (.NOT. GotIt) THEN
            CALL Fatal( 'ComputeDevStressSSA','Gravity must be vertical and prescribed in Flow BodyForce 2')
        END IF
     ELSE
        Gravity = GetConstReal(BodyForce, 'Flow BodyForce 3', GotIt)
        IF (.NOT. GotIt) THEN
            CALL Fatal( 'ComputeDevStressSSA','Gravity must be vertical and prescribed in Flow BodyForce 3')
        END IF
     END IF
  
     g = abs(Gravity)

     IF ( AllocationsDone ) THEN
        DEALLOCATE( ElementNodes % x,     &
             ElementNodes % y,     &
             ElementNodes % z,     &
             LocalVelo, LocalP,    &                      
             Density,              &
             Basis, ddBasisddx,    &
             dBasisdx,             &
             LocalMassMatrix,      &
             LocalStiffMatrix,     &
             LocalForce )
     END IF

     ALLOCATE( ElementNodes % x( N ), &
          ElementNodes % y( N ), &
          ElementNodes % z( N ), &
          LocalVelo( 3,N ), LocalP( N ), &
          Density( N ), &                                     
          Basis( 2*N ),ddBasisddx(1,1,1), dBasisdx( 2*N,3 ), &
          LocalMassMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
          LocalStiffMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
          LocalForce( 2*STDOFs*N ),  &
          STAT=istat)

     StiffMatrixAssemblyDone=.FALSE.

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'ComputeDevStressSSA', 'Memory allocation error.' )
     END IF
!------------------------------------------------------------------------------

     AllocationsDone = .TRUE.
  END IF



!------------------------------------------------------------------------------
!!!!!Loop for each stress
     DO COMP = 1, 2*dim

        WRITE(Message,'(a,i3)' ) ' Component : ', COMP  
        CALL Info( 'ComputeDevStressSSA', Message, Level=5 )

        StiffMatrixAssemblyDone=StiffMatrixAssemblyDone.AND.ASSOCIATED(Solver % Matrix % BulkValues)
        IF (StiffMatrixAssemblyDone) THEN
           Solver % Matrix % Values = Solver % Matrix % BulkValues
           ForceVector = 0._dp
        ELSE
   CALL DefaultInitialize()
END IF
!------------------------------------------------------------------------------
        DO t=1,Solver % NumberOFActiveElements

           CurrentElement => GetActiveElement(t)
           n = GetElementNOFNodes()
           NodeIndexes => CurrentElement % NodeIndexes

           ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
           ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
           ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

!------------------------------------------------------------------------------
!    Read in material constants from Material section
!------------------------------------------------------------------------------
           Material => GetMaterial(CurrentElement)
           
           Density = 0._dp
           Density(1:n) = ListGetReal (Material, 'SSA Mean Density', n, NodeIndexes, GotIt)
 
           Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )
           IF (.NOT.Gotit) THEN
              Cauchy = .FALSE.
              WRITE(Message,'(A)') 'Cauchy set to False'
              CALL INFO('ComputeDevStressSSA', Message, Level = 20)
           END IF
           
           LocalVelo = 0.0_dp
           DO i=1, dim
              LocalVelo(i,1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + i)
           END DO

!!! Loop on nodes to calculate P cryostatic (/!\ Different from SSA pressure)
           DO j=1,n 
              rho = ListGetRealAtNode( Material, 'Density', NodeIndexes(j), GotIt)
              IF (.NOT.GotIt) rho=Density(j)

              IF (dim==2) THEN
                 alti = Model % Nodes % y (NodeIndexes(j))
                 altiTop = Model % Nodes % y (TopPointer(NodeIndexes(j)))
              ELSE
                 alti = Model % Nodes % z (NodeIndexes(j))
                 altiTop = Model % Nodes % z (TopPointer(NodeIndexes(j)))
              END IF

              LocalP(j) = rho * g * (altiTop-alti) ! Cryostatic pressure not the real one  
            END DO 
          
           CALL LocalNSMatrix(COMP, LocalMassMatrix, LocalStiffMatrix, &
                LocalForce,  LocalVelo, LocalP, &
                CurrentElement, n, &
                ElementNodes, Cauchy)
!------------------------------------------------------------------------------
!        Update global stiffness matrix from local stiffness matrices 
!------------------------------------------------------------------------------
!This is done only once because the stiffmatrix is constant 
           IF ( .NOT. StiffMatrixAssemblyDone ) THEN
             CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce(1:n) )
           ELSE
             CALL DefaultUpdateForce( LocalForce(1:n) )
           END IF          

        END DO
        
        IF (.NOT.StiffMatrixAssemblyDone) THEN
           CALL DefaultFinishBulkAssembly()
        END IF
       
        CALL Info( 'ComputeDevStressSSA', 'Assembly done', Level=4 )
        CALL DefaultFinishAssembly()

        StiffMatrixAssemblyDone = .TRUE.
!------------------------------------------------------------------------------

        CALL Info( 'ComputeDevStressSSA', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------

        UNorm = DefaultSolve()

        DO t=1,Solver % NumberOfActiveElements
           CurrentElement => GetActiveElement(t) 
           n = GetElementNOFNodes()
           DO i=1,n
              k = CurrentElement % NodeIndexes(i)
              Stress( StressDOFs*(StressPerm(k)-1) + COMP ) =    & 
                   Solver % Variable % Values( Solver % Variable % Perm(k) )
           END DO
        END DO

     END DO ! End DO Comp

!------------------------------------------------------------------------------
!  Fill the empty space of the SSABulkVelocity variable with pressure
!------------------------------------------------------------------------------
    
!!! Loop on the nodes
     DO t=1,Model % Mesh % NumberOfNodes
        rhoM = ListGetRealAtNode( Material, 'Density', t, GotIt)
        IF (.NOT. GotIt) rhoM = ListGetRealAtNode( Material, 'SSA Mean Density',t, GotIt )

        IF (dim==2) THEN
           altiM = Model % Nodes % y (t)
           altiTopM = Model % Nodes % y (TopPointer(t))
        ELSE
           altiM = Model % Nodes % z (t)
           altiTopM = Model % Nodes % z (TopPointer(t))
        END IF
        Pcryo = g*rho*(altiTopM-altiM)

        IF (.NOT. Cauchy) THEN
           IF (dim==2) THEN
              FlowValues((dim+1)*FlowPerm(t))=Pcryo-Stress(StressDOFs*(StressPerm(t)-1) + 1)
           ELSE
              FlowValues((dim+1)*FlowPerm(t))=Pcryo-Stress(StressDOFs*(StressPerm(t)-1) + 1) &         
                                                   -Stress(StressDOFs*(StressPerm(t)-1)+2)
           END IF
        ELSE
           IF (dim==2) THEN
              FlowValues((dim+1)*FlowPerm(t))=0.5*(Pcryo-Stress(StressDOFs*(StressPerm(t)-1) + 1))
           ELSE
              FlowValues((dim+1)*FlowPerm(t))=(1/3)*(Pcryo-Stress(StressDOFs*(StressPerm(t)-1) &
                                        + 1)-Stress(StressDOFs*(StressPerm(t)-1)+2))
           END IF
        END IF
     END DO
  
!~~~~~~~~~~~~~~~~~
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalNSMatrix(COMP, MassMatrix, StiffMatrix, ForceVector, &
       NodalVelo, NodalP, &
       Element, n, Nodes, Cauchy )
!------------------------------------------------------------------------------
    
    REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
    REAL(KIND=dp) ::  NodalVelo(:,:)
    REAL(KIND=dp), DIMENSION(:) :: ForceVector, NodalP
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    LOGICAL ::  Cauchy
    INTEGER :: n, COMP
!------------------------------------------------------------------------------
!
    REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
    REAL(KIND=dp) :: dBasisdx(2*n,3),detJ, pBasis(n)
    
    REAL(KIND=dp) :: Stress, Sxx, Syy, epsi

    REAL(KIND=dp) :: Pressure
    REAL(KIND=dp) :: dVelodx(3,3), SR(3,3)
    
    INTEGER :: i, j, k, p, q, t, dim, cc, NBasis,  LinearBasis
    
    REAL(KIND=dp) :: s, u, v, w, eta
    REAL(KIND=dp) :: SecInv, stmp

    REAL(KIND=dp) :: dDispldx(3,3), Viscosity
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER :: N_Integ, nd
    INTEGER, DIMENSION(6), PARAMETER :: indx = (/1, 2, 3, 1, 2, 3/), &
                                         indy = (/1, 2, 3, 2, 3, 1/)
    
    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
    
    LOGICAL :: stat
    
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    cc=2*dim
    
    ForceVector = 0.0_dp
    StiffMatrix = 0.0_dp
    MassMatrix  = 0.0_dp
    
    IntegStuff = GaussPoints( Element )
     
    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n
!
!   Now we start integrating
!
    DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo(Element,Nodes,u,v,w,detJ, &
            Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)
       
       s = detJ * S_Integ(t)
       
     
       dVelodx = 0.0_dp  !Initialisation :dvx/dz=dvy/dz=0 (SSA) dvz/dx et dvz/dy négligeables
       IF (dim ==2) THEN
          dVelodx(1,1)= SUM( NodalVelo(1,1:n)*dBasisdx(1:n,1) )
          dVelodx(2,2)= - dVelodx(1,1)
       ELSE
         DO j=1,2
            dVelodx(1,j) = SUM( NodalVelo(1,1:n)*dBasisdx(1:n,j) ) !dvx/dx puis dvx/dy
            dVelodx(2,j) = SUM( NodalVelo(2,1:n)*dBasisdx(1:n,j) ) ! dvy/dx puis dvy/dy
         END DO
         dVelodx(3,3) = - dVelodx(1,1) - dVelodx(2,2)  !!!!dvz/dz  Conservation Masse 
       END IF


!  Calcul du second invariant (ATTENTION coordonnées cartésiennes requises
       
       SecInv = 0.0D0  

       DO i=1,3
         DO j=1,3
          stmp = 0.5 * (dVelodx(i,j) + dVelodx(j,i))
          SecInv = SecInv + stmp*stmp
         END DO
       END DO

       SecInv = SQRT(2.0_dp * SecInv) !!SecInv est le gama de la formulation power law (gama^2=2DijDij etIe^2=0.5DijDij)
       
       Viscosity = EffectiveViscositySSA( 1.0_dp, SecInv, Element, Nodes, &
             n, n, u, v, w )

!
! Strain-Rate
!
      
       SR = 0.5 * ( dVelodx + TRANSPOSE(dVelodx) )
       
!
!    Compute deviatoric stresses or Cauchy stresses: 
!    ----------------------------
      
       Stress = 2.0 * Viscosity * SR(indx(COMP),indy(COMP))
       Sxx = 2.0 * Viscosity * SR(1,1)
       Syy = 2.0 * Viscosity * SR(2,2)
       
       IF ((Cauchy).AND.(COMP.LE.3)) THEN
          IF (dim==2) THEN
             Pressure = SUM( NodalP(1:n)*Basis(1:n) ) - Sxx
          ELSE
             Pressure = SUM( NodalP(1:n)*Basis(1:n) ) - Sxx -Syy
          END IF
          Stress = Stress - Pressure
       END IF
       
       IF (.NOT.StiffMatrixAssemblyDone) THEN
        DO p=1,n         
          DO q=1,n        
             StiffMatrix(p,q) =  &
                  StiffMatrix(p,q) + s*Basis(q)*Basis(p)
          END DO
        END DO
       END IF

       DO p=1,n
          ForceVector(p) =  &
               ForceVector(p) + s*Stress*Basis(p) 
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalNSMatrix
!------------------------------------------------------------------------------

     
!CONTAINS

!------------------------------------------------------------------------------
  FUNCTION EffectiveViscositySSA(Density,gama,Element,Nodes, &
       n,nd,u,v,w) RESULT(mu)
    !------------------------------------------

     USE ModelDescription

     REAL(KIND=dp)  :: Density,u,v,w,mu,gama
!     REAL(KIND=dp), OPTIONAL :: muder
     TYPE(Nodes_t)  :: Nodes
     INTEGER :: n,nd
     TYPE(Element_t),POINTER :: Element

     !------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3)
     REAL(KIND=dp) :: Ie,SqrtMetric,SqrtElementMetric,Velo(3)
     REAL(KIND=dp) :: Metric(3,3), dVelodx(3,3), CtrMetric(3,3), &
          Symb(3,3,3), dSymb(3,3,3,3)

     INTEGER :: i,j,k
     LOGICAL :: stat,GotIt

     CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityFlag, TemperatureName
     TYPE(ValueList_t), POINTER :: Material
     REAL(KIND=dp) :: x, y, z, c1n(n), c2n(n), c3n(n), c4n(n), &
          c1, c2, c3, c4, c5, c6, c7, Temp, NodalTemperature(n), eta0, NodalViscosity(n), Tlimit, TempCoeff, &
          h, A1, A2, Q1, Q2, R, NodalEhF(n), EhF, ArrheniusFactor

     ! Temperature is needed for thermal models
     TYPE(Variable_t), POINTER :: TempSol 
     REAL(KIND=dp), POINTER :: Temperature(:)
     INTEGER, POINTER :: TempPerm(:)

     INTEGER(KIND=AddrInt) :: Fnc

     TYPE(Variable_t), POINTER :: Var

     REAL(KIND=dp) :: dist,F2,F3
     REAL(KIND=dp) :: KE_K, KE_E, KE_Z, CT, TimeScale,Clip, Cmu, Vals(n)

     CHARACTER(LEN=MAX_NAME_LEN) :: str

     LOGICAL :: SetArrheniusFactor=.FALSE.
     
     mu = 0.0_dp

     k = ListGetInteger( CurrentModel % Bodies(Element % BodyId) % Values, 'Material', &
          minv=1, maxv=CurrentModel % NumberOFMaterials )

     Material => CurrentModel % Materials(k) % Values

     ViscosityFlag = ListGetString( Material,'Viscosity Model', GotIt)

     IF(.NOT. gotIt) RETURN
     !------------------------------------------------------------------------------
     !    Basis function values & derivatives at the calculation point
     !------------------------------------------------------------------------------
     stat = ElementInfo( Element,Nodes,u,v,w, &
          SqrtElementMetric, Basis,dBasisdx )

     SELECT CASE( ViscosityFlag )


     CASE('glen')
        c2n = ListGetReal( Material, 'Glen Exponent', n, Element % NodeIndexes, GotIt ) ! this is the real exponent, n, not 1/n
        IF (.NOT.GotIt) c2n(1:n) = 3.0_dp
        c2 = SUM( Basis(1:n) * c2n(1:n) )
        Ie = gama/2.0_dp ! Ie is the real second invariant used in the Glen's flow law formulation (gama^2=4*Ie^2)        
        
        SetArrheniusFactor = GetLogical(Material, 'Set Arrhenius Factor', GotIt)
        IF ( (.NOT.GotIt) .OR. .NOT.(SetArrheniusFactor)) THEN
           NodalTemperature(1:n) = ListGetReal(Material, 'Constant Temperature', n, Element % NodeIndexes, GotIt) !we are happy as is
           IF(.NOT.GotIt) THEN !we have to find a temperature field

              TemperatureName = GetString(Material, 'Temperature Field Variable', GotIt)
              IF (.NOT.GotIt) WRITE(TemperatureName,'(A)') 'Temperature'
              TempSol => VariableGet( CurrentModel % Variables,TRIM(TemperatureName))
              IF ( ASSOCIATED( TempSol) ) THEN
                 TempPerm    => TempSol % Perm
                 Temperature => TempSol % Values   
                 Temp =  SUM(Basis(1:n) * Temperature(TempPerm(Element % NodeIndexes(1:n))))
              ELSE
                 WRITE(Message, '(A,A,A)') 'Could not find variable ',&
                      TRIM(TemperatureName),' to inquire temperatur field for Glen'
                 CALL FATAL('ComputeDevStressSSA',Message)
              END IF
           
           ELSE
              Temp = SUM(Basis(1:n) * NodalTemperature(1:n))
           END IF
        
           R = GetConstReal( CurrentModel % Constants,'Gas Constant',GotIt)
           IF (.NOT.GotIt) R = 8.314_dp
           ! lets for the time being have this hardcoded
           Tlimit = GetConstReal(Material, 'Limit Temperature', GotIt)
           IF (.NOT.GotIt) THEN
              Tlimit = -10.0_dp
              CALL INFO('ComputeDevStressSSA','Limit Temperature not found. Setting to -10', Level=5)
           END IF
           A1 = GetConstReal(Material, 'Rate Factor 1', GotIt)
           IF (.NOT.GotIt) THEN
              A1 = 3.985d-13
              CALL INFO('ComputeDevStressSSA','Rate Factor 1 not found. Setting to 3.985e-13', Level=5)
           END IF
           A2 = GetConstReal(Material, 'Rate Factor 2', GotIt)
           IF (.NOT.GotIt) THEN
              A2 = 1.916d03
              CALL INFO('ComputeDevStressSSA','Rate Factor 2 not found. Setting to 1.916E03', Level=5)
           END IF
           Q1 = GetConstReal(Material, 'Activation Energy 1', GotIt)
           IF (.NOT.GotIt) THEN
              Q1 = 60.0d03
              CALL INFO('ComputeDevStessSSA','Activation Energy 1 not found. Setting to 60.0E03', Level=5)
           END IF
           Q2 = GetConstReal(Material, 'Activation Energy 2', GotIt)
           IF (.NOT.GotIt) THEN
              Q2 = 139.0d03
              CALL INFO('ComputeDevStressSSA','Activation Energy 2 not found. Setting to 139.0d03', Level=5)
           END IF
        
           IF (Temp.LE. Tlimit) THEN
              ArrheniusFactor = A1 * EXP( -Q1/(R * (273.15 + Temp)))
           ELSE IF((Tlimit<Temp) .AND. (Temp .LE. 0.0_dp)) THEN
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15 + Temp)))
           ELSE
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15)))
              CALL INFO('ComputeDevStressSSA','Positive Temperature detected in Glen - limiting to zero!', Level = 5)
           END IF
        ELSE
           ArrheniusFactor = GetConstReal(Material,'Arrhenius Factor', GotIt)
           IF (.NOT.(GotIt)) THEN 
              CALL FATAL('ComputeDevStressSSA','<Set Arrhenius Factor> is TRUE, but no value <Arrhenius Factor> found')
           END IF
        END IF

        NodalEhF(1:n) =  ListGetReal( Material, 'Glen Enhancement Factor', n, Element % NodeIndexes, GotIt )
        IF (.NOT.GotIt) NodalEhF(1:n) = 1.0_dp
        EhF = SUM(Basis(1:n) * NodalEhF(1:n))

        c3n = ListGetReal( Material, 'Critical Shear Rate',n, Element % NodeIndexes,GotIt )
        IF (GotIt) THEN
           c3 = SUM( Basis(1:n) * c3n(1:n) )
           IF(Ie < c3) THEN
              Ie = c3
           END IF
        END IF

        ! compute the effective viscosity
        mu = 0.5_dp * (EhF * ArrheniusFactor)**(-1.0_dp/c2) * Ie**((1.0_dp/c2)-1.0_dp);

     CASE('power law')
           
        NodalViscosity(1:n) = ListGetReal(Material, 'Viscosity', n, Element % NodeIndexes, GotIt) 
        IF(.NOT.GotIt) THEN 
           WRITE(Message,'(A)')'Variable Viscosity not found. Setting to 1.0'
           CALL INFO('ComputeDevStressSSA',Message, Level = 20)
           NodalViscosity(1:n) = 1.0_dp
        END IF
        eta0 = SUM ( NodalViscosity(1:n)*Basis(1:n) )
    
        c2n = ListGetReal( Material, 'Viscosity Exponent', n, Element % NodeIndexes )
        c2 = SUM( Basis(1:n) * c2n(1:n) )

        c3n = ListGetReal( Material, 'Critical Shear Rate',n, Element % NodeIndexes,gotIt )
        IF (GotIt) THEN
           c3 = SUM( Basis(1:n) * c3n(1:n) )
           IF((gama/2.0_dp) < c3) THEN
              gama = 2.0_dp * c3
           END IF
        END IF
        mu = eta0 * gama**(c2-1)

        c4n = ListGetReal( Material, 'Nominal Shear Rate',n, Element % NodeIndexes,gotIt )
        IF (GotIt) THEN
           c4 = SUM( Basis(1:n) * c4n(1:n) )
           mu = mu / c4**(c2-1)
        END IF

     CASE('power law too')
           
        NodalViscosity(1:n) = ListGetReal(Material, 'Viscosity', n, Element % NodeIndexes, GotIt) 
        IF(.NOT.GotIt) THEN 
           WRITE(Message,'(A)')'Variable Viscosity not found. Setting to 1.0'
           CALL INFO('ComputeDevStressSSA',Message, Level = 20)
           NodalViscosity(1:n) = 1.0_dp
        END IF
        eta0 = SUM ( NodalViscosity(1:n)*Basis(1:n) )
       
        c2n = ListGetReal( Material, 'Viscosity Exponent', n, Element % NodeIndexes )
        c2 = SUM( Basis(1:n) * c2n(1:n) )
        mu = eta0 **(-1/c2)* gama**(-(c2-1)/(2*c2)) 

     CASE DEFAULT
        CALL WARN('ComputeDevStressSSA','Unknown material model')

     END SELECT

  END FUNCTION EffectiveViscositySSA
!------------------------------------------------------------------------------
END SUBROUTINE ComputeDevStressSSA
!------------------------------------------------------------------------------
