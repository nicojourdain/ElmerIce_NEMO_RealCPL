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
!> SSolver to inquire the velocity from the SSA solution            
SUBROUTINE SSABasalSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the in-plane basal velocity with the SSA solution !
!  To be computed only at the base. Use then the SSASolver to export verticaly 
!  the basal velocity and compute the vertical velocity and pressure (if needed)
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
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: PointerToVariable, ZsSol, ZbSol, VeloSol, Nsol
  TYPE(Variable_t), POINTER :: GMSol,BedrockSol

  LOGICAL, SAVE :: AllocationsDone = .FALSE.
  LOGICAL :: Found, GotIt
  LOGICAL :: Newton

  INTEGER :: i, n, m, t, istat, DIM, p, STDOFs
  INTEGER :: NonlinearIter, NewtonIter,NewtonMinIter, iter
  INTEGER :: NewtonConvergence,NewtonConvergenceTol
  INTEGER :: GLnIP ! number of Integ. Points for GL Sub-element parametrization
  LOGICAL :: SEP ! Sub-element parametrization for Grounding line

          
  INTEGER, POINTER :: Permutation(:), ZsPerm(:), ZbPerm(:), NPerm(:)

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), Zs(:), Zb(:), Nval(:)
                            
  REAL(KIND=dp) :: UNorm, NonlinearTol, NewtonTol,PrevUNorm, relativeChange, minv

  REAL(KIND=dp),SAVE :: rhow, sealevel

  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: SolverName
  REAL(KIND=dp) :: st, st0
  REAL(KIND=dp),SAVE :: totat0, totst0
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at, at0
#else
  REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif
  LOGICAL,SAVE :: Firsttime=.True.


!! Global vars for SSAStar
  LOGICAL,SAVE :: SolveSSAStar=.False.,InternalIntegrate=.False.,GotMueff=.False.,Slip0=.False.
  LOGICAL,SAVE :: Timer=.False.
  INTEGER,SAVE :: nlevels
  REAL(kind=dp),dimension(:),allocatable,save :: zr
  REAL(kind=dp) :: mr,h1,hn
  LOGICAL :: GotRatio

  TYPE(Solver_t), POINTER :: PSolver
  INTEGER, POINTER,SAVE :: TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:)
  INTEGER :: nlayers

!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  STDOFs = PointerToVariable % DOFs 
  WRITE(SolverName, '(A)') 'SSASolver-SSABasalSolver'
  IF (.NOT.((STDOFS.EQ.1).OR.(STDOFS.EQ.2))) &
     CALL FATAL(SolverName,'Var DOFs has to be 1 or 2')

  IF (Firsttime) then
     totat0=0._dp
     totst0=0._dp
     Firsttime=.False.
  END IF
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
        DIM = CoordinateSystemDimension()

 ! IF DIM = STDOFs+1  Normal-Tangential can not be used => trick temporary set Model Dimension to STDOFs
       IF (DIM.eq.(STDOFs+1)) CurrentModel % Dimension = STDOFs
 
 ! Get Pointer to auxiliary variables requires to solve the SSA
        ZbSol => VariableGet( Solver % Mesh % Variables, 'Zb' )
        IF (ASSOCIATED(ZbSol)) THEN
           Zb => ZbSol % Values
           ZbPerm => ZbSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Zb<')
        END IF

        ZsSol => VariableGet( Solver % Mesh % Variables, 'Zs' )
        IF (ASSOCIATED(ZsSol)) THEN
           Zs => ZsSol % Values
           ZsPerm => ZsSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Zs<')
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

      Timer = GetLogical(Solver % Values, 'compute solver time',Found)

      Slip0 = GetLogical(Solver % Values,'Free slip by element',Found)

      SolveSSAStar=GetLogical(Solver % Values,'Solve SSA Star',Found)
      IF (SolveSSAStar) then
        CALL INFO(SolverName, 'Soving SSA Star equation', level=10)
        InternalIntegrate =GetLogical(Solver % Values,'SSA Star Internal viscosity integration',Found)
        if (InternalIntegrate) then
           CALL INFO(SolverName, 'SSA Star viscosity z-independant', level=10)
           nlevels=GetInteger(Solver % Values,'SSA Star integrations levels',Found)
           If (.NOT.Found) then
             nlevels=20
           endif
           write(message,'(A,A,I0,A)') trim(SolverName),' using ',nlevels,' integration levels'
           CALL INFO(SolverName, message,level=10)
           IF (Solver % Mesh % Changed) deallocate(zr)
           allocate(zr(nlevels))
           mr = ListGetConstReal(Solver % Values,'SSA Star integration levels ratio',GotRatio)
           if (GotRatio.AND.(mr.NE.1.0)) then
             h1 = (1-mr**(1.0_dp/(nlevels-1)))/(1-mr)
             zr(1) = 0.0_dp
             hn = h1;
             DO i=2,nlevels-1
               zr(i) = zr(i-1) + hn;
               hn = hn * ( mr**(1.0_dp/(nlevels-1)) )
             END DO
             zr(nlevels) = 1.0_dp
           else
             DO i=1,nlevels     
               zr(i) = (i-1)/(1._dp * (nlevels-1))
             END DO
           endif
           CALL Info(trim(SolverName),'SSA Star z-integration levels ready',Level=11)
           DO i=1,nlevels
               WRITE( Message, '(A,I0,A,ES12.4)') 'zr(',i,') : ',zr(i)
               CALL Info(trim(SolverName), Message, Level=11 )
           END DO
       Else
          GotMueff = GetLogical(Solver % Values,'SSA Star user defined viscosity',Found)
          If (.NOT.GotMueff) then ! viscosity has nodal values => detect Extruded Structure
            
           PSolver => Solver
           CALL DetectExtrudedStructure( Solver%Mesh, PSolver, &
                 TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
                 UpNodePointer = UpPointer, DownNodePointer = DownPointer, &
                 NumberOfLayers = nlayers )
           nlevels = nlayers + 1
           allocate(zr(nlevels))
          End if  
       Endif
    END IF

     ! Allocate

     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

     ! Get sea level and water density from the constant section
     sealevel = GetCReal( Model % Constants, 'Sea Level', Found )
     If (.NOT.Found) THEN
       WRITE(Message,'(A)') 'Constant >Sea Level< not found. &
         &Setting to 0.0'
       CALL INFO(SolverName, Message, level=3)
       sealevel = 0._dp
     End if
     rhow = GetCReal( Model % Constants, 'water density', Found )
     If (.NOT.Found) THEN
       WRITE(Message,'(A)') 'Constant Water Density not found. &
          &Setting to 1.03225e-18'
       CALL INFO(SolverName, Message, level=3)
       rhow = 1.03225d-18
     END IF
!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
     SEP=GetLogical( Solver % Values, 'Sub-Element GL parameterization',Found)
     IF (.NOT.Found) SEP=.False.
     IF (SEP) THEN
        GLnIP=GetInteger( Solver % Values, &
              'GL integration points number',GotIt )
        IF ( .NOT.GotIt ) &
          CALL FATAL(SolverName,'<GL integration points number> not found')

        GMSol => VariableGet( Solver % Mesh % Variables, 'GroundedMask' )
        IF (.NOT.ASSOCIATED(GMSol)) &
           CALL FATAL(SolverName,'Variable <GoundedMask> not found')
        BedrockSol => VariableGet( Solver % Mesh % Variables, 'bedrock' )
        IF (.NOT.ASSOCIATED(BedrockSol)) &
           CALL FATAL(SolverName,'Variable <bedrock> not found')
     END IF

!------------------------------------------------------------------------------
      NonlinearTol = GetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )

      NonlinearIter = GetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

      IF ( .NOT.GotIt ) NonlinearIter = 1

      NewtonTol = ListGetConstReal( Solver % Values, &
              'Nonlinear System Newton After Tolerance', minv=0.0d0 )

      NewtonIter = ListGetInteger( Solver % Values, &
              'Nonlinear System Newton After Iterations', GotIt )
      if (.NOT.Gotit) NewtonIter = NonlinearIter + 1

      NewtonMinIter = ListGetInteger( Solver % Values, &
              'Nonlinear System Newton minimum Iterations', GotIt )
      if (.NOT.Gotit) NewtonMinIter =  1

      NewtonConvergence=0
      NewtonConvergenceTol = ListGetInteger( Solver % Values, &
            'Nonlinear System Newton Max divergent iterations', GotIt )
      if (.NOT.Gotit)  NewtonConvergenceTol=NonlinearIter
    
      Newton=.False.

!------------------------------------------------------------------------------
      DO iter=1,NonlinearIter


       CALL Info( SolverName, ' ', Level=4 )
       CALL Info( SolverName, ' ', Level=4 )
       CALL Info( SolverName, &
                   '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'SSA BASAL VELOCITY NON-LINEAR ITERATION', iter
       CALL Info( SolverName, Message, Level=4 )
       If (Newton) Then
           WRITE( Message, * ) 'Newton linearisation is used'
           CALL Info( SolverName, Message, Level=4 )
       Endif
       CALL Info( SolverName, ' ', Level=4 )
       CALL Info( SolverName, &
                   '-------------------------------------',Level=4 )
       CALL Info( SolverName, ' ', Level=4 )

       IF (Timer) then
         at  = CPUTime()
         at0 = RealTime()
       END IF

       CALL DefaultInitialize()
       
       CALL BulkAssembly()
       CALL DefaultFinishBulkAssembly()
  
       CALL BoundaryAssembly()
       CALL DefaultFinishAssembly()
       ! Dirichlet 
       CALL DefaultDirichletBCs()
  
       IF (Timer) then
         at=CPUTime()-at
         at0 = RealTime() - at0
       END IF

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm

      IF (Timer) then
        st=CPUTime()
        st0 = RealTime()
      END IF

      UNorm = DefaultSolve()

      IF (Timer) then
        st=CPUTime()-st
        st0 = RealTime() - st0

        totat0 = totat0 + at0
        totst0 = totst0 + st0

        WRITE( Message, '(A,F8.2,F8.2,F8.2)') 'Assembly time (CPU,REAL,TOT REAL)',at,at0,totat0
        CALL Info(SolverName, Message, Level=3 )
        WRITE( Message, '(A,F8.2,F8.2,F8.2)' ) 'solve time (CPU,REAL,TOT REAL)',st,st0,totst0
        CALL Info(SolverName, Message, Level=3 )
      END IF

      if (Newton) then
          IF (Solver % Variable % NonlinChange.GT.RelativeChange) THEN
             NewtonConvergence=NewtonConvergence+1
          ELSE
             NewtonConvergence=0
          ENDIF
      endif

      RelativeChange = Solver % Variable % NonlinChange

      WRITE( Message, * ) 'Result Norm   : ', UNorm, PrevUNorm
      CALL Info(SolverName, Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ', RelativeChange
      CALL Info(SolverName, Message, Level=4 )

  
     IF (( RelativeChange < NewtonTol .OR. &
           iter > NewtonIter ).AND.(iter>NewtonMinIter)) Newton =.TRUE.

     IF (NewtonConvergence.GT.NewtonConvergenceTol)  Newton = .False.

!------------------------------------------------------------------------------
      IF ( RelativeChange < NonLinearTol ) EXIT
!------------------------------------------------------------------------------

  END DO ! Loop Non-Linear Iterations

  !!! reset Model Dimension to dim
  IF (DIM.eq.(STDOFs+1)) CurrentModel % Dimension = DIM

CONTAINS

!------------------------------------------------------------------------------
SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Element
      INTEGER :: t,n,nd

      !$omp parallel do private(Element,n,nd)
      DO t=1,GetNOFActive()
           Element => GetActiveElement(t)
           n  = GetElementNOFNodes(Element)
           nd = GetElementNOFDOFs(Element)
           CALL LocalMatrixUVSSA( Element, n, nd, STDOFs )
      END DO
      !$omp end parallel do
!------------------------------------------------------------------------------
END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE LocalMatrixUVSSA( Element, n, nd , STDOFs)
!------------------------------------------------------------------------------

TYPE(Element_t), POINTER :: Element
INTEGER :: n !ElementNOFNodes
INTEGER :: nd !ElementNOFDOFs
INTEGER :: STDOFs !Variable % DOFs

REAL(KIND=dp) :: STIFF(STDOFs*nd,STDOFs*nd),Jac(STDOFs*nd,STDOFs*nd)
REAL(KIND=dp) :: FORCE(STDOFs*nd),SOL(STDOFs*nd)
REAL(KIND=dp) :: gravity(n), Density(n),Viscosity(n),LocalZb(n),LocalZs(n)
REAL(KIND=dp) :: LocalU(n),LocalV(n),LocalBeta(n),LocalLinVelo(n)
REAL(KIND=dp) :: LocalC(n),LocalN(n), Localf(n)
REAL(KIND=dp) :: LocalGM(n), LocalBedrock(n) !Local Grounded Mask, bedrock
REAL(KIND=dp) :: cm, fm, PostPeak
INTEGER :: iFriction
CHARACTER(LEN=MAX_NAME_LEN) :: Friction

INTEGER, POINTER :: NodeIndexes(:)

LOGICAL :: Found
LOGICAL :: NewtonLin, fNewtonLin
LOGICAL :: PartlyGroundedElement
LOGICAL, SAVE :: AllocationDone=.False.
!$omp threadprivate(AllocationDone)

TYPE(ValueList_t), POINTER ::  BodyForce, Material

TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)

LOGICAL :: stat
TYPE(GaussIntegrationPoints_t) :: IP
REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ

REAL(KIND=dp) :: g, rho, eta,eta0, h, bottom, dhdx, dhdy , muder,DEe
REAL(KIND=dp) :: bedrock,Hf
REAL(KIND=dp) :: beta, LinVelo, fC, fN, fT, Velo(2), ub, alpha, fB
REAL(KIND=dp) :: gradS(2), A(2,2), StrainA(2,2), StrainB(2,2), Exx, Eyy, Exy, Ezz, Ee, MinSRInv ,MinH                           
REAL(KIND=dp) :: mueff1,mueff2,SIAInv,rgs2,dmueff1,dmueff2
REAL(KIND=dp) :: Slip, Slip2, MinN, SlipC, SlipW
REAL(KIND=dp) :: LevelEta(n),eta1,zz
INTEGER :: LevelNodeIndexes(n)

INTEGER :: i, j, t, p, q , dim

if (.NOT.AllocationDone) then
   allocate(Nodes%x(Model % MaxElementNodes),&
            Nodes%y(Model % MaxElementNodes),&
            Nodes%z(Model % MaxElementNodes))
   AllocationDone=.True.
end if

! set coords of highest occuring dimension to zero (to get correct path element)
NodeIndexes => Element%NodeIndexes

!PRINT *, 'In Bulk'

Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
IF (STDOFs == 1) THEN !1D SSA
   Nodes % y(1:n) = 0.0_dp
   Nodes % z(1:n) = 0.0_dp
ELSE !2D SSA
   Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
   Nodes % z(1:n) = 0.0_dp
END IF

! Get Body Forces and Material properties
Material => GetMaterial(Element)
BodyForce=> GetBodyForce(Element)


gravity = 0.0_dp
IF ( ASSOCIATED( BodyForce ) ) THEN
   IF (STDOFs==1) THEN
      gravity(1:n) = GetReal(BodyForce, 'Flow BodyForce 2', Found, Element)
   ELSE
      gravity(1:n) = GetReal(BodyForce, 'Flow BodyForce 3', Found, Element)
   END IF
END IF

Density(1:n) = GetReal( Material, 'SSA Mean Density', Found, Element)
IF (.NOT.Found) &
   CALL FATAL(SolverName,'Could not find Material prop.  >SSA Mean Density<')


cm = ListGetConstReal( Material, 'Viscosity Exponent',Found)
If (.NOT.Found) &
   CALL FATAL(SolverName,'Could not find Material prop.  >Viscosity Exponent<')
MinSRInv = ListGetConstReal( Material, 'Critical Shear Rate',Found)
If (.NOT.Found) MinSRInv = AEPS

MinH = ListGetConstReal( Material, 'SSA Critical Thickness',Found)
If (.NOT.Found) MinH=EPSILON(MinH)

IF (SolveSSAStar) then
   If (InternalIntegrate) then
      Viscosity(1:n) = GetReal( Material, 'Viscosity', Found, Element)
      IF (.NOT.Found) &
         CALL FATAL(SolverName,'Could not find Material prop. >Viscosity<')
   Else if (GotMueff) then
      Viscosity(1:n) = GetReal( Material, 'SSAStar Integrated Viscosity', Found, Element)
      IF (.NOT.Found) &
         CALL FATAL(SolverName,'Could not find Material prop. >SSAStar Integrated Viscosity<')
   End if
ELSE
   Viscosity(1:n) = GetReal( Material, 'SSA Mean Viscosity',Found, Element)
   IF (.NOT.Found) &
      CALL FATAL(SolverName,'Could not find Material prop. >SSA Mean Viscosity<')
END IF

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
   CASE('tsai')
     iFriction = 4
   CASE DEFAULT
     CALL FATAL(SolverName,'Friction should be linear, Weertman or Coulomb')
END SELECT

! for all friction law
LocalBeta = 0.0_dp
LocalBeta(1:n) = GetReal( Material, 'SSA Friction Parameter', Found, Element)
IF (.NOT.Found) &
   CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Parameter<')
IF (iFriction > 1) THEN
   fm = ListGetConstReal( Material, 'SSA Friction Exponent', Found )
   IF (.NOT.Found) &
      CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Exponent<')
   LocalLinVelo = 0.0_dp
   LocalLinVelo(1:n) = GetReal(Material, 'SSA Friction Linear Velocity', Found, Element)
   IF (.NOT.Found) &
      CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Linear Velocity<')
END IF

! for Coulomb and Tsai friction
IF (iFriction > 2) THEN
   MinN = ListGetConstReal( Material, 'SSA Min Effective Pressure', Found )
   IF (.NOT.Found) THEN
      MinN = 1.0e-6_dp
      WRITE( Message, * ) 'Parameter >SSA Min Effective Pressure< not found. Set to 10E-6'
      CALL Info(SolverName, Message, Level=20 )
   END IF
   ! Get the effective pressure
   IF (ASSOCIATED(NSol)) THEN
      LocalN(1:n) = Nval(NPerm(NodeIndexes(1:n)))
   ELSE
      CALL FATAL(SolverName,'Could not find variable >Effective Pressure<')
   END IF
END IF

! for Coulomb only 
IF (iFriction == 3) THEN  
   PostPeak = ListGetConstReal( Material, 'SSA Friction Post-Peak', Found )
   IF (.NOT.Found) &
      CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Post-Peak<')
   LocalC = 0.0_dp
   LocalC(1:n) = GetReal(Material, 'SSA Friction Maximum Value', Found, Element)
   IF (.NOT.Found) &
      CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Maximum Value<')
END IF

! for Tsai only
IF (iFriction == 4) THEN
   Localf = 0.0_dp
   Localf(1:n) = GetReal(Material, 'SSA Tsai Coulomb Friction Parameter', Found, Element)
   IF (.NOT.Found) &
      CALL FATAL(SolverName,'Could not find Material prop. >SSA Tsai Coulomb Friction Parameter<')
END IF

! Get the Nodal value of Zb and Zs
LocalZb(1:n) = Zb(ZbPerm(NodeIndexes(1:n)))
LocalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))


! Previous Velocity
LocalU(1:n) = VariableValues(STDOFs*(Permutation(NodeIndexes(1:n))-1)+1)
LocalV = 0.0_dp
IF (STDOFs.EQ.2) LocalV(1:n) = VariableValues(STDOFs*(Permutation(NodeIndexes(1:n))-1)+2)


! Use Newton Linearisation
NewtonLin = (Newton.AND.(cm.NE.1.0_dp))
fNewtonLin = (Newton.AND.(fm.NE.1.0_dp))


! Do Local Assembly
STIFF = 0.0_dp
FORCE = 0.0_dp
Jac=0.0_dp

IF (SEP) THEN
  ! Get Nodal Values of Grounded Mask
  LocalGM(1:n)=GMSol%Values(GMSol%Perm(NodeIndexes(1:n)))
  LocalBedrock(1:n)=BedrockSol%Values(BedrockSol%Perm(NodeIndexes(1:n)))
  PartlyGroundedElement=(ANY(LocalGM(1:n).GE.0._dp).AND.ANY(LocalGM(1:n).LT.0._dp))
  IF (PartlyGroundedElement) THEN
     IP = GaussPoints( Element , np=GLnIP )
  ELSE
     IP = GaussPoints( Element )
  ENDIF
ELSE
  IP = GaussPoints( Element )
ENDIF

DO t=1,IP % n
   stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
   IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

! Needed Intergration Point value

   g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
   rho = SUM( Density(1:n) * Basis(1:n) )
   eta0 = SUM( Viscosity(1:n) * Basis(1:n) )
   gradS = 0.0_dp
   gradS(1) = SUM( LocalZs(1:n) * dBasisdx(1:n,1) )
   IF (STDOFs == 2) gradS(2) = SUM( LocalZs(1:n) * dBasisdx(1:n,2) )
   h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n) )
   h=max(h,MinH)

   bottom=SUM( LocalZb(1:n) * Basis(1:n))

   beta = SUM( LocalBeta(1:n) * Basis(1:n) )

   if (Slip0) then
     if (any(LocalBeta(1:n).lt.1.0e-16)) then
        beta=0._dp
     endif
   End if

   IF (SEP) THEN
     IF (ALL(LocalGM(1:n).LT.0._dp)) THEN
        beta=0._dp
     ELSE IF (PartlyGroundedElement) THEN
        bedrock = SUM( LocalBedrock(1:n) * Basis(1:n) )
        Hf= rhow * (sealevel-bedrock) / rho
        if (h.lt.Hf) beta=0._dp
     END IF
   END IF

   IF (iFriction > 1) THEN
      LinVelo = SUM( LocalLinVelo(1:n) * Basis(1:n) )
      IF ((iFriction == 2).AND.(fm==1.0_dp)) iFriction=1
         Velo = 0.0_dp
         Velo(1) = SUM(LocalU(1:n) * Basis(1:n))
         IF (STDOFs == 2) Velo(2) = SUM(LocalV(1:n) * Basis(1:n))
         ub = SQRT(Velo(1)*Velo(1)+Velo(2)*Velo(2))
         Slip2=1.0_dp
         IF (ub < LinVelo) then 
            ub = LinVelo
            Slip2=0.0_dp
         ENDIF
   END IF

   IF (iFriction==3) THEN
      fC = SUM( LocalC(1:n) * Basis(1:n) )
      fN = MAX(SUM( LocalN(1:n) * Basis(1:n) ) , MinN)
   END IF

   IF (iFriction==4) THEN
      fT = SUM( Localf(1:n) * Basis(1:n) )
      fN = MAX(SUM( LocalN(1:n) * Basis(1:n) ) , MinN)
   END IF
       
   IF (iFriction==1) THEN
      Slip = beta
      fNewtonLin = .FALSE.
   ELSE IF (iFriction==2) THEN
      Slip = beta * ub**(fm-1.0_dp) 
      Slip2 = Slip2*Slip*(fm-1.0_dp)/(ub*ub)
   ELSE IF (iFriction==3) THEN
      IF (PostPeak.NE.1.0_dp) THEN
         alpha = (PostPeak-1.0_dp)**(PostPeak-1.0_dp) / PostPeak**PostPeak
      ELSE
         alpha = 1.0_dp
      END IF
      fB = alpha * (beta / (fC*fN))**(PostPeak/fm)
      Slip = beta * ub**(fm-1.0_dp) / (1.0_dp + fB * ub**PostPeak)**fm
      Slip2 = Slip2 * Slip * ((fm-1.0_dp) / (ub*ub) - &
           fm*PostPeak*fB*ub**(PostPeak-2.0_dp)/(1.0_dp+fB*ub**PostPeak))
   ELSE IF (iFriction==4) THEN
      SlipW = beta * ub**(fm-1.0_dp)
      SlipC = fT * fN * ub**(-1.0_dp)
      IF (SlipW .LT.  SlipC) THEN
         Slip = SlipW
         Slip2 = Slip2*Slip*(fm-1.0_dp)/(ub*ub)
      ELSE 
         Slip = SlipC
         Slip2 = - Slip2*Slip/(ub*ub)
      END IF
   END IF
!------------------------------------------------------------------------------
! In the non-linear case, effective viscosity       
   eta = 0._dp
   IF (cm.NE.1.0_dp) THEN
      Exx = SUM(LocalU(1:n)*dBasisdx(1:n,1))
      Eyy = 0.0_dp
      Exy = 0.0_dp
      IF (STDOFs.EQ.2) THEN
         Eyy = SUM(LocalV(1:n)*dBasisdx(1:n,2))
         Ezz = -Exx - Eyy
         Exy = SUM(LocalU(1:n)*dBasisdx(1:n,2))
         Exy = 0.5_dp*(Exy + SUM(LocalV(1:n)*dBasisdx(1:n,1)))
         Ee = 0.5_dp*(Exx*Exx + Eyy*Eyy + Ezz*Ezz) + Exy*Exy
      ELSE
         Ee = Exx*Exx
      END IF
      muder = eta0 * 0.5_dp * (2.0_dp**cm) * ((cm-1.0_dp)/2.0_dp) *  Ee**((cm-1.0_dp)/2.0_dp - 1.0_dp)
      dEe = 1.0_dp
      IF (sqrt(Ee) < MinSRInv) THEN
         Ee = MinSRInv*MinSRInv
         dEe = 0.0_dp
      END IF
         IF (SolveSSAStar) then
            If (GotMueFF) then
               ! effective viscosity prescribed by the user
               eta = eta0
               muder=0._dp
            ! need to integrate effective viscosity
            Else 
              ! eta = int_zb^zs mueff dz
              eta=0._dp
              ! muder = deta/dEe
              muder=0._dp
            
              ! constant terms for the SIA stress invaraint
              rgs2=(gradS(1)*gradS(1)+gradS(2)*gradS(2))*(rho*g)**2

              ! passe en coordonnées xsi pour l'integration verticale
              ! pour l'instant suppose eta0=f(x,y)
              If (InternalIntegrate) then ! constant viscosity with z
                  eta1=eta0
              else !need to get eta at ecah level
                 LevelNodeIndexes(1:n)=Element % NodeIndexes(1:n)
                 LevelEta(1:n) = ListGetReal( Material, 'Viscosity',n, LevelNodeIndexes,Found)
                 IF (.NOT.Found) &
                      CALL FATAL(SolverName,'Could not find Material prop. >Viscosity<')
                 eta1=SUM( LevelEta(1:n) * Basis(1:n) )
                 zr(1)=0._dp
              endif

              !mueff=sol(mueff*eta0**-n (4 * mueff**2 * Ee + phi2)**(n-1)/2 - 1 = 0) (1)
              !phi2=rgs2*depth^2
              ! Get 1rst value
              SIAInv=rgs2*h*h*(1.0-zr(1))*(1.0-zr(1))
              mueff1=MueffRoot(eta1,cm,Ee,SIAInv,dmueff1)

              Do i=2,nlevels

                If (.NOT.InternalIntegrate) then !Nedd to get viscosity from
                                                 !nodes in the layer above
                  LevelNodeIndexes(1:n)=UpPointer(LevelNodeIndexes(1:n))
                  LevelEta(1:n) = ListGetReal( Material, 'Viscosity',n, LevelNodeIndexes,Found)
                  IF (.NOT.Found) &
                       CALL FATAL(SolverName,'Could not find Material prop. >Viscosity<')
                  eta1=SUM( LevelEta(1:n) * Basis(1:n) )

                  If (STDOFs.eq.1) then
                    zz=SUM(Model%Nodes%y(LevelNodeIndexes(1:n))* Basis(1:n) )
                  else
                    zz=SUM(Model%Nodes%z(LevelNodeIndexes(1:n))* Basis(1:n) )
                  End if
                   zr(i)=(zz-bottom)/h
               End if

                SIAInv=rgs2*h*h*(1.0-zr(i))*(1.0-zr(i))
                mueff2=MueffRoot(eta1,cm,Ee,SIAInv,dmueff2)
            
                eta=eta+0.5*(mueff1+mueff2)*h*(zr(i)-zr(i-1))

                muder=muder+0.5*(zr(i)-zr(i-1))*(dmueff1+dmueff2)
                mueff1=mueff2
                dmueff1=dmueff2
            Enddo
           Endif
         ELSE
            ! normal SSA viso-cosity H * mu_bar
            eta = h * eta0 * 0.5_dp * (2.0_dp**cm) * Ee**((cm-1.0_dp)/2.0_dp)
         ENDIF
   ELSE !linear flow law
            eta = h * eta0
   END IF 

   StrainA=0.0_dp
   StrainB=0.0_dp
   IF (NewtonLin) THEN
      StrainA(1,1)=SUM(2.0_dp*dBasisdx(1:n,1)*LocalU(1:n))

      IF (STDOFs.EQ.2) THEN
         StrainB(1,1)=SUM(0.5_dp*dBasisdx(1:n,2)*LocalU(1:n))

         StrainA(1,2)=SUM(dBasisdx(1:n,2)*LocalV(1:n))
         StrainB(1,2)=SUM(0.5_dp*dBasisdx(1:n,1)*LocalV(1:n))

         StrainA(2,1)=SUM(dBasisdx(1:n,1)*LocalU(1:n))
         StrainB(2,1)=SUM(0.5_dp*dBasisdx(1:n,2)*LocalU(1:n))

         StrainA(2,2)=SUM(2.0_dp*dBasisdx(1:n,2)*LocalV(1:n))
         StrainB(2,2)=SUM(0.5_dp*dBasisdx(1:n,1)*LocalV(1:n))
      END IF
   END IF

   A = 0.0_dp
   DO p=1,n
      DO q=1,n
         A(1,1) = 2.0_dp*dBasisdx(q,1)*dBasisdx(p,1)  
         IF (STDOFs.EQ.2) THEN
            A(1,1) = A(1,1) + 0.5_dp*dBasisdx(q,2)*dBasisdx(p,2)
            A(1,2) = dBasisdx(q,2)*dBasisdx(p,1) + &
                             0.5_dp*dBasisdx(q,1)*dBasisdx(p,2)
            A(2,1) = dBasisdx(q,1)*dBasisdx(p,2) + &
                             0.5_dp*dBasisdx(q,2)*dBasisdx(p,1)
            A(2,2) = 2.0*dBasisdx(q,2)*dBasisdx(p,2) +&
                             0.5_dp*dBasisdx(q,1)*dBasisdx(p,1)  
         END IF
         A = 2.0_dp * eta * A
         DO i=1,STDOFs
            STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) +&
                  Slip * Basis(q) * Basis(p) * IP % S(t) * detJ
            DO j=1,STDOFs
               STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +& 
                      A(i,j) * IP % S(t) * detJ 
            END DO 
         END DO

         IF ((fNewtonLin).AND.(iFriction > 1)) THEN
            DO i=1,STDOFs
              Do j=1,STDOFs
               STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +&
                  Slip2 * Velo(i) * Velo(j) * Basis(q) * Basis(p) * IP % S(t) * detJ
              End do
            END DO
         END IF

         IF (NewtonLin) then
            muder=muder*dEe
            IF (STDOFs.EQ.1) THEN
               Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                   IP % S(t) * detJ * 2.0_dp * h * StrainA(1,1)*dBasisdx(p,1) * &
                   muder * 2.0_dp * Exx*dBasisdx(q,1) 
            ELSE IF (STDOFs.EQ.2) THEN
               Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                 IP % S(t) * detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ &
                 (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*Exx+Eyy)*dBasisdx(q,1)+Exy*dBasisdx(q,2)) 

               Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) +&
                 IP % S(t) * detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ & 
                 (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*Eyy+Exx)*dBasisdx(q,2)+Exy*dBasisdx(q,1)) 

               Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) +&
                 IP % S(t) * detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ & 
                 (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*Exx+Eyy)*dBasisdx(q,1)+Exy*dBasisdx(q,2)) 

               Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) +&
                 IP % S(t) * detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ &
                 (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*Eyy+Exx)*dBasisdx(q,2)+Exy*dBasisdx(q,1)) 
            END IF
         END IF

      END DO

      DO i=1,STDOFs
         FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) - &   
               rho*g*h*gradS(i) * IP % s(t) * detJ * Basis(p) 
      END DO
         
      IF ((fNewtonLin).AND.(iFriction>1)) THEN
         DO i=1,STDOFs
            FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) + &   
               Slip2 * Velo(i) * ub * ub * IP % s(t) * detJ * Basis(p) 
         END DO
      END IF
          
   END DO
END DO

IF (NewtonLin) THEN
   SOL(1:STDOFs*n:STDOFs)=LocalU(1:n)
   IF (STDOFs.EQ.2) SOL(2:STDOFs*n:STDOFs)=LocalV(1:n)

   STIFF(1:STDOFs*n,1:STDOFs*n) = STIFF(1:STDOFs*n,1:STDOFs*n) + &
       Jac(1:STDOFs*n,1:STDOFs*n)
   FORCE(1:STDOFs*n) = FORCE(1:STDOFs*n) + &
       MATMUL(Jac(1:STDOFs*n,1:STDOFs*n),SOL(1:STDOFs*n))
END IF


CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element,USolver=Solver)

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUVSSA
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION MueffRoot(eta0,cm,Ee,phi2,dmu) RESULT(mu)
 ! Solve for the effective viscosity mueff; according to Eq. 13 
 ! (Adaptive mesh, finite volume modeling of marine ice sheets, Cornford et al.,2013)
 !mueff=sol(mueff*eta0**-n (4 * mueff**2 * Ee + phi2)**(n-1)/2 - 1 = 0) (1)
 !
 ! use a Newton scheme from the initial guess given by SSA only effective viscosity
 !  -Is there only 1 solution whatever n??
 !  - Maybe could be usefull to implement analytical solution for n=3 to use
 !  newton linearistaion of the SSA system ??
  implicit none
  REAL(KIND=dp) :: eta0,cm,Ee,phi2
  REAL(KIND=dp),OPTIONAL :: dmu
  REAL(KIND=dp) :: mu,F,dF,dF_a,dF_b
  INTEGER :: compt

  REAL(KIND=dp),parameter :: eps=1.0e-12
  INTEGER,parameter :: MaxIter=150
  
  !initial guess phi2=0
  mu = eta0 * 2**(cm-1) * Ee**((cm-1.0)/2.0)
  F = mu * eta0**(-1.0/cm) * (4*mu*mu*Ee+phi2)**((1.0/cm -1)/2.0) - 1.0_dp
  dF_a=eta0**(-1.0/cm) * (4*mu*mu*Ee+phi2)**((1.0/cm -1)/2.0)
  dF_b=mu * eta0**(-1.0/cm) * ((1.0/cm -1)/2.0) * (4*mu*mu*Ee+phi2)**((1.0/cm-1)/2.0 - 1.0) 
  dF = dF_a + 8.0 * Ee * mu * dF_b

  compt=0
  Do while ((abs(F).gt.eps).AND.(compt.LT.MaxIter))
     mu=mu-F/dF
     F = mu * eta0**(-1.0/cm) * (4*mu*mu*Ee+phi2)**((1.0/cm -1)/2.0) - 1.0_dp
     dF_a=eta0**(-1.0/cm) * (4*mu*mu*Ee+phi2)**((1.0/cm -1)/2.0)
     dF_b=mu * eta0**(-1.0/cm) * ((1.0/cm -1)/2.0) * (4*mu*mu*Ee+phi2)**((1.0/cm-1)/2.0 - 1.0) 
     dF = dF_a + 8.0 * Ee * mu * dF_b
     compt=compt+1
  end do
     IF (compt.GE.MaxIter) CALL FATAL(SolverName,'MueffRoot: unable to find effective viscosity')
     IF (mu.LT.0._dp) CALL FATAL(SolverName,'MueffRoot: get negative effective viscosity')

     !dmueff1=d(mueff1)/d(Ee) from (1)
     IF (PRESENT(dmu)) dmu=-4.0*dF_b*mu*mu/dF

  END FUNCTION MueffRoot


!  
!------------------------------------------------------------------------------
  SUBROUTINE BoundaryAssembly()
!------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    INTEGER :: t,n,nd

!$omp parallel do private(Element,n,nd)
    DO t=1,GetNOFBoundaryElements()
      Element => GetBoundaryElement(t)
      IF (STDOFS.NE.1) Then
          IF(.NOT.ActiveBoundaryElement(Element)) CYCLE
      END IF
      IF ( GetElementFamily() == 1 ) CYCLE
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE

      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      CALL LocalMatrixBCSSA( Element, n, nd , STDOFs )
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryAssembly
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCSSA( BoundaryElement, n, nd, STDOFs)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: BoundaryElement,ParentElement
    INTEGER :: n, nd , STDOFs

    TYPE(ValueList_t), POINTER ::  BodyForce, Material, BC

    REAL(KIND=dp) :: STIFF(STDOFs*nd,STDOFs*nd), FORCE(STDOFs*nd)

    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL,SAVE :: AllocationDone=.False.
    !$omp threadprivate(AllocationDone)

    TYPE(Nodes_t), SAVE :: Nodes
    !$omp threadprivate(Nodes)

!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3), rhoi, g, alpha, h, h_im,norm

    REAL(KIND=dp) :: LocalZb(n),LocalZs(n),gravity(n),density(n)
    REAL(KIND=dp) :: MinH
    REAL(KIND=dp) :: Butt
    LOGICAL :: Stat
    INTEGER :: t, i, p
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: CalvingFront,Found
    INTEGER :: other_body_id 

!------------------------------------------------------------------------------

     BC => GetBC(BoundaryElement)
     IF (.NOT.ASSOCIATED( BC ) ) Return

     NodeIndexes => BoundaryElement % NodeIndexes

     if (.NOT.AllocationDone) then
        allocate(Nodes%x(Model % MaxElementNodes),&
                 Nodes%y(Model % MaxElementNodes),&
                 Nodes%z(Model % MaxElementNodes))
     AllocationDone=.True.
     end if

     ! set coords of highest occuring dimension to zero (to get correct path element)
     NodeIndexes => BoundaryElement%NodeIndexes

     Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
     IF (STDOFs == 1) THEN !1D SSA
        Nodes % y(1:n) = 0.0_dp
        Nodes % z(1:n) = 0.0_dp
     ELSE !2D SSA
        Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        Nodes % z(1:n) = 0.0_dp
     END IF


! Find the nodes for which 'Calving Front' = True             
     CalvingFront=.False. 
     CalvingFront = ListGetLogical( BC, 'Calving Front', Found )
     
     IF (.NOT.CalvingFront) RETURN

! Is there any buttressing prescribed as a force at the front 
     Butt = GetConstReal ( BC, 'Buttressing Factor', Found )
     IF (.NOT.Found) THEN
       WRITE(Message,'(A)')'Buttressing Factor Not Found. Default value is 1, i.e. no buttressing at the calving front'
       CALL INFO(SolverName, Message, Level =20)
       Butt=1.0_dp
     END IF

     LocalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))
     LocalZb(1:n) = Zb(ZbPerm(NodeIndexes(1:n)))
     
     ! Need to access Parent Element to get Material properties
     other_body_id = BoundaryElement % BoundaryInfo % outbody
     IF (other_body_id < 1) THEN ! only one body in calculation
          ParentElement => BoundaryElement % BoundaryInfo % Right
          IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
     ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
          ParentElement => BoundaryElement %  BoundaryInfo % Right
          IF (ParentElement % BodyId == other_body_id) ParentElement =>  BoundaryElement % BoundaryInfo % Left
     END IF

     ! Read Density in the Material Section
     Material => GetMaterial(ParentElement)
     BodyForce => GetBodyForce(ParentElement)

     Density=0.0_dp
     Density(1:n) = ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,Found)
     IF (.NOT.Found) &
           CALL FATAL(SolverName,'Could not find Material prop.  >SSA Mean Density<')

     MinH = ListGetConstReal( Material, 'SSA Critical Thickness',Found)
     If (.NOT.Found) MinH=EPSILON(MinH)

    ! Read the gravity in the Body Force Section 
     Gravity = 0.0_dp
     IF ( ASSOCIATED( BodyForce ) ) THEN
        IF (STDOFs==1) THEN 
           Gravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 2', n, NodeIndexes, Found)
        ELSE 
           Gravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 3', n, NodeIndexes, Found)
        END IF
     END IF

    STIFF = 0.0_dp
    FORCE = 0.0_dp

! The front force is a concentrated nodal force in 1D-SSA and
! a force distributed along a line in 2D-SSA    

! 1D-SSA Case : concentrated force at each nodes
    IF (STDOFs==1) THEN  !1D SSA but should be 2D problem (does elmer work in 1D?)
      DO i = 1, n
         g = ABS( Gravity(i) )
         rhoi = Density(i)
         h = LocalZs(i)-LocalZb(i) 
         h = max(h,MinH)
         h_im = max(0.0_dp,sealevel-LocalZb(i))
         alpha= Butt*(0.5_dp * g * (rhoi * h*h - rhow * h_im*h_im))
         FORCE(i) = FORCE(i) + alpha
      END DO

! 2D-SSA Case : force distributed along the line       
! This will work in DIM=3D only if working with Extruded Mesh and Preserve
! Baseline as been set to True to keep the 1D-BC 
    ELSE IF (STDOFs==2) THEN

          IP = GaussPoints( BoundaryElement )
          DO t=1,IP % n
             stat = ElementInfo( BoundaryElement, Nodes, IP % U(t), IP % V(t), &
                 IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
 
             g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
             rhoi = SUM( Density(1:n) * Basis(1:n) )
             h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n))
             h_im = max(0.0_dp , SUM( (sealevel-LocalZb(1:n)) * Basis(1:n)) )
             alpha= Butt*(0.5_dp * g * (rhoi * h*h - rhow * h_im*h_im))

! Normal in the (x,y) plane
             Normal = NormalVector( BoundaryElement, Nodes, IP % U(t), IP % V(t), .TRUE.)
             norm=SQRT(normal(1)*normal(1) +normal(2)*normal(2))
             Normal(1) = Normal(1)/norm
             Normal(2) = Normal(2)/norm

             DO p=1,n
                DO i=1,STDOFs
                   FORCE(STDOFs*(p-1)+i) =   FORCE(STDOFs*(p-1)+i) +&   
                    alpha * Normal(i) * IP % s(t) * detJ * Basis(p) 
                END DO
             END DO
          END DO

    ELSE   

      CALL FATAL('SSASolver-SSABasalSolver','Do not work for STDOFs <> 1 or 2')

    END IF

    CALL DefaultUpdateEquations( STIFF, FORCE,UElement=BoundaryElement )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCSSA
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE SSABasalSolver
!------------------------------------------------------------------------------










! *****************************************************************************
!>   Compute the depth integrated viscosity = sum_zb^zs eta dz
!>     and the depth integrated density = sum_zb^zs rho dz
SUBROUTINE GetMeanValueSolver( Model,Solver,dt,TransientSimulation )
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
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement, &
                             BoundaryElement
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable, IntViscoSol, IntDensSol,&
                                DepthSol  

  LOGICAL :: AllocationsDone = .FALSE., Found

  INTEGER :: i, n, m, t, istat, DIM, COMP, other_body_id   
  INTEGER, POINTER :: Permutation(:), NodeIndexes(:), IntViscoPerm(:),&
                      IntDensPerm(:), DepthPerm(:) 
       
  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), IntVisco(:), IntDens(:), Depth(:)
  REAL(KIND=dp) :: Norm, cn, dd 

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
           NodalVar(:) 

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName
  SAVE NodalVar 
!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  WRITE(SolverName, '(A)') 'SSASolver-IntValue'

  IntViscoSol => VariableGet( Solver % Mesh % Variables, 'Mean Viscosity' )
  IF (ASSOCIATED(IntViscoSol)) THEN
     IntVisco => IntViscoSol % Values
     IntViscoPerm => IntViscoSol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find variable >Mean Viscosity<')
  END IF
  IntDensSol => VariableGet( Solver % Mesh % Variables, 'Mean Density' )
  IF (ASSOCIATED(IntDensSol)) THEN
     IntDens => IntDensSol % Values
     IntDensPerm => IntDensSol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find variable >Mean Density<')
  END IF
  DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
  IF (ASSOCIATED(DepthSol)) THEN
     Depth => DepthSol % Values
     DepthPerm => DepthSol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find variable >Depth<')
  END IF
  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalVar) 

     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), NodalVar(N), &
                          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

! Loop for viscosity and density
DO COMP=1, 2
! No non-linear iteration, no time dependency  
  VariableValues = 0.0d0
  Norm = Solver % Variable % Norm

  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()
  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()

     NodeIndexes => Element % NodeIndexes
     Material => GetMaterial(Element)

     IF (COMP==1) THEN
     ! Read the Viscosity eta, 
     ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
     NodalVar = 0.0D0
     NodalVar(1:n) = ListGetReal( &
         Material, 'Viscosity', n, NodeIndexes, Found )
     ELSE IF (COMP==2) THEN
     NodalVar = 0.0D0
     NodalVar(1:n) = ListGetReal( &
         Material, 'Density', n, NodeIndexes, Found )
     END IF

     CALL LocalMatrix (  STIFF, FORCE, Element, n, NodalVar )
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  
  ! Neumann conditions 
  DO t=1,Solver % Mesh % NUmberOfBoundaryElements
     BoundaryElement => GetBoundaryElement(t)
     IF ( GetElementFamily() == 1 ) CYCLE
     NodeIndexes => BoundaryElement % NodeIndexes
     IF (ParEnv % myPe .NE. BoundaryElement % partIndex) CYCLE
     n = GetElementNOFNodes()

! Find the Parent element     
     other_body_id = BoundaryElement % BoundaryInfo % outbody
     IF (other_body_id < 1) THEN ! only one body in calculation
         ParentElement => BoundaryElement % BoundaryInfo % Right
         IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
         ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
             ParentElement => BoundaryElement % BoundaryInfo % Right
             IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
         END IF

     Material => GetMaterial(ParentElement)

     IF (COMP==1) THEN
     ! Read the Viscosity eta, 
     ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
     NodalVar = 0.0D0
     NodalVar(1:n) = ListGetReal( &
         Material, 'Viscosity', n, NodeIndexes, Found )
     ELSE IF (COMP==2) THEN
     NodalVar = 0.0D0
     NodalVar(1:n) = ListGetReal( &
         Material, 'Density', n, NodeIndexes, Found )
     END IF
     CALL LocalMatrixBC(  STIFF, FORCE, BoundaryElement, n, NodalVar)
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO

  CALL DefaultFinishAssembly()
  ! Dirichlet 
  IF (COMP==1) THEN
     CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
          'Mean Viscosity', 1,1, Permutation )
  ELSE
     CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
          'Mean Density', 1,1, Permutation )
  END IF
  Norm = DefaultSolve()

  ! Save the solution on the right variable
  IF (COMP==1) THEN
     DO i = 1, Model % Mesh % NumberOfNodes
        IF (IntViscoPerm(i)>0) THEN
            IntVisco(IntViscoPerm(i)) = VariableValues(Permutation(i)) 
            IF (Depth(DepthPerm(i))>0.0_dp) IntVisco(IntViscoPerm(i)) = &
                 IntVisco(IntViscoPerm(i)) / Depth(DepthPerm(i))
        END IF
     END DO
  ELSE IF (COMP==2) THEN
     DO i = 1, Model % Mesh % NumberOfNodes
        IF (IntDensPerm(i)>0) THEN
            IntDens(IntDensPerm(i)) = VariableValues(Permutation(i)) 
            IF (Depth(DepthPerm(i))>0.0_dp) IntDens(IntDensPerm(i)) = &
                IntDens(IntDensPerm(i)) / Depth(DepthPerm(i))
                                                 
                                                 
        END IF
     END DO
  END IF
  
END DO !COMP


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n, var)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), var(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ, grad
    LOGICAL :: Stat
    INTEGER :: t, p,q ,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
         
        grad  = SUM( var(1:n) * dBasisdx(1:n,dim) )
        FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ  * Basis(1:n)
       
       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(  STIFF, FORCE, Element, n, var ) 
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), var(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3), eta, grad 
    LOGICAL :: Stat
    INTEGER :: t, dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       grad  = SUM( var(1:n) * Basis(1:n) )

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) - grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------
END SUBROUTINE GetMeanValueSolver
!------------------------------------------------------------------------------


! *****************************************************************************
SUBROUTINE SSASolver( Model,Solver,dt,TransientSimulation )
!DEC$ATTRIBUTES DLLEXPORT :: SSASolver
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Export vertically the SSABasal Velocity (given as a Dirichlet Boundary condition) 
!  Compute also the vertical velocity and the pressure
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
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable, Grad1Sol, Grad2Sol, &
                               DepthSol, VeloSol

  LOGICAL :: AllocationsDone = .FALSE., Found

  INTEGER :: i, n, m, t, istat, DIM, p, Indexes(128), COMP 
  INTEGER, POINTER :: Permutation(:), VeloPerm(:), &
       DepthPerm(:), GradSurface1Perm(:), GradSurface2Perm(:), &
       NodeIndexes(:)

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), Depth(:), GradSurface1(:), &
                            GradSurface2(:), Velocity(:), PrevVelo(:,:)
  REAL(KIND=dp) :: Norm, cn, dd 

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
           NodalGravity(:), NodalDensity(:), &
           NodalDepth(:), NodalSurfGrad1(:), NodalSurfGrad2(:), &
           NodalU(:), NodalV(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
       

  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName
  SAVE NodalGravity, NodalDensity, &
           NodalDepth, NodalSurfGrad1, NodalSurfGrad2, &
           NodalU, NodalV
!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  WRITE(SolverName, '(A)') 'SSASolver'

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
        DIM = CoordinateSystemDimension()

        VeloSol => VariableGet( Solver % Mesh % Variables, 'SSAFlow' )
        IF (ASSOCIATED(veloSol)) THEN
           Velocity => VeloSol % Values
           VeloPerm => VeloSol % Perm
           PrevVelo => veloSol % PrevValues
        ELSE
           CALL FATAL(SolverName,'Could not find variable >SSAFlow<')
        END IF
        DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
        IF (ASSOCIATED(DepthSol)) THEN
           Depth => DepthSol % Values
           DepthPerm => DepthSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Depth<')
        END IF
        Grad1Sol => VariableGet( Solver % Mesh % Variables, 'FreeSurfGrad1')
        IF (ASSOCIATED(Grad1Sol)) THEN
           GradSurface1 => Grad1Sol % Values
           GradSurface1Perm => Grad1Sol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >FreeSurfGrad1<')
        END IF
        IF (dim > 2) THEN
           Grad2Sol => VariableGet( Solver % Mesh % Variables, 'FreeSurfGrad2')
           IF (ASSOCIATED(Grad2Sol)) THEN
              GradSurface2 => Grad2Sol % Values
              GradSurface2Perm => Grad2Sol % Perm
           ELSE
              CALL FATAL(SolverName,'Could not find variable >FreeSurfGrad2<')
           END IF
        END IF

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalGravity, &
                       NodalDensity, NodalDepth, &
                       NodalSurfGrad1, NodalSurfGrad2, NodalU, NodalV )

     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), &
          NodalGravity(N), NodalDensity(N), &
          NodalDepth(N), NodalSurfGrad1(N), NodalSurfGrad2(N), &
          NodalU(N), NodalV(N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF


     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

  ! Loop over the velocity components and pressure 
  ! If DIM = 2 u, w, p
  ! If DIM = 3 u, v, w, p
  !-----------------------------------------------
  DO  COMP = 1, DIM+1

! No non-linear iteration, no time dependency  
  VariableValues = 0.0d0
  Norm = Solver % Variable % Norm


  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()
  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()

     NodeIndexes => Element % NodeIndexes

     ! Read the gravity in the Body Force Section 
     BodyForce => GetBodyForce()
     NodalGravity = 0.0_dp
     IF ( ASSOCIATED( BodyForce ) ) THEN
           IF (DIM==2) THEN 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 2', n, NodeIndexes, Found)
           ELSE 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 3', n, NodeIndexes, Found)
           END IF
     END IF
     
     ! Read the Viscosity eta, density, and exponent m in Material Section
     ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
     Material => GetMaterial()

     NodalDensity = 0.0D0
     NodalDensity(1:n) = ListGetReal( &
         Material, 'Density', n, NodeIndexes, Found )

     ! Get the Nodal value of Depth, FreeSurfGrad1 and FreeSurfGrad2
     NodalDepth(1:n) = Depth(DepthPerm(NodeIndexes(1:n)))
     NodalSurfGrad1(1:n) = GradSurface1(GradSurface1Perm(NodeIndexes(1:n)))
     NodalSurfGrad2 = 0.0D0
     IF (DIM==3) NodalSurfGrad2(1:n) = GradSurface2(GradSurface2Perm(NodeIndexes(1:n)))

     IF (COMP==1) THEN     ! u
        CALL LocalMatrixUV (  STIFF, FORCE, Element, n ) 

     ELSE IF (COMP==DIM) THEN  ! w
        NodalU(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+1)
        NodalV = 0.0D0
        IF (DIM==3) NodalV(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+2)
        CALL LocalMatrixW (  STIFF, FORCE, Element, n, NodalU, NodalV ) 

     ELSE IF (COMP==DIM+1) THEN ! p
        CALL LocalMatrixP (  STIFF, FORCE, Element, n )

     ELSE               ! v if dim=3
        CALL LocalMatrixUV (  STIFF, FORCE, Element, n )

     END IF

     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  
  ! Neumann conditions only for w and p
  IF (COMP .GE. DIM) THEN
  DO t=1,Solver % Mesh % NUmberOfBoundaryElements
     Element => GetBoundaryElement(t)
     IF ( GetElementFamily() == 1 ) CYCLE
     NodeIndexes => Element % NodeIndexes
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()
     STIFF = 0.0D00
     FORCE = 0.0D00

     IF (COMP==DIM) THEN
     ! only for the surface nodes
        dd = SUM(ABS(Depth(Depthperm(NodeIndexes(1:n)))))
        IF (dd < 1.0e-6) THEN
           NodalU(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+1)
           NodalV = 0.0D0
           IF (DIM==3) NodalV(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+2)
           CALL LocalMatrixBCW (  STIFF, FORCE, Element, n, NodalU, NodalV ) 
        END IF
     ELSE IF (COMP==DIM+1) THEN
            CALL LocalMatrixBCP(  STIFF, FORCE, Element, n, NodalDensity, &
                    NodalGravity )
     END IF
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  END IF

  CALL DefaultFinishAssembly()

  ! Dirichlet 
     CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
          ComponentName('SSAFlow',COMP), 1,1, Permutation )
  
  !Solve the system
  Norm = DefaultSolve()

  ! Save the solution on the right variable
         DO i = 1, Model % Mesh % NumberOfNodes
           IF (VeloPerm(i)>0) THEN
           Velocity ((DIM+1)*(VeloPerm(i)-1) + COMP) = VariableValues(Permutation(i)) 
           END IF
         END DO 

  END DO ! Loop p

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixUV(  STIFF, FORCE, Element, n ) 
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ 
    LOGICAL :: Stat
    INTEGER :: t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()


    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO

    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUV
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixW(  STIFF, FORCE, Element, n, VeloU, VeloV)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), VeloU(:), VeloV(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ, &
                     dU2dxz, dV2dyz
    LOGICAL :: Stat
    INTEGER :: t, p,q , DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    DIM = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .TRUE. )

       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO

       dU2dxz = SUM(VeloU(1:n)*ddBasisddx(1:n,1,dim))
       dV2dyz = 0.0d0
       IF (DIM==3) dV2dyz = SUM(VeloV(1:n)*ddBasisddx(1:n,2,3))
       

       FORCE(1:n) = FORCE(1:n) + (dU2dxz + dV2dyz) * IP % s(t) * detJ * Basis(1:n) 

    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixW

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixP(  STIFF, FORCE, Element, n)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ
    LOGICAL :: Stat
    INTEGER :: t, p,q ,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixP
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCW(  STIFF, FORCE, Element, n, VeloU, VeloV )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), veloU(:), veloV(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ, Normal(3), grad, dUdx, dVdy  
    LOGICAL :: Stat
    INTEGER :: t, DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    DIM = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       dUdx = SUM( VeloU(1:n) * dBasisdx(1:n,1) )
       dVdy = 0.0e0
       IF (DIM==3) dVdy = SUM( VeloV(1:n) * dBasisdx(1:n,2) )

       grad = - (dUdx + dVdy) 

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCW
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCP(  STIFF, FORCE, Element, n, Density, & 
                      Gravity)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), density(:), Gravity(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3), rho, g, grad
    LOGICAL :: Stat
    INTEGER :: t, dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
       rho = SUM( Density(1:n) * Basis(1:n) )

       grad = - rho * g 

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCP
!------------------------------------------------------------------------------
END SUBROUTINE SSASolver
!------------------------------------------------------------------------------


