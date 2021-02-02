!-----------------------------------------------------------------------------
!> 
!> 
!> 
SUBROUTINE SSATEMPICE_Init(Model,Solver,dt,TransientSimulation )
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t),TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh3D
  TYPE(Valuelist_t),POINTER :: SolverParams
  LOGICAL :: Found
  CHARACTER(LEN=MAX_NAME_LEN) :: MeshName,DefaultName="3DGrid"
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='SSATEMPICE'

  SolverParams => Solver % Values

  MeshName=ListGetString(SolverParams,'Extruded Mesh Name',Found)
  IF (.NOT.Found) MeshName=DefaultName

  Mesh3D => CurrentModel % Meshes
  IF (TRIM(Mesh3D % Name) .NE. TRIM(MeshName)) THEN
     Mesh3D => CurrentModel % Meshes % Next
     DO WHILE (ASSOCIATED(Mesh3D))
       IF (TRIM(Mesh3D%Name) == TRIM(MeshName)) EXIT
     END DO
  ENDIF

  IF (.NOT.ASSOCIATED(Mesh3D)) &
      CALL FATAL(SolverName,'Extruded Mesh Not found')

  PRINT *,'Solver Mesh assoxciated with 3D:',ASSOCIATED(Solver%Mesh,Mesh3D)
  IF (.NOT.ASSOCIATED(Solver%Mesh,Mesh3D)) &
     Solver%Mesh => Mesh3D


END SUBROUTINE SSATEMPICE_Init

!------------------------------------------------------------------------------
SUBROUTINE SSATEMPICE( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
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
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter
  LOGICAL :: Found
!------------------------------------------------------------------------------
  TYPE(Mesh_t),POINTER :: Mesh2D,Mesh3D 
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t),POINTER :: HSol,BMBSol,DHDTSol,SSASol
  REAL(KIND=dp),POINTER :: zcoord(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: MeshName,DefaultName="3DGrid"
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='SSATEMPICE'

!!!! Variables for DetectExtrudedStructure
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Variable_t), POINTER :: Var
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),Upointer(:),DownPointer(:)

  INTEGER :: DIM


!!!!!!!!
  SolverParams => Solver % Values

  MeshName=ListGetString(SolverParams,'Extruded Mesh Name',Found)
  IF (.NOT.Found) MeshName=DefaultName

  Mesh2D => CurrentModel % Meshes

  Mesh3D => CurrentModel % Meshes
  IF (TRIM(Mesh3D % Name) .NE. TRIM(MeshName)) THEN
     Mesh3D => CurrentModel % Meshes % Next
     DO WHILE (ASSOCIATED(Mesh3D))
       IF (TRIM(Mesh3D%Name) == TRIM(MeshName)) EXIT
     END DO
  ENDIF

  IF (.NOT.ASSOCIATED(Mesh3D)) &
      CALL FATAL(SolverName,'Extruded Mesh Not found')
  CurrentModel % Mesh => Mesh3D
  CurrentModel % Variables => Mesh3D % Variables

  PSolver => Solver
  CALL DetectExtrudedStructure( Mesh3D , PSolver, Var, &
     TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
     UpNodePointer = Upointer , DownNodePointer = DownPointer)

!!!!!!!!!
  DIM = Mesh3D%MeshDim
  IF (DIM.EQ.2) THEN
     zcoord=> Mesh3D%Nodes%y
  ELSE IF (DIM.EQ.3) THEN
     zcoord=> Mesh3D%Nodes%z
  ELSE
    CALL FATAL(SolverName,'Invalid mesh dimension')
  ENDIF

  IF ((MINVAL(zcoord).NE.0.d0).OR.(MAXVAL(zcoord).NE.1.d0)) &
     CALL FATAL(SolverName,'This solver is only for reduced vertical coordinates')

!!!!!!!!!
  SSASol => VariableGet(Mesh2D%Variables,'SSAVelocity',UnFoundFatal=.TRUE.)
  IF (SSASol%DOFs.NE.(DIM-1)) THEN
    write(Message,'(a,i0,a,i0)') 'Pb with dimension SSA:',SSASol%DOFs,' should be:',DIM-1
    CALL FATAL(SolverName,trim(Message))
  END IF
  HSol   => VariableGet(Mesh2D%Variables,'H',UnFoundFatal=.TRUE.)
  DHDTSol=> VariableGet(Mesh2D%Variables,'dhdt',UnFoundFatal=.TRUE.)
  BMBSol => VariableGet(Mesh2D%Variables,'bmb',UnFoundFatal=.TRUE.)
!!

  CALL DefaultStart(Solver)
  
  maxiter = ListGetInteger( SolverParams,&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    
    CALL ResetTimer(trim(SolverName))

    ! System assembly:
    !----------------
    CALL DefaultInitialize(Solver)

    Active = GetNOFActive(Solver)
    DO t=1,Active
      Element => GetActiveElement(t,Solver)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)
      !PRINT *,n,nd,nb
      CALL LocalMatrix(  Element, n, nd+nb )
    END DO

    CALL DefaultFinishBulkAssembly(Solver)


!!!!!!
    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBC(  Element, n, nd+nb )
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
 

    CALL CheckTimer(trim(SolverName),Level=3)
    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    CALL CheckTimer(trim(SolverName),Level=3,Delete=.TRUE.)

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO

  CALL DefaultFinish()
 
!!!
  CurrentModel % Mesh => Mesh2D
  CurrentModel % Variables => Mesh2D % Variables
!!!
CONTAINS

! Assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Kappa(n), Rho(n),Cond(n)
    REAL(KIND=dp) :: Velo(3,n),NodalH(n),NodalDHDT(n),NodalBMB(n)
    REAL(KIND=dp) :: rC,H,Vgauss(3),D,LoadAtIP,xsi
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight
    REAL(KIND=dp) :: dNodalBasisdxsi(nd,n)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: M,A,F
    REAL(KIND=dp) :: U,V,W
    REAL(KIND=dp) :: hK,mK,VNorm,Pe
    REAL(KIND=dp) :: SW(nd)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t),SAVE :: Nodes
    INTEGER,POINTER :: NodeIndexes(:)

    LOGICAL :: Stabilize=.TRUE.
    REAL(KIND=dp) :: Tau
!------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    SW=0._dp

    NodeIndexes => Element % NodeIndexes

    Material => GetMaterial(Element)
    Kappa(1:n)=ListGetReal(Material,'Heat Conductivity',n, NodeIndexes,UnFoundFatal=.TRUE.)
    Cond(1:n)=ListGetReal(Material,'Heat Capacity',n, NodeIndexes,UnFoundFatal=.TRUE.)
    Rho(1:n)=ListGetReal(Material,'Ice Density',n, NodeIndexes,UnFoundFatal=.TRUE.)
    
    
    Velo = 0._dp
    DO i=1,dim-1
      Velo(i,1:n)=SSASol%Values(SSASol%DOFs*(SSASol%Perm(BotPointer(NodeIndexes(1:n)))-1)+i) 
    END DO
    NodalH(1:n)=HSol%Values(HSol%Perm(BotPointer(NodeIndexes(1:n))))
    NodalDHDT(1:n)=DHDTSol%Values(DHDTSol%Perm(BotPointer(NodeIndexes(1:n))))
    NodalBMB(1:n)=BMBSol%Values(BMBSol%Perm(BotPointer(NodeIndexes(1:n))))

    IF (Stabilize) THEN
      hK = element % hK
      mK = element % StabilizationMK

      dNodalBasisdxsi = 0._dp
      DO p=1,n
        u = Element % TYPE % NodeU(p)
        v = Element % TYPE % NodeV(p)
        w = Element % TYPE % NodeW(p)
        stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
        dNodalBasisdxsi(1:nd,p) = dBasisdx(1:nd,dim)
      END DO
    END IF

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      LoadAtIP=0._dp    

      H=SUM(Basis(1:n)*NodalH(1:n))
      rC=SUM(Basis(1:n)*Rho(1:n)*Cond(1:n))

      ! The vertical diffusion term:
      ! -----------------------------------
      D=SUM(Basis(1:n)*Kappa(1:n))/(H*H)

      !! convection velocity
      !! Horizontal
      Vgauss(1:dim-1) = MATMUL(Velo(1:dim-1,1:n),Basis(1:n))
      !! vertical
      xsi=SUM(Basis(1:n)*zcoord(NodeIndexes(1:n)))
      Vgauss(dim)=0._dp
      DO i=1,dim-1
        Vgauss(dim)=Vgauss(dim)+H*SUM(Velo(i,1:n)*dBasisdx(1:n,i))+Vgauss(i)*SUM(NodalH(1:n)*dBasisdx(1:n,i))
      END DO
      Vgauss(dim)=Vgauss(dim)+SUM(Basis(1:n)*NodalDHDT(1:n))
      Vgauss(dim)=-xsi*Vgauss(dim)
      Vgauss(dim)=Vgauss(dim)+SUM(NodalBMB(1:n)*Basis(1:n))    !have to decide the sign for bmb!!
      Vgauss(dim)=Vgauss(dim)/H

      ! Stabilization terms
      !------------------------
      IF (Stabilize) THEN
         DO p=1,nd
            SW(p) = 0.0_dp
            DO i=1,dim
               SW(p) = SW(p) + Vgauss(i) * dBasisdx(p,i)
            END DO
         END DO
         VNorm = SQRT( SUM(Vgauss(1:dim)**2) )
         Pe  = MIN( 1.0D0, mK*hK*rc*VNorm/(2*D)) 
         IF ( VNorm /= 0.0 ) THEN
           Tau = hK * Pe / (2 * rc * VNorm)
         END IF
      END IF


      DO p=1,nd
        DO q=1,nd

          ! diffusion term ( D * du/dxsi,dv/dxsi)
          ! -----------------------------------
          A =  D * dBasisdx(p,dim) * dBasisdx(q,dim)

          ! advection term (C*grad(u),v)
          ! -----------------------------------
          A = A + rC * SUM(Vgauss(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! time derivative (rho*du/dt,v):
          ! ------------------------------
          M = rC * Basis(q) * Basis(p)

          
          IF (Stabilize) THEN
             ! time derivative
             M = M + Tau * rC * Basis(q) * SW(p)
             ! diffusion: 
             !  -  -dD/dxsi dT/dxsi
             A = A - Tau * (SUM(Kappa(1:n)*dBasisdx(1:n,dim))/(H*H)) * dBasisdx(q,dim) * SW(p)
             ! - -D d/dxsi(dT/dxsi)
             A = A - Tau * D * SUM(dNodalBasisdxsi(q,1:n)*dBasisdx(1:n,dim)) * SW (p)
             ! advection 
             A = A + Tau * rC * SW(q) * SW (p)
          END IF

          MASS(p,q) = MASS(p,q) + Weight * M
          STIFF(p,q) = STIFF(p,q) + Weight * A
        END DO
      END DO

      ! right-hand side : 
      ! ------------------------------
      DO p=1,nd
         F = LoadAtIP * Basis(p)
         IF (Stabilize) F = F + Tau * LoadAtIP * SW (p)

         FORCE(p) = FORCE(p) + Weight * F
      END DO
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n),F,Weight
    REAL(KIND=dp) :: NodalH(n),H
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC
    INTEGER,POINTER :: NodeIndexes(:)

    TYPE(Nodes_t),SAVE :: Nodes
!------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    NodeIndexes => Element % NodeIndexes(:)
    NodalH(1:n)=HSol%Values(HSol%Perm(BotPointer(NodeIndexes(1:n))))
    Flux(1:n)  = GetReal( BC,'Heat flux', Found )

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux -kdT/dxsi=-(1/H)kdT/dz :
      ! -----------
      H = SUM(Basis(1:n)*NodalH(1:n))
      F = SUM(Basis(1:n)*flux(1:n))/H

      FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
    END DO
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

! Perform static condensation in case bubble dofs are present
!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE SSATEMPICE
!------------------------------------------------------------------------------
