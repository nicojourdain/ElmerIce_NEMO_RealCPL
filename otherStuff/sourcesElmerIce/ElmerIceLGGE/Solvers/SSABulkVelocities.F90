!----------------------- BE CAREFUL ---------------------------------------------------------
!
! The solver was modified in December 2015 by Julien B.
! Now the solver compute the 3D (resp. 2D) Velocity field only from the 2D (resp. 1D) SSA Solution
! and leave an empty space ((u,v,w,free space) in 3D) for pressure which is filled by the new ComputeDevStressSSA
!
! This is because the expression of pressure in SSA is not only the hydrostatic pressure but the deviatoric stresses Sxx and Syy must be withdrawn (see SSA documentation)

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~Description of the solver before modification~~~~~~~~~~~~~~~~~~~~~~~

!!! Compute the 3D (resp. 2D) Velocity/Pressure field from the 2D (resp. 1D) SSA Solution
!
!  Execute this solver on the 2D bottom boundary where the SSA solution is computed
!  REQUIRE an ExtrudedMesh along the vertical direction
!  There is a consistency check that z varies between Zb and Zs.
!
!  OUTPUT Variable:
!    SSABulkVelocity (DOFS=dim+1 => (u,v,w,p)) (EXPORT THIS VARIABLE in a solver executed on the Bulk)
!  
!  BE CAREFUL : Pressure is not calculated in this solver but in the solver
!  ComputeDevStressSSA
!
!    OPTIONAL: if the varaibles are found
!      DZbDx,DZsDx (dofs=dim) : the horizontal derivatives of Zb and Zs
!
!  INPUT Variables:
!    SSAVelocity     (SSADOFs=DOFs-2)
!    Zb 
!    Zs
!    DZbDt (if Transient simulation)
!    DZsDt (if Transient simulation)
!
!  INPUT Body Forces: 
!    "Flow Body Force 2" (if dim=2);"Flow Body Force 3" (if dim=3)
!    "Top Surface Accumulation" (0 by default)
!    "Bottom Surface Accumulation" (0 by default)
!
!  INPUT Material Properties:
!    "Density"  use "SSA Mean Density" if not found => used to compute pressure
!       WARNING "Density" HAS TO BE DEFINED in the SAME MATERIAL SECTION than "SSA Mean Density"
! 
!  INPUT Solver Parameters:
!     "Active Coorinate = Integer " 2 or 3 the direction of extrusion
!
!------------------------------------------------------------------------------
SUBROUTINE SSABulkVelocities( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE MeshUtils
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: BodyForce,Material
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: Var,SSAVar,Zb,Zs,DZbDt,DZsDt,BulkVar
  TYPE(Element_t), POINTER :: Element

  REAL(KIND=DP),DIMENSION(:),POINTER :: Values
  REAL(KIND=DP),DIMENSION(:),ALLOCATABLE,SAVE :: DZb,DZs
  REAL(KIND=DP),DIMENSION(:),ALLOCATABLE,SAVE :: as,ab,Gravity,Density

  REAL(KIND=DP) :: s,b,bmesh,smesh,z,H
  REAL(KIND=DP) :: epss,epsb
  REAL(KIND=DP) :: rho,g
  REAL(KIND=DP) :: wb,ws,p
  REAL(KIND=DP),PARAMETER :: TOL=1.0d-6

  INTEGER,DIMENSION(:),POINTER :: Perm
  INTEGER, POINTER,SAVE :: BotPointer(:),TopPointer(:),UpPointer(:)
  INTEGER, POINTER,SAVE :: NodeIndexes(:)

  INTEGER :: M,N
  INTEGER :: DOFs,SSADOFs
  INTEGER :: node,botnode,topnode
  INTEGER :: dim
  INTEGER :: i,k,t

  LOGICAL,SAVE :: Initialized = .FALSE.
  LOGICAL :: Found

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='SSABulkVelocities'

!------------------------------------------------------------------------------
  Mesh => Model % Mesh


!!! get the 3D SSABulkVelocity variable
  BulkVar => VariableGet( Model % Mesh % Variables, 'SSABulkVelocity' )
  IF (ASSOCIATED(BulkVar)) THEN
     Values => BulkVar % Values
     Perm => BulkVar % Perm
     DOFs = BulkVar % Dofs
  ELSE
     Message='SSABulkVelocity not found'
     CALL FATAL(SolverName,Message)
  END IF
!!! get dimension of the mesh
!!!   should be consistent with DOFs
  dim=CoordinateSystemDimension()
  Message='CoordinateSystemDimension not consistent with SSABulkVelocity%Dofs'
  IF (dim.NE.(DOFs-1)) CALL FATAL(SolverName,Message)

!!! Do some initialisation/allocation
  IF ((.NOT.Initialized).OR.Mesh%Changed) THEN
    ! Choose active direction coordinate and set corresponding unit vector
    !---------------------------------------------------------------------
    PSolver => Solver
    CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, BotNodePointer = BotPointer , &
                                TopNodePointer = TopPointer, UpNodePointer = UpPointer)

    !! ALLOCATE
    IF (Initialized) DEALLOCATE(DZb,DZs,&
                                as,ab,&
                                Gravity,Density)

    M=Model%Mesh%NumberOfNodes
    N=Model % MaxElementNodes
    ALLOCATE(DZb(DOFS*M),DZs(DOFS*M),&
             as(N),ab(N), &
             Gravity(N),Density(N))

    Initialized = .TRUE.
  END IF

!!! get 2D SSAVelocity
  SSAVar => VariableGet( Model % Mesh % Variables, 'SSAVelocity' )
  IF (ASSOCIATED(SSAVar)) THEN
     SSADOFs=SSAVar%Dofs
     ! SSAVar%Dofs should be DOFs-2
     IF (SSADOFs.NE.(DOFs-2)) THEN
           Message='SSAVar%Dofs inconstistent with SSABulkVelocity%Dofs'
           CALL FATAL(SolverName,Message)
     END IF
  ELSE
     Message='SSAVelocity  not found'
     CALL FATAL(SolverName,Message)
  END IF

!!! get Zb and Zs
  zb => VariableGet( Model % Mesh % Variables, 'Zb')
  IF (.NOT.ASSOCIATED(zb)) THEN
     Message='Zb  not found'
     CALL FATAL(SolverName,Message)
  END IF
  zs => VariableGet( Model % Mesh % Variables, 'Zs')
  IF (.NOT.ASSOCIATED(zs)) THEN
     Message='Zs  not found'
     CALL FATAL(SolverName,Message)
  END IF
!!
  IF (Transient) THEN
    DzbDt => VariableGet( Model % Mesh % Variables, 'DZbDt')
    IF (.NOT.ASSOCIATED(DzbDt)) THEN
       Message='DZbDt  not found'
       CALL FATAL(SolverName,Message)
    END IF
    DzsDt => VariableGet( Model % Mesh % Variables, 'DZsDt')
    IF (.NOT.ASSOCIATED(DzsDt)) THEN
       Message='DZsDt  not found'
       CALL FATAL(SolverName,Message)
    END IF
  END IF

!! compute Dzb and DZs (required to compute bottom and top vertical velocity)
  CALL ComputeNodalSlope(Model,Solver,Zb,SSADOFs,DZb)
  CALL ComputeNodalSlope(Model,Solver,Zs,SSADOFs,DZs)

!!
 Do t=1,Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes

    BodyForce => GetBodyForce(Element)
    IF (.NOT.ASSOCIATED(BodyForce)) THEN
    END IF
    Material => GetMaterial(Element)
    IF (.NOT.ASSOCIATED(BodyForce)) THEN
    END IF

    Gravity=0._dp
    IF (SSADOFs==1) THEN
       Gravity(1:n) = ListGetReal(BodyForce, 'Flow BodyForce 2', n, NodeIndexes, Found)
    ELSE
       Gravity(1:n) = ListGetReal(BodyForce, 'Flow BodyForce 3', n, NodeIndexes, Found)
    END IF
    as=0._dp
    as(1:n) = ListGetReal(BodyForce, 'Top surface accumulation', n, NodeIndexes, Found)
    ab=0._dp
    ab(1:n) = ListGetReal(BodyForce, 'Bottom surface accumulation', n, NodeIndexes, Found)
    Density = 0._dp
    Density(1:n) = ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,Found)

    Do i=1,n 
       g=abs(Gravity(i))

       botnode=BotPointer(NodeIndexes(i))
       topnode=TopPointer(NodeIndexes(i))
       node=BotPointer(NodeIndexes(i))

       b=zb%Values(zb%Perm(botnode))
       s=zs%Values(zs%Perm(botnode))
       H=s-b

       !Consistency check
       IF (dim.eq.2) then
          bmesh=Mesh%nodes%y(botnode)
          smesh=Mesh%nodes%y(topnode)
       ELSE
          bmesh=Mesh%nodes%z(botnode)
          smesh=Mesh%nodes%z(topnode)
       END IF
       epsb=abs(b-bmesh)
       epss=abs(s-smesh)
       IF ((epsb.GT.TOL).OR.(epss.GT.TOL)) THEN
          WRITE(Message,*) 'Topography inconsistency detected',b,bmesh,s,smesh,epsb,epss
          CALL FATAL(SolverName,Message)
       END IF
       
       ws=-as(i)
       wb=ab(i)
       Do k=1,SSADOFs
         ws=ws+SSAVar%Values(SSADOFs*(SSAVar%Perm(botnode)-1)+k)*DZs(SSADOFS*(botnode-1)+k)
         wb=wb+SSAVar%Values(SSADOFs*(SSAVar%Perm(botnode)-1)+k)*DZb(SSADOFS*(botnode-1)+k)
       End do
       IF (Transient) THEN
         ws=ws+DZsDt%Values(DZsDt%Perm(botnode))
         wb=wb+DZbDt%Values(DZbDt%Perm(botnode))
       END IF
       p=(ws-wb)/H

       Do while (.True.)
         IF (dim.eq.2) then
          z=Mesh%nodes%y(node)
         else
          z=Mesh%nodes%z(node)
         endif
         Do k=1,SSAVar%Dofs
           Values(DOFs*(Perm(node)-1)+k)=SSAVar%Values(SSADOFs*(SSAVar%Perm(botnode)-1)+k) !u,v
         End Do
         Values(DOFs*(Perm(node)-1)+DOFs-1)= wb +&
            (z-b)*p !w
         if (node.eq.topnode) exit
         node=UpPointer(node)
       End do
    End do

 End Do

CONTAINS

! *****************************************************************************
! Compute DVar/Dx(1:STDOFs), 
! IF The varaible <'D',Var%Name,'Dx'> is found save the results
SUBROUTINE ComputeNodalSlope(Model,Solver,Var,STDOFs,DVar)
!------------------------------------------------------------------------------
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  TYPE(Variable_t), POINTER  :: Var,
!     INPUT: The variable
!
!  INTEGER :: STDOFs
!     INPUT: dimension of the probleme (1D (flowline SSA) or 2D (plane view SSA))
!
!  REAL(KIND=dp), DIMENSION(:) :: DVar  
!     OUTPUT: Nodal spatial derivative of var
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t),TARGET :: Solver
  TYPE(Model_t) :: Model
  TYPE(Variable_t), POINTER :: Var
  INTEGER :: STDOFs
  REAL(KIND=dp), DIMENSION(:) :: DVar
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: DVarDx
  LOGICAL,SAVE :: FEConsistent=.True. ! Do FE consistent averages
!!  
  TYPE(Nodes_t),SAVE   :: ElementNodes
  TYPE(Element_t),POINTER ::  Element
  TYPE(GaussIntegrationPoints_t) :: IP

  LOGICAL,SAVE :: AllocationsDone = .FALSE.
  LOGICAL :: Found
  LOGICAL :: stat

  LOGICAL,ALLOCATABLE,SAVE :: ActiveNode(:)

  INTEGER :: i,j,k,l,t
  INTEGER :: n,m
  INTEGER :: istat
  INTEGER :: proc,ierr

  INTEGER, POINTER :: Perm(:)
  INTEGER, POINTER,SAVE :: NodeIndexes(:)

  REAL(KIND=dp), POINTER :: Values(:)
  REAL(KIND=dp) :: U,V,W
  REAL(KIND=dp) :: dsdx(2)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: grad(:),weight(:)
  REAL(KIND=dp),allocatable,save :: Basis(:),dBasisdx(:,:),ddBasisddx(:,:,:)
  REAL(KIND=dp) :: detJ

  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: SolverName
  CHARACTER(LEN=MAX_NAME_LEN) :: DVarName

  TYPE lbuff_t
        INTEGER, ALLOCATABLE :: buff(:)
        REAL(KIND=dp), ALLOCATABLE :: values(:)
  END TYPE lbuff_t
  INTEGER, POINTER :: nlist(:)
  TYPE(lbuff_t), ALLOCATABLE :: n_index(:)
  REAL(KIND=dp), ALLOCATABLE :: nbuff(:)
  INTEGER, ALLOCATABLE :: n_count(:), gbuff(:), n_comp(:)
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

!------------------------------------------------------------------------------
  WRITE(SolverName, '(A)') 'ComputeNodalSlope'

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
        IF (ASSOCIATED(Var)) THEN
           Values => Var % Values
           Perm => Var % Perm
        ELSE
           WRITE(Message, '(A)') 'Input Variable not associated'
           CALL FATAL(SolverName,Message)
        END IF

        write(DVarName,'(A,A,A)') 'D',trim(Var%Name),'Dx'
        DVarDx => VariableGet( Model % Mesh % Variables, trim(DVarName))
        IF (ASSOCIATED(DVarDx)) THEN
           IF (DVarDx%DOFs.NE.STDOFS) THEN
             WRITE(Message, '(A,A,A)') 'Found Variable ',trim(DVarName),' but Dofs &
               &inconsistent'
             CALL FATAL(SolverName,Message)
           END IF
         END IF

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN

     N = Model % MaxElementNodes
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(ElementNodes % x, &
                       ElementNodes % y, ElementNodes % z, &
                       Basis,dBasisdx,ddBasisddx,&
                       grad, &
                       ActiveNode, &
                       weight)

     ALLOCATE(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N), &
             Basis(N),dBasisdx(N,3),ddBasisddx(N,3,3),&
              grad(STDOFs*M), &
              ActiveNode(M), &
              weight(M),&
           STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

          ActiveNode=.False.
          grad=0._dp
          weight=0.0

   !!!!!! Compute surface slope element wise
          Do t=1,Solver % NumberOfActiveElements
             Element => GetActiveElement(t)
             IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
             n = GetElementNOFNodes()
             NodeIndexes => Element % NodeIndexes

             ! set coords of highest occuring dimension to zero (to get correct path element)
             !-------------------------------------------------------------------------------
             ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
             IF (STDOFs == 1) THEN !1D 
               ElementNodes % y(1:n) = 0.0_dp
               ElementNodes % z(1:n) = 0.0_dp
             ELSE IF (STDOFs == 2) THEN !2D 
               ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
               ElementNodes % z(1:n) = 0.0_dp
             ELSE
                WRITE(Message,'(a,i1,a)')&
                  'It is not possible to compute slope with DOFs=',&
                  STDOFs, ' . Aborting'
                CALL Fatal( SolverName, Message)
                STOP
             END IF

             Do i=1,n
                  k=NodeIndexes(i)
                  ActiveNode(k)=.True.

                  IF (FEConsistent) THEN
                     IP = GaussPoints(Element)
                     DO j=1,IP % n
                        stat = ElementInfo(Element, ElementNodes, IP % U(j), &
                                           IP % v(j), IP % W(j), detJ, &
                                           Basis,dBasisdx, ddBasisddx,.FALSE. )

                        dsdx(1)=SUM( Values(Perm(NodeIndexes(1:n))) * dBasisdx(1:n,1) )
                        IF (STDOFs == 2) &
                           dsdx(2)=SUM( Values(Perm(NodeIndexes(1:n))) * dBasisdx(1:n,2) )

                        weight(k)=weight(k)+IP%s(j)*detJ*Basis(i)

                        grad(STDOFs*(k-1)+1) = grad(STDOFs*(k-1)+1) + &
                                               dsdx(1)*IP%s(j)*detJ*Basis(i)
                        IF (STDOFs == 2) grad(STDOFs*(k-1)+2) =grad(STDOFs*(k-1)+2) + &
                                                dsdx(2)*IP%s(j)*detJ*Basis(i)
                      END DO !IPs
                  ELSE

                     U=Element % TYPE % NodeU(i)
                     V=Element % TYPE % NodeV(i)
                     W=0.0

                     stat = ElementInfo( Element, ElementNodes, U, V, &
                         W,  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

                     weight(k)=weight(k)+1.0_dp


                     grad(STDOFs*(k-1)+1) = grad(STDOFs*(k-1)+1) + &
                             SUM( Values(Perm(NodeIndexes(1:n))) * dBasisdx(1:n,1) )

                     IF (STDOFs == 2) grad(STDOFs*(k-1)+2) =grad(STDOFs*(k-1)+2) + &
                             SUM( Values(Perm(NodeIndexes(1:n))) * dBasisdx(1:n,2) )
                 END IF
              End do !elements nodes

          End do !elements

          IF (ParEnv % PEs>1 ) THEN  !! if parallel need to sum values at interfaces
                                     !! here is a copy of what is done in
                                     !SolverUtils.src to average boundary normals
             ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
             n_count = 0

             DO i=1,Solver%Mesh % NumberOfNodes
                IF (.NOT.ActiveNode(i)) CYCLE
                IF (.NOT.Solver % Mesh % ParallelInfo % INTERFACE(i) ) CYCLE
  
                nlist => Solver%Mesh % ParallelInfo % NeighbourList(i) % Neighbours
                DO j=1,SIZE(nlist)
                   k = nlist(j)+1
                   IF ( k-1 == ParEnv % myPE ) CYCLE
                   n_count(k) = n_count(k)+1
                END DO
             END DO
             DO i=1,ParEnv % PEs
                IF ( n_count(i)>0 ) &
                  ALLOCATE( n_index(i) % buff(n_count(i)), &
                            n_index(i) % values((STDOFs+1)*n_count(i)) )
             END DO

             n_count = 0
             DO i=1,Model % NumberOfNodes
               IF (.NOT.ActiveNode(i)) CYCLE
               IF (.NOT.Solver % Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

               nlist =>Solver% Mesh % ParallelInfo % NeighbourList(i) % Neighbours
               DO j=1,SIZE(nlist)
                 k = nlist(j)+1
                 IF ( k-1 == ParEnv % myPE ) CYCLE
                 n_count(k) = n_count(k)+1
                 n_index(k) % buff(n_count(k)) = Solver%Mesh % Parallelinfo % &
                 GlobalDOFs(i)
                 n_index(k) % values((STDOFs+1)*(n_count(k)-1)+1)=grad(STDOFs*(i-1)+1)
                 if (STDOFS.EQ.2) &
                  n_index(k) % Values((STDOFs+1)*(n_count(k)-1)+2)=grad(STDOFs*(i-1)+2)
                 n_index(k) % values((STDOFs+1)*(n_count(k)-1)+STDOFs+1)=weight(i)
               END DO
             END DO

             DO i=1,ParEnv % PEs
               IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
                CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                  900, MPI_COMM_WORLD, ierr )
                IF ( n_count(i)>0 ) THEN
                 CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                  901, MPI_COMM_WORLD, ierr )
                 CALL MPI_BSEND( n_index(i) % values, (STDOFs+1)*n_count(i), MPI_DOUBLE_PRECISION, &
                    i-1,  902, MPI_COMM_WORLD, ierr )
               END IF
              END IF
            END DO
            DO i=1,ParEnv % PEs
               IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff, n_index(i) % values)

               IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
                  CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     900, MPI_COMM_WORLD, status, ierr )
                 IF ( n>0 ) THEN
                   proc = status(MPI_SOURCE)
                   ALLOCATE( gbuff(n), nbuff((STDOFs+1)*n) )
                   CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                     901, MPI_COMM_WORLD, status, ierr )

                   CALL MPI_RECV( nbuff, (STDOFs+1)*n, MPI_DOUBLE_PRECISION, proc, &
                     902, MPI_COMM_WORLD, status, ierr )

                   DO j=1,n
                     k = SearchNodeL(Solver% Mesh % ParallelInfo, gbuff(j), Solver%Mesh % NumberOfNodes )

                     IF ( k>0 ) THEN
                       !n_comp(k) = n_comp(k)+1
                      ! IF ( l>0 ) THEN
                         grad(STDOFs*(k-1)+1)=grad(STDOFs*(k-1)+1)+nbuff((STDOFs+1)*(j-1)+1)
                         IF (STDOFs.EQ.2) &
                            grad(STDOFs*(k-1)+2)=grad(STDOFs*(k-1)+2)+nbuff((STDOFs+1)*(j-1)+2)
                         weight(k)=weight(k)+nbuff((STDOFs+1)*(j-1)+STDOFs+1)
                      ! END IF
                     END IF
                   END DO
                   DEALLOCATE(gbuff, nbuff)
                 END IF
               END IF
           END DO
           DEALLOCATE( n_index, n_count )
           !DEALLOCATE(n_comp)

       END IF !end do parallel reduction

        Do t=1,Solver % Mesh % NumberOfNodes
             IF (.NOT.ActiveNode(t)) CYCLE 
             DVar(STDOFs*(t-1)+1) = grad(STDOFs*(t-1)+1) / weight(t)
             IF (ASSOCIATED(DVarDx)) THEN
                 DVarDx%Values(STDOFs*(DVarDx%Perm(t)-1)+1)=DVar(STDOFs*(t-1)+1)
             END IF
             IF (STDOFs == 2) THEN
                DVar(STDOFs*(t-1)+2)=grad(STDOFs*(t-1)+2) / weight(t)
                IF (ASSOCIATED(DVarDx)) THEN
                    DVarDx%Values(STDOFs*(DVarDx%Perm(t)-1)+2)=DVar(STDOFs*(t-1)+2)
                END IF
             END IF
          END DO

!------------------------------------------------------------------------------
END SUBROUTINE ComputeNodalSlope
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE SSABulkVelocities
!------------------------------------------------------------------------------
