
!------------------------------------------------------------------------------
!> Map results from mesh to mesh. The from-Mesh is stored in an octree from 
!> which it is relatively fast to find the to-nodes. When the node is found
!> interpolation is performed. Optionally there may be an existing projector
!> that speeds up the interpolation.
!------------------------------------------------------------------------------
     SUBROUTINE InterpolateMeshToMeshR( OldMesh, NewMesh, OldVariables, &
            NewVariables, UseQuadrantTree, Projector, MaskName )
!------------------------------------------------------------------------------
       USE Lists
       USE SParIterComm
       USE Interpolation
       USE CoordinateSystems
!-------------------------------------------------------------------------------
       TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
       TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
       LOGICAL, OPTIONAL :: UseQuadrantTree
       TYPE(Projector_t), POINTER, OPTIONAL :: Projector
       CHARACTER(LEN=*),OPTIONAL :: MaskName
!-------------------------------------------------------------------------------
       INTEGER, ALLOCATABLE :: perm(:), vperm(:)
       INTEGER, POINTER :: nperm(:)
       LOGICAL, ALLOCATABLE :: FoundNodes(:)
       TYPE(Mesh_t), POINTER :: nMesh
       TYPE(VAriable_t), POINTER :: Var, nVar
       INTEGER :: i,j,k,l,nfound,n,ierr,nvars,npart,proc,status(MPI_STATUS_SIZE)
       REAL(KIND=dp) :: myBB(6), eps2, dn
       REAL(KIND=dp), POINTER :: store(:)
       REAL(KIND=dp), ALLOCATABLE, TARGET :: astore(:),vstore(:,:), BB(:,:), &
             nodes_x(:),nodes_y(:),nodes_z(:), xpart(:), ypart(:), zpart(:)

       TYPE ProcRecv_t
         INTEGER :: n = 0
         REAL(KIND=dp), ALLOCATABLE :: nodes_x(:),nodes_y(:),nodes_z(:)
       END TYPE ProcRecv_t
       TYPE(ProcRecv_t),  ALLOCATABLE, TARGET :: ProcRecv(:)

       TYPE ProcSend_t
         INTEGER :: n = 0
         INTEGER, ALLOCATABLE :: perm(:)
       END TYPE ProcSend_t
       TYPE(ProcSend_t),  ALLOCATABLE :: ProcSend(:)

!-------------------------------------------------------------------------------
       INTERFACE
         SUBROUTINE InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, &
             NewVariables, UseQuadrantTree, Projector, MaskName, FoundNodes )
           USE Types
           TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
           TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
           LOGICAL, OPTIONAL :: UseQuadrantTree,FoundNodes(:)
           CHARACTER(LEN=*),OPTIONAL :: MaskName
           TYPE(Projector_t), POINTER, OPTIONAL :: Projector
         END SUBROUTINE InterpolateMeshToMeshQ
       END INTERFACE
!-------------------------------------------------------------------------------

      IF ( ParEnv % PEs<=1 ) THEN
!         CALL InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, &
!            NewVariables, UseQuadrantTree, Projector, MaskName )!
         CALL InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, &
            NewVariables, UseQuadrantTree )
         RETURN
      END IF

      ! Interpolate within our own partition, flag the points
      ! we found:
      ! -----------------------------------------------------
      ALLOCATE( FoundNodes(NewMesh % NumberOfNodes) ); FoundNodes=.FALSE.
      CALL InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, &
         NewVariables, UseQuadrantTree, MaskName=MaskName, FoundNodes=FoundNodes )

      ! special case "all found":
      !--------------------------
      n = COUNT(.NOT.FoundNodes); dn = n
      CALL SParActiveSUM(dn,2)
      IF ( dn==0 ) RETURN

      ! Exchange partition bounding boxes:
      ! ----------------------------------
      myBB(1) = MINVAL(OldMesh % Nodes % x)
      myBB(2) = MINVAL(OldMesh % Nodes % y)
      myBB(3) = MINVAL(OldMesh % Nodes % z)
      myBB(4) = MAXVAL(OldMesh % Nodes % x)
      myBB(5) = MAXVAL(OldMesh % Nodes % y)
      myBB(6) = MAXVAL(OldMesh % Nodes % z)

      eps2 = 0.1_dp * MAXVAL(myBB(4:6)-myBB(1:3))
      myBB(1:3) = myBB(1:3) - eps2
      myBB(4:6) = myBB(4:6) + eps2

      ALLOCATE(BB(6,ParEnv % PEs))
      DO i=1,ParEnv % PEs
        IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
        proc = i-1
        CALL MPI_BSEND( myBB, 6, MPI_DOUBLE_PRECISION, proc, &
                 999, MPI_COMM_WORLD, ierr )
      END DO
      DO i=1,COUNT(ParEnv % Active)-1
        CALL MPI_RECV( myBB, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                 999, MPI_COMM_WORLD, status, ierr )
        proc = status(MPI_SOURCE)
        BB(:,proc+1) = myBB
      END DO

      IF ( n==0 ) THEN
        DEALLOCATE(FoundNodes, BB)
        DO i=1,ParEnv % PEs
          IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
          proc = i-1
          CALL MPI_BSEND( n, 1, MPI_INTEGER, proc, &
                1001, MPI_COMM_WORLD, ierr )
        END DO
      ELSE
        ! Extract nodes that we didn't find from our own partition...
        ! ------------------------------------------------------------
        ALLOCATE( Perm(n), nodes_x(n), nodes_y(n),nodes_z(n) ); Perm=0
        j = 0
        DO i=1,NewMesh % NumberOfNodes
          IF ( FoundNodes(i) ) CYCLE
          j = j + 1
          perm(j) = i
          nodes_x(j) = NewMesh % Nodes % x(i)
          nodes_y(j) = NewMesh % Nodes % y(i)
          nodes_z(j) = NewMesh % Nodes % z(i)
        END DO
        DEALLOCATE(FoundNodes)

        ! ...and ask those from others
        ! -------------------------------
        ALLOCATE(ProcSend(ParEnv % PEs))
        DO i=1,ParEnv % PEs
          IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
          proc = i-1

          ! extract those of the missing nodes that are within the other
          ! partions bounding box:
          ! ------------------------------------------------------------
          myBB = BB(:,i)
          npart = 0
          DO j=1,n
            IF ( nodes_x(j)<myBB(1) .OR. nodes_x(j)>myBB(4) .OR. &
                 nodes_y(j)<myBB(2) .OR. nodes_y(j)>myBB(5) .OR. &
                 nodes_z(j)<myBB(3) .OR. nodes_z(j)>myBB(6) ) CYCLE
            npart = npart+1
          END DO
          ProcSend(proc+1) % n = npart
          IF ( npart>0 ) THEN
            ALLOCATE( xpart(npart),ypart(npart),zpart(npart),ProcSend(proc+1) % perm(npart) )
            npart = 0
            DO j=1,n
              IF ( nodes_x(j)<myBB(1) .OR. nodes_x(j)>myBB(4) .OR. &
                   nodes_y(j)<myBB(2) .OR. nodes_y(j)>myBB(5) .OR. &
                   nodes_z(j)<myBB(3) .OR. nodes_z(j)>myBB(6) ) CYCLE
              npart=npart+1
              ProcSend(proc+1) % perm(npart)=j
              xpart(npart) = Nodes_x(j)
              ypart(npart) = Nodes_y(j)
              zpart(npart) = Nodes_z(j)
            END DO
          END IF

          ! send count...
          ! -------------
          CALL MPI_BSEND( npart, 1, MPI_INTEGER, proc, &
                  1001, MPI_COMM_WORLD, ierr )

          IF ( npart==0 ) CYCLE

          ! ...and points
          ! -------------
          CALL MPI_BSEND( xpart, npart, MPI_DOUBLE_PRECISION, proc, &
                  1002, MPI_COMM_WORLD, ierr )
          CALL MPI_BSEND( ypart, npart, MPI_DOUBLE_PRECISION, proc, &
                  1003, MPI_COMM_WORLD, ierr )
          CALL MPI_BSEND( zpart, npart, MPI_DOUBLE_PRECISION, proc, &
                  1004, MPI_COMM_WORLD, ierr )

          DEALLOCATE(xpart,ypart,zpart)
        END DO
        DEALLOCATE(nodes_x,nodes_y,nodes_z,BB)
      END IF
       

      ! receive points from others:
      ! ----------------------------
      ALLOCATE(ProcRecv(Parenv % Pes))
      DO i=1,COUNT(ParEnv % Active)-1
        CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
              1001, MPI_COMM_WORLD, status, ierr )

        proc = status(MPI_SOURCE)
        ProcRecv(proc+1) % n = n

        IF ( n<=0 ) CYCLE

        ALLOCATE(ProcRecv(proc+1) % Nodes_x(n), &
              ProcRecv(proc+1) % Nodes_y(n),ProcRecv(proc+1) % Nodes_z(n))

        CALL MPI_RECV( ProcRecv(proc+1) % nodes_x, n, MPI_DOUBLE_PRECISION, proc, &
               1002, MPI_COMM_WORLD, status, ierr )
        CALL MPI_RECV( ProcRecv(proc+1) % nodes_y, n, MPI_DOUBLE_PRECISION, proc, &
               1003, MPI_COMM_WORLD, status, ierr )
        CALL MPI_RECV( ProcRecv(proc+1) % nodes_z, n, MPI_DOUBLE_PRECISION, proc, &
               1004, MPI_COMM_WORLD, status, ierr )
      END DO

      ! Check the received points and extract values for the to-be-interpolated-
      ! variables, if we have the points within our domain: 
      ! ------------------------------------------------------------------------
      DO i=1,ParEnv % PEs
        IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE

        proc = i-1
        n = ProcRecv(i) % n

        IF ( n==0 ) THEN
          CALL MPI_BSEND( n, 1, MPI_INTEGER, proc, &
                2001, MPI_COMM_WORLD, ierr )
          CYCLE
        END IF
      
        ! Construct temporary mesh structure for the received points:
        ! -----------------------------------------------------------
        Nmesh => AllocateMesh()
        Nmesh % Nodes % x => ProcRecv(i) % nodes_x
        Nmesh % Nodes % y => ProcRecv(i) % nodes_y
        Nmesh % Nodes % z => ProcRecv(i) % nodes_z
        Nmesh % NumberOfNodes = n

        ALLOCATE(nperm(n))
        DO j=1,n
          nPerm(j)=j
        END DO

        Var => OldVariables
        nvars = 0
        DO WHILE(ASSOCIATED(Var))
          IF ( Var % DOFs==1 .AND. ASSOCIATED(Var % Perm).AND..NOT.Var % Secondary ) THEN
            ALLOCATE(store(n)); store=0
            nvars = nvars+1
            CALL VariableAdd(nMesh % Variables,nMesh,Var % Solver, &
                     Var % Name,1,store,nperm )
            IF ( ASSOCIATED(Var % PrevValues) ) THEN
              j = SIZE(Var % PrevValues,2)
              nvars = nvars+j
              Nvar => VariableGet( Nmesh % Variables,Var % Name,ThisOnly=.TRUE.)
              ALLOCATE(Nvar % PrevValues(n,j))
            END IF
          END IF
          Var => Var % Next
        END DO

        ! try interpolating values for the points:
        ! ----------------------------------------
        ALLOCATE( FoundNodes(n) ); FoundNodes=.FALSE.
        CALL InterpolateMeshToMeshQ( OldMesh, nMesh, OldVariables, &
           nMesh % Variables, UseQuadrantTree, MaskName=MaskName, FoundNodes=FoundNodes )

        nfound = COUNT(FoundNodes)

        CALL MPI_BSEND( nfound, 1, MPI_INTEGER, proc, &
                2001, MPI_COMM_WORLD, ierr )

        ! send interpolated values back to the owner:
        ! -------------------------------------------
        IF ( nfound>0 ) THEN
          ALLOCATE(vstore(nfound,nvars), vperm(nfound)); vstore=0
          k = 0
          DO j=1,n
            IF ( .NOT.FoundNodes(j)) CYCLE   
            k = k + 1
            vperm(k) = j
            Var => OldVariables
            nvars = 0
            DO WHILE(ASSOCIATED(Var))
              IF ( Var % DOFs==1  .AND. ASSOCIATED(Var % Perm).AND..NOT.Var % Secondary) THEN
                Nvar => VariableGet( Nmesh % Variables,Var % Name,ThisOnly=.TRUE.)
                nvars=nvars+1
                vstore(k,nvars)=Nvar % Values(j)
                IF ( ASSOCIATED(Var % PrevValues) ) THEN
                  DO l=1,SIZE(Var % PrevValues,2)
                    nvars = nvars+1
                    vstore(k,nvars)=Nvar % PrevValues(j,l)
                  END DO
                END IF
              END IF
              Var => Var % Next
            END DO
          END DO

          CALL MPI_BSEND( vperm, nfound, MPI_INTEGER, proc, &
                2002, MPI_COMM_WORLD, status, ierr )

          DO j=1,nvars
            CALL MPI_BSEND( vstore(:,j), nfound,MPI_DOUBLE_PRECISION, proc, &
                       2002+j, MPI_COMM_WORLD,ierr )
          END DO

          DEALLOCATE(vstore, vperm)
        END IF

        DEALLOCATE(ProcRecv(i) % Nodes_x, ProcRecv(i) % Nodes_y,&
                   ProcRecv(i) % Nodes_z )

        Var => Nmesh % Variables
        DO WHILE(ASSOCIATED(Var))
          DEALLOCATE(Var % Values)
          IF (ASSOCIATED(Var % PrevValues)) DEALLOCATE(Var % PrevValues)
          nVar => Var
          Var => Var % Next
          DEALLOCATE(nVar)
        END DO
        DEALLOCATE(nperm,foundnodes, Nmesh)
      END DO
      DEALLOCATE(ProcRecv)

     ! Receive interpolated values:
     ! ----------------------------
      DO i=1,COUNT(ParEnv % Active)-1

        ! recv count:
        ! -----------
        CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
              2001, MPI_COMM_WORLD, status, ierr )

        proc = status(MPI_SOURCE)
        IF ( n<=0 ) THEN
          IF ( ALLOCATED(ProcSend) ) THEN
            IF ( ALLOCATED(ProcSend(proc+1) % Perm)) &
                       DEALLOCATE(ProcSend(proc+1) % Perm)
          END IF
          CYCLE
        END IF

        ALLOCATE(astore(n),vperm(n))

        ! recv permutation (where in the original array the
        ! points the partition found are):
        ! --------------------------------------------------
        CALL MPI_RECV( vperm, n, MPI_INTEGER, proc, &
              2002, MPI_COMM_WORLD, status, ierr )

        ! recv values and store:
        ! ----------------------
        Var => OldVariables
        nvars=0
        DO WHILE(ASSOCIATED(Var))
          IF ( Var % DOFs==1 .AND. ASSOCIATED(Var % Perm) .AND..NOT.Var % Secondary ) THEN

            nvars=nvars+1
            CALL MPI_RECV( astore, n, MPI_DOUBLE_PRECISION, proc, &
                2002+nvars, MPI_COMM_WORLD, status, ierr )

            Nvar => VariableGet( NewMesh % Variables,Var % Name,ThisOnly=.TRUE.)

            IF ( ASSOCIATED(Nvar) ) THEN
              DO j=1,n
                k=perm(ProcSend(proc+1) % Perm(vperm(j)))
                IF ( Nvar % perm(k)>0 ) &
                  Nvar % Values(Nvar % Perm(k)) = astore(j)
              END DO
            END IF

            IF ( ASSOCIATED(Var % PrevValues) ) THEN
              DO l=1,SIZE(Var % PrevValues,2)
                nvars=nvars+1
                CALL MPI_RECV( astore, n, MPI_DOUBLE_PRECISION, proc, &
                    2002+nvars, MPI_COMM_WORLD, status, ierr )

                IF ( ASSOCIATED(Nvar) ) THEN
                  DO j=1,n
                    k=perm(ProcSend(proc+1) % Perm(vperm(j)))
                    IF ( Nvar % perm(k)>0 ) &
                      Nvar % PrevValues(Nvar % Perm(k),l) = astore(j)
                  END DO
                END IF
              END DO
            END IF
          END IF
          Var => Var % Next
        END DO
        DEALLOCATE(astore,vperm,ProcSend(proc+1) % perm)
      END DO
      IF ( ALLOCATED(Perm) ) DEALLOCATE(Perm,ProcSend)

      CALL MPI_BARRIER(ParEnv % ActiveComm,ierr)

CONTAINS

!------------------------------------------------------------------------------
   FUNCTION AllocateMesh() RESULT(Mesh)
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
     INTEGER :: istat

     ALLOCATE( Mesh, STAT=istat )
     IF ( istat /= 0 ) &
        CALL Fatal( 'AllocateMesh', 'Unable to allocate a few bytes of memory?' )

!    Nothing computed on this mesh yet!
!    ----------------------------------
     Mesh % SavesDone    = 0
     Mesh % OutputActive = .FALSE.

     Mesh % AdaptiveDepth = 0
     Mesh % Changed   = .FALSE. !  TODO: Change this sometime

     Mesh % Stabilize = .FALSE.

     Mesh % Variables => NULL()
     Mesh % Parent => NULL()
     Mesh % Child => NULL()
     Mesh % Next => NULL()
     Mesh % RootQuadrant => NULL()
     Mesh % Elements => NULL()
     Mesh % Edges => NULL()
     Mesh % Faces => NULL()
     Mesh % Projector => NULL()
     Mesh % NumberOfEdges = 0
     Mesh % NumberOfFaces = 0
     Mesh % NumberOfNodes = 0
     Mesh % NumberOfBulkElements = 0
     Mesh % NumberOfBoundaryElements = 0

     Mesh % MaxFaceDOFs = 0
     Mesh % MaxEdgeDOFs = 0
     Mesh % MaxBDOFs = 0
     Mesh % MaxElementDOFs  = 0
     Mesh % MaxElementNodes = 0

     Mesh % ViewFactors => NULL()

     ALLOCATE( Mesh % Nodes, STAT=istat )
     IF ( istat /= 0 ) &
        CALL Fatal( 'AllocateMesh', 'Unable to allocate a few bytes of memory?' )
     NULLIFY( Mesh % Nodes % x )
     NULLIFY( Mesh % Nodes % y )
     NULLIFY( Mesh % Nodes % z )
     Mesh % Nodes % NumberOfNodes = 0

     Mesh % ParallelInfo % NumberOfIfDOFs =  0
     NULLIFY( Mesh % ParallelInfo % GlobalDOFs )
     NULLIFY( Mesh % ParallelInfo % INTERFACE )
     NULLIFY( Mesh % ParallelInfo % NeighbourList )

  END FUNCTION AllocateMesh
!-------------------------------------------------------------------------------
     END SUBROUTINE InterpolateMeshToMeshR
!-------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Interpolates values of all variables from a mesh associated with
!>    the old model to the mesh of the new model.
!------------------------------------------------------------------------------
     SUBROUTINE InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, &
            NewVariables, UseQuadrantTree, Projector, MaskName, FoundNodes )
!------------------------------------------------------------------------------
       USE Interpolation
       USE CRSMatrix
       USE CoordinateSystems
!-------------------------------------------------------------------------------
       TYPE(Mesh_t), TARGET  :: OldMesh   !< Old mesh structure
       TYPE(Mesh_t), TARGET  :: NewMesh   !< New mesh structure
       TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables  !< Old model variable structure
       TYPE(Variable_t), POINTER, OPTIONAL :: NewVariables  !< New model variable structure. NewVariables defines the variables to be interpolated
       LOGICAL, OPTIONAL :: UseQuadrantTree  !< If true use the RootQuadrant of the old mesh in interpolation.
       TYPE(Projector_t), POINTER, OPTIONAL :: Projector  !< Use projector between meshes for interpolation, if available
       CHARACTER(LEN=*),OPTIONAL :: MaskName  !< Mask the old variable set by the given MaskName when trying to define the interpolation.
       LOGICAL, OPTIONAL :: FoundNodes(:)     !< List of nodes where the interpolation was a success
!------------------------------------------------------------------------------
       INTEGER :: dim
       TYPE(Nodes_t) :: ElementNodes
       INTEGER :: nBulk, i, j, k, l, n, bf_id, QTreeFails, TotFails
       REAL(KIND=dp), DIMENSION(3) :: Point
       INTEGER, POINTER :: NodeIndexes(:)
       REAL(KIND=dp), DIMENSION(3) :: LocalCoordinates
       TYPE(Variable_t), POINTER :: OldSol, NewSol, Var
       INTEGER, POINTER :: OldPerm(:)
       REAL(KIND=dp), POINTER :: OldValue(:), NewValue(:), ElementValues(:)
       TYPE(Quadrant_t), POINTER :: LeafQuadrant
       TYPE(Element_t),POINTER :: Element
       
       REAL(KIND=dp), ALLOCATABLE :: Basis(:),Vals(:),dBasisdx(:,:)
       REAL(KIND=dp) :: BoundingBox(6), detJ, u,v,w,s,val,rowsum
       
       LOGICAL :: UseQTree, TryQTree, Stat, UseProjector
       TYPE(Quadrant_t), POINTER :: RootQuadrant
       
       INTEGER, POINTER   :: Rows(:), Cols(:)
       
       TYPE Epntr_t
         TYPE(Element_t), POINTER :: Element
       END TYPE Epntr_t
       
       TYPE(Epntr_t), ALLOCATABLE :: ElemPtrs(:)
       
       INTEGER, ALLOCATABLE:: RInd(:)
       LOGICAL :: Found, EpsAbsGiven,EpsRelGiven, MaskExists, ProjectorAllocated
       INTEGER :: eps_tries, nrow
       REAL(KIND=dp) :: eps1 = 0.1, eps2, eps_global, eps_local, eps_basis,eps_numeric
       REAL(KIND=dp), POINTER :: Values(:), LocalU(:), LocalV(:), LocalW(:)
       !$OMP THREADPRIVATE(eps1)
!------------------------------------------------------------------------------

!
!      If projector argument given, search for existing
!      projector matrix, or generate new projector, if
!      not already there:
!      ------------------------------------------------
       IF ( PRESENT(Projector) ) THEN
         Projector => NewMesh % Projector
         
         DO WHILE( ASSOCIATED( Projector ) )
           IF ( ASSOCIATED(Projector % Mesh, OldMesh) ) THEN
             IF ( PRESENT(OldVariables) ) CALL ApplyProjector()
             RETURN
           END IF
           Projector => Projector % Next
         END DO

         n = NewMesh % NumberOfNodes
         ALLOCATE( LocalU(n), LocalV(n), LocalW(n), ElemPtrs(n) )
         DO i=1,n
           NULLIFY( ElemPtrs(i) % Element )
         END DO
       END IF
!
!      Check if using the spatial division hierarchy for the search:
!      -------------------------------------------------------------

       RootQuadrant => OldMesh % RootQuadrant
       dim = CoordinateSystemDimension()
       
       IF ( .NOT. PRESENT( UseQuadrantTree ) ) THEN
         UseQTree = .TRUE.
       ELSE
         UseQTree = UseQuadrantTree
       ENDIF
    
       IF ( UseQTree ) THEN
         IF ( .NOT.ASSOCIATED( RootQuadrant ) ) THEN
           BoundingBox(1) = MINVAL(OldMesh % Nodes % x)
           BoundingBox(2) = MINVAL(OldMesh % Nodes % y)
           BoundingBox(3) = MINVAL(OldMesh % Nodes % z)
           BoundingBox(4) = MAXVAL(OldMesh % Nodes % x)
           BoundingBox(5) = MAXVAL(OldMesh % Nodes % y)
           BoundingBox(6) = MAXVAL(OldMesh % Nodes % z)
           
           eps2 = 0.1_dp * MAXVAL(BoundingBox(4:6)-BoundingBox(1:3))
           BoundingBox(1:3) = BoundingBox(1:3) - eps2
           BoundingBox(4:6) = BoundingBox(4:6) + eps2
           
           CALL BuildQuadrantTree( OldMesh,BoundingBox,OldMesh % RootQuadrant)
           RootQuadrant => OldMesh % RootQuadrant
         END IF
       END IF
       
! Use mask or not
!---------------------------------------
       MaskExists = PRESENT( MaskName )

!------------------------------------------------------------------------------

       n = OldMesh % MaxElementNodes
       ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), &
           ElementNodes % z(n), ElementValues(n) )
       
       eps_global = ListGetConstReal( CurrentModel % Simulation,  &
           'Interpolation Global Epsilon', Stat)
       IF(.NOT. Stat) eps_global = 2.0d-10
       
       eps_local = ListGetConstReal( CurrentModel % Simulation,  &
           'Interpolation Local Epsilon', Stat )
       IF(.NOT. Stat) eps_local = 1.0d-10

       eps_tries = ListGetInteger( CurrentModel % Simulation,  &
           'Interpolation Max Iterations', Stat )
       IF(.NOT. Stat) eps_tries = 12

       eps_numeric = ListGetConstReal( CurrentModel % Simulation, &
           'Interpolation Numeric Epsilon', Stat)
       IF(.NOT. Stat) eps_numeric = 1.0e-10

       QTreeFails = 0
       TotFails = 0

!------------------------------------------------------------------------------
! Loop over all nodes in the new mesh
!------------------------------------------------------------------------------
       DO i=1,NewMesh % NumberOfNodes
!------------------------------------------------------------------------------
         Point(1) = NewMesh % Nodes % x(i)
         Point(2) = NewMesh % Nodes % y(i)
         Point(3) = NewMesh % Nodes % z(i)

!------------------------------------------------------------------------------
! Find in which old mesh bulk element the point belongs to
!------------------------------------------------------------------------------
         Found = .FALSE.
         TryQTree = ASSOCIATED(RootQuadrant) .AND. UseQTree 

         IF( TryQTree ) THEN
!------------------------------------------------------------------------------
! Find the last existing quadrant that the point belongs to
!------------------------------------------------------------------------------
           Element => NULL()
           CALL FindLeafElements(Point, dim, RootQuadrant, LeafQuadrant)
           
           IF ( ASSOCIATED(LeafQuadrant) ) THEN
             ! Go through the bulk elements in the last ChildQuadrant
             ! only.  Try to find matching element with progressively
             ! sloppier tests. Allow at most 100 % of slack:
             ! -------------------------------------------------------
             Eps1 = eps_global
             Eps2 = eps_local
             
             DO j=1,eps_tries
               DO k=1, LeafQuadrant % NElemsInQuadrant
                 Element => OldMesh % Elements(LeafQuadrant % Elements(k))
                 
                 IF( MaskExists ) THEN
                   bf_id = ListGetInteger( CurrentModel % Bodies(Element % BodyId) % Values, &
                       'Body Force', Found )
                   IF( .NOT. Found ) CYCLE
                   IF(.NOT. ListCheckPresent( &
                       CurrentModel % BodyForces(bf_id) % Values,MaskName) ) CYCLE
                 END IF
                  
                 NodeIndexes => Element % NodeIndexes
                 n = Element % TYPE % NumberOfNodes
                 
                 ElementNodes % x(1:n) = OldMesh % Nodes % x(NodeIndexes)
                 ElementNodes % y(1:n) = OldMesh % Nodes % y(NodeIndexes)
                 ElementNodes % z(1:n) = OldMesh % Nodes % z(NodeIndexes)
                 
                 Found = PointInElement( Element, ElementNodes, &
                     Point, LocalCoordinates, Eps1, Eps2, NumericEps=eps_numeric )
                 IF ( Found ) EXIT
               END DO
               IF ( Found ) EXIT  
               
               Eps1 = 10 * Eps1
               Eps2 = 10 * Eps2               
               IF( Eps1 > 1.0_dp ) EXIT
             END DO
           END IF
         END IF

         IF( .NOT. TryQTree .OR. &
             (.NOT. Found .AND. .NOT. PRESENT( FoundNodes) ) ) THEN
           !------------------------------------------------------------------------------
           ! Go through all old mesh bulk elements
           !------------------------------------------------------------------------------
           DO k=1,OldMesh % NumberOfBulkElements
             Element => OldMesh % Elements(k)
             
             n = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes
             
             ElementNodes % x(1:n) = OldMesh % Nodes % x(NodeIndexes)
             ElementNodes % y(1:n) = OldMesh % Nodes % y(NodeIndexes)
             ElementNodes % z(1:n) = OldMesh % Nodes % z(NodeIndexes)
             
             Found =  PointInElement( Element, ElementNodes, &
                 Point, LocalCoordinates  ) 
             IF( Found ) THEN
               IF( TryQTree ) QTreeFails = QtreeFails + 1
               EXIT
             END IF
           END DO
         END IF
         
         IF (.NOT.Found) THEN
           Element => NULL()
           IF (.NOT.PRESENT(FoundNodes) ) THEN
             WRITE( Message,'(A,I0,A)' ) 'Point ',i,' was not found in any of the elements!'
             CALL Info( 'InterpolateMeshToMesh', Message, Level=20 )
             TotFails = TotFails + 1
           END IF
           CYCLE
         END IF
         IF ( PRESENT(FoundNodes) ) FoundNodes(i) = .TRUE.

!------------------------------------------------------------------------------
!
!         Found Element in OldModel:
!         ---------------------------------
          IF ( PRESENT(Projector) ) THEN
             ElemPtrs(i) % Element => Element
             LocalU(i) = LocalCoordinates(1)
             LocalV(i) = LocalCoordinates(2)
             LocalW(i) = LocalCoordinates(3)
          END IF

          IF ( .NOT.PRESENT(OldVariables) .OR. PRESENT(Projector) ) CYCLE
!------------------------------------------------------------------------------
!
!         Go through all variables to be interpolated:
!         --------------------------------------------
          Var => OldVariables
          DO WHILE( ASSOCIATED( Var ) )

             IF( SIZE( Var % Values ) == Var % DOFs ) THEN
               Var => Var % Next
               CYCLE
             END IF          

             IF( Var % Secondary ) THEN
               Var => Var % Next
               CYCLE
             END IF

             IF ( Var % DOFs == 1 .AND. &
                 Var % Name(1:10) /= 'coordinate') THEN

!------------------------------------------------------------------------------
!
!               Interpolate variable at Point in Element:
!               ------------------------------------------------

                NewSol => VariableGet( NewMesh % Variables, Var % Name, .TRUE. )
                IF ( .NOT. ASSOCIATED( NewSol ) ) THEN
                   Var => Var % Next
                   CYCLE
                END IF
                OldSol => VariableGet( OldMesh % Variables, Var % Name, .TRUE. )


                ! Check that the node was found in the old mesh:
                ! ----------------------------------------------
                IF ( ASSOCIATED (Element) ) THEN
!------------------------------------------------------------------------------
!
!                  Check for rounding errors:
!                  --------------------------
                   NodeIndexes => Element % NodeIndexes
                   IF ( ALL(OldSol % Perm(NodeIndexes)>0) ) THEN
                     IF ( NewSol % Perm(i) /= 0 ) THEN
                       ElementValues(1:n) = & 
                              OldSol % Values(OldSol % Perm(NodeIndexes))
                       NewSol % Values(NewSol % Perm(i)) = InterpolateInElement( &
                            Element, ElementValues, LocalCoordinates(1), &
                                LocalCoordinates(2), LocalCoordinates(3) )

                       IF ( ASSOCIATED( OldSol % PrevValues ) ) THEN
                         DO j=1,SIZE(OldSol % PrevValues,2)
                           ElementValues(1:n) = &
                               OldSol % PrevValues(OldSol % Perm(NodeIndexes),j)
                           NewSol % PrevValues(NewSol % Perm(i),j) = &
                             InterpolateInElement( Element, ElementValues, &
                               LocalCoordinates(1), &
                                  LocalCoordinates(2), LocalCoordinates(3) )
                         END DO
                       END IF
                     END IF
                   END IF
                ELSE
                   IF ( NewSol % Perm(i)/=0 ) NewValue(NewSol % Perm(i))=0.0_dp
                END IF

!------------------------------------------------------------------------------
             END IF
             Var => Var % Next
           END DO
!------------------------------------------------------------------------------
         END DO

         IF( .NOT. PRESENT( FoundNodes ) ) THEN
           IF( QtreeFails > 0 ) THEN
             WRITE( Message,'(A,I0)' ) 'Number of points not found in quadtree: ',QtreeFails
             CALL Info( 'InterpolateMeshToMesh', Message )
             IF( TotFails == 0 ) THEN
               CALL Info( 'InterpolateMeshToMesh','All nodes still found by N^2 dummy search!' )               
             END IF
           END IF
           IF( TotFails == 0 ) THEN
             CALL Info('InterpolateMeshToMesh','Found all nodes in the target mesh',Level=6)
           ELSE
             WRITE( Message,'(A,I0,A,I0,A)') 'Points not found: ',TotFails,' (found ',&
                 NewMesh % NumberOfNodes - TotFails,')'
             CALL Warn( 'InterpolateMeshToMesh', Message )
           END IF
         END IF

!------------------------------------------------------------------------------
!
!      Construct mesh projector, if requested. Next time around
!      will use the existing projector to interpolate values:
!      ---------------------------------------------------------
       IF ( PRESENT(Projector) ) THEN

          n = NewMesh % NumberOfNodes
          ALLOCATE( Basis(OldMesh % MaxElementNodes), Vals(OldMesh % MaxElementNodes) )

          ! The critical value of basis function that is accepted to the 
          ! projector. Note that the sum of weights is one, so this
          ! we know the scale for this one. 
          eps_basis = ListGetConstReal( CurrentModel % Simulation,  &
              'Interpolation Basis Epsilon', Stat )
          IF(.NOT. Stat) eps_basis = 0.0d-12

          ALLOCATE( Rows(n+1) )
          Rows(1) = 1
          ProjectorAllocated = .FALSE.

100       nrow = 1

          DO i=1,n

            Element => ElemPtrs(i) % Element
            Found = ASSOCIATED( Element ) 
            
            IF( .NOT. Found ) THEN
             ! It seems unnecessary to make a matrix entry in case no target element is found!
              IF(.FALSE.) THEN
                IF( ProjectorAllocated ) THEN
                  Cols(nrow) = 1
                  Values(nrow) = 0._dp
                END IF
                nrow = nrow + 1
              END IF
            ELSE
              k = Element % TYPE % NumberOfNodes

              NodeIndexes => Element % NodeIndexes
              
              u = LocalU(i)
              v = LocalV(i)
              w = LocalW(i)
 
              Basis(1:k) = 0.0d0
              DO j=1,k
                Basis(j) = 1.0d0
                Vals(j) = InterpolateInElement(Element,Basis,u,v,w)
                Basis(j) = 0.0d0
              END DO
              
              rowsum = 0.0_dp
              DO j=1,k
                IF( ABS( vals(j) ) < eps_basis ) CYCLE
                rowsum = rowsum + vals(j)
                IF( .NOT. ProjectorAllocated ) nrow = nrow + 1
              END DO
              

              IF( ProjectorAllocated ) THEN
                DO j=1,k
                  IF( ABS( vals(j) ) < eps_basis ) CYCLE

                  Rind(NodeIndexes(j)) = Rind(NodeIndexes(j)) + 1
                  Cols(nrow) = NodeIndexes(j)

                  ! Always normalize the weights to one!
                  Values(nrow) = vals(j) / rowsum
                  nrow = nrow + 1                  
                END DO
              END IF                
            END IF

            Rows(i+1) = nrow
         END DO

         IF( .NOT. ProjectorAllocated ) THEN
            ALLOCATE( Cols(Rows(n+1)-1), Values(Rows(n+1)-1) )
            Cols   = 0
            Values = 0
            
            ALLOCATE( Projector )
            Projector % Matrix => AllocateMatrix()
            Projector % Matrix % NumberOfRows = n
            Projector % Matrix % Rows   => Rows
            Projector % Matrix % Cols   => Cols 
            Projector % Matrix % Values => Values
            
            Projector % Next => NewMesh % Projector
            NewMesh % Projector => Projector
            NewMesh % Projector % Mesh => OldMesh
            
            ALLOCATE( RInd(OldMesh % NumberOfNodes) )
            RInd = 0

            ProjectorAllocated = .TRUE.

            GOTO 100
          END IF


          DEALLOCATE( Basis, Vals, ElemPtrs, LocalU, LocalV, LocalW )

!         Store also the transpose of the projector:
!         ------------------------------------------ 
          Projector % TMatrix => NULL()
          IF ( Found ) THEN
            n = OldMesh % NumberOfNodes
            ! Needed for some matrices
            n = MAX( n, MAXVAL( Projector % Matrix % Cols ) )

            ALLOCATE( Rows(n+1) )
            Rows(1) = 1
            DO i=2,n+1
               Rows(i) = Rows(i-1)+RInd(i-1)
            END DO

            ALLOCATE( Cols(Rows(n+1)-1), Values(Rows(n+1)-1) )
            Projector % TMatrix => AllocateMatrix()
            Projector % TMatrix % NumberOfRows = n
            Projector % TMatrix % Rows   => Rows
            Projector % TMatrix % Cols   => Cols 
            Projector % TMatrix % Values => Values

            RInd = 0
            DO i=1,Projector % Matrix % NumberOfRows
              DO j=Projector % Matrix % Rows(i), Projector % Matrix % Rows(i+1)-1
                 k = Projector % Matrix % Cols(j)
                 l = Rows(k)+RInd(k)
                 RInd(k) = RInd(k)+1
                 Cols(l) = i
                 Values(l) = Projector % Matrix % Values(j)
              END DO
            END DO
          END IF

          DEALLOCATE(RInd)

          IF ( PRESENT(OldVariables) ) CALL ApplyProjector
       END IF

       DEALLOCATE( ElementNodes % x, ElementNodes % y, &
                   ElementNodes % z, ElementValues )

CONTAINS

!------------------------------------------------------------------------------
     SUBROUTINE ApplyProjector
!------------------------------------------------------------------------------
        INTEGER :: i
        TYPE(Variable_t), POINTER :: Var
!------------------------------------------------------------------------------
        Var => OldVariables
        DO WHILE( ASSOCIATED(Var) )
           IF( SIZE( Var % Values ) == Var % DOFs ) THEN   
             Var => Var % Next
             CYCLE
           END IF 

           IF( Var % Secondary ) THEN
             Var => Var % Next
             CYCLE
           END IF 

           IF ( Var % DOFs == 1 .AND. &
             Var % Name(1:10) /= 'coordinate') THEN

              OldSol => VariableGet( OldMesh % Variables, Var % Name, .TRUE. )
              NewSol => VariableGet( NewMesh % Variables, Var % Name, .TRUE. )
              IF ( .NOT. (ASSOCIATED (NewSol) ) ) THEN
                 Var => Var % Next
                 CYCLE
              END IF

              CALL CRS_ApplyProjector( Projector % Matrix, &
                   OldSol % Values, OldSol % Perm,         &
                   NewSol % Values, NewSol % Perm )

              IF ( ASSOCIATED( OldSol % PrevValues ) ) THEN
                 DO i=1,SIZE(OldSol % PrevValues,2)
                    CALL CRS_ApplyProjector( Projector % Matrix,  &
                         OldSol % PrevValues(:,i), OldSol % Perm, &
                         NewSol % PrevValues(:,i), NewSol % Perm )
                 END DO
              END IF
           END IF
           Var => Var % Next
        END DO
!------------------------------------------------------------------------------
     END SUBROUTINE ApplyProjector
!------------------------------------------------------------------------------
  END SUBROUTINE InterpolateMeshToMeshQ
!------------------------------------------------------------------------------

