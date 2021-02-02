      SUBROUTINE hRefinementSolver( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
      USE DefUtils
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t), TARGET :: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
      TYPE(Solver_t), POINTER :: PSolver
      Type(Mesh_t),POINTER :: NewMesh,PrevMesh
      Type(Variable_t),POINTER :: ESVar,Var,NewVar
      TYPE(ValueList_t), POINTER :: SolverParams

      INTEGER :: i,j,k,t

      LOGICAL :: Found,UseProjector
      LOGICAL :: DoIT=.True.

      INTEGER :: MaxRefinementDepth

      CHARACTER(LEN=MAX_NAME_LEN) :: method
      CHARACTER(LEN=MAX_NAME_LEN),parameter :: &
                       SolverName='MeshAdaptationSolver'
!------------------------------------------------------------------------------

     ! Check for keyword "Do hRefinement"
     !  If Found => another solver (or function) decide when to refine
      DoIT=ListGetLogical(CurrentModel % Simulation,'Do hRefinement',Found)
     ! keyword do not exist => do refinement
      If (.NOT.Found) DoIt=.True. 
     ! if no need to refine => return
      If (.NOT.DoIt) RETURN ! 
     ! refinement will be done set 'Do hRefinement' and wait for new
     ! instruction
      IF (Found) CALL ListAddLogical(CurrentModel % Simulation,&
                                'Do hRefinement',.FALSE.)


      CALL ResetTimer(trim(SolverName))

      !! A variable that correspond to the required Element Size
      ESVar => VariableGet(Solver%Mesh%Variables,'ElementSize',.FALSE.)
      IF (.NOT.ASSOCIATED(ESVar)) Then
         CALL FATAL(trim(SolverName),'Variable <ElementSize> not found')
      END IF

      !!! Create Edge Tables if not done
      IF (Solver % Mesh % NumberOfEdges==0) then
         call SetMeshEdgeFaceDOFs(Solver%Mesh)
      End if

      ! Get some global parameters
      SolverParams => GetSolverParams()

      method=GetString(SolverParams,'h-refinement method',Found)
      IF (.NOT.Found) CALL FATAL(trim(SolverName),&
                      'Keyword <h-refinement method> not found')

      !! Max number of refinement
      MaxRefinementDepth = &
         GetInteger(SolverParams,'Maximum Refinement Depth',Found)
          If (.NOT.Found) CALL FATAL(trim(SolverName),&
                'Solver parameter <Maximum Refinement Depth> not found')

      ! Refine Solver%Mesh
      SELECT CASE(method)
          CASE ('lepp')
             !! SplitMesh
             CALL INFO(trim(SolverName),&
                  'Refine mesh using Lepp-bisection algorithm')
             NewMesh => LeppAdaptivity(Solver%Mesh)
             !!!
          CASE ('rgb')
             CALL INFO(trim(SolverName),&
                  'Refine mesh using RGB algorithm')
             NewMesh => SplitMesh(Solver%Mesh)

          CASE DEFAULT
          write(Message,'(A,A,A)') 'refinement method ',&
                           trim(method),' not implemented'
          CALL FATAL(trim(SolverName),trim(Message))
      END SELECT


     !!! Get Previous child mesh to do variable interpolations
      IF (ASSOCIATED(Solver%Mesh%Child)) then 
        PrevMesh => Solver%Mesh%Child
      Else
        PrevMesh => Solver%Mesh
      Endif

     !!!!!!! Temporary trick to interoplate vectors in the master mesh; 
     !otherwise seg.fault. when resultouputsolver tries to
     !interpolate compenents....????
      Var => Solver % Mesh % Variables
      DO WHILE( ASSOCIATED( Var ) )
         IF ( Var % DOFs > 1 ) THEN
            NewVar => VariableGet( Solver% Mesh % Variables,Var % Name,.FALSE. )
            k = SIZE( NewVar % Values )
            IF ( ASSOCIATED( NewVar % Perm ) ) THEN
              k = COUNT( NewVar % Perm > 0 )
            END IF
            IF ( GetVarName( NewVar ) == 'flow solution' ) THEN
               NewVar % Norm = 0.0d0
               DO i=1,NewMesh % NumberOfNodes
                  DO j=1,NewVar % DOFs-1
                     NewVar % Norm = NewVar % Norm + &
                          NewVar % Values( NewVar % DOFs*(i-1)+j )**2
                  END DO
               END DO
               NewVar % Norm = SQRT( NewVar % Norm / k )
            ELSE
               NewVar % Norm = SQRT( SUM(NewVar % Values**2) / k )
            END IF
         END IF
         Var => Var % Next
      END DO
     !!!!

     !! WARNING !!!
     !!!!!!!!!!!!!!
     !If Use Mesh Projector is true (default) the same mesh projector
     !will be used for interpolation which is not good if mesh has changed
     !!!!
      UseProjector=.FALSE.
      UseProjector=&
         ListGetLogical(CurrentModel % Simulation,'Use Mesh Projector',Found)
      IF ((.NOT.Found).OR.(UseProjector)) &
      CALL FATAL('TestSolver',&
            'Set Use Mesh Projector to False in Simulation')

     !!! The new mesh becomes the computation mesh
      CALL MakeMeshGlobal(NewMesh,PrevMesh)  

     !   Add the new mesh to the global list of meshes:
     !   ----------------------------------------------
      NewMesh % Next   => Solver%Mesh
      Model % Meshes   => NewMesh
      NewMesh % Parent => Solver%Mesh
      NewMesh % Child => NULL()

     !   Update Solver Meshes to the new Mesh:
     !   ----------------------------------------------
      Do i=1,Model % NumberOfSolvers
       PSolver => Model % Solvers(i)
       IF (.NOT.ASSOCIATED(PSolver)) CYCLE
       if (PSolver % Variable % Name == Solver % Variable % name) cycle

       !PRINT *,PSolver % Variable % Name,ASSOCIATED(PSolver % Matrix)
       !!! if Matrix not associated no need to update solver Mesh
       IF (.NOT.ASSOCIATED(PSolver % Matrix))  CYCLE
       call MyUpdateSolverMesh( PSolver, NewMesh )

       write(Message,'(A,A)') '<Solver % Mesh> Updated for variable',&
                                        trim(PSolver % Variable % Name)
       CALL INFO(trim(SolverName),trim(Message),level=6)
      End do

      !! Release previous child
      IF (ASSOCIATED(Solver%Mesh%Child)) THEN
        !! BoundaryInfo is not deallocated by Remlease Mesh??
        Do i=1,PrevMesh%NumberOfBoundaryElements
           k=i+PrevMesh%NumberOfBulkElements
           if (ASSOCIATED(PrevMesh%Elements(k)%BoundaryInfo)) then
              deallocate(PrevMesh%Elements(k)%BoundaryInfo)
              PrevMesh%Elements(k)%BoundaryInfo=>NULL()
           endif
        End do
        CALL ReleaseMesh(PrevMesh)
        Deallocate(PrevMesh)
      ENDIF
      !! And set Child to New Mesh
      Solver%Mesh%Child => NewMesh
      NewMesh % Changed = .true.

      write(Message,'(A)') 'New mesh ready'
       CALL INFO(trim(SolverName),trim(Message),level=5)
      write(Message,'(A,I0)') '    NumberOfNodes: ',NewMesh%NumberOfNodes
       CALL INFO(trim(SolverName),trim(Message),level=5)
      write(Message,'(A,I0)') '    NumberOfBulkElements: ',NewMesh%NumberOfBulkElements
       CALL INFO(trim(SolverName),trim(Message),level=5)
      write(Message,'(A,I0)') '    NumberOfEdges: ',NewMesh%NumberOfEdges
       CALL INFO(trim(SolverName),trim(Message),level=5)
      write(Message,'(A,I0)') '    NumberOfBoundaryElements: ',NewMesh%NumberOfBoundaryElements
       CALL INFO(trim(SolverName),trim(Message),level=5)

      CALL CheckTimer(trim(SolverName),Delete=.TRUE.)

     CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LeppAdaptivity(RefMesh) RESULT(NewMesh)
      !
      TYPE(Mesh_t), POINTER :: NewMesh, RefMesh
      !
      Type(Element_t), POINTER :: Edge,Element,NeighbourElement,ParentElement
      Type(Mesh_t),POINTER :: tmpMesh,BoundaryMesh
      Type(Nodes_t),SAVE :: ElementNodes
      TYPE(ValueList_t), POINTER :: SolverParams

      REAL(KIND=dp) :: x1,x2,y1,y2
      REAL(KIND=dp) :: Dx,ESize,hk

      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: LeppTable
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: Parent ! Element Index of the parent in the initial mesh
      INTEGER,DIMENSION(:,:),ALLOCATABLE,SAVE :: Child ! Element Index of the childs
      INTEGER :: Node1,Node2
      INTEGER :: LongestEdge,NeighbourLongestEdge,LeppSize
      INTEGER :: NONodes,NOEdges
      INTEGER :: NOBulkElements,NOBoundaryElements
      INTEGER :: FacePoint
      INTEGER,SAVE :: MaxNumberOfNodes,MaxNumberOfElements,MaxNumberOfEdges
      INTEGER,SAVE :: MaxNumberOfBoundaryElements
      INTEGER,SAVE :: MaxNumberOfChilds
      INTEGER :: i,j,k,t
      INTEGER :: NodeIndexe,NofBulk
      INTEGER :: ParentIndex
      INTEGER :: Counter

      INTEGER,SAVE  :: MeshSizeIncrease

      LOGICAL :: Done
      LOGICAL :: Found
      LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: MarkedElement
      LOGICAL,SAVE :: AllocationDone=.False.

      CHARACTER(LEN=MAX_NAME_LEN),parameter :: &
                       SolverName='LeppAdaptationSolver'

      NULLIFY(NewMesh)

      !! Do some allocation
      if ((.NOT.AllocationDone).OR.RefMesh%Changed) then
          If (AllocationDone) then
             deallocate(ElementNodes%x,&
                        ElementNodes%y,&
                        ElementNodes%z,&
                        LeppTable,&
                        Parent,&
                        Child,&
                        MarkedElement)
          Endif

          !! Read some solver parameters
          SolverParams => GetSolverParams()

          !! an interger coeff. to allocate space for the construction of the new mesh
          MeshSizeIncrease=GetInteger(SolverParams,'Mesh Size Increase',Found)
          If (.NOT.Found) CALL FATAL(trim(SolverName),&
                      'Solver parameter <Mesh Size Increase> not found')

          MaxNumberOfElements=MeshSizeIncrease*RefMesh%NumberOfBulkElements
          MaxNumberOfBoundaryElements=MeshSizeIncrease*RefMesh%NumberOfBoundaryElements
          MaxNumberOfNodes=MeshSizeIncrease*RefMesh%NumberOfNodes
          MaxNumberOfEdges = 2*MeshSizeIncrease*RefMesh%NumberOfEdges
          !MaxNumberOfChilds = 2**MeshSizeIncrease
          MaxNumberOfChilds = 2**(MaxRefinementDepth+2)
          Allocate(ElementNodes%x(3),ElementNodes%y(3),ElementNodes%z(3))
          Allocate(LeppTable(RefMesh%NumberOfBulkElements))
          ALLOCATE(Parent(MaxNumberOfElements))
          ALLOCATE(Child(RefMesh % NumberOfBulkElements,MaxNumberOfChilds))
          ALLOCATE(MarkedElement(MaxNumberOfElements))


          AllocationDone=.True.
      endif


      NONodes=RefMesh%NumberOfNodes
      NOEdges=RefMesh % NumberOfEdges
      NOBulkElements = RefMesh % NumberOfBulkElements
      NOBoundaryElements = RefMesh % NumberOfBoundaryElements

      !Parent and child arrays
      Parent=-1
      Child=-1

CALL ResetTimer('MeshCopy')

      !temporary Mesh: Do the allocations and copy the initial mesh
      tmpMesh => AllocateMesh()
      CALL AllocateVector( tmpMesh % Elements, MaxNumberOfElements )
      CALL AllocateVector( tmpMesh % Edges , MaxNumberOfEdges)
      CALL AllocateVector( tmpMesh % Nodes % x, MaxNumberOfNodes )
      CALL AllocateVector( tmpMesh % Nodes % y, MaxNumberOfNodes )
      CALL AllocateVector( tmpMesh % Nodes % z, MaxNumberOfNodes )

      tmpMesh % NumberOfBulkElements = RefMesh % NumberOfBulkElements
      tmpMesh % NumberOfBoundaryElements = 0
      tmpMesh % NumberOfNodes = RefMesh % NumberOfNodes
      tmpMesh % NumberOfEdges = RefMesh % NumberOfEdges

      tmpMesh % Nodes % x(1:RefMesh%NumberOfNodes)=RefMesh%Nodes % x(1:RefMesh%NumberOfNodes)
      tmpMesh % Nodes % y(1:RefMesh%NumberOfNodes)=RefMesh%Nodes % y(1:RefMesh%NumberOfNodes)
      tmpMesh % Nodes % z(1:RefMesh%NumberOfNodes)=RefMesh%Nodes % z(1:RefMesh%NumberOfNodes)

      ! Copy Bulk Elements
      Do i=1,RefMesh % NumberOfBulkElements
         tmpMesh % Elements(i)=RefMesh%Elements(i)
         CALL AllocateVector(tmpMesh %Elements(i) % NodeIndexes, 3)
         CALL AllocateVector(tmpMesh %Elements(i) % EdgeIndexes, 3)
         tmpMesh %Elements(i) % NodeIndexes(1:3)=RefMesh % Elements(i) % NodeIndexes(1:3)
         tmpMesh %Elements(i) % EdgeIndexes(1:3)=RefMesh % Elements(i) % EdgeIndexes(1:3)
         Parent(i) = RefMesh % Elements(i)%ElementIndex
         Child(i,1) = RefMesh % Elements(i)%ElementIndex
      End do

      ! Copy Edges NodeIndexes and BoundaryInfo
      Do i=1,RefMesh % NumberOfEdges
         tmpMesh % Edges(i)=RefMesh%Edges(i)
         CALL AllocateVector(tmpMesh % Edges(i) % NodeIndexes,2)
         ALLOCATE(tmpMesh % Edges(i) % BoundaryInfo)
         tmpMesh % Edges(i) % NodeIndexes(1:2)=RefMesh % Edges(i) % NodeIndexes(1:2)
         tmpMesh % Edges(i) % BoundaryInfo = RefMesh % Edges(i) % BoundaryInfo
      End do

      !Copy Boundary Element in another structure
      BoundaryMesh => AllocateMesh()
      CALL AllocateVector( BoundaryMesh % Elements,MaxNumberOfBoundaryElements )

      BoundaryMesh % NumberOfBulkElements = 0
      BoundaryMesh % NumberOfBoundaryElements = RefMesh%NumberOfBoundaryElements
      BoundaryMesh % NumberOfNodes = 0

      Do i=1 ,RefMesh % NumberOfBoundaryElements
          k=RefMesh % NumberOfBulkElements+i
          BoundaryMesh % Elements(i)=RefMesh%Elements(k)
          CALL AllocateVector(BoundaryMesh % Elements(i) % NodeIndexes, 2)
          CALL AllocateVector(BoundaryMesh % Elements(i) % EdgeIndexes, 1) 
          ALLOCATE(BoundaryMesh%Elements(i)%BoundaryInfo)
          BoundaryMesh % Elements(i) % NodeIndexes(1:2)=&
              RefMesh % Elements(k) % NodeIndexes(1:2)

          ! Edge Indexes not created in the initial mesh!
          ! need to find them
          IF (.NOT.ASSOCIATED(RefMesh % Elements(k)%BoundaryInfo%Left)) &
             CALL FATAL(trim(SolverName),&
                   'BC Element BoundaryInfo%Left NOT ASSOCIATED??')
          IF (ASSOCIATED(RefMesh % Elements(k)%BoundaryInfo%Right)) &
             CALL FATAL(trim(SolverName),&
                   'BC Element BoundaryInfo%Right ASSOCIATED??')

          !!! Keep the parent in the initial mesh  
          BoundaryMesh % Elements(i) % BoundaryInfo = RefMesh % Elements(k) % BoundaryInfo
          BoundaryMesh % Elements(i) %BoundaryInfo % Left => &
              RefMesh % Elements(k)%BoundaryInfo%Left
          BoundaryMesh % Elements(i) %BoundaryInfo % Constraint = &
              RefMesh % Elements(k)%BoundaryInfo%Constraint
          NULLIFY( Edge )
          DO j=1,3
             Edge => RefMesh%Edges(RefMesh%Elements(k)%BoundaryInfo%Left%EdgeIndexes(j))

             IF ( (Edge%NodeIndexes(1) == BoundaryMesh%Elements(i)%NodeIndexes(1) .AND. &
                  Edge%NodeIndexes(2) ==  BoundaryMesh%Elements(i)%NodeIndexes(2)) .OR. &
                  (Edge%NodeIndexes(1) == BoundaryMesh%Elements(i)%NodeIndexes(2) .AND. &
                  Edge%NodeIndexes(2) == BoundaryMesh%Elements(i)%NodeIndexes(1))) EXIT
          END DO
          IF (j.GT.3) CALL FATAL(trim(SolverName),'BC Dont found Edge')
          BoundaryMesh % Elements(i) % EdgeIndexes(1)=&
            RefMesh % Elements(k)%BoundaryInfo%Left% EdgeIndexes(j)
      End do
      
CALL CheckTimer('MeshCopy',Delete=.True.)

      Do t=1,MaxRefinementDepth

         write(Message,'(A,I0)') 'Refinement depth :',t
         CALL INFO(trim(SolverName),trim(Message),Level=5)
      
         NofBulk=tmpMesh%NumberOfBulkElements


        !! besoin de trouver la taille de l'element sur le maillage maitre
        !=> comme on conait le parent on peut le faire 

        !!!! MARK elements to split !!!
        !!! Compare Required Element Size (on the master mesh) with
        !local element diameter
        CALL ResetTimer('MarkElement')
        MarkedElement=.False.
        counter=0
        Do i=1,tmpMesh%NumberOfBulkElements
           ParentElement => RefMesh %elements(Parent(i))
           
           !local element diameter
           ElementNodes % x(1:3)=tmpMesh % Nodes % x(tmpMesh%Elements(i) % NodeIndexes(1:3))
           ElementNodes % y(1:3)=tmpMesh % Nodes % y(tmpMesh%Elements(i) % NodeIndexes(1:3))
           ElementNodes % z=0.0
           hk=ElementDiameter( tmpMesh%Elements(i), ElementNodes )

           !Required Element Size (on the master mesh)
           ESize=MinVal(ESVar%Values(ESVar%Perm(ParentElement%NodeIndexes(1:3))))

           !Mark element if hk>ESize
           if (hk.GT.ESize)  then
             MarkedElement(i)=.True.
             counter=counter+1
           endif
        End Do

        CALL CheckTimer('MarkElement',Delete=.True.)
        write(Message,'(A,I0,A,I0)') 'Number of marked elements :',counter,&
                                '/',tmpMesh%NumberOfBulkElements
        CALL INFO(trim(SolverName),trim(Message),Level=6)

        IF (counter.LT.1) exit

        CALL ResetTimer('Split')
        !! ******* Do the splitting *******
        Do i=1,NofBulk

          !!!! Element not in the list
          if (.NOT.MarkedElement(i)) cycle

          write(Message,'(A,I0,A,I0)') &
               'Remaining number of elements to split',COUNT(MarkedElement),&
                                '/',counter
          CALL INFO(trim(SolverName),trim(Message),Level=15)

          Element => tmpMesh % Elements(i)
          j=Element % Type % NumberOfEdges  !TODO: it has to be 3; work only for linear triangles

          !Create Lepp table
          LeppSize=0
          LeppTable=-1
          ! Find element longest edge
          LongestEdge = FindLongestEdge(Element,tmpMesh)

          Done=.False.
          Do While (.NOT.Done)
       
         !*** first time=> create initial Lepp Table
          if (LeppSize<1) then
            Do While(.true.) 
              LeppSize=LeppSize+1
              LeppTable(LeppSize)=Element%ElementIndex

              Edge => tmpMesh % Edges (LongestEdge) 
              NULLIFY(NeighbourElement)
              IF (ASSOCIATED(Edge%BoundaryInfo%Left)) then
                IF (Edge%BoundaryInfo%Left%ElementIndex.NE.Element%ElementIndex) &
                   NeighbourElement => tmpMesh % Elements(Edge%BoundaryInfo%Left%ElementIndex)
              ENDIF
              IF (ASSOCIATED(Edge%BoundaryInfo%Right)) then
                IF (Edge%BoundaryInfo%Right%ElementIndex.NE.Element%ElementIndex) &
                   NeighbourElement => tmpMesh %  Elements(Edge%BoundaryInfo%Right%ElementIndex)
              ENDIF

              !has a neighbour element=> find longestEdge
              IF (ASSOCIATED(NeighbourElement)) THEN 
                NeighbourLongestEdge = FindLongestEdge(NeighbourElement,tmpMesh)

                ! share longest edge -> End Lepp
                If (NeighbourLongestEdge.EQ.LongestEdge) then 
                   LeppTable(LeppSize+1)=NeighbourElement%ElementIndex
                   exit
                ! not same longest edge continue   
                else 
                  Element => NeighbourElement
                  LongestEdge = NeighbourLongestEdge
                Endif
              ELSE !has no neighbour;i.e. boundary Element -> End Lepp
                exit
              ENDIF
            End do

         !*** new round: update lepp table
         !***  remove tN and tN+1 from Lepp table and find new tN+1
          Else
             LeppSize=LeppSize-1
             Element => tmpMesh % Elements(LeppTable(LeppSize))
             LongestEdge = FindLongestEdge(Element,tmpMesh)
             Edge => tmpMesh % Edges (LongestEdge)
             NULLIFY(NeighbourElement)
             IF (ASSOCIATED(Edge%BoundaryInfo%Left)) then
                IF (Edge%BoundaryInfo%Left%ElementIndex.NE.Element%ElementIndex) &
                   NeighbourElement => tmpMesh % Elements(Edge%BoundaryInfo%Left%ElementIndex)
             ENDIF
             IF (ASSOCIATED(Edge%BoundaryInfo%Right)) then
                IF (Edge%BoundaryInfo%Right%ElementIndex.NE.Element%ElementIndex) &
                   NeighbourElement => tmpMesh %  Elements(Edge%BoundaryInfo%Right%ElementIndex)
             ENDIF
             IF (ASSOCIATED(NeighbourElement)) THEN
                LeppTable(LeppSize+1)=NeighbourElement%ElementIndex
             ELSE
                CALL FATAL(trim(SolverName),'Dont find Neighbour')
             ENDIF
           Endif

          !!!! Lepp Size = 1 ; we will split the initially marked element
          ! then quit
           Done=(LeppSize.EQ.1)

          !Insert longest edge mid-point
           NONodes=NONodes+1
           tmpMesh%NumberOfNodes=NONodes
           if (tmpMesh%NumberOfNodes.GT.MaxNumberOfNodes) &
                 CALL FATAL(trim(SolverName),'MaxNumberOfNodes reached')

           Node1 = tmpMesh % Edges(LongestEdge) % NodeIndexes(1)
           Node2 = tmpMesh % Edges(LongestEdge) % NodeIndexes(2)
           x1 = tmpMesh % Nodes % x( Node1 )
           x2 = tmpMesh % Nodes % x( Node2 )
           y1 = tmpMesh % Nodes % y( Node1 )
           y2 = tmpMesh % Nodes % y( Node2 )
           tmpMesh % Nodes % x(NONodes) = (x1+x2) / 2.0d0
           tmpMesh % Nodes % y(NONodes) = (y1+y2) / 2.0d0
           tmpMesh % Nodes % z(NONodes) = 0.0d0

          ! A new Edge will be created by splitting the longest edge
          !   allocate space
           NOEdges=NOEdges+1
           tmpMesh%NumberOfEdges = NOEdges
           if (tmpMesh%NumberOfEdges.GT.MaxNumberOfEdges) &
                 CALL FATAL(trim(SolverName),'MaxNumberOfEdges reached')
           tmpMesh % Edges(NOEdges) = tmpMesh % Edges(LongestEdge)      
           call AllocateVector(tmpMesh % Edges(NOEdges) % NodeIndexes,2) !! WARNING 2 hard coded
           tmpMesh % Edges(NOEdges) % NodeIndexes(1:2)=Edge % NodeIndexes(1:2)
           ALLOCATE(tmpMesh % Edges(NOEdges)% BoundaryInfo)
           tmpMesh % Edges(NOEdges) % BoundaryInfo % Left => NULL()
           tmpMesh % Edges(NOEdges) % BoundaryInfo % Right => NULL()

          !! Split Element LeppTable(Lepp)
           Element => tmpMesh % Elements(LeppTable(LeppSize))
           IF ((NOBulkElements+1).GT.MaxNumberOfElements) &
              CALL FATAL(trim(SolverName),'MaxNumberOfElements reached')
           call SplitElement(Element,tmpMesh,&
              LongestEdge,NOBulkElements,NONodes,NOEdges,1,FacePoint,&
              Parent,Child)

          !! The Updated element has already been splitted in 2 => no more splitting
           MarkedElement(Element%ElementIndex)=.FALSE.

          !! Update left and right Info of newest Edge
           NULLIFY(tmpMesh % Edges(LongestEdge) % BoundaryInfo % Left)
           NULLIFY(tmpMesh % Edges(LongestEdge) % BoundaryInfo % Right)
           NULLIFY(tmpMesh % Edges(NOEdges) % BoundaryInfo %Left)
           NULLIFY(tmpMesh % Edges(NOEdges) % BoundaryInfo %Right)
           tmpMesh % Edges(LongestEdge) % BoundaryInfo % Left => & 
                 tmpMesh % Elements(LeppTable(LeppSize))
           tmpMesh % Edges(NOEdges) % BoundaryInfo % Left => & 
                 tmpMesh % Elements(NOBulkElements)

          !New Edge between the 2 elements           
           if ((NOEdges+1).GT.MaxNumberOfEdges) &
                 CALL FATAL(trim(SolverName),'MaxNumberOfEdges reached')
           tmpMesh % Edges(NOEdges+1) = tmpMesh % Edges(longestEdge)
           call AllocateVector(tmpMesh%Edges(NOEdges+1)%NodeIndexes,2)
           tmpMesh%Edges(NOEdges+1)%NodeIndexes(1:2)=Edge%NodeIndexes(1:2)
           ALLOCATE(tmpMesh%Edges(NOEdges+1)%BoundaryInfo)
           tmpMesh%Edges(NOEdges+1)%BoundaryInfo%Left => NULL()
           tmpMesh%Edges(NOEdges+1)%BoundaryInfo%Right => NULL()

           tmpMesh%Edges(NOEdges+1)%BoundaryInfo % Left => &
                       tmpMesh % Elements(LeppTable(LeppSize))
           tmpMesh%Edges(NOEdges+1) % BoundaryInfo % Right => &
                       tmpMesh % Elements(NOBulkElements)
           tmpMesh % Edges(NOEdges+1)%NodeIndexes(1) = NONodes
           tmpMesh % Edges(NOEdges+1)%NodeIndexes(2) = FacePoint

          !LeppTable(LeppSize+1) is an element to split
           If (LeppTable(LeppSize+1) > 0) then 
            Element => tmpMesh % Elements(LeppTable(LeppSize+1))
            IF ((NOBulkElements+1).GT.MaxNumberOfElements) &
              CALL FATAL(trim(SolverName),'MaxNumberOfElements reached')
            call SplitElement(Element,tmpMesh,&
              LongestEdge,NOBulkElements,NONodes,NOEdges,2,FacePoint,&
              Parent,Child)
           
           !! The Update element has already been splitted in 2 => no more splitting
            MarkedElement(Element%ElementIndex)=.FALSE.

           !! Update left and right Info of newest Edge
            tmpMesh % Edges(LongestEdge) % BoundaryInfo % Right => & 
                 tmpMesh % Elements(LeppTable(LeppSize+1))
            tmpMesh % Edges(NOEdges) % BoundaryInfo % Right => & 
                 tmpMesh % Elements(NOBulkElements)

            !New Edge between the 2 elements           
            if ((NOEdges+2).GT.MaxNumberOfEdges) &
                 CALL FATAL(trim(SolverName),'MaxNumberOfEdges reached')
             tmpMesh % Edges(NOEdges+2) = tmpMesh % Edges(longestEdge)
             call AllocateVector(tmpMesh % Edges(NOEdges+2) % NodeIndexes,2)
             ALLOCATE(tmpMesh%Edges(NOEdges+2)%BoundaryInfo)
             tmpMesh%Edges(NOEdges+2)%BoundaryInfo%Left => NULL()
             tmpMesh%Edges(NOEdges+2)%BoundaryInfo%Right => NULL()
             tmpMesh%Edges(NOEdges+2)%BoundaryInfo%Left => &
                       tmpMesh % Elements(LeppTable(LeppSize+1))
             tmpMesh%Edges(NOEdges+2)%BoundaryInfo%Right => &
                       tmpMesh%Elements(NOBulkElements)

             tmpMesh%Edges(NOEdges+2)%NodeIndexes(1) = NONodes
             tmpMesh%Edges(NOEdges+2)%NodeIndexes(2) = FacePoint
          else ! Lepp finish on a boundary = create new BC element
             Do j=1,BoundaryMesh%NumberOfBoundaryElements
               If (BoundaryMesh%Elements(j)%EdgeIndexes(1).EQ.LongestEdge) exit
            End Do
            If (j.GT.BoundaryMesh%NumberOfBoundaryElements) &
                 CALL FATAL(trim(SolverName),&
                   'Don t find BC element corresponding to longest edge')

             NOBoundaryElements=NOBoundaryElements+1
             BoundaryMesh%NumberOfBoundaryElements=NOBoundaryElements
             IF (NOBoundaryElements.GT.MaxNumberOfBoundaryElements) &
                  CALL FATAL(trim(SolverName),&
                             'MaxNumberOfBoundaryElements reached')

             BoundaryMesh % Elements(NOBoundaryElements) = BoundaryMesh%Elements(j)
             CALL AllocateVector(BoundaryMesh%Elements(NOBoundaryElements)%NodeIndexes,2)
             CALL AllocateVector(BoundaryMesh%Elements(NOBoundaryElements)%EdgeIndexes,1)
             Allocate(BoundaryMesh%Elements(NOBoundaryElements)%BoundaryInfo)

             BoundaryMesh % Elements(NOBoundaryElements)%NodeIndexes(1)=NONodes
             BoundaryMesh % Elements(NOBoundaryElements)%NodeIndexes(2)=tmpMesh%Edges(LongestEdge)%NodeIndexes(2)
             BoundaryMesh % Elements(NOBoundaryElements)%EdgeIndexes(1)=NOEdges
             BoundaryMesh%Elements(j)%NodeIndexes(1)=tmpMesh%Edges(LongestEdge)%NodeIndexes(1) !!HEUH???
             BoundaryMesh%Elements(j)%NodeIndexes(2)=NONodes 

             !IF (ASSOCIATED(BoundaryMesh%Elements(j)%BoundaryInfo%Left)) Then
             !END IF
             !IF (ASSOCIATED(BoundaryMesh%Elements(j)%BoundaryInfo%Right)) Then
             !END IF
            !!! Keep pointing on the initial mesh parent
             BoundaryMesh%Elements(NOBoundaryElements)%BoundaryInfo%Left => &
                 RefMesh%Elements(BoundaryMesh%Elements(j)%BoundaryInfo%Left%ElementIndex)
             
           
             BoundaryMesh%Elements(NOBoundaryElements)%BoundaryInfo%Constraint=&
                 BoundaryMesh%Elements(j)%BoundaryInfo%Constraint
                 
             

          endif

         !! Finally Update Edge NodeIndexes 
         tmpMesh%Edges(NOEdges)%NodeIndexes(1)=NONodes
         tmpMesh%Edges(NOEdges)%NodeIndexes(2) = &
                tmpMesh%Edges(LongestEdge)%NodeIndexes(2)

         tmpMesh%Edges(LongestEdge)%NodeIndexes(2)=NONodes
     
         !! For the new Edges created by splitting the Elements
         If (LeppTable(LeppSize+1) > 0) then
           NOEdges=NOEdges+2
           tmpMesh%NumberOfEdges = NOEdges
         Else
           NOEdges=NOEdges+1
           tmpMesh%NumberOfEdges = NOEdges
         Endif

     !!!!!! End Do while t(i) has to be splitted
        End do

     !! End do Elements to split
       End Do
       CALL CheckTimer('Split',Delete=.True.)

     !!!!! End do refinement depth
      End Do

      CALL INFO(trim(SolverName),&
               'Refinement Done, copy tmpMesh to NewMesh')

      !   Initialize the new mesh:
      NewMesh => AllocateMesh()
      NewMesh%MaxElementNodes=3
      NewMesh%MaxElementDOFs=3
      NewMesh%MeshDim=RefMesh%MeshDim

      !Copy Nodes
      NewMesh % NumberOfNodes = tmpMesh % NumberOfNodes
      CALL AllocateVector( NewMesh % Nodes % x, NewMesh % NumberOfNodes)
      CALL AllocateVector( NewMesh % Nodes % y, NewMesh % NumberOfNodes)
      CALL AllocateVector( NewMesh % Nodes % z, NewMesh % NumberOfNodes)
      NewMesh % Nodes % x(1:tmpMesh % NumberOfNodes) = &
         tmpMesh % Nodes % x(1:tmpMesh % NumberOfNodes)
      NewMesh % Nodes % y(1:tmpMesh % NumberOfNodes) = &
         tmpMesh % Nodes % y(1:tmpMesh % NumberOfNodes)
      NewMesh % Nodes % z(1:tmpMesh % NumberOfNodes) = &
         tmpMesh % Nodes % z(1:tmpMesh % NumberOfNodes)

      !Copy Elements 
      NewMesh % NumberOfBulkElements = tmpMesh % NumberOfBulkElements
      NewMesh % NumberOfBoundaryElements = BoundaryMesh % NumberOfBoundaryElements

      CALL AllocateVector( NewMesh % Elements, NewMesh % NumberOfBulkElements + &
               NewMesh % NumberOfBoundaryElements )

      Do i=1,NewMesh % NumberOfBulkElements
         NewMesh % Elements(i)=tmpMesh % Elements(i)
         NewMesh % Elements(i) % ElementIndex = i
         CALL AllocateVector(NewMesh % Elements(i) % NodeIndexes, 3)
         CALL AllocateVector(NewMesh % Elements(i) % EdgeIndexes, 3)
         NewMesh % Elements(i) % NodeIndexes(1:3)=tmpMesh % Elements(i) % NodeIndexes(1:3)
         NewMesh % Elements(i) % EdgeIndexes(1:3)=tmpMesh % Elements(i) % EdgeIndexes(1:3)
      End do

      !Copy Edges
      NewMesh % NumberOfEdges = tmpMesh % NumberOfEdges
      CALL AllocateVector( NewMesh % Edges , NewMesh % NumberOfEdges)
      Do i=1,NewMesh % NumberOfEdges
         NewMesh % Edges(i)=tmpMesh % Edges(i)
         CALL AllocateVector(NewMesh % Edges(i) % NodeIndexes,2)
         NewMesh % Edges(i) % NodeIndexes(1:2)=tmpMesh % Edges(i) % NodeIndexes(1:2)
         ALLOCATE(NewMesh % Edges(i) % BoundaryInfo)
         NewMesh % Edges(i ) % BoundaryInfo=tmpMesh % Edges(i) % BoundaryInfo
      End Do

      ! Copy BoundaryElements
      Do i=1,NewMesh % NumberOfBoundaryElements
         k=NewMesh % NumberOfBulkElements+i
         deallocate(BoundaryMesh % Elements(i) % EdgeIndexes)
         BoundaryMesh % Elements(i) % EdgeIndexes => NULL()
         NewMesh % Elements(k)=BoundaryMesh % Elements(i)
         NewMesh % Elements(k) % ElementIndex = k
         CALL AllocateVector(NewMesh % Elements(k) % NodeIndexes, 2)
         NewMesh % Elements(k) % NodeIndexes(1:2)=BoundaryMesh % Elements(i) % NodeIndexes(1:2)
         ALLOCATE(NewMesh % Elements(k) % BoundaryInfo)
         !! As BoundaryMesh%Elements(i)%BoundaryInfo still assoiciated
         !      with Solver%mesh%BoundaryInfo => need to find parent
         ParentIndex=BoundaryMesh % Elements(i) % BoundaryInfo % Left % ElementIndex
         CALL SetParent(NewMesh,NewMesh % Elements(k),ParentIndex,Child)
         NewMesh % Elements(k) % BoundaryInfo % Constraint = BoundaryMesh % Elements(i) % BoundaryInfo % Constraint
         deallocate(BoundaryMesh % Elements(i) % BoundaryInfo)
         BoundaryMesh % Elements(i) % BoundaryInfo=>NULL()
      End Do


      ! Release and deallocate meshes
      CALL ReleaseMesh(tmpMesh)
      Deallocate(tmpMesh)
      CALL ReleaseMesh(BoundaryMesh)
      Deallocate(BoundaryMesh)

      END FUNCTION LeppAdaptivity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION FindLongestEdge(Element,Mesh) RESULT(LongestEdge)
       !! Find longest edge of the element
        Type(Mesh_t), POINTER :: Mesh
        Type(Element_t), POINTER :: Element
        INTEGER :: LongestEdge

        REAL(KIND=dp) :: MaxLength,EdgeLength
        REAL(KIND=dp) :: x1,x2,y1,y2
        INTEGER :: EdgeNumber,Node1,Node2
        INTEGER :: i,j

        MaxLength   = 0.0D0
        LongestEdge = 0
        DO j = 1,Element % Type % NumberOfEdges
           EdgeNumber = Element % EdgeIndexes(j)
           Node1 = Mesh % Edges( EdgeNumber ) % NodeIndexes(1)
           Node2 = Mesh % Edges( EdgeNumber ) % NodeIndexes(2)
           x1 = Mesh % Nodes % x( Node1 )
           x2 = Mesh % Nodes % x( Node2 )
           y1 = Mesh % Nodes % y( Node1 )
           y2 = Mesh % Nodes % y( Node2 )
           EdgeLength = SQRT((x2-x1)**2+(y2-y1)**2)
           IF (EdgeLength >= MaxLength) THEN
               MaxLength = EdgeLength
               LongestEdge = EdgeNumber 
           END IF
        END DO
       END FUNCTION FindLongestEdge

       SUBROUTINE SplitElement(Element,Mesh,LongestEdge,NumberOfBulkElements,NewNode,NewEdge,MiddleEdge,FacePoint,Parent,Child)
       !! Split Linear Triangle element in 2 using middle point of the
       !! longest edge
       !!   Update the element and create the second
         Type(Element_t), POINTER :: Element
         INTEGER, Dimension(:) :: Parent
         INTEGER, Dimension(:,:) :: Child
         Type(Mesh_t), POINTER :: Mesh
         INTEGER :: LongestEdge,NumberOfBulkElements,NewNode,NewEdge,MiddleEdge,FacePoint
         INTEGER :: k,n
         INTEGER :: EdgeIndex
         INTEGER :: MaxNumberOfChilds


         DO k = 1,3
           IF ( Element % NodeIndexes(k) /= &
             Mesh % Edges(LongestEdge) % NodeIndexes(1) .AND. &
             Element % NodeIndexes(k) /= &
             Mesh % Edges(LongestEdge) % NodeIndexes(2) ) EXIT
         END DO
         FacePoint=Element % NodeIndexes(k)

         !Update 1rst Element
         Mesh % Elements(Element%ElementIndex) % NodeIndexes(1) = &
              FacePoint
         Mesh % Elements(Element%ElementIndex) % NodeIndexes(2) = &
              Mesh % Edges(LongestEdge) % NodeIndexes(1)
         Mesh % Elements(Element%ElementIndex) % NodeIndexes(3) = &
              NewNode

         !Create New Element
         NumberOfBulkElements=NumberOfBulkElements+1
         Mesh%NumberOfBulkElements = NumberOfBulkElements
         Mesh % Elements(NumberOfBulkElements) = Mesh % Elements(Element%ElementIndex)
         call AllocateVector(Mesh %Elements(NumberOfBulkElements) % NodeIndexes, 3) !! WARNING 3 hard coded
         call AllocateVector(Mesh %Elements(NumberOfBulkElements) % EdgeIndexes, 3)
         Mesh % Elements(NumberOfBulkElements) % EdgeIndexes(1:3) = &
              Element % EdgeIndexes(1:3)
         Mesh % Elements(NumberOfBulkElements) % ElementIndex = NumberOfBulkElements
         Mesh % Elements(NumberOfBulkElements) % BodyId = Element % BodyId

         Mesh % Elements(NumberOfBulkElements) % NodeIndexes(1) = &
               FacePoint
         Mesh % Elements(NumberOfBulkElements) % NodeIndexes(2) = &
               NewNode
         Mesh % Elements(NumberOfBulkElements) % NodeIndexes(3) = &
            Mesh % Edges(LongestEdge) % NodeIndexes(2)

         Parent(NumberOfBulkElements) = Parent(Element%ElementIndex)
         MaxNumberOfChilds=SIZE(Child,2)

         n=1
         Do while((Child(Parent(NumberOfBulkElements),n).GT.0).AND.(n.LE.MaxNumberOfChilds))
           n=n+1
         End do
         IF (n.GT.MaxNumberOfChilds) Then
             CALL FATAL(trim(SolverName),'MaxNumberOfChilds reached')
         END IF
         Child(Parent(NumberOfBulkElements),n)=NumberOfBulkElements

         !! Update EdgesIndexes of Element
         k = FindSameEdge(Mesh,NumberOfBulkElements,3,1)
         Mesh % Elements(NumberOfBulkElements) % EdgeIndexes(3) = &
              Mesh % Elements(NumberOfBulkElements) % EdgeIndexes(k)
         Mesh % Elements(NumberOfBulkElements) % EdgeIndexes(1) = &
              NewEdge+MiddleEdge
         Mesh % Elements(NumberOfBulkElements) % EdgeIndexes(2) = &
              NewEdge

         EdgeIndex=Mesh % Elements(NumberOfBulkElements) % EdgeIndexes(3)
         IF (ASSOCIATED(Mesh%Edges(EdgeIndex)%BoundaryInfo%Left)) then
            IF (Mesh%Edges(EdgeIndex)%BoundaryInfo%Left%ElementIndex.EQ.Element%ElementIndex) then
               Mesh%Edges(EdgeIndex)%BoundaryInfo%Left => Mesh % Elements(NumberOfBulkElements)
             endif
         ENDIF
         IF (ASSOCIATED(Mesh%Edges(EdgeIndex)%BoundaryInfo%Right)) then
            IF (Mesh%Edges(EdgeIndex)%BoundaryInfo%Right%ElementIndex.EQ.Element%ElementIndex) then
               Mesh%Edges(EdgeIndex)%BoundaryInfo%Right => Mesh % Elements(NumberOfBulkElements)
            endif
         ENDIF

         k = FindSameEdge(Mesh,Element%ElementIndex,1,2)
         Mesh % Elements(Element%ElementIndex) % EdgeIndexes(1) = &
              Mesh % Elements(Element%ElementIndex) % EdgeIndexes(k)
         Mesh % Elements(Element%ElementIndex) % EdgeIndexes(2) = &
              LongestEdge
         Mesh % Elements(Element%ElementIndex) % EdgeIndexes(3) = &
              NewEdge+MiddleEdge
              
      END SUBROUTINE SplitElement

      FUNCTION FindSameEdge(Mesh,EIndex,Node1,Node2) RESULT(k)
      ! Find the Number k of the Edge between Node1 and Node2
      TYPE(Mesh_t),POINTER :: Mesh
      INTEGER :: EIndex,Node1,Node2
      INTEGER :: k
       ! go through the edges and check the one between node1 and node2
        DO k = 1,3 
            if (((Mesh % Elements(EIndex) % NodeIndexes(Node1) == &
                Mesh%Edges(Mesh % Elements(EIndex) % EdgeIndexes(k)) % NodeIndexes(1)) &
                 .AND. &
                (Mesh % Elements(EIndex) % NodeIndexes(Node2) == &
                Mesh%Edges(Mesh % Elements(EIndex) % EdgeIndexes(k)) % NodeIndexes(2))) &
                 .OR. &
               ((Mesh % Elements(EIndex) % NodeIndexes(Node1) == &
                Mesh%Edges(Mesh % Elements(EIndex) % EdgeIndexes(k)) % NodeIndexes(2)) &
                 .AND. &
                (Mesh % Elements(EIndex) % NodeIndexes(Node2) == &
                Mesh%Edges(Mesh % Elements(EIndex) % EdgeIndexes(k)) % NodeIndexes(1)))) &
                    exit
         END DO

         IF (k.GT.3) CALL FATAL(trim(SolverName),'FindSameEdge FAILED')
        ! IF (k.GT.3) PRINT *,'************************************'
        ! IF (k.GT.3) PRINT *,'           FindSameEdge FAILED      '
        ! IF (k.GT.3) PRINT *, EIndex,Node1,Node2
        ! IF (k.GT.3) PRINT *, Mesh % Elements(EIndex) %NodeIndexes(Node1),Mesh % Elements(EIndex) % NodeIndexes(Node2)
        ! IF (k.GT.3) PRINT *, Mesh%Edges(Mesh % Elements(EIndex) % EdgeIndexes(1)) % NodeIndexes(1:2)
        ! IF (k.GT.3) PRINT *, Mesh%Edges(Mesh % Elements(EIndex) % EdgeIndexes(2)) % NodeIndexes(1:2)
        ! IF (k.GT.3) PRINT *, Mesh%Edges(Mesh % Elements(EIndex) % EdgeIndexes(3)) % NodeIndexes(1:2)
        ! IF (k.GT.3) PRINT *,'************************************'
      END FUNCTION

      SUBROUTINE SetParent(Mesh,BCElement,ParentIndex,Child)
       Type(Mesh_t) :: Mesh
       Type(Element_t) :: BCElement
       Type(Element_t), POINTER :: Element,Edge
       INTEGER :: ParentIndex
       INTEGER, DIMENSION(:,:) :: Child
       INTEGER :: n,i
       INTEGER :: MaxNumberOfChilds
      
        BCElement%BoundaryInfo%Left=>NULL()
        BCElement%BoundaryInfo%Right=>NULL()

        n=1
        NULLIFY(Element)
        NULLIFY(Edge)
        MaxNumberOfChilds=SIZE(Child,2)
        Do while((Child(ParentIndex,n).GT.0).AND.(n.LE.MaxNumberOfChilds))
          Element => Mesh % Elements(Child(ParentIndex,n)) 
          Do i=1,3
             Edge => Mesh % Edges(Element%EdgeIndexes(i))
             if ( ((BCElement%NodeIndexes(1).EQ.Edge%NodeIndexes(1)).AND.&
                  (BCElement%NodeIndexes(2).EQ.Edge%NodeIndexes(2))).OR.&
                  ((BCElement%NodeIndexes(1).EQ.Edge%NodeIndexes(2)).AND.&
                  (BCElement%NodeIndexes(2).EQ.Edge%NodeIndexes(1)))) then
                BCElement%BoundaryInfo%Left=>Mesh % Elements(Child(ParentIndex,n))
                Return
              endif
          End do
          n=n+1
        End do
        CALL FATAL(trim(SolverName),' SetParent: Parent Not found')
 
      END SUBROUTINE SetParent

      SUBROUTINE MakeMeshGlobal(NewMesh,RefMesh)
      TYPE(Mesh_t),POINTER :: NewMesh,RefMesh
      TYPE( Variable_t ), POINTER :: Var,NewVar
      CHARACTER(LEN=1024) :: Path
      LOGICAL :: Found

        NewMesh % MaxBDOFs = RefMesh % MaxBDOFs

        NewMesh % Name = ListGetString( Solver % Values, &
           'Adaptive Mesh Name', Found )
       IF ( .NOT. Found ) NewMesh % Name = 'MyRefinedMesh'
       Path = TRIM(NewMesh % Name)
         CALL MakeDirectory( TRIM(path) // CHAR(0) )
       IF ( ListGetLogical( Solver % Values, 'Adaptive Save Mesh', Found ) ) &
         CALL WriteMeshToDisk( NewMesh, Path )

     !
     !   Initialize local variables for the new mesh:
     !   --------------------------------------------
      ALLOCATE( NewMesh % Variables )
      NULLIFY( NewMesh % Variables )

      CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
         'Coordinate 1', 1, NewMesh % Nodes % x )

      CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
         'Coordinate 2', 1, NewMesh % Nodes % y )

      CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
         'Coordinate 3', 1, NewMesh % Nodes % z )

     !   Time must always be there:
     !   --------------------------
      Var => VariableGet( RefMesh % Variables,'Time',ThisOnly=.TRUE. )
      CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                     'Time', 1, Var % Values )

      Var => VariableGet( RefMesh % Variables,'Timestep',ThisOnly=.TRUE. )
       CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                     'Timestep', 1, Var % Values )

      Var => VariableGet( RefMesh % Variables,'Timestep size',ThisOnly=.TRUE. )
       CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                     'Timestep size', 1, Var % Values )

      NewVar => VariableGet( NewMesh % Variables,'Timestep size',ThisOnly=.TRUE. )
      NewVar % PrevValues => Var % PrevValues

      Var => VariableGet( RefMesh % Variables,'Timestep interval',ThisOnly=.TRUE. )
      CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                 'Timestep interval', 1, Var % Values )

      Var => VariableGet( RefMesh % Variables,'Coupled iter',ThisOnly=.TRUE. )
      CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                 'Coupled iter', 1, Var % Values )

      Var => VariableGet( RefMesh % Variables,'Nonlin iter',ThisOnly=.TRUE. )
      CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                 'Nonlin iter', 1, Var % Values )

      ! Initialize the field variables for the new mesh. These are
      ! interpolated from the old meshes variables. Vector variables
      ! are in the variable lists in two ways: as vectors and as
      ! vector components. We MUST update the vectors (i.e. DOFs>1)
      ! first!!!!!
      ! -----------------------------------------------------------
      CALL SetCurrentMesh( Model, NewMesh )
      Var => RefMesh % Variables
      DO WHILE( ASSOCIATED( Var ) )
         IF ( Var % DOFs > 1 ) THEN
            NewVar => VariableGet( NewMesh % Variables,Var % Name,.FALSE. )
            k = SIZE( NewVar % Values )
            IF ( ASSOCIATED( NewVar % Perm ) ) THEN
              k = COUNT( NewVar % Perm > 0 )
            END IF
            IF ( GetVarName( NewVar ) == 'flow solution' ) THEN
               NewVar % Norm = 0.0d0
               DO i=1,NewMesh % NumberOfNodes
                  DO j=1,NewVar % DOFs-1
                     NewVar % Norm = NewVar % Norm + &
                          NewVar % Values( NewVar % DOFs*(i-1)+j )**2
                  END DO
               END DO
               NewVar % Norm = SQRT( NewVar % Norm / k )
            ELSE
               NewVar % Norm = SQRT( SUM(NewVar % Values**2) / k )
            END IF
         END IF
         Var => Var % Next
      END DO

     !   Second time around, update scalar variables and
     !   vector components:
     !   -----------------------------------------------
      Var => RefMesh % Variables
      DO WHILE( ASSOCIATED( Var ) )

        IF( SIZE( Var % Values ) == Var % DOFs ) THEN
           Var => Var % Next
           CYCLE
         END IF

         SELECT CASE( Var % Name )
         CASE( 'coordinate 1', 'coordinate 2', 'coordinate 3', 'time' , &
             'timestep', 'timestep size', 'timestep interval', &
             'coupled iter', 'nonlin iter' )
         CASE DEFAULT
            IF ( Var % DOFs == 1 ) THEN
               Found = .FALSE.
               IF ( Found ) THEN
                  k = Solver % Variable % NameLen
                  IF ( Var % Name(1:k) /= Solver % Variable % Name(1:k) ) THEN
                     Var => Var % Next
                     CYCLE
                  END IF
               END IF

               NewVar => VariableGet( NewMesh % Variables, Var % Name, .FALSE. )
               NewVar % PrimaryMesh => NewMesh
               k = SIZE( NewVar % Values )
               IF ( ASSOCIATED( NewVar % Perm ) ) THEN
                  k = COUNT( NewVar % Perm > 0 )
               END IF
               NewVar % Norm = SQRT( SUM(NewVar % Values**2) / k )
            END IF
         END SELECT
         Var => Var % Next
      END DO
!
!   Update Solver structure to use the new mesh:
!   ---------------------------------------------    
     CALL MeshStabParams( NewMesh )
!
!   Nothing computed on this mesh yet:
!   ----------------------------------
      NewMesh % SavesDone    = 0  ! start new output file
      NewMesh % OutputActive = .True.

      NewMesh % Changed   = .TRUE.

      END SUBROUTINE MakeMeshGlobal


!------------------------------------------------------------------------------
      SUBROUTINE MyUpdateSolverMesh( Solver, Mesh )
      !!!! A copy of the subroutine in ElementUtils to allow
      !deallocation of previous Matrix structure
!------------------------------------------------------------------------------
      TYPE( Mesh_t ), POINTER :: Mesh
      TYPE( Solver_t ), TARGET :: Solver
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,n,n1,n2,DOFs
      LOGICAL :: Found, OptimizeBandwidth
      TYPE(Matrix_t), POINTER   :: Matrix
      REAL(KIND=dp), POINTER :: Work(:)
      INTEGER, POINTER :: Permutation(:)
      TYPE(Variable_t), POINTER :: TimeVar, SaveVar
      TYPE(Matrix_t), POINTER   :: PARENT
      integer :: nParents
!------------------------------------------------------------------------------
      SaveVar => Solver % Variable
      DOFs = SaveVar % DOFs

      Solver % Mesh => Mesh
      CALL SetCurrentMesh( CurrentModel, Mesh )
!
!    Create matrix and variable structures for
!    current equation on the new mesh:
!    -----------------------------------------
      Solver % Variable => VariableGet( Mesh % Variables, &
        Solver % Variable % Name, ThisOnly = .FALSE. )

      CALL AllocateVector( Permutation, SIZE(Solver % Variable % Perm) )

      OptimizeBandwidth = ListGetLogical( Solver % Values, &
                                         'Optimize Bandwidth', Found )
      IF ( .NOT. Found ) OptimizeBandwidth = .TRUE.

      Matrix => CreateMatrix( CurrentModel, Solver, &
        Mesh, Permutation, DOFs, MATRIX_CRS, OptimizeBandwidth, &
        ListGetString( Solver % Values, 'Equation' ) )

      Matrix % Symmetric = ListGetLogical( Solver % Values, &
             'Linear System Symmetric', Found )

      Matrix % Lumped = ListGetLogical( Solver % Values, &
             'Lumped Mass Matrix', Found )

      ALLOCATE( Work(SIZE(Solver % Variable % Values)) )
      Work = Solver % Variable % Values
      DO k=0,DOFs-1
        DO i=1,SIZE(Permutation)
           IF ( Permutation(i) > 0 ) THEN
              Solver % Variable % Values( DOFs*Permutation(i)-k ) = &
                 Work( DOFs*Solver % Variable % Perm(i)-k )
           END IF
        END DO
      END DO

      IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        DO j=1,SIZE(Solver % Variable % PrevValues,2)
           Work = Solver % Variable % PrevValues(:,j)
           DO k=0,DOFs-1
              DO i=1,SIZE(Permutation)
                 IF ( Permutation(i) > 0 ) THEN
                    Solver % Variable % PrevValues( DOFs*Permutation(i) - k,j ) =  &
                        Work( DOFs * Solver % Variable % Perm(i) - k )
                  END IF
              END DO
           END DO
        END DO
      END IF
      DEALLOCATE( Work )

      Solver % Variable % Perm = Permutation
      Solver % Variable % Solver => Solver

     DEALLOCATE( Permutation )
     CALL AllocateVector( Matrix % RHS, Matrix % NumberOfRows )

      IF ( ASSOCIATED(SaveVar % EigenValues) ) THEN
        n = SIZE(SaveVar % EigenValues)

        IF ( n > 0 ) THEN
           Solver % NOFEigenValues = n
           CALL AllocateVector( Solver % Variable % EigenValues,n )
           CALL AllocateArray( Solver % Variable % EigenVectors, n, &
                    SIZE(Solver % Variable % Values) ) 

           Solver % Variable % EigenValues  = 0.0d0
           Solver % Variable % EigenVectors = 0.0d0

           CALL AllocateVector( Matrix % MassValues, SIZE(Matrix % Values) )
           Matrix % MassValues = 0.0d0
        END IF
      ELSE IF ( ASSOCIATED( Solver % Matrix % Force ) ) THEN
        n1 = Matrix % NumberOFRows
        n2 = SIZE(Solver % Matrix % Force,2)
        ALLOCATE(Matrix % Force(n1,n2))
        Matrix % Force = 0.0d0
      END IF
 
      !!!! Fab lines Added
      !!!! IN CASE WE USE MULTIGRID PARENT MATRIX ARE CREATED AND SHOULD
      !BE RECOMPUTED => DEALLOCATE PARENTS
      !!!   SHOULD WE CHECK FOR CHILDS ??
      IF ( ASSOCIATED( Solver % Matrix ) ) THEN
        nParents=0
        PARENT => Solver % Matrix % PARENT
        IF (ASSOCIATED(PARENT)) THEN
          DO WHILE (ASSOCIATED(PARENT))
             PARENT => PARENT %  PARENT
             nParents=nParents+1
          END DO
          Do While (nParents.GT.0)
             PARENT => Solver % Matrix % PARENT
             Do j=1,nParents-1
                PARENT => PARENT %  PARENT
             End do
             IF ( ASSOCIATED( PARENT % EMatrix ) ) CALL FreeMatrix( PARENT % EMatrix ) 
             IF ( ASSOCIATED( PARENT ) ) CALL FreeMatrix( PARENT ) 
             nParents=nParents-1
          End do
        END IF
      !!!  FINALLY DEALLOCATE Solver % Matrix
        CALL FreeMatrix( Solver % Matrix )
        NULLIFY( Solver % Matrix )
      END IF
      !!!! Fab End: lines Added

      Solver % Matrix => Matrix
      Solver % Mesh % Changed = .TRUE.

!------------------------------------------------------------------------------
      END SUBROUTINE MyUpdateSolverMesh
!------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FUNCTIONS to split Mesh using RGB algorithm
!!!! This is a copy adaptation of the material found in the "Adaptive" module under elmerfem/fem
!------------------------------------------------------------------------------
      FUNCTION SplitMesh( RefMesh ) RESULT(NewMesh)
!------------------------------------------------------------------------------
      TYPE(Mesh_t), POINTER :: NewMesh, RefMesh
!------------------------------------------------------------------------------
      TYPE(Mesh_t), POINTER :: NewMesh1
      INTEGER :: i,j,k,n,MarkedElements
      INTEGER :: Lambda
      TYPE(Element_t), POINTER :: RefElement
      REAL(Kind=dp) :: hk,ESize
      Type(Nodes_t),SAVE :: ElementNodes
      LOGICAL,SAVE :: AllocationDone=.False.
!------------------------------------------------------------------------------

      NULLIFY( NewMesh )


      IF (.NOT.AllocationDone) THEN
         allocate(ElementNodes%x(3),&
                  ElementNodes%y(3),&
                  ElementNodes%z(3))
         AllocationDone=.true.
      END IF
    
!   Determine the marked elements:
!   ------------------------------
      MarkedElements = 0

      DO i = 1,RefMesh % NumberOfBulkElements
       RefElement => RefMesh % Elements(i)

       IF ( RefElement % TYPE % ElementCode /= 303 ) THEN
          CALL Fatal( 'SplitMesh',&
          'Internal splitting implemented only for linear triangles.' )
       END IF

       ! Compute Element Diameter    
       ElementNodes % x(1:3)=RefMesh % Nodes % x(RefElement % NodeIndexes(1:3))
       ElementNodes % y(1:3)=RefMesh % Nodes % y(RefElement % NodeIndexes(1:3))
       ElementNodes % z=0.0
       hk=ElementDiameter( RefElement , ElementNodes )
 
       ! Required Element Size
       ESize=MinVal(ESVar%Values(ESVar%Perm(RefElement % NodeIndexes(1:3))))

         
       ! Element has to be splitted hk/ESize times
       Lambda=hk/ESize
       RefElement % Splitted = 0
       RefElement % Splitted = MIN( MaxRefinementDepth, Lambda )

       IF ( RefElement % Splitted > 0 ) MarkedElements = MarkedElements  + 1
      END DO

!   PRINT*,MarkedElements,' marked elements'

      IF ( MarkedElements == 0 ) THEN
       RefMesh % Changed = .FALSE.
       RETURN
      END IF

!   Refine until all elements splitted specified times:
!   ---------------------------------------------------
    NewMesh => SplitOneLevel( RefMesh )
    DO WHILE( .TRUE. )
       MarkedElements = 0
       DO i=1,NewMesh % NumberOfBulkElements
          IF ( NewMesh % Elements(i) % Splitted > 0 ) THEN
             MarkedElements = MarkedElements + 1
          END IF
       END DO

       IF ( MarkedElements == 0 ) EXIT

       NewMesh1 => SplitOneLevel( NewMesh )
       CALL ReleaseMesh( NewMesh )
       DEALLOCATE( NewMesh )

       NewMesh => NewMesh1
    END DO

!------------------------------------------------------------------------------
  END FUNCTION SplitMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION SplitOneLevel( RefMesh ) RESULT( NewMesh )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE( Mesh_t ), POINTER :: RefMesh, NewMesh
!------------------------------------------------------------------------------
     !REAL(KIND=dp) :: CPUTime,t
    REAL(KIND=dp) :: t

    INTEGER :: EdgeNumber,LongestEdge,Node1,Node2
    INTEGER :: i,j,k,l,n,NewElCnt,NewNodeCnt,MarkedEdges

    TYPE(Element_t), POINTER :: RefElement,Parent,Child,Edge

    LOGICAL, POINTER :: EdgeSplitted(:)
    INTEGER, POINTER :: MarkedOrder(:), Children(:,:)

    TYPE(PElementDefs_t), POINTER :: PD
    REAL(KIND=dp) :: x1, x2, y1, y2, EdgeLength, MaxLength
!------------------------------------------------------------------------------

    t = CPUTime()
    CALL FindMeshEdges( RefMesh )
    WRITE( Message, * ) 'Find mesh edges time (cpu-secs):                 ',CPUTime()-t
    CALL Info( 'SplitOneLevel', Message, Level=6 )

!   RGB Refinement:
!   ---------------
    t = CPUTime()
    CALL AllocateVector( EdgeSplitted, RefMesh % NumberOfEdges )
    MarkedEdges = RGBRefinement( EdgeSplitted,RefMesh )
    WRITE( Message, * ) 'RGB Refinement time (cpu-secs):                  ',CPUTime()-t
    CALL Info( 'SplitOneLevel', Message, Level=6 )

!   Initialize the new mesh:
!   ------------------------
    NewMesh => AllocateMesh()
    NewMesh % MaxElementNodes = 3
    NewMesh % MaxElementDOFs  = 3
    NewMesh % MeshDim = RefMesh % MeshDim

!   Create node tables for the new mesh:
!   ------------------------------------    
    t = CPUTime()
    NewMesh % NumberOfNodes = RefMesh % NumberOfNodes + MarkedEdges
    CALL AllocateVector( NewMesh % Nodes % x, NewMesh % NumberOfNodes )
    CALL AllocateVector( NewMesh % Nodes % y, NewMesh % NumberOfNodes )
    CALL AllocateVector( NewMesh % Nodes % z, NewMesh % NumberOfNodes )

!   Add old nodes to the new mesh:
!   ------------------------------    
    NewMesh % Nodes % x(1:RefMesh % NumberOfNodes) = &
               RefMesh % Nodes % x(1:RefMesh % NumberOfNodes)
    NewMesh % Nodes % y(1:RefMesh % NumberOfNodes) = &
               RefMesh % Nodes % y(1:RefMesh % NumberOfNodes)
    NewMesh % Nodes % z(1:RefMesh % NumberOfNodes) = &
               RefMesh % Nodes % z(1:RefMesh % NumberOfNodes)

!   Add new nodes to the new mesh:
!   ------------------------------    
    NewNodeCnt = RefMesh % NumberOfNodes
    DO i = 1,RefMesh % NumberOfEdges
       IF ( EdgeSplitted(i) ) THEN
          Node1 = RefMesh % Edges(i) % NodeIndexes(1)
          Node2 = RefMesh % Edges(i) % NodeIndexes(2)
          x1 = RefMesh % Nodes % x(Node1)
          x2 = RefMesh % Nodes % x(Node2)
          y1 = RefMesh % Nodes % y(Node1)
          y2 = RefMesh % Nodes % y(Node2)
          NewNodeCnt = NewNodeCnt + 1
          NewMesh % Nodes % x(NewNodeCnt) = (x1+x2) / 2.0d0
          NewMesh % Nodes % y(NewNodeCnt) = (y1+y2) / 2.0d0
          NewMesh % Nodes % z(NewNodeCnt) = 0.0d0
       END IF
    END DO
    WRITE( Message, * ) 'Node tables generation time (cpu-secs):          ',CPUTime()-t
    CALL Info( 'SplitOneLevel', Message, Level=6 )

!   Count the new number of bulk elements:
!   --------------------------------------
    CALL AllocateVector( MarkedOrder, RefMesh % NumberOfEdges )
    MarkedOrder = 0

    k = 0
    NewElCnt = 0
    DO i = 1,RefMesh % NumberOfBulkElements
       MarkedEdges = 0
       DO j = 1,3
          EdgeNumber = RefMesh % Elements(i) % EdgeIndexes(j)
          IF( EdgeSplitted(EdgeNumber) ) THEN
             MarkedEdges = MarkedEdges + 1
             IF ( MarkedOrder(EdgeNumber) == 0 ) THEN
                k = k + 1
                MarkedOrder(EdgeNumber) = k + RefMesh % NumberOfNodes
             END IF
          END IF
       END DO
       NewElCnt = NewElCnt + MarkedEdges + 1
    END DO
    NewMesh % NumberOfBulkElements = NewElCnt
!
!   Count the new number of boundary elements:
!   ------------------------------------------
    NewElCnt = 0
    DO i = RefMesh % NumberOfBulkElements+1,RefMesh % NumberOfBulkElements+&
         RefMesh % NumberOfBoundaryElements

       RefElement => RefMesh % Elements(i) % BoundaryInfo % Left
       IF ( .NOT.ASSOCIATED( RefElement) ) &
            RefElement => RefMesh % Elements(i) % BoundaryInfo % Right

       IF ( ASSOCIATED( RefElement ) ) THEN
          NULLIFY( Edge )

          DO j=1,3
             Edge => RefMesh % Edges(RefElement % EdgeIndexes(j))

             IF ( Edge % NodeIndexes(1) == RefMesh % Elements(i) % NodeIndexes(1) .AND. &
                  Edge % NodeIndexes(2) == RefMesh % Elements(i) % NodeIndexes(2) .OR.  &
                  Edge % NodeIndexes(2) == RefMesh % Elements(i) % NodeIndexes(1) .AND. &
                  Edge % NodeIndexes(1) == RefMesh % Elements(i) % NodeIndexes(2) ) EXIT
          END DO
   
          IF ( EdgeSplitted( RefElement % EdgeIndexes(j) ) ) THEN
             NewElCnt = NewElCnt + 2
          ELSE
             NewElCnt = NewElCnt + 1
          END IF
       ELSE
          NewElCnt = NewElCnt + 1
       END IF
    END DO

    NewMesh % NumberOfBoundaryElements = NewElCnt

!   Allocate element tables:
!   ------------------------
    t = CPUTime()
    CALL AllocateVector( NewMesh % Elements, NewMesh % NumberOfBulkElements + &
         NewMesh % NumberOfBoundaryElements )

    CALL AllocateArray( Children, RefMesh % NumberOfBulkElements + &
             RefMesh % NumberOfBoundaryElements, 4 )
    Children = 0

!   Find the new bulk elements:
!   ---------------------------
    NewElCnt    = 0
    DO i = 1,RefMesh % NumberOfBulkElements
       RefElement => RefMesh % Elements(i)
       n = RefElement % TYPE % NumberOfNodes

       MarkedEdges = 0
       DO j = 1,3
          EdgeNumber = RefElement % EdgeIndexes(j)
          IF ( EdgeSplitted(EdgeNumber) ) THEN
             MarkedEdges = MarkedEdges + 1
          END IF
       END DO

!      Make elements for the new mesh:
!      --------------------------------
       SELECT CASE(MarkedEdges)
       CASE(0)
!         Just copy of the old one:
!         -------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
          NewMesh % Elements(NewElCnt) % NodeIndexes(1:n) = &
               RefElement % NodeIndexes(1:n)

          Children(i,1) = NewElCnt
          
!-------------------------------------------------------------------------
       CASE(1)
!         Bisect the longest edge to give two triangles:
!         ----------------------------------------------
          DO j = 1,3
             EdgeNumber = RefElement % EdgeIndexes(j)
             IF ( EdgeSplitted( EdgeNumber ) ) EXIT
          END DO
            
!         Find node (k) opposite to the splitted edge:
!         --------------------------------------------
          DO k = 1,3
             IF ( RefElement % NodeIndexes(k) /= &
                  RefMesh % Edges(EdgeNumber) % NodeIndexes(1) .AND. &
                  RefElement % NodeIndexes(k) /= &
                  RefMesh % Edges(EdgeNumber) % NodeIndexes(2) ) EXIT
          END DO

!         New element 1
!         -------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
               RefElement % NodeIndexes(k)
          NewMesh % Elements(NewElCnt) % NodeIndexes(2) = &
               RefMesh % Edges(EdgeNumber) % NodeIndexes(1)

          NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
               MarkedOrder(RefElement % EdgeIndexes(j))

          Children(i,1) = NewElCnt

!         New element 2
!----------------------------------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
               RefElement % NodeIndexes(k)

          NewMesh % Elements(NewElCnt) % NodeIndexes(2) = &
               MarkedOrder(RefElement % EdgeIndexes(j))

          NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
               RefMesh % Edges(EdgeNumber) % NodeIndexes(2)

          Children(i,2) = NewElCnt

!-------------------------------------------------------------------------
       CASE(2)
!         Bisect two of the edges to give three new elements:
!         ---------------------------------------------------

!         Find the edge NOT splitted:
!         ---------------------------
          DO j = 1,3
             EdgeNumber = RefElement % EdgeIndexes(j)
             IF ( .NOT.EdgeSplitted( EdgeNumber ) ) EXIT             
          END DO

!         Find node (k) opposite to the edge NOT splitted:
!         ------------------------------------------------
          DO k = 1,3
             IF (RefElement % NodeIndexes(k) /= &
                  RefMesh % Edges(EdgeNumber) % NodeIndexes(1) .AND. &
                  RefElement % NodeIndexes(k) /= &
                  RefMesh % Edges(EdgeNumber) % NodeIndexes(2) ) EXIT
          END DO

!         New element 1
!----------------------------------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
               RefElement % NodeIndexes(k)

          l = 1
          DO k = 1,3
             IF ( k /= j ) THEN
                l = l + 1
                NewMesh % Elements(NewElCnt) % NodeIndexes(l) = &
                     MarkedOrder(RefElement % EdgeIndexes(k))
             END IF
          END DO

          Children(i,1) = NewElCnt

!         New element 2
!----------------------------------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          l = 0
          DO k = 1,3
             IF ( k /= j ) THEN
                l = l + 1
                NewMesh % Elements(NewElCnt) % NodeIndexes(l) = &
                     MarkedOrder(RefElement % EdgeIndexes(k))
             END IF
          END DO

          MaxLength = 0.0d0
          DO k = 1,3
             IF ( k /= j ) THEN
                EdgeNumber = RefElement % EdgeIndexes(k)
                Node1 = RefMesh % Edges( EdgeNumber ) % NodeIndexes(1)
                Node2 = RefMesh % Edges( EdgeNumber ) % NodeIndexes(2)
                x1 = RefMesh % Nodes % x( Node1 )
                x2 = RefMesh % Nodes % x( Node2 )
                y1 = RefMesh % Nodes % y( Node1 )
                y2 = RefMesh % Nodes % y( Node2 )
                EdgeLength = SQRT((x2-x1)**2+(y2-y1)**2)
                IF (EdgeLength >= MaxLength) THEN
                   MaxLength = EdgeLength
                   LongestEdge = k
                END IF
             END IF
          END DO
          k = LongestEdge
IF ( k <= 0 .OR. k > 3 ) PRINT*,k

          IF ( RefMesh % Edges(RefElement % EdgeIndexes(j)) % NodeIndexes(1) ==  &
               RefMesh % Edges(RefElement % EdgeIndexes(k)) % NodeIndexes(1) .OR.&
               RefMesh % Edges(RefElement % EdgeIndexes(j)) % NodeIndexes(1) ==  &
               RefMesh % Edges(RefElement % EdgeIndexes(k)) % NodeIndexes(2) ) THEN
             NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
                  RefMesh % Edges(RefElement % EdgeIndexes(j)) % NodeIndexes(2)
          ELSE
             NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
                  RefMesh % Edges(RefElement % EdgeIndexes(j)) % NodeIndexes(1)
          END IF

          Children(i,2) = NewElCnt

!         New element 3
!----------------------------------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          DO j = 1,3
             EdgeNumber = RefElement % EdgeIndexes(j)
             IF ( .NOT.EdgeSplitted( EdgeNumber ) ) EXIT             
          END DO

          DO k = 1,2
             NewMesh % Elements(NewElCnt) % NodeIndexes(k) = &
                  RefMesh % Edges(EdgeNumber) % NodeIndexes(k)
          END DO

          NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
               MarkedOrder(RefElement % EdgeIndexes(LongestEdge))

          Children(i,3) = NewElCnt

!-------------------------------------------------------------------------
       CASE(3)
!         Bisect all the edges to give four new elements:
!         -----------------------------------------------

!         New element 1
!----------------------------------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
               RefElement % NodeIndexes(1)

          j = RefElement % EdgeIndexes(1)
          NewMesh % Elements(NewElCnt) % NodeIndexes(2) = MarkedOrder(j)

          j = RefElement % EdgeIndexes(3)
          NewMesh % Elements(NewElCnt) % NodeIndexes(3) = MarkedOrder(j)

          Children(i,1) = NewElCnt

!         New element 2
!----------------------------------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
               RefElement % NodeIndexes(2)

          j = RefElement % EdgeIndexes(2)
          NewMesh % Elements(NewElCnt) % NodeIndexes(2) = MarkedOrder(j)

          j = RefElement % EdgeIndexes(1)
          NewMesh % Elements(NewElCnt) % NodeIndexes(3) = MarkedOrder(j)

          Children(i,2) = NewElCnt

!         New element 3
!----------------------------------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
               RefElement % NodeIndexes(3)

          j = RefElement % EdgeIndexes(3)
          NewMesh % Elements(NewElCnt) % NodeIndexes(2) = MarkedOrder(j)

          j = RefElement % EdgeIndexes(2)
          NewMesh % Elements(NewElCnt) % NodeIndexes(3) = MarkedOrder(j)

          Children(i,3) = NewElCnt

!         New element 4
!----------------------------------------------------
          NewElCnt = NewElCnt + 1
          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

          DO j=1,n
             NewMesh % Elements(NewElCnt) % NodeIndexes(j) = &
                  MarkedOrder( RefElement % EdgeIndexes(j) )
          END DO

          Children(i,4) = NewElCnt
!----------------------------------------------------
       END SELECT

!----------------------------------------------------
       DO j=1,4
          k = Children(i,j)
          IF ( k > 0 ) THEN
             NewMesh % Elements(k) % Splitted = RefElement % Splitted-1
          END IF
       END DO
    END DO


    WRITE( Message, * ) 'Bulk element tables generation time (cpu-secs):  ',CPUTime()-t
    CALL Info( 'SplitOneLevel', Message, Level=6 )
    
!
!   Update boundary elements:
!   -------------------------
    t = CPUTime()
    NewElCnt = NewMesh % NumberOfBulkElements
    DO j = RefMesh % NumberOfBulkElements + 1, &
       RefMesh % NumberOfBulkElements + &
          RefMesh % NumberOfBoundaryElements

       RefElement => RefMesh % Elements(j) % BoundaryInfo % Left
       IF ( .NOT.ASSOCIATED( RefElement) ) &
            RefElement => RefMesh % Elements(j) % BoundaryInfo % Right

       IF ( ASSOCIATED( RefElement ) ) THEN
          NULLIFY( Edge )
          DO i=1,3
             Edge => RefMesh % Edges(RefElement % EdgeIndexes(i))
             IF ( Edge % NodeIndexes(1) == RefMesh % Elements(j) % NodeIndexes(1) .AND. &
                  Edge % NodeIndexes(2) == RefMesh % Elements(j) % NodeIndexes(2) .OR.  &
                  Edge % NodeIndexes(2) == RefMesh % Elements(j) % NodeIndexes(1) .AND. &
                  Edge % NodeIndexes(1) == RefMesh % Elements(j) % NodeIndexes(2) ) EXIT
          END DO
          EdgeNumber = RefElement % EdgeIndexes(i)

          RefElement => RefMesh % Elements(j)
          n = RefElement % TYPE % NumberOfNodes
            
          IF ( EdgeSplitted(EdgeNumber) ) THEN
!
!            New element 1:
!            --------------
             NewElCnt = NewElCnt + 1
             NewMesh % Elements(NewElCnt) = RefElement
             NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
             CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
             NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
                  RefElement % NodeIndexes(1)
             NewMesh % Elements(NewElCnt) % NodeIndexes(2) = &
                  MarkedOrder(EdgeNumber)

             ALLOCATE( NewMesh % Elements(NewElCnt) % BoundaryInfo )
             NewMesh % Elements(NewElCnt) % BoundaryInfo = &
                  RefElement % BoundaryInfo

             NULLIFY( NewMesh % Elements(NewElCnt) % &
                Boundaryinfo % GebhardtFactors )
               
             CALL SetParents( NewMesh % Elements(NewElCnt), &
                  NewMesh, Children, Edge )

             Children(j,1) = NewElCnt
               
!
!            New element 2:
!            --------------
             NewElCnt = NewElCnt + 1
             NewMesh % Elements(NewElCnt) = RefElement
             NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
             CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
             NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
                  MarkedOrder(EdgeNumber)
             NewMesh % Elements(NewElCnt) % NodeIndexes(2) = &
                  RefElement % NodeIndexes(2)

             ALLOCATE( NewMesh % Elements(NewElCnt) % BoundaryInfo )
             NewMesh % Elements(NewElCnt) % BoundaryInfo = &
                  RefElement % BoundaryInfo

             NULLIFY( NewMesh % Elements(NewElCnt) % &
                Boundaryinfo % GebhardtFactors )

             CALL SetParents( NewMesh % Elements(NewElCnt), &
                  NewMesh, Children, Edge )

             Children(j,2) = NewElCnt
          ELSE
!
!            New element 1:
!            --------------
             NewElCnt = NewElCnt + 1
             NewMesh % Elements(NewElCnt) = RefElement
             NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
             CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
             NewMesh % Elements(NewElCnt) % NodeIndexes = &
                  RefElement % NodeIndexes

             ALLOCATE( NewMesh % Elements(NewElCnt) % BoundaryInfo )

             NewMesh % Elements(NewElCnt) % BoundaryInfo = &
                  RefElement % BoundaryInfo

             NULLIFY( NewMesh % Elements(NewElCnt) % &
                Boundaryinfo % GebhardtFactors )
            
             CALL SetParents( NewMesh % Elements(NewElCnt), &
                  NewMesh, Children, Edge )

             Children(j,1) = NewElCnt
          END IF
       ELSE
!
!         New element 1, this is point element:
!         -------------------------------------
          NewElCnt = NewElCnt + 1
          RefElement => RefMesh % Elements(j)
          n = RefElement % TYPE % NumberOfNodes

          NewMesh % Elements(NewElCnt) = RefElement
          NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
          CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
          NewMesh % Elements(NewElCnt) % NodeIndexes = &
               RefElement % NodeIndexes
               
          ALLOCATE( NewMesh % Elements(NewElCnt) % BoundaryInfo )

          NewMesh % Elements(NewElCnt) % BoundaryInfo = &
               RefElement % BoundaryInfo
 
          NULLIFY( NewMesh % Elements(NewElCnt) % &
             Boundaryinfo % GebhardtFactors )

          NULLIFY( NewMesh % Elements(NewElCnt) % BoundaryInfo % Left )
          NULLIFY( NewMesh % Elements(NewElCnt) % BoundaryInfo % Right )

          Children(j,1) = NewElCnt
       END IF
    END DO

    NewMesh % MaxBDOFs = RefMesh % MaxBDOFs
    DO i = 1,NewMesh % NumberOfBulkElements+NewMesh % NumberOfBoundaryElements
      RefElement => NewMesh % Elements(i)
      NULLIFY( RefElement % PDefs )
      NULLIFY( RefElement % DGIndexes )
      NULLIFY( RefElement % EdgeIndexes )
      NULLIFY( RefElement % FaceIndexes )
      NULLIFY( RefElement % BubbleIndexes )
      IF ( RefElement % BDOFs > 0 ) THEN
        ALLOCATE( RefElement % BubbleIndexes(RefElement % BDOFs) )
        DO j=1,RefElement % BDOFs
          RefElement % BubbleIndexes(j) = NewMesh % MaxBDOFs*(i-1)+j
        END DO
      END IF

      IF ( ASSOCIATED(RefElement % PDefs) ) THEN
        PD => RefElement % PDefs
        CALL AllocatePDefinitions(RefElement)
        RefElement % PDefs = PD
      END IF
    END DO

!
!   Update Gebhardt factors, if present and the current solver
!   is a heat equation solver:
!   ------------------------------------------------------------
 !   IF ( ListGetString( Solver % Values, 'Equation' ) == 'heat equation' ) &
 !        CALL UpdateGebhardtFactors( RefMesh, NewMesh, Children )

    WRITE( Message, * ) 'Bndry element tables generation time (cpu-secs): ',CPUTime()-t
    CALL Info( 'SplitOneLevel', Message, Level=6 )

    DEALLOCATE( EdgeSplitted, MarkedOrder, Children )
!------------------------------------------------------------------------------
  END FUNCTION SplitOneLevel
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION RGBRefinement(  EdgeSplitted,RefMesh ) RESULT(MarkedEdges)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    LOGICAL :: EdgeSplitted(:)
    INTEGER :: MarkedEdges
    TYPE(Mesh_t), POINTER :: RefMesh
!------------------------------------------------------------------------------
    LOGICAL :: MarkedEdgesFound
    INTEGER :: i,j,EdgeNumber,HangingNodes,RGBIterations,Node1,Node2,&
         LongestEdge
    REAL(KIND=dp) :: x1,y1,x2,y2,EdgeLength,MaxLength
!------------------------------------------------------------------------------
    EdgeSplitted = .FALSE.

!   Mark all three edges of the marked elements (RED refinement):
!   -------------------------------------------------------------
!     DO i = 1,RefMesh % NumberOfBulkElements
!        IF ( RefMesh % Elements(i) % Splitted > 0 ) THEN
!           DO j = 1,3
!              EdgeNumber = RefMesh % Elements(i) % EdgeIndexes(j)
!              EdgeSplitted( EdgeNumber ) = .TRUE.
!           END DO
!        END IF
!     END DO

!   Mark the longest edges of the marked elements (GREEN refinement):
!   -----------------------------------------------------------------
    DO i = 1,RefMesh % NumberOfBulkElements
       IF ( RefMesh % Elements(i) % Splitted > 0 ) THEN
          MaxLength   = 0.0D0
          LongestEdge = 0
          DO j = 1,3
             EdgeNumber = RefMesh % Elements(i) % EdgeIndexes(j)
             Node1 = RefMesh % Edges( EdgeNumber ) % NodeIndexes(1)
             Node2 = RefMesh % Edges( EdgeNumber ) % NodeIndexes(2)
             x1 = RefMesh % Nodes % x( Node1 )
             x2 = RefMesh % Nodes % x( Node2 )
             y1 = RefMesh % Nodes % y( Node1 )
             y2 = RefMesh % Nodes % y( Node2 )
             EdgeLength = SQRT((x2-x1)**2+(y2-y1)**2)
             IF (EdgeLength >= MaxLength) THEN
                MaxLength = EdgeLength
                LongestEdge = EdgeNumber
             END IF
          END DO
          EdgeSplitted( LongestEdge ) = .TRUE.
       END IF
    END DO

    MarkedEdges = 0
    DO i = 1,RefMesh % NumberOfEdges
       IF ( EdgeSplitted(i) ) THEN
          MarkedEdges = MarkedEdges + 1
       END IF
    END DO
!   PRINT*,MarkedEdges,' marked edges'

!   Mark longest edges until we have a RGB-refinement:
!   --------------------------------------------------
    RGBiterations = 0
    DO WHILE( .TRUE. )
       HangingNodes = 0
       RGBiterations = RGBiterations+1
       DO i = 1,RefMesh % NumberOfBulkElements
            
!         Check for marked edges and find the longest edge:
!         -------------------------------------------------
          MarkedEdgesFound = .FALSE.
          LongestEdge      = 0
          MaxLength        = 0.0d0
          DO j = 1,3
             EdgeNumber = RefMesh % Elements(i) % EdgeIndexes(j)
             MarkedEdgesFound = MarkedEdgesFound.OR.EdgeSplitted(EdgeNumber)
             Node1 = RefMesh % Edges(EdgeNumber) % NodeIndexes(1)
             Node2 = RefMesh % Edges(EdgeNumber) % NodeIndexes(2)
             x1 = RefMesh % Nodes % x( Node1 )
             x2 = RefMesh % Nodes % x( Node2 )
             y1 = RefMesh % Nodes % y( Node1 )
             y2 = RefMesh % Nodes % y( Node2 )
             EdgeLength = SQRT((x2-x1)**2+(y2-y1)**2)
             IF (EdgeLength >= MaxLength) THEN
                MaxLength = EdgeLength
                LongestEdge = EdgeNumber
             END IF
          END DO
          
!         If there are marked edges, the longest edge must be one of them:
!         ----------------------------------------------------------------
          IF ( MarkedEdgesFound.AND.(.NOT.EdgeSplitted(LongestEdge)) ) THEN
             HangingNodes = HangingNodes + 1
             EdgeSplitted( LongestEdge ) = .TRUE.
          END IF
       END DO

       IF( HangingNodes > 0) THEN
          WRITE( Message, * ) 'RGB ',RGBiterations,' : ',HangingNodes,' new nodes'
          CALL Info( 'RGBRefinement', Message, Level=6 )
          MarkedEdges = MarkedEdges + HangingNodes
       ELSE
          EXIT
       END IF
    END DO
!------------------------------------------------------------------------------
  END FUNCTION RGBRefinement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Find the parent elements to the splitted boundary element
! among the children of the original parent element:
! ---------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE SetParents( Element, Mesh, Children, Edge )
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    TYPE(Element_t), POINTER :: Edge

    INTEGER :: Children(:,:)
    TYPE(Mesh_t), POINTER :: Mesh

    INTEGER j,k,l,n,i0,j0,k0

    TYPE(Element_t), POINTER :: Child

    n = Element % TYPE % NumberOfNodes

    k = Edge % BoundaryInfo % Left % ElementIndex
    NULLIFY( Child )
    DO l=1,4
       IF ( Children(k,l)>0 ) THEN
          Child => Mesh % Elements( Children(k,l) )
          i0 = 0
          DO j0=1,n
             DO k0=1,Child % TYPE % NumberOfNodes
                IF ( Child % NodeIndexes(k0) == Element % NodeIndexes(j0) ) THEN
                   i0 = i0 + 1 
                   EXIT
                END IF
             END DO
          END DO
          IF ( i0 == n ) EXIT
       END IF
    END DO

    IF ( l > 4 ) STOP 'Adaptive: parent 1 not found'
        
    Element % BoundaryInfo % Left  => Child
    NULLIFY( Element % BoundaryInfo % Right )
        
    NULLIFY( Child )
    IF ( ASSOCIATED(Edge % BoundaryInfo % Right) ) THEN
       k = Edge % BoundaryInfo % Right % ElementIndex
       DO l=1,4
          IF ( Children(k,l)>0 ) THEN
             Child => Mesh % Elements( Children(k,l) )
             i0 = 0
             DO j0=1,n
                DO k0=1,Child % TYPE % NumberOfNodes
                   IF ( Child % NodeIndexes(k0) == Element % NodeIndexes(j0) ) THEN
                      i0 = i0 + 1 
                      EXIT
                   END IF
                END DO
             END DO
             IF ( i0 == n ) EXIT
          END IF
       END DO
           
       Element % BoundaryInfo % Right => Child
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SetParents
!------------------------------------------------------------------------------


      END SUBROUTINE hRefinementSolver



