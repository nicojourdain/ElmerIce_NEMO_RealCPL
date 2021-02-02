      SUBROUTINE Activate3DGrid_init( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      IMPLICIT NONE
      !------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
      !--------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: MeshName,DefaultName="3DGrid"
      TYPE(ValueList_t), POINTER :: SolverParams
      LOGICAL :: GotIt,Found
      REAL(KIND=dp) :: q
      INTEGER :: ExtrudeLayers
      TYPE(Mesh_t),POINTER :: ExtrudedMesh
      CHARACTER(LEN=1024) :: Path
      LOGICAL :: SAVE_MESH=.FALSE.

      SolverParams => Solver % Values

      ExtrudeLayers=ListGetInteger(SolverParams,'Extruded Mesh Levels',UnFoundFatal=.TRUE.)-1
      q=ListGetConstReal(SolverParams,'Extruded Mesh Ratio',Found)
      IF (Found) THEN
        CALL ListAddConstReal(Model%Simulation,'Extruded Mesh Ratio',q)
      END IF

      ExtrudedMesh => MeshExtrude(CurrentModel % Meshes,ExtrudeLayers-1)
      CALL MeshStabParams(ExtrudedMesh)
      ExtrudedMesh % OutputActive=.TRUE.

      MeshName=ListGetString(SolverParams,'Extruded Mesh Name',Found)
      IF (.NOT.Found) MeshName=DefaultName
      ExtrudedMesh % Name = trim(MeshName)

      Path = TRIM(ExtrudedMesh % Name)
      CALL MakeDirectory( TRIM(path) // CHAR(0) )

      CurrentModel % Meshes % Next => ExtrudedMesh

      Solver % Mesh => ExtrudedMesh

      SAVE_MESH=ListGetLogical(SolverParams,'Save Extruded Mesh',Found)
      IF (SAVE_MESH) THEN 
        IF (ParEnv % Pes > 1) THEN
           CALL WriteMeshToDisk2(Model, ExtrudedMesh , Path ,ParEnv%MyPE)
        ELSE
           CALL WriteMeshToDisk2(Model, ExtrudedMesh , Path)
        END IF
      ENDIF

      END SUBROUTINE Activate3DGrid_init
      
      SUBROUTINE Activate3DGrid(Model,Solver,dt,TransientSimulation )
      USE DefUtils
      IMPLICIT NONE
      !------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
      !--------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: Name,VarName
      TYPE(ValueList_t), POINTER :: SolverParams
      LOGICAL :: GotIt
      Type(Mesh_t),POINTER :: Mesh2D,Mesh3D
      Type(Variable_t),POINTER :: MuSol,VarSol,Var3DSol
      INTEGER :: i,t
      INTEGER :: n
      LOGICAL :: stat
      REAL(KIND=dp) :: x,y,x0,y0,dist
      TYPE(Element_t), POINTER :: Element
      TYPE(GaussIntegrationPoints_t) :: IP
      REAL(KIND=dp) :: u,v,w
      REAL(KIND=dp) :: detJ
      REAL(KIND=dp),ALLOCATABLE :: Basis(:),dBasisdx(:,:)
      TYPE(Nodes_t),SAVE :: Nodes
      !!!! Variables for DetectExtrudedStructure
      TYPE(Solver_t), POINTER :: PSolver
      TYPE(Variable_t), POINTER :: Var
      INTEGER, POINTER :: TopPointer(:),BotPointer(:),Upointer(:),DownPointer(:)
      INTEGER :: Up
      LOGICAL :: GotVar
      INTEGER :: j
#if 1
      SolverParams => Solver % Values

      Mesh2D => CurrentModel % Meshes
      Mesh3D => Solver % Mesh

      PSolver => Solver
      CALL DetectExtrudedStructure( Mesh3D , PSolver, Var, &
         TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
         UpNodePointer = Upointer , DownNodePointer = DownPointer)
      
      j=1
      DO WHILE(.TRUE.) 
        WRITE (Name,'(A,I0)') 'Exported Variable ',j
        VarName = ListGetString(SolverParams,TRIM(Name),GotVar)
        IF (.NOT.GotVar) EXIT
        Var3DSol => VariableGet(Mesh3D%Variables,TRIM(VarName),UnFoundFatal=.TRUE.)
        WRITE (Name,'(A,I0)') 'Target Variable ',j
        VarName = ListGetString(SolverParams,TRIM(Name),UnFoundFatal=.TRUE.)
        VarSol => VariableGet(Mesh2D%Variables,TRIM(VarName),UnFoundFatal=.TRUE.)


        Do i=1,Mesh2D%NumberOfNodes
         Up = i
         Var3DSol%Values(Var3DSol%Perm(Up))=VarSol%Values(VarSol%Perm(i))
         DO WHILE (Up /= TopPointer(i))
              Up = Upointer(Up)
              Var3DSol%Values(Var3DSol%Perm(Up))=VarSol%Values(VarSol%Perm(i))
          END DO

        End do

        j=j+1
      END DO
#else
      SolverParams => Solver % Values

      PRINT *,ParEnv%MyPe,'HERE WE ARE!!'

      Mesh2D => CurrentModel % Meshes
      Mesh3D => Solver % Mesh

      ALLOCATE(Basis(Mesh3D%MaxElementNodes),dBasisdx(Mesh3D%MaxElementNodes,3))

      PRINT *,Solver%Mesh%NumberOfNodes,Mesh2D%NumberOfNodes

      MuSol => VariableGet(Mesh2D%Variables,'Mu',UnFoundFatal=.TRUE.)

      Do i=1,Mesh2D%NumberOfNodes
         x0=Mesh2D%Nodes%x(i)
         y0=Mesh2D%Nodes%y(i)
         x=Mesh3D%Nodes%x(i)
         y=Mesh3D%Nodes%y(i)
         dist=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))
         IF (dist.GT.1.0e-6) THEN
            PRINT *,ParEnv%MyPe,dist
            CALL FATAL('TEST3D','Distance check failed')
         END IF
         MuSol % Values(MuSol %Perm(i))=0.5
      End DO

      t=1
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      IP = GaussPoints( Element)
      CALL GetElementNodes(Nodes,Element,Solver,Mesh3D)

      PRINT *,'Element: ',Element%TYPE%ElementCode
      PRINT *,'nIP: ',IP%n
      PRINT *,'i, U, V, W '
      Do i=1,IP%n
        u=IP % U(i)
        v=IP % V(i)
        w=IP % W(i)
        stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis,dBasisdx )
        PRINT *,'-------------------------------------------'
        PRINT *,i,IP % U(i),IP % V(i),IP % W(i)
        PRINT *,Basis(1:n)
        PRINT *,dBasisdx(1:n,1)
        PRINT *,dBasisdx(1:n,2)
        PRINT *,dBasisdx(1:n,3)
      End do

      DEALLOCATE(Basis,dBasisdx)
#endif

      END SUBROUTINE Activate3DGrid
      
