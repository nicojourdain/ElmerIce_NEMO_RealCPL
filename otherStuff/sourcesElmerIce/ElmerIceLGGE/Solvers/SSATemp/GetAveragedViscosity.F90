      !! Potential Pb with functions in Mat properties as by defaul Model%Mesh is mesh2D?
      !! e.g. enhencement factor function of z??? 
      !! Constant Temp??
      !!=> ok si on met CurrentModel % Mesh and Variables => Mesh3D equivalents temporairement
      !! mais du coup Pb si on fait dependre de variables qui vivent sur le maillage 2D!!
      
      SUBROUTINE GetAveragedViscosity(Model,Solver,dt,TransientSimulation )
      USE DefUtils
      IMPLICIT NONE
      !------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
      !--------------------------------------------------------------------------
      Type(Mesh_t),POINTER :: Mesh2D,Mesh3D
      TYPE(Element_t), POINTER :: Element
      TYPE(ValueList_t), POINTER :: SolverParams
      TYPE(ValueList_t), POINTER :: Material
      Type(Variable_t),POINTER :: MuSol

      REAL(KIND=dp) :: zi,ziup
      REAL(KIND=dp) :: Intval
      REAL(KIND=dp),POINTER :: zcoord(:)

      INTEGER :: i,t
      INTEGER :: Active,n
      INTEGER :: DIM
      INTEGER, POINTER :: NodeIndexes(:)
      INTEGER :: Up

      LOGICAL :: Found
      LOGICAL,ALLOCATABLE :: VisistedNode(:)

      CHARACTER(LEN=MAX_NAME_LEN) :: MeshName,DefaultName="3DGrid"
      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='GetAveragedViscosity'

      !!!! Variables for DetectExtrudedStructure
      TYPE(Solver_t), POINTER :: PSolver
      TYPE(Variable_t), POINTER :: Var
      INTEGER, POINTER :: TopPointer(:),BotPointer(:),Upointer(:),DownPointer(:)

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
      CurrentModel % Mesh => Mesh3D
      CurrentModel % Variables => Mesh3D % Variables

      PSolver => Solver
      CALL DetectExtrudedStructure( Mesh3D , PSolver, Var, &
         TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
         UpNodePointer = Upointer , DownNodePointer = DownPointer)

      DIM = Mesh3D%MeshDim
      IF (DIM.EQ.2) THEN
         zcoord=> Mesh3D%Nodes%y
      ELSE IF (DIM.EQ.3) THEN 
         zcoord=> Mesh3D%Nodes%z
      ELSE
        CALL FATAL(SolverName,'Invalid mesh dimension')
      ENDIF
      
      Mesh2D => Solver % Mesh

      ALLOCATE(VisistedNode(Mesh2D%NumberOfNodes))
      VisistedNode=.FALSE.

      MuSol => VariableGet(Mesh2D%Variables,'Mu',UnFoundFatal=.TRUE.)
      
      Active=GetNOFActive( USolver=Solver )
      Do t=1,Active
         Element => GetActiveElement(t,USolver=Solver)
         NodeIndexes => Element % NodeIndexes
         n  = GetElementNOFNodes(Element)
         Material => GetMaterial(Element)

         Do i=1,n
           IF (VisistedNode(NodeIndexes(i))) CYCLE
           CALL consistency(Mesh2D,Mesh3D,NodeIndexes(i),DIM)
           
           IntVal=0._dp
           Up = NodeIndexes(i)
           DO WHILE (Up /= TopPointer(NodeIndexes(i)))
              zi=zcoord(Up)
              ziup=zcoord(UPointer(Up))
              IntVal=IntVal+0.5*(visco(Up)+visco(UPointer(Up)))*(ziup-zi)
              Up = Upointer(Up)
           END DO
           !! Averaged value
           IntVal=IntVal/(zcoord(TopPointer(NodeIndexes(i)))-zcoord(NodeIndexes(i)))
           
           MuSol%Values(MuSol%Perm(NodeIndexes(i)))=IntVal
           VisistedNode(NodeIndexes(i))=.TRUE.
         End Do
      End do

      DEALLOCATE(VisistedNode)

      CurrentModel % Mesh => Mesh2D
      CurrentModel % Variables => Mesh2D % Variables
      RETURN

      CONTAINS
       ! check that node numbering is consistent between Mesh2D and
       ! Mesh3D so that we will be able to exchange variables
       SUBROUTINE consistency(Mesh2D,Mesh3D,n,DIM) 
       IMPLICIT NONE
       TYPE(Mesh_t),POINTER :: Mesh2D,Mesh3D
       INTEGER :: n,DIM
       LOGICAL :: check

       REAL(KIND=dp) :: x0,y0,x,y
       REAL(KIND=dp) :: dist

         x0=Mesh2D%Nodes%x(n)
         y0=Mesh2D%Nodes%y(n)

         x=Mesh3D%Nodes%x(n)
         y=Mesh3D%Nodes%y(n)

         dist=(x-x0)*(x-x0)
         IF (DIM.EQ.3) dist=dist+(y-y0)*(y-y0)
         dist=sqrt(dist)

         IF (dist.GT.1.0e-6) THEN
            PRINT *,ParEnv%MyPe,dist
            CALL FATAL(SolverName,'Mesh consitency check failed')
         END IF

       END SUBROUTINE  consistency
   
       FUNCTION Visco(nodenumber) RESULT(eta0)
    !------------------------------------------
     IMPLICIT NONE
    !------------------------------------------
     INTEGER :: nodenumber
     REAL(KIND=dp)  :: eta0
    !------------------------------------------
     LOGICAL :: GotIt

     CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityFlag, TemperatureName
     REAL(KIND=dp) :: c1,Temp, Tlimit, &
          A1, A2, Q1, Q2, R, EhF, ArrheniusFactor

     ! Temperature is needed for thermal models
     TYPE(Variable_t), POINTER :: TempSol 
     REAL(KIND=dp), POINTER :: Temperature(:)
     INTEGER, POINTER :: TempPerm(:)

     LOGICAL :: SetArrheniusFactor=.FALSE.
     
     eta0 = 0.0_dp

     ViscosityFlag = ListGetString( Material,'Viscosity Model', UnFoundFatal=.TRUE.)

     SELECT CASE( ViscosityFlag )

     CASE('glen')
        c1 = ListGetConstReal( Material, 'Glen Exponent', UnFoundFatal=.TRUE. ) ! this is the real exponent, n, not 1/n
        
        SetArrheniusFactor = GetLogical(Material, 'Set Arrhenius Factor', GotIt)
        IF ( (.NOT.GotIt) .OR. .NOT.(SetArrheniusFactor)) THEN
           Temp = ListGetRealAtNode(Material, 'Constant Temperature', nodenumber, GotIt) !we are happy as is
           IF(.NOT.GotIt) THEN !we have to find a temperature field
              TemperatureName = GetString(Material, 'Temperature Field Variable', GotIt)
              IF (.NOT.GotIt) WRITE(TemperatureName,'(A)') 'Temperature'
              TempSol => VariableGet( Mesh3D % Variables,TRIM(TemperatureName),UnFoundFatal=.TRUE.)
              TempPerm    => TempSol % Perm
              Temperature => TempSol % Values   
              Temp =  Temperature(TempPerm(nodenumber))
           END IF
        
           R = ListGetConstReal( Model % Constants,'Gas Constant',UnFoundFatal=.TRUE.)
           Tlimit = ListGetConstReal(Material, 'Limit Temperature', UnFoundFatal=.TRUE.)
           A1 = ListGetConstReal(Material, 'Rate Factor 1', UnFoundFatal=.TRUE.)
           A2 = ListGetConstReal(Material, 'Rate Factor 2', UnFoundFatal=.TRUE.)
           Q1 = ListGetConstReal(Material, 'Activation Energy 1', UnFoundFatal=.TRUE.)
           Q2 = ListGetConstReal(Material, 'Activation Energy 2', UnFoundFatal=.TRUE.)
        
           IF (Temp.LE. Tlimit) THEN
              ArrheniusFactor = A1 * EXP( -Q1/(R * (273.15 + Temp)))
           ELSE IF((Tlimit<Temp) .AND. (Temp .LE. 0.0_dp)) THEN
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15 + Temp)))
           ELSE
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15)))
              CALL INFO('IntegrateViscositySSA','Positive Temperature detected in Glen - limiting to zero!', Level = 5)
           END IF
        ELSE
           ArrheniusFactor = ListGetConstReal(Material,'Arrhenius Factor',UnFoundFatal=.TRUE.)
        END IF

        !PRINT *,'In SUB',nodenumber
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
      END SUBROUTINE GetAveragedViscosity
      
