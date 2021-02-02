!Compute the Cost function of the Arhtern/Gudmundsson inverse Problem
!      as Sum_Surface (vn-vd).(sigma_n-sigma_d).n
!      with a regularization as Sum_bedrock 0.5 Lambda (dBeta/dx)^2
!
!   Serial/Parallel    2D/3D
!
! Need : 
!   - Name of the Cost Variable
!   - Solutions of Neumann and Dirchlet problem
!       (Velocities and Stresses)
!   - Lambda and Beta for regularization
!   - define in the sif Name='surface' and Name='bed' in appropriate BC.
!
! *****************************************************************************
SUBROUTINE CostSolver_Robin_Stresses( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,NeumannSolName,DirichletSolName,CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: NeumannStressName,DirichletStressName,BCName
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName, VarSolName 
  TYPE(Solver_t), POINTER :: ParSolver
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar,VeloSolN,VeloSolD,StressSolD,StressSolN,BetaSol
  TYPE(ValueList_t), POINTER :: BC,SolverParams
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp), POINTER :: VelocityN(:),VelocityD(:),StressD(:),StressN(:),Beta(:)
  INTEGER, POINTER :: VeloNPerm(:),VeloDPerm(:),StressNPerm(:),StressDPerm(:),NodeIndexes(:), BetaPerm(:)
  Logical :: Firsttime=.true.,Found,Parallel,stat
  integer :: i,j,k,l,t,n,NMAX,NActiveNodes,DIM,ierr
  real(kind=dp) :: Unorm,Cost,Cost_surf,Cost_bed,Cost_S,Cost_surf_S,Cost_bed_S,Normal(3),vn(3),vd(3),Sn(3,3),Sd(3,3),sTimesN,Lambda
  real(kind=dp) :: Bu,Bv,u,v,w,s,coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: NodeCost(Model % MaxElementNodes),Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  CHARACTER*10 :: date,temps

  save Firsttime,Parallel,CostFile,DIM,ElementNodes
  save SolverName,NeumannSolName,DirichletSolName,VarSolname,CostSolName
  save NeumannStressName,DirichletStressName,Lambda


  If (Firsttime) then

     DIM = CoordinateSystemDimension()
     WRITE(SolverName, '(A)') 'CostSolver_Robin_Z'

!!!!!!!!!!! get Solver Variables
  SolverParams => GetSolverParams()

  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile
    OPEN (12, FILE=CostFile)
                    CALL DATE_AND_TIME(date,temps)
                    write(12,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
    CLOSE(12)


   NeumannSolName =  GetString( SolverParams,'Neumann Solution Name', Found)
         IF(.NOT.Found) THEN
                 CALL WARN(SolverName,'Keyword >Neumann Solution Name< not found  in section >Equation<')
                 CALL WARN(SolverName,'Taking default value >Flow Solution<')
                 WRITE(NeumannSolName,'(A)') 'Flow Solution'
         END IF
   NeumannStressName =  GetString( SolverParams,'Neumann Stress Name', Found)
         IF(.NOT.Found) THEN
                 CALL WARN(SolverName,'Keyword >Neumann Stress Name< not found  in section >Equation<')
                 CALL WARN(SolverName,'Taking default value >StressN<')
                 WRITE(NeumannStressName,'(A)') 'StressN'
         END IF

   DirichletSolName =  GetString( SolverParams,'Dirichlet Solution Name', Found)
       IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Dirichlet Solution Name< not found  in section >Equation<')
           CALL WARN(SolverName,'Taking default value >Flow Solution<')
           WRITE(NeumannSolName,'(A)') 'Flow Solution'
       End if
   DirichletStressName =  GetString( SolverParams,'Dirichlet Stress Name', Found)
       IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Dirichlet Stress Name< not found  in section >Equation<')
           CALL WARN(SolverName,'Taking default value >StressD<')
           WRITE(DirichletStressName,'(A)') 'StressD'
       End if
   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
       IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Equation<')
           CALL WARN(SolverName,'Taking default value Lambda=0.0')
           Lambda = 0.0
       End if
   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF
   VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >Beta<')
                    WRITE(VarSolName,'(A)') 'Beta'
          END IF


!!!!!!! Check for parallel run 
    Parallel = .FALSE.
      IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
        IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
          Parallel = .TRUE.
        END IF
      END IF
  
  !!! End of First visit
    Firsttime=.false.
  Endif


    VeloSolN => VariableGet( Solver % Mesh % Variables, NeumannSolName)
        IF ( ASSOCIATED( VeloSolN ) ) THEN
            VelocityN => VeloSolN % Values
            VeloNPerm => VeloSolN % Perm
        ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',NeumannSolName,'< found'
            CALL FATAL(SolverName,Message)
        END IF            
    StressSolN => VariableGet( Solver % Mesh % Variables, NeumannStressName)
        IF ( ASSOCIATED( StressSolN ) ) THEN
            StressN => StressSolN % Values
            StressNPerm => StressSolN % Perm
        ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',NeumannStressName,'< found'
            CALL FATAL(SolverName,Message)
        END IF            
    VeloSolD => VariableGet( Solver % Mesh % Variables, DirichletSolName )
        IF ( ASSOCIATED( VeloSolD ) ) THEN
            VelocityD => VeloSolD % Values
            VeloDPerm => VeloSolD % Perm
        ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',DirichletSolName,'< found'
            CALL FATAL(SolverName,Message)
        END IF
    StressSolD => VariableGet( Solver % Mesh % Variables, DirichletStressName )
        IF ( ASSOCIATED( StressSolD ) ) THEN
            StressD => StressSolD % Values
            StressDPerm => StressSolD % Perm
        ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',DirichletStressName,'< found'
            CALL FATAL(SolverName,Message)
        END IF      

    BetaSol => VariableGet( Solver % Mesh % Variables, VarSolName )
        IF ( ASSOCIATED( BetaSol ) ) THEN
            Beta => BetaSol % Values
            BetaPerm => BetaSol % Perm
        ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >', VarSolName, '< found'
            CALL FATAL(SolverName,Message)
        END IF            


    Cost=0._dp
    Cost_surf=0.0_dp
    Cost_bed=0.0_dp
    DO t=1,Solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)

      BC => GetBC()
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE

      BCName =  ListGetString( BC,'Name', Found)
      IF((BCName /= 'surface').AND.(BCName /= 'bed')) CYCLE


      CALL GetElementNodes( ElementNodes )
      n = GetElementNOFNodes()
      NodeIndexes => Element % NodeIndexes

      NodeCost=0.0_dp
      Do i=1,n
           j =  NodeIndexes( i )
           
       IF (BCName == 'surface') THEN
           Bu = Element % Type % NodeU(i)
           IF ( Element % Type % Dimension > 1 ) THEN
              Bv = Element % Type % NodeV(i)
           ELSE
              Bv = 0.0D0
           END IF
           Normal = NormalVector(Element, ElementNodes, Bu, Bv, .TRUE.)
           s = SQRT( SUM( Normal(1:DIM)**2 ) )
           Normal = Normal/s

           vn=0.0_dp
           vd=0.0_dp
           vn(1:DIM)=VelocityN((DIM+1)*(VeloNPerm(j)-1)+1:DIM)
           vd(1:DIM)=VelocityD((DIM+1)*(VeloDPerm(j)-1)+1:DIM)

           sd=0.0_dp
           sd(1,1)=StressD(2*DIM*(StressDPerm(j)-1)+1)-VelocityD((DIM+1)*(VeloDPerm(j)-1)+DIM+1)
           sd(2,2)=StressD(2*DIM*(StressDPerm(j)-1)+2)-VelocityD((DIM+1)*(VeloDPerm(j)-1)+DIM+1)
           sd(3,3)=StressD(2*DIM*(StressDPerm(j)-1)+3)-VelocityD((DIM+1)*(VeloDPerm(j)-1)+DIM+1)
           sd(1,2)=StressD(2*DIM*(StressDPerm(j)-1)+4)
           sd(2,1)=sd(1,2)

           sn=0.0_dp
           sn(1,1)=StressN(2*DIM*(StressNPerm(j)-1)+1)-VelocityN((DIM+1)*(VeloNPerm(j)-1)+DIM+1)
           sn(2,2)=StressN(2*DIM*(StressNPerm(j)-1)+2)-VelocityN((DIM+1)*(VeloNPerm(j)-1)+DIM+1)
           sn(3,3)=StressN(2*DIM*(StressNPerm(j)-1)+3)-VelocityN((DIM+1)*(VeloNPerm(j)-1)+DIM+1)
           sn(1,2)=StressN(2*DIM*(StressNPerm(j)-1)+4)
           sn(2,1)=sn(1,2)

          if (DIM.eq.3) then
              sd(2,3)=StressD(2*DIM*(StressDPerm(j)-1)+5)
              sd(3,2)=sd(2,3)
              sd(1,3)=StressD(2*DIM*(StressDPerm(j)-1)+6)
              sd(3,1)=sd(1,3)

              sn(2,3)=StressN(2*DIM*(StressNPerm(j)-1)+5)
              sn(3,2)=sn(2,3)
              sn(1,3)=StressN(2*DIM*(StressNPerm(j)-1)+6)
              sn(3,1)=sd(1,3)

           end if
               
           Do k=1,3
               sTimesN=0.0_dp
               Do l=1,3
                 sTimesN=sTimesN+(sn(k,l)-sd(k,l))*Normal(l)
               End do
               NodeCost(i)=NodeCost(i)+(vn(k)-vd(k))*sTimesN
          End do

          Else IF (BCName == 'bed') Then
              NodeCost(i)=Beta(BetaPerm(j))
          End if
         End do

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------
        IntegStuff = GaussPoints( Element )

        DO i=1,IntegStuff % n
          U = IntegStuff % u(i)
          V = IntegStuff % v(i)
          W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )

          x = SUM( ElementNodes % x(1:n) * Basis(1:n) )
          s = 1.0d0

            IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
              s = 2.0d0 * PI * x 
            END IF
            s = s * SqrtElementMetric * IntegStuff % s(i)
          
           IF (BCName == 'surface') Then   
            coeff = SUM(NodeCost(1:n)  * Basis(1:n))
            Cost_surf=Cost_surf+coeff*s
           else IF (BCName == 'bed') Then
            coeff = SUM(NodeCost(1:n) * dBasisdx(1:n,1))
            coeff = coeff*coeff
            IF (DIM.eq.3) then
                    coeff=coeff+ & 
                    SUM(NodeCost(1:n)*dBasisdx(1:n,2))*SUM(NodeCost(1:n) * dBasisdx(1:n,2))
            END IF
            !coeff = 0.5*Lambda * coeff 
            Cost_bed=Cost_bed+coeff*s
           else 
            coeff = 0.0
           End if

            !Cost=Cost+coeff*s
        End do
    End do

   Cost=Cost_surf+0.5*Lambda*Cost_bed

   TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

   IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  
           CALL MPI_ALLREDUCE(Cost_surf,Cost_surf_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  
           CALL MPI_ALLREDUCE(Cost_bed,Cost_bed_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
                 CostVar % Values(1)=Cost_S
          END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (12, FILE=CostFile,POSITION='APPEND')
                 write(12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost_S,Cost_surf_S,Cost_bed_S
                 CLOSE(12)
         End if
   ELSE
              CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
              IF (ASSOCIATED(CostVar)) THEN
                   CostVar % Values(1)=Cost
              END IF
               OPEN (10, FILE=CostFile,POSITION='APPEND')
               write(10,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost,Cost_surf,Cost_bed
               close(10)
   END IF

   Return

!------------------------------------------------------------------------------
END SUBROUTINE CostSolver_Robin_Stresses
!------------------------------------------------------------------------------
! *****************************************************************************
