! Calcul du gradient de la fonction cout du probleme inverse  de Robin 
!  (Arthern & Gudmundsson, J. Glaciol., 2010)
!   par rapport a un enhencement de la viscosité (i.e mu=E mu_0)
!   en supposant une perturbation independante de Z
!
!   !! Le maillage doit etre extrudé verticalement, i.e. les points alignes selon la verticale !!
!
!  !!!!!!  Regularisation ajoutee comme Lambda * int_bed (DMu/Dx)^2  !!!!!!!
!
!  Serial/Parallel   and 2D/3D
!
! Paramètres dans le sif:
!    Dans la section du Solver:
!     - Active Coordinate = integer ;  la direction de l'extrusion
!     - Neumann Solution Name = String ; le nom de la variable du probleme de Neumann ("Flow Solution" par defaut)
!     - Dirichlet Solution Name = String ; le nom de la variable du probleme de Dirichlet ("VeloD" par defaut)
!     - mu0 = String ; le nom de la variable mu_0 ("mu0" par defaut)
!     - Optimized Variable Name = String ; le nom de la variable a optimiser ("E" par default)
!     - Gradient Variable Name = String ; le nom de la variable dans lequel est mis le gradient ("DJDE" par default)
!     - SquareFormulation = Logical ; True si la viscosite definie comme alpha^2 et optimisation sur alpha pour assurer une viscosite positive
!     - Lambda = Real ; valeur du parametre de regularisation
!
!    Solver a executer dans le body
! 
!   Dans la section Material:
!     Si SquareFormulation = False:
!          Viscosity = Variable E , mu0
!              Real MATC "tx(0)*tx(1)"
!     Si SquareFormulation = True:
!         Viscosity = Variable E , mu0
!           Real MATC "tx(0)*tx(0)*tx(1)"
!
!    Dans le section BC:
!     Pour la condition limite inferieure:
!       Name = String "bed"  ! Pour connaitre la BC inferieure et se reperer par rapport a l'Extrusion
!
!  Solver a executer dans le body; en conjonction avec :
!       -CostSolver_Robin.f90: pour le calcul de la fonction cout;
!       -Optimize_m1qn3[Serial/Parallel]_Z.f90: pour l'optimisation
!

! *****************************************************************************
SUBROUTINE DJDE_Robin_Z( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: SolverParams,Material,BC ! pointeurs pour avoir acces aux infos du sif

!!!!  Variables utiles pour les elements et les fonctions de base
  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  real(kind=dp),allocatable :: Basis(:),dBasisdx(:,:)
  real(kind=dp) :: u,v,w,SqrtElementMetric
  INTEGER, POINTER :: NodeIndexes(:)

!!!!! variables Elmer
  TYPE(Variable_t), POINTER :: Variable, GradVariable, VeloSolN,VeloSolD,Mu0Sol
  REAL(KIND=dp), POINTER :: Values(:),GradValues(:),VelocityN(:),VelocityD(:),Mu0(:)
  INTEGER, POINTER :: Perm(:), VeloNPerm(:),VeloDPerm(:),GradPerm(:),Mu0Perm(:)

!!!! Variables pour DetectExtrudedStructure
  TYPE(Solver_t), POINTER :: PSolver  
  TYPE(Variable_t), POINTER :: Var
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),Upointer(:),DownPointer(:)

!!!! Autres variables
  real(kind=dp),allocatable :: VisitedNode(:),db(:)
  real(kind=dp),allocatable :: NodeDJ(:),NodalVeloN(:,:),NodalVeloD(:,:),NodalMu0(:)
  real(kind=dp),allocatable :: m(:),cs(:)
  real(kind=dp) :: Lambda
  real(kind=dp) :: vn(3),vd(3),LGradN(3,3),LGradD(3,3),SRD(3,3),SRN(3,3)
  real(kind=dp) :: IPGrad
  real(kind=dp) :: Sgrad,SecInv


  integer :: i,j,k,t,n,NMAX,NActiveNodes,DIM
  integer :: up

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  CHARACTER(LEN=MAX_NAME_LEN) :: VarSolName,GradSolName,NeumannSolName,DirichletSolName,Mu0SolName
  CHARACTER(LEN=MAX_NAME_LEN) :: BCName
  
  Logical ::  Firsttime=.true.
  Logical ::  Found,stat
  logical ::  SquareFormulation


  save Firsttime,DIM
  save ElementNodes
  save SolverName
  save NeumannSolName,DirichletSolName,VarSolName,GradSolName,Mu0SolName
  SAVE Lambda
  save SquareFormulation
  save VisitedNode,db,NodeDJ,Basis,dBasisdx,NodalVeloN,NodalVeloD,NodalMu0,m,cs
  SAVE TopPointer,BotPointer,Upointer,DownPointer

  If (Firsttime) then
      WRITE(SolverName, '(A)') 'DJDE_Robin_Z'
      DIM = CoordinateSystemDimension()

      !!!!  
      PSolver => Solver
      CALL DetectExtrudedStructure(Solver % Mesh , PSolver, Var, &
           TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
           UpNodePointer = Upointer , DownNodePointer = DownPointer)


      NMAX=Solver % Mesh % NumberOfNodes
      allocate(VisitedNode(NMAX),db(NMAX), NodeDJ(Model %  MaxElementNodes), &
               Basis(Model % MaxElementNodes),  &
               dBasisdx(Model % MaxElementNodes,3), &
               NodalVeloN(3,Model % MaxElementNodes),NodalVeloD(3,Model % MaxElementNodes), &
               NodalMu0(Model % MaxElementNodes), &
               m(Model % MaxElementNodes),cs(Model % MaxElementNodes))

!!!!!!!!!!! get Solver Variables
      SolverParams => GetSolverParams()

      NeumannSolName =  GetString( SolverParams,'Neumann Solution Name', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >Neumann Solution Name< not found in section >Solver<')
               CALL WARN(SolverName,'Taking default value >Flow Solution<')
               WRITE(NeumannSolName,'(A)') 'Flow Solution'
          END IF
      DirichletSolName =  GetString( SolverParams,'Dirichlet Solution Name', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >Dirichlet Solution Name< not found in section >Solver<')
               CALL WARN(SolverName,'Taking default value >VeloD<')
               WRITE(DirichletSolName,'(A)') 'VeloD'
          END IF
      Mu0SolName = GetString( SolverParams,'Mu0', Found)
          IF(.NOT.Found) THEN
                  CALL WARN(SolverName,'Keyword >Mu0< not found in section  >Solver<')
                  CALL WARN(SolverName,'Taking default value >Mu0<')
                  WRITE(Mu0SolName,'(A)') 'Mu0'
          ENDIF
      VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
             IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >E<')
                    WRITE(VarSolName,'(A)') 'E'
              END IF
      GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
             IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >DJDE<')
                    WRITE(GradSolName,'(A)') 'DJDE'
             END IF

       Lambda =  GetConstReal( SolverParams,'Lambda', Found)
            IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value Lambda=0.0')
                    Lambda = 0.0_dp
            END IF

       SquareFormulation=GetLogical( SolverParams, 'SquareFormulation', Found)
           IF(.NOT.Found) THEN
                   CALL WARN(SolverName,'Keyword >SquareFormulation< not found  in section >Solver<')
                   CALL WARN(SolverName,'Taking default value >FALSE<')
                   SquareFormulation=.FALSE.
           END IF
  
  !!! End of First visit
    Firsttime=.false.
  Endif

 ! Get variables needed by the Solver

        GradVariable => VariableGet( Solver % Mesh % Variables, GradSolName )
           IF (ASSOCIATED(GradVariable)) THEN
              GradValues => GradVariable % Values
              GradPerm => GradVariable % Perm
           ELSE
              WRITE(Message,'(A,A,A)') 'No variable >',GradSolName,'< found'
             CALL FATAL(SolverName,Message)
           END IF
        Variable => VariableGet( Solver % Mesh % Variables, VarSolName )
           IF (ASSOCIATED(Variable)) THEN
                Values => Variable % Values
                Perm => Variable % Perm
           ELSE
                WRITE(Message,'(A,A,A)') 'No variable >',VarSolName,'< found'
                CALL FATAL(SolverName,Message)
           END IF
        VeloSolN => VariableGet( Solver % Mesh % Variables, NeumannSolName )
           IF ( ASSOCIATED( VeloSolN ) ) THEN
             VelocityN => VeloSolN % Values
             VeloNPerm => VeloSolN % Perm
           ELSE
              WRITE(Message,'(A,A,A)') &
                   'No variable >',NeumannSolName,'< found'
              CALL FATAL(SolverName,Message)              
           END IF
         VeloSolD => VariableGet( Solver % Mesh % Variables, DirichletSolName )
          IF (ASSOCIATED(veloSolD)) THEN
             VelocityD => VeloSolD % Values
             VeloDPerm => VeloSolD % Perm
          ELSE
              WRITE(Message,'(A,A,A)') &
                   'No variable >',DirichletSolName,'< found'
              CALL FATAL(SolverName,Message)              
           END IF

        Mu0Sol=> VariableGet( Solver % Mesh % Variables, Mu0SolName )
           IF (ASSOCIATED(Mu0Sol)) THEN
                   Mu0 => Mu0Sol % Values
                   Mu0Perm => Mu0Sol % Perm
            ELSE
                    WRITE(Message,'(A,A,A)') &
                         'No variable >',Mu0SolName,'< found'
                    CALL FATAL(SolverName,Message)
            ENDIF

    !!!!! calcul les gradient nodaux       
    VisitedNode=0.0_dp
    db=0.0_dp

    DO t=1,Solver % NumberOfActiveElements

          Element => GetActiveElement(t)
          Material => GetMaterial()
          CALL GetElementNodes( ElementNodes )
          n = GetElementNOFNodes()
          NodeIndexes => Element % NodeIndexes

          NodalMu0=0._dp
          NodalMu0(1:n)=Mu0(Mu0Perm(NodeIndexes(1:n)))
          !!! vitesses nodales
          NodalVeloN = 0.0d0
          NodalVeloD = 0.0d0
          DO i=1, dim
             NodalVeloN(i,1:n) = VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(1:n))-1)+i)
             NodalVeloD(i,1:n) = VelocityD((DIM+1)*(VeloDPerm(NodeIndexes(1:n))-1)+i)
          END DO

          !! exposant nodal
          m=ListGetReal(Material, 'Viscosity Exponent', n, NodeIndexes)
          cs=ListGetReal(Material, 'Critical Shear Rate', n, NodeIndexes)

          ! Compute Nodal Value of DJDE=DJDmu*DmuDE
          Do i=1,n
             VisitedNode(NodeIndexes(i))=VisitedNode(NodeIndexes(i))+1.0_dp

             !!! infos sur les fonctions de bases pour calculer les vitesses de deformation
             u=Element % Type % NodeU(i)
             v=Element % Type % NodeV(i)
             w=Element % Type % NodeW(i)
             stat=ElementInfo(Element,ElementNodes,u,v,w, &
                                SqrtElementMetric,Basis,dBasisdx)

             LGradN=0.0_dp
             LGradD=0.0_dp
             LGradN = MATMUL( NodalVeloN(:,1:n), dBasisdx(1:n,:) )
             SRN = 0.5 * ( LGradN + TRANSPOSE(LGradN) )
             LGradD = MATMUL( NodalVeloD(:,1:n), dBasisdx(1:n,:) )
             SRD = 0.5 * ( LGradD + TRANSPOSE(LGradD) )

             !!! DJ/DV_eff = (2*(|e^D|^2-e^N|^2)
             NodeDJ(i)=2._dp*(calcNorm2(SRD)-calcNorm2(SRN))

             !!!! si m != 1; DJ/DV=DJ/DV_eff * DV_eff/DV; V_eff=V*gamma^(m-1)
             if (m(i).ne.1.0_dp) then
                     SecInv=2.0_dp*calcNorm2(SRN) ! le carre du Second invariant de D
                     If (SecInv.lt.cs(i)*cs(i)) SecInv=cs(i)*cs(i)
                     NodeDJ(i)=NodeDJ(i)*(SecInv)**(0.5_dp*(m(i)-1._dp))
             endif

             !!!!! DJ/DE_eff=DJ/DV*DV/DE_eff; V=E_eff*Mu0
             NodeDJ(i)=NodeDJ(i)*NodalMu0(i)

             !!! si E_eff=E**2;
             IF (SquareFormulation) then
                      NodeDJ(i)=NodeDJ(i)*2.0_dp*Values(Perm(NodeIndexes(i)))
              ENDIF
           End do

           ! Compute Integrated Nodal Value of DJDE
           IntegStuff = GaussPoints( Element )
           DO j=1,IntegStuff % n
              U = IntegStuff % u(j)
              V = IntegStuff % v(j)
              W = IntegStuff % w(j)
              stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                             Basis,dBasisdx )
              Do i=1,n
                    IPGrad=NodeDJ(i)*Basis(i)
                    db(NodeIndexes(i)) = db(NodeIndexes(i)) + &
                                   SqrtElementMetric*IntegStuff % s(j)*IPGrad
              End do
            End Do
    End do


    !!!!! Derivee du terme de regularisation int_bed (DValues/Dx)2 db
    If (Lambda /= 0.0) then
      DoBoundaryElements: DO k=1,Solver % Mesh % NumberOfBoundaryElements
                  Element => GetBoundaryElement(k)
                  BC => GetBC()
                  IF ( .NOT. ASSOCIATED(BC) ) CYCLE
                  BCName =  ListGetString( BC,'Name', Found)
                  IF(BCName /= 'bed') CYCLE
                  CALL GetElementNodes( ElementNodes )
                  n = GetElementNOFNodes()
                  NodeIndexes => Element % NodeIndexes

                  IntegStuff = GaussPoints( Element )
                  DO j=1,IntegStuff % n
                     U = IntegStuff % u(j)
                     V = IntegStuff % v(j)
                     W = IntegStuff % w(j)
                     stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                                                Basis,dBasisdx )
                     Do i=1,n
                        IPGrad=Lambda*(SUM(dBasisdx(1:n,1)*Values(Perm(NodeIndexes(1:n))))*dBasisdx(i,1))
                        IF (DIM.eq.3) then
                            IPGrad=IPGrad+Lambda*(SUM(dBasisdx(1:n,2)*Values(Perm(NodeIndexes(1:n))))*dBasisdx(i,2))
                        End if
                        db(NodeIndexes(i)) = db(NodeIndexes(i)) + &
                                             SqrtElementMetric*IntegStuff % s(j)*IPGrad
                     End do
                  End Do
     END DO DoBoundaryElements
   End if

   !!!!!!! Suppose une pertubation independante de Z=>Gradient=Sum_z Gradients(x,z)
   !!!!!!! On calcul Gradient et on l'affecte de maniere homogene suivant z a GradValues
   DoBoundaryElements2: DO k=1,Solver % Mesh % NumberOfBoundaryElements
                  Element => GetBoundaryElement(k)
                  BC => GetBC()
                  IF ( .NOT. ASSOCIATED(BC) ) CYCLE
                  BCName =  ListGetString( BC,'Name', Found)
                  IF(BCName /= 'bed') CYCLE
                  CALL GetElementNodes( ElementNodes )
                  n = GetElementNOFNodes()
                  NodeIndexes => Element % NodeIndexes


                  DoNodes: Do i=1,n
                        up=NodeIndexes(i)
                        Sgrad=db(TopPointer(NodeIndexes(i)))
                        Do while (up.ne.TopPointer(NodeIndexes(i)))
                           Sgrad=Sgrad+db(up)
                           up=Upointer(up)
                        End do

                        up=NodeIndexes(i)
                        Do while (up.ne.TopPointer(NodeIndexes(i)))
                           GradValues(GradPerm(up))=Sgrad
                           up=Upointer(up)
                        End do
                           GradValues(GradPerm(TopPointer(NodeIndexes(i))))=Sgrad

                  END DO DoNodes

   END DO DoBoundaryElements2



   Return

   CONTAINS
          ! Calcul le second invariant de vij comme v2=vij*vji
           function calcNorm2(v) result(v2)
             implicit none
             real(kind=dp) :: v(3,3),v2
             integer :: i,j
             v2=0._dp
             Do i=1,3
               Do j=1,3
                 v2=v2+v(i,j)*v(j,i)
               End do
             End do
           end function calcNorm2
!------------------------------------------------------------------------------
END SUBROUTINE DJDE_Robin_Z
!------------------------------------------------------------------------------


