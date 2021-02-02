!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CalvingSolver.f90
! 7 October : cleaning
! 18 July 2014 : Interpolation of Calving absciss along x-direction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The shape of crevasse depth is given by the level of damage. The deepest nodes where
! the damage exceed a threshol value "Dc" defines a damage enveloppe on which
! the LEFM calculation are done (see below)
!
! Calving Criterion based on LEFM : If stress intensity factor KI > ice toughness "KIc" 
! at depth given by damage shape AND KI > arrest criterion "KIa" at sea
! level, then, fracture occurs along a vertical line.
!
! Horizontal interpolation : When a nodes validates above criteria, they are
! check again at the upstream which belong to the damage enveloppe. LEFM
! criteria are computed at the previous node, and the new calving absciss is
! interpolated between these two nodes.
! CAUTION : because of infinite value, this is not possible when the "previous"
! node of the damage enveloppe is on the upper surface.
!
!________________________________________________________________
! 
! Variable : FrontRetrat = How far the node has to retreat 
! Variable : DContour = shape of the damage enveloppe, such as D > Dc
! Variable : KI = stress intensity factor at each node
!
! Requires :
! USF_Damage.f90 for damage computation
! FlowDepth.f90  solver for depth
! ComputeDevStressNS for stress computation : Cauchy Stresses are recalculated if
! necessary.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CalvingSolver( Model, Solver, dt, TransientSimulation )

   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   USE GeneralUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver 
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
   
   TYPE(ValueList_t), POINTER :: Simulation, SolverParams, Material
   TYPE(Solver_t), POINTER :: PSolver
   TYPE(Element_t),POINTER :: Element
   TYPE(Mesh_t),POINTER :: Mesh, MeshInit
   TYPE(Variable_t), POINTER :: DamageVariable, FrontRetreatVariable, &
                                DepthVariable, StressVariable, ExtVar, &
                                TimeVar, DContourVariable, &
                                KIVariable
   REAL(KIND=dp), POINTER :: DamageValues(:), FrontRetreatValues(:),&
                             DepthValues(:), StressValues(:), &
                             DContourValues(:), KIValues(:)
   INTEGER, POINTER :: DamagePerm(:), FrontRetreatPerm(:), DepthPerm(:), &
                       StressPerm(:), TopNodePointer(:), BotNodePointer(:), &
                       UpNodePointer(:), DownNodePointer(:), &
                       DContourPerm(:), KIPerm(:)
   INTEGER :: NONodes, NodeNumber, t, i, ni, Ind(3,3), DIM, NodeIndex, &
              NodeIndex2, PreviousNode, NodeIndexPrevious, CalvingNode,&
              compteur, uu
   REAL (KIND=dp) :: tampon, CriticalValue, KIc, KIa, alpha, Nb_Nodes, &
                     Eps = 0.1_dp, MinDist, Dist, Mask2, Time, d, w, &
                     M1, M2, M3, KI, xa, ya, da, &
                     sealevel, PreviousD, Previousdc, Previousw, Previousxa, &
                     Previousya, PreviousKI, PreviousKISL, KISL, &
                     fa, fb, xb, xcalv, xcalvSL
   REAL (KIND=dp), ALLOCATABLE :: Mask(:), Xcriterion(:), todelete(:)
   LOGICAL :: FirstTime = .TRUE., CalvingOccurs = .FALSE., Found, Cauchy, &
              NoPreviousNode 
   CHARACTER*20 :: SolverName='CalvingSolver'

 
   SAVE FirstTime, Ind, DIM, KI, Mask, Xcriterion, todelete, Cauchy
   SAVE TopNodePointer, BotNodePointer, UpNodePointer, DownNodePointer, sealevel


  Mesh => Solver % Mesh

  CALL INFO(SolverName,"Entering the calving solver", level=3)

  IF (FirstTime) THEN
     FirstTime = .FALSE.
     NONodes = Model % NumberOfNodes
     ALLOCATE(Xcriterion(NONodes))
     ALLOCATE(Mask(NONodes))
     ALLOCATE(todelete(NONodes))
 
     DIM = CoordinateSystemDimension()

     MeshInit => LoadMesh(Model,'./',TRIM(Model % Mesh % Name),.FALSE.,1,0)
     IF( .NOT. ASSOCIATED(MeshInit) ) THEN
       CALL FATAL('CalvingSolver','Cannot load mesh directory ')
     ENDIF
     PRINT *,Trim(Model % Mesh % Name)        
     PRINT *,MeshInit % NumberOfNodes

     CALL VariableAdd(MeshInit % Variables,MeshInit,Solver, &
                  'Coordinate 1',1,MeshInit % Nodes % x )
     CALL VariableAdd(MeshInit % Variables,MeshInit,Solver, &
                  'Coordinate 2',1,MeshInit % Nodes % y )
     CALL VariableAdd(MeshInit % Variables,MeshInit,Solver, &
                  'Coordinate 3',1,MeshInit % Nodes % z )
     PSolver => Solver

     CALL Info(SolverName,'---------------------------------------',Level=4 )
     CALL Info(SolverName,'Performing mapping on a structured mesh',Level=4 )
     CALL Info(SolverName,'---------------------------------------',Level=4 ) 
   
       CALL DetectExtrudedStructure( MeshInit, PSolver, ExtVar, &
           TopNodePointer, BotNodePointer, UpNodePointer, DownNodePointer )
           NULLIFY(PSolver,MeshInit)

    ! Be sure that the Bottom-Left corner is node number 1 
     IF ( (BotNodePointer(1) .NE. 1) .OR. & 
     ( Model % Nodes % x(1) > MINVAL(Model % Nodes % x(:)) + Eps) ) THEN
       CALL FATAL('CalvingSolver','Bottom-Left corner MUST be node number 1')
     ENDIF
   
    ! Cauchy or Deviatoric Stress
     Material => GetMaterial()
     Cauchy = ListGetLogical( Material , 'Cauchy', Found )
     WRITE(*,*)'Cauchy',Cauchy
     IF (.NOT. Cauchy) THEN
       CALL FATAL('CalvingSolver', 'This Solver requires the Cauchy Stress tensor, &
                                    and not the deviatoric !!!!')
       RETURN
     ENDIF 
 
     sealevel = GetConstReal( Material, 'Sea level', Found) 
     IF (.NOT.Found) THEN
       CALL FATAL('CalvingSolver','No "Sea level" given')
     END IF
  
  ENDIF ! FirstTime

  TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
  Time = TimeVar % Values(1)

  FrontRetreatVariable => Solver % Variable
  FrontRetreatPerm  => FrontRetreatVariable % Perm
  FrontRetreatValues => FrontRetreatVariable % Values
  
  WRITE(SolverName, '(A)') 'CalvingSolver'
  
  NONodes = Model % NumberOfNodes

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get the different variables needed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! DContourVariable
  DContourVariable => VariableGet( Model % Mesh % Variables, 'DContour')
  DContourValues => DContourVariable % Values
  DContourPerm => DContourVariable % Perm
   
   ! KIVariable
  KIVariable => VariableGet( Model % Mesh % Variables, 'KI')
  KIValues => KIVariable % Values
  KIPerm => KIVariable % Perm

   !Get Depth
   DepthVariable => VariableGet( Model % Mesh % Variables, 'Depth')
   IF ( ASSOCIATED( DepthVariable ) ) THEN
     DepthValues => DepthVariable % Values
     DepthPerm => DepthVariable % Perm
   ELSE
     CALL FATAL( SolverName, 'need to get Flowdepth Solver for the use of depth' )
   END IF
   
   !Get Damage
   DamageVariable => VariableGet( Model % Mesh % Variables, 'Damage')
   IF ( ASSOCIATED( DamageVariable ) ) THEN
     DamageValues => DamageVariable % Values
     DamagePerm => DamageVariable % Perm
   ELSE
     CALL FATAL( SolverName, 'need to get Damage Solver for the use of damage' )
   END IF
   
   ! Get the Stress                     
   StressVariable => VariableGet( Model % Variables, 'Sxx' )
   IF ( ASSOCIATED( StressVariable ) ) THEN
     StressPerm    => StressVariable % Perm
     StressValues  => StressVariable % Values
   ELSE
     CALL FATAL( SolverName, 'need to get ComputeDevStress Solver for the use of Sxx' )
   END IF

   SolverParams => GetSolverParams()
   ! Get the damage and the LEFM parameters
   CriticalValue = GetConstReal( SolverParams, 'Critical Damage',  Found )
   KIc =  GetConstReal( SolverParams, 'Ice Toughness',  Found )
   alpha =  GetConstReal( SolverParams, 'Crack Arrest Parameter',  Found )
   
   KIa = alpha*KIc


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 1st step : get the shape of damage threshold 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Mask(1:NONodes) = -1.0_dp
   ! Loop over all the nodes in order to see where
   ! the criterion on damage is satisfied
   ! This loop is necessary to initialize the mask
   DO t = 1, NONodes
      IF (DamageValues(DamagePerm(t)) > CriticalValue) THEN 
         Mask(t) = 1.0_dp      
      ELSE
         Mask(t) = -1.0_dp
      ENDIF
   ENDDO 

   todelete(1:NONodes) = -1.0_dp
   DO t = 1, NONodes
      IF ( (Mask(t) > 0.0_dp) .AND. (Mask(DownNodePointer(t)) > 0.0_dp)  ) THEN
         todelete(t) = 1.0_dp ! Mark the nodes whose mask must be switched to -1 
      ENDIF
   ENDDO 

   DO t = 1, NONodes
      IF ( (todelete(t) > 0.0_dp) .OR. (t == BotNodePointer(t)) ) THEN
         Mask(t) = -1.0_dp ! Turn the mask to -1for identified nodes 
      ENDIF
   ENDDO 

   DContourValues(DContourPerm(:)) = Mask(:)
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 2nd step : check the first point where 
   ! the LEFM criterion is reached 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   Mask2 = -1.0_dp
   Xcriterion = 9999999.0_dp !faut-il vraiment l'initialiser ??
   KIValues(KIPerm(:)) = 0.0_dp
     
 
   DO t = 1, NONodes
      IF ( Mask(t) < 0.0_dp ) CYCLE

      ! Here we need to compute the weight function m(y,d) associated with the 
      ! shape of the crevasse. This requires sxx*m(y,d) integrated over d
      d = DepthValues(DepthPerm(t))
      w = DepthValues(DepthPerm(BotNodePointer(t))) !wrong for highly-deformed mesh
      !Compute parameter for the weight function
      CALL ComputeWeightParameters(d,w,M1,M2,M3)
      !Integrate sxx*m(y,d)
      CALL ComputeStressIntensityFactor(M1,M2,M3,d,t,KI)
      KIValues(KIPerm(t)) = KI
      IF ( KI > KIc ) THEN
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! 3rd step : check if the LEFM arrest
         ! criterion is satisfied at the sea level 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         MinDist = 1000.0_dp
         NodeIndex = -1

         ! Here we need to compute the weight function m(y,d) associated with the 
         ! shape of the crevasse. 
         !depth where KI have to be calculated
         xa = Model % Nodes % x(t)
         ya = sealevel
         DO i = 1, NONodes
            Dist = SQRT((Model % Nodes % x(i) - xa)**2 &
                       + (Model % Nodes % y(i) - ya)**2 )
            IF ( Dist < MinDist ) THEN
               MinDist = Dist
               NodeIndex = i 
            ENDIF 
         ENDDO
         d = DepthValues(DepthPerm(NodeIndex))
         w = DepthValues(DepthPerm(BotNodePointer(t))) !may be wrong for highly-deformed mesh
         !Compute parameter for the weight function
         CALL ComputeWeightParameters(d,w,M1,M2,M3)
         !Integrate sxx*m(y,d)
         CALL ComputeStressIntensityFactor(M1,M2,M3,d,NodeIndex,KI)
         KIValues(KIPerm(NodeIndex)) = KI
         IF (KI > KIa ) THEN
         PRINT*,"NodeIndex",NodeIndex 
         PRINT*,"x(NodeIndex)",Model % Nodes % x(NodeIndex) 
         PRINT*,"y(NodeIndex)",Model % Nodes % y(NodeIndex) 
         PRINT*,"----------" 
             Xcriterion(t) = Model % Nodes % x (t)
             Mask2 = 1.0_dp !Calving definitely occurs
         ENDIF 
      ENDIF ! end loop KI>KIc
   ENDDO ! end loop on all nodes


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 4th step : The calving point may be somewhere between 
   ! the NodeIndex Node and the previous node
   ! Thus, we need to interpolate KI and KI at Sea Level 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! find the nodenumber before the one where calving should occur
   InDContour: IF (Mask2 > 0.5_dp) THEN
      DO t = 1,NoNodes
         IF (DContourValues(DContourPerm(t)) < 0.0_dp) CYCLE
         IF ( Xcriterion(t) < MINVAL(Xcriterion)+Eps) THEN
            PreviousNode = t-1
            CalvingNode = t
            PRINT*,"x(PreviousNode)=",Model % Nodes % x(PreviousNode)
            PRINT*,"x(t)=", Model % Nodes % x(t)
            KI = KIValues(KIPerm(t))
            EXIT
         ENDIF
      ENDDO

      ! -------------------------------------------------------------
      ! Get the node on the same vertical line which satisfies D>Dc
      uu = PreviousNode
      NoPreviousNode = .FALSE.
      IF (DamageValues(DamagePerm(PreviousNode)) > CriticalValue) THEN
         print*, "Previous node validate D>Dc"
         PreviousNode = uu
      ELSE
         compteur = 1 ! Compteur initialisation
         DO WHILE (UpNodePointer(uu) .NE. TopNodePointer(uu))
            uu = UpNodePointer(uu)   
            compteur = compteur + 1
         ENDDO
         uu = PreviousNode
         DO i=1,compteur+1
            IF (DamageValues(DamagePerm(uu)) < CriticalValue) THEN
               IF (uu == TopNodePointer(uu)) THEN
                  PRINT*,"There is no previous node validating D > Dc"
                  NoPreviousNode = .TRUE.
               ENDIF
               uu = UpNodePointer(uu)
            ELSE
               PreviousNode = uu
               print*, "The node which validates D>Dc is n°", PreviousNode
               IF (uu == TopNodePointer(uu)) THEN
                  PRINT*,"The node validating D > Dc is at surface, no KI &
                          computation possible"
                  NoPreviousNode = .TRUE.
               ENDIF
            ENDIF
         ENDDO
      ENDIF 

      Interpolation_IF: IF (.NOT. NoPreviousNode) THEN ! We interplate 
         PRINT*,"PreviousNode=", PreviousNode

         ! -------------------------------------------------------------
         ! Recompute KI at sea level for the original calving node
   
         MinDist = 1000.0_dp
         xa = Model % Nodes % x(CalvingNode)
         ya = sealevel
   
         DO i = 1, NONodes
            Dist = SQRT((Model % Nodes % x(i) - xa)**2 &
                       + (Model % Nodes % y(i) - ya)**2 )
            IF ( Dist < MinDist ) THEN
               MinDist = Dist
               NodeIndex2 = i 
            ENDIF 
         ENDDO
         d = DepthValues(DepthPerm(NodeIndex2))
         w = DepthValues(DepthPerm(BotNodePointer(NodeIndex2))) ! may be wrong for highly-deformed mesh
         !Compute parameter for the weight function
         CALL ComputeWeightParameters(d,w,M1,M2,M3)
         !Integrate sxx*m(y,d)
         CALL ComputeStressIntensityFactor(M1,M2,M3,d,NodeIndex2,KISL)



        ! -------------------------------------------------------------
        ! Infer KI at PreviousNode
         Previousdc = DepthValues(DepthPerm(PreviousNode)) 
         Previousw = DepthValues(DepthPerm(BotNodePointer(PreviousNode)))  ! may be wrong for highly-deformed mesh
         !Compute parameter for the weight function
         CALL ComputeWeightParameters(Previousdc,Previousw,M1,M2,M3)
         !Integrate sxx*m(y,d)
         CALL ComputeStressIntensityFactor(M1,M2,M3,Previousdc,PreviousNode,PreviousKI)


        ! -------------------------------------------------------------
        ! Infer KI(Sea level) at PreviousNode

         MinDist = 1000.0_dp
         Previousxa = Model % Nodes % x(PreviousNode)
         Previousya = sealevel

         DO i = 1, NONodes
            Dist = SQRT((Model % Nodes % x(i) - Previousxa)**2 &
                       + (Model % Nodes % y(i) - Previousya)**2 )
            IF ( Dist < MinDist ) THEN
               MinDist = Dist
               NodeIndexPrevious = i 
            ENDIF 
         ENDDO
         d = DepthValues(DepthPerm(NodeIndexPrevious))
         w = DepthValues(DepthPerm(BotNodePointer(NodeIndexPrevious))) !wrong for highly-deformed mesh
         !Compute parameter for the weight function
         CALL ComputeWeightParameters(d,w,M1,M2,M3)
         !Integrate sxx*m(y,d)
         CALL ComputeStressIntensityFactor(M1,M2,M3,d,NodeIndexPrevious,PreviousKISL)

         ! -------------------------------------------------------------
         ! Interpolate KI and KISL between two nodes 

            PRINT*,"KI=",KI
            PRINT*,"KISL=",KISL

            PRINT*,"PreviousKI=",PreviousKI
            PRINT*,"PreviousKISL=",PreviousKISL


         ! Interpolate KI
         xa =  Model % Nodes % x(PreviousNode)
         xb =  Model % Nodes % x(CalvingNode)

         fa = PreviousKI
         fb = KI

         xcalv = (KIc-fa)/(fb-fa)*(xb-xa) + xa 
    PRINT*,"xcalv =",xcalv

         ! Interpolate KISL
         xa =  Model % Nodes % x(PreviousNode)
         xb =  Model % Nodes % x(CalvingNode)

         fa = PreviousKISL
         fb = KISL

         xcalvSL = (KIa-fa)/(fb-fa)*(xb-xa) + xa 
    PRINT*,"xcalvSL =",xcalvSL

         tampon = MAX(xcalv, xcalvSL)
         IF (tampon > xb) tampon = xb
         IF (tampon < xa) tampon = xa

         PRINT*,"xa =",xa   
         PRINT*,"xb =",xb   
         PRINT*,"tampon =",tampon   

      ELSE ! Previous Node does not satisfy D>Dc, and no interpolation is possible
           ! In this case, we keep the original absciss for calving

         PRINT*,"Previous line does not validate D>Dc, or at surface only"
         PRINT*,"==> no possible interpolation, calving at original coord"
         tampon = MINVAL(Xcriterion) 
   
      ENDIF Interpolation_IF
   ENDIF InDContour 
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 5th step : 
   ! Compute the size of the calving event 
   ! can be improved and applied on front nodes only
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !Compute KI for each node (FOR VIZUALIZATION)
   DO t = 1,NoNodes
      d = DepthValues(DepthPerm(t))
      w = DepthValues(DepthPerm(BotNodePointer(t))) ! may be wrong for highly-deformed mesh
      IF (w-d < 0.2 * w) CYCLE
      CALL ComputeWeightParameters(d,w,M1,M2,M3)
      !Integrate sxx*m(y,d)
      CALL ComputeStressIntensityFactor(M1,M2,M3,d,t,KI)
      KIValues(KIPerm(t)) = KI
   ENDDO


   IF ( Mask2 > 0.0_dp) THEN
      CalvingOccurs = .TRUE.
      DO t = 1, NONodes
         IF (Model % Nodes % x ( t ) > tampon) THEN
            FrontRetreatValues(FrontRetreatPerm(t)) = - (Model % Nodes % x ( t ) - tampon)
         ELSE
            FrontRetreatValues(FrontRetreatPerm(t)) = 0.0
         ENDIF
      ENDDO
      CALL INFO(SolverName,"Calving Occurs", level=3)
   ELSE
      CalvingOccurs = .FALSE.
      FrontRetreatValues(FrontRetreatPerm(1:NONodes)) = 0.0_dp
      CALL INFO(SolverName,"No calving event", level=3)
   ENDIF

   ! Used in subsequent routines (FrontDisplacement, TwoMeshes, InterpBedrock)
   CALL ListAddLogical( Model % Simulation, 'CalvingOccurs', CalvingOccurs ) 
   ! Used in subsequent routines (FrontDisplacement, TwoMeshes)
   CALL ListAddConstReal( Model % Simulation, 'XCalving', tampon )

CONTAINS

!-----------------------------------------------------------------------------
  SUBROUTINE ComputeWeightParameters(d,w,M1,M2,M3)

  IMPLICIT NONE  

  REAL(KIND=dp), INTENT(IN) :: d,w
  REAL(KIND=dp), INTENT(OUT) :: M1,M2,M3

  M1 = 0.0719768 - 1.51346*(d/w) - 61.1001*(d/w)**2 + 1554.95*(d/w)**3 &
       - 14583.8*(d/w)**4 + 71590.7*(d/w)**5 - 205384.0*(d/w)**6 & 
       + 356469.0*(d/w)**7 -368270.0*(d/w)**8 + 208233.0*(d/w)**9 &
       - 49544.0*(d/w)**10

  M2 = 0.246984 + 6.47583*(d/w) + 176.456*(d/w)**2 - 4058.76*(d/w)**3 &
       + 37303.8*(d/w)**4 - 181755*(d/w)**5 + 520551*(d/w)**6 &
       - 904370*(d/w)**7 +936863*(d/w)**8 - 531940*(d/w)**9 &
       + 127291*(d/w)**10

   M3 = 0.529659 - 22.3235*(d/w) + 532.074*(d/w)**2 - 5479.53*(d/w)**3 &
        + 28592.2*(d/w)**4 - 81388.6*(d/w)**5 + 128746*(d/w)**6 &
        - 106246*(d/w)**7 + 35780.7*(d/w)**8 

  END SUBROUTINE ComputeWeightParameters
  
!-----------------------------------------------------------------------------

  SUBROUTINE ComputeStressIntensityFactor(M1,M2,M3,d,t,KI)
  
  IMPLICIT NONE  

  REAL(KIND=dp), INTENT(IN) :: M1, M2, M3, d
  INTEGER, INTENT(IN) :: t 
  REAL(KIND=dp), INTENT(OUT) :: KI

  INTEGER :: n, i, compteur
  REAL(KIND=dp), ALLOCATABLE :: y(:), sol(:)

  !Loop to find the number of nodes in the crevasse
  compteur = 1
  n = t ! node number of the actual node
  DO WHILE (UpNodePointer(n) .NE. TopNodePointer(n))
     n = UpNodePointer(n)   
     compteur = compteur + 1
  ENDDO
  compteur = compteur + 1 !Add the top layer node to the compteur  
  ALLOCATE(sol(compteur))
  ALLOCATE(y(compteur))
  
  n = t
  DO i = 1, compteur

    y(i) = DepthValues(DepthPerm(n))
   ! In order to avoid division by zero
    IF (i == 1)  y(1) = DepthValues(DepthPerm(n))-0.64 ! 05
    sol(i) = 2.0/sqrt(2.0*Pi*(d-y(i))) *    &
             ( 1.0 + M1*(1.0-y(i)/d)**0.5   & 
                   + M2*(1.0-y(i)/d)        &
                   + M3*(1.0-y(i)/d)**1.5)  &
             * StressValues(StressPerm(n))
    n = UpNodePointer(n)
  ENDDO
  
  ! Integrate over d
  KI = 0.0
  DO i = 1, compteur-1
     KI = KI + 0.5 * (sol(i)+sol(i+1)) * abs((y(i)-y(i+1)))
  ENDDO

  END SUBROUTINE ComputeStressIntensityFactor

!-----------------------------------------------------------------------------

END SUBROUTINE CalvingSolver 
