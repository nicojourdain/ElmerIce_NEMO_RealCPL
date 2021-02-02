!  Fait l'optimisation d'une fonction cout
!    using quasi-Newton M1QN3 Routine in Reverse Communication
!    Using Euclidian inner product
!*
!    !!! Suppose une perturbation independante de z!!!!
!
! !!! Le maillage doit etre extrudé verticalement, i.e. les points alignes selon la verticale !!
!
!  Fait l'optimisation uniquement pour les points de la BC ou Name=String 'bed'
!    et applique la meme perturbation suivant la verticale
!
!  Serial Only   2D/3D
!
!
!Paramètres dans le sif:
!    Dans la section du Solver:
!       - Cost Variable Name = String  !! Nom de la variable dans lequel il y a le cout a minimiser (default "Costvalue")
!       - Optimized Variable Name = String  !! Nom de la variable a optimiser (default "Mu")
!       - Gradient Variable Name = String   !! Nom de la variable contenant le gradient (default "DJDMu")
!       - Optimisation Mask Variable (Optional): name of a mask variable. If
!         mask.lt.0 the variable is considered fixed
!
!       - Les parametres de m1qn3 (ont tous une valeur par default)    
!
!       - Active Coordinate = integer   !!!! La direction de l'extrusion
!
! *****************************************************************************
SUBROUTINE Optimize_m1qn3Serial_Z( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: SolverParams,BC
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  CHARACTER(LEN=MAX_NAME_LEN) :: BCName

!!!! Variables pour DetectExtrudedStructure
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Variable_t), POINTER :: Var
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),Upointer(:),DownPointer(:)
  
  TYPE(Element_t),POINTER ::  Element
  INTEGER, POINTER :: NodeIndexes(:)

  TYPE(Variable_t), POINTER :: Variable,CostVariable,GradVariable,MaskVar,TimeVar
  REAL(KIND=dp), POINTER :: Values(:),CostValues(:),GradValues(:),MaskValues(:)
  INTEGER, POINTER :: Perm(:),GradPerm(:),MaskPerm(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName,GradSolName,MaskVarName,NormFile


  REAL(KIND=dp),allocatable :: x(:),g(:)
  REAL(KIND=dp) :: f,NormG

  integer :: i,j,t,n,NMAX,NActiveNodes
  integer :: up
  integer,allocatable :: ActiveNodes(:)
  integer,allocatable :: NewNode(:)

  Logical :: FirstVisit=.true.,Found,UseMask,ComputeNormG=.False.
  logical,allocatable :: VisitedNode(:)

  CHARACTER*10 :: date,temps

!Variables for m1qn3
  external simul_rc,euclid,ctonbe,ctcabe
  character*3 normtype
  REAL(KIND=dp) :: dxmin,df1,epsrel,dzs(1)
  real(kind=dp), allocatable :: dz(:)
  REAL :: rzs(1)
  integer :: imp,io,imode(3),omode=-1,niter,nsim,iz(5),ndz,reverse,indic,izs(1)
  integer :: ierr,Npes,ntot
  CHARACTER(LEN=MAX_NAME_LEN) :: IOM1QN3,NormM1QN3
  logical :: DISbool


!
  save NActiveNodes
  save x,g
  save ActiveNodes
  save normtype,dxmin,df1,epsrel,dz,dzs,rzs,imp,io,imode,omode,niter,nsim,iz,ndz,reverse,indic,izs
  save FirstVisit
  save ComputeNormG,NormFile
  save SolverName
  save CostSolName,VarSolName,GradSolName,IOM1QN3
  SAVE TopPointer,BotPointer,Upointer,DownPointer


!  Read Constant from sif solver section
      IF(FirstVisit) Then
            FirstVisit=.FALSE.
            WRITE(SolverName, '(A)') 'Optimize_m1qn3Serial_Z'

           ! Check we don't have a parallel run
          IF(ASSOCIATED(Solver %  Matrix % ParMatrix)) Then
             CALL FATAL(SolverName,'ParMatrix associated! Ths solver for serial only!!')
          End if
!!!!!
          PSolver => Solver
          CALL DetectExtrudedStructure(Solver % Mesh , PSolver, Var, &
               TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
               UpNodePointer = Upointer , DownNodePointer = DownPointer)
!!!!!!!!!
          NMAX=Solver % Mesh % NumberOfNodes
          allocate(VisitedNode(NMAX),NewNode(NMAX))

            SolverParams => GetSolverParams()
          MaskVarName = GetString( SolverParams,'Optimisation Mask Variable',UseMask)
            IF (UseMask) Then
                MaskVar => VariableGet( Solver % Mesh % Variables, MaskVarName ) 
                IF (ASSOCIATED(MaskVar)) THEN 
                   MaskValues => MaskVar % Values 
                   MaskPerm => MaskVar % Perm 
               ELSE 
                   WRITE(Message,'(A,A,A)') 'No variable >',MaskVarName,'< found' 
                   CALL FATAL(SolverName,Message) 
               ENDIF
            ENDIF


!!!!!!!!!!!!find active nodes ! Uniquements ceux sur le Bed comme independant de z
           VisitedNode=.false.  
           NewNode=-1
           NActiveNodes=0 
           Do t=1,Solver % Mesh % NumberOfBoundaryElements
              Element => GetBoundaryElement(t)
              BC => GetBC()
              IF ( .NOT. ASSOCIATED(BC) ) CYCLE
              BCName =  ListGetString( BC,'Name', Found)
              IF(BCName /= 'bed') CYCLE

              n = GetElementNOFNodes()
              NodeIndexes => Element % NodeIndexes
              Do i=1,n
                 if (VisitedNode(NodeIndexes(i))) then
                     cycle
                 else
                     VisitedNode(NodeIndexes(i))=.true.
                     IF (UseMask) Then
                       IF (MaskValues(MaskPerm(NodeIndexes(i))).lt.0) cycle
                     END IF
                     NActiveNodes=NActiveNodes+1
                     NewNode(NActiveNodes)=NodeIndexes(i)
                 endif
             End do
           End do

           if (NActiveNodes.eq.0) THEN
              WRITE(Message,'(A)') 'NActiveNodes = 0 !! Require a boundary with <BCName = bed> !!'
              CALL FATAL(SolverName,Message)
           End if

           allocate(ActiveNodes(NActiveNodes),x(NActiveNodes),g(NActiveNodes))
           ActiveNodes(1:NActiveNodes)=NewNode(1:NActiveNodes)

           deallocate(VisitedNode,NewNode)

!!!!!!!  Solver Params

            CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
                END IF
            VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >Mu<')
                    WRITE(VarSolName,'(A)') 'Mu'
                END IF
            GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >DJDmu<')
                    WRITE(GradSolName,'(A)') 'DJDmu'
                END IF

             NormFile=GetString( SolverParams,'gradient Norm File',Found)
             IF(Found)  Then
                 ComputeNormG=.True.
                 open(io,file=trim(NormFile))
                    CALL DATE_AND_TIME(date,temps)
                    write(io,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)')'#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
                 close(io)
             END IF

!!  initialization of m1qn3 variables
            dxmin=GetConstReal( SolverParams,'M1QN3 dxmin', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 dxmin< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >1.e-10<')
                    dxmin=1.e-10
                END IF
            epsrel=GetConstReal( SolverParams,'M1QN3 epsg', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 epsg< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >1.e-06<')
                    epsrel=1.e-6
                END IF
            niter=GetInteger(SolverParams,'M1QN3 niter', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 niter< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >200<')
                    niter=200
                END IF
            nsim=GetInteger(SolverParams,'M1QN3 nsim', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 nsim< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >200<')
                    nsim=200
                END IF
            imp=GetInteger(SolverParams,'M1QN3 impres', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 impres< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >5<')
                    imp=5
                END IF
              ndz=GetInteger( SolverParams,'M1QN3 ndz', Found)
                  IF(.NOT.Found) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 ndz< not found  in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >5< update')
                       ndz=5
                   END IF
            DISbool=GetLogical( SolverParams, 'M1QN3 DIS Mode', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 DIS Mode< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >FALSE<')
                    DISbool=.False.
                END IF
                if(DISbool) then
                    imode(1)=0 !DIS Mode
                else
                    imode(1)=1 !SIS Mode
                End if
                IF (DISbool) then
                   ndz=4*NActiveNodes+ndz*(2*NActiveNodes+1)+10
                else
                   ndz=3*NActiveNodes+ndz*(2*NActiveNodes+1)+10
               end if
           allocate(dz(ndz))
            df1=GetConstReal( SolverParams,'M1QN3 df1', Found)
                IF(.NOT.Found) THEN
                   CALL WARN(SolverName,'Keyword >M1QN3 df1< not found  in section >Solver<')
                   CALL WARN(SolverName,'Taking default value >0.2<')
                   df1=0.2
                End if
                CostVariable => VariableGet( Solver % Mesh % Variables, CostSolName )
                IF (ASSOCIATED(CostVariable)) THEN
                    CostValues => CostVariable % Values
                 ELSE
                     WRITE(Message,'(A,A,A)') 'No variable >',CostSolName,'< found'
                     CALL FATAL(SolverName,Message)
                 ENDIF
                 df1=CostValues(1)*df1
             NormM1QN3 = GetString( SolverParams,'M1QN3 normtype', Found)
                 IF((.NOT.Found).AND.((NormM1QN3(1:3).ne.'dfn').OR.(NormM1QN3(1:3).ne.'sup') &
                     .OR.(NormM1QN3(1:3).ne.'two'))) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 normtype< not good in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >dfn<')
                       PRINT *,'M1QN3 normtype  ',NormM1QN3(1:3)
                       normtype = 'dfn'
                  ELSE
                       PRINT *,'M1QN3 normtype  ',NormM1QN3(1:3)
                       normtype = NormM1QN3(1:3)
                  END IF
              IOM1QN3 = GetString( SolverParams,'M1QN3 OutputFile', Found)
                 IF(.NOT.Found) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 OutputFile< not found  in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >M1QN3.out<')
                       WRITE(IOM1QN3,'(A)') 'M1QN3.out'
                 END IF
                 io=20
                 open(io,file=trim(IOM1QN3))
                    CALL DATE_AND_TIME(date,temps)
                    write(io,*) '******** M1QN3 Output file ************'
                    write(io,'(a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
                    write(io,*) '*****************************************'
                 close(io)

                    imode(2)=0 
                    imode(3)=0 
                    reverse=1 
                    omode=-1 
                    dzs=0.0 
                    rzs=0.0
                    izs=0

        End if


! Omode from previous iter; if > 0 m1qn3 has terminated => return 
     IF (omode.gt.0) then 
             WRITE(Message,'(a,I3)') 'm1qn3 finished; omode=',omode 
             CALL Info(SolverName, Message, Level=1) 
             return  
     End if

!  Get Variables CostValue, Mu and DJDmu
     CostVariable => VariableGet( Solver % Mesh % Variables, CostSolName )
     IF (ASSOCIATED(CostVariable)) THEN 
             CostValues => CostVariable % Values 
     ELSE
            WRITE(Message,'(A,A,A)') 'No variable >',CostSolName,'< found' 
            CALL FATAL(SolverName,Message) 
    ENDIF 
    f=CostValues(1)

     Variable => VariableGet( Solver % Mesh % Variables, VarSolName ) 
     IF (ASSOCIATED(Variable)) THEN 
             Values => Variable % Values 
             Perm => Variable % Perm 
     ELSE 
             WRITE(Message,'(A,A,A)') 'No variable >',VarSolName,'< found' 
             CALL FATAL(SolverName,Message) 
     ENDIF

     GradVariable => VariableGet( Solver % Mesh % Variables, GradSolName) 
     IF (ASSOCIATED(GradVariable)) THEN 
             GradValues   => GradVariable % Values 
             GradPerm => GradVariable % Perm 
     ELSE 
             WRITE(Message,'(A,A,A)') 'No variable >',GradSolName,'< found' 
             CALL FATAL(SolverName,Message)    
     END IF


     x(1:NActiveNodes)=Values(Perm(ActiveNodes(1:NActiveNodes)))
     g(1:NActiveNodes)=GradValues(GradPerm(ActiveNodes(1:NActiveNodes)))

     If (ComputeNormG) then
             TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
             NormG=0.0_dp
             Do i=1,NActiveNodes
                NormG=NormG+g(i)*g(i)
             End do
             open(io,file=trim(NormFile),position='append')
                write(io,'(e13.5,2x,e15.8)') TimeVar % Values(1),sqrt(NormG)
             close(io)
     End if      

     ! go to minimization
      open(io,file=trim(IOM1QN3),position='append')
       call m1qn3 (simul_rc,Euclid,ctonbe,ctcabe,NActiveNodes,x,f,g,dxmin,df1, &
                        epsrel,normtype,imp,io,imode,omode,niter,nsim,iz, &
                        dz,ndz,reverse,indic,izs,rzs,dzs)

     close(io)
     WRITE(Message,'(a,F15.8,x,I2)') 'm1qn3: Cost,omode= ',f,omode
     CALL Info(SolverName, Message, Level=3)

     ! Update Variable Values ! en conservant le profile Vertical
     Do i=1,NActiveNodes
        up=Upointer(ActiveNodes(i))
        Do while (up.ne.TopPointer(ActiveNodes(i)))
           Values(Perm(up))=x(i)+Values(Perm(up))-Values(Perm(ActiveNodes(i)))
           up=Upointer(up)
        End do
        Values(Perm(TopPointer(ActiveNodes(i))))=x(i)+ &
                    Values(Perm(TopPointer(ActiveNodes(i))))-Values(Perm(ActiveNodes(i)))

        Values(Perm(ActiveNodes(i)))=x(i)
     End do

   Return
!------------------------------------------------------------------------------
END SUBROUTINE Optimize_m1qn3Serial_Z
!------------------------------------------------------------------------------


