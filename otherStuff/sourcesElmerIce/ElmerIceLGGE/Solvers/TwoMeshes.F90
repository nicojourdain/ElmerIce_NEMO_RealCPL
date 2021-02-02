!*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  This solver may be used to map the field variables of a strecthed mesh
! *  to the unstreched initial configuarion of the same mesh. In addition
! *  the vertical size of the mesh is accomodated so that the height of 
! *  the mesh stays fixed. All-in-all two different mapping stages and 
! *  one solution of Laplace equation is needed for the mapping.
! *  The user may also optionally nullify some field variables such as the 
! *  Mesh Update which should freshly start from the new configuration. 
! *
! *  Note: If running in parallel, this solver requires the -halo switch
! *  when using ElmerGrid
! ******************************************************************************
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 19.10.2010
! *
! ****************************************************************************/ 

SUBROUTINE TwoMeshes( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE CRSMatrix
  USE GeneralUtils
  USE ElementDescription
  USE MeshUtils  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
   TYPE(Solver_t), POINTER :: PSolver

  TYPE(Mesh_t), POINTER :: MeshA => NULL(), MeshB => NULL(), Mesh0 => Null(),&
       MeshTop => Null(), MeshBot => Null(), MeshInit
  TYPE(Variable_t), POINTER :: Var, Var2, VarNew, VarTop, VarBot, VarDisplace1, VarDisplace2, &
                               CalvingMaskVariable, ExtVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), POINTER :: NodesA, NodesB, Nodes0, NodesTop, NodesBot
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName, HeightName, &
      TopMaskName, BotMaskName, FrontMaskName
  REAL(KIND=dp) :: BoundingBox(6),eps1,eps2, Norm, ExtensionCoefficient, &
      ExtensionAmplitude, rnd, TopDisplacement, BotDisplacement, H, B, T, y, p,&
      XCalving, homo, xold, xnew, b0, x0, L, L0
  INTEGER :: i,j,n,dim,istat,HeightDim, NoNodes, TopNodes, BotNodes, &
      FrontNodes, active, MapStep, TopCornerIndex, BotCornerIndex, ii, nn, jj, &
      NoLayers, kk
  INTEGER, POINTER :: TopPerm(:), BotPerm(:), FrontPerm(:), HeightPerm(:), FieldPerm(:), CalvingMaskPerm(:) ,&
                      TopNodePointer(:),BotNodePointer(:),UpNodePointer(:),DownNodePointer(:)
  REAL(KIND=dp), POINTER :: TopValues(:), BotValues(:), FrontValues(:), Height(:), ForceVector(:), Field(:), &
                            CalvingMaskValues(:)
  TYPE(Matrix_t), POINTER  :: StiffMatrix
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), NVNodes(:), n0(:), &
                                homov(:), NVNodesTemp(:), homoh(:)
  LOGICAL :: Found,RebuiltQuadrantTree,FSTop,FSBot,NeedInit=.FALSE., MapCondition
  INTEGER :: VisitedTimes = 0
  CHARACTER*20 :: SolverName='TwoMeshes'
  SAVE MeshA, MeshB, Mesh0, MeshTop, MeshBot, NodesA, NodesB, Nodes0, &
      NodesTop, NodesBot, dim, HeightDim, TopCornerIndex, BotCornerIndex, &
      TopPerm, BotPerm, FrontPerm, TopValues, BotValues, FrontValues, NoNodes,&
      TopNodes, BotNodes, TopMaskName, BotMaskName, Height, HeightPerm, &
      STIFF, FORCE, VisitedTimes, FSTop, FSBot
  SAVE TopNodePointer, BotNodePointer, UpNodePointer, DownNodePointer, homoh, &
       homov, b0, x0, H 

!------------------------------------------------------------------------------
  
  CALL Info( 'TwoMeshes', '-----------------------------------' )
  CALL Info( 'TwoMeshes', ' Projecting solution to itself' )
  CALL Info( 'TwoMeshes', '-----------------------------------' )

  VisitedTimes = VisitedTimes + 1
FirstTime:IF ( VisitedTimes == 1) THEN

  MeshInit => LoadMesh(Model,'./',TRIM(Model % Mesh % Name),.FALSE.,1,0)
  IF( .NOT. ASSOCIATED(MeshInit) ) THEN
    CALL FATAL('TwoMeshes','Cannot load mesh directory ')
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
 
   !If both top and bottom surfaces are free, we simply interpolate 
   !them from the previous (t-1) timestep.  However, repeated interpolation
   !on a fixed surface will lead to artificial smoothing, so we must store
   !the t=0 surface and interpolate from that each time.
    FSTop = ListGetLogical( Solver % Values, 'FS Top', Found)
    IF(.NOT. Found) THEN
       Call Warn('TwoMeshes','No value "FS Top" found, assuming Free Surface at surface.')
       FSTop = .TRUE.
    END IF

    FSBot = ListGetLogical( Solver % Values, 'FS Bottom', Found)
    IF(.NOT. Found) THEN
       Call Warn('TwoMeshes','No value "FS Bottom" found, assuming no Free Surface at bed.')
       FSBot = .FALSE.
    END IF

    IF((.NOT. FSBot).OR.(.NOT. FSTop))THEN
       NeedInit = .TRUE.
    END IF
    MeshA => Solver % Mesh
    NodesA => MeshA % Nodes

    MeshB => AllocateMesh()
    MeshB = MeshA
    MeshB % Name = TRIM(MeshA % Name)//'_copy'

    CALL Info('TwoMeshes','Allocated a new copy of the mesh')


    IF(NeedInit)THEN
       Mesh0 => AllocateMesh()
       Mesh0 = MeshA
       Mesh0 % Name = TRIM(Solver % Mesh % Name)//'_initial'
       CALL Info('TwoMeshes','Created initial reference mesh because at least one surface is fixed')
    END IF

    NoNodes = SIZE( MeshA % Nodes % x )
    NULLIFY( MeshB % Nodes )
    ALLOCATE( NodesB )
    MeshB % Nodes => NodesB

    ALLOCATE( NodesB % x(NoNodes), NodesB % y(NoNodes), NodesB % z(NoNodes) )
    NodesB % x = NodesA % x
    NodesB % y = NodesA % y
    NodesB % z = NodesA % z

    ALLOCATE( Nodes0 ) 
    ALLOCATE( Nodes0 % x(NoNodes), Nodes0 % y(NoNodes), Nodes0 % z(NoNodes) )
    Nodes0 % x = NodesA % x
    Nodes0 % y = NodesA % y
    Nodes0 % z = NodesA % z
    IF(NeedInit) Mesh0 % Nodes => Nodes0

    CALL Info('TwoMeshes','Created two new sets of coordinates')

    HeightPerm => Solver % Variable % Perm
    Height => Solver % Variable % Values

    !The convenience variables MeshTop/Bot and NodesTop/Bot remove the
    !need for multiple IF statements later on.  
    !------------------------------------------------------------
    IF(FSTop) THEN
       MeshTop => MeshA
       NodesTop => NodesA
    ELSE
       MeshTop => Mesh0
       NodesTop => Nodes0
    END IF

    IF(FSBot) THEN
       MeshBot => MeshA
       NodesBot => NodesA
    ELSE
       MeshBot => Mesh0
       NodesBot => Nodes0
    END IF

    ! Nullify the list of variables to take use of the automatic
    ! features of the mapping
    NULLIFY( MeshB % Variables )

    ! Create the variables in MeshB
    ! These need to be there for successful mapping
    !--------------------------------------------------------------
    DO i=1,99
      IF( i < 10 ) THEN
        WRITE( Name,'(A,I2)') 'Variable',i
      ELSE
        WRITE( Name,'(A,I3)') 'Variable',i
      END IF
      VarName = ListGetString( Solver % Values,TRIM(Name),Found)
      IF(.NOT. Found ) EXIT
      Var => VariableGet( MeshA % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 
      
      IF( ASSOCIATED( Var) ) THEN
        NULLIFY( Field )
        NULLIFY( FieldPerm ) 
        ALLOCATE( Field( SIZE( Var % Values ) ) ) 
    !This is a bit dirty.  Any node which isn't found by InterpolateMeshToMesh
    !will retain its 'previous value'. If you instead use Field = 0.0_dp, the
    !unfound nodes will take a value of zero.  This *shouldn't* be a problem 
    !because all nodes should be found...
        Field = Var % Values
        !Field = 0.0_dp
        ALLOCATE( FieldPerm( SIZE( Var % Perm ) ) )
        FieldPerm = Var % Perm
        CALL VariableAdd( MeshB % Variables, MeshB, Solver, TRIM( VarName), 1, Field, FieldPerm )

!The next two statements are required because VariableAdd nullifies the PrevValues field (why?) which is required by the InterpolateMeshToMeshQ subroutine.
!The dimension of the PrevValues array must be at least Nodes, 2, for some reason.
        VarNew => VariableGet( MeshB % Variables, TRIM( VarName), ThisOnly = .TRUE. )
        ALLOCATE( VarNew % PrevValues ( SIZE( VarNew % Values), 2 ) )
        CALL Info('TwoMeshes','Created variable: '//TRIM( VarName ) )
      END IF
    END DO

    DO i=1,99
      !Variables created here (listed in the sif Solver section) are those 
      !which are only defined on the upper and lower surface, and whose value
      !is equal to the height coordinate (y in 2D).  They cannot be 
      !'interpolated', because the non-boundary nodes will have a value of 0.
      !Instead, we simply equate them to the height coordinate.  
      !(Not yet generalized to 3D because I'm lazy!)

      IF( i < 10 ) THEN
        WRITE( Name,'(A,I2)') 'CoordVariable',i
      ELSE
        WRITE( Name,'(A,I3)') 'CoordVariable',i
      END IF
      VarName = ListGetString( Solver % Values,TRIM(Name),Found)
      IF(.NOT. Found ) EXIT
      Var => VariableGet( MeshA % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 
      
      IF( ASSOCIATED( Var) ) THEN
        NULLIFY( Field )
        NULLIFY( FieldPerm ) 
        ALLOCATE( Field( SIZE( Var % Values ) ) ) 
        Field = Var % Values
        !Field = 0.0_dp
        ALLOCATE( FieldPerm( SIZE( Var % Perm ) ) )
        FieldPerm = Var % Perm
        CALL VariableAdd( MeshB % Variables, MeshB, Solver, TRIM( VarName), 1, Field, FieldPerm )

!The next two statements are required because VariableAdd nullifies the PrevValues field (why?) which is required by the InterpolateMeshToMeshQ subroutine.
!The dimension of the PrevValues array must be at least Nodes, 2, for some reason.
        VarNew => VariableGet( MeshB % Variables, TRIM( VarName), ThisOnly = .TRUE. )
        ALLOCATE( VarNew % PrevValues ( SIZE( VarNew % Values), 2 ) )
        CALL Info('TwoMeshes','Created coordinate variable: '//TRIM( VarName ) )
      END IF
    END DO


    dim = CoordinateSystemDimension()
    HeightDim = dim

    ! Create variable for the top surface
    !---------------------------------------------------
    TopMaskName = 'Interp Top Surface'
    ALLOCATE( TopPerm(NoNodes) )
    CALL MakePermUsingMask( Model,Solver,MeshTop,TopMaskName, &
        .FALSE., TopPerm, TopNodes )
    ALLOCATE( TopValues(TopNodes) )
    CALL VariableAdd( MeshTop % Variables, MeshTop, Solver, &
        TopMaskName, 1, TopValues, TopPerm )

    ! Create variable for the bottom surface
    !---------------------------------------------------
    BotMaskName = 'Interp Bottom Surface'
    ALLOCATE( BotPerm(NoNodes) )
    CALL MakePermUsingMask( Model,Solver,MeshBot,BotMaskName, &
        .FALSE., BotPerm, BotNodes )
    ALLOCATE( BotValues(BotNodes) )
    CALL VariableAdd( MeshBot % Variables, MeshBot, Solver, &
        BotMaskName, 1, BotValues, BotPerm )

    ! Create variable for the calving front
    !---------------------------------------------------
    FrontMaskName = 'Calving Front'
    ALLOCATE( FrontPerm(NoNodes) )
    CALL MakePermUsingMask( Model,Solver,MeshBot,FrontMaskName, &
        .FALSE., FrontPerm, FrontNodes )
    ALLOCATE( FrontValues(FrontNodes) )
    CALL VariableAdd( MeshA % Variables, MeshA, Solver, &
        FrontMaskName, 1, FrontValues, FrontPerm )

    ! Find two corner nodes at calving front
    !---------------------------------------------------

    DO i=1,NoNodes
      IF(TopPerm(i) > 0 .AND. FrontPerm(i) > 0) THEN
         TopCornerIndex = i
      END IF
      IF(BotPerm(i) > 0 .AND. FrontPerm(i) > 0) THEN
         BotCornerIndex = i
      END IF
    END DO

    n = MeshA % MaxElementNodes 
    ALLOCATE( FORCE(n), STIFF(n,n), STAT=istat )

! Here we save the horizontal scaling of the nodes, that we re-use later for
! remeshing

! First, we save the horizontal scaling
    x0 = MINVAL(Model % Nodes % x(:))
    WRITE(*,*),'Mesh start at:', x0
    L0 = NodesA % x(BotCornerIndex) - x0 !Length of the mesh in the x-direction
    WRITE(*,*),'Length of the mesh in the x direction:', L0
    ALLOCATE(homoh(SIZE(BotValues(:))))
    j = 1
    DO ii = 1, MeshA % NumberOfNodes
       IF (DownNodePointer(ii) == ii) THEN ! We stay on the bottom line
          ! The horizontal scaling is saved in homoh vector (for "homothety
          ! horizontal")
          homoh(j) = (NodesA % x(ii) - x0) / L0
         ! WRITE(*,*),'NodesA % x(ii)=', NodesA % x(ii)
         ! WRITE(*,*),'homoh(j)=', homoh(j)
          j = j+1
       ENDIF
    ENDDO
    IF (MAXVAL(homoh(:)) .NE. 1.0_dp) THEN
       CALL FATAL(SolverName,"Bad detection of the horizontal scaling")
    ENDIF 

  ! Then, we save the vertical scaling
  ! Ice Thickness at the inlet (used later)
    DO i=1,NoNodes
       IF (NodesA % x(i) < x0 + 0.001) THEN
          H = NodesA % y(TopNodePointer(i)) - NodesA % y(BotNodePointer(i))  
       ENDIF
    END DO
    WRITE(*,*),'Ice Thickness at inlet boundary =', H
    ! May be improved using extruded mesh properties (TopNodePointer,...) 
    b0 = 999999.9
    DO ii = 1,NoNodes
       ! Looking for the node above, still on the inlet boundary
       IF (NodesA % x(ii) < x0 + 0.001) THEN
         IF (NodesA % y(ii) < b0) b0 = NodesA % y(ii)
       END IF
    END DO
    !WRITE(*,*),'b0= ',b0
    ALLOCATE(homov(SIZE(FrontValues(:))))
    j = 1
    DO ii = 1,NoNodes
     IF (NodesA % x(ii) < x0 + 0.001) THEN ! We stay on the left vertical
        ! boundary. The vertical scaling is save in the homov vector.
         homov(j) =  (NodesA % y(ii) - b0) / H
        ! WRITE(*,*),'NodesA % y(ii)=', NodesA % y(ii)
        ! WRITE(*,*),'homov(j)=', homov(j)
         j=j+1
     END IF
    END DO
    IF (MAXVAL(homov(:)) .NE. 1.0_dp) THEN
       CALL FATAL(SolverName,"Bad detection of the vertical scaling")
    ENDIF 

  END IF FirstTime

  ! Add a suitable condition here for calving
  !---------------------------------------------
  ! This has been linked to the CalvingSolver.f90 solver by means of the "CalvingOccurs"
  ! logical in the list.
  MapCondition = ListGetLogical( Model % Simulation, 'CalvingOccurs', Found)
  IF(.NOT. Found) THEN
     CALL WARN("Two Meshes","Can't find CalvingOccurs Logical, exiting... ")
     RETURN
  END IF
  
  IF ( Found .AND. .NOT. MapCondition ) THEN 
       CALL INFO("TwoMeshes","No Calving, leaving TwoMeshes Solver", level=3)
  RETURN
  ENDIF

  MapStep = ListGetInteger( Solver % Values,'Map Interval',Found)
  IF( Found .AND. MODULO( VisitedTimes, MapStep) /= 0 ) RETURN

  ExtensionAmplitude = ListGetCReal( Solver % values,'Extension Amplitude',Found )
  IF( Found ) THEN 
    CALL RANDOM_NUMBER( rnd )
    ExtensionCoefficient = 1.0_dp + rnd * ExtensionAmplitude 
    PRINT *,'rnd',rnd,ExtensionCoefficient
  END IF

  PRINT *,'Initial ranges'
  PRINT *,'X0:',MINVAL( Nodes0 % x), MAXVAL( Nodes0 % x)
  PRINT *,'XA:',MINVAL( NodesA % x), MAXVAL( NodesA % x)
  PRINT *,'XB:',MINVAL( NodesB % x), MAXVAL( NodesB % x)

  ! First copy the top and bottom height to variables
  !----------------------------------------------------------
  DO i=1,NoNodes
    j = TopPerm(i)
    IF( j > 0 ) THEN
      IF( HeightDim == 1 ) THEN
        TopValues(j) = NodesTop % x(i)
      ELSE IF( HeightDim == 2 ) THEN
        TopValues(j) = NodesTop % y(i)
      ELSE
        TopValues(j) = NodesTop % z(i)
      END IF
    END IF
    j = BotPerm(i)
    IF( j > 0 ) THEN
      IF( HeightDim == 1 ) THEN
        BotValues(j) = NodesBot % x(i)
      ELSE IF( HeightDim == 2 ) THEN
        BotValues(j) = NodesBot % y(i)
      ELSE
        BotValues(j) = NodesBot % z(i)
      END IF
    END IF
!TEST - is this ideal?
! GAG should it be if HeightDim=2 FrontValues(j) = x ?
    j = FrontPerm(i)
    IF( j > 0 ) THEN
      IF( HeightDim == 1 ) THEN
        FrontValues(j) = NodesA % x(i)
      ELSE IF( HeightDim == 2 ) THEN
        FrontValues(j) = NodesA % y(i)
      ELSE
        FrontValues(j) = NodesA % z(i)
      END IF
    END IF
  END DO
  

  ! Map the top and bottom variables to the reference mesh
  ! For that purpose set the coordinates of MeshB to initial ones
  ! Stretching may be accounted for here.
  !--------------------------------------------------------------
  NodesB % x = NodesA % x
  NodesB % y = NodesA % y
  NodesB % z = NodesA % z

  !At this point, our nodes are displaced due to calving.  If we deal only with
  !vertical calving faces and uniform calving retreat (i.e. a perfect rectangle
  !falls off) then we can simply say NodesB % x = NodesB % x - CalvingEventSize
  !However, I use a more complex method because I have non-vertical calving 
  !faces and events.

  VarDisplace1 => VariableGet(MeshA % Variables, 'Front Displacement 1', ThisOnly=.TRUE.)
  VarDisplace2 => VariableGet(MeshA % Variables, 'Front Displacement 2', ThisOnly=.TRUE.)

  XCalving = ListGetConstReal( Model % Simulation, 'XCalving', Found)
  IF(.NOT. Found) THEN
     CALL WARN("Two Meshes","Can't find XCalving real, exiting... ")
     RETURN
  END IF

  PRINT*, "========================"
  PRINT*, "Calving occurs at abscis :", XCalving
  PRINT*, "X-coord ot the bottom front node :", NodesA % x(BotCornerIndex)
  PRINT*, "========================"

! Displacement of the BotCornerNode :
! The opperation is different if the calving should happen upstream the
! BotCornerIndex Node or downstream: 
!
! 1 - Calving happens downstream (i.e. between BotCornerIndex and TopCornerIndex):
! The calving face is not necessarily vertical. We cannot properly align the
! nodes. Thus, we simply use the Variable VarDisplace from FrontDisplacement
! Solver to move all the mesh back. This may lead to progressive mesh
! degeneration.
!
! 2 - Calving happens upstream:
! The subsequent calving face is vertical : the bottom line is made shorter
! using the homoh vector, and all the nodes are then vertically aligned on the
! bottom layer
! 

! 1 -
  IF (XCalving > NodesA % x(BotCornerIndex)) THEN
     NodesB % x = NodesA % x + VarDisplace1 % Values(VarDisplace1 % Perm)
  ELSE
! 2 -   
     j = 1
     DO ii = 1, MeshA % NumberOfNodes
        IF (DownNodePointer(ii) == ii) THEN ! We are on the bottom line
           L = XCalving - x0
           NodesB % x(ii) = homoh(j) * L + x0  ! We applied the homothety
           j = j+1 
           nn=ii
           DO WHILE (UpNodePointer(nn) .NE. TopNodePointer(ii))
              NodesB % x(UpNodePointer(nn)) = NodesB % x(ii)
              nn = UpNodePointer(nn)   
           ENDDO
           NodesB % x(TopNodePointer(ii)) = NodesB % x(ii)
        ENDIF
     ENDDO 
  ENDIF

  NodesB % y = Nodes0 % y 
  
  ! Note that the stretching here assumes that (x=0) stays fixed. 
  ! This could be made more generic


  CALL InterpolateVarToVarReduced(MeshTop, MeshB, TopMaskName, HeightDim )
  CALL Info('TwoMeshes','Interpolated the top values to the original mesh')

  CALL InterpolateVarToVarReduced(MeshBot, MeshB, BotMaskName, HeightDim )
  CALL Info('TwoMeshes','Interpolated the bottom values to the original mesh')


  ! First copy the top and bottom height to variables
  !----------------------------------------------------------

  VarTop => VariableGet(MeshB % Variables, TopMaskName, ThisOnly = .TRUE.)
  VarBot => VariableGet(MeshB % Variables, BotMaskName, ThisOnly = .TRUE.)
  DO i=1,NoNodes
    j = TopPerm(i)
    IF( j > 0 ) THEN
        TopValues(j) = VarTop % Values(VarTop % Perm(i)) 
    END IF
    j = BotPerm(i) 
    IF( j > 0 ) THEN
        BotValues(j) = VarBot % Values(VarBot % Perm(i)) 
    END IF
  END DO

  ! Find the y-shift in the calving front nodes
  ! using the two end values from MeshA and B to guide the deformation
  ! This is necessary for situations where the calving front is not 
  ! vertical, as the nodes will tend towards the outwardmost end otherwise
  ! NB: This bit won't work in 3D!
  !----------------------------------------------------------

  H = NodesA % y(TopCornerIndex) - NodesA % y(BotCornerIndex)
  T = NodesA % y(TopCornerIndex)
  B = NodesA % y(BotCornerIndex)
  TopDisplacement = VarTop % Values(VarTop % Perm(TopCornerIndex)) - T 
  BotDisplacement = VarBot % Values(VarBot % Perm(BotCornerIndex)) - B 
  DO i=1,NoNodes
     j = FrontPerm(i)
     IF( j > 0 ) THEN
        y = Frontvalues(j)
        p = ((y-B)/H)
        FrontValues(j) = y + (p*TopDisplacement) + ((1-p)*BotDisplacement) 
     END IF
  END DO


  ! Then map the mesh using Laplace equation in height direction
  !----------------------------------------------------------
  CALL Info('TwoMeshes','Solving the height from poisson equation')
  Solver % Mesh % Nodes => NodesB

  CALL DefaultInitialize()

  active = GetNOFActive()
  DO i=1,active
    Element => GetActiveElement(i)
    n = GetElementNOFNodes()
    CALL LocalMatrix(  STIFF, FORCE, Element, n )
    CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO

  CALL DefaultFinishAssembly()

  ! Set top and bottom Dirichlet BCs using the ones mapped from the
  ! deformed mesh. This eliminates the need for external BCs. 
  !-----------------------------------------------------------------  
  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  DO i=1, NoNodes
    j =  TopPerm(i)
    IF( j > 0 ) THEN
      CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
          HeightPerm, i, TopValues(j) )                
    END IF
    j = BotPerm(i)
    IF( j > 0 ) THEN
      CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
          HeightPerm, i, BotValues(j) )                
    END IF
    j = FrontPerm(i)
    IF( j > 0 ) THEN
      CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
          HeightPerm, i, FrontValues(j) )                
    END IF
  END DO

  Norm = DefaultSolve()

  ! Return the nodes to the original ones in which all other variables 
  ! have been computed in order not to case problems later on...
  Solver % Mesh % Nodes => NodesA

  ! Copy the height now to the original mesh deformed to comply
  ! with the deformed top and bottom surfaces
  !------------------------------------------------------------
  DO i=1,NoNodes
    j = HeightPerm(i)
    IF( j == 0 ) CYCLE
    IF( HeightDim == 1 ) THEN
      NodesB % x(i) = Height(j)
    ELSE IF( HeightDim == 2 ) THEN
      NodesB % y(i) = Height(j)
    ELSE
      NodesB % z(i) = Height(j)
    END IF
  END DO


  ! Map the mesh to the real positions
  ! Always built a new tree as MeshA has changed
  ! This would ideally be built inside the interpolation...
  !--------------------------------------------------------------
  RebuiltQuadrantTree = .TRUE.
  IF( RebuiltQuadrantTree ) THEN
    BoundingBox(1) = MINVAL( NodesA % x )
    BoundingBox(2) = MINVAL( NodesA % y )
    BoundingBox(3) = MINVAL( NodesA % z )
    BoundingBox(4) = MAXVAL( NodesA % x )
    BoundingBox(5) = MAXVAL( NodesA % y )
    BoundingBox(6) = MAXVAL( NodesA % z )

    eps1 = 0.1_dp
    eps2 = eps1 * MAXVAL( BoundingBox(4:6) - BoundingBox(1:3) )
    BoundingBox(1:3) = BoundingBox(1:3) - eps2
    BoundingBox(4:6) = BoundingBox(4:6) + eps2

    CALL BuildQuadrantTree( MeshA,BoundingBox,MeshA % RootQuadrant)
  END IF


!C++ Jean
! Vertical adjustement for Nodes = stay at the x=0 absciss

    DO j = 1,NoNodes
        IF (DownNodePointer(j) == j) THEN
        nn = j
        kk=2
          DO WHILE (UpNodePointer(nn) .NE. TopNodePointer(j))
            H = NodesB % y(TopNodePointer(j)) - NodesB % y(j)  
            NodesB % y(UpNodePointer(nn)) = homov(kk)*H + NodesB % y(j)
            nn = UpNodePointer(nn)
            kk=kk+1
          ENDDO    
        ENDIF
    ENDDO
!C-- Jean
    
  ! This one uses the standard interpolation routines with two meshes
     CALL InterpolateMeshToMeshR( MeshA, MeshB, MeshA % Variables,MeshB % Variables, .TRUE. ) 
     CALL Info('TwoMesh','Interpolation done')

  ! Now we still have the problem that the interpolated values sit on MeshB while
  ! we would like to work with MeshA. It seems easiest to copy all the relevent
  ! stuff back to MeshA. 
  !------------------------------------------------------------------------------


  NodesA % x = NodesB % x
  NodesA % y = NodesB % y
  NodesA % z = NodesB % z
  DO i=1,99
    IF( i < 10 ) THEN
      WRITE( Name,'(A,I2)') 'Variable',i
    ELSE
      WRITE( Name,'(A,I3)') 'Variable',i
    END IF
    VarName = ListGetString( Solver % Values,TRIM(Name),Found)
    IF(.NOT. Found ) EXIT
    Var => VariableGet( MeshB % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 
    Var2 => VariableGet( MeshA % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 

    IF( ASSOCIATED( Var) .AND. ASSOCIATED( Var2 ) ) THEN
      CALL Info('TwoMeshes','Copying back variable: '//TRIM( VarName ) )
      Var2 % Values = Var % Values
      PRINT *,'Range ', TRIM(VarName),MINVAL( Var2 % Values), MAXVAL( Var2 % Values )
      !C++ Jean
    !VarName = ListGetString( Solver % Values,TRIM(Name),Found)
    !  PRINT *,'RangeOLD ', TRIM(VarName),MINVAL( Var % Values), MAXVAL( Var % Values )
    END IF
  END DO
  
  DO i=1,99
    IF( i < 10 ) THEN
      WRITE( Name,'(A,I2)') 'Nullify',i
    ELSE
      WRITE( Name,'(A,I3)') 'Nullify',i
    END IF
    VarName = ListGetString( Solver % Values,TRIM(Name),Found)
    IF(.NOT. Found ) EXIT
    Var2 => VariableGet( MeshA % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 
    IF( ASSOCIATED( Var2) ) THEN
      CALL Info('TwoMeshes','Zeroing variable: '//TRIM( VarName ) )
      Var2 % Values = 0.0_dp
    END IF
  END DO

!This bit not ready for 3D:
  DO i=1,99
    IF( i < 10 ) THEN
      WRITE( Name,'(A,I2)') 'CoordVariable',i
    ELSE
      WRITE( Name,'(A,I3)') 'CoordVariable',i
    END IF
    VarName = ListGetString( Solver % Values,TRIM(Name),Found)
    IF(.NOT. Found ) EXIT

    Var2 => VariableGet( MeshA % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 
    IF( ASSOCIATED( Var2 ) ) THEN
       CALL Info('TwoMeshes','Copying back variable: '//TRIM( VarName ) )
       ! Added by GAG  
       ! we treat the x variable first and the y after
       IF (TRIM(VarName) == 'xf') THEN
          DO j=1,SIZE(Var2 % Perm)
            IF(Var2 % Perm(j) == 0) CYCLE
            Var2 % Values(Var2 % Perm(j)) = NodesA % x(j) 
          END DO
          PRINT *,'Range ', TRIM(VarName),MINVAL( Var2 % Values), MAXVAL( Var2 % Values )
       ELSEIF (TRIM(VarName) == 'referencexf') THEN
          DO j=1,SIZE(Var2 % Perm)
            IF(Var2 % Perm(j) == 0) CYCLE
            Var2 % Values(Var2 % Perm(j)) = NodesA % x(j) 
          END DO
          PRINT *,'Range ', TRIM(VarName),MINVAL( Var2 % Values), MAXVAL( Var2 % Values )
       ELSE
          DO j=1,SIZE(Var2 % Perm)
             IF(Var2 % Perm(j) == 0) CYCLE
             Var2 % Values(Var2 % Perm(j)) = NodesA % y(j)
          END DO
          PRINT *,'Range ', TRIM(VarName),MINVAL( Var2 % Values), MAXVAL( Var2 % Values )
       ENDIF
    ELSE
      CALL FATAL('TwoMeshes','Coordinate variable listed but not associated: '//TRIM( VarName))
    END IF
  END DO


  PRINT *,'Final ranges'
  PRINT *,'X0:',MINVAL( Nodes0 % x), MAXVAL( Nodes0 % x)
  PRINT *,'XA:',MINVAL( NodesA % x), MAXVAL( NodesA % x)
  PRINT *,'XB:',MINVAL( NodesB % x), MAXVAL( NodesB % x)

  CALL Info('TwoMeshes','All done')



CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    
    FORCE = 0.0_dp
    STIFF = 0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t),  detJ, Basis, dBasisdx )
      STIFF(1:n,1:n) = STIFF(1:n,1:n) + IP % s(t) * DetJ * &
          MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix


!------------------------------------------------------------------------------
  SUBROUTINE SetDirichtletPoint( StiffMatrix, ForceVector,DOF, NDOFs, &
      Perm, NodeIndex, NodeValue) 
!------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(Matrix_t), POINTER :: StiffMatrix
    REAL(KIND=dp) :: ForceVector(:), NodeValue
    INTEGER :: DOF, NDOFs, Perm(:), NodeIndex
!------------------------------------------------------------------------------

    INTEGER :: PermIndex
    REAL(KIND=dp) :: s

!------------------------------------------------------------------------------

    PermIndex = Perm(NodeIndex)
    
    IF ( PermIndex > 0 ) THEN
      PermIndex = NDOFs * (PermIndex-1) + DOF
      
      IF ( StiffMatrix % FORMAT == MATRIX_SBAND ) THEN        
        CALL SBand_SetDirichlet( StiffMatrix,ForceVector,PermIndex,NodeValue )        
      ELSE IF ( StiffMatrix % FORMAT == MATRIX_CRS .AND. &
          StiffMatrix % Symmetric ) THEN        
        CALL CRS_SetSymmDirichlet(StiffMatrix,ForceVector,PermIndex,NodeValue)        
      ELSE                          
        s = StiffMatrix % Values(StiffMatrix % Diag(PermIndex))
        ForceVector(PermIndex) = NodeValue * s
        CALL ZeroRow( StiffMatrix,PermIndex )
        CALL SetMatrixElement( StiffMatrix,PermIndex,PermIndex,1.0d0*s )        
      END IF
    END IF
    
!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichtletPoint
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE InterpolateVartoVarReduced( OldMesh, NewMesh,HeightName,HeightDimension)

    USE SParIterComm
    USE Interpolation
    USE CoordinateSystems
!-------------------------------------------------------------------------------
    TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
    CHARACTER(LEN=*) :: HeightName
    INTEGER :: HeightDimension
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, nVar
    TYPE(Mesh_t), POINTER :: nMesh
    REAL(KIND=dp), POINTER :: OldHeight(:), NewHeight(:), nHeight(:)
    INTEGER, POINTER :: OldPerm(:), NewPerm(:), nperm(:)
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: i, j, k, l, n, ierr, npart, nfound, proc, status(MPI_STATUS_SIZE), unfound
    REAL(KIND=dp), DIMENSION(3) :: Point
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER, ALLOCATABLE :: perm(:), vperm(:)
    REAL(KIND=dp), DIMENSION(3) :: LocalCoordinates
    REAL(KIND=dp), POINTER :: ElementValues(:)
    TYPE(Element_t),POINTER :: Element       
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE, TARGET :: BB(:,:), nodes_x(:),nodes_y(:),nodes_z(:),vstore(:),&
         astore(:), xpart(:), ypart(:), zpart(:)
    REAL(KIND=dp) :: detJ, u,v,w,s, dn
    LOGICAL :: Found
    REAL(KIND=dp) :: eps1, eps2, eps_global, eps_local, myBB(6)
    LOGICAL, ALLOCATABLE :: FoundNodes(:)

!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------


    IF ( ParEnv % PEs<=1 ) THEN
       CALL InterpolateVarToVarReducedQ( OldMesh, NewMesh, HeightName, HeightDimension )
       RETURN
    END IF 

    !Need to add 'FoundNodes' array to InterpoateVarToVarReducedQ
    ALLOCATE( FoundNodes(NewMesh % NumberOfNodes) ); FoundNodes=.TRUE.
    CALL InterpolateVarToVarReducedQ( OldMesh, NewMesh, HeightName, HeightDimension, FoundNodes=FoundNodes )

    !Special case: all found
    !----------------------
    n = COUNT(.NOT.FoundNodes); dn = n
    CALL SParActiveSUM(dn,2)
    IF ( dn==0 ) RETURN
    !TEST
    PRINT *, 'Partition ',ParEnv % MyPE,' couldnt find ',n,'points!'

    ! Exchange partition bounding boxes:                             
    ! ----------------------------------                              
    myBB(1) = MINVAL(OldMesh % Nodes % x)
    myBB(2) = MINVAL(OldMesh % Nodes % y)
    myBB(3) = MINVAL(OldMesh % Nodes % z)
    myBB(4) = MAXVAL(OldMesh % Nodes % x)
    myBB(5) = MAXVAL(OldMesh % Nodes % y)
    myBB(6) = MAXVAL(OldMesh % Nodes % z)

!Possibly need to adjust this - is it necessary to be extending the bounding box by a 
!factor of 0
    eps2 = 0.1_dp * MAXVAL(myBB(4:6)-myBB(1:3))
    myBB(1:3) = myBB(1:3) - eps2
    myBB(4:6) = myBB(4:6) + eps2

    ALLOCATE(BB(6,ParEnv % PEs))
    DO i=1,ParEnv % PEs
      IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
      proc = i-1
      CALL MPI_BSEND( myBB, 6, MPI_DOUBLE_PRECISION, proc, &
               1099, MPI_COMM_WORLD, ierr )
    END DO
    DO i=1,COUNT(ParEnv % Active)-1
      CALL MPI_RECV( myBB, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
               1099, MPI_COMM_WORLD, status, ierr )
      proc = status(MPI_SOURCE)
      BB(:,proc+1) = myBB
    END DO

Sending:IF ( n==0 ) THEN
      DEALLOCATE(FoundNodes, BB)
      DO i=1,ParEnv % PEs
        IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
        proc = i-1
        CALL MPI_BSEND( n, 1, MPI_INTEGER, proc, &
              1101, MPI_COMM_WORLD, ierr )
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
        ! --------------------------------------------------------  
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
        !TEST
        PRINT *,'Partition ',ParEnv % MyPE,' sending points: ',xpart(1),ypart(1),' and ',xpart(2),ypart(2)
        ! send count...                                               
        ! -------------                                               
        CALL MPI_BSEND( npart, 1, MPI_INTEGER, proc, &
                1101, MPI_COMM_WORLD, ierr )
        IF ( npart==0 ) CYCLE
        ! ...and points                                               
        ! -------------                                               
        CALL MPI_BSEND( xpart, npart, MPI_DOUBLE_PRECISION, proc, &
                1102, MPI_COMM_WORLD, ierr )
        CALL MPI_BSEND( ypart, npart, MPI_DOUBLE_PRECISION, proc, &
                1103, MPI_COMM_WORLD, ierr )
        CALL MPI_BSEND( zpart, npart, MPI_DOUBLE_PRECISION, proc, &
                1104, MPI_COMM_WORLD, ierr )

        DEALLOCATE(xpart,ypart,zpart)
      END DO
      DEALLOCATE(nodes_x,nodes_y,nodes_z,BB)
    END IF Sending

    ! receive points from others:                                     
    ! ----------------------------                                    
    ALLOCATE(ProcRecv(Parenv % Pes))
    DO i=1,COUNT(ParEnv % Active)-1
      CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
            1101, MPI_COMM_WORLD, status, ierr )
 
      proc = status(MPI_SOURCE)
      ProcRecv(proc+1) % n = n

      IF ( n<=0 ) CYCLE

      ALLOCATE(ProcRecv(proc+1) % Nodes_x(n), &
            ProcRecv(proc+1) % Nodes_y(n),ProcRecv(proc+1) % Nodes_z(n))

      CALL MPI_RECV( ProcRecv(proc+1) % nodes_x, n, MPI_DOUBLE_PRECISION, proc, &
             1102, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV( ProcRecv(proc+1) % nodes_y, n, MPI_DOUBLE_PRECISION, proc, &
             1103, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV( ProcRecv(proc+1) % nodes_z, n, MPI_DOUBLE_PRECISION, proc, &
             1104, MPI_COMM_WORLD, status, ierr )
    END DO

    DO i=1,ParEnv % PEs
      IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE

      proc = i-1
      n = ProcRecv(i) % n

      IF ( n==0 ) THEN
        CALL MPI_BSEND( n, 1, MPI_INTEGER, proc, &
              2101, MPI_COMM_WORLD, ierr )
        CYCLE
      END IF

      ! Construct temporary mesh structure for the received points:   
      ! -----------------------------------------------------------   
      Nmesh => AllocateMesh()
      Nmesh % Nodes % x => ProcRecv(i) % nodes_x
      Nmesh % Nodes % y => ProcRecv(i) % nodes_y
      Nmesh % Nodes % z => ProcRecv(i) % nodes_z
      Nmesh % NumberOfNodes = n
!JOE TEST - This is the problem - the new variable created in the InterpolateVarToVarReducedQ for the NMesh gives a
!permutation table with values 0 1, for some reason...
      ALLOCATE(nperm(n))
      DO j=1,n
        nPerm(j)=j
      END DO
      ALLOCATE(nHeight(SIZE(nPerm)))
      nHeight = 0.0_dp
      CALL VariableAdd( NMesh % Variables, NMesh, CurrentModel % Solver, &
          HeightName, 1, nHeight, nPerm )

      ! try interpolating values for the points:                      
      ! ----------------------------------------                      
      ALLOCATE( FoundNodes(n) ); FoundNodes=.TRUE.
!TEST - remove last logical variable on interpvartovar
      CALL InterpolateVarToVarReducedQ( OldMesh, nMesh, HeightName, HeightDimension, &
           FoundNodes=FoundNodes )

      nfound = COUNT(FoundNodes)
      CALL MPI_BSEND( nfound, 1, MPI_INTEGER, proc, &
              2101, MPI_COMM_WORLD, ierr )

      !TEST
      unfound = n - nfound
      IF(unfound > 0) THEN
         PRINT *, 'TwoMeshes','Parallel: Found ',nfound,' nodes but still cant find ',unfound,' nodes!'
      END IF

      ! send interpolated values back to the owner:                   
      ! -------------------------------------------                   
      IF ( nfound>0 ) THEN
        ALLOCATE(vstore(nfound), vperm(nfound)); vstore=0
        k = 0
        DO j=1,n
          IF ( .NOT.FoundNodes(j)) CYCLE
          k = k + 1
          vperm(k) = j
          Nvar => VariableGet( Nmesh % Variables,HeightName,ThisOnly=.TRUE.)
          vstore(k)=Nvar % Values(j)
          !TEST
          PRINT *,'Partition ',ParEnv % MyPE,' found point with value: ',vstore(k)
        END DO
 
        CALL MPI_BSEND( vperm, nfound, MPI_INTEGER, proc, &
              2102, MPI_COMM_WORLD, status, ierr )

          CALL MPI_BSEND( vstore(:), nfound,MPI_DOUBLE_PRECISION, proc, &
                     2103, MPI_COMM_WORLD,ierr )

        DEALLOCATE(vstore, vperm)
      END IF

      DEALLOCATE(ProcRecv(i) % Nodes_x, ProcRecv(i) % Nodes_y,&
                 ProcRecv(i) % Nodes_z)
!      IF(ALLOCATED(Nvar)) DEALLOCATE(Nvar)
      DEALLOCATE(nperm,foundnodes, Nmesh)
    END DO
    DEALLOCATE(ProcRecv)

    ! Receive interpolated values:                                     
    ! ----------------------------                                     
     DO i=1,COUNT(ParEnv % Active)-1

       ! recv count:                                                   
       ! -----------                                                   
       CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
             2101, MPI_COMM_WORLD, status, ierr )

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
             2102, MPI_COMM_WORLD, status, ierr )

       ! recv values and store:                                        
       ! ----------------------                                        

       CALL MPI_RECV( astore, n, MPI_DOUBLE_PRECISION, proc, &
           2103, MPI_COMM_WORLD, status, ierr )

       Nvar => VariableGet( NewMesh % Variables,HeightName,ThisOnly=.TRUE.)

       IF ( ASSOCIATED(Nvar) ) THEN
         DO j=1,n
           k=perm(ProcSend(proc+1) % Perm(vperm(j)))
           IF ( Nvar % perm(k)>0 ) &
             Nvar % Values(Nvar % Perm(k)) = astore(j)
         END DO
       END IF

       DEALLOCATE(astore,vperm,ProcSend(proc+1) % perm)

     END DO
     IF ( ALLOCATED(Perm) ) DEALLOCATE(Perm,ProcSend)

     CALL MPI_BARRIER(ParEnv % ActiveComm,ierr)



!------------------------------------------------------------------------------
  END SUBROUTINE InterpolateVarToVarReduced

!------------------------------------------------------------------------------
  SUBROUTINE InterpolateVarToVarReducedQ( OldMesh, NewMesh,HeightName,HeightDimension, FoundNodes)
!This subroutine takes each boundary node on the specified boundary of the new mesh and finds its height (y coord in 2D) by performing (DIM - 1) interpolaton through boundary elements of the old mesh.
!------------------------------------------------------------------------------
    USE Interpolation
    USE CoordinateSystems
!-------------------------------------------------------------------------------
    TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
    CHARACTER(LEN=*) :: HeightName
    INTEGER :: HeightDimension
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp), POINTER :: OldHeight(:), NewHeight(:)
    INTEGER, POINTER :: OldPerm(:), NewPerm(:)
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: i, j, k, l, n
    REAL(KIND=dp), DIMENSION(3) :: Point
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp), DIMENSION(3) :: LocalCoordinates
    REAL(KIND=dp), POINTER :: ElementValues(:)
    TYPE(Element_t),POINTER :: Element       
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:)
    REAL(KIND=dp) :: detJ, u,v,w,s
    LOGICAL :: Found
    LOGICAL, OPTIONAL :: FoundNodes(:)
    REAL(KIND=dp) :: eps1, eps2, eps_global, eps_local
!------------------------------------------------------------------------------

    
    ! Get the height variable
    !---------------------------------------    
    Var => VariableGet( OldMesh % Variables, HeightName, ThisOnly = .TRUE. )
    OldHeight => Var % Values
    OldPerm => Var % Perm
    
    ! If the target variable does not exist, create it
    !----------------------------------------------------------
    Var => VariableGet( NewMesh % Variables, HeightName, ThisOnly = .TRUE. )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      ALLOCATE( NewHeight(SIZE(OldHeight) ) )
      NewHeight = 0.0_dp
      ALLOCATE( NewPerm(SIZE(OldPerm) ) )
      NewPerm = OldPerm
      CALL VariableAdd( NewMesh % Variables, NewMesh, CurrentModel % Solver, &
          HeightName, 1, NewHeight, NewPerm )
      Var => VariableGet( NewMesh % Variables, HeightName, ThisOnly = .TRUE. )
    END IF
    NewHeight => Var % Values
    NewPerm => Var % Perm

    !------------------------------------------------------------------------------
    n = OldMesh % MaxElementNodes
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), &
        ElementNodes % z(n), ElementValues(n) )
    ElementNodes % x = 0.0_dp
    ElementNodes % y = 0.0_dp
    ElementNodes % z = 0.0_dp
    
    eps_global = ListGetConstReal( CurrentModel % Simulation,  &
        'Interpolation Global Epsilon', Found)
    IF(.NOT. Found) eps_global = 2.0e-10
    
    eps_local = ListGetConstReal( CurrentModel % Simulation,  &
        'Interpolation Local Epsilon', Found )
    IF(.NOT. Found) eps_local = 1.0e-10

    !------------------------------------------------------------------------------
    ! Loop over all nodes in the new mesh
    !------------------------------------------------------------------------------
    DO i=1,NewMesh % NumberOfNodes
      !------------------------------------------------------------------------------


      IF( NewPerm(i) == 0 ) CYCLE
      
      Point(1) = NewMesh % Nodes % x(i)
      Point(2) = NewMesh % Nodes % y(i)
      Point(3) = NewMesh % Nodes % z(i)
      
      Point(HeightDim) = 0.0_dp
      !------------------------------------------------------------------------------
      ! Go through all old mesh bulk elements
      !------------------------------------------------------------------------------
      DO k=OldMesh % NumberOfBulkElements+1,&
          OldMesh % NumberOfBulkElements + OldMesh % NumberOfBoundaryElements

        Element => OldMesh % Elements(k)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes
        !TEST
        IF( ANY( OldPerm( NodeIndexes ) == 0 ) ) CYCLE
        
        IF( HeightDim /= 1 ) &
            ElementNodes % x(1:n) = OldMesh % Nodes % x(NodeIndexes)
        
        IF( HeightDim /= 2 ) &
            ElementNodes % y(1:n) = OldMesh % Nodes % y(NodeIndexes)
        
        IF( HeightDim /= 3 ) &
            ElementNodes % z(1:n) = OldMesh % Nodes % z(NodeIndexes)
        
        Found =  PointInElement( Element, ElementNodes, &
            Point, LocalCoordinates ) 
!TEST
        IF( i == 1 .AND. Found) THEN
           PRINT *, 'Point ',i,' with coords x ', Point(1),' y ',Point(2),' z ',Point(3)
           PRINT *, 'interpolated in element ',ElementNodes % x(1:n),' y ',ElementNodes % y(1:n)
        END IF
        IF( Found ) EXIT
      END DO
      
      IF (.NOT.Found) THEN
        NULLIFY(Element)
        IF(PRESENT(FoundNodes)) FoundNodes(i) = .FALSE.
        WRITE( Message, * ) 'Point was not found in any of the elements!',i
        CALL Warn( 'InterpolateVarToVarReduced', Message )
        !TEST
        NewHeight(NewPerm(i)) = OldHeight( OldPerm(i) )
        CYCLE
      END IF
      ElementValues(1:n) = OldHeight( OldPerm(NodeIndexes) )
      
      NewHeight(NewPerm(i)) = InterpolateInElement( &
          Element, ElementValues, LocalCoordinates(1), &
          LocalCoordinates(2), LocalCoordinates(3) )

    END DO

    DEALLOCATE( ElementNodes % x, ElementNodes % y, &
        ElementNodes % z, ElementValues )
    
!------------------------------------------------------------------------------
  END SUBROUTINE InterpolateVarToVarReducedQ
  
  
!------------------------------------------------------------------------------
END SUBROUTINE TwoMeshes
