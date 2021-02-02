!/*****************************************************************************/
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
! *****************************************************************************
!
!******************************************************************************
! *
! *  Authors: Olivier Gagliardini 
! *  Email:   olivier.gagliardini@ujf-grenoble.fr
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 19 June 2013
! *
! ****************************************************************************

!------------------------------------------------------------------------------
!>  Thickness equation solved with the discontinous Galerkin method.
!>  Work only in 3D (Thickness solved in the 2D horizontal plane then)
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE ThicknessSolverDG( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model            !< All model information (mesh, materials, BCs, etc...)
     TYPE(Solver_t), TARGET :: Solver  !< Linear & nonlinear equation solver options
     REAL(KIND=dp) :: dt               !< Timestep size for time dependent simulations
     LOGICAL :: Transient    !< Steady state or transient simulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: BC, BodyForce, Material,&
          Equation, Constants, SolverParams

     TYPE(Element_t), POINTER :: Element, Face, &
       ParentElement, LeftParent, RightParent

     LOGICAL :: AllocationsDone = .FALSE., Found, Stat, &
          LimitSolution, FoundLowerLimit
     INTEGER :: Active, DIM, NonLinearIterMin, NonlinearIterMax, iter,&
          CorrectedLowerLimit
     INTEGER :: n1, n2, k, n, t, istat, i, j, dummyInt, NumberOfFAces, Indexes(128)
     REAL(KIND=dp) :: Norm, RelativeChange, totat, st, totst,&
          OriginalValue 
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), &
              FORCE(:), Velo(:,:), MeshVelo(:,:), Ref(:), &
              LowerLimit(:)
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Variable_t), POINTER ::  Var
     CHARACTER(LEN=MAX_NAME_LEN) ::  VariableName, SolverName, ExpVariableName

     SAVE MASS, STIFF, LOAD, FORCE, Velo, MeshVelo, &
          AllocationsDone, DIM, VariableName, SolverName, &
          LowerLimit
!*******************************************************************************
     TYPE( Element_t ), POINTER :: Faces(:)
     TYPE(Nodes_t) :: ElementNodes

     INTEGER :: DOFs
!*******************************************************************************
     SolverName = 'ThicknessSolverDG ('// TRIM(Solver % Variable % Name) // ')'
     VariableName = TRIM(Solver % Variable % Name)
     WRITE(Message,'(A,A)')&
          'ThicknessSolverDG for variable ', VariableName
     CALL INFO(SolverName,Message,Level=1)

     Mesh => GetMesh()

     DIM = CoordinateSystemDimension() 
     IF (DIM /= 3) CALL FATAL( SolverName, 'Only work for 3D problem' )
     Faces => Mesh % Edges
     ! Only the edge for which ThickPerm > 0 (TODO)
     NumberOfFaces = Mesh % NumberOfEdges

     ! The Perm of the variable is false as we are running it on a BC
     ! Make it correct here
     Active = GetNOFActive()
     Solver % Variable % Perm = 0
     j = 0
     DO t=1, Active
        Element => GetActiveElement( t )
        n = GetElementNOfNodes( Element )
        dummyInt = GetElementDOFs( Indexes )
        DO i=1, n
           j = j + 1
           Solver % Variable % Perm(Indexes(i)) = j 
           write(*,*)Indexes(i),j
        END DO 
     END DO 
     write(*,*)'j',j,Active,Active*n
     ALLOCATE(Solver % Variable % Values(j))
     Solver % Variable % Values = 10.0


     ! Initialize & allocate some permanent storage, this is done first time only:
     !----------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
        N = 2 * MAX(Mesh % MaxElementDOFs, Mesh % MaxElementNodes )
        ALLOCATE( FORCE(N), MASS(n,n), STIFF(N,N), LOAD(N),  &
                  Velo(3,N), MeshVelo( 3,N ), &
                  LowerLimit(n), STAT = istat )
        
       IF ( istat /= 0 ) THEN
          CALL FATAL(SolverName,'Memory allocation error.' )
       ELSE
          CALL INFO(SolverName,'Memory allocation done',Level=1 )
       END IF
       DIM = CoordinateSystemDimension()
       AllocationsDone = .TRUE.
     END IF

     !------------------------------------------------------------------------------
     !    Read physical and numerical constants and initialize 
     !------------------------------------------------------------------------------
     Constants => GetConstants()
     SolverParams => GetSolverParams()

     NonlinearIterMax = GetInteger(   SolverParams, &
                     'Nonlinear System Max Iterations', Found )
     IF ( .NOT.Found ) THEN
        CALL WARN(SolverName,'No > Nonlinear System Max Iterations < found. Setting 1')
        NonlinearIterMax = 1
     END IF

     NonlinearIterMin = GetInteger(   SolverParams, &
                     'Nonlinear System Min Iterations', Found )
     IF ( .NOT.Found ) THEN        
        CALL WARN(SolverName,'No >Nonlinear System Min Iterations< found. Setting 1')
        NonlinearIterMin = 1
     ELSE IF (NonlinearIterMin > NonlinearIterMax) THEN
        CALL WARN(SolverName,&
             '>Nonlinear System Min Iterations< is exceeding >Nonlinear System Max Iterations<.')
         CALL WARN(SolverName,&
             'First is being reset to the latter.')
        NonlinearIterMin = NonlinearIterMax
     END IF

     ExpVariableName = GetString(SolverParams , 'Thickness Variable Name', Found )
     IF (.NOT.Found) &
        CALL FATAL(SolverName,'No value > Thickness Variable Name < found in Solver')

     ! Thickness is always limited to a minimum thickness
     LimitSolution = .True. 

     !------------------------------------------------------------------------------
     !       non-linear system iteration loop
     !------------------------------------------------------------------------------
     totat = 0; totst = 0
     DO iter=1,NonlinearIterMax
        ! Assembly of the bulk elements:
        !-------------------------------
        !------------------------------------------------------------------------------
        ! print out some information
        !------------------------------------------------------------------------------
        at  = CPUTime()
        at0 = RealTime()

        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, '-------------------------------------',Level=4 )
        WRITE( Message,'(A,I3,A,I3)') &
             'Nonlinear iteration no.', iter,' of max',NonlinearIterMax
        CALL Info( SolverName, Message, Level=4 )
        CALL Info( SolverName, '-------------------------------------',Level=4 )
        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, 'Starting Assembly...', Level=4 )


        CALL DefaultInitialize()
        Active = GetNOFActive()
        CorrectedLowerLimit = 0
        write(*,*)'SIZE',SIZE(Solver % Variable % Perm )
        DO i=1, SIZE(Solver % Variable % Perm )
          IF (Solver % Variable % Perm(i)/=0) &
          write(*,*)i,Solver % Variable % Perm(i)
        END DO 
        write(*,*)'SIZE Values',SIZE(Solver % Variable % Values)

        DO t = 1, Active  
           write(*,*)t,Active
           !------------------------------------------------------------------------------
           ! write some info on status of assembly
           !------------------------------------------------------------------------------
           IF ( RealTime() - at0 > 1.0 ) THEN
              WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                   (Active-t) / &
                   (1.0*Active)), ' % done'

              CALL Info( SolverName, Message, Level=5 )

              at0 = RealTime()
           END IF
           !------------------------------------------------------------------------------
           ! assign pointers and get number of nodes in element
           !------------------------------------------------------------------------------  
           Element => GetActiveElement( t )
           n = GetElementNOfNodes( Element )
           write(*,*)'n',n
           Material => GetMaterial()
           BodyForce => GetBodyForce( Element )
           Equation => GetEquation()

           dummyInt = GetElementDOFs( Indexes )
           write(*,*)'Indexes',Indexes(1:n)
           write(*,*)'Perm',Solver % Variable % Perm(Indexes(1:n))
           IF (ANY(Solver % Variable % Perm(Indexes(1:n))==0)) CYCLE

           LowerLimit(1:n) = GetReal(Material,'Minimum Ice Thickness', FoundLowerLimit)
           IF (.Not.FoundLowerLimit)  &
              CALL FATAL( SolverName, '> Minimum Ice Thickness < not given in Material' )
           
           DO i=1,n
              OriginalValue = Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
              write(*,*)OriginalValue, LowerLimit(i)
              IF (OriginalValue < LowerLimit(i)) THEN
                 Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) ) = LowerLimit(i)
                 CorrectedLowerLimit = CorrectedLowerLimit + 1
              END IF
           END DO
           !------------------------------------------------------------------------------
           ! the body force (r.h.s) = source
           !------------------------------------------------------------------------------         
           LOAD(1:n) = GetReal( BodyForce, 'Top Surface Accumulation ', Found )    
           IF (.NOT.Found) THEN
              WRITE(Message,'(A,A,A)') '> Top Surface Accumulation < not found'
              CALL INFO(SolverName, Message, Level=42)
              LOAD(1:n)  = 0.0d0
           END IF
           LOAD(1:n) = LOAD(1:n) + GetReal( BodyForce, 'Bottom Surface Accumulation ', Found )    
           IF (.NOT.Found) THEN
              WRITE(Message,'(A,A,A)') '> Bottom Surface Accumulation < not found'
              CALL INFO(SolverName, Message, Level=42)
           END IF
           !------------------------------------------------------------------------------
           ! Get convection and mesh velocity
           CALL GetLocalALEVelocity(Velo, MeshVelo, SolverName, Material,&
                Equation, Solver, Model, Element)
           !-----------------------
           ! get reaction constant
           !-----------------------

           CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, Velo, Element, n ) 
           IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
           CALL DefaultUpdateEquations( STIFF, FORCE )
        END DO

        ! Assembly of the face terms:
        !----------------------------
        FORCE = 0.0d0
        DO t=1,NumberOfFaces
           Face => Faces(t)
           IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

           LeftParent  => Face % BoundaryInfo % Left
           RightParent => Face % BoundaryInfo % Right
           IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
              n  = GetElementNOFNodes( Face )
              n1 = GetElementNOFNodes( LeftParent )
              n2 = GetElementNOFNodes( RightParent )

              !------------------------------------------------------------------------------
              ! Get convection and mesh velocity
              !------------------------------------------------------------------------------ 
              Material => GetMaterial( LeftParent )
              Equation => GetEquation( LeftParent )          
              CALL GetLocalALEVelocity(Velo, MeshVelo, SolverName, Material,&
                   Equation, Solver, Model, Face)

              CALL LocalJumps( STIFF, Face, n, LeftParent, n1, RightParent, n2, Velo )
              CALL DefaultUpdateEquations( STIFF, FORCE, Face )
           END IF
        END DO

        CALL DefaultFinishBulkAssembly()

        ! Loop over the boundary elements:
        !---------------------------------
        DO t=1,Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(t)
           IF( .NOT. ActiveBoundaryElement() )  CYCLE

           !------------------------------------------------------------------------------
           ! Check that the dimension of element is suitable for BCs
           !------------------------------------------------------------------------------
           IF( .NOT. PossibleFluxElement(Element) ) CYCLE

           ParentElement => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( ParentElement ) ) &
                ParentElement => Element % BoundaryInfo % Right

           n  = GetElementNOFNodes( Element )
           n1 = GetElementNOFnodes( ParentElement )

           !------------------------------------------------------------------------------
           ! Get convection and mesh velocity
           !------------------------------------------------------------------------------ 
           Material => GetMaterial( ParentElement )
           Equation => GetEquation( ParentElement )        
           CALL GetLocalALEVelocity(Velo, MeshVelo, SolverName, Material,&
                Equation, Solver, Model, Element)

           BC => GetBC()
           LOAD = 0.0d0
           Found = .FALSE.
           IF ( ASSOCIATED(BC) ) THEN
              LOAD(1:n) = GetReal( BC, TRIM(ExpVariableName), Found )
           END IF

           MASS = 0.0d0
           CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD, &
                Element, n, ParentElement, n1, Velo, Found )

           IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
           CALL DefaultUpdateEquations( STIFF, FORCE )
        END DO

        CALL DefaultFinishAssembly()
        CALL Info( SolverName, 'Assembly done', Level=4 )

        !------------------------------------------------------------------------------
        !     Solve the system and check for convergence
        !------------------------------------------------------------------------------
        at = CPUTime() - at
        st = CPUTime()

        ! Solve the system:
        !------------------
        Norm = DefaultSolve()

        st = CPUTIme()-st
        totat = totat + at
        totst = totst + st
        WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
        CALL Info( SolverName, Message, Level=4 )
        WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
        CALL Info( SolverName, Message, Level=4 )

        RelativeChange = Solver % Variable % NonlinChange 
        
        IF ( Solver % Variable % NonlinConverged == 1 )  THEN 
           WRITE(Message,'(A,I6,A,I6,A)') &
                'Nonlinear iteration converged after ', iter, &
                ' out of max ',NonlinearIterMax,' iterations'
           CALL INFO(SolverName,Message)
           EXIT
        ELSE IF ((iter .EQ. NonlinearIterMax) .AND. (NonlinearIterMax > 1)) THEN
           CALL WARN(SolverName,'Maximum nonlinear iterations reached, but system not converged')
        END IF

     !-----------------------------------------------
     END DO ! End of nonlinear iteration loop
     !----------------------------------------------

     !------------------------------------------------------------------------------
     ! limit solution
     !------------------------------------------------------------------------------
     CorrectedLowerLimit = 0
     Active = GetNOFActive()
     DO t = 1, Active
        Element => GetActiveElement( t )
        n = GetElementNOfNodes( Element )
        Material => GetMaterial()
        dummyInt = GetElementDOFs( Indexes )
        IF (ANY(Solver % Variable % Perm(Indexes(1:n))==0)) CYCLE
        LowerLimit(1:n) = GetReal(Material,'Minimum Ice Thickness', FoundLowerLimit)
        IF (.Not.FoundLowerLimit)  &
              CALL FATAL( SolverName, '> Minimum Ice Thickness < not given in Material' )

        DO i=1,n
           OriginalValue = Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
           IF (OriginalValue < LowerLimit(i)) THEN
              Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) ) = &
                         MAX(LowerLimit(i), OriginalValue )
              CorrectedLowerLimit = CorrectedLowerLimit + 1
           END IF
        END DO
     END DO
     WRITE(Message,'(a,i10)') 'Limited values for lower limit: ', CorrectedLowerLimit
     CALL Info( SolverName, Message, Level=3 )


     ! Average the elemental results to nodal values:
     !-----------------------------------------------
     Var => VariableGet( Mesh % Variables,TRIM(ExpVariableName))
     IF (.Not.ASSOCIATED(Var)) & 
        CALL FATAL( SolverName, TRIM(ExpVariableName) //' is not associated !' ) 

     IF  (Var % TYPE /= 0) & 
        CALL FATAL( SolverName, TRIM(ExpVariableName) //' should be a nodal variable' ) 
     
     n1 = Mesh % NumberOfNodes
     ALLOCATE( Ref(n1) )
     Ref = 0

     Var % Values = 0.0d0
     DO t=1,Active
        Element => GetActiveElement(t) 
        n = GetElementDOFs( Indexes )
        IF (ANY(Solver % Variable % Perm(Indexes(1:n))==0)) CYCLE
        n = GetElementNOFNodes()
        DO i=1,n
           k = Element % NodeIndexes(i)
           Var % Values(Var % Perm(k)) = Var % Values(Var % Perm(k)) + &
               Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
           Ref(k) = Ref(k) + 1
        END DO
     END DO

     WHERE( Ref > 0 )
        Var % Values(Var % Perm(1:n1)) = Var % Values(Var % Perm(1:n1)) / Ref
     END WHERE
     DEALLOCATE( Ref )

   CONTAINS

!------------------------------------------------------------------------------      
     SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, Velo, Element, n)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:), &
                  LOAD(:), Velo(:,:)
       INTEGER :: n
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n), dBasisdx(n,3)
       REAL(KIND=dp) :: detJ, U, V, W, S, A, L, cu(3)
       LOGICAL :: Stat
       INTEGER :: i, p, q, t, dim
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       TYPE(Nodes_t) :: Nodes
       SAVE Nodes
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       FORCE = 0.0d0
       STIFF = 0.0d0
       MASS  = 0.0d0
       CALL GetElementNodes( Nodes, Element )
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
                 Basis, dBasisdx )

         S = S * detJ
         L = SUM( LOAD(1:n) *  Basis(1:n) )

         cu   = 0.0d0
         DO i=1,dim-1
           cu(i) = SUM( Basis(1:n) * Velo(i,1:n) )
         END DO
!------------------------------------------------------------------------------
!        The advection-reaction equation: dc/dt + grad(u . c) + gamma c = s
!------------------------------------------------------------------------------
         DO p=1,n
            DO q=1,n
              MASS(p,q)  = MASS(p,q)  + s * Basis(q) * Basis(p)
              DO i=1,dim-1
                STIFF(p,q) = STIFF(p,q) - s * cu(i) * Basis(q) * dBasisdx(p,i)
              END DO
            END DO
         END DO
         FORCE(1:n) = FORCE(1:n) + s*L*Basis(1:n)
!------------------------------------------------------------------------------
      END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW(Face, nFace, Parent, nParent, U, V, W, Basis)
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(Element_t) :: Face, Parent
      INTEGER :: nFace, nParent
      REAL( KIND=dp ) :: U, V, W, Basis(:)
!------------------------------------------------------------------------------
      INTEGER :: i,j
      REAL(KIND=dp) :: ParentU(nFace), ParentV(nFace), ParentW(nFace)
!------------------------------------------------------------------------------
      DO i = 1,nFace
        DO j = 1,nParent
          IF ( Face % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
            ParentU(i) = Parent % TYPE % NodeU(j)
            ParentV(i) = Parent % TYPE % NodeV(j)
            ParentW(i) = Parent % TYPE % NodeW(j)
            EXIT
          END IF
        END DO
      END DO
      U = SUM( Basis(1:nFace) * ParentU(1:nFace) )
      V = SUM( Basis(1:nFace) * ParentV(1:nFace) )
      W = 0.0_dp 
!------------------------------------------------------------------------------      
    END SUBROUTINE FindParentUVW
!------------------------------------------------------------------------------      


!------------------------------------------------------------------------------
    SUBROUTINE LocalJumps( STIFF, Face, n, LeftParent, n1, RightParent, n2, Velo )
!------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=dp) :: STIFF(:,:), Velo(:,:)
      INTEGER :: n,n1,n2
      TYPE(Element_t), POINTER :: Face, LeftParent, RightParent
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: FaceBasis(n), FacedBasisdx(n,3)
      REAL(KIND=dp) :: LeftBasis(n1), LeftdBasisdx(n1,3)
      REAL(KIND=dp) :: RightBasis(n2), RightdBasisdx(n2,3)
      REAL(KIND=dp) :: Jump(n1+n2), Average(n1+n2)
      REAL(KIND=dp) :: detJ, U, V, W, S, Udotn, xx, yy
      LOGICAL :: Stat
      INTEGER :: i, j, p, q, dim, t, nFace, nParent
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: hE, Normal(3), cu(3), LeftOut(3)

      TYPE(Nodes_t) :: FaceNodes, LeftParentNodes, RightParentNodes
      SAVE FaceNodes, LeftParentNodes, RightParentNodes
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      STIFF = 0.0d0

      CALL GetElementNodes( FaceNodes, Face )
      CALL GetElementNodes( LeftParentNodes,  LeftParent )
      CALL GetElementNodes( RightParentNodes, RightParent )
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Face )

      LeftOut(1) = SUM(LeftParentNodes % x(1:n1)) / n1
      LeftOut(2) = SUM(LeftParentNodes % y(1:n1)) / n1
      LeftOut(3) = 0.0_dp
      LeftOut(1) = SUM(FaceNodes % x(1:n)) / n - LeftOut(1)
      LeftOut(2) = SUM(FaceNodes % y(1:n)) / n - LeftOut(2)
      LeftOut(3) = 0.0_dp 

      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Face, FaceNodes, U, V, W, detJ, &
             FaceBasis, FacedBasisdx )

        S = S * detJ

        Normal = NormalVector( Face, FaceNodes, U, V, .FALSE. )
        ! we are in 2D 
        Normal(3) = 0.0_dp
        Normal = Normal / SQRT(SUM(Normal*Normal))
        IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

        ! Find basis functions for the parent elements:
        ! ---------------------------------------------
        CALL FindParentUVW( Face, n, LeftParent, n1, U, V, W, FaceBasis )
        stat = ElementInfo( LeftParent, LeftParentNodes, U, V, W, detJ, &
                LeftBasis, LeftdBasisdx )

        CALL FindParentUVW( Face, n, RightParent, n2, U, V, W, FaceBasis )
        stat = ElementInfo( RightParent, RightParentNodes, U, V, W, detJ, &
              RightBasis, RightdBasisdx )

        ! Integrate jump terms:
        ! ---------------------
        Jump(1:n1) = LeftBasis(1:n1)
        Jump(n1+1:n1+n2) = -RightBasis(1:n2)

        Average(1:n1) = LeftBasis(1:n1) / 2
        Average(n1+1:n1+n2) = RightBasis(1:n2) / 2

        cu = 0.0d0
        DO i=1,dim-1
          cu(i) = SUM( Velo(i,1:n) * FaceBasis(1:n) )
        END DO
        Udotn = SUM( Normal * cu )

        DO p=1,n1+n2
          DO q=1,n1+n2
            STIFF(p,q) = STIFF(p,q) + s * Udotn * Average(q) * Jump(p)
            STIFF(p,q) = STIFF(p,q) + s * ABS(Udotn)/2 * Jump(q) * Jump(p)
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
        Element, n, ParentElement, np, Velo, InFlowBC )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), LOAD(:), Velo(:,:)
     INTEGER :: n, np
     LOGICAL :: InFlowBC
     TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3)
     REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3)
     INTEGER :: i,j,p,q,t,dim

     REAL(KIND=dp) :: Normal(3), g, L, Udotn, cu(3), cu1(3), detJ,U,V,W,S
     LOGICAL :: Stat, Inflow
     TYPE(GaussIntegrationPoints_t) :: IntegStuff

     TYPE(Nodes_t) :: Nodes, ParentNodes
     SAVE Nodes, ParentNodes
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     FORCE = 0.0d0
     STIFF = 0.0d0
 
     CALL GetElementNodes( Nodes, Element )
     CALL GetElementNodes( ParentNodes, ParentElement )

     Normal = NormalVector( Element, Nodes, 0.0d0, 0.0d0, .TRUE. ) 
     ! we are in 2D 
     Normal(3) = 0.0_dp
     Normal = Normal / SQRT(SUM(Normal*Normal))
     DO i=1,dim-1
       cu(i) = SUM( Velo(i,1:n) ) / n
     END DO
     Inflow = InFlowBC .AND. SUM( Normal * cu ) < 0.0d0

     ! Numerical integration:
     !-----------------------
     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       Normal = NormalVector( Element, Nodes, U, V, .TRUE. ) 

       ! Basis function values & derivatives at the integration point:
       ! -------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx )
       S = S * detJ

       CALL FindParentUVW( Element, n, ParentElement, np, U, V, W, Basis )
       stat = ElementInfo( ParentElement, ParentNodes, U, V, W, &
            detJ, ParentBasis, ParentdBasisdx )

       L = SUM( LOAD(1:n) * Basis(1:n) )
       cu = 0.0d0
       DO i=1,dim-1
          cu(i)  = SUM( Velo(i,1:n) * Basis(1:n) )
       END DO
       Udotn = SUM( Normal * cu )

       DO p = 1,np
         IF ( Inflow ) THEN
            FORCE(p) = FORCE(p) - s * Udotn*L*ParentBasis(p)
         ELSE
           DO q=1,np
             STIFF(p,q) = STIFF(p,q) + s*Udotn*ParentBasis(q)*ParentBasis(p)
           END DO
         END IF
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
   SUBROUTINE GetLocalALEVelocity(Velo, MeshVelo, SolverName, Material,&
        Equation, Solver, Model, Element)
     IMPLICIT NONE
     !------------------------------------------------------------------------------
     REAL (KIND=dp) :: Velo(:,:), MeshVelo(:,:)
     CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
     TYPE(ValueList_t), POINTER :: Material, Equation
     TYPE(Solver_t), TARGET :: Solver
     TYPE(Model_t), TARGET :: Model
     TYPE(Element_t),POINTER :: Element
     !------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, FlowSolName
     INTEGER :: i, j, k, N, FlowDOFs
     INTEGER, POINTER :: FlowPerm(:)
     REAL(KIND=dp), POINTER :: FlowSolution(:)
     TYPE(Variable_t), POINTER :: FlowSol
     LOGICAL :: Found
     !------------------------------------------------------------------------------
     n  = GetElementNOFNodes( Element )
     
     ! Velocity are given by the integrated velocity
     FlowSol => VariableGet( Model % Mesh % Variables, 'Integrated Velocity' )
     IF (ASSOCIATED(FlowSol)) THEN
        FlowPerm => FlowSol % Perm
        FlowDOFs =  FlowSol % DOFs
        FlowSolution => FlowSol % Values
     ELSE
        WRITE(Message,'(A)') '> Integrated Velocity < not associated'
        CALL FATAL( SolverName, Message )
     END IF

     DO i=1,n
        k = FlowPerm(Element % NodeIndexes(i))
        IF ( k > 0 ) THEN
           ! Only two component here
           Velo(1,i) = FlowSolution( FlowDOFs*k-1 )
           Velo(2,i) = FlowSolution( FlowDOFs*k )
           Velo(3,i) = 0.0d0
        END IF
     END DO

     !-------------------------------------------------
     ! Add mesh velocity (if any) for ALE formulation
     !-------------------------------------------------
     MeshVelo = 0.0d0
     CALL GetVectorLocalSolution( MeshVelo, 'Mesh Velocity', Element)
     IF (ANY( MeshVelo /= 0.0d0 )) THEN
        DO i=1,FlowDOFs
           Velo(i,1:N) = Velo(i,1:N) - MeshVelo(i,1:N)
        END DO
     END IF
   END SUBROUTINE GetLocalALEVelocity

!------------------------------------------------------------------------------
 END SUBROUTINE ThicknessSolverDG
!------------------------------------------------------------------------------
