!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! *****************************************************************************
SUBROUTINE SIAFluxSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t),TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(ValueList_t), POINTER,SAVE :: Material,BodyForce
  TYPE(Variable_t), POINTER :: PointerToVariable
  TYPE(Variable_t), POINTER :: ZsSol,DZsSol,ZbSol,SIABulkSol,HeatSol
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL, SAVE :: AllocationsDone = .FALSE.
  LOGICAL :: Found, GotIt
          
  INTEGER, POINTER :: Permutation(:)
  INTEGER, POINTER :: ZsPerm(:),DZsPerm(:),ZbPerm(:),SIAPerm(:),HeatPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp), POINTER :: Zs(:),DZs(:),Zb(:),SIA(:),Heat(:)

  REAL(KIND=dp) :: g,rho,rateA,NGlen
  REAL(KIND=dp) :: H,S,B,z
  REAL(KIND=dp) :: flux,grads
  REAL(KIND=dp) :: Hmesh,bottom,top
  REAL(KIND=dp) :: f1,f2,g1,g2
  REAL(KIND=dp) :: INTEGRAL,INTEGRAL2
  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: zr

  REAL(KIND=dp) :: U,V,W,detJ

  REAL(KIND=dp),DIMENSION(:),ALLOCATABLE,SAVE :: div,weight
  REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE,SAVE :: UV
  REAL(KIND=dp),ALLOCATABLE,SAVE  :: Basis(:),dBasisdx(:,:)
  LOGICAL :: stat

  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: SolverName
  CHARACTER(LEN=MAX_NAME_LEN) :: GName

  INTEGER,SAVE :: NActive
  INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: ActiveNode
  LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: IsActiveNode

  LOGICAL,SAVE :: VIntegrate=.FALSE.

  INTEGER :: DIM,STDOFs,NIndex
  INTEGER :: M,n
  INTEGER :: i,j,k,l

  TYPE(Solver_t), POINTER :: PSolver
  INTEGER, POINTER,SAVE :: TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:)
  INTEGER,SAVE :: nlayers,nlevels
  INTEGER :: CActive

  !Stuff for the parallel reduction
  TYPE lbuff_t
        INTEGER, ALLOCATABLE :: buff(:)
        REAL(KIND=dp), ALLOCATABLE :: values(:)
  END TYPE lbuff_t
  INTEGER, POINTER :: nlist(:)
  TYPE(lbuff_t), ALLOCATABLE :: n_index(:)
  REAL(KIND=dp), ALLOCATABLE :: nbuff(:)
  INTEGER, ALLOCATABLE :: n_count(:), gbuff(:), n_comp(:)
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
  INTEGER :: proc,ierr

!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  STDOFs = PointerToVariable % DOFs 
  WRITE(SolverName, '(A)') 'SIAFluxSolver'
  IF (.NOT.((STDOFS.EQ.1).OR.(STDOFS.EQ.2))) &
     CALL FATAL(SolverName,'Var DOFs has to be 1 or 2')

  SolverParams => GetSolverParams()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
        DIM = CoordinateSystemDimension()

        IF (DIM.GT.STDOFs) THEN
           ! for the special case where the mesh is of higher DIM than the pb
           ! but A uniform in the vertical => solution is analytical
           VIntegrate = ListGetLogical(SolverParams,'Do Vertical Integration',Found)
           IF (.NOT.Found) VIntegrate=.True.
           IF (VIntegrate) THEN
              CActive=ListGetInteger(SolverParams,'Active Coordinate',Found)
              IF (.NOT.Found) CALL FATAL(SolverName,'Keyword <Active Coordinate> not Found')
           END IF
        ENDIF
 
 ! Get Pointer to auxiliary variables requires to solve the SIA
        ZsSol => VariableGet( Solver % Mesh % Variables, 'Zs' )
        IF (ASSOCIATED(ZsSol)) THEN
           Zs => ZsSol % Values
           ZsPerm => ZsSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Zs<')
        END IF
        DZsSol => VariableGet( Solver % Mesh % Variables, 'DZs' )
        IF (ASSOCIATED(DZsSol)) THEN
           DZs => DZsSol % Values
           DZsPerm => DZsSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >DZs<')
        END IF
        ZbSol => VariableGet( Solver % Mesh % Variables, 'Zb' )
        IF (ASSOCIATED(ZbSol)) THEN
           Zb => ZbSol % Values
           ZbPerm => ZbSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Zb<')
        END IF

        SIABulkSol => VariableGet( Solver % Mesh % Variables, 'SIABulkVelocity' )
        IF (ASSOCIATED(SIABulkSol)) THEN
           SIA => SIABulkSol % Values
           SIAPerm => SIABulkSol % Perm
           IF (SIABulkSol % DOFs.NE.(STDOFs+2)) &
              CALL FATAL(SolverName,'SIABulkVelocity DOFs inconsitent with STDOFs')
        END IF
        HeatSol => VariableGet( Solver % Mesh % Variables, 'SIA Heat Source' )
        IF (ASSOCIATED(HeatSol)) THEN
           Heat => HeatSol % Values
           HeatPerm => HeatSol % Perm
        END IF


  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
       
      IF (VIntegrate) THEN
        PSolver => Solver
        CALL DetectExtrudedStructure( Solver%Mesh, PSolver, &
                 TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
                 UpNodePointer = UpPointer, DownNodePointer = DownPointer, &
                 NumberOfLayers = nlayers )
        nlevels = nlayers + 1
        IF (AllocationsDone) deallocate(zr)
        allocate(zr(nlevels))
      ENDIF

   !!! Get Material and BF associated with the first active element (should be
   !the same for all elements!!!)
     Element => GetActiveElement(1)
     Material => GetMaterial(Element)
     BodyForce=> GetBodyForce(Element)
     IF (.NOT.ASSOCIATED(Material)) CALL FATAL(SolverName,'No <Material> Found')
     IF (.NOT.ASSOCIATED(BodyForce)) CALL FATAL(SolverName,'No <Body Force> Found')

  !! Get Active Nodes
     M = Solver % Mesh % NumberOfNodes
     N = Model % Mesh % MaxElementNodes
     if (AllocationsDone) deallocate(IsActiveNode,ActiveNode,weight,div)
     if (AllocationsDone) deallocate(UV,Basis,dBasisdx)
     allocate(IsActiveNode(M),ActiveNode(M),weight(M),div(M))
     allocate(UV(N,2),Basis(N),dBasisdx(N,3))
     IsActiveNode=.False.
     NActive=0
     Do i=1,GetNOFActive(Solver)
       Element => GetActiveElement(i)
       n  = GetElementNOFNodes(Element)
       Do k=1,n
          NIndex=Element%NodeIndexes(k)
          IF (.NOT.IsActiveNode(NIndex)) THEN
             IsActiveNode(NIndex)=.TRUE.
             NActive=NActive+1
             ActiveNode(NActive)=NIndex
           ENDIF
        END DO
     END DO

     ! Allocate
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF


  Do i=1, NActive
     NIndex = ActiveNode(i)

     ! get gravity
     write(GName,'(A,I1)') 'Flow BodyForce ',STDOFs+1
     g = ListGetRealAtNode(BodyForce,trim(GName),NIndex,Found)
     write(Message,'(A,A)') trim(GName),' not found'
     IF (.NOT.FOUND) CALL FATAL(SolverName,trim(Message))
     !get Density
     rho = ListGetRealAtNode(Material,'Density',NIndex,Found)  !! LA FORMULATION
     IF (.NOT.FOUND) CALL FATAL(SolverName,'<Density> not found') !! SUPOSE rho uniform sur la verticale

     !get Glen rate factor                                                          
     rateA = ListGetRealAtNode(Material,'Rate Factor',NIndex,Found)
     IF (.NOT.FOUND) CALL FATAL(SolverName,'<Rate Factor> not found')
     NGlen = ListGetRealAtNode(Material,'Glen Exponent',NIndex,Found)
     IF (.NOT.FOUND) CALL FATAL(SolverName,'<Glen Exponent> not found')

     !get ice thickness
     S=Zs(ZsPerm(NIndex))
     B=Zb(ZbPerm(NIndex))
     H=S-B
     
     !get surface slope
     grads=0._dp
     Do j=1,STDOFs
        grads=grads+DZs(STDOFs*(DZsPerm(NIndex)-1)+j)*DZs(STDOFs*(DZsPerm(NIndex)-1)+j)
     End do

     ! the constant term
     flux=-2.0_dp*(rho*g)**NGlen
     flux=flux*grads**((NGlen-1.0)/2.0)

     !Do int_b_S A*(S-z) dz   
     IF (VIntegrate) THEN ! numerical integration
        
        !!! Use reduced coordinates, need to get Mesh thickness
        IF (STDOFS.eq.1) THEN
           bottom=Solver%Mesh%Nodes%y(BotPointer(NIndex))
           top=Solver%Mesh%Nodes%y(TopPointer(NIndex))
        ELSE
           bottom=Solver%Mesh%Nodes%z(BotPointer(NIndex))
           top=Solver%Mesh%Nodes%z(TopPointer(NIndex))
        ENDIF
        Hmesh=top-bottom

        zr(1)=0._dp
        INTEGRAL=0._dp
        INTEGRAL2=0._dp
        f1=rateA*(H*(1.0-zr(1)))**(NGlen+1.0)
        g1=rateA*(1.0-zr(1))**NGlen

        IF (ASSOCIATED(SIABulkSol)) THEN
           DO j=1,STDOFs
              SIA((STDOFS+2)*(SIAPerm(NIndex)-1)+j)=0._dp
           END DO
           SIA((STDOFS+2)*(SIAPerm(NIndex)))=rho*g*(1.0-zr(1))*H
        END IF
        IF (ASSOCIATED(HeatSol)) THEN
           Heat(HeatPerm(NIndex))=2.0*f1*(sqrt(grads)*rho*g)**(NGlen+1.0)
        END IF

        Do k=2,nlevels
          NIndex = UpPointer(NIndex)
          rateA = ListGetRealAtNode(Material,'Rate Factor',NIndex,Found)
          IF (STDOFS.eq.1) THEN
             z=Solver%Mesh%Nodes%y(NIndex)
          ELSE
             z=Solver%Mesh%Nodes%z(NIndex)
          ENDIF
          zr(k)=(z-bottom)/Hmesh
          f2=rateA*(H*(1.0-zr(k)))**(NGlen+1.0)
          g2=rateA*(1.0-zr(k))**NGlen
          INTEGRAL=INTEGRAL+0.5*(f1+f2)*(zr(k)-zr(k-1))
          INTEGRAL2=INTEGRAL2+0.5*(g1+g2)*(zr(k)-zr(k-1))*H**(NGlen+1.0)
          IF (ASSOCIATED(SIABulkSol)) THEN
             DO j=1,STDOFs
                SIA((STDOFS+2)*(SIAPerm(NIndex)-1)+j)=flux*INTEGRAL2*DZs(STDOFs*(DZsPerm(ActiveNode(i))-1)+j)
             END DO
             SIA((STDOFS+2)*(SIAPerm(NIndex)))=rho*g*(1.0-zr(k))*H
          END IF
          IF (ASSOCIATED(HeatSol)) THEN
             Heat(HeatPerm(NIndex))=2.0*f2*(sqrt(grads)*rho*g)**(NGlen+1.0)
          END IF
          f1=f2
          g1=g2
        END DO
     ELSE ! A unifrom in the vertical
       INTEGRAL=rateA*H**(NGlen+1.0)
       INTEGRAL=INTEGRAL/(NGlen+2.0)
     ENDIF
     flux=flux*INTEGRAL

     NIndex=ActiveNode(i)
     Do j=1,STDOFs
        VariableValues(STDOFs*(Permutation(NIndex)-1)+j)=flux*DZs(STDOFs*(DZsPerm(NIndex)-1)+j)
     End Do
  End Do
!------------------------------------------------------------------------------
 
 !! Compute vertical Bulk velocity
 IF (VIntegrate.AND.(ASSOCIATED(SIABulkSol))) THEN

 !!! NEED TO COMPUTE THE HORIZONTAL DIVERGENCE FIRST
     div=0.0
     weight = 0.0

     Do i=1,Model%Mesh%NumberOfBulkElements
       Element => Model%Elements(i)
       n  = GetElementNOFNodes(Element)
       CALL GetElementNodes( Nodes,Element )
      
       Do l=1,STDOFs
            UV(1:n,l)=SIA((STDOFS+2)*(SIAPerm(Element%NodeIndexes(1:n))-1)+l)
       End do
       
       !! compute div at the nodes; 
       !!   TODO: FE consistent average.
       Do j=1,n

          k = Element%NodeIndexes(j)

          U=Element % TYPE % NodeU(j)
          V=Element % TYPE % NodeV(j)
          W=Element % TYPE % NodeW(j)

          stat = ElementInfo( Element, Nodes, U, V, &
                  W,  detJ, Basis, dBasisdx)

          weight(k)=weight(k)+1.0_dp
          Do l=1,STDOFS
             div(k)=div(k)-SUM(UV(1:n,l)*dBasisdx(1:n,l))
          End Do
        End do
     END DO

     !!TODO: Parallel reduction
          IF (ParEnv % PEs>1 ) THEN  !! if parallel need to sum values at interfaces
                                     !! here is a copy of what is done in
                                     !SolverUtils.src to average boundary normals
             ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
             n_count = 0

             DO i=1,Solver%Mesh % NumberOfNodes
                IF (.NOT.Solver % Mesh % ParallelInfo % INTERFACE(i) ) CYCLE
  
                nlist => Solver%Mesh % ParallelInfo % NeighbourList(i) % Neighbours
                DO j=1,SIZE(nlist)
                   k = nlist(j)+1
                   IF ( k-1 == ParEnv % myPE ) CYCLE
                   n_count(k) = n_count(k)+1
                END DO
             END DO
             DO i=1,ParEnv % PEs
                IF ( n_count(i)>0 ) &
                  ALLOCATE( n_index(i) % buff(n_count(i)), &
                            n_index(i) % values(2*n_count(i)) )
             END DO

             n_count = 0
             DO i=1,Model % NumberOfNodes
               IF (.NOT.Solver % Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

               nlist =>Solver% Mesh % ParallelInfo % NeighbourList(i) % Neighbours
               DO j=1,SIZE(nlist)
                 k = nlist(j)+1
                 IF ( k-1 == ParEnv % myPE ) CYCLE
                 n_count(k) = n_count(k)+1
                 n_index(k) % buff(n_count(k)) = Solver%Mesh % Parallelinfo % &
                 GlobalDOFs(i)
                 n_index(k) % values(2*n_count(k)-1)=div(i)
                 n_index(k) % values(2*n_count(k))=weight(i)
               END DO
             END DO

             DO i=1,ParEnv % PEs
               IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
                CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                  900, MPI_COMM_WORLD, ierr )
                IF ( n_count(i)>0 ) THEN
                 CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                  901, MPI_COMM_WORLD, ierr )
                 CALL MPI_BSEND( n_index(i) % values, 2*n_count(i), MPI_DOUBLE_PRECISION, &
                    i-1,  902, MPI_COMM_WORLD, ierr )
               END IF
              END IF
            END DO
            DO i=1,ParEnv % PEs
               IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff, n_index(i) % values)

               IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
                  CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     900, MPI_COMM_WORLD, status, ierr )
                 IF ( n>0 ) THEN
                   proc = status(MPI_SOURCE)
                   ALLOCATE( gbuff(n), nbuff(2*n) )
                   CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                     901, MPI_COMM_WORLD, status, ierr )

                   CALL MPI_RECV( nbuff, 2*n, MPI_DOUBLE_PRECISION, proc, &
                     902, MPI_COMM_WORLD, status, ierr )

                   DO j=1,n
                     k = SearchNodeL(Solver% Mesh % ParallelInfo, gbuff(j), Solver%Mesh % NumberOfNodes )

                     IF ( k>0 ) THEN
                         div(k)=div(k)+nbuff(2*j-1)
                         weight(k)=weight(k)+nbuff(2*j)
                     END IF
                   END DO
                   DEALLOCATE(gbuff, nbuff)
                 END IF
               END IF
           END DO
           DEALLOCATE( n_index, n_count )
       END IF !end do parallel reduction

     !! Do the vertical integration to get w
     Do i=1, NActive
        NIndex = ActiveNode(i)

        S=Zs(ZsPerm(NIndex))
        B=Zb(ZbPerm(NIndex))
        H=S-B

        IF (STDOFS.eq.1) THEN
           bottom=Solver%Mesh%Nodes%y(BotPointer(NIndex))
           top=Solver%Mesh%Nodes%y(TopPointer(NIndex))
        ELSE
           bottom=Solver%Mesh%Nodes%z(BotPointer(NIndex))
           top=Solver%Mesh%Nodes%z(TopPointer(NIndex))
        ENDIF
        Hmesh=top-bottom

        INTEGRAL=0._dp
        SIA((STDOFS+2)*(SIAPerm(NIndex)-1)+STDOFs+1)=INTEGRAL

        zr(1)=0.0
        f1=div(NIndex)/weight(NIndex)
        Do k=2,nlevels
          NIndex = UpPointer(NIndex)
          IF (STDOFS.eq.1) THEN
             z=Solver%Mesh%Nodes%y(NIndex)
          ELSE
             z=Solver%Mesh%Nodes%z(NIndex)
          ENDIF
          zr(k)=(z-bottom)/Hmesh
          f2=div(NIndex)/weight(NIndex)
          INTEGRAL=INTEGRAL+0.5*(f1+f2)*(zr(k)-zr(k-1))*H
          SIA((STDOFS+2)*(SIAPerm(NIndex)-1)+STDOFs+1)=INTEGRAL
       !   PRINT *,k,f1,f2
          f1=f2
        END DO
      END DO  
  ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE SIAFluxSolver
!------------------------------------------------------------------------------

