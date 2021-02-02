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
! ******************************************************************************
! *
! *  Authors: Lionel Favier, Nicolas Jourdain
! *  Email:   lionel.favier@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 21 March 2019
! * 
! *****************************************************************************
!
SUBROUTINE MISMIP_Melt_Param_ForIsmip6( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
USE CoordinateSystems
USE MeshUtils
USE DefUtils

IMPLICIT NONE
!------------------------------------------------------------------------------
TYPE(Model_t)  :: Model
TYPE(Solver_t), TARGET :: Solver
LOGICAL ::  Transient
REAL(KIND=dp) :: dt

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
TYPE(Mesh_t),POINTER :: Mesh
TYPE(Solver_t),POINTER :: PSolver
TYPE(Variable_t),POINTER :: MeltVar=>NULL(), GMVar=>NULL(), DepthVar=>NULL()
TYPE(Variable_t),POINTER :: tforcingVar=>NULL(), dtbasinVar=>NULL(), basinVar=>NULL()
TYPE(Variable_t),POINTER :: TimeVar=>NULL(), shelvesVar=>NULL()
TYPE(Nodes_t) :: ElementNodes
TYPE(GaussIntegrationPoints_t) :: IntegStuff
TYPE(Element_t),POINTER ::  Element

REAL(kind=dp),allocatable :: VisitedNode(:),Basis(:),dBasisdx(:,:)
REAL(kind=dp) :: u,v,w,SqrtElementMetric,s

INTEGER , POINTER :: MeltPerm(:), GMPerm(:), DepthPerm(:), NodeIndexes(:)
INTEGER , POINTER :: tforcingPerm(:), basinPerm(:), dtbasinPerm(:), shelvesPerm(:)
REAL(KIND=dp) , POINTER :: Melt(:), GM(:), Depth(:), tforcing(:), basin(:), dtbasin(:), shelves(:)

LOGICAL :: Found, stat, Parallel, llGL, addit

CHARACTER(len=MAX_NAME_LEN) :: para_name
CHARACTER(len = 200) :: meltValue
CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='MELT_MISMIP'

INTEGER :: Nmax, node, e, t, n, i, j, ii, ff, gg, ierr
INTEGER :: basinnb
INTEGER :: status(MPI_STATUS_SIZE)

REAL(KIND=dp) :: localInteg, Integ, Integ_Reduced
REAL(KIND=dp) :: Tf, time00, sealevel, meltfac, mskcrit, time, epsz, gamma0

REAL(KIND=dp), ALLOCATABLE :: basinarea(:), tfmean(:), localbasinarea(:), localtfmean(:)
REAL(KIND=dp), ALLOCATABLE :: basinarea_reduced(:), tfmean_reduced(:), AvgTmTf(:)

!------------------------------------------------------------------------------
! 1- Read constants and parameters of the simulation :
!------------------------------------------------------------------------------

!- Simulation parameters

time00 = GetCReal( Model % Simulation, 'Experiment initial time',Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Experiment initial time> not found')

para_name = GetString( Model % Simulation, 'Melt Parameterization',Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Melt Parameterization> not found')

llGL = GetLogical( Model % Simulation, 'Grounding Line Melt',Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Grounding Line Melt> not found')

!- General constants

sealevel = GetCReal( Model % Constants, 'Sea Level', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Sea Level> not found')

meltfac = GetCReal( Model % Constants, 'Melt factor',Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Melt factor> not found')

gamma0 = GetCReal( Model % Constants, 'gamma0', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<gamma0> not found')

basinnb = GetInteger( Model % Constants, 'basinnb', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<basinnb> not found')

epsz = GetCReal( Model % Constants, 'Minimum thickness to Melt', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<No melt layer thickness beneath sea level> not found')

!- Time varying variables :
TimeVar => VariableGet( Model % Variables,'Time')
time = TimeVar % Values(1)

!------------------------------------------------------------------------------
! 2- Read variables 
!------------------------------------------------------------------------------

Nmax = Solver % Mesh % NumberOfNodes

ALLOCATE(VisitedNode(Nmax),  &
        Basis(Model % MaxElementNodes),  &
        dBasisdx(Model % MaxElementNodes,3))

Mesh => Model % Mesh

Parallel = .FALSE.
IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
  IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
    Parallel = .TRUE.
  END IF
END IF

tforcingVar => VariableGet( Model % Mesh % Variables, 'thermal_forcing')
IF (.NOT.ASSOCIATED(tforcingVar)) &
&    CALL FATAL(SolverName,'thermal_forcing not found')

tforcingPerm => tforcingVar % Perm
tforcing => tforcingVar % Values

dtbasinVar => VariableGet( Model % Mesh % Variables, 'deltaT_basin')
IF (.NOT.ASSOCIATED(dtbasinVar)) &
&    CALL FATAL(SolverName,'deltaT_basin not found')

dtbasinPerm => dtbasinVar % Perm
dtbasin => dtbasinVar % Values

basinVar => VariableGet( Model % Mesh % Variables, 'basinNumber')
IF (.NOT.ASSOCIATED(basinVar)) &
&    CALL FATAL(SolverName,'basinNumber not found')

basinPerm => basinVar % Perm
basin => basinVar % Values

shelvesVar => VariableGet( Model % Mesh % Variables, 'shelves')
IF (.NOT.ASSOCIATED(shelvesVar)) &
&    CALL FATAL(SolverName,'shelves not found')

shelvesPerm => shelvesVar % Perm
shelves => shelvesVar % Values

MeltVar => VariableGet( Model % Mesh % Variables, 'Melt')
IF (.NOT.ASSOCIATED(MeltVar)) &
&    CALL FATAL(SolverName,'Melt not found')

MeltPerm => MeltVar % Perm
Melt => MeltVar % Values

GMVar => VariableGet( Model % Mesh % Variables, 'GroundedMask')
IF (.NOT.ASSOCIATED(GMVar)) &
&    CALL FATAL(SolverName,'GroundedMask not found')

GMPerm => GMVar % Perm
GM => GMVar % Values

DepthVar => VariableGet( Model % Mesh % Variables, 'Zb')
IF (.NOT.ASSOCIATED(DepthVar)) &
&    CALL FATAL(SolverName,'Zb not found')

DepthPerm => DepthVar % Perm
Depth => DepthVar % Values

if ( llGL ) then
  mskcrit =  0.5 ! Melt is at the Grounding Line and floating points
else
  mskcrit = -0.5 ! No melt at the Grounding Line, only floating points
endif

SELECT CASE (para_name)

CASE('quad')

  DO node = 1, Nmax

    IF ( GM(GMPerm(node)) .LT. mskcrit .AND. Depth(DepthPerm(node)) .LT. -epsz ) THEN
    !IF ( GM(GMPerm(node)) .LT. mskcrit ) THEN
      Melt(MeltPerm(node)) = -gamma0 * meltfac**2 * max(tforcing(tforcingperm(node))+ &
&                            dtbasin(dtbasinperm(node)),0.0)**2
    ELSE
      Melt(MeltPerm(node)) = 0.0_dp
    ENDIF

  ENDDO

CASE('plus')

  ALLOCATE(basinarea(basinnb),tfmean(basinnb))
  ALLOCATE(localbasinarea(basinnb),localtfmean(basinnb))
  ALLOCATE(basinarea_reduced(basinnb),tfmean_reduced(basinnb))
  ALLOCATE(AvgTmTf(basinnb))

  basinarea = 0.0_dp
  tfmean = 0.0_dp
  basinarea_reduced = 0.0_dp
  tfmean_reduced = 0.0_dp

  AvgTmTf = 0.0_dp

  !- we first store T0-Tf in the melt variable :
  DO node = 1, Nmax

    IF ( GM(GMPerm(node)) .LT. mskcrit  .AND. Depth(DepthPerm(node)) .LT. -epsz ) THEN
      ! NB: we first store TF in the melt variable :
      Melt(MeltPerm(node)) = tforcing(tforcingperm(node))
    ELSE
      Melt(MeltPerm(node)) = 0.0_dp
    ENDIF

  ENDDO

  !- we now calculate the average thermal forcing :
  ! Integrate melt rate over the Elmer grid for conservative interpolation purpose
  DO e=1,Solver % NumberOfActiveElements
    Element => GetActiveElement(e)

    CALL GetElementNodes( ElementNodes )
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes
    VisitedNode(NodeIndexes(1:n))=VisitedNode(NodeIndexes(1:n))+1.0_dp
    ! leave the loop if grounded point in the element
    IF ( ANY( GM(GMPerm(NodeIndexes(:))) .GE. mskcrit )  .OR. ANY( Depth(DepthPerm(NodeIndexes(:))) .GE. -epsz ) ) CYCLE

    localbasinarea = 0.0_dp
    localtfmean = 0.0_dp

    DO ff=1,basinnb
      addit = .FALSE.
      DO gg=1,n
        IF (NINT(basin(basinPerm(NodeIndexes(gg)))) .EQ. ff) THEN
          addit = .TRUE.
          EXIT
        ENDIF
      ENDDO

      IF (addit) THEN
        IntegStuff = GaussPoints( Element )
        DO t=1,IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric,Basis,dBasisdx )
          s = SqrtElementMetric * IntegStuff % s(t)

          localbasinarea(ff) = localbasinarea(ff) + s * SUM(Basis(:))
          localtfmean(ff) = localtfmean(ff) + s * SUM(Basis(:) * Melt(MeltPerm(NodeIndexes(:))))
        ENDDO

        tfmean(ff) = tfmean(ff) + (localtfmean(ff)/localbasinarea(ff)) * localbasinarea(ff)
        basinarea(ff) = basinarea(ff) + localbasinarea(ff)

      ENDIF
    ENDDO

  END DO

  IF (Parallel) THEN
    DO ii=1,basinnb
      CALL MPI_ALLREDUCE(basinarea(ii),basinarea_reduced(ii),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(tfmean(ii),tfmean_reduced(ii),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    ENDDO
  ELSE
    DO ii=1,basinnb
      basinarea_reduced(ii) = basinarea(ii)
      tfmean_reduced(ii) = tfmean(ii)
    ENDDO
  END IF

  DO ii=1,basinnb
    AvgTmTf(ii) = tfmean_reduced(ii) / basinarea_reduced(ii)
  ENDDO

  DO node = 1, Nmax

    IF ( GM(GMPerm(node)) .LT. mskcrit  .AND. Depth(DepthPerm(node)) .LT. -epsz ) THEN
    !IF ( GM(GMPerm(node)) .LT. mskcrit ) THEN
      Melt(MeltPerm(node)) = -gamma0 * meltfac**2 * (tforcing(tforcingperm(node))+dtbasin(dtbasinperm(node))) * &
&                                                  ABS(AvgTmTf(NINT(basin(basinPerm(node)))) + dtbasin(dtbasinperm(node)))
    ELSE
      Melt(MeltPerm(node)) = 0.0_dp
    ENDIF
  ENDDO

END SELECT

!---------------------------------------------

Integ = 0.0_dp  
! Integrate melt rate over the Elmer grid for conservative interpolation purpose
DO e = 1, Solver % NumberOfActiveElements
  Element => GetActiveElement(e)
  CALL GetElementNodes( ElementNodes )
  n = GetElementNOFNodes()
  NodeIndexes => Element % NodeIndexes
  VisitedNode(NodeIndexes(1:n)) = VisitedNode(NodeIndexes(1:n)) + 1.0_dp
  localInteg = 0.0_dp
  IntegStuff = GaussPoints( Element )
  DO t = 1, IntegStuff % n
    U = IntegStuff % u(t)
    V = IntegStuff % v(t)
    W = IntegStuff % w(t)
    stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
         Basis,dBasisdx )
    s = SqrtElementMetric * IntegStuff % s(t)
    localInteg = localInteg + s * SUM(Basis(:) * Melt(MeltPerm(NodeIndexes(:))))
  ENDDO
  Integ = Integ + localInteg
END DO

IF (Parallel) THEN
  CALL MPI_ALLREDUCE(Integ,Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
ELSE
  Integ_Reduced = Integ
ENDIF

IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
  WRITE(meltValue,'(F20.2)') Integ_Reduced 
  Message='TOTAL_MELT_RATE: '//meltValue
  CALL INFO(SolverName,Message,Level=1)
END IF

END SUBROUTINE MISMIP_Melt_Param_ForIsmip6


