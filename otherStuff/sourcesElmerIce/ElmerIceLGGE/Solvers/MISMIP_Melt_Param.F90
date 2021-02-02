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
! *  Email:   nicolas.jourdain@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 15 January 2019
! * 
! *****************************************************************************
! Computes the sub-shelf melt rates using ocean temperature and salinity based parameterizations
! Linked to Favier et al., 2019, GMD (doi:...)
! 
!
SUBROUTINE MISMIP_Melt_Param( Model,Solver,dt,Transient )
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
TYPE(Variable_t),POINTER :: isfslopeVar=>NULL(), distGLVar=>NULL(), distIFVar=>NULL()
TYPE(Variable_t),POINTER :: TimeVar=>NULL()
TYPE(Nodes_t) :: ElementNodes
TYPE(GaussIntegrationPoints_t) :: IntegStuff
TYPE(Element_t),POINTER ::  Element

REAL(kind=dp),allocatable :: VisitedNode(:),db(:),Basis(:),dBasisdx(:,:)
REAL(kind=dp),allocatable :: xgl(:),ygl(:),dgl(:),xglloc(:),yglloc(:),dglloc(:),xgltot(:),ygltot(:),dgltot(:)
REAL(kind=dp) :: u,v,w,SqrtElementMetric,s

INTEGER , POINTER :: MeltPerm(:), GMPerm(:), DepthPerm(:),NodeIndexes(:)
INTEGER , POINTER :: isfslopePerm(:), distGLPerm(:), distIFPerm(:)
REAL(KIND=dp) , POINTER :: Melt(:),GM(:), Depth(:), isfslope(:)
REAL(KIND=dp) , POINTER :: DATAPointer(:,:), distGL(:), distIF(:)
REAL(KIND=dp), POINTER :: pp(:,:)

LOGICAL :: Found, Got, stat, Parallel, llGL, cond1, cond2, cond3, ZaUniform, tmp

CHARACTER(len=MAX_NAME_LEN) :: para_name
CHARACTER(len = 200) :: meltValue
CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='MELT_MISMIP'

INTEGER :: Nmax, node, e, t, n, i, j, kk, nD, ii, LazerType, ierr
INTEGER :: ngl, nglloc, ngltot, cpt, nn, debilos, nbdir, cpttmp, jj, hwidth, kkstart, halfw, cptglob
INTEGER :: status(MPI_STATUS_SIZE)

REAL(KIND=dp) :: localInteg, Integ, Integ_Reduced, zzz, T0, S0,  &
&                z1_0, z2_0, z1, z2, Tsrf0, Tsrf, Tbot0, Tbot, Ssrf0, Ssrf, Sbot0, Sbot, dz1, dz2, dTsrf, dTbot, timsc, zb, &
&                Tf, xP , yP, time00, sealevel, lbd1, lbd2, lbd3, meltfac, K, gT, &
&                zGL, wn, div, x0, E0, Cd, GamT, GefT, GTS, MM, M_hat, X_hat, ll, alpha, AvgTmTf, localunity,  &
&                Area, Area_Reduced, Za, rhostar, CC, beta, g1, g2, rr, xbox, Tstar, tmp1, qqq,  &
&                mskcrit, time, epsz, sn, epsy, disty, zGLmin, yGLmin, xGLmin, radius, &
&                angle, angle1, angle2, xp2, yp2, xp3, yp3, xp12, yp12, xp23, yp23, xp31, yp31, zGLtmp, ltmp, det, &
&                zGLloc, mindist, dist, xGLfin, yGLfin, alphax, alphay, locslope, locslope1, locslope2, dy, ymax

REAL(KIND=dp), DIMENSION(:), ALLOCATABLE ::  Zbox, Abox, Tbox, Sbox, Mbox, xgltot2, ygltot2, dgltot2

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

IF (para_name .EQ. 'lin' .OR. para_name .EQ. 'quad' .OR. para_name .EQ. 'plus') THEN
  ZaUniform = GetLogical( Model % Simulation, 'Depth Ambiant Uniform', Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Depth Ambiant Uniform> not found')
ENDIF

!- General constants

sealevel = GetCReal( Model % Constants, 'Sea Level', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Sea Level> not found')

IF (para_name .EQ. 'lin' .OR. para_name .EQ. 'quad' .OR. para_name .EQ. 'plus' .OR. para_name .EQ. 'pico') THEN
  meltfac = GetCReal( Model % Constants, 'Melt factor',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Melt factor> not found')
ENDIF

IF (para_name .EQ. 'pico') THEN
  ZaUniform = .TRUE.
ENDIF

IF (para_name .EQ. 'lazer') THEN
  ZaUniform = .FALSE.
ENDIF

IF (ZaUniform) THEN
  Za = GetCReal( Model % Constants, 'Depth Ambiant',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Depth Ambiant> not found')
ENDIF

lbd1 = GetCReal( Model % Constants, 'Liquidus slope', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Liquidus slope> not found')

lbd2 = GetCReal( Model % Constants, 'Liquidus intercept', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Liquidus intercept> not found')

lbd3 = GetCReal( Model % Constants, 'Liquidus pressure coeff', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<Liquidus pressure coeff> not found')

epsz = GetCReal( Model % Constants, 'Minimum thickness to Melt', Found )
IF (.NOT.Found) &
&  CALL FATAL(SolverName,'<No melt layer thickness beneath sea level> not found')

!- Parameters according to chosen parameterization:

SELECT CASE (para_name)

! Linear dependency to Thermal Forcingi (Beckmann and Goose, 2002, OM)
CASE('lin')

  gT = GetCReal( Model % Constants, 'Multiplying Factor LIN', Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Multiplying Factor LIN> not found')

! Quadratic dependency to Thermal Forcing
! 'quad' uses local thermal forcing (Pollard & DeConto, 2012, GMD)
! 'plus' uses local/nonlocal thermal forcing

CASE('quad','plus')

  gT = GetCReal( Model % Constants, 'Multiplying Factor QUAD', Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Multiplying Factor QUAD> not found')

! 2D plume emulator model (Lazeroms et al., 2018, TC)
CASE('lazer')

  ! 4 different flavours
  ! 1 -> from the original article (Lazeroms et al., 2018, TC)
  ! 2 -> from the discussion paper (Lazeroms et al., 2017, TCD, Appendix B)
  ! 3 -> more simple law, the effective grounding line depth is taken at deepest Grounding line (only for ideal simulations so far)
  ! 4 -> Jenkins law, find the closest GL point, take the deepest GL counterclockwise (only working for ideal geometry of MISOMIP)
  LazerType = GetInteger( Model % Constants, 'Type LAZER', Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Type LAZER> not found')

  pp => ListGetConstRealArray( Solver % Values, 'Polynomial coefficients lazer', Found)
  IF (.NOT. Found) &
  & CALL FATAL(SolverName,'<Polynomial coefficients lazer> not found, need to be in Solver')

  if ( llGL ) CALL INFO(SolverName, &
  & 'WARNING : no melt at the grounding line with this parameterization',Level=1)

  E0 = GetCReal( Model % Constants, 'Entrainment Coeff LAZER',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Entrainment Coeff LAZER> not found')

  GamT = GetCReal( Model % Constants, 'Heat Exchange Coeff LAZER',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Heat Exchange Coeff LAZER> not found')

  GefT = GetCReal( Model % Constants, 'Effective Exchange Coeff LAZER',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Effective Exchange Coeff LAZER> not found')

  K = GetCReal( Model % Constants, 'K Coeff LAZER',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<K Coeff LAZER> not found')

  x0 = GetCReal( Model % Constants, 'x0 LAZER',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<x0 LAZER> not found')

  Cd = GetCReal( Model % Constants, 'Cd LAZER',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Cd LAZER> not found')

! The PICO box model for vertical overturning, from Reese et al. 2018
CASE('pico')

  nD = GetInteger( Model % Constants, 'Number Boxes PICO',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Number Boxes PICO> not found')

  CC = GetCReal( Model % Constants, 'Circulation Parameter PICO',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Circulation Parameter PICO> not found')

  gT = GetCReal( Model % Constants, 'Effective Exchange Velocity PICO',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Effective Exchange Velocity PICO> not found')

  alpha = GetCReal( Model % Constants, 'Thermal expansion coeff PICO',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Thermal expansion coeff PICO> not found')

  beta = GetCReal( Model % Constants, 'Salinity contraction coeff PICO',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<Salinity contraction coeff PICO> not found')

  rhostar = GetCReal( Model % Constants, 'EOS ref Density PICO',Found )
  IF (.NOT.Found) &
  &  CALL FATAL(SolverName,'<EOS ref Density PICO> not found')

CASE DEFAULT

  CALL FATAL(SolverName,'Unknown Melt Parameterization')

END SELECT

!------------------------------------------------------------------------------
! 2- Define temperature and salinity profiles 
! IN THE FUTURE, MAY HAVE TO BE CHANGED TO IMPORT PROFILES FROM EXTERNAL FILES
!------------------------------------------------------------------------------

!-- ISOMIP EXP3 (Warm0)
z1_0  =   0.0 ! initial depth of upper thermocline
z2_0  = 720.0 ! initial depth of lower thermocline
Tsrf0 =  -1.9 ! initial surface temperature
Tbot0 =   1.0 ! initial bottom temperature
Ssrf0 =  33.8 ! initial surface salinity 
Sbot0 =  34.7 ! initial bottom salinity (an increase of 0.2psu/700m is added below)
dz1   =   0.0 ! vertical displacement of the upper thermocline  (over time scale timsc)
dz2   =   0.0 ! vertical displacement of the lower thermocline  (over time scale timsc)
dTsrf =   0.0 ! warming  (over time scale timsc)
dTbot =   0.0 !    "            "
timsc =   1.e9! time scale for changes dz1, dz2, dTsrf, dTbot (in years)

!- Time varying variables :
TimeVar => VariableGet( Model % Variables,'Time')
time = TimeVar % Values(1)
z1 = z1_0 + (time-time00) * dz1 / timsc
z2 = z2_0 + (time-time00) * dz2 / timsc
Tsrf = Tsrf0 + (time-time00) * dTsrf / timsc
Tbot = Tbot0 + (time-time00) * dTbot / timsc
Ssrf = Ssrf0  ! no trend on salinity except if z1 and z2 vary
Sbot = Sbot0  ! no trend on salinity except if z1 and z2 vary

!------------------------------------------------------------------------------
! 2- Calculate melt rates at each floating  node
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

MeltVar => VariableGet( Model % Mesh % Variables, 'Melt')
IF (.NOT.ASSOCIATED(MeltVar)) &
&    CALL FATAL(SolverName,'Melt not found')

GMVar => VariableGet( Model % Mesh % Variables, 'GroundedMask')
IF (.NOT.ASSOCIATED(GMVar)) &
&    CALL FATAL(SolverName,'GroundedMask not found')

DepthVar => VariableGet( Model % Mesh % Variables, 'Zb')
IF (.NOT.ASSOCIATED(DepthVar)) &
&    CALL FATAL(SolverName,'Zb not found')

MeltPerm => MeltVar % Perm
Melt => MeltVar % Values

GMPerm => GMVar % Perm
GM => GMVar % Values

! THIS IS TO BE MODIFIED TO HAVE ZB NEGATIVE
DepthPerm => DepthVar % Perm
Depth => DepthVar % Values
Depth = sealevel - Depth     ! Depth < 0 under sea level

if ( llGL ) then
  mskcrit =  0.5 ! Melt is at the Grounding Line and floating points
else
  mskcrit = -0.5 ! No melt at the Grounding Line, only floating points
endif

SELECT CASE (para_name)

CASE('lin')

  DO node = 1, Nmax

    IF (ZaUniform) THEN
      zzz = Za
    ELSE
      zzz = Depth(DepthPerm(node))
    ENDIF

    zb = Depth(DepthPerm(node))

    IF ( GM(GMPerm(node)) .LT. mskcrit .AND. zzz .GE. epsz ) THEN
      IF ( zzz .le. z1 ) THEN
        T0 = Tsrf
        S0 = Ssrf
      ELSEIF ( zzz .GE. z2 ) THEN
        T0 = Tbot
        S0 = Sbot
      ELSE
        T0 = (Tbot-Tsrf)*zzz/(z2-z1) + Tsrf - (Tbot-Tsrf)*z1/(z2-z1)
        S0 = (Sbot-Ssrf)*zzz/(z2-z1) + Ssrf - (Sbot-Ssrf)*z1/(z2-z1)
      ENDIF
      Tf = lbd1*S0 + lbd2 + lbd3*zb  ! Sea water freezing temperature
      ! Melt in m/yr (meters of ice per year), positive if ice ablation
      Melt(MeltPerm(node)) = -gT * meltfac * (T0-Tf) ! 3 Equations with uniform velocity
    ELSE
      Melt(MeltPerm(node)) = 0.0_dp
    ENDIF
  ENDDO

CASE('quad')

  DO node = 1, Nmax

    IF (ZaUniform) THEN
      zzz = Za
    ELSE
      zzz = Depth(DepthPerm(node))
    ENDIF

    zb = Depth(DepthPerm(node))

    IF ( GM(GMPerm(node)) .LT. mskcrit .AND. zzz .GE. epsz ) THEN
      IF ( zzz .le. z1 ) THEN
        T0 = Tsrf
        S0 = Ssrf
      ELSEIF ( zzz .GE. z2 ) THEN
        T0 = Tbot
        S0 = Sbot
      ELSE
        T0 = (Tbot-Tsrf)*zzz/(z2-z1) + Tsrf - (Tbot-Tsrf)*z1/(z2-z1)
        S0 = (Sbot-Ssrf)*zzz/(z2-z1) + Ssrf - (Sbot-Ssrf)*z1/(z2-z1)
      ENDIF
      Tf = lbd1*S0 + lbd2 + lbd3*zb  ! Sea water freezing temperature
      Melt(MeltPerm(node)) = - gT * meltfac * abs(T0-Tf) * (T0-Tf)
    ELSE
      Melt(MeltPerm(node)) = 0.0_dp
    ENDIF
  ENDDO

CASE('plus')

   !- we first store T0-Tf in the melt variable :
  DO node = 1, Nmax

    IF (ZaUniform) THEN
      zzz = Za
    ELSE
      zzz = Depth(DepthPerm(node))
    ENDIF

    zb = Depth(DepthPerm(node))

    IF ( GM(GMPerm(node)) .LT. mskcrit .AND. zzz .GE. epsz ) THEN
      IF ( zzz .le. z1 ) THEN
        T0 = Tsrf
        S0 = Ssrf
      ELSEIF ( zzz .GE. z2 ) THEN
        T0 = Tbot
        S0 = Sbot
      ELSE
        T0 = (Tbot-Tsrf)*zzz/(z2-z1) + Tsrf - (Tbot-Tsrf)*z1/(z2-z1)
        S0 = (Sbot-Ssrf)*zzz/(z2-z1) + Ssrf - (Sbot-Ssrf)*z1/(z2-z1)
      ENDIF
      Tf = lbd1*S0 + lbd2 + lbd3*zb  ! Sea water freezing temperature
      ! NB: we first store (T0-Tf) in the melt variable :
      Melt(MeltPerm(node)) = T0-Tf
    ELSE
      Melt(MeltPerm(node)) = 0.0_dp
    ENDIF
  ENDDO

  !- we now calculate the average thermal forcing :
  Integ = 0.0_dp
  Area  = 0.0_dp
  ! Integrate melt rate over the Elmer grid for conservative interpolation purpose
  DO e=1,Solver % NumberOfActiveElements
    Element => GetActiveElement(e)
    CALL GetElementNodes( ElementNodes )
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes
    VisitedNode(NodeIndexes(1:n))=VisitedNode(NodeIndexes(1:n))+1.0_dp
    ! leave the loop if grounded point in the element
    IF ( ANY( GM(GMPerm(NodeIndexes(:))) .GE. mskcrit ) .OR. ANY(Depth(DepthPerm(NodeIndexes(:))) .LT. epsz ) ) CYCLE
    localInteg = 0.0_dp
    localunity = 0.0_dp
    IntegStuff = GaussPoints( Element )
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric,Basis,dBasisdx )
      s = SqrtElementMetric * IntegStuff % s(t)
      localInteg = localInteg + s * SUM(Basis(:) * Melt(MeltPerm(NodeIndexes(:))))
      localunity = localunity + s * SUM(Basis(:))
    END DO
    Integ = Integ + localInteg
    Area  = Area + localunity
  END DO

  IF (Parallel) THEN
    CALL MPI_ALLREDUCE(Integ,Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(Area,Area_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  ELSE
    Integ_Reduced = Integ
    Area_Reduced = Area
  END IF

  AvgTmTf = Integ_Reduced / Area_Reduced

  DO node = 1, Nmax

    IF (ZaUniform) THEN
      zzz = Za
    ELSE
      zzz = Depth(DepthPerm(node))
    ENDIF

    zb = Depth(DepthPerm(node))

    IF ( GM(GMPerm(node)) .LT. mskcrit .AND. zzz .GE. epsz ) THEN
      IF ( zzz .le. z1 ) THEN
        T0 = Tsrf
        S0 = Ssrf
      ELSEIF ( zzz .GE. z2 ) THEN
        T0 = Tbot
        S0 = Sbot
      ELSE
        T0 = (Tbot-Tsrf)*zzz/(z2-z1) + Tsrf - (Tbot-Tsrf)*z1/(z2-z1)
        S0 = (Sbot-Ssrf)*zzz/(z2-z1) + Ssrf - (Sbot-Ssrf)*z1/(z2-z1)
      ENDIF
      Tf = lbd1*S0 + lbd2 + lbd3*zb  ! Sea water freezing temperature
      ! NB: we first store (T0-Tf) in the melt variable :
      Melt(MeltPerm(node)) = - gT * meltfac * AvgTmTf * (T0-Tf)
    ELSE
      Melt(MeltPerm(node)) = 0.0_dp
    ENDIF
  ENDDO

CASE('lazer')

  ! take the local slope
  IF (LazerType .EQ. 1 .OR. LazerType .EQ. 2 .OR. LazerType .EQ. 4) THEN

    isfslopeVar => VariableGet( Model % Mesh % Variables, 'isfslope')
    IF (.NOT.ASSOCIATED(isfslopeVar)) &
    &    CALL FATAL(SolverName,'isfslope not found')

    isfslopePerm => isfslopeVar % Perm
    isfslope => isfslopeVar % Values

  ENDIF

  IF (LazerType .EQ. 1 .OR. LazerType .EQ. 2) THEN
    nbdir = GetInteger( Model % Constants, 'nb dir LAZER',Found )
    IF (.NOT.Found) &
    &  CALL FATAL(SolverName,'<nb dir LAZER> not found')
  ENDIF

  ! Parallel crux
  ! Exchanging informations between partitions
  ! including: Number of GL nodes, GLx and GLy, GL depth
  if (Parallel) then
    ! count the number of GL nodes in this partition 
    ngl=0
    DO node=1,Nmax
      IF (GM(GMPerm(node)) .EQ. 0) THEN
        ngl = ngl+1
      ENDIF
    ENDDO

    debilos = -9999

    ! if nodes in GL, fill up xgl,ygl and dgl
    IF (ngl .GT. 0) then
      ALLOCATE(xgl(ngl),ygl(ngl),dgl(ngl))
      cpt = 1
      DO node = 1, Nmax
        IF (GM(GMPerm(node)) .EQ. 0) THEN
          xgl(cpt) = Mesh % Nodes % x(node)
          ygl(cpt) = Mesh % Nodes % y(node)
          dgl(cpt) = Depth(DepthPerm(node))
          cpt = cpt + 1
        ENDIF
      ENDDO
    ELSE
      ngl = debilos
      ALLOCATE(xgl(1),ygl(1),dgl(1))
      xgl(1) = debilos
      ygl(1) = debilos
      dgl(1) = debilos
    ENDIF

    ! first, send and receive the number of GL nodes between all partitions
    DO i = 1, ParEnv % PEs
      ! need to send info to all partitions but this one
      IF ( ParEnv % MyPE .NE. (i-1) ) THEN
        CALL MPI_SEND(ngl,1,MPI_Integer,i-1,1004,ELMER_COMM_WORLD,ierr)
      ENDIF
    ENDDO

    ngltot=0
    DO i = 1, ParEnv % PEs
      ! need to send info to all partitions but this one
      if ( ParEnv % MyPE .NE. (i-1) ) THEN
        CALL MPI_RECV(nglloc,1,MPI_Integer,i-1,1004,ELMER_COMM_WORLD,status,ierr)
        IF (nglloc .NE. debilos) THEN
          ngltot = ngltot + nglloc
        ENDIF
      ELSE
        IF (ngl .NE. debilos) THEN
          ngltot=ngltot+ngl
        ENDIF
      ENDIF
    ENDDO

    ! ngltot should always be positive whenever you have a grounding line
    ! means ngltot > 0
    ALLOCATE(xgltot(ngltot),ygltot(ngltot),dgltot(ngltot))

    DO i=1,ParEnv % PEs
      IF ( ParEnv % MyPE .NE. (i-1) ) THEN
        CALL MPI_SEND(ngl,1,MPI_Integer,i-1,1006,ELMER_COMM_WORLD,ierr)
        IF (ngl .NE. debilos) THEN
          CALL MPI_SEND(xgl,ngl,MPI_DOUBLE_PRECISION,i-1,1001,ELMER_COMM_WORLD,ierr)
          CALL MPI_SEND(ygl,ngl,MPI_DOUBLE_PRECISION,i-1,1002,ELMER_COMM_WORLD,ierr)
          CALL MPI_SEND(dgl,ngl,MPI_DOUBLE_PRECISION,i-1,1003,ELMER_COMM_WORLD,ierr)
        ELSE
          CALL MPI_SEND(xgl,1,MPI_DOUBLE_PRECISION,i-1,1001,ELMER_COMM_WORLD,ierr)
          CALL MPI_SEND(ygl,1,MPI_DOUBLE_PRECISION,i-1,1002,ELMER_COMM_WORLD,ierr)
          CALL MPI_SEND(dgl,1,MPI_DOUBLE_PRECISION,i-1,1003,ELMER_COMM_WORLD,ierr)
        ENDIF
      ENDIF
    ENDDO

    cpt=1
    DO i=1,ParEnv % PEs
      IF ( ParEnv % MyPE .NE. (i-1) ) THEN

        CALL MPI_RECV(nglloc,1,MPI_Integer,i-1,1006,ELMER_COMM_WORLD,status,ierr)

        IF (nglloc .NE. debilos) THEN

          ALLOCATE(xglloc(nglloc),yglloc(nglloc),dglloc(nglloc))
          CALL MPI_RECV(xglloc,nglloc,MPI_DOUBLE_PRECISION,i-1,1001,ELMER_COMM_WORLD,status,ierr)
          CALL MPI_RECV(yglloc,nglloc,MPI_DOUBLE_PRECISION,i-1,1002,ELMER_COMM_WORLD,status,ierr)
          CALL MPI_RECV(dglloc,nglloc,MPI_DOUBLE_PRECISION,i-1,1003,ELMER_COMM_WORLD,status,ierr)

          xgltot(cpt:cpt+nglloc-1)=xglloc
          ygltot(cpt:cpt+nglloc-1)=yglloc
          dgltot(cpt:cpt+nglloc-1)=dglloc

          cpt=cpt+nglloc

          DEALLOCATE(xglloc,yglloc,dglloc)
        ELSE
          ALLOCATE(xglloc(1),yglloc(1),dglloc(1))
          CALL MPI_RECV(xglloc,1,MPI_DOUBLE_PRECISION,i-1,1001,ELMER_COMM_WORLD,status,ierr)
          CALL MPI_RECV(yglloc,1,MPI_DOUBLE_PRECISION,i-1,1002,ELMER_COMM_WORLD,status,ierr)
          CALL MPI_RECV(dglloc,1,MPI_DOUBLE_PRECISION,i-1,1003,ELMER_COMM_WORLD,status,ierr)
          DEALLOCATE(xglloc,yglloc,dglloc)
        ENDIF
      ELSE
        IF (ngl .NE. debilos) THEN
          xgltot(cpt:cpt+ngl-1)=xgl
          ygltot(cpt:cpt+ngl-1)=ygl
          dgltot(cpt:cpt+ngl-1)=dgl
          cpt=cpt+ngl
        ENDIF
      ENDIF
    ENDDO

  ELSE !if not parallel
    WRITE(Message,'(A)') 'Lazeroms not coded in serial simulations'
    WRITE(Message,'(A)') 'PLease try in parallel'
    CALL FATAL(SolverName,Message)
  ENDIF
  ! End of exhanging informations between partitions

  ! For simple law, get a group of points around the central flowline
  ! that should correspond to the deepest points at the GL
  IF (LazerType .EQ. 3) THEN
    zGLmin = 0.0_dp
    xGLmin = 0.0_dp
    !yGLmin = 40000.0_dp ! CHANGE EVENTUALLY
    yGLmin=(minval(Mesh % Nodes % y)+maxval(Mesh % Nodes % y))/2.0_dp
    cpt = 0
    DO ii = 1, ngltot
      IF (dgltot(ii) .GE. zGLmin) THEN
        zGLmin = dgltot(ii)
        xGlmin = xgltot(ii)
      ENDIF
    ENDDO
  ENDIF

  ! For asymmetric Jenkins law
  IF (LazerType .EQ. 4) THEN

    dy = 1000.0_dp
    ymax = maxval(Mesh % Nodes % y)
    nn = INT(ymax/dy)

    ALLOCATE(xgltot2(nn),ygltot2(nn),dgltot2(nn))

    DO ii=1,nn
      ygltot2(ii) = ymax-(ii-0.5)*dy
      cpt = 0
      xgltot2(ii) = 0.0_dp
      dgltot2(ii) = 0.0_dp
      DO jj=1,ngltot
        IF (ygltot(jj) .LE. ymax-(ii-1)*dy .AND. ygltot(jj) .GE. ymax-ii*dy) THEN
          xgltot2(ii) = xgltot2(ii) + xgltot(jj)
          dgltot2(ii) = dgltot2(ii) + dgltot(jj)
          cpt = cpt + 1
        ENDIF
      ENDDO
      xgltot2(ii) = xgltot2(ii)/cpt
      dgltot2(ii) = dgltot2(ii)/cpt
    ENDDO
  ENDIF

  DO node = 1, Nmax

    zzz = Depth(DepthPerm(node))

    ! NB: this param does not allow melt at the GL
    IF ( GM(GMPerm(node)) .LT. -0.5 .AND. zzz .GE. epsz ) THEN    
      ! Effective slope angle calculated through Solver Compute2DNodalGradient
      ! tan(alpha) = sqrt( (dZb/dx)^2 + (dZb/dy)^2 ) (Lazeroms et al. 2017, Appendix B) :
      ! maximum local slope only accounted for in the AppenB type
      IF (LazerType .EQ. 2) THEN
        alpha = ATAN( SQRT( isfslope(2*(isfslopePerm(node)-1)+1)**2 + isfslope(2*(isfslopePerm(node)-1)+2)**2 ) )
      ENDIF

      ! Calculate the effective grounding line depth (Lazeroms et al. 2017, Appendix B) :
      xP  =  Mesh % Nodes % x(node)
      yP  =  Mesh % Nodes % y(node)

      IF (LazerType .EQ. 3) THEN
        IF (zGLmin .GT. zzz) THEN
          zGL=zGLmin
          alpha = atan(abs(zGL-zzz)/sqrt((xGLmin-xP)**2+(yGLmin-yP)**2))
        ELSE
          zGL=0.0_dp
          alpha=0.0_dp
        ENDIF

      ENDIF

      IF (LazerType .EQ. 1) THEN
        zGL = 0.0_dp
        sn = 0.0_dp
        cptglob = 0
        ! the two local slopes along x and y
        alphax = isfslope(2*(isfslopePerm(node)-1)+1)
        alphay = isfslope(2*(isfslopePerm(node)-1)+2)
        ! define a radius, bigger than the domain size
        radius = (maxval(Mesh % Nodes % x)-minval(Mesh % Nodes % x))*(maxval(Mesh % Nodes % y)-minval(Mesh % Nodes % y))
        angle = 360.0_dp / nbdir * PI / 180.0_dp

        DO ii=1,nbdir

          ! calculate the positions of the two other points in every triangle
          angle1 = (ii-1) * angle
          angle2 = ii * angle

          locslope1 = ( COS(ATAN(alphax)) * SIN(ATAN(alphay)) * SIN(angle1) + COS(ATAN(alphay)) * SIN(ATAN(alphax)) &
                       * COS(angle1) ) / ( COS(ATAN(alphax)) * COS(ATAN(alphay)) )
          locslope2 = ( COS(ATAN(alphax)) * SIN(ATAN(alphay)) * SIN(angle2) + COS(ATAN(alphay)) * SIN(ATAN(alphax)) &
                       * COS(angle2) ) / ( COS(ATAN(alphax)) * COS(ATAN(alphay)) )

          locslope = (locslope1 + locslope2) / 2.0_dp
          if (locslope .GE. 0.0_dp) CYCLE

          xP2 = xP + radius * COS(angle1)
          yP2 = yP + radius * SIN(angle1)
          xP3 = xP + radius * COS(angle2)
          yP3 = yP + radius * SIN(angle2)

          xP12 = xP - xP2
          yP12 = yP - yP2
          xP23 = xP2 - xP3
          yP23 = yP2 - yP3
          xP31 = xP3 - xP
          yP31 = yP3 - yP

          ! now check if GL points lie inside the triangles
          cpttmp = 0
          zGLtmp = 0.0_dp
          DO jj = 1, ngltot
            cond1 = SIGN(1.0_dp,determinant(xP31,yP31,xP23,yP23)) &
                     * SIGN(1.0_dp,determinant(xP3-xgltot(jj),yP3-ygltot(jj),xP23,yP23)) .GE. 0
            cond2 = SIGN(1.0_dp,determinant(xP12,yP12,xP31,yP31)) &
                     * SIGN(1.0_dp,determinant(xP-xgltot(jj),yP-ygltot(jj),xP31,yP31)) .GE. 0
            cond3 = SIGN(1.0_dp,determinant(xP23,yP23,xP12,yP12)) &
                     * SIGN(1.0_dp,determinant(xP2-xgltot(jj),yP2-ygltot(jj),xP12,yP12)) .GE. 0
            IF (cond1 .AND. cond2 .AND. cond3) THEN
              zGLtmp = zGLtmp + dgltot(jj)
              cpttmp = cpttmp + 1
            ENDIF
          ENDDO

          IF (cpttmp .EQ. 0) CYCLE

          zGLloc = zGLtmp / cpttmp

          if (zGLloc .GT. zzz) THEN
            zGL = zGL + zGLloc
            sn = sn + abs(locslope)
            cptglob = cptglob + 1
          ENDIF

        ENDDO

        IF (cptglob .GT. 0) THEN
          zGL = zGL / cptglob
          alpha = ATAN(sn/cptglob)
        ELSE
          zGL = zzz
          alpha = 0.0_dp
        ENDIF

      ENDIF

      IF (LazerType .EQ. 2) THEN
        zGL = 0.0_dp
        div = 0.0_dp
        ! define a radius, bigger than the domain size
        radius = (maxval(Mesh % Nodes % x)-minval(Mesh % Nodes % x))*(maxval(Mesh % Nodes % y)-minval(Mesh % Nodes % y))
        angle = 360.0_dp / nbdir * PI / 180.0_dp

        DO ii = 1, nbdir

          ! calculate the positions of the two other points in every triangle
          angle1 = (ii-1) * angle
          angle2 = ii * angle

          xP2 = xP + radius * COS(angle1)
          yP2 = yP + radius * SIN(angle1)
          xP3 = xP + radius * COS(angle2)
          yP3 = yP + radius * SIN(angle2)

          xP12 = xP - xP2
          yP12 = yP - yP2
          xP23 = xP2 - xP3
          yP23 = yP2 - yP3
          xP31 = xP3 - xP
          yP31 = yP3 - yP

          ! now check if GL points lie inside the triangles
          cpttmp = 0
          zGLtmp = 0.0_dp
          ltmp = 0.0_dp
          DO jj=1,ngltot
            cond1 = SIGN(1.0_dp,determinant(xP31,yP31,xP23,yP23)) &
                     * SIGN(1.0_dp,determinant(xP3-xgltot(jj),yP3-ygltot(jj),xP23,yP23)) .GE. 0
            cond2 = SIGN(1.0_dp,determinant(xP12,yP12,xP31,yP31)) &
                     * SIGN(1.0_dp,determinant(xP-xgltot(jj),yP-ygltot(jj),xP31,yP31)) .GE. 0
            cond3 = SIGN(1.0_dp,determinant(xP23,yP23,xP12,yP12)) &
                     * SIGN(1.0_dp,determinant(xP2-xgltot(jj),yP2-ygltot(jj),xP12,yP12)) .GE. 0
            IF (cond1 .AND. cond2 .AND. cond3) THEN
              zGLtmp = zGLtmp + dgltot(jj)
              ltmp = ltmp + SQRT((xP-xgltot(jj))**2+(yP-ygltot(jj))**2)
              cpttmp = cpttmp + 1
            ENDIF
          ENDDO

          IF (cpttmp .EQ. 0) CYCLE

          zGLloc = zGLtmp / cpttmp

          IF (zGLloc .GT. zzz) THEN
            ltmp = ltmp / cpttmp
            wn = (zGLloc - zzz) / ltmp
            zGL = zGL + wn * zGLloc
            div = div + wn
          ENDIF

        ENDDO

        IF (div .GT. 0.0_dp) THEN
          zGL = zGL / div
        ELSE
          zGL = zzz
          alpha = 0.0_dp
        ENDIF
      ENDIF

      IF (LazerType .EQ. 4) THEN
        ! find the closest GL point
        mindist = (maxval(Mesh % Nodes % x)-minval(Mesh % Nodes % x))*(maxval(Mesh % Nodes % y)-minval(Mesh % Nodes % y))
        kkstart = -9999
        DO kk=1,SIZE(xgltot2)
          dist = SQRT((xP-xgltot2(kk))**2+(yP-ygltot2(kk))**2)
          IF (dist .eq. 0.0) CYCLE
          IF (dgltot2(kk) .GT. zzz .AND. dist .LT. mindist) THEN
            mindist = dist
            kkstart = kk
          ENDIF
        ENDDO

        IF (kkstart .LE. 0) THEN
          zGL = zzz
          alpha = 0.0
        ELSE
          ! go find the GL points from where will start the plume
          halfw = 2
          DO kk=kkstart,SIZE(xgltot2)-halfw
            IF (dgltot2(kk) .GE. dgltot2(kk+1) .AND. dgltot2(kk) .GE. dgltot2(kk+2)) EXIT
          ENDDO
          zGL = dgltot2(kk)
          xGLfin = xgltot2(kk)
          yGLfin = ygltot2(kk)
        ENDIF

        !still need to check if zGL>zzz
        IF (zGL .GE. zzz) THEN
          dist = sqrt((xP-xGLfin)**2+(yP-yGLfin)**2)
          alpha = atan(abs(zzz-zGL)/dist)
        ELSE
          zGL=zzz
          alpha=0.0
        ENDIF

      ENDIF

      ! CALCULATE T0 AND S0 FROM INPUT DATA : TODO

      ! Sea water freezing temperature at the effective grounding line depth :
      Tf = lbd1*S0 + lbd2 + lbd3*zGL  
      ! Effective heat exchange coefficient :
      GTS = GamT * ( 0.545 + 3.5e-5 * (T0-Tf)/lbd3 * E0*sin(alpha)/(GefT+E0*sin(alpha)) )

      ! Melt scale :
      MM = 10. * (T0-Tf)**2 * sqrt( sin(alpha)/(Cd+E0*sin(alpha)) ) * sqrt( GTS/(GTS+E0*sin(alpha)) ) &
               * E0*sin(alpha)/(GTS+E0*sin(alpha))
      ! Length scale :
      ll = (T0-Tf)/lbd3 * (x0*GTS+E0*sin(alpha)) / (x0*(GTS+E0*sin(alpha)))
      ! Dimensionless coordinate :
      X_hat = MAX( 0.0_dp, MIN( 1.0_dp, ( zzz - zGL ) / ll ) )
      ! Dimensionless melt curve :
      M_hat = 0.0_dp

      print *, Size(pp)

      DO kk=1,12
        M_hat = M_hat + pp(kk,1) * X_hat**(kk-1)
      ENDDO

      ! Melt rate in m/yr:
      Melt(MeltPerm(node)) = - K * MM * M_hat

      IF (zzz + Melt(MeltPerm(node)) * dt .LE. epsz) THEN
        Melt(MeltPerm(node)) = 0.0_dp
      ENDIF

    ELSE
      Melt(MeltPerm(node)) = 0.0_dp
    ENDIF

  ENDDO

CASE('pico')

  distGLVar => VariableGet( Model % Mesh % Variables, 'distGL')
  IF (.NOT.ASSOCIATED(distGLVar)) &
  &    CALL FATAL(SolverName,'distGL not found')

  distIFVar => VariableGet( Model % Mesh % Variables, 'distIF')
  IF (.NOT.ASSOCIATED(distIFVar)) &
  &    CALL FATAL(SolverName,'distIF not found')

  distGLPerm => distGLVar % Perm
  distGL => distGLVar % Values

  distIFPerm => distIFVar % Perm
  distIF => distIFVar % Values

  ALLOCATE( Zbox(nD), Abox(nD), Tbox(nD), Sbox(nD), Mbox(nD) )

  !- Ambiant temperature and salinity (at depth Za)
  IF ( Za .LE. z1 ) THEN
    T0 = Tsrf
    S0 = Ssrf
  ELSEIF ( Za .GE. z2 ) THEN
    T0 = Tbot
    S0 = Sbot
  ELSE
    T0 = (Tbot-Tsrf)*Za/(z2-z1) + Tsrf - (Tbot-Tsrf)*z1/(z2-z1)
    S0 = (Sbot-Ssrf)*Za/(z2-z1) + Ssrf - (Sbot-Ssrf)*z1/(z2-z1)
  ENDIF

  !- Calculate total area and mean depth of each box :
  Abox(:) = 0.0_dp
  Zbox(:) = 0.0_dp
  DO e = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(e)
    CALL GetElementNodes( ElementNodes )
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes
    VisitedNode(NodeIndexes(1:n)) = VisitedNode(NodeIndexes(1:n)) + 1.0_dp
    ! leave the loop if grounded point in the element

    IF ( ANY( GM(GMPerm(NodeIndexes(:))) .GE. mskcrit ) ) CYCLE

    rr = SUM(distGL(distGLPerm(NodeIndexes(:)))    &
    &             / ( distGL(distGLPerm(NodeIndexes(:))) + distIF(distIFPerm(NodeIndexes(:))) )) &
    &    / MAX(1,SIZE(distGL(distGLPerm(NodeIndexes(:)))))
    localInteg = 0.0_dp
    localunity = 0.0_dp
    IntegStuff = GaussPoints( Element )
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
           Basis,dBasisdx )
      s = SqrtElementMetric * IntegStuff % s(t)
      localInteg = localInteg + s * SUM(Basis(:)*Depth(DepthPerm(NodeIndexes(:))))
      localunity = localunity + s * SUM(Basis(:))
    ENDDO
    DO kk=1,nD
      IF ( rr .ge. 1.0-sqrt(1.0*(nD-kk+1)/nD) .and. rr .le. 1.0-sqrt(1.0*(nD-kk)/nD) ) THEN
        Zbox(kk) = Zbox(kk) + localInteg
        Abox(kk) = Abox(kk) + localunity
      ENDIF
    ENDDO
  END DO
  DO kk=1,nD
    IF (Parallel) THEN
      CALL MPI_ALLREDUCE(Zbox(kk),Integ_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(Abox(kk),Area_Reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      Zbox(kk) = Integ_Reduced
      Abox(kk) = Area_Reduced
    END IF
  ENDDO

  Zbox(:) = Zbox(:) / Abox(:) ! mean depth of each box

  !- Temperature and salinity in Box #1 :
  Tstar = lbd1*S0 + lbd2 + lbd3*Zbox(1) - T0
  g1 = Abox(1) * gT 
  tmp1 = g1 / (CC*rhostar*(beta*S0*meltfac-alpha))
  xbox = - 0.5*tmp1 + sqrt( (0.5*tmp1)**2 - tmp1*Tstar )
  Tbox(1) = T0 - xbox
  Sbox(1) = S0 - xbox*S0*meltfac
  qqq = CC*rhostar*(beta*(S0-Sbox(1))-alpha*(T0-Tbox(1)))
  Mbox(1) = - gT * meltfac * ( lbd1*Sbox(1) + lbd2 + lbd3*Zbox(1) - Tbox(1) )

  !- Temperature and salinity in possible other boxes :
  DO kk = 2, nD
    Tstar = lbd1*Sbox(kk-1) + lbd2 + lbd3*Zbox(kk) - Tbox(kk-1)
    g1  = Abox(kk) * gT
    g2  = g1 * meltfac
    xbox = - g1 * Tstar / ( qqq + g1 - g2*lbd1*Sbox(kk-1) )
    Tbox(kk) = Tbox(kk-1) - xbox
    Sbox(kk) = Sbox(kk-1) - xbox*Sbox(kk-1)*meltfac
    Mbox(kk) = - gT * meltfac * ( lbd1*Sbox(kk) + lbd2 + lbd3*Zbox(kk) - Tbox(kk) )
  ENDDO

  !- Attribute melt at each node to a box value :
  DO node = 1, Nmax

    zzz=Depth(DepthPerm(node))

    IF ( GM(GMPerm(node)) .LT. mskcrit .AND. zzz .GE. epsz ) THEN

      rr = distGL(distGLPerm(node)) / ( distGL(distGLPerm(node)) + distIF(distIFPerm(node)) )

      DO kk = 1, nD

        IF ( rr .GE. 1.0-sqrt(1.0*(nD-kk+1)/nD) .AND. rr .LE. 1.0-sqrt(1.0*(nD-kk)/nD) ) THEN

          Melt(MeltPerm(node)) = - Mbox(kk)
          IF (zzz + Melt(MeltPerm(node)) * dt .LE. epsz) THEN
            Melt(MeltPerm(node)) = 0.0_dp
          ENDIF

        ENDIF
      ENDDO
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

!!!
CONTAINS

  ! For Lazeroms params
  FUNCTION determinant(x1,y1,x2,y2) RESULT(det)
    REAL(KIND=dp), INTENT(IN) :: x1, y1, x2, y2
    REAL(KIND=dp) :: det

    det=x1*y2-x2*y1

  END FUNCTION determinant

END SUBROUTINE MISMIP_Melt_Param


