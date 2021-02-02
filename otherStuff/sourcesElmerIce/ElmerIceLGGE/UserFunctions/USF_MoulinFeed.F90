!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Made up for the HydroSolver, compute a recharge !
!due to a background flux and a number of 
!discrete recharge fluxes
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION MoulinFeed ( Model, nodenumber, x ) RESULT(InFlow)

  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE
  !-------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(KIND=dp) :: x, InFlow 
  !-------------------------------------------------
  TYPE(Variable_t), POINTER:: TimeVar
  TYPE(ValueList_t), POINTER :: BodyForce
  TYPE(Element_t), POINTER :: Element, BoundElement, CurElt
  TYPE(Nodes_t) :: ElementNodes

  REAL(KIND=DP), ALLOCATABLE :: MoulinPosition(:,:), & !position of each moulin read in file(NM*DIM)
       Distance(:,:), & !Distance of each moulin to the points of an element (NM*DIM)
       MoulinIndex(:,:), & !Nodenumber of each moulin and their distance to the given coordinates(NM*2)
       MoulinFlux(:,:,:), & !Flux injected through each moulin at each given time (NM*NtFmax*2)
       FuncTime(:,:), & !time at and after which a given function is applied
       AuxReal(:), LocalArea(:)
  REAL(KIND=DP) :: tps, ratio, MoulinFlow

  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: bf_id, DIM, N
  INTEGER :: NM ! Number of Moulins
  INTEGER :: NtFmax, NtFmin, NtFlux!Max of time for flux
  INTEGER :: i, j, k

  LOGICAL :: FirstTime = .TRUE., AllocationsDone = .FALSE., Found

  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: MoulinFluxFile(:), &! name of the files containig the flux (one file per moulin)
       MoulinFunc(:,:) ! tab containing the flux function for each moulin which needs it
  CHARACTER(LEN=MAX_NAME_LEN) :: MoulinFile, & !name of the file containing moulin position and moulin flux file names
       InputDir, &!name of the input directory
       cmd, tmp_str ! used for the matc Call
 
  SAVE MoulinIndex, MoulinFlux
  SAVE FuncTime, MoulinFunc, NM
  SAVE NtFmax, NtFmin
  SAVE FirstTime, AllocationsDone, LocalArea


  Timevar => VariableGet( Model % Variables,'Time')
  tps = TimeVar % Values(1)
  Element => Model % CurrentElement
  NodeIndexes => Element % NodeIndexes
  N = GetElementNOFNodes()
  IF (ALLOCATED(AuxReal)) DEALLOCATE(AuxReal)
  ALLOCATE(AuxReal(N))

  BodyForce => GetBodyForce(Element) 

  IF ( ASSOCIATED( BodyForce ) ) THEN
     bf_id = GetBodyForceId()
     AuxReal(1:N) = ListGetReal (Model % BodyForces(bf_id) % Values,&
          'Water Recharge Background', N, Element % NodeIndexes, Found )
     DO i=1, N
        IF (NodeNumber .EQ. Element % NodeIndexes(i)) EXIT 
     END DO
     InFlow = AuxReal(i)
     MoulinFlow = 0.0
  END IF

  IF (FirstTime) THEN
     FirstTime = .FALSE.
     !-------------------------------------------------
     !Get the name of the moulin position file
     !-------------------------------------------------

     DIM = CoordinateSystemDimension()
     N = Model % MaxElementNodes
     CurElt => Model % CurrentElement

     BodyForce => GetBodyForce() 
     IF ( ASSOCIATED( BodyForce ) ) THEN
        bf_id = GetBodyForceId()

        InputDir = ListGetString(Model % BodyForces(bf_id) % Values,&
             'Moulin Input Directory', Found)
        IF ( .NOT.Found )CALL FATAL('USF_MoulinFeed', 'No Moulin Input Directory found')

        MoulinFile = ListGetString(Model % BodyForces(bf_id) % Values,&
             'Moulin Position File', Found)
        IF ( .NOT.Found )CALL FATAL('USF_MoulinFeed', 'No Moulin Position File found')

     END IF

     !-------------------------------------------------
     !Get the number of moulins, their position and 
     !flux file
     !-------------------------------------------------
     OPEN(12,file=TRIM(InputDir)//"/"//TRIM(MoulinFile))
     READ(12,fmt="(i6)")NM

     IF ( AllocationsDone ) THEN

        DEALLOCATE( MoulinPosition, &
             MoulinIndex             , &
             MoulinFluxFile          , &
             Distance                , &
             LocalArea)
     END IF

     ALLOCATE(MoulinPosition(NM,DIM), &
          MoulinIndex (NM,3)        , &
          MoulinFluxFile(NM)        , &
          Distance (NM,N)           , &
          LocalArea(N))

     AllocationsDone = .TRUE.

     !Reading position of the moulins depending of the problem dimension
     !------------------------------------------------------------------

     SELECT CASE (DIM)
     CASE (1)
        READ(12,*)(MoulinPosition(i,1), MoulinFluxFile(i), i=1,NM)
     CASE (2)
        READ(12,*)(MoulinPosition(i,1), MoulinPosition(i,2), &
             MoulinFluxFile(i), i=1,NM)
     CASE (3)
        READ(12,*)(MoulinPosition(i,1), MoulinPosition(i,2), &
             MoulinPosition(i,3),MoulinFluxFile(i), i=1,NM)
     END SELECT

     CLOSE(12)

     !-------------------------------------------------
     !Get nodenumbers where to inject the well input, 
     !nodenumbers that are the closest to the given
     !coordinates
     !-------------------------------------------------

     DO i=1, Model % NumberOfBoundaryElements

        BoundElement => GetBoundaryElement(i)
        N = GetElementNOFNodes(BoundElement)
        NodeIndexes => BoundElement % NodeIndexes
        CALL GetElementNodes(ElementNodes,BoundElement) 

        CALL GetScalarLocalSolution(LocalArea,'hwater Boundary Weights')
        write(*,*)LocalArea(1:N)

        !Compute distance between each node and moulin
        !---------------------------------------------

        SELECT CASE (DIM)
        CASE (1)
           DO j=1,N
              Distance (:,j) = ((MoulinPosition(:,1)-ElementNodes % x(j))**2.0)**0.5  

              !keep indexes of the closer points to the given moulins
              !------------------------------------------------------
              IF (i+j.EQ.2)THEN
                 !CALL GetScalarLocalSolution(LocalArea,'Area')
                 MoulinIndex(:,2) = Distance(:,1)
                 MoulinIndex(:,1) = BoundElement % NodeIndexes(1)
                 MoulinIndex(:,3) = LocalArea(1)
                 
              ELSE
              
                 WHERE(Distance(:,j).LT.MoulinIndex(:,2))
                    !CALL GetScalarLocalSolution(LocalArea,'Area')
                    MoulinIndex(:,1) = BoundElement % NodeIndexes(j)
                    MoulinIndex(:,2) = Distance(:,j)
                    MoulinIndex(:,3) = LocalArea(j)
                 END WHERE
              END IF
           END DO
           
        CASE (2)
           DO j=1,N
              Distance (:,j) = ((MoulinPosition(:,1)-ElementNodes % x(j))**2.0 &
                   + (MoulinPosition(:,2)-ElementNodes % y(j))**2.0)**0.5  
              
              !keep indexes of the closer points to the given moulins
              !------------------------------------------------------
              IF (i+j.EQ.2)THEN
                 !CALL GetScalarLocalSolution(LocalArea,'Area')
                 MoulinIndex(:,2) = Distance(:,1)
                 MoulinIndex(:,1) = BoundElement % NodeIndexes(1)
                 MoulinIndex(:,3) = LocalArea(1)
                 
              ELSE
              
                 WHERE(Distance(:,j).LT.MoulinIndex(:,2))
                    !CALL GetScalarLocalSolution(LocalArea,'Area')
                    MoulinIndex(:,1) = BoundElement % NodeIndexes(j)
                    MoulinIndex(:,2) = Distance(:,j)
                    MoulinIndex(:,3) = LocalArea(j)
                 END WHERE
              END IF
           END DO

        CASE (3)
           DO j=1,N
              Distance (:,j) =((MoulinPosition(:,1)-ElementNodes % x(j))**2.0 &
                   + (MoulinPosition(:,2)-ElementNodes % y(j))**2.0 &
                   + (MoulinPosition(:,3)-ElementNodes % z(j))**2.0)**0.5

              !keep indexes of the closer points to the given moulins
              !------------------------------------------------------
              IF (i+j.EQ.2)THEN
                 !CALL GetScalarLocalSolution(LocalArea,'Area')
                 MoulinIndex(:,2) = Distance(:,1)
                 MoulinIndex(:,1) = BoundElement % NodeIndexes(1)
                 MoulinIndex(:,3) = LocalArea(1)
                 
              ELSE
              
                 WHERE(Distance(:,j).LT.MoulinIndex(:,2))
                    !CALL GetScalarLocalSolution(LocalArea,'Area')
                    MoulinIndex(:,1) = BoundElement % NodeIndexes(j)
                    MoulinIndex(:,2) = Distance(:,j)
                    MoulinIndex(:,3) = LocalArea(j)
                 END WHERE
              END IF
           END DO
        END SELECT
  
     END DO

     !-------------------------------------------------
     !Open Moulin injection files and store the 
     !injection volumes in an array
     !-------------------------------------------------
     DO i=1, NM
        !NtFlux is the number of timesteps given for the injection
        !NtFmax is the max of NtFluxes
        !---------------------------------------------------------

        OPEN(11,file=TRIM(InputDir)//"/"//TRIM(MoulinFluxFile(i)))
        READ(11,fmt="(i6)")NtFlux
        IF (i.EQ.1)THEN
           NtFmax = INT(NtFlux)
           NtFmin = NtFlux
        ELSE
           IF (NtFlux.GT.NtFmax)NtFmax = NtFlux
           IF (NtFlux.LT.NtFmin)NtFmin = NtFlux
        END IF
        CLOSE(11)
     END DO
     

     IF ((NtFmax.GT.0).AND.(NtFmin.GE.0))THEN
        NtFmin = -1
        IF (ALLOCATED(MoulinFlux))DEALLOCATE(MoulinFlux)  
        ALLOCATE(MoulinFlux(NM,NtFmax,2)) !number of moulin * Max number of timesteps * 2 (time and input)
        IF (ALLOCATED(MoulinFunc))DEALLOCATE(MoulinFunc, FuncTime)
        ALLOCATE(MoulinFunc(NM,ABS(NtFmin)), &
             FuncTime(NM,ABS(NtFmin)))

     ELSEIF((NtFmax.GT.0).AND.(NtFmin.LT.0))THEN
        IF (ALLOCATED(MoulinFlux))DEALLOCATE(MoulinFlux)  
        ALLOCATE(MoulinFlux(NM,NtFmax,2)) 
        IF (ALLOCATED(MoulinFunc))DEALLOCATE(MoulinFunc, FuncTime)  
        ALLOCATE(MoulinFunc(NM,ABS(NtFmin)), & !number of moulin * Max number of function 
             FuncTime(NM,ABS(NtFmin)))        !number of moulin * Max number of function (at and after given time execute given function)

     ELSEIF((NtFmax.LE.0).AND.(NtFmin.LT.0))THEN
        NtFmax = 1
        IF (ALLOCATED(MoulinFunc))DEALLOCATE(MoulinFunc, FuncTime)  
        ALLOCATE(MoulinFunc(NM,ABS(NtFmin)), & 
             FuncTime(NM,ABS(NtFmin)))
        IF (ALLOCATED(MoulinFlux))DEALLOCATE(MoulinFlux)
        ALLOCATE(MoulinFlux(NM,NtFmax,2))  
     
     END IF 

     DO i=1, NM
        OPEN(10,file=TRIM(InputDir)//"/"//TRIM(MoulinFluxFile(i)))
        READ(10,fmt="(i6)")NtFlux
        
        !test pour lecture de fonction (matc)
        IF (NtFlux.LE.0)THEN
           !we have to read matc function
           READ(10,fmt="(e12.5,a)")(FuncTime(i,j), MoulinFunc(i,j), j=1,ABS(NtFlux))
           CLOSE(10)

           !complete the empty tab with null value (-9999.9999)
           IF (NtFlux.GT.NtFmin)FuncTime(i,ABS(NtFlux)+1:ABS(NtFmin)) = -9999.9999
           !Moulin Flux is Dummy fill with dummy code (9999.9999)
           MoulinFlux(i,:,:) = 9999.9999
           
        ELSE
           READ(10,*)(MoulinFlux(i,j,1), MoulinFlux(i,j,2), j=1,NtFlux)
           CLOSE(10)
           
           !complete the empty tab with null value (-9999.9999)
           IF (NtFlux.LT.Ntfmax)MoulinFlux(i,NtFlux+1:NtFmax,:) = -9999.9999
           !FuncTime is Dummy fill with dummy code (9999.9999)
           FuncTime(i,:) = 9999.9999
        END IF
     END DO
     Model % CurrentElement=>CurElt

  END IF !end of the initialisation

  !-------------------------------------------------
  !If we are on a moulin location
  !-------------------------------------------------
  DO i=1,NM
     IF(nodenumber.EQ.MoulinIndex(i,1))THEN

        !-------------------------------------------------
        !Compute flux at current Time
        !-------------------------------------------------

        !If the flux is constant through time
        !------------------------------------
        IF (MoulinFlux(i,2,1).EQ.-9999.9999)THEN 
           MoulinFlow = MoulinFlux(i,1,2)
           MoulinFlow = MoulinFlow / MoulinIndex(i,3)
        !If the flux is given by a function
        !------------------------------------
        ELSEIF (MoulinFlux(i,1,1).EQ.9999.9999)THEN 
           DO j=1,ABS(NtFmin)-1


              IF ((tps.GE.FuncTime(i,j)).AND.(tps.LE.FuncTime(i,j+1)))THEN
                 
                 WRITE( cmd, '(a,e15.8)' ) 'tx = ', tps
                 k = LEN_TRIM(cmd)
                 CALL matc( cmd, tmp_str, k )
  
                 cmd = MoulinFunc(i,j)
                 k = LEN_TRIM(cmd)
                 CALL matc( cmd, tmp_str, k )
                 READ( tmp_str(1:k), * ) MoulinFlow
              END IF
           END DO
           MoulinFlow = MoulinFlow / MoulinIndex(i,3)

        ELSEIF((tps.GT.Maxval(MoulinFlux(i,:,1))).OR.(tps.LT.Minval(MoulinFlux(i,:,1))))THEN 
           write(*,*)"Current time (", tps, ") out of bound of given entering fluxes times (",&
                MoulinFlux(i,1,1),"-",MoulinFlux(i,NtFmax,1), ") for Moulin", i
           STOP

       
           !iteration on times given by the user
           !------------------------------------
        ELSE
           DO j=1,NtFmax-1

              IF ((tps.GE.MoulinFlux(i,j,1)).AND.(tps.LE.MoulinFlux(i,j+1,1)))THEN
                 !linear interpolation of the flux
                 !--------------------------------
                 ratio = (tps - MoulinFlux(i,j,1))/(MoulinFlux(i,j+1,1) - MoulinFlux(i,j,1))
                 MoulinFlow = ((1.0 - ratio) * MoulinFlux(i,j,2) + (ratio) * MoulinFlux(i,j+1,2)) 
                 MoulinFlow = MoulinFlow / MoulinIndex(i,3)
              END IF
           END DO
        END IF
        
        write(42,*)tps, i, MoulinFlow
     END IF
  END DO
  

InFlow = InFlow + MoulinFlow

END FUNCTION MoulinFeed
