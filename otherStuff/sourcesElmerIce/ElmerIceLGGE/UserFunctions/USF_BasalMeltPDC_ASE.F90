! Pollard, D., DeConto, R.M., 2012. Description of a hybrid ice
! sheet-shelf model, and application to Antarctica. Geoscientific Model
! Development 5, 1273–1295. doi:10.5194/gmd-5-1273-2012
!
! Pas de modification en fn de l'angle d'open ocean (Eq.18)
! Basé sur PDC_BASALMELT.F90 de Julien Brondex
! Considérant seulement la mer d'Amudsen
! The formula doesn't calculate the Tf
! But gives dt = To-Tf as a function of depth
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION BasalMeltPDC_ASE(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(2),VarOut
  REAL(kind=dp) :: zb,mask
  REAL(kind=dp) :: depth,Tf,T0,K,dt
  REAL(kind=dp) :: melt
  LOGICAL :: Found
  REAL(kind=dp) :: sealevel,rw,cw,Lf,ri,Kt
  CHARACTER(LEN=MAX_NAME_LEN) :: USFName='USF_BasalMeltPDC_ASE'

  mask = VarIn(2)

  IF (mask .GT. -0.5) THEN 
    VarOut = 0._dp
  ELSE
    sealevel = GetCReal( Model % Constants, 'Sea Level', Found )
    IF (.NOT.Found) &
      CALL FATAL(USFName,'<Sea Level> not found')

    rw=GetCReal( Model % Constants, 'water density', Found )
    IF (.NOT.Found) &
      CALL FATAL(USFName,'<water density> not found')

    cw=GetCReal( Model % Constants, 'Sea Water Specific heat',Found )
    IF (.NOT.Found) &
      CALL FATAL(USFName,'<Sea Water Specific heat> not found')

    Lf=GetCReal( Model % Constants, 'Ice fusion latent heat',Found )
    IF (.NOT.Found) &
      CALL FATAL(USFName,'<Ice fusion latent heat> not found')

    ri=GetCReal( Model % Constants, 'Ice density',Found )
    IF (.NOT.Found) &
      CALL FATAL(USFName,'<Ice density> not found')

    kt=GetCReal( Model % Constants, 'Transfer Factor',Found )
    IF (.NOT.Found) &
      CALL FATAL(USFName,'<Transfer Factor> not found')

    zb = VarIn(1)
    depth = sealevel-zb
    if (depth .LT. 0._dp) THEN 
      VarOut=0._dp
      RETURN
    ENDIF

    IF (depth .LT. 170.0) THEN
      dt = 0.5
    ELSEIF (depth .GT. 680.0) THEN
      dt = 3.5
    ELSE
      dt = 0.5 + (depth-170.0)*3.0/510.0
    ENDIF
    K = 8.0

    melt = K*Kt*rw*cw*abs(dt)*dt/(ri*Lf)

    VarOut = -melt

  END IF

END FUNCTION BasalMeltPDC_ASE 

