!
! Clean a contour file by deleting multiple points with same coordinates
!
   PROGRAM  CleanContour
       USE types
       USE DefUtils
!
IMPLICIT NONE
!
REAL(KIND=dp), ALLOCATABLE  :: x(:), y(:)                       
REAL(KIND=dp)  :: ux, uy, vx, vy, kx, ky                        
REAL(KIND=dp)  :: Eps = 1.0e-3_dp, big = 1.0e20_dp                          
CHARACTER :: NameFile*20
INTEGER :: Ni, Nf, i, j, k 
INTEGER, ALLOCATABLE :: Pt(:) 
!
!
!
      WRITE(*,*)'Name of the contour file ?'
      READ(*,*)NameFile

      OPEN(10,file=TRIM(NameFile)//".dat")
      READ(10,1000)Ni
      ALLOCATE(x(Ni), y(Ni))
      READ(10,*)(x(i), y(i), i=1,Ni)
      CLOSE(10)
      

DO k = 1, 10
      ALLOCATE(Pt(Ni))
      Pt = 1
      DO i=2,Ni-1
        ux = (x(i)-x(i-1))
        uy = (y(i)-y(i-1))
        vx = (x(i+1)-x(i))
        vy = (y(i+1)-y(i))
        
      
        IF (ABS(vx)>Eps) THEN
          kx = ux / vx
        ELSE IF (ABS(ux) > Eps) THEN
          kx = big * 1.2_dp
        ELSE
          kx = -big * 1.2_dp 
        END IF

        IF (ABS(vy)>Eps) THEN
          ky = uy / vy
        ELSE IF (ABS(uy)> Eps) THEN
          ky = big * 1.2_dp
        ELSE
          ky = -big * 1.2_dp
        END IF
        
        IF (ky < -big) ky=kx
        IF (kx < -big) kx=ky
        
        IF (ABS(kx-ky)>Eps) THEN 
           CYCLE
        ELSE
  ! kx=ky -> points are aligned
    ! i-1, i, i+1 are ordered
           IF (kx > 0.) CYCLE
  ! Pts i-1 and i are identical, delete i-1
           IF (ABS(kx)<Eps) THEN 
              Pt(i-1) = 0
              CYCLE
  ! Pts i and i+1 are identical, delete i
           ELSE IF (ABS(kx)>big) THEN
              Pt(i)=0
              CYCLE
           END IF
           IF (kx < -1.) THEN
              Pt(i)=0
           ELSE IF (kx > -1.) THEN
              Pt(i)=0
           ELSE IF (ABS(kx+1.)<Eps) THEN
              Pt(i-1)=0
              Pt(i)=0
           END IF
        END IF
      END DO

      Nf = SUM(Pt)
      j = 0
      DO i=1, Ni
         IF (Pt(i)==1) THEN
           j = j + 1
           x(j) = x(i) 
           y(j) = y(i)
         END IF
      END DO

      WRITE(*,*)Ni-Nf,' Points have been deleted for loop ', k 
      IF ((Ni-Nf)==0) EXIT 
      Ni = Nf
      DEALLOCATE(Pt)
END DO 

      OPEN(10,file=TRIM(NameFile)//"_new.dat")
      WRITE(10,1000)Nf
      DO i=1, Nf
         WRITE(10,*)x(i),y(i)
      END DO
      CLOSE(10)

      WRITE(*,*)'END WITH NO TROUBLE ...'
!
1000 FORMAT(I6)

End Program CleanContour

