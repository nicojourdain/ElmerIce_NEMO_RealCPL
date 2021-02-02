!
! Create a 2D mesh given by its bed and surface topography 
! By deforming a square 1x1 mesh
!
! Input :
!    NameMesh of the grd mesh (a square 1 x 1)
!    x_start end x_end 
!    yi=surf(xi) and yi=bed(xi) are in NameMesh_surf.dat and NameMesh_bed.dat
!
!  Maillage type ELMER 
!
! Output :
!    NameMesh/mesh.nodes
!
!
   PROGRAM  MSH_Glacier 
       USE types
       USE DefUtils
!
IMPLICIT NONE
!
REAL(KIND=dp)  ::   x0, x1, x, y, z, xnew, ynew, zb, zs, Xconstraint
REAL(KIND=dp), ALLOCATABLE :: xbed(:), ybed(:), xsurf(:), & 
                       ysurf(:), y2s(:), y2b(:)
REAL(KIND=dp), ALLOCATABLE  :: xnode(:), ynode(:)                       
CHARACTER :: NameMsh*20
INTEGER :: NtN, i, j, NptS, NptB, n, NtNx, xi 
!
       INTERFACE 
         SUBROUTINE spline(x,y,yp1,ypn,y2) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x,y 
         REAL(KIND=dp), INTENT(IN) :: yp1,ypn 
         REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: y2 
         END SUBROUTINE spline
         FUNCTION splint(xa,ya,y2a,x) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xa,ya,y2a 
         REAL(KIND=dp), INTENT(IN) :: x 
         REAL(KIND=dp) :: splint
         END FUNCTION splint
       END INTERFACE         
!
!
      WRITE(*,*)'Name of the Elmer mesh directory ?'
      READ(*,*)NameMsh
      WRITE(*,*)'The mesh x-coordinates correspond to the data? (Yes=1, No=0)'
      READ(*,*)XConstraint
      IF (XConstraint < 1.0) THEN
         WRITE(*,*)'X start :'
         READ(*,*)x0
         WRITE(*,*)'X end :'
         READ(*,*)x1
      END IF

      OPEN(10,file=TRIM(NameMsh)//"/mesh.header")
        READ(10,1000)NtN
      CLOSE(10)
      ALLOCATE(xnode(NtN), ynode(NtN))

      OPEN(10,file=TRIM(NameMsh)//"_surf.dat")
        READ(10,1000)NptS
        ALLOCATE(xsurf(NptS), ysurf(NptS), y2s(NptS))
        READ(10,*)(xsurf(i), ysurf(i), i=1,NptS)
      CLOSE(10)

      OPEN(10,file=TRIM(NameMsh)//"_bed.dat")
        READ(10,1000)NptB
        ALLOCATE(xbed(NptB), ybed(NptB), y2b(NptB))
        READ(10,*)(xbed(i), ybed(i), i=1,NptB)
      CLOSE(10)
      
      IF ((XConstraint < 1.0) .AND. ((MINVAL(xbed)>x0) .OR. (MAXVAL(xbed)<x1)) ) THEN
         WRITE(*,*)'MUST BE : x0 > MIN(xbed) AND x1 < MAX(xbed)'
         STOP
      END IF

    
      CALL spline(xsurf, ysurf, 0.0_dp, 0.0_dp, y2s) 
      CALL spline(xbed, ybed, 0.0_dp, 0.0_dp, y2b) 
      
      OPEN(11,file=TRIM(NameMsh)//"/mesh.nodes0")
      OPEN(12,file=TRIM(NameMsh)//"/mesh.nodes")
      READ(12,*)(N, j, xnode(i), ynode(i), z, i=1,NtN)
      write(*,*)z
!     WRITE(11,1200)(i, j, xnode(i), ynode(i), z, i=1,NtN)
      CLOSE(11)
      REWIND(12)


      IF (XConstraint > 0.0) THEN
              WRITE(*,*)'Will have to check if number of nodes in x direction &
                    &equals number of data points'
! First tests
          IF (NptS/=NptB) THEN
                  WRITE(*,*)'Surface and Bed data must have the same number of &
                  &points'
                  STOP
          END IF
          NtNx = 0
          DO n=1,NtN
             IF (ABS(ynode(n))<1.0e-4) NtNx = NtNx + 1 
          END DO
          IF (NptS/=NtNx) THEN
                  WRITE(*,*)'Mesh must have',NptS,' nodes in x direction' 
                  STOP
          END IF
          x0 = MINVAL(xbed)
          x1 = MAXVAL(xbed)
      END IF

      DO n=1, NtN

        x = xnode(n)
        y = ynode(n)


        IF (Xconstraint > 0.0) THEN
           xi = ANINT(x / (1.0_dp / (NptS-1.0_dp))) + 1
           x = (xsurf(xi)-x0)/(x1-x0)
        END IF

        xnew = x0 + x * (x1 - x0) 
        zs = splint(xsurf, ysurf, y2s, xnew) 
        zb = splint(xbed, ybed, y2b, xnew) 
        ynew = zb + y * (zs - zb)
        
        WRITE(12,1200)N,j,xnew,ynew,z
      END DO
      WRITE(*,*)'END WITH NO TROUBLE ...'
!
1000 FORMAT(I6)
1200 FORMAT(i5,2x,i5,3(2x,e22.15)) 

End Program MSH_Glacier

! --------------------------------------------------------
FUNCTION splint(xa,ya,y2a,x) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xa,ya,y2a 
REAL(KIND=dp), INTENT(IN) :: x 
REAL(KIND=dp) :: splint
INTEGER :: khi,klo,n
REAL(KIND=dp) :: a,b,h 
n = size(xa)
klo=max(min(locate(xa,x),n-1),1)
khi=klo+1 
h=xa(khi)-xa(klo) 
if (h == 0.0) THEN
    Write(*,*)' bad xa input in splint'  
    Write(*,*)khi,klo,xa(khi),xa(klo),x
End IF
a=(xa(khi)-x)/h 
b=(x-xa(klo))/h 
splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp 
CONTAINS
! ----------------------------------------------------------------------------

FUNCTION locate(xx,x) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xx 
REAL(KIND=dp), INTENT(IN) :: x 
INTEGER :: locate 
INTEGER :: n,jl,jm,ju 
LOGICAL :: ascnd
n=size(xx) 
ascnd = (xx(n) >= xx(1)) 
jl=0 
ju=n+1  
do 
if (ju-jl <= 1) exit 
jm=(ju+jl)/2 
if (ascnd .eqv. (x >= xx(jm))) then 
jl=jm 
else 
ju=jm 
end if 
end do 
if (x == xx(1)) then 
locate=1 
else if (x == xx(n)) then 
locate=n-1 
else 
locate=jl 
end if 
END FUNCTION locate

END FUNCTION splint
! ----------------------------------------------------------------------------

SUBROUTINE spline(x,y,yp1,ypn,y2) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x,y 
REAL(KIND=dp), INTENT(IN) :: yp1,ypn 
REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: y2 
INTEGER :: n 
REAL(KIND=dp), DIMENSION(size(x)) :: a,b,c,r 
INTERFACE
SUBROUTINE tridag(a,b,c,r,u) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: a,b,c,r 
REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: u 
END SUBROUTINE tridag
END INTERFACE
n = size(x)
c(1:n-1)=x(2:n)-x(1:n-1) 
r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1)) 
r(2:n-1)=r(2:n-1)-r(1:n-2) 
a(2:n-1)=c(1:n-2) 
b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1)) 
b(1)=1.0 
b(n)=1.0 
if (yp1 > 0.99e30_dp) then 
r(1)=0.0 
c(1)=0.0 
else 
r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
c(1)=0.5 
end if 
if (ypn > 0.99e30_dp) then 
r(n)=0.0 
a(n)=0.0 
else 
r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn) 
a(n)=0.5 
end if 
call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n)) 
END SUBROUTINE spline

! ----------------------------------------------------------------------------

SUBROUTINE tridag(a,b,c,r,u) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: a,b,c,r 
REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: u 
REAL(KIND=dp), DIMENSION(size(b)) :: gam
INTEGER :: n,j 
REAL(KIND=dp) :: bet 
n = size(b)
bet=b(1) 
if (bet == 0.0) Write(*,*)' tridag_ser: Error at code stage 1 ' 
u(1)=r(1)/bet 
do j=2,n 
gam(j)=c(j-1)/bet 
bet=b(j)-a(j-1)*gam(j) 
if (bet == 0.0) Write(*,*)' tridag_ser: Error at code stage 2 ' 
u(j)=(r(j)-a(j-1)*u(j-1))/bet 
end do 
do j=n-1,1,-1 
u(j)=u(j)-gam(j+1)*u(j+1) 
end do 
END SUBROUTINE tridag

