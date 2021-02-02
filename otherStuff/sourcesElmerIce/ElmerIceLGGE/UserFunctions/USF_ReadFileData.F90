!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Read Data from a file      
!  
!  work only in 2D at that time
! OG 2009/04/01
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION ReadFileData ( Model, nodenumber, x ) RESULT(value)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: value, x 

   REAL(KIND=dp), ALLOCATABLE :: data_value(:), data_value2(:)
   REAL(KIND=dp), ALLOCATABLE :: Coordinate(:)
   INTEGER, ALLOCATABLE :: ind(:)
   INTEGER :: Ns, i 
   LOGICAL :: FirstTime = .TRUE.  
   CHARACTER(LEN=MAX_NAME_LEN) :: FileDataName

   SAVE FirstTime, Ns
   SAVE data_value, data_value2, Coordinate

       INTERFACE 
         SUBROUTINE indexx(arr,ind) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp) :: arr(:) 
         INTEGER :: ind(:) 
         END SUBROUTINE indexx

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

   IF (FirstTime) THEN
      FirstTime = .FALSE.

!------------------------------------
! Get the name of the input data file         
!------------------------------------
!TODO
      FileDataName = 'aFactor_variegated.dat'

      
!------------------------------------
! Which coordinate are given x, y, z ?          
!------------------------------------
!TODO


!------------------------
! Read and Store the data                     
!------------------------
      
      OPEN(1,file=TRIM(FileDataName))
      READ(1,*)Ns

      ALLOCATE (Coordinate(Ns), data_value(Ns), data_value2(Ns))

      READ(1,*)(Coordinate(i), data_value(i), i=1,Ns)
      CLOSE(1) 

!------------------------------------
! Be sure the coordinate are sorted         
!------------------------------------
      ALLOCATE (ind(Ns))
      ind = 0
      CALL indexx(Coordinate,ind)
      data_value = data_value(ind) 
      DEALLOCATE(ind)
     

!-----------------------------------
! Construct the spline interpolation 
!-----------------------------------
      CALL spline(coordinate,data_value,0.0_dp,0.0_dp,data_value2)
      
   ENDIF ! End First Time


!-----------------------
! Interpolate for that x
!-----------------------
    x = Model % Nodes % x ( nodenumber)
    value = splint(coordinate,data_value,data_value2,x)
    write(*,*)x,value


END FUNCTION ReadFileData       




! -------------------------------------------------------------
! Numerical recipes subroutines  and functions          -------
! -------------------------------------------------------------
! --------------------------------------------------------
       SUBROUTINE indexx(arr,ind) 
       USE types
       USE DefUtils
       IMPLICIT NONE 
       REAL(KIND=dp) :: arr(:) 
       INTEGER :: ind(:) 
       INTEGER, PARAMETER :: NN=15, NSTACK=50
       REAL(KIND=dp) :: a 
       INTEGER :: n,k,i,j,indext,jstack,l,r 
       INTEGER, DIMENSION(NSTACK) :: istack 

        n = size(arr)
        Do i = 1, n
        ind( i ) = i 
        END DO
jstack=0 
l=1 
r=n 
do 
if (r-l < NN) then 
do j=l+1,r 
indext=ind(j) 
a=arr(indext) 
do i=j-1,l,-1 
if (arr(ind(i)) <= a) exit 
ind(i+1)=ind(i) 
end do 
ind(i+1)=indext 
end do 
if (jstack == 0) RETURN 
r=istack(jstack) 
l=istack(jstack-1) 
jstack=jstack-2 
else 
k=(l+r)/2 
call swap(ind(k),ind(l+1)) 
call icomp_xchg(ind(l),ind(r)) 
call icomp_xchg(ind(l+1),ind(r)) 
call icomp_xchg(ind(l),ind(l+1)) 
i=l+1 
j=r 
indext=ind(l+1) 
a=arr(indext) 
do 
do 
i=i+1 
if (arr(ind(i)) >= a) exit 
end do 
do 
j=j-1 
if (arr(ind(j)) <= a) exit 
end do 
if (j < i) exit 
call swap(ind(i),ind(j)) 
end do 
ind(l+1)=ind(j) 
ind(j)=indext 
jstack=jstack+2 
if (jstack > NSTACK) Write(*,*)' indexx: NSTACK too small'
if (r-i+1 >= j-l) then 
istack(jstack)=r
istack(jstack-1)=i 
r=j-1 
else 
istack(jstack)=j-1 
istack(jstack-1)=l 
l=i 
end if 
end if 
end do 
CONTAINS 
SUBROUTINE icomp_xchg(i,j) 
INTEGER :: i,j 
INTEGER :: swp 
if (arr(j) < arr(i)) then 
swp=i 
i=j 
j=swp 
end if 
END SUBROUTINE icomp_xchg 
SUBROUTINE swap(a,b) 
INTEGER :: a, b
INTEGER :: dum
dum = a
a = b
b = dum
END SUBROUTINE swap
END SUBROUTINE indexx
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

! ----------------------------------------------------------------------------
