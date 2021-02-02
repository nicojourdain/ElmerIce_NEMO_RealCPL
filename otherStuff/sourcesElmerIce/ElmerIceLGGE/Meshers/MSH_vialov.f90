!
! Create Vialov  mesh 
! by deforming a square 1x1 mesh
!
! Input :
!    Name of the grd mesh (a square 1 x 1)
!    (h/h0)^{2+2/n} + (x/L)^{1+1/n} = 1
!    (h/h0)^{nh} + (x/L)^{nx} = 1
!    Read h0, n, Eps = h0/L and the mesh length Lm
!
! Output :
!    mesh/mesh.nodes
!
!
   PROGRAM  MSH_vialov  
!
IMPLICIT NONE
!
REAL(KIND=8)  ::   x, y, z, xnew, ynew, zb, zs, L, acc
REAL(KIND=8)  ::   h0, Eps, ng, Lm, nx, nh, rhog, Bg, K                   
REAL(KIND=8), ALLOCATABLE  :: xnode(:), ynode(:)                       
CHARACTER :: NameMsh*20
INTEGER :: NtN, i, j, N
!
!
      rhog = 0.00893   ! MPa.m-1 (g=9.81 and rho = 910kg/m3)
      Bg = 200.0       ! MPa-3 a-1
!
      WRITE(*,*)'Name of the Elmer mesh directory ?'
      READ(*,*)NameMsh
      WRITE(*,*)'Heigth at the dome H0 [m] :'
      READ(*,*)h0
      WRITE(*,*)'Accumulation [m/a] :'
      READ(*,*)acc
      WRITE(*,*)'Glen law exponent n :'
      READ(*,*)ng 
      nx = (ng+1.0)/ng
      nh = 2.0*nx
      K = 2.0 / rhog * ((ng+2.0) * acc / Bg )**(1.0/ng)
      L = (H0**nh / K)**(1.0/nx)
      Eps = H0 / L
      WRITE(*,*)'** Ice Sheet length L = ', L
      WRITE(*,*)'** Ice sheet aspect ratio Eps = ', Eps 
      WRITE(*,*)'Lenght of the mesh Lm [m] < L :'
      READ(*,*)Lm 

!     IF (Lm >= L ) THEN
!             WRITE(*,*)'Mesh length Lm must be smaller than L = ', L
!             STOP
!     END IF

      OPEN(10,file=TRIM(NameMsh)//"/mesh.header")
        READ(10,1000)NtN
      CLOSE(10)

      ALLOCATE(xnode(NtN), ynode(NtN))

      
      OPEN(12,file=TRIM(NameMsh)//"/mesh.nodes")
      READ(12,*)(N, j, xnode(i), ynode(i), z, i=1,NtN)
      REWIND(12)

      DO N=1, NtN

        x = xnode(N)
        y = ynode(N)
        xnew = x * Lm 
        IF (xnew < L) THEN 
           zs = H0 * (1.0 - (xnew/L)**nx)**(1.0/nh) 
        ELSE
           zs = zb + 1.0
        END IF
        zb = 0.0 
        ynew = zb + y * (zs - zb)
        
        WRITE(12,1200)N,j,xnew,ynew,z
      END DO
!
      DEALLOCATE(xnode, ynode)

1000 FORMAT(I6)
1200 FORMAT(i6,2x,i5,3(2x,e22.15)) 

End Program MSH_vialov 

