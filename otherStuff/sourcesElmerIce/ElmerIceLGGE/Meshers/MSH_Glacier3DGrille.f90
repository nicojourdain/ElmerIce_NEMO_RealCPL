!
! Create a 3D mesh given by its bed and surface topography 
!
!
!
   PROGRAM  MSH_Glacier3DGrille 
       USE types
       USE DefUtils
!
IMPLICIT NONE
!
REAL(KIND=dp)  ::  x, y, z 
REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:), xs(:), & 
                       ys(:), zs(:)
REAL(KIND=dp), ALLOCATABLE  :: xnode(:), ynode(:), znode(:)                       
CHARACTER :: NameMsh*20, Rien*1, NameSurf*20, NameBed*20
REAL(KIND=dp)  :: dsx, dsy, dbx, dby, x1, x2, y1, y2, zi(2,2)
REAL(KIND=dp)  :: xs0, ys0, xb0, yb0, zbed, zsurf, znew, Rmin, R
REAL(KIND=dp)  :: lsx, lsy, lbx, lby 
INTEGER :: NtN, i, j, Ns, Nsx, Nsy, Nb, Nbx, Nby, n, Npt, ix, iy, imin, is, ib
!  
!
!  Read data input from mesh_input.dat
!
!
      OPEN(10,file="mesh_input.dat")
      READ(10,*)Rien
      READ(10,*)NameMsh
      READ(10,*)Rien
      READ(10,*)NameSurf
      READ(10,*)Rien
      READ(10,*)Nsx, Nsy
      READ(10,*)Rien
      READ(10,*)xs0, ys0
      READ(10,*)Rien
      READ(10,*)lsx, lsy
      READ(10,*)Rien
      READ(10,*)NameBed
      READ(10,*)Rien
      READ(10,*)Nbx, Nby
      READ(10,*)Rien
      READ(10,*)xb0, yb0
      READ(10,*)Rien
      READ(10,*)lbx, lby
      CLOSE(10)
      Ns = Nsx*Nsy
      Nb = Nbx*Nby
      
      dsx = lsx / (Nsx-1.0)
      dsy = lsy / (Nsy-1.0)
      dbx = lbx / (Nbx-1.0)
      dby = lby / (Nby-1.0)
      
      write(*,*)Ns,Nb
      write(*,*)dsx,dsy,dbx,dby

      OPEN(10,file=TRIM(NameMsh)//"/mesh.header")
        READ(10,1000)NtN
      CLOSE(10)
      ALLOCATE(xnode(NtN), ynode(NtN), znode(NtN))

      OPEN(10,file=TRIM(NameSurf))
        ALLOCATE(xs(Ns), ys(Ns), zs(Ns))
        READ(10,*)(xs(i), ys(i), zs(i), i=1,Ns)
      CLOSE(10)

      OPEN(10,file=TRIM(NameBed))
        ALLOCATE(xb(Nb), yb(Nb), zb(Nb))
        READ(10,*)(xb(i), yb(i), zb(i), i=1,Nb)
      CLOSE(10)
      
      
      OPEN(12,file=TRIM(NameMsh)//"/mesh.nodes")
      READ(12,*)(N, j, xnode(i), ynode(i), znode(i), i=1,NtN)

      REWIND(12)


      DO n=1, NtN

        x = xnode(n)
        y = ynode(n)
        z = znode(n)

! Find zbed for that point from the Bedrock MNT 

        ix = INT((x-xb0)/dbx)+1
        iy = INT((y-yb0)/dby)+1
        ib = Nbx * (iy - 1) + ix
        
        x1 = xb(ib)
        x2 = xb(ib+1)
        y1 = yb(ib)
        y2 = yb(ib + Nbx)
        
        IF ((x<x1).OR.(x>x2).OR.(y<y1).OR.(y>y2)) WRITE(*,*)'Pb bedrock',i
        
        zi(1,1) = zb(ib)
        zi(2,1) = zb(ib+1)
        zi(2,2) = zb(ib + Nbx + 1)
        zi(1,2) = zb(ib + Nbx)
        
        
        DO i=1, 2
          DO J=1, 2
            IF (zi(i,j) < 0.0) write(*,*)'Bedrock zi(1,1)',n,i,j, zi(i,j)
          END DO
        END DO
        
        
        IF ((zi(1,1)<-9990.0).OR.(zi(1,2)<-9990.0).OR.(zi(2,1)<-9990.0).OR.(zi(2,2)<-9990.0)) THEN
           IF ((zi(1,1)<-9990.0).AND.(zi(1,2)<-9990.0).AND.(zi(2,1)<-9990.0).AND.(zi(2,2)<-9990.0)) THEN
           ! Find the nearest point avalable
             Rmin = 9999.0
             DO i=1, Nb
               IF (zb(i)>0.0) THEN
                 R = SQRT((x-xb(i))**2.0+(y-yb(i))**2.0)
                 IF (R<Rmin) THEN
                   Rmin = R
                   imin = i
                 END IF
               END IF
             END DO
            zbed = zb(imin)
            write(*,*)'all null', zbed,imin
            
           ELSE
            ! Mean value over the avalable data
             zbed = 0.0
             Npt = 0
             DO i=1, 2
               DO J=1, 2
                  IF (zi(i,j) > 0.0) THEN 
                     zbed = zbed + zi(i,j)
                     Npt = Npt + 1
                  END IF   
               END DO
             END DO
             zbed = zbed / Npt
             write(*,*)'not all null',zbed,Npt
           END IF
        ELSE
          zbed = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dbx*dby)
        END IF

! Find zsurf for that point from the Surface MNT
        
        ix = INT((x-xs0)/dsx)+1
        iy = INT((y-ys0)/dsy)+1
        is = Nsx * (iy - 1) + ix
        
        x1 = xs(is)
        x2 = xs(is+1)
        y1 = ys(is)
        y2 = ys(is + Nsx)
        
        IF ((x<x1).OR.(x>x2).OR.(y<y1).OR.(y>y2)) WRITE(*,*)'Pb surface',is,x,x1,x2,y,y1,y2
        
        zi(1,1) = zs(is)
        zi(2,1) = zs(is+1)
        zi(2,2) = zs(is + Nsx + 1)
        zi(1,2) = zs(is + Nsx)
        
        DO i=1, 2
          DO J=1, 2
            IF ((zi(i,j) < 0.0).OR.(zi(i,j)>4000.0)) write(*,*)'surface zi(1,1)',n,i,j, zi(i,j)
          END DO
        END DO
        
        IF ((zi(1,1)<-9990.0).OR.(zi(1,2)<-9990.0).OR.(zi(2,1)<-9990.0).OR.(zi(2,2)<-9990.0)) THEN
           IF ((zi(1,1)<-9990.0).AND.(zi(1,2)<-9990.0).AND.(zi(2,1)<-9990.0).AND.(zi(2,2)<-9990.0)) THEN
           ! Find the nearest point avalable
             Rmin = 9999.0
             DO i=1, Ns
               IF (zs(i)>0.0) THEN
                 R = SQRT((x-xs(i))**2.0+(y-ys(i))**2.0)
                 IF (R<Rmin) THEN
                   Rmin = R
                   imin = i
                 END IF
               END IF
             END DO
            zsurf = zs(imin)
            write(*,*)'all null', zsurf,imin
            
           ELSE
            ! Mean value over the avalable data
             zsurf = 0.0
             Npt = 0
             DO i=1, 2
               DO J=1, 2
                  IF (zi(i,j) > 0.0) THEN 
                     zsurf = zsurf + zi(i,j)
                     Npt = Npt + 1
                  END IF   
               END DO
             END DO
             zsurf = zsurf / Npt
             write(*,*)'not all null',zsurf,Npt
           END IF
        ELSE
           zsurf = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dsx*dsy)
        END IF      
        
        IF ((zsurf + 1.0) .LE. zbed) THEN
           !write(*,*)n,zsurf,zbed,z
           zsurf = zbed + 1.0
        END IF
        
        znew = zbed + z * (zsurf - zbed) 
        IF (znew > 4000.0) write(*,*)zbed,zsurf,z,n
        
        WRITE(12,1200)N,j,x,y,znew
      END DO
      WRITE(*,*)'END WITH NO TROUBLE ...'
!
1000 FORMAT(4x,I6)
1200 FORMAT(i5,2x,i5,3(2x,e22.15)) 

End Program MSH_Glacier3DGrille

! --------------------------------------------------------

