! Model program to compress given Gausssians
! Kernal Density Compression Algorithm
!
! Authour: Anji Babu Kapakayala
!          Dept. of Chemistry
!          IIT Kanpur, India.
!          anjibabu480@gmail.com
!
! OUTPUTS: Original_Gaus.dat
!          Compressed_Gaus.dat 
!
PROGRAM kde
IMPLICIT NONE
REAL*8::x,dx0,x_min,x_max,x_width,dum,dum1,error
INTEGER::i,j,nbin=10001,ng=15000,CNG
REAL*8,ALLOCATABLE::grid(:),G(:),ht(:),sigma(:),x0(:)
REAL*8,ALLOCATABLE::C_hill(:),C_width(:),C_height(:)
!Initializing grid values
x_min=-100.d0;x_max=100.d0
x_width=(x_max-x_min)/(nbin-1)

open(10,file="Original_Gaus.dat",status='unknown')
open(11,file="Compressed_Gaus.dat",status='unknown')

ALLOCATE(grid(nbin),G(nbin))
ALLOCATE(ht(ng),sigma(ng),x0(ng))
ALLOCATE(C_hill(ng),C_width(ng),C_height(ng))

! Gridding data
do i=1,nbin
grid(i)=DFLOAT(i-1)*x_width+x_min
enddo

!Initializing Gaussian Parameters
x0(1)=-35.0;sigma(1)=0.050
ht(1)=0.20;dx0=0.0050

! Generate Gaussians centers 
DO j=1,ng
ht(j)=0.50
sigma(j)=0.05
!write(*,*) x0(j),sigma(j),ht(j)
x0(j+1)=x0(j)+dx0
ENDDO

!Summing the original Gaussians
DO i=1,nbin
   dum=0.0
   DO j=1,ng
     dum=dum+ht(j)*EXP(-(((grid(i)-x0(j))**2)/2*(sigma(j)*sigma(j))))
   ENDDO
write(10,*) grid(i),dum
ENDDO

! Compressing the given Gaussians using Kernel Density Compression Algorithm
CALL Compress(ng,x0,ht,sigma,CNG,C_hill,C_width,C_height)

!Summing the Compressed Gaussians
DO i=1,nbin
   dum1=0.0
   DO j=1,CNG
    if(ng.eq.CNG)then
      dum1=dum1+C_height(j)*EXP(-(((grid(i)-C_hill(j))**2)/2*(C_width(j)*C_width(j))))
    else
      dum1=dum1+C_height(j)*EXP(-(((grid(i)-C_hill(j))**2)/2*(C_width(j))))
    endif
   ENDDO
write(11,*) grid(i),dum1
ENDDO


WRITE(*,*) "No of Original Gassians    : ",ng,"No. of Compressed Gaussians: ",CNG
WRITE(*,*) "Sum of original Gaussians  : ",dum
WRITE(*,*) "Sum of compressed Gaussians: ",dum1

! Calculating relative error
error=(abs(dum1-dum)/dum)*100
WRITE(*,*) "Relative error(%)          : ",error
close(10)
close(11)
ENDPROGRAM kde
!------------------------------------------------------------------------------------------!
!Subroutine to compress Given gasuuians using Kernel Density Compression Algorithm
SUBROUTINE Compress(NG,hill,height,width,CNG,C_hill,C_width,C_height)
IMPLICIT NONE
 INTEGER, INTENT(IN)::NG
 REAL*8,  INTENT(IN)::hill(NG),width(NG),height(NG)
 REAL*8, INTENT(OUT)::C_hill(NG),C_width(NG),C_height(NG)
 INTEGER,INTENT(OUT)::CNG
!
 integer :: i,j,n,min_j
 real*8  :: distance,cutoff=1.0,dum,new_height,new_width,new_hill,min_d
!
! Initialization of Compressed Gaussians
 C_hill(1) =hill(1)
 C_width(1)=width(1)
 C_height(1)=height(1)
 n=1
 
 DO i=2,NG
    min_d=sqrt(ABS(((hill(i))-C_hill(1))/C_width(1))**2)
    min_j=1
      do j=1,n
         distance=sqrt(ABS(((hill(i))-C_hill(j))/C_width(j))**2)
         if(distance.lt.min_d) then
            min_d=distance
            min_j=j
         endif
        ! print*,"Distance: ",distance
     enddo
     !print*,"Minimum Distance: ",min_d,"LOCATION: ",min_j
     !If minimun distance is less than given cutoff, it merges the both gaussinas
     ! Add as new compressed kernal otherwise.   
     if(min_d.lt.cutoff)then
           new_height=height(i)+C_height(min_j)
           new_hill  =(1.0/new_height)*((height(i)*hill(i))+(C_height(min_j)*C_hill(min_j)))
           dum=(height(i)*(width(i)**2+hill(i)**2)+C_height(min_j)*(C_width(min_j)**2+C_hill(min_j)**2))
           new_width=(1.0/new_height)*dum-new_hill**2
           ! Store new variables to min_jth kernel
           C_height(min_j)=new_height           
           C_width(min_j)=new_width
           C_hill(min_j)=new_hill
        else
           n=n+1
           C_height(n)=height(i)  
           C_width(n)=width(i)
           C_hill(n)=hill(i)
        endif   

 ENDDO
! Store n as No. of Compressed Gaussians
 CNG=n

END SUBROUTINE Compress
!------------------------------------------------------------------------------------------!
!                           Anji Babu Kapakayala, Jan, 2021			           !
!------------------------------------------------------------------------------------------!
