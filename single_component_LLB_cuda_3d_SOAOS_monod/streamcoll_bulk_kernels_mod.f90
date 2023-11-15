#include "defines.h"
 module streamcoll_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine streamcoll_shared()
	
	implicit none  
	
    integer :: i,j,k,gi,gj,gk,myblock,iidblock
    integer :: gii,gjj,gkk,ii,jj,kk
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
    integer :: li,lj,lk
	real(kind=db), shared :: f00(1:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f01(0:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f02(1:TILE_DIMx_d+1,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f03(1:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f04(1:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f05(1:TILE_DIMx_d,1:TILE_DIMy_d,0:TILE_DIMz_d)
	real(kind=db), shared :: f06(1:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
	real(kind=db), shared :: f07(0:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f08(1:TILE_DIMx_d+1,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f09(0:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f10(1:TILE_DIMx_d+1,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f11(1:TILE_DIMx_d,0:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f12(1:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f13(1:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f14(1:TILE_DIMx_d,1:TILE_DIMy_d+1,0:TILE_DIMz_d)
    real(kind=db), shared :: f15(0:TILE_DIMx_d,1:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f16(1:TILE_DIMx_d+1,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f17(1:TILE_DIMx_d+1,1:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f18(0:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
    
		  
		  
	!if(isfluid(i,j,k).ne.1)return
    	
    
    gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	li=threadIdx%x
	lj=threadIdx%y
	lk=threadIdx%z
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
	
	myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
!	 myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
!    zblock=(idblock-1)/nxyblock +1
!    yblock=((idblock-1)-(zblock-1)*nxyblock)/nxblock +1
!    xblock=(idblock-1)-(zblock-1)*nxyblock-(yblock-1)*nxblock +1
	
	!myblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
	iidblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	!iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	if(myblock.ne.iidblock)write(*,*)'cazzo amaro streamcoll_shared',gi,gj,gk,myblock,iidblock
    
        
    uu=halfonecssq*(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10)) + &
     hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)) + &
     hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))
    
    !0
	f00(li,lj,lk)=p0*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))))
	
    
	!1 -1  0  0
	udotc=hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))&
	 +hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_7_8*hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_7_8*hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_9_10*hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_9_10*hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_15_16*hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_15_16*hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
 	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_17_18*hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_17_18*hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_11_12*hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_11_12*hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))-&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_13_14*hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_13_14*hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    ! Halo Faces
    if(li==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) +&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) +&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !1 -1  0  0
	  udotc=hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f01(li-1,lj,lk)=p1*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10)))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f07(li-1,lj,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f09(li-1,lj,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f15(li-1,lj,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f18(li-1,lj,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) +&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) +&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !2 +1  0  0
	  udotc=hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f02(li+1,lj,lk)=p1*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10)))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f08(li+1,lj,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f10(li+1,lj,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f16(li+1,lj,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f17(li+1,lj,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) +&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) +& 
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !3 0 -1  0
	  udotc=hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f03(li,lj-1,lk)=p1*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10)))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f07(li,lj-1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f10(li,lj-1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f11(li,lj-1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	  hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f13(li,lj-1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !4  0 +1  0
	  udotc=hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f04(li,lj+1,lk)=p1*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10)))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f08(li,lj+1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f09(li,lj+1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f12(li,lj+1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f14(li,lj+1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      gii=gi
      gjj=gj
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !5  0  0 -1
	  udotc=hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f05(li,lj,lk-1)=p1*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10)))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f15(li,lj,lk-1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f17(li,lj,lk-1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f11(li,lj,lk-1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f14(li,lj,lk-1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !6  0  0  +1
	  udotc=hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f06(li,lj,lk+1)=p1*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10)))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f16(li,lj,lk+1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f18(li,lj,lk+1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f12(li,lj,lk+1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f13(li,lj,lk+1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      gii=gi-1
      gjj=gj-1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !7 -1 -1  0
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f07(li-1,lj-1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      gii=gi-1
      gjj=gj+1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !9  -1 +1 0
      udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f09(li-1,lj+1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !15  -1  0 -1
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f15(li-1,lj,lk-1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi-1
      gjj=gj
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !18   -1   0  +1
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f18(li-1,lj,lk+1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      gii=gi+1
      gjj=gj-1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !10   +1 -1  0
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f10(li+1,lj-1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      gii=gi+1
      gjj=gj+1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !8 +1 +1  0
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f08(li+1,lj+1,lk)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfields(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      gii=gi+1
      gjj=gj
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !17  +1  0 -1
	  udotc=(-hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f17(li+1,lj,lk-1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !16  +1  0 +1
	  udotc=(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f16(li+1,lj,lk+1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfields(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !11  0  -1  -1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f11(li,lj-1,lk-1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj-1
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !13  0  -1   +1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f13(li,lj-1,lk+1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      gii=gi
      gjj=gj+1
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !14  0  +1  -1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f14(li,lj+1,lk-1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !12  0  +1  +1
	  udotc=(hfields(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfields(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f12(li,lj+1,lk+1)=p2*(hfields(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfields(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfields(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfields(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
#ifdef IFBC    
    if(abs(isfluid(gi,gj,gk)).ne.1)return
#endif
    udotc=f00(li,lj,lk)+f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))=udotc
	!if(gk==1 .and. gi==3 .and. gj==3)write(*,*)gi,gj,gk,mystep,udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))=udotc
     
#ifdef PRESSCORR
	
	uu=halfonecssq*(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*&
	 hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10)) + &
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)) + &
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))
    
	!1 -1  0  0
	udotc=hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))-&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	    
    return
#endif	

  end subroutine streamcoll_shared
  
  attributes(global) subroutine streamcoll_shared_flop()
	
	implicit none  
	
	
	integer :: i,j,k,gi,gj,gk,myblock,iidblock
    integer :: gii,gjj,gkk,ii,jj,kk
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
    integer :: li,lj,lk
	real(kind=db), shared :: f00(1:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f01(0:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f02(1:TILE_DIMx_d+1,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f03(1:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f04(1:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f05(1:TILE_DIMx_d,1:TILE_DIMy_d,0:TILE_DIMz_d)
	real(kind=db), shared :: f06(1:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
	real(kind=db), shared :: f07(0:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f08(1:TILE_DIMx_d+1,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f09(0:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f10(1:TILE_DIMx_d+1,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f11(1:TILE_DIMx_d,0:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f12(1:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f13(1:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f14(1:TILE_DIMx_d,1:TILE_DIMy_d+1,0:TILE_DIMz_d)
    real(kind=db), shared :: f15(0:TILE_DIMx_d,1:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f16(1:TILE_DIMx_d+1,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f17(1:TILE_DIMx_d+1,1:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f18(0:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
    
		  
		  
	!if(isfluid(i,j,k).ne.1)return
    
    
    gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	li=threadIdx%x
	lj=threadIdx%y
	lk=threadIdx%z
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
	
	
	
	 myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
!    zblock=(idblock-1)/nxyblock +1
!    yblock=((idblock-1)-(zblock-1)*nxyblock)/nxblock +1
!    xblock=(idblock-1)-(zblock-1)*nxyblock-(yblock-1)*nxblock +1
	
	!myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	!iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	!if(myblock.ne.iidblock)write(*,*)'cazzo2',gi,gj,gk,myblock,iidblock
    
    
    
    uu=halfonecssq*(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*&
     hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10)) + &
     hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*&
     hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)) + &
     hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*&
     hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))
    
    !0
	f00(li,lj,lk)=p0*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))))
	
    
	!1 -1  0  0
	udotc=hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_7_8*hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_7_8*hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_9_10*hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_9_10*hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_15_16*hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_15_16*hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_17_18*hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_17_18*hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_11_12*hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_11_12*hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))-&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_13_14*hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_13_14*hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    ! Halo Faces
    if(li==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !1 -1  0  0
	  udotc=hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f01(li-1,lj,lk)=p1*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10)))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f07(li-1,lj,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f09(li-1,lj,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f15(li-1,lj,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f18(li-1,lj,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !2 +1  0  0
	  udotc=hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f02(li+1,lj,lk)=p1*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10)))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f08(li+1,lj,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f10(li+1,lj,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+two*qxy_9_10*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f16(li+1,lj,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f17(li+1,lj,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !3 0 -1  0
	  udotc=hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f03(li,lj-1,lk)=p1*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10)))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f07(li,lj-1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f10(li,lj-1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f11(li,lj-1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f13(li,lj-1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !4  0 +1  0
	  udotc=hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f04(li,lj+1,lk)=p1*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10)))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f08(li,lj+1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f09(li,lj+1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f12(li,lj+1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f14(li,lj+1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      gii=gi
      gjj=gj
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !5  0  0 -1
	  udotc=hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f05(li,lj,lk-1)=p1*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10)))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f15(li,lj,lk-1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f17(li,lj,lk-1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f11(li,lj,lk-1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f14(li,lj,lk-1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !6  0  0  +1
	  udotc=hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*onecssq
	  f06(li,lj,lk+1)=p1*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*(hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10)))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f16(li,lj,lk+1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f18(li,lj,lk+1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f12(li,lj,lk+1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f13(li,lj,lk+1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      gii=gi-1
      gjj=gj-1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
      hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
      hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
      hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
      hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
      hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !7 -1 -1  0
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f07(li-1,lj-1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      gii=gi-1
      gjj=gj+1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !9  -1 +1 0
      udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f09(li-1,lj+1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !15  -1  0 -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f15(li-1,lj,lk-1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi-1
      gjj=gj
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !18   -1   0  +1
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f18(li-1,lj,lk+1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      gii=gi+1
      gjj=gj-1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !10   +1 -1  0
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f10(li+1,lj-1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_9_10*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      gii=gi+1
      gjj=gj+1
      gkk=gk
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !8 +1 +1  0
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f08(li+1,lj+1,lk)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxy_7_8*hfieldsh(idx5d(ii,jj,kk,8,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      gii=gi+1
      gjj=gj
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !17  +1  0 -1
	  udotc=(-hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f17(li+1,lj,lk-1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_17_18*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !16  +1  0 +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f16(li+1,lj,lk+1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qxz_15_16*hfieldsh(idx5d(ii,jj,kk,9,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !11  0  -1  -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f11(li,lj-1,lk-1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj-1
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !13  0  -1   +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f13(li,lj-1,lk+1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      gii=gi
      gjj=gj+1
      gkk=gk-1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !14  0  +1  -1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))-&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f14(li,lj+1,lk-1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_13_14*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk+1
      xblock=(gii+2*TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,2,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10)) + &
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10))*&
       hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))
      
      !12  0  +1  +1
	  udotc=(hfieldsh(idx5d(ii,jj,kk,3,iidblock,nx_d,ny_d,nz_d,10))+&
	   hfieldsh(idx5d(ii,jj,kk,4,iidblock,nx_d,ny_d,nz_d,10)))*onecssq
	  f12(li,lj+1,lk+1)=p2*(hfieldsh(idx5d(ii,jj,kk,1,iidblock,nx_d,ny_d,nz_d,10))+&
	   (-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(ii,jj,kk,6,iidblock,nx_d,ny_d,nz_d,10))+&
	   qzz*hfieldsh(idx5d(ii,jj,kk,7,iidblock,nx_d,ny_d,nz_d,10))-&
	   cssq*hfieldsh(idx5d(ii,jj,kk,5,iidblock,nx_d,ny_d,nz_d,10))+&
	   two*qyz_11_12*hfieldsh(idx5d(ii,jj,kk,10,iidblock,nx_d,ny_d,nz_d,10))) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
#ifdef IFBC    
    if(abs(isfluid(gi,gj,gk)).ne.1)return
#endif
    udotc=f00(li,lj,lk)+f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))=udotc
     
#ifdef PRESSCORR

	
	uu=halfonecssq*(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*&
	 hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10)) + &
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)) + &
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))
    
	!1 -1  0  0
	udotc=hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))-&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))-udotc
	
#endif
    
    
    return
 
  end subroutine streamcoll_shared_flop
  
    attributes(global) subroutine streamcoll_shared_halo()
	
	implicit none  
	  
    integer :: i,j,k,gi,gj,gk,myblock,iidblock
    integer :: gii,gjj,gkk,ii,jj,kk
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
    integer :: li,lj,lk
	real(kind=db), shared :: f00(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f01(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f02(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f03(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f04(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f05(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f06(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f07(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f08(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f09(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f10(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f11(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f12(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f13(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f14(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f15(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f16(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f17(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f18(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    
		  
		  
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x-1
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y-1
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z-1
	
	li = threadIdx%x-1
    lj = threadIdx%y-1
    lk = threadIdx%z-1

    xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    
    i=gi-xblock*TILE_DIMx_d+2*TILE_DIMx_d
	j=gj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
	
	myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	!if(gi==6 .and. isfluid(gi,gj,gk).ne.3)write(*,*)'test',gi,gj,gk,myblock,isfluid(gi,gj,gk)
    
    
    
    uu=halfonecssq*(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10)) + &
     hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)) + &
     hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))
    
    !0
	f00(li,lj,lk)=p0*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))))
	
    
	!1 -1  0  0
	udotc=hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_7_8*hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_7_8*hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_9_10*hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_9_10*hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_15_16*hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_15_16*hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_17_18*hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_17_18*hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_11_12*hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_11_12*hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))-&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_13_14*hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_13_14*hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    call syncthreads
    
    if(li<1 .or. lj<1 .or. lk<1)return
    if(li>TILE_DIMx_d .or. lj>TILE_DIMy_d .or. lk>TILE_DIMz_d)return
#ifdef IFBC    
    if(abs(isfluid(gi,gj,gk)).ne.1)return
#endif
    udotc=f00(li,lj,lk)+f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))=udotc
     
#ifdef PRESSCORR

	!call syncthreads
	
	uu=halfonecssq*(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*&
	 hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10)) &
     + hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*&
     hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)) &
     + hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*&
     hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))
    
	!1 -1  0  0
	udotc=hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))-&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	    
    return
#endif	
  end subroutine streamcoll_shared_halo
  
  attributes(global) subroutine streamcoll_shared_halo_flop()
	
	implicit none  
	
	
	integer :: i,j,k,gi,gj,gk,myblock,iidblock
    integer :: gii,gjj,gkk,ii,jj,kk
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
    integer :: li,lj,lk
	real(kind=db), shared :: f00(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f01(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f02(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f03(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f04(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f05(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f06(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f07(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f08(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f09(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f10(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
	real(kind=db), shared :: f11(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f12(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f13(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f14(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f15(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f16(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f17(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    real(kind=db), shared :: f18(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
    
		  
		  
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x -1
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y -1
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z -1
	
	li = threadIdx%x-1
    lj = threadIdx%y-1
    lk = threadIdx%z-1

    xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    
    i=gi-xblock*TILE_DIMx_d+2*TILE_DIMx_d
	j=gj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
	
	myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	!if(gi==6 .and. isfluid(gi,gj,gk).ne.3)write(*,*)'test',gi,gj,gk,myblock,isfluid(gi,gj,gk)
    
    
    
    uu=halfonecssq*(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*&
     hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10)) &
     + hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*&
     hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)) &
     + hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*&
     hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))
    
    !0
	f00(li,lj,lk)=p0*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))))
	
    
	!1 -1  0  0
	udotc=hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-&
	 cssq*(hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*(hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_7_8*hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_7_8*hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_9_10*hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxy_9_10*hfieldsh(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_15_16*hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_15_16*hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfieldsh(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_17_18*hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qxz_17_18*hfieldsh(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_11_12*hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_11_12*hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfieldsh(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))-&
	 hfieldsh(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_13_14*hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfieldsh(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))+&
	 qzz*hfieldsh(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10)) &
	 -cssq*hfieldsh(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))+&
	 two*qyz_13_14*hfieldsh(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    
    call syncthreads
    
    if(li<1 .or. lj<1 .or. lk<1)return
    if(li>TILE_DIMx_d .or. lj>TILE_DIMy_d .or. lk>TILE_DIMz_d)return
#ifdef IFBC    
    if(abs(isfluid(gi,gj,gk)).ne.1)return
#endif
    udotc=f00(li,lj,lk)+f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))=udotc
     
#ifdef PRESSCORR

	
	uu=halfonecssq*(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*&
	 hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10)) &
	 + hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)) &
	 + hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))
    
	!1 -1  0  0
	udotc=hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10))*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(idx5d(i,j,k,2,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(idx5d(i,j,k,3,myblock,nx_d,ny_d,nz_d,10))-&
	 hfields(idx5d(i,j,k,4,myblock,nx_d,ny_d,nz_d,10)))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(idx5d(i,j,k,1,myblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,5,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,6,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,7,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,8,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,9,myblock,nx_d,ny_d,nz_d,10))-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,10,myblock,nx_d,ny_d,nz_d,10))-udotc
	
#endif
    
    return
 
  end subroutine streamcoll_shared_halo_flop
  
 end module streamcoll_kernels
