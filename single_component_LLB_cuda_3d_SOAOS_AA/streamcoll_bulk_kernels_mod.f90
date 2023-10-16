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
    
        
    uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(hfields(i,j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock)-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock)-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
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
      
      uu=halfonecssq*(xshell1(jj,kk,2,iidblock)*xshell1(jj,kk,2,iidblock) + xshell1(jj,kk,3,iidblock)*xshell1(jj,kk,3,iidblock) + xshell1(jj,kk,4,iidblock)*xshell1(jj,kk,4,iidblock))
      
      !1 -1  0  0
	  udotc=xshell1(jj,kk,2,iidblock)*onecssq
	  f01(li-1,lj,lk)=p1*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*xshell1(jj,kk,5,iidblock)-cssq*(xshell1(jj,kk,6,iidblock)+xshell1(jj,kk,7,iidblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(xshell1(jj,kk,2,iidblock)+xshell1(jj,kk,3,iidblock))*onecssq
	  f07(li-1,lj,lk)=p2*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(jj,kk,5,iidblock)+qyy*xshell1(jj,kk,6,iidblock)-cssq*xshell1(jj,kk,7,iidblock)+two*qxy_7_8*xshell1(jj,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-xshell1(jj,kk,2,iidblock)+xshell1(jj,kk,3,iidblock))*onecssq
	  f09(li-1,lj,lk)=p2*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(jj,kk,5,iidblock)+qyy*xshell1(jj,kk,6,iidblock)-cssq*xshell1(jj,kk,7,iidblock)+two*qxy_9_10*xshell1(jj,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(xshell1(jj,kk,2,iidblock)+xshell1(jj,kk,4,iidblock))*onecssq
	  f15(li-1,lj,lk)=p2*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(jj,kk,5,iidblock)+qzz*xshell1(jj,kk,7,iidblock)-cssq*xshell1(jj,kk,6,iidblock)+two*qxz_15_16*xshell1(jj,kk,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-xshell1(jj,kk,2,iidblock)+xshell1(jj,kk,4,iidblock))*onecssq
	  f18(li-1,lj,lk)=p2*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(jj,kk,5,iidblock)+qzz*xshell1(jj,kk,7,iidblock)-cssq*xshell1(jj,kk,6,iidblock)+two*qxz_17_18*xshell1(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(xshell0(jj,kk,2,iidblock)*xshell0(jj,kk,2,iidblock) + xshell0(jj,kk,3,iidblock)*xshell0(jj,kk,3,iidblock) + xshell0(jj,kk,4,iidblock)*xshell0(jj,kk,4,iidblock))
      
      !2 +1  0  0
	  udotc=xshell0(jj,kk,2,iidblock)*onecssq
	  f02(li+1,lj,lk)=p1*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*xshell0(jj,kk,5,iidblock)-cssq*(xshell0(jj,kk,6,iidblock)+xshell0(jj,kk,7,iidblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(xshell0(jj,kk,2,iidblock)+xshell0(jj,kk,3,iidblock))*onecssq
	  f08(li+1,lj,lk)=p2*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(jj,kk,5,iidblock)+qyy*xshell0(jj,kk,6,iidblock)-cssq*xshell0(jj,kk,7,iidblock)+two*qxy_7_8*xshell0(jj,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-xshell0(jj,kk,2,iidblock)+xshell0(jj,kk,3,iidblock))*onecssq
	  f10(li+1,lj,lk)=p2*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(jj,kk,5,iidblock)+qyy*xshell0(jj,kk,6,iidblock)-cssq*xshell0(jj,kk,7,iidblock)+two*qxy_9_10*xshell0(jj,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(xshell0(jj,kk,2,iidblock)+xshell0(jj,kk,4,iidblock))*onecssq
	  f16(li+1,lj,lk)=p2*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(jj,kk,5,iidblock)+qzz*xshell0(jj,kk,7,iidblock)-cssq*xshell0(jj,kk,6,iidblock)+two*qxz_15_16*xshell0(jj,kk,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-xshell0(jj,kk,2,iidblock)+xshell0(jj,kk,4,iidblock))*onecssq
	  f17(li+1,lj,lk)=p2*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*xshell0(jj,kk,5,iidblock)+qzz*xshell0(jj,kk,7,iidblock)-cssq*xshell0(jj,kk,6,iidblock)+two*qxz_17_18*xshell0(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(yshell1(ii,kk,2,iidblock)*yshell1(ii,kk,2,iidblock) + yshell1(ii,kk,3,iidblock)*yshell1(ii,kk,3,iidblock) + yshell1(ii,kk,4,iidblock)*yshell1(ii,kk,4,iidblock))
      
      !3 0 -1  0
	  udotc=yshell1(ii,kk,3,iidblock)*onecssq
	  f03(li,lj-1,lk)=p1*(yshell1(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*yshell1(ii,kk,6,iidblock)-cssq*(yshell1(ii,kk,5,iidblock)+yshell1(ii,kk,7,iidblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(yshell1(ii,kk,2,iidblock)+yshell1(ii,kk,3,iidblock))*onecssq
	  f07(li,lj-1,lk)=p2*(yshell1(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell1(ii,kk,5,iidblock)+qyy*yshell1(ii,kk,6,iidblock)-cssq*yshell1(ii,kk,7,iidblock)+two*qxy_7_8*yshell1(ii,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-yshell1(ii,kk,2,iidblock)+yshell1(ii,kk,3,iidblock))*onecssq
	  f10(li,lj-1,lk)=p2*(yshell1(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell1(ii,kk,5,iidblock)+qyy*yshell1(ii,kk,6,iidblock)-cssq*yshell1(ii,kk,7,iidblock)+two*qxy_9_10*yshell1(ii,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(yshell1(ii,kk,3,iidblock)+yshell1(ii,kk,4,iidblock))*onecssq
	  f11(li,lj-1,lk)=p2*(yshell1(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1(ii,kk,6,iidblock)+qzz*yshell1(ii,kk,7,iidblock)-cssq*yshell1(ii,kk,5,iidblock)+two*qyz_11_12*yshell1(ii,kk,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(yshell1(ii,kk,3,iidblock)-yshell1(ii,kk,4,iidblock))*onecssq
	  f13(li,lj-1,lk)=p2*(yshell1(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1(ii,kk,6,iidblock)+qzz*yshell1(ii,kk,7,iidblock)-cssq*yshell1(ii,kk,5,iidblock)+two*qyz_13_14*yshell1(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(yshell0(ii,kk,2,iidblock)*yshell0(ii,kk,2,iidblock) + yshell0(ii,kk,3,iidblock)*yshell0(ii,kk,3,iidblock) + yshell0(ii,kk,4,iidblock)*yshell0(ii,kk,4,iidblock))
      
      !4  0 +1  0
	  udotc=yshell0(ii,kk,3,iidblock)*onecssq
	  f04(li,lj+1,lk)=p1*(yshell0(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*yshell0(ii,kk,6,iidblock)-cssq*(yshell0(ii,kk,5,iidblock)+yshell0(ii,kk,7,iidblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(yshell0(ii,kk,2,iidblock)+yshell0(ii,kk,3,iidblock))*onecssq
	  f08(li,lj+1,lk)=p2*(yshell0(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell0(ii,kk,5,iidblock)+qyy*yshell0(ii,kk,6,iidblock)-cssq*yshell0(ii,kk,7,iidblock)+two*qxy_7_8*yshell0(ii,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-yshell0(ii,kk,2,iidblock)+yshell0(ii,kk,3,iidblock))*onecssq
	  f09(li,lj+1,lk)=p2*(yshell0(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell0(ii,kk,5,iidblock)+qyy*yshell0(ii,kk,6,iidblock)-cssq*yshell0(ii,kk,7,iidblock)+two*qxy_9_10*yshell0(ii,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(yshell0(ii,kk,3,iidblock)+yshell0(ii,kk,4,iidblock))*onecssq
	  f12(li,lj+1,lk)=p2*(yshell0(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0(ii,kk,6,iidblock)+qzz*yshell0(ii,kk,7,iidblock)-cssq*yshell0(ii,kk,5,iidblock)+two*qyz_11_12*yshell0(ii,kk,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(yshell0(ii,kk,3,iidblock)-yshell0(ii,kk,4,iidblock))*onecssq
	  f14(li,lj+1,lk)=p2*(yshell0(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0(ii,kk,6,iidblock)+qzz*yshell0(ii,kk,7,iidblock)-cssq*yshell0(ii,kk,5,iidblock)+two*qyz_13_14*yshell0(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(zshell1(ii,jj,2,iidblock)*zshell1(ii,jj,2,iidblock) + zshell1(ii,jj,3,iidblock)*zshell1(ii,jj,3,iidblock) + zshell1(ii,jj,4,iidblock)*zshell1(ii,jj,4,iidblock))
      
      !5  0  0 -1
	  udotc=zshell1(ii,jj,4,iidblock)*onecssq
	  f05(li,lj,lk-1)=p1*(zshell1(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*zshell1(ii,jj,7,iidblock)-cssq*(zshell1(ii,jj,5,iidblock)+zshell1(ii,jj,6,iidblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(zshell1(ii,jj,2,iidblock)+zshell1(ii,jj,4,iidblock))*onecssq
	  f15(li,lj,lk-1)=p2*(zshell1(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell1(ii,jj,5,iidblock)+qzz*zshell1(ii,jj,7,iidblock)-cssq*zshell1(ii,jj,6,iidblock)+two*qxz_15_16*zshell1(ii,jj,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-zshell1(ii,jj,2,iidblock)+zshell1(ii,jj,4,iidblock))*onecssq
	  f17(li,lj,lk-1)=p2*(zshell1(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell1(ii,jj,5,iidblock)+qzz*zshell1(ii,jj,7,iidblock)-cssq*zshell1(ii,jj,6,iidblock)+two*qxz_17_18*zshell1(ii,jj,9,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(zshell1(ii,jj,3,iidblock)+zshell1(ii,jj,4,iidblock))*onecssq
	  f11(li,lj,lk-1)=p2*(zshell1(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell1(ii,jj,6,iidblock)+qzz*zshell1(ii,jj,7,iidblock)-cssq*zshell1(ii,jj,5,iidblock)+two*qyz_11_12*zshell1(ii,jj,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(zshell1(ii,jj,3,iidblock)-zshell1(ii,jj,4,iidblock))*onecssq
	  f14(li,lj,lk-1)=p2*(zshell1(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell1(ii,jj,6,iidblock)+qzz*zshell1(ii,jj,7,iidblock)-cssq*zshell1(ii,jj,5,iidblock)+two*qyz_13_14*zshell1(ii,jj,10,iidblock)) &
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
      
      uu=halfonecssq*(zshell0(ii,jj,2,iidblock)*zshell0(ii,jj,2,iidblock) + zshell0(ii,jj,3,iidblock)*zshell0(ii,jj,3,iidblock) + zshell0(ii,jj,4,iidblock)*zshell0(ii,jj,4,iidblock))
      
      !6  0  0  +1
	  udotc=zshell0(ii,jj,4,iidblock)*onecssq
	  f06(li,lj,lk+1)=p1*(zshell0(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*zshell0(ii,jj,7,iidblock)-cssq*(zshell0(ii,jj,5,iidblock)+zshell0(ii,jj,6,iidblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(zshell0(ii,jj,2,iidblock)+zshell0(ii,jj,4,iidblock))*onecssq
	  f16(li,lj,lk+1)=p2*(zshell0(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell0(ii,jj,5,iidblock)+qzz*zshell0(ii,jj,7,iidblock)-cssq*zshell0(ii,jj,6,iidblock)+two*qxz_15_16*zshell0(ii,jj,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-zshell0(ii,jj,2,iidblock)+zshell0(ii,jj,4,iidblock))*onecssq
	  f18(li,lj,lk+1)=p2*(zshell0(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell0(ii,jj,5,iidblock)+qzz*zshell0(ii,jj,7,iidblock)-cssq*zshell0(ii,jj,6,iidblock)+two*qxz_17_18*zshell0(ii,jj,9,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(zshell0(ii,jj,3,iidblock)+zshell0(ii,jj,4,iidblock))*onecssq
	  f12(li,lj,lk+1)=p2*(zshell0(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell0(ii,jj,6,iidblock)+qzz*zshell0(ii,jj,7,iidblock)-cssq*zshell0(ii,jj,5,iidblock)+two*qyz_11_12*zshell0(ii,jj,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(zshell0(ii,jj,3,iidblock)-zshell0(ii,jj,4,iidblock))*onecssq
	  f13(li,lj,lk+1)=p2*(zshell0(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell0(ii,jj,6,iidblock)+qzz*zshell0(ii,jj,7,iidblock)-cssq*zshell0(ii,jj,5,iidblock)+two*qyz_13_14*zshell0(ii,jj,10,iidblock)) &
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
      
      uu=halfonecssq*(xshell1(jj,kk,2,iidblock)*xshell1(jj,kk,2,iidblock) + xshell1(jj,kk,3,iidblock)*xshell1(jj,kk,3,iidblock) + xshell1(jj,kk,4,iidblock)*xshell1(jj,kk,4,iidblock))
      
      !7 -1 -1  0
	  udotc=(xshell1(jj,kk,2,iidblock)+xshell1(jj,kk,3,iidblock))*onecssq
	  f07(li-1,lj-1,lk)=p2*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(jj,kk,5,iidblock)+qyy*xshell1(jj,kk,6,iidblock)-cssq*xshell1(jj,kk,7,iidblock)+two*qxy_7_8*xshell1(jj,kk,8,iidblock)) &
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
      
      uu=halfonecssq*(xshell1(jj,kk,2,iidblock)*xshell1(jj,kk,2,iidblock) + xshell1(jj,kk,3,iidblock)*xshell1(jj,kk,3,iidblock) + xshell1(jj,kk,4,iidblock)*xshell1(jj,kk,4,iidblock))
      
      !9  -1 +1 0
      udotc=(-xshell1(jj,kk,2,iidblock)+xshell1(jj,kk,3,iidblock))*onecssq
	  f09(li-1,lj+1,lk)=p2*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(jj,kk,5,iidblock)+qyy*xshell1(jj,kk,6,iidblock)-cssq*xshell1(jj,kk,7,iidblock)+two*qxy_9_10*xshell1(jj,kk,8,iidblock)) &
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
      
      uu=halfonecssq*(xshell1(jj,kk,2,iidblock)*xshell1(jj,kk,2,iidblock) + xshell1(jj,kk,3,iidblock)*xshell1(jj,kk,3,iidblock) + xshell1(jj,kk,4,iidblock)*xshell1(jj,kk,4,iidblock))
      
      !15  -1  0 -1
	  udotc=(xshell1(jj,kk,2,iidblock)+xshell1(jj,kk,4,iidblock))*onecssq
	  f15(li-1,lj,lk-1)=p2*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(jj,kk,5,iidblock)+qzz*xshell1(jj,kk,7,iidblock)-cssq*xshell1(jj,kk,6,iidblock)+two*qxz_15_16*xshell1(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(xshell1(jj,kk,2,iidblock)*xshell1(jj,kk,2,iidblock) + xshell1(jj,kk,3,iidblock)*xshell1(jj,kk,3,iidblock) + xshell1(jj,kk,4,iidblock)*xshell1(jj,kk,4,iidblock))
      
      !18   -1   0  +1
	  udotc=(-xshell1(jj,kk,2,iidblock)+xshell1(jj,kk,4,iidblock))*onecssq
	  f18(li-1,lj,lk+1)=p2*(xshell1(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(jj,kk,5,iidblock)+qzz*xshell1(jj,kk,7,iidblock)-cssq*xshell1(jj,kk,6,iidblock)+two*qxz_17_18*xshell1(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(xshell0(jj,kk,2,iidblock)*xshell0(jj,kk,2,iidblock) + xshell0(jj,kk,3,iidblock)*xshell0(jj,kk,3,iidblock) + xshell0(jj,kk,4,iidblock)*xshell0(jj,kk,4,iidblock))
      
      !10   +1 -1  0
	  udotc=(-xshell0(jj,kk,2,iidblock)+xshell0(jj,kk,3,iidblock))*onecssq
	  f10(li+1,lj-1,lk)=p2*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(jj,kk,5,iidblock)+qyy*xshell0(jj,kk,6,iidblock)-cssq*xshell0(jj,kk,7,iidblock)+two*qxy_9_10*xshell0(jj,kk,8,iidblock)) &
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
      
      uu=halfonecssq*(xshell0(jj,kk,2,iidblock)*xshell0(jj,kk,2,iidblock) + xshell0(jj,kk,3,iidblock)*xshell0(jj,kk,3,iidblock) + xshell0(jj,kk,4,iidblock)*xshell0(jj,kk,4,iidblock))
      
      !8 +1 +1  0
	  udotc=(xshell0(jj,kk,2,iidblock)+xshell0(jj,kk,3,iidblock))*onecssq
	  f08(li+1,lj+1,lk)=p2*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(jj,kk,5,iidblock)+qyy*xshell0(jj,kk,6,iidblock)-cssq*xshell0(jj,kk,7,iidblock)+two*qxy_7_8*xshell0(jj,kk,8,iidblock)) &
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
      
      uu=halfonecssq*(xshell0(jj,kk,2,iidblock)*xshell0(jj,kk,2,iidblock) + xshell0(jj,kk,3,iidblock)*xshell0(jj,kk,3,iidblock) + xshell0(jj,kk,4,iidblock)*xshell0(jj,kk,4,iidblock))
      
      !17  +1  0 -1
	  udotc=(-xshell0(jj,kk,2,iidblock)+xshell0(jj,kk,4,iidblock))*onecssq
	  f17(li+1,lj,lk-1)=p2*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*xshell0(jj,kk,5,iidblock)+qzz*xshell0(jj,kk,7,iidblock)-cssq*xshell0(jj,kk,6,iidblock)+two*qxz_17_18*xshell0(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(xshell0(jj,kk,2,iidblock)*xshell0(jj,kk,2,iidblock) + xshell0(jj,kk,3,iidblock)*xshell0(jj,kk,3,iidblock) + xshell0(jj,kk,4,iidblock)*xshell0(jj,kk,4,iidblock))
      
      !16  +1  0 +1
	  udotc=(xshell0(jj,kk,2,iidblock)+xshell0(jj,kk,4,iidblock))*onecssq
	  f16(li+1,lj,lk+1)=p2*(xshell0(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(jj,kk,5,iidblock)+qzz*xshell0(jj,kk,7,iidblock)-cssq*xshell0(jj,kk,6,iidblock)+two*qxz_15_16*xshell0(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(yshell1(ii,kk,2,iidblock)*yshell1(ii,kk,2,iidblock) + yshell1(ii,kk,3,iidblock)*yshell1(ii,kk,3,iidblock) + yshell1(ii,kk,4,iidblock)*yshell1(ii,kk,4,iidblock))
      
      !11  0  -1  -1
	  udotc=(yshell1(ii,kk,3,iidblock)+yshell1(ii,kk,4,iidblock))*onecssq
	  f11(li,lj-1,lk-1)=p2*(yshell1(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1(ii,kk,6,iidblock)+qzz*yshell1(ii,kk,7,iidblock)-cssq*yshell1(ii,kk,5,iidblock)+two*qyz_11_12*yshell1(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(yshell1(ii,kk,2,iidblock)*yshell1(ii,kk,2,iidblock) + yshell1(ii,kk,3,iidblock)*yshell1(ii,kk,3,iidblock) + yshell1(ii,kk,4,iidblock)*yshell1(ii,kk,4,iidblock))
      
      !13  0  -1   +1
	  udotc=(yshell1(ii,kk,3,iidblock)-yshell1(ii,kk,4,iidblock))*onecssq
	  f13(li,lj-1,lk+1)=p2*(yshell1(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1(ii,kk,6,iidblock)+qzz*yshell1(ii,kk,7,iidblock)-cssq*yshell1(ii,kk,5,iidblock)+two*qyz_13_14*yshell1(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(yshell0(ii,kk,2,iidblock)*yshell0(ii,kk,2,iidblock) + yshell0(ii,kk,3,iidblock)*yshell0(ii,kk,3,iidblock) + yshell0(ii,kk,4,iidblock)*yshell0(ii,kk,4,iidblock))
      
      !14  0  +1  -1
	  udotc=(yshell0(ii,kk,3,iidblock)-yshell0(ii,kk,4,iidblock))*onecssq
	  f14(li,lj+1,lk-1)=p2*(yshell0(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0(ii,kk,6,iidblock)+qzz*yshell0(ii,kk,7,iidblock)-cssq*yshell0(ii,kk,5,iidblock)+two*qyz_13_14*yshell0(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(yshell0(ii,kk,2,iidblock)*yshell0(ii,kk,2,iidblock) + yshell0(ii,kk,3,iidblock)*yshell0(ii,kk,3,iidblock) + yshell0(ii,kk,4,iidblock)*yshell0(ii,kk,4,iidblock))
      
      !12  0  +1  +1
	  udotc=(yshell0(ii,kk,3,iidblock)+yshell0(ii,kk,4,iidblock))*onecssq
	  f12(li,lj+1,lk+1)=p2*(yshell0(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0(ii,kk,6,iidblock)+qzz*yshell0(ii,kk,7,iidblock)-cssq*yshell0(ii,kk,5,iidblock)+two*qyz_11_12*yshell0(ii,kk,10,iidblock)) &
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
	hfields(i,j,k,1,myblock)=udotc
	!if(gk==1 .and. gi==3 .and. gj==3)write(*,*)gi,gj,gk,mystep,udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,2,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,3,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,4,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,5,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(i,j,k,6,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,7,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(i,j,k,8,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,9,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,10,myblock)=udotc
     
#ifdef PRESSCORR
	
	uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,5,myblock)=hfields(i,j,k,5,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(i,j,k,6,myblock)=hfields(i,j,k,6,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,7,myblock)=hfields(i,j,k,7,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(i,j,k,8,myblock)=hfields(i,j,k,8,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(i,j,k,9,myblock)=hfields(i,j,k,9,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(i,j,k,10,myblock)=hfields(i,j,k,10,myblock)-udotc
	
	    
    
#endif	
  
  ! Halo Faces
    if(li==1)then
      
      xshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(li==TILE_DIMx_d)then
      
      xshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
	  
    endif

    if(lj==1)then
      
      yshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lj==TILE_DIMy_d)then
      
      yshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    if(lk==1)then
      
      zshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lk==TILE_DIMz_d)then
      
      zshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
  
  return

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
	
!	 myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
!    zblock=(idblock-1)/nxyblock +1
!    yblock=((idblock-1)-(zblock-1)*nxyblock)/nxblock +1
!    xblock=(idblock-1)-(zblock-1)*nxyblock-(yblock-1)*nxblock +1
	
	!myblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
	iidblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	!iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	if(myblock.ne.iidblock)write(*,*)'cazzo amaro streamcoll_shared',gi,gj,gk,myblock,iidblock
    
        
    uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(hfields(i,j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock)-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock)-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
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
      
      uu=halfonecssq*(xshell1h(jj,kk,2,iidblock)*xshell1h(jj,kk,2,iidblock) + xshell1h(jj,kk,3,iidblock)*xshell1h(jj,kk,3,iidblock) + xshell1h(jj,kk,4,iidblock)*xshell1h(jj,kk,4,iidblock))
      
      !1 -1  0  0
	  udotc=xshell1h(jj,kk,2,iidblock)*onecssq
	  f01(li-1,lj,lk)=p1*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*xshell1h(jj,kk,5,iidblock)-cssq*(xshell1h(jj,kk,6,iidblock)+xshell1h(jj,kk,7,iidblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(xshell1h(jj,kk,2,iidblock)+xshell1h(jj,kk,3,iidblock))*onecssq
	  f07(li-1,lj,lk)=p2*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(jj,kk,5,iidblock)+qyy*xshell1h(jj,kk,6,iidblock)-cssq*xshell1h(jj,kk,7,iidblock)+two*qxy_7_8*xshell1h(jj,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-xshell1h(jj,kk,2,iidblock)+xshell1h(jj,kk,3,iidblock))*onecssq
	  f09(li-1,lj,lk)=p2*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(jj,kk,5,iidblock)+qyy*xshell1h(jj,kk,6,iidblock)-cssq*xshell1h(jj,kk,7,iidblock)+two*qxy_9_10*xshell1h(jj,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(xshell1h(jj,kk,2,iidblock)+xshell1h(jj,kk,4,iidblock))*onecssq
	  f15(li-1,lj,lk)=p2*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(jj,kk,5,iidblock)+qzz*xshell1h(jj,kk,7,iidblock)-cssq*xshell1h(jj,kk,6,iidblock)+two*qxz_15_16*xshell1h(jj,kk,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-xshell1h(jj,kk,2,iidblock)+xshell1h(jj,kk,4,iidblock))*onecssq
	  f18(li-1,lj,lk)=p2*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(jj,kk,5,iidblock)+qzz*xshell1h(jj,kk,7,iidblock)-cssq*xshell1h(jj,kk,6,iidblock)+two*qxz_17_18*xshell1h(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(xshell0h(jj,kk,2,iidblock)*xshell0h(jj,kk,2,iidblock) + xshell0h(jj,kk,3,iidblock)*xshell0h(jj,kk,3,iidblock) + xshell0h(jj,kk,4,iidblock)*xshell0h(jj,kk,4,iidblock))
      
      !2 +1  0  0
	  udotc=xshell0h(jj,kk,2,iidblock)*onecssq
	  f02(li+1,lj,lk)=p1*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*xshell0h(jj,kk,5,iidblock)-cssq*(xshell0h(jj,kk,6,iidblock)+xshell0h(jj,kk,7,iidblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(xshell0h(jj,kk,2,iidblock)+xshell0h(jj,kk,3,iidblock))*onecssq
	  f08(li+1,lj,lk)=p2*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(jj,kk,5,iidblock)+qyy*xshell0h(jj,kk,6,iidblock)-cssq*xshell0h(jj,kk,7,iidblock)+two*qxy_7_8*xshell0h(jj,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-xshell0h(jj,kk,2,iidblock)+xshell0h(jj,kk,3,iidblock))*onecssq
	  f10(li+1,lj,lk)=p2*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(jj,kk,5,iidblock)+qyy*xshell0h(jj,kk,6,iidblock)-cssq*xshell0h(jj,kk,7,iidblock)+two*qxy_9_10*xshell0h(jj,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(xshell0h(jj,kk,2,iidblock)+xshell0h(jj,kk,4,iidblock))*onecssq
	  f16(li+1,lj,lk)=p2*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(jj,kk,5,iidblock)+qzz*xshell0h(jj,kk,7,iidblock)-cssq*xshell0h(jj,kk,6,iidblock)+two*qxz_15_16*xshell0h(jj,kk,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-xshell0h(jj,kk,2,iidblock)+xshell0h(jj,kk,4,iidblock))*onecssq
	  f17(li+1,lj,lk)=p2*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*xshell0h(jj,kk,5,iidblock)+qzz*xshell0h(jj,kk,7,iidblock)-cssq*xshell0h(jj,kk,6,iidblock)+two*qxz_17_18*xshell0h(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(yshell1h(ii,kk,2,iidblock)*yshell1h(ii,kk,2,iidblock) + yshell1h(ii,kk,3,iidblock)*yshell1h(ii,kk,3,iidblock) + yshell1h(ii,kk,4,iidblock)*yshell1h(ii,kk,4,iidblock))
      
      !3 0 -1  0
	  udotc=yshell1h(ii,kk,3,iidblock)*onecssq
	  f03(li,lj-1,lk)=p1*(yshell1h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*yshell1h(ii,kk,6,iidblock)-cssq*(yshell1h(ii,kk,5,iidblock)+yshell1h(ii,kk,7,iidblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(yshell1h(ii,kk,2,iidblock)+yshell1h(ii,kk,3,iidblock))*onecssq
	  f07(li,lj-1,lk)=p2*(yshell1h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell1h(ii,kk,5,iidblock)+qyy*yshell1h(ii,kk,6,iidblock)-cssq*yshell1h(ii,kk,7,iidblock)+two*qxy_7_8*yshell1h(ii,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-yshell1h(ii,kk,2,iidblock)+yshell1h(ii,kk,3,iidblock))*onecssq
	  f10(li,lj-1,lk)=p2*(yshell1h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell1h(ii,kk,5,iidblock)+qyy*yshell1h(ii,kk,6,iidblock)-cssq*yshell1h(ii,kk,7,iidblock)+two*qxy_9_10*yshell1h(ii,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(yshell1h(ii,kk,3,iidblock)+yshell1h(ii,kk,4,iidblock))*onecssq
	  f11(li,lj-1,lk)=p2*(yshell1h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1h(ii,kk,6,iidblock)+qzz*yshell1h(ii,kk,7,iidblock)-cssq*yshell1h(ii,kk,5,iidblock)+two*qyz_11_12*yshell1h(ii,kk,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(yshell1h(ii,kk,3,iidblock)-yshell1h(ii,kk,4,iidblock))*onecssq
	  f13(li,lj-1,lk)=p2*(yshell1h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1h(ii,kk,6,iidblock)+qzz*yshell1h(ii,kk,7,iidblock)-cssq*yshell1h(ii,kk,5,iidblock)+two*qyz_13_14*yshell1h(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(yshell0h(ii,kk,2,iidblock)*yshell0h(ii,kk,2,iidblock) + yshell0h(ii,kk,3,iidblock)*yshell0h(ii,kk,3,iidblock) + yshell0h(ii,kk,4,iidblock)*yshell0h(ii,kk,4,iidblock))
      
      !4  0 +1  0
	  udotc=yshell0h(ii,kk,3,iidblock)*onecssq
	  f04(li,lj+1,lk)=p1*(yshell0h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*yshell0h(ii,kk,6,iidblock)-cssq*(yshell0h(ii,kk,5,iidblock)+yshell0h(ii,kk,7,iidblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(yshell0h(ii,kk,2,iidblock)+yshell0h(ii,kk,3,iidblock))*onecssq
	  f08(li,lj+1,lk)=p2*(yshell0h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell0h(ii,kk,5,iidblock)+qyy*yshell0h(ii,kk,6,iidblock)-cssq*yshell0h(ii,kk,7,iidblock)+two*qxy_7_8*yshell0h(ii,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-yshell0h(ii,kk,2,iidblock)+yshell0h(ii,kk,3,iidblock))*onecssq
	  f09(li,lj+1,lk)=p2*(yshell0h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell0h(ii,kk,5,iidblock)+qyy*yshell0h(ii,kk,6,iidblock)-cssq*yshell0h(ii,kk,7,iidblock)+two*qxy_9_10*yshell0h(ii,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(yshell0h(ii,kk,3,iidblock)+yshell0h(ii,kk,4,iidblock))*onecssq
	  f12(li,lj+1,lk)=p2*(yshell0h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0h(ii,kk,6,iidblock)+qzz*yshell0h(ii,kk,7,iidblock)-cssq*yshell0h(ii,kk,5,iidblock)+two*qyz_11_12*yshell0h(ii,kk,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(yshell0h(ii,kk,3,iidblock)-yshell0h(ii,kk,4,iidblock))*onecssq
	  f14(li,lj+1,lk)=p2*(yshell0h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0h(ii,kk,6,iidblock)+qzz*yshell0h(ii,kk,7,iidblock)-cssq*yshell0h(ii,kk,5,iidblock)+two*qyz_13_14*yshell0h(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(zshell1h(ii,jj,2,iidblock)*zshell1h(ii,jj,2,iidblock) + zshell1h(ii,jj,3,iidblock)*zshell1h(ii,jj,3,iidblock) + zshell1h(ii,jj,4,iidblock)*zshell1h(ii,jj,4,iidblock))
      
      !5  0  0 -1
	  udotc=zshell1h(ii,jj,4,iidblock)*onecssq
	  f05(li,lj,lk-1)=p1*(zshell1h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*zshell1h(ii,jj,7,iidblock)-cssq*(zshell1h(ii,jj,5,iidblock)+zshell1h(ii,jj,6,iidblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(zshell1h(ii,jj,2,iidblock)+zshell1h(ii,jj,4,iidblock))*onecssq
	  f15(li,lj,lk-1)=p2*(zshell1h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell1h(ii,jj,5,iidblock)+qzz*zshell1h(ii,jj,7,iidblock)-cssq*zshell1h(ii,jj,6,iidblock)+two*qxz_15_16*zshell1h(ii,jj,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-zshell1h(ii,jj,2,iidblock)+zshell1h(ii,jj,4,iidblock))*onecssq
	  f17(li,lj,lk-1)=p2*(zshell1h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell1h(ii,jj,5,iidblock)+qzz*zshell1h(ii,jj,7,iidblock)-cssq*zshell1h(ii,jj,6,iidblock)+two*qxz_17_18*zshell1h(ii,jj,9,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(zshell1h(ii,jj,3,iidblock)+zshell1h(ii,jj,4,iidblock))*onecssq
	  f11(li,lj,lk-1)=p2*(zshell1h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell1h(ii,jj,6,iidblock)+qzz*zshell1h(ii,jj,7,iidblock)-cssq*zshell1h(ii,jj,5,iidblock)+two*qyz_11_12*zshell1h(ii,jj,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(zshell1h(ii,jj,3,iidblock)-zshell1h(ii,jj,4,iidblock))*onecssq
	  f14(li,lj,lk-1)=p2*(zshell1h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell1h(ii,jj,6,iidblock)+qzz*zshell1h(ii,jj,7,iidblock)-cssq*zshell1h(ii,jj,5,iidblock)+two*qyz_13_14*zshell1h(ii,jj,10,iidblock)) &
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
      
      uu=halfonecssq*(zshell0h(ii,jj,2,iidblock)*zshell0h(ii,jj,2,iidblock) + zshell0h(ii,jj,3,iidblock)*zshell0h(ii,jj,3,iidblock) + zshell0h(ii,jj,4,iidblock)*zshell0h(ii,jj,4,iidblock))
      
      !6  0  0  +1
	  udotc=zshell0h(ii,jj,4,iidblock)*onecssq
	  f06(li,lj,lk+1)=p1*(zshell0h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*zshell0h(ii,jj,7,iidblock)-cssq*(zshell0h(ii,jj,5,iidblock)+zshell0h(ii,jj,6,iidblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(zshell0h(ii,jj,2,iidblock)+zshell0h(ii,jj,4,iidblock))*onecssq
	  f16(li,lj,lk+1)=p2*(zshell0h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell0h(ii,jj,5,iidblock)+qzz*zshell0h(ii,jj,7,iidblock)-cssq*zshell0h(ii,jj,6,iidblock)+two*qxz_15_16*zshell0h(ii,jj,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-zshell0h(ii,jj,2,iidblock)+zshell0h(ii,jj,4,iidblock))*onecssq
	  f18(li,lj,lk+1)=p2*(zshell0h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell0h(ii,jj,5,iidblock)+qzz*zshell0h(ii,jj,7,iidblock)-cssq*zshell0h(ii,jj,6,iidblock)+two*qxz_17_18*zshell0h(ii,jj,9,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(zshell0h(ii,jj,3,iidblock)+zshell0h(ii,jj,4,iidblock))*onecssq
	  f12(li,lj,lk+1)=p2*(zshell0h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell0h(ii,jj,6,iidblock)+qzz*zshell0h(ii,jj,7,iidblock)-cssq*zshell0h(ii,jj,5,iidblock)+two*qyz_11_12*zshell0h(ii,jj,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(zshell0h(ii,jj,3,iidblock)-zshell0h(ii,jj,4,iidblock))*onecssq
	  f13(li,lj,lk+1)=p2*(zshell0h(ii,jj,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell0h(ii,jj,6,iidblock)+qzz*zshell0h(ii,jj,7,iidblock)-cssq*zshell0h(ii,jj,5,iidblock)+two*qyz_13_14*zshell0h(ii,jj,10,iidblock)) &
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
      
      uu=halfonecssq*(xshell1h(jj,kk,2,iidblock)*xshell1h(jj,kk,2,iidblock) + xshell1h(jj,kk,3,iidblock)*xshell1h(jj,kk,3,iidblock) + xshell1h(jj,kk,4,iidblock)*xshell1h(jj,kk,4,iidblock))
      
      !7 -1 -1  0
	  udotc=(xshell1h(jj,kk,2,iidblock)+xshell1h(jj,kk,3,iidblock))*onecssq
	  f07(li-1,lj-1,lk)=p2*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(jj,kk,5,iidblock)+qyy*xshell1h(jj,kk,6,iidblock)-cssq*xshell1h(jj,kk,7,iidblock)+two*qxy_7_8*xshell1h(jj,kk,8,iidblock)) &
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
      
      uu=halfonecssq*(xshell1h(jj,kk,2,iidblock)*xshell1h(jj,kk,2,iidblock) + xshell1h(jj,kk,3,iidblock)*xshell1h(jj,kk,3,iidblock) + xshell1h(jj,kk,4,iidblock)*xshell1h(jj,kk,4,iidblock))
      
      !9  -1 +1 0
      udotc=(-xshell1h(jj,kk,2,iidblock)+xshell1h(jj,kk,3,iidblock))*onecssq
	  f09(li-1,lj+1,lk)=p2*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(jj,kk,5,iidblock)+qyy*xshell1h(jj,kk,6,iidblock)-cssq*xshell1h(jj,kk,7,iidblock)+two*qxy_9_10*xshell1h(jj,kk,8,iidblock)) &
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
      
      uu=halfonecssq*(xshell1h(jj,kk,2,iidblock)*xshell1h(jj,kk,2,iidblock) + xshell1h(jj,kk,3,iidblock)*xshell1h(jj,kk,3,iidblock) + xshell1h(jj,kk,4,iidblock)*xshell1h(jj,kk,4,iidblock))
      
      !15  -1  0 -1
	  udotc=(xshell1h(jj,kk,2,iidblock)+xshell1h(jj,kk,4,iidblock))*onecssq
	  f15(li-1,lj,lk-1)=p2*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(jj,kk,5,iidblock)+qzz*xshell1h(jj,kk,7,iidblock)-cssq*xshell1h(jj,kk,6,iidblock)+two*qxz_15_16*xshell1h(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(xshell1h(jj,kk,2,iidblock)*xshell1h(jj,kk,2,iidblock) + xshell1h(jj,kk,3,iidblock)*xshell1h(jj,kk,3,iidblock) + xshell1h(jj,kk,4,iidblock)*xshell1h(jj,kk,4,iidblock))
      
      !18   -1   0  +1
	  udotc=(-xshell1h(jj,kk,2,iidblock)+xshell1h(jj,kk,4,iidblock))*onecssq
	  f18(li-1,lj,lk+1)=p2*(xshell1h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(jj,kk,5,iidblock)+qzz*xshell1h(jj,kk,7,iidblock)-cssq*xshell1h(jj,kk,6,iidblock)+two*qxz_17_18*xshell1h(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(xshell0h(jj,kk,2,iidblock)*xshell0h(jj,kk,2,iidblock) + xshell0h(jj,kk,3,iidblock)*xshell0h(jj,kk,3,iidblock) + xshell0h(jj,kk,4,iidblock)*xshell0h(jj,kk,4,iidblock))
      
      !10   +1 -1  0
	  udotc=(-xshell0h(jj,kk,2,iidblock)+xshell0h(jj,kk,3,iidblock))*onecssq
	  f10(li+1,lj-1,lk)=p2*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(jj,kk,5,iidblock)+qyy*xshell0h(jj,kk,6,iidblock)-cssq*xshell0h(jj,kk,7,iidblock)+two*qxy_9_10*xshell0h(jj,kk,8,iidblock)) &
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
      
      uu=halfonecssq*(xshell0h(jj,kk,2,iidblock)*xshell0h(jj,kk,2,iidblock) + xshell0h(jj,kk,3,iidblock)*xshell0h(jj,kk,3,iidblock) + xshell0h(jj,kk,4,iidblock)*xshell0h(jj,kk,4,iidblock))
      
      !8 +1 +1  0
	  udotc=(xshell0h(jj,kk,2,iidblock)+xshell0h(jj,kk,3,iidblock))*onecssq
	  f08(li+1,lj+1,lk)=p2*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(jj,kk,5,iidblock)+qyy*xshell0h(jj,kk,6,iidblock)-cssq*xshell0h(jj,kk,7,iidblock)+two*qxy_7_8*xshell0h(jj,kk,8,iidblock)) &
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
      
      uu=halfonecssq*(xshell0h(jj,kk,2,iidblock)*xshell0h(jj,kk,2,iidblock) + xshell0h(jj,kk,3,iidblock)*xshell0h(jj,kk,3,iidblock) + xshell0h(jj,kk,4,iidblock)*xshell0h(jj,kk,4,iidblock))
      
      !17  +1  0 -1
	  udotc=(-xshell0h(jj,kk,2,iidblock)+xshell0h(jj,kk,4,iidblock))*onecssq
	  f17(li+1,lj,lk-1)=p2*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*xshell0h(jj,kk,5,iidblock)+qzz*xshell0h(jj,kk,7,iidblock)-cssq*xshell0h(jj,kk,6,iidblock)+two*qxz_17_18*xshell0h(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(xshell0h(jj,kk,2,iidblock)*xshell0h(jj,kk,2,iidblock) + xshell0h(jj,kk,3,iidblock)*xshell0h(jj,kk,3,iidblock) + xshell0h(jj,kk,4,iidblock)*xshell0h(jj,kk,4,iidblock))
      
      !16  +1  0 +1
	  udotc=(xshell0h(jj,kk,2,iidblock)+xshell0h(jj,kk,4,iidblock))*onecssq
	  f16(li+1,lj,lk+1)=p2*(xshell0h(jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(jj,kk,5,iidblock)+qzz*xshell0h(jj,kk,7,iidblock)-cssq*xshell0h(jj,kk,6,iidblock)+two*qxz_15_16*xshell0h(jj,kk,9,iidblock)) &
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
      
      uu=halfonecssq*(yshell1h(ii,kk,2,iidblock)*yshell1h(ii,kk,2,iidblock) + yshell1h(ii,kk,3,iidblock)*yshell1h(ii,kk,3,iidblock) + yshell1h(ii,kk,4,iidblock)*yshell1h(ii,kk,4,iidblock))
      
      !11  0  -1  -1
	  udotc=(yshell1h(ii,kk,3,iidblock)+yshell1h(ii,kk,4,iidblock))*onecssq
	  f11(li,lj-1,lk-1)=p2*(yshell1h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1h(ii,kk,6,iidblock)+qzz*yshell1h(ii,kk,7,iidblock)-cssq*yshell1h(ii,kk,5,iidblock)+two*qyz_11_12*yshell1h(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(yshell1h(ii,kk,2,iidblock)*yshell1h(ii,kk,2,iidblock) + yshell1h(ii,kk,3,iidblock)*yshell1h(ii,kk,3,iidblock) + yshell1h(ii,kk,4,iidblock)*yshell1h(ii,kk,4,iidblock))
      
      !13  0  -1   +1
	  udotc=(yshell1h(ii,kk,3,iidblock)-yshell1h(ii,kk,4,iidblock))*onecssq
	  f13(li,lj-1,lk+1)=p2*(yshell1h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1h(ii,kk,6,iidblock)+qzz*yshell1h(ii,kk,7,iidblock)-cssq*yshell1h(ii,kk,5,iidblock)+two*qyz_13_14*yshell1h(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(yshell0h(ii,kk,2,iidblock)*yshell0h(ii,kk,2,iidblock) + yshell0h(ii,kk,3,iidblock)*yshell0h(ii,kk,3,iidblock) + yshell0h(ii,kk,4,iidblock)*yshell0h(ii,kk,4,iidblock))
      
      !14  0  +1  -1
	  udotc=(yshell0h(ii,kk,3,iidblock)-yshell0h(ii,kk,4,iidblock))*onecssq
	  f14(li,lj+1,lk-1)=p2*(yshell0h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0h(ii,kk,6,iidblock)+qzz*yshell0h(ii,kk,7,iidblock)-cssq*yshell0h(ii,kk,5,iidblock)+two*qyz_13_14*yshell0h(ii,kk,10,iidblock)) &
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
      
      uu=halfonecssq*(yshell0h(ii,kk,2,iidblock)*yshell0h(ii,kk,2,iidblock) + yshell0h(ii,kk,3,iidblock)*yshell0h(ii,kk,3,iidblock) + yshell0h(ii,kk,4,iidblock)*yshell0h(ii,kk,4,iidblock))
      
      !12  0  +1  +1
	  udotc=(yshell0h(ii,kk,3,iidblock)+yshell0h(ii,kk,4,iidblock))*onecssq
	  f12(li,lj+1,lk+1)=p2*(yshell0h(ii,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0h(ii,kk,6,iidblock)+qzz*yshell0h(ii,kk,7,iidblock)-cssq*yshell0h(ii,kk,5,iidblock)+two*qyz_11_12*yshell0h(ii,kk,10,iidblock)) &
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
	hfields(i,j,k,1,myblock)=udotc
	!if(gk==1 .and. gi==3 .and. gj==3)write(*,*)gi,gj,gk,mystep,udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,2,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,3,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,4,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,5,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(i,j,k,6,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,7,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(i,j,k,8,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,9,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,10,myblock)=udotc
     
#ifdef PRESSCORR
	
	uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,5,myblock)=hfields(i,j,k,5,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(i,j,k,6,myblock)=hfields(i,j,k,6,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,7,myblock)=hfields(i,j,k,7,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(i,j,k,8,myblock)=hfields(i,j,k,8,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(i,j,k,9,myblock)=hfields(i,j,k,9,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(i,j,k,10,myblock)=hfields(i,j,k,10,myblock)-udotc
	
	    
    
#endif	
  
  ! Halo Faces
    if(li==1)then
      
      xshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(li==TILE_DIMx_d)then
      
      xshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
	  
    endif

    if(lj==1)then
      
      yshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lj==TILE_DIMy_d)then
      
      yshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    if(lk==1)then
      
      zshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lk==TILE_DIMz_d)then
      
      zshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
  
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
    
    !inner cube
    if(li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d .and. &
     lk>=1 .and. lk<=TILE_DIMz_d)then
     
    uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
     + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(hfields(i,j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
	
    endif
    
    ! Halo Faces
    if(li==0)then
      
    uu=halfonecssq*(xshell1(j,k,2,myblock)*xshell1(j,k,2,myblock) &
     + xshell1(j,k,3,myblock)*xshell1(j,k,3,myblock) + xshell1(j,k,4,myblock)*xshell1(j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(xshell1(j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(xshell1(j,k,6,myblock)+xshell1(j,k,5,myblock)+xshell1(j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=xshell1(j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*xshell1(j,k,5,myblock) &
	 -cssq*(xshell1(j,k,6,myblock)+xshell1(j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*xshell1(j,k,5,myblock) &
	 -cssq*(xshell1(j,k,6,myblock)+xshell1(j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=xshell1(j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*xshell1(j,k,6,myblock) &
	 -cssq*(xshell1(j,k,5,myblock)+xshell1(j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*xshell1(j,k,6,myblock) &
	 -cssq*(xshell1(j,k,5,myblock)+xshell1(j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=xshell1(j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*xshell1(j,k,7,myblock) &
	 -cssq*(xshell1(j,k,5,myblock)+xshell1(j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*xshell1(j,k,7,myblock) &
	 -cssq*(xshell1(j,k,5,myblock)+xshell1(j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(xshell1(j,k,2,myblock)+xshell1(j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qyy*xshell1(j,k,6,myblock) &
	 -cssq*xshell1(j,k,7,myblock)+two*qxy_7_8*xshell1(j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qyy*xshell1(j,k,6,myblock) &
	 -cssq*xshell1(j,k,7,myblock)+two*qxy_7_8*xshell1(j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-xshell1(j,k,2,myblock)+xshell1(j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qyy*xshell1(j,k,6,myblock) &
	 -cssq*xshell1(j,k,7,myblock)+two*qxy_9_10*xshell1(j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qyy*xshell1(j,k,6,myblock) &
	 -cssq*xshell1(j,k,7,myblock)+two*qxy_9_10*xshell1(j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(xshell1(j,k,2,myblock)+xshell1(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qzz*xshell1(j,k,7,myblock) &
	 -cssq*xshell1(j,k,6,myblock)+two*qxz_15_16*xshell1(j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qzz*xshell1(j,k,7,myblock) &
	 -cssq*xshell1(j,k,6,myblock)+two*qxz_15_16*xshell1(j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-xshell1(j,k,2,myblock)+xshell1(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qzz*xshell1(j,k,7,myblock) &
	 -cssq*xshell1(j,k,6,myblock)+two*qxz_17_18*xshell1(j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qzz*xshell1(j,k,7,myblock) &
	 -cssq*xshell1(j,k,6,myblock)+two*qxz_17_18*xshell1(j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(xshell1(j,k,3,myblock)+xshell1(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell1(j,k,6,myblock)+qzz*xshell1(j,k,7,myblock) &
	 -cssq*xshell1(j,k,5,myblock)+two*qyz_11_12*xshell1(j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell1(j,k,6,myblock)+qzz*xshell1(j,k,7,myblock) &
	 -cssq*xshell1(j,k,5,myblock)+two*qyz_11_12*xshell1(j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(xshell1(j,k,3,myblock)-xshell1(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell1(j,k,6,myblock)+qzz*xshell1(j,k,7,myblock) &
	 -cssq*xshell1(j,k,5,myblock)+two*qyz_13_14*xshell1(j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell1(j,k,6,myblock)+qzz*xshell1(j,k,7,myblock) &
	 -cssq*xshell1(j,k,5,myblock)+two*qyz_13_14*xshell1(j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
    endif
    if(li==TILE_DIMx_d+1)then
      
    uu=halfonecssq*(xshell0(j,k,2,myblock)*xshell0(j,k,2,myblock) &
     + xshell0(j,k,3,myblock)*xshell0(j,k,3,myblock) + xshell0(j,k,4,myblock)*xshell0(j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(xshell0(j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(xshell0(j,k,6,myblock)+xshell0(j,k,5,myblock)+xshell0(j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=xshell0(j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*xshell0(j,k,5,myblock) &
	 -cssq*(xshell0(j,k,6,myblock)+xshell0(j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*xshell0(j,k,5,myblock) &
	 -cssq*(xshell0(j,k,6,myblock)+xshell0(j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=xshell0(j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*xshell0(j,k,6,myblock) &
	 -cssq*(xshell0(j,k,5,myblock)+xshell0(j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*xshell0(j,k,6,myblock) &
	 -cssq*(xshell0(j,k,5,myblock)+xshell0(j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=xshell0(j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*xshell0(j,k,7,myblock) &
	 -cssq*(xshell0(j,k,5,myblock)+xshell0(j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*xshell0(j,k,7,myblock) &
	 -cssq*(xshell0(j,k,5,myblock)+xshell0(j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(xshell0(j,k,2,myblock)+xshell0(j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qyy*xshell0(j,k,6,myblock) &
	 -cssq*xshell0(j,k,7,myblock)+two*qxy_7_8*xshell0(j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qyy*xshell0(j,k,6,myblock) &
	 -cssq*xshell0(j,k,7,myblock)+two*qxy_7_8*xshell0(j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-xshell0(j,k,2,myblock)+xshell0(j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qyy*xshell0(j,k,6,myblock) &
	 -cssq*xshell0(j,k,7,myblock)+two*qxy_9_10*xshell0(j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qyy*xshell0(j,k,6,myblock) &
	 -cssq*xshell0(j,k,7,myblock)+two*qxy_9_10*xshell0(j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(xshell0(j,k,2,myblock)+xshell0(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qzz*xshell0(j,k,7,myblock) &
	 -cssq*xshell0(j,k,6,myblock)+two*qxz_15_16*xshell0(j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qzz*xshell0(j,k,7,myblock) &
	 -cssq*xshell0(j,k,6,myblock)+two*qxz_15_16*xshell0(j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-xshell0(j,k,2,myblock)+xshell0(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qzz*xshell0(j,k,7,myblock) &
	 -cssq*xshell0(j,k,6,myblock)+two*qxz_17_18*xshell0(j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qzz*xshell0(j,k,7,myblock) &
	 -cssq*xshell0(j,k,6,myblock)+two*qxz_17_18*xshell0(j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(xshell0(j,k,3,myblock)+xshell0(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell0(j,k,6,myblock)+qzz*xshell0(j,k,7,myblock) &
	 -cssq*xshell0(j,k,5,myblock)+two*qyz_11_12*xshell0(j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell0(j,k,6,myblock)+qzz*xshell0(j,k,7,myblock) &
	 -cssq*xshell0(j,k,5,myblock)+two*qyz_11_12*xshell0(j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(xshell0(j,k,3,myblock)-xshell0(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell0(j,k,6,myblock)+qzz*xshell0(j,k,7,myblock) &
	 -cssq*xshell0(j,k,5,myblock)+two*qyz_13_14*xshell0(j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell0(j,k,6,myblock)+qzz*xshell0(j,k,7,myblock) &
	 -cssq*xshell0(j,k,5,myblock)+two*qyz_13_14*xshell0(j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
	  
    endif

    if(lj==0 .and. li>=1 .and. li<=TILE_DIMx_d)then
      
    uu=halfonecssq*(yshell1(i,k,2,myblock)*yshell1(i,k,2,myblock) &
     + yshell1(i,k,3,myblock)*yshell1(i,k,3,myblock) + yshell1(i,k,4,myblock)*yshell1(i,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(yshell1(i,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(yshell1(i,k,6,myblock)+yshell1(i,k,5,myblock)+yshell1(i,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=yshell1(i,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*yshell1(i,k,5,myblock) &
	 -cssq*(yshell1(i,k,6,myblock)+yshell1(i,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*yshell1(i,k,5,myblock) &
	 -cssq*(yshell1(i,k,6,myblock)+yshell1(i,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=yshell1(i,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*yshell1(i,k,6,myblock) &
	 -cssq*(yshell1(i,k,5,myblock)+yshell1(i,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*yshell1(i,k,6,myblock) &
	 -cssq*(yshell1(i,k,5,myblock)+yshell1(i,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=yshell1(i,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*yshell1(i,k,7,myblock) &
	 -cssq*(yshell1(i,k,5,myblock)+yshell1(i,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*yshell1(i,k,7,myblock) &
	 -cssq*(yshell1(i,k,5,myblock)+yshell1(i,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(yshell1(i,k,2,myblock)+yshell1(i,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qyy*yshell1(i,k,6,myblock) &
	 -cssq*yshell1(i,k,7,myblock)+two*qxy_7_8*yshell1(i,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qyy*yshell1(i,k,6,myblock) &
	 -cssq*yshell1(i,k,7,myblock)+two*qxy_7_8*yshell1(i,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-yshell1(i,k,2,myblock)+yshell1(i,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qyy*yshell1(i,k,6,myblock) &
	 -cssq*yshell1(i,k,7,myblock)+two*qxy_9_10*yshell1(i,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qyy*yshell1(i,k,6,myblock) &
	 -cssq*yshell1(i,k,7,myblock)+two*qxy_9_10*yshell1(i,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(yshell1(i,k,2,myblock)+yshell1(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qzz*yshell1(i,k,7,myblock) &
	 -cssq*yshell1(i,k,6,myblock)+two*qxz_15_16*yshell1(i,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qzz*yshell1(i,k,7,myblock) &
	 -cssq*yshell1(i,k,6,myblock)+two*qxz_15_16*yshell1(i,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-yshell1(i,k,2,myblock)+yshell1(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qzz*yshell1(i,k,7,myblock) &
	 -cssq*yshell1(i,k,6,myblock)+two*qxz_17_18*yshell1(i,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qzz*yshell1(i,k,7,myblock) &
	 -cssq*yshell1(i,k,6,myblock)+two*qxz_17_18*yshell1(i,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(yshell1(i,k,3,myblock)+yshell1(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell1(i,k,6,myblock)+qzz*yshell1(i,k,7,myblock) &
	 -cssq*yshell1(i,k,5,myblock)+two*qyz_11_12*yshell1(i,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell1(i,k,6,myblock)+qzz*yshell1(i,k,7,myblock) &
	 -cssq*yshell1(i,k,5,myblock)+two*qyz_11_12*yshell1(i,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(yshell1(i,k,3,myblock)-yshell1(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell1(i,k,6,myblock)+qzz*yshell1(i,k,7,myblock) &
	 -cssq*yshell1(i,k,5,myblock)+two*qyz_13_14*yshell1(i,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell1(i,k,6,myblock)+qzz*yshell1(i,k,7,myblock) &
	 -cssq*yshell1(i,k,5,myblock)+two*qyz_13_14*yshell1(i,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
    endif
    
    if(lj==TILE_DIMy_d+1  .and. li>=1 .and. li<=TILE_DIMx_d)then
      
    uu=halfonecssq*(yshell0(i,k,2,myblock)*yshell0(i,k,2,myblock) &
     + yshell0(i,k,3,myblock)*yshell0(i,k,3,myblock) + yshell0(i,k,4,myblock)*yshell0(i,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(yshell0(i,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(yshell0(i,k,6,myblock)+yshell0(i,k,5,myblock)+yshell0(i,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=yshell0(i,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*yshell0(i,k,5,myblock) &
	 -cssq*(yshell0(i,k,6,myblock)+yshell0(i,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*yshell0(i,k,5,myblock) &
	 -cssq*(yshell0(i,k,6,myblock)+yshell0(i,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=yshell0(i,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*yshell0(i,k,6,myblock) &
	 -cssq*(yshell0(i,k,5,myblock)+yshell0(i,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*yshell0(i,k,6,myblock) &
	 -cssq*(yshell0(i,k,5,myblock)+yshell0(i,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=yshell0(i,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*yshell0(i,k,7,myblock) &
	 -cssq*(yshell0(i,k,5,myblock)+yshell0(i,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*yshell0(i,k,7,myblock) &
	 -cssq*(yshell0(i,k,5,myblock)+yshell0(i,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(yshell0(i,k,2,myblock)+yshell0(i,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qyy*yshell0(i,k,6,myblock) &
	 -cssq*yshell0(i,k,7,myblock)+two*qxy_7_8*yshell0(i,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qyy*yshell0(i,k,6,myblock) &
	 -cssq*yshell0(i,k,7,myblock)+two*qxy_7_8*yshell0(i,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-yshell0(i,k,2,myblock)+yshell0(i,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qyy*yshell0(i,k,6,myblock) &
	 -cssq*yshell0(i,k,7,myblock)+two*qxy_9_10*yshell0(i,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qyy*yshell0(i,k,6,myblock) &
	 -cssq*yshell0(i,k,7,myblock)+two*qxy_9_10*yshell0(i,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(yshell0(i,k,2,myblock)+yshell0(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qzz*yshell0(i,k,7,myblock) &
	 -cssq*yshell0(i,k,6,myblock)+two*qxz_15_16*yshell0(i,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qzz*yshell0(i,k,7,myblock) &
	 -cssq*yshell0(i,k,6,myblock)+two*qxz_15_16*yshell0(i,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-yshell0(i,k,2,myblock)+yshell0(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qzz*yshell0(i,k,7,myblock) &
	 -cssq*yshell0(i,k,6,myblock)+two*qxz_17_18*yshell0(i,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qzz*yshell0(i,k,7,myblock) &
	 -cssq*yshell0(i,k,6,myblock)+two*qxz_17_18*yshell0(i,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(yshell0(i,k,3,myblock)+yshell0(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell0(i,k,6,myblock)+qzz*yshell0(i,k,7,myblock) &
	 -cssq*yshell0(i,k,5,myblock)+two*qyz_11_12*yshell0(i,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell0(i,k,6,myblock)+qzz*yshell0(i,k,7,myblock) &
	 -cssq*yshell0(i,k,5,myblock)+two*qyz_11_12*yshell0(i,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(yshell0(i,k,3,myblock)-yshell0(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell0(i,k,6,myblock)+qzz*yshell0(i,k,7,myblock) &
	 -cssq*yshell0(i,k,5,myblock)+two*qyz_13_14*yshell0(i,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell0(i,k,6,myblock)+qzz*yshell0(i,k,7,myblock) &
	 -cssq*yshell0(i,k,5,myblock)+two*qyz_13_14*yshell0(i,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
    endif

    if(lk==0  .and. li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d)then
      
      uu=halfonecssq*(zshell1(i,j,2,myblock)*zshell1(i,j,2,myblock) &
     + zshell1(i,j,3,myblock)*zshell1(i,j,3,myblock) + zshell1(i,j,4,myblock)*zshell1(i,j,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(zshell1(i,j,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(zshell1(i,j,6,myblock)+zshell1(i,j,5,myblock)+zshell1(i,j,7,myblock)))
	
    
	!1 -1  0  0
	udotc=zshell1(i,j,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*zshell1(i,j,5,myblock) &
	 -cssq*(zshell1(i,j,6,myblock)+zshell1(i,j,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*zshell1(i,j,5,myblock) &
	 -cssq*(zshell1(i,j,6,myblock)+zshell1(i,j,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=zshell1(i,j,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*zshell1(i,j,6,myblock) &
	 -cssq*(zshell1(i,j,5,myblock)+zshell1(i,j,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*zshell1(i,j,6,myblock) &
	 -cssq*(zshell1(i,j,5,myblock)+zshell1(i,j,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=zshell1(i,j,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*zshell1(i,j,7,myblock) &
	 -cssq*(zshell1(i,j,5,myblock)+zshell1(i,j,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*zshell1(i,j,7,myblock) &
	 -cssq*(zshell1(i,j,5,myblock)+zshell1(i,j,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(zshell1(i,j,2,myblock)+zshell1(i,j,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qyy*zshell1(i,j,6,myblock) &
	 -cssq*zshell1(i,j,7,myblock)+two*qxy_7_8*zshell1(i,j,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qyy*zshell1(i,j,6,myblock) &
	 -cssq*zshell1(i,j,7,myblock)+two*qxy_7_8*zshell1(i,j,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-zshell1(i,j,2,myblock)+zshell1(i,j,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qyy*zshell1(i,j,6,myblock) &
	 -cssq*zshell1(i,j,7,myblock)+two*qxy_9_10*zshell1(i,j,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qyy*zshell1(i,j,6,myblock) &
	 -cssq*zshell1(i,j,7,myblock)+two*qxy_9_10*zshell1(i,j,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(zshell1(i,j,2,myblock)+zshell1(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qzz*zshell1(i,j,7,myblock) &
	 -cssq*zshell1(i,j,6,myblock)+two*qxz_15_16*zshell1(i,j,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qzz*zshell1(i,j,7,myblock) &
	 -cssq*zshell1(i,j,6,myblock)+two*qxz_15_16*zshell1(i,j,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-zshell1(i,j,2,myblock)+zshell1(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qzz*zshell1(i,j,7,myblock) &
	 -cssq*zshell1(i,j,6,myblock)+two*qxz_17_18*zshell1(i,j,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qzz*zshell1(i,j,7,myblock) &
	 -cssq*zshell1(i,j,6,myblock)+two*qxz_17_18*zshell1(i,j,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(zshell1(i,j,3,myblock)+zshell1(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell1(i,j,6,myblock)+qzz*zshell1(i,j,7,myblock) &
	 -cssq*zshell1(i,j,5,myblock)+two*qyz_11_12*zshell1(i,j,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell1(i,j,6,myblock)+qzz*zshell1(i,j,7,myblock) &
	 -cssq*zshell1(i,j,5,myblock)+two*qyz_11_12*zshell1(i,j,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(zshell1(i,j,3,myblock)-zshell1(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell1(i,j,6,myblock)+qzz*zshell1(i,j,7,myblock) &
	 -cssq*zshell1(i,j,5,myblock)+two*qyz_13_14*zshell1(i,j,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell1(i,j,6,myblock)+qzz*zshell1(i,j,7,myblock) &
	 -cssq*zshell1(i,j,5,myblock)+two*qyz_13_14*zshell1(i,j,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
    endif
    
    if(lk==TILE_DIMz_d+1 .and. li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d)then
      
      uu=halfonecssq*(zshell0(i,j,2,myblock)*zshell0(i,j,2,myblock) &
     + zshell0(i,j,3,myblock)*zshell0(i,j,3,myblock) + zshell0(i,j,4,myblock)*zshell0(i,j,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(zshell0(i,j,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(zshell0(i,j,6,myblock)+zshell0(i,j,5,myblock)+zshell0(i,j,7,myblock)))
	
    
	!1 -1  0  0
	udotc=zshell0(i,j,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*zshell0(i,j,5,myblock) &
	 -cssq*(zshell0(i,j,6,myblock)+zshell0(i,j,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*zshell0(i,j,5,myblock) &
	 -cssq*(zshell0(i,j,6,myblock)+zshell0(i,j,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=zshell0(i,j,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*zshell0(i,j,6,myblock) &
	 -cssq*(zshell0(i,j,5,myblock)+zshell0(i,j,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*zshell0(i,j,6,myblock) &
	 -cssq*(zshell0(i,j,5,myblock)+zshell0(i,j,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=zshell0(i,j,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*zshell0(i,j,7,myblock) &
	 -cssq*(zshell0(i,j,5,myblock)+zshell0(i,j,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*zshell0(i,j,7,myblock) &
	 -cssq*(zshell0(i,j,5,myblock)+zshell0(i,j,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(zshell0(i,j,2,myblock)+zshell0(i,j,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qyy*zshell0(i,j,6,myblock) &
	 -cssq*zshell0(i,j,7,myblock)+two*qxy_7_8*zshell0(i,j,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qyy*zshell0(i,j,6,myblock) &
	 -cssq*zshell0(i,j,7,myblock)+two*qxy_7_8*zshell0(i,j,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-zshell0(i,j,2,myblock)+zshell0(i,j,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qyy*zshell0(i,j,6,myblock) &
	 -cssq*zshell0(i,j,7,myblock)+two*qxy_9_10*zshell0(i,j,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qyy*zshell0(i,j,6,myblock) &
	 -cssq*zshell0(i,j,7,myblock)+two*qxy_9_10*zshell0(i,j,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(zshell0(i,j,2,myblock)+zshell0(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qzz*zshell0(i,j,7,myblock) &
	 -cssq*zshell0(i,j,6,myblock)+two*qxz_15_16*zshell0(i,j,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qzz*zshell0(i,j,7,myblock) &
	 -cssq*zshell0(i,j,6,myblock)+two*qxz_15_16*zshell0(i,j,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-zshell0(i,j,2,myblock)+zshell0(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qzz*zshell0(i,j,7,myblock) &
	 -cssq*zshell0(i,j,6,myblock)+two*qxz_17_18*zshell0(i,j,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qzz*zshell0(i,j,7,myblock) &
	 -cssq*zshell0(i,j,6,myblock)+two*qxz_17_18*zshell0(i,j,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(zshell0(i,j,3,myblock)+zshell0(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell0(i,j,6,myblock)+qzz*zshell0(i,j,7,myblock) &
	 -cssq*zshell0(i,j,5,myblock)+two*qyz_11_12*zshell0(i,j,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell0(i,j,6,myblock)+qzz*zshell0(i,j,7,myblock) &
	 -cssq*zshell0(i,j,5,myblock)+two*qyz_11_12*zshell0(i,j,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(zshell0(i,j,3,myblock)-zshell0(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell0(i,j,6,myblock)+qzz*zshell0(i,j,7,myblock) &
	 -cssq*zshell0(i,j,5,myblock)+two*qyz_13_14*zshell0(i,j,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell0(i,j,6,myblock)+qzz*zshell0(i,j,7,myblock) &
	 -cssq*zshell0(i,j,5,myblock)+two*qyz_13_14*zshell0(i,j,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
	
    endif  
    
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
	hfields(i,j,k,1,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,2,myblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,3,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,4,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,5,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(i,j,k,6,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,7,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(i,j,k,8,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,9,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,10,myblock)=udotc
     
#ifdef PRESSCORR

	!call syncthreads
	
	uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
	+ hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,5,myblock)=hfields(i,j,k,5,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(i,j,k,6,myblock)=hfields(i,j,k,6,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,7,myblock)=hfields(i,j,k,7,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(i,j,k,8,myblock)=hfields(i,j,k,8,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(i,j,k,9,myblock)=hfields(i,j,k,9,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(i,j,k,10,myblock)=hfields(i,j,k,10,myblock)-udotc
	
#endif	
    
    ! Halo Faces
    if(li==1)then
      
      xshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(li==TILE_DIMx_d)then
      
      xshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
	  
    endif

    if(lj==1)then
      
      yshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lj==TILE_DIMy_d)then
      
      yshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    if(lk==1)then
      
      zshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lk==TILE_DIMz_d)then
      
      zshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    return
    
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
    
    !inner cube
    if(li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d .and. &
     lk>=1 .and. lk<=TILE_DIMz_d)then
     
    uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
     + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(hfields(i,j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
	
    endif
    
    ! Halo Faces
    if(li==0)then
      
    uu=halfonecssq*(xshell1h(j,k,2,myblock)*xshell1h(j,k,2,myblock) &
     + xshell1h(j,k,3,myblock)*xshell1h(j,k,3,myblock) + xshell1h(j,k,4,myblock)*xshell1h(j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(xshell1h(j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(xshell1h(j,k,6,myblock)+xshell1h(j,k,5,myblock)+xshell1h(j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=xshell1h(j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*xshell1h(j,k,5,myblock) &
	 -cssq*(xshell1h(j,k,6,myblock)+xshell1h(j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*xshell1h(j,k,5,myblock) &
	 -cssq*(xshell1h(j,k,6,myblock)+xshell1h(j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=xshell1h(j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*xshell1h(j,k,6,myblock) &
	 -cssq*(xshell1h(j,k,5,myblock)+xshell1h(j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*xshell1h(j,k,6,myblock) &
	 -cssq*(xshell1h(j,k,5,myblock)+xshell1h(j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=xshell1h(j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*xshell1h(j,k,7,myblock) &
	 -cssq*(xshell1h(j,k,5,myblock)+xshell1h(j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*xshell1h(j,k,7,myblock) &
	 -cssq*(xshell1h(j,k,5,myblock)+xshell1h(j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(xshell1h(j,k,2,myblock)+xshell1h(j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qyy*xshell1h(j,k,6,myblock) &
	 -cssq*xshell1h(j,k,7,myblock)+two*qxy_7_8*xshell1h(j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qyy*xshell1h(j,k,6,myblock) &
	 -cssq*xshell1h(j,k,7,myblock)+two*qxy_7_8*xshell1h(j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-xshell1h(j,k,2,myblock)+xshell1h(j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qyy*xshell1h(j,k,6,myblock) &
	 -cssq*xshell1h(j,k,7,myblock)+two*qxy_9_10*xshell1h(j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qyy*xshell1h(j,k,6,myblock) &
	 -cssq*xshell1h(j,k,7,myblock)+two*qxy_9_10*xshell1h(j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(xshell1h(j,k,2,myblock)+xshell1h(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qzz*xshell1h(j,k,7,myblock) &
	 -cssq*xshell1h(j,k,6,myblock)+two*qxz_15_16*xshell1h(j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qzz*xshell1h(j,k,7,myblock) &
	 -cssq*xshell1h(j,k,6,myblock)+two*qxz_15_16*xshell1h(j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-xshell1h(j,k,2,myblock)+xshell1h(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qzz*xshell1h(j,k,7,myblock) &
	 -cssq*xshell1h(j,k,6,myblock)+two*qxz_17_18*xshell1h(j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qzz*xshell1h(j,k,7,myblock) &
	 -cssq*xshell1h(j,k,6,myblock)+two*qxz_17_18*xshell1h(j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(xshell1h(j,k,3,myblock)+xshell1h(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell1h(j,k,6,myblock)+qzz*xshell1h(j,k,7,myblock) &
	 -cssq*xshell1h(j,k,5,myblock)+two*qyz_11_12*xshell1h(j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell1h(j,k,6,myblock)+qzz*xshell1h(j,k,7,myblock) &
	 -cssq*xshell1h(j,k,5,myblock)+two*qyz_11_12*xshell1h(j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(xshell1h(j,k,3,myblock)-xshell1h(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell1h(j,k,6,myblock)+qzz*xshell1h(j,k,7,myblock) &
	 -cssq*xshell1h(j,k,5,myblock)+two*qyz_13_14*xshell1h(j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell1h(j,k,6,myblock)+qzz*xshell1h(j,k,7,myblock) &
	 -cssq*xshell1h(j,k,5,myblock)+two*qyz_13_14*xshell1h(j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
    endif
    if(li==TILE_DIMx_d+1)then
      
    uu=halfonecssq*(xshell0h(j,k,2,myblock)*xshell0h(j,k,2,myblock) &
     + xshell0h(j,k,3,myblock)*xshell0h(j,k,3,myblock) + xshell0h(j,k,4,myblock)*xshell0h(j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(xshell0h(j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(xshell0h(j,k,6,myblock)+xshell0h(j,k,5,myblock)+xshell0h(j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=xshell0h(j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*xshell0h(j,k,5,myblock) &
	 -cssq*(xshell0h(j,k,6,myblock)+xshell0h(j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*xshell0h(j,k,5,myblock) &
	 -cssq*(xshell0h(j,k,6,myblock)+xshell0h(j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=xshell0h(j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*xshell0h(j,k,6,myblock) &
	 -cssq*(xshell0h(j,k,5,myblock)+xshell0h(j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*xshell0h(j,k,6,myblock) &
	 -cssq*(xshell0h(j,k,5,myblock)+xshell0h(j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=xshell0h(j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*xshell0h(j,k,7,myblock) &
	 -cssq*(xshell0h(j,k,5,myblock)+xshell0h(j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*xshell0h(j,k,7,myblock) &
	 -cssq*(xshell0h(j,k,5,myblock)+xshell0h(j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(xshell0h(j,k,2,myblock)+xshell0h(j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qyy*xshell0h(j,k,6,myblock) &
	 -cssq*xshell0h(j,k,7,myblock)+two*qxy_7_8*xshell0h(j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qyy*xshell0h(j,k,6,myblock) &
	 -cssq*xshell0h(j,k,7,myblock)+two*qxy_7_8*xshell0h(j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-xshell0h(j,k,2,myblock)+xshell0h(j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qyy*xshell0h(j,k,6,myblock) &
	 -cssq*xshell0h(j,k,7,myblock)+two*qxy_9_10*xshell0h(j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qyy*xshell0h(j,k,6,myblock) &
	 -cssq*xshell0h(j,k,7,myblock)+two*qxy_9_10*xshell0h(j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(xshell0h(j,k,2,myblock)+xshell0h(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qzz*xshell0h(j,k,7,myblock) &
	 -cssq*xshell0h(j,k,6,myblock)+two*qxz_15_16*xshell0h(j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qzz*xshell0h(j,k,7,myblock) &
	 -cssq*xshell0h(j,k,6,myblock)+two*qxz_15_16*xshell0h(j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-xshell0h(j,k,2,myblock)+xshell0h(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qzz*xshell0h(j,k,7,myblock) &
	 -cssq*xshell0h(j,k,6,myblock)+two*qxz_17_18*xshell0h(j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qzz*xshell0h(j,k,7,myblock) &
	 -cssq*xshell0h(j,k,6,myblock)+two*qxz_17_18*xshell0h(j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(xshell0h(j,k,3,myblock)+xshell0h(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell0h(j,k,6,myblock)+qzz*xshell0h(j,k,7,myblock) &
	 -cssq*xshell0h(j,k,5,myblock)+two*qyz_11_12*xshell0h(j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell0h(j,k,6,myblock)+qzz*xshell0h(j,k,7,myblock) &
	 -cssq*xshell0h(j,k,5,myblock)+two*qyz_11_12*xshell0h(j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(xshell0h(j,k,3,myblock)-xshell0h(j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell0h(j,k,6,myblock)+qzz*xshell0h(j,k,7,myblock) &
	 -cssq*xshell0h(j,k,5,myblock)+two*qyz_13_14*xshell0h(j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*xshell0h(j,k,6,myblock)+qzz*xshell0h(j,k,7,myblock) &
	 -cssq*xshell0h(j,k,5,myblock)+two*qyz_13_14*xshell0h(j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
	  
    endif

    if(lj==0 .and. li>=1 .and. li<=TILE_DIMx_d)then
      
    uu=halfonecssq*(yshell1h(i,k,2,myblock)*yshell1h(i,k,2,myblock) &
     + yshell1h(i,k,3,myblock)*yshell1h(i,k,3,myblock) + yshell1h(i,k,4,myblock)*yshell1h(i,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(yshell1h(i,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(yshell1h(i,k,6,myblock)+yshell1h(i,k,5,myblock)+yshell1h(i,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=yshell1h(i,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*yshell1h(i,k,5,myblock) &
	 -cssq*(yshell1h(i,k,6,myblock)+yshell1h(i,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*yshell1h(i,k,5,myblock) &
	 -cssq*(yshell1h(i,k,6,myblock)+yshell1h(i,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=yshell1h(i,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*yshell1h(i,k,6,myblock) &
	 -cssq*(yshell1h(i,k,5,myblock)+yshell1h(i,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*yshell1h(i,k,6,myblock) &
	 -cssq*(yshell1h(i,k,5,myblock)+yshell1h(i,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=yshell1h(i,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*yshell1h(i,k,7,myblock) &
	 -cssq*(yshell1h(i,k,5,myblock)+yshell1h(i,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*yshell1h(i,k,7,myblock) &
	 -cssq*(yshell1h(i,k,5,myblock)+yshell1h(i,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(yshell1h(i,k,2,myblock)+yshell1h(i,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qyy*yshell1h(i,k,6,myblock) &
	 -cssq*yshell1h(i,k,7,myblock)+two*qxy_7_8*yshell1h(i,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qyy*yshell1h(i,k,6,myblock) &
	 -cssq*yshell1h(i,k,7,myblock)+two*qxy_7_8*yshell1h(i,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-yshell1h(i,k,2,myblock)+yshell1h(i,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qyy*yshell1h(i,k,6,myblock) &
	 -cssq*yshell1h(i,k,7,myblock)+two*qxy_9_10*yshell1h(i,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qyy*yshell1h(i,k,6,myblock) &
	 -cssq*yshell1h(i,k,7,myblock)+two*qxy_9_10*yshell1h(i,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(yshell1h(i,k,2,myblock)+yshell1h(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qzz*yshell1h(i,k,7,myblock) &
	 -cssq*yshell1h(i,k,6,myblock)+two*qxz_15_16*yshell1h(i,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qzz*yshell1h(i,k,7,myblock) &
	 -cssq*yshell1h(i,k,6,myblock)+two*qxz_15_16*yshell1h(i,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-yshell1h(i,k,2,myblock)+yshell1h(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qzz*yshell1h(i,k,7,myblock) &
	 -cssq*yshell1h(i,k,6,myblock)+two*qxz_17_18*yshell1h(i,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qzz*yshell1h(i,k,7,myblock) &
	 -cssq*yshell1h(i,k,6,myblock)+two*qxz_17_18*yshell1h(i,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(yshell1h(i,k,3,myblock)+yshell1h(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell1h(i,k,6,myblock)+qzz*yshell1h(i,k,7,myblock) &
	 -cssq*yshell1h(i,k,5,myblock)+two*qyz_11_12*yshell1h(i,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell1h(i,k,6,myblock)+qzz*yshell1h(i,k,7,myblock) &
	 -cssq*yshell1h(i,k,5,myblock)+two*qyz_11_12*yshell1h(i,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(yshell1h(i,k,3,myblock)-yshell1h(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell1h(i,k,6,myblock)+qzz*yshell1h(i,k,7,myblock) &
	 -cssq*yshell1h(i,k,5,myblock)+two*qyz_13_14*yshell1h(i,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell1h(i,k,6,myblock)+qzz*yshell1h(i,k,7,myblock) &
	 -cssq*yshell1h(i,k,5,myblock)+two*qyz_13_14*yshell1h(i,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
    endif
    
    if(lj==TILE_DIMy_d+1  .and. li>=1 .and. li<=TILE_DIMx_d)then
      
    uu=halfonecssq*(yshell0h(i,k,2,myblock)*yshell0h(i,k,2,myblock) &
     + yshell0h(i,k,3,myblock)*yshell0h(i,k,3,myblock) + yshell0h(i,k,4,myblock)*yshell0h(i,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(yshell0h(i,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(yshell0h(i,k,6,myblock)+yshell0h(i,k,5,myblock)+yshell0h(i,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=yshell0h(i,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*yshell0h(i,k,5,myblock) &
	 -cssq*(yshell0h(i,k,6,myblock)+yshell0h(i,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*yshell0h(i,k,5,myblock) &
	 -cssq*(yshell0h(i,k,6,myblock)+yshell0h(i,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=yshell0h(i,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*yshell0h(i,k,6,myblock) &
	 -cssq*(yshell0h(i,k,5,myblock)+yshell0h(i,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*yshell0h(i,k,6,myblock) &
	 -cssq*(yshell0h(i,k,5,myblock)+yshell0h(i,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=yshell0h(i,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*yshell0h(i,k,7,myblock) &
	 -cssq*(yshell0h(i,k,5,myblock)+yshell0h(i,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*yshell0h(i,k,7,myblock) &
	 -cssq*(yshell0h(i,k,5,myblock)+yshell0h(i,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(yshell0h(i,k,2,myblock)+yshell0h(i,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qyy*yshell0h(i,k,6,myblock) &
	 -cssq*yshell0h(i,k,7,myblock)+two*qxy_7_8*yshell0h(i,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qyy*yshell0h(i,k,6,myblock) &
	 -cssq*yshell0h(i,k,7,myblock)+two*qxy_7_8*yshell0h(i,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-yshell0h(i,k,2,myblock)+yshell0h(i,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qyy*yshell0h(i,k,6,myblock) &
	 -cssq*yshell0h(i,k,7,myblock)+two*qxy_9_10*yshell0h(i,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qyy*yshell0h(i,k,6,myblock) &
	 -cssq*yshell0h(i,k,7,myblock)+two*qxy_9_10*yshell0h(i,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(yshell0h(i,k,2,myblock)+yshell0h(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qzz*yshell0h(i,k,7,myblock) &
	 -cssq*yshell0h(i,k,6,myblock)+two*qxz_15_16*yshell0h(i,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qzz*yshell0h(i,k,7,myblock) &
	 -cssq*yshell0h(i,k,6,myblock)+two*qxz_15_16*yshell0h(i,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-yshell0h(i,k,2,myblock)+yshell0h(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qzz*yshell0h(i,k,7,myblock) &
	 -cssq*yshell0h(i,k,6,myblock)+two*qxz_17_18*yshell0h(i,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qzz*yshell0h(i,k,7,myblock) &
	 -cssq*yshell0h(i,k,6,myblock)+two*qxz_17_18*yshell0h(i,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(yshell0h(i,k,3,myblock)+yshell0h(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell0h(i,k,6,myblock)+qzz*yshell0h(i,k,7,myblock) &
	 -cssq*yshell0h(i,k,5,myblock)+two*qyz_11_12*yshell0h(i,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell0h(i,k,6,myblock)+qzz*yshell0h(i,k,7,myblock) &
	 -cssq*yshell0h(i,k,5,myblock)+two*qyz_11_12*yshell0h(i,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(yshell0h(i,k,3,myblock)-yshell0h(i,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell0h(i,k,6,myblock)+qzz*yshell0h(i,k,7,myblock) &
	 -cssq*yshell0h(i,k,5,myblock)+two*qyz_13_14*yshell0h(i,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*yshell0h(i,k,6,myblock)+qzz*yshell0h(i,k,7,myblock) &
	 -cssq*yshell0h(i,k,5,myblock)+two*qyz_13_14*yshell0h(i,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
    endif

    if(lk==0  .and. li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d)then
      
      uu=halfonecssq*(zshell1h(i,j,2,myblock)*zshell1h(i,j,2,myblock) &
     + zshell1h(i,j,3,myblock)*zshell1h(i,j,3,myblock) + zshell1h(i,j,4,myblock)*zshell1h(i,j,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(zshell1h(i,j,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(zshell1h(i,j,6,myblock)+zshell1h(i,j,5,myblock)+zshell1h(i,j,7,myblock)))
	
    
	!1 -1  0  0
	udotc=zshell1h(i,j,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*zshell1h(i,j,5,myblock) &
	 -cssq*(zshell1h(i,j,6,myblock)+zshell1h(i,j,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*zshell1h(i,j,5,myblock) &
	 -cssq*(zshell1h(i,j,6,myblock)+zshell1h(i,j,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=zshell1h(i,j,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*zshell1h(i,j,6,myblock) &
	 -cssq*(zshell1h(i,j,5,myblock)+zshell1h(i,j,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*zshell1h(i,j,6,myblock) &
	 -cssq*(zshell1h(i,j,5,myblock)+zshell1h(i,j,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=zshell1h(i,j,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*zshell1h(i,j,7,myblock) &
	 -cssq*(zshell1h(i,j,5,myblock)+zshell1h(i,j,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*zshell1h(i,j,7,myblock) &
	 -cssq*(zshell1h(i,j,5,myblock)+zshell1h(i,j,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(zshell1h(i,j,2,myblock)+zshell1h(i,j,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qyy*zshell1h(i,j,6,myblock) &
	 -cssq*zshell1h(i,j,7,myblock)+two*qxy_7_8*zshell1h(i,j,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qyy*zshell1h(i,j,6,myblock) &
	 -cssq*zshell1h(i,j,7,myblock)+two*qxy_7_8*zshell1h(i,j,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-zshell1h(i,j,2,myblock)+zshell1h(i,j,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qyy*zshell1h(i,j,6,myblock) &
	 -cssq*zshell1h(i,j,7,myblock)+two*qxy_9_10*zshell1h(i,j,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qyy*zshell1h(i,j,6,myblock) &
	 -cssq*zshell1h(i,j,7,myblock)+two*qxy_9_10*zshell1h(i,j,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(zshell1h(i,j,2,myblock)+zshell1h(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qzz*zshell1h(i,j,7,myblock) &
	 -cssq*zshell1h(i,j,6,myblock)+two*qxz_15_16*zshell1h(i,j,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qzz*zshell1h(i,j,7,myblock) &
	 -cssq*zshell1h(i,j,6,myblock)+two*qxz_15_16*zshell1h(i,j,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-zshell1h(i,j,2,myblock)+zshell1h(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qzz*zshell1h(i,j,7,myblock) &
	 -cssq*zshell1h(i,j,6,myblock)+two*qxz_17_18*zshell1h(i,j,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qzz*zshell1h(i,j,7,myblock) &
	 -cssq*zshell1h(i,j,6,myblock)+two*qxz_17_18*zshell1h(i,j,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(zshell1h(i,j,3,myblock)+zshell1h(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell1h(i,j,6,myblock)+qzz*zshell1h(i,j,7,myblock) &
	 -cssq*zshell1h(i,j,5,myblock)+two*qyz_11_12*zshell1h(i,j,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell1h(i,j,6,myblock)+qzz*zshell1h(i,j,7,myblock) &
	 -cssq*zshell1h(i,j,5,myblock)+two*qyz_11_12*zshell1h(i,j,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(zshell1h(i,j,3,myblock)-zshell1h(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell1h(i,j,6,myblock)+qzz*zshell1h(i,j,7,myblock) &
	 -cssq*zshell1h(i,j,5,myblock)+two*qyz_13_14*zshell1h(i,j,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell1h(i,j,6,myblock)+qzz*zshell1h(i,j,7,myblock) &
	 -cssq*zshell1h(i,j,5,myblock)+two*qyz_13_14*zshell1h(i,j,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
      
    endif
    
    if(lk==TILE_DIMz_d+1 .and. li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d)then
      
      uu=halfonecssq*(zshell0h(i,j,2,myblock)*zshell0h(i,j,2,myblock) &
     + zshell0h(i,j,3,myblock)*zshell0h(i,j,3,myblock) + zshell0h(i,j,4,myblock)*zshell0h(i,j,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(zshell0h(i,j,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(zshell0h(i,j,6,myblock)+zshell0h(i,j,5,myblock)+zshell0h(i,j,7,myblock)))
	
    
	!1 -1  0  0
	udotc=zshell0h(i,j,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*zshell0h(i,j,5,myblock) &
	 -cssq*(zshell0h(i,j,6,myblock)+zshell0h(i,j,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*zshell0h(i,j,5,myblock) &
	 -cssq*(zshell0h(i,j,6,myblock)+zshell0h(i,j,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=zshell0h(i,j,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*zshell0h(i,j,6,myblock) &
	 -cssq*(zshell0h(i,j,5,myblock)+zshell0h(i,j,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*zshell0h(i,j,6,myblock) &
	 -cssq*(zshell0h(i,j,5,myblock)+zshell0h(i,j,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=zshell0h(i,j,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*zshell0h(i,j,7,myblock) &
	 -cssq*(zshell0h(i,j,5,myblock)+zshell0h(i,j,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*zshell0h(i,j,7,myblock) &
	 -cssq*(zshell0h(i,j,5,myblock)+zshell0h(i,j,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(zshell0h(i,j,2,myblock)+zshell0h(i,j,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qyy*zshell0h(i,j,6,myblock) &
	 -cssq*zshell0h(i,j,7,myblock)+two*qxy_7_8*zshell0h(i,j,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qyy*zshell0h(i,j,6,myblock) &
	 -cssq*zshell0h(i,j,7,myblock)+two*qxy_7_8*zshell0h(i,j,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-zshell0h(i,j,2,myblock)+zshell0h(i,j,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qyy*zshell0h(i,j,6,myblock) &
	 -cssq*zshell0h(i,j,7,myblock)+two*qxy_9_10*zshell0h(i,j,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qyy*zshell0h(i,j,6,myblock) &
	 -cssq*zshell0h(i,j,7,myblock)+two*qxy_9_10*zshell0h(i,j,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(zshell0h(i,j,2,myblock)+zshell0h(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qzz*zshell0h(i,j,7,myblock) &
	 -cssq*zshell0h(i,j,6,myblock)+two*qxz_15_16*zshell0h(i,j,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qzz*zshell0h(i,j,7,myblock) &
	 -cssq*zshell0h(i,j,6,myblock)+two*qxz_15_16*zshell0h(i,j,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-zshell0h(i,j,2,myblock)+zshell0h(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qzz*zshell0h(i,j,7,myblock) &
	 -cssq*zshell0h(i,j,6,myblock)+two*qxz_17_18*zshell0h(i,j,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qzz*zshell0h(i,j,7,myblock) &
	 -cssq*zshell0h(i,j,6,myblock)+two*qxz_17_18*zshell0h(i,j,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(zshell0h(i,j,3,myblock)+zshell0h(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell0h(i,j,6,myblock)+qzz*zshell0h(i,j,7,myblock) &
	 -cssq*zshell0h(i,j,5,myblock)+two*qyz_11_12*zshell0h(i,j,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell0h(i,j,6,myblock)+qzz*zshell0h(i,j,7,myblock) &
	 -cssq*zshell0h(i,j,5,myblock)+two*qyz_11_12*zshell0h(i,j,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(zshell0h(i,j,3,myblock)-zshell0h(i,j,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell0h(i,j,6,myblock)+qzz*zshell0h(i,j,7,myblock) &
	 -cssq*zshell0h(i,j,5,myblock)+two*qyz_13_14*zshell0h(i,j,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*zshell0h(i,j,6,myblock)+qzz*zshell0h(i,j,7,myblock) &
	 -cssq*zshell0h(i,j,5,myblock)+two*qyz_13_14*zshell0h(i,j,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
	
    endif  
    
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
	hfields(i,j,k,1,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,2,myblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,3,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,4,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,5,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(i,j,k,6,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,7,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(i,j,k,8,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,9,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,10,myblock)=udotc
     
#ifdef PRESSCORR

	!call syncthreads
	
	uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
	+ hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,5,myblock)=hfields(i,j,k,5,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(i,j,k,6,myblock)=hfields(i,j,k,6,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,7,myblock)=hfields(i,j,k,7,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(i,j,k,8,myblock)=hfields(i,j,k,8,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(i,j,k,9,myblock)=hfields(i,j,k,9,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(i,j,k,10,myblock)=hfields(i,j,k,10,myblock)-udotc
	
#endif	
    
    ! Halo Faces
    if(li==1)then
      
      xshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(li==TILE_DIMx_d)then
      
      xshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
	  
    endif

    if(lj==1)then
      
      yshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lj==TILE_DIMy_d)then
      
      yshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    if(lk==1)then
      
      zshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lk==TILE_DIMz_d)then
      
      zshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    return
    
 
  end subroutine streamcoll_shared_halo_flop
  
  attributes(global) subroutine streamcoll_minshared_halo()
	
	implicit none  
	  
    integer :: i,j,k,gi,gj,gk,myblock
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
    
    !inner cube
    if(li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d .and. &
     lk>=1 .and. lk<=TILE_DIMz_d)then
     
    uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
     + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(hfields(i,j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
	
    endif
    
    ! Halo Faces
    if(li==0 .and. lj>=1 .and. lj<=TILE_DIMy_d .and. lk>=1 .and. lk<=TILE_DIMz_d)then
      
      
      uu=halfonecssq*(xshell1(j,k,2,myblock)*xshell1(j,k,2,myblock) + xshell1(j,k,3,myblock)*xshell1(j,k,3,myblock) + xshell1(j,k,4,myblock)*xshell1(j,k,4,myblock))
      
      !1 -1  0  0
	  udotc=xshell1(j,k,2,myblock)*onecssq
	  f01(li,lj,lk)=p1*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*xshell1(j,k,5,myblock)-cssq*(xshell1(j,k,6,myblock)+xshell1(j,k,7,myblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(xshell1(j,k,2,myblock)+xshell1(j,k,3,myblock))*onecssq
	  f07(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qyy*xshell1(j,k,6,myblock)-cssq*xshell1(j,k,7,myblock)+two*qxy_7_8*xshell1(j,k,8,myblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-xshell1(j,k,2,myblock)+xshell1(j,k,3,myblock))*onecssq
	  f09(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qyy*xshell1(j,k,6,myblock)-cssq*xshell1(j,k,7,myblock)+two*qxy_9_10*xshell1(j,k,8,myblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(xshell1(j,k,2,myblock)+xshell1(j,k,4,myblock))*onecssq
	  f15(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qzz*xshell1(j,k,7,myblock)-cssq*xshell1(j,k,6,myblock)+two*qxz_15_16*xshell1(j,k,9,myblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-xshell1(j,k,2,myblock)+xshell1(j,k,4,myblock))*onecssq
	  f18(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qzz*xshell1(j,k,7,myblock)-cssq*xshell1(j,k,6,myblock)+two*qxz_17_18*xshell1(j,k,9,myblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    
    if(li==TILE_DIMx_d+1 .and. lj>=1 .and. lj<=TILE_DIMy_d .and. lk>=1 .and. lk<=TILE_DIMz_d)then
      
      
      uu=halfonecssq*(xshell0(j,k,2,myblock)*xshell0(j,k,2,myblock) + xshell0(j,k,3,myblock)*xshell0(j,k,3,myblock) + xshell0(j,k,4,myblock)*xshell0(j,k,4,myblock))
      
      !2 +1  0  0
	  udotc=xshell0(j,k,2,myblock)*onecssq
	  f02(li,lj,lk)=p1*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*xshell0(j,k,5,myblock)-cssq*(xshell0(j,k,6,myblock)+xshell0(j,k,7,myblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(xshell0(j,k,2,myblock)+xshell0(j,k,3,myblock))*onecssq
	  f08(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qyy*xshell0(j,k,6,myblock)-cssq*xshell0(j,k,7,myblock)+two*qxy_7_8*xshell0(j,k,8,myblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-xshell0(j,k,2,myblock)+xshell0(j,k,3,myblock))*onecssq
	  f10(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qyy*xshell0(j,k,6,myblock)-cssq*xshell0(j,k,7,myblock)+two*qxy_9_10*xshell0(j,k,8,myblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(xshell0(j,k,2,myblock)+xshell0(j,k,4,myblock))*onecssq
	  f16(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qzz*xshell0(j,k,7,myblock)-cssq*xshell0(j,k,6,myblock)+two*qxz_15_16*xshell0(j,k,9,myblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-xshell0(j,k,2,myblock)+xshell0(j,k,4,myblock))*onecssq
	  f17(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qzz*xshell0(j,k,7,myblock)-cssq*xshell0(j,k,6,myblock)+two*qxz_17_18*xshell0(j,k,9,myblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==0 .and. li>=1 .and. li<=TILE_DIMx_d .and. lk>=1 .and. lk<=TILE_DIMz_d)then
      
      
      uu=halfonecssq*(yshell1(i,k,2,myblock)*yshell1(i,k,2,myblock) + yshell1(i,k,3,myblock)*yshell1(i,k,3,myblock) + yshell1(i,k,4,myblock)*yshell1(i,k,4,myblock))
      
      !3 0 -1  0
	  udotc=yshell1(i,k,3,myblock)*onecssq
	  f03(li,lj,lk)=p1*(yshell1(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*yshell1(i,k,6,myblock)-cssq*(yshell1(i,k,5,myblock)+yshell1(i,k,7,myblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(yshell1(i,k,2,myblock)+yshell1(i,k,3,myblock))*onecssq
	  f07(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qyy*yshell1(i,k,6,myblock)-cssq*yshell1(i,k,7,myblock)+two*qxy_7_8*yshell1(i,k,8,myblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-yshell1(i,k,2,myblock)+yshell1(i,k,3,myblock))*onecssq
	  f10(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell1(i,k,5,myblock)+qyy*yshell1(i,k,6,myblock)-cssq*yshell1(i,k,7,myblock)+two*qxy_9_10*yshell1(i,k,8,myblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(yshell1(i,k,3,myblock)+yshell1(i,k,4,myblock))*onecssq
	  f11(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1(i,k,6,myblock)+qzz*yshell1(i,k,7,myblock)-cssq*yshell1(i,k,5,myblock)+two*qyz_11_12*yshell1(i,k,10,myblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(yshell1(i,k,3,myblock)-yshell1(i,k,4,myblock))*onecssq
	  f13(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1(i,k,6,myblock)+qzz*yshell1(i,k,7,myblock)-cssq*yshell1(i,k,5,myblock)+two*qyz_13_14*yshell1(i,k,10,myblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d+1 .and. li>=1 .and. li<=TILE_DIMx_d .and. lk>=1 .and. lk<=TILE_DIMz_d)then
      
      
      uu=halfonecssq*(yshell0(i,k,2,myblock)*yshell0(i,k,2,myblock) + yshell0(i,k,3,myblock)*yshell0(i,k,3,myblock) + yshell0(i,k,4,myblock)*yshell0(i,k,4,myblock))
      
      !4  0 +1  0
	  udotc=yshell0(i,k,3,myblock)*onecssq
	  f04(li,lj,lk)=p1*(yshell0(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*yshell0(i,k,6,myblock)-cssq*(yshell0(i,k,5,myblock)+yshell0(i,k,7,myblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(yshell0(i,k,2,myblock)+yshell0(i,k,3,myblock))*onecssq
	  f08(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qyy*yshell0(i,k,6,myblock)-cssq*yshell0(i,k,7,myblock)+two*qxy_7_8*yshell0(i,k,8,myblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-yshell0(i,k,2,myblock)+yshell0(i,k,3,myblock))*onecssq
	  f09(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell0(i,k,5,myblock)+qyy*yshell0(i,k,6,myblock)-cssq*yshell0(i,k,7,myblock)+two*qxy_9_10*yshell0(i,k,8,myblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(yshell0(i,k,3,myblock)+yshell0(i,k,4,myblock))*onecssq
	  f12(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0(i,k,6,myblock)+qzz*yshell0(i,k,7,myblock)-cssq*yshell0(i,k,5,myblock)+two*qyz_11_12*yshell0(i,k,10,myblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(yshell0(i,k,3,myblock)-yshell0(i,k,4,myblock))*onecssq
	  f14(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0(i,k,6,myblock)+qzz*yshell0(i,k,7,myblock)-cssq*yshell0(i,k,5,myblock)+two*qyz_13_14*yshell0(i,k,10,myblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==0 .and. li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d)then
      
      
      uu=halfonecssq*(zshell1(i,j,2,myblock)*zshell1(i,j,2,myblock) + zshell1(i,j,3,myblock)*zshell1(i,j,3,myblock) + zshell1(i,j,4,myblock)*zshell1(i,j,4,myblock))
      
      !5  0  0 -1
	  udotc=zshell1(i,j,4,myblock)*onecssq
	  f05(li,lj,lk)=p1*(zshell1(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*zshell1(i,j,7,myblock)-cssq*(zshell1(i,j,5,myblock)+zshell1(i,j,6,myblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(zshell1(i,j,2,myblock)+zshell1(i,j,4,myblock))*onecssq
	  f15(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qzz*zshell1(i,j,7,myblock)-cssq*zshell1(i,j,6,myblock)+two*qxz_15_16*zshell1(i,j,9,myblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-zshell1(i,j,2,myblock)+zshell1(i,j,4,myblock))*onecssq
	  f17(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell1(i,j,5,myblock)+qzz*zshell1(i,j,7,myblock)-cssq*zshell1(i,j,6,myblock)+two*qxz_17_18*zshell1(i,j,9,myblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(zshell1(i,j,3,myblock)+zshell1(i,j,4,myblock))*onecssq
	  f11(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell1(i,j,6,myblock)+qzz*zshell1(i,j,7,myblock)-cssq*zshell1(i,j,5,myblock)+two*qyz_11_12*zshell1(i,j,10,myblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(zshell1(i,j,3,myblock)-zshell1(i,j,4,myblock))*onecssq
	  f14(li,lj,lk)=p2*(zshell1(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell1(i,j,6,myblock)+qzz*zshell1(i,j,7,myblock)-cssq*zshell1(i,j,5,myblock)+two*qyz_13_14*zshell1(i,j,10,myblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d+1 .and. li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d)then
      
      
      uu=halfonecssq*(zshell0(i,j,2,myblock)*zshell0(i,j,2,myblock) + zshell0(i,j,3,myblock)*zshell0(i,j,3,myblock) + zshell0(i,j,4,myblock)*zshell0(i,j,4,myblock))
      
      !6  0  0  +1
	  udotc=zshell0(i,j,4,myblock)*onecssq
	  f06(li,lj,lk)=p1*(zshell0(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*zshell0(i,j,7,myblock)-cssq*(zshell0(i,j,5,myblock)+zshell0(i,j,6,myblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(zshell0(i,j,2,myblock)+zshell0(i,j,4,myblock))*onecssq
	  f16(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qzz*zshell0(i,j,7,myblock)-cssq*zshell0(i,j,6,myblock)+two*qxz_15_16*zshell0(i,j,9,myblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-zshell0(i,j,2,myblock)+zshell0(i,j,4,myblock))*onecssq
	  f18(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell0(i,j,5,myblock)+qzz*zshell0(i,j,7,myblock)-cssq*zshell0(i,j,6,myblock)+two*qxz_17_18*zshell0(i,j,9,myblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(zshell0(i,j,3,myblock)+zshell0(i,j,4,myblock))*onecssq
	  f12(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell0(i,j,6,myblock)+qzz*zshell0(i,j,7,myblock)-cssq*zshell0(i,j,5,myblock)+two*qyz_11_12*zshell0(i,j,10,myblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(zshell0(i,j,3,myblock)-zshell0(i,j,4,myblock))*onecssq
	  f13(li,lj,lk)=p2*(zshell0(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell0(i,j,6,myblock)+qzz*zshell0(i,j,7,myblock)-cssq*zshell0(i,j,5,myblock)+two*qyz_13_14*zshell0(i,j,10,myblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==0 .and. lj==0)then
      
      
      uu=halfonecssq*(xshell1(j,k,2,myblock)*xshell1(j,k,2,myblock) + xshell1(j,k,3,myblock)*xshell1(j,k,3,myblock) + xshell1(j,k,4,myblock)*xshell1(j,k,4,myblock))
      
      !7 -1 -1  0
	  udotc=(xshell1(j,k,2,myblock)+xshell1(j,k,3,myblock))*onecssq
	  f07(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qyy*xshell1(j,k,6,myblock)-cssq*xshell1(j,k,7,myblock)+two*qxy_7_8*xshell1(j,k,8,myblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==0 .and. lj==TILE_DIMy_d+1)then
      
      
      uu=halfonecssq*(xshell1(j,k,2,myblock)*xshell1(j,k,2,myblock) + xshell1(j,k,3,myblock)*xshell1(j,k,3,myblock) + xshell1(j,k,4,myblock)*xshell1(j,k,4,myblock))
      
      !9  -1 +1 0
      udotc=(-xshell1(j,k,2,myblock)+xshell1(j,k,3,myblock))*onecssq
	  f09(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qyy*xshell1(j,k,6,myblock)-cssq*xshell1(j,k,7,myblock)+two*qxy_9_10*xshell1(j,k,8,myblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==0 .and. lk==0)then
      
      
      uu=halfonecssq*(xshell1(j,k,2,myblock)*xshell1(j,k,2,myblock) + xshell1(j,k,3,myblock)*xshell1(j,k,3,myblock) + xshell1(j,k,4,myblock)*xshell1(j,k,4,myblock))
      
      !15  -1  0 -1
	  udotc=(xshell1(j,k,2,myblock)+xshell1(j,k,4,myblock))*onecssq
	  f15(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qzz*xshell1(j,k,7,myblock)-cssq*xshell1(j,k,6,myblock)+two*qxz_15_16*xshell1(j,k,9,myblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==0 .and. lk==TILE_DIMz_d+1)then
      
      
      uu=halfonecssq*(xshell1(j,k,2,myblock)*xshell1(j,k,2,myblock) + xshell1(j,k,3,myblock)*xshell1(j,k,3,myblock) + xshell1(j,k,4,myblock)*xshell1(j,k,4,myblock))
      
      !18   -1   0  +1
	  udotc=(-xshell1(j,k,2,myblock)+xshell1(j,k,4,myblock))*onecssq
	  f18(li,lj,lk)=p2*(xshell1(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1(j,k,5,myblock)+qzz*xshell1(j,k,7,myblock)-cssq*xshell1(j,k,6,myblock)+two*qxz_17_18*xshell1(j,k,9,myblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d+1 .and. lj==0)then

      
      uu=halfonecssq*(xshell0(j,k,2,myblock)*xshell0(j,k,2,myblock) + xshell0(j,k,3,myblock)*xshell0(j,k,3,myblock) + xshell0(j,k,4,myblock)*xshell0(j,k,4,myblock))
      
      !10   +1 -1  0
	  udotc=(-xshell0(j,k,2,myblock)+xshell0(j,k,3,myblock))*onecssq
	  f10(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qyy*xshell0(j,k,6,myblock)-cssq*xshell0(j,k,7,myblock)+two*qxy_9_10*xshell0(j,k,8,myblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d+1 .and. lj==TILE_DIMy_d+1)then

      
      uu=halfonecssq*(xshell0(j,k,2,myblock)*xshell0(j,k,2,myblock) + xshell0(j,k,3,myblock)*xshell0(j,k,3,myblock) + xshell0(j,k,4,myblock)*xshell0(j,k,4,myblock))
      
      !8 +1 +1  0
	  udotc=(xshell0(j,k,2,myblock)+xshell0(j,k,3,myblock))*onecssq
	  f08(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qyy*xshell0(j,k,6,myblock)-cssq*xshell0(j,k,7,myblock)+two*qxy_7_8*xshell0(j,k,8,myblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d+1 .and. lk==0)then
      
      
      uu=halfonecssq*(xshell0(j,k,2,myblock)*xshell0(j,k,2,myblock) + xshell0(j,k,3,myblock)*xshell0(j,k,3,myblock) + xshell0(j,k,4,myblock)*xshell0(j,k,4,myblock))
      
      !17  +1  0 -1
	  udotc=(-xshell0(j,k,2,myblock)+xshell0(j,k,4,myblock))*onecssq
	  f17(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qzz*xshell0(j,k,7,myblock)-cssq*xshell0(j,k,6,myblock)+two*qxz_17_18*xshell0(j,k,9,myblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d+1 .and. lk==TILE_DIMz_d+1)then

      
      uu=halfonecssq*(xshell0(j,k,2,myblock)*xshell0(j,k,2,myblock) + xshell0(j,k,3,myblock)*xshell0(j,k,3,myblock) + xshell0(j,k,4,myblock)*xshell0(j,k,4,myblock))
      
      !16  +1  0 +1
	  udotc=(xshell0(j,k,2,myblock)+xshell0(j,k,4,myblock))*onecssq
	  f16(li,lj,lk)=p2*(xshell0(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0(j,k,5,myblock)+qzz*xshell0(j,k,7,myblock)-cssq*xshell0(j,k,6,myblock)+two*qxz_15_16*xshell0(j,k,9,myblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==0 .and. lk==0)then

      
      uu=halfonecssq*(yshell1(i,k,2,myblock)*yshell1(i,k,2,myblock) + yshell1(i,k,3,myblock)*yshell1(i,k,3,myblock) + yshell1(i,k,4,myblock)*yshell1(i,k,4,myblock))
      
      !11  0  -1  -1
	  udotc=(yshell1(i,k,3,myblock)+yshell1(i,k,4,myblock))*onecssq
	  f11(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1(i,k,6,myblock)+qzz*yshell1(i,k,7,myblock)-cssq*yshell1(i,k,5,myblock)+two*qyz_11_12*yshell1(i,k,10,myblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==0 .and. lk==TILE_DIMz_d+1)then

      
      uu=halfonecssq*(yshell1(i,k,2,myblock)*yshell1(i,k,2,myblock) + yshell1(i,k,3,myblock)*yshell1(i,k,3,myblock) + yshell1(i,k,4,myblock)*yshell1(i,k,4,myblock))
      
      !13  0  -1   +1
	  udotc=(yshell1(i,k,3,myblock)-yshell1(i,k,4,myblock))*onecssq
	  f13(li,lj,lk)=p2*(yshell1(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1(i,k,6,myblock)+qzz*yshell1(i,k,7,myblock)-cssq*yshell1(i,k,5,myblock)+two*qyz_13_14*yshell1(i,k,10,myblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d+1 .and. lk==0)then

      
      uu=halfonecssq*(yshell0(i,k,2,myblock)*yshell0(i,k,2,myblock) + yshell0(i,k,3,myblock)*yshell0(i,k,3,myblock) + yshell0(i,k,4,myblock)*yshell0(i,k,4,myblock))
      
      !14  0  +1  -1
	  udotc=(yshell0(i,k,3,myblock)-yshell0(i,k,4,myblock))*onecssq
	  f14(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0(i,k,6,myblock)+qzz*yshell0(i,k,7,myblock)-cssq*yshell0(i,k,5,myblock)+two*qyz_13_14*yshell0(i,k,10,myblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d+1 .and. lk==TILE_DIMz_d+1)then

      
      uu=halfonecssq*(yshell0(i,k,2,myblock)*yshell0(i,k,2,myblock) + yshell0(i,k,3,myblock)*yshell0(i,k,3,myblock) + yshell0(i,k,4,myblock)*yshell0(i,k,4,myblock))
      
      !12  0  +1  +1
	  udotc=(yshell0(i,k,3,myblock)+yshell0(i,k,4,myblock))*onecssq
	  f12(li,lj,lk)=p2*(yshell0(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0(i,k,6,myblock)+qzz*yshell0(i,k,7,myblock)-cssq*yshell0(i,k,5,myblock)+two*qyz_11_12*yshell0(i,k,10,myblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
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
	hfields(i,j,k,1,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,2,myblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,3,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,4,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,5,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(i,j,k,6,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,7,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(i,j,k,8,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,9,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,10,myblock)=udotc
     
#ifdef PRESSCORR

	!call syncthreads
	
	uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
	+ hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,5,myblock)=hfields(i,j,k,5,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(i,j,k,6,myblock)=hfields(i,j,k,6,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,7,myblock)=hfields(i,j,k,7,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(i,j,k,8,myblock)=hfields(i,j,k,8,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(i,j,k,9,myblock)=hfields(i,j,k,9,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(i,j,k,10,myblock)=hfields(i,j,k,10,myblock)-udotc
	
#endif	
    
    ! Halo Faces
    if(li==1)then
      
      xshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(li==TILE_DIMx_d)then
      
      xshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
	  
    endif

    if(lj==1)then
      
      yshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lj==TILE_DIMy_d)then
      
      yshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    if(lk==1)then
      
      zshell0h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell0h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell0h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell0h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell0h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell0h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell0h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell0h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell0h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell0h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lk==TILE_DIMz_d)then
      
      zshell1h(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell1h(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell1h(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell1h(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell1h(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell1h(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell1h(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell1h(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell1h(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell1h(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    return
    
  end subroutine streamcoll_minshared_halo
  
    attributes(global) subroutine streamcoll_minshared_halo_flop()
	
	implicit none  
	  
    integer :: i,j,k,gi,gj,gk,myblock
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
    
    !inner cube
    if(li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d .and. &
     lk>=1 .and. lk<=TILE_DIMz_d)then
     
    uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
     + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(hfields(i,j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
	
    endif
    
    ! Halo Faces
    if(li==0 .and. lj>=1 .and. lj<=TILE_DIMy_d .and. lk>=1 .and. lk<=TILE_DIMz_d)then
      
      
      uu=halfonecssq*(xshell1h(j,k,2,myblock)*xshell1h(j,k,2,myblock) + xshell1h(j,k,3,myblock)*xshell1h(j,k,3,myblock) + xshell1h(j,k,4,myblock)*xshell1h(j,k,4,myblock))
      
      !1 -1  0  0
	  udotc=xshell1h(j,k,2,myblock)*onecssq
	  f01(li,lj,lk)=p1*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*xshell1h(j,k,5,myblock)-cssq*(xshell1h(j,k,6,myblock)+xshell1h(j,k,7,myblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(xshell1h(j,k,2,myblock)+xshell1h(j,k,3,myblock))*onecssq
	  f07(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qyy*xshell1h(j,k,6,myblock)-cssq*xshell1h(j,k,7,myblock)+two*qxy_7_8*xshell1h(j,k,8,myblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-xshell1h(j,k,2,myblock)+xshell1h(j,k,3,myblock))*onecssq
	  f09(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qyy*xshell1h(j,k,6,myblock)-cssq*xshell1h(j,k,7,myblock)+two*qxy_9_10*xshell1h(j,k,8,myblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(xshell1h(j,k,2,myblock)+xshell1h(j,k,4,myblock))*onecssq
	  f15(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qzz*xshell1h(j,k,7,myblock)-cssq*xshell1h(j,k,6,myblock)+two*qxz_15_16*xshell1h(j,k,9,myblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-xshell1h(j,k,2,myblock)+xshell1h(j,k,4,myblock))*onecssq
	  f18(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qzz*xshell1h(j,k,7,myblock)-cssq*xshell1h(j,k,6,myblock)+two*qxz_17_18*xshell1h(j,k,9,myblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    
    if(li==TILE_DIMx_d+1 .and. lj>=1 .and. lj<=TILE_DIMy_d .and. lk>=1 .and. lk<=TILE_DIMz_d)then
      
      
      uu=halfonecssq*(xshell0h(j,k,2,myblock)*xshell0h(j,k,2,myblock) + xshell0h(j,k,3,myblock)*xshell0h(j,k,3,myblock) + xshell0h(j,k,4,myblock)*xshell0h(j,k,4,myblock))
      
      !2 +1  0  0
	  udotc=xshell0h(j,k,2,myblock)*onecssq
	  f02(li,lj,lk)=p1*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*xshell0h(j,k,5,myblock)-cssq*(xshell0h(j,k,6,myblock)+xshell0h(j,k,7,myblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(xshell0h(j,k,2,myblock)+xshell0h(j,k,3,myblock))*onecssq
	  f08(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qyy*xshell0h(j,k,6,myblock)-cssq*xshell0h(j,k,7,myblock)+two*qxy_7_8*xshell0h(j,k,8,myblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-xshell0h(j,k,2,myblock)+xshell0h(j,k,3,myblock))*onecssq
	  f10(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qyy*xshell0h(j,k,6,myblock)-cssq*xshell0h(j,k,7,myblock)+two*qxy_9_10*xshell0h(j,k,8,myblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(xshell0h(j,k,2,myblock)+xshell0h(j,k,4,myblock))*onecssq
	  f16(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qzz*xshell0h(j,k,7,myblock)-cssq*xshell0h(j,k,6,myblock)+two*qxz_15_16*xshell0h(j,k,9,myblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-xshell0h(j,k,2,myblock)+xshell0h(j,k,4,myblock))*onecssq
	  f17(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qzz*xshell0h(j,k,7,myblock)-cssq*xshell0h(j,k,6,myblock)+two*qxz_17_18*xshell0h(j,k,9,myblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==0 .and. li>=1 .and. li<=TILE_DIMx_d .and. lk>=1 .and. lk<=TILE_DIMz_d)then
      
      
      uu=halfonecssq*(yshell1h(i,k,2,myblock)*yshell1h(i,k,2,myblock) + yshell1h(i,k,3,myblock)*yshell1h(i,k,3,myblock) + yshell1h(i,k,4,myblock)*yshell1h(i,k,4,myblock))
      
      !3 0 -1  0
	  udotc=yshell1h(i,k,3,myblock)*onecssq
	  f03(li,lj,lk)=p1*(yshell1h(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*yshell1h(i,k,6,myblock)-cssq*(yshell1h(i,k,5,myblock)+yshell1h(i,k,7,myblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(yshell1h(i,k,2,myblock)+yshell1h(i,k,3,myblock))*onecssq
	  f07(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qyy*yshell1h(i,k,6,myblock)-cssq*yshell1h(i,k,7,myblock)+two*qxy_7_8*yshell1h(i,k,8,myblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-yshell1h(i,k,2,myblock)+yshell1h(i,k,3,myblock))*onecssq
	  f10(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell1h(i,k,5,myblock)+qyy*yshell1h(i,k,6,myblock)-cssq*yshell1h(i,k,7,myblock)+two*qxy_9_10*yshell1h(i,k,8,myblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(yshell1h(i,k,3,myblock)+yshell1h(i,k,4,myblock))*onecssq
	  f11(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1h(i,k,6,myblock)+qzz*yshell1h(i,k,7,myblock)-cssq*yshell1h(i,k,5,myblock)+two*qyz_11_12*yshell1h(i,k,10,myblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(yshell1h(i,k,3,myblock)-yshell1h(i,k,4,myblock))*onecssq
	  f13(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1h(i,k,6,myblock)+qzz*yshell1h(i,k,7,myblock)-cssq*yshell1h(i,k,5,myblock)+two*qyz_13_14*yshell1h(i,k,10,myblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d+1 .and. li>=1 .and. li<=TILE_DIMx_d .and. lk>=1 .and. lk<=TILE_DIMz_d)then
      
      
      uu=halfonecssq*(yshell0h(i,k,2,myblock)*yshell0h(i,k,2,myblock) + yshell0h(i,k,3,myblock)*yshell0h(i,k,3,myblock) + yshell0h(i,k,4,myblock)*yshell0h(i,k,4,myblock))
      
      !4  0 +1  0
	  udotc=yshell0h(i,k,3,myblock)*onecssq
	  f04(li,lj,lk)=p1*(yshell0h(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*yshell0h(i,k,6,myblock)-cssq*(yshell0h(i,k,5,myblock)+yshell0h(i,k,7,myblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(yshell0h(i,k,2,myblock)+yshell0h(i,k,3,myblock))*onecssq
	  f08(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qyy*yshell0h(i,k,6,myblock)-cssq*yshell0h(i,k,7,myblock)+two*qxy_7_8*yshell0h(i,k,8,myblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-yshell0h(i,k,2,myblock)+yshell0h(i,k,3,myblock))*onecssq
	  f09(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*yshell0h(i,k,5,myblock)+qyy*yshell0h(i,k,6,myblock)-cssq*yshell0h(i,k,7,myblock)+two*qxy_9_10*yshell0h(i,k,8,myblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(yshell0h(i,k,3,myblock)+yshell0h(i,k,4,myblock))*onecssq
	  f12(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0h(i,k,6,myblock)+qzz*yshell0h(i,k,7,myblock)-cssq*yshell0h(i,k,5,myblock)+two*qyz_11_12*yshell0h(i,k,10,myblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(yshell0h(i,k,3,myblock)-yshell0h(i,k,4,myblock))*onecssq
	  f14(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0h(i,k,6,myblock)+qzz*yshell0h(i,k,7,myblock)-cssq*yshell0h(i,k,5,myblock)+two*qyz_13_14*yshell0h(i,k,10,myblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==0 .and. li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d)then
      
      
      uu=halfonecssq*(zshell1h(i,j,2,myblock)*zshell1h(i,j,2,myblock) + zshell1h(i,j,3,myblock)*zshell1h(i,j,3,myblock) + zshell1h(i,j,4,myblock)*zshell1h(i,j,4,myblock))
      
      !5  0  0 -1
	  udotc=zshell1h(i,j,4,myblock)*onecssq
	  f05(li,lj,lk)=p1*(zshell1h(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*zshell1h(i,j,7,myblock)-cssq*(zshell1h(i,j,5,myblock)+zshell1h(i,j,6,myblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(zshell1h(i,j,2,myblock)+zshell1h(i,j,4,myblock))*onecssq
	  f15(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qzz*zshell1h(i,j,7,myblock)-cssq*zshell1h(i,j,6,myblock)+two*qxz_15_16*zshell1h(i,j,9,myblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-zshell1h(i,j,2,myblock)+zshell1h(i,j,4,myblock))*onecssq
	  f17(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell1h(i,j,5,myblock)+qzz*zshell1h(i,j,7,myblock)-cssq*zshell1h(i,j,6,myblock)+two*qxz_17_18*zshell1h(i,j,9,myblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(zshell1h(i,j,3,myblock)+zshell1h(i,j,4,myblock))*onecssq
	  f11(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell1h(i,j,6,myblock)+qzz*zshell1h(i,j,7,myblock)-cssq*zshell1h(i,j,5,myblock)+two*qyz_11_12*zshell1h(i,j,10,myblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(zshell1h(i,j,3,myblock)-zshell1h(i,j,4,myblock))*onecssq
	  f14(li,lj,lk)=p2*(zshell1h(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell1h(i,j,6,myblock)+qzz*zshell1h(i,j,7,myblock)-cssq*zshell1h(i,j,5,myblock)+two*qyz_13_14*zshell1h(i,j,10,myblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d+1 .and. li>=1 .and. li<=TILE_DIMx_d .and. lj>=1 .and. lj<=TILE_DIMy_d)then
      
      
      uu=halfonecssq*(zshell0h(i,j,2,myblock)*zshell0h(i,j,2,myblock) + zshell0h(i,j,3,myblock)*zshell0h(i,j,3,myblock) + zshell0h(i,j,4,myblock)*zshell0h(i,j,4,myblock))
      
      !6  0  0  +1
	  udotc=zshell0h(i,j,4,myblock)*onecssq
	  f06(li,lj,lk)=p1*(zshell0h(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*zshell0h(i,j,7,myblock)-cssq*(zshell0h(i,j,5,myblock)+zshell0h(i,j,6,myblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(zshell0h(i,j,2,myblock)+zshell0h(i,j,4,myblock))*onecssq
	  f16(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qzz*zshell0h(i,j,7,myblock)-cssq*zshell0h(i,j,6,myblock)+two*qxz_15_16*zshell0h(i,j,9,myblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-zshell0h(i,j,2,myblock)+zshell0h(i,j,4,myblock))*onecssq
	  f18(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*zshell0h(i,j,5,myblock)+qzz*zshell0h(i,j,7,myblock)-cssq*zshell0h(i,j,6,myblock)+two*qxz_17_18*zshell0h(i,j,9,myblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(zshell0h(i,j,3,myblock)+zshell0h(i,j,4,myblock))*onecssq
	  f12(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell0h(i,j,6,myblock)+qzz*zshell0h(i,j,7,myblock)-cssq*zshell0h(i,j,5,myblock)+two*qyz_11_12*zshell0h(i,j,10,myblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(zshell0h(i,j,3,myblock)-zshell0h(i,j,4,myblock))*onecssq
	  f13(li,lj,lk)=p2*(zshell0h(i,j,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*zshell0h(i,j,6,myblock)+qzz*zshell0h(i,j,7,myblock)-cssq*zshell0h(i,j,5,myblock)+two*qyz_13_14*zshell0h(i,j,10,myblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==0 .and. lj==0)then
      
      
      uu=halfonecssq*(xshell1h(j,k,2,myblock)*xshell1h(j,k,2,myblock) + xshell1h(j,k,3,myblock)*xshell1h(j,k,3,myblock) + xshell1h(j,k,4,myblock)*xshell1h(j,k,4,myblock))
      
      !7 -1 -1  0
	  udotc=(xshell1h(j,k,2,myblock)+xshell1h(j,k,3,myblock))*onecssq
	  f07(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qyy*xshell1h(j,k,6,myblock)-cssq*xshell1h(j,k,7,myblock)+two*qxy_7_8*xshell1h(j,k,8,myblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==0 .and. lj==TILE_DIMy_d+1)then
      
      
      uu=halfonecssq*(xshell1h(j,k,2,myblock)*xshell1h(j,k,2,myblock) + xshell1h(j,k,3,myblock)*xshell1h(j,k,3,myblock) + xshell1h(j,k,4,myblock)*xshell1h(j,k,4,myblock))
      
      !9  -1 +1 0
      udotc=(-xshell1h(j,k,2,myblock)+xshell1h(j,k,3,myblock))*onecssq
	  f09(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qyy*xshell1h(j,k,6,myblock)-cssq*xshell1h(j,k,7,myblock)+two*qxy_9_10*xshell1h(j,k,8,myblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==0 .and. lk==0)then
      
      
      uu=halfonecssq*(xshell1h(j,k,2,myblock)*xshell1h(j,k,2,myblock) + xshell1h(j,k,3,myblock)*xshell1h(j,k,3,myblock) + xshell1h(j,k,4,myblock)*xshell1h(j,k,4,myblock))
      
      !15  -1  0 -1
	  udotc=(xshell1h(j,k,2,myblock)+xshell1h(j,k,4,myblock))*onecssq
	  f15(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qzz*xshell1h(j,k,7,myblock)-cssq*xshell1h(j,k,6,myblock)+two*qxz_15_16*xshell1h(j,k,9,myblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==0 .and. lk==TILE_DIMz_d+1)then
      
      
      uu=halfonecssq*(xshell1h(j,k,2,myblock)*xshell1h(j,k,2,myblock) + xshell1h(j,k,3,myblock)*xshell1h(j,k,3,myblock) + xshell1h(j,k,4,myblock)*xshell1h(j,k,4,myblock))
      
      !18   -1   0  +1
	  udotc=(-xshell1h(j,k,2,myblock)+xshell1h(j,k,4,myblock))*onecssq
	  f18(li,lj,lk)=p2*(xshell1h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell1h(j,k,5,myblock)+qzz*xshell1h(j,k,7,myblock)-cssq*xshell1h(j,k,6,myblock)+two*qxz_17_18*xshell1h(j,k,9,myblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d+1 .and. lj==0)then

      
      uu=halfonecssq*(xshell0h(j,k,2,myblock)*xshell0h(j,k,2,myblock) + xshell0h(j,k,3,myblock)*xshell0h(j,k,3,myblock) + xshell0h(j,k,4,myblock)*xshell0h(j,k,4,myblock))
      
      !10   +1 -1  0
	  udotc=(-xshell0h(j,k,2,myblock)+xshell0h(j,k,3,myblock))*onecssq
	  f10(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qyy*xshell0h(j,k,6,myblock)-cssq*xshell0h(j,k,7,myblock)+two*qxy_9_10*xshell0h(j,k,8,myblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d+1 .and. lj==TILE_DIMy_d+1)then

      
      uu=halfonecssq*(xshell0h(j,k,2,myblock)*xshell0h(j,k,2,myblock) + xshell0h(j,k,3,myblock)*xshell0h(j,k,3,myblock) + xshell0h(j,k,4,myblock)*xshell0h(j,k,4,myblock))
      
      !8 +1 +1  0
	  udotc=(xshell0h(j,k,2,myblock)+xshell0h(j,k,3,myblock))*onecssq
	  f08(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qyy*xshell0h(j,k,6,myblock)-cssq*xshell0h(j,k,7,myblock)+two*qxy_7_8*xshell0h(j,k,8,myblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d+1 .and. lk==0)then
      
      
      uu=halfonecssq*(xshell0h(j,k,2,myblock)*xshell0h(j,k,2,myblock) + xshell0h(j,k,3,myblock)*xshell0h(j,k,3,myblock) + xshell0h(j,k,4,myblock)*xshell0h(j,k,4,myblock))
      
      !17  +1  0 -1
	  udotc=(-xshell0h(j,k,2,myblock)+xshell0h(j,k,4,myblock))*onecssq
	  f17(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qzz*xshell0h(j,k,7,myblock)-cssq*xshell0h(j,k,6,myblock)+two*qxz_17_18*xshell0h(j,k,9,myblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d+1 .and. lk==TILE_DIMz_d+1)then

      
      uu=halfonecssq*(xshell0h(j,k,2,myblock)*xshell0h(j,k,2,myblock) + xshell0h(j,k,3,myblock)*xshell0h(j,k,3,myblock) + xshell0h(j,k,4,myblock)*xshell0h(j,k,4,myblock))
      
      !16  +1  0 +1
	  udotc=(xshell0h(j,k,2,myblock)+xshell0h(j,k,4,myblock))*onecssq
	  f16(li,lj,lk)=p2*(xshell0h(j,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*xshell0h(j,k,5,myblock)+qzz*xshell0h(j,k,7,myblock)-cssq*xshell0h(j,k,6,myblock)+two*qxz_15_16*xshell0h(j,k,9,myblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==0 .and. lk==0)then

      
      uu=halfonecssq*(yshell1h(i,k,2,myblock)*yshell1h(i,k,2,myblock) + yshell1h(i,k,3,myblock)*yshell1h(i,k,3,myblock) + yshell1h(i,k,4,myblock)*yshell1h(i,k,4,myblock))
      
      !11  0  -1  -1
	  udotc=(yshell1h(i,k,3,myblock)+yshell1h(i,k,4,myblock))*onecssq
	  f11(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1h(i,k,6,myblock)+qzz*yshell1h(i,k,7,myblock)-cssq*yshell1h(i,k,5,myblock)+two*qyz_11_12*yshell1h(i,k,10,myblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==0 .and. lk==TILE_DIMz_d+1)then

      
      uu=halfonecssq*(yshell1h(i,k,2,myblock)*yshell1h(i,k,2,myblock) + yshell1h(i,k,3,myblock)*yshell1h(i,k,3,myblock) + yshell1h(i,k,4,myblock)*yshell1h(i,k,4,myblock))
      
      !13  0  -1   +1
	  udotc=(yshell1h(i,k,3,myblock)-yshell1h(i,k,4,myblock))*onecssq
	  f13(li,lj,lk)=p2*(yshell1h(i,k,1,myblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell1h(i,k,6,myblock)+qzz*yshell1h(i,k,7,myblock)-cssq*yshell1h(i,k,5,myblock)+two*qyz_13_14*yshell1h(i,k,10,myblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d+1 .and. lk==0)then

      
      uu=halfonecssq*(yshell0h(i,k,2,myblock)*yshell0h(i,k,2,myblock) + yshell0h(i,k,3,myblock)*yshell0h(i,k,3,myblock) + yshell0h(i,k,4,myblock)*yshell0h(i,k,4,myblock))
      
      !14  0  +1  -1
	  udotc=(yshell0h(i,k,3,myblock)-yshell0h(i,k,4,myblock))*onecssq
	  f14(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0h(i,k,6,myblock)+qzz*yshell0h(i,k,7,myblock)-cssq*yshell0h(i,k,5,myblock)+two*qyz_13_14*yshell0h(i,k,10,myblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d+1 .and. lk==TILE_DIMz_d+1)then

      
      uu=halfonecssq*(yshell0h(i,k,2,myblock)*yshell0h(i,k,2,myblock) + yshell0h(i,k,3,myblock)*yshell0h(i,k,3,myblock) + yshell0h(i,k,4,myblock)*yshell0h(i,k,4,myblock))
      
      !12  0  +1  +1
	  udotc=(yshell0h(i,k,3,myblock)+yshell0h(i,k,4,myblock))*onecssq
	  f12(li,lj,lk)=p2*(yshell0h(i,k,1,myblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*yshell0h(i,k,6,myblock)+qzz*yshell0h(i,k,7,myblock)-cssq*yshell0h(i,k,5,myblock)+two*qyz_11_12*yshell0h(i,k,10,myblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
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
	hfields(i,j,k,1,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,2,myblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,3,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,4,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,5,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(i,j,k,6,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,7,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(i,j,k,8,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,9,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,10,myblock)=udotc
     
#ifdef PRESSCORR

	!call syncthreads
	
	uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
	+ hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,5,myblock)=hfields(i,j,k,5,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(i,j,k,6,myblock)=hfields(i,j,k,6,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,7,myblock)=hfields(i,j,k,7,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(i,j,k,8,myblock)=hfields(i,j,k,8,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(i,j,k,9,myblock)=hfields(i,j,k,9,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(i,j,k,10,myblock)=hfields(i,j,k,10,myblock)-udotc
	
#endif	
    
    ! Halo Faces
    if(li==1)then
      
      xshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(li==TILE_DIMx_d)then
      
      xshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      xshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      xshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      xshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      xshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      xshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      xshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      xshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      xshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      xshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
	  
    endif

    if(lj==1)then
      
      yshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lj==TILE_DIMy_d)then
      
      yshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      yshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      yshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      yshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      yshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      yshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      yshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      yshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      yshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      yshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    if(lk==1)then
      
      zshell0(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell0(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell0(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell0(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell0(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell0(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell0(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell0(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell0(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell0(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif
    if(lk==TILE_DIMz_d)then
      
      zshell1(j,k,1,myblock)=hfields(i,j,k,1,myblock)
      zshell1(j,k,2,myblock)=hfields(i,j,k,2,myblock)
      zshell1(j,k,3,myblock)=hfields(i,j,k,3,myblock)
      zshell1(j,k,4,myblock)=hfields(i,j,k,4,myblock)
      zshell1(j,k,5,myblock)=hfields(i,j,k,5,myblock)
      zshell1(j,k,6,myblock)=hfields(i,j,k,6,myblock)
      zshell1(j,k,7,myblock)=hfields(i,j,k,7,myblock)
      zshell1(j,k,8,myblock)=hfields(i,j,k,8,myblock)
      zshell1(j,k,9,myblock)=hfields(i,j,k,9,myblock)
      zshell1(j,k,10,myblock)=hfields(i,j,k,10,myblock)
      
    endif

    return
    
  end subroutine streamcoll_minshared_halo_flop
  
 end module streamcoll_kernels
