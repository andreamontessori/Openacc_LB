#include "defines.h"
 module streamcoll_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine streamcoll_shared()
	
	implicit none  
	
	integer :: i,j,k
    integer :: gi,gj,gk,myblock,iidblock
    integer :: gii,gjj,gkk,ii,jj,kk
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
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
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
	
	
	
!	 myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
!    zblock=(idblock-1)/nxyblock +1
!    yblock=((idblock-1)-(zblock-1)*nxyblock)/nxblock +1
!    xblock=(idblock-1)-(zblock-1)*nxyblock-(yblock-1)*nxblock +1
	
	!myblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	!iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	!if(myblock.ne.iidblock)write(*,*)'cazzo1',gi,gj,gk,myblock,iidblock
    
    uu=halfonecssq*(u(i,j,k,myblock)*u(i,j,k,myblock) + v(i,j,k,myblock)*v(i,j,k,myblock) + w(i,j,k,myblock)*w(i,j,k,myblock))
    
    !0
	f00(i,j,k)=p0*(rho(i,j,k,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(pyy(i,j,k,myblock)+pxx(i,j,k,myblock)+pzz(i,j,k,myblock)))
	
    
	!1 -1  0  0
	udotc=u(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(i,j,k)=p1*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxx(i,j,k,myblock)-cssq*(pyy(i,j,k,myblock)+pzz(i,j,k,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(i,j,k)=p1*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxx(i,j,k,myblock)-cssq*(pyy(i,j,k,myblock)+pzz(i,j,k,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=v(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(i,j,k)=p1*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyy(i,j,k,myblock)-cssq*(pxx(i,j,k,myblock)+pzz(i,j,k,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(i,j,k)=p1*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyy(i,j,k,myblock)-cssq*(pxx(i,j,k,myblock)+pzz(i,j,k,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=w(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(i,j,k)=p1*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k,myblock)-cssq*(pxx(i,j,k,myblock)+pyy(i,j,k,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(i,j,k)=p1*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k,myblock)-cssq*(pxx(i,j,k,myblock)+pyy(i,j,k,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(u(i,j,k,myblock)+v(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qyy*pyy(i,j,k,myblock)-cssq*pzz(i,j,k,myblock)+two*qxy_7_8*pxy(i,j,k,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qyy*pyy(i,j,k,myblock)-cssq*pzz(i,j,k,myblock)+two*qxy_7_8*pxy(i,j,k,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-u(i,j,k,myblock)+v(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qyy*pyy(i,j,k,myblock)-cssq*pzz(i,j,k,myblock)+two*qxy_9_10*pxy(i,j,k,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qyy*pyy(i,j,k,myblock)-cssq*pzz(i,j,k,myblock)+two*qxy_9_10*pxy(i,j,k,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(u(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pyy(i,j,k,myblock)+two*qxz_15_16*pxz(i,j,k,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pyy(i,j,k,myblock)+two*qxz_15_16*pxz(i,j,k,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-u(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pyy(i,j,k,myblock)+two*qxz_17_18*pxz(i,j,k,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pyy(i,j,k,myblock)+two*qxz_17_18*pxz(i,j,k,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(v(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pxx(i,j,k,myblock)+two*qyz_11_12*pyz(i,j,k,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pxx(i,j,k,myblock)+two*qyz_11_12*pyz(i,j,k,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(v(i,j,k,myblock)-w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pxx(i,j,k,myblock)+two*qyz_13_14*pyz(i,j,k,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pxx(i,j,k,myblock)+two*qyz_13_14*pyz(i,j,k,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    ! Halo Faces
    if(i==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !1 -1  0  0
	  udotc=u(ii,jj,kk,iidblock)*onecssq
	  f01(i-1,j,k)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxx(ii,jj,kk,iidblock)-cssq*(pyy(ii,jj,kk,iidblock)+pzz(ii,jj,kk,iidblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f07(i-1,j,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f09(i-1,j,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f15(i-1,j,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f18(i-1,j,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(i==TILE_DIMx_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !2 +1  0  0
	  udotc=u(ii,jj,kk,iidblock)*onecssq
	  f02(i+1,j,k)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxx(ii,jj,kk,iidblock)-cssq*(pyy(ii,jj,kk,iidblock)+pzz(ii,jj,kk,iidblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f08(i+1,j,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f10(i+1,j,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f16(i+1,j,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f17(i+1,j,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(j==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !3 0 -1  0
	  udotc=v(ii,jj,kk,iidblock)*onecssq
	  f03(i,j-1,k)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyy(ii,jj,kk,iidblock)-cssq*(pxx(ii,jj,kk,iidblock)+pzz(ii,jj,kk,iidblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f07(i,j-1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f10(i,j-1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f11(i,j-1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f13(i,j-1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(j==TILE_DIMy_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !4  0 +1  0
	  udotc=v(ii,jj,kk,iidblock)*onecssq
	  f04(i,j+1,k)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyy(ii,jj,kk,iidblock)-cssq*(pxx(ii,jj,kk,iidblock)+pzz(ii,jj,kk,iidblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f08(i,j+1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f09(i,j+1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f12(i,j+1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f14(i,j+1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(k==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !5  0  0 -1
	  udotc=w(ii,jj,kk,iidblock)*onecssq
	  f05(i,j,k-1)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzz(ii,jj,kk,iidblock)-cssq*(pxx(ii,jj,kk,iidblock)+pyy(ii,jj,kk,iidblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f15(i,j,k-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f17(i,j,k-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f11(i,j,k-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f14(i,j,k-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !6  0  0  +1
	  udotc=w(ii,jj,kk,iidblock)*onecssq
	  f06(i,j,k+1)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzz(ii,jj,kk,iidblock)-cssq*(pxx(ii,jj,kk,iidblock)+pyy(ii,jj,kk,iidblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f16(i,j,k+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f18(i,j,k+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f12(i,j,k+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f13(i,j,k+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(i==1 .and. j==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !7 -1 -1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f07(i-1,j-1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(i==1 .and. j==TILE_DIMy_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !9  -1 +1 0
      udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f09(i-1,j+1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(i==1 .and. k==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !15  -1  0 -1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f15(i-1,j,k-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(i==1 .and. k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !18   -1   0  +1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f18(i-1,j,k+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(i==TILE_DIMx_d .and. j==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !10   +1 -1  0
	  udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f10(i+1,j-1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(i==TILE_DIMx_d .and. j==TILE_DIMy_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !8 +1 +1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f08(i+1,j+1,k)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(i==TILE_DIMx_d .and. k==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !17  +1  0 -1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f17(i+1,j,k-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(i==TILE_DIMx_d .and. k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !16  +1  0 +1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f16(i+1,j,k+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(j==1 .and. k==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !11  0  -1  -1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f11(i,j-1,k-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(j==1 .and. k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !13  0  -1   +1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f13(i,j-1,k+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(j==TILE_DIMy_d .and. k==1)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !14  0  +1  -1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f14(i,j+1,k-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(j==TILE_DIMy_d .and. k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !12  0  +1  +1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f12(i,j+1,k+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
    
    if(abs(isfluid(gi,gj,gk)).ne.1)return
    
    udotc=f00(i,j,k)+f01(i-1,j,k)+f02(i+1,j,k)+  &
     f03(i,j-1,k)+f04(i,j+1,k)+  &
     f05(i,j,k-1)+f06(i,j,k+1)+  &
     f07(i-1,j-1,k)+f08(i+1,j+1,k)+ &
     f09(i-1,j+1,k)+f10(i+1,j-1,k)+ &
     f11(i,j-1,k-1)+f12(i,j+1,k+1)+ &
     f13(i,j-1,k+1)+f14(i,j+1,k-1)+ &
     f15(i-1,j,k-1)+f16(i+1,j,k+1)+ &
     f17(i+1,j,k-1)+f18(i-1,j,k+1)
	rhoh(i,j,k,myblock)=udotc
	
	udotc=f01(i-1,j,k)-f02(i+1,j,k)+  &
     f07(i-1,j-1,k)-f08(i+1,j+1,k)+ &
     f09(i-1,j+1,k)-f10(i+1,j-1,k)+ &
     f15(i-1,j,k-1)-f16(i+1,j,k+1)- &
     f17(i+1,j,k-1)+f18(i-1,j,k+1)
	uh(i,j,k,myblock)=udotc
	
	
	udotc=f03(i,j-1,k)-f04(i,j+1,k)+ &
     f07(i-1,j-1,k)-f08(i+1,j+1,k)- &
     f09(i-1,j+1,k)+f10(i+1,j-1,k)+ &
     f11(i,j-1,k-1)-f12(i,j+1,k+1)+ &
     f13(i,j-1,k+1)-f14(i,j+1,k-1)
	vh(i,j,k,myblock)=udotc
	
	udotc=f05(i,j,k-1)-f06(i,j,k+1)+  &
     f11(i,j-1,k-1)-f12(i,j+1,k+1)- &
     f13(i,j-1,k+1)+f14(i,j+1,k-1)+ &
     f15(i-1,j,k-1)-f16(i+1,j,k+1)+ &
     f17(i+1,j,k-1)-f18(i-1,j,k+1)
	wh(i,j,k,myblock)=udotc
	
	udotc=f01(i-1,j,k)+f02(i+1,j,k)+  &
     f07(i-1,j-1,k)+f08(i+1,j+1,k)+ &
     f09(i-1,j+1,k)+f10(i+1,j-1,k)+ &
     f15(i-1,j,k-1)+f16(i+1,j,k+1)+ &
     f17(i+1,j,k-1)+f18(i-1,j,k+1)
	pxxh(i,j,k,myblock)=udotc
	
	udotc=f03(i,j-1,k)+f04(i,j+1,k)+  &
     f07(i-1,j-1,k)+f08(i+1,j+1,k)+ &
     f09(i-1,j+1,k)+f10(i+1,j-1,k)+ &
     f11(i,j-1,k-1)+f12(i,j+1,k+1)+ &
     f13(i,j-1,k+1)+f14(i,j+1,k-1)
	pyyh(i,j,k,myblock)=udotc
	
	udotc=f05(i,j,k-1)+f06(i,j,k+1)+  &
     f11(i,j-1,k-1)+f12(i,j+1,k+1)+ &
     f13(i,j-1,k+1)+f14(i,j+1,k-1)+ &
     f15(i-1,j,k-1)+f16(i+1,j,k+1)+ &
     f17(i+1,j,k-1)+f18(i-1,j,k+1)
	pzzh(i,j,k,myblock)=udotc
	
	udotc=f07(i-1,j-1,k)+f08(i+1,j+1,k)- &
     f09(i-1,j+1,k)-f10(i+1,j-1,k)
	pxyh(i,j,k,myblock)=udotc
	
	udotc=f15(i-1,j,k-1)+f16(i+1,j,k+1)- &
     f17(i+1,j,k-1)-f18(i-1,j,k+1)
	pxzh(i,j,k,myblock)=udotc
	
	udotc=f11(i,j-1,k-1)+f12(i,j+1,k+1)- &
     f13(i,j-1,k+1)-f14(i,j+1,k-1)
	pyzh(i,j,k,myblock)=udotc
     
#ifdef PRESSCORR

	call syncthreads
	
	uu=halfonecssq*(uh(i,j,k,myblock)*uh(i,j,k,myblock) + vh(i,j,k,myblock)*vh(i,j,k,myblock) + wh(i,j,k,myblock)*wh(i,j,k,myblock))
    
	!1 -1  0  0
	udotc=uh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=vh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=wh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(uh(i,j,k,myblock)+vh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-uh(i,j,k,myblock)+vh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(uh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-uh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(vh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(vh(i,j,k,myblock)-wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(i,j,k)+f02(i,j,k)+  &
     f07(i,j,k)+f08(i,j,k)+ &
     f09(i,j,k)+f10(i,j,k)+ &
     f15(i,j,k)+f16(i,j,k)+ &
     f17(i,j,k)+f18(i,j,k)
	pxxh(i,j,k,myblock)=pxxh(i,j,k,myblock)-udotc
	
	udotc=f03(i,j,k)+f04(i,j,k)+  &
     f07(i,j,k)+f08(i,j,k)+ &
     f09(i,j,k)+f10(i,j,k)+ &
     f11(i,j,k)+f12(i,j,k)+ &
     f13(i,j,k)+f14(i,j,k)
	pyyh(i,j,k,myblock)=pyyh(i,j,k,myblock)-udotc
	
	udotc=f05(i,j,k)+f06(i,j,k)+  &
     f11(i,j,k)+f12(i,j,k)+ &
     f13(i,j,k)+f14(i,j,k)+ &
     f15(i,j,k)+f16(i,j,k)+ &
     f17(i,j,k)+f18(i,j,k)
	pzzh(i,j,k,myblock)=pzzh(i,j,k,myblock)-udotc
	
	udotc=f07(i,j,k)+f08(i,j,k)- &
     f09(i,j,k)-f10(i,j,k)
	pxyh(i,j,k,myblock)=pxyh(i,j,k,myblock)-udotc
	
	udotc=f15(i,j,k)+f16(i,j,k)- &
     f17(i,j,k)-f18(i,j,k)
	pxzh(i,j,k,myblock)=pxzh(i,j,k,myblock)-udotc
	
	udotc=f11(i,j,k)+f12(i,j,k)- &
     f13(i,j,k)-f14(i,j,k)
	pyzh(i,j,k,myblock)=pyzh(i,j,k,myblock)-udotc
	
	    
    return
#endif	
  end subroutine streamcoll_shared
  
  attributes(global) subroutine streamcoll_shared_flop()
	
	implicit none  
	
	integer :: i,j,k
	integer :: gi,gj,gk,myblock,iidblock
    integer :: gii,gjj,gkk,ii,jj,kk
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
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
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
	
	
	
!	 myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
!    zblock=(idblock-1)/nxyblock +1
!    yblock=((idblock-1)-(zblock-1)*nxyblock)/nxblock +1
!    xblock=(idblock-1)-(zblock-1)*nxyblock-(yblock-1)*nxblock +1
	
	!myblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	!iidblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	!if(myblock.ne.iidblock)write(*,*)'cazzo2',gi,gj,gk,myblock,iidblock
    
    uu=halfonecssq*(uh(i,j,k,myblock)*uh(i,j,k,myblock) + vh(i,j,k,myblock)*vh(i,j,k,myblock) + wh(i,j,k,myblock)*wh(i,j,k,myblock))
    
    !0
	f00(i,j,k)=p0*(rhoh(i,j,k,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(pyyh(i,j,k,myblock)+pxxh(i,j,k,myblock)+pzzh(i,j,k,myblock)))
	
    
	!1 -1  0  0
	udotc=uh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k,myblock)-cssq*(pyyh(i,j,k,myblock)+pzzh(i,j,k,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k,myblock)-cssq*(pyyh(i,j,k,myblock)+pzzh(i,j,k,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=vh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k,myblock)-cssq*(pxxh(i,j,k,myblock)+pzzh(i,j,k,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k,myblock)-cssq*(pxxh(i,j,k,myblock)+pzzh(i,j,k,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=wh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k,myblock)-cssq*(pxxh(i,j,k,myblock)+pyyh(i,j,k,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(i,j,k)=p1*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k,myblock)-cssq*(pxxh(i,j,k,myblock)+pyyh(i,j,k,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(uh(i,j,k,myblock)+vh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qyy*pyyh(i,j,k,myblock)-cssq*pzzh(i,j,k,myblock)+two*qxy_7_8*pxyh(i,j,k,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qyy*pyyh(i,j,k,myblock)-cssq*pzzh(i,j,k,myblock)+two*qxy_7_8*pxyh(i,j,k,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-uh(i,j,k,myblock)+vh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qyy*pyyh(i,j,k,myblock)-cssq*pzzh(i,j,k,myblock)+two*qxy_9_10*pxyh(i,j,k,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qyy*pyyh(i,j,k,myblock)-cssq*pzzh(i,j,k,myblock)+two*qxy_9_10*pxyh(i,j,k,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(uh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pyyh(i,j,k,myblock)+two*qxz_15_16*pxzh(i,j,k,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pyyh(i,j,k,myblock)+two*qxz_15_16*pxzh(i,j,k,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-uh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pyyh(i,j,k,myblock)+two*qxz_17_18*pxzh(i,j,k,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pyyh(i,j,k,myblock)+two*qxz_17_18*pxzh(i,j,k,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(vh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pxxh(i,j,k,myblock)+two*qyz_11_12*pyzh(i,j,k,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pxxh(i,j,k,myblock)+two*qyz_11_12*pyzh(i,j,k,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(vh(i,j,k,myblock)-wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pxxh(i,j,k,myblock)+two*qyz_13_14*pyzh(i,j,k,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(i,j,k)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pxxh(i,j,k,myblock)+two*qyz_13_14*pyzh(i,j,k,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    ! Halo Faces
    if(i==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !1 -1  0  0
	  udotc=uh(ii,jj,kk,iidblock)*onecssq
	  f01(i-1,j,k)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxxh(ii,jj,kk,iidblock)-cssq*(pyyh(ii,jj,kk,iidblock)+pzzh(ii,jj,kk,iidblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f07(i-1,j,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f09(i-1,j,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f15(i-1,j,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f18(i-1,j,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(i==TILE_DIMx_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !2 +1  0  0
	  udotc=uh(ii,jj,kk,iidblock)*onecssq
	  f02(i+1,j,k)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxxh(ii,jj,kk,iidblock)-cssq*(pyyh(ii,jj,kk,iidblock)+pzzh(ii,jj,kk,iidblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f08(i+1,j,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f10(i+1,j,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f16(i+1,j,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f17(i+1,j,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(j==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !3 0 -1  0
	  udotc=vh(ii,jj,kk,iidblock)*onecssq
	  f03(i,j-1,k)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyyh(ii,jj,kk,iidblock)-cssq*(pxxh(ii,jj,kk,iidblock)+pzzh(ii,jj,kk,iidblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f07(i,j-1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f10(i,j-1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f11(i,j-1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f13(i,j-1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(j==TILE_DIMy_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !4  0 +1  0
	  udotc=vh(ii,jj,kk,iidblock)*onecssq
	  f04(i,j+1,k)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyyh(ii,jj,kk,iidblock)-cssq*(pxxh(ii,jj,kk,iidblock)+pzzh(ii,jj,kk,iidblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f08(i,j+1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f09(i,j+1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f12(i,j+1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f14(i,j+1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(k==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !5  0  0 -1
	  udotc=wh(ii,jj,kk,iidblock)*onecssq
	  f05(i,j,k-1)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzzh(ii,jj,kk,iidblock)-cssq*(pxxh(ii,jj,kk,iidblock)+pyyh(ii,jj,kk,iidblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f15(i,j,k-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f17(i,j,k-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f11(i,j,k-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f14(i,j,k-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !6  0  0  +1
	  udotc=wh(ii,jj,kk,iidblock)*onecssq
	  f06(i,j,k+1)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzzh(ii,jj,kk,iidblock)-cssq*(pxxh(ii,jj,kk,iidblock)+pyyh(ii,jj,kk,iidblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f16(i,j,k+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f18(i,j,k+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f12(i,j,k+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f13(i,j,k+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(i==1 .and. j==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !7 -1 -1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f07(i-1,j-1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(i==1 .and. j==TILE_DIMy_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !9  -1 +1 0
      udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f09(i-1,j+1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(i==1 .and. k==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !15  -1  0 -1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f15(i-1,j,k-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(i==1 .and. k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !18   -1   0  +1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f18(i-1,j,k+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(i==TILE_DIMx_d .and. j==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !10   +1 -1  0
	  udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f10(i+1,j-1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(i==TILE_DIMx_d .and. j==TILE_DIMy_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !8 +1 +1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f08(i+1,j+1,k)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(i==TILE_DIMx_d .and. k==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !17  +1  0 -1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f17(i+1,j,k-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(i==TILE_DIMx_d .and. k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !16  +1  0 +1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f16(i+1,j,k+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(j==1 .and. k==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !11  0  -1  -1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f11(i,j-1,k-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(j==1 .and. k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !13  0  -1   +1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f13(i,j-1,k+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(j==TILE_DIMy_d .and. k==1)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !14  0  +1  -1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f14(i,j+1,k-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(j==TILE_DIMy_d .and. k==TILE_DIMz_d)then
      
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
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !12  0  +1  +1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f12(i,j+1,k+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
    
    if(abs(isfluid(gi,gj,gk)).ne.1)return
    
    udotc=f00(i,j,k)+f01(i-1,j,k)+f02(i+1,j,k)+  &
     f03(i,j-1,k)+f04(i,j+1,k)+  &
     f05(i,j,k-1)+f06(i,j,k+1)+  &
     f07(i-1,j-1,k)+f08(i+1,j+1,k)+ &
     f09(i-1,j+1,k)+f10(i+1,j-1,k)+ &
     f11(i,j-1,k-1)+f12(i,j+1,k+1)+ &
     f13(i,j-1,k+1)+f14(i,j+1,k-1)+ &
     f15(i-1,j,k-1)+f16(i+1,j,k+1)+ &
     f17(i+1,j,k-1)+f18(i-1,j,k+1)
	rho(i,j,k,myblock)=udotc
	
	udotc=f01(i-1,j,k)-f02(i+1,j,k)+  &
     f07(i-1,j-1,k)-f08(i+1,j+1,k)+ &
     f09(i-1,j+1,k)-f10(i+1,j-1,k)+ &
     f15(i-1,j,k-1)-f16(i+1,j,k+1)- &
     f17(i+1,j,k-1)+f18(i-1,j,k+1)
	u(i,j,k,myblock)=udotc
	
	
	udotc=f03(i,j-1,k)-f04(i,j+1,k)+ &
     f07(i-1,j-1,k)-f08(i+1,j+1,k)- &
     f09(i-1,j+1,k)+f10(i+1,j-1,k)+ &
     f11(i,j-1,k-1)-f12(i,j+1,k+1)+ &
     f13(i,j-1,k+1)-f14(i,j+1,k-1)
	v(i,j,k,myblock)=udotc
	
	udotc=f05(i,j,k-1)-f06(i,j,k+1)+  &
     f11(i,j-1,k-1)-f12(i,j+1,k+1)- &
     f13(i,j-1,k+1)+f14(i,j+1,k-1)+ &
     f15(i-1,j,k-1)-f16(i+1,j,k+1)+ &
     f17(i+1,j,k-1)-f18(i-1,j,k+1)
	w(i,j,k,myblock)=udotc
	
	udotc=f01(i-1,j,k)+f02(i+1,j,k)+  &
     f07(i-1,j-1,k)+f08(i+1,j+1,k)+ &
     f09(i-1,j+1,k)+f10(i+1,j-1,k)+ &
     f15(i-1,j,k-1)+f16(i+1,j,k+1)+ &
     f17(i+1,j,k-1)+f18(i-1,j,k+1)
	pxx(i,j,k,myblock)=udotc
	
	udotc=f03(i,j-1,k)+f04(i,j+1,k)+  &
     f07(i-1,j-1,k)+f08(i+1,j+1,k)+ &
     f09(i-1,j+1,k)+f10(i+1,j-1,k)+ &
     f11(i,j-1,k-1)+f12(i,j+1,k+1)+ &
     f13(i,j-1,k+1)+f14(i,j+1,k-1)
	pyy(i,j,k,myblock)=udotc
	
	udotc=f05(i,j,k-1)+f06(i,j,k+1)+  &
     f11(i,j-1,k-1)+f12(i,j+1,k+1)+ &
     f13(i,j-1,k+1)+f14(i,j+1,k-1)+ &
     f15(i-1,j,k-1)+f16(i+1,j,k+1)+ &
     f17(i+1,j,k-1)+f18(i-1,j,k+1)
	pzz(i,j,k,myblock)=udotc
	
	udotc=f07(i-1,j-1,k)+f08(i+1,j+1,k)- &
     f09(i-1,j+1,k)-f10(i+1,j-1,k)
	pxy(i,j,k,myblock)=udotc
	
	udotc=f15(i-1,j,k-1)+f16(i+1,j,k+1)- &
     f17(i+1,j,k-1)-f18(i-1,j,k+1)
	pxz(i,j,k,myblock)=udotc
	
	udotc=f11(i,j-1,k-1)+f12(i,j+1,k+1)- &
     f13(i,j-1,k+1)-f14(i,j+1,k-1)
	pyz(i,j,k,myblock)=udotc
     
#ifdef PRESSCORR

	call syncthreads
	
	uu=halfonecssq*(u(i,j,k,myblock)*u(i,j,k,myblock) + v(i,j,k,myblock)*v(i,j,k,myblock) + w(i,j,k,myblock)*w(i,j,k,myblock))
    
	!1 -1  0  0
	udotc=u(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(i,j,k)=p1*(rho(i,j,k,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(i,j,k)=p1*(rho(i,j,k,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=v(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(i,j,k)=p1*(rho(i,j,k,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(i,j,k)=p1*(rho(i,j,k,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=w(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(i,j,k)=p1*(rho(i,j,k,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(i,j,k)=p1*(rho(i,j,k,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(u(i,j,k,myblock)+v(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-u(i,j,k,myblock)+v(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(u(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-u(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(v(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(v(i,j,k,myblock)-w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(i,j,k)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(i,j,k)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(i,j,k)+f02(i,j,k)+  &
     f07(i,j,k)+f08(i,j,k)+ &
     f09(i,j,k)+f10(i,j,k)+ &
     f15(i,j,k)+f16(i,j,k)+ &
     f17(i,j,k)+f18(i,j,k)
	pxx(i,j,k,myblock)=pxx(i,j,k,myblock)-udotc
	
	udotc=f03(i,j,k)+f04(i,j,k)+  &
     f07(i,j,k)+f08(i,j,k)+ &
     f09(i,j,k)+f10(i,j,k)+ &
     f11(i,j,k)+f12(i,j,k)+ &
     f13(i,j,k)+f14(i,j,k)
	pyy(i,j,k,myblock)=pyy(i,j,k,myblock)-udotc
	
	udotc=f05(i,j,k)+f06(i,j,k)+  &
     f11(i,j,k)+f12(i,j,k)+ &
     f13(i,j,k)+f14(i,j,k)+ &
     f15(i,j,k)+f16(i,j,k)+ &
     f17(i,j,k)+f18(i,j,k)
	pzz(i,j,k,myblock)=pzz(i,j,k,myblock)-udotc
	
	udotc=f07(i,j,k)+f08(i,j,k)- &
     f09(i,j,k)-f10(i,j,k)
	pxy(i,j,k,myblock)=pxy(i,j,k,myblock)-udotc
	
	udotc=f15(i,j,k)+f16(i,j,k)- &
     f17(i,j,k)-f18(i,j,k)
	pxz(i,j,k,myblock)=pxz(i,j,k,myblock)-udotc
	
	udotc=f11(i,j,k)+f12(i,j,k)- &
     f13(i,j,k)-f14(i,j,k)
	pyz(i,j,k,myblock)=pyz(i,j,k,myblock)-udotc
	
#endif
    
    return
 
  end subroutine streamcoll_shared_flop
  
    attributes(global) subroutine streamcoll_shared_halo()
	
	implicit none   
	
	integer :: i,j,k
	integer :: li,lj,lk
	integer :: gi,gj,gk,myblock
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
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
    
		  
		  
	!if(isfluid(i,j,k).ne.1)return
    
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
    
!    zblock=(idblock-1)/nxyblock +1
!    yblock=((idblock-1)-(zblock-1)*nxyblock)/nxblock +1
!    xblock=(idblock-1)-(zblock-1)*nxyblock-(yblock-1)*nxblock +1
    
    uu=halfonecssq*(u(i,j,k,myblock)*u(i,j,k,myblock) + v(i,j,k,myblock)*v(i,j,k,myblock) + w(i,j,k,myblock)*w(i,j,k,myblock))
    
    !0
	f00(li,lj,lk)=p0*(rho(i,j,k,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(pyy(i,j,k,myblock)+pxx(i,j,k,myblock)+pzz(i,j,k,myblock)))
	
    
	!1 -1  0  0
	udotc=u(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxx(i,j,k,myblock)-cssq*(pyy(i,j,k,myblock)+pzz(i,j,k,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxx(i,j,k,myblock)-cssq*(pyy(i,j,k,myblock)+pzz(i,j,k,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=v(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyy(i,j,k,myblock)-cssq*(pxx(i,j,k,myblock)+pzz(i,j,k,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyy(i,j,k,myblock)-cssq*(pxx(i,j,k,myblock)+pzz(i,j,k,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=w(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k,myblock)-cssq*(pxx(i,j,k,myblock)+pyy(i,j,k,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k,myblock)-cssq*(pxx(i,j,k,myblock)+pyy(i,j,k,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(u(i,j,k,myblock)+v(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qyy*pyy(i,j,k,myblock)-cssq*pzz(i,j,k,myblock)+two*qxy_7_8*pxy(i,j,k,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qyy*pyy(i,j,k,myblock)-cssq*pzz(i,j,k,myblock)+two*qxy_7_8*pxy(i,j,k,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-u(i,j,k,myblock)+v(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qyy*pyy(i,j,k,myblock)-cssq*pzz(i,j,k,myblock)+two*qxy_9_10*pxy(i,j,k,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qyy*pyy(i,j,k,myblock)-cssq*pzz(i,j,k,myblock)+two*qxy_9_10*pxy(i,j,k,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(u(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pyy(i,j,k,myblock)+two*qxz_15_16*pxz(i,j,k,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pyy(i,j,k,myblock)+two*qxz_15_16*pxz(i,j,k,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-u(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pyy(i,j,k,myblock)+two*qxz_17_18*pxz(i,j,k,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pyy(i,j,k,myblock)+two*qxz_17_18*pxz(i,j,k,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(v(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pxx(i,j,k,myblock)+two*qyz_11_12*pyz(i,j,k,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pxx(i,j,k,myblock)+two*qyz_11_12*pyz(i,j,k,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(v(i,j,k,myblock)-w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pxx(i,j,k,myblock)+two*qyz_13_14*pyz(i,j,k,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,myblock)+qzz*pzz(i,j,k,myblock)-cssq*pxx(i,j,k,myblock)+two*qyz_13_14*pyz(i,j,k,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    call syncthreads
    
    if(li<1 .or. lj<1 .or. lk<1)return
    if(li>TILE_DIMx_d .or. lj>TILE_DIMy_d .or. lk>TILE_DIMz_d)return
    
    if(abs(isfluid(gi,gj,gk)).ne.1)return
    
    udotc=f00(li,lj,lk)+f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	rhoh(i,j,k,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	uh(i,j,k,myblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	vh(i,j,k,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	wh(i,j,k,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pxxh(i,j,k,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	pyyh(i,j,k,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pzzh(i,j,k,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	pxyh(i,j,k,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	pxzh(i,j,k,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	pyzh(i,j,k,myblock)=udotc
     
#ifdef PRESSCORR
	
	uu=halfonecssq*(uh(i,j,k,myblock)*uh(i,j,k,myblock) + &
	 vh(i,j,k,myblock)*vh(i,j,k,myblock) + wh(i,j,k,myblock)*wh(i,j,k,myblock))
    
	!1 -1  0  0
	udotc=uh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=vh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=wh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(uh(i,j,k,myblock)+vh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-uh(i,j,k,myblock)+vh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(uh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-uh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(vh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(vh(i,j,k,myblock)-wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pxxh(i,j,k,myblock)=pxxh(i,j,k,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	pyyh(i,j,k,myblock)=pyyh(i,j,k,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pzzh(i,j,k,myblock)=pzzh(i,j,k,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	pxyh(i,j,k,myblock)=pxyh(i,j,k,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	pxzh(i,j,k,myblock)=pxzh(i,j,k,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	pyzh(i,j,k,myblock)=pyzh(i,j,k,myblock)-udotc
	
#endif	
	    
    return
    
  end subroutine streamcoll_shared_halo
  
  attributes(global) subroutine streamcoll_shared_halo_flop()
	
	implicit none  
	
    integer :: i,j,k
    integer :: li,lj,lk
	integer :: gi,gj,gk,myblock
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
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
    
		  
		  
	!if(isfluid(i,j,k).ne.1)return
    
	
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
    
!    zblock=(idblock-1)/nxyblock +1
!    yblock=((idblock-1)-(zblock-1)*nxyblock)/nxblock +1
!    xblock=(idblock-1)-(zblock-1)*nxyblock-(yblock-1)*nxblock +1
    
    uu=halfonecssq*(uh(i,j,k,myblock)*uh(i,j,k,myblock) + vh(i,j,k,myblock)*vh(i,j,k,myblock) + wh(i,j,k,myblock)*wh(i,j,k,myblock))
    
    !0
	f00(li,lj,lk)=p0*(rhoh(i,j,k,myblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(pyyh(i,j,k,myblock)+pxxh(i,j,k,myblock)+pzzh(i,j,k,myblock)))
	
    
	!1 -1  0  0
	udotc=uh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k,myblock)-cssq*(pyyh(i,j,k,myblock)+pzzh(i,j,k,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k,myblock)-cssq*(pyyh(i,j,k,myblock)+pzzh(i,j,k,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=vh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k,myblock)-cssq*(pxxh(i,j,k,myblock)+pzzh(i,j,k,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k,myblock)-cssq*(pxxh(i,j,k,myblock)+pzzh(i,j,k,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=wh(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k,myblock)-cssq*(pxxh(i,j,k,myblock)+pyyh(i,j,k,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k,myblock)-cssq*(pxxh(i,j,k,myblock)+pyyh(i,j,k,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(uh(i,j,k,myblock)+vh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qyy*pyyh(i,j,k,myblock)-cssq*pzzh(i,j,k,myblock)+two*qxy_7_8*pxyh(i,j,k,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qyy*pyyh(i,j,k,myblock)-cssq*pzzh(i,j,k,myblock)+two*qxy_7_8*pxyh(i,j,k,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-uh(i,j,k,myblock)+vh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qyy*pyyh(i,j,k,myblock)-cssq*pzzh(i,j,k,myblock)+two*qxy_9_10*pxyh(i,j,k,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qyy*pyyh(i,j,k,myblock)-cssq*pzzh(i,j,k,myblock)+two*qxy_9_10*pxyh(i,j,k,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(uh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pyyh(i,j,k,myblock)+two*qxz_15_16*pxzh(i,j,k,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pyyh(i,j,k,myblock)+two*qxz_15_16*pxzh(i,j,k,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-uh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pyyh(i,j,k,myblock)+two*qxz_17_18*pxzh(i,j,k,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pyyh(i,j,k,myblock)+two*qxz_17_18*pxzh(i,j,k,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(vh(i,j,k,myblock)+wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pxxh(i,j,k,myblock)+two*qyz_11_12*pyzh(i,j,k,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pxxh(i,j,k,myblock)+two*qyz_11_12*pyzh(i,j,k,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(vh(i,j,k,myblock)-wh(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pxxh(i,j,k,myblock)+two*qyz_13_14*pyzh(i,j,k,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(rhoh(i,j,k,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,myblock)+qzz*pzzh(i,j,k,myblock)-cssq*pxxh(i,j,k,myblock)+two*qyz_13_14*pyzh(i,j,k,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    
    call syncthreads
    
    if(li<1 .or. lj<1 .or. lk<1)return
    if(li>TILE_DIMx_d .or. lj>TILE_DIMy_d .or. lk>TILE_DIMz_d)return
    
    if(abs(isfluid(gi,gj,gk)).ne.1)return
    
    
   udotc=f00(li,lj,lk)+f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	rho(i,j,k,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	u(i,j,k,myblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	v(i,j,k,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	w(i,j,k,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pxx(i,j,k,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	pyy(i,j,k,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pzz(i,j,k,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	pxy(i,j,k,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	pxz(i,j,k,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	pyz(i,j,k,myblock)=udotc
    
#ifdef PRESSCORR
	
	uu=halfonecssq*(u(i,j,k,myblock)*u(i,j,k,myblock) + &
	 v(i,j,k,myblock)*v(i,j,k,myblock) + w(i,j,k,myblock)*w(i,j,k,myblock))
    
	!1 -1  0  0
	udotc=u(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=v(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=w(i,j,k,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(rho(i,j,k,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(u(i,j,k,myblock)+v(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-u(i,j,k,myblock)+v(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(u(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-u(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(v(i,j,k,myblock)+w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(v(i,j,k,myblock)-w(i,j,k,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(rho(i,j,k,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pxx(i,j,k,myblock)=pxx(i,j,k,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	pyy(i,j,k,myblock)=pyy(i,j,k,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pzz(i,j,k,myblock)=pzz(i,j,k,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	pxy(i,j,k,myblock)=pxy(i,j,k,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	pxz(i,j,k,myblock)=pxz(i,j,k,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	pyz(i,j,k,myblock)=pyz(i,j,k,myblock)-udotc
	
	
#endif
    
    return
 
  end subroutine streamcoll_shared_halo_flop
  
 end module streamcoll_kernels
