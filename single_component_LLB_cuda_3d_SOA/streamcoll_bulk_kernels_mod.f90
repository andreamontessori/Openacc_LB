#include "defines.h"
 module streamcoll_bulk_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine streamcoll_bulk()
	
	implicit none  
	  
    integer :: i,j,k
    integer :: idblock=1
	real(kind=db) :: uu,udotc,temp
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
#ifdef USESHARE 
    integer :: li,lj,lk
	real(kind=db), shared :: loc_rho(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
#endif
		  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
		  
	!if(isfluid(i,j,k).ne.1)return
#ifdef USESHARE  	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    loc_rho(li,lj,lk) = rho(i,j,k,idblock)

    ! Halo Faces
    if(li==1)then
      loc_rho(li-1,lj,lk) = rho(i-1,j,k,idblock)
    endif
    if(li==TILE_DIMx_d)then
      loc_rho(li+1,lj,lk) = rho(i+1,j,k,idblock)
    endif

    if(lj==1)then
      loc_rho(li,lj-1,lk) = rho(i,j-1,k,idblock)
    endif
    if(lj==TILE_DIMy_d)then
      loc_rho(li,lj+1,lk) = rho(i,j+1,k,idblock)
    endif

    if(lk==1)then
      loc_rho(li,lj,lk-1) = rho(i,j,k-1,idblock)
    endif
    if(lk==TILE_DIMz_d)then
      loc_rho(li,lj,lk+1) = rho(i,j,k+1,idblock)
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      loc_rho(li-1,lj-1,lk) = rho(i-1,j-1,k,idblock)
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      loc_rho(li-1,lj+1,lk) = rho(i-1,j+1,k,idblock)
    endif
    if(li==1 .and. lk==1)then
      loc_rho(li-1,lj,lk-1) = rho(i-1,j,k-1,idblock)
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      loc_rho(li-1,lj,lk+1) = rho(i-1,j,k+1,idblock)
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      loc_rho(li+1,lj-1,lk) = rho(i+1,j-1,k,idblock)
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      loc_rho(li+1,lj+1,lk) = rho(i+1,j+1,k,idblock)
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      loc_rho(li+1,lj,lk-1) = rho(i+1,j,k-1,idblock)
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      loc_rho(li+1,lj,lk+1) = rho(i+1,j,k+1,idblock)
    endif
    if(lj==1 .and. lk==1)then
      loc_rho(li,lj-1,lk-1) = rho(i,j-1,k-1,idblock)
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      loc_rho(li,lj-1,lk+1) = rho(i,j-1,k+1,idblock)
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      loc_rho(li,lj+1,lk-1) = rho(i,j+1,k-1,idblock)
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      loc_rho(li,lj+1,lk+1) = rho(i,j+1,k+1,idblock)
    endif      

    call syncthreads
#endif
    
    uu=halfonecssq*(u(i,j,k,idblock)*u(i,j,k,idblock) + v(i,j,k,idblock)*v(i,j,k,idblock) + w(i,j,k,idblock)*w(i,j,k,idblock))
	!0
#ifdef USESHARE 
    temp_pop=p0*(loc_rho(li,lj,lk)-uu)
#else
	temp_pop=p0*(rho(i,j,k,idblock)-uu)
#endif
	temp_rho=temp_pop + oneminusomega*pi2cssq0*(-cssq*(pyy(i,j,k,idblock)+pxx(i,j,k,idblock)+pzz(i,j,k,idblock)))
    
	!1
	uu=halfonecssq*(u(i-1,j,k,idblock)*u(i-1,j,k,idblock) + v(i-1,j,k,idblock)*v(i-1,j,k,idblock) + w(i-1,j,k,idblock)*w(i-1,j,k,idblock))
	udotc=u(i-1,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li-1,lj,lk)+(temp + udotc))
#else
	temp_pop=p1*(rho(i-1,j,k,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i-1,j,k,idblock)-cssq*(pyy(i-1,j,k,idblock)+pzz(i-1,j,k,idblock)))
	temp_pop=temp_pop + fx*p1dcssq
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop

	!2
	uu=halfonecssq*(u(i+1,j,k,idblock)*u(i+1,j,k,idblock) + v(i+1,j,k,idblock)*v(i+1,j,k,idblock) + w(i+1,j,k,idblock)*w(i+1,j,k,idblock))
	udotc=u(i+1,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li+1,lj,lk)+(temp - udotc))
#else
	temp_pop=p1*(rho(i+1,j,k,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i+1,j,k,idblock)-cssq*(pyy(i+1,j,k,idblock)+pzz(i+1,j,k,idblock)))
	temp_pop=temp_pop - fx*p1dcssq
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
    		
	!3
	uu=halfonecssq*(u(i,j-1,k,idblock)*u(i,j-1,k,idblock) + v(i,j-1,k,idblock)*v(i,j-1,k,idblock) + w(i,j-1,k,idblock)*w(i,j-1,k,idblock))
	udotc=v(i,j-1,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p1*(rho(i,j-1,k,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j-1,k,idblock)-cssq*(pxx(i,j-1,k,idblock)+pzz(i,j-1,k,idblock)))
	temp_pop=temp_pop + fy*p1dcssq
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	uu=halfonecssq*(u(i,j+1,k,idblock)*u(i,j+1,k,idblock) + v(i,j+1,k,idblock)*v(i,j+1,k,idblock) + w(i,j+1,k,idblock)*w(i,j+1,k,idblock))
	udotc=v(i,j+1,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p1*(rho(i,j+1,k,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j+1,k,idblock)-cssq*(pxx(i,j+1,k,idblock)+pzz(i,j+1,k,idblock)))
	temp_pop=temp_pop - fy*p1dcssq
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	uu=halfonecssq*(u(i,j,k-1,idblock)*u(i,j,k-1,idblock) + v(i,j,k-1,idblock)*v(i,j,k-1,idblock) + w(i,j,k-1,idblock)*w(i,j,k-1,idblock))
	udotc=w(i,j,k-1,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p1*(rho(i,j,k-1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k-1,idblock)-cssq*(pxx(i,j,k-1,idblock)+pyy(i,j,k-1,idblock)))
	temp_pop=temp_pop + fz*p1dcssq
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop

	!6 +1
	uu=halfonecssq*(u(i,j,k+1,idblock)*u(i,j,k+1,idblock) + v(i,j,k+1,idblock)*v(i,j,k+1,idblock) + w(i,j,k+1,idblock)*w(i,j,k+1,idblock))
	udotc=w(i,j,k+1,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p1*(rho(i,j,k+1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k+1,idblock)-cssq*(pxx(i,j,k+1,idblock)+pyy(i,j,k+1,idblock)))
	temp_pop=temp_pop - fz*p1dcssq
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
    	
	!7
	uu=halfonecssq*(u(i-1,j-1,k,idblock)*u(i-1,j-1,k,idblock) + v(i-1,j-1,k,idblock)*v(i-1,j-1,k,idblock) + w(i-1,j-1,k,idblock)*w(i-1,j-1,k,idblock))
	udotc=(u(i-1,j-1,k,idblock)+v(i-1,j-1,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li-1,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p2*(rho(i-1,j-1,k,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j-1,k,idblock)+qyy*pyy(i-1,j-1,k,idblock)-cssq*pzz(i-1,j-1,k,idblock)+two*qxy_7_8*pxy(i-1,j-1,k,idblock))
	temp_pop=temp_pop + (fx+fy)*p2dcssq 
	!+1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pop
	
	!8
	uu=halfonecssq*(u(i+1,j+1,k,idblock)*u(i+1,j+1,k,idblock) + v(i+1,j+1,k,idblock)*v(i+1,j+1,k,idblock) + w(i+1,j+1,k,idblock)*w(i+1,j+1,k,idblock))
	udotc=(u(i+1,j+1,k,idblock)+v(i+1,j+1,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li+1,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p2*(rho(i+1,j+1,k,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j+1,k,idblock)+qyy*pyy(i+1,j+1,k,idblock)-cssq*pzz(i+1,j+1,k,idblock)+two*qxy_7_8*pxy(i+1,j+1,k,idblock))
	temp_pop=temp_pop - (fx+fy)*p2dcssq
	!-1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy+temp_pop
	
	!10   +1 -1
	uu=halfonecssq*(u(i+1,j-1,k,idblock)*u(i+1,j-1,k,idblock) + v(i+1,j-1,k,idblock)*v(i+1,j-1,k,idblock) + w(i+1,j-1,k,idblock)*w(i+1,j-1,k,idblock))
	udotc=(-u(i+1,j-1,k,idblock)+v(i+1,j-1,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li+1,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p2*(rho(i+1,j-1,k,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j-1,k,idblock)+qyy*pyy(i+1,j-1,k,idblock)-cssq*pzz(i+1,j-1,k,idblock)+two*qxy_9_10*pxy(i+1,j-1,k,idblock))
	temp_pop=temp_pop +(fy-fx)*p2dcssq
	!-1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
	
	!9  -1 +1
	uu=halfonecssq*(u(i-1,j+1,k,idblock)*u(i-1,j+1,k,idblock) + v(i-1,j+1,k,idblock)*v(i-1,j+1,k,idblock) + w(i-1,j+1,k,idblock)*w(i-1,j+1,k,idblock))
	udotc=(-u(i-1,j+1,k,idblock)+v(i-1,j+1,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li-1,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p2*(rho(i-1,j+1,k,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j+1,k,idblock)+qyy*pyy(i-1,j+1,k,idblock)-cssq*pzz(i-1,j+1,k,idblock)+two*qxy_9_10*pxy(i-1,j+1,k,idblock))
	temp_pop=temp_pop + (fx-fy)*p2dcssq
	!+1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
		

	!15  -1  -1
	uu=halfonecssq*(u(i-1,j,k-1,idblock)*u(i-1,j,k-1,idblock) + v(i-1,j,k-1,idblock)*v(i-1,j,k-1,idblock) + w(i-1,j,k-1,idblock)*w(i-1,j,k-1,idblock))
	udotc=(u(i-1,j,k-1,idblock)+w(i-1,j,k-1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li-1,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rho(i-1,j,k-1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k-1,idblock)+qzz*pzz(i-1,j,k-1,idblock)-cssq*pyy(i-1,j,k-1,idblock)+two*qxz_15_16*pxz(i-1,j,k-1,idblock))
	temp_pop=temp_pop + (fx+fz)*p2dcssq 

	!+1  0  +1
	
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pop

	!16  +1  +1
	uu=halfonecssq*(u(i+1,j,k+1,idblock)*u(i+1,j,k+1,idblock) + v(i+1,j,k+1,idblock)*v(i+1,j,k+1,idblock) + w(i+1,j,k+1,idblock)*w(i+1,j,k+1,idblock))
	udotc=(u(i+1,j,k+1,idblock)+w(i+1,j,k+1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li+1,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rho(i+1,j,k+1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k+1,idblock)+qzz*pzz(i+1,j,k+1,idblock)-cssq*pyy(i+1,j,k+1,idblock)+two*qxz_15_16*pxz(i+1,j,k+1,idblock))
	temp_pop=temp_pop - (fx+fz)*p2dcssq
	!-1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz+temp_pop

	!17  +1   -1
	uu=halfonecssq*(u(i+1,j,k-1,idblock)*u(i+1,j,k-1,idblock) + v(i+1,j,k-1,idblock)*v(i+1,j,k-1,idblock) + w(i+1,j,k-1,idblock)*w(i+1,j,k-1,idblock))
	udotc=(-u(i+1,j,k-1,idblock)+w(i+1,j,k-1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li+1,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rho(i+1,j,k-1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k-1,idblock)+qzz*pzz(i+1,j,k-1,idblock)-cssq*pyy(i+1,j,k-1,idblock)+two*qxz_17_18*pxz(i+1,j,k-1,idblock))
	temp_pop=temp_pop +(fz-fx)*p2dcssq
	!-1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!18   -1   +1
	uu=halfonecssq*(u(i-1,j,k+1,idblock)*u(i-1,j,k+1,idblock) + v(i-1,j,k+1,idblock)*v(i-1,j,k+1,idblock) + w(i-1,j,k+1,idblock)*w(i-1,j,k+1,idblock))
	udotc=(-u(i-1,j,k+1,idblock)+w(i-1,j,k+1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li-1,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rho(i-1,j,k+1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k+1,idblock)+qzz*pzz(i-1,j,k+1,idblock)-cssq*pyy(i-1,j,k+1,idblock)+two*qxz_17_18*pxz(i-1,j,k+1,idblock))
	temp_pop=temp_pop + (fx-fz)*p2dcssq
	!+1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!11  -1  -1
	uu=halfonecssq*(u(i,j-1,k-1,idblock)*u(i,j-1,k-1,idblock) + v(i,j-1,k-1,idblock)*v(i,j-1,k-1,idblock) + w(i,j-1,k-1,idblock)*w(i,j-1,k-1,idblock))
	udotc=(v(i,j-1,k-1,idblock)+w(i,j-1,k-1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li,lj-1,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rho(i,j-1,k-1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k-1,idblock)+qzz*pzz(i,j-1,k-1,idblock)-cssq*pxx(i,j-1,k-1,idblock)+two*qyz_11_12*pyz(i,j-1,k-1,idblock))
	temp_pop=temp_pop + (fy+fz)*p2dcssq
	! 0 +1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pop
	
	!12   +1  +1
	uu=halfonecssq*(u(i,j+1,k+1,idblock)*u(i,j+1,k+1,idblock) + v(i,j+1,k+1,idblock)*v(i,j+1,k+1,idblock) + w(i,j+1,k+1,idblock)*w(i,j+1,k+1,idblock))
	udotc=(v(i,j+1,k+1,idblock)+w(i,j+1,k+1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li,lj+1,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rho(i,j+1,k+1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k+1,idblock)+qzz*pzz(i,j+1,k+1,idblock)-cssq*pxx(i,j+1,k+1,idblock)+two*qyz_11_12*pyz(i,j+1,k+1,idblock))
	temp_pop=temp_pop - (fy+fz)*p2dcssq
	! 0 -1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz+temp_pop

	!13   -1   +1
	uu=halfonecssq*(u(i,j-1,k+1,idblock)*u(i,j-1,k+1,idblock) + v(i,j-1,k+1,idblock)*v(i,j-1,k+1,idblock) + w(i,j-1,k+1,idblock)*w(i,j-1,k+1,idblock))
	udotc=(v(i,j-1,k+1,idblock)-w(i,j-1,k+1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li,lj-1,lk+1)+(temp + udotc))
#else
	temp_pop=p2*(rho(i,j-1,k+1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k+1,idblock)+qzz*pzz(i,j-1,k+1,idblock)-cssq*pxx(i,j-1,k+1,idblock)+two*qyz_13_14*pyz(i,j-1,k+1,idblock))
	temp_pop=temp_pop + (fy-fz)*p2dcssq
	! 0 +1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	!14  +1 -1
	uu=halfonecssq*(u(i,j+1,k-1,idblock)*u(i,j+1,k-1,idblock) + v(i,j+1,k-1,idblock)*v(i,j+1,k-1,idblock) + w(i,j+1,k-1,idblock)*w(i,j+1,k-1,idblock))
	udotc=(v(i,j+1,k-1,idblock)-w(i,j+1,k-1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li,lj+1,lk-1)+(temp - udotc))
#else
	temp_pop=p2*(rho(i,j+1,k-1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k-1,idblock)+qzz*pzz(i,j+1,k-1,idblock)-cssq*pxx(i,j+1,k-1,idblock)+two*qyz_13_14*pyz(i,j+1,k-1,idblock))
	temp_pop=temp_pop + (fz-fy)*p2dcssq
	! 0 -1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
     	
	rhoh(i,j,k,idblock)=temp_rho
	
	uh(i,j,k,idblock)=temp_u
	vh(i,j,k,idblock)=temp_v
	wh(i,j,k,idblock)=temp_w
	
	pxxh(i,j,k,idblock)=temp_pxx
	pyyh(i,j,k,idblock)=temp_pyy
	pzzh(i,j,k,idblock)=temp_pzz
	pxyh(i,j,k,idblock)=temp_pxy
	pxzh(i,j,k,idblock)=temp_pxz
	pyzh(i,j,k,idblock)=temp_pyz
    
    return
	
  end subroutine streamcoll_bulk
  
  attributes(global) subroutine streamcoll_bulk_flop()
	
	implicit none  
	  
    integer :: i,j,k
    integer :: idblock=1
	real(kind=db) :: uu,udotc,temp
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
#ifdef USESHARE 
    integer :: li,lj,lk
	real(kind=db), shared :: loc_rhoh(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
#endif
		  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
		  
	!if(isfluid(i,j,k).ne.1)return
#ifdef USESHARE	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    loc_rhoh(li,lj,lk) = rhoh(i,j,k,idblock)

    ! Halo Faces
    if(li==1)then
      loc_rhoh(li-1,lj,lk) = rhoh(i-1,j,k,idblock)
    endif
    if(li==TILE_DIMx_d)then
      loc_rhoh(li+1,lj,lk) = rhoh(i+1,j,k,idblock)
    endif

    if(lj==1)then
      loc_rhoh(li,lj-1,lk) = rhoh(i,j-1,k,idblock)
    endif
    if(lj==TILE_DIMy_d)then
      loc_rhoh(li,lj+1,lk) = rhoh(i,j+1,k,idblock)
    endif

    if(lk==1)then
      loc_rhoh(li,lj,lk-1) = rhoh(i,j,k-1,idblock)
    endif
    if(lk==TILE_DIMz_d)then
      loc_rhoh(li,lj,lk+1) = rhoh(i,j,k+1,idblock)
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      loc_rhoh(li-1,lj-1,lk) = rhoh(i-1,j-1,k,idblock)
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      loc_rhoh(li-1,lj+1,lk) = rhoh(i-1,j+1,k,idblock)
    endif
    if(li==1 .and. lk==1)then
      loc_rhoh(li-1,lj,lk-1) = rhoh(i-1,j,k-1,idblock)
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      loc_rhoh(li-1,lj,lk+1) = rhoh(i-1,j,k+1,idblock)
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      loc_rhoh(li+1,lj-1,lk) = rhoh(i+1,j-1,k,idblock)
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      loc_rhoh(li+1,lj+1,lk) = rhoh(i+1,j+1,k,idblock)
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      loc_rhoh(li+1,lj,lk-1) = rhoh(i+1,j,k-1,idblock)
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      loc_rhoh(li+1,lj,lk+1) = rhoh(i+1,j,k+1,idblock)
    endif
    if(lj==1 .and. lk==1)then
      loc_rhoh(li,lj-1,lk-1) = rhoh(i,j-1,k-1,idblock)
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      loc_rhoh(li,lj-1,lk+1) = rhoh(i,j-1,k+1,idblock)
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      loc_rhoh(li,lj+1,lk-1) = rhoh(i,j+1,k-1,idblock)
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      loc_rhoh(li,lj+1,lk+1) = rhoh(i,j+1,k+1,idblock)
    endif      

    call syncthreads
#endif
    
    uu=halfonecssq*(uh(i,j,k,idblock)*uh(i,j,k,idblock) + vh(i,j,k,idblock)*vh(i,j,k,idblock) + wh(i,j,k,idblock)*wh(i,j,k,idblock))
	!0
#ifdef USESHARE 
    temp_pop=p0*(loc_rhoh(li,lj,lk)-uu)
#else
	temp_pop=p0*(rhoh(i,j,k,idblock)-uu)
#endif
	temp_rho=temp_pop + oneminusomega*pi2cssq0*(-cssq*(pyyh(i,j,k,idblock)+pxxh(i,j,k,idblock)+pzzh(i,j,k,idblock)))
    
	!1
	uu=halfonecssq*(uh(i-1,j,k,idblock)*uh(i-1,j,k,idblock) + vh(i-1,j,k,idblock)*vh(i-1,j,k,idblock) + wh(i-1,j,k,idblock)*wh(i-1,j,k,idblock))
	udotc=uh(i-1,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li-1,lj,lk)+(temp + udotc))
#else
	temp_pop=p1*(rhoh(i-1,j,k,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i-1,j,k,idblock)-cssq*(pyyh(i-1,j,k,idblock)+pzzh(i-1,j,k,idblock)))
	temp_pop=temp_pop + fx*p1dcssq
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop

	!2
	uu=halfonecssq*(uh(i+1,j,k,idblock)*uh(i+1,j,k,idblock) + vh(i+1,j,k,idblock)*vh(i+1,j,k,idblock) + wh(i+1,j,k,idblock)*wh(i+1,j,k,idblock))
	udotc=uh(i+1,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li+1,lj,lk)+(temp - udotc))
#else
	temp_pop=p1*(rhoh(i+1,j,k,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i+1,j,k,idblock)-cssq*(pyyh(i+1,j,k,idblock)+pzzh(i+1,j,k,idblock)))
	temp_pop=temp_pop - fx*p1dcssq
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
    		
	!3
	uu=halfonecssq*(uh(i,j-1,k,idblock)*uh(i,j-1,k,idblock) + vh(i,j-1,k,idblock)*vh(i,j-1,k,idblock) + wh(i,j-1,k,idblock)*wh(i,j-1,k,idblock))
	udotc=vh(i,j-1,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p1*(rhoh(i,j-1,k,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j-1,k,idblock)-cssq*(pxxh(i,j-1,k,idblock)+pzzh(i,j-1,k,idblock)))
	temp_pop=temp_pop + fy*p1dcssq
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	uu=halfonecssq*(uh(i,j+1,k,idblock)*uh(i,j+1,k,idblock) + vh(i,j+1,k,idblock)*vh(i,j+1,k,idblock) + wh(i,j+1,k,idblock)*wh(i,j+1,k,idblock))
	udotc=vh(i,j+1,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p1*(rhoh(i,j+1,k,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j+1,k,idblock)-cssq*(pxxh(i,j+1,k,idblock)+pzzh(i,j+1,k,idblock)))
	temp_pop=temp_pop - fy*p1dcssq
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	uu=halfonecssq*(uh(i,j,k-1,idblock)*uh(i,j,k-1,idblock) + vh(i,j,k-1,idblock)*vh(i,j,k-1,idblock) + wh(i,j,k-1,idblock)*wh(i,j,k-1,idblock))
	udotc=wh(i,j,k-1,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p1*(rhoh(i,j,k-1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k-1,idblock)-cssq*(pxxh(i,j,k-1,idblock)+pyyh(i,j,k-1,idblock)))
	temp_pop=temp_pop + fz*p1dcssq
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop
	
	!6 +1
	uu=halfonecssq*(uh(i,j,k+1,idblock)*uh(i,j,k+1,idblock) + vh(i,j,k+1,idblock)*vh(i,j,k+1,idblock) + wh(i,j,k+1,idblock)*wh(i,j,k+1,idblock))
	udotc=wh(i,j,k+1,idblock)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p1*(rhoh(i,j,k+1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k+1,idblock)-cssq*(pxxh(i,j,k+1,idblock)+pyyh(i,j,k+1,idblock)))
	temp_pop=temp_pop - fz*p1dcssq
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
	
	!7
	uu=halfonecssq*(uh(i-1,j-1,k,idblock)*uh(i-1,j-1,k,idblock) + vh(i-1,j-1,k,idblock)*vh(i-1,j-1,k,idblock) + wh(i-1,j-1,k,idblock)*wh(i-1,j-1,k,idblock))
	udotc=(uh(i-1,j-1,k,idblock)+vh(i-1,j-1,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li-1,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i-1,j-1,k,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j-1,k,idblock)+qyy*pyyh(i-1,j-1,k,idblock)-cssq*pzzh(i-1,j-1,k,idblock)+two*qxy_7_8*pxyh(i-1,j-1,k,idblock))
	temp_pop=temp_pop + (fx+fy)*p2dcssq 
	!+1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pop
	
	!8
	uu=halfonecssq*(uh(i+1,j+1,k,idblock)*uh(i+1,j+1,k,idblock) + vh(i+1,j+1,k,idblock)*vh(i+1,j+1,k,idblock) + wh(i+1,j+1,k,idblock)*wh(i+1,j+1,k,idblock))
	udotc=(uh(i+1,j+1,k,idblock)+vh(i+1,j+1,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li+1,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i+1,j+1,k,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j+1,k,idblock)+qyy*pyyh(i+1,j+1,k,idblock)-cssq*pzzh(i+1,j+1,k,idblock)+two*qxy_7_8*pxyh(i+1,j+1,k,idblock))
	temp_pop=temp_pop - (fx+fy)*p2dcssq
	!-1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy+temp_pop
	
	!10   +1 -1
	uu=halfonecssq*(uh(i+1,j-1,k,idblock)*uh(i+1,j-1,k,idblock) + vh(i+1,j-1,k,idblock)*vh(i+1,j-1,k,idblock) + wh(i+1,j-1,k,idblock)*wh(i+1,j-1,k,idblock))
	udotc=(-uh(i+1,j-1,k,idblock)+vh(i+1,j-1,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li+1,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i+1,j-1,k,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j-1,k,idblock)+qyy*pyyh(i+1,j-1,k,idblock)-cssq*pzzh(i+1,j-1,k,idblock)+two*qxy_9_10*pxyh(i+1,j-1,k,idblock))
	temp_pop=temp_pop +(fy-fx)*p2dcssq
	!-1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
	
	!9  -1 +1
	uu=halfonecssq*(uh(i-1,j+1,k,idblock)*uh(i-1,j+1,k,idblock) + vh(i-1,j+1,k,idblock)*vh(i-1,j+1,k,idblock) + wh(i-1,j+1,k,idblock)*wh(i-1,j+1,k,idblock))
	udotc=(-uh(i-1,j+1,k,idblock)+vh(i-1,j+1,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li-1,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i-1,j+1,k,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j+1,k,idblock)+qyy*pyyh(i-1,j+1,k,idblock)-cssq*pzzh(i-1,j+1,k,idblock)+two*qxy_9_10*pxyh(i-1,j+1,k,idblock))
	temp_pop=temp_pop + (fx-fy)*p2dcssq
	!+1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
		

	!15  -1  -1
	uu=halfonecssq*(uh(i-1,j,k-1,idblock)*uh(i-1,j,k-1,idblock) + vh(i-1,j,k-1,idblock)*vh(i-1,j,k-1,idblock) + wh(i-1,j,k-1,idblock)*wh(i-1,j,k-1,idblock))
	udotc=(uh(i-1,j,k-1,idblock)+wh(i-1,j,k-1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li-1,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i-1,j,k-1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k-1,idblock)+qzz*pzzh(i-1,j,k-1,idblock)-cssq*pyyh(i-1,j,k-1,idblock)+two*qxz_15_16*pxzh(i-1,j,k-1,idblock))
	temp_pop=temp_pop + (fx+fz)*p2dcssq 

	!+1  0  +1
	
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pop

	!16  +1  +1
	uu=halfonecssq*(uh(i+1,j,k+1,idblock)*uh(i+1,j,k+1,idblock) + vh(i+1,j,k+1,idblock)*vh(i+1,j,k+1,idblock) + wh(i+1,j,k+1,idblock)*wh(i+1,j,k+1,idblock))
	udotc=(uh(i+1,j,k+1,idblock)+wh(i+1,j,k+1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li+1,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i+1,j,k+1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k+1,idblock)+qzz*pzzh(i+1,j,k+1,idblock)-cssq*pyyh(i+1,j,k+1,idblock)+two*qxz_15_16*pxzh(i+1,j,k+1,idblock))
	temp_pop=temp_pop - (fx+fz)*p2dcssq
	!-1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz+temp_pop

	!17  +1   -1
	uu=halfonecssq*(uh(i+1,j,k-1,idblock)*uh(i+1,j,k-1,idblock) + vh(i+1,j,k-1,idblock)*vh(i+1,j,k-1,idblock) + wh(i+1,j,k-1,idblock)*wh(i+1,j,k-1,idblock))
	udotc=(-uh(i+1,j,k-1,idblock)+wh(i+1,j,k-1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li+1,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i+1,j,k-1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k-1,idblock)+qzz*pzzh(i+1,j,k-1,idblock)-cssq*pyyh(i+1,j,k-1,idblock)+two*qxz_17_18*pxzh(i+1,j,k-1,idblock))
	temp_pop=temp_pop +(fz-fx)*p2dcssq
	!-1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!18   -1   +1
	uu=halfonecssq*(uh(i-1,j,k+1,idblock)*uh(i-1,j,k+1,idblock) + vh(i-1,j,k+1,idblock)*vh(i-1,j,k+1,idblock) + wh(i-1,j,k+1,idblock)*wh(i-1,j,k+1,idblock))
	udotc=(-uh(i-1,j,k+1,idblock)+wh(i-1,j,k+1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li-1,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i-1,j,k+1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k+1,idblock)+qzz*pzzh(i-1,j,k+1,idblock)-cssq*pyyh(i-1,j,k+1,idblock)+two*qxz_17_18*pxzh(i-1,j,k+1,idblock))
	temp_pop=temp_pop + (fx-fz)*p2dcssq
	!+1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!11  -1  -1
	uu=halfonecssq*(uh(i,j-1,k-1,idblock)*uh(i,j-1,k-1,idblock) + vh(i,j-1,k-1,idblock)*vh(i,j-1,k-1,idblock) + wh(i,j-1,k-1,idblock)*wh(i,j-1,k-1,idblock))
	udotc=(vh(i,j-1,k-1,idblock)+wh(i,j-1,k-1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li,lj-1,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i,j-1,k-1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k-1,idblock)+qzz*pzzh(i,j-1,k-1,idblock)-cssq*pxxh(i,j-1,k-1,idblock)+two*qyz_11_12*pyzh(i,j-1,k-1,idblock))
	temp_pop=temp_pop + (fy+fz)*p2dcssq
	! 0 +1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pop
	
	!12   +1  +1
	uu=halfonecssq*(uh(i,j+1,k+1,idblock)*uh(i,j+1,k+1,idblock) + vh(i,j+1,k+1,idblock)*vh(i,j+1,k+1,idblock) + wh(i,j+1,k+1,idblock)*wh(i,j+1,k+1,idblock))
	udotc=(vh(i,j+1,k+1,idblock)+wh(i,j+1,k+1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li,lj+1,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i,j+1,k+1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k+1,idblock)+qzz*pzzh(i,j+1,k+1,idblock)-cssq*pxxh(i,j+1,k+1,idblock)+two*qyz_11_12*pyzh(i,j+1,k+1,idblock))
	temp_pop=temp_pop - (fy+fz)*p2dcssq
	! 0 -1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz+temp_pop

	!13   -1   +1
	uu=halfonecssq*(uh(i,j-1,k+1,idblock)*uh(i,j-1,k+1,idblock) + vh(i,j-1,k+1,idblock)*vh(i,j-1,k+1,idblock) + wh(i,j-1,k+1,idblock)*wh(i,j-1,k+1,idblock))
	udotc=(vh(i,j-1,k+1,idblock)-wh(i,j-1,k+1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li,lj-1,lk+1)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i,j-1,k+1,idblock)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k+1,idblock)+qzz*pzzh(i,j-1,k+1,idblock)-cssq*pxxh(i,j-1,k+1,idblock)+two*qyz_13_14*pyzh(i,j-1,k+1,idblock))
	temp_pop=temp_pop + (fy-fz)*p2dcssq
	! 0 +1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	!14  +1 -1
	uu=halfonecssq*(uh(i,j+1,k-1,idblock)*uh(i,j+1,k-1,idblock) + vh(i,j+1,k-1,idblock)*vh(i,j+1,k-1,idblock) + wh(i,j+1,k-1,idblock)*wh(i,j+1,k-1,idblock))
	udotc=(vh(i,j+1,k-1,idblock)-wh(i,j+1,k-1,idblock))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li,lj+1,lk-1)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i,j+1,k-1,idblock)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k-1,idblock)+qzz*pzzh(i,j+1,k-1,idblock)-cssq*pxxh(i,j+1,k-1,idblock)+two*qyz_13_14*pyzh(i,j+1,k-1,idblock))
	temp_pop=temp_pop + (fz-fy)*p2dcssq
	! 0 -1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
    	
	rho(i,j,k,idblock)=temp_rho
	
	u(i,j,k,idblock)=temp_u
	v(i,j,k,idblock)=temp_v
	w(i,j,k,idblock)=temp_w
	
	pxx(i,j,k,idblock)=temp_pxx
	pyy(i,j,k,idblock)=temp_pyy
	pzz(i,j,k,idblock)=temp_pzz
	pxy(i,j,k,idblock)=temp_pxy
	pxz(i,j,k,idblock)=temp_pxz
	pyz(i,j,k,idblock)=temp_pyz
    
    return
	
  end subroutine streamcoll_bulk_flop
  
  attributes(global) subroutine streamcoll_bulk_shared()
	
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
    	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    
    
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
    
    ! Halo Faces
    if(li==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !1 -1  0  0
	  udotc=u(ii,jj,kk,iidblock)*onecssq
	  f01(li-1,lj,lk)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxx(ii,jj,kk,iidblock)-cssq*(pyy(ii,jj,kk,iidblock)+pzz(ii,jj,kk,iidblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f07(li-1,lj,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f09(li-1,lj,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f15(li-1,lj,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f18(li-1,lj,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !2 +1  0  0
	  udotc=u(ii,jj,kk,iidblock)*onecssq
	  f02(li+1,lj,lk)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxx(ii,jj,kk,iidblock)-cssq*(pyy(ii,jj,kk,iidblock)+pzz(ii,jj,kk,iidblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f08(li+1,lj,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f10(li+1,lj,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f16(li+1,lj,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f17(li+1,lj,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !3 0 -1  0
	  udotc=v(ii,jj,kk,iidblock)*onecssq
	  f03(li,lj-1,lk)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyy(ii,jj,kk,iidblock)-cssq*(pxx(ii,jj,kk,iidblock)+pzz(ii,jj,kk,iidblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f07(li,lj-1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f10(li,lj-1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f11(li,lj-1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f13(li,lj-1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !4  0 +1  0
	  udotc=v(ii,jj,kk,iidblock)*onecssq
	  f04(li,lj+1,lk)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyy(ii,jj,kk,iidblock)-cssq*(pxx(ii,jj,kk,iidblock)+pzz(ii,jj,kk,iidblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f08(li,lj+1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f09(li,lj+1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f12(li,lj+1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f14(li,lj+1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      gii=gi
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !5  0  0 -1
	  udotc=w(ii,jj,kk,iidblock)*onecssq
	  f05(li,lj,lk-1)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzz(ii,jj,kk,iidblock)-cssq*(pxx(ii,jj,kk,iidblock)+pyy(ii,jj,kk,iidblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f15(li,lj,lk-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f17(li,lj,lk-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f11(li,lj,lk-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f14(li,lj,lk-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !6  0  0  +1
	  udotc=w(ii,jj,kk,iidblock)*onecssq
	  f06(li,lj,lk+1)=p1*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzz(ii,jj,kk,iidblock)-cssq*(pxx(ii,jj,kk,iidblock)+pyy(ii,jj,kk,iidblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f16(li,lj,lk+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f18(li,lj,lk+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f12(li,lj,lk+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f13(li,lj,lk+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      gii=gi-1
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !7 -1 -1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f07(li-1,lj-1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      gii=gi-1
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !9  -1 +1 0
      udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f09(li-1,lj+1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !15  -1  0 -1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f15(li-1,lj,lk-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi-1
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !18   -1   0  +1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f18(li-1,lj,lk+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      gii=gi+1
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !10   +1 -1  0
	  udotc=(-u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f10(li+1,lj-1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_9_10*pxy(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      gii=gi+1
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !8 +1 +1  0
	  udotc=(u(ii,jj,kk,iidblock)+v(ii,jj,kk,iidblock))*onecssq
	  f08(li+1,lj+1,lk)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qyy*pyy(ii,jj,kk,iidblock)-cssq*pzz(ii,jj,kk,iidblock)+two*qxy_7_8*pxy(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      gii=gi+1
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !17  +1  0 -1
	  udotc=(-u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f17(li+1,lj,lk-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_17_18*pxz(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !16  +1  0 +1
	  udotc=(u(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f16(li+1,lj,lk+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pyy(ii,jj,kk,iidblock)+two*qxz_15_16*pxz(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !11  0  -1  -1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f11(li,lj-1,lk-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj-1
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !13  0  -1   +1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f13(li,lj-1,lk+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      gii=gi
      gjj=gj+1
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !14  0  +1  -1
	  udotc=(v(ii,jj,kk,iidblock)-w(ii,jj,kk,iidblock))*onecssq
	  f14(li,lj+1,lk-1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_13_14*pyz(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(u(ii,jj,kk,iidblock)*u(ii,jj,kk,iidblock) + v(ii,jj,kk,iidblock)*v(ii,jj,kk,iidblock) + w(ii,jj,kk,iidblock)*w(ii,jj,kk,iidblock))
      
      !12  0  +1  +1
	  udotc=(v(ii,jj,kk,iidblock)+w(ii,jj,kk,iidblock))*onecssq
	  f12(li,lj+1,lk+1)=p2*(rho(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(ii,jj,kk,iidblock)+qzz*pzz(ii,jj,kk,iidblock)-cssq*pxx(ii,jj,kk,iidblock)+two*qyz_11_12*pyz(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
    
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

	call syncthreads
	
	uu=halfonecssq*(uh(i,j,k,myblock)*uh(i,j,k,myblock) + vh(i,j,k,myblock)*vh(i,j,k,myblock) + wh(i,j,k,myblock)*wh(i,j,k,myblock))
    
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
	
	    
    return
#endif	
  end subroutine streamcoll_bulk_shared
  
  attributes(global) subroutine streamcoll_bulk_shared_flop()
	
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
    	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    
    
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
    
    ! Halo Faces
    if(li==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !1 -1  0  0
	  udotc=uh(ii,jj,kk,iidblock)*onecssq
	  f01(li-1,lj,lk)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxxh(ii,jj,kk,iidblock)-cssq*(pyyh(ii,jj,kk,iidblock)+pzzh(ii,jj,kk,iidblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f07(li-1,lj,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f09(li-1,lj,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f15(li-1,lj,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f18(li-1,lj,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !2 +1  0  0
	  udotc=uh(ii,jj,kk,iidblock)*onecssq
	  f02(li+1,lj,lk)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxxh(ii,jj,kk,iidblock)-cssq*(pyyh(ii,jj,kk,iidblock)+pzzh(ii,jj,kk,iidblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f08(li+1,lj,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f10(li+1,lj,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f16(li+1,lj,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f17(li+1,lj,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !3 0 -1  0
	  udotc=vh(ii,jj,kk,iidblock)*onecssq
	  f03(li,lj-1,lk)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyyh(ii,jj,kk,iidblock)-cssq*(pxxh(ii,jj,kk,iidblock)+pzzh(ii,jj,kk,iidblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f07(li,lj-1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f10(li,lj-1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f11(li,lj-1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f13(li,lj-1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !4  0 +1  0
	  udotc=vh(ii,jj,kk,iidblock)*onecssq
	  f04(li,lj+1,lk)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyyh(ii,jj,kk,iidblock)-cssq*(pxxh(ii,jj,kk,iidblock)+pzzh(ii,jj,kk,iidblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f08(li,lj+1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f09(li,lj+1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f12(li,lj+1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f14(li,lj+1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      gii=gi
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !5  0  0 -1
	  udotc=wh(ii,jj,kk,iidblock)*onecssq
	  f05(li,lj,lk-1)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzzh(ii,jj,kk,iidblock)-cssq*(pxxh(ii,jj,kk,iidblock)+pyyh(ii,jj,kk,iidblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f15(li,lj,lk-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f17(li,lj,lk-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f11(li,lj,lk-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f14(li,lj,lk-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !6  0  0  +1
	  udotc=wh(ii,jj,kk,iidblock)*onecssq
	  f06(li,lj,lk+1)=p1*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzzh(ii,jj,kk,iidblock)-cssq*(pxxh(ii,jj,kk,iidblock)+pyyh(ii,jj,kk,iidblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f16(li,lj,lk+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f18(li,lj,lk+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f12(li,lj,lk+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f13(li,lj,lk+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      gii=gi-1
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !7 -1 -1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f07(li-1,lj-1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      gii=gi-1
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !9  -1 +1 0
      udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f09(li-1,lj+1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !15  -1  0 -1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f15(li-1,lj,lk-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi-1
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !18   -1   0  +1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f18(li-1,lj,lk+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      gii=gi+1
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !10   +1 -1  0
	  udotc=(-uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f10(li+1,lj-1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_9_10*pxyh(ii,jj,kk,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      gii=gi+1
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !8 +1 +1  0
	  udotc=(uh(ii,jj,kk,iidblock)+vh(ii,jj,kk,iidblock))*onecssq
	  f08(li+1,lj+1,lk)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qyy*pyyh(ii,jj,kk,iidblock)-cssq*pzzh(ii,jj,kk,iidblock)+two*qxy_7_8*pxyh(ii,jj,kk,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      gii=gi+1
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !17  +1  0 -1
	  udotc=(-uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f17(li+1,lj,lk-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_17_18*pxzh(ii,jj,kk,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !16  +1  0 +1
	  udotc=(uh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f16(li+1,lj,lk+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pyyh(ii,jj,kk,iidblock)+two*qxz_15_16*pxzh(ii,jj,kk,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !11  0  -1  -1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f11(li,lj-1,lk-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj-1
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !13  0  -1   +1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f13(li,lj-1,lk+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      gii=gi
      gjj=gj+1
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !14  0  +1  -1
	  udotc=(vh(ii,jj,kk,iidblock)-wh(ii,jj,kk,iidblock))*onecssq
	  f14(li,lj+1,lk-1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_13_14*pyzh(ii,jj,kk,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(uh(ii,jj,kk,iidblock)*uh(ii,jj,kk,iidblock) + vh(ii,jj,kk,iidblock)*vh(ii,jj,kk,iidblock) + wh(ii,jj,kk,iidblock)*wh(ii,jj,kk,iidblock))
      
      !12  0  +1  +1
	  udotc=(vh(ii,jj,kk,iidblock)+wh(ii,jj,kk,iidblock))*onecssq
	  f12(li,lj+1,lk+1)=p2*(rhoh(ii,jj,kk,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(ii,jj,kk,iidblock)+qzz*pzzh(ii,jj,kk,iidblock)-cssq*pxxh(ii,jj,kk,iidblock)+two*qyz_11_12*pyzh(ii,jj,kk,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
    
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

	call syncthreads
	
	uu=halfonecssq*(u(i,j,k,myblock)*u(i,j,k,myblock) + v(i,j,k,myblock)*v(i,j,k,myblock) + w(i,j,k,myblock)*w(i,j,k,myblock))
    
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
 
  end subroutine streamcoll_bulk_shared_flop
  
    attributes(global) subroutine streamcoll_bulk_shared_halo()
	
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
    
		  
		  
	!if(isfluid(i,j,k).ne.1)return
    	
	li = threadIdx%x-1
    lj = threadIdx%y-1
    lk = threadIdx%z-1
    
    gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x-1
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y-1
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z-1
	
	i=threadIdx%x-1
	j=threadIdx%y-1
	k=threadIdx%z-1
	
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    
    
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

	!call syncthreads
	
	uu=halfonecssq*(uh(i,j,k,myblock)*uh(i,j,k,myblock) + vh(i,j,k,myblock)*vh(i,j,k,myblock) + wh(i,j,k,myblock)*wh(i,j,k,myblock))
    
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
	
	    
    return
#endif	
  end subroutine streamcoll_bulk_shared_halo
  
  attributes(global) subroutine streamcoll_bulk_shared_halo_flop()
	
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
    
		  
		  
	!if(isfluid(i,j,k).ne.1)return
    	
	li = threadIdx%x-1
    lj = threadIdx%y-1
    lk = threadIdx%z-1
    
    gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x-1
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y-1
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z-1
	
	i=threadIdx%x-1
	j=threadIdx%y-1
	k=threadIdx%z-1
	
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    
    
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

	!call syncthreads
	
	uu=halfonecssq*(u(i,j,k,myblock)*u(i,j,k,myblock) + v(i,j,k,myblock)*v(i,j,k,myblock) + w(i,j,k,myblock)*w(i,j,k,myblock))
    
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
 
  end subroutine streamcoll_bulk_shared_halo_flop
  
 end module streamcoll_bulk_kernels
