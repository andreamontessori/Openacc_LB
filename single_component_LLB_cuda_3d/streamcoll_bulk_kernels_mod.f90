#include "defines.h"
 module streamcoll_bulk_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine streamcoll_bulk()
	
	implicit none  
	  
    integer :: i,j,k
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
    
    loc_rho(li,lj,lk) = rho(i,j,k)

    ! Halo Faces
    if(li==1)then
      loc_rho(li-1,lj,lk) = rho(i-1,j,k)
    endif
    if(li==TILE_DIMx_d)then
      loc_rho(li+1,lj,lk) = rho(i+1,j,k)
    endif

    if(lj==1)then
      loc_rho(li,lj-1,lk) = rho(i,j-1,k)
    endif
    if(lj==TILE_DIMy_d)then
      loc_rho(li,lj+1,lk) = rho(i,j+1,k)
    endif

    if(lk==1)then
      loc_rho(li,lj,lk-1) = rho(i,j,k-1)
    endif
    if(lk==TILE_DIMz_d)then
      loc_rho(li,lj,lk+1) = rho(i,j,k+1)
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      loc_rho(li-1,lj-1,lk) = rho(i-1,j-1,k)
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      loc_rho(li-1,lj+1,lk) = rho(i-1,j+1,k)
    endif
    if(li==1 .and. lk==1)then
      loc_rho(li-1,lj,lk-1) = rho(i-1,j,k-1)
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      loc_rho(li-1,lj,lk+1) = rho(i-1,j,k+1)
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      loc_rho(li+1,lj-1,lk) = rho(i+1,j-1,k)
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      loc_rho(li+1,lj+1,lk) = rho(i+1,j+1,k)
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      loc_rho(li+1,lj,lk-1) = rho(i+1,j,k-1)
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      loc_rho(li+1,lj,lk+1) = rho(i+1,j,k+1)
    endif
    if(lj==1 .and. lk==1)then
      loc_rho(li,lj-1,lk-1) = rho(i,j-1,k-1)
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      loc_rho(li,lj-1,lk+1) = rho(i,j-1,k+1)
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      loc_rho(li,lj+1,lk-1) = rho(i,j+1,k-1)
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      loc_rho(li,lj+1,lk+1) = rho(i,j+1,k+1)
    endif      

    call syncthreads
#endif
    
    uu=halfonecssq*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))
	!0
#ifdef USESHARE 
    temp_pop=p0*(loc_rho(li,lj,lk)-uu)
#else
	temp_pop=p0*(rho(i,j,k)-uu)
#endif
	temp_rho=temp_pop + oneminusomega*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k)))
    
	!1
	uu=halfonecssq*(u(i-1,j,k)*u(i-1,j,k) + v(i-1,j,k)*v(i-1,j,k) + w(i-1,j,k)*w(i-1,j,k))
	udotc=u(i-1,j,k)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li-1,lj,lk)+(temp + udotc))
#else
	temp_pop=p1*(rho(i-1,j,k)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i-1,j,k)-cssq*(pyy(i-1,j,k)+pzz(i-1,j,k)))
	temp_pop=temp_pop + fx*p1dcssq
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop

	!2
	uu=halfonecssq*(u(i+1,j,k)*u(i+1,j,k) + v(i+1,j,k)*v(i+1,j,k) + w(i+1,j,k)*w(i+1,j,k))
	udotc=u(i+1,j,k)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li+1,lj,lk)+(temp - udotc))
#else
	temp_pop=p1*(rho(i+1,j,k)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i+1,j,k)-cssq*(pyy(i+1,j,k)+pzz(i+1,j,k)))
	temp_pop=temp_pop - fx*p1dcssq
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
    		
	!3
	uu=halfonecssq*(u(i,j-1,k)*u(i,j-1,k) + v(i,j-1,k)*v(i,j-1,k) + w(i,j-1,k)*w(i,j-1,k))
	udotc=v(i,j-1,k)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p1*(rho(i,j-1,k)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j-1,k)-cssq*(pxx(i,j-1,k)+pzz(i,j-1,k)))
	temp_pop=temp_pop + fy*p1dcssq
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	uu=halfonecssq*(u(i,j+1,k)*u(i,j+1,k) + v(i,j+1,k)*v(i,j+1,k) + w(i,j+1,k)*w(i,j+1,k))
	udotc=v(i,j+1,k)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p1*(rho(i,j+1,k)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j+1,k)-cssq*(pxx(i,j+1,k)+pzz(i,j+1,k)))
	temp_pop=temp_pop - fy*p1dcssq
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	uu=halfonecssq*(u(i,j,k-1)*u(i,j,k-1) + v(i,j,k-1)*v(i,j,k-1) + w(i,j,k-1)*w(i,j,k-1))
	udotc=w(i,j,k-1)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p1*(rho(i,j,k-1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k-1)-cssq*(pxx(i,j,k-1)+pyy(i,j,k-1)))
	temp_pop=temp_pop + fz*p1dcssq
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop

	!6 +1
	uu=halfonecssq*(u(i,j,k+1)*u(i,j,k+1) + v(i,j,k+1)*v(i,j,k+1) + w(i,j,k+1)*w(i,j,k+1))
	udotc=w(i,j,k+1)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rho(li,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p1*(rho(i,j,k+1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k+1)-cssq*(pxx(i,j,k+1)+pyy(i,j,k+1)))
	temp_pop=temp_pop - fz*p1dcssq
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
    	
	!7
	uu=halfonecssq*(u(i-1,j-1,k)*u(i-1,j-1,k) + v(i-1,j-1,k)*v(i-1,j-1,k) + w(i-1,j-1,k)*w(i-1,j-1,k))
	udotc=(u(i-1,j-1,k)+v(i-1,j-1,k))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li-1,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p2*(rho(i-1,j-1,k)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j-1,k)+qyy*pyy(i-1,j-1,k)-cssq*pzz(i-1,j-1,k)+two*qxy_7_8*pxy(i-1,j-1,k))
	temp_pop=temp_pop + (fx+fy)*p2dcssq 
	!+1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pop
	
	!8
	uu=halfonecssq*(u(i+1,j+1,k)*u(i+1,j+1,k) + v(i+1,j+1,k)*v(i+1,j+1,k) + w(i+1,j+1,k)*w(i+1,j+1,k))
	udotc=(u(i+1,j+1,k)+v(i+1,j+1,k))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li+1,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p2*(rho(i+1,j+1,k)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j+1,k)+qyy*pyy(i+1,j+1,k)-cssq*pzz(i+1,j+1,k)+two*qxy_7_8*pxy(i+1,j+1,k))
	temp_pop=temp_pop - (fx+fy)*p2dcssq
	!-1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy+temp_pop
	
	!10   +1 -1
	uu=halfonecssq*(u(i+1,j-1,k)*u(i+1,j-1,k) + v(i+1,j-1,k)*v(i+1,j-1,k) + w(i+1,j-1,k)*w(i+1,j-1,k))
	udotc=(-u(i+1,j-1,k)+v(i+1,j-1,k))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li+1,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p2*(rho(i+1,j-1,k)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j-1,k)+qyy*pyy(i+1,j-1,k)-cssq*pzz(i+1,j-1,k)+two*qxy_9_10*pxy(i+1,j-1,k))
	temp_pop=temp_pop +(fy-fx)*p2dcssq
	!-1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
	
	!9  -1 +1
	uu=halfonecssq*(u(i-1,j+1,k)*u(i-1,j+1,k) + v(i-1,j+1,k)*v(i-1,j+1,k) + w(i-1,j+1,k)*w(i-1,j+1,k))
	udotc=(-u(i-1,j+1,k)+v(i-1,j+1,k))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li-1,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p2*(rho(i-1,j+1,k)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j+1,k)+qyy*pyy(i-1,j+1,k)-cssq*pzz(i-1,j+1,k)+two*qxy_9_10*pxy(i-1,j+1,k))
	temp_pop=temp_pop + (fx-fy)*p2dcssq
	!+1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
		

	!15  -1  -1
	uu=halfonecssq*(u(i-1,j,k-1)*u(i-1,j,k-1) + v(i-1,j,k-1)*v(i-1,j,k-1) + w(i-1,j,k-1)*w(i-1,j,k-1))
	udotc=(u(i-1,j,k-1)+w(i-1,j,k-1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li-1,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rho(i-1,j,k-1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k-1)+qzz*pzz(i-1,j,k-1)-cssq*pyy(i-1,j,k-1)+two*qxz_15_16*pxz(i-1,j,k-1))
	temp_pop=temp_pop + (fx+fz)*p2dcssq 

	!+1  0  +1
	
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pop

	!16  +1  +1
	uu=halfonecssq*(u(i+1,j,k+1)*u(i+1,j,k+1) + v(i+1,j,k+1)*v(i+1,j,k+1) + w(i+1,j,k+1)*w(i+1,j,k+1))
	udotc=(u(i+1,j,k+1)+w(i+1,j,k+1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li+1,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rho(i+1,j,k+1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k+1)+qzz*pzz(i+1,j,k+1)-cssq*pyy(i+1,j,k+1)+two*qxz_15_16*pxz(i+1,j,k+1))
	temp_pop=temp_pop - (fx+fz)*p2dcssq
	!-1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz+temp_pop

	!17  +1   -1
	uu=halfonecssq*(u(i+1,j,k-1)*u(i+1,j,k-1) + v(i+1,j,k-1)*v(i+1,j,k-1) + w(i+1,j,k-1)*w(i+1,j,k-1))
	udotc=(-u(i+1,j,k-1)+w(i+1,j,k-1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li+1,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rho(i+1,j,k-1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k-1)+qzz*pzz(i+1,j,k-1)-cssq*pyy(i+1,j,k-1)+two*qxz_17_18*pxz(i+1,j,k-1))
	temp_pop=temp_pop +(fz-fx)*p2dcssq
	!-1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!18   -1   +1
	uu=halfonecssq*(u(i-1,j,k+1)*u(i-1,j,k+1) + v(i-1,j,k+1)*v(i-1,j,k+1) + w(i-1,j,k+1)*w(i-1,j,k+1))
	udotc=(-u(i-1,j,k+1)+w(i-1,j,k+1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li-1,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rho(i-1,j,k+1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k+1)+qzz*pzz(i-1,j,k+1)-cssq*pyy(i-1,j,k+1)+two*qxz_17_18*pxz(i-1,j,k+1))
	temp_pop=temp_pop + (fx-fz)*p2dcssq
	!+1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!11  -1  -1
	uu=halfonecssq*(u(i,j-1,k-1)*u(i,j-1,k-1) + v(i,j-1,k-1)*v(i,j-1,k-1) + w(i,j-1,k-1)*w(i,j-1,k-1))
	udotc=(v(i,j-1,k-1)+w(i,j-1,k-1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li,lj-1,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rho(i,j-1,k-1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k-1)+qzz*pzz(i,j-1,k-1)-cssq*pxx(i,j-1,k-1)+two*qyz_11_12*pyz(i,j-1,k-1))
	temp_pop=temp_pop + (fy+fz)*p2dcssq
	! 0 +1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pop
	
	!12   +1  +1
	uu=halfonecssq*(u(i,j+1,k+1)*u(i,j+1,k+1) + v(i,j+1,k+1)*v(i,j+1,k+1) + w(i,j+1,k+1)*w(i,j+1,k+1))
	udotc=(v(i,j+1,k+1)+w(i,j+1,k+1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li,lj+1,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rho(i,j+1,k+1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k+1)+qzz*pzz(i,j+1,k+1)-cssq*pxx(i,j+1,k+1)+two*qyz_11_12*pyz(i,j+1,k+1))
	temp_pop=temp_pop - (fy+fz)*p2dcssq
	! 0 -1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz+temp_pop

	!13   -1   +1
	uu=halfonecssq*(u(i,j-1,k+1)*u(i,j-1,k+1) + v(i,j-1,k+1)*v(i,j-1,k+1) + w(i,j-1,k+1)*w(i,j-1,k+1))
	udotc=(v(i,j-1,k+1)-w(i,j-1,k+1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li,lj-1,lk+1)+(temp + udotc))
#else
	temp_pop=p2*(rho(i,j-1,k+1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k+1)+qzz*pzz(i,j-1,k+1)-cssq*pxx(i,j-1,k+1)+two*qyz_13_14*pyz(i,j-1,k+1))
	temp_pop=temp_pop + (fy-fz)*p2dcssq
	! 0 +1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	!14  +1 -1
	uu=halfonecssq*(u(i,j+1,k-1)*u(i,j+1,k-1) + v(i,j+1,k-1)*v(i,j+1,k-1) + w(i,j+1,k-1)*w(i,j+1,k-1))
	udotc=(v(i,j+1,k-1)-w(i,j+1,k-1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rho(li,lj+1,lk-1)+(temp - udotc))
#else
	temp_pop=p2*(rho(i,j+1,k-1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k-1)+qzz*pzz(i,j+1,k-1)-cssq*pxx(i,j+1,k-1)+two*qyz_13_14*pyz(i,j+1,k-1))
	temp_pop=temp_pop + (fz-fy)*p2dcssq
	! 0 -1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
     	
	rhoh(i,j,k)=temp_rho
	
	uh(i,j,k)=temp_u
	vh(i,j,k)=temp_v
	wh(i,j,k)=temp_w
	
	pxxh(i,j,k)=temp_pxx
	pyyh(i,j,k)=temp_pyy
	pzzh(i,j,k)=temp_pzz
	pxyh(i,j,k)=temp_pxy
	pxzh(i,j,k)=temp_pxz
	pyzh(i,j,k)=temp_pyz
    
    return
	
  end subroutine streamcoll_bulk
  
  attributes(global) subroutine streamcoll_bulk_flop()
	
	implicit none  
	  
    integer :: i,j,k
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
    
    loc_rhoh(li,lj,lk) = rhoh(i,j,k)

    ! Halo Faces
    if(li==1)then
      loc_rhoh(li-1,lj,lk) = rhoh(i-1,j,k)
    endif
    if(li==TILE_DIMx_d)then
      loc_rhoh(li+1,lj,lk) = rhoh(i+1,j,k)
    endif

    if(lj==1)then
      loc_rhoh(li,lj-1,lk) = rhoh(i,j-1,k)
    endif
    if(lj==TILE_DIMy_d)then
      loc_rhoh(li,lj+1,lk) = rhoh(i,j+1,k)
    endif

    if(lk==1)then
      loc_rhoh(li,lj,lk-1) = rhoh(i,j,k-1)
    endif
    if(lk==TILE_DIMz_d)then
      loc_rhoh(li,lj,lk+1) = rhoh(i,j,k+1)
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      loc_rhoh(li-1,lj-1,lk) = rhoh(i-1,j-1,k)
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      loc_rhoh(li-1,lj+1,lk) = rhoh(i-1,j+1,k)
    endif
    if(li==1 .and. lk==1)then
      loc_rhoh(li-1,lj,lk-1) = rhoh(i-1,j,k-1)
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      loc_rhoh(li-1,lj,lk+1) = rhoh(i-1,j,k+1)
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      loc_rhoh(li+1,lj-1,lk) = rhoh(i+1,j-1,k)
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      loc_rhoh(li+1,lj+1,lk) = rhoh(i+1,j+1,k)
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      loc_rhoh(li+1,lj,lk-1) = rhoh(i+1,j,k-1)
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      loc_rhoh(li+1,lj,lk+1) = rhoh(i+1,j,k+1)
    endif
    if(lj==1 .and. lk==1)then
      loc_rhoh(li,lj-1,lk-1) = rhoh(i,j-1,k-1)
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      loc_rhoh(li,lj-1,lk+1) = rhoh(i,j-1,k+1)
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      loc_rhoh(li,lj+1,lk-1) = rhoh(i,j+1,k-1)
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      loc_rhoh(li,lj+1,lk+1) = rhoh(i,j+1,k+1)
    endif      

    call syncthreads
#endif
    
    uu=halfonecssq*(uh(i,j,k)*uh(i,j,k) + vh(i,j,k)*vh(i,j,k) + wh(i,j,k)*wh(i,j,k))
	!0
#ifdef USESHARE 
    temp_pop=p0*(loc_rhoh(li,lj,lk)-uu)
#else
	temp_pop=p0*(rhoh(i,j,k)-uu)
#endif
	temp_rho=temp_pop + oneminusomega*pi2cssq0*(-cssq*(pyyh(i,j,k)+pxxh(i,j,k)+pzzh(i,j,k)))
    
	!1
	uu=halfonecssq*(uh(i-1,j,k)*uh(i-1,j,k) + vh(i-1,j,k)*vh(i-1,j,k) + wh(i-1,j,k)*wh(i-1,j,k))
	udotc=uh(i-1,j,k)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li-1,lj,lk)+(temp + udotc))
#else
	temp_pop=p1*(rhoh(i-1,j,k)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i-1,j,k)-cssq*(pyyh(i-1,j,k)+pzzh(i-1,j,k)))
	temp_pop=temp_pop + fx*p1dcssq
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop

	!2
	uu=halfonecssq*(uh(i+1,j,k)*uh(i+1,j,k) + vh(i+1,j,k)*vh(i+1,j,k) + wh(i+1,j,k)*wh(i+1,j,k))
	udotc=uh(i+1,j,k)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li+1,lj,lk)+(temp - udotc))
#else
	temp_pop=p1*(rhoh(i+1,j,k)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i+1,j,k)-cssq*(pyyh(i+1,j,k)+pzzh(i+1,j,k)))
	temp_pop=temp_pop - fx*p1dcssq
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
    		
	!3
	uu=halfonecssq*(uh(i,j-1,k)*uh(i,j-1,k) + vh(i,j-1,k)*vh(i,j-1,k) + wh(i,j-1,k)*wh(i,j-1,k))
	udotc=vh(i,j-1,k)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p1*(rhoh(i,j-1,k)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j-1,k)-cssq*(pxxh(i,j-1,k)+pzzh(i,j-1,k)))
	temp_pop=temp_pop + fy*p1dcssq
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	uu=halfonecssq*(uh(i,j+1,k)*uh(i,j+1,k) + vh(i,j+1,k)*vh(i,j+1,k) + wh(i,j+1,k)*wh(i,j+1,k))
	udotc=vh(i,j+1,k)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p1*(rhoh(i,j+1,k)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j+1,k)-cssq*(pxxh(i,j+1,k)+pzzh(i,j+1,k)))
	temp_pop=temp_pop - fy*p1dcssq
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	uu=halfonecssq*(uh(i,j,k-1)*uh(i,j,k-1) + vh(i,j,k-1)*vh(i,j,k-1) + wh(i,j,k-1)*wh(i,j,k-1))
	udotc=wh(i,j,k-1)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p1*(rhoh(i,j,k-1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k-1)-cssq*(pxxh(i,j,k-1)+pyyh(i,j,k-1)))
	temp_pop=temp_pop + fz*p1dcssq
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop
	
	!6 +1
	uu=halfonecssq*(uh(i,j,k+1)*uh(i,j,k+1) + vh(i,j,k+1)*vh(i,j,k+1) + wh(i,j,k+1)*wh(i,j,k+1))
	udotc=wh(i,j,k+1)*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p1*(loc_rhoh(li,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p1*(rhoh(i,j,k+1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k+1)-cssq*(pxxh(i,j,k+1)+pyyh(i,j,k+1)))
	temp_pop=temp_pop - fz*p1dcssq
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
	
	!7
	uu=halfonecssq*(uh(i-1,j-1,k)*uh(i-1,j-1,k) + vh(i-1,j-1,k)*vh(i-1,j-1,k) + wh(i-1,j-1,k)*wh(i-1,j-1,k))
	udotc=(uh(i-1,j-1,k)+vh(i-1,j-1,k))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li-1,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i-1,j-1,k)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j-1,k)+qyy*pyyh(i-1,j-1,k)-cssq*pzzh(i-1,j-1,k)+two*qxy_7_8*pxyh(i-1,j-1,k))
	temp_pop=temp_pop + (fx+fy)*p2dcssq 
	!+1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pop
	
	!8
	uu=halfonecssq*(uh(i+1,j+1,k)*uh(i+1,j+1,k) + vh(i+1,j+1,k)*vh(i+1,j+1,k) + wh(i+1,j+1,k)*wh(i+1,j+1,k))
	udotc=(uh(i+1,j+1,k)+vh(i+1,j+1,k))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li+1,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i+1,j+1,k)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j+1,k)+qyy*pyyh(i+1,j+1,k)-cssq*pzzh(i+1,j+1,k)+two*qxy_7_8*pxyh(i+1,j+1,k))
	temp_pop=temp_pop - (fx+fy)*p2dcssq
	!-1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy+temp_pop
	
	!10   +1 -1
	uu=halfonecssq*(uh(i+1,j-1,k)*uh(i+1,j-1,k) + vh(i+1,j-1,k)*vh(i+1,j-1,k) + wh(i+1,j-1,k)*wh(i+1,j-1,k))
	udotc=(-uh(i+1,j-1,k)+vh(i+1,j-1,k))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li+1,lj-1,lk)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i+1,j-1,k)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j-1,k)+qyy*pyyh(i+1,j-1,k)-cssq*pzzh(i+1,j-1,k)+two*qxy_9_10*pxyh(i+1,j-1,k))
	temp_pop=temp_pop +(fy-fx)*p2dcssq
	!-1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
	
	!9  -1 +1
	uu=halfonecssq*(uh(i-1,j+1,k)*uh(i-1,j+1,k) + vh(i-1,j+1,k)*vh(i-1,j+1,k) + wh(i-1,j+1,k)*wh(i-1,j+1,k))
	udotc=(-uh(i-1,j+1,k)+vh(i-1,j+1,k))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li-1,lj+1,lk)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i-1,j+1,k)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j+1,k)+qyy*pyyh(i-1,j+1,k)-cssq*pzzh(i-1,j+1,k)+two*qxy_9_10*pxyh(i-1,j+1,k))
	temp_pop=temp_pop + (fx-fy)*p2dcssq
	!+1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
		

	!15  -1  -1
	uu=halfonecssq*(uh(i-1,j,k-1)*uh(i-1,j,k-1) + vh(i-1,j,k-1)*vh(i-1,j,k-1) + wh(i-1,j,k-1)*wh(i-1,j,k-1))
	udotc=(uh(i-1,j,k-1)+wh(i-1,j,k-1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li-1,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i-1,j,k-1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k-1)+qzz*pzzh(i-1,j,k-1)-cssq*pyyh(i-1,j,k-1)+two*qxz_15_16*pxzh(i-1,j,k-1))
	temp_pop=temp_pop + (fx+fz)*p2dcssq 

	!+1  0  +1
	
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pop

	!16  +1  +1
	uu=halfonecssq*(uh(i+1,j,k+1)*uh(i+1,j,k+1) + vh(i+1,j,k+1)*vh(i+1,j,k+1) + wh(i+1,j,k+1)*wh(i+1,j,k+1))
	udotc=(uh(i+1,j,k+1)+wh(i+1,j,k+1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li+1,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i+1,j,k+1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k+1)+qzz*pzzh(i+1,j,k+1)-cssq*pyyh(i+1,j,k+1)+two*qxz_15_16*pxzh(i+1,j,k+1))
	temp_pop=temp_pop - (fx+fz)*p2dcssq
	!-1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz+temp_pop

	!17  +1   -1
	uu=halfonecssq*(uh(i+1,j,k-1)*uh(i+1,j,k-1) + vh(i+1,j,k-1)*vh(i+1,j,k-1) + wh(i+1,j,k-1)*wh(i+1,j,k-1))
	udotc=(-uh(i+1,j,k-1)+wh(i+1,j,k-1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li+1,lj,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i+1,j,k-1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k-1)+qzz*pzzh(i+1,j,k-1)-cssq*pyyh(i+1,j,k-1)+two*qxz_17_18*pxzh(i+1,j,k-1))
	temp_pop=temp_pop +(fz-fx)*p2dcssq
	!-1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!18   -1   +1
	uu=halfonecssq*(uh(i-1,j,k+1)*uh(i-1,j,k+1) + vh(i-1,j,k+1)*vh(i-1,j,k+1) + wh(i-1,j,k+1)*wh(i-1,j,k+1))
	udotc=(-uh(i-1,j,k+1)+wh(i-1,j,k+1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li-1,lj,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i-1,j,k+1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k+1)+qzz*pzzh(i-1,j,k+1)-cssq*pyyh(i-1,j,k+1)+two*qxz_17_18*pxzh(i-1,j,k+1))
	temp_pop=temp_pop + (fx-fz)*p2dcssq
	!+1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!11  -1  -1
	uu=halfonecssq*(uh(i,j-1,k-1)*uh(i,j-1,k-1) + vh(i,j-1,k-1)*vh(i,j-1,k-1) + wh(i,j-1,k-1)*wh(i,j-1,k-1))
	udotc=(vh(i,j-1,k-1)+wh(i,j-1,k-1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li,lj-1,lk-1)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i,j-1,k-1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k-1)+qzz*pzzh(i,j-1,k-1)-cssq*pxxh(i,j-1,k-1)+two*qyz_11_12*pyzh(i,j-1,k-1))
	temp_pop=temp_pop + (fy+fz)*p2dcssq
	! 0 +1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pop
	
	!12   +1  +1
	uu=halfonecssq*(uh(i,j+1,k+1)*uh(i,j+1,k+1) + vh(i,j+1,k+1)*vh(i,j+1,k+1) + wh(i,j+1,k+1)*wh(i,j+1,k+1))
	udotc=(vh(i,j+1,k+1)+wh(i,j+1,k+1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li,lj+1,lk+1)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i,j+1,k+1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k+1)+qzz*pzzh(i,j+1,k+1)-cssq*pxxh(i,j+1,k+1)+two*qyz_11_12*pyzh(i,j+1,k+1))
	temp_pop=temp_pop - (fy+fz)*p2dcssq
	! 0 -1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz+temp_pop

	!13   -1   +1
	uu=halfonecssq*(uh(i,j-1,k+1)*uh(i,j-1,k+1) + vh(i,j-1,k+1)*vh(i,j-1,k+1) + wh(i,j-1,k+1)*wh(i,j-1,k+1))
	udotc=(vh(i,j-1,k+1)-wh(i,j-1,k+1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li,lj-1,lk+1)+(temp + udotc))
#else
	temp_pop=p2*(rhoh(i,j-1,k+1)+(temp + udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k+1)+qzz*pzzh(i,j-1,k+1)-cssq*pxxh(i,j-1,k+1)+two*qyz_13_14*pyzh(i,j-1,k+1))
	temp_pop=temp_pop + (fy-fz)*p2dcssq
	! 0 +1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	!14  +1 -1
	uu=halfonecssq*(uh(i,j+1,k-1)*uh(i,j+1,k-1) + vh(i,j+1,k-1)*vh(i,j+1,k-1) + wh(i,j+1,k-1)*wh(i,j+1,k-1))
	udotc=(vh(i,j+1,k-1)-wh(i,j+1,k-1))*onecssq
	temp = -uu + half*udotc*udotc
#ifdef USESHARE 
    temp_pop=p2*(loc_rhoh(li,lj+1,lk-1)+(temp - udotc))
#else
	temp_pop=p2*(rhoh(i,j+1,k-1)+(temp - udotc))
#endif
	temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k-1)+qzz*pzzh(i,j+1,k-1)-cssq*pxxh(i,j+1,k-1)+two*qyz_13_14*pyzh(i,j+1,k-1))
	temp_pop=temp_pop + (fz-fy)*p2dcssq
	! 0 -1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
    	
	rho(i,j,k)=temp_rho
	
	u(i,j,k)=temp_u
	v(i,j,k)=temp_v
	w(i,j,k)=temp_w
	
	pxx(i,j,k)=temp_pxx
	pyy(i,j,k)=temp_pyy
	pzz(i,j,k)=temp_pzz
	pxy(i,j,k)=temp_pxy
	pxz(i,j,k)=temp_pxz
	pyz(i,j,k)=temp_pyz
    
    return
	
  end subroutine streamcoll_bulk_flop
  
  attributes(global) subroutine streamcoll_bulk_shared()
	
	implicit none  
	  
    integer :: i,j,k
	real(kind=db) :: udotc,uu
    
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
    
		  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
		  
	!if(isfluid(i,j,k).ne.1)return
    	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    uu=halfonecssq*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))
    
    !0
	f00(li,lj,lk)=p0*(rho(i,j,k)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k)))
	
    
	!1 -1  0  0
	udotc=u(i,j,k)*onecssq
	f01(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	udotc=u(i,j,k)*onecssq
	f02(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=v(i,j,k)*onecssq
	f03(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	udotc=v(i,j,k)*onecssq
	f04(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=w(i,j,k)*onecssq
	f05(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	udotc=w(i,j,k)*onecssq
	f06(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(u(i,j,k)+v(i,j,k))*onecssq
	f07(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+two*qxy_7_8*pxy(i,j,k)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	udotc=(u(i,j,k)+v(i,j,k))*onecssq
	f08(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+two*qxy_7_8*pxy(i,j,k)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-u(i,j,k)+v(i,j,k))*onecssq
	f10(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+two*qxy_9_10*pxy(i,j,k)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	udotc=(-u(i,j,k)+v(i,j,k))*onecssq
	f09(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+two*qxy_9_10*pxy(i,j,k)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(u(i,j,k)+w(i,j,k))*onecssq
	f15(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+two*qxz_15_16*pxz(i,j,k)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	udotc=(u(i,j,k)+w(i,j,k))*onecssq
	f16(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+two*qxz_15_16*pxz(i,j,k)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-u(i,j,k)+w(i,j,k))*onecssq
	f17(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+two*qxz_17_18*pxz(i,j,k)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	udotc=(-u(i,j,k)+w(i,j,k))*onecssq
	f18(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+two*qxz_17_18*pxz(i,j,k)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(v(i,j,k)+w(i,j,k))*onecssq
	f11(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+two*qyz_11_12*pyz(i,j,k)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	udotc=(v(i,j,k)+w(i,j,k))*onecssq
	f12(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+two*qyz_11_12*pyz(i,j,k)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(v(i,j,k)-w(i,j,k))*onecssq
	f13(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+two*qyz_13_14*pyz(i,j,k)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	udotc=(v(i,j,k)-w(i,j,k))*onecssq
	f14(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+two*qyz_13_14*pyz(i,j,k)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
        ! Halo Faces
    if(li==1)then
      
      uu=halfonecssq*(u(i-1,j,k)*u(i-1,j,k) + v(i-1,j,k)*v(i-1,j,k) + w(i-1,j,k)*w(i-1,j,k))
      
      !1 -1  0  0
	  udotc=u(i-1,j,k)*onecssq
	  f01(li-1,lj,lk)=p1*(rho(i-1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxx(i-1,j,k)-cssq*(pyy(i-1,j,k)+pzz(i-1,j,k))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(u(i-1,j,k)+v(i-1,j,k))*onecssq
	  f07(li-1,lj,lk)=p2*(rho(i-1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k)+qyy*pyy(i-1,j,k)-cssq*pzz(i-1,j,k)+two*qxy_7_8*pxy(i-1,j,k)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-u(i-1,j,k)+v(i-1,j,k))*onecssq
	  f09(li-1,lj,lk)=p2*(rho(i-1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k)+qyy*pyy(i-1,j,k)-cssq*pzz(i-1,j,k)+two*qxy_9_10*pxy(i-1,j,k)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(u(i-1,j,k)+w(i-1,j,k))*onecssq
	  f15(li-1,lj,lk)=p2*(rho(i-1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k)+qzz*pzz(i-1,j,k)-cssq*pyy(i-1,j,k)+two*qxz_15_16*pxz(i-1,j,k)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-u(i-1,j,k)+w(i-1,j,k))*onecssq
	  f18(li-1,lj,lk)=p2*(rho(i-1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k)+qzz*pzz(i-1,j,k)-cssq*pyy(i-1,j,k)+two*qxz_17_18*pxz(i-1,j,k)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      uu=halfonecssq*(u(i+1,j,k)*u(i+1,j,k) + v(i+1,j,k)*v(i+1,j,k) + w(i+1,j,k)*w(i+1,j,k))
      
      !2 +1  0  0
	  udotc=u(i+1,j,k)*onecssq
	  f02(li+1,lj,lk)=p1*(rho(i+1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxx(i+1,j,k)-cssq*(pyy(i+1,j,k)+pzz(i+1,j,k))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(u(i+1,j,k)+v(i+1,j,k))*onecssq
	  f08(li+1,lj,lk)=p2*(rho(i+1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k)+qyy*pyy(i+1,j,k)-cssq*pzz(i+1,j,k)+two*qxy_7_8*pxy(i+1,j,k)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-u(i+1,j,k)+v(i+1,j,k))*onecssq
	  f10(li+1,lj,lk)=p2*(rho(i+1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k)+qyy*pyy(i+1,j,k)-cssq*pzz(i+1,j,k)+two*qxy_9_10*pxy(i+1,j,k)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(u(i+1,j,k)+w(i+1,j,k))*onecssq
	  f16(li+1,lj,lk)=p2*(rho(i+1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k)+qzz*pzz(i+1,j,k)-cssq*pyy(i+1,j,k)+two*qxz_15_16*pxz(i+1,j,k)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-u(i+1,j,k)+w(i+1,j,k))*onecssq
	  f17(li+1,lj,lk)=p2*(rho(i+1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k)+qzz*pzz(i+1,j,k)-cssq*pyy(i+1,j,k)+two*qxz_17_18*pxz(i+1,j,k)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      uu=halfonecssq*(u(i,j-1,k)*u(i,j-1,k) + v(i,j-1,k)*v(i,j-1,k) + w(i,j-1,k)*w(i,j-1,k))
      
            !3 0 -1  0
	  udotc=v(i,j-1,k)*onecssq
	  f03(li,lj-1,lk)=p1*(rho(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyy(i,j-1,k)-cssq*(pxx(i,j-1,k)+pzz(i,j-1,k))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(u(i,j-1,k)+v(i,j-1,k))*onecssq
	  f07(li,lj-1,lk)=p2*(rho(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j-1,k)+qyy*pyy(i,j-1,k)-cssq*pzz(i,j-1,k)+two*qxy_7_8*pxy(i,j-1,k)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-u(i,j-1,k)+v(i,j-1,k))*onecssq
	  f10(li,lj-1,lk)=p2*(rho(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j-1,k)+qyy*pyy(i,j-1,k)-cssq*pzz(i,j-1,k)+two*qxy_9_10*pxy(i,j-1,k)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(v(i,j-1,k)+w(i,j-1,k))*onecssq
	  f11(li,lj-1,lk)=p2*(rho(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k)+qzz*pzz(i,j-1,k)-cssq*pxx(i,j-1,k)+two*qyz_11_12*pyz(i,j-1,k)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(v(i,j-1,k)-w(i,j-1,k))*onecssq
	  f13(li,lj-1,lk)=p2*(rho(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k)+qzz*pzz(i,j-1,k)-cssq*pxx(i,j-1,k)+two*qyz_13_14*pyz(i,j-1,k)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(u(i,j+1,k)*u(i,j+1,k) + v(i,j+1,k)*v(i,j+1,k) + w(i,j+1,k)*w(i,j+1,k))
      
      !4  0 +1  0
	  udotc=v(i,j+1,k)*onecssq
	  f04(li,lj+1,lk)=p1*(rho(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyy(i,j+1,k)-cssq*(pxx(i,j+1,k)+pzz(i,j+1,k))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(u(i,j+1,k)+v(i,j+1,k))*onecssq
	  f08(li,lj+1,lk)=p2*(rho(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j+1,k)+qyy*pyy(i,j+1,k)-cssq*pzz(i,j+1,k)+two*qxy_7_8*pxy(i,j+1,k)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-u(i,j+1,k)+v(i,j+1,k))*onecssq
	  f09(li,lj+1,lk)=p2*(rho(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j+1,k)+qyy*pyy(i,j+1,k)-cssq*pzz(i,j+1,k)+two*qxy_9_10*pxy(i,j+1,k)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(v(i,j+1,k)+w(i,j+1,k))*onecssq
	  f12(li,lj+1,lk)=p2*(rho(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k)+qzz*pzz(i,j+1,k)-cssq*pxx(i,j+1,k)+two*qyz_11_12*pyz(i,j+1,k)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(v(i,j+1,k)-w(i,j+1,k))*onecssq
	  f14(li,lj+1,lk)=p2*(rho(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k)+qzz*pzz(i,j+1,k)-cssq*pxx(i,j+1,k)+two*qyz_13_14*pyz(i,j+1,k)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      uu=halfonecssq*(u(i,j,k-1)*u(i,j,k-1) + v(i,j,k-1)*v(i,j,k-1) + w(i,j,k-1)*w(i,j,k-1))
      
      !5  0  0 -1
	  udotc=w(i,j,k-1)*onecssq
	  f05(li,lj,lk-1)=p1*(rho(i,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k-1)-cssq*(pxx(i,j,k-1)+pyy(i,j,k-1))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(u(i,j,k-1)+w(i,j,k-1))*onecssq
	  f15(li,lj,lk-1)=p2*(rho(i,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k-1)+qzz*pzz(i,j,k-1)-cssq*pyy(i,j,k-1)+two*qxz_15_16*pxz(i,j,k-1)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-u(i,j,k-1)+w(i,j,k-1))*onecssq
	  f17(li,lj,lk-1)=p2*(rho(i,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k-1)+qzz*pzz(i,j,k-1)-cssq*pyy(i,j,k-1)+two*qxz_17_18*pxz(i,j,k-1)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(v(i,j,k-1)+w(i,j,k-1))*onecssq
	  f11(li,lj,lk-1)=p2*(rho(i,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k-1)+qzz*pzz(i,j,k-1)-cssq*pxx(i,j,k-1)+two*qyz_11_12*pyz(i,j,k-1)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(v(i,j,k-1)-w(i,j,k-1))*onecssq
	  f14(li,lj,lk-1)=p2*(rho(i,j,k-1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k-1)+qzz*pzz(i,j,k-1)-cssq*pxx(i,j,k-1)+two*qyz_13_14*pyz(i,j,k-1)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i,j,k+1)*u(i,j,k+1) + v(i,j,k+1)*v(i,j,k+1) + w(i,j,k+1)*w(i,j,k+1))
      
      !6  0  0  +1
	  udotc=w(i,j,k+1)*onecssq
	  f06(li,lj,lk+1)=p1*(rho(i,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k+1)-cssq*(pxx(i,j,k+1)+pyy(i,j,k+1))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(u(i,j,k+1)+w(i,j,k+1))*onecssq
	  f16(li,lj,lk+1)=p2*(rho(i,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k+1)+qzz*pzz(i,j,k+1)-cssq*pyy(i,j,k+1)+two*qxz_15_16*pxz(i,j,k+1)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-u(i,j,k+1)+w(i,j,k+1))*onecssq
	  f18(li,lj,lk+1)=p2*(rho(i,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k+1)+qzz*pzz(i,j,k+1)-cssq*pyy(i,j,k+1)+two*qxz_17_18*pxz(i,j,k+1)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(v(i,j,k+1)+w(i,j,k+1))*onecssq
	  f12(li,lj,lk+1)=p2*(rho(i,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k+1)+qzz*pzz(i,j,k+1)-cssq*pxx(i,j,k+1)+two*qyz_11_12*pyz(i,j,k+1)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(v(i,j,k+1)-w(i,j,k+1))*onecssq
	  f13(li,lj,lk+1)=p2*(rho(i,j,k+1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k+1)+qzz*pzz(i,j,k+1)-cssq*pxx(i,j,k+1)+two*qyz_13_14*pyz(i,j,k+1)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      uu=halfonecssq*(u(i-1,j-1,k)*u(i-1,j-1,k) + v(i-1,j-1,k)*v(i-1,j-1,k) + w(i-1,j-1,k)*w(i-1,j-1,k))
      
      !7 -1 -1  0
	  udotc=(u(i-1,j-1,k)+v(i-1,j-1,k))*onecssq
	  f07(li-1,lj-1,lk)=p2*(rho(i-1,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j-1,k)+qyy*pyy(i-1,j-1,k)-cssq*pzz(i-1,j-1,k)+two*qxy_7_8*pxy(i-1,j-1,k)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(u(i-1,j+1,k)*u(i-1,j+1,k) + v(i-1,j+1,k)*v(i-1,j+1,k) + w(i-1,j+1,k)*w(i-1,j+1,k))
      
      !9  -1 +1 0
      udotc=(-u(i-1,j+1,k)+v(i-1,j+1,k))*onecssq
	  f09(li-1,lj+1,lk)=p2*(rho(i-1,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j+1,k)+qyy*pyy(i-1,j+1,k)-cssq*pzz(i-1,j+1,k)+two*qxy_9_10*pxy(i-1,j+1,k)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      uu=halfonecssq*(u(i-1,j,k-1)*u(i-1,j,k-1) + v(i-1,j,k-1)*v(i-1,j,k-1) + w(i-1,j,k-1)*w(i-1,j,k-1))
      
      !15  -1  0 -1
	  udotc=(u(i-1,j,k-1)+w(i-1,j,k-1))*onecssq
	  f15(li-1,lj,lk-1)=p2*(rho(i-1,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k-1)+qzz*pzz(i-1,j,k-1)-cssq*pyy(i-1,j,k-1)+two*qxz_15_16*pxz(i-1,j,k-1)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i-1,j,k+1)*u(i-1,j,k+1) + v(i-1,j,k+1)*v(i-1,j,k+1) + w(i-1,j,k+1)*w(i-1,j,k+1))
      
      !18   -1   0  +1
	  udotc=(-u(i-1,j,k+1)+w(i-1,j,k+1))*onecssq
	  f18(li-1,lj,lk+1)=p2*(rho(i-1,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k+1)+qzz*pzz(i-1,j,k+1)-cssq*pyy(i-1,j,k+1)+two*qxz_17_18*pxz(i-1,j,k+1)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      uu=halfonecssq*(u(i+1,j-1,k)*u(i+1,j-1,k) + v(i+1,j-1,k)*v(i+1,j-1,k) + w(i+1,j-1,k)*w(i+1,j-1,k))
      
      !10   +1 -1  0
	  udotc=(-u(i+1,j-1,k)+v(i+1,j-1,k))*onecssq
	  f10(li+1,lj-1,lk)=p2*(rho(i+1,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j-1,k)+qyy*pyy(i+1,j-1,k)-cssq*pzz(i+1,j-1,k)+two*qxy_9_10*pxy(i+1,j-1,k)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(u(i+1,j+1,k)*u(i+1,j+1,k) + v(i+1,j+1,k)*v(i+1,j+1,k) + w(i+1,j+1,k)*w(i+1,j+1,k))
      
      !8 +1 +1  0
	  udotc=(u(i+1,j+1,k)+v(i+1,j+1,k))*onecssq
	  f08(li+1,lj+1,lk)=p2*(rho(i+1,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j+1,k)+qyy*pyy(i+1,j+1,k)-cssq*pzz(i+1,j+1,k)+two*qxy_7_8*pxy(i+1,j+1,k)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      uu=halfonecssq*(u(i+1,j,k-1)*u(i+1,j,k-1) + v(i+1,j,k-1)*v(i+1,j,k-1) + w(i+1,j,k-1)*w(i+1,j,k-1))
      
      !17  +1  0 -1
	  udotc=(-u(i+1,j,k-1)+w(i+1,j,k-1))*onecssq
	  f17(li+1,lj,lk-1)=p2*(rho(i+1,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k-1)+qzz*pzz(i+1,j,k-1)-cssq*pyy(i+1,j,k-1)+two*qxz_17_18*pxz(i+1,j,k-1)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i+1,j,k+1)*u(i+1,j,k+1) + v(i+1,j,k+1)*v(i+1,j,k+1) + w(i+1,j,k+1)*w(i+1,j,k+1))
      
      !16  +1  0 +1
	  udotc=(u(i+1,j,k+1)+w(i+1,j,k+1))*onecssq
	  f16(li+1,lj,lk+1)=p2*(rho(i+1,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k+1)+qzz*pzz(i+1,j,k+1)-cssq*pyy(i+1,j,k+1)+two*qxz_15_16*pxz(i+1,j,k+1)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      uu=halfonecssq*(u(i,j-1,k-1)*u(i,j-1,k-1) + v(i,j-1,k-1)*v(i,j-1,k-1) + w(i,j-1,k-1)*w(i,j-1,k-1))
      
      !11  0  -1  -1
	  udotc=(v(i,j-1,k-1)+w(i,j-1,k-1))*onecssq
	  f11(li,lj-1,lk-1)=p2*(rho(i,j-1,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k-1)+qzz*pzz(i,j-1,k-1)-cssq*pxx(i,j-1,k-1)+two*qyz_11_12*pyz(i,j-1,k-1)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i,j-1,k+1)*u(i,j-1,k+1) + v(i,j-1,k+1)*v(i,j-1,k+1) + w(i,j-1,k+1)*w(i,j-1,k+1))
      
      !13  0  -1   +1
	  udotc=(v(i,j-1,k+1)-w(i,j-1,k+1))*onecssq
	  f13(li,lj-1,lk+1)=p2*(rho(i,j-1,k+1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k+1)+qzz*pzz(i,j-1,k+1)-cssq*pxx(i,j-1,k+1)+two*qyz_13_14*pyz(i,j-1,k+1)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      uu=halfonecssq*(u(i,j+1,k-1)*u(i,j+1,k-1) + v(i,j+1,k-1)*v(i,j+1,k-1) + w(i,j+1,k-1)*w(i,j+1,k-1))
      
      !14  0  +1  -1
	  udotc=(v(i,j+1,k-1)-w(i,j+1,k-1))*onecssq
	  f14(li,lj+1,lk-1)=p2*(rho(i,j+1,k-1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k-1)+qzz*pzz(i,j+1,k-1)-cssq*pxx(i,j+1,k-1)+two*qyz_13_14*pyz(i,j+1,k-1)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i,j+1,k+1)*u(i,j+1,k+1) + v(i,j+1,k+1)*v(i,j+1,k+1) + w(i,j+1,k+1)*w(i,j+1,k+1))
      
      !12  0  +1  +1
	  udotc=(v(i,j+1,k+1)+w(i,j+1,k+1))*onecssq
	  f12(li,lj+1,lk+1)=p2*(rho(i,j+1,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k+1)+qzz*pzz(i,j+1,k+1)-cssq*pxx(i,j+1,k+1)+two*qyz_11_12*pyz(i,j+1,k+1)) &
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
	rhoh(i,j,k)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	uh(i,j,k)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	vh(i,j,k)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	wh(i,j,k)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pxxh(i,j,k)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	pyyh(i,j,k)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pzzh(i,j,k)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	pxyh(i,j,k)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	pxzh(i,j,k)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	pyzh(i,j,k)=udotc
     return
#if 0
	call syncthreads
	
	uu=halfonecssq*(uh(i,j,k)*uh(i,j,k) + vh(i,j,k)*vh(i,j,k) + wh(i,j,k)*wh(i,j,k))
    
	!1 -1  0  0
	udotc=uh(i,j,k)*onecssq
	f01(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!+1  0  0


	!2 +1  0  0
	udotc=uh(i,j,k)*onecssq
	f02(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=vh(i,j,k)*onecssq
	f03(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	udotc=vh(i,j,k)*onecssq
	f04(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=wh(i,j,k)*onecssq
	f05(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	! 0  0 +1


	!6  0  0  +1
	udotc=wh(i,j,k)*onecssq
	f06(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(uh(i,j,k)+vh(i,j,k))*onecssq
	f07(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	udotc=(uh(i,j,k)+vh(i,j,k))*onecssq
	f08(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-uh(i,j,k)+vh(i,j,k))*onecssq
	f10(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	udotc=(-uh(i,j,k)+vh(i,j,k))*onecssq
	f09(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(uh(i,j,k)+wh(i,j,k))*onecssq
	f15(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!+1  0  +1


	!16  +1  0 +1
	udotc=(uh(i,j,k)+wh(i,j,k))*onecssq
	f16(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-uh(i,j,k)+wh(i,j,k))*onecssq
	f17(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!-1  0  +1


	!18   -1   0  +1
	udotc=(-uh(i,j,k)+wh(i,j,k))*onecssq
	f18(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(vh(i,j,k)+wh(i,j,k))*onecssq
	f11(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	udotc=(vh(i,j,k)+wh(i,j,k))*onecssq
	f12(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(vh(i,j,k)-wh(i,j,k))*onecssq
	f13(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	udotc=(vh(i,j,k)-wh(i,j,k))*onecssq
	f14(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pxxh(i,j,k)=pxxh(i,j,k)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	pyyh(i,j,k)=pyyh(i,j,k)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pzzh(i,j,k)=pzzh(i,j,k)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	pxyh(i,j,k)=pxyh(i,j,k)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	pxzh(i,j,k)=pxzh(i,j,k)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	pyzh(i,j,k)=pyzh(i,j,k)-udotc
	
	    
    return
#endif	
  end subroutine streamcoll_bulk_shared
  
  attributes(global) subroutine streamcoll_bulk_shared_flop()
	
	implicit none  
	  
    integer :: i,j,k
	real(kind=db) :: udotc,uu
    
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
    
		  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
		  
	!if(isfluid(i,j,k).ne.1)return
    	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    uu=halfonecssq*(uh(i,j,k)*uh(i,j,k) + vh(i,j,k)*vh(i,j,k) + wh(i,j,k)*wh(i,j,k))
    
    !0
	f00(li,lj,lk)=p0*(rhoh(i,j,k)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(pyyh(i,j,k)+pxxh(i,j,k)+pzzh(i,j,k)))
	
    
	!1 -1  0  0
	udotc=uh(i,j,k)*onecssq
	f01(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k)-cssq*(pyyh(i,j,k)+pzzh(i,j,k))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	udotc=uh(i,j,k)*onecssq
	f02(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k)-cssq*(pyyh(i,j,k)+pzzh(i,j,k))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=vh(i,j,k)*onecssq
	f03(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k)-cssq*(pxxh(i,j,k)+pzzh(i,j,k))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	udotc=vh(i,j,k)*onecssq
	f04(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k)-cssq*(pxxh(i,j,k)+pzzh(i,j,k))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=wh(i,j,k)*onecssq
	f05(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k)-cssq*(pxxh(i,j,k)+pyyh(i,j,k))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	udotc=wh(i,j,k)*onecssq
	f06(li,lj,lk)=p1*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k)-cssq*(pxxh(i,j,k)+pyyh(i,j,k))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(uh(i,j,k)+vh(i,j,k))*onecssq
	f07(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qyy*pyyh(i,j,k)-cssq*pzzh(i,j,k)+two*qxy_7_8*pxyh(i,j,k)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	udotc=(uh(i,j,k)+vh(i,j,k))*onecssq
	f08(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qyy*pyyh(i,j,k)-cssq*pzzh(i,j,k)+two*qxy_7_8*pxyh(i,j,k)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-uh(i,j,k)+vh(i,j,k))*onecssq
	f10(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qyy*pyyh(i,j,k)-cssq*pzzh(i,j,k)+two*qxy_9_10*pxyh(i,j,k)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	udotc=(-uh(i,j,k)+vh(i,j,k))*onecssq
	f09(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qyy*pyyh(i,j,k)-cssq*pzzh(i,j,k)+two*qxy_9_10*pxyh(i,j,k)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(uh(i,j,k)+wh(i,j,k))*onecssq
	f15(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pyyh(i,j,k)+two*qxz_15_16*pxzh(i,j,k)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	udotc=(uh(i,j,k)+wh(i,j,k))*onecssq
	f16(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pyyh(i,j,k)+two*qxz_15_16*pxzh(i,j,k)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-uh(i,j,k)+wh(i,j,k))*onecssq
	f17(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pyyh(i,j,k)+two*qxz_17_18*pxzh(i,j,k)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	udotc=(-uh(i,j,k)+wh(i,j,k))*onecssq
	f18(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pyyh(i,j,k)+two*qxz_17_18*pxzh(i,j,k)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(vh(i,j,k)+wh(i,j,k))*onecssq
	f11(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pxxh(i,j,k)+two*qyz_11_12*pyzh(i,j,k)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	udotc=(vh(i,j,k)+wh(i,j,k))*onecssq
	f12(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pxxh(i,j,k)+two*qyz_11_12*pyzh(i,j,k)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(vh(i,j,k)-wh(i,j,k))*onecssq
	f13(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pxxh(i,j,k)+two*qyz_13_14*pyzh(i,j,k)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	udotc=(vh(i,j,k)-wh(i,j,k))*onecssq
	f14(li,lj,lk)=p2*(rhoh(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pxxh(i,j,k)+two*qyz_13_14*pyzh(i,j,k)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    ! Halo Faces
    if(li==1)then
      
      uu=halfonecssq*(uh(i-1,j,k)*uh(i-1,j,k) + vh(i-1,j,k)*vh(i-1,j,k) + wh(i-1,j,k)*wh(i-1,j,k))
      
      !1 -1  0  0
	  udotc=uh(i-1,j,k)*onecssq
	  f01(li-1,lj,lk)=p1*(rhoh(i-1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxxh(i-1,j,k)-cssq*(pyyh(i-1,j,k)+pzzh(i-1,j,k))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(uh(i-1,j,k)+vh(i-1,j,k))*onecssq
	  f07(li-1,lj,lk)=p2*(rhoh(i-1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k)+qyy*pyyh(i-1,j,k)-cssq*pzzh(i-1,j,k)+two*qxy_7_8*pxyh(i-1,j,k)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-uh(i-1,j,k)+vh(i-1,j,k))*onecssq
	  f09(li-1,lj,lk)=p2*(rhoh(i-1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k)+qyy*pyyh(i-1,j,k)-cssq*pzzh(i-1,j,k)+two*qxy_9_10*pxyh(i-1,j,k)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(uh(i-1,j,k)+wh(i-1,j,k))*onecssq
	  f15(li-1,lj,lk)=p2*(rhoh(i-1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k)+qzz*pzzh(i-1,j,k)-cssq*pyyh(i-1,j,k)+two*qxz_15_16*pxzh(i-1,j,k)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-uh(i-1,j,k)+wh(i-1,j,k))*onecssq
	  f18(li-1,lj,lk)=p2*(rhoh(i-1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k)+qzz*pzzh(i-1,j,k)-cssq*pyyh(i-1,j,k)+two*qxz_17_18*pxzh(i-1,j,k)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      uu=halfonecssq*(uh(i+1,j,k)*uh(i+1,j,k) + vh(i+1,j,k)*vh(i+1,j,k) + wh(i+1,j,k)*wh(i+1,j,k))
      
      !2 +1  0  0
	  udotc=uh(i+1,j,k)*onecssq
	  f02(li+1,lj,lk)=p1*(rhoh(i+1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxxh(i+1,j,k)-cssq*(pyyh(i+1,j,k)+pzzh(i+1,j,k))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(uh(i+1,j,k)+vh(i+1,j,k))*onecssq
	  f08(li+1,lj,lk)=p2*(rhoh(i+1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k)+qyy*pyyh(i+1,j,k)-cssq*pzzh(i+1,j,k)+two*qxy_7_8*pxyh(i+1,j,k)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-uh(i+1,j,k)+vh(i+1,j,k))*onecssq
	  f10(li+1,lj,lk)=p2*(rhoh(i+1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k)+qyy*pyyh(i+1,j,k)-cssq*pzzh(i+1,j,k)+two*qxy_9_10*pxyh(i+1,j,k)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(uh(i+1,j,k)+wh(i+1,j,k))*onecssq
	  f16(li+1,lj,lk)=p2*(rhoh(i+1,j,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k)+qzz*pzzh(i+1,j,k)-cssq*pyyh(i+1,j,k)+two*qxz_15_16*pxzh(i+1,j,k)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-uh(i+1,j,k)+wh(i+1,j,k))*onecssq
	  f17(li+1,lj,lk)=p2*(rhoh(i+1,j,k)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k)+qzz*pzzh(i+1,j,k)-cssq*pyyh(i+1,j,k)+two*qxz_17_18*pxzh(i+1,j,k)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      uu=halfonecssq*(uh(i,j-1,k)*uh(i,j-1,k) + vh(i,j-1,k)*vh(i,j-1,k) + wh(i,j-1,k)*wh(i,j-1,k))
      
            !3 0 -1  0
	  udotc=vh(i,j-1,k)*onecssq
	  f03(li,lj-1,lk)=p1*(rhoh(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyyh(i,j-1,k)-cssq*(pxxh(i,j-1,k)+pzzh(i,j-1,k))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(uh(i,j-1,k)+vh(i,j-1,k))*onecssq
	  f07(li,lj-1,lk)=p2*(rhoh(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j-1,k)+qyy*pyyh(i,j-1,k)-cssq*pzzh(i,j-1,k)+two*qxy_7_8*pxyh(i,j-1,k)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-uh(i,j-1,k)+vh(i,j-1,k))*onecssq
	  f10(li,lj-1,lk)=p2*(rhoh(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j-1,k)+qyy*pyyh(i,j-1,k)-cssq*pzzh(i,j-1,k)+two*qxy_9_10*pxyh(i,j-1,k)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(vh(i,j-1,k)+wh(i,j-1,k))*onecssq
	  f11(li,lj-1,lk)=p2*(rhoh(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k)+qzz*pzzh(i,j-1,k)-cssq*pxxh(i,j-1,k)+two*qyz_11_12*pyzh(i,j-1,k)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(vh(i,j-1,k)-wh(i,j-1,k))*onecssq
	  f13(li,lj-1,lk)=p2*(rhoh(i,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k)+qzz*pzzh(i,j-1,k)-cssq*pxxh(i,j-1,k)+two*qyz_13_14*pyzh(i,j-1,k)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(uh(i,j+1,k)*uh(i,j+1,k) + vh(i,j+1,k)*vh(i,j+1,k) + wh(i,j+1,k)*wh(i,j+1,k))
      
      !4  0 +1  0
	  udotc=vh(i,j+1,k)*onecssq
	  f04(li,lj+1,lk)=p1*(rhoh(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyyh(i,j+1,k)-cssq*(pxxh(i,j+1,k)+pzzh(i,j+1,k))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(uh(i,j+1,k)+vh(i,j+1,k))*onecssq
	  f08(li,lj+1,lk)=p2*(rhoh(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j+1,k)+qyy*pyyh(i,j+1,k)-cssq*pzzh(i,j+1,k)+two*qxy_7_8*pxyh(i,j+1,k)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-uh(i,j+1,k)+vh(i,j+1,k))*onecssq
	  f09(li,lj+1,lk)=p2*(rhoh(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j+1,k)+qyy*pyyh(i,j+1,k)-cssq*pzzh(i,j+1,k)+two*qxy_9_10*pxyh(i,j+1,k)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(vh(i,j+1,k)+wh(i,j+1,k))*onecssq
	  f12(li,lj+1,lk)=p2*(rhoh(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k)+qzz*pzzh(i,j+1,k)-cssq*pxxh(i,j+1,k)+two*qyz_11_12*pyzh(i,j+1,k)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(vh(i,j+1,k)-wh(i,j+1,k))*onecssq
	  f14(li,lj+1,lk)=p2*(rhoh(i,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k)+qzz*pzzh(i,j+1,k)-cssq*pxxh(i,j+1,k)+two*qyz_13_14*pyzh(i,j+1,k)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      uu=halfonecssq*(uh(i,j,k-1)*uh(i,j,k-1) + vh(i,j,k-1)*vh(i,j,k-1) + wh(i,j,k-1)*wh(i,j,k-1))
      
      !5  0  0 -1
	  udotc=wh(i,j,k-1)*onecssq
	  f05(li,lj,lk-1)=p1*(rhoh(i,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k-1)-cssq*(pxxh(i,j,k-1)+pyyh(i,j,k-1))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(uh(i,j,k-1)+wh(i,j,k-1))*onecssq
	  f15(li,lj,lk-1)=p2*(rhoh(i,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k-1)+qzz*pzzh(i,j,k-1)-cssq*pyyh(i,j,k-1)+two*qxz_15_16*pxzh(i,j,k-1)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-uh(i,j,k-1)+wh(i,j,k-1))*onecssq
	  f17(li,lj,lk-1)=p2*(rhoh(i,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k-1)+qzz*pzzh(i,j,k-1)-cssq*pyyh(i,j,k-1)+two*qxz_17_18*pxzh(i,j,k-1)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(vh(i,j,k-1)+wh(i,j,k-1))*onecssq
	  f11(li,lj,lk-1)=p2*(rhoh(i,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k-1)+qzz*pzzh(i,j,k-1)-cssq*pxxh(i,j,k-1)+two*qyz_11_12*pyzh(i,j,k-1)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(vh(i,j,k-1)-wh(i,j,k-1))*onecssq
	  f14(li,lj,lk-1)=p2*(rhoh(i,j,k-1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k-1)+qzz*pzzh(i,j,k-1)-cssq*pxxh(i,j,k-1)+two*qyz_13_14*pyzh(i,j,k-1)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i,j,k+1)*uh(i,j,k+1) + vh(i,j,k+1)*vh(i,j,k+1) + wh(i,j,k+1)*wh(i,j,k+1))
      
      !6  0  0  +1
	  udotc=wh(i,j,k+1)*onecssq
	  f06(li,lj,lk+1)=p1*(rhoh(i,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k+1)-cssq*(pxxh(i,j,k+1)+pyyh(i,j,k+1))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(uh(i,j,k+1)+wh(i,j,k+1))*onecssq
	  f16(li,lj,lk+1)=p2*(rhoh(i,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k+1)+qzz*pzzh(i,j,k+1)-cssq*pyyh(i,j,k+1)+two*qxz_15_16*pxzh(i,j,k+1)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-uh(i,j,k+1)+wh(i,j,k+1))*onecssq
	  f18(li,lj,lk+1)=p2*(rhoh(i,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k+1)+qzz*pzzh(i,j,k+1)-cssq*pyyh(i,j,k+1)+two*qxz_17_18*pxzh(i,j,k+1)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(vh(i,j,k+1)+wh(i,j,k+1))*onecssq
	  f12(li,lj,lk+1)=p2*(rhoh(i,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k+1)+qzz*pzzh(i,j,k+1)-cssq*pxxh(i,j,k+1)+two*qyz_11_12*pyzh(i,j,k+1)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(vh(i,j,k+1)-wh(i,j,k+1))*onecssq
	  f13(li,lj,lk+1)=p2*(rhoh(i,j,k+1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k+1)+qzz*pzzh(i,j,k+1)-cssq*pxxh(i,j,k+1)+two*qyz_13_14*pyzh(i,j,k+1)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      uu=halfonecssq*(uh(i-1,j-1,k)*uh(i-1,j-1,k) + vh(i-1,j-1,k)*vh(i-1,j-1,k) + wh(i-1,j-1,k)*wh(i-1,j-1,k))
      
      !7 -1 -1  0
	  udotc=(uh(i-1,j-1,k)+vh(i-1,j-1,k))*onecssq
	  f07(li-1,lj-1,lk)=p2*(rhoh(i-1,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j-1,k)+qyy*pyyh(i-1,j-1,k)-cssq*pzzh(i-1,j-1,k)+two*qxy_7_8*pxyh(i-1,j-1,k)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(uh(i-1,j+1,k)*uh(i-1,j+1,k) + vh(i-1,j+1,k)*vh(i-1,j+1,k) + wh(i-1,j+1,k)*wh(i-1,j+1,k))
      
      !9  -1 +1 0
      udotc=(-uh(i-1,j+1,k)+vh(i-1,j+1,k))*onecssq
	  f09(li-1,lj+1,lk)=p2*(rhoh(i-1,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j+1,k)+qyy*pyyh(i-1,j+1,k)-cssq*pzzh(i-1,j+1,k)+two*qxy_9_10*pxyh(i-1,j+1,k)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      uu=halfonecssq*(uh(i-1,j,k-1)*uh(i-1,j,k-1) + vh(i-1,j,k-1)*vh(i-1,j,k-1) + wh(i-1,j,k-1)*wh(i-1,j,k-1))
      
      !15  -1  0 -1
	  udotc=(uh(i-1,j,k-1)+wh(i-1,j,k-1))*onecssq
	  f15(li-1,lj,lk-1)=p2*(rhoh(i-1,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k-1)+qzz*pzzh(i-1,j,k-1)-cssq*pyyh(i-1,j,k-1)+two*qxz_15_16*pxzh(i-1,j,k-1)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i-1,j,k+1)*uh(i-1,j,k+1) + vh(i-1,j,k+1)*vh(i-1,j,k+1) + wh(i-1,j,k+1)*wh(i-1,j,k+1))
      
      !18   -1   0  +1
	  udotc=(-uh(i-1,j,k+1)+wh(i-1,j,k+1))*onecssq
	  f18(li-1,lj,lk+1)=p2*(rhoh(i-1,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k+1)+qzz*pzzh(i-1,j,k+1)-cssq*pyyh(i-1,j,k+1)+two*qxz_17_18*pxzh(i-1,j,k+1)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      uu=halfonecssq*(uh(i+1,j-1,k)*uh(i+1,j-1,k) + vh(i+1,j-1,k)*vh(i+1,j-1,k) + wh(i+1,j-1,k)*wh(i+1,j-1,k))
      
      !10   +1 -1  0
	  udotc=(-uh(i+1,j-1,k)+vh(i+1,j-1,k))*onecssq
	  f10(li+1,lj-1,lk)=p2*(rhoh(i+1,j-1,k)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j-1,k)+qyy*pyyh(i+1,j-1,k)-cssq*pzzh(i+1,j-1,k)+two*qxy_9_10*pxyh(i+1,j-1,k)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(uh(i+1,j+1,k)*uh(i+1,j+1,k) + vh(i+1,j+1,k)*vh(i+1,j+1,k) + wh(i+1,j+1,k)*wh(i+1,j+1,k))
      
      !8 +1 +1  0
	  udotc=(uh(i+1,j+1,k)+vh(i+1,j+1,k))*onecssq
	  f08(li+1,lj+1,lk)=p2*(rhoh(i+1,j+1,k)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j+1,k)+qyy*pyyh(i+1,j+1,k)-cssq*pzzh(i+1,j+1,k)+two*qxy_7_8*pxyh(i+1,j+1,k)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      uu=halfonecssq*(uh(i+1,j,k-1)*uh(i+1,j,k-1) + vh(i+1,j,k-1)*vh(i+1,j,k-1) + wh(i+1,j,k-1)*wh(i+1,j,k-1))
      
      !17  +1  0 -1
	  udotc=(-uh(i+1,j,k-1)+wh(i+1,j,k-1))*onecssq
	  f17(li+1,lj,lk-1)=p2*(rhoh(i+1,j,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k-1)+qzz*pzzh(i+1,j,k-1)-cssq*pyyh(i+1,j,k-1)+two*qxz_17_18*pxzh(i+1,j,k-1)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i+1,j,k+1)*uh(i+1,j,k+1) + vh(i+1,j,k+1)*vh(i+1,j,k+1) + wh(i+1,j,k+1)*wh(i+1,j,k+1))
      
      !16  +1  0 +1
	  udotc=(uh(i+1,j,k+1)+wh(i+1,j,k+1))*onecssq
	  f16(li+1,lj,lk+1)=p2*(rhoh(i+1,j,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k+1)+qzz*pzzh(i+1,j,k+1)-cssq*pyyh(i+1,j,k+1)+two*qxz_15_16*pxzh(i+1,j,k+1)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      uu=halfonecssq*(uh(i,j-1,k-1)*uh(i,j-1,k-1) + vh(i,j-1,k-1)*vh(i,j-1,k-1) + wh(i,j-1,k-1)*wh(i,j-1,k-1))
      
      !11  0  -1  -1
	  udotc=(vh(i,j-1,k-1)+wh(i,j-1,k-1))*onecssq
	  f11(li,lj-1,lk-1)=p2*(rhoh(i,j-1,k-1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k-1)+qzz*pzzh(i,j-1,k-1)-cssq*pxxh(i,j-1,k-1)+two*qyz_11_12*pyzh(i,j-1,k-1)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i,j-1,k+1)*uh(i,j-1,k+1) + vh(i,j-1,k+1)*vh(i,j-1,k+1) + wh(i,j-1,k+1)*wh(i,j-1,k+1))
      
      !13  0  -1   +1
	  udotc=(vh(i,j-1,k+1)-wh(i,j-1,k+1))*onecssq
	  f13(li,lj-1,lk+1)=p2*(rhoh(i,j-1,k+1)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k+1)+qzz*pzzh(i,j-1,k+1)-cssq*pxxh(i,j-1,k+1)+two*qyz_13_14*pyzh(i,j-1,k+1)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      uu=halfonecssq*(uh(i,j+1,k-1)*uh(i,j+1,k-1) + vh(i,j+1,k-1)*vh(i,j+1,k-1) + wh(i,j+1,k-1)*wh(i,j+1,k-1))
      
      !14  0  +1  -1
	  udotc=(vh(i,j+1,k-1)-wh(i,j+1,k-1))*onecssq
	  f14(li,lj+1,lk-1)=p2*(rhoh(i,j+1,k-1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k-1)+qzz*pzzh(i,j+1,k-1)-cssq*pxxh(i,j+1,k-1)+two*qyz_13_14*pyzh(i,j+1,k-1)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i,j+1,k+1)*uh(i,j+1,k+1) + vh(i,j+1,k+1)*vh(i,j+1,k+1) + wh(i,j+1,k+1)*wh(i,j+1,k+1))
      
      !12  0  +1  +1
	  udotc=(vh(i,j+1,k+1)+wh(i,j+1,k+1))*onecssq
	  f12(li,lj+1,lk+1)=p2*(rhoh(i,j+1,k+1)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k+1)+qzz*pzzh(i,j+1,k+1)-cssq*pxxh(i,j+1,k+1)+two*qyz_11_12*pyzh(i,j+1,k+1)) &
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
	rho(i,j,k)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	u(i,j,k)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	v(i,j,k)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	w(i,j,k)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pxx(i,j,k)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	pyy(i,j,k)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pzz(i,j,k)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	pxy(i,j,k)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	pxz(i,j,k)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	pyz(i,j,k)=udotc
     return
#if 0     
	call syncthreads
	
	uu=halfonecssq*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))
    
	!1 -1  0  0
	udotc=u(i,j,k)*onecssq
	f01(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!+1  0  0


	!2 +1  0  0
	udotc=u(i,j,k)*onecssq
	f02(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=v(i,j,k)*onecssq
	f03(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	udotc=v(i,j,k)*onecssq
	f04(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=w(i,j,k)*onecssq
	f05(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	! 0  0 +1


	!6  0  0  +1
	udotc=w(i,j,k)*onecssq
	f06(li,lj,lk)=p1*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(u(i,j,k)+v(i,j,k))*onecssq
	f07(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	udotc=(u(i,j,k)+v(i,j,k))*onecssq
	f08(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-u(i,j,k)+v(i,j,k))*onecssq
	f10(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	udotc=(-u(i,j,k)+v(i,j,k))*onecssq
	f09(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(u(i,j,k)+w(i,j,k))*onecssq
	f15(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!+1  0  +1


	!16  +1  0 +1
	udotc=(u(i,j,k)+w(i,j,k))*onecssq
	f16(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-u(i,j,k)+w(i,j,k))*onecssq
	f17(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	!-1  0  +1


	!18   -1   0  +1
	udotc=(-u(i,j,k)+w(i,j,k))*onecssq
	f18(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(v(i,j,k)+w(i,j,k))*onecssq
	f11(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	udotc=(v(i,j,k)+w(i,j,k))*onecssq
	f12(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(v(i,j,k)-w(i,j,k))*onecssq
	f13(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	udotc=(v(i,j,k)-w(i,j,k))*onecssq
	f14(li,lj,lk)=p2*(rho(i,j,k)+(-uu &
	 + half*udotc*udotc - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pxx(i,j,k)=pxx(i,j,k)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	pyy(i,j,k)=pyy(i,j,k)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pzz(i,j,k)=pzz(i,j,k)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	pxy(i,j,k)=pxy(i,j,k)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	pxz(i,j,k)=pxz(i,j,k)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	pyz(i,j,k)=pyz(i,j,k)-udotc
	
	    
    return
#endif	
  end subroutine streamcoll_bulk_shared_flop
  
 end module streamcoll_bulk_kernels
