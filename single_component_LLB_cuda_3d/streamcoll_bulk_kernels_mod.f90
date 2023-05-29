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
	temp_u=temp_pop
	temp_v=temp_pop
	temp_w=temp_pop
	temp_pxx=temp_pop
	temp_pyy=temp_pop
	temp_pzz=temp_pop
	temp_pxy=temp_pop
	temp_pxz=temp_pop
	temp_pyz=temp_pop

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
	temp_u=temp_pop
	temp_v=temp_pop
	temp_w=temp_pop
	temp_pxx=temp_pop
	temp_pyy=temp_pop
	temp_pzz=temp_pop
	temp_pxy=temp_pop
	temp_pxz=temp_pop
	temp_pyz=temp_pop

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
  
 end module streamcoll_bulk_kernels
