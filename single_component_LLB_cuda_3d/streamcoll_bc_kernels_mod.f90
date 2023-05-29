#include "defines.h"
 module streamcoll_bc_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine streamcoll_bc()
	
	implicit none  
	  
    integer :: i,j,k
    !logical :: alltrue
	real(kind=db) :: uu,udotc,temp,uu0
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
		  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
    
    !alltrue = .false.
    !if( allthreads(isfluid(i,j,k).ne.-1) )alltrue = .true.
    !call syncthreads
    !if(alltrue)return
    
    if(isfluid(i,j,k).ne.-1)return
    
	uu0=halfonecssq*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))
	!0
	temp_rho=p0*(rho(i,j,k)-uu0)
	temp_rho=temp_rho + oneminusomega*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k)))
	
	!1
	if(isfluid(i-1,j,k).eq.0)then
	  udotc=u(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
	  temp_pop=temp_pop - fx*p1dcssq
	else
	  uu=halfonecssq*(u(i-1,j,k)*u(i-1,j,k) + v(i-1,j,k)*v(i-1,j,k) + w(i-1,j,k)*w(i-1,j,k))
	  udotc=u(i-1,j,k)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i-1,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i-1,j,k)-cssq*(pyy(i-1,j,k)+pzz(i-1,j,k)))
	  temp_pop=temp_pop + fx*p1dcssq
	endif
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop
	
	!2
	if(isfluid(i+1,j,k).eq.0)then
	  udotc=u(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
	  temp_pop=temp_pop + fx*p1dcssq
	else
	  uu=halfonecssq*(u(i+1,j,k)*u(i+1,j,k) + v(i+1,j,k)*v(i+1,j,k) + w(i+1,j,k)*w(i+1,j,k))
	  udotc=u(i+1,j,k)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i+1,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i+1,j,k)-cssq*(pyy(i+1,j,k)+pzz(i+1,j,k)))
	  temp_pop=temp_pop - fx*p1dcssq
	endif
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
	
	!3
	if(isfluid(i,j-1,k).eq.0)then
	  udotc=v(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
	  temp_pop=temp_pop - fy*p1dcssq
	else
	  uu=halfonecssq*(u(i,j-1,k)*u(i,j-1,k) + v(i,j-1,k)*v(i,j-1,k) + w(i,j-1,k)*w(i,j-1,k))
	  udotc=v(i,j-1,k)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i,j-1,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j-1,k)-cssq*(pxx(i,j-1,k)+pzz(i,j-1,k)))
	  temp_pop=temp_pop + fy*p1dcssq
	endif
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	if(isfluid(i,j+1,k).eq.0)then
	  udotc=v(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
	  temp_pop=temp_pop + fy*p1dcssq
	else
	  uu=halfonecssq*(u(i,j+1,k)*u(i,j+1,k) + v(i,j+1,k)*v(i,j+1,k) + w(i,j+1,k)*w(i,j+1,k))
	  udotc=v(i,j+1,k)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i,j+1,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j+1,k)-cssq*(pxx(i,j+1,k)+pzz(i,j+1,k)))
	  temp_pop=temp_pop - fy*p1dcssq
	endif
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	if(isfluid(i,j,k-1).eq.0)then
	  udotc=w(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(u(i,j,k-1)*u(i,j,k-1) + v(i,j,k-1)*v(i,j,k-1) + w(i,j,k-1)*w(i,j,k-1))
	  udotc=w(i,j,k-1)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k-1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k-1)-cssq*(pxx(i,j,k-1)+pyy(i,j,k-1)))
	  temp_pop=temp_pop + fz*p1dcssq
	endif
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop
	
	!6 +1
	if(isfluid(i,j,k+1).eq.0)then
	  udotc=w(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(u(i,j,k+1)*u(i,j,k+1) + v(i,j,k+1)*v(i,j,k+1) + w(i,j,k+1)*w(i,j,k+1))
	  udotc=w(i,j,k+1)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k+1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k+1)-cssq*(pxx(i,j,k+1)+pyy(i,j,k+1)))
	  temp_pop=temp_pop - fz*p1dcssq
	endif
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
	
	!7 -1 -1
	if(isfluid(i-1,j-1,k).eq.0)then
	  udotc=(u(i,j,k)+v(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+two*qxy_7_8*pxy(i,j,k))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	else
	  uu=halfonecssq*(u(i-1,j-1,k)*u(i-1,j-1,k) + v(i-1,j-1,k)*v(i-1,j-1,k) + w(i-1,j-1,k)*w(i-1,j-1,k))
	  udotc=(u(i-1,j-1,k)+v(i-1,j-1,k))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i-1,j-1,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j-1,k)+qyy*pyy(i-1,j-1,k)-cssq*pzz(i-1,j-1,k)+two*qxy_7_8*pxy(i-1,j-1,k))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	endif
	!+1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pop
	
	!8 +1 +1
	if(isfluid(i+1,j+1,k).eq.0)then
	  udotc=(u(i,j,k)+v(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+two*qxy_7_8*pxy(i,j,k))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	else
	  uu=halfonecssq*(u(i+1,j+1,k)*u(i+1,j+1,k) + v(i+1,j+1,k)*v(i+1,j+1,k) + w(i+1,j+1,k)*w(i+1,j+1,k))
	  udotc=(u(i+1,j+1,k)+v(i+1,j+1,k))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i+1,j+1,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j+1,k)+qyy*pyy(i+1,j+1,k)-cssq*pzz(i+1,j+1,k)+two*qxy_7_8*pxy(i+1,j+1,k))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	endif
	!-1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy+temp_pop
	
	!10   +1 -1
	if(isfluid(i+1,j-1,k).eq.0)then
	  udotc=(-u(i,j,k)+v(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+two*qxy_9_10*pxy(i,j,k))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	else
	  uu=halfonecssq*(u(i+1,j-1,k)*u(i+1,j-1,k) + v(i+1,j-1,k)*v(i+1,j-1,k) + w(i+1,j-1,k)*w(i+1,j-1,k))
	  udotc=(-u(i+1,j-1,k)+v(i+1,j-1,k))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i+1,j-1,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j-1,k)+qyy*pyy(i+1,j-1,k)-cssq*pzz(i+1,j-1,k)+two*qxy_9_10*pxy(i+1,j-1,k))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	endif
	!-1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
	
	!9  -1 +1
	if(isfluid(i-1,j+1,k).eq.0)then
	  udotc=(-u(i,j,k)+v(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+two*qxy_9_10*pxy(i,j,k))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	else
	  uu=halfonecssq*(u(i-1,j+1,k)*u(i-1,j+1,k) + v(i-1,j+1,k)*v(i-1,j+1,k) + w(i-1,j+1,k)*w(i-1,j+1,k))
	  udotc=(-u(i-1,j+1,k)+v(i-1,j+1,k))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i-1,j+1,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j+1,k)+qyy*pyy(i-1,j+1,k)-cssq*pzz(i-1,j+1,k)+two*qxy_9_10*pxy(i-1,j+1,k))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	endif
	!+1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop

	!15  -1  -1
	if(isfluid(i-1,j,k-1).eq.0)then
	  udotc=(u(i,j,k)+w(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+two*qxz_15_16*pxz(i,j,k))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i-1,j,k-1)*u(i-1,j,k-1) + v(i-1,j,k-1)*v(i-1,j,k-1) + w(i-1,j,k-1)*w(i-1,j,k-1))
	  udotc=(u(i-1,j,k-1)+w(i-1,j,k-1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i-1,j,k-1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k-1)+qzz*pzz(i-1,j,k-1)-cssq*pyy(i-1,j,k-1)+two*qxz_15_16*pxz(i-1,j,k-1))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq 
	endif
	!+1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pop
	
	!16  +1  +1
	if(isfluid(i+1,j,k+1).eq.0)then
	  udotc=(u(i,j,k)+w(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+two*qxz_15_16*pxz(i,j,k))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i+1,j,k+1)*u(i+1,j,k+1) + v(i+1,j,k+1)*v(i+1,j,k+1) + w(i+1,j,k+1)*w(i+1,j,k+1))
	  udotc=(u(i+1,j,k+1)+w(i+1,j,k+1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i+1,j,k+1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k+1)+qzz*pzz(i+1,j,k+1)-cssq*pyy(i+1,j,k+1)+two*qxz_15_16*pxz(i+1,j,k+1))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	endif
	!-1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz+temp_pop

	!17  +1   -1
	if(isfluid(i+1,j,k-1).eq.0)then
	  udotc=(-u(i,j,k)+w(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+two*qxz_17_18*pxz(i,j,k))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i+1,j,k-1)*u(i+1,j,k-1) + v(i+1,j,k-1)*v(i+1,j,k-1) + w(i+1,j,k-1)*w(i+1,j,k-1))
	  udotc=(-u(i+1,j,k-1)+w(i+1,j,k-1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i+1,j,k-1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k-1)+qzz*pzz(i+1,j,k-1)-cssq*pyy(i+1,j,k-1)+two*qxz_17_18*pxz(i+1,j,k-1))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	endif
	!-1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop
	
	!18   -1   +1
	if(isfluid(i-1,j,k+1).eq.0)then
	  udotc=(-u(i,j,k)+w(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+two*qxz_17_18*pxz(i,j,k))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	else
	  uu=halfonecssq*(u(i-1,j,k+1)*u(i-1,j,k+1) + v(i-1,j,k+1)*v(i-1,j,k+1) + w(i-1,j,k+1)*w(i-1,j,k+1))
	  udotc=(-u(i-1,j,k+1)+w(i-1,j,k+1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i-1,j,k+1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k+1)+qzz*pzz(i-1,j,k+1)-cssq*pyy(i-1,j,k+1)+two*qxz_17_18*pxz(i-1,j,k+1))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	endif
	!+1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!11  -1  -1
	if(isfluid(i,j-1,k-1).eq.0)then
	  udotc=(v(i,j,k)+w(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+two*qyz_11_12*pyz(i,j,k))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i,j-1,k-1)*u(i,j-1,k-1) + v(i,j-1,k-1)*v(i,j-1,k-1) + w(i,j-1,k-1)*w(i,j-1,k-1))
	  udotc=(v(i,j-1,k-1)+w(i,j-1,k-1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i,j-1,k-1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k-1)+qzz*pzz(i,j-1,k-1)-cssq*pxx(i,j-1,k-1)+two*qyz_11_12*pyz(i,j-1,k-1))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	endif
	! 0 +1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pop
	
	!12   +1  +1
	if(isfluid(i,j+1,k+1).eq.0)then
	  udotc=(v(i,j,k)+w(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+two*qyz_11_12*pyz(i,j,k))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i,j+1,k+1)*u(i,j+1,k+1) + v(i,j+1,k+1)*v(i,j+1,k+1) + w(i,j+1,k+1)*w(i,j+1,k+1))
	  udotc=(v(i,j+1,k+1)+w(i,j+1,k+1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i,j+1,k+1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k+1)+qzz*pzz(i,j+1,k+1)-cssq*pxx(i,j+1,k+1)+two*qyz_11_12*pyz(i,j+1,k+1))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	endif
	! 0 -1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz+temp_pop

	!13   -1   +1
	if(isfluid(i,j-1,k+1).eq.0)then
	  udotc=(v(i,j,k)-w(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+two*qyz_13_14*pyz(i,j,k))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	else
	  uu=halfonecssq*(u(i,j-1,k+1)*u(i,j-1,k+1) + v(i,j-1,k+1)*v(i,j-1,k+1) + w(i,j-1,k+1)*w(i,j-1,k+1))
	  udotc=(v(i,j-1,k+1)-w(i,j-1,k+1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i,j-1,k+1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k+1)+qzz*pzz(i,j-1,k+1)-cssq*pxx(i,j-1,k+1)+two*qyz_13_14*pyz(i,j-1,k+1))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	endif
	! 0 +1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	!14  +1 -1
	if(isfluid(i,j+1,k-1).eq.0)then
	  udotc=(v(i,j,k)-w(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+two*qyz_13_14*pyz(i,j,k))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i,j+1,k-1)*u(i,j+1,k-1) + v(i,j+1,k-1)*v(i,j+1,k-1) + w(i,j+1,k-1)*w(i,j+1,k-1))
	  udotc=(v(i,j+1,k-1)-w(i,j+1,k-1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i,j+1,k-1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k-1)+qzz*pzz(i,j+1,k-1)-cssq*pxx(i,j+1,k-1)+two*qyz_13_14*pyz(i,j+1,k-1))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	endif
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
    
  end subroutine streamcoll_bc
  
  attributes(global) subroutine streamcoll_bc_flop()
	
	implicit none  
	  
    integer :: i,j,k
    !logical :: alltrue
	real(kind=db) :: uu,udotc,temp,feq,uu0,fneq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
		  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
    
    !alltrue = .false.
    !if( allthreads(isfluid(i,j,k).ne.-1) )alltrue = .true.
    !call syncthreads
    !if(alltrue)return
    
    if(isfluid(i,j,k).ne.-1)return
    
	uu0=halfonecssq*(uh(i,j,k)*uh(i,j,k) + vh(i,j,k)*vh(i,j,k) + wh(i,j,k)*wh(i,j,k))
	!0
	temp_rho=p0*(rhoh(i,j,k)-uu0)
	temp_rho=temp_rho + oneminusomega*pi2cssq0*(-cssq*(pyyh(i,j,k)+pxxh(i,j,k)+pzzh(i,j,k)))
	
	!1
	if(isfluid(i-1,j,k).eq.0)then
	  udotc=uh(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k)-cssq*(pyyh(i,j,k)+pzzh(i,j,k)))
	  temp_pop=temp_pop - fx*p1dcssq
	else
	  uu=halfonecssq*(uh(i-1,j,k)*uh(i-1,j,k) + vh(i-1,j,k)*vh(i-1,j,k) + wh(i-1,j,k)*wh(i-1,j,k))
	  udotc=uh(i-1,j,k)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i-1,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i-1,j,k)-cssq*(pyyh(i-1,j,k)+pzzh(i-1,j,k)))
	  temp_pop=temp_pop + fx*p1dcssq
	endif
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop
	
	!2
	if(isfluid(i+1,j,k).eq.0)then
	  udotc=uh(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k)-cssq*(pyyh(i,j,k)+pzzh(i,j,k)))
	  temp_pop=temp_pop + fx*p1dcssq
	else
	  uu=halfonecssq*(uh(i+1,j,k)*uh(i+1,j,k) + vh(i+1,j,k)*vh(i+1,j,k) + wh(i+1,j,k)*wh(i+1,j,k))
	  udotc=uh(i+1,j,k)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i+1,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i+1,j,k)-cssq*(pyyh(i+1,j,k)+pzzh(i+1,j,k)))
	  temp_pop=temp_pop - fx*p1dcssq
	endif
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
	
	!3
	if(isfluid(i,j-1,k).eq.0)then
	  udotc=vh(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k)-cssq*(pxxh(i,j,k)+pzzh(i,j,k)))
	  temp_pop=temp_pop - fy*p1dcssq
	else
	  uu=halfonecssq*(uh(i,j-1,k)*uh(i,j-1,k) + vh(i,j-1,k)*vh(i,j-1,k) + wh(i,j-1,k)*wh(i,j-1,k))
	  udotc=vh(i,j-1,k)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j-1,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j-1,k)-cssq*(pxxh(i,j-1,k)+pzzh(i,j-1,k)))
	  temp_pop=temp_pop + fy*p1dcssq
	endif
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	if(isfluid(i,j+1,k).eq.0)then
	  udotc=vh(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k)-cssq*(pxxh(i,j,k)+pzzh(i,j,k)))
	  temp_pop=temp_pop + fy*p1dcssq
	else
	  uu=halfonecssq*(uh(i,j+1,k)*uh(i,j+1,k) + vh(i,j+1,k)*vh(i,j+1,k) + wh(i,j+1,k)*wh(i,j+1,k))
	  udotc=vh(i,j+1,k)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j+1,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j+1,k)-cssq*(pxxh(i,j+1,k)+pzzh(i,j+1,k)))
	  temp_pop=temp_pop - fy*p1dcssq
	endif
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	if(isfluid(i,j,k-1).eq.0)then
	  udotc=wh(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k)-cssq*(pxxh(i,j,k)+pyyh(i,j,k)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(uh(i,j,k-1)*uh(i,j,k-1) + vh(i,j,k-1)*vh(i,j,k-1) + wh(i,j,k-1)*wh(i,j,k-1))
	  udotc=wh(i,j,k-1)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k-1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k-1)-cssq*(pxxh(i,j,k-1)+pyyh(i,j,k-1)))
	  temp_pop=temp_pop + fz*p1dcssq
	endif
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop
	
	!6 +1
	if(isfluid(i,j,k+1).eq.0)then
	  udotc=wh(i,j,k)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k)-cssq*(pxxh(i,j,k)+pyyh(i,j,k)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(uh(i,j,k+1)*uh(i,j,k+1) + vh(i,j,k+1)*vh(i,j,k+1) + wh(i,j,k+1)*wh(i,j,k+1))
	  udotc=wh(i,j,k+1)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k+1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k+1)-cssq*(pxxh(i,j,k+1)+pyyh(i,j,k+1)))
	  temp_pop=temp_pop - fz*p1dcssq
	endif
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
	
	!7 -1 -1
	if(isfluid(i-1,j-1,k).eq.0)then
	  udotc=(uh(i,j,k)+vh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qyy*pyyh(i,j,k)-cssq*pzzh(i,j,k)+two*qxy_7_8*pxyh(i,j,k))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	else
	  uu=halfonecssq*(uh(i-1,j-1,k)*uh(i-1,j-1,k) + vh(i-1,j-1,k)*vh(i-1,j-1,k) + wh(i-1,j-1,k)*wh(i-1,j-1,k))
	  udotc=(uh(i-1,j-1,k)+vh(i-1,j-1,k))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i-1,j-1,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j-1,k)+qyy*pyyh(i-1,j-1,k)-cssq*pzzh(i-1,j-1,k)+two*qxy_7_8*pxyh(i-1,j-1,k))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	endif
	!+1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pop
	
	!8 +1 +1
	if(isfluid(i+1,j+1,k).eq.0)then
	  udotc=(uh(i,j,k)+vh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qyy*pyyh(i,j,k)-cssq*pzzh(i,j,k)+two*qxy_7_8*pxyh(i,j,k))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	else
	  uu=halfonecssq*(uh(i+1,j+1,k)*uh(i+1,j+1,k) + vh(i+1,j+1,k)*vh(i+1,j+1,k) + wh(i+1,j+1,k)*wh(i+1,j+1,k))
	  udotc=(uh(i+1,j+1,k)+vh(i+1,j+1,k))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i+1,j+1,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j+1,k)+qyy*pyyh(i+1,j+1,k)-cssq*pzzh(i+1,j+1,k)+two*qxy_7_8*pxyh(i+1,j+1,k))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	endif
	!-1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy+temp_pop
	
	!10   +1 -1
	if(isfluid(i+1,j-1,k).eq.0)then
	  udotc=(-uh(i,j,k)+vh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qyy*pyyh(i,j,k)-cssq*pzzh(i,j,k)+two*qxy_9_10*pxyh(i,j,k))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	else
	  uu=halfonecssq*(uh(i+1,j-1,k)*uh(i+1,j-1,k) + vh(i+1,j-1,k)*vh(i+1,j-1,k) + wh(i+1,j-1,k)*wh(i+1,j-1,k))
	  udotc=(-uh(i+1,j-1,k)+vh(i+1,j-1,k))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i+1,j-1,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j-1,k)+qyy*pyyh(i+1,j-1,k)-cssq*pzzh(i+1,j-1,k)+two*qxy_9_10*pxyh(i+1,j-1,k))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	endif
	!-1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
	
	!9  -1 +1
	if(isfluid(i-1,j+1,k).eq.0)then
	  udotc=(-uh(i,j,k)+vh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qyy*pyyh(i,j,k)-cssq*pzzh(i,j,k)+two*qxy_9_10*pxyh(i,j,k))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	else
	  uu=halfonecssq*(uh(i-1,j+1,k)*uh(i-1,j+1,k) + vh(i-1,j+1,k)*vh(i-1,j+1,k) + wh(i-1,j+1,k)*wh(i-1,j+1,k))
	  udotc=(-uh(i-1,j+1,k)+vh(i-1,j+1,k))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i-1,j+1,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j+1,k)+qyy*pyyh(i-1,j+1,k)-cssq*pzzh(i-1,j+1,k)+two*qxy_9_10*pxyh(i-1,j+1,k))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	endif
	!+1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop

	!15  -1  -1
	if(isfluid(i-1,j,k-1).eq.0)then
	  udotc=(uh(i,j,k)+wh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pyyh(i,j,k)+two*qxz_15_16*pxzh(i,j,k))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i-1,j,k-1)*uh(i-1,j,k-1) + vh(i-1,j,k-1)*vh(i-1,j,k-1) + wh(i-1,j,k-1)*wh(i-1,j,k-1))
	  udotc=(uh(i-1,j,k-1)+wh(i-1,j,k-1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i-1,j,k-1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k-1)+qzz*pzzh(i-1,j,k-1)-cssq*pyyh(i-1,j,k-1)+two*qxz_15_16*pxzh(i-1,j,k-1))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq 
	endif
	!+1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pop
	
	!16  +1  +1
	if(isfluid(i+1,j,k+1).eq.0)then
	  udotc=(uh(i,j,k)+wh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pyyh(i,j,k)+two*qxz_15_16*pxzh(i,j,k))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i+1,j,k+1)*uh(i+1,j,k+1) + vh(i+1,j,k+1)*vh(i+1,j,k+1) + wh(i+1,j,k+1)*wh(i+1,j,k+1))
	  udotc=(uh(i+1,j,k+1)+wh(i+1,j,k+1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i+1,j,k+1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k+1)+qzz*pzzh(i+1,j,k+1)-cssq*pyyh(i+1,j,k+1)+two*qxz_15_16*pxzh(i+1,j,k+1))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	endif
	!-1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz+temp_pop

	!17  +1   -1
	if(isfluid(i+1,j,k-1).eq.0)then
	  udotc=(-uh(i,j,k)+wh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pyyh(i,j,k)+two*qxz_17_18*pxzh(i,j,k))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i+1,j,k-1)*uh(i+1,j,k-1) + vh(i+1,j,k-1)*vh(i+1,j,k-1) + wh(i+1,j,k-1)*wh(i+1,j,k-1))
	  udotc=(-uh(i+1,j,k-1)+wh(i+1,j,k-1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i+1,j,k-1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k-1)+qzz*pzzh(i+1,j,k-1)-cssq*pyyh(i+1,j,k-1)+two*qxz_17_18*pxzh(i+1,j,k-1))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	endif
	!-1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop
	
	!18   -1   +1
	if(isfluid(i-1,j,k+1).eq.0)then
	  udotc=(-uh(i,j,k)+wh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pyyh(i,j,k)+two*qxz_17_18*pxzh(i,j,k))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	else
	  uu=halfonecssq*(uh(i-1,j,k+1)*uh(i-1,j,k+1) + vh(i-1,j,k+1)*vh(i-1,j,k+1) + wh(i-1,j,k+1)*wh(i-1,j,k+1))
	  udotc=(-uh(i-1,j,k+1)+wh(i-1,j,k+1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i-1,j,k+1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k+1)+qzz*pzzh(i-1,j,k+1)-cssq*pyyh(i-1,j,k+1)+two*qxz_17_18*pxzh(i-1,j,k+1))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	endif
	!+1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!11  -1  -1
	if(isfluid(i,j-1,k-1).eq.0)then
	  udotc=(vh(i,j,k)+wh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pxxh(i,j,k)+two*qyz_11_12*pyzh(i,j,k))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i,j-1,k-1)*uh(i,j-1,k-1) + vh(i,j-1,k-1)*vh(i,j-1,k-1) + wh(i,j-1,k-1)*wh(i,j-1,k-1))
	  udotc=(vh(i,j-1,k-1)+wh(i,j-1,k-1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j-1,k-1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k-1)+qzz*pzzh(i,j-1,k-1)-cssq*pxxh(i,j-1,k-1)+two*qyz_11_12*pyzh(i,j-1,k-1))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	endif
	! 0 +1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pop
	
	!12   +1  +1
	if(isfluid(i,j+1,k+1).eq.0)then
	  udotc=(vh(i,j,k)+wh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pxxh(i,j,k)+two*qyz_11_12*pyzh(i,j,k))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i,j+1,k+1)*uh(i,j+1,k+1) + vh(i,j+1,k+1)*vh(i,j+1,k+1) + wh(i,j+1,k+1)*wh(i,j+1,k+1))
	  udotc=(vh(i,j+1,k+1)+wh(i,j+1,k+1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j+1,k+1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k+1)+qzz*pzzh(i,j+1,k+1)-cssq*pxxh(i,j+1,k+1)+two*qyz_11_12*pyzh(i,j+1,k+1))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	endif
	! 0 -1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz+temp_pop

	!13   -1   +1
	if(isfluid(i,j-1,k+1).eq.0)then
	  udotc=(vh(i,j,k)-wh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pxxh(i,j,k)+two*qyz_13_14*pyzh(i,j,k))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	else
	  uu=halfonecssq*(uh(i,j-1,k+1)*uh(i,j-1,k+1) + vh(i,j-1,k+1)*vh(i,j-1,k+1) + wh(i,j-1,k+1)*wh(i,j-1,k+1))
	  udotc=(vh(i,j-1,k+1)-wh(i,j-1,k+1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j-1,k+1)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k+1)+qzz*pzzh(i,j-1,k+1)-cssq*pxxh(i,j-1,k+1)+two*qyz_13_14*pyzh(i,j-1,k+1))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	endif
	! 0 +1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	!14  +1 -1
	if(isfluid(i,j+1,k-1).eq.0)then
	  udotc=(vh(i,j,k)-wh(i,j,k))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k)+qzz*pzzh(i,j,k)-cssq*pxxh(i,j,k)+two*qyz_13_14*pyzh(i,j,k))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i,j+1,k-1)*uh(i,j+1,k-1) + vh(i,j+1,k-1)*vh(i,j+1,k-1) + wh(i,j+1,k-1)*wh(i,j+1,k-1))
	  udotc=(vh(i,j+1,k-1)-wh(i,j+1,k-1))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j+1,k-1)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k-1)+qzz*pzzh(i,j+1,k-1)-cssq*pxxh(i,j+1,k-1)+two*qyz_13_14*pyzh(i,j+1,k-1))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	endif
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
  
  end subroutine streamcoll_bc_flop
  
 end module streamcoll_bc_kernels
