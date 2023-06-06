#include "defines.h"
 module streamcoll_bc_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine streamcoll_bc()
	
	implicit none  
	  
    integer :: i,j,k,gi,gj,gk,idblock
    !logical :: alltrue
	real(kind=db) :: uu,udotc,temp,uu0
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
    
    integer :: idblock=1
    
	
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    !alltrue = .false.
    !if( allthreads(isfluid(gi,gj,gk).ne.-1) )alltrue = .true.
    !call syncthreads
    !if(alltrue)return
    
    if(isfluid(gi,gj,gk).ne.-1)return
    
	uu0=halfonecssq*(u(i,j,k,idblock)*u(i,j,k,idblock) + v(i,j,k,idblock)*v(i,j,k,idblock) + w(i,j,k,idblock)*w(i,j,k,idblock))
	!0
	temp_rho=p0*(rho(i,j,k,idblock)-uu0)
	temp_rho=temp_rho + oneminusomega*pi2cssq0*(-cssq*(pyy(i,j,k,idblock)+pxx(i,j,k,idblock)+pzz(i,j,k,idblock)))
	
	!1
	if(isfluid(gi-1,gj,gk).eq.0)then
	  udotc=u(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i,j,k,idblock)-cssq*(pyy(i,j,k,idblock)+pzz(i,j,k,idblock)))
	  temp_pop=temp_pop - fx*p1dcssq
	else
	  uu=halfonecssq*(u(i-1,j,k,idblock)*u(i-1,j,k,idblock) + v(i-1,j,k,idblock)*v(i-1,j,k,idblock) + w(i-1,j,k,idblock)*w(i-1,j,k,idblock))
	  udotc=u(i-1,j,k,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i-1,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i-1,j,k,idblock)-cssq*(pyy(i-1,j,k,idblock)+pzz(i-1,j,k,idblock)))
	  temp_pop=temp_pop + fx*p1dcssq
	endif
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop
	
	!2
	if(isfluid(gi+1,gj,gk).eq.0)then
	  udotc=u(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i,j,k,idblock)-cssq*(pyy(i,j,k,idblock)+pzz(i,j,k,idblock)))
	  temp_pop=temp_pop + fx*p1dcssq
	else
	  uu=halfonecssq*(u(i+1,j,k,idblock)*u(i+1,j,k,idblock) + v(i+1,j,k,idblock)*v(i+1,j,k,idblock) + w(i+1,j,k,idblock)*w(i+1,j,k,idblock))
	  udotc=u(i+1,j,k,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i+1,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxx(i+1,j,k,idblock)-cssq*(pyy(i+1,j,k,idblock)+pzz(i+1,j,k,idblock)))
	  temp_pop=temp_pop - fx*p1dcssq
	endif
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
	
	!3
	if(isfluid(gi,gj-1,gk).eq.0)then
	  udotc=v(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j,k,idblock)-cssq*(pxx(i,j,k,idblock)+pzz(i,j,k,idblock)))
	  temp_pop=temp_pop - fy*p1dcssq
	else
	  uu=halfonecssq*(u(i,j-1,k,idblock)*u(i,j-1,k,idblock) + v(i,j-1,k,idblock)*v(i,j-1,k,idblock) + w(i,j-1,k,idblock)*w(i,j-1,k,idblock))
	  udotc=v(i,j-1,k,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i,j-1,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j-1,k,idblock)-cssq*(pxx(i,j-1,k,idblock)+pzz(i,j-1,k,idblock)))
	  temp_pop=temp_pop + fy*p1dcssq
	endif
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	if(isfluid(i,j+1,k).eq.0)then
	  udotc=v(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j,k,idblock)-cssq*(pxx(i,j,k,idblock)+pzz(i,j,k,idblock)))
	  temp_pop=temp_pop + fy*p1dcssq
	else
	  uu=halfonecssq*(u(i,j+1,k,idblock)*u(i,j+1,k,idblock) + v(i,j+1,k,idblock)*v(i,j+1,k,idblock) + w(i,j+1,k,idblock)*w(i,j+1,k,idblock))
	  udotc=v(i,j+1,k,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i,j+1,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyy(i,j+1,k,idblock)-cssq*(pxx(i,j+1,k,idblock)+pzz(i,j+1,k,idblock)))
	  temp_pop=temp_pop - fy*p1dcssq
	endif
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	if(isfluid(i,j,k-1).eq.0)then
	  udotc=w(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k,idblock)-cssq*(pxx(i,j,k,idblock)+pyy(i,j,k,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(u(i,j,k-1,idblock)*u(i,j,k-1,idblock) + v(i,j,k-1,idblock)*v(i,j,k-1,idblock) + w(i,j,k-1,idblock)*w(i,j,k-1,idblock))
	  udotc=w(i,j,k-1,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k-1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k-1,idblock)-cssq*(pxx(i,j,k-1,idblock)+pyy(i,j,k-1,idblock)))
	  temp_pop=temp_pop + fz*p1dcssq
	endif
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop
	
	!6 +1
	if(isfluid(i,j,k+1).eq.0)then
	  udotc=w(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k,idblock)-cssq*(pxx(i,j,k,idblock)+pyy(i,j,k,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(u(i,j,k+1,idblock)*u(i,j,k+1,idblock) + v(i,j,k+1,idblock)*v(i,j,k+1,idblock) + w(i,j,k+1,idblock)*w(i,j,k+1,idblock))
	  udotc=w(i,j,k+1,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rho(i,j,k+1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzz(i,j,k+1,idblock)-cssq*(pxx(i,j,k+1,idblock)+pyy(i,j,k+1,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	endif
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
	
	!7 -1 -1
	if(isfluid(i-1,j-1,k).eq.0)then
	  udotc=(u(i,j,k,idblock)+v(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qyy*pyy(i,j,k,idblock)-cssq*pzz(i,j,k,idblock)+two*qxy_7_8*pxy(i,j,k,idblock))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	else
	  uu=halfonecssq*(u(i-1,j-1,k,idblock)*u(i-1,j-1,k,idblock) + v(i-1,j-1,k,idblock)*v(i-1,j-1,k,idblock) + w(i-1,j-1,k,idblock)*w(i-1,j-1,k,idblock))
	  udotc=(u(i-1,j-1,k,idblock)+v(i-1,j-1,k,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i-1,j-1,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j-1,k,idblock)+qyy*pyy(i-1,j-1,k,idblock)-cssq*pzz(i-1,j-1,k,idblock)+two*qxy_7_8*pxy(i-1,j-1,k,idblock))
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
	  udotc=(u(i,j,k,idblock)+v(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qyy*pyy(i,j,k,idblock)-cssq*pzz(i,j,k,idblock)+two*qxy_7_8*pxy(i,j,k,idblock))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	else
	  uu=halfonecssq*(u(i+1,j+1,k,idblock)*u(i+1,j+1,k,idblock) + v(i+1,j+1,k,idblock)*v(i+1,j+1,k,idblock) + w(i+1,j+1,k,idblock)*w(i+1,j+1,k,idblock))
	  udotc=(u(i+1,j+1,k,idblock)+v(i+1,j+1,k,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i+1,j+1,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j+1,k,idblock)+qyy*pyy(i+1,j+1,k,idblock)-cssq*pzz(i+1,j+1,k,idblock)+two*qxy_7_8*pxy(i+1,j+1,k,idblock))
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
	  udotc=(-u(i,j,k,idblock)+v(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qyy*pyy(i,j,k,idblock)-cssq*pzz(i,j,k,idblock)+two*qxy_9_10*pxy(i,j,k,idblock))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	else
	  uu=halfonecssq*(u(i+1,j-1,k,idblock)*u(i+1,j-1,k,idblock) + v(i+1,j-1,k,idblock)*v(i+1,j-1,k,idblock) + w(i+1,j-1,k,idblock)*w(i+1,j-1,k,idblock))
	  udotc=(-u(i+1,j-1,k,idblock)+v(i+1,j-1,k,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i+1,j-1,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j-1,k,idblock)+qyy*pyy(i+1,j-1,k,idblock)-cssq*pzz(i+1,j-1,k,idblock)+two*qxy_9_10*pxy(i+1,j-1,k,idblock))
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
	  udotc=(-u(i,j,k,idblock)+v(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qyy*pyy(i,j,k,idblock)-cssq*pzz(i,j,k,idblock)+two*qxy_9_10*pxy(i,j,k,idblock))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	else
	  uu=halfonecssq*(u(i-1,j+1,k,idblock)*u(i-1,j+1,k,idblock) + v(i-1,j+1,k,idblock)*v(i-1,j+1,k,idblock) + w(i-1,j+1,k,idblock)*w(i-1,j+1,k,idblock))
	  udotc=(-u(i-1,j+1,k,idblock)+v(i-1,j+1,k,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i-1,j+1,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j+1,k,idblock)+qyy*pyy(i-1,j+1,k,idblock)-cssq*pzz(i-1,j+1,k,idblock)+two*qxy_9_10*pxy(i-1,j+1,k,idblock))
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
	  udotc=(u(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pyy(i,j,k,idblock)+two*qxz_15_16*pxz(i,j,k,idblock))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i-1,j,k-1,idblock)*u(i-1,j,k-1,idblock) + v(i-1,j,k-1,idblock)*v(i-1,j,k-1,idblock) + w(i-1,j,k-1,idblock)*w(i-1,j,k-1,idblock))
	  udotc=(u(i-1,j,k-1,idblock)+w(i-1,j,k-1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i-1,j,k-1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k-1,idblock)+qzz*pzz(i-1,j,k-1,idblock)-cssq*pyy(i-1,j,k-1,idblock)+two*qxz_15_16*pxz(i-1,j,k-1,idblock))
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
	  udotc=(u(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pyy(i,j,k,idblock)+two*qxz_15_16*pxz(i,j,k,idblock))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i+1,j,k+1,idblock)*u(i+1,j,k+1,idblock) + v(i+1,j,k+1,idblock)*v(i+1,j,k+1,idblock) + w(i+1,j,k+1,idblock)*w(i+1,j,k+1,idblock))
	  udotc=(u(i+1,j,k+1,idblock)+w(i+1,j,k+1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i+1,j,k+1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k+1,idblock)+qzz*pzz(i+1,j,k+1,idblock)-cssq*pyy(i+1,j,k+1,idblock)+two*qxz_15_16*pxz(i+1,j,k+1,idblock))
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
	  udotc=(-u(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pyy(i,j,k,idblock)+two*qxz_17_18*pxz(i,j,k,idblock))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i+1,j,k-1,idblock)*u(i+1,j,k-1,idblock) + v(i+1,j,k-1,idblock)*v(i+1,j,k-1,idblock) + w(i+1,j,k-1,idblock)*w(i+1,j,k-1,idblock))
	  udotc=(-u(i+1,j,k-1,idblock)+w(i+1,j,k-1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i+1,j,k-1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k-1,idblock)+qzz*pzz(i+1,j,k-1,idblock)-cssq*pyy(i+1,j,k-1,idblock)+two*qxz_17_18*pxz(i+1,j,k-1,idblock))
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
	  udotc=(-u(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pyy(i,j,k,idblock)+two*qxz_17_18*pxz(i,j,k,idblock))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	else
	  uu=halfonecssq*(u(i-1,j,k+1,idblock)*u(i-1,j,k+1,idblock) + v(i-1,j,k+1,idblock)*v(i-1,j,k+1,idblock) + w(i-1,j,k+1,idblock)*w(i-1,j,k+1,idblock))
	  udotc=(-u(i-1,j,k+1,idblock)+w(i-1,j,k+1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i-1,j,k+1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k+1,idblock)+qzz*pzz(i-1,j,k+1,idblock)-cssq*pyy(i-1,j,k+1,idblock)+two*qxz_17_18*pxz(i-1,j,k+1,idblock))
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
	  udotc=(v(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pxx(i,j,k,idblock)+two*qyz_11_12*pyz(i,j,k,idblock))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i,j-1,k-1,idblock)*u(i,j-1,k-1,idblock) + v(i,j-1,k-1,idblock)*v(i,j-1,k-1,idblock) + w(i,j-1,k-1,idblock)*w(i,j-1,k-1,idblock))
	  udotc=(v(i,j-1,k-1,idblock)+w(i,j-1,k-1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i,j-1,k-1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k-1,idblock)+qzz*pzz(i,j-1,k-1,idblock)-cssq*pxx(i,j-1,k-1,idblock)+two*qyz_11_12*pyz(i,j-1,k-1,idblock))
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
	  udotc=(v(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pxx(i,j,k,idblock)+two*qyz_11_12*pyz(i,j,k,idblock))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i,j+1,k+1,idblock)*u(i,j+1,k+1,idblock) + v(i,j+1,k+1,idblock)*v(i,j+1,k+1,idblock) + w(i,j+1,k+1,idblock)*w(i,j+1,k+1,idblock))
	  udotc=(v(i,j+1,k+1,idblock)+w(i,j+1,k+1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i,j+1,k+1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k+1,idblock)+qzz*pzz(i,j+1,k+1,idblock)-cssq*pxx(i,j+1,k+1,idblock)+two*qyz_11_12*pyz(i,j+1,k+1,idblock))
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
	  udotc=(v(i,j,k,idblock)-w(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pxx(i,j,k,idblock)+two*qyz_13_14*pyz(i,j,k,idblock))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	else
	  uu=halfonecssq*(u(i,j-1,k+1,idblock)*u(i,j-1,k+1,idblock) + v(i,j-1,k+1,idblock)*v(i,j-1,k+1,idblock) + w(i,j-1,k+1,idblock)*w(i,j-1,k+1,idblock))
	  udotc=(v(i,j-1,k+1,idblock)-w(i,j-1,k+1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i,j-1,k+1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k+1,idblock)+qzz*pzz(i,j-1,k+1,idblock)-cssq*pxx(i,j-1,k+1,idblock)+two*qyz_13_14*pyz(i,j-1,k+1,idblock))
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
	  udotc=(v(i,j,k,idblock)-w(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rho(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pxx(i,j,k,idblock)+two*qyz_13_14*pyz(i,j,k,idblock))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	else
	  uu=halfonecssq*(u(i,j+1,k-1,idblock)*u(i,j+1,k-1,idblock) + v(i,j+1,k-1,idblock)*v(i,j+1,k-1,idblock) + w(i,j+1,k-1,idblock)*w(i,j+1,k-1,idblock))
	  udotc=(v(i,j+1,k-1,idblock)-w(i,j+1,k-1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rho(i,j+1,k-1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k-1,idblock)+qzz*pzz(i,j+1,k-1,idblock)-cssq*pxx(i,j+1,k-1,idblock)+two*qyz_13_14*pyz(i,j+1,k-1,idblock))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	endif
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
    
  end subroutine streamcoll_bc
  
  attributes(global) subroutine streamcoll_bc_flop()
	
	implicit none  
	  
    integer :: i,j,k
    !logical :: alltrue
	real(kind=db) :: uu,udotc,temp,feq,uu0,fneq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
	integer :: idblock=1
	
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
    
    !alltrue = .false.
    !if( allthreads(isfluid(i,j,k).ne.-1) )alltrue = .true.
    !call syncthreads
    !if(alltrue)return
    
    if(isfluid(i,j,k).ne.-1)return
    
	uu0=halfonecssq*(uh(i,j,k,idblock)*uh(i,j,k,idblock) + vh(i,j,k,idblock)*vh(i,j,k,idblock) + wh(i,j,k,idblock)*wh(i,j,k,idblock))
	!0
	temp_rho=p0*(rhoh(i,j,k,idblock)-uu0)
	temp_rho=temp_rho + oneminusomega*pi2cssq0*(-cssq*(pyyh(i,j,k,idblock)+pxxh(i,j,k,idblock)+pzzh(i,j,k,idblock)))
	
	!1
	if(isfluid(i-1,j,k).eq.0)then
	  udotc=uh(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k,idblock)-cssq*(pyyh(i,j,k,idblock)+pzzh(i,j,k,idblock)))
	  temp_pop=temp_pop - fx*p1dcssq
	else
	  uu=halfonecssq*(uh(i-1,j,k,idblock)*uh(i-1,j,k,idblock) + vh(i-1,j,k,idblock)*vh(i-1,j,k,idblock) + wh(i-1,j,k,idblock)*wh(i-1,j,k,idblock))
	  udotc=uh(i-1,j,k,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i-1,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i-1,j,k,idblock)-cssq*(pyyh(i-1,j,k,idblock)+pzzh(i-1,j,k,idblock)))
	  temp_pop=temp_pop + fx*p1dcssq
	endif
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop
	
	!2
	if(isfluid(i+1,j,k).eq.0)then
	  udotc=uh(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k,idblock)-cssq*(pyyh(i,j,k,idblock)+pzzh(i,j,k,idblock)))
	  temp_pop=temp_pop + fx*p1dcssq
	else
	  uu=halfonecssq*(uh(i+1,j,k,idblock)*uh(i+1,j,k,idblock) + vh(i+1,j,k,idblock)*vh(i+1,j,k,idblock) + wh(i+1,j,k,idblock)*wh(i+1,j,k,idblock))
	  udotc=uh(i+1,j,k,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i+1,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*pxxh(i+1,j,k,idblock)-cssq*(pyyh(i+1,j,k,idblock)+pzzh(i+1,j,k,idblock)))
	  temp_pop=temp_pop - fx*p1dcssq
	endif
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
	
	!3
	if(isfluid(i,j-1,k).eq.0)then
	  udotc=vh(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k,idblock)-cssq*(pxxh(i,j,k,idblock)+pzzh(i,j,k,idblock)))
	  temp_pop=temp_pop - fy*p1dcssq
	else
	  uu=halfonecssq*(uh(i,j-1,k,idblock)*uh(i,j-1,k,idblock) + vh(i,j-1,k,idblock)*vh(i,j-1,k,idblock) + wh(i,j-1,k,idblock)*wh(i,j-1,k,idblock))
	  udotc=vh(i,j-1,k,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j-1,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j-1,k,idblock)-cssq*(pxxh(i,j-1,k,idblock)+pzzh(i,j-1,k,idblock)))
	  temp_pop=temp_pop + fy*p1dcssq
	endif
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	if(isfluid(i,j+1,k).eq.0)then
	  udotc=vh(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k,idblock)-cssq*(pxxh(i,j,k,idblock)+pzzh(i,j,k,idblock)))
	  temp_pop=temp_pop + fy*p1dcssq
	else
	  uu=halfonecssq*(uh(i,j+1,k,idblock)*uh(i,j+1,k,idblock) + vh(i,j+1,k,idblock)*vh(i,j+1,k,idblock) + wh(i,j+1,k,idblock)*wh(i,j+1,k,idblock))
	  udotc=vh(i,j+1,k,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j+1,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*pyyh(i,j+1,k,idblock)-cssq*(pxxh(i,j+1,k,idblock)+pzzh(i,j+1,k,idblock)))
	  temp_pop=temp_pop - fy*p1dcssq
	endif
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	if(isfluid(i,j,k-1).eq.0)then
	  udotc=wh(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k,idblock)-cssq*(pxxh(i,j,k,idblock)+pyyh(i,j,k,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(uh(i,j,k-1,idblock)*uh(i,j,k-1,idblock) + vh(i,j,k-1,idblock)*vh(i,j,k-1,idblock) + wh(i,j,k-1,idblock)*wh(i,j,k-1,idblock))
	  udotc=wh(i,j,k-1,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k-1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k-1,idblock)-cssq*(pxxh(i,j,k-1,idblock)+pyyh(i,j,k-1,idblock)))
	  temp_pop=temp_pop + fz*p1dcssq
	endif
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop
	
	!6 +1
	if(isfluid(i,j,k+1).eq.0)then
	  udotc=wh(i,j,k,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k,idblock)-cssq*(pxxh(i,j,k,idblock)+pyyh(i,j,k,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(uh(i,j,k+1,idblock)*uh(i,j,k+1,idblock) + vh(i,j,k+1,idblock)*vh(i,j,k+1,idblock) + wh(i,j,k+1,idblock)*wh(i,j,k+1,idblock))
	  udotc=wh(i,j,k+1,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(rhoh(i,j,k+1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k+1,idblock)-cssq*(pxxh(i,j,k+1,idblock)+pyyh(i,j,k+1,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	endif
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
	
	!7 -1 -1
	if(isfluid(i-1,j-1,k).eq.0)then
	  udotc=(uh(i,j,k,idblock)+vh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qyy*pyyh(i,j,k,idblock)-cssq*pzzh(i,j,k,idblock)+two*qxy_7_8*pxyh(i,j,k,idblock))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	else
	  uu=halfonecssq*(uh(i-1,j-1,k,idblock)*uh(i-1,j-1,k,idblock) + vh(i-1,j-1,k,idblock)*vh(i-1,j-1,k,idblock) + wh(i-1,j-1,k,idblock)*wh(i-1,j-1,k,idblock))
	  udotc=(uh(i-1,j-1,k,idblock)+vh(i-1,j-1,k,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i-1,j-1,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j-1,k,idblock)+qyy*pyyh(i-1,j-1,k,idblock)-cssq*pzzh(i-1,j-1,k,idblock)+two*qxy_7_8*pxyh(i-1,j-1,k,idblock))
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
	  udotc=(uh(i,j,k,idblock)+vh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qyy*pyyh(i,j,k,idblock)-cssq*pzzh(i,j,k,idblock)+two*qxy_7_8*pxyh(i,j,k,idblock))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	else
	  uu=halfonecssq*(uh(i+1,j+1,k,idblock)*uh(i+1,j+1,k,idblock) + vh(i+1,j+1,k,idblock)*vh(i+1,j+1,k,idblock) + wh(i+1,j+1,k,idblock)*wh(i+1,j+1,k,idblock))
	  udotc=(uh(i+1,j+1,k,idblock)+vh(i+1,j+1,k,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i+1,j+1,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j+1,k,idblock)+qyy*pyyh(i+1,j+1,k,idblock)-cssq*pzzh(i+1,j+1,k,idblock)+two*qxy_7_8*pxyh(i+1,j+1,k,idblock))
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
	  udotc=(-uh(i,j,k,idblock)+vh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qyy*pyyh(i,j,k,idblock)-cssq*pzzh(i,j,k,idblock)+two*qxy_9_10*pxyh(i,j,k,idblock))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	else
	  uu=halfonecssq*(uh(i+1,j-1,k,idblock)*uh(i+1,j-1,k,idblock) + vh(i+1,j-1,k,idblock)*vh(i+1,j-1,k,idblock) + wh(i+1,j-1,k,idblock)*wh(i+1,j-1,k,idblock))
	  udotc=(-uh(i+1,j-1,k,idblock)+vh(i+1,j-1,k,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i+1,j-1,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j-1,k,idblock)+qyy*pyyh(i+1,j-1,k,idblock)-cssq*pzzh(i+1,j-1,k,idblock)+two*qxy_9_10*pxyh(i+1,j-1,k,idblock))
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
	  udotc=(-uh(i,j,k,idblock)+vh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qyy*pyyh(i,j,k,idblock)-cssq*pzzh(i,j,k,idblock)+two*qxy_9_10*pxyh(i,j,k,idblock))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	else
	  uu=halfonecssq*(uh(i-1,j+1,k,idblock)*uh(i-1,j+1,k,idblock) + vh(i-1,j+1,k,idblock)*vh(i-1,j+1,k,idblock) + wh(i-1,j+1,k,idblock)*wh(i-1,j+1,k,idblock))
	  udotc=(-uh(i-1,j+1,k,idblock)+vh(i-1,j+1,k,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i-1,j+1,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j+1,k,idblock)+qyy*pyyh(i-1,j+1,k,idblock)-cssq*pzzh(i-1,j+1,k,idblock)+two*qxy_9_10*pxyh(i-1,j+1,k,idblock))
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
	  udotc=(uh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pyyh(i,j,k,idblock)+two*qxz_15_16*pxzh(i,j,k,idblock))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i-1,j,k-1,idblock)*uh(i-1,j,k-1,idblock) + vh(i-1,j,k-1,idblock)*vh(i-1,j,k-1,idblock) + wh(i-1,j,k-1,idblock)*wh(i-1,j,k-1,idblock))
	  udotc=(uh(i-1,j,k-1,idblock)+wh(i-1,j,k-1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i-1,j,k-1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k-1,idblock)+qzz*pzzh(i-1,j,k-1,idblock)-cssq*pyyh(i-1,j,k-1,idblock)+two*qxz_15_16*pxzh(i-1,j,k-1,idblock))
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
	  udotc=(uh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pyyh(i,j,k,idblock)+two*qxz_15_16*pxzh(i,j,k,idblock))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i+1,j,k+1,idblock)*uh(i+1,j,k+1,idblock) + vh(i+1,j,k+1,idblock)*vh(i+1,j,k+1,idblock) + wh(i+1,j,k+1,idblock)*wh(i+1,j,k+1,idblock))
	  udotc=(uh(i+1,j,k+1,idblock)+wh(i+1,j,k+1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i+1,j,k+1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k+1,idblock)+qzz*pzzh(i+1,j,k+1,idblock)-cssq*pyyh(i+1,j,k+1,idblock)+two*qxz_15_16*pxzh(i+1,j,k+1,idblock))
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
	  udotc=(-uh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pyyh(i,j,k,idblock)+two*qxz_17_18*pxzh(i,j,k,idblock))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i+1,j,k-1,idblock)*uh(i+1,j,k-1,idblock) + vh(i+1,j,k-1,idblock)*vh(i+1,j,k-1,idblock) + wh(i+1,j,k-1,idblock)*wh(i+1,j,k-1,idblock))
	  udotc=(-uh(i+1,j,k-1,idblock)+wh(i+1,j,k-1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i+1,j,k-1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k-1,idblock)+qzz*pzzh(i+1,j,k-1,idblock)-cssq*pyyh(i+1,j,k-1,idblock)+two*qxz_17_18*pxzh(i+1,j,k-1,idblock))
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
	  udotc=(-uh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pyyh(i,j,k,idblock)+two*qxz_17_18*pxzh(i,j,k,idblock))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	else
	  uu=halfonecssq*(uh(i-1,j,k+1,idblock)*uh(i-1,j,k+1,idblock) + vh(i-1,j,k+1,idblock)*vh(i-1,j,k+1,idblock) + wh(i-1,j,k+1,idblock)*wh(i-1,j,k+1,idblock))
	  udotc=(-uh(i-1,j,k+1,idblock)+wh(i-1,j,k+1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i-1,j,k+1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k+1,idblock)+qzz*pzzh(i-1,j,k+1,idblock)-cssq*pyyh(i-1,j,k+1,idblock)+two*qxz_17_18*pxzh(i-1,j,k+1,idblock))
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
	  udotc=(vh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pxxh(i,j,k,idblock)+two*qyz_11_12*pyzh(i,j,k,idblock))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i,j-1,k-1,idblock)*uh(i,j-1,k-1,idblock) + vh(i,j-1,k-1,idblock)*vh(i,j-1,k-1,idblock) + wh(i,j-1,k-1,idblock)*wh(i,j-1,k-1,idblock))
	  udotc=(vh(i,j-1,k-1,idblock)+wh(i,j-1,k-1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j-1,k-1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k-1,idblock)+qzz*pzzh(i,j-1,k-1,idblock)-cssq*pxxh(i,j-1,k-1,idblock)+two*qyz_11_12*pyzh(i,j-1,k-1,idblock))
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
	  udotc=(vh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pxxh(i,j,k,idblock)+two*qyz_11_12*pyzh(i,j,k,idblock))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i,j+1,k+1,idblock)*uh(i,j+1,k+1,idblock) + vh(i,j+1,k+1,idblock)*vh(i,j+1,k+1,idblock) + wh(i,j+1,k+1,idblock)*wh(i,j+1,k+1,idblock))
	  udotc=(vh(i,j+1,k+1,idblock)+wh(i,j+1,k+1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j+1,k+1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k+1,idblock)+qzz*pzzh(i,j+1,k+1,idblock)-cssq*pxxh(i,j+1,k+1,idblock)+two*qyz_11_12*pyzh(i,j+1,k+1,idblock))
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
	  udotc=(vh(i,j,k,idblock)-wh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pxxh(i,j,k,idblock)+two*qyz_13_14*pyzh(i,j,k,idblock))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	else
	  uu=halfonecssq*(uh(i,j-1,k+1,idblock)*uh(i,j-1,k+1,idblock) + vh(i,j-1,k+1,idblock)*vh(i,j-1,k+1,idblock) + wh(i,j-1,k+1,idblock)*wh(i,j-1,k+1,idblock))
	  udotc=(vh(i,j-1,k+1,idblock)-wh(i,j-1,k+1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j-1,k+1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k+1,idblock)+qzz*pzzh(i,j-1,k+1,idblock)-cssq*pxxh(i,j-1,k+1,idblock)+two*qyz_13_14*pyzh(i,j-1,k+1,idblock))
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
	  udotc=(vh(i,j,k,idblock)-wh(i,j,k,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pxxh(i,j,k,idblock)+two*qyz_13_14*pyzh(i,j,k,idblock))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	else
	  uu=halfonecssq*(uh(i,j+1,k-1,idblock)*uh(i,j+1,k-1,idblock) + vh(i,j+1,k-1,idblock)*vh(i,j+1,k-1,idblock) + wh(i,j+1,k-1,idblock)*wh(i,j+1,k-1,idblock))
	  udotc=(vh(i,j+1,k-1,idblock)-wh(i,j+1,k-1,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(rhoh(i,j+1,k-1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k-1,idblock)+qzz*pzzh(i,j+1,k-1,idblock)-cssq*pxxh(i,j+1,k-1,idblock)+two*qyz_13_14*pyzh(i,j+1,k-1,idblock))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	endif
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
  
  end subroutine streamcoll_bc_flop
  
  attributes(global) subroutine streamcoll_bc_shared()
	
	implicit none  
	  
    integer :: i,j,k,gi,gj,gk,idblock
	real(kind=db) :: udotc,uu,temp
    integer :: idblock=1
    
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
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    
    
    uu=halfonecssq*(u(i,j,k,idblock)*u(i,j,k,idblock) + v(i,j,k,idblock)*v(i,j,k,idblock) + w(i,j,k,idblock)*w(i,j,k,idblock))
    
    !0
	f00(li,lj,lk)=p0*(rho(i,j,k,idblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(pyy(i,j,k,idblock)+pxx(i,j,k,idblock)+pzz(i,j,k,idblock)))
	
    
	!1 -1  0  0
	udotc=u(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxx(i,j,k,idblock)-cssq*(pyy(i,j,k,idblock)+pzz(i,j,k,idblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxx(i,j,k,idblock)-cssq*(pyy(i,j,k,idblock)+pzz(i,j,k,idblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=v(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyy(i,j,k,idblock)-cssq*(pxx(i,j,k,idblock)+pzz(i,j,k,idblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyy(i,j,k,idblock)-cssq*(pxx(i,j,k,idblock)+pzz(i,j,k,idblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=w(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k,idblock)-cssq*(pxx(i,j,k,idblock)+pyy(i,j,k,idblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k,idblock)-cssq*(pxx(i,j,k,idblock)+pyy(i,j,k,idblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(u(i,j,k,idblock)+v(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qyy*pyy(i,j,k,idblock)-cssq*pzz(i,j,k,idblock)+two*qxy_7_8*pxy(i,j,k,idblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qyy*pyy(i,j,k,idblock)-cssq*pzz(i,j,k,idblock)+two*qxy_7_8*pxy(i,j,k,idblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-u(i,j,k,idblock)+v(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qyy*pyy(i,j,k,idblock)-cssq*pzz(i,j,k,idblock)+two*qxy_9_10*pxy(i,j,k,idblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qyy*pyy(i,j,k,idblock)-cssq*pzz(i,j,k,idblock)+two*qxy_9_10*pxy(i,j,k,idblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(u(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pyy(i,j,k,idblock)+two*qxz_15_16*pxz(i,j,k,idblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pyy(i,j,k,idblock)+two*qxz_15_16*pxz(i,j,k,idblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-u(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pyy(i,j,k,idblock)+two*qxz_17_18*pxz(i,j,k,idblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pyy(i,j,k,idblock)+two*qxz_17_18*pxz(i,j,k,idblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(v(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pxx(i,j,k,idblock)+two*qyz_11_12*pyz(i,j,k,idblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pxx(i,j,k,idblock)+two*qyz_11_12*pyz(i,j,k,idblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(v(i,j,k,idblock)-w(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pxx(i,j,k,idblock)+two*qyz_13_14*pyz(i,j,k,idblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k,idblock)+qzz*pzz(i,j,k,idblock)-cssq*pxx(i,j,k,idblock)+two*qyz_13_14*pyz(i,j,k,idblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
        ! Halo Faces
    if(li==1)then
      
      uu=halfonecssq*(u(i-1,j,k,idblock)*u(i-1,j,k,idblock) + v(i-1,j,k,idblock)*v(i-1,j,k,idblock) + w(i-1,j,k,idblock)*w(i-1,j,k,idblock))
      
      !1 -1  0  0
	  udotc=u(i-1,j,k,idblock)*onecssq
	  f01(li-1,lj,lk)=p1*(rho(i-1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxx(i-1,j,k,idblock)-cssq*(pyy(i-1,j,k,idblock)+pzz(i-1,j,k,idblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(u(i-1,j,k,idblock)+v(i-1,j,k,idblock))*onecssq
	  f07(li-1,lj,lk)=p2*(rho(i-1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k,idblock)+qyy*pyy(i-1,j,k,idblock)-cssq*pzz(i-1,j,k,idblock)+two*qxy_7_8*pxy(i-1,j,k,idblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-u(i-1,j,k,idblock)+v(i-1,j,k,idblock))*onecssq
	  f09(li-1,lj,lk)=p2*(rho(i-1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k,idblock)+qyy*pyy(i-1,j,k,idblock)-cssq*pzz(i-1,j,k,idblock)+two*qxy_9_10*pxy(i-1,j,k,idblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(u(i-1,j,k,idblock)+w(i-1,j,k,idblock))*onecssq
	  f15(li-1,lj,lk)=p2*(rho(i-1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k,idblock)+qzz*pzz(i-1,j,k,idblock)-cssq*pyy(i-1,j,k,idblock)+two*qxz_15_16*pxz(i-1,j,k,idblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-u(i-1,j,k,idblock)+w(i-1,j,k,idblock))*onecssq
	  f18(li-1,lj,lk)=p2*(rho(i-1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k,idblock)+qzz*pzz(i-1,j,k,idblock)-cssq*pyy(i-1,j,k,idblock)+two*qxz_17_18*pxz(i-1,j,k,idblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      uu=halfonecssq*(u(i+1,j,k,idblock)*u(i+1,j,k,idblock) + v(i+1,j,k,idblock)*v(i+1,j,k,idblock) + w(i+1,j,k,idblock)*w(i+1,j,k,idblock))
      
      !2 +1  0  0
	  udotc=u(i+1,j,k,idblock)*onecssq
	  f02(li+1,lj,lk)=p1*(rho(i+1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxx(i+1,j,k,idblock)-cssq*(pyy(i+1,j,k,idblock)+pzz(i+1,j,k,idblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(u(i+1,j,k,idblock)+v(i+1,j,k,idblock))*onecssq
	  f08(li+1,lj,lk)=p2*(rho(i+1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k,idblock)+qyy*pyy(i+1,j,k,idblock)-cssq*pzz(i+1,j,k,idblock)+two*qxy_7_8*pxy(i+1,j,k,idblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-u(i+1,j,k,idblock)+v(i+1,j,k,idblock))*onecssq
	  f10(li+1,lj,lk)=p2*(rho(i+1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k,idblock)+qyy*pyy(i+1,j,k,idblock)-cssq*pzz(i+1,j,k,idblock)+two*qxy_9_10*pxy(i+1,j,k,idblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(u(i+1,j,k,idblock)+w(i+1,j,k,idblock))*onecssq
	  f16(li+1,lj,lk)=p2*(rho(i+1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k,idblock)+qzz*pzz(i+1,j,k,idblock)-cssq*pyy(i+1,j,k,idblock)+two*qxz_15_16*pxz(i+1,j,k,idblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-u(i+1,j,k,idblock)+w(i+1,j,k,idblock))*onecssq
	  f17(li+1,lj,lk)=p2*(rho(i+1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k,idblock)+qzz*pzz(i+1,j,k,idblock)-cssq*pyy(i+1,j,k,idblock)+two*qxz_17_18*pxz(i+1,j,k,idblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      uu=halfonecssq*(u(i,j-1,k,idblock)*u(i,j-1,k,idblock) + v(i,j-1,k,idblock)*v(i,j-1,k,idblock) + w(i,j-1,k,idblock)*w(i,j-1,k,idblock))
      
            !3 0 -1  0
	  udotc=v(i,j-1,k,idblock)*onecssq
	  f03(li,lj-1,lk)=p1*(rho(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyy(i,j-1,k,idblock)-cssq*(pxx(i,j-1,k,idblock)+pzz(i,j-1,k,idblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(u(i,j-1,k,idblock)+v(i,j-1,k,idblock))*onecssq
	  f07(li,lj-1,lk)=p2*(rho(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j-1,k,idblock)+qyy*pyy(i,j-1,k,idblock)-cssq*pzz(i,j-1,k,idblock)+two*qxy_7_8*pxy(i,j-1,k,idblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-u(i,j-1,k,idblock)+v(i,j-1,k,idblock))*onecssq
	  f10(li,lj-1,lk)=p2*(rho(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j-1,k,idblock)+qyy*pyy(i,j-1,k,idblock)-cssq*pzz(i,j-1,k,idblock)+two*qxy_9_10*pxy(i,j-1,k,idblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(v(i,j-1,k,idblock)+w(i,j-1,k,idblock))*onecssq
	  f11(li,lj-1,lk)=p2*(rho(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k,idblock)+qzz*pzz(i,j-1,k,idblock)-cssq*pxx(i,j-1,k,idblock)+two*qyz_11_12*pyz(i,j-1,k,idblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(v(i,j-1,k,idblock)-w(i,j-1,k,idblock))*onecssq
	  f13(li,lj-1,lk)=p2*(rho(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k,idblock)+qzz*pzz(i,j-1,k,idblock)-cssq*pxx(i,j-1,k,idblock)+two*qyz_13_14*pyz(i,j-1,k,idblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(u(i,j+1,k,idblock)*u(i,j+1,k,idblock) + v(i,j+1,k,idblock)*v(i,j+1,k,idblock) + w(i,j+1,k,idblock)*w(i,j+1,k,idblock))
      
      !4  0 +1  0
	  udotc=v(i,j+1,k,idblock)*onecssq
	  f04(li,lj+1,lk)=p1*(rho(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyy(i,j+1,k,idblock)-cssq*(pxx(i,j+1,k,idblock)+pzz(i,j+1,k,idblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(u(i,j+1,k,idblock)+v(i,j+1,k,idblock))*onecssq
	  f08(li,lj+1,lk)=p2*(rho(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j+1,k,idblock)+qyy*pyy(i,j+1,k,idblock)-cssq*pzz(i,j+1,k,idblock)+two*qxy_7_8*pxy(i,j+1,k,idblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-u(i,j+1,k,idblock)+v(i,j+1,k,idblock))*onecssq
	  f09(li,lj+1,lk)=p2*(rho(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j+1,k,idblock)+qyy*pyy(i,j+1,k,idblock)-cssq*pzz(i,j+1,k,idblock)+two*qxy_9_10*pxy(i,j+1,k,idblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(v(i,j+1,k,idblock)+w(i,j+1,k,idblock))*onecssq
	  f12(li,lj+1,lk)=p2*(rho(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k,idblock)+qzz*pzz(i,j+1,k,idblock)-cssq*pxx(i,j+1,k,idblock)+two*qyz_11_12*pyz(i,j+1,k,idblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(v(i,j+1,k,idblock)-w(i,j+1,k,idblock))*onecssq
	  f14(li,lj+1,lk)=p2*(rho(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k,idblock)+qzz*pzz(i,j+1,k,idblock)-cssq*pxx(i,j+1,k,idblock)+two*qyz_13_14*pyz(i,j+1,k,idblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      uu=halfonecssq*(u(i,j,k-1,idblock)*u(i,j,k-1,idblock) + v(i,j,k-1,idblock)*v(i,j,k-1,idblock) + w(i,j,k-1,idblock)*w(i,j,k-1,idblock))
      
      !5  0  0 -1
	  udotc=w(i,j,k-1,idblock)*onecssq
	  f05(li,lj,lk-1)=p1*(rho(i,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k-1,idblock)-cssq*(pxx(i,j,k-1,idblock)+pyy(i,j,k-1,idblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(u(i,j,k-1,idblock)+w(i,j,k-1,idblock))*onecssq
	  f15(li,lj,lk-1)=p2*(rho(i,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k-1,idblock)+qzz*pzz(i,j,k-1,idblock)-cssq*pyy(i,j,k-1,idblock)+two*qxz_15_16*pxz(i,j,k-1,idblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-u(i,j,k-1,idblock)+w(i,j,k-1,idblock))*onecssq
	  f17(li,lj,lk-1)=p2*(rho(i,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k-1,idblock)+qzz*pzz(i,j,k-1,idblock)-cssq*pyy(i,j,k-1,idblock)+two*qxz_17_18*pxz(i,j,k-1,idblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(v(i,j,k-1,idblock)+w(i,j,k-1,idblock))*onecssq
	  f11(li,lj,lk-1)=p2*(rho(i,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k-1,idblock)+qzz*pzz(i,j,k-1,idblock)-cssq*pxx(i,j,k-1,idblock)+two*qyz_11_12*pyz(i,j,k-1,idblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(v(i,j,k-1,idblock)-w(i,j,k-1,idblock))*onecssq
	  f14(li,lj,lk-1)=p2*(rho(i,j,k-1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k-1,idblock)+qzz*pzz(i,j,k-1,idblock)-cssq*pxx(i,j,k-1,idblock)+two*qyz_13_14*pyz(i,j,k-1,idblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i,j,k+1,idblock)*u(i,j,k+1,idblock) + v(i,j,k+1,idblock)*v(i,j,k+1,idblock) + w(i,j,k+1,idblock)*w(i,j,k+1,idblock))
      
      !6  0  0  +1
	  udotc=w(i,j,k+1,idblock)*onecssq
	  f06(li,lj,lk+1)=p1*(rho(i,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzz(i,j,k+1,idblock)-cssq*(pxx(i,j,k+1,idblock)+pyy(i,j,k+1,idblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(u(i,j,k+1,idblock)+w(i,j,k+1,idblock))*onecssq
	  f16(li,lj,lk+1)=p2*(rho(i,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k+1,idblock)+qzz*pzz(i,j,k+1,idblock)-cssq*pyy(i,j,k+1,idblock)+two*qxz_15_16*pxz(i,j,k+1,idblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-u(i,j,k+1,idblock)+w(i,j,k+1,idblock))*onecssq
	  f18(li,lj,lk+1)=p2*(rho(i,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i,j,k+1,idblock)+qzz*pzz(i,j,k+1,idblock)-cssq*pyy(i,j,k+1,idblock)+two*qxz_17_18*pxz(i,j,k+1,idblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(v(i,j,k+1,idblock)+w(i,j,k+1,idblock))*onecssq
	  f12(li,lj,lk+1)=p2*(rho(i,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k+1,idblock)+qzz*pzz(i,j,k+1,idblock)-cssq*pxx(i,j,k+1,idblock)+two*qyz_11_12*pyz(i,j,k+1,idblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(v(i,j,k+1,idblock)-w(i,j,k+1,idblock))*onecssq
	  f13(li,lj,lk+1)=p2*(rho(i,j,k+1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j,k+1,idblock)+qzz*pzz(i,j,k+1,idblock)-cssq*pxx(i,j,k+1,idblock)+two*qyz_13_14*pyz(i,j,k+1,idblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      uu=halfonecssq*(u(i-1,j-1,k,idblock)*u(i-1,j-1,k,idblock) + v(i-1,j-1,k,idblock)*v(i-1,j-1,k,idblock) + w(i-1,j-1,k,idblock)*w(i-1,j-1,k,idblock))
      
      !7 -1 -1  0
	  udotc=(u(i-1,j-1,k,idblock)+v(i-1,j-1,k,idblock))*onecssq
	  f07(li-1,lj-1,lk)=p2*(rho(i-1,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j-1,k,idblock)+qyy*pyy(i-1,j-1,k,idblock)-cssq*pzz(i-1,j-1,k,idblock)+two*qxy_7_8*pxy(i-1,j-1,k,idblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(u(i-1,j+1,k,idblock)*u(i-1,j+1,k,idblock) + v(i-1,j+1,k,idblock)*v(i-1,j+1,k,idblock) + w(i-1,j+1,k,idblock)*w(i-1,j+1,k,idblock))
      
      !9  -1 +1 0
      udotc=(-u(i-1,j+1,k,idblock)+v(i-1,j+1,k,idblock))*onecssq
	  f09(li-1,lj+1,lk)=p2*(rho(i-1,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j+1,k,idblock)+qyy*pyy(i-1,j+1,k,idblock)-cssq*pzz(i-1,j+1,k,idblock)+two*qxy_9_10*pxy(i-1,j+1,k,idblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      uu=halfonecssq*(u(i-1,j,k-1,idblock)*u(i-1,j,k-1,idblock) + v(i-1,j,k-1,idblock)*v(i-1,j,k-1,idblock) + w(i-1,j,k-1,idblock)*w(i-1,j,k-1,idblock))
      
      !15  -1  0 -1
	  udotc=(u(i-1,j,k-1,idblock)+w(i-1,j,k-1,idblock))*onecssq
	  f15(li-1,lj,lk-1)=p2*(rho(i-1,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k-1,idblock)+qzz*pzz(i-1,j,k-1,idblock)-cssq*pyy(i-1,j,k-1,idblock)+two*qxz_15_16*pxz(i-1,j,k-1,idblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i-1,j,k+1,idblock)*u(i-1,j,k+1,idblock) + v(i-1,j,k+1,idblock)*v(i-1,j,k+1,idblock) + w(i-1,j,k+1,idblock)*w(i-1,j,k+1,idblock))
      
      !18   -1   0  +1
	  udotc=(-u(i-1,j,k+1,idblock)+w(i-1,j,k+1,idblock))*onecssq
	  f18(li-1,lj,lk+1)=p2*(rho(i-1,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i-1,j,k+1,idblock)+qzz*pzz(i-1,j,k+1,idblock)-cssq*pyy(i-1,j,k+1,idblock)+two*qxz_17_18*pxz(i-1,j,k+1,idblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      uu=halfonecssq*(u(i+1,j-1,k,idblock)*u(i+1,j-1,k,idblock) + v(i+1,j-1,k,idblock)*v(i+1,j-1,k,idblock) + w(i+1,j-1,k,idblock)*w(i+1,j-1,k,idblock))
      
      !10   +1 -1  0
	  udotc=(-u(i+1,j-1,k,idblock)+v(i+1,j-1,k,idblock))*onecssq
	  f10(li+1,lj-1,lk)=p2*(rho(i+1,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j-1,k,idblock)+qyy*pyy(i+1,j-1,k,idblock)-cssq*pzz(i+1,j-1,k,idblock)+two*qxy_9_10*pxy(i+1,j-1,k,idblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(u(i+1,j+1,k,idblock)*u(i+1,j+1,k,idblock) + v(i+1,j+1,k,idblock)*v(i+1,j+1,k,idblock) + w(i+1,j+1,k,idblock)*w(i+1,j+1,k,idblock))
      
      !8 +1 +1  0
	  udotc=(u(i+1,j+1,k,idblock)+v(i+1,j+1,k,idblock))*onecssq
	  f08(li+1,lj+1,lk)=p2*(rho(i+1,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j+1,k,idblock)+qyy*pyy(i+1,j+1,k,idblock)-cssq*pzz(i+1,j+1,k,idblock)+two*qxy_7_8*pxy(i+1,j+1,k,idblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      uu=halfonecssq*(u(i+1,j,k-1,idblock)*u(i+1,j,k-1,idblock) + v(i+1,j,k-1,idblock)*v(i+1,j,k-1,idblock) + w(i+1,j,k-1,idblock)*w(i+1,j,k-1,idblock))
      
      !17  +1  0 -1
	  udotc=(-u(i+1,j,k-1,idblock)+w(i+1,j,k-1,idblock))*onecssq
	  f17(li+1,lj,lk-1)=p2*(rho(i+1,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k-1,idblock)+qzz*pzz(i+1,j,k-1,idblock)-cssq*pyy(i+1,j,k-1,idblock)+two*qxz_17_18*pxz(i+1,j,k-1,idblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i+1,j,k+1,idblock)*u(i+1,j,k+1,idblock) + v(i+1,j,k+1,idblock)*v(i+1,j,k+1,idblock) + w(i+1,j,k+1,idblock)*w(i+1,j,k+1,idblock))
      
      !16  +1  0 +1
	  udotc=(u(i+1,j,k+1,idblock)+w(i+1,j,k+1,idblock))*onecssq
	  f16(li+1,lj,lk+1)=p2*(rho(i+1,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxx(i+1,j,k+1,idblock)+qzz*pzz(i+1,j,k+1,idblock)-cssq*pyy(i+1,j,k+1,idblock)+two*qxz_15_16*pxz(i+1,j,k+1,idblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      uu=halfonecssq*(u(i,j-1,k-1,idblock)*u(i,j-1,k-1,idblock) + v(i,j-1,k-1,idblock)*v(i,j-1,k-1,idblock) + w(i,j-1,k-1,idblock)*w(i,j-1,k-1,idblock))
      
      !11  0  -1  -1
	  udotc=(v(i,j-1,k-1,idblock)+w(i,j-1,k-1,idblock))*onecssq
	  f11(li,lj-1,lk-1)=p2*(rho(i,j-1,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k-1,idblock)+qzz*pzz(i,j-1,k-1,idblock)-cssq*pxx(i,j-1,k-1,idblock)+two*qyz_11_12*pyz(i,j-1,k-1,idblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i,j-1,k+1,idblock)*u(i,j-1,k+1,idblock) + v(i,j-1,k+1,idblock)*v(i,j-1,k+1,idblock) + w(i,j-1,k+1,idblock)*w(i,j-1,k+1,idblock))
      
      !13  0  -1   +1
	  udotc=(v(i,j-1,k+1,idblock)-w(i,j-1,k+1,idblock))*onecssq
	  f13(li,lj-1,lk+1)=p2*(rho(i,j-1,k+1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j-1,k+1,idblock)+qzz*pzz(i,j-1,k+1,idblock)-cssq*pxx(i,j-1,k+1,idblock)+two*qyz_13_14*pyz(i,j-1,k+1,idblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      uu=halfonecssq*(u(i,j+1,k-1,idblock)*u(i,j+1,k-1,idblock) + v(i,j+1,k-1,idblock)*v(i,j+1,k-1,idblock) + w(i,j+1,k-1,idblock)*w(i,j+1,k-1,idblock))
      
      !14  0  +1  -1
	  udotc=(v(i,j+1,k-1,idblock)-w(i,j+1,k-1,idblock))*onecssq
	  f14(li,lj+1,lk-1)=p2*(rho(i,j+1,k-1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k-1,idblock)+qzz*pzz(i,j+1,k-1,idblock)-cssq*pxx(i,j+1,k-1,idblock)+two*qyz_13_14*pyz(i,j+1,k-1,idblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(u(i,j+1,k+1,idblock)*u(i,j+1,k+1,idblock) + v(i,j+1,k+1,idblock)*v(i,j+1,k+1,idblock) + w(i,j+1,k+1,idblock)*w(i,j+1,k+1,idblock))
      
      !12  0  +1  +1
	  udotc=(v(i,j+1,k+1,idblock)+w(i,j+1,k+1,idblock))*onecssq
	  f12(li,lj+1,lk+1)=p2*(rho(i,j+1,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyy(i,j+1,k+1,idblock)+qzz*pzz(i,j+1,k+1,idblock)-cssq*pxx(i,j+1,k+1,idblock)+two*qyz_11_12*pyz(i,j+1,k+1,idblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
    
    if(isfluid(gi,gj,gk).eq.0)then
#ifdef IFBC
          
          if(isfluid(gi+1,gj,gk-1).eq.-1)f18(li,lj,lk)=f17(li+1,lj,lk-1) !gpc 
          if(isfluid(gi-1,gj,gk+1).eq.-1)f17(li,lj,lk)=f18(li-1,lj,lk+1) !hpc

          if(isfluid(gi-1,gj,gk-1).eq.-1)f16(li,lj,lk)=f15(li-1,lj,lk-1) !gpc 
          if(isfluid(gi+1,gj,gk+1).eq.-1)f15(li,lj,lk)=f16(li+1,lj,lk+1) !hpc

          if(isfluid(gi,gj-1,gk+1).eq.-1)f14(li,lj,lk)=f13(li,lj-1,lk+1)!gpc 
          if(isfluid(gi,gj+1,gk-1).eq.-1)f13(li,lj,lk)=f14(li,lj+1,lk-1)!hpc
        
          if(isfluid(i,j-1,k-1).eq.-1)f12(li,lj,lk)=f11(li,lj-1,lk-1)!gpc
          if(isfluid(i,j+1,k+1).eq.-1)f11(li,lj,lk)=f12(li,lj+1,lk+1)!hpc

          if(isfluid(i-1,j+1,k).eq.-1)f10(li,lj,lk)=f09(li-1,lj+1,lk)!gpc  
          if(isfluid(i+1,j-1,k).eq.-1)f09(li,lj,lk)=f10(li+1,lj-1,lk)!hpc

          if(isfluid(i-1,j-1,k).eq.-1)f08(li,lj,lk)=f07(li-1,lj-1,lk)!gpc 
          if(isfluid(i+1,j+1,k).eq.-1)f07(li,lj,lk)=f08(li+1,lj+1,lk)!hpc

          if(isfluid(i,j,k-1).eq.-1)f06(li,lj,lk)=f05(li,lj,lk-1)!gpc
          if(isfluid(i,j,k+1).eq.-1)f05(li,lj,lk)=f06(li,lj,lk+1)!hpc 

          if(isfluid(i,j-1,k).eq.-1)f04(li,lj,lk)=f03(li,lj-1,lk)!gpc 
          if(isfluid(i,j+1,k).eq.-1)f03(li,lj,lk)=f04(li,lj+1,lk)!hpc 

          if(isfluid(i-1,j,k).eq.-1)f02(li,lj,lk)=f01(li-1,lj,lk)!gpc
          if(isfluid(i+1,j,k).eq.-1)f01(li,lj,lk)=f02(li+1,lj,lk)!hpc 
          
#else
          
          f18(li,lj,lk)=f17(li+1,lj,lk-1) !gpc 
          f17(li,lj,lk)=f18(li-1,lj,lk+1) !hpc

          f16(li,lj,lk)=f15(li-1,lj,lk-1) !gpc 
          f15(li,lj,lk)=f16(li+1,lj,lk+1) !hpc

          f14(li,lj,lk)=f13(li,lj-1,lk+1)!gpc 
          f13(li,lj,lk)=f14(li,lj+1,lk-1)!hpc
        
          f12(li,lj,lk)=f11(li,lj-1,lk-1)!gpc 
          f11(li,lj,lk)=f12(li,lj+1,lk+1)!hpc

          f10(li,lj,lk)=f09(li-1,lj+1,lk)!gpc 
          f09(li,lj,lk)=f10(li+1,lj-1,lk)!hpc

          f08(li,lj,lk)=f07(li-1,lj-1,lk)!gpc 
          f07(li,lj,lk)=f08(li+1,lj+1,lk)!hpc

          f06(li,lj,lk)=f05(li,lj,lk-1)!gpc 
          f05(li,lj,lk)=f06(li,lj,lk+1)!hpc 

          f04(li,lj,lk)=f03(li,lj-1,lk)!gpc 
          f03(li,lj,lk)=f04(li,lj+1,lk)!hpc 

          f02(li,lj,lk)=f01(li-1,lj,lk)!gpc 
          f01(li,lj,lk)=f02(li+1,lj,lk)!hpc          

#endif
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
	rhoh(i,j,k,idblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	uh(i,j,k,idblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	vh(i,j,k,idblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	wh(i,j,k,idblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pxxh(i,j,k,idblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	pyyh(i,j,k,idblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pzzh(i,j,k,idblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	pxyh(i,j,k,idblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	pxzh(i,j,k,idblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	pyzh(i,j,k,idblock)=udotc
     
#ifdef PRESSCORR

	call syncthreads
	
	uu=halfonecssq*(uh(i,j,k,idblock)*uh(i,j,k,idblock) + vh(i,j,k,idblock)*vh(i,j,k,idblock) + wh(i,j,k,idblock)*wh(i,j,k,idblock))
    
	!1 -1  0  0
	udotc=uh(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=vh(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=wh(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(uh(i,j,k,idblock)+vh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-uh(i,j,k,idblock)+vh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(uh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-uh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(vh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(vh(i,j,k,idblock)-wh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pxxh(i,j,k,idblock)=pxxh(i,j,k,idblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	pyyh(i,j,k,idblock)=pyyh(i,j,k,idblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pzzh(i,j,k,idblock)=pzzh(i,j,k,idblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	pxyh(i,j,k,idblock)=pxyh(i,j,k,idblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	pxzh(i,j,k,idblock)=pxzh(i,j,k,idblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	pyzh(i,j,k,idblock)=pyzh(i,j,k,idblock)-udotc
	
	    
    return
#endif	
  end subroutine streamcoll_bc_shared
  
    attributes(global) subroutine streamcoll_bc_shared_flop()
	
	implicit none  
	  
    integer :: i,j,k,gi,gj,gk,idblock
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
    integer :: idblock=1
		  
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
		  
	!if(isfluid(i,j,k).ne.1)return
    	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    uu=halfonecssq*(uh(i,j,k,idblock)*uh(i,j,k,idblock) + vh(i,j,k,idblock)*vh(i,j,k,idblock) + wh(i,j,k,idblock)*wh(i,j,k,idblock))
    
    !0
	f00(li,lj,lk)=p0*(rhoh(i,j,k,idblock)-uu) &
	 + oneminusomega*pi2cssq0*(-cssq*(pyyh(i,j,k,idblock)+pxxh(i,j,k,idblock)+pzzh(i,j,k,idblock)))
	
    
	!1 -1  0  0
	udotc=uh(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k,idblock)-cssq*(pyyh(i,j,k,idblock)+pzzh(i,j,k,idblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*pxxh(i,j,k,idblock)-cssq*(pyyh(i,j,k,idblock)+pzzh(i,j,k,idblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=vh(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k,idblock)-cssq*(pxxh(i,j,k,idblock)+pzzh(i,j,k,idblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*pyyh(i,j,k,idblock)-cssq*(pxxh(i,j,k,idblock)+pzzh(i,j,k,idblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=wh(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k,idblock)-cssq*(pxxh(i,j,k,idblock)+pyyh(i,j,k,idblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k,idblock)-cssq*(pxxh(i,j,k,idblock)+pyyh(i,j,k,idblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(uh(i,j,k,idblock)+vh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qyy*pyyh(i,j,k,idblock)-cssq*pzzh(i,j,k,idblock)+two*qxy_7_8*pxyh(i,j,k,idblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qyy*pyyh(i,j,k,idblock)-cssq*pzzh(i,j,k,idblock)+two*qxy_7_8*pxyh(i,j,k,idblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-uh(i,j,k,idblock)+vh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qyy*pyyh(i,j,k,idblock)-cssq*pzzh(i,j,k,idblock)+two*qxy_9_10*pxyh(i,j,k,idblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qyy*pyyh(i,j,k,idblock)-cssq*pzzh(i,j,k,idblock)+two*qxy_9_10*pxyh(i,j,k,idblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(uh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pyyh(i,j,k,idblock)+two*qxz_15_16*pxzh(i,j,k,idblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pyyh(i,j,k,idblock)+two*qxz_15_16*pxzh(i,j,k,idblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-uh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pyyh(i,j,k,idblock)+two*qxz_17_18*pxzh(i,j,k,idblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pyyh(i,j,k,idblock)+two*qxz_17_18*pxzh(i,j,k,idblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(vh(i,j,k,idblock)+wh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pxxh(i,j,k,idblock)+two*qyz_11_12*pyzh(i,j,k,idblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pxxh(i,j,k,idblock)+two*qyz_11_12*pyzh(i,j,k,idblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(vh(i,j,k,idblock)-wh(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pxxh(i,j,k,idblock)+two*qyz_13_14*pyzh(i,j,k,idblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(rhoh(i,j,k,idblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k,idblock)+qzz*pzzh(i,j,k,idblock)-cssq*pxxh(i,j,k,idblock)+two*qyz_13_14*pyzh(i,j,k,idblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    ! Halo Faces
    if(li==1)then
      
      uu=halfonecssq*(uh(i-1,j,k,idblock)*uh(i-1,j,k,idblock) + vh(i-1,j,k,idblock)*vh(i-1,j,k,idblock) + wh(i-1,j,k,idblock)*wh(i-1,j,k,idblock))
      
      !1 -1  0  0
	  udotc=uh(i-1,j,k,idblock)*onecssq
	  f01(li-1,lj,lk)=p1*(rhoh(i-1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxxh(i-1,j,k,idblock)-cssq*(pyyh(i-1,j,k,idblock)+pzzh(i-1,j,k,idblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(uh(i-1,j,k,idblock)+vh(i-1,j,k,idblock))*onecssq
	  f07(li-1,lj,lk)=p2*(rhoh(i-1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k,idblock)+qyy*pyyh(i-1,j,k,idblock)-cssq*pzzh(i-1,j,k,idblock)+two*qxy_7_8*pxyh(i-1,j,k,idblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-uh(i-1,j,k,idblock)+vh(i-1,j,k,idblock))*onecssq
	  f09(li-1,lj,lk)=p2*(rhoh(i-1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k,idblock)+qyy*pyyh(i-1,j,k,idblock)-cssq*pzzh(i-1,j,k,idblock)+two*qxy_9_10*pxyh(i-1,j,k,idblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(uh(i-1,j,k,idblock)+wh(i-1,j,k,idblock))*onecssq
	  f15(li-1,lj,lk)=p2*(rhoh(i-1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k,idblock)+qzz*pzzh(i-1,j,k,idblock)-cssq*pyyh(i-1,j,k,idblock)+two*qxz_15_16*pxzh(i-1,j,k,idblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-uh(i-1,j,k,idblock)+wh(i-1,j,k,idblock))*onecssq
	  f18(li-1,lj,lk)=p2*(rhoh(i-1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k,idblock)+qzz*pzzh(i-1,j,k,idblock)-cssq*pyyh(i-1,j,k,idblock)+two*qxz_17_18*pxzh(i-1,j,k,idblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      uu=halfonecssq*(uh(i+1,j,k,idblock)*uh(i+1,j,k,idblock) + vh(i+1,j,k,idblock)*vh(i+1,j,k,idblock) + wh(i+1,j,k,idblock)*wh(i+1,j,k,idblock))
      
      !2 +1  0  0
	  udotc=uh(i+1,j,k,idblock)*onecssq
	  f02(li+1,lj,lk)=p1*(rhoh(i+1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*pxxh(i+1,j,k,idblock)-cssq*(pyyh(i+1,j,k,idblock)+pzzh(i+1,j,k,idblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(uh(i+1,j,k,idblock)+vh(i+1,j,k,idblock))*onecssq
	  f08(li+1,lj,lk)=p2*(rhoh(i+1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k,idblock)+qyy*pyyh(i+1,j,k,idblock)-cssq*pzzh(i+1,j,k,idblock)+two*qxy_7_8*pxyh(i+1,j,k,idblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-uh(i+1,j,k,idblock)+vh(i+1,j,k,idblock))*onecssq
	  f10(li+1,lj,lk)=p2*(rhoh(i+1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k,idblock)+qyy*pyyh(i+1,j,k,idblock)-cssq*pzzh(i+1,j,k,idblock)+two*qxy_9_10*pxyh(i+1,j,k,idblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(uh(i+1,j,k,idblock)+wh(i+1,j,k,idblock))*onecssq
	  f16(li+1,lj,lk)=p2*(rhoh(i+1,j,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k,idblock)+qzz*pzzh(i+1,j,k,idblock)-cssq*pyyh(i+1,j,k,idblock)+two*qxz_15_16*pxzh(i+1,j,k,idblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-uh(i+1,j,k,idblock)+wh(i+1,j,k,idblock))*onecssq
	  f17(li+1,lj,lk)=p2*(rhoh(i+1,j,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k,idblock)+qzz*pzzh(i+1,j,k,idblock)-cssq*pyyh(i+1,j,k,idblock)+two*qxz_17_18*pxzh(i+1,j,k,idblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      uu=halfonecssq*(uh(i,j-1,k,idblock)*uh(i,j-1,k,idblock) + vh(i,j-1,k,idblock)*vh(i,j-1,k,idblock) + wh(i,j-1,k,idblock)*wh(i,j-1,k,idblock))
      
            !3 0 -1  0
	  udotc=vh(i,j-1,k,idblock)*onecssq
	  f03(li,lj-1,lk)=p1*(rhoh(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyyh(i,j-1,k,idblock)-cssq*(pxxh(i,j-1,k,idblock)+pzzh(i,j-1,k,idblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(uh(i,j-1,k,idblock)+vh(i,j-1,k,idblock))*onecssq
	  f07(li,lj-1,lk)=p2*(rhoh(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j-1,k,idblock)+qyy*pyyh(i,j-1,k,idblock)-cssq*pzzh(i,j-1,k,idblock)+two*qxy_7_8*pxyh(i,j-1,k,idblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-uh(i,j-1,k,idblock)+vh(i,j-1,k,idblock))*onecssq
	  f10(li,lj-1,lk)=p2*(rhoh(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j-1,k,idblock)+qyy*pyyh(i,j-1,k,idblock)-cssq*pzzh(i,j-1,k,idblock)+two*qxy_9_10*pxyh(i,j-1,k,idblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(vh(i,j-1,k,idblock)+wh(i,j-1,k,idblock))*onecssq
	  f11(li,lj-1,lk)=p2*(rhoh(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k,idblock)+qzz*pzzh(i,j-1,k,idblock)-cssq*pxxh(i,j-1,k,idblock)+two*qyz_11_12*pyzh(i,j-1,k,idblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(vh(i,j-1,k,idblock)-wh(i,j-1,k,idblock))*onecssq
	  f13(li,lj-1,lk)=p2*(rhoh(i,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k,idblock)+qzz*pzzh(i,j-1,k,idblock)-cssq*pxxh(i,j-1,k,idblock)+two*qyz_13_14*pyzh(i,j-1,k,idblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(uh(i,j+1,k,idblock)*uh(i,j+1,k,idblock) + vh(i,j+1,k,idblock)*vh(i,j+1,k,idblock) + wh(i,j+1,k,idblock)*wh(i,j+1,k,idblock))
      
      !4  0 +1  0
	  udotc=vh(i,j+1,k,idblock)*onecssq
	  f04(li,lj+1,lk)=p1*(rhoh(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*pyyh(i,j+1,k,idblock)-cssq*(pxxh(i,j+1,k,idblock)+pzzh(i,j+1,k,idblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(uh(i,j+1,k,idblock)+vh(i,j+1,k,idblock))*onecssq
	  f08(li,lj+1,lk)=p2*(rhoh(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j+1,k,idblock)+qyy*pyyh(i,j+1,k,idblock)-cssq*pzzh(i,j+1,k,idblock)+two*qxy_7_8*pxyh(i,j+1,k,idblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-uh(i,j+1,k,idblock)+vh(i,j+1,k,idblock))*onecssq
	  f09(li,lj+1,lk)=p2*(rhoh(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j+1,k,idblock)+qyy*pyyh(i,j+1,k,idblock)-cssq*pzzh(i,j+1,k,idblock)+two*qxy_9_10*pxyh(i,j+1,k,idblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(vh(i,j+1,k,idblock)+wh(i,j+1,k,idblock))*onecssq
	  f12(li,lj+1,lk)=p2*(rhoh(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k,idblock)+qzz*pzzh(i,j+1,k,idblock)-cssq*pxxh(i,j+1,k,idblock)+two*qyz_11_12*pyzh(i,j+1,k,idblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(vh(i,j+1,k,idblock)-wh(i,j+1,k,idblock))*onecssq
	  f14(li,lj+1,lk)=p2*(rhoh(i,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k,idblock)+qzz*pzzh(i,j+1,k,idblock)-cssq*pxxh(i,j+1,k,idblock)+two*qyz_13_14*pyzh(i,j+1,k,idblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      uu=halfonecssq*(uh(i,j,k-1,idblock)*uh(i,j,k-1,idblock) + vh(i,j,k-1,idblock)*vh(i,j,k-1,idblock) + wh(i,j,k-1,idblock)*wh(i,j,k-1,idblock))
      
      !5  0  0 -1
	  udotc=wh(i,j,k-1,idblock)*onecssq
	  f05(li,lj,lk-1)=p1*(rhoh(i,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k-1,idblock)-cssq*(pxxh(i,j,k-1,idblock)+pyyh(i,j,k-1,idblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(uh(i,j,k-1,idblock)+wh(i,j,k-1,idblock))*onecssq
	  f15(li,lj,lk-1)=p2*(rhoh(i,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k-1,idblock)+qzz*pzzh(i,j,k-1,idblock)-cssq*pyyh(i,j,k-1,idblock)+two*qxz_15_16*pxzh(i,j,k-1,idblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-uh(i,j,k-1,idblock)+wh(i,j,k-1,idblock))*onecssq
	  f17(li,lj,lk-1)=p2*(rhoh(i,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k-1,idblock)+qzz*pzzh(i,j,k-1,idblock)-cssq*pyyh(i,j,k-1,idblock)+two*qxz_17_18*pxzh(i,j,k-1,idblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(vh(i,j,k-1,idblock)+wh(i,j,k-1,idblock))*onecssq
	  f11(li,lj,lk-1)=p2*(rhoh(i,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k-1,idblock)+qzz*pzzh(i,j,k-1,idblock)-cssq*pxxh(i,j,k-1,idblock)+two*qyz_11_12*pyzh(i,j,k-1,idblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(vh(i,j,k-1,idblock)-wh(i,j,k-1,idblock))*onecssq
	  f14(li,lj,lk-1)=p2*(rhoh(i,j,k-1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k-1,idblock)+qzz*pzzh(i,j,k-1,idblock)-cssq*pxxh(i,j,k-1,idblock)+two*qyz_13_14*pyzh(i,j,k-1,idblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i,j,k+1,idblock)*uh(i,j,k+1,idblock) + vh(i,j,k+1,idblock)*vh(i,j,k+1,idblock) + wh(i,j,k+1,idblock)*wh(i,j,k+1,idblock))
      
      !6  0  0  +1
	  udotc=wh(i,j,k+1,idblock)*onecssq
	  f06(li,lj,lk+1)=p1*(rhoh(i,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*pzzh(i,j,k+1,idblock)-cssq*(pxxh(i,j,k+1,idblock)+pyyh(i,j,k+1,idblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(uh(i,j,k+1,idblock)+wh(i,j,k+1,idblock))*onecssq
	  f16(li,lj,lk+1)=p2*(rhoh(i,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k+1,idblock)+qzz*pzzh(i,j,k+1,idblock)-cssq*pyyh(i,j,k+1,idblock)+two*qxz_15_16*pxzh(i,j,k+1,idblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-uh(i,j,k+1,idblock)+wh(i,j,k+1,idblock))*onecssq
	  f18(li,lj,lk+1)=p2*(rhoh(i,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i,j,k+1,idblock)+qzz*pzzh(i,j,k+1,idblock)-cssq*pyyh(i,j,k+1,idblock)+two*qxz_17_18*pxzh(i,j,k+1,idblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(vh(i,j,k+1,idblock)+wh(i,j,k+1,idblock))*onecssq
	  f12(li,lj,lk+1)=p2*(rhoh(i,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k+1,idblock)+qzz*pzzh(i,j,k+1,idblock)-cssq*pxxh(i,j,k+1,idblock)+two*qyz_11_12*pyzh(i,j,k+1,idblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(vh(i,j,k+1,idblock)-wh(i,j,k+1,idblock))*onecssq
	  f13(li,lj,lk+1)=p2*(rhoh(i,j,k+1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j,k+1,idblock)+qzz*pzzh(i,j,k+1,idblock)-cssq*pxxh(i,j,k+1,idblock)+two*qyz_13_14*pyzh(i,j,k+1,idblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      uu=halfonecssq*(uh(i-1,j-1,k,idblock)*uh(i-1,j-1,k,idblock) + vh(i-1,j-1,k,idblock)*vh(i-1,j-1,k,idblock) + wh(i-1,j-1,k,idblock)*wh(i-1,j-1,k,idblock))
      
      !7 -1 -1  0
	  udotc=(uh(i-1,j-1,k,idblock)+vh(i-1,j-1,k,idblock))*onecssq
	  f07(li-1,lj-1,lk)=p2*(rhoh(i-1,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j-1,k,idblock)+qyy*pyyh(i-1,j-1,k,idblock)-cssq*pzzh(i-1,j-1,k,idblock)+two*qxy_7_8*pxyh(i-1,j-1,k,idblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(uh(i-1,j+1,k,idblock)*uh(i-1,j+1,k,idblock) + vh(i-1,j+1,k,idblock)*vh(i-1,j+1,k,idblock) + wh(i-1,j+1,k,idblock)*wh(i-1,j+1,k,idblock))
      
      !9  -1 +1 0
      udotc=(-uh(i-1,j+1,k,idblock)+vh(i-1,j+1,k,idblock))*onecssq
	  f09(li-1,lj+1,lk)=p2*(rhoh(i-1,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j+1,k,idblock)+qyy*pyyh(i-1,j+1,k,idblock)-cssq*pzzh(i-1,j+1,k,idblock)+two*qxy_9_10*pxyh(i-1,j+1,k,idblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      uu=halfonecssq*(uh(i-1,j,k-1,idblock)*uh(i-1,j,k-1,idblock) + vh(i-1,j,k-1,idblock)*vh(i-1,j,k-1,idblock) + wh(i-1,j,k-1,idblock)*wh(i-1,j,k-1,idblock))
      
      !15  -1  0 -1
	  udotc=(uh(i-1,j,k-1,idblock)+wh(i-1,j,k-1,idblock))*onecssq
	  f15(li-1,lj,lk-1)=p2*(rhoh(i-1,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k-1,idblock)+qzz*pzzh(i-1,j,k-1,idblock)-cssq*pyyh(i-1,j,k-1,idblock)+two*qxz_15_16*pxzh(i-1,j,k-1,idblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i-1,j,k+1,idblock)*uh(i-1,j,k+1,idblock) + vh(i-1,j,k+1,idblock)*vh(i-1,j,k+1,idblock) + wh(i-1,j,k+1,idblock)*wh(i-1,j,k+1,idblock))
      
      !18   -1   0  +1
	  udotc=(-uh(i-1,j,k+1,idblock)+wh(i-1,j,k+1,idblock))*onecssq
	  f18(li-1,lj,lk+1)=p2*(rhoh(i-1,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i-1,j,k+1,idblock)+qzz*pzzh(i-1,j,k+1,idblock)-cssq*pyyh(i-1,j,k+1,idblock)+two*qxz_17_18*pxzh(i-1,j,k+1,idblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      uu=halfonecssq*(uh(i+1,j-1,k,idblock)*uh(i+1,j-1,k,idblock) + vh(i+1,j-1,k,idblock)*vh(i+1,j-1,k,idblock) + wh(i+1,j-1,k,idblock)*wh(i+1,j-1,k,idblock))
      
      !10   +1 -1  0
	  udotc=(-uh(i+1,j-1,k,idblock)+vh(i+1,j-1,k,idblock))*onecssq
	  f10(li+1,lj-1,lk)=p2*(rhoh(i+1,j-1,k,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j-1,k,idblock)+qyy*pyyh(i+1,j-1,k,idblock)-cssq*pzzh(i+1,j-1,k,idblock)+two*qxy_9_10*pxyh(i+1,j-1,k,idblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      uu=halfonecssq*(uh(i+1,j+1,k,idblock)*uh(i+1,j+1,k,idblock) + vh(i+1,j+1,k,idblock)*vh(i+1,j+1,k,idblock) + wh(i+1,j+1,k,idblock)*wh(i+1,j+1,k,idblock))
      
      !8 +1 +1  0
	  udotc=(uh(i+1,j+1,k,idblock)+vh(i+1,j+1,k,idblock))*onecssq
	  f08(li+1,lj+1,lk)=p2*(rhoh(i+1,j+1,k,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j+1,k,idblock)+qyy*pyyh(i+1,j+1,k,idblock)-cssq*pzzh(i+1,j+1,k,idblock)+two*qxy_7_8*pxyh(i+1,j+1,k,idblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      uu=halfonecssq*(uh(i+1,j,k-1,idblock)*uh(i+1,j,k-1,idblock) + vh(i+1,j,k-1,idblock)*vh(i+1,j,k-1,idblock) + wh(i+1,j,k-1,idblock)*wh(i+1,j,k-1,idblock))
      
      !17  +1  0 -1
	  udotc=(-uh(i+1,j,k-1,idblock)+wh(i+1,j,k-1,idblock))*onecssq
	  f17(li+1,lj,lk-1)=p2*(rhoh(i+1,j,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k-1,idblock)+qzz*pzzh(i+1,j,k-1,idblock)-cssq*pyyh(i+1,j,k-1,idblock)+two*qxz_17_18*pxzh(i+1,j,k-1,idblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i+1,j,k+1,idblock)*uh(i+1,j,k+1,idblock) + vh(i+1,j,k+1,idblock)*vh(i+1,j,k+1,idblock) + wh(i+1,j,k+1,idblock)*wh(i+1,j,k+1,idblock))
      
      !16  +1  0 +1
	  udotc=(uh(i+1,j,k+1,idblock)+wh(i+1,j,k+1,idblock))*onecssq
	  f16(li+1,lj,lk+1)=p2*(rhoh(i+1,j,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*pxxh(i+1,j,k+1,idblock)+qzz*pzzh(i+1,j,k+1,idblock)-cssq*pyyh(i+1,j,k+1,idblock)+two*qxz_15_16*pxzh(i+1,j,k+1,idblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      uu=halfonecssq*(uh(i,j-1,k-1,idblock)*uh(i,j-1,k-1,idblock) + vh(i,j-1,k-1,idblock)*vh(i,j-1,k-1,idblock) + wh(i,j-1,k-1,idblock)*wh(i,j-1,k-1,idblock))
      
      !11  0  -1  -1
	  udotc=(vh(i,j-1,k-1,idblock)+wh(i,j-1,k-1,idblock))*onecssq
	  f11(li,lj-1,lk-1)=p2*(rhoh(i,j-1,k-1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k-1,idblock)+qzz*pzzh(i,j-1,k-1,idblock)-cssq*pxxh(i,j-1,k-1,idblock)+two*qyz_11_12*pyzh(i,j-1,k-1,idblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i,j-1,k+1,idblock)*uh(i,j-1,k+1,idblock) + vh(i,j-1,k+1,idblock)*vh(i,j-1,k+1,idblock) + wh(i,j-1,k+1,idblock)*wh(i,j-1,k+1,idblock))
      
      !13  0  -1   +1
	  udotc=(vh(i,j-1,k+1,idblock)-wh(i,j-1,k+1,idblock))*onecssq
	  f13(li,lj-1,lk+1)=p2*(rhoh(i,j-1,k+1,idblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j-1,k+1,idblock)+qzz*pzzh(i,j-1,k+1,idblock)-cssq*pxxh(i,j-1,k+1,idblock)+two*qyz_13_14*pyzh(i,j-1,k+1,idblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      uu=halfonecssq*(uh(i,j+1,k-1,idblock)*uh(i,j+1,k-1,idblock) + vh(i,j+1,k-1,idblock)*vh(i,j+1,k-1,idblock) + wh(i,j+1,k-1,idblock)*wh(i,j+1,k-1,idblock))
      
      !14  0  +1  -1
	  udotc=(vh(i,j+1,k-1,idblock)-wh(i,j+1,k-1,idblock))*onecssq
	  f14(li,lj+1,lk-1)=p2*(rhoh(i,j+1,k-1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k-1,idblock)+qzz*pzzh(i,j+1,k-1,idblock)-cssq*pxxh(i,j+1,k-1,idblock)+two*qyz_13_14*pyzh(i,j+1,k-1,idblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      uu=halfonecssq*(uh(i,j+1,k+1,idblock)*uh(i,j+1,k+1,idblock) + vh(i,j+1,k+1,idblock)*vh(i,j+1,k+1,idblock) + wh(i,j+1,k+1,idblock)*wh(i,j+1,k+1,idblock))
      
      !12  0  +1  +1
	  udotc=(vh(i,j+1,k+1,idblock)+wh(i,j+1,k+1,idblock))*onecssq
	  f12(li,lj+1,lk+1)=p2*(rhoh(i,j+1,k+1,idblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*pyyh(i,j+1,k+1,idblock)+qzz*pzzh(i,j+1,k+1,idblock)-cssq*pxxh(i,j+1,k+1,idblock)+two*qyz_11_12*pyzh(i,j+1,k+1,idblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
    
    if(isfluid(gi,gj,gk).eq.0)then
#ifdef IFBC
          
          if(isfluid(i+1,j,k-1).eq.-1)f18(li,lj,lk)=f17(li+1,lj,lk-1) !gpc 
          if(isfluid(i-1,j,k+1).eq.-1)f17(li,lj,lk)=f18(li-1,lj,lk+1) !hpc

          if(isfluid(i-1,j,k-1).eq.-1)f16(li,lj,lk)=f15(li-1,lj,lk-1) !gpc 
          if(isfluid(i+1,j,k+1).eq.-1)f15(li,lj,lk)=f16(li+1,lj,lk+1) !hpc

          if(isfluid(i,j-1,k+1).eq.-1)f14(li,lj,lk)=f13(li,lj-1,lk+1)!gpc 
          if(isfluid(i,j+1,k-1).eq.-1)f13(li,lj,lk)=f14(li,lj+1,lk-1)!hpc
        
          if(isfluid(i,j-1,k-1).eq.-1)f12(li,lj,lk)=f11(li,lj-1,lk-1)!gpc
          if(isfluid(i,j+1,k+1).eq.-1)f11(li,lj,lk)=f12(li,lj+1,lk+1)!hpc

          if(isfluid(i-1,j+1,k).eq.-1)f10(li,lj,lk)=f09(li-1,lj+1,lk)!gpc  
          if(isfluid(i+1,j-1,k).eq.-1)f09(li,lj,lk)=f10(li+1,lj-1,lk)!hpc

          if(isfluid(i-1,j-1,k).eq.-1)f08(li,lj,lk)=f07(li-1,lj-1,lk)!gpc 
          if(isfluid(i+1,j+1,k).eq.-1)f07(li,lj,lk)=f08(li+1,lj+1,lk)!hpc

          if(isfluid(i,j,k-1).eq.-1)f06(li,lj,lk)=f05(li,lj,lk-1)!gpc
          if(isfluid(i,j,k+1).eq.-1)f05(li,lj,lk)=f06(li,lj,lk+1)!hpc 

          if(isfluid(i,j-1,k).eq.-1)f04(li,lj,lk)=f03(li,lj-1,lk)!gpc 
          if(isfluid(i,j+1,k).eq.-1)f03(li,lj,lk)=f04(li,lj+1,lk)!hpc 

          if(isfluid(i-1,j,k).eq.-1)f02(li,lj,lk)=f01(li-1,lj,lk)!gpc
          if(isfluid(i+1,j,k).eq.-1)f01(li,lj,lk)=f02(li+1,lj,lk)!hpc 
          
#else
          
          f18(li,lj,lk)=f17(li+1,lj,lk-1) !gpc 
          f17(li,lj,lk)=f18(li-1,lj,lk+1) !hpc

          f16(li,lj,lk)=f15(li-1,lj,lk-1) !gpc 
          f15(li,lj,lk)=f16(li+1,lj,lk+1) !hpc

          f14(li,lj,lk)=f13(li,lj-1,lk+1)!gpc 
          f13(li,lj,lk)=f14(li,lj+1,lk-1)!hpc
        
          f12(li,lj,lk)=f11(li,lj-1,lk-1)!gpc 
          f11(li,lj,lk)=f12(li,lj+1,lk+1)!hpc

          f10(li,lj,lk)=f09(li-1,lj+1,lk)!gpc 
          f09(li,lj,lk)=f10(li+1,lj-1,lk)!hpc

          f08(li,lj,lk)=f07(li-1,lj-1,lk)!gpc 
          f07(li,lj,lk)=f08(li+1,lj+1,lk)!hpc

          f06(li,lj,lk)=f05(li,lj,lk-1)!gpc 
          f05(li,lj,lk)=f06(li,lj,lk+1)!hpc 

          f04(li,lj,lk)=f03(li,lj-1,lk)!gpc 
          f03(li,lj,lk)=f04(li,lj+1,lk)!hpc 

          f02(li,lj,lk)=f01(li-1,lj,lk)!gpc 
          f01(li,lj,lk)=f02(li+1,lj,lk)!hpc          

#endif
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
	rho(i,j,k,idblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	u(i,j,k,idblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	v(i,j,k,idblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	w(i,j,k,idblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pxx(i,j,k,idblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	pyy(i,j,k,idblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	pzz(i,j,k,idblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	pxy(i,j,k,idblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	pxz(i,j,k,idblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	pyz(i,j,k,idblock)=udotc
    
#ifdef PRESSCORR
	call syncthreads
	
	uu=halfonecssq*(u(i,j,k,idblock)*u(i,j,k,idblock) + v(i,j,k,idblock)*v(i,j,k,idblock) + w(i,j,k,idblock)*w(i,j,k,idblock))
    
	!1 -1  0  0
	udotc=u(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=v(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=w(i,j,k,idblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(rho(i,j,k,idblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(u(i,j,k,idblock)+v(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-u(i,j,k,idblock)+v(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(u(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-u(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(v(i,j,k,idblock)+w(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(v(i,j,k,idblock)-w(i,j,k,idblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(rho(i,j,k,idblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pxx(i,j,k,idblock)=pxx(i,j,k,idblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	pyy(i,j,k,idblock)=pyy(i,j,k,idblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	pzz(i,j,k,idblock)=pzz(i,j,k,idblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	pxy(i,j,k,idblock)=pxy(i,j,k,idblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	pxz(i,j,k,idblock)=pxz(i,j,k,idblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	pyz(i,j,k,idblock)=pyz(i,j,k,idblock)-udotc
	
	    
    return
#endif	
  end subroutine streamcoll_bc_shared_flop
  
 end module streamcoll_bc_kernels
