#include "defines.h"
 module correct_press_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine correct_pressure
    
    implicit none
    
    integer :: i,j,k,gi,gj,gk,idblock
     
	real(kind=db) :: uu,udotc,temp,feq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
		  
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	
	if(abs(isfluid(gi,gj,gk)).ne.1)return
                        
	uu=half*(uh(i,j,k,idblock)*uh(i,j,k,idblock) + vh(i,j,k,idblock)*vh(i,j,k,idblock) + wh(i,j,k,idblock)*wh(i,j,k,idblock))/cssq
	!1-2
	udotc=uh(i,j,k,idblock)/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pxx=feq
	feq=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq

	!3-4
	udotc=vh(i,j,k,idblock)/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pyy=feq
	feq=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	!5-6
	udotc=wh(i,j,k,idblock)/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pzz=feq
	feq=p1*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pzz=temp_pzz+feq
	!7-8
	udotc=(uh(i,j,k,idblock)+vh(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=feq
	feq=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy+feq
	!10-9
	udotc=(-uh(i,j,k,idblock)+vh(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	feq=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	!11-12
	udotc=(vh(i,j,k,idblock)+wh(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=feq
	feq=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz+feq
	!13-14
	udotc=(vh(i,j,k,idblock)-wh(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	feq=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	!15-16
	udotc=(uh(i,j,k,idblock)+wh(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=feq
	feq=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz+feq
	!17-18
	udotc=(-uh(i,j,k,idblock)+wh(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rhoh(i,j,k,idblock)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	feq=p2*(rhoh(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	!ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
	
	pxxh(i,j,k,idblock)=pxxh(i,j,k,idblock)-temp_pxx
	pyyh(i,j,k,idblock)=pyyh(i,j,k,idblock)-temp_pyy
	pzzh(i,j,k,idblock)=pzzh(i,j,k,idblock)-temp_pzz
	pxyh(i,j,k,idblock)=pxyh(i,j,k,idblock)-temp_pxy
	pxzh(i,j,k,idblock)=pxzh(i,j,k,idblock)-temp_pxz
	pyzh(i,j,k,idblock)=pyzh(i,j,k,idblock)-temp_pyz
    
    return
    
  end subroutine correct_pressure
  
  attributes(global) subroutine correct_pressure_flop
    
    implicit none
    
    integer :: i,j,k,gi,gj,gk,idblock
    
	real(kind=db) :: uu,udotc,temp,feq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
		  
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	
	if(abs(isfluid(gi,gj,gk)).ne.1)return
	
    uu=half*(u(i,j,k,idblock)*u(i,j,k,idblock) + v(i,j,k,idblock)*v(i,j,k,idblock) + w(i,j,k,idblock)*w(i,j,k,idblock))/cssq
	!1-2
	udotc=u(i,j,k,idblock)/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pxx=feq
	feq=p1*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq

	!3-4
	udotc=v(i,j,k,idblock)/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pyy=feq
	feq=p1*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	!5-6
	udotc=w(i,j,k,idblock)/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pzz=feq
	feq=p1*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pzz=temp_pzz+feq
	!7-8
	udotc=(u(i,j,k,idblock)+v(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=feq
	feq=p2*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy+feq
	!10-9
	udotc=(-u(i,j,k,idblock)+v(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	feq=p2*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	!11-12
	udotc=(v(i,j,k,idblock)+w(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=feq
	feq=p2*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz+feq
	!13-14
	udotc=(v(i,j,k,idblock)-w(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	feq=p2*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	!15-16
	udotc=(u(i,j,k,idblock)+w(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=feq
	feq=p2*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz+feq
	!17-18
	udotc=(-u(i,j,k,idblock)+w(i,j,k,idblock))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(rho(i,j,k,idblock)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	feq=p2*(rho(i,j,k,idblock)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	!ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
	
	pxx(i,j,k,idblock)=pxx(i,j,k,idblock)-temp_pxx
	pyy(i,j,k,idblock)=pyy(i,j,k,idblock)-temp_pyy
	pzz(i,j,k,idblock)=pzz(i,j,k,idblock)-temp_pzz
	pxy(i,j,k,idblock)=pxy(i,j,k,idblock)-temp_pxy
	pxz(i,j,k,idblock)=pxz(i,j,k,idblock)-temp_pxz
	pyz(i,j,k,idblock)=pyz(i,j,k,idblock)-temp_pyz
	
	return
    
  end subroutine correct_pressure_flop
  
 end module correct_press_kernels
