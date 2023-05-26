#include "defines.h"
 module correct_press_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine correct_pressure
    
    implicit none
    
    integer :: i,j,k
    
	real(kind=db) :: uu,udotc,temp,feq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
		  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	if(abs(isfluid(i,j,k)).ne.1)return
                        
	uu=0.5_db*(uh(i,j,k)*uh(i,j,k) + vh(i,j,k)*vh(i,j,k) + wh(i,j,k)*wh(i,j,k))/cssq
	!1-2
	udotc=uh(i,j,k)/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p1*(rhoh(i,j,k)+(temp + udotc))
	temp_pxx=feq
	feq=p1*(rhoh(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq

	!3-4
	udotc=vh(i,j,k)/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p1*(rhoh(i,j,k)+(temp + udotc))
	temp_pyy=feq
	feq=p1*(rhoh(i,j,k)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	!5-6
	udotc=wh(i,j,k)/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p1*(rhoh(i,j,k)+(temp + udotc))
	temp_pzz=feq
	feq=p1*(rhoh(i,j,k)+(temp - udotc))
	temp_pzz=temp_pzz+feq
	!7-8
	udotc=(uh(i,j,k)+vh(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rhoh(i,j,k)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=feq
	feq=p2*(rhoh(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy+feq
	!10-9
	udotc=(-uh(i,j,k)+vh(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rhoh(i,j,k)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	feq=p2*(rhoh(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	!11-12
	udotc=(vh(i,j,k)+wh(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rhoh(i,j,k)+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=feq
	feq=p2*(rhoh(i,j,k)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz+feq
	!13-14
	udotc=(vh(i,j,k)-wh(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rhoh(i,j,k)+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	feq=p2*(rhoh(i,j,k)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	!15-16
	udotc=(uh(i,j,k)+wh(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rhoh(i,j,k)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=feq
	feq=p2*(rhoh(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz+feq
	!17-18
	udotc=(-uh(i,j,k)+wh(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rhoh(i,j,k)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	feq=p2*(rhoh(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	!ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
	
	pxxh(i,j,k)=pxxh(i,j,k)-temp_pxx
	pyyh(i,j,k)=pyyh(i,j,k)-temp_pyy
	pzzh(i,j,k)=pzzh(i,j,k)-temp_pzz
	pxyh(i,j,k)=pxyh(i,j,k)-temp_pxy
	pxzh(i,j,k)=pxzh(i,j,k)-temp_pxz
	pyzh(i,j,k)=pyzh(i,j,k)-temp_pyz
    
    return
    
  end subroutine correct_pressure
  
  attributes(global) subroutine correct_pressure_flop
    
    implicit none
    
    integer :: i,j,k
    
	real(kind=db) :: uu,udotc,temp,feq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
		  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	if(abs(isfluid(i,j,k)).ne.1)return
	
    uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
	!1-2
	udotc=u(i,j,k)/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p1*(rho(i,j,k)+(temp + udotc))
	temp_pxx=feq
	feq=p1*(rho(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq

	!3-4
	udotc=v(i,j,k)/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p1*(rho(i,j,k)+(temp + udotc))
	temp_pyy=feq
	feq=p1*(rho(i,j,k)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	!5-6
	udotc=w(i,j,k)/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p1*(rho(i,j,k)+(temp + udotc))
	temp_pzz=feq
	feq=p1*(rho(i,j,k)+(temp - udotc))
	temp_pzz=temp_pzz+feq
	!7-8
	udotc=(u(i,j,k)+v(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rho(i,j,k)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=feq
	feq=p2*(rho(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy+feq
	!10-9
	udotc=(-u(i,j,k)+v(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rho(i,j,k)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	feq=p2*(rho(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	!11-12
	udotc=(v(i,j,k)+w(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rho(i,j,k)+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=feq
	feq=p2*(rho(i,j,k)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz+feq
	!13-14
	udotc=(v(i,j,k)-w(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rho(i,j,k)+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	feq=p2*(rho(i,j,k)+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	!15-16
	udotc=(u(i,j,k)+w(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rho(i,j,k)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=feq
	feq=p2*(rho(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz+feq
	!17-18
	udotc=(-u(i,j,k)+w(i,j,k))/cssq
	temp = -uu + 0.5_db*udotc*udotc
	feq=p2*(rho(i,j,k)+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	feq=p2*(rho(i,j,k)+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	!ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
	
	pxx(i,j,k)=pxx(i,j,k)-temp_pxx
	pyy(i,j,k)=pyy(i,j,k)-temp_pyy
	pzz(i,j,k)=pzz(i,j,k)-temp_pzz
	pxy(i,j,k)=pxy(i,j,k)-temp_pxy
	pxz(i,j,k)=pxz(i,j,k)-temp_pxz
	pyz(i,j,k)=pyz(i,j,k)-temp_pyz
	
	return
    
  end subroutine correct_pressure_flop
  
 end module correct_press_kernels
