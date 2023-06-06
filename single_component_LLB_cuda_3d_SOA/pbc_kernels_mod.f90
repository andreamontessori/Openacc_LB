#include "defines.h"
 module pbc_kernels
  
  use cudavars
  
  implicit none
  
 contains
 
 attributes(global) subroutine pbc_edge_x()
	  
	integer :: j,k
	
	  
	j = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	if (j>ny_d .or. k>nz_d)return
	  
	rho(0,j,k)=rho(nx_d,j,k)
	u(0,j,k)=u(nx_d,j,k)
	v(0,j,k)=v(nx_d,j,k)
	w(0,j,k)=w(nx_d,j,k)
	pxx(0,j,k)=pxx(nx_d,j,k)
	pyy(0,j,k)=pyy(nx_d,j,k)
	pzz(0,j,k)=pzz(nx_d,j,k)
	pxy(0,j,k)=pxy(nx_d,j,k)
	pxz(0,j,k)=pxz(nx_d,j,k)
	pyz(0,j,k)=pyz(nx_d,j,k)
	
	rho(nx_d+1,j,k)=rho(1,j,k)
	u(nx_d+1,j,k)=u(1,j,k)
	v(nx_d+1,j,k)=v(1,j,k)
	w(nx_d+1,j,k)=w(1,j,k)
	pxx(nx_d+1,j,k)=pxx(1,j,k)
	pyy(nx_d+1,j,k)=pyy(1,j,k)
	pzz(nx_d+1,j,k)=pzz(1,j,k)
	pxy(nx_d+1,j,k)=pxy(1,j,k)
	pxz(nx_d+1,j,k)=pxz(1,j,k)
	pyz(nx_d+1,j,k)=pyz(1,j,k)
	  
    return
		  
  end subroutine pbc_edge_x
  
  attributes(global) subroutine pbc_edge_x_flop()
	  
	integer :: j,k
	
	  
	j = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	if (j>ny_d .or. k>nz_d)return
	  
	rhoh(0,j,k)=rhoh(nx_d,j,k)
	uh(0,j,k)=uh(nx_d,j,k)
	vh(0,j,k)=vh(nx_d,j,k)
	wh(0,j,k)=wh(nx_d,j,k)
	pxxh(0,j,k)=pxxh(nx_d,j,k)
	pyyh(0,j,k)=pyyh(nx_d,j,k)
	pzzh(0,j,k)=pzzh(nx_d,j,k)
	pxyh(0,j,k)=pxyh(nx_d,j,k)
	pxzh(0,j,k)=pxzh(nx_d,j,k)
	pyzh(0,j,k)=pyzh(nx_d,j,k)
	
	rhoh(nx_d+1,j,k)=rhoh(1,j,k)
	uh(nx_d+1,j,k)=uh(1,j,k)
	vh(nx_d+1,j,k)=vh(1,j,k)
	wh(nx_d+1,j,k)=wh(1,j,k)
	pxxh(nx_d+1,j,k)=pxxh(1,j,k)
	pyyh(nx_d+1,j,k)=pyyh(1,j,k)
	pzzh(nx_d+1,j,k)=pzzh(1,j,k)
	pxyh(nx_d+1,j,k)=pxyh(1,j,k)
	pxzh(nx_d+1,j,k)=pxzh(1,j,k)
	pyzh(nx_d+1,j,k)=pyzh(1,j,k)
	  
    return
		  
  end subroutine pbc_edge_x_flop

  attributes(global) subroutine pbc_edge_y()
	
	implicit none
	integer :: i,k
	
	  
	i = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

	if (i>nx_d .or. k>nz_d)return
	
	rho(i,0,k)=rho(i,ny_d,k)
	u(i,0,k)=u(i,ny_d,k)
	v(i,0,k)=v(i,ny_d,k)
	w(i,0,k)=w(i,ny_d,k)
	pxx(i,0,k)=pxx(i,ny_d,k)
	pyy(i,0,k)=pyy(i,ny_d,k)
	pzz(i,0,k)=pzz(i,ny_d,k)
	pxy(i,0,k)=pxy(i,ny_d,k)
	pxz(i,0,k)=pxz(i,ny_d,k)
	pyz(i,0,k)=pyz(i,ny_d,k)
	
	rho(i,ny_d+1,k)=rho(i,1,k)
	u(i,ny_d+1,k)=u(i,1,k)
	v(i,ny_d+1,k)=v(i,1,k)
	w(i,ny_d+1,k)=w(i,1,k)
	pxx(i,ny_d+1,k)=pxx(i,1,k)
	pyy(i,ny_d+1,k)=pyy(i,1,k)
	pzz(i,ny_d+1,k)=pzz(i,1,k)
	pxy(i,ny_d+1,k)=pxy(i,1,k)
	pxz(i,ny_d+1,k)=pxz(i,1,k)
	pyz(i,ny_d+1,k)=pyz(i,1,k)
	  
    return

  end subroutine pbc_edge_y
  
  attributes(global) subroutine pbc_edge_y_flop()
	
	implicit none
	
	integer :: i,k
	
	  
	i = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

	if (i>nx_d .or. k>nz_d)return
	
	rhoh(i,0,k)=rhoh(i,ny_d,k)
	uh(i,0,k)=uh(i,ny_d,k)
	vh(i,0,k)=vh(i,ny_d,k)
	wh(i,0,k)=wh(i,ny_d,k)
	pxxh(i,0,k)=pxxh(i,ny_d,k)
	pyyh(i,0,k)=pyyh(i,ny_d,k)
	pzzh(i,0,k)=pzzh(i,ny_d,k)
	pxyh(i,0,k)=pxyh(i,ny_d,k)
	pxzh(i,0,k)=pxzh(i,ny_d,k)
	pyzh(i,0,k)=pyzh(i,ny_d,k)
	
	rhoh(i,ny_d+1,k)=rhoh(i,1,k)
	uh(i,ny_d+1,k)=uh(i,1,k)
	vh(i,ny_d+1,k)=vh(i,1,k)
	wh(i,ny_d+1,k)=wh(i,1,k)
	pxxh(i,ny_d+1,k)=pxxh(i,1,k)
	pyyh(i,ny_d+1,k)=pyyh(i,1,k)
	pzzh(i,ny_d+1,k)=pzzh(i,1,k)
	pxyh(i,ny_d+1,k)=pxyh(i,1,k)
	pxzh(i,ny_d+1,k)=pxzh(i,1,k)
	pyzh(i,ny_d+1,k)=pyzh(i,1,k)
	  
    return

  end subroutine pbc_edge_y_flop
 
 end module pbc_kernels
