#include "defines.h"
 module pbc_kernels
  
  use cudavars
  
  implicit none
  
 contains
 
 attributes(global) subroutine pbc_edge_x()
	  
	integer :: i,j,k,ii,jj,kk,xblock,yblock,zblock,idblock
	
	j = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	if (j>ny_d .or. k>nz_d)return
	
	i=1
	xblock=(i+TILE_DIMx_d-1)/TILE_DIMx_d
    yblock=(j+TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(k+TILE_DIMz_d-1)/TILE_DIMz_d
	idblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
	ii=i-xblock*TILE_DIMx_d+TILE_DIMx_d
    jj=j-yblock*TILE_DIMy_d+TILE_DIMy_d
    kk=k-zblock*TILE_DIMz_d+TILE_DIMz_d
	  
	rho(0,j,k,idblock)=rho(nx_d,j,k,idblock)
	u(0,j,k,idblock)=u(nx_d,j,k,idblock)
	v(0,j,k,idblock)=v(nx_d,j,k,idblock)
	w(0,j,k,idblock)=w(nx_d,j,k,idblock)
	pxx(0,j,k,idblock)=pxx(nx_d,j,k,idblock)
	pyy(0,j,k,idblock)=pyy(nx_d,j,k,idblock)
	pzz(0,j,k,idblock)=pzz(nx_d,j,k,idblock)
	pxy(0,j,k,idblock)=pxy(nx_d,j,k,idblock)
	pxz(0,j,k,idblock)=pxz(nx_d,j,k,idblock)
	pyz(0,j,k,idblock)=pyz(nx_d,j,k,idblock)
	
	rho(nx_d+1,j,k,idblock)=rho(1,j,k,idblock)
	u(nx_d+1,j,k,idblock)=u(1,j,k,idblock)
	v(nx_d+1,j,k,idblock)=v(1,j,k,idblock)
	w(nx_d+1,j,k,idblock)=w(1,j,k,idblock)
	pxx(nx_d+1,j,k,idblock)=pxx(1,j,k,idblock)
	pyy(nx_d+1,j,k,idblock)=pyy(1,j,k,idblock)
	pzz(nx_d+1,j,k,idblock)=pzz(1,j,k,idblock)
	pxy(nx_d+1,j,k,idblock)=pxy(1,j,k,idblock)
	pxz(nx_d+1,j,k,idblock)=pxz(1,j,k,idblock)
	pyz(nx_d+1,j,k,idblock)=pyz(1,j,k,idblock)
	  
    return
		  
  end subroutine pbc_edge_x
  
  attributes(global) subroutine pbc_edge_x_flop()
	  
	integer :: j,k,ii,jj,kk,xblock,yblock,zblock,idblock
	  
	j = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	if (j>ny_d .or. k>nz_d)return
	  
	rhoh(0,j,k,idblock)=rhoh(nx_d,j,k,idblock)
	uh(0,j,k,idblock)=uh(nx_d,j,k,idblock)
	vh(0,j,k,idblock)=vh(nx_d,j,k,idblock)
	wh(0,j,k,idblock)=wh(nx_d,j,k,idblock)
	pxxh(0,j,k,idblock)=pxxh(nx_d,j,k,idblock)
	pyyh(0,j,k,idblock)=pyyh(nx_d,j,k,idblock)
	pzzh(0,j,k,idblock)=pzzh(nx_d,j,k,idblock)
	pxyh(0,j,k,idblock)=pxyh(nx_d,j,k,idblock)
	pxzh(0,j,k,idblock)=pxzh(nx_d,j,k,idblock)
	pyzh(0,j,k,idblock)=pyzh(nx_d,j,k,idblock)
	
	rhoh(nx_d+1,j,k,idblock)=rhoh(1,j,k,idblock)
	uh(nx_d+1,j,k,idblock)=uh(1,j,k,idblock)
	vh(nx_d+1,j,k,idblock)=vh(1,j,k,idblock)
	wh(nx_d+1,j,k,idblock)=wh(1,j,k,idblock)
	pxxh(nx_d+1,j,k,idblock)=pxxh(1,j,k,idblock)
	pyyh(nx_d+1,j,k,idblock)=pyyh(1,j,k,idblock)
	pzzh(nx_d+1,j,k,idblock)=pzzh(1,j,k,idblock)
	pxyh(nx_d+1,j,k,idblock)=pxyh(1,j,k,idblock)
	pxzh(nx_d+1,j,k,idblock)=pxzh(1,j,k,idblock)
	pyzh(nx_d+1,j,k,idblock)=pyzh(1,j,k,idblock)
	  
    return
		  
  end subroutine pbc_edge_x_flop

  attributes(global) subroutine pbc_edge_y()
	
	implicit none
	integer :: i,k,ii,jj,kk,xblock,yblock,zblock,idblock
	  
	i = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

	if (i>nx_d .or. k>nz_d)return
	
	rho(i,0,k,idblock)=rho(i,ny_d,k,idblock)
	u(i,0,k,idblock)=u(i,ny_d,k,idblock)
	v(i,0,k,idblock)=v(i,ny_d,k,idblock)
	w(i,0,k,idblock)=w(i,ny_d,k,idblock)
	pxx(i,0,k,idblock)=pxx(i,ny_d,k,idblock)
	pyy(i,0,k,idblock)=pyy(i,ny_d,k,idblock)
	pzz(i,0,k,idblock)=pzz(i,ny_d,k,idblock)
	pxy(i,0,k,idblock)=pxy(i,ny_d,k,idblock)
	pxz(i,0,k,idblock)=pxz(i,ny_d,k,idblock)
	pyz(i,0,k,idblock)=pyz(i,ny_d,k,idblock)
	
	rho(i,ny_d+1,k,idblock)=rho(i,1,k,idblock)
	u(i,ny_d+1,k,idblock)=u(i,1,k,idblock)
	v(i,ny_d+1,k,idblock)=v(i,1,k,idblock)
	w(i,ny_d+1,k,idblock)=w(i,1,k,idblock)
	pxx(i,ny_d+1,k,idblock)=pxx(i,1,k,idblock)
	pyy(i,ny_d+1,k,idblock)=pyy(i,1,k,idblock)
	pzz(i,ny_d+1,k,idblock)=pzz(i,1,k,idblock)
	pxy(i,ny_d+1,k,idblock)=pxy(i,1,k,idblock)
	pxz(i,ny_d+1,k,idblock)=pxz(i,1,k,idblock)
	pyz(i,ny_d+1,k,idblock)=pyz(i,1,k,idblock)
	  
    return

  end subroutine pbc_edge_y
  
  attributes(global) subroutine pbc_edge_y_flop()
	
	implicit none
	
	integer :: i,k,ii,jj,kk,xblock,yblock,zblock,idblock
	  
	i = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

	if (i>nx_d .or. k>nz_d)return
	
	rhoh(i,0,k,idblock)=rhoh(i,ny_d,k,idblock)
	uh(i,0,k,idblock)=uh(i,ny_d,k,idblock)
	vh(i,0,k,idblock)=vh(i,ny_d,k,idblock)
	wh(i,0,k,idblock)=wh(i,ny_d,k,idblock)
	pxxh(i,0,k,idblock)=pxxh(i,ny_d,k,idblock)
	pyyh(i,0,k,idblock)=pyyh(i,ny_d,k,idblock)
	pzzh(i,0,k,idblock)=pzzh(i,ny_d,k,idblock)
	pxyh(i,0,k,idblock)=pxyh(i,ny_d,k,idblock)
	pxzh(i,0,k,idblock)=pxzh(i,ny_d,k,idblock)
	pyzh(i,0,k,idblock)=pyzh(i,ny_d,k,idblock)
	
	rhoh(i,ny_d+1,k,idblock)=rhoh(i,1,k,idblock)
	uh(i,ny_d+1,k,idblock)=uh(i,1,k,idblock)
	vh(i,ny_d+1,k,idblock)=vh(i,1,k,idblock)
	wh(i,ny_d+1,k,idblock)=wh(i,1,k,idblock)
	pxxh(i,ny_d+1,k,idblock)=pxxh(i,1,k,idblock)
	pyyh(i,ny_d+1,k,idblock)=pyyh(i,1,k,idblock)
	pzzh(i,ny_d+1,k,idblock)=pzzh(i,1,k,idblock)
	pxyh(i,ny_d+1,k,idblock)=pxyh(i,1,k,idblock)
	pxzh(i,ny_d+1,k,idblock)=pxzh(i,1,k,idblock)
	pyzh(i,ny_d+1,k,idblock)=pyzh(i,1,k,idblock)
	  
    return

  end subroutine pbc_edge_y_flop
 
 end module pbc_kernels
