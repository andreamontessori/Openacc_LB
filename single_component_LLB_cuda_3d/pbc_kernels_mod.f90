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
  
  attributes(global) subroutine pbc_edge_z()
	  
    integer :: k
	
	
	k = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (k>nz_d)return
	
	
!	gio=nx_d
!	gie=0
!	gjo=ny_d
!	gje=0
	
    rho(0,0,k)=rho(nx_d,ny_d,k)
	u(0,0,k)=u(nx_d,ny_d,k)
	v(0,0,k)=v(nx_d,ny_d,k)
	w(0,0,k)=w(nx_d,ny_d,k)
	pxx(0,0,k)=pxx(nx_d,ny_d,k)
	pyy(0,0,k)=pyy(nx_d,ny_d,k)
	pzz(0,0,k)=pzz(nx_d,ny_d,k)
	pxy(0,0,k)=pxy(nx_d,ny_d,k)
	pxz(0,0,k)=pxz(nx_d,ny_d,k)
	pyz(0,0,k)=pyz(nx_d,ny_d,k)
	
!	gio=1
!	gie=nx_d+1
!	gjo=1
!	gje=ny_d+1
	
	rho(nx_d+1,ny_d+1,k)=rho(1,1,k)
	u(nx_d+1,ny_d+1,k)=u(1,1,k)
	v(nx_d+1,ny_d+1,k)=v(1,1,k)
	w(nx_d+1,ny_d+1,k)=w(1,1,k)
	pxx(nx_d+1,ny_d+1,k)=pxx(1,1,k)
	pyy(nx_d+1,ny_d+1,k)=pyy(1,1,k)
	pzz(nx_d+1,ny_d+1,k)=pzz(1,1,k)
	pxy(nx_d+1,ny_d+1,k)=pxy(1,1,k)
	pxz(nx_d+1,ny_d+1,k)=pxz(1,1,k)
	pyz(nx_d+1,ny_d+1,k)=pyz(1,1,k)
	
!	gio=1
!	gie=nx_d+1
!	gjo=ny_d
!	gje=0
	
	rho(nx_d+1,0,k)=rho(1,ny_d,k)
	u(nx_d+1,0,k)=u(1,ny_d,k)
	v(nx_d+1,0,k)=v(1,ny_d,k)
	w(nx_d+1,0,k)=w(1,ny_d,k)
	pxx(nx_d+1,0,k)=pxx(1,ny_d,k)
	pyy(nx_d+1,0,k)=pyy(1,ny_d,k)
	pzz(nx_d+1,0,k)=pzz(1,ny_d,k)
	pxy(nx_d+1,0,k)=pxy(1,ny_d,k)
	pxz(nx_d+1,0,k)=pxz(1,ny_d,k)
	pyz(nx_d+1,0,k)=pyz(1,ny_d,k)
	
!	gio=nx_d
!	gie=0
!	gjo=1
!	gje=ny_d+1
	
	rho(0,ny_d+1,k)=rho(nx_d,1,k)
	u(0,ny_d+1,k)=u(nx_d,1,k)
	v(0,ny_d+1,k)=v(nx_d,1,k)
	w(0,ny_d+1,k)=w(nx_d,1,k)
	pxx(0,ny_d+1,k)=pxx(nx_d,1,k)
	pyy(0,ny_d+1,k)=pyy(nx_d,1,k)
	pzz(0,ny_d+1,k)=pzz(nx_d,1,k)
	pxy(0,ny_d+1,k)=pxy(nx_d,1,k)
	pxz(0,ny_d+1,k)=pxz(nx_d,1,k)
	pyz(0,ny_d+1,k)=pyz(nx_d,1,k)
	
	
    return
		  
  end subroutine pbc_edge_z
  
  attributes(global) subroutine pbc_edge_z_flop()
	  
    integer :: k
	
	
	k = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (k>nz_d)return
	
	
!	gio=nx_d
!	gie=0
!	gjo=ny_d
!	gje=0
	
    rhoh(0,0,k)=rhoh(nx_d,ny_d,k)
	uh(0,0,k)=uh(nx_d,ny_d,k)
	vh(0,0,k)=vh(nx_d,ny_d,k)
	wh(0,0,k)=wh(nx_d,ny_d,k)
	pxxh(0,0,k)=pxxh(nx_d,ny_d,k)
	pyyh(0,0,k)=pyyh(nx_d,ny_d,k)
	pzzh(0,0,k)=pzzh(nx_d,ny_d,k)
	pxyh(0,0,k)=pxyh(nx_d,ny_d,k)
	pxzh(0,0,k)=pxzh(nx_d,ny_d,k)
	pyzh(0,0,k)=pyzh(nx_d,ny_d,k)
	
!	gio=1
!	gie=nx_d+1
!	gjo=1
!	gje=ny_d+1
	
	rhoh(nx_d+1,ny_d+1,k)=rhoh(1,1,k)
	uh(nx_d+1,ny_d+1,k)=uh(1,1,k)
	vh(nx_d+1,ny_d+1,k)=vh(1,1,k)
	wh(nx_d+1,ny_d+1,k)=wh(1,1,k)
	pxxh(nx_d+1,ny_d+1,k)=pxxh(1,1,k)
	pyyh(nx_d+1,ny_d+1,k)=pyyh(1,1,k)
	pzzh(nx_d+1,ny_d+1,k)=pzzh(1,1,k)
	pxyh(nx_d+1,ny_d+1,k)=pxyh(1,1,k)
	pxzh(nx_d+1,ny_d+1,k)=pxzh(1,1,k)
	pyzh(nx_d+1,ny_d+1,k)=pyzh(1,1,k)
	
!	gio=1
!	gie=nx_d+1
!	gjo=ny_d
!	gje=0
	
	rhoh(nx_d+1,0,k)=rhoh(1,ny_d,k)
	uh(nx_d+1,0,k)=uh(1,ny_d,k)
	vh(nx_d+1,0,k)=vh(1,ny_d,k)
	wh(nx_d+1,0,k)=wh(1,ny_d,k)
	pxxh(nx_d+1,0,k)=pxxh(1,ny_d,k)
	pyyh(nx_d+1,0,k)=pyyh(1,ny_d,k)
	pzzh(nx_d+1,0,k)=pzzh(1,ny_d,k)
	pxyh(nx_d+1,0,k)=pxyh(1,ny_d,k)
	pxzh(nx_d+1,0,k)=pxzh(1,ny_d,k)
	pyzh(nx_d+1,0,k)=pyzh(1,ny_d,k)
	
!	gio=nx_d
!	gie=0
!	gjo=1
!	gje=ny_d+1
	
	rhoh(0,ny_d+1,k)=rhoh(nx_d,1,k)
	uh(0,ny_d+1,k)=uh(nx_d,1,k)
	vh(0,ny_d+1,k)=vh(nx_d,1,k)
	wh(0,ny_d+1,k)=wh(nx_d,1,k)
	pxxh(0,ny_d+1,k)=pxxh(nx_d,1,k)
	pyyh(0,ny_d+1,k)=pyyh(nx_d,1,k)
	pzzh(0,ny_d+1,k)=pzzh(nx_d,1,k)
	pxyh(0,ny_d+1,k)=pxyh(nx_d,1,k)
	pxzh(0,ny_d+1,k)=pxzh(nx_d,1,k)
	pyzh(0,ny_d+1,k)=pyzh(nx_d,1,k)
	
	
    return
		  
  end subroutine pbc_edge_z_flop
  
  attributes(global) subroutine bc_per_x_hvar(step)
  
      integer, value :: step
      
      integer :: j,k
  
      j = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	  if (j>ny_d .or. k>nz_d)return
	  
	  rho(nx_d,j,k)=rho(2,j,k)
	  u(nx_d,j,k)=u(2,j,k)
	  v(nx_d,j,k)=v(2,j,k)
	  w(nx_d,j,k)=w(2,j,k)
	  pxx(nx_d,j,k)=pxx(2,j,k)
	  pyy(nx_d,j,k)=pyy(2,j,k)
	  pzz(nx_d,j,k)=pzz(2,j,k)
	  pxy(nx_d,j,k)=pxy(2,j,k)
	  pxz(nx_d,j,k)=pxz(2,j,k)
	  pyz(nx_d,j,k)=pyz(2,j,k)
	  
	  rho(1,j,k)=rho(nx_d-1,j,k)
	  u(1,j,k)=u(nx_d-1,j,k)
	  v(1,j,k)=v(nx_d-1,j,k)
  	  w(1,j,k)=w(nx_d-1,j,k)
	  pxx(1,j,k)=pxx(nx_d-1,j,k)
	  pyy(1,j,k)=pyy(nx_d-1,j,k)
	  pzz(1,j,k)=pzz(nx_d-1,j,k)
	  pxy(1,j,k)=pxy(nx_d-1,j,k)
	  pxz(1,j,k)=pxz(nx_d-1,j,k)
	  pyz(1,j,k)=pyz(nx_d-1,j,k)
	  
	  
      
      
  end subroutine bc_per_x_hvar
  
  attributes(global) subroutine bc_per_x_hvar_flop(step)
  
      integer, value :: step
      
      integer :: j,k
  
      j = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	  if (j>ny_d .or. k>nz_d)return
	  
	  rhoh(nx_d,j,k)=rhoh(2,j,k)
	  uh(nx_d,j,k)=uh(2,j,k)
	  vh(nx_d,j,k)=vh(2,j,k)
	  wh(nx_d,j,k)=wh(2,j,k)
	  pxxh(nx_d,j,k)=pxxh(2,j,k)
	  pyyh(nx_d,j,k)=pyyh(2,j,k)
	  pzzh(nx_d,j,k)=pzzh(2,j,k)
	  pxyh(nx_d,j,k)=pxyh(2,j,k)
	  pxzh(nx_d,j,k)=pxzh(2,j,k)
	  pyzh(nx_d,j,k)=pyzh(2,j,k)
	  
	  rhoh(1,j,k)=rhoh(nx_d-1,j,k)
	  uh(1,j,k)=uh(nx_d-1,j,k)
	  vh(1,j,k)=vh(nx_d-1,j,k)
  	  wh(1,j,k)=wh(nx_d-1,j,k)
	  pxxh(1,j,k)=pxxh(nx_d-1,j,k)
	  pyyh(1,j,k)=pyyh(nx_d-1,j,k)
	  pzzh(1,j,k)=pzzh(nx_d-1,j,k)
	  pxyh(1,j,k)=pxyh(nx_d-1,j,k)
	  pxzh(1,j,k)=pxzh(nx_d-1,j,k)
	  pyzh(1,j,k)=pyzh(nx_d-1,j,k)
	  
	  
      
      
  end subroutine bc_per_x_hvar_flop
 
 end module pbc_kernels
