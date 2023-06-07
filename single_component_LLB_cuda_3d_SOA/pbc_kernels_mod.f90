#include "defines.h"
 module pbc_kernels
  
  use cudavars
  
  implicit none
  
 contains
 
 attributes(global) subroutine pbc_side_x()
	  
	integer :: io,ie,j,k,gio,gie,gj,gk,xblocko,xblocke,yblock,zblock,idblocko,idblocke
	
	
	gj = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	if (gj>ny_d .or. gk>nz_d)return
	
	yblock=(gj+TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
    j=gj-yblock*TILE_DIMy_d+TILE_DIMy_d
    k=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
	
	gio=nx_d
	gie=0
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=xblocko+yblock*nxblock_d+zblock*nxyblock_d+1
    
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=xblocke+yblock*nxblock_d+zblock*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    
	  
	rho(ie,j,k,idblocke)=rho(io,j,k,idblocko)
	u(ie,j,k,idblocke)=u(io,j,k,idblocko)
	v(ie,j,k,idblocke)=v(io,j,k,idblocko)
	w(ie,j,k,idblocke)=w(io,j,k,idblocko)
	pxx(ie,j,k,idblocke)=pxx(io,j,k,idblocko)
	pyy(ie,j,k,idblocke)=pyy(io,j,k,idblocko)
	pzz(ie,j,k,idblocke)=pzz(io,j,k,idblocko)
	pxy(ie,j,k,idblocke)=pxy(io,j,k,idblocko)
	pxz(ie,j,k,idblocke)=pxz(io,j,k,idblocko)
	pyz(ie,j,k,idblocke)=pyz(io,j,k,idblocko)
	
	gio=1
	gie=nx_d+1
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=xblocko+yblock*nxblock_d+zblock*nxyblock_d+1
    
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=xblocke+yblock*nxblock_d+zblock*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
	
	rho(ie,j,k,idblocke)=rho(io,j,k,idblocko)
	u(ie,j,k,idblocke)=u(io,j,k,idblocko)
	v(ie,j,k,idblocke)=v(io,j,k,idblocko)
	w(ie,j,k,idblocke)=w(io,j,k,idblocko)
	pxx(ie,j,k,idblocke)=pxx(io,j,k,idblocko)
	pyy(ie,j,k,idblocke)=pyy(io,j,k,idblocko)
	pzz(ie,j,k,idblocke)=pzz(io,j,k,idblocko)
	pxy(ie,j,k,idblocke)=pxy(io,j,k,idblocko)
	pxz(ie,j,k,idblocke)=pxz(io,j,k,idblocko)
	pyz(ie,j,k,idblocke)=pyz(io,j,k,idblocko)
	  
    return
		  
  end subroutine pbc_side_x
  
  attributes(global) subroutine pbc_side_x_flop()
	  
	integer :: io,ie,j,k,gio,gie,gj,gk,xblocko,xblocke,yblock,zblock,idblocko,idblocke
	
	
	gj = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	if (gj>ny_d .or. gk>nz_d)return
	
	yblock=(gj+TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
    j=gj-yblock*TILE_DIMy_d+TILE_DIMy_d
    k=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
	
	gio=nx_d
	gie=0
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=xblocko+yblock*nxblock_d+zblock*nxyblock_d+1
    
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=xblocke+yblock*nxblock_d+zblock*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    
	  
	rhoh(ie,j,k,idblocke)=rhoh(io,j,k,idblocko)
	uh(ie,j,k,idblocke)=uh(io,j,k,idblocko)
	vh(ie,j,k,idblocke)=vh(io,j,k,idblocko)
	wh(ie,j,k,idblocke)=wh(io,j,k,idblocko)
	pxxh(ie,j,k,idblocke)=pxxh(io,j,k,idblocko)
	pyyh(ie,j,k,idblocke)=pyyh(io,j,k,idblocko)
	pzzh(ie,j,k,idblocke)=pzzh(io,j,k,idblocko)
	pxyh(ie,j,k,idblocke)=pxyh(io,j,k,idblocko)
	pxzh(ie,j,k,idblocke)=pxzh(io,j,k,idblocko)
	pyzh(ie,j,k,idblocke)=pyzh(io,j,k,idblocko)
	
	gio=1
	gie=nx_d+1
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=xblocko+yblock*nxblock_d+zblock*nxyblock_d+1
    
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=xblocke+yblock*nxblock_d+zblock*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
	
	rhoh(ie,j,k,idblocke)=rhoh(io,j,k,idblocko)
	uh(ie,j,k,idblocke)=uh(io,j,k,idblocko)
	vh(ie,j,k,idblocke)=vh(io,j,k,idblocko)
	wh(ie,j,k,idblocke)=wh(io,j,k,idblocko)
	pxxh(ie,j,k,idblocke)=pxxh(io,j,k,idblocko)
	pyyh(ie,j,k,idblocke)=pyyh(io,j,k,idblocko)
	pzzh(ie,j,k,idblocke)=pzzh(io,j,k,idblocko)
	pxyh(ie,j,k,idblocke)=pxyh(io,j,k,idblocko)
	pxzh(ie,j,k,idblocke)=pxzh(io,j,k,idblocko)
	pyzh(ie,j,k,idblocke)=pyzh(io,j,k,idblocko)
	  
    return
		  
  end subroutine pbc_side_x_flop
  
  attributes(global) subroutine pbc_edge_x()
	  
    integer :: i,jo,je,ko,ke,gi,gjo,gje,gko,gke,xblock,yblocko,yblocke,zblocko,zblocke,idblocko,idblocke
	
	
	gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (gi>nx_d)return
	
	xblock=(gi+TILE_DIMx_d-1)/TILE_DIMx_d
	i=gi-xblock*TILE_DIMx_d+TILE_DIMx_d
	
	gjo=ny_d
	gje=0
	gko=nz_d
	gke=0
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	ko=gko-zblocko*TILE_DIMz_d+TILE_DIMz_d
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	ke=gke-zblocke*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblocko*nxblock_d+zblocko*nxyblock_d+1
	idblocke=xblock+yblocke*nxblock_d+zblocke*nxyblock_d+1
	
	rho(i,je,ke,idblocke)=rho(i,jo,ko,idblocko)
	u(i,je,ke,idblocke)=u(i,jo,ko,idblocko)
	v(i,je,ke,idblocke)=v(i,jo,ko,idblocko)
	w(i,je,ke,idblocke)=w(i,jo,ko,idblocko)
	pxx(i,je,ke,idblocke)=pxx(i,jo,ko,idblocko)
	pyy(i,je,ke,idblocke)=pyy(i,jo,ko,idblocko)
	pzz(i,je,ke,idblocke)=pzz(i,jo,ko,idblocko)
	pxy(i,je,ke,idblocke)=pxy(i,jo,ko,idblocko)
	pxz(i,je,ke,idblocke)=pxz(i,jo,ko,idblocko)
	pyz(i,je,ke,idblocke)=pyz(i,jo,ko,idblocko)
	
	gjo=1
	gje=ny_d+1
	gko=1
	gke=nz_d+1
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	ko=gko-zblocko*TILE_DIMz_d+TILE_DIMz_d
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	ke=gke-zblocke*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblocko*nxblock_d+zblocko*nxyblock_d+1
	idblocke=xblock+yblocke*nxblock_d+zblocke*nxyblock_d+1
	
	rho(i,je,ke,idblocke)=rho(i,jo,ko,idblocko)
	u(i,je,ke,idblocke)=u(i,jo,ko,idblocko)
	v(i,je,ke,idblocke)=v(i,jo,ko,idblocko)
	w(i,je,ke,idblocke)=w(i,jo,ko,idblocko)
	pxx(i,je,ke,idblocke)=pxx(i,jo,ko,idblocko)
	pyy(i,je,ke,idblocke)=pyy(i,jo,ko,idblocko)
	pzz(i,je,ke,idblocke)=pzz(i,jo,ko,idblocko)
	pxy(i,je,ke,idblocke)=pxy(i,jo,ko,idblocko)
	pxz(i,je,ke,idblocke)=pxz(i,jo,ko,idblocko)
	pyz(i,je,ke,idblocke)=pyz(i,jo,ko,idblocko)
	
	gjo=ny_d
	gje=0
	gko=1
	gke=nz_d+1
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	ko=gko-zblocko*TILE_DIMz_d+TILE_DIMz_d
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	ke=gke-zblocke*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblocko*nxblock_d+zblocko*nxyblock_d+1
	idblocke=xblock+yblocke*nxblock_d+zblocke*nxyblock_d+1
	
	rho(i,je,ke,idblocke)=rho(i,jo,ko,idblocko)
	u(i,je,ke,idblocke)=u(i,jo,ko,idblocko)
	v(i,je,ke,idblocke)=v(i,jo,ko,idblocko)
	w(i,je,ke,idblocke)=w(i,jo,ko,idblocko)
	pxx(i,je,ke,idblocke)=pxx(i,jo,ko,idblocko)
	pyy(i,je,ke,idblocke)=pyy(i,jo,ko,idblocko)
	pzz(i,je,ke,idblocke)=pzz(i,jo,ko,idblocko)
	pxy(i,je,ke,idblocke)=pxy(i,jo,ko,idblocko)
	pxz(i,je,ke,idblocke)=pxz(i,jo,ko,idblocko)
	pyz(i,je,ke,idblocke)=pyz(i,jo,ko,idblocko)
	
	gjo=1
	gje=ny_d+1
	gko=nz_d
	gke=0
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	ko=gko-zblocko*TILE_DIMz_d+TILE_DIMz_d
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	ke=gke-zblocke*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblocko*nxblock_d+zblocko*nxyblock_d+1
	idblocke=xblock+yblocke*nxblock_d+zblocke*nxyblock_d+1
	
	rho(i,je,ke,idblocke)=rho(i,jo,ko,idblocko)
	u(i,je,ke,idblocke)=u(i,jo,ko,idblocko)
	v(i,je,ke,idblocke)=v(i,jo,ko,idblocko)
	w(i,je,ke,idblocke)=w(i,jo,ko,idblocko)
	pxx(i,je,ke,idblocke)=pxx(i,jo,ko,idblocko)
	pyy(i,je,ke,idblocke)=pyy(i,jo,ko,idblocko)
	pzz(i,je,ke,idblocke)=pzz(i,jo,ko,idblocko)
	pxy(i,je,ke,idblocke)=pxy(i,jo,ko,idblocko)
	pxz(i,je,ke,idblocke)=pxz(i,jo,ko,idblocko)
	pyz(i,je,ke,idblocke)=pyz(i,jo,ko,idblocko)
	
	
    return
		  
  end subroutine pbc_edge_x
  
  attributes(global) subroutine pbc_edge_x_flop()
	  
    integer :: i,jo,je,ko,ke,gi,gjo,gje,gko,gke,xblock,yblocko,yblocke,zblocko,zblocke,idblocko,idblocke
	
	
	gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (gi>nx_d)return
	
	xblock=(gi+TILE_DIMx_d-1)/TILE_DIMx_d
	i=gi-xblock*TILE_DIMx_d+TILE_DIMx_d
	
	gjo=ny_d
	gje=0
	gko=nz_d
	gke=0
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	ko=gko-zblocko*TILE_DIMz_d+TILE_DIMz_d
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	ke=gke-zblocke*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblocko*nxblock_d+zblocko*nxyblock_d+1
	idblocke=xblock+yblocke*nxblock_d+zblocke*nxyblock_d+1
	
	rhoh(i,je,ke,idblocke)=rhoh(i,jo,ko,idblocko)
	uh(i,je,ke,idblocke)=uh(i,jo,ko,idblocko)
	vh(i,je,ke,idblocke)=vh(i,jo,ko,idblocko)
	wh(i,je,ke,idblocke)=wh(i,jo,ko,idblocko)
	pxxh(i,je,ke,idblocke)=pxxh(i,jo,ko,idblocko)
	pyyh(i,je,ke,idblocke)=pyyh(i,jo,ko,idblocko)
	pzzh(i,je,ke,idblocke)=pzzh(i,jo,ko,idblocko)
	pxyh(i,je,ke,idblocke)=pxyh(i,jo,ko,idblocko)
	pxzh(i,je,ke,idblocke)=pxzh(i,jo,ko,idblocko)
	pyzh(i,je,ke,idblocke)=pyzh(i,jo,ko,idblocko)
	
	gjo=1
	gje=ny_d+1
	gko=1
	gke=nz_d+1
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	ko=gko-zblocko*TILE_DIMz_d+TILE_DIMz_d
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	ke=gke-zblocke*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblocko*nxblock_d+zblocko*nxyblock_d+1
	idblocke=xblock+yblocke*nxblock_d+zblocke*nxyblock_d+1
	
	rhoh(i,je,ke,idblocke)=rhoh(i,jo,ko,idblocko)
	uh(i,je,ke,idblocke)=uh(i,jo,ko,idblocko)
	vh(i,je,ke,idblocke)=vh(i,jo,ko,idblocko)
	wh(i,je,ke,idblocke)=wh(i,jo,ko,idblocko)
	pxxh(i,je,ke,idblocke)=pxxh(i,jo,ko,idblocko)
	pyyh(i,je,ke,idblocke)=pyyh(i,jo,ko,idblocko)
	pzzh(i,je,ke,idblocke)=pzzh(i,jo,ko,idblocko)
	pxyh(i,je,ke,idblocke)=pxyh(i,jo,ko,idblocko)
	pxzh(i,je,ke,idblocke)=pxzh(i,jo,ko,idblocko)
	pyzh(i,je,ke,idblocke)=pyzh(i,jo,ko,idblocko)
	
	gjo=ny_d
	gje=0
	gko=1
	gke=nz_d+1
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	ko=gko-zblocko*TILE_DIMz_d+TILE_DIMz_d
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	ke=gke-zblocke*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblocko*nxblock_d+zblocko*nxyblock_d+1
	idblocke=xblock+yblocke*nxblock_d+zblocke*nxyblock_d+1
	
	rhoh(i,je,ke,idblocke)=rhoh(i,jo,ko,idblocko)
	uh(i,je,ke,idblocke)=uh(i,jo,ko,idblocko)
	vh(i,je,ke,idblocke)=vh(i,jo,ko,idblocko)
	wh(i,je,ke,idblocke)=wh(i,jo,ko,idblocko)
	pxxh(i,je,ke,idblocke)=pxxh(i,jo,ko,idblocko)
	pyyh(i,je,ke,idblocke)=pyyh(i,jo,ko,idblocko)
	pzzh(i,je,ke,idblocke)=pzzh(i,jo,ko,idblocko)
	pxyh(i,je,ke,idblocke)=pxyh(i,jo,ko,idblocko)
	pxzh(i,je,ke,idblocke)=pxzh(i,jo,ko,idblocko)
	pyzh(i,je,ke,idblocke)=pyzh(i,jo,ko,idblocko)
	
	gjo=1
	gje=ny_d+1
	gko=nz_d
	gke=0
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	ko=gko-zblocko*TILE_DIMz_d+TILE_DIMz_d
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	ke=gke-zblocke*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblocko*nxblock_d+zblocko*nxyblock_d+1
	idblocke=xblock+yblocke*nxblock_d+zblocke*nxyblock_d+1
	
	rhoh(i,je,ke,idblocke)=rhoh(i,jo,ko,idblocko)
	uh(i,je,ke,idblocke)=uh(i,jo,ko,idblocko)
	vh(i,je,ke,idblocke)=vh(i,jo,ko,idblocko)
	wh(i,je,ke,idblocke)=wh(i,jo,ko,idblocko)
	pxxh(i,je,ke,idblocke)=pxxh(i,jo,ko,idblocko)
	pyyh(i,je,ke,idblocke)=pyyh(i,jo,ko,idblocko)
	pzzh(i,je,ke,idblocke)=pzzh(i,jo,ko,idblocko)
	pxyh(i,je,ke,idblocke)=pxyh(i,jo,ko,idblocko)
	pxzh(i,je,ke,idblocke)=pxzh(i,jo,ko,idblocko)
	pyzh(i,je,ke,idblocke)=pyzh(i,jo,ko,idblocko)
	
	
    return
		  
  end subroutine pbc_edge_x_flop

  attributes(global) subroutine pbc_side_y()
	
	implicit none
	integer :: i,jo,je,k,gi,gjo,gje,gk,xblock,yblocko,yblocke,zblock,idblocko,idblocke
	  
	gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

	if (gi>nx_d .or. gk>nz_d)return
	
	xblock=(gi+TILE_DIMx_d-1)/TILE_DIMx_d
    zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
    i=gi-xblock*TILE_DIMx_d+TILE_DIMx_d
    k=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
    
    gjo=ny_d
	gje=0
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=xblock+yblocko*nxblock_d+zblock*nxyblock_d+1
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=xblock+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	
	rho(i,je,k,idblocke)=rho(i,jo,k,idblocko)
	u(i,je,k,idblocke)=u(i,jo,k,idblocko)
	v(i,je,k,idblocke)=v(i,jo,k,idblocko)
	w(i,je,k,idblocke)=w(i,jo,k,idblocko)
	pxx(i,je,k,idblocke)=pxx(i,jo,k,idblocko)
	pyy(i,je,k,idblocke)=pyy(i,jo,k,idblocko)
	pzz(i,je,k,idblocke)=pzz(i,jo,k,idblocko)
	pxy(i,je,k,idblocke)=pxy(i,jo,k,idblocko)
	pxz(i,je,k,idblocke)=pxz(i,jo,k,idblocko)
	pyz(i,je,k,idblocke)=pyz(i,jo,k,idblocko)
	
	gjo=1
	gje=ny_d+1
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=xblock+yblocko*nxblock_d+zblock*nxyblock_d+1
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=xblock+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	
	rho(i,je,k,idblocke)=rho(i,jo,k,idblocko)
	u(i,je,k,idblocke)=u(i,jo,k,idblocko)
	v(i,je,k,idblocke)=v(i,jo,k,idblocko)
	w(i,je,k,idblocke)=w(i,jo,k,idblocko)
	pxx(i,je,k,idblocke)=pxx(i,jo,k,idblocko)
	pyy(i,je,k,idblocke)=pyy(i,jo,k,idblocko)
	pzz(i,je,k,idblocke)=pzz(i,jo,k,idblocko)
	pxy(i,je,k,idblocke)=pxy(i,jo,k,idblocko)
	pxz(i,je,k,idblocke)=pxz(i,jo,k,idblocko)
	pyz(i,je,k,idblocke)=pyz(i,jo,k,idblocko)
	  
    return

  end subroutine pbc_side_y
  
  attributes(global) subroutine pbc_side_y_flop()
	
	implicit none
	
	integer :: i,jo,je,k,gi,gjo,gje,gk,xblock,yblocko,yblocke,zblock,idblocko,idblocke
	  
	gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

	if (gi>nx_d .or. gk>nz_d)return
	
	xblock=(gi+TILE_DIMx_d-1)/TILE_DIMx_d
    zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
    i=gi-xblock*TILE_DIMx_d+TILE_DIMx_d
    k=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
    
    gjo=ny_d
	gje=0
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=xblock+yblocko*nxblock_d+zblock*nxyblock_d+1
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=xblock+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	
	rhoh(i,je,k,idblocke)=rhoh(i,jo,k,idblocko)
	uh(i,je,k,idblocke)=uh(i,jo,k,idblocko)
	vh(i,je,k,idblocke)=vh(i,jo,k,idblocko)
	wh(i,je,k,idblocke)=wh(i,jo,k,idblocko)
	pxxh(i,je,k,idblocke)=pxxh(i,jo,k,idblocko)
	pyyh(i,je,k,idblocke)=pyyh(i,jo,k,idblocko)
	pzzh(i,je,k,idblocke)=pzzh(i,jo,k,idblocko)
	pxyh(i,je,k,idblocke)=pxyh(i,jo,k,idblocko)
	pxzh(i,je,k,idblocke)=pxzh(i,jo,k,idblocko)
	pyzh(i,je,k,idblocke)=pyzh(i,jo,k,idblocko)
	
	gjo=1
	gje=ny_d+1
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=xblock+yblocko*nxblock_d+zblock*nxyblock_d+1
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=xblock+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	
	rhoh(i,je,k,idblocke)=rhoh(i,jo,k,idblocko)
	uh(i,je,k,idblocke)=uh(i,jo,k,idblocko)
	vh(i,je,k,idblocke)=vh(i,jo,k,idblocko)
	wh(i,je,k,idblocke)=wh(i,jo,k,idblocko)
	pxxh(i,je,k,idblocke)=pxxh(i,jo,k,idblocko)
	pyyh(i,je,k,idblocke)=pyyh(i,jo,k,idblocko)
	pzzh(i,je,k,idblocke)=pzzh(i,jo,k,idblocko)
	pxyh(i,je,k,idblocke)=pxyh(i,jo,k,idblocko)
	pxzh(i,je,k,idblocke)=pxzh(i,jo,k,idblocko)
	pyzh(i,je,k,idblocke)=pyzh(i,jo,k,idblocko)
	  
    return

  end subroutine pbc_side_y_flop
  
  attributes(global) subroutine pbc_edge_z()
	  
    integer :: io,ie,jo,je,k,gio,gie,gjo,gje,gk,xblocko,xblocke,yblocko,yblocke,zblock,idblocko,idblocke
	
	
	gk = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (gk>nz_d)return
	
	zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
	k=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
	
	gio=nx_d
	gie=0
	gjo=ny_d
	gje=0
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMx_d
	
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    
    ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMx_d
	
    
    idblocko=xblocko+yblocko*nxblock_d+zblock*nxyblock_d+1
	idblocke=xblocke+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	rho(ie,je,k,idblocke)=rho(io,jo,k,idblocko)
	u(ie,je,k,idblocke)=u(io,jo,k,idblocko)
	v(ie,je,k,idblocke)=v(io,jo,k,idblocko)
	w(ie,je,k,idblocke)=w(io,jo,k,idblocko)
	pxx(ie,je,k,idblocke)=pxx(io,jo,k,idblocko)
	pyy(ie,je,k,idblocke)=pyy(io,jo,k,idblocko)
	pzz(ie,je,k,idblocke)=pzz(io,jo,k,idblocko)
	pxy(ie,je,k,idblocke)=pxy(io,jo,k,idblocko)
	pxz(ie,je,k,idblocke)=pxz(io,jo,k,idblocko)
	pyz(ie,je,k,idblocke)=pyz(io,jo,k,idblocko)
	
	gio=1
	gie=nx_d+1
	gjo=1
	gje=ny_d+1
	
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMx_d
	
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    
    ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMx_d
	
    
    idblocko=xblocko+yblocko*nxblock_d+zblock*nxyblock_d+1
	idblocke=xblocke+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	rho(ie,je,k,idblocke)=rho(io,jo,k,idblocko)
	u(ie,je,k,idblocke)=u(io,jo,k,idblocko)
	v(ie,je,k,idblocke)=v(io,jo,k,idblocko)
	w(ie,je,k,idblocke)=w(io,jo,k,idblocko)
	pxx(ie,je,k,idblocke)=pxx(io,jo,k,idblocko)
	pyy(ie,je,k,idblocke)=pyy(io,jo,k,idblocko)
	pzz(ie,je,k,idblocke)=pzz(io,jo,k,idblocko)
	pxy(ie,je,k,idblocke)=pxy(io,jo,k,idblocko)
	pxz(ie,je,k,idblocke)=pxz(io,jo,k,idblocko)
	pyz(ie,je,k,idblocke)=pyz(io,jo,k,idblocko)
	
	gio=1
	gie=nx_d+1
	gjo=ny_d
	gje=0
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMx_d
	
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    
    ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMx_d
	
    
    idblocko=xblocko+yblocko*nxblock_d+zblock*nxyblock_d+1
	idblocke=xblocke+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	rho(ie,je,k,idblocke)=rho(io,jo,k,idblocko)
	u(ie,je,k,idblocke)=u(io,jo,k,idblocko)
	v(ie,je,k,idblocke)=v(io,jo,k,idblocko)
	w(ie,je,k,idblocke)=w(io,jo,k,idblocko)
	pxx(ie,je,k,idblocke)=pxx(io,jo,k,idblocko)
	pyy(ie,je,k,idblocke)=pyy(io,jo,k,idblocko)
	pzz(ie,je,k,idblocke)=pzz(io,jo,k,idblocko)
	pxy(ie,je,k,idblocke)=pxy(io,jo,k,idblocko)
	pxz(ie,je,k,idblocke)=pxz(io,jo,k,idblocko)
	pyz(ie,je,k,idblocke)=pyz(io,jo,k,idblocko)
	
	gio=nx_d
	gie=0
	gjo=1
	gje=ny_d+1
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMx_d
	
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    
    ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMx_d
	
    
    idblocko=xblocko+yblocko*nxblock_d+zblock*nxyblock_d+1
	idblocke=xblocke+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	rho(ie,je,k,idblocke)=rho(io,jo,k,idblocko)
	u(ie,je,k,idblocke)=u(io,jo,k,idblocko)
	v(ie,je,k,idblocke)=v(io,jo,k,idblocko)
	w(ie,je,k,idblocke)=w(io,jo,k,idblocko)
	pxx(ie,je,k,idblocke)=pxx(io,jo,k,idblocko)
	pyy(ie,je,k,idblocke)=pyy(io,jo,k,idblocko)
	pzz(ie,je,k,idblocke)=pzz(io,jo,k,idblocko)
	pxy(ie,je,k,idblocke)=pxy(io,jo,k,idblocko)
	pxz(ie,je,k,idblocke)=pxz(io,jo,k,idblocko)
	pyz(ie,je,k,idblocke)=pyz(io,jo,k,idblocko)
	
	
    return
		  
  end subroutine pbc_edge_z
  
  attributes(global) subroutine pbc_edge_z_flop()
	  
    integer :: io,ie,jo,je,k,gio,gie,gjo,gje,gk,xblocko,xblocke,yblocko,yblocke,zblock,idblocko,idblocke
	
	
	gk = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (gk>nz_d)return
	
	zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
	k=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
	
	gio=nx_d
	gie=0
	gjo=ny_d
	gje=0
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMx_d
	
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    
    ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMx_d
	
    
    idblocko=xblocko+yblocko*nxblock_d+zblock*nxyblock_d+1
	idblocke=xblocke+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	rhoh(ie,je,k,idblocke)=rhoh(io,jo,k,idblocko)
	uh(ie,je,k,idblocke)=uh(io,jo,k,idblocko)
	vh(ie,je,k,idblocke)=vh(io,jo,k,idblocko)
	wh(ie,je,k,idblocke)=wh(io,jo,k,idblocko)
	pxxh(ie,je,k,idblocke)=pxxh(io,jo,k,idblocko)
	pyyh(ie,je,k,idblocke)=pyyh(io,jo,k,idblocko)
	pzzh(ie,je,k,idblocke)=pzzh(io,jo,k,idblocko)
	pxyh(ie,je,k,idblocke)=pxyh(io,jo,k,idblocko)
	pxzh(ie,je,k,idblocke)=pxzh(io,jo,k,idblocko)
	pyzh(ie,je,k,idblocke)=pyzh(io,jo,k,idblocko)
	
	gio=1
	gie=nx_d+1
	gjo=1
	gje=ny_d+1
	
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMx_d
	
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    
    ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMx_d
	
    
    idblocko=xblocko+yblocko*nxblock_d+zblock*nxyblock_d+1
	idblocke=xblocke+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	rhoh(ie,je,k,idblocke)=rhoh(io,jo,k,idblocko)
	uh(ie,je,k,idblocke)=uh(io,jo,k,idblocko)
	vh(ie,je,k,idblocke)=vh(io,jo,k,idblocko)
	wh(ie,je,k,idblocke)=wh(io,jo,k,idblocko)
	pxxh(ie,je,k,idblocke)=pxxh(io,jo,k,idblocko)
	pyyh(ie,je,k,idblocke)=pyyh(io,jo,k,idblocko)
	pzzh(ie,je,k,idblocke)=pzzh(io,jo,k,idblocko)
	pxyh(ie,je,k,idblocke)=pxyh(io,jo,k,idblocko)
	pxzh(ie,je,k,idblocke)=pxzh(io,jo,k,idblocko)
	pyzh(ie,je,k,idblocke)=pyzh(io,jo,k,idblocko)
	
	gio=1
	gie=nx_d+1
	gjo=ny_d
	gje=0
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMx_d
	
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    
    ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMx_d
	
    
    idblocko=xblocko+yblocko*nxblock_d+zblock*nxyblock_d+1
	idblocke=xblocke+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	rhoh(ie,je,k,idblocke)=rhoh(io,jo,k,idblocko)
	uh(ie,je,k,idblocke)=uh(io,jo,k,idblocko)
	vh(ie,je,k,idblocke)=vh(io,jo,k,idblocko)
	wh(ie,je,k,idblocke)=wh(io,jo,k,idblocko)
	pxxh(ie,je,k,idblocke)=pxxh(io,jo,k,idblocko)
	pyyh(ie,je,k,idblocke)=pyyh(io,jo,k,idblocko)
	pzzh(ie,je,k,idblocke)=pzzh(io,jo,k,idblocko)
	pxyh(ie,je,k,idblocke)=pxyh(io,jo,k,idblocko)
	pxzh(ie,je,k,idblocke)=pxzh(io,jo,k,idblocko)
	pyzh(ie,je,k,idblocke)=pyzh(io,jo,k,idblocko)
	
	gio=nx_d
	gie=0
	gjo=1
	gje=ny_d+1
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMx_d
	
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
    
    ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+TILE_DIMx_d
	
    
    idblocko=xblocko+yblocko*nxblock_d+zblock*nxyblock_d+1
	idblocke=xblocke+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	rhoh(ie,je,k,idblocke)=rhoh(io,jo,k,idblocko)
	uh(ie,je,k,idblocke)=uh(io,jo,k,idblocko)
	vh(ie,je,k,idblocke)=vh(io,jo,k,idblocko)
	wh(ie,je,k,idblocke)=wh(io,jo,k,idblocko)
	pxxh(ie,je,k,idblocke)=pxxh(io,jo,k,idblocko)
	pyyh(ie,je,k,idblocke)=pyyh(io,jo,k,idblocko)
	pzzh(ie,je,k,idblocke)=pzzh(io,jo,k,idblocko)
	pxyh(ie,je,k,idblocke)=pxyh(io,jo,k,idblocko)
	pxzh(ie,je,k,idblocke)=pxzh(io,jo,k,idblocko)
	pyzh(ie,je,k,idblocke)=pyzh(io,jo,k,idblocko)
	
	
    return
		  
  end subroutine pbc_edge_z_flop
 
 end module pbc_kernels
