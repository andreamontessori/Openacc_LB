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
    
	  
	hfields(ie,j,k,1,idblocke)=hfields(io,j,k,1,idblocko)
	hfields(ie,j,k,2,idblocke)=hfields(io,j,k,2,idblocko)
	hfields(ie,j,k,3,idblocke)=hfields(io,j,k,3,idblocko)
	hfields(ie,j,k,4,idblocke)=hfields(io,j,k,4,idblocko)
	hfields(ie,j,k,5,idblocke)=hfields(io,j,k,5,idblocko)
	hfields(ie,j,k,6,idblocke)=hfields(io,j,k,6,idblocko)
	hfields(ie,j,k,7,idblocke)=hfields(io,j,k,7,idblocko)
	hfields(ie,j,k,8,idblocke)=hfields(io,j,k,8,idblocko)
	hfields(ie,j,k,9,idblocke)=hfields(io,j,k,9,idblocko)
	hfields(ie,j,k,10,idblocke)=hfields(io,j,k,10,idblocko)
	
	gio=1
	gie=nx_d+1
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=xblocko+yblock*nxblock_d+zblock*nxyblock_d+1
    
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=xblocke+yblock*nxblock_d+zblock*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
	
	hfields(ie,j,k,1,idblocke)=hfields(io,j,k,1,idblocko)
	hfields(ie,j,k,2,idblocke)=hfields(io,j,k,2,idblocko)
	hfields(ie,j,k,3,idblocke)=hfields(io,j,k,3,idblocko)
	hfields(ie,j,k,4,idblocke)=hfields(io,j,k,4,idblocko)
	hfields(ie,j,k,5,idblocke)=hfields(io,j,k,5,idblocko)
	hfields(ie,j,k,6,idblocke)=hfields(io,j,k,6,idblocko)
	hfields(ie,j,k,7,idblocke)=hfields(io,j,k,7,idblocko)
	hfields(ie,j,k,8,idblocke)=hfields(io,j,k,8,idblocko)
	hfields(ie,j,k,9,idblocke)=hfields(io,j,k,9,idblocko)
	hfields(ie,j,k,10,idblocke)=hfields(io,j,k,10,idblocko)
	  
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
    
	  
	hfieldsh(ie,j,k,1,idblocke)=hfieldsh(io,j,k,1,idblocko)
	hfieldsh(ie,j,k,2,idblocke)=hfieldsh(io,j,k,2,idblocko)
	hfieldsh(ie,j,k,3,idblocke)=hfieldsh(io,j,k,3,idblocko)
	hfieldsh(ie,j,k,4,idblocke)=hfieldsh(io,j,k,4,idblocko)
	hfieldsh(ie,j,k,5,idblocke)=hfieldsh(io,j,k,5,idblocko)
	hfieldsh(ie,j,k,6,idblocke)=hfieldsh(io,j,k,6,idblocko)
	hfieldsh(ie,j,k,7,idblocke)=hfieldsh(io,j,k,7,idblocko)
	hfieldsh(ie,j,k,8,idblocke)=hfieldsh(io,j,k,8,idblocko)
	hfieldsh(ie,j,k,9,idblocke)=hfieldsh(io,j,k,9,idblocko)
	hfieldsh(ie,j,k,10,idblocke)=hfieldsh(io,j,k,10,idblocko)
	
	gio=1
	gie=nx_d+1
	
	xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=xblocko+yblock*nxblock_d+zblock*nxyblock_d+1
    
    xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=xblocke+yblock*nxblock_d+zblock*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
	
	hfieldsh(ie,j,k,1,idblocke)=hfieldsh(io,j,k,1,idblocko)
	hfieldsh(ie,j,k,2,idblocke)=hfieldsh(io,j,k,2,idblocko)
	hfieldsh(ie,j,k,3,idblocke)=hfieldsh(io,j,k,3,idblocko)
	hfieldsh(ie,j,k,4,idblocke)=hfieldsh(io,j,k,4,idblocko)
	hfieldsh(ie,j,k,5,idblocke)=hfieldsh(io,j,k,5,idblocko)
	hfieldsh(ie,j,k,6,idblocke)=hfieldsh(io,j,k,6,idblocko)
	hfieldsh(ie,j,k,7,idblocke)=hfieldsh(io,j,k,7,idblocko)
	hfieldsh(ie,j,k,8,idblocke)=hfieldsh(io,j,k,8,idblocko)
	hfieldsh(ie,j,k,9,idblocke)=hfieldsh(io,j,k,9,idblocko)
	hfieldsh(ie,j,k,10,idblocke)=hfieldsh(io,j,k,10,idblocko)
	  
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
	
	hfields(i,je,ke,1,idblocke)=hfields(i,jo,ko,1,idblocko)
	hfields(i,je,ke,2,idblocke)=hfields(i,jo,ko,2,idblocko)
	hfields(i,je,ke,3,idblocke)=hfields(i,jo,ko,3,idblocko)
	hfields(i,je,ke,4,idblocke)=hfields(i,jo,ko,4,idblocko)
	hfields(i,je,ke,5,idblocke)=hfields(i,jo,ko,5,idblocko)
	hfields(i,je,ke,6,idblocke)=hfields(i,jo,ko,6,idblocko)
	hfields(i,je,ke,7,idblocke)=hfields(i,jo,ko,7,idblocko)
	hfields(i,je,ke,8,idblocke)=hfields(i,jo,ko,8,idblocko)
	hfields(i,je,ke,9,idblocke)=hfields(i,jo,ko,9,idblocko)
	hfields(i,je,ke,10,idblocke)=hfields(i,jo,ko,10,idblocko)
	
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
	
	hfields(i,je,ke,1,idblocke)=hfields(i,jo,ko,1,idblocko)
	hfields(i,je,ke,2,idblocke)=hfields(i,jo,ko,2,idblocko)
	hfields(i,je,ke,3,idblocke)=hfields(i,jo,ko,3,idblocko)
	hfields(i,je,ke,4,idblocke)=hfields(i,jo,ko,4,idblocko)
	hfields(i,je,ke,5,idblocke)=hfields(i,jo,ko,5,idblocko)
	hfields(i,je,ke,6,idblocke)=hfields(i,jo,ko,6,idblocko)
	hfields(i,je,ke,7,idblocke)=hfields(i,jo,ko,7,idblocko)
	hfields(i,je,ke,8,idblocke)=hfields(i,jo,ko,8,idblocko)
	hfields(i,je,ke,9,idblocke)=hfields(i,jo,ko,9,idblocko)
	hfields(i,je,ke,10,idblocke)=hfields(i,jo,ko,10,idblocko)
	
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
	
	hfields(i,je,ke,1,idblocke)=hfields(i,jo,ko,1,idblocko)
	hfields(i,je,ke,2,idblocke)=hfields(i,jo,ko,2,idblocko)
	hfields(i,je,ke,3,idblocke)=hfields(i,jo,ko,3,idblocko)
	hfields(i,je,ke,4,idblocke)=hfields(i,jo,ko,4,idblocko)
	hfields(i,je,ke,5,idblocke)=hfields(i,jo,ko,5,idblocko)
	hfields(i,je,ke,6,idblocke)=hfields(i,jo,ko,6,idblocko)
	hfields(i,je,ke,7,idblocke)=hfields(i,jo,ko,7,idblocko)
	hfields(i,je,ke,8,idblocke)=hfields(i,jo,ko,8,idblocko)
	hfields(i,je,ke,9,idblocke)=hfields(i,jo,ko,9,idblocko)
	hfields(i,je,ke,10,idblocke)=hfields(i,jo,ko,10,idblocko)
	
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
	
	hfields(i,je,ke,1,idblocke)=hfields(i,jo,ko,1,idblocko)
	hfields(i,je,ke,2,idblocke)=hfields(i,jo,ko,2,idblocko)
	hfields(i,je,ke,3,idblocke)=hfields(i,jo,ko,3,idblocko)
	hfields(i,je,ke,4,idblocke)=hfields(i,jo,ko,4,idblocko)
	hfields(i,je,ke,5,idblocke)=hfields(i,jo,ko,5,idblocko)
	hfields(i,je,ke,6,idblocke)=hfields(i,jo,ko,6,idblocko)
	hfields(i,je,ke,7,idblocke)=hfields(i,jo,ko,7,idblocko)
	hfields(i,je,ke,8,idblocke)=hfields(i,jo,ko,8,idblocko)
	hfields(i,je,ke,9,idblocke)=hfields(i,jo,ko,9,idblocko)
	hfields(i,je,ke,10,idblocke)=hfields(i,jo,ko,10,idblocko)
	
	
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
	
	hfieldsh(i,je,ke,1,idblocke)=hfieldsh(i,jo,ko,1,idblocko)
	hfieldsh(i,je,ke,2,idblocke)=hfieldsh(i,jo,ko,2,idblocko)
	hfieldsh(i,je,ke,3,idblocke)=hfieldsh(i,jo,ko,3,idblocko)
	hfieldsh(i,je,ke,4,idblocke)=hfieldsh(i,jo,ko,4,idblocko)
	hfieldsh(i,je,ke,5,idblocke)=hfieldsh(i,jo,ko,5,idblocko)
	hfieldsh(i,je,ke,6,idblocke)=hfieldsh(i,jo,ko,6,idblocko)
	hfieldsh(i,je,ke,7,idblocke)=hfieldsh(i,jo,ko,7,idblocko)
	hfieldsh(i,je,ke,8,idblocke)=hfieldsh(i,jo,ko,8,idblocko)
	hfieldsh(i,je,ke,9,idblocke)=hfieldsh(i,jo,ko,9,idblocko)
	hfieldsh(i,je,ke,10,idblocke)=hfieldsh(i,jo,ko,10,idblocko)
	
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
	
	hfieldsh(i,je,ke,1,idblocke)=hfieldsh(i,jo,ko,1,idblocko)
	hfieldsh(i,je,ke,2,idblocke)=hfieldsh(i,jo,ko,2,idblocko)
	hfieldsh(i,je,ke,3,idblocke)=hfieldsh(i,jo,ko,3,idblocko)
	hfieldsh(i,je,ke,4,idblocke)=hfieldsh(i,jo,ko,4,idblocko)
	hfieldsh(i,je,ke,5,idblocke)=hfieldsh(i,jo,ko,5,idblocko)
	hfieldsh(i,je,ke,6,idblocke)=hfieldsh(i,jo,ko,6,idblocko)
	hfieldsh(i,je,ke,7,idblocke)=hfieldsh(i,jo,ko,7,idblocko)
	hfieldsh(i,je,ke,8,idblocke)=hfieldsh(i,jo,ko,8,idblocko)
	hfieldsh(i,je,ke,9,idblocke)=hfieldsh(i,jo,ko,9,idblocko)
	hfieldsh(i,je,ke,10,idblocke)=hfieldsh(i,jo,ko,10,idblocko)
	
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
	
	hfieldsh(i,je,ke,1,idblocke)=hfieldsh(i,jo,ko,1,idblocko)
	hfieldsh(i,je,ke,2,idblocke)=hfieldsh(i,jo,ko,2,idblocko)
	hfieldsh(i,je,ke,3,idblocke)=hfieldsh(i,jo,ko,3,idblocko)
	hfieldsh(i,je,ke,4,idblocke)=hfieldsh(i,jo,ko,4,idblocko)
	hfieldsh(i,je,ke,5,idblocke)=hfieldsh(i,jo,ko,5,idblocko)
	hfieldsh(i,je,ke,6,idblocke)=hfieldsh(i,jo,ko,6,idblocko)
	hfieldsh(i,je,ke,7,idblocke)=hfieldsh(i,jo,ko,7,idblocko)
	hfieldsh(i,je,ke,8,idblocke)=hfieldsh(i,jo,ko,8,idblocko)
	hfieldsh(i,je,ke,9,idblocke)=hfieldsh(i,jo,ko,9,idblocko)
	hfieldsh(i,je,ke,10,idblocke)=hfieldsh(i,jo,ko,10,idblocko)
	
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
	
	hfieldsh(i,je,ke,1,idblocke)=hfieldsh(i,jo,ko,1,idblocko)
	hfieldsh(i,je,ke,2,idblocke)=hfieldsh(i,jo,ko,2,idblocko)
	hfieldsh(i,je,ke,3,idblocke)=hfieldsh(i,jo,ko,3,idblocko)
	hfieldsh(i,je,ke,4,idblocke)=hfieldsh(i,jo,ko,4,idblocko)
	hfieldsh(i,je,ke,5,idblocke)=hfieldsh(i,jo,ko,5,idblocko)
	hfieldsh(i,je,ke,6,idblocke)=hfieldsh(i,jo,ko,6,idblocko)
	hfieldsh(i,je,ke,7,idblocke)=hfieldsh(i,jo,ko,7,idblocko)
	hfieldsh(i,je,ke,8,idblocke)=hfieldsh(i,jo,ko,8,idblocko)
	hfieldsh(i,je,ke,9,idblocke)=hfieldsh(i,jo,ko,9,idblocko)
	hfieldsh(i,je,ke,10,idblocke)=hfieldsh(i,jo,ko,10,idblocko)
	
	
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
	
	hfields(i,je,k,1,idblocke)=hfields(i,jo,k,1,idblocko)
	hfields(i,je,k,2,idblocke)=hfields(i,jo,k,2,idblocko)
	hfields(i,je,k,3,idblocke)=hfields(i,jo,k,3,idblocko)
	hfields(i,je,k,4,idblocke)=hfields(i,jo,k,4,idblocko)
	hfields(i,je,k,5,idblocke)=hfields(i,jo,k,5,idblocko)
	hfields(i,je,k,6,idblocke)=hfields(i,jo,k,6,idblocko)
	hfields(i,je,k,7,idblocke)=hfields(i,jo,k,7,idblocko)
	hfields(i,je,k,8,idblocke)=hfields(i,jo,k,8,idblocko)
	hfields(i,je,k,9,idblocke)=hfields(i,jo,k,9,idblocko)
	hfields(i,je,k,10,idblocke)=hfields(i,jo,k,10,idblocko)
	
	gjo=1
	gje=ny_d+1
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=xblock+yblocko*nxblock_d+zblock*nxyblock_d+1
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=xblock+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	
	hfields(i,je,k,1,idblocke)=hfields(i,jo,k,1,idblocko)
	hfields(i,je,k,2,idblocke)=hfields(i,jo,k,2,idblocko)
	hfields(i,je,k,3,idblocke)=hfields(i,jo,k,3,idblocko)
	hfields(i,je,k,4,idblocke)=hfields(i,jo,k,4,idblocko)
	hfields(i,je,k,5,idblocke)=hfields(i,jo,k,5,idblocko)
	hfields(i,je,k,6,idblocke)=hfields(i,jo,k,6,idblocko)
	hfields(i,je,k,7,idblocke)=hfields(i,jo,k,7,idblocko)
	hfields(i,je,k,8,idblocke)=hfields(i,jo,k,8,idblocko)
	hfields(i,je,k,9,idblocke)=hfields(i,jo,k,9,idblocko)
	hfields(i,je,k,10,idblocke)=hfields(i,jo,k,10,idblocko)
	  
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
	
	hfieldsh(i,je,k,1,idblocke)=hfieldsh(i,jo,k,1,idblocko)
	hfieldsh(i,je,k,2,idblocke)=hfieldsh(i,jo,k,2,idblocko)
	hfieldsh(i,je,k,3,idblocke)=hfieldsh(i,jo,k,3,idblocko)
	hfieldsh(i,je,k,4,idblocke)=hfieldsh(i,jo,k,4,idblocko)
	hfieldsh(i,je,k,5,idblocke)=hfieldsh(i,jo,k,5,idblocko)
	hfieldsh(i,je,k,6,idblocke)=hfieldsh(i,jo,k,6,idblocko)
	hfieldsh(i,je,k,7,idblocke)=hfieldsh(i,jo,k,7,idblocko)
	hfieldsh(i,je,k,8,idblocke)=hfieldsh(i,jo,k,8,idblocko)
	hfieldsh(i,je,k,9,idblocke)=hfieldsh(i,jo,k,9,idblocko)
	hfieldsh(i,je,k,10,idblocke)=hfieldsh(i,jo,k,10,idblocko)
	
	gjo=1
	gje=ny_d+1
	
	yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=xblock+yblocko*nxblock_d+zblock*nxyblock_d+1
    
    yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=xblock+yblocke*nxblock_d+zblock*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	
	hfieldsh(i,je,k,1,idblocke)=hfieldsh(i,jo,k,1,idblocko)
	hfieldsh(i,je,k,2,idblocke)=hfieldsh(i,jo,k,2,idblocko)
	hfieldsh(i,je,k,3,idblocke)=hfieldsh(i,jo,k,3,idblocko)
	hfieldsh(i,je,k,4,idblocke)=hfieldsh(i,jo,k,4,idblocko)
	hfieldsh(i,je,k,5,idblocke)=hfieldsh(i,jo,k,5,idblocko)
	hfieldsh(i,je,k,6,idblocke)=hfieldsh(i,jo,k,6,idblocko)
	hfieldsh(i,je,k,7,idblocke)=hfieldsh(i,jo,k,7,idblocko)
	hfieldsh(i,je,k,8,idblocke)=hfieldsh(i,jo,k,8,idblocko)
	hfieldsh(i,je,k,9,idblocke)=hfieldsh(i,jo,k,9,idblocko)
	hfieldsh(i,je,k,10,idblocke)=hfieldsh(i,jo,k,10,idblocko)
	  
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
	
	hfields(ie,je,k,1,idblocke)=hfields(io,jo,k,1,idblocko)
	hfields(ie,je,k,2,idblocke)=hfields(io,jo,k,2,idblocko)
	hfields(ie,je,k,3,idblocke)=hfields(io,jo,k,3,idblocko)
	hfields(ie,je,k,4,idblocke)=hfields(io,jo,k,4,idblocko)
	hfields(ie,je,k,5,idblocke)=hfields(io,jo,k,5,idblocko)
	hfields(ie,je,k,6,idblocke)=hfields(io,jo,k,6,idblocko)
	hfields(ie,je,k,7,idblocke)=hfields(io,jo,k,7,idblocko)
	hfields(ie,je,k,8,idblocke)=hfields(io,jo,k,8,idblocko)
	hfields(ie,je,k,9,idblocke)=hfields(io,jo,k,9,idblocko)
	hfields(ie,je,k,10,idblocke)=hfields(io,jo,k,10,idblocko)
	
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
	
	hfields(ie,je,k,1,idblocke)=hfields(io,jo,k,1,idblocko)
	hfields(ie,je,k,2,idblocke)=hfields(io,jo,k,2,idblocko)
	hfields(ie,je,k,3,idblocke)=hfields(io,jo,k,3,idblocko)
	hfields(ie,je,k,4,idblocke)=hfields(io,jo,k,4,idblocko)
	hfields(ie,je,k,5,idblocke)=hfields(io,jo,k,5,idblocko)
	hfields(ie,je,k,6,idblocke)=hfields(io,jo,k,6,idblocko)
	hfields(ie,je,k,7,idblocke)=hfields(io,jo,k,7,idblocko)
	hfields(ie,je,k,8,idblocke)=hfields(io,jo,k,8,idblocko)
	hfields(ie,je,k,9,idblocke)=hfields(io,jo,k,9,idblocko)
	hfields(ie,je,k,10,idblocke)=hfields(io,jo,k,10,idblocko)
	
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
	
	hfields(ie,je,k,1,idblocke)=hfields(io,jo,k,1,idblocko)
	hfields(ie,je,k,2,idblocke)=hfields(io,jo,k,2,idblocko)
	hfields(ie,je,k,3,idblocke)=hfields(io,jo,k,3,idblocko)
	hfields(ie,je,k,4,idblocke)=hfields(io,jo,k,4,idblocko)
	hfields(ie,je,k,5,idblocke)=hfields(io,jo,k,5,idblocko)
	hfields(ie,je,k,6,idblocke)=hfields(io,jo,k,6,idblocko)
	hfields(ie,je,k,7,idblocke)=hfields(io,jo,k,7,idblocko)
	hfields(ie,je,k,8,idblocke)=hfields(io,jo,k,8,idblocko)
	hfields(ie,je,k,9,idblocke)=hfields(io,jo,k,9,idblocko)
	hfields(ie,je,k,10,idblocke)=hfields(io,jo,k,10,idblocko)
	
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
	
	hfields(ie,je,k,1,idblocke)=hfields(io,jo,k,1,idblocko)
	hfields(ie,je,k,2,idblocke)=hfields(io,jo,k,2,idblocko)
	hfields(ie,je,k,3,idblocke)=hfields(io,jo,k,3,idblocko)
	hfields(ie,je,k,4,idblocke)=hfields(io,jo,k,4,idblocko)
	hfields(ie,je,k,5,idblocke)=hfields(io,jo,k,5,idblocko)
	hfields(ie,je,k,6,idblocke)=hfields(io,jo,k,6,idblocko)
	hfields(ie,je,k,7,idblocke)=hfields(io,jo,k,7,idblocko)
	hfields(ie,je,k,8,idblocke)=hfields(io,jo,k,8,idblocko)
	hfields(ie,je,k,9,idblocke)=hfields(io,jo,k,9,idblocko)
	hfields(ie,je,k,10,idblocke)=hfields(io,jo,k,10,idblocko)
	
	
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
	
	hfieldsh(ie,je,k,1,idblocke)=hfieldsh(io,jo,k,1,idblocko)
	hfieldsh(ie,je,k,2,idblocke)=hfieldsh(io,jo,k,2,idblocko)
	hfieldsh(ie,je,k,3,idblocke)=hfieldsh(io,jo,k,3,idblocko)
	hfieldsh(ie,je,k,4,idblocke)=hfieldsh(io,jo,k,4,idblocko)
	hfieldsh(ie,je,k,5,idblocke)=hfieldsh(io,jo,k,5,idblocko)
	hfieldsh(ie,je,k,6,idblocke)=hfieldsh(io,jo,k,6,idblocko)
	hfieldsh(ie,je,k,7,idblocke)=hfieldsh(io,jo,k,7,idblocko)
	hfieldsh(ie,je,k,8,idblocke)=hfieldsh(io,jo,k,8,idblocko)
	hfieldsh(ie,je,k,9,idblocke)=hfieldsh(io,jo,k,9,idblocko)
	hfieldsh(ie,je,k,10,idblocke)=hfieldsh(io,jo,k,10,idblocko)
	
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
	
	hfieldsh(ie,je,k,1,idblocke)=hfieldsh(io,jo,k,1,idblocko)
	hfieldsh(ie,je,k,2,idblocke)=hfieldsh(io,jo,k,2,idblocko)
	hfieldsh(ie,je,k,3,idblocke)=hfieldsh(io,jo,k,3,idblocko)
	hfieldsh(ie,je,k,4,idblocke)=hfieldsh(io,jo,k,4,idblocko)
	hfieldsh(ie,je,k,5,idblocke)=hfieldsh(io,jo,k,5,idblocko)
	hfieldsh(ie,je,k,6,idblocke)=hfieldsh(io,jo,k,6,idblocko)
	hfieldsh(ie,je,k,7,idblocke)=hfieldsh(io,jo,k,7,idblocko)
	hfieldsh(ie,je,k,8,idblocke)=hfieldsh(io,jo,k,8,idblocko)
	hfieldsh(ie,je,k,9,idblocke)=hfieldsh(io,jo,k,9,idblocko)
	hfieldsh(ie,je,k,10,idblocke)=hfieldsh(io,jo,k,10,idblocko)
	
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
	
	hfieldsh(ie,je,k,1,idblocke)=hfieldsh(io,jo,k,1,idblocko)
	hfieldsh(ie,je,k,2,idblocke)=hfieldsh(io,jo,k,2,idblocko)
	hfieldsh(ie,je,k,3,idblocke)=hfieldsh(io,jo,k,3,idblocko)
	hfieldsh(ie,je,k,4,idblocke)=hfieldsh(io,jo,k,4,idblocko)
	hfieldsh(ie,je,k,5,idblocke)=hfieldsh(io,jo,k,5,idblocko)
	hfieldsh(ie,je,k,6,idblocke)=hfieldsh(io,jo,k,6,idblocko)
	hfieldsh(ie,je,k,7,idblocke)=hfieldsh(io,jo,k,7,idblocko)
	hfieldsh(ie,je,k,8,idblocke)=hfieldsh(io,jo,k,8,idblocko)
	hfieldsh(ie,je,k,9,idblocke)=hfieldsh(io,jo,k,9,idblocko)
	hfieldsh(ie,je,k,10,idblocke)=hfieldsh(io,jo,k,10,idblocko)
	
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
	
	hfieldsh(ie,je,k,1,idblocke)=hfieldsh(io,jo,k,1,idblocko)
	hfieldsh(ie,je,k,2,idblocke)=hfieldsh(io,jo,k,2,idblocko)
	hfieldsh(ie,je,k,3,idblocke)=hfieldsh(io,jo,k,3,idblocko)
	hfieldsh(ie,je,k,4,idblocke)=hfieldsh(io,jo,k,4,idblocko)
	hfieldsh(ie,je,k,5,idblocke)=hfieldsh(io,jo,k,5,idblocko)
	hfieldsh(ie,je,k,6,idblocke)=hfieldsh(io,jo,k,6,idblocko)
	hfieldsh(ie,je,k,7,idblocke)=hfieldsh(io,jo,k,7,idblocko)
	hfieldsh(ie,je,k,8,idblocke)=hfieldsh(io,jo,k,8,idblocko)
	hfieldsh(ie,je,k,9,idblocke)=hfieldsh(io,jo,k,9,idblocko)
	hfieldsh(ie,je,k,10,idblocke)=hfieldsh(io,jo,k,10,idblocko)
	
	
    return
		  
  end subroutine pbc_edge_z_flop
 
 end module pbc_kernels
