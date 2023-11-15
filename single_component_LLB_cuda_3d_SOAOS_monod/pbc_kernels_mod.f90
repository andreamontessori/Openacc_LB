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
	
	yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    
    j=gj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
	
	gio=nx_d
	gie=0
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=(xblocko-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=(xblocke-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
	  
	hfields(idx5d(ie,j,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gio=1
	gie=nx_d+1
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=(xblocko-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=(xblocke-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
	
	hfields(idx5d(ie,j,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,j,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,10,idblocko,nx_d,ny_d,nz_d,10))
	  
    return
		  
  end subroutine pbc_side_x
  
  attributes(global) subroutine pbc_side_x_flop()
	  
	integer :: io,ie,j,k,gio,gie,gj,gk,xblocko,xblocke,yblock,zblock,idblocko,idblocke
	
	
	gj = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	  
	if (gj>ny_d .or. gk>nz_d)return
	
	yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    
    j=gj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
	
	gio=nx_d
	gie=0
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=(xblocko-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=(xblocke-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    
	  
	hfieldsh(idx5d(ie,j,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gio=1
	gie=nx_d+1
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
    idblocko=(xblocko-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
	idblocke=(xblocke-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
	ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
	
	hfieldsh(idx5d(ie,j,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,j,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,10,idblocko,nx_d,ny_d,nz_d,10))
	  
    return
		  
  end subroutine pbc_side_x_flop
  
  attributes(global) subroutine pbc_edge_x()
	  
    integer :: i,jo,je,ko,ke,gi,gjo,gje,gko,gke,xblock,yblocko,yblocke,zblocko,zblocke,idblocko,idblocke
	
	
	gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (gi>nx_d)return
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
	i=gi-xblock*TILE_DIMx_d+2*TILE_DIMx_d
	
	gjo=ny_d
	gje=0
	gko=nz_d
	gke=0
	
    yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+2*TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
    ko=gko-zblocko*TILE_DIMz_d+2*TILE_DIMz_d
	
	yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+2*TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
    ke=gke-zblocke*TILE_DIMz_d+2*TILE_DIMz_d
    
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblocko-1)*nxyblock_d+1
    idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblocke-1)*nxyblock_d+1
	
	hfields(idx5d(i,je,ke,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gjo=1
	gje=ny_d+1
	gko=1
	gke=nz_d+1
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+2*TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
    ko=gko-zblocko*TILE_DIMz_d+2*TILE_DIMz_d
	
	yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+2*TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
    ke=gke-zblocke*TILE_DIMz_d+2*TILE_DIMz_d
    
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblocko-1)*nxyblock_d+1
    idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblocke-1)*nxyblock_d+1
	
	hfields(idx5d(i,je,ke,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gjo=ny_d
	gje=0
	gko=1
	gke=nz_d+1
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+2*TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
    ko=gko-zblocko*TILE_DIMz_d+2*TILE_DIMz_d
	
	yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+2*TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
    ke=gke-zblocke*TILE_DIMz_d+2*TILE_DIMz_d
    
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblocko-1)*nxyblock_d+1
    idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblocke-1)*nxyblock_d+1
	
	hfields(idx5d(i,je,ke,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gjo=1
	gje=ny_d+1
	gko=nz_d
	gke=0
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+2*TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
    ko=gko-zblocko*TILE_DIMz_d+2*TILE_DIMz_d
	
	yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+2*TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
    ke=gke-zblocke*TILE_DIMz_d+2*TILE_DIMz_d
    
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblocko-1)*nxyblock_d+1
    idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblocke-1)*nxyblock_d+1
	
	hfields(idx5d(i,je,ke,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,ke,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,ko,10,idblocko,nx_d,ny_d,nz_d,10))
	
	
    return
		  
  end subroutine pbc_edge_x
  
  attributes(global) subroutine pbc_edge_x_flop()
	  
    integer :: i,jo,je,ko,ke,gi,gjo,gje,gko,gke,xblock,yblocko,yblocke,zblocko,zblocke,idblocko,idblocke
	
	
	gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (gi>nx_d)return
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
	i=gi-xblock*TILE_DIMx_d+2*TILE_DIMx_d
	
	gjo=ny_d
	gje=0
	gko=nz_d
	gke=0
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+2*TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
    ko=gko-zblocko*TILE_DIMz_d+2*TILE_DIMz_d
	
	yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+2*TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
    ke=gke-zblocke*TILE_DIMz_d+2*TILE_DIMz_d
    
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblocko-1)*nxyblock_d+1
    idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblocke-1)*nxyblock_d+1
	
	hfieldsh(idx5d(i,je,ke,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gjo=1
	gje=ny_d+1
	gko=1
	gke=nz_d+1
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+2*TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
    ko=gko-zblocko*TILE_DIMz_d+2*TILE_DIMz_d
	
	yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+2*TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
    ke=gke-zblocke*TILE_DIMz_d+2*TILE_DIMz_d
    
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblocko-1)*nxyblock_d+1
    idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblocke-1)*nxyblock_d+1
	
	hfieldsh(idx5d(i,je,ke,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gjo=ny_d
	gje=0
	gko=1
	gke=nz_d+1
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+2*TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
    ko=gko-zblocko*TILE_DIMz_d+2*TILE_DIMz_d
	
	yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+2*TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
    ke=gke-zblocke*TILE_DIMz_d+2*TILE_DIMz_d
    
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblocko-1)*nxyblock_d+1
    idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblocke-1)*nxyblock_d+1
	
	hfieldsh(idx5d(i,je,ke,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gjo=1
	gje=ny_d+1
	gko=nz_d
	gke=0
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocko=(gko+2*TILE_DIMz_d-1)/TILE_DIMz_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
    ko=gko-zblocko*TILE_DIMz_d+2*TILE_DIMz_d
	
	yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblocke=(gke+2*TILE_DIMz_d-1)/TILE_DIMz_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
    ke=gke-zblocke*TILE_DIMz_d+2*TILE_DIMz_d
    
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblocko-1)*nxyblock_d+1
    idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblocke-1)*nxyblock_d+1
	
	hfieldsh(idx5d(i,je,ke,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,ke,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,ko,10,idblocko,nx_d,ny_d,nz_d,10))
	
	
    return
		  
  end subroutine pbc_edge_x_flop

  attributes(global) subroutine pbc_side_y()
	
	implicit none
	integer :: i,jo,je,k,gi,gjo,gje,gk,xblock,yblocko,yblocke,zblock,idblocko,idblocke
	  
	gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

	if (gi>nx_d .or. gk>nz_d)return
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    i=gi-xblock*TILE_DIMx_d+2*TILE_DIMx_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
    
    gjo=ny_d
	gje=0
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
	
	hfields(idx5d(i,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gjo=1
	gje=ny_d+1
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
	
	hfields(idx5d(i,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(i,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(i,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	  
    return

  end subroutine pbc_side_y
  
  attributes(global) subroutine pbc_side_y_flop()
	
	implicit none
	
	integer :: i,jo,je,k,gi,gjo,gje,gk,xblock,yblocko,yblocke,zblock,idblocko,idblocke
	  
	gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

	if (gi>nx_d .or. gk>nz_d)return
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    i=gi-xblock*TILE_DIMx_d+2*TILE_DIMx_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
    
    gjo=ny_d
	gje=0
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
	
	hfieldsh(idx5d(i,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gjo=1
	gje=ny_d+1
	
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    idblocko=(xblock-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
    
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	idblocke=(xblock-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMy_d
	je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMy_d
	
	hfieldsh(idx5d(i,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(i,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(i,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	  
    return

  end subroutine pbc_side_y_flop
  
  attributes(global) subroutine pbc_edge_z()
	  
    integer :: io,ie,jo,je,k,gio,gie,gjo,gje,gk,xblocko,xblocke,yblocko,yblocke,zblock,idblocko,idblocke
	
	
	gk = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (gk>nz_d)return
	
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
	
	
	gio=nx_d
	gie=0
	gjo=ny_d
	gje=0
	
    xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMx_d
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	
    ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMx_d
	
    
	idblocko=(xblocko-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	idblocke=(xblocke-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	hfields(idx5d(ie,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gio=1
	gie=nx_d+1
	gjo=1
	gje=ny_d+1
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMx_d
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	
    ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMx_d
	
    
	idblocko=(xblocko-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	idblocke=(xblocke-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	hfields(idx5d(ie,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gio=1
	gie=nx_d+1
	gjo=ny_d
	gje=0
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMx_d
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	
    ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMx_d
	
    
	idblocko=(xblocko-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	idblocke=(xblocke-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	hfields(idx5d(ie,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gio=nx_d
	gie=0
	gjo=1
	gje=ny_d+1
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMx_d
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	
    ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMx_d
	
    
	idblocko=(xblocko-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	idblocke=(xblocke-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	hfields(idx5d(ie,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfields(idx5d(ie,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	
    return
		  
  end subroutine pbc_edge_z
  
  attributes(global) subroutine pbc_edge_z_flop()
	  
    integer :: io,ie,jo,je,k,gio,gie,gjo,gje,gk,xblocko,xblocke,yblocko,yblocke,zblock,idblocko,idblocke
	
	
	gk = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	  
	if (gk>nz_d)return
	
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
	
	gio=nx_d
	gie=0
	gjo=ny_d
	gje=0
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMx_d
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	
    ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMx_d
	
    
	idblocko=(xblocko-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	idblocke=(xblocke-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	hfieldsh(idx5d(ie,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gio=1
	gie=nx_d+1
	gjo=1
	gje=ny_d+1
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMx_d
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	
    ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMx_d
	
    
	idblocko=(xblocko-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	idblocke=(xblocke-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	hfieldsh(idx5d(ie,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gio=1
	gie=nx_d+1
	gjo=ny_d
	gje=0
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMx_d
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	
    ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMx_d
	
    
	idblocko=(xblocko-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	idblocke=(xblocke-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	hfieldsh(idx5d(ie,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	gio=nx_d
	gie=0
	gjo=1
	gje=ny_d+1
	
	xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblocko=(gjo+2*TILE_DIMy_d-1)/TILE_DIMy_d
    
    io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
    jo=gjo-yblocko*TILE_DIMy_d+2*TILE_DIMx_d
    
    xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblocke=(gje+2*TILE_DIMy_d-1)/TILE_DIMy_d
	
    ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
    je=gje-yblocke*TILE_DIMy_d+2*TILE_DIMx_d
	
    
	idblocko=(xblocko-1)+(yblocko-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	idblocke=(xblocke-1)+(yblocke-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	hfieldsh(idx5d(ie,je,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,1,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,2,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,3,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,4,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,5,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,6,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,7,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,8,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,9,idblocko,nx_d,ny_d,nz_d,10))
	hfieldsh(idx5d(ie,je,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,jo,k,10,idblocko,nx_d,ny_d,nz_d,10))
	
	
    return
		  
  end subroutine pbc_edge_z_flop
  
  attributes(global) subroutine bc_per_x_hvar(step)
  
      integer, value :: step
      
      integer :: j,k,gj,gk
      integer :: yblock,zblock
      integer :: xblocko,xblocke
      integer :: idblocko,idblocke
      integer :: io,ie,gio,gie
  
      gj = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
      gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
      
      if (gj>ny_d .or. gk>nz_d) return
      
      yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      j=gj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      gio=2
	  gie=nx_d
	  
	  xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
      idblocko=(xblocko-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	  
	  xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  idblocke=(xblocke-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	  
	  io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
	  ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
	  
	  hfields(idx5d(ie,j,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,1,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,2,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,3,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,4,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,5,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,6,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,7,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,8,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,9,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,10,idblocko,nx_d,ny_d,nz_d,10))
	  
	  gio=nx_d-1
	  gie=1
	  
	  xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
      idblocko=(xblocko-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	  
	  xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  idblocke=(xblocke-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	  
	  io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
	  ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
	  
      hfields(idx5d(ie,j,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,1,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,2,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,3,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,4,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,5,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,6,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,7,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,8,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,9,idblocko,nx_d,ny_d,nz_d,10))
	  hfields(idx5d(ie,j,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfields(idx5d(io,j,k,10,idblocko,nx_d,ny_d,nz_d,10))
      
      
    end subroutine bc_per_x_hvar
    
    attributes(global) subroutine bc_per_x_hvar_flop(step)
  
      integer, value :: step
      
      integer :: j,k,gj,gk
      integer :: yblock,zblock
      integer :: xblocko,xblocke
      integer :: idblocko,idblocke
      integer :: io,ie,gio,gie
  
      gj = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
      gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
      
      if (gj>ny_d .or. gk>nz_d) return
      
      yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      j=gj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      gio=2
	  gie=nx_d
	  
	  xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
      idblocko=(xblocko-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	  
	  xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  idblocke=(xblocke-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	  
	  io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
	  ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
	  
	  hfieldsh(idx5d(ie,j,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,1,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,2,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,3,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,4,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,5,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,6,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,7,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,8,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,9,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,10,idblocko,nx_d,ny_d,nz_d,10))
	  
	  gio=nx_d-1
	  gie=1
	  
	  xblocko=(gio+2*TILE_DIMx_d-1)/TILE_DIMx_d
      idblocko=(xblocko-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	  
	  xblocke=(gie+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  idblocke=(xblocke-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	  
	  io=gio-xblocko*TILE_DIMx_d+2*TILE_DIMx_d
	  ie=gie-xblocke*TILE_DIMx_d+2*TILE_DIMx_d
	  
      hfieldsh(idx5d(ie,j,k,1,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,1,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,2,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,2,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,3,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,3,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,4,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,4,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,5,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,5,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,6,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,6,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,7,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,7,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,8,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,8,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,9,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,9,idblocko,nx_d,ny_d,nz_d,10))
	  hfieldsh(idx5d(ie,j,k,10,idblocke,nx_d,ny_d,nz_d,10))=hfieldsh(idx5d(io,j,k,10,idblocko,nx_d,ny_d,nz_d,10))
      
      
    end subroutine bc_per_x_hvar_flop
 
 end module pbc_kernels
