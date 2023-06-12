#include "defines.h"
 module cudavars
   
    use cudafor
    
    implicit none
    
    integer, parameter :: db=kind(1.e0)
    
    real(kind=db),parameter :: zero=real(0.d0,kind=db)
    real(kind=db),parameter :: one=real(1.d0,kind=db)
    real(kind=db),parameter :: half=real(0.5d0,kind=db)
    real(kind=db),parameter :: halfthree=real(1.5d0,kind=db)
    real(kind=db),parameter :: two=real(2.d0,kind=db)
    real(kind=db),parameter :: three=real(3.d0,kind=db)
    real(kind=db),parameter :: five=real(5.d0,kind=db)
    real(kind=db),parameter :: ten=real(10.d0,kind=db)
    real(kind=db),parameter :: eighteen=real(18.d0,kind=db)
    real(kind=db),parameter :: thirtysix=real(36.d0,kind=db)
    
    real(kind=db),parameter :: pi_greek=real(3.141592653589793238462643383279502884d0,kind=db)
    real(kind=db),parameter :: p0 = (one/three)
    real(kind=db),parameter :: p1 = (one/eighteen)
    real(kind=db),parameter :: p2 = (one/thirtysix)
    real(kind=db),parameter :: cssq = one/three
    real(kind=db),parameter :: onecssq = three
    real(kind=db),parameter :: halfonecssq = halfthree
    real(kind=db),parameter :: p1dcssq=p1/cssq
    real(kind=db),parameter :: p2dcssq=p2/cssq
    
    real(kind=db),parameter :: pi2cssq0=p0/(two*cssq**two)
    real(kind=db),parameter :: pi2cssq1=p1/(two*cssq**two)
    real(kind=db),parameter :: pi2cssq2=p2/(two*cssq**two)

    real(kind=db),parameter :: qxx=one-cssq
    real(kind=db),parameter :: qyy=one-cssq
    real(kind=db),parameter :: qzz=one-cssq
    real(kind=db),parameter :: qxy_7_8=one
    real(kind=db),parameter :: qxy_9_10=-one
    real(kind=db),parameter :: qxz_15_16=one
    real(kind=db),parameter :: qxz_17_18=-one
    real(kind=db),parameter :: qyz_11_12=one
    real(kind=db),parameter :: qyz_13_14=-one
    
    integer, constant :: nx_d,ny_d,nz_d
    
    real(kind=db), constant :: omega,fx,fy,fz,oneminusomega
    integer(kind=cuda_Stream_Kind) :: stream1,stream2
    type (cudaDeviceProp) :: prop
    type (cudaEvent) :: startEvent, stopEvent, dummyEvent, dummyEvent1, dummyEvent2
    type (dim3) :: dimGrid,dimBlock,dimGridx,dimGridy,dimBlock2, &
     dimGridhalo,dimBlockhalo,dimGridshared,dimBlockshared
    
    integer, constant :: TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
    
    integer, constant :: nxblock_d,nxyblock_d,nblocks_d
    
    !real(kind=db), allocatable, dimension(:,:,:,:), device :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    !real(kind=db), allocatable, dimension(:,:,:,:), device :: rhoh,uh,vh,wh,pxxh,pxyh,pxzh,pyyh,pyzh,pzzh
    real(kind=db), allocatable, dimension(:,:,:,:,:), device :: hfields,hfieldsh
    real(kind=4), allocatable, dimension(:,:,:), device :: rhoprint_d
    real(kind=4), allocatable, dimension(:,:,:,:), device :: velprint_d
    integer(kind=1), allocatable, dimension(:,:,:), device   :: isfluid
    
    contains
    
    attributes(global) subroutine setup_system(rhos,vxs,vys,vzs)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    real :: mytest
    
    !integer :: i,j,k,ii,jj,kk,xblock,yblock,zblock,idblock
    integer :: i,j,k,gi,gj,gk,idblock       
            
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	!xblock=(i+TILE_DIMx_d-1)/TILE_DIMx_d
    !yblock=(j+TILE_DIMy_d-1)/TILE_DIMy_d
    !zblock=(k+TILE_DIMz_d-1)/TILE_DIMz_d
	!idblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
	!ii=i-xblock*TILE_DIMx_d+TILE_DIMx_d
    !jj=j-yblock*TILE_DIMy_d+TILE_DIMy_d
    !kk=k-zblock*TILE_DIMz_d+TILE_DIMz_d
    !if(ii/=threadIdx%x .or. jj/=threadIdx%y .or. kk/=threadIdx%z)write(*,*)'cazzo1'
    !if(xblock/=(blockIdx%x) .or. yblock/=(blockIdx%y) .or. zblock/=(blockIdx%z))write(*,*)'cazzo2'
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	
	mytest=rhos!real(gi**3+gj**3+gk**3)
	
	hfields(i,j,k,1,idblock)=mytest!rhos
	hfields(i,j,k,2,idblock)=vxs
	hfields(i,j,k,3,idblock)=vys
	hfields(i,j,k,4,idblock)=vzs
	hfields(i,j,k,5,idblock)=zero
	hfields(i,j,k,6,idblock)=zero
	hfields(i,j,k,7,idblock)=zero
	hfields(i,j,k,8,idblock)=zero
	hfields(i,j,k,9,idblock)=zero
	hfields(i,j,k,10,idblock)=zero
	
	hfieldsh(i,j,k,1,idblock)=mytest!rhos
	hfieldsh(i,j,k,2,idblock)=vxs
	hfieldsh(i,j,k,3,idblock)=vys
	hfieldsh(i,j,k,4,idblock)=vzs
	hfieldsh(i,j,k,5,idblock)=zero
	hfieldsh(i,j,k,6,idblock)=zero
	hfieldsh(i,j,k,7,idblock)=zero
	hfieldsh(i,j,k,8,idblock)=zero
	hfieldsh(i,j,k,9,idblock)=zero
	hfieldsh(i,j,k,10,idblock)=zero
	
    
    return

 end subroutine setup_system
 
 attributes(global) subroutine setup_system_halo(rhos,vxs,vys,vzs)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    real :: mytest
    
    !integer :: i,j,k,ii,jj,kk,xblock,yblock,zblock,idblock
    integer :: i,j,k,gi,gj,gk,idblock       
            
	gi = (blockIdx%x-2) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-2) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-2) * TILE_DIMz_d + threadIdx%z
	
	
	if(gi>nx_d .or. gj>ny_d .or. gk>nz_d)return
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	!xblock=(i+TILE_DIMx_d-1)/TILE_DIMx_d
    !yblock=(j+TILE_DIMy_d-1)/TILE_DIMy_d
    !zblock=(k+TILE_DIMz_d-1)/TILE_DIMz_d
	!idblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
	!ii=i-xblock*TILE_DIMx_d+TILE_DIMx_d
    !jj=j-yblock*TILE_DIMy_d+TILE_DIMy_d
    !kk=k-zblock*TILE_DIMz_d+TILE_DIMz_d
    !if(ii/=threadIdx%x .or. jj/=threadIdx%y .or. kk/=threadIdx%z)write(*,*)'cazzo1'
    !if(xblock/=(blockIdx%x) .or. yblock/=(blockIdx%y) .or. zblock/=(blockIdx%z))write(*,*)'cazzo2'
	idblock=(blockIdx%x-1)+(blockIdx%y-1)*nxblock_d+(blockIdx%z-1)*nxyblock_d+1
	
	mytest=rhos!real(gi**3+gj**3+gk**3)
	
	
	hfields(i,j,k,1,idblock)=mytest!rhos
	hfields(i,j,k,2,idblock)=vxs
	hfields(i,j,k,3,idblock)=vys
	hfields(i,j,k,4,idblock)=vzs
	hfields(i,j,k,5,idblock)=zero
	hfields(i,j,k,6,idblock)=zero
	hfields(i,j,k,7,idblock)=zero
	hfields(i,j,k,8,idblock)=zero
	hfields(i,j,k,9,idblock)=zero
	hfields(i,j,k,10,idblock)=zero
	
	!if(gi==3 .and. gj==3 .and. gk==0)write(*,*)'CAZZONE ',gi,gj,gk
	
	hfieldsh(i,j,k,1,idblock)=mytest!rhos
	hfieldsh(i,j,k,2,idblock)=vxs
	hfieldsh(i,j,k,3,idblock)=vys
	hfieldsh(i,j,k,4,idblock)=vzs
	hfieldsh(i,j,k,5,idblock)=zero
	hfieldsh(i,j,k,6,idblock)=zero
	hfieldsh(i,j,k,7,idblock)=zero
	hfieldsh(i,j,k,8,idblock)=zero
	hfieldsh(i,j,k,9,idblock)=zero
	hfieldsh(i,j,k,10,idblock)=zero
    
    return

 end subroutine setup_system_halo
 
 attributes(global) subroutine store_print()
	
	integer :: i,j,k,gi,gj,gk,idblock,ii,jj,kk  ,xblock,yblock,zblock,idblocko
    real :: mytest
    
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	
	xblock=(gi+TILE_DIMx_d-1)/TILE_DIMx_d
	yblock=(gj+TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
    ii=gi-xblock*TILE_DIMx_d+TILE_DIMx_d
    jj=gj-yblock*TILE_DIMy_d+TILE_DIMy_d
    kk=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
    
    idblocko=xblock+yblock*nxblock_d+zblock*nxyblock_d-1
	
	mytest=real(gi**3+gj**3+gk**3)
	
	!if(rho(ii,jj,kk,idblocko).ne. mytest)write(*,*)'cazzo'
	  
	!write(*,*)i,j,p_d(0)*myrho_d
	if(abs(isfluid(gi,gj,gk)).eq.1)then
      rhoprint_d(gi,gj,gk)=hfields(i,j,k,1,idblock)
	  velprint_d(1,gi,gj,gk)=hfields(i,j,k,2,idblock)
	  velprint_d(2,gi,gj,gk)=hfields(i,j,k,3,idblock)
	  velprint_d(3,gi,gj,gk)=hfields(i,j,k,4,idblock)
		
	else
	  
	  rhoprint_d(gi,gj,gk)=zero
	  velprint_d(1,gi,gj,gk)=zero
	  velprint_d(2,gi,gj,gk)=zero
	  velprint_d(3,gi,gj,gk)=zero
	  
	endif
	  
	  return

 end subroutine store_print
 
 attributes(global) subroutine store_print_flop()
	  
	integer :: i,j,k,gi,gj,gk,idblock       
            
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	  
	  
	!write(*,*)i,j,p_d(0)*myrho_d
	if(abs(isfluid(gi,gj,gk)).eq.1)then
      rhoprint_d(gi,gj,gk)=hfieldsh(i,j,k,1,idblock)
	  velprint_d(1,gi,gj,gk)=hfieldsh(i,j,k,2,idblock)
	  velprint_d(2,gi,gj,gk)=hfieldsh(i,j,k,3,idblock)
	  velprint_d(3,gi,gj,gk)=hfieldsh(i,j,k,4,idblock)
		
	else
	  
	  rhoprint_d(gi,gj,gk)=zero
	  velprint_d(1,gi,gj,gk)=zero
	  velprint_d(2,gi,gj,gk)=zero
	  velprint_d(3,gi,gj,gk)=zero
	  
	endif
	  
	return

 end subroutine store_print_flop
  
 subroutine abortOnLastErrorAndSync(msg, step)
    implicit none
    integer, intent(in) :: step
    character(len=*), intent(in) :: msg
    integer :: istat0
    
    istat0 = cudaGetLastError()
    
    if (istat0/=0) then
      write(*,*) 'status after ',msg,':', cudaGetErrorString(istat0)
      write(*,*) 'Exiting at step:', step
      stop
    endif
    
    return
    
 end subroutine  abortOnLastErrorAndSync

 end module cudavars
