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
    type (dim3) :: dimGrid,dimBlock,dimGridx,dimGridy,dimBlock2
    
    integer, constant :: TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
    
    integer, constant :: nxblock_d,nxyblock_d
    
    real(kind=db), allocatable, dimension(:,:,:,:), device :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:,:), device :: rhoh,uh,vh,wh,pxxh,pxyh,pxzh,pyyh,pyzh,pzzh
    real(kind=4), allocatable, dimension(:,:,:), device :: rhoprint_d
    real(kind=4), allocatable, dimension(:,:,:,:), device :: velprint_d
    integer(kind=1), allocatable, dimension(:,:,:), device   :: isfluid
    
    contains
    
    attributes(global) subroutine setup_system(rhos,vxs,vys,vzs)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    
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
	
	
	u(i,j,k,idblock)=vxs
	v(i,j,k,idblock)=vys
	w(i,j,k,idblock)=vzs
	rho(i,j,k,idblock)=rhos
	pxx(i,j,k,idblock)=zero
	pxy(i,j,k,idblock)=zero
	pxz(i,j,k,idblock)=zero
	pyy(i,j,k,idblock)=zero
	pyz(i,j,k,idblock)=zero
	pzz(i,j,k,idblock)=zero
	
	uh(i,j,k,idblock)=vxs
	vh(i,j,k,idblock)=vys
	wh(i,j,k,idblock)=vzs
	rhoh(i,j,k,idblock)=rhos
	pxxh(i,j,k,idblock)=zero
	pxyh(i,j,k,idblock)=zero
	pxzh(i,j,k,idblock)=zero
	pyyh(i,j,k,idblock)=zero
	pyzh(i,j,k,idblock)=zero
	pzzh(i,j,k,idblock)=zero
    
    return

 end subroutine setup_system
 
 attributes(global) subroutine store_print()
	
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
      rhoprint_d(gi,gj,gk)=rho(i,j,k,idblock)
	  velprint_d(1,gi,gj,gk)=u(i,j,k,idblock)
	  velprint_d(2,gi,gj,gk)=v(i,j,k,idblock)
	  velprint_d(3,gi,gj,gk)=w(i,j,k,idblock)
		
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
      rhoprint_d(gi,gj,gk)=rhoh(i,j,k,idblock)
	  velprint_d(1,gi,gj,gk)=uh(i,j,k,idblock)
	  velprint_d(2,gi,gj,gk)=vh(i,j,k,idblock)
	  velprint_d(3,gi,gj,gk)=wh(i,j,k,idblock)
		
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
