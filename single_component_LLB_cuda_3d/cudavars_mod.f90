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
    type (dim3) :: dimGrid,dimBlock,dimGridx,dimGridy,dimGridz,dimBlock2,dimBlockshared,&
     dimGridhalo,dimBlockhalo
    
    integer, constant :: TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
    
    real(kind=db), allocatable, dimension(:,:,:), device :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:), device :: rhoh,uh,vh,wh,pxxh,pxyh,pxzh,pyyh,pyzh,pzzh
    real(kind=4), allocatable, dimension(:,:,:), device :: rhoprint_d
    real(kind=4), allocatable, dimension(:,:,:,:), device :: velprint_d
    integer(kind=1), allocatable, dimension(:,:,:), device   :: isfluid
    
    contains
    
    attributes(global) subroutine setup_system(rhos,vxs,vys,vzs)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    
    integer :: i,j,k
          
            
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x -1
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y -1
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z -1
	
	if (i>nx_d+1 .or. j>ny_d+1 .or. k>nz_d+1)return
	
	u(i,j,k)=vxs
	v(i,j,k)=vys
	w(i,j,k)=vzs
	rho(i,j,k)=rhos
	pxx(i,j,k)=zero
	pxy(i,j,k)=zero
	pxz(i,j,k)=zero
	pyy(i,j,k)=zero
	pyz(i,j,k)=zero
	pzz(i,j,k)=zero
	
	uh(i,j,k)=vxs
	vh(i,j,k)=vys
	wh(i,j,k)=vzs
	rhoh(i,j,k)=rhos
	pxxh(i,j,k)=zero
	pxyh(i,j,k)=zero
	pxzh(i,j,k)=zero
	pyyh(i,j,k)=zero
	pyzh(i,j,k)=zero
	pzzh(i,j,k)=zero
    
    return

 end subroutine setup_system
 
 attributes(global) subroutine store_print()
	  
	integer :: i,j,k
	
	  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	  
	  
	!write(*,*)i,j,p_d(0)*myrho_d
	if(abs(isfluid(i,j,k)).eq.1)then
      rhoprint_d(i,j,k)=rho(i,j,k)
	  velprint_d(1,i,j,k)=u(i,j,k)
	  velprint_d(2,i,j,k)=v(i,j,k)
	  velprint_d(3,i,j,k)=w(i,j,k)
		
	else
	  
	  rhoprint_d(i,j,k)=zero
	  velprint_d(1,i,j,k)=zero
	  velprint_d(2,i,j,k)=zero
	  velprint_d(3,i,j,k)=zero
	  
	endif
	  
	  return

 end subroutine store_print
 
 attributes(global) subroutine store_print_flop()
	  
	integer :: i,j,k
	
	  
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	  
	  
	!write(*,*)i,j,p_d(0)*myrho_d
	if(abs(isfluid(i,j,k)).eq.1)then
      rhoprint_d(i,j,k)=rhoh(i,j,k)
	  velprint_d(1,i,j,k)=uh(i,j,k)
	  velprint_d(2,i,j,k)=vh(i,j,k)
	  velprint_d(3,i,j,k)=wh(i,j,k)
		
	else
	  
	  rhoprint_d(i,j,k)=zero
	  velprint_d(1,i,j,k)=zero
	  velprint_d(2,i,j,k)=zero
	  velprint_d(3,i,j,k)=zero
	  
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
