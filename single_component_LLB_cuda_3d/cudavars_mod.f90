#include "defines.h"
 module cudavars
   
    use cudafor
    
    implicit none
    
    integer, parameter :: db=kind(1.e0)
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    real(kind=db),parameter :: p0 = (1.0_db/3.0_db)
    real(kind=db),parameter :: p1 = (1.0_db/18.0_db)
    real(kind=db),parameter :: p2 = (1.0_db/36.0_db)
    real(kind=db),parameter :: cssq = 1.0_db/3.0_db
    real(kind=db),parameter :: onecssq = 3.0_db
    real(kind=db),parameter :: halfonecssq = 1.5_db
    real(kind=db),parameter :: p1dcssq=p1/cssq
    real(kind=db),parameter :: p2dcssq=p2/cssq
    
    real(kind=db),parameter :: pi2cssq0=p0/(2.0_db*cssq**2.0_db)
    real(kind=db),parameter :: pi2cssq1=p1/(2.0_db*cssq**2.0_db)
    real(kind=db),parameter :: pi2cssq2=p2/(2.0_db*cssq**2.0_db)

    real(kind=db),parameter :: qxx=1.0_db-cssq
    real(kind=db),parameter :: qyy=1.0_db-cssq
    real(kind=db),parameter :: qzz=1.0_db-cssq
    real(kind=db),parameter :: qxy_7_8=1.0_db
    real(kind=db),parameter :: qxy_9_10=-1.0_db
    real(kind=db),parameter :: qxz_15_16=1.0_db
    real(kind=db),parameter :: qxz_17_18=-1.0_db
    real(kind=db),parameter :: qyz_11_12=1.0_db
    real(kind=db),parameter :: qyz_13_14=-1.0_db
    
    integer, constant :: nx_d,ny_d,nz_d
    
    real(kind=db), constant :: omega,fx,fy,fz,oneminusomega
    integer(kind=cuda_Stream_Kind) :: stream1,stream2
    type (cudaDeviceProp) :: prop
    type (cudaEvent) :: startEvent, stopEvent, dummyEvent, dummyEvent1, dummyEvent2
    type (dim3) :: dimGrid,dimBlock,dimGridx,dimGridy,dimBlock2
    
    integer, constant :: TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
    
    real(kind=db), allocatable, dimension(:,:,:), device :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:), device :: rhoh,uh,vh,wh,pxxh,pxyh,pxzh,pyyh,pyzh,pzzh
    real(kind=4), allocatable, dimension(:,:,:), device :: rhoprint_d
    real(kind=4), allocatable, dimension(:,:,:,:), device :: velprint_d
    integer(kind=4), allocatable, dimension(:,:,:), device   :: isfluid
    
    contains
    
    attributes(global) subroutine setup_system(rhos,vxs,vys,vzs)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    
    integer :: i,j,k
          
            
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	
	u(i,j,k)=vxs
	v(i,j,k)=vys
	w(i,j,k)=vzs
	rho(i,j,k)=rhos
	pxx(i,j,k)=0.0_db
	pxy(i,j,k)=0.0_db
	pxz(i,j,k)=0.0_db
	pyy(i,j,k)=0.0_db
	pyz(i,j,k)=0.0_db
	pzz(i,j,k)=0.0_db
	
	uh(i,j,k)=vxs
	vh(i,j,k)=vys
	wh(i,j,k)=vzs
	rhoh(i,j,k)=rhos
	pxxh(i,j,k)=0.0_db
	pxyh(i,j,k)=0.0_db
	pxzh(i,j,k)=0.0_db
	pyyh(i,j,k)=0.0_db
	pyzh(i,j,k)=0.0_db
	pzzh(i,j,k)=0.0_db
    
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
	  
	  rhoprint_d(i,j,k)=0.0_db
	  velprint_d(1,i,j,k)=0.0_db
	  velprint_d(2,i,j,k)=0.0_db
	  velprint_d(3,i,j,k)=0.0_db
	  
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
	  
	  rhoprint_d(i,j,k)=0.0_db
	  velprint_d(1,i,j,k)=0.0_db
	  velprint_d(2,i,j,k)=0.0_db
	  velprint_d(3,i,j,k)=0.0_db
	  
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
