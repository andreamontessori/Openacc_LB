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
    type (dim3) :: dimGrid,dimBlock,dimGridx,dimGridy,dimGridz,dimBlock2, &
     dimGridhalo,dimBlockhalo,dimGridshared,dimBlockshared
    
    integer, constant :: TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
    
    integer, constant :: nxblock_d,nxyblock_d,nblocks_d
    
    real(kind=db), allocatable, dimension(:,:,:,:), device :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:,:), device :: rhoh,uh,vh,wh,pxxh,pxyh,pxzh,pyyh,pyzh,pzzh
    real(kind=4), allocatable, dimension(:,:,:), device :: rhoprint_d
    real(kind=4), allocatable, dimension(:,:,:,:), device :: velprint_d
    integer(kind=1), allocatable, dimension(:,:,:), device   :: isfluid
    
    contains
    
    attributes(global) subroutine setup_system(rhos,vxs,vys,vzs)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    real :: mytest
    
    integer :: i,j,k,gi,gj,gk,idblock       
            
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	
	mytest=rhos!real(gi**3+gj**3+gk**3)
	
	u(i,j,k,idblock)=vxs
	v(i,j,k,idblock)=vys
	w(i,j,k,idblock)=vzs
	rho(i,j,k,idblock)=mytest!rhos
	pxx(i,j,k,idblock)=zero
	pxy(i,j,k,idblock)=zero
	pxz(i,j,k,idblock)=zero
	pyy(i,j,k,idblock)=zero
	pyz(i,j,k,idblock)=zero
	pzz(i,j,k,idblock)=zero
	
	uh(i,j,k,idblock)=vxs
	vh(i,j,k,idblock)=vys
	wh(i,j,k,idblock)=vzs
	rhoh(i,j,k,idblock)=mytest!rhos
	pxxh(i,j,k,idblock)=zero
	pxyh(i,j,k,idblock)=zero
	pxzh(i,j,k,idblock)=zero
	pyyh(i,j,k,idblock)=zero
	pyzh(i,j,k,idblock)=zero
	pzzh(i,j,k,idblock)=zero
    
    return

 end subroutine setup_system
 
 attributes(global) subroutine setup_system_halo(rhos,vxs,vys,vzs)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    real :: mytest
    
    integer :: ii,jj,kk,xblock,yblock,zblock,idblock
    integer :: i,j,k,gi,gj,gk       
            
	gi = (blockIdx%x-2) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-2) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-2) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
    
    xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
	zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
	
	idblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	!if(gi==6)write(*,*)'cazzo',gi,gj,gk,idblock
	
	mytest=rhos!real(gi**3+gj**3+gk**3)
	
	u(i,j,k,idblock)=vxs
	v(i,j,k,idblock)=vys
	w(i,j,k,idblock)=vzs
	rho(i,j,k,idblock)=mytest!rhos
	pxx(i,j,k,idblock)=zero
	pxy(i,j,k,idblock)=zero
	pxz(i,j,k,idblock)=zero
	pyy(i,j,k,idblock)=zero
	pyz(i,j,k,idblock)=zero
	pzz(i,j,k,idblock)=zero
	
	uh(i,j,k,idblock)=vxs
	vh(i,j,k,idblock)=vys
	wh(i,j,k,idblock)=vzs
	rhoh(i,j,k,idblock)=mytest!rhos
	pxxh(i,j,k,idblock)=zero
	pxyh(i,j,k,idblock)=zero
	pxzh(i,j,k,idblock)=zero
	pyyh(i,j,k,idblock)=zero
	pyzh(i,j,k,idblock)=zero
	pzzh(i,j,k,idblock)=zero
    
    return

 end subroutine setup_system_halo
 
  attributes(global) subroutine setup_system_halo2(rhos,vxs,vys,vzs,idblock,ii,jj,kk,iii,jjj,kkk)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    integer, value ::idblock,ii,jj,kk,iii,jjj,kkk
    real :: mytest
    
    integer :: i,j,k,gi,gj,gk,myblock,coordblock_d(3),xblock,yblock,zblock
            
	gi = (blockIdx%x-2) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-2) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-2) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
	yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
	zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
	
	myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	
	!if(gi==1 .and. gj==1 .and. gk==1)
	!write(*,*)gi,gj,gk,blockIdx%x,blockIdx%y,blockIdx%z
	
	if(gi==ii .and. gj==jj .and. gk==kk)then
	coordblock_d(3)=(idblock-1)/nxyblock_d+1
    coordblock_d(2)=((idblock-1)-(coordblock_d(3)-1)*nxyblock_d)/nxblock_d +1
    coordblock_d(1)=(idblock-1)-(coordblock_d(3)-1)*nxyblock_d-(coordblock_d(2)-1)*nxblock_d+1
	!write(*,*)ii,jj,kk
	  if(idblock .ne. myblock)write(*,*)'SONO CAZZo1',idblock,myblock
  	  if(i .ne. iii)write(*,*)'SONO CAZZI2',idblock,myblock
  	  if(j .ne. jjj)write(*,*)'SONO CAZZI3',idblock,myblock
  	  if(k .ne. kkk)write(*,*)'SONO CAZZI4',idblock,myblock
	endif
	
  end subroutine setup_system_halo2
 
  attributes(global) subroutine setup_system_bulk2(rhos,vxs,vys,vzs,idblock,ii,jj,kk,iii,jjj,kkk)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    integer, value ::idblock,ii,jj,kk,iii,jjj,kkk
    real :: mytest
    
    integer :: i,j,k,gi,gj,gk,myblock,coordblock_d(3)       
            
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	
	if(gi==ii .and. gj==jj .and. gk==kk)then
	coordblock_d(3)=(myblock-1)/nxyblock_d+1
    coordblock_d(2)=((myblock-1)-(coordblock_d(3)-1)*nxyblock_d)/nxblock_d +1
    coordblock_d(1)=(myblock-1)-(coordblock_d(3)-1)*nxyblock_d-(coordblock_d(2)-1)*nxblock_d+1
	!write(*,*)ii,jj,kk
	  if(idblock .ne. myblock)write(*,*)'SONO CAZZI2',idblock,myblock
  	  if(i .ne. iii)write(*,*)'SONO CAZZI2',idblock,myblock
  	  if(j .ne. jjj)write(*,*)'SONO CAZZI2',idblock,myblock
  	  if(k .ne. kkk)write(*,*)'SONO CAZZI2',idblock,myblock
	endif
	
 end subroutine setup_system_bulk2
 
 attributes(global) subroutine setup_system_bulk3(rhos,vxs,vys,vzs,idblock,ii,jj,kk,iii,jjj,kkk)
    
    real(kind=db), value :: rhos,vxs,vys,vzs
    integer, value ::idblock,ii,jj,kk,iii,jjj,kkk,myxb,myyb,myzb
    real :: mytest
    
    integer :: i,j,k,gi,gj,gk,myblock,coordblock_d(3),xblock,yblock,zblock,iidblock 
            
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x -1
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y -1
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z -1
	
	ii=threadIdx%x-1
	jj=threadIdx%y-1
	kk=threadIdx%z-1
	
	xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
	
	i=gi-xblock*TILE_DIMx_d+2*TILE_DIMx_d
	j=gj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
    k=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
	
	myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	!iidblock=(blockIdx%x)+(blockIdx%y)*nxblock_d+(blockIdx%z)*nxyblock_d+1
	
	myxb=blockIdx%x+1
	myyb=blockIdx%y+1
	myzb=blockIdx%z+1
	
	if(ii>=1 .and. jj>=1 .and. kk>=1 .and. ii<=TILE_DIMx_d .and. jj<=TILE_DIMy_d .and. kk<=TILE_DIMz_d)then
	  if(myxb==2 .and. myyb==2 .and. myzb==2)write(*,*)'dentro',gi,gj,gk,myblock
	else
	  if(myxb==2 .and. myyb==2 .and. myzb==2)write(*,*)'fuori',gi,gj,gk,myblock
	endif
	
	if(ii<1 .or. jj<1 .or. kk<1)return
	if(ii>TILE_DIMx_d .or. jj>TILE_DIMy_d .or. kk>TILE_DIMz_d)return
	
	!write(*,*)gi,gj,gk,myblock
	
!	if(gi==ii .and. gj==jj .and. gk==kk)then
!	coordblock_d(3)=(myblock-1)/nxyblock_d+1
!    coordblock_d(2)=((myblock-1)-(coordblock_d(3)-1)*nxyblock_d)/nxblock_d +1
!    coordblock_d(1)=(myblock-1)-(coordblock_d(3)-1)*nxyblock_d-(coordblock_d(2)-1)*nxblock_d+1
!	!write(*,*)ii,jj,kk
!	  if(idblock .ne. iidblock)write(*,*)'SONO CAZZI2',idblock,myblock,iidblock
!  	  if(i .ne. iii)write(*,*)'SONO CAZZI3',i,iii,iidblock
!  	  if(j .ne. jjj)write(*,*)'SONO CAZZI4',j,jjj,iidblock
!  	  if(k .ne. kkk)write(*,*)'SONO CAZZI5',k,kkk,iidblock
!	endif
	
 end subroutine setup_system_bulk3
 
 attributes(global) subroutine store_print()
	
	integer :: i,j,k,gi,gj,gk,idblock,ii,jj,kk,xblock,yblock,zblock,idblocko
    real :: mytest
    
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    !xblock=(gi+2*TILE_DIMx_d-1)/TILE_DIMx_d
    !yblock=(gj+2*TILE_DIMy_d-1)/TILE_DIMy_d
    !zblock=(gk+2*TILE_DIMz_d-1)/TILE_DIMz_d
    
    !ii=gi-xblock*TILE_DIMx_d+2*TILE_DIMx_d
    !jj=gj-yblock*TILE_DIMy_d+2*TILE_DIMy_d
    !kk=gk-zblock*TILE_DIMz_d+2*TILE_DIMz_d
    
    !idblocko=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
	!if(idblocko.ne. idblock)write(*,*)'cazzo amaro!'
	!mytest=real(gi**3+gj**3+gk**3)
	
	!if(rho(ii,jj,kk,idblocko).ne. mytest)write(*,*)'cazzo'
	  
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
