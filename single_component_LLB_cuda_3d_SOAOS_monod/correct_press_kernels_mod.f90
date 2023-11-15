#include "defines.h"
 module correct_press_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine correct_pressure
    
    implicit none
    
    integer :: i,j,k,gi,gj,gk,idblock
     
	real(kind=db) :: uu,udotc,temp,feq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
		  
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	
	if(abs(isfluid(gi,gj,gk)).ne.1)return
                        
	uu=half*(hfieldsh(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))*&
	 hfieldsh(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10)) + &
	 hfieldsh(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10))*&
	 hfieldsh(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10)) + &
	 hfieldsh(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10))*&
	 hfieldsh(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	!1-2
	udotc=hfieldsh(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=feq
	feq=p1*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq

	!3-4
	udotc=hfieldsh(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10))/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pyy=feq
	feq=p1*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pyy=temp_pyy+feq
	!5-6
	udotc=hfieldsh(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10))/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pzz=feq
	feq=p1*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pzz=temp_pzz+feq
	!7-8
	udotc=(hfieldsh(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=feq
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy+feq
	!10-9
	udotc=(-hfieldsh(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	!11-12
	udotc=(hfieldsh(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=feq
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz+feq
	!13-14
	udotc=(hfieldsh(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10))-&
	 hfieldsh(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	!15-16
	udotc=(hfieldsh(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=feq
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz+feq
	!17-18
	udotc=(-hfieldsh(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))+&
	 hfieldsh(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	feq=p2*(hfieldsh(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	!ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
	
	hfieldsh(idx5d(i,j,k,5,idblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,5,idblock,nx_d,ny_d,nz_d,10))-temp_pxx
	hfieldsh(idx5d(i,j,k,6,idblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,6,idblock,nx_d,ny_d,nz_d,10))-temp_pyy
	hfieldsh(idx5d(i,j,k,7,idblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,7,idblock,nx_d,ny_d,nz_d,10))-temp_pzz
	hfieldsh(idx5d(i,j,k,8,idblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,8,idblock,nx_d,ny_d,nz_d,10))-temp_pxy
	hfieldsh(idx5d(i,j,k,9,idblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,9,idblock,nx_d,ny_d,nz_d,10))-temp_pxz
	hfieldsh(idx5d(i,j,k,10,idblock,nx_d,ny_d,nz_d,10))=&
	 hfieldsh(idx5d(i,j,k,10,idblock,nx_d,ny_d,nz_d,10))-temp_pyz
    
    return
    
  end subroutine correct_pressure
  
  attributes(global) subroutine correct_pressure_flop
    
    implicit none
    
    integer :: i,j,k,gi,gj,gk,idblock
    
	real(kind=db) :: uu,udotc,temp,feq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
		  
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	
	if(abs(isfluid(gi,gj,gk)).ne.1)return
	
    uu=half*(hfields(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10)) + &
     hfields(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10)) + &
     hfields(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10))*&
     hfields(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
    
	!1-2
	udotc=hfields(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=feq
	feq=p1*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq

	!3-4
	udotc=hfields(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10))/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pyy=feq
	feq=p1*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pyy=temp_pyy+feq
	!5-6
	udotc=hfields(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10))/cssq
	temp = -uu + half*udotc*udotc
	feq=p1*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pzz=feq
	feq=p1*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pzz=temp_pzz+feq
	!7-8
	udotc=(hfields(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=feq
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy+feq
	!10-9
	udotc=(-hfields(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pyy=temp_pyy+feq
	temp_pxy=temp_pxy-feq
	!11-12
	udotc=(hfields(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=feq
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz+feq
	!13-14
	udotc=(hfields(idx5d(i,j,k,3,idblock,nx_d,ny_d,nz_d,10))-&
	 hfields(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pyy=temp_pyy+feq
	temp_pzz=temp_pzz+feq
	temp_pyz=temp_pyz-feq
	!15-16
	udotc=(hfields(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=feq
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz+feq
	!17-18
	udotc=(-hfields(idx5d(i,j,k,2,idblock,nx_d,ny_d,nz_d,10))+&
	 hfields(idx5d(i,j,k,4,idblock,nx_d,ny_d,nz_d,10)))/cssq
	temp = -uu + half*udotc*udotc
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp + udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	feq=p2*(hfields(idx5d(i,j,k,1,idblock,nx_d,ny_d,nz_d,10))+(temp - udotc))
	temp_pxx=temp_pxx+feq
	temp_pzz=temp_pzz+feq
	temp_pxz=temp_pxz-feq
	!ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
	
	hfields(idx5d(i,j,k,5,idblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,5,idblock,nx_d,ny_d,nz_d,10))-temp_pxx
	hfields(idx5d(i,j,k,6,idblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,6,idblock,nx_d,ny_d,nz_d,10))-temp_pyy
	hfields(idx5d(i,j,k,7,idblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,7,idblock,nx_d,ny_d,nz_d,10))-temp_pzz
	hfields(idx5d(i,j,k,8,idblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,8,idblock,nx_d,ny_d,nz_d,10))-temp_pxy
	hfields(idx5d(i,j,k,9,idblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,9,idblock,nx_d,ny_d,nz_d,10))-temp_pxz
	hfields(idx5d(i,j,k,10,idblock,nx_d,ny_d,nz_d,10))=&
	 hfields(idx5d(i,j,k,10,idblock,nx_d,ny_d,nz_d,10))-temp_pyz
	
	return
    
  end subroutine correct_pressure_flop
  
 end module correct_press_kernels
