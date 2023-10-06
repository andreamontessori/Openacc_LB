#include "defines.h"
 module streamcoll_bc_kernels
 
  use cudavars
  
  implicit none
  
  contains
  
  attributes(global) subroutine streamcoll_bc()
	
	implicit none  
	  
    integer :: i,j,k,gi,gj,gk,idblock
    !logical :: alltrue
	real(kind=db) :: uu,udotc,temp,uu0
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
    
    
    
	
	gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    !alltrue = .false.
    !if( allthreads(isfluid(gi,gj,gk).ne.-1) )alltrue = .true.
    !call syncthreads
    !if(alltrue)return
    
    if(isfluid(gi,gj,gk).ne.-1)return
    
	uu0=halfonecssq*(hfields(i,j,k,2,idblock)*hfields(i,j,k,2,idblock) &
	+ hfields(i,j,k,3,idblock)*hfields(i,j,k,3,idblock) + hfields(i,j,k,4,idblock)*hfields(i,j,k,4,idblock))
	!0
	temp_rho=p0*(hfields(i,j,k,1,idblock)-uu0)
	temp_rho=temp_rho + oneminusomega*pi2cssq0*( &
	-cssq*(hfields(i,j,k,6,idblock)+hfields(i,j,k,4,idblock)+hfields(i,j,k,7,idblock)))
	
	!1
	if(isfluid(gi-1,gj,gk).eq.0)then
	  udotc=hfields(i,j,k,2,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,4,idblock) &
	  -cssq*(hfields(i,j,k,6,idblock)+hfields(i,j,k,7,idblock)))
	  temp_pop=temp_pop - fx*p1dcssq
	else
	  uu=halfonecssq*(hfields(i-1,j,k,2,idblock)*hfields(i-1,j,k,2,idblock) &
	  + hfields(i-1,j,k,3,idblock)*hfields(i-1,j,k,3,idblock) + hfields(i-1,j,k,4,idblock)*hfields(i-1,j,k,4,idblock))
	  udotc=hfields(i-1,j,k,2,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfields(i-1,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*hfields(i-1,j,k,5,idblock) &
	  -cssq*(hfields(i-1,j,k,6,idblock)+hfields(i-1,j,k,7,idblock)))
	  temp_pop=temp_pop + fx*p1dcssq
	endif
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop
	
	!2
	if(isfluid(gi+1,gj,gk).eq.0)then
	  udotc=hfields(i,j,k,2,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,4,idblock) &
	   -cssq*(hfields(i,j,k,6,idblock)+hfields(i,j,k,7,idblock)))
	  temp_pop=temp_pop + fx*p1dcssq
	else
	  uu=halfonecssq*(hfields(i+1,j,k,2,idblock)*hfields(i+1,j,k,2,idblock) &
	   + hfields(i+1,j,k,3,idblock)*hfields(i+1,j,k,3,idblock) + hfields(i+1,j,k,4,idblock)*hfields(i+1,j,k,4,idblock))
	  udotc=hfields(i+1,j,k,2,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfields(i+1,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*hfields(i+1,j,k,5,idblock) &
	   -cssq*(hfields(i+1,j,k,6,idblock)+hfields(i+1,j,k,7,idblock)))
	  temp_pop=temp_pop - fx*p1dcssq
	endif
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
	
	!3
	if(isfluid(gi,gj-1,gk).eq.0)then
	  udotc=hfields(i,j,k,3,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,idblock) &
	   -cssq*(hfields(i,j,k,4,idblock)+hfields(i,j,k,7,idblock)))
	  temp_pop=temp_pop - fy*p1dcssq
	else
	  uu=halfonecssq*(hfields(i,j-1,k,2,idblock)*hfields(i,j-1,k,2,idblock) &
	   + hfields(i,j-1,k,3,idblock)*hfields(i,j-1,k,3,idblock) + hfields(i,j-1,k,4,idblock)*hfields(i,j-1,k,4,idblock))
	  udotc=hfields(i,j-1,k,3,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j-1,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*hfields(i,j-1,k,6,idblock) &
	   -cssq*(hfields(i,j-1,k,5,idblock)+hfields(i,j-1,k,7,idblock)))
	  temp_pop=temp_pop + fy*p1dcssq
	endif
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	if(isfluid(i,j+1,k).eq.0)then
	  udotc=hfields(i,j,k,3,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,idblock) &
	   -cssq*(hfields(i,j,k,4,idblock)+hfields(i,j,k,7,idblock)))
	  temp_pop=temp_pop + fy*p1dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j+1,k,2,idblock)*hfieldsh(i,j+1,k,2,idblock) &
	   + hfieldsh(i,j+1,k,3,idblock)*hfieldsh(i,j+1,k,3,idblock) + hfieldsh(i,j+1,k,4,idblock)*hfieldsh(i,j+1,k,4,idblock))
	  udotc=hfieldsh(i,j+1,k,3,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j+1,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*hfieldsh(i,j+1,k,6,idblock) &
	   -cssq*(hfieldsh(i,j+1,k,5,idblock)+hfieldsh(i,j+1,k,7,idblock)))
	  temp_pop=temp_pop - fy*p1dcssq
	endif
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	if(isfluid(i,j,k-1).eq.0)then
	  udotc=hfields(i,j,k,4,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,idblock) &
	   -cssq*(hfields(i,j,k,4,idblock)+hfields(i,j,k,6,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(hfields(i,j,k-1,2,idblock)*hfields(i,j,k-1,2,idblock) &
	   + hfields(i,j,k-1,3,idblock)*hfields(i,j,k-1,3,idblock) + hfields(i,j,k-1,4,idblock)*hfields(i,j,k-1,4,idblock))
	  udotc=hfields(i,j,k-1,4,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j,k-1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*hfields(i,j,k-1,7,idblock) &
	   -cssq*(hfields(i,j,k-1,5,idblock)+hfields(i,j,k-1,6,idblock)))
	  temp_pop=temp_pop + fz*p1dcssq
	endif
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop
	
	!6 +1
	if(isfluid(i,j,k+1).eq.0)then
	  udotc=hfields(i,j,k,4,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,idblock) &
	   -cssq*(hfields(i,j,k,4,idblock)+hfields(i,j,k,6,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(hfields(i,j,k+1,2,idblock)*hfields(i,j,k+1,2,idblock) &
	   + hfields(i,j,k+1,3,idblock)*hfields(i,j,k+1,3,idblock) + hfields(i,j,k+1,4,idblock)*hfields(i,j,k+1,4,idblock))
	  udotc=hfields(i,j,k+1,4,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfields(i,j,k+1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*hfields(i,j,k+1,7,idblock) &
	   -cssq*(hfields(i,j,k+1,5,idblock)+hfields(i,j,k+1,6,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	endif
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
	
	!7 -1 -1
	if(isfluid(i-1,j-1,k).eq.0)then
	  udotc=(hfields(i,j,k,2,idblock)+hfields(i,j,k,3,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,4,idblock)+qyy*hfields(i,j,k,6,idblock) &
	   -cssq*hfields(i,j,k,7,idblock)+two*qxy_7_8*hfields(i,j,k,8,idblock))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i-1,j-1,k,2,idblock)*hfields(i-1,j-1,k,2,idblock) &
	   + hfields(i-1,j-1,k,3,idblock)*hfields(i-1,j-1,k,3,idblock) + hfields(i-1,j-1,k,4,idblock)*hfields(i-1,j-1,k,4,idblock))
	  udotc=(hfields(i-1,j-1,k,2,idblock)+hfields(i-1,j-1,k,3,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i-1,j-1,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i-1,j-1,k,5,idblock)+qyy*hfields(i-1,j-1,k,6,idblock) &
	   -cssq*hfields(i-1,j-1,k,7,idblock)+two*qxy_7_8*hfields(i-1,j-1,k,8,idblock))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	endif
	!+1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pop
	
	!8 +1 +1
	if(isfluid(i+1,j+1,k).eq.0)then
	  udotc=(hfields(i,j,k,2,idblock)+hfields(i,j,k,3,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,4,idblock)+qyy*hfields(i,j,k,6,idblock) &
	   -cssq*hfields(i,j,k,7,idblock)+two*qxy_7_8*hfields(i,j,k,8,idblock))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	else
	  uu=halfonecssq*(hfields(i+1,j+1,k,2,idblock)*hfields(i+1,j+1,k,2,idblock) &
	   + hfields(i+1,j+1,k,3,idblock)*hfields(i+1,j+1,k,3,idblock) + hfields(i+1,j+1,k,4,idblock)*hfields(i+1,j+1,k,4,idblock))
	  udotc=(hfields(i+1,j+1,k,2,idblock)+hfields(i+1,j+1,k,3,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i+1,j+1,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i+1,j+1,k,5,idblock)+qyy*hfields(i+1,j+1,k,6,idblock) &
	   -cssq*hfields(i+1,j+1,k,7,idblock)+two*qxy_7_8*hfields(i+1,j+1,k,8,idblock))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	endif
	!-1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy+temp_pop
	
	!10   +1 -1
	if(isfluid(i+1,j-1,k).eq.0)then
	  udotc=(-hfields(i,j,k,2,idblock)+hfields(i,j,k,3,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,4,idblock)+qyy*hfields(i,j,k,6,idblock) &
	   -cssq*hfields(i,j,k,7,idblock)+two*qxy_9_10*hfields(i,j,k,8,idblock))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i+1,j-1,k,2,idblock)*hfields(i+1,j-1,k,2,idblock) &
	   + hfields(i+1,j-1,k,3,idblock)*hfields(i+1,j-1,k,3,idblock) + hfields(i+1,j-1,k,4,idblock)*hfields(i+1,j-1,k,4,idblock))
	  udotc=(-hfields(i+1,j-1,k,2,idblock)+hfields(i+1,j-1,k,3,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i+1,j-1,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i+1,j-1,k,5,idblock)+qyy*hfields(i+1,j-1,k,6,idblock) &
	   -cssq*hfields(i+1,j-1,k,7,idblock)+two*qxy_9_10*hfields(i+1,j-1,k,8,idblock))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	endif
	!-1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
	
	!9  -1 +1
	if(isfluid(i-1,j+1,k).eq.0)then
	  udotc=(-hfields(i,j,k,2,idblock)+hfields(i,j,k,3,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,4,idblock)+qyy*hfields(i,j,k,6,idblock) &
	   -cssq*hfields(i,j,k,7,idblock)+two*qxy_9_10*hfields(i,j,k,8,idblock))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i-1,j+1,k,2,idblock)*hfields(i-1,j+1,k,2,idblock) &
	   + hfields(i-1,j+1,k,3,idblock)*hfields(i-1,j+1,k,3,idblock) + hfields(i-1,j+1,k,4,idblock)*hfields(i-1,j+1,k,4,idblock))
	  udotc=(-hfields(i-1,j+1,k,2,idblock)+hfields(i-1,j+1,k,3,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i-1,j+1,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i-1,j+1,k,5,idblock)+qyy*hfields(i-1,j+1,k,6,idblock) &
	   -cssq*hfields(i-1,j+1,k,7,idblock)+two*qxy_9_10*hfields(i-1,j+1,k,8,idblock))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	endif
	!+1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop

	!15  -1  -1
	if(isfluid(i-1,j,k-1).eq.0)then
	  udotc=(hfields(i,j,k,2,idblock)+hfields(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,4,idblock)+qzz*hfields(i,j,k,7,idblock) &
	   -cssq*hfields(i,j,k,6,idblock)+two*qxz_15_16*hfields(i,j,k,9,idblock))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i-1,j,k-1,2,idblock)*hfields(i-1,j,k-1,2,idblock) &
	   + hfields(i-1,j,k-1,3,idblock)*hfields(i-1,j,k-1,3,idblock) + hfields(i-1,j,k-1,4,idblock)*hfields(i-1,j,k-1,4,idblock))
	  udotc=(hfields(i-1,j,k-1,2,idblock)+hfields(i-1,j,k-1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i-1,j,k-1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i-1,j,k-1,5,idblock)+qzz*hfields(i-1,j,k-1,7,idblock) &
	   -cssq*hfields(i-1,j,k-1,6,idblock)+two*qxz_15_16*hfields(i-1,j,k-1,9,idblock))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq 
	endif
	!+1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pop
	
	!16  +1  +1
	if(isfluid(i+1,j,k+1).eq.0)then
	  udotc=(hfields(i,j,k,2,idblock)+hfields(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,4,idblock)+qzz*hfields(i,j,k,7,idblock) &
	   -cssq*hfields(i,j,k,6,idblock)+two*qxz_15_16*hfields(i,j,k,9,idblock))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i+1,j,k+1,2,idblock)*hfields(i+1,j,k+1,2,idblock) &
	   + hfields(i+1,j,k+1,3,idblock)*hfields(i+1,j,k+1,3,idblock) + hfields(i+1,j,k+1,4,idblock)*hfields(i+1,j,k+1,4,idblock))
	  udotc=(hfields(i+1,j,k+1,2,idblock)+hfields(i+1,j,k+1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i+1,j,k+1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i+1,j,k+1,5,idblock)+qzz*hfields(i+1,j,k+1,7,idblock) &
	   -cssq*hfields(i+1,j,k+1,6,idblock)+two*qxz_15_16*hfields(i+1,j,k+1,9,idblock))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	endif
	!-1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz+temp_pop

	!17  +1   -1
	if(isfluid(i+1,j,k-1).eq.0)then
	  udotc=(-hfields(i,j,k,2,idblock)+hfields(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,4,idblock)+qzz*hfields(i,j,k,7,idblock) &
	   -cssq*hfields(i,j,k,6,idblock)+two*qxz_17_18*hfields(i,j,k,9,idblock))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i+1,j,k-1,2,idblock)*hfields(i+1,j,k-1,2,idblock) &
	   + hfields(i+1,j,k-1,3,idblock)*hfields(i+1,j,k-1,3,idblock) + hfields(i+1,j,k-1,4,idblock)*hfields(i+1,j,k-1,4,idblock))
	  udotc=(-hfields(i+1,j,k-1,2,idblock)+hfields(i+1,j,k-1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i+1,j,k-1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i+1,j,k-1,5,idblock)+qzz*hfields(i+1,j,k-1,7,idblock) &
	   -cssq*hfields(i+1,j,k-1,6,idblock)+two*qxz_17_18*hfields(i+1,j,k-1,9,idblock))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	endif
	!-1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop
	
	!18   -1   +1
	if(isfluid(i-1,j,k+1).eq.0)then
	  udotc=(-hfields(i,j,k,2,idblock)+hfields(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,4,idblock)+qzz*hfields(i,j,k,7,idblock) &
	   -cssq*hfields(i,j,k,6,idblock)+two*qxz_17_18*hfields(i,j,k,9,idblock))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i-1,j,k+1,2,idblock)*hfields(i-1,j,k+1,2,idblock) &
	   + hfields(i-1,j,k+1,3,idblock)*hfields(i-1,j,k+1,3,idblock) + hfields(i-1,j,k+1,4,idblock)*hfields(i-1,j,k+1,4,idblock))
	  udotc=(-hfields(i-1,j,k+1,2,idblock)+hfields(i-1,j,k+1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i-1,j,k+1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfields(i-1,j,k+1,5,idblock)+qzz*hfields(i-1,j,k+1,7,idblock) &
	   -cssq*hfields(i-1,j,k+1,6,idblock)+two*qxz_17_18*hfields(i-1,j,k+1,9,idblock))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	endif
	!+1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!11  -1  -1
	if(isfluid(i,j-1,k-1).eq.0)then
	  udotc=(hfields(i,j,k,3,idblock)+hfields(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,idblock)+qzz*hfields(i,j,k,7,idblock) &
	   -cssq*hfields(i,j,k,4,idblock)+two*qyz_11_12*hfields(i,j,k,10,idblock))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i,j-1,k-1,2,idblock)*hfields(i,j-1,k-1,2,idblock) &
	   + hfields(i,j-1,k-1,3,idblock)*hfields(i,j-1,k-1,3,idblock) + hfields(i,j-1,k-1,4,idblock)*hfields(i,j-1,k-1,4,idblock))
	  udotc=(hfields(i,j-1,k-1,3,idblock)+hfields(i,j-1,k-1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j-1,k-1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfields(i,j-1,k-1,6,idblock)+qzz*hfields(i,j-1,k-1,7,idblock) &
	   -cssq*hfields(i,j-1,k-1,5,idblock)+two*qyz_11_12*hfields(i,j-1,k-1,10,idblock))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	endif
	! 0 +1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pop
	
	!12   +1  +1
	if(isfluid(i,j+1,k+1).eq.0)then
	  udotc=(hfields(i,j,k,3,idblock)+hfields(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,idblock)+qzz*hfields(i,j,k,7,idblock) &
	   -cssq*hfields(i,j,k,4,idblock)+two*qyz_11_12*hfields(i,j,k,10,idblock))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i,j+1,k+1,2,idblock)*hfields(i,j+1,k+1,2,idblock) &
	   + hfields(i,j+1,k+1,3,idblock)*hfields(i,j+1,k+1,3,idblock) + hfields(i,j+1,k+1,4,idblock)*hfields(i,j+1,k+1,4,idblock))
	  udotc=(hfields(i,j+1,k+1,3,idblock)+hfields(i,j+1,k+1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j+1,k+1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfields(i,j+1,k+1,6,idblock)+qzz*hfields(i,j+1,k+1,7,idblock) &
	   -cssq*hfields(i,j+1,k+1,5,idblock)+two*qyz_11_12*hfields(i,j+1,k+1,10,idblock))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	endif
	! 0 -1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz+temp_pop

	!13   -1   +1
	if(isfluid(i,j-1,k+1).eq.0)then
	  udotc=(hfields(i,j,k,3,idblock)-hfields(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,idblock)+qzz*hfields(i,j,k,7,idblock) &
	   -cssq*hfields(i,j,k,4,idblock)+two*qyz_13_14*hfields(i,j,k,10,idblock))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i,j-1,k+1,2,idblock)*hfields(i,j-1,k+1,2,idblock) &
	   + hfields(i,j-1,k+1,3,idblock)*hfields(i,j-1,k+1,3,idblock) + hfields(i,j-1,k+1,4,idblock)*hfields(i,j-1,k+1,4,idblock))
	  udotc=(hfields(i,j-1,k+1,3,idblock)-hfields(i,j-1,k+1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j-1,k+1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfields(i,j-1,k+1,6,idblock)+qzz*hfields(i,j-1,k+1,7,idblock) &
	   -cssq*hfields(i,j-1,k+1,5,idblock)+two*qyz_13_14*hfields(i,j-1,k+1,10,idblock))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	endif
	! 0 +1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	!14  +1 -1
	if(isfluid(i,j+1,k-1).eq.0)then
	  udotc=(hfields(i,j,k,3,idblock)-hfields(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,idblock)+qzz*hfields(i,j,k,7,idblock) &
	   -cssq*hfields(i,j,k,4,idblock)+two*qyz_13_14*hfields(i,j,k,10,idblock))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	else
	  uu=halfonecssq*(hfields(i,j+1,k-1,2,idblock)*hfields(i,j+1,k-1,2,idblock) &
	   + hfields(i,j+1,k-1,3,idblock)*hfields(i,j+1,k-1,3,idblock) + hfields(i,j+1,k-1,4,idblock)*hfields(i,j+1,k-1,4,idblock))
	  udotc=(hfields(i,j+1,k-1,3,idblock)-hfields(i,j+1,k-1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfields(i,j+1,k-1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfields(i,j+1,k-1,6,idblock)+qzz*hfields(i,j+1,k-1,7,idblock) &
	   -cssq*hfields(i,j+1,k-1,5,idblock)+two*qyz_13_14*hfields(i,j+1,k-1,10,idblock))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	endif
	! 0 -1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	hfieldsh(i,j,k,1,idblock)=temp_rho
	
	hfieldsh(i,j,k,2,idblock)=temp_u
	hfieldsh(i,j,k,3,idblock)=temp_v
	hfieldsh(i,j,k,4,idblock)=temp_w
	
	hfieldsh(i,j,k,5,idblock)=temp_pxx
	hfieldsh(i,j,k,6,idblock)=temp_pyy
	hfieldsh(i,j,k,7,idblock)=temp_pzz
	hfieldsh(i,j,k,8,idblock)=temp_pxy
	hfieldsh(i,j,k,9,idblock)=temp_pxz
	hfieldsh(i,j,k,10,idblock)=temp_pyz

    
    return
    
  end subroutine streamcoll_bc
  
  attributes(global) subroutine streamcoll_bc_flop()
	
	implicit none  
	  
    integer :: i,j,k
    !logical :: alltrue
	real(kind=db) :: uu,udotc,temp,feq,uu0,fneq
	real(kind=db) :: temp_pop,temp_rho,temp_u,temp_v,temp_w
	real(kind=db) :: temp_pxx,temp_pyy,temp_pzz,temp_pxy,temp_pxz,temp_pyz
	integer :: idblock
	
	i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
    
    !alltrue = .false.
    !if( allthreads(isfluid(i,j,k).ne.-1) )alltrue = .true.
    !call syncthreads
    !if(alltrue)return
    
    if(isfluid(i,j,k).ne.-1)return
    
	uu0=halfonecssq*(hfieldsh(i,j,k,2,idblock)*hfieldsh(i,j,k,2,idblock) &
	 + hfieldsh(i,j,k,3,idblock)*hfieldsh(i,j,k,3,idblock) + hfieldsh(i,j,k,4,idblock)*hfieldsh(i,j,k,4,idblock))
	!0
	temp_rho=p0*(hfieldsh(i,j,k,1,idblock)-uu0)
	temp_rho=temp_rho + oneminusomega*pi2cssq0*( &
	 -cssq*(hfieldsh(i,j,k,6,idblock)+hfieldsh(i,j,k,5,idblock)+hfieldsh(i,j,k,7,idblock)))
	
	!1
	if(isfluid(i-1,j,k).eq.0)then
	  udotc=hfieldsh(i,j,k,2,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*hfieldsh(i,j,k,5,idblock) &
	   -cssq*(hfieldsh(i,j,k,6,idblock)+hfieldsh(i,j,k,7,idblock)))
	  temp_pop=temp_pop - fx*p1dcssq
	else
	  uu=halfonecssq*(hfieldsh(i-1,j,k,2,idblock)*hfieldsh(i-1,j,k,2,idblock) &
	   + hfieldsh(i-1,j,k,3,idblock)*hfieldsh(i-1,j,k,3,idblock) + hfieldsh(i-1,j,k,4,idblock)*hfieldsh(i-1,j,k,4,idblock))
	  udotc=hfieldsh(i-1,j,k,2,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i-1,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*hfieldsh(i-1,j,k,5,idblock) &
	   -cssq*(hfieldsh(i-1,j,k,6,idblock)+hfieldsh(i-1,j,k,7,idblock)))
	  temp_pop=temp_pop + fx*p1dcssq
	endif
	!+1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_pop
	temp_pxx=temp_pop
	
	!2
	if(isfluid(i+1,j,k).eq.0)then
	  udotc=hfieldsh(i,j,k,2,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*hfieldsh(i,j,k,5,idblock) &
	   -cssq*(hfieldsh(i,j,k,6,idblock)+hfieldsh(i,j,k,7,idblock)))
	  temp_pop=temp_pop + fx*p1dcssq
	else
	  uu=halfonecssq*(hfieldsh(i+1,j,k,2,idblock)*hfieldsh(i+1,j,k,2,idblock) &
	   + hfieldsh(i+1,j,k,3,idblock)*hfieldsh(i+1,j,k,3,idblock) + hfieldsh(i+1,j,k,4,idblock)*hfieldsh(i+1,j,k,4,idblock))
	  udotc=hfieldsh(i+1,j,k,2,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i+1,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qxx*hfieldsh(i+1,j,k,5,idblock) &
	   -cssq*(hfieldsh(i+1,j,k,6,idblock)+hfieldsh(i+1,j,k,7,idblock)))
	  temp_pop=temp_pop - fx*p1dcssq
	endif
	!-1  0  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_pxx=temp_pxx+temp_pop
	
	!3
	if(isfluid(i,j-1,k).eq.0)then
	  udotc=hfieldsh(i,j,k,3,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*hfieldsh(i,j,k,6,idblock) &
	   -cssq*(hfieldsh(i,j,k,5,idblock)+hfieldsh(i,j,k,7,idblock)))
	  temp_pop=temp_pop - fy*p1dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j-1,k,2,idblock)*hfieldsh(i,j-1,k,2,idblock) &
	   + hfieldsh(i,j-1,k,3,idblock)*hfieldsh(i,j-1,k,3,idblock) + hfieldsh(i,j-1,k,4,idblock)*hfieldsh(i,j-1,k,4,idblock))
	  udotc=hfieldsh(i,j-1,k,3,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j-1,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*hfieldsh(i,j-1,k,6,idblock) &
	   -cssq*(hfieldsh(i,j-1,k,5,idblock)+hfieldsh(i,j-1,k,7,idblock)))
	  temp_pop=temp_pop + fy*p1dcssq
	endif
	! 0 +1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_pop
	temp_pyy=temp_pop
	
	!4
	if(isfluid(i,j+1,k).eq.0)then
	  udotc=hfieldsh(i,j,k,3,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*hfieldsh(i,j,k,6,idblock) &
	   -cssq*(hfieldsh(i,j,k,5,idblock)+hfieldsh(i,j,k,7,idblock)))
	  temp_pop=temp_pop + fy*p1dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j+1,k,2,idblock)*hfieldsh(i,j+1,k,2,idblock) &
	   + hfieldsh(i,j+1,k,3,idblock)*hfieldsh(i,j+1,k,3,idblock) + hfieldsh(i,j+1,k,4,idblock)*hfieldsh(i,j+1,k,4,idblock))
	  udotc=hfieldsh(i,j+1,k,3,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j+1,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qyy*hfieldsh(i,j+1,k,6,idblock) &
	   -cssq*(hfieldsh(i,j+1,k,5,idblock)+hfieldsh(i,j+1,k,7,idblock)))
	  temp_pop=temp_pop - fy*p1dcssq
	endif
	! 0 -1  0
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_pyy=temp_pyy+temp_pop
	
	!5 -1
	if(isfluid(i,j,k-1).eq.0)then
	  udotc=hfieldsh(i,j,k,4,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*(hfieldsh(i,j,k,5,idblock)+hfieldsh(i,j,k,6,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j,k-1,2,idblock)*hfieldsh(i,j,k-1,2,idblock) &
	   + hfieldsh(i,j,k-1,3,idblock)*hfieldsh(i,j,k-1,3,idblock) + hfieldsh(i,j,k-1,4,idblock)*hfieldsh(i,j,k-1,4,idblock))
	  udotc=hfieldsh(i,j,k-1,4,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j,k-1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*hfieldsh(i,j,k-1,7,idblock) &
	   -cssq*(hfieldsh(i,j,k-1,5,idblock)+hfieldsh(i,j,k-1,6,idblock)))
	  temp_pop=temp_pop + fz*p1dcssq
	endif
	! 0  0 +1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_pop
	temp_pzz=temp_pop
	
	!6 +1
	if(isfluid(i,j,k+1).eq.0)then
	  udotc=hfieldsh(i,j,k,4,idblock)*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*(hfieldsh(i,j,k,5,idblock)+hfieldsh(i,j,k,6,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j,k+1,2,idblock)*hfieldsh(i,j,k+1,2,idblock) &
	   + hfieldsh(i,j,k+1,3,idblock)*hfieldsh(i,j,k+1,3,idblock) + hfieldsh(i,j,k+1,4,idblock)*hfieldsh(i,j,k+1,4,idblock))
	  udotc=hfieldsh(i,j,k+1,4,idblock)*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p1*(hfieldsh(i,j,k+1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq1*(qzz*hfieldsh(i,j,k+1,7,idblock) &
	   -cssq*(hfieldsh(i,j,k+1,5,idblock)+hfieldsh(i,j,k+1,6,idblock)))
	  temp_pop=temp_pop - fz*p1dcssq
	endif
	! 0  0 -1
	temp_rho=temp_rho+temp_pop
	temp_w=temp_w-temp_pop
	temp_pzz=temp_pzz+temp_pop
	
	!7 -1 -1
	if(isfluid(i-1,j-1,k).eq.0)then
	  udotc=(hfieldsh(i,j,k,2,idblock)+hfieldsh(i,j,k,3,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,idblock)+qyy*hfieldsh(i,j,k,6,idblock) &
	   -cssq*hfieldsh(i,j,k,7,idblock)+two*qxy_7_8*hfieldsh(i,j,k,8,idblock))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i-1,j-1,k,2,idblock)*hfieldsh(i-1,j-1,k,2,idblock) &
	   + hfieldsh(i-1,j-1,k,3,idblock)*hfieldsh(i-1,j-1,k,3,idblock) + hfieldsh(i-1,j-1,k,4,idblock)*hfieldsh(i-1,j-1,k,4,idblock))
	  udotc=(hfieldsh(i-1,j-1,k,2,idblock)+hfieldsh(i-1,j-1,k,3,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i-1,j-1,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i-1,j-1,k,5,idblock)+qyy*hfieldsh(i-1,j-1,k,6,idblock) &
	   -cssq*hfieldsh(i-1,j-1,k,7,idblock)+two*qxy_7_8*hfieldsh(i-1,j-1,k,8,idblock))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	endif
	!+1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pop
	
	!8 +1 +1
	if(isfluid(i+1,j+1,k).eq.0)then
	  udotc=(hfieldsh(i,j,k,2,idblock)+hfieldsh(i,j,k,3,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,idblock)+qyy*hfieldsh(i,j,k,6,idblock) &
	   -cssq*hfieldsh(i,j,k,7,idblock)+two*qxy_7_8*hfieldsh(i,j,k,8,idblock))
	  temp_pop=temp_pop + (fx+fy)*p2dcssq 
	else
	  uu=halfonecssq*(hfieldsh(i+1,j+1,k,2,idblock)*hfieldsh(i+1,j+1,k,2,idblock) &
	   + hfieldsh(i+1,j+1,k,3,idblock)*hfieldsh(i+1,j+1,k,3,idblock) + hfieldsh(i+1,j+1,k,4,idblock)*hfieldsh(i+1,j+1,k,4,idblock))
	  udotc=(hfieldsh(i+1,j+1,k,2,idblock)+hfieldsh(i+1,j+1,k,3,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i+1,j+1,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i+1,j+1,k,5,idblock)+qyy*hfieldsh(i+1,j+1,k,6,idblock) &
	   -cssq*hfieldsh(i+1,j+1,k,7,idblock)+two*qxy_7_8*hfieldsh(i+1,j+1,k,8,idblock))
	  temp_pop=temp_pop - (fx+fy)*p2dcssq
	endif
	!-1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy+temp_pop
	
	!10   +1 -1
	if(isfluid(i+1,j-1,k).eq.0)then
	  udotc=(-hfieldsh(i,j,k,2,idblock)+hfieldsh(i,j,k,3,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,idblock)+qyy*hfieldsh(i,j,k,6,idblock) &
	   -cssq*hfieldsh(i,j,k,7,idblock)+two*qxy_9_10*hfieldsh(i,j,k,8,idblock))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i+1,j-1,k,2,idblock)*hfieldsh(i+1,j-1,k,2,idblock) &
	   + hfieldsh(i+1,j-1,k,3,idblock)*hfieldsh(i+1,j-1,k,3,idblock) + hfieldsh(i+1,j-1,k,4,idblock)*hfieldsh(i+1,j-1,k,4,idblock))
	  udotc=(-hfieldsh(i+1,j-1,k,2,idblock)+hfieldsh(i+1,j-1,k,3,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i+1,j-1,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i+1,j-1,k,5,idblock)+qyy*hfieldsh(i+1,j-1,k,6,idblock) &
	   -cssq*hfieldsh(i+1,j-1,k,7,idblock)+two*qxy_9_10*hfieldsh(i+1,j-1,k,8,idblock))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	endif
	!-1 +1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_v=temp_v+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop
	
	!9  -1 +1
	if(isfluid(i-1,j+1,k).eq.0)then
	  udotc=(-hfieldsh(i,j,k,2,idblock)+hfieldsh(i,j,k,3,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,idblock)+qyy*hfieldsh(i,j,k,6,idblock) &
	   -cssq*hfieldsh(i,j,k,7,idblock)+two*qxy_9_10*hfieldsh(i,j,k,8,idblock))
	  temp_pop=temp_pop +(fy-fx)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i-1,j+1,k,2,idblock)*hfieldsh(i-1,j+1,k,2,idblock) &
	   + hfieldsh(i-1,j+1,k,3,idblock)*hfieldsh(i-1,j+1,k,3,idblock) + hfieldsh(i-1,j+1,k,4,idblock)*hfieldsh(i-1,j+1,k,4,idblock))
	  udotc=(-hfieldsh(i-1,j+1,k,2,idblock)+hfieldsh(i-1,j+1,k,3,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i-1,j+1,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i-1,j+1,k,5,idblock)+qyy*hfieldsh(i-1,j+1,k,6,idblock) &
	   -cssq*hfieldsh(i-1,j+1,k,7,idblock)+two*qxy_9_10*hfieldsh(i-1,j+1,k,8,idblock))
	  temp_pop=temp_pop + (fx-fy)*p2dcssq
	endif
	!+1 -1  0
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_v=temp_v-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pxy=temp_pxy-temp_pop

	!15  -1  -1
	if(isfluid(i-1,j,k-1).eq.0)then
	  udotc=(hfieldsh(i,j,k,2,idblock)+hfieldsh(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,idblock)+qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*hfieldsh(i,j,k,6,idblock)+two*qxz_15_16*hfieldsh(i,j,k,9,idblock))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i-1,j,k-1,2,idblock)*hfieldsh(i-1,j,k-1,2,idblock) &
	   + hfieldsh(i-1,j,k-1,3,idblock)*hfieldsh(i-1,j,k-1,3,idblock) + hfieldsh(i-1,j,k-1,4,idblock)*hfieldsh(i-1,j,k-1,4,idblock))
	  udotc=(hfieldsh(i-1,j,k-1,2,idblock)+hfieldsh(i-1,j,k-1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i-1,j,k-1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i-1,j,k-1,5,idblock)+qzz*hfieldsh(i-1,j,k-1,7,idblock) &
	   -cssq*hfieldsh(i-1,j,k-1,6,idblock)+two*qxz_15_16*hfieldsh(i-1,j,k-1,9,idblock))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq 
	endif
	!+1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pop
	
	!16  +1  +1
	if(isfluid(i+1,j,k+1).eq.0)then
	  udotc=(hfieldsh(i,j,k,2,idblock)+hfieldsh(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,idblock)+qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*hfieldsh(i,j,k,6,idblock)+two*qxz_15_16*hfieldsh(i,j,k,9,idblock))
	  temp_pop=temp_pop + (fx+fz)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i+1,j,k+1,2,idblock)*hfieldsh(i+1,j,k+1,2,idblock) &
	   + hfieldsh(i+1,j,k+1,3,idblock)*hfieldsh(i+1,j,k+1,3,idblock) + hfieldsh(i+1,j,k+1,4,idblock)*hfieldsh(i+1,j,k+1,4,idblock))
	  udotc=(hfieldsh(i+1,j,k+1,2,idblock)+hfieldsh(i+1,j,k+1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i+1,j,k+1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i+1,j,k+1,5,idblock)+qzz*hfieldsh(i+1,j,k+1,7,idblock) &
	   -cssq*hfieldsh(i+1,j,k+1,6,idblock)+two*qxz_15_16*hfieldsh(i+1,j,k+1,9,idblock))
	  temp_pop=temp_pop - (fx+fz)*p2dcssq
	endif
	!-1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz+temp_pop

	!17  +1   -1
	if(isfluid(i+1,j,k-1).eq.0)then
	  udotc=(-hfieldsh(i,j,k,2,idblock)+hfieldsh(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,idblock)+qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*hfieldsh(i,j,k,6,idblock)+two*qxz_17_18*hfieldsh(i,j,k,9,idblock))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i+1,j,k-1,2,idblock)*hfieldsh(i+1,j,k-1,2,idblock) &
	   + hfieldsh(i+1,j,k-1,3,idblock)*hfieldsh(i+1,j,k-1,3,idblock) + hfieldsh(i+1,j,k-1,4,idblock)*hfieldsh(i+1,j,k-1,4,idblock))
	  udotc=(-hfieldsh(i+1,j,k-1,2,idblock)+hfieldsh(i+1,j,k-1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i+1,j,k-1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i+1,j,k-1,5,idblock)+qzz*hfieldsh(i+1,j,k-1,7,idblock) &
	   -cssq*hfieldsh(i+1,j,k-1,6,idblock)+two*qxz_17_18*hfieldsh(i+1,j,k-1,9,idblock))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	endif
	!-1  0  +1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u-temp_pop
	temp_w=temp_w+temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop
	
	!18   -1   +1
	if(isfluid(i-1,j,k+1).eq.0)then
	  udotc=(-hfieldsh(i,j,k,2,idblock)+hfieldsh(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,idblock)+qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*hfieldsh(i,j,k,6,idblock)+two*qxz_17_18*hfieldsh(i,j,k,9,idblock))
	  temp_pop=temp_pop +(fz-fx)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i-1,j,k+1,2,idblock)*hfieldsh(i-1,j,k+1,2,idblock) &
	   + hfieldsh(i-1,j,k+1,3,idblock)*hfieldsh(i-1,j,k+1,3,idblock) + hfieldsh(i-1,j,k+1,4,idblock)*hfieldsh(i-1,j,k+1,4,idblock))
	  udotc=(-hfieldsh(i-1,j,k+1,2,idblock)+hfieldsh(i-1,j,k+1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i-1,j,k+1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qxx*hfieldsh(i-1,j,k+1,5,idblock)+qzz*hfieldsh(i-1,j,k+1,7,idblock) &
	   -cssq*hfieldsh(i-1,j,k+1,6,idblock)+two*qxz_17_18*hfieldsh(i-1,j,k+1,9,idblock))
	  temp_pop=temp_pop + (fx-fz)*p2dcssq
	endif
	!+1  0  -1
	temp_rho=temp_rho+temp_pop
	temp_u=temp_u+temp_pop
	temp_w=temp_w-temp_pop
	temp_pxx=temp_pxx+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pxz=temp_pxz-temp_pop

	!11  -1  -1
	if(isfluid(i,j-1,k-1).eq.0)then
	  udotc=(hfieldsh(i,j,k,3,idblock)+hfieldsh(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j,k,6,idblock)+qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*hfieldsh(i,j,k,5,idblock)+two*qyz_11_12*hfieldsh(i,j,k,10,idblock))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j-1,k-1,2,idblock)*hfieldsh(i,j-1,k-1,2,idblock) &
	   + hfieldsh(i,j-1,k-1,3,idblock)*hfieldsh(i,j-1,k-1,3,idblock) + hfieldsh(i,j-1,k-1,4,idblock)*hfieldsh(i,j-1,k-1,4,idblock))
	  udotc=(hfieldsh(i,j-1,k-1,3,idblock)+hfieldsh(i,j-1,k-1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j-1,k-1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j-1,k-1,6,idblock)+qzz*hfieldsh(i,j-1,k-1,7,idblock) &
	   -cssq*hfieldsh(i,j-1,k-1,5,idblock)+two*qyz_11_12*hfieldsh(i,j-1,k-1,10,idblock))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	endif
	! 0 +1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pop
	
	!12   +1  +1
	if(isfluid(i,j+1,k+1).eq.0)then
	  udotc=(hfieldsh(i,j,k,3,idblock)+hfieldsh(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j,k,6,idblock)+qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*hfieldsh(i,j,k,5,idblock)+two*qyz_11_12*hfieldsh(i,j,k,10,idblock))
	  temp_pop=temp_pop + (fy+fz)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j+1,k+1,2,idblock)*hfieldsh(i,j+1,k+1,2,idblock) &
	   + hfieldsh(i,j+1,k+1,3,idblock)*hfieldsh(i,j+1,k+1,3,idblock) + hfieldsh(i,j+1,k+1,4,idblock)*hfieldsh(i,j+1,k+1,4,idblock))
	  udotc=(hfieldsh(i,j+1,k+1,3,idblock)+hfieldsh(i,j+1,k+1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j+1,k+1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j+1,k+1,6,idblock)+qzz*hfieldsh(i,j+1,k+1,7,idblock) &
	   -cssq*hfieldsh(i,j+1,k+1,5,idblock)+two*qyz_11_12*hfieldsh(i,j+1,k+1,10,idblock))
	  temp_pop=temp_pop - (fy+fz)*p2dcssq
	endif
	! 0 -1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz+temp_pop

	!13   -1   +1
	if(isfluid(i,j-1,k+1).eq.0)then
	  udotc=(hfieldsh(i,j,k,3,idblock)-hfieldsh(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j,k,6,idblock)+qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*hfieldsh(i,j,k,5,idblock)+two*qyz_13_14*hfieldsh(i,j,k,10,idblock))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j-1,k+1,2,idblock)*hfieldsh(i,j-1,k+1,2,idblock) &
	   + hfieldsh(i,j-1,k+1,3,idblock)*hfieldsh(i,j-1,k+1,3,idblock) + hfieldsh(i,j-1,k+1,4,idblock)*hfieldsh(i,j-1,k+1,4,idblock))
	  udotc=(hfieldsh(i,j-1,k+1,3,idblock)-hfieldsh(i,j-1,k+1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j-1,k+1,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j-1,k+1,6,idblock)+qzz*hfieldsh(i,j-1,k+1,7,idblock) &
	   -cssq*hfieldsh(i,j-1,k+1,5,idblock)+two*qyz_13_14*hfieldsh(i,j-1,k+1,10,idblock))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	endif
	! 0 +1 -1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v+temp_pop
	temp_w=temp_w-temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	!14  +1 -1
	if(isfluid(i,j+1,k-1).eq.0)then
	  udotc=(hfieldsh(i,j,k,3,idblock)-hfieldsh(i,j,k,4,idblock))*onecssq
	  temp = -uu0 + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j,k,1,idblock)+(temp + udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j,k,6,idblock)+qzz*hfieldsh(i,j,k,7,idblock) &
	   -cssq*hfieldsh(i,j,k,5,idblock)+two*qyz_13_14*hfieldsh(i,j,k,10,idblock))
	  temp_pop=temp_pop + (fy-fz)*p2dcssq
	else
	  uu=halfonecssq*(hfieldsh(i,j+1,k-1,2,idblock)*hfieldsh(i,j+1,k-1,2,idblock) &
	   + hfieldsh(i,j+1,k-1,3,idblock)*hfieldsh(i,j+1,k-1,3,idblock) + hfieldsh(i,j+1,k-1,4,idblock)*hfieldsh(i,j+1,k-1,4,idblock))
	  udotc=(hfieldsh(i,j+1,k-1,3,idblock)-hfieldsh(i,j+1,k-1,4,idblock))*onecssq
	  temp = -uu + half*udotc*udotc
	  temp_pop=p2*(hfieldsh(i,j+1,k-1,1,idblock)+(temp - udotc))
	  temp_pop=temp_pop+oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j+1,k-1,6,idblock)+qzz*hfieldsh(i,j+1,k-1,7,idblock) &
	   -cssq*hfieldsh(i,j+1,k-1,5,idblock)+two*qyz_13_14*hfieldsh(i,j+1,k-1,10,idblock))
	  temp_pop=temp_pop + (fz-fy)*p2dcssq
	endif
	! 0 -1 +1
	temp_rho=temp_rho+temp_pop
	temp_v=temp_v-temp_pop
	temp_w=temp_w+temp_pop
	temp_pyy=temp_pyy+temp_pop
	temp_pzz=temp_pzz+temp_pop
	temp_pyz=temp_pyz-temp_pop
	
	hfields(i,j,k,1,idblock)=temp_rho
	
	hfields(i,j,k,2,idblock)=temp_u
	hfields(i,j,k,3,idblock)=temp_v
	hfields(i,j,k,4,idblock)=temp_w
	
	hfields(i,j,k,4,idblock)=temp_pxx
	hfields(i,j,k,6,idblock)=temp_pyy
	hfields(i,j,k,7,idblock)=temp_pzz
	hfields(i,j,k,8,idblock)=temp_pxy
	hfields(i,j,k,9,idblock)=temp_pxz
	hfields(i,j,k,10,idblock)=temp_pyz
	
	return
  
  end subroutine streamcoll_bc_flop
  
  attributes(global) subroutine streamcoll_bc_shared()
	
	implicit none  
	  
    integer :: i,j,k,gi,gj,gk,myblock,iidblock
    integer :: gii,gjj,gkk,ii,jj,kk
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
    integer :: li,lj,lk
	real(kind=db), shared :: f00(1:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f01(0:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f02(1:TILE_DIMx_d+1,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f03(1:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f04(1:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f05(1:TILE_DIMx_d,1:TILE_DIMy_d,0:TILE_DIMz_d)
	real(kind=db), shared :: f06(1:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
	real(kind=db), shared :: f07(0:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f08(1:TILE_DIMx_d+1,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f09(0:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f10(1:TILE_DIMx_d+1,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f11(1:TILE_DIMx_d,0:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f12(1:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f13(1:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f14(1:TILE_DIMx_d,1:TILE_DIMy_d+1,0:TILE_DIMz_d)
    real(kind=db), shared :: f15(0:TILE_DIMx_d,1:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f16(1:TILE_DIMx_d+1,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f17(1:TILE_DIMx_d+1,1:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f18(0:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
    
		  
		  
	!if(isfluid(i,j,k).ne.1)return
    	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    
    
    uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
     + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(hfields(i,j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*( &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock) &
	 -cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock) &
	 -cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock) &
	 -cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock) &
	 -cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    ! Halo Faces
    if(li==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !1 -1  0  0
	  udotc=hfields(ii,jj,kk,2,iidblock)*onecssq
	  f01(li-1,lj,lk)=p1*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*hfields(ii,jj,kk,5,iidblock) &
	   -cssq*(hfields(ii,jj,kk,6,iidblock)+hfields(ii,jj,kk,7,iidblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f07(li-1,lj,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfields(ii,jj,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f09(li-1,lj,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfields(ii,jj,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f15(li-1,lj,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfields(ii,jj,kk,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f18(li-1,lj,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfields(ii,jj,kk,9,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !2 +1  0  0
	  udotc=hfields(ii,jj,kk,2,iidblock)*onecssq
	  f02(li+1,lj,lk)=p1*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*hfields(ii,jj,kk,5,iidblock) &
	   -cssq*(hfields(ii,jj,kk,6,iidblock)+hfields(ii,jj,kk,7,iidblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f08(li+1,lj,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfields(ii,jj,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f10(li+1,lj,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfields(ii,jj,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f16(li+1,lj,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfields(ii,jj,kk,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f17(li+1,lj,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfields(ii,jj,kk,9,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !3 0 -1  0
	  udotc=hfields(ii,jj,kk,3,iidblock)*onecssq
	  f03(li,lj-1,lk)=p1*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*(hfields(ii,jj,kk,5,iidblock)+hfields(ii,jj,kk,7,iidblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f07(li,lj-1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfields(ii,jj,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f10(li,lj-1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfields(ii,jj,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(hfields(ii,jj,kk,3,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f11(li,lj-1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfields(ii,jj,kk,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(hfields(ii,jj,kk,3,iidblock)-hfields(ii,jj,kk,4,iidblock))*onecssq
	  f13(li,lj-1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfields(ii,jj,kk,10,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !4  0 +1  0
	  udotc=hfields(ii,jj,kk,3,iidblock)*onecssq
	  f04(li,lj+1,lk)=p1*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*(hfields(ii,jj,kk,5,iidblock)+hfields(ii,jj,kk,7,iidblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f08(li,lj+1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfields(ii,jj,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f09(li,lj+1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfields(ii,jj,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(hfields(ii,jj,kk,3,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f12(li,lj+1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfields(ii,jj,kk,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(hfields(ii,jj,kk,3,iidblock)-hfields(ii,jj,kk,4,iidblock))*onecssq
	  f14(li,lj+1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfields(ii,jj,kk,10,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      gii=gi
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !5  0  0 -1
	  udotc=hfields(ii,jj,kk,4,iidblock)*onecssq
	  f05(li,lj,lk-1)=p1*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*(hfields(ii,jj,kk,5,iidblock)+hfields(ii,jj,kk,6,iidblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f15(li,lj,lk-1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfields(ii,jj,kk,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f17(li,lj,lk-1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfields(ii,jj,kk,9,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(hfields(ii,jj,kk,3,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f11(li,lj,lk-1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfields(ii,jj,kk,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(hfields(ii,jj,kk,3,iidblock)-hfields(ii,jj,kk,4,iidblock))*onecssq
	  f14(li,lj,lk-1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfields(ii,jj,kk,10,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !6  0  0  +1
	  udotc=hfields(ii,jj,kk,4,iidblock)*onecssq
	  f06(li,lj,lk+1)=p1*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*(hfields(ii,jj,kk,5,iidblock)+hfields(ii,jj,kk,6,iidblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f16(li,lj,lk+1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfields(ii,jj,kk,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f18(li,lj,lk+1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfields(ii,jj,kk,9,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(hfields(ii,jj,kk,3,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f12(li,lj,lk+1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfields(ii,jj,kk,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(hfields(ii,jj,kk,3,iidblock)-hfields(ii,jj,kk,4,iidblock))*onecssq
	  f13(li,lj,lk+1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfields(ii,jj,kk,10,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      gii=gi-1
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !7 -1 -1  0
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f07(li-1,lj-1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfields(ii,jj,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      gii=gi-1
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !9  -1 +1 0
      udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f09(li-1,lj+1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfields(ii,jj,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !15  -1  0 -1
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f15(li-1,lj,lk-1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfields(ii,jj,kk,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi-1
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !18   -1   0  +1
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f18(li-1,lj,lk+1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfields(ii,jj,kk,9,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      gii=gi+1
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !10   +1 -1  0
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f10(li+1,lj-1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfields(ii,jj,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      gii=gi+1
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !8 +1 +1  0
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,3,iidblock))*onecssq
	  f08(li+1,lj+1,lk)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qyy*hfields(ii,jj,kk,6,iidblock) &
	   -cssq*hfields(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfields(ii,jj,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      gii=gi+1
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !17  +1  0 -1
	  udotc=(-hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f17(li+1,lj,lk-1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock)&
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfields(ii,jj,kk,9,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !16  +1  0 +1
	  udotc=(hfields(ii,jj,kk,2,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f16(li+1,lj,lk+1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfields(ii,jj,kk,5,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfields(ii,jj,kk,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !11  0  -1  -1
	  udotc=(hfields(ii,jj,kk,3,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f11(li,lj-1,lk-1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfields(ii,jj,kk,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj-1
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !13  0  -1   +1
	  udotc=(hfields(ii,jj,kk,3,iidblock)-hfields(ii,jj,kk,4,iidblock))*onecssq
	  f13(li,lj-1,lk+1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfields(ii,jj,kk,10,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      gii=gi
      gjj=gj+1
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !14  0  +1  -1
	  udotc=(hfields(ii,jj,kk,3,iidblock)-hfields(ii,jj,kk,4,iidblock))*onecssq
	  f14(li,lj+1,lk-1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfields(ii,jj,kk,10,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfields(ii,jj,kk,2,iidblock)*hfields(ii,jj,kk,2,iidblock) &
      + hfields(ii,jj,kk,3,iidblock)*hfields(ii,jj,kk,3,iidblock) + hfields(ii,jj,kk,4,iidblock)*hfields(ii,jj,kk,4,iidblock))
      
      !12  0  +1  +1
	  udotc=(hfields(ii,jj,kk,3,iidblock)+hfields(ii,jj,kk,4,iidblock))*onecssq
	  f12(li,lj+1,lk+1)=p2*(hfields(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfields(ii,jj,kk,6,iidblock)+qzz*hfields(ii,jj,kk,7,iidblock) &
	   -cssq*hfields(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfields(ii,jj,kk,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
    
    if(isfluid(gi,gj,gk).eq.0)then
#ifdef IFBC
          
          if(isfluid(gi+1,gj,gk-1).eq.-1)f18(li,lj,lk)=f17(li+1,lj,lk-1) !gpc 
          if(isfluid(gi-1,gj,gk+1).eq.-1)f17(li,lj,lk)=f18(li-1,lj,lk+1) !hpc

          if(isfluid(gi-1,gj,gk-1).eq.-1)f16(li,lj,lk)=f15(li-1,lj,lk-1) !gpc 
          if(isfluid(gi+1,gj,gk+1).eq.-1)f15(li,lj,lk)=f16(li+1,lj,lk+1) !hpc

          if(isfluid(gi,gj-1,gk+1).eq.-1)f14(li,lj,lk)=f13(li,lj-1,lk+1)!gpc 
          if(isfluid(gi,gj+1,gk-1).eq.-1)f13(li,lj,lk)=f14(li,lj+1,lk-1)!hpc
        
          if(isfluid(gi,gj-1,gk-1).eq.-1)f12(li,lj,lk)=f11(li,lj-1,lk-1)!gpc
          if(isfluid(gi,gj+1,gk+1).eq.-1)f11(li,lj,lk)=f12(li,lj+1,lk+1)!hpc

          if(isfluid(gi-1,gj+1,gk).eq.-1)f10(li,lj,lk)=f09(li-1,lj+1,lk)!gpc  
          if(isfluid(gi+1,gj-1,gk).eq.-1)f09(li,lj,lk)=f10(li+1,lj-1,lk)!hpc

          if(isfluid(gi-1,gj-1,gk).eq.-1)f08(li,lj,lk)=f07(li-1,lj-1,lk)!gpc 
          if(isfluid(gi+1,gj+1,gk).eq.-1)f07(li,lj,lk)=f08(li+1,lj+1,lk)!hpc

          if(isfluid(gi,gj,gk-1).eq.-1)f06(li,lj,lk)=f05(li,lj,lk-1)!gpc
          if(isfluid(gi,gj,gk+1).eq.-1)f05(li,lj,lk)=f06(li,lj,lk+1)!hpc 

          if(isfluid(gi,gj-1,gk).eq.-1)f04(li,lj,lk)=f03(li,lj-1,lk)!gpc 
          if(isfluid(gi,gj+1,gk).eq.-1)f03(li,lj,lk)=f04(li,lj+1,lk)!hpc 

          if(isfluid(gi-1,gj,gk).eq.-1)f02(li,lj,lk)=f01(li-1,lj,lk)!gpc
          if(isfluid(gi+1,gj,gk).eq.-1)f01(li,lj,lk)=f02(li+1,lj,lk)!hpc 
          
#else
          
        f01(li,lj,lk)=f02(li+1,lj,lk)!hpc 
        f02(li,lj,lk)=f01(li-1,lj,lk)!gpc 
        
        f03(li,lj,lk)=f04(li,lj+1,lk)!hpc 
        f04(li,lj,lk)=f03(li,lj-1,lk)!gpc 
        
        f05(li,lj,lk)=f06(li,lj,lk+1)!hpc 
        f06(li,lj,lk)=f05(li,lj,lk-1)!gpc 
        
        f07(li,lj,lk)=f08(li+1,lj+1,lk)!hpc
        f08(li,lj,lk)=f07(li-1,lj-1,lk)!gpc 
        
        f09(li,lj,lk)=f10(li+1,lj-1,lk)!hpc
        f10(li,lj,lk)=f09(li-1,lj+1,lk)!gpc 
        
        f11(li,lj,lk)=f12(li,lj+1,lk+1)!hpc
        f12(li,lj,lk)=f11(li,lj-1,lk-1)!gpc 
        
        f13(li,lj,lk)=f14(li,lj+1,lk-1)!hpc
        f14(li,lj,lk)=f13(li,lj-1,lk+1)!gpc 
        
        f15(li,lj,lk)=f16(li+1,lj,lk+1) !hpc
        f16(li,lj,lk)=f15(li-1,lj,lk-1) !gpc 
        
        f17(li,lj,lk)=f18(li-1,lj,lk+1) !hpc
        f18(li,lj,lk)=f17(li+1,lj,lk-1) !gpc       

#endif
    endif
    
    ! Halo Faces
    if(li==1)then
      if(isfluid(gi-1,gj,gk).eq.0)then
        f01(li-1,lj,lk)=f02(li,lj,lk)!hpc 
        
        f07(li-1,lj,lk)=f08(li,lj+1,lk)!hpc
        
        f09(li-1,lj,lk)=f10(li,lj-1,lk)!hpc
        
        f15(li-1,lj,lk)=f16(li,lj,lk+1) !hpc

        f18(li-1,lj,lk)=f17(li,lj,lk-1) !gpc 
      endif
    endif
    if(li==TILE_DIMx_d)then
      if(isfluid(gi+1,gj,gk).eq.0)then
        f02(li+1,lj,lk)=f01(li,lj,lk)!gpc 
        
        f08(li+1,lj,lk)=f07(li,lj-1,lk)!gpc 
        
        f10(li+1,lj,lk)=f09(li,lj+1,lk)!gpc 
        
        f16(li+1,lj,lk)=f15(li,lj,lk-1) !gpc 
        
        f17(li+1,lj,lk)=f18(li,lj,lk+1) !hpc 
      endif
    endif

    if(lj==1)then
      if(isfluid(gi,gj-1,gk).eq.0)then
        f03(li,lj-1,lk)=f04(li,lj,lk)!hpc 
        
        f07(li,lj-1,lk)=f08(li+1,lj,lk)!hpc
        
        f10(li,lj-1,lk)=f09(li-1,lj,lk)!gpc 
        
        f11(li,lj-1,lk)=f12(li,lj,lk+1)!hpc
        
        f13(li,lj-1,lk)=f14(li,lj,lk-1)!hpc    
      endif
    endif
    if(lj==TILE_DIMy_d)then
      if(isfluid(gi,gj+1,gk).eq.0)then
        f04(li,lj+1,lk)=f03(li,lj,lk)!gpc 
        
        f08(li,lj+1,lk)=f07(li-1,lj,lk)!gpc 
        
        f09(li,lj+1,lk)=f10(li+1,lj,lk)!hpc
        
        f12(li,lj+1,lk)=f11(li,lj,lk-1)!gpc 
        
        f14(li,lj+1,lk)=f13(li,lj,lk+1)!gpc   
      endif
    endif

    if(lk==1)then
      if(isfluid(gi,gj,gk-1).eq.0)then
        f05(li,lj,lk-1)=f06(li,lj,lk)!hpc 
        
        f11(li,lj,lk-1)=f12(li,lj+1,lk)!hpc
        
        f14(li,lj,lk-1)=f13(li,lj-1,lk)!gpc 
        
        f15(li,lj,lk-1)=f16(li+1,lj,lk) !hpc
        
        f17(li,lj,lk-1)=f18(li-1,lj,lk) !hpc
      endif
    endif
    if(lk==TILE_DIMz_d)then
      if(isfluid(gi,gj,gk+1).eq.0)then
        f06(li,lj,lk+1)=f05(li,lj,lk)!gpc 
        
        f12(li,lj,lk+1)=f11(li,lj-1,lk)!gpc 
        
        f13(li,lj,lk+1)=f14(li,lj+1,lk)!hpc
        
        f16(li,lj,lk+1)=f15(li-1,lj,lk) !gpc 
        
        f18(li,lj,lk+1)=f17(li+1,lj,lk) !gpc 
      endif
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      if(isfluid(gi-1,gj-1,gk).eq.0)then
        f07(li-1,lj-1,lk)=f08(li,lj,lk)!hpc
      endif
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      if(isfluid(gi-1,gj+1,gk).eq.0)then
        f09(li-1,lj+1,lk)=f10(li,lj,lk)!hpc
      endif
    endif
    if(li==1 .and. lk==1)then
      if(isfluid(gi-1,gj,gk-1).eq.0)then
        f15(li-1,lj,lk-1)=f16(li,lj,lk) !hpc
      endif
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      if(isfluid(gi-1,gj,gk+1).eq.0)then
        f18(li-1,lj,lk+1)=f17(li,lj,lk) !gpc    
      endif
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      if(isfluid(gi+1,gj-1,gk).eq.0)then
        f10(li+1,lj-1,lk)=f09(li,lj,lk)!gpc 
      endif
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      if(isfluid(gi+1,gj+1,gk).eq.0)then
        f08(li+1,lj+1,lk)=f07(li,lj,lk)!gpc 
      endif
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      if(isfluid(gi+1,gj,gk-1).eq.0)then
        f17(li+1,lj,lk-1)=f18(li,lj,lk) !hpc
      endif
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      if(isfluid(gi+1,gj,gk+1).eq.0)then
        f16(li+1,lj,lk+1)=f15(li,lj,lk) !gpc 
      endif
    endif
    if(lj==1 .and. lk==1)then
      if(isfluid(gi,gj-1,gk-1).eq.0)then
        f11(li,lj-1,lk-1)=f12(li,lj,lk)!hpc
      endif
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      if(isfluid(gi,gj-1,gk+1).eq.0)then
        f13(li,lj-1,lk+1)=f14(li,lj,lk)!hpc
      endif
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      if(isfluid(gi,gj+1,gk-1).eq.0)then
        f14(li,lj+1,lk-1)=f13(li,lj,lk)!gpc 
      endif
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      if(isfluid(gi,gj+1,gk+1).eq.0)then
        f12(li,lj+1,lk+1)=f11(li,lj,lk)!gpc 
      endif
    endif      
    
    call syncthreads
    
    udotc=f00(li,lj,lk)+f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(i,j,k,1,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(i,j,k,2,myblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfieldsh(i,j,k,3,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfieldsh(i,j,k,4,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(i,j,k,5,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfieldsh(i,j,k,6,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfieldsh(i,j,k,7,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfieldsh(i,j,k,8,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfieldsh(i,j,k,9,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfieldsh(i,j,k,10,myblock)=udotc
     
#ifdef PRESSCORR

	call syncthreads
	
	uu=halfonecssq*(hfieldsh(i,j,k,2,myblock)*hfieldsh(i,j,k,2,myblock) &
	 + hfieldsh(i,j,k,3,myblock)*hfieldsh(i,j,k,3,myblock) + hfieldsh(i,j,k,4,myblock)*hfieldsh(i,j,k,4,myblock))
    
	!1 -1  0  0
	udotc=hfieldsh(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfieldsh(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfieldsh(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfieldsh(i,j,k,2,myblock)+hfieldsh(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfieldsh(i,j,k,2,myblock)+hfieldsh(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfieldsh(i,j,k,2,myblock)+hfieldsh(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfieldsh(i,j,k,2,myblock)+hfieldsh(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfieldsh(i,j,k,3,myblock)+hfieldsh(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfieldsh(i,j,k,3,myblock)-hfieldsh(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfieldsh(i,j,k,5,myblock)=hfieldsh(i,j,k,5,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfieldsh(i,j,k,6,myblock)=hfieldsh(i,j,k,6,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfieldsh(i,j,k,7,myblock)=hfieldsh(i,j,k,7,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfieldsh(i,j,k,8,myblock)=hfieldsh(i,j,k,8,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfieldsh(i,j,k,9,myblock)=hfieldsh(i,j,k,9,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfieldsh(i,j,k,10,myblock)=hfieldsh(i,j,k,10,myblock)-udotc
	
	    
    return
#endif	
  end subroutine streamcoll_bc_shared
  
  attributes(global) subroutine streamcoll_bc_shared_flop()
	
	implicit none  
	
	
	integer :: i,j,k,gi,gj,gk,myblock,iidblock
    integer :: gii,gjj,gkk,ii,jj,kk
    integer :: xblock,yblock,zblock
	real(kind=db) :: udotc,uu,temp
    
    
    integer :: li,lj,lk
	real(kind=db), shared :: f00(1:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f01(0:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f02(1:TILE_DIMx_d+1,1:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f03(1:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f04(1:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f05(1:TILE_DIMx_d,1:TILE_DIMy_d,0:TILE_DIMz_d)
	real(kind=db), shared :: f06(1:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
	real(kind=db), shared :: f07(0:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f08(1:TILE_DIMx_d+1,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f09(0:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d)
	real(kind=db), shared :: f10(1:TILE_DIMx_d+1,0:TILE_DIMy_d,1:TILE_DIMz_d)
	real(kind=db), shared :: f11(1:TILE_DIMx_d,0:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f12(1:TILE_DIMx_d,1:TILE_DIMy_d+1,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f13(1:TILE_DIMx_d,0:TILE_DIMy_d,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f14(1:TILE_DIMx_d,1:TILE_DIMy_d+1,0:TILE_DIMz_d)
    real(kind=db), shared :: f15(0:TILE_DIMx_d,1:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f16(1:TILE_DIMx_d+1,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
    real(kind=db), shared :: f17(1:TILE_DIMx_d+1,1:TILE_DIMy_d,0:TILE_DIMz_d)
    real(kind=db), shared :: f18(0:TILE_DIMx_d,1:TILE_DIMy_d,1:TILE_DIMz_d+1)
    
		  
		  
	!if(isfluid(i,j,k).ne.1)return
    	
	li = threadIdx%x
    lj = threadIdx%y
    lk = threadIdx%z
    
    gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	i=threadIdx%x
	j=threadIdx%y
	k=threadIdx%z
	
	myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
    
    
    
    uu=halfonecssq*(hfieldsh(i,j,k,2,myblock)*hfieldsh(i,j,k,2,myblock) &
     + hfieldsh(i,j,k,3,myblock)*hfieldsh(i,j,k,3,myblock) + hfieldsh(i,j,k,4,myblock)*hfieldsh(i,j,k,4,myblock))
    
    !0
	f00(li,lj,lk)=p0*(hfieldsh(i,j,k,1,myblock)-uu) &
	 + oneminusomega*pi2cssq0*( &
	 -cssq*(hfieldsh(i,j,k,6,myblock)+hfieldsh(i,j,k,5,myblock)+hfieldsh(i,j,k,7,myblock)))
	
    
	!1 -1  0  0
	udotc=hfieldsh(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfieldsh(i,j,k,5,myblock) &
	 -cssq*(hfieldsh(i,j,k,6,myblock)+hfieldsh(i,j,k,7,myblock))) &
	 + fx*p1dcssq
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qxx*hfieldsh(i,j,k,5,myblock) &
	 -cssq*(hfieldsh(i,j,k,6,myblock)+hfieldsh(i,j,k,7,myblock))) &
	 - fx*p1dcssq
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfieldsh(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfieldsh(i,j,k,6,myblock) &
	 -cssq*(hfieldsh(i,j,k,5,myblock)+hfieldsh(i,j,k,7,myblock))) &
	 + fy*p1dcssq
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qyy*hfieldsh(i,j,k,6,myblock) &
	 -cssq*(hfieldsh(i,j,k,5,myblock)+hfieldsh(i,j,k,7,myblock))) &
	 - fy*p1dcssq
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfieldsh(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*(hfieldsh(i,j,k,5,myblock)+hfieldsh(i,j,k,6,myblock))) &
	 + fz*p1dcssq
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq1*(qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*(hfieldsh(i,j,k,5,myblock)+hfieldsh(i,j,k,6,myblock))) &
	 - fz*p1dcssq
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfieldsh(i,j,k,2,myblock)+hfieldsh(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,myblock)+qyy*hfieldsh(i,j,k,6,myblock) &
	 -cssq*hfieldsh(i,j,k,7,myblock)+two*qxy_7_8*hfieldsh(i,j,k,8,myblock)) &
	 + (fx+fy)*p2dcssq 
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,myblock)+qyy*hfieldsh(i,j,k,6,myblock) &
	 -cssq*hfieldsh(i,j,k,7,myblock)+two*qxy_7_8*hfieldsh(i,j,k,8,myblock)) &
	 - (fx+fy)*p2dcssq
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfieldsh(i,j,k,2,myblock)+hfieldsh(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,myblock)+qyy*hfieldsh(i,j,k,6,myblock) &
	 -cssq*hfieldsh(i,j,k,7,myblock)+two*qxy_9_10*hfieldsh(i,j,k,8,myblock)) &
	 +(fy-fx)*p2dcssq
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,myblock)+qyy*hfieldsh(i,j,k,6,myblock) &
	 -cssq*hfieldsh(i,j,k,7,myblock)+two*qxy_9_10*hfieldsh(i,j,k,8,myblock)) &
	 + (fx-fy)*p2dcssq
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfieldsh(i,j,k,2,myblock)+hfieldsh(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,myblock)+qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*hfieldsh(i,j,k,6,myblock)+two*qxz_15_16*hfieldsh(i,j,k,9,myblock)) &
	 + (fx+fz)*p2dcssq 
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,myblock)+qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*hfieldsh(i,j,k,6,myblock)+two*qxz_15_16*hfieldsh(i,j,k,9,myblock)) &
	 - (fx+fz)*p2dcssq
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfieldsh(i,j,k,2,myblock)+hfieldsh(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,myblock)+qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*hfieldsh(i,j,k,6,myblock)+two*qxz_17_18*hfieldsh(i,j,k,9,myblock)) &
	 +(fz-fx)*p2dcssq
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qxx*hfieldsh(i,j,k,5,myblock)+qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*hfieldsh(i,j,k,6,myblock)+two*qxz_17_18*hfieldsh(i,j,k,9,myblock)) &
	 + (fx-fz)*p2dcssq
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfieldsh(i,j,k,3,myblock)+hfieldsh(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j,k,6,myblock)+qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*hfieldsh(i,j,k,5,myblock)+two*qyz_11_12*hfieldsh(i,j,k,10,myblock)) &
	 + (fy+fz)*p2dcssq
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j,k,6,myblock)+qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*hfieldsh(i,j,k,5,myblock)+two*qyz_11_12*hfieldsh(i,j,k,10,myblock)) &
	 - (fy+fz)*p2dcssq
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfieldsh(i,j,k,3,myblock)-hfieldsh(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp + udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j,k,6,myblock)+qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*hfieldsh(i,j,k,5,myblock)+two*qyz_13_14*hfieldsh(i,j,k,10,myblock)) &
	 + (fy-fz)*p2dcssq
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfieldsh(i,j,k,1,myblock)+(temp - udotc)) &
	 +oneminusomega*pi2cssq2*(qyy*hfieldsh(i,j,k,6,myblock)+qzz*hfieldsh(i,j,k,7,myblock) &
	 -cssq*hfieldsh(i,j,k,5,myblock)+two*qyz_13_14*hfieldsh(i,j,k,10,myblock)) &
	 + (fz-fy)*p2dcssq
	! 0 -1 +1
    
    ! Halo Faces
    if(li==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !1 -1  0  0
	  udotc=hfieldsh(ii,jj,kk,2,iidblock)*onecssq
	  f01(li-1,lj,lk)=p1*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*hfieldsh(ii,jj,kk,5,iidblock) &
	   -cssq*(hfieldsh(ii,jj,kk,6,iidblock)+hfieldsh(ii,jj,kk,7,iidblock))) &
	   + fx*p1dcssq
	  !+1  0  0
	  
	  !7 -1 -1  0
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f07(li-1,lj,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfieldsh(ii,jj,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !9  -1 +1 0
      udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f09(li-1,lj,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfieldsh(ii,jj,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !15  -1  0 -1
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f15(li-1,lj,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfieldsh(ii,jj,kk,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
	  
	  !18   -1   0  +1
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f18(li-1,lj,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfieldsh(ii,jj,kk,9,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !2 +1  0  0
	  udotc=hfieldsh(ii,jj,kk,2,iidblock)*onecssq
	  f02(li+1,lj,lk)=p1*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qxx*hfieldsh(ii,jj,kk,5,iidblock) &
	   -cssq*(hfieldsh(ii,jj,kk,6,iidblock)+hfieldsh(ii,jj,kk,7,iidblock))) &
	   - fx*p1dcssq
	  !-1  0  0
	  
	  !8 +1 +1  0
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f08(li+1,lj,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfieldsh(ii,jj,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
	  
	  !10   +1 -1  0
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f10(li+1,lj,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfieldsh(ii,jj,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !16  +1  0 +1
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f16(li+1,lj,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfieldsh(ii,jj,kk,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
	  
	  !17  +1  0 -1
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f17(li+1,lj,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfieldsh(ii,jj,kk,9,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
    endif

    if(lj==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !3 0 -1  0
	  udotc=hfieldsh(ii,jj,kk,3,iidblock)*onecssq
	  f03(li,lj-1,lk)=p1*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*(hfieldsh(ii,jj,kk,5,iidblock)+hfieldsh(ii,jj,kk,7,iidblock))) &
	   + fy*p1dcssq
	  ! 0 +1  0
	  
	  !7 -1 -1  0
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f07(li,lj-1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfieldsh(ii,jj,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
	  
	  !10   +1 -1  0
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f10(li,lj-1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfieldsh(ii,jj,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
	  
	  !11  0  -1  -1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f11(li,lj-1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !13  0  -1   +1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)-hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f13(li,lj-1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !4  0 +1  0
	  udotc=hfieldsh(ii,jj,kk,3,iidblock)*onecssq
	  f04(li,lj+1,lk)=p1*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*(hfieldsh(ii,jj,kk,5,iidblock)+hfieldsh(ii,jj,kk,7,iidblock))) &
	   - fy*p1dcssq
	  ! 0 -1  0
      
      !8 +1 +1  0
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f08(li,lj+1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfieldsh(ii,jj,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
      !9  -1 +1 0
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f09(li,lj+1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfieldsh(ii,jj,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
	  
	  !12  0  +1  +1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f12(li,lj+1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfieldsh(ii,jj,kk,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
	  
	  !14  0  +1  -1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)-hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f14(li,lj+1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif

    if(lk==1)then
      
      gii=gi
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !5  0  0 -1
	  udotc=hfieldsh(ii,jj,kk,4,iidblock)*onecssq
	  f05(li,lj,lk-1)=p1*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*(hfieldsh(ii,jj,kk,5,iidblock)+hfieldsh(ii,jj,kk,6,iidblock))) &
	   + fz*p1dcssq
	  ! 0  0 +1
      
      !15  -1  0 -1
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f15(li,lj,lk-1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfieldsh(ii,jj,kk,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
      !17  +1  0 -1
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f17(li,lj,lk-1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfieldsh(ii,jj,kk,9,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
	  
	  !11  0  -1  -1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f11(li,lj,lk-1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
	  
	  !14  0  +1  -1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)-hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f14(li,lj,lk-1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
	  
      
    endif
    if(lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !6  0  0  +1
	  udotc=hfieldsh(ii,jj,kk,4,iidblock)*onecssq
	  f06(li,lj,lk+1)=p1*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq1*(qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*(hfieldsh(ii,jj,kk,5,iidblock)+hfieldsh(ii,jj,kk,6,iidblock))) &
	   - fz*p1dcssq
	  ! 0  0 -1
	  
	  !16  +1  0 +1
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f16(li,lj,lk+1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfieldsh(ii,jj,kk,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
      !18   -1   0  +1
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f18(li,lj,lk+1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfieldsh(ii,jj,kk,9,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
	  
	  !12  0  +1  +1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f12(li,lj,lk+1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfieldsh(ii,jj,kk,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1


	  !13  0  -1   +1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)-hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f13(li,lj,lk+1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      
      gii=gi-1
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !7 -1 -1  0
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f07(li-1,lj-1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfieldsh(ii,jj,kk,8,iidblock)) &
	   + (fx+fy)*p2dcssq 
	  !+1 +1  0
      
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      
      gii=gi-1
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !9  -1 +1 0
      udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f09(li-1,lj+1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfieldsh(ii,jj,kk,8,iidblock)) &
	   + (fx-fy)*p2dcssq
	  !+1 -1  0
      
    endif
    if(li==1 .and. lk==1)then
      
      gii=gi-1
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !15  -1  0 -1
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f15(li-1,lj,lk-1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfieldsh(ii,jj,kk,9,iidblock)) &
	   + (fx+fz)*p2dcssq 
	  !+1  0  +1
      
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi-1
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !18   -1   0  +1
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f18(li-1,lj,lk+1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfieldsh(ii,jj,kk,9,iidblock)) &
	   + (fx-fz)*p2dcssq
	  !+1  0  -1
      
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      
      gii=gi+1
      gjj=gj-1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !10   +1 -1  0
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f10(li+1,lj-1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_9_10*hfieldsh(ii,jj,kk,8,iidblock)) &
	   +(fy-fx)*p2dcssq
	  !-1 +1  0
      
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      
      gii=gi+1
      gjj=gj+1
      gkk=gk
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !8 +1 +1  0
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,3,iidblock))*onecssq
	  f08(li+1,lj+1,lk)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qyy*hfieldsh(ii,jj,kk,6,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,7,iidblock)+two*qxy_7_8*hfieldsh(ii,jj,kk,8,iidblock)) &
	   - (fx+fy)*p2dcssq
	  !-1 -1  0
      
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      
      gii=gi+1
      gjj=gj
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !17  +1  0 -1
	  udotc=(-hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f17(li+1,lj,lk-1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   + oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_17_18*hfieldsh(ii,jj,kk,9,iidblock)) &
	   +(fz-fx)*p2dcssq
	  !-1  0  +1
      
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      
      gii=gi+1
      gjj=gj
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !16  +1  0 +1
	  udotc=(hfieldsh(ii,jj,kk,2,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f16(li+1,lj,lk+1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qxx*hfieldsh(ii,jj,kk,5,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,6,iidblock)+two*qxz_15_16*hfieldsh(ii,jj,kk,9,iidblock)) &
	   - (fx+fz)*p2dcssq
	  !-1  0  -1
      
    endif
    if(lj==1 .and. lk==1)then
      
      gii=gi
      gjj=gj-1
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !11  0  -1  -1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f11(li,lj-1,lk-1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fy+fz)*p2dcssq
	  ! 0 +1 +1
      
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj-1
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !13  0  -1   +1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)-hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f13(li,lj-1,lk+1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc + udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fy-fz)*p2dcssq
	  ! 0 +1 -1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      
      gii=gi
      gjj=gj+1
      gkk=gk-1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !14  0  +1  -1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)-hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f14(li,lj+1,lk-1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_13_14*hfieldsh(ii,jj,kk,10,iidblock)) &
	   + (fz-fy)*p2dcssq
	  ! 0 -1 +1
      
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      
      gii=gi
      gjj=gj+1
      gkk=gk+1
      xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
      yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
      zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
      iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
      
      ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
      jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
      kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
      
      uu=halfonecssq*(hfieldsh(ii,jj,kk,2,iidblock)*hfieldsh(ii,jj,kk,2,iidblock) &
      + hfieldsh(ii,jj,kk,3,iidblock)*hfieldsh(ii,jj,kk,3,iidblock) + hfieldsh(ii,jj,kk,4,iidblock)*hfieldsh(ii,jj,kk,4,iidblock))
      
      !12  0  +1  +1
	  udotc=(hfieldsh(ii,jj,kk,3,iidblock)+hfieldsh(ii,jj,kk,4,iidblock))*onecssq
	  f12(li,lj+1,lk+1)=p2*(hfieldsh(ii,jj,kk,1,iidblock)+(-uu + half*udotc*udotc - udotc)) &
	   +oneminusomega*pi2cssq2*(qyy*hfieldsh(ii,jj,kk,6,iidblock)+qzz*hfieldsh(ii,jj,kk,7,iidblock) &
	   -cssq*hfieldsh(ii,jj,kk,5,iidblock)+two*qyz_11_12*hfieldsh(ii,jj,kk,10,iidblock)) &
	   - (fy+fz)*p2dcssq
	  ! 0 -1 -1
      
    endif      
    
    call syncthreads
    
    if(isfluid(gi,gj,gk).eq.0)then
#ifdef IFBC
          
          if(isfluid(gi+1,gj,gk-1).eq.-1)f18(li,lj,lk)=f17(li+1,lj,lk-1) !gpc 
          if(isfluid(gi-1,gj,gk+1).eq.-1)f17(li,lj,lk)=f18(li-1,lj,lk+1) !hpc

          if(isfluid(gi-1,gj,gk-1).eq.-1)f16(li,lj,lk)=f15(li-1,lj,lk-1) !gpc 
          if(isfluid(gi+1,gj,gk+1).eq.-1)f15(li,lj,lk)=f16(li+1,lj,lk+1) !hpc

          if(isfluid(gi,gj-1,gk+1).eq.-1)f14(li,lj,lk)=f13(li,lj-1,lk+1)!gpc 
          if(isfluid(gi,gj+1,gk-1).eq.-1)f13(li,lj,lk)=f14(li,lj+1,lk-1)!hpc
        
          if(isfluid(gi,gj-1,gk-1).eq.-1)f12(li,lj,lk)=f11(li,lj-1,lk-1)!gpc
          if(isfluid(gi,gj+1,gk+1).eq.-1)f11(li,lj,lk)=f12(li,lj+1,lk+1)!hpc

          if(isfluid(gi-1,gj+1,gk).eq.-1)f10(li,lj,lk)=f09(li-1,lj+1,lk)!gpc  
          if(isfluid(gi+1,gj-1,gk).eq.-1)f09(li,lj,lk)=f10(li+1,lj-1,lk)!hpc

          if(isfluid(gi-1,gj-1,gk).eq.-1)f08(li,lj,lk)=f07(li-1,lj-1,lk)!gpc 
          if(isfluid(gi+1,gj+1,gk).eq.-1)f07(li,lj,lk)=f08(li+1,lj+1,lk)!hpc

          if(isfluid(gi,gj,gk-1).eq.-1)f06(li,lj,lk)=f05(li,lj,lk-1)!gpc
          if(isfluid(gi,gj,gk+1).eq.-1)f05(li,lj,lk)=f06(li,lj,lk+1)!hpc 

          if(isfluid(gi,gj-1,gk).eq.-1)f04(li,lj,lk)=f03(li,lj-1,lk)!gpc 
          if(isfluid(gi,gj+1,gk).eq.-1)f03(li,lj,lk)=f04(li,lj+1,lk)!hpc 

          if(isfluid(gi-1,gj,gk).eq.-1)f02(li,lj,lk)=f01(li-1,lj,lk)!gpc
          if(isfluid(gi+1,gj,gk).eq.-1)f01(li,lj,lk)=f02(li+1,lj,lk)!hpc 
          
#else
          
        f01(li,lj,lk)=f02(li+1,lj,lk)!hpc 
        f02(li,lj,lk)=f01(li-1,lj,lk)!gpc 
        
        f03(li,lj,lk)=f04(li,lj+1,lk)!hpc 
        f04(li,lj,lk)=f03(li,lj-1,lk)!gpc 
        
        f05(li,lj,lk)=f06(li,lj,lk+1)!hpc 
        f06(li,lj,lk)=f05(li,lj,lk-1)!gpc 
        
        f07(li,lj,lk)=f08(li+1,lj+1,lk)!hpc
        f08(li,lj,lk)=f07(li-1,lj-1,lk)!gpc 
        
        f09(li,lj,lk)=f10(li+1,lj-1,lk)!hpc
        f10(li,lj,lk)=f09(li-1,lj+1,lk)!gpc 
        
        f11(li,lj,lk)=f12(li,lj+1,lk+1)!hpc
        f12(li,lj,lk)=f11(li,lj-1,lk-1)!gpc 
        
        f13(li,lj,lk)=f14(li,lj+1,lk-1)!hpc
        f14(li,lj,lk)=f13(li,lj-1,lk+1)!gpc 
        
        f15(li,lj,lk)=f16(li+1,lj,lk+1) !hpc
        f16(li,lj,lk)=f15(li-1,lj,lk-1) !gpc 
        
        f17(li,lj,lk)=f18(li-1,lj,lk+1) !hpc
        f18(li,lj,lk)=f17(li+1,lj,lk-1) !gpc       

#endif
    endif
    
    ! Halo Faces
    if(li==1)then
      if(isfluid(gi-1,gj,gk).eq.0)then
        f01(li-1,lj,lk)=f02(li,lj,lk)!hpc 
        
        f07(li-1,lj,lk)=f08(li,lj+1,lk)!hpc
        
        f09(li-1,lj,lk)=f10(li,lj-1,lk)!hpc
        
        f15(li-1,lj,lk)=f16(li,lj,lk+1) !hpc

        f18(li-1,lj,lk)=f17(li,lj,lk-1) !gpc 
      endif
    endif
    if(li==TILE_DIMx_d)then
      if(isfluid(gi+1,gj,gk).eq.0)then
        f02(li+1,lj,lk)=f01(li,lj,lk)!gpc 
        
        f08(li+1,lj,lk)=f07(li,lj-1,lk)!gpc 
        
        f10(li+1,lj,lk)=f09(li,lj+1,lk)!gpc 
        
        f16(li+1,lj,lk)=f15(li,lj,lk-1) !gpc 
        
        f17(li+1,lj,lk)=f18(li,lj,lk+1) !hpc 
      endif
    endif

    if(lj==1)then
      if(isfluid(gi,gj-1,gk).eq.0)then
        f03(li,lj-1,lk)=f04(li,lj,lk)!hpc 
        
        f07(li,lj-1,lk)=f08(li+1,lj,lk)!hpc
        
        f10(li,lj-1,lk)=f09(li-1,lj,lk)!gpc 
        
        f11(li,lj-1,lk)=f12(li,lj,lk+1)!hpc
        
        f13(li,lj-1,lk)=f14(li,lj,lk-1)!hpc    
      endif
    endif
    if(lj==TILE_DIMy_d)then
      if(isfluid(gi,gj+1,gk).eq.0)then
        f04(li,lj+1,lk)=f03(li,lj,lk)!gpc 
        
        f08(li,lj+1,lk)=f07(li-1,lj,lk)!gpc 
        
        f09(li,lj+1,lk)=f10(li+1,lj,lk)!hpc
        
        f12(li,lj+1,lk)=f11(li,lj,lk-1)!gpc 
        
        f14(li,lj+1,lk)=f13(li,lj,lk+1)!gpc   
      endif
    endif

    if(lk==1)then
      if(isfluid(gi,gj,gk-1).eq.0)then
        f05(li,lj,lk-1)=f06(li,lj,lk)!hpc 
        
        f11(li,lj,lk-1)=f12(li,lj+1,lk)!hpc
        
        f14(li,lj,lk-1)=f13(li,lj-1,lk)!gpc 
        
        f15(li,lj,lk-1)=f16(li+1,lj,lk) !hpc
        
        f17(li,lj,lk-1)=f18(li-1,lj,lk) !hpc
      endif
    endif
    if(lk==TILE_DIMz_d)then
      if(isfluid(gi,gj,gk+1).eq.0)then
        f06(li,lj,lk+1)=f05(li,lj,lk)!gpc 
        
        f12(li,lj,lk+1)=f11(li,lj-1,lk)!gpc 
        
        f13(li,lj,lk+1)=f14(li,lj+1,lk)!hpc
        
        f16(li,lj,lk+1)=f15(li-1,lj,lk) !gpc 
        
        f18(li,lj,lk+1)=f17(li+1,lj,lk) !gpc 
      endif
    endif

    ! Halo edges
    if(li==1 .and. lj==1)then
      if(isfluid(gi-1,gj-1,gk).eq.0)then
        f07(li-1,lj-1,lk)=f08(li,lj,lk)!hpc
      endif
    endif
    if(li==1 .and. lj==TILE_DIMy_d)then
      if(isfluid(gi-1,gj+1,gk).eq.0)then
        f09(li-1,lj+1,lk)=f10(li,lj,lk)!hpc
      endif
    endif
    if(li==1 .and. lk==1)then
      if(isfluid(gi-1,gj,gk-1).eq.0)then
        f15(li-1,lj,lk-1)=f16(li,lj,lk) !hpc
      endif
    endif
    if(li==1 .and. lk==TILE_DIMz_d)then
      if(isfluid(gi-1,gj,gk+1).eq.0)then
        f18(li-1,lj,lk+1)=f17(li,lj,lk) !gpc    
      endif
    endif
    if(li==TILE_DIMx_d .and. lj==1)then
      if(isfluid(gi+1,gj-1,gk).eq.0)then
        f10(li+1,lj-1,lk)=f09(li,lj,lk)!gpc 
      endif
    endif
    if(li==TILE_DIMx_d .and. lj==TILE_DIMy_d)then
      if(isfluid(gi+1,gj+1,gk).eq.0)then
        f08(li+1,lj+1,lk)=f07(li,lj,lk)!gpc 
      endif
    endif
    if(li==TILE_DIMx_d .and. lk==1)then
      if(isfluid(gi+1,gj,gk-1).eq.0)then
        f17(li+1,lj,lk-1)=f18(li,lj,lk) !hpc
      endif
    endif
    if(li==TILE_DIMx_d .and. lk==TILE_DIMz_d)then
      if(isfluid(gi+1,gj,gk+1).eq.0)then
        f16(li+1,lj,lk+1)=f15(li,lj,lk) !gpc 
      endif
    endif
    if(lj==1 .and. lk==1)then
      if(isfluid(gi,gj-1,gk-1).eq.0)then
        f11(li,lj-1,lk-1)=f12(li,lj,lk)!hpc
      endif
    endif
    if(lj==1 .and. lk==TILE_DIMz_d)then
      if(isfluid(gi,gj-1,gk+1).eq.0)then
        f13(li,lj-1,lk+1)=f14(li,lj,lk)!hpc
      endif
    endif
    if(lj==TILE_DIMy_d .and. lk==1)then
      if(isfluid(gi,gj+1,gk-1).eq.0)then
        f14(li,lj+1,lk-1)=f13(li,lj,lk)!gpc 
      endif
    endif
    if(lj==TILE_DIMy_d .and. lk==TILE_DIMz_d)then
      if(isfluid(gi,gj+1,gk+1).eq.0)then
        f12(li,lj+1,lk+1)=f11(li,lj,lk)!gpc 
      endif
    endif      
    
    call syncthreads
    
    udotc=f00(li,lj,lk)+f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,1,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)-f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,2,myblock)=udotc
	
	
	udotc=f03(li,lj-1,lk)-f04(li,lj+1,lk)+ &
     f07(li-1,lj-1,lk)-f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,3,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)-f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)-f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)-f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,4,myblock)=udotc
	
	udotc=f01(li-1,lj,lk)+f02(li+1,lj,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,5,myblock)=udotc
	
	udotc=f03(li,lj-1,lk)+f04(li,lj+1,lk)+  &
     f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)+ &
     f09(li-1,lj+1,lk)+f10(li+1,lj-1,lk)+ &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)
	hfields(i,j,k,6,myblock)=udotc
	
	udotc=f05(li,lj,lk-1)+f06(li,lj,lk+1)+  &
     f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)+ &
     f13(li,lj-1,lk+1)+f14(li,lj+1,lk-1)+ &
     f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)+ &
     f17(li+1,lj,lk-1)+f18(li-1,lj,lk+1)
	hfields(i,j,k,7,myblock)=udotc
	
	udotc=f07(li-1,lj-1,lk)+f08(li+1,lj+1,lk)- &
     f09(li-1,lj+1,lk)-f10(li+1,lj-1,lk)
	hfields(i,j,k,8,myblock)=udotc
	
	udotc=f15(li-1,lj,lk-1)+f16(li+1,lj,lk+1)- &
     f17(li+1,lj,lk-1)-f18(li-1,lj,lk+1)
	hfields(i,j,k,9,myblock)=udotc
	
	udotc=f11(li,lj-1,lk-1)+f12(li,lj+1,lk+1)- &
     f13(li,lj-1,lk+1)-f14(li,lj+1,lk-1)
	hfields(i,j,k,10,myblock)=udotc
     
#ifdef PRESSCORR

	call syncthreads
	
	uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) &
	+ hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
    
	!1 -1  0  0
	udotc=hfields(i,j,k,2,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f01(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  0


	!2 +1  0  0
	f02(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  0

    		
	!3 0 -1  0
	udotc=hfields(i,j,k,3,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f03(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1  0

	
	!4  0 +1  0
	f04(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1  0

	
	!5  0  0 -1
	udotc=hfields(i,j,k,4,myblock)*onecssq
	temp = -uu + half*udotc*udotc
	f05(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0  0 +1


	!6  0  0  +1
	f06(li,lj,lk)=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0  0 -1

    	
	!7 -1 -1  0
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f07(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1 +1  0

	
	!8 +1 +1  0
	f08(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1 -1  0

	
	!10   +1 -1  0
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f10(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1 +1  0

	
	!9  -1 +1 0
	f09(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1 -1  0

		

	!15  -1  0 -1
	udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f15(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!+1  0  +1


	!16  +1  0 +1
	f16(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!-1  0  -1


	!17  +1  0 -1
	udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f17(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	!-1  0  +1


	!18   -1   0  +1
	f18(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	!+1  0  -1


	!11  0  -1  -1
	udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f11(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 +1

	
	!12  0  +1  +1
	f12(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 -1


	!13  0  -1   +1
	udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
	temp = -uu + half*udotc*udotc
	f13(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
	! 0 +1 -1

	
	!14  0  +1  -1
	f14(li,lj,lk)=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
	! 0 -1 +1
	
	udotc=f01(li,lj,lk)+f02(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,5,myblock)=hfields(i,j,k,5,myblock)-udotc
	
	udotc=f03(li,lj,lk)+f04(li,lj,lk)+  &
     f07(li,lj,lk)+f08(li,lj,lk)+ &
     f09(li,lj,lk)+f10(li,lj,lk)+ &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)
	hfields(i,j,k,6,myblock)=hfields(i,j,k,6,myblock)-udotc
	
	udotc=f05(li,lj,lk)+f06(li,lj,lk)+  &
     f11(li,lj,lk)+f12(li,lj,lk)+ &
     f13(li,lj,lk)+f14(li,lj,lk)+ &
     f15(li,lj,lk)+f16(li,lj,lk)+ &
     f17(li,lj,lk)+f18(li,lj,lk)
	hfields(i,j,k,7,myblock)=hfields(i,j,k,7,myblock)-udotc
	
	udotc=f07(li,lj,lk)+f08(li,lj,lk)- &
     f09(li,lj,lk)-f10(li,lj,lk)
	hfields(i,j,k,8,myblock)=hfields(i,j,k,8,myblock)-udotc
	
	udotc=f15(li,lj,lk)+f16(li,lj,lk)- &
     f17(li,lj,lk)-f18(li,lj,lk)
	hfields(i,j,k,9,myblock)=hfields(i,j,k,9,myblock)-udotc
	
	udotc=f11(li,lj,lk)+f12(li,lj,lk)- &
     f13(li,lj,lk)-f14(li,lj,lk)
	hfields(i,j,k,10,myblock)=hfields(i,j,k,10,myblock)-udotc
	
#endif
    
    return
 
  end subroutine streamcoll_bc_shared_flop
  
 end module streamcoll_bc_kernels
