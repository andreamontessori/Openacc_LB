#include "defines.h"

program lb_openacc
    
    use cudavars
    use pbc_kernels
    use streamcoll_bulk_kernels
    use streamcoll_bc_kernels
    use correct_press_kernels
    use prints
    
    implicit none
    
    
    integer :: i,j,k,ii,jj,kk,ll
    integer :: nx,ny,nz,step,istep,stamp,nlinks,nsteps,ngpus,devNum
    integer :: TILE_DIMx,TILE_DIMy,TILE_DIMz,TILE_DIM,istat,iframe
    
    logical :: lpbc
    
    
    
    real(kind=4)  :: ts1,ts2
    real(kind=db) :: visc_LB,tau,one_ov_nu,h_fx,h_fy,h_fz,h_omega,h_oneminusomega
    

    integer, parameter :: npops=18
                                                  
    integer, parameter, dimension(0:npops) :: &
     ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
         !0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18
    integer, parameter, dimension(0:npops) :: &
     ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
    integer, parameter, dimension(0:npops) :: &
     ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
    integer, parameter, dimension(0:npops) :: &
    opp=(/0, 2,  1, 4,  3,  6,  5,  8,   7, 10,   9, 12,  11, 14,  13, 16,  15,  18,  17/)
    integer(kind=4), allocatable, dimension(:,:,:) :: h_isfluid
    
    real(kind=db) :: mymemory,totmemory
    real(kind=db) :: time
    
    integer :: mshared
    
    
       
    nlinks=18 !pari!
    tau=real(1.5d0,kind=db)
    visc_LB=cssq*(tau-half)
    one_ov_nu=one/visc_LB


    !*******************************user parameters and allocations**************************m
        nx=256
        ny=256
        nz=256
        nsteps=10
        stamp=100
        h_fx=one*ten**(-real(5.d0,kind=db))
        h_fy=zero*ten**(-real(5.d0,kind=db))
        h_fz=zero*ten**(-real(5.d0,kind=db))
        lpbc=.true.
        lprint=.false.
        lvtk=.false.
        lasync=.false.
        
        TILE_DIMx=256
        TILE_DIMy=1
        TILE_DIMz=1
        TILE_DIM=16
        if (mod(nx, TILE_DIMx)/= 0) then
          write(*,*) 'nx must be a multiple of TILE_DIM'
          stop
        end if
        if (mod(ny, TILE_DIMy) /= 0) then
          write(*,*) 'ny must be a multiple of TILE_DIMy'
          stop
        end if
        if (mod(nz, TILE_DIMz) /= 0) then
          write(*,*) 'nz must be a multiple of TILE_DIMz'
          stop
        end if
        dimGrid  = dim3(nx/TILE_DIMx, ny/TILE_DIMy, nz/TILE_DIMz)
        dimBlock = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
        
        dimGridx  = dim3((ny+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
        dimGridy  = dim3((nx+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
        dimBlock2 = dim3(TILE_DIM, TILE_DIM, 1)
        
        allocate(rho(0:nx+1,0:ny+1,0:nz+1),u(0:nx+1,0:ny+1,0:nz+1),v(0:nx+1,0:ny+1,0:nz+1),w(0:nx+1,0:ny+1,0:nz+1))
        allocate(pxx(0:nx+1,0:ny+1,0:nz+1),pxy(0:nx+1,0:ny+1,0:nz+1),pxz(0:nx+1,0:ny+1,0:nz+1),pyy(0:nx+1,0:ny+1,0:nz+1))
        allocate(pyz(0:nx+1,0:ny+1,0:nz+1),pzz(0:nx+1,0:ny+1,0:nz+1))
        allocate(rhoh(0:nx+1,0:ny+1,0:nz+1),uh(0:nx+1,0:ny+1,0:nz+1),vh(0:nx+1,0:ny+1,0:nz+1),wh(0:nx+1,0:ny+1,0:nz+1))
        allocate(pxxh(0:nx+1,0:ny+1,0:nz+1),pxyh(0:nx+1,0:ny+1,0:nz+1),pxzh(0:nx+1,0:ny+1,0:nz+1),pyyh(0:nx+1,0:ny+1,0:nz+1))
        allocate(pyzh(0:nx+1,0:ny+1,0:nz+1),pzzh(0:nx+1,0:ny+1,0:nz+1))
        allocate(h_isfluid(0:nx+1,0:ny+1,0:nz+1))
        allocate(isfluid(0:nx+1,0:ny+1,0:nz+1)) 
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz))
          allocate(velprint(1:3,1:nx,1:ny,1:nz))
          rhoprint(1:nx,1:ny,1:nz)=0.0
          velprint(1:3,1:nx,1:ny,1:nz)=0.0
          
          allocate(rhoprint_d(1:nx,1:ny,1:nz))
          allocate(velprint_d(1:3,1:nx,1:ny,1:nz))
          rhoprint_d(1:nx,1:ny,1:nz)=0.0
          velprint_d(1:3,1:nx,1:ny,1:nz)=0.0
        endif
        
       
        h_omega=one/tau
        h_oneminusomega=one-h_omega
        fx=h_fx
        fy=h_fy
        fz=h_fz
        omega=h_omega
        oneminusomega=h_oneminusomega
        nx_d=nx
        ny_d=ny
        nz_d=nz
        TILE_DIMx_d=TILE_DIMx
        TILE_DIMy_d=TILE_DIMy
        TILE_DIMz_d=TILE_DIMz
        TILE_DIM_d=TILE_DIM
    !*****************************************geometry************************
        h_isfluid=0
        h_isfluid(1:nx,1:ny,1:nz)=1
!        h_isfluid=1
!        h_isfluid(1,:,:)=0 !left
!        h_isfluid(nx,:,:)=0 !right
!        h_isfluid(:,1,:)=0 !front 
!        h_isfluid(:,ny,:)=0 !rear
!        h_isfluid(:,:,1)=0 !bottom
!        h_isfluid(:,:,nz)=0 !top
        if(lpbc)then
          h_isfluid=1
          h_isfluid(:,:,0)=0 !bottom
          h_isfluid(:,:,nz+1)=0 !top
        endif
        do k=1,nz
	      do j=1,ny
		    do i=1,nx
			  if(h_isfluid(i,j,k).eq.1)then
			    do ll=1,npops
				  ii=i+ex(ll)
			      jj=j+ey(ll)
				  kk=k+ez(ll) 
				  if(ii.gt.0 .and. ii.lt.nx+1 .and. jj.gt.0 .and. jj.lt.ny+1 .and. kk.gt.0 .and. kk.lt.nz+1)then
				    if(h_isfluid(ii,jj,kk).eq.0)then
					  h_isfluid(i,j,k)=-1
					endif
				  endif
				enddo
			  endif
		    enddo
		  enddo
	    enddo
	    
	    istat = cudaDeviceSynchronize
        istat = cudaMemcpy(isfluid,h_isfluid,(nx+2)*(ny+2)*(nz+2) )
        istat = cudaDeviceSynchronize
    !****************************************hermite projection vars**********
        

    !*************************************initial conditions ************************    
    
    call setup_system<<<dimGrid,dimBlock>>>(one,zero,zero,zero)
        
        
    !*************************************check data ************************ 
	istat = cudaGetDeviceCount(ngpus)
	istat = cudaGetDevice(devNum)
	
	write(6,*) '*******************LB data*****************'
	write(6,*) 'tau',tau
	write(6,*) 'omega',h_omega
	write(6,*) 'visc',visc_LB
	write(6,*) 'fx',h_fx
	write(6,*) 'fy',h_fy
	write(6,*) 'fz',h_fz
	write(6,*) 'cssq',cssq
	write(6,*) '*******************INPUT data*****************'
	write(6,*) 'nx',nx
	write(6,*) 'ny',ny
	write(6,*) 'ny',nz
	write(6,*) 'lpbc',lpbc
	write(6,*) 'lprint',lprint
	write(6,*) 'lvtk',lvtk
	write(6,*) 'lasync',lasync
	write(6,*) 'nsteps',nsteps
	write(6,*) 'stamp',stamp
	write(6,*) 'max fx',huge(fx)
	write(6,*) 'max fx',huge(fy)
	write(6,*) 'max fx',huge(fz)
	write(6,*) 'available gpus',ngpus
	write(6,*) 'my ID gpu',devNum
	write(6,*) '*******************************************'
	write(6,*) 'TILE_DIMx',TILE_DIMx
    write(6,*) 'TILE_DIMy',TILE_DIMy
    write(6,*) 'TILE_DIMz',TILE_DIMz
    write(6,*) 'TILE_DIM ',TILE_DIM
	write(6,*) '*******************************************'
	
	
	istat = cudaGetDeviceProperties(prop, devNum)
#ifdef USESHARE 
    mshared = prop%sharedMemPerBlock
#else
    mshared = 0
#endif
	call printDeviceProperties(prop,6, devNum)
        
	! create events and streams
	istat = cudaStreamCreate(stream1)
	istat = cudaStreamCreate(stream2)
	istat = cudaforSetDefaultstream(stream1)
	istat = cudaDeviceSynchronize


	istat = cudaEventCreate(startEvent)
	istat = cudaEventCreate(stopEvent)  
	istat = cudaEventCreate(dummyEvent)
	istat = cudaEventCreate(dummyEvent1)
	istat = cudaEventCreate(dummyEvent2)
        
    step = 0
    
    
    if(lprint)then  
      call init_output(nx,ny,nz,1,lvtk)
      call string_char(head1,nheadervtk(1),headervtk(1))
      call string_char(head2,nheadervtk(2),headervtk(2))
    endif
    
    iframe=0
    write(6,'(a,i8,a,i8,3f16.4)')'start step : ',0,' frame ',iframe
    
    istat = cudaDeviceSynchronize
    
    if(lprint)then
      call store_print<<<dimGrid,dimBlock,0,stream1>>>()
      istat = cudaEventRecord(dummyEvent1, stream1)
      istat = cudaEventSynchronize(dummyEvent1)
      if(lasync)then
        istat = cudaMemcpyAsync(rhoprint,rhoprint_d,nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
        istat = cudaMemcpyAsync(velprint,velprint_d,3*nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
      else
        istat = cudaMemcpy(rhoprint,rhoprint_d,nx*ny*nz )
        istat = cudaMemcpy(velprint,velprint_d,3*nx*ny*nz )
        istat = cudaEventRecord(dummyEvent, 0)
        istat = cudaEventSynchronize(dummyEvent)
        if(lvtk)then
          call print_vtk_sync(iframe)
        else
          call print_raw_sync(iframe)
        endif
      endif
    endif
      
    !*************************************time loop************************  
    call cpu_time(ts1)
    istat = cudaEventRecord(startEvent,0)
    do step=1,nsteps,2 
        !flip
        istep=step
        !******************************************call other bcs************************
        if(lpbc)then      
            !periodic along x 
            call pbc_edge_x<<<dimGridx,dimBlock2,0,stream1>>>()
            !periodic along x 
            call pbc_edge_y<<<dimGridy,dimBlock2,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
	    endif
        
        !***********************************PRINT************************
        if(mod(istep,stamp).eq.0)write(6,'(a,i8)')'step : ',istep
        if(lprint)then
          if(mod(istep,stamp).eq.0)then
            iframe=iframe+1
            call store_print<<<dimGrid,dimBlock,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
            if(lasync)then
              call close_print_async
              istat = cudaMemcpyAsync(rhoprint,rhoprint_d,nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
              istat = cudaMemcpyAsync(velprint,velprint_d,3*nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
            else
              istat = cudaMemcpy(rhoprint,rhoprint_d,nx*ny*nz )
              istat = cudaMemcpy(velprint,velprint_d,3*nx*ny*nz )
              istat = cudaEventRecord(dummyEvent, 0)
              istat = cudaEventSynchronize(dummyEvent)
              if(lvtk)then
                call print_vtk_sync(iframe)
              else
                call print_raw_sync(iframe)
              endif
            endif
          endif
          if(mod(istep-stamp/4,stamp).eq.0 .and. lasync)then
            !write(6,*)'ciao 2',step,iframe
            istat = cudaEventRecord(dummyEvent2, stream2)
            istat = cudaEventSynchronize(dummyEvent2)
            if(lvtk)then
              call print_vtk_async(iframe)
            else
              call print_raw_async(iframe)
            endif
          endif
        endif
        
        !***********************************collision + no slip + forcing: fused implementation*********
        call streamcoll_bulk<<<dimGrid,dimBlock,mshared,stream1>>>()
        !call streamcoll_bulk01<<<dimGrid,dimBlock,mshared,stream1>>>()
        !call streamcoll_bulk012xy<<<dimGrid,dimBlock,mshared,stream1>>>()
        !call streamcoll_bulk012xz<<<dimGrid,dimBlock,mshared,stream1>>>()
        !call streamcoll_bulk012yz<<<dimGrid,dimBlock,mshared,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after streamcoll_bulk', istep)

        !********************************close to boundary conditions no slip everywhere********************************!
        
        call streamcoll_bc<<<dimGrid,dimBlock,0,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after streamcoll_bc', istep)
        
        !***********************************correct pressor*********
        call correct_pressure<<<dimGrid,dimBlock,0,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after correct_pressure', istep)
        
        !flop
        istep=step+1
        !******************************************call other bcs************************
        if(lpbc)then      
            !periodic along x 
            call pbc_edge_x_flop<<<dimGridx,dimBlock2,0,stream1>>>()
            !periodic along x 
            call pbc_edge_y_flop<<<dimGridy,dimBlock2,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
	    endif
        
        !***********************************PRINT************************
        if(mod(istep,stamp).eq.0)write(6,'(a,i8)')'step : ',istep
        if(lprint)then
          if(mod(istep,stamp).eq.0)then
            iframe=iframe+1
            call store_print_flop<<<dimGrid,dimBlock,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
            if(lasync)then
              call close_print_async
              istat = cudaMemcpyAsync(rhoprint,rhoprint_d,nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
              istat = cudaMemcpyAsync(velprint,velprint_d,3*nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
            else
              istat = cudaMemcpy(rhoprint,rhoprint_d,nx*ny*nz )
              istat = cudaMemcpy(velprint,velprint_d,3*nx*ny*nz )
              istat = cudaEventRecord(dummyEvent, 0)
              istat = cudaEventSynchronize(dummyEvent)
              if(lvtk)then
                call print_vtk_sync(iframe)
              else
                call print_raw_sync(iframe)
              endif
            endif
          endif
          if(mod(istep-stamp/4,stamp).eq.0 .and. lasync)then
            !write(6,*)'ciao 2',step,iframe
            istat = cudaEventRecord(dummyEvent2, stream2)
            istat = cudaEventSynchronize(dummyEvent2)
            if(lvtk)then
              call print_vtk_async(iframe)
            else
              call print_raw_async(iframe)
            endif
          endif
        endif
        
        
        !***********************************collision + no slip + forcing: fused implementation*********
        
        call streamcoll_bulk_flop<<<dimGrid,dimBlock,mshared,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after streamcoll_bulk_flop', istep)
        
        !********************************close to boundary conditions no slip everywhere********************************!
        
        call streamcoll_bc_flop<<<dimGrid,dimBlock,0,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after streamcoll_bc_flop', istep)
        
        !***********************************correct pressor*********
        call correct_pressure_flop<<<dimGrid,dimBlock,0,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after correct_pressure_flop', istep)
        
    enddo 
    
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time, startEvent, stopEvent)
    if(lasync)then
      istat = cudaEventRecord(dummyEvent2, stream2)
      istat = cudaEventSynchronize(dummyEvent2)
      if(lvtk)then
        call print_vtk_sync(iframe)
      else
        call print_raw_sync(iframe)
      endif
    endif
    call cpu_time(ts2)
    
    

    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time'
    write(6,*) 'cuda time elapsed: ', time, ' s of your life time'
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)*real(1.d-9,kind=db)/(ts2-ts1)
    write(6,*) 'glups cuda: ',  real(nx)*real(ny)*real(nz)*real(nsteps)*real(1.d-9,kind=db)*1000/time
    
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)


  

  
    
end program
