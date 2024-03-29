#include "defines.h"

program lb_openacc
    
    use cudavars
    use pbc_kernels
    use streamcoll_kernels
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
    integer(kind=1), allocatable, dimension(:,:,:) :: h_isfluid
    
    real(kind=db) :: mymemory,totmemory
    real(kind=db) :: time
    
    integer :: mshared
    integer :: nxblock,nyblock,nzblock,nblocks,idblock,xblock,yblock,zblock,nxyblock
    integer :: ciao(3)
    
    logical :: xperiodic,yperiodic,zperiodic
    
    nlinks=18 !pari!
    tau=real(1.5d0,kind=db)
    visc_LB=cssq*(tau-half)
    one_ov_nu=one/visc_LB


    !*******************************user parameters and allocations**************************m
        nx=32
        ny=32
        nz=32
        nsteps=5
        stamp=1
        h_fx=one*ten**(-real(4.d0,kind=db))
        h_fy=zero*ten**(-real(5.d0,kind=db))
        h_fz=zero*ten**(-real(5.d0,kind=db))
        lpbc=.true.
        lprint=.true.
        lvtk=.true.
        lasync=.false.
        
        TILE_DIMx=4
        TILE_DIMy=4
        TILE_DIMz=4
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
        
        
        dimGridhalo  = dim3((nx+TILE_DIMx-1)/TILE_DIMx+2,(ny+TILE_DIMy-1)/TILE_DIMy+2,(nz+TILE_DIMz-1)/TILE_DIMz+2)
        dimBlockhalo = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
        
        dimBlockshared = dim3(TILE_DIMx+2, TILE_DIMy+2, TILE_DIMz+2)
        
        
        dimGridx  = dim3((ny+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
        dimGridy  = dim3((nx+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
        dimGridz  = dim3((nx+TILE_DIM-1)/TILE_DIM, (ny+TILE_DIM-1)/TILE_DIM, 1)
        
        dimBlock2 = dim3(TILE_DIM, TILE_DIM, 1)
        !plus 2 for the halo forward and backward
        nxblock=nx/TILE_DIMx +2
        nyblock=ny/TILE_DIMy +2
        nzblock=nz/TILE_DIMz +2
        
        nxyblock=nxblock*nyblock
        
        nblocks=nxblock*nyblock*nzblock
        
        write(6,*)'nx,ny,nz',nx,ny,nz
        write(6,*)'TILE_DIMx,TILE_DIMy,TILE_DIMz',TILE_DIMx,TILE_DIMy,TILE_DIMz
        write(6,*)'nxblock,nyblock,nzblock',nxblock,nyblock,nzblock
        write(6,*)'nblocks',nblocks
        
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
        nxblock_d=nxblock
        nxyblock_d=nxyblock
        nblocks_d=nblocks
        
#if 0   
       k=1
       j=1   
        do k=1-TILE_DIMz,nz+TILE_DIMz
        do j=1-TILE_DIMy,ny+TILE_DIMy
        do i=1-TILE_DIMx,nx+TILE_DIMx
          xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
          yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
          zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
          !idblock start from 1
          idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*(nxblock*nyblock)+1
          !return coord of block from its id
          ciao=coordblock(idblock)
          ii=i-xblock*TILE_DIMx+2*TILE_DIMx
          jj=j-yblock*TILE_DIMy+2*TILE_DIMy
          kk=k-zblock*TILE_DIMz+2*TILE_DIMz
          write(6,'(10i4)')i,j,k,xblock,yblock,zblock,i-xblock*TILE_DIMx+TILE_DIMx, &
           j-yblock*TILE_DIMy+TILE_DIMy,k-zblock*TILE_DIMz+TILE_DIMz
          if(xblock/=ciao(1) .or. yblock/=ciao(2) .or. zblock/=ciao(3) .or. idblock>nblocks)then
            write(6,*)'cazzo',ciao(1),xblock
            stop
          endif
        enddo
        enddo
        enddo
        
#endif    

#if 0  
        do k=1-TILE_DIMz,nz+TILE_DIMz
        do j=1-TILE_DIMy,ny+TILE_DIMy
        do i=1-TILE_DIMx,nx+TILE_DIMx
          xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
          yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
          zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
          !idblock start from 1
          idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*(nxblock*nyblock)+1
          !return coord of block from its id
          ciao=coordblock(idblock)
          ii=i-xblock*TILE_DIMx+2*TILE_DIMx
          jj=j-yblock*TILE_DIMy+2*TILE_DIMy
          kk=k-zblock*TILE_DIMz+2*TILE_DIMz
          if(abs(xblock-2)<=1 .and. abs(yblock-2)<=1 .and. abs(zblock-2)<=1)then
            if(xblock==2 .and. yblock==2 .and. zblock==2)then
              write(6,'(a,7i5)')'azzo dentro',i,j,k,idblock!,xblock,yblock,zblock

            else
              if(i>=0 .and. i<=5 .and. j>=0 .and. j<=3 .and. k>=0 .and. k<=3)then
              write(6,'(a,7i5)')'azzo intorno',i,j,k,idblock!,xblock,yblock,zblock
              endif
            endif
          endif
          
          !if(idblock==30)write(6,'(10i4)')i,j,k,ii,jj,kk,idblock,xblock,yblock,zblock
          if(xblock/=ciao(1) .or. yblock/=ciao(2) .or. zblock/=ciao(3) .or. idblock>nblocks)then
            write(6,*)'cazzo',ciao(1),xblock
            stop
          endif
           call setup_system_halo2<<<dimGridhalo,dimBlockhalo>>>(1.0,0.0,0.0,0.0,idblock,i,j,k,ii,jj,kk)
           istat = cudaDeviceSynchronize
           
           if(i<0)cycle
           if(i>nx)cycle
           if(j<0)cycle
           if(j>ny)cycle
           if(k<0)cycle
           if(k>nz)cycle
           call setup_system_bulk2<<<dimGrid,dimBlock>>>(1.0,0.0,0.0,0.0,idblock,i,j,k,ii,jj,kk)
           istat = cudaDeviceSynchronize
        enddo
        enddo
        enddo
        
#endif   
        
#if 0
 
        do k=1,nz
        do j=1,ny
        do i=1,nx
          xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
          yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
          zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
          !idblock start from 1
          idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*(nxblock*nyblock)+1
          !return coord of block from its id
          ciao=coordblock(idblock)
          ii=i-xblock*TILE_DIMx+2*TILE_DIMx
          jj=j-yblock*TILE_DIMy+2*TILE_DIMy
          kk=k-zblock*TILE_DIMz+2*TILE_DIMz
          !write(6,'(10i4)')i,j,k,xblock,yblock,zblock,idblock
!          if(xblock/=ciao(1) .or. yblock/=ciao(2) .or. zblock/=ciao(3) .or. idblock>nblocks)then
!            write(6,*)'cazzo',ciao(1),xblock
!            stop
!          endif
!          if(xblock==2 .and. yblock==2 .and. zblock==2)then
!            write(6,'(a,7i5)')'azzo',i,j,k,idblock,xblock,yblock,zblock

!          endif

          !if(abs(xblock-2)<=1 .and. abs(yblock-2)<=1 .and. abs(zblock-2)<=1)then
            if(xblock==2 .and. yblock==2 .and. zblock==2)then
            write(6,'(a,7i5)')'azzo dentro2',i,j,k,idblock,xblock,yblock,zblock

          !  else
          !    write(6,'(a,7i5)')'azzo intorno',i,j,k,idblock,xblock,yblock,zblock
            endif
          !endif
        enddo
        enddo
        enddo
        call flush(6)
#endif 
#if 0         
         call setup_system_bulk3<<<dimGrid,dimBlockshared>>>(1.0,0.0,0.0,0.0,idblock,i,j,k,ii,jj,kk)
           istat = cudaDeviceSynchronize
           stop
#endif        
            
    
  
    
        
!        allocate(rho(0:nx+1,0:ny+1,0:nz+1),u(0:nx+1,0:ny+1,0:nz+1),v(0:nx+1,0:ny+1,0:nz+1),w(0:nx+1,0:ny+1,0:nz+1))
!        allocate(pxx(0:nx+1,0:ny+1,0:nz+1),pxy(0:nx+1,0:ny+1,0:nz+1),pxz(0:nx+1,0:ny+1,0:nz+1),pyy(0:nx+1,0:ny+1,0:nz+1))
!        allocate(pyz(0:nx+1,0:ny+1,0:nz+1),pzz(0:nx+1,0:ny+1,0:nz+1))
!        allocate(rhoh(0:nx+1,0:ny+1,0:nz+1),uh(0:nx+1,0:ny+1,0:nz+1),vh(0:nx+1,0:ny+1,0:nz+1),wh(0:nx+1,0:ny+1,0:nz+1))
!        allocate(pxxh(0:nx+1,0:ny+1,0:nz+1),pxyh(0:nx+1,0:ny+1,0:nz+1),pxzh(0:nx+1,0:ny+1,0:nz+1),pyyh(0:nx+1,0:ny+1,0:nz+1))
!        allocate(pyzh(0:nx+1,0:ny+1,0:nz+1),pzzh(0:nx+1,0:ny+1,0:nz+1))
        
        
!        allocate(rho(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),u(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks), &
!         v(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),w(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(pxx(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pxy(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),&
!         pxz(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pyy(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(pyz(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pzz(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(rhoh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),uh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks), &
!         vh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),wh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(pxxh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pxyh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks), &
!         pxzh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pyyh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(pyzh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pzzh(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
        
        !allocate(hfields(TILE_DIMx,TILE_DIMy,TILE_DIMz,10,nblocks),hfieldsh(TILE_DIMx,TILE_DIMy,TILE_DIMz,10,nblocks))
        allocate(hfields(TILE_DIMx,TILE_DIMy,TILE_DIMz,10,nblocks))
        
        allocate(xshell0(TILE_DIMy,TILE_DIMz,10,nblocks),xshell1(TILE_DIMy,TILE_DIMz,10,nblocks))
        allocate(xshell0h(TILE_DIMy,TILE_DIMz,10,nblocks),xshell1h(TILE_DIMy,TILE_DIMz,10,nblocks))
        
        allocate(yshell0(TILE_DIMx,TILE_DIMz,10,nblocks),yshell1(TILE_DIMx,TILE_DIMz,10,nblocks))
        allocate(yshell0h(TILE_DIMx,TILE_DIMz,10,nblocks),yshell1h(TILE_DIMx,TILE_DIMz,10,nblocks))
        
        allocate(zshell0(TILE_DIMx,TILE_DIMy,10,nblocks),zshell1(TILE_DIMx,TILE_DIMy,10,nblocks))
        allocate(zshell0h(TILE_DIMx,TILE_DIMy,10,nblocks),zshell1h(TILE_DIMx,TILE_DIMy,10,nblocks))
        
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
        
    !*****************************************geometry************************
        h_isfluid=0
        h_isfluid(1:nx,1:ny,1:nz)=1
!        h_isfluid=1
        h_isfluid(1,:,:)=0 !left
        h_isfluid(nx,:,:)=0 !right
        h_isfluid(:,1,:)=0 !front 
        h_isfluid(:,ny,:)=0 !rear
        h_isfluid(:,:,1)=0 !bottom
        h_isfluid(:,:,nz)=0 !top
        if(lpbc)then
          h_isfluid=1
          h_isfluid(:,:,0)=3 !bottom
          h_isfluid(:,:,nz+1)=3 !top
        endif
        do k=0,nz+1
	      do j=0,ny+1
		    do i=0,nx+1
			  if(h_isfluid(i,j,k).eq.1)then
			    do ll=1,npops
				  ii=i+ex(ll)
			      jj=j+ey(ll)
				  kk=k+ez(ll) 
				  if(ii.ge.0 .and. ii.le.nx+1 .and. jj.ge.0 .and. jj.le.ny+1 .and. kk.ge.0 .and. kk.le.nz+1)then
				    if(h_isfluid(ii,jj,kk).eq.3)then
					  h_isfluid(i,j,k)=0
					endif
				  endif
				enddo
			  endif
			  
		    enddo
		  enddo
	    enddo
	    
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
     xperiodic=.false.
    !*************************************initial conditions ************************    

    call setup_system_halo<<<dimGridhalo,dimBlockhalo>>>(one,zero,zero,zero)

    istat = cudaDeviceSynchronize
    call abortOnLastErrorAndSync('after setup_system', 0)
    
        
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
    write(6,*)'nxblock,nyblock,nzblock',nxblock,nyblock,nzblock
    write(6,*)'nblocks',nblocks
	write(6,*) '*******************************************'
	
	
	istat = cudaGetDeviceProperties(prop, devNum)
	
#ifdef USESHARE 
    mshared = prop%sharedMemPerBlock
#else
    mshared = 0
#endif
    mshared = prop%sharedMemPerBlock
    
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
            call pbc_side_x<<<dimGridx,dimBlock2,0,stream1>>>()
            !periodic along x 
            call pbc_side_y<<<dimGridy,dimBlock2,0,stream1>>>()
            
            call pbc_edge_z<<<(nz+TILE_DIM-1)/TILE_DIM, TILE_DIM,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
	    endif
        call fixPeriodic_hvar(istep)
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
        !********************************close to boundary conditions no slip everywhere********************************!
#ifdef HALOSHARED  
        call streamcoll_shared_halo<<<dimGrid,dimBlockshared,mshared,stream1>>>()
#else
        call streamcoll_shared<<<dimGrid,dimBlock,mshared,stream1>>>()
#endif
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after streamcoll_bulk', istep)

 
        !***********************************correct pressor*********
#ifndef PRESSCORR
        call correct_pressure<<<dimGrid,dimBlock,0,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after correct_pressure', istep)
#endif        
        !flop
        istep=step+1
        !******************************************call other bcs************************
        if(lpbc)then      
            !periodic along x 
            call pbc_side_x<<<dimGridx,dimBlock2,0,stream1>>>()
            !periodic along x 
            call pbc_side_y<<<dimGridy,dimBlock2,0,stream1>>>()
            call pbc_edge_z<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
	    endif
        call fixPeriodic_hvar(istep)
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
        !********************************close to boundary conditions no slip everywhere********************************!   
#ifdef HALOSHARED  
        call streamcoll_shared_halo_flop<<<dimGrid,dimBlockshared,mshared,stream1>>>()
#else
        call streamcoll_shared_flop<<<dimGrid,dimBlock,mshared,stream1>>>()
#endif
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after streamcoll_bulk_flop', istep)


        !***********************************correct pressor*********
#ifndef PRESSCORR
        call correct_pressure<<<dimGrid,dimBlock,0,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after correct_pressure_flop', istep)
#endif
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

    contains
    
    function coordblock(idblock)
    !return block coordinate from id block (idblock start from 1 so we apply minus 1)
     implicit none
     integer, intent(in) :: idblock
     integer, dimension(3) :: coordblock
   
     coordblock(3)=(idblock-1)/nxyblock +1
     coordblock(2)=((idblock-1)-(coordblock(3)-1)*nxyblock)/nxblock +1
     coordblock(1)=(idblock-1)-(coordblock(3)-1)*nxyblock-(coordblock(2)-1)*nxblock +1
      
    end function coordblock
    
    subroutine fixPeriodic_hvar(mystep)

      implicit none
      
      integer, intent(in) :: mystep
    
   
         if(xperiodic)then
          call bc_per_x_hvar<<<dimGridx, dimBlock2>>>(mystep)
          call abortOnLastErrorAndSync('after bc_per_x_hvar', mystep)
         endif
!         if(yperiodic)then
!          call bc_per_y_hvar<<<dimGridy, dimBlock2>>>(mystep)
!          call abortOnLastErrorAndSync('after bc_per_y_hvar', mystep)
!         if(zperiodic)then
!           call bc_per_z_hvar<<<dimGridz, dimBlock2>>>(mystep)
!           call abortOnLastErrorAndSync('after bc_per_z_hvar', mystep)
!         endif
         
!         if(xperiodic)then
!           call bc_edge_x_hvar<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(mystep)
!           call abortOnLastErrorAndSync('after bc_edge_x_hvar', mystep)
!         endif
!         if(yperiodic)then
!           call bc_edge_y_hvar<<<(ny+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(mystep)
!           call abortOnLastErrorAndSync('after bc_edge_y_hvar', mystep)
!         endif
!         if(zperiodic)then
!           call bc_edge_z_hvar<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(mystep)
!           call abortOnLastErrorAndSync('after bc_edge_z_hvar', mystep)
!         endif

    end subroutine fixPeriodic_hvar
    
    
end program
