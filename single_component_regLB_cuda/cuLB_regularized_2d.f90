 
 module mysubs
   
   use cudafor
   
   integer, parameter :: db=4 !kind(1.0)
   ! device arrays
    integer(kind=4), allocatable,  dimension(:,:), device   :: isfluid_d
    integer, constant :: nx_d,ny_d,TILE_DIMx_d,TILE_DIMy_d
    integer, parameter :: nz=1,nz_d=1
    real(kind=db), dimension(0:8), constant :: p_d
    real(kind=db), constant :: fx_d,fy_d,omega_d,qxx_d,qyy_d,qxy5_7_d,qxy6_8_d, &
     pi2cssq0_d,pi2cssq1_d,pi2cssq2_d,myrho_d,myu_d,myv_d
    real(kind=db), allocatable, dimension(:,:), device  :: rho_d,u_d,v_d,pxx_d,pyy_d,pxy_d
    real(kind=db), allocatable, dimension(:,:), device  :: f0_d,f1_d,f2_d,f3_d,f4_d,f5_d,f6_d,f7_d,f8_d
    real(kind=db), allocatable, dimension(:,:,:), device :: rhoprint_d
    real(kind=db), allocatable, dimension(:,:,:,:), device :: velprint_d
    type (dim3) :: dimGrid,dimBlock
   
  contains
  
    attributes(global) subroutine setup_pops()
      
      integer :: i,j
    
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      
      !write(*,*)i,j,p_d(0)*myrho_d
      
      f0_d(i,j)=p_d(0)*myrho_d
      f1_d(i,j)=p_d(1)*myrho_d
      f2_d(i,j)=p_d(2)*myrho_d
      f3_d(i,j)=p_d(3)*myrho_d
      f4_d(i,j)=p_d(4)*myrho_d
      f5_d(i,j)=p_d(5)*myrho_d
      f6_d(i,j)=p_d(6)*myrho_d
      f7_d(i,j)=p_d(7)*myrho_d
      f8_d(i,j)=p_d(8)*myrho_d
     

  end subroutine setup_pops
  
  attributes(global) subroutine moments()
      
      integer :: i,j
      real(kind=db) ::uu,udotc,temp,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8
    
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      
      if(isfluid_d(i,j).ne.1)return
      
       rho_d(i,j) = f0_d(i,j)+f1_d(i,j)+f2_d(i,j)+f3_d(i,j)+f4_d(i,j)+f5_d(i,j)+f6_d(i,j)+f7_d(i,j)+f8_d(i,j)
                    u_d(i,j) = (f1_d(i,j) +f5_d(i,j) +f8_d(i,j)-f3_d(i,j) -f6_d(i,j) -f7_d(i,j)) !/rho_d(i,j)
                    v_d(i,j) = (f5_d(i,j) +f2_d(i,j) +f6_d(i,j)-f7_d(i,j) -f4_d(i,j) -f8_d(i,j))
                    ! non equilibrium pressor components
                    uu=0.5_db*(u_d(i,j)*u_d(i,j) + v_d(i,j)*v_d(i,j))/cssq
                    !1-3
                    udotc=u_d(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq1=f1_d(i,j)-p_d(1)*(rho_d(i,j)+(temp + udotc))
                    fneq3=f3_d(i,j)-p_d(3)*(rho_d(i,j)+(temp - udotc))
                    !2-4
                    udotc=v_d(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq2=f2_d(i,j)-p_d(2)*(rho_d(i,j)+(temp + udotc))
                    fneq4=f4_d(i,j)-p_d(4)*(rho_d(i,j)+(temp - udotc))
                    !5-7
                    udotc=(u_d(i,j)+v_d(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq5=f5_d(i,j)-p_d(5)*(rho_d(i,j)+(temp + udotc))
                    fneq7=f7_d(i,j)-p_d(7)*(rho_d(i,j)+(temp - udotc))
                    !6-8
                    udotc=(-u_d(i,j)+v_d(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq6=f6_d(i,j)-p_d(6)*(rho_d(i,j)+(temp + udotc))
                    fneq8=f8_d(i,j)-p_d(8)*(rho_d(i,j)+(temp - udotc))

                    pxx_d(i,j)= fneq1 + fneq3 + fneq5 + fneq6 + fneq7 + fneq8
                    pyy_d(i,j)= fneq2 + fneq4 + fneq5 + fneq6 + fneq7 + fneq8
                    pxy_d(i,j)= fneq5 - fneq6 + fneq7 - fneq8
     

  end subroutine moments
  
  attributes(global) subroutine streamcoll()
      
      integer :: i,j
      real(kind=db) ::uu,udotc,temp,feq
    
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      
      if(isfluid_d(i,j).ne.1)return
      
      uu=0.5_db*(u_d(i,j)*u_d(i,j) + v_d(i,j)*v_d(i,j))/cssq
      !oneminusuu= -uu !1.0_db - uu
      !0
      feq=p_d(0)*(rho_d(i,j)-uu)
      f0_d(i,j)=feq + (1.0_db-omega_d)*pi2cssq0_d*(-cssq*pxx_d(i,j)-cssq*pyy_d(i,j))
      !1
      udotc=u_d(i,j)/cssq
      temp = -uu + 0.5_db*udotc*udotc
      feq=p_d(1)*(rho_d(i,j)+(temp + udotc))
      f1_d(i+1,j)= feq + (1.0_db-omega_d)*pi2cssq1_d*(qxx_d*pxx_d(i,j)-cssq*pyy_d(i,j)) + fx_d*p_d(1)/cssq !f1(i-1,j,nsp) + omega_d*(feq - f1(i-1,j,nsp)) + fx*p(1)/cssq
      !3
      feq=p_d(3)*(rho_d(i,j)+(temp - udotc))
      f3_d(i-1,j)= feq + (1.0_db-omega_d)*pi2cssq1_d*(qxx_d*pxx_d(i,j)-cssq*pyy_d(i,j))  - fx_d*p_d(3)/cssq !f3(i+1,j,nsp) + omega_d*(feq - f3(i+1,j,nsp)) - fx*p(3)/cssq
      !2
      udotc=v_d(i,j)/cssq
      temp = -uu + 0.5_db*udotc*udotc
      feq=p_d(2)*(rho_d(i,j)+(temp + udotc))
      f2_d(i,j+1)= feq + (1.0_db-omega_d)*pi2cssq1_d*(qyy_d*pyy_d(i,j)-cssq*pxx_d(i,j))  + fy_d*p_d(2)/cssq !f2(i,j-1,nsp) + omega_d*(feq - f2(i,j-1,nsp)) + fy*p(2)/cssq
      !4
      feq=p_d(4)*(rho_d(i,j)+(temp - udotc))
      f4_d(i,j-1)= feq + (1.0_db-omega_d)*pi2cssq1_d*(qyy_d*pyy_d(i,j)-cssq*pxx_d(i,j))  - fy_d*p_d(4)/cssq !f4(i,j+1,nsp) + omega_d*(feq - f4(i,j+1,nsp)) - fy*p(4)/cssq
      !5
      udotc=(u_d(i,j)+v_d(i,j))/cssq
      temp = -uu + 0.5_db*udotc*udotc
      feq=p_d(5)*(rho_d(i,j)+(temp + udotc))
      f5_d(i+1,j+1)= feq + (1.0_db-omega_d)*pi2cssq2_d*(qxx_d*pxx_d(i,j)+qyy_d*pyy_d(i,j)+2.0_db*qxy5_7_d*pxy_d(i,j)) + fx_d*p_d(5)/cssq + fy_d*p_d(5)/cssq!f5(i-1,j-1,nsp) + omega_d*(feq - f5(i-1,j-1,nsp)) + fx*p(5)/cssq + fy*p(5)/cssq 
      !7
      feq=p_d(7)*(rho_d(i,j)+(temp - udotc))
      f7_d(i-1,j-1)=feq + (1.0_db-omega_d)*pi2cssq2_d*(qxx_d*pxx_d(i,j)+qyy_d*pyy_d(i,j)+2.0_db*qxy5_7_d*pxy_d(i,j)) - fx_d*p_d(7)/cssq - fy_d*p_d(7)/cssq !f7(i+1,j+1,nsp) + omega_d*(feq - f7(i+1,j+1,nsp)) - fx*p(7)/cssq - fy*p(7)/cssq
      !6
      udotc=(-u_d(i,j)+v_d(i,j))/cssq
      temp = -uu + 0.5_db*udotc*udotc
      feq=p_d(6)*(rho_d(i,j)+(temp + udotc))
      f6_d(i-1,j+1)= feq + (1.0_db-omega_d)*pi2cssq2_d*(qxx_d*pxx_d(i,j)+qyy_d*pyy_d(i,j)+2.0_db*qxy6_8_d*pxy_d(i,j)) - fx_d*p_d(6)/cssq + fy_d*p_d(6)/cssq !f6(i+1,j-1,nsp) + omega_d*(feq - f6(i+1,j-1,nsp)) - fx*p(6)/cssq + fy*p(6)/cssq
      !8
      feq=p_d(8)*(rho_d(i,j)+(temp - udotc))
      f8_d(i+1,j-1)=feq + (1.0_db-omega_d)*pi2cssq2_d*(qxx_d*pxx_d(i,j)+qyy_d*pyy_d(i,j)+2.0_db*qxy6_8_d*pxy_d(i,j)) + fx_d*p_d(8)/cssq - fy_d*p_d(8)/cssq !f8(i-1


  end subroutine streamcoll
  
  attributes(global) subroutine bcs_no_slip()
      
      integer :: i,j
    
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      
      !write(*,*)i,j,p_d(0)*myrho_d
      if(isfluid_d(i,j).ne.0)return
      
      

	  f8_d(i+1,j-1)=f6_d(i,j)!gpc 
	  f7_d(i-1,j-1)=f5_d(i,j)!hpc

	  f6_d(i-1,j+1)=f8_d(i,j)!gpc 
	  f5_d(i+1,j+1)=f7_d(i,j)!hpc 


	  f4_d(i,j-1)=f2_d(i,j)!gpc 
	  f3_d(i-1,j)=f1_d(i,j)!hpc 

	  f2_d(i,j+1)=f4_d(i,j)!gpc 
	  f1_d(i+1,j)=f3_d(i,j)!hpc 
     

  end subroutine bcs_no_slip
  
  attributes(global) subroutine pbc_edge_y()
      
      integer :: i,j
    
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      
      if(i<2 .or. i>nx_d-1)return
      
      f5_d(i,2)=f5_d(i,ny_d)
      f2_d(i,2)=f2_d(i,ny_d)
      f6_d(i,2)=f6_d(i,ny_d)
      f8_d(i,ny_d-1)=f8_d(i,1)
      f4_d(i,ny_d-1)=f4_d(i,1)
      f7_d(i,ny_d-1)=f7_d(i,1)
     

  end subroutine pbc_edge_y
  
  attributes(global) subroutine store_print()
      
      integer :: i,j
    
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      
      !write(*,*)i,j,p_d(0)*myrho_d
      if(isfluid_d(i,j).eq.1)then
        rhoprint_d(i,j,1)=rho_d(i,j)
        velprint_d(1,i,j,1)=u_d(i,j)
        velprint_d(2,i,j,1)=v_d(i,j)
        return
      endif
      
      rhoprint_d(i,j,1)=0.0
      velprint_d(1,i,j,1)=0.0
      velprint_d(2,i,j,1)=0.0
      
      return

  end subroutine store_print

 end module mysubs

program lb_openacc
  
    use cudafor
    use mysubs
    
    implicit none
    
    
    integer(kind=4) :: i,j,ll,l,dumm
    integer(kind=4) :: nx,ny,step,stamp,nlinks,nsteps,ngpus
    integer :: TILE_DIMx,TILE_DIMy,istat,iframe
    real(kind=db),parameter :: pi_greek=3.141592653589793238462643383279502884_db
    logical :: lprint=.false.
    
    real(kind=4)  :: ts1,ts2 
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8
    real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq1,pi2cssq2,pi2cssq0
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,dummy,myrho,myu,myv
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p
    !real(kind=db), allocatable, dimension(:,:) :: rho,u,v,pxx,pyy,pxy
    !real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8
    real(kind=db), allocatable, dimension(:,:,:) :: rhoprint
    real(kind=db), allocatable, dimension(:,:,:,:) :: velprint
    
    integer, parameter :: mxln=120
    character(len=mxln) :: sevt1,sevt2
    
    sevt1=repeat(' ',mxln)
    sevt2=repeat(' ',mxln)
    
       
    nlinks=8 !pari!
    tau=1.0_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif

    !*******************************user parameters**************************
    nx=256
    ny=512
    TILE_DIMx=64
    TILE_DIMy=1
    if (mod(nx, TILE_DIMx)/= 0) then
        write(*,*) 'nx must be a multiple of TILE_DIM'
        stop
    end if
    if (mod(ny, TILE_DIMy) /= 0) then
        write(*,*) 'ny must be a multiple of TILE_DIMy'
        stop
    end if
    dimGrid  = dim3(nx/TILE_DIMx, ny/TILE_DIMy, 1)
    dimBlock = dim3(TILE_DIMx, TILE_DIMy, 1)
    
    nsteps=1000
    stamp=50
    lprint=.true.
    fx=0.0_db*10.0_db**(-7.0_db)
    fy=1.0_db*10.0_db**(-8.0_db)
    allocate(p(0:nlinks))
    !allocate(f0(0:nx+1,0:ny+1),f1(0:nx+1,0:ny+1),f2(0:nx+1,0:ny+1),f3(0:nx+1,0:ny+1),f4(0:nx+1,0:ny+1))
    !allocate(f5(0:nx+1,0:ny+1),f6(0:nx+1,0:ny+1),f7(0:nx+1,0:ny+1),f8(0:nx+1,0:ny+1))
    !allocate(rho(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny),pxx(1:nx,1:ny),pyy(1:nx,1:ny),pxy(1:nx,1:ny))
    allocate(isfluid(1:nx,1:ny)) !,omega_2d(1:nx,1:ny)) 
    
    
    !ex=(/0,1,0,-1,0,1,-1,-1,1/)
    !ey=(/0,0,1,0,-1,1,1,-1,-1/)

    p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
    1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
    
   
    omega=1.0_db/tau
    
    ! regularized: hermite 
    qxx=1.0_db-cssq
    qyy=1.0_db-cssq
    qxy5_7=1.0_db
    qxy6_8=-1.0_db
    pi2cssq0=p(0)/(2.0_db*cssq**2)
    pi2cssq1=p(1)/(2.0_db*cssq**2)
    pi2cssq2=p(5)/(2.0_db*cssq**2)
    
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:)=0 !EAST
    isfluid(nx,:)=0 !WEST
    isfluid(:,1)=0 !SOUTH 
    isfluid(:,ny)=0 !NORTH
    !*************************************initial conditions ************************    
    myu=0.0_db
    myv=0.0_db
    myrho=1.0_db     !rho!
    !do ll=0,nlinks
!    f0(1:nx,1:ny)=p(0)*rho(:,:)!0.0_db
!    f1(1:nx,1:ny)=p(1)*rho(:,:)
!    f2(1:nx,1:ny)=p(2)*rho(:,:)
!    f3(1:nx,1:ny)=p(3)*rho(:,:)
!    f4(1:nx,1:ny)=p(4)*rho(:,:)
!    f5(1:nx,1:ny)=p(5)*rho(:,:)
!    f6(1:nx,1:ny)=p(6)*rho(:,:)
!    f7(1:nx,1:ny)=p(7)*rho(:,:)
!    f8(1:nx,1:ny)=p(8)*rho(:,:)
    !enddo
    !*************************************check data ************************ 
    write(6,*) '*******************LB data*****************'
    write(6,*) 'tau',tau
    write(6,*) 'omega',omega
    write(6,*) 'visc',visc_LB
    write(6,*) 'fx',fx
    write(6,*) 'cssq',cssq
    write(6,*) '*******************INPUT data*****************'
    write(6,*) 'nx',nx
    write(6,*) 'ny',ny
    write(6,*) 'nsteps',nsteps
    write(6,*) 'stamp',stamp
    write(6,*) 'max fx',huge(fx)
    write(6,*) 'available gpus',ngpus
    write(6,*) '*******************************************'

!!!!!!!!!!!!!from host to device!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nx_d=nx
    ny_d=ny
    TILE_DIMx_d=TILE_DIMx
    TILE_DIMy_d=TILE_DIMy
    myrho_d=myrho
    myu_d=myu
    myv_d=myv
    fx_d=fx
    fy_d=fy
    p_d=p
    omega_d=omega
    qxx_d=qxx
    qyy_d=qyy
    qxy5_7_d=qxy5_7
    qxy6_8_d=qxy6_8
    pi2cssq0_d= pi2cssq0
    pi2cssq1_d=pi2cssq1
    pi2cssq2_d=pi2cssq2   
    allocate(isfluid_d(1:nx_d,1:ny_d))
    istat = cudaMemcpy(isfluid_d,isfluid,nx*ny*nz )
    if (istat/=0) write(*,*) 'status after copy isfluid:',istat
    allocate(rho_d(1:nx_d,1:ny_d),u_d(1:nx_d,1:ny_d),v_d(1:nx_d,1:ny_d),pxx_d(1:nx_d,1:ny_d),pyy_d(1:nx_d,1:ny_d),pxy_d(1:nx_d,1:ny_d))
    allocate(f0_d(0:nx_d+1,0:ny_d+1),f1_d(0:nx_d+1,0:ny_d+1),f2_d(0:nx_d+1,0:ny_d+1),f3_d(0:nx_d+1,0:ny_d+1),f4_d(0:nx_d+1,0:ny_d+1))
    allocate(f5_d(0:nx_d+1,0:ny_d+1),f6_d(0:nx_d+1,0:ny_d+1),f7_d(0:nx_d+1,0:ny_d+1),f8_d(0:nx_d+1,0:ny_d+1))
    
    call setup_pops<<<dimGrid,dimBlock>>>()
    
    if(lprint)then
      allocate(rhoprint(1:nx,1:ny,1),velprint(3,1:nx,1:ny,1))
      allocate(rhoprint_d(1:nx_d,1:ny_d,1),velprint_d(3,1:nx_d,1:ny_d,1))
    endif
    
    istat = cudaDeviceSynchronize
    iframe=0
    
    

    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moment + neq pressor*********
        
        call moments<<<dimGrid,dimBlock>>>()
        
        !***********************************PRINT************************
        if(mod(step,stamp).eq.0)then
          istat = cudaDeviceSynchronize
          call store_print<<<dimGrid,dimBlock>>>()
          istat = cudaDeviceSynchronize
          istat = cudaMemcpy(rhoprint,rhoprint_d,nx*ny*nz )
          istat = cudaMemcpy(velprint,velprint_d,3*nx*ny*nz )
          istat = cudaDeviceSynchronize
          iframe=iframe+1
          write(6,'(a,2i8)')'stamp frame : ',step,iframe
          sevt1 = 'rho'//write_fmtnumb(iframe)//'.out'
          sevt2 = 'vel'//write_fmtnumb(iframe)//'.out'
          open(unit=345,file=trim(sevt1), &
           status='replace',action='write',access='stream',form='unformatted')
          write(345)rhoprint
          close(345)
          open(unit=346,file=trim(sevt2), &
           status='replace',action='write',access='stream',form='unformatted')
          write(346)velprint
          close(346)
        endif
        
        
        !***********************************collision + no slip + forcing: fused implementation*********
        call  streamcoll<<<dimGrid,dimBlock>>>()
        
        
          
        !********************************************bcs no slip*****************************************!
        
        call bcs_no_slip<<<dimGrid,dimBlock>>>()
        
        
     
        !!$acc end kernels
        !******************************************call periodic bcs: always after fused************************
        !periodic along y
        !!$acc kernels 
        
        call pbc_edge_y<<<(nx+TILE_DIMx-1)/TILE_DIMx, TILE_DIMx>>>()
        
        !istat = cudaDeviceSynchronize
        

    enddo 
    call cpu_time(ts2)
    


    !************************************************test points**********************************************!
!    write(6,*) 'u=',u(nx/2,ny/2) ,'v=',v(nx/2,ny/2),'rho',rho(nx/2,ny/2) !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
!    write(6,*) 'u=',u(2,ny/2) ,'v=',v(2,ny/2),'rho',rho(2,ny/2)
!    write(6,*) 'u=',u(1,ny/2) ,'v=',v(1,ny/2),'rho',rho(1,ny/2)
    
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  nx*ny*nsteps/10.0_db**9/ts2-ts1

    open(101, file = 'v.out', status = 'replace')
    do j=1,ny
        do i=1,nx
            write(101,*) velprint(2,i,j,1) 
        enddo
    enddo
    close(101) 
    
    
  contains 
  !*****************************************************functions********************************************************!
    function dimenumb(inum)

    !***********************************************************************
    !    
    !     LBsoft function for returning the number of digits
    !     of an integer number
    !     originally written in JETSPIN by M. Lauricella et al.
    !    
    !     licensed under the 3-Clause BSD License (BSD-3-Clause)
    !     author: M. Lauricella
    !     last modification July 2018
    !    
    !***********************************************************************

        implicit none

        integer,intent(in) :: inum
        integer :: dimenumb
        integer :: i
        real*8 :: tmp

        i=1
        tmp=real(inum,kind=8)
        do
        if(tmp< 10.d0 )exit
        i=i+1
        tmp=tmp/ 10.0d0
        enddo

        dimenumb=i

        return

        end function dimenumb

    function write_fmtnumb(inum)

    !***********************************************************************
    !    
    !     LBsoft function for returning the string of six characters
    !     with integer digits and leading zeros to the left
    !     originally written in JETSPIN by M. Lauricella et al.
    !    
    !     licensed under the 3-Clause BSD License (BSD-3-Clause)
    !     author: M. Lauricella
    !     last modification July 2018
    !    
    !***********************************************************************

    implicit none

    integer,intent(in) :: inum
    character(len=6) :: write_fmtnumb
    integer :: numdigit,irest
    !real*8 :: tmp
    character(len=22) :: cnumberlabel
    numdigit=dimenumb(inum)
    irest=6-numdigit
    if(irest>0)then
        write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
        write(write_fmtnumb,fmt=cnumberlabel)repeat('0',irest),inum
    else
        write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
        write(write_fmtnumb,fmt=cnumberlabel)inum
    endif

    return
    end function write_fmtnumb   
    
end program

 
