program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    integer, parameter :: db=4 !kind(1.0)
    integer(kind=8) :: i,j,k
    integer(kind=8) :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus,nsp,nsk,dum
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp,dummy,dummy2
    
    integer(kind=4), allocatable,  dimension(:,:,:)   :: isfluid
    
    real(kind=db) :: rho,u,v,w
     real(kind=db), allocatable, dimension(:,:,:) :: f0
    real(kind=db), allocatable, dimension(:,:,:,:) :: f1,f2,f3,f4,f5,f6,f7,f8,f9
    real(kind=db), allocatable, dimension(:,:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18

       
    nlinks=18 !pari!
    tau=1.5_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif

    !*******************************user parameters**************************
    nx=300
    ny=300
    nz=300
    nsteps=1000
    stamp=1000
    fx=1.0_db*10.0**(-7)
    fy=0.0_db*10.0**(-5)
    fz=0.0_db*10.0**(-5)

    allocate(f0(0:nx+1,0:ny+1,0:nz+1),f1(0:nx+1,0:ny+1,0:nz+1,2),f2(0:nx+1,0:ny+1,0:nz+1,2),f3(0:nx+1,0:ny+1,0:nz+1,2))
    allocate(f4(0:nx+1,0:ny+1,0:nz+1,2),f5(0:nx+1,0:ny+1,0:nz+1,2),f6(0:nx+1,0:ny+1,0:nz+1,2),f7(0:nx+1,0:ny+1,0:nz+1,2))
    allocate(f8(0:nx+1,0:ny+1,0:nz+1,2),f9(0:nx+1,0:ny+1,0:nz+1,2),f10(0:nx+1,0:ny+1,0:nz+1,2),f11(0:nx+1,0:ny+1,0:nz+1,2))
    allocate(f12(0:nx+1,0:ny+1,0:nz+1,2),f13(0:nx+1,0:ny+1,0:nz+1,2),f14(0:nx+1,0:ny+1,0:nz+1,2),f15(0:nx+1,0:ny+1,0:nz+1,2))
    allocate(f16(0:nx+1,0:ny+1,0:nz+1,2),f17(0:nx+1,0:ny+1,0:nz+1,2),f18(0:nx+1,0:ny+1,0:nz+1,2))
    allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
    
    !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
    nsp=1 ! flip-flop
    nsk=2

    p0=(1.0_db/3.0_db)
	p1=(1.0_db/18.0_db)
	p2=(1.0_db/36.0_db)
    p1dcssq=p1/cssq
    p2dcssq=p2/cssq
    omega=1.0_db/tau
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:,:)=0 !left
    isfluid(nx,:,:)=0 !right
    isfluid(:,1,:)=0 !front 
    isfluid(:,ny,:)=0 !rear
    isfluid(:,:,1)=0 !bottom
    isfluid(:,:,nz)=0 !top
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    w=0.0_db
    rho=0.0_db     !is to be intended as a delta rho
    !do ll=0,nlinks
    f0(1:nx,1:ny,1:nz)=0.0_db
    f1(1:nx,1:ny,1:nz,:)=0.0_db
    f2(1:nx,1:ny,1:nz,:)=0.0_db
    f3(1:nx,1:ny,1:nz,:)=0.0_db
    f4(1:nx,1:ny,1:nz,:)=0.0_db
    f5(1:nx,1:ny,1:nz,:)=0.0_db
    f6(1:nx,1:ny,1:nz,:)=0.0_db
    f7(1:nx,1:ny,1:nz,:)=0.0_db
    f8(1:nx,1:ny,1:nz,:)=0.0_db
    f9(1:nx,1:ny,1:nz,:)=0.0_db
    f10(1:nx,1:ny,1:nz,:)=0.0_db
    f11(1:nx,1:ny,1:nz,:)=0.0_db
    f12(1:nx,1:ny,1:nz,:)=0.0_db
    f13(1:nx,1:ny,1:nz,:)=0.0_db
    f14(1:nx,1:ny,1:nz,:)=0.0_db
    f15(1:nx,1:ny,1:nz,:)=0.0_db
    f16(1:nx,1:ny,1:nz,:)=0.0_db
    f17(1:nx,1:ny,1:nz,:)=0.0_db
    f18(1:nx,1:ny,1:nz,:)=0.0_db
    !enddo
    !*************************************check data ************************ 
    write(6,*) '*******************LB data*****************'
    write(6,*) 'tau',tau
    write(6,*) 'omega',omega
    write(6,*) 'visc',visc_LB
    write(6,*) 'fx',fx
    write(6,*) 'fy',fy
    write(6,*) 'fz',fz
    write(6,*) 'cssq',cssq
    write(6,*) '*******************INPUT data*****************'
    write(6,*) 'nx',nx
    write(6,*) 'ny',ny
    write(6,*) 'ny',nz
    write(6,*) 'nsteps',nsteps
    write(6,*) 'stamp',stamp
    write(6,*) 'max fx',huge(fx)
    write(6,*) 'max fx',huge(fy)
    write(6,*) 'max fx',huge(fz)
    write(6,*) '*******************************************'

    !$acc data copy(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,isfluid,nsp,nsk,p0,p1,p2)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moments collision&forcing************************ 
        !!$acc update host(rho,u,v)
        !$acc update device(nsp,nsk)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) !!private(uu,temp,udotc,feq,dummy)
        !$acc loop private(uu,temp,udotc,feq,dummy,u,v,w,rho)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
                        rho = f0(i,j,k)+f1(i,j,k,nsp)+f2(i,j,k,nsp)+f3(i,j,k,nsp)+f4(i,j,k,nsp)+f5(i,j,k,nsp) &
                            +f6(i,j,k,nsp)+f7(i,j,k,nsp)+f8(i,j,k,nsp)+f9(i,j,k,nsp)+f10(i,j,k,nsp)+f11(i,j,k,nsp) &
                            +f12(i,j,k,nsp)+f13(i,j,k,nsp)+f14(i,j,k,nsp)+f15(i,j,k,nsp)+f16(i,j,k,nsp)+f17(i,j,k,nsp) &
                            +f18(i,j,k,nsp)

                        u = (f1(i,j,k,nsp)+f7(i,j,k,nsp)+f9(i,j,k,nsp)+f15(i,j,k,nsp)+f18(i,j,k,nsp)) &
                             -(f2(i,j,k,nsp)+f8(i,j,k,nsp)+f10(i,j,k,nsp)+f16(i,j,k,nsp)+f17(i,j,k,nsp)) 
                        
                        v = (f3(i,j,k,nsp)+f7(i,j,k,nsp)+f10(i,j,k,nsp)+f11(i,j,k,nsp)+f13(i,j,k,nsp)) &
                            -(f4(i,j,k,nsp)+f8(i,j,k,nsp)+f9(i,j,k,nsp)+f12(i,j,k,nsp)+f14(i,j,k,nsp))

                        w = (f5(i,j,k,nsp)+f11(i,j,k,nsp)+f14(i,j,k,nsp)+f15(i,j,k,nsp)+f17(i,j,k,nsp)) &
                            -(f6(i,j,k,nsp)+f12(i,j,k,nsp)+f13(i,j,k,nsp)+f16(i,j,k,nsp)+f18(i,j,k,nsp))
                        !
                        uu=0.5_db*(u*u + v*v + w*w)/cssq
                        !0
                        feq=p0*(rho-uu)
                        f0(i,j,k)=f0(i,j,k) + omega*(feq - f0(i,j,k)) 
                        
                        !1
                        udotc=u/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho+(temp + udotc))
                        f1(i,j,k,nsp)=f1(i,j,k,nsp) + omega*(feq - f1(i,j,k,nsp)) + fx*p1dcssq
                        
                        !2
                        feq=p1*(rho+(temp - udotc))
                        f2(i,j,k,nsp)=f2(i,j,k,nsp) + omega*(feq - f2(i,j,k,nsp)) - fx*p1dcssq
                        
                        !3
                        udotc=v/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho+(temp + udotc))
                        f3(i,j,k,nsp)=f3(i,j,k,nsp) + omega*(feq - f3(i,j,k,nsp)) + fy*p1dcssq
                        
                        !4
                        feq=p1*(rho+(temp - udotc))
                        f4(i,j,k,nsp)=f4(i,j,k,nsp) + omega*(feq - f4(i,j,k,nsp)) - fy*p1dcssq
                        
                        !7
                        udotc=(u+v)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f7(i,j,k,nsp)=f7(i,j,k,nsp) + omega*(feq - f7(i,j,k,nsp)) + (fx+fy)*p2dcssq 
                        
                        !8
                        feq=p2*(rho+(temp - udotc))
                        f8(i,j,k,nsp)=f8(i,j,k,nsp) + omega*(feq - f8(i,j,k,nsp)) - (fx+fy)*p2dcssq
                        
                        !10
                        udotc=(-u+v)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f10(i,j,k,nsp)=f10(i,j,k,nsp) + omega*(feq - f10(i,j,k,nsp)) +(fy-fx)*p2dcssq
                        
                        !9
                        feq=p2*(rho+(temp - udotc))
                        f9(i,j,k,nsp)=f9(i,j,k,nsp) + omega*(feq - f9(i,j,k,nsp)) + (fx-fy)*p2dcssq

                        !5
                        udotc=w/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho+(temp + udotc))
                        f5(i,j,k,nsp)=f5(i,j,k,nsp) + omega*(feq - f5(i,j,k,nsp)) + fz*p1dcssq
                        
                        !6
                        feq=p1*(rho+(temp - udotc))
                        f6(i,j,k,nsp)=f6(i,j,k,nsp) + omega*(feq - f6(i,j,k,nsp)) - fz*p1dcssq

                        !15
                        udotc=(u+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f15(i,j,k,nsp)=f15(i,j,k,nsp) + omega*(feq - f15(i,j,k,nsp)) + (fx+fz)*p2dcssq 
                        
                        !16
                        feq=p2*(rho+(temp - udotc))
                        f16(i,j,k,nsp)=f16(i,j,k,nsp) + omega*(feq - f16(i,j,k,nsp)) - (fx+fz)*p2dcssq

                        !17
                        udotc=(-u+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f17(i,j,k,nsp)=f17(i,j,k,nsp) + omega*(feq - f17(i,j,k,nsp)) +(fz-fx)*p2dcssq
                        
                        !18
                        feq=p2*(rho+(temp - udotc))
                        f18(i,j,k,nsp)=f18(i,j,k,nsp) + omega*(feq - f18(i,j,k,nsp)) + (fx-fz)*p2dcssq

                        !11
                        udotc=(v+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f11(i,j,k,nsp)=f11(i,j,k,nsp) + omega*(feq - f11(i,j,k,nsp)) +(fy+fz)*p2dcssq
                        
                        !12
                        feq=p2*(rho+(temp - udotc))
                        f12(i,j,k,nsp)=f12(i,j,k,nsp) + omega*(feq - f12(i,j,k,nsp)) - (fy+fz)*p2dcssq

                        !13
                        udotc=(v-w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f13(i,j,k,nsp)=f13(i,j,k,nsp) + omega*(feq - f13(i,j,k,nsp)) + (fy-fz)*p2dcssq
                        
                        !14
                        feq=p2*(rho+(temp - udotc))
                        f14(i,j,k,nsp)=f14(i,j,k,nsp) + omega*(feq - f14(i,j,k,nsp)) + (fz-fy)*p2dcssq

                    elseif(isfluid(i,j,k).eq.0)then
                        !
                        dummy=f1(i,j,k,nsp)
                        f1(i,j,k,nsp)=f2(i,j,k,nsp)
                        f2(i,j,k,nsp)=dummy
                        !
                        dummy=f3(i,j,k,nsp)
                        f3(i,j,k,nsp)=f4(i,j,k,nsp)
                        f4(i,j,k,nsp)=dummy
                        !
                        dummy=f5(i,j,k,nsp)
                        f5(i,j,k,nsp)=f6(i,j,k,nsp)
                        f6(i,j,k,nsp)=dummy
                        !
                        dummy=f7(i,j,k,nsp)
                        f7(i,j,k,nsp)=f8(i,j,k,nsp)
                        f8(i,j,k,nsp)=dummy
                        !
                        dummy=f9(i,j,k,nsp)
                        f9(i,j,k,nsp)=f10(i,j,k,nsp)
                        f10(i,j,k,nsp)=dummy
                        !
                        dummy=f11(i,j,k,nsp)
                        f11(i,j,k,nsp)=f12(i,j,k,nsp)
                        f12(i,j,k,nsp)=dummy
                        !
                        dummy=f13(i,j,k,nsp)
                        f13(i,j,k,nsp)=f14(i,j,k,nsp)
                        f14(i,j,k,nsp)=dummy
                        !
                        dummy=f15(i,j,k,nsp)
                        f15(i,j,k,nsp)=f16(i,j,k,nsp)
                        f16(i,j,k,nsp)=dummy
                        !
                        dummy=f17(i,j,k,nsp)
                        f17(i,j,k,nsp)=f18(i,j,k,nsp)
                        f18(i,j,k,nsp)=dummy
                        
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        !******************************************call bcs************************
        !periodic along y
        !x=1     
        !$acc kernels 
        f1(1,:,:,:)=f1(nx-1,:,:,:)
        !!$acc end kernels

        !!$acc kernels 
        f7(1,:,:,:)=f7(nx-1,:,:,:)
        !!$acc end kernels

        !!$acc kernels 
        f9(1,:,:,:)=f9(nx-1,:,:,:)
        !!$acc end kernels

        !!$acc kernels 
        f15(1,:,:,:)=f15(nx-1,:,:,:)
        !!$acc end kernels
        !!$acc kernels 
        f18(1,:,:,:)=f18(nx-1,:,:,:)
        !!$acc end kernels
        !
        !x=nx 
        !!$acc kernels
        f2(nx,:,:,:)=f2(2,:,:,:)
        !!$acc end kernels

        !!$acc kernels
        f8(nx,:,:,:)=f8(2,:,:,:)
        !!$acc end kernels

        !!$acc kernels
        f10(nx,:,:,:)=f10(2,:,:,:)
        !!$acc end kernels

        !!$acc kernels
        f16(nx,:,:,:)=f16(2,:,:,:)
        !!$acc end kernels

        !!$acc kernels
        f17(nx,:,:,:)=f17(2,:,:,:)
        !!$acc end kernels

        !y=1
        !!$acc kernels
        f3(:,1,:,:)=f3(:,ny-1,:,:)
        !!$acc end kernels

        !!$acc kernels
        f7(:,1,:,:)=f7(:,ny-1,:,:)
        !!$acc end kernels
        !!$acc kernels
        f10(:,1,:,:)=f10(:,ny-1,:,:)
        !!$acc end kernels

        !!$acc kernels
        f11(:,1,:,:)=f11(:,ny-1,:,:)
        !!$acc end kernels

        !!$acc kernels
        f13(:,1,:,:)=f13(:,ny-1,:,:)
        !!$acc end kernels

        !y=ny
        !!$acc kernels
        f4(:,ny,:,:)=f4(:,2,:,:)
        !!$acc end kernels

        !!$acc kernels
        f8(:,ny,:,:)=f8(:,2,:,:)
        !!$acc end kernels

        !!$acc kernels
        f9(:,ny,:,:)=f9(:,2,:,:)
        !!$acc end kernels

        !!$acc kernels
        f12(:,ny,:,:)=f12(:,2,:,:)
        !!$acc end kernels

        !!$acc kernels
        f14(:,ny,:,:)=f14(:,2,:,:)
        !$acc end kernels
        ! !******************************************streaming***************************  
        !
        !!$acc update host(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !!$acc update device(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) async(3) 
        !$acc loop independent collapse(3)  !independent  
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    f1(i,j,k,nsk)=f1(i-1,j,k,nsp)
                    f2(i,j,k,nsk)=f2(i+1,j,k,nsp)
                    f3(i,j,k,nsk)=f3(i,j-1,k,nsp)
                    f4(i,j,k,nsk)=f4(i,j+1,k,nsp)
                    f5(i,j,k,nsk)=f5(i,j,k-1,nsp)
                    f6(i,j,k,nsk)=f6(i,j,k+1,nsp)
                    f7(i,j,k,nsk)=f7(i-1,j-1,k,nsp)
                    f8(i,j,k,nsk)=f8(i+1,j+1,k,nsp)
                    f9(i,j,k,nsk)=f9(i-1,j+1,k,nsp)
                    f10(i,j,k,nsk)=f10(i+1,j-1,k,nsp)
                    f11(i,j,k,nsk)=f11(i,j-1,k-1,nsp)
                    f12(i,j,k,nsk)=f12(i,j+1,k+1,nsp)
                    f13(i,j,k,nsk)=f13(i,j-1,k+1,nsp)
                    f14(i,j,k,nsk)=f14(i,j+1,k-1,nsp)
                    f15(i,j,k,nsk)=f15(i-1,j,k-1,nsp)
                    f16(i,j,k,nsk)=f16(i+1,j,k+1,nsp)
                    f17(i,j,k,nsk)=f17(i+1,j,k-1,nsp)
                    f18(i,j,k,nsk)=f18(i-1,j,k+1,nsp)
                enddo
            enddo
        enddo
        !$acc end kernels
        !flip-flop
        dum=nsp
        nsp=nsk
        nsk=dum
        !
    enddo 
    call cpu_time(ts2)
    !$acc update host(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,nsk,nsp)
    !$acc end data
    rho = f0(nx/2,ny/2,nz/2) &
                            +f1(nx/2,ny/2,nz/2,nsk) &
                            +f2(nx/2,ny/2,nz/2,nsk) &
                            +f3(nx/2,ny/2,nz/2,nsk) &
                            +f4(nx/2,ny/2,nz/2,nsk) &
                            +f5(nx/2,ny/2,nz/2,nsk) &
                            +f6(nx/2,ny/2,nz/2,nsk) &
                            +f7(nx/2,ny/2,nz/2,nsk) &
                            +f8(nx/2,ny/2,nz/2,nsk) &
                            +f9(nx/2,ny/2,nz/2,nsk) &
                            +f10(nx/2,ny/2,nz/2,nsk) &
                            +f11(nx/2,ny/2,nz/2,nsk) &
                            +f12(nx/2,ny/2,nz/2,nsk) &
                            +f13(nx/2,ny/2,nz/2,nsk) &
                            +f14(nx/2,ny/2,nz/2,nsk) &
                            +f15(nx/2,ny/2,nz/2,nsk) &
                            +f16(nx/2,ny/2,nz/2,nsk) &
                            +f17(nx/2,ny/2,nz/2,nsk) &
                            +f18(nx/2,ny/2,nz/2,nsk)

    u = (f1(nx/2,ny/2,nz/2,nsk)+f7(nx/2,ny/2,nz/2,nsk)+f9(nx/2,ny/2,nz/2,nsk)+f15(nx/2,ny/2,nz/2,nsk)+f18(nx/2,ny/2,nz/2,nsk)) &
            -(f2(nx/2,ny/2,nz/2,nsk)+f8(nx/2,ny/2,nz/2,nsk)+f10(nx/2,ny/2,nz/2,nsk)+f16(nx/2,ny/2,nz/2,nsk)+f17(nx/2,ny/2,nz/2,nsk)) 
    
    v = (f3(nx/2,ny/2,nz/2,nsk)+f7(nx/2,ny/2,nz/2,nsk)+f10(nx/2,ny/2,nz/2,nsk)+f11(nx/2,ny/2,nz/2,nsk)+f13(nx/2,ny/2,nz/2,nsk)) &
        -(f4(nx/2,ny/2,nz/2,nsk)+f8(nx/2,ny/2,nz/2,nsk)+f9(nx/2,ny/2,nz/2,nsk)+f12(nx/2,ny/2,nz/2,nsk)+f14(nx/2,ny/2,nz/2,nsk))

    w = (f5(nx/2,ny/2,nz/2,nsk)+f11(nx/2,ny/2,nz/2,nsk)+f14(nx/2,ny/2,nz/2,nsk)+f15(nx/2,ny/2,nz/2,nsk)+f17(nx/2,ny/2,nz/2,nsk)) &
        -(f6(nx/2,ny/2,nz/2,nsk)+f12(nx/2,ny/2,nz/2,nsk)+f13(nx/2,ny/2,nz/2,nsk)+f16(nx/2,ny/2,nz/2,nsk)+f18(nx/2,ny/2,nz/2,nsk))

    write(6,*) 'u=',u,'v=',v,'w=',w,'rho=',rho,nsk,nsp

    write(6,*) 'You''ve just wasted ', ts2-ts1, ' s of your life time' 

    
end program