program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    !***************************************var block********************************!
        integer, parameter :: db=4 !kind(1.0)
        integer(kind=8) :: i,j,k
        integer(kind=8) :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus,nsp,nsk,dum
        
        real(kind=db),parameter :: pi_greek=3.14159265359793234626433
        
        real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
        real(kind=db) :: visc_LB,uu,udotc,omega,feq
        real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp
        
        integer(kind=4), allocatable,  dimension(:,:,:)   :: isfluid
        
        real(kind=db) :: rho,u,v,w
        real(kind=db), allocatable, dimension(:,:,:) :: f0
        real(kind=db), allocatable, dimension(:,:,:,:) :: f1,f2,f3,f4,f5,f6,f7,f8,f9
        real(kind=db), allocatable, dimension(:,:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18

    !***********************************sim pars***************************************! 
        nlinks=18 !pari!
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
        ny=256
        nz=256
        nsteps=10
        stamp=1000
        fx=0.0_db*10.0**(-7)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
    !***********************************************allocation and lattice vars************************
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
        rho=1.0_db  !not to be intended as a delta rho
    !init distros
        f0(1:nx,1:ny,1:nz)=p0
        f1(1:nx,1:ny,1:nz,:)=p1
        f2(1:nx,1:ny,1:nz,:)=p1
        f3(1:nx,1:ny,1:nz,:)=p1
        f4(1:nx,1:ny,1:nz,:)=p1
        f5(1:nx,1:ny,1:nz,:)=p1
        f6(1:nx,1:ny,1:nz,:)=p1
        f7(1:nx,1:ny,1:nz,:)=p2
        f8(1:nx,1:ny,1:nz,:)=p2
        f9(1:nx,1:ny,1:nz,:)=p2
        f10(1:nx,1:ny,1:nz,:)=p2
        f11(1:nx,1:ny,1:nz,:)=p2
        f12(1:nx,1:ny,1:nz,:)=p2
        f13(1:nx,1:ny,1:nz,:)=p2
        f14(1:nx,1:ny,1:nz,:)=p2
        f15(1:nx,1:ny,1:nz,:)=p2
        f16(1:nx,1:ny,1:nz,:)=p2
        f17(1:nx,1:ny,1:nz,:)=p2
        f18(1:nx,1:ny,1:nz,:)=p2
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
        !***********************************moments collision bbck + forcing************************ 
        !$acc update device(nsp,nsk)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) 
        !$acc loop independent collapse (3) private(uu,temp,udotc,u,v,w,rho)
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
                        f1(i+1,j,k,nsk)=f1(i,j,k,nsp) + omega*(feq - f1(i,j,k,nsp)) + fx*p1dcssq
                        
                        !2
                        feq=p1*(rho+(temp - udotc))
                        f2(i-1,j,k,nsk)=f2(i,j,k,nsp) + omega*(feq - f2(i,j,k,nsp)) - fx*p1dcssq
                        
                        !3
                        udotc=v/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho+(temp + udotc))
                        f3(i,j+1,k,nsk)=f3(i,j,k,nsp) + omega*(feq - f3(i,j,k,nsp)) + fy*p1dcssq
                        
                        !4
                        feq=p1*(rho+(temp - udotc))
                        f4(i,j-1,k,nsk)=f4(i,j,k,nsp) + omega*(feq - f4(i,j,k,nsp)) - fy*p1dcssq
                        
                        !7
                        udotc=(u+v)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f7(i+1,j+1,k,nsk)=f7(i,j,k,nsp) + omega*(feq - f7(i,j,k,nsp)) + (fx+fy)*p2dcssq 
                        
                        !8
                        feq=p2*(rho+(temp - udotc))
                        f8(i-1,j-1,k,nsk)=f8(i,j,k,nsp) + omega*(feq - f8(i,j,k,nsp)) - (fx+fy)*p2dcssq
                        
                        !10
                        udotc=(-u+v)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f10(i-1,j+1,k,nsk)=f10(i,j,k,nsp) + omega*(feq - f10(i,j,k,nsp)) +(fy-fx)*p2dcssq
                        
                        !9
                        feq=p2*(rho+(temp - udotc))
                        f9(i+1,j-1,k,nsk)=f9(i,j,k,nsp) + omega*(feq - f9(i,j,k,nsp)) + (fx-fy)*p2dcssq

                        !5
                        udotc=w/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho+(temp + udotc))
                        f5(i,j,k+1,nsk)=f5(i,j,k,nsp) + omega*(feq - f5(i,j,k,nsp)) + fz*p1dcssq
                        
                        !6
                        feq=p1*(rho+(temp - udotc))
                        f6(i,j,k-1,nsk)=f6(i,j,k,nsp) + omega*(feq - f6(i,j,k,nsp)) - fz*p1dcssq

                        !15
                        udotc=(u+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f15(i+1,j,k+1,nsk)=f15(i,j,k,nsp) + omega*(feq - f15(i,j,k,nsp)) + (fx+fz)*p2dcssq 
                        
                        !16
                        feq=p2*(rho+(temp - udotc))
                        f16(i-1,j,k-1,nsk)=f16(i,j,k,nsp) + omega*(feq - f16(i,j,k,nsp)) - (fx+fz)*p2dcssq

                        !17
                        udotc=(-u+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f17(i-1,j,k+1,nsk)=f17(i,j,k,nsp) + omega*(feq - f17(i,j,k,nsp)) +(fz-fx)*p2dcssq
                        
                        !18
                        feq=p2*(rho+(temp - udotc))
                        f18(i+1,j,k-1,nsk)=f18(i,j,k,nsp) + omega*(feq - f18(i,j,k,nsp)) + (fx-fz)*p2dcssq

                        !11
                        udotc=(v+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f11(i,j+1,k+1,nsk)=f11(i,j,k,nsp) + omega*(feq - f11(i,j,k,nsp)) +(fy+fz)*p2dcssq
                        
                        !12
                        feq=p2*(rho+(temp - udotc))
                        f12(i,j-1,k-1,nsk)=f12(i,j,k,nsp) + omega*(feq - f12(i,j,k,nsp)) - (fy+fz)*p2dcssq

                        !13
                        udotc=(v-w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho+(temp + udotc))
                        f13(i,j+1,k-1,nsk)=f13(i,j,k,nsp) + omega*(feq - f13(i,j,k,nsp)) + (fy-fz)*p2dcssq
                        
                        !14
                        feq=p2*(rho+(temp - udotc))
                        f14(i,j-1,k+1,nsk)=f14(i,j,k,nsp) + omega*(feq - f14(i,j,k,nsp)) + (fz-fy)*p2dcssq

                    elseif(isfluid(i,j,k).eq.0)then
                        !
                        f1(i+1,j,k,nsk)=f2(i,j,k,nsp)
                        f2(i-1,j,k,nsk)=f1(i,j,k,nsp)
                        f3(i,j+1,k,nsk)=f4(i,j,k,nsp)
                        f4(i,j-1,k,nsk)=f3(i,j,k,nsp)
                        f5(i,j,k+1,nsk)=f6(i,j,k,nsp)
                        f6(i,j,k-1,nsk)=f5(i,j,k,nsp)
                        f7(i+1,j+1,k,nsk)=f8(i,j,k,nsp)
                        f8(i-1,j-1,k,nsk)=f7(i,j,k,nsp)
                        f9(i+1,j-1,k,nsk)=f10(i,j,k,nsp)
                        f10(i-1,j+1,k,nsk)=f9(i,j,k,nsp)
                        f11(i,j+1,k+1,nsk)=f12(i,j,k,nsp)
                        f12(i,j-1,k-1,nsk)=f11(i,j,k,nsp)
                        f13(i,j+1,k-1,nsk)=f14(i,j,k,nsp)
                        f14(i,j-1,k+1,nsk)=f13(i,j,k,nsp)
                        f15(i+1,j,k+1,nsk)=f16(i,j,k,nsp)
                        f16(i-1,j,k-1,nsk)=f15(i,j,k,nsp)
                        f17(i-1,j,k+1,nsk)=f18(i,j,k,nsp)
                        f18(i+1,j,k-1,nsk)=f17(i,j,k,nsp)
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        !******************************************call periodic (or other) bcs************************
            !periodic along y
            !x=1     
            !$acc kernels 
            f1(2,:,:,:)=f1(nx,:,:,:)
        
            f7(2,:,:,:)=f7(nx,:,:,:)
        
            f9(2,:,:,:)=f9(nx,:,:,:)
        
            f15(2,:,:,:)=f15(nx,:,:,:)
        
            f18(2,:,:,:)=f18(nx,:,:,:)

            !x=nx 
            f2(nx-1,:,:,:)=f2(1,:,:,:)
        
            f8(nx-1,:,:,:)=f8(1,:,:,:)
            
            f10(nx-1,:,:,:)=f10(1,:,:,:)
        
            f16(nx-1,:,:,:)=f16(1,:,:,:)
        
            f17(nx-1,:,:,:)=f17(1,:,:,:)

            !y=1
            f3(:,2,:,:)=f3(:,ny,:,:)
        
            f7(:,2,:,:)=f7(:,ny,:,:)
        
            f10(:,2,:,:)=f10(:,ny,:,:)
            
            f11(:,2,:,:)=f11(:,ny,:,:)
        
            f13(:,2,:,:)=f13(:,ny,:,:)
        
            !y=ny
            f4(:,ny-1,:,:)=f4(:,1,:,:)

            f8(:,ny-1,:,:)=f8(:,1,:,:)
    
            f9(:,ny-1,:,:)=f9(:,1,:,:)
        
            f12(:,ny-1,:,:)=f12(:,1,:,:)
        
            f14(:,ny-1,:,:)=f14(:,1,:,:)
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

    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  nx*ny*nz*nsteps/10.0_db**9/ts2-ts1

    
end program