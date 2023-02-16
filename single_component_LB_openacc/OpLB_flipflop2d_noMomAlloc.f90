program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    integer, parameter :: db=4 !kind(1.0)
    integer(kind=8) :: i,j,ll,l,nsp,nsk,dumm
    integer(kind=8) :: nx,ny,step,stamp,nlinks,nsteps,ngpus
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2 
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,dummy
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p
    real(kind=db)  :: rho,u,v
    real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8

       
    nlinks=8 !pari!
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
    nx=1024
    ny=2048
    nsteps=1000
    stamp=1000
    fx=0.0_db*10.0**(-7)
    fy=1.0_db*10.0**(-4)
    allocate(p(0:nlinks))
    allocate(f0(0:nx+1,0:ny+1,2),f1(0:nx+1,0:ny+1,2),f2(0:nx+1,0:ny+1,2),f3(0:nx+1,0:ny+1,2),f4(0:nx+1,0:ny+1,2))
    allocate(f5(0:nx+1,0:ny+1,2),f6(0:nx+1,0:ny+1,2),f7(0:nx+1,0:ny+1,2),f8(0:nx+1,0:ny+1,2))
    allocate(isfluid(1:nx,1:ny)) !,omega_2d(1:nx,1:ny)) 
    
    !ex=(/0,1,0,-1,0,1,-1,-1,1/)
    !ey=(/0,0,1,0,-1,1,1,-1,-1/)

    p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
    1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
    omega=1.0_db/tau
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:)=0 !EAST
    isfluid(nx,:)=0 !WEST
    isfluid(:,1)=0 !SOUTH 
    isfluid(:,ny)=0 !NORTH
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    rho=1.0_db     !not to be intended as a delta rho
    !do ll=0,nlinks
    f0(1:nx,1:ny,1:2)=p(0)
    f1(1:nx,1:ny,1:2)=p(1)
    f2(1:nx,1:ny,1:2)=p(2)
    f3(1:nx,1:ny,1:2)=p(3)
    f4(1:nx,1:ny,1:2)=p(4)
    f5(1:nx,1:ny,1:2)=p(5)
    f6(1:nx,1:ny,1:2)=p(6)
    f7(1:nx,1:ny,1:2)=p(7)
    f8(1:nx,1:ny,1:2)=p(8)
    nsp=1
    nsk=2
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
    write(6,*) 'numero di gpu disponibili',ngpus
    write(6,*) '*******************************************'


    !$acc data copy(p,f0,f1,f2,f3,f4,f5,f6,f7,f8,isfluid,nsp,nsk)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************collision&forcing************************ 
        !$acc update device(nsp,nsk)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8) !async(2)
        !$acc loop private(uu,temp,udotc,feq,dummy) 
        do j=1,ny
           !$acc loop private(rho,u,v,uu,temp,udotc,feq,dummy)
            do i=1,nx
                if(isfluid(i,j).eq.1)then
                    rho =  f0(i,j,nsp) &
                          +f1(i,j,nsp) &
                          +f2(i,j,nsp) &
                          +f3(i,j,nsp) &
                          +f4(i,j,nsp) &
                          +f5(i,j,nsp) &
                          +f6(i,j,nsp) &
                          +f7(i,j,nsp) &
                          +f8(i,j,nsp)

                    u = (f1(i,j,nsp) +f5(i,j,nsp) +f8(i,j,nsp) &
                            -f3(i,j,nsp) -f6(i,j,nsp) -f7(i,j,nsp)) !/rho(i,j)
                        
                    v = (f5(i,j,nsp) +f2(i,j,nsp) +f6(i,j,nsp) &
                            -f7(i,j,nsp) -f4(i,j,nsp) -f8(i,j,nsp)) !/rho(i,j)
                    
                    uu=0.5_db*(u*u + v*v)/cssq
                    !oneminusuu= -uu !1.0_db - uu
                    
                    !0
                    feq=p(0)*(rho-uu)
                    f0(i,j,nsp)=f0(i,j,nsp) + omega*(feq - f0(i,j,nsp)) 
                    
                    !1
                    udotc=u/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p(1)*(rho+(temp + udotc))
                    f1(i,j,nsp)=f1(i,j,nsp) + omega*(feq - f1(i,j,nsp)) + fx*p(1)/cssq
                    
                    !3
                    feq=p(3)*(rho+(temp - udotc))
                    f3(i,j,nsp)=f3(i,j,nsp) + omega*(feq - f3(i,j,nsp)) - fx*p(3)/cssq
                    
                    !2
                    udotc=v/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p(2)*(rho+(temp + udotc))
                    f2(i,j,nsp)=f2(i,j,nsp) + omega*(feq - f2(i,j,nsp)) + fy*p(2)/cssq
                    
                    !4
                    feq=p(4)*(rho+(temp - udotc))
                    f4(i,j,nsp)=f4(i,j,nsp) + omega*(feq - f4(i,j,nsp)) - fy*p(4)/cssq
                    
                    !5
                    udotc=(u+v)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p(5)*(rho+(temp + udotc))
                    f5(i,j,nsp)=f5(i,j,nsp) + omega*(feq - f5(i,j,nsp)) + fx*p(5)/cssq + fy*p(5)/cssq 
                    
                    !7
                    feq=p(7)*(rho+(temp - udotc))
                    f7(i,j,nsp)=f7(i,j,nsp) + omega*(feq - f7(i,j,nsp)) - fx*p(7)/cssq - fy*p(7)/cssq
                    
                    !6
                    udotc=(-u+v)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p(6)*(rho+(temp + udotc))
                    f6(i,j,nsp)=f6(i,j,nsp) + omega*(feq - f6(i,j,nsp)) - fx*p(6)/cssq + fy*p(6)/cssq

                    !8
                    feq=p(8)*(rho+(temp - udotc))
                    f8(i,j,nsp)=f8(i,j,nsp) + omega*(feq - f8(i,j,nsp)) + fx*p(8)/cssq - fy*p(8)/cssq
                    !
                elseif(isfluid(i,j).eq.0)then
                    !1-3
                    dummy=f1(i,j,nsp)
					f1(i,j,nsp)=f3(i,j,nsp)
					f3(i,j,nsp)=dummy
                    !2-4
                    dummy=f2(i,j,nsp)
					f2(i,j,nsp)=f4(i,j,nsp)
					f4(i,j,nsp)=dummy
                    !5-7
                    dummy=f5(i,j,nsp)
					f5(i,j,nsp)=f7(i,j,nsp)
					f7(i,j,nsp)=dummy
                    !6-8
                    dummy=f6(i,j,nsp)
					f6(i,j,nsp)=f8(i,j,nsp)
					f8(i,j,nsp)=dummy
                endif
            enddo
        enddo
        !$acc end kernels
        !******************************************call bcs************************
        !periodic along y
        !$acc kernels 
        f0(2:nx-1,1,nsp)=f0(2:nx-1,ny-1,nsp)

        f1(2:nx-1,1,nsp)=f1(2:nx-1,ny-1,nsp)

        f2(2:nx-1,1,nsp)=f2(2:nx-1,ny-1,nsp)

        f3(2:nx-1,1,nsp)=f3(2:nx-1,ny-1,nsp)

        f4(2:nx-1,1,nsp)=f4(2:nx-1,ny-1,nsp)

        f5(2:nx-1,1,nsp)=f5(2:nx-1,ny-1,nsp)

        f6(2:nx-1,1,nsp)=f6(2:nx-1,ny-1,nsp)

        f7(2:nx-1,1,nsp)=f7(2:nx-1,ny-1,nsp)

        f8(2:nx-1,1,nsp)=f8(2:nx-1,ny-1,nsp)

        f0(2:nx-1,ny,nsp)=f0(2:nx-1,2,nsp)

        f1(2:nx-1,ny,nsp)=f1(2:nx-1,2,nsp)

        f2(2:nx-1,ny,nsp)=f2(2:nx-1,2,nsp)

        f3(2:nx-1,ny,nsp)=f3(2:nx-1,2,nsp)

        f4(2:nx-1,ny,nsp)=f4(2:nx-1,2,nsp)

        f5(2:nx-1,ny,nsp)=f5(2:nx-1,2,nsp)

        f6(2:nx-1,ny,nsp)=f6(2:nx-1,2,nsp)

        f7(2:nx-1,ny,nsp)=f7(2:nx-1,2,nsp)

        f8(2:nx-1,ny,nsp)=f8(2:nx-1,2,nsp)
        !$acc end kernels
        ! !******************************************streaming***************************  
        !
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8) !async(1)
        !$acc loop independent collapse(2)
        do j=1,ny
            !!$acc loop independent
            do i=1,nx
                    f1(i,j,nsk)=f1(i-1,j,nsp)
                    f3(i,j,nsk)=f3(i+1,j,nsp)
                    f2(i,j,nsk)=f2(i,j-1,nsp)
                    f4(i,j,nsk)=f4(i,j+1,nsp)
                    f5(i,j,nsk)=f5(i-1,j-1,nsp)
                    f6(i,j,nsk)=f6(i+1,j-1,nsp)
                    f7(i,j,nsk)=f7(i+1,j+1,nsp)
                    f8(i,j,nsk)=f8(i-1,j+1,nsp)
            enddo
        enddo
        !$acc end kernels
        dumm=nsp
        nsp=nsk
        nsk=dumm
        
    enddo 
    call cpu_time(ts2)
    !$acc update host(f0,f1,f2,f3,f4,f5,f6,f7,f8)
    !$acc end data
    rho = (f0(nx/2,1+(ny-1)/2,nsp) +f1(nx/2,1+(ny-1)/2,nsp) +f2(nx/2,1+(ny-1)/2,nsp) +f3(nx/2,1+(ny-1)/2,nsp) &
    +f4(nx/2,1+(ny-1)/2,nsp) +f5(nx/2,1+(ny-1)/2,nsp) +f6(nx/2,1+(ny-1)/2,nsp) +f7(nx/2,1+(ny-1)/2,nsp)+f8(nx/2,1+(ny-1)/2,nsp)  )
    u = (f1(nx/2,1+(ny-1)/2,nsp) +f5(nx/2,1+(ny-1)/2,nsp) +f8(nx/2,1+(ny-1)/2,nsp) &
                         -f3(nx/2,1+(ny-1)/2,nsp) -f6(nx/2,1+(ny-1)/2,nsp) -f7(nx/2,1+(ny-1)/2,nsp))
    v = (f5(nx/2,1+(ny-1)/2,nsp) +f2(nx/2,1+(ny-1)/2,nsp) +f6(nx/2,1+(ny-1)/2,nsp) &
                         -f7(nx/2,1+(ny-1)/2,nsp) -f4(nx/2,1+(ny-1)/2,nsp) -f8(nx/2,1+(ny-1)/2,nsp))
    write(6,*) 'u=',u ,'v=',v,'rho',rho !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'You''ve just wasted ', ts2-ts1, ' s of your life time' 

    
end program