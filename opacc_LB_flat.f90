program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    integer, parameter :: db=4 !kind(1.0)
    integer(kind=8) :: i,j,ll,l
    integer(kind=8) :: nx,ny,step,stamp,nlinks,nsteps,ngpus
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,dummy
    
    integer(kind=4), allocatable,  dimension(:)     :: ex,ey,opp
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    integer(kind=4), allocatable,  dimension(:,:,:) :: opposite
    
    real(kind=db), allocatable, dimension(:)     :: p,dex,dey,fdum
    real(kind=db), allocatable, dimension(:,:)   :: rho,u,v
    real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f1dummy

       
    nlinks=8 !pari!
    tau=0.8_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif

    !*******************************user parameters**************************
    nx=16384
    ny=16384
    nsteps=100
    stamp=1000
    fx=0.0_db*10.0**(-5)
    fy=5.0_db*10.0**(-5)
    allocate(p(0:nlinks),ex(0:nlinks),ey(0:nlinks),dex(0:nlinks))
    allocate(dey(0:nlinks),opp(0:nlinks),fdum(0:nlinks))
    allocate(f0(1:nx,1:ny),f1(1:nx,1:ny),f2(1:nx,1:ny),f3(1:nx,1:ny),f4(1:nx,1:ny))
    allocate(f5(1:nx,1:ny),f6(1:nx,1:ny),f7(1:nx,1:ny),f8(1:nx,1:ny))
    allocate(u(1:nx,1:ny), v(1:nx,1:ny), rho(1:nx,1:ny))
    allocate(isfluid(1:nx,1:ny), opposite(1:nx,1:ny,1:nlinks)) !,omega_2d(1:nx,1:ny)) 
    allocate(f1dummy(1:nx,1:ny))
    
    ex=(/0,1,0,-1,0,1,-1,-1,1/)
    ey=(/0,0,1,0,-1,1,1,-1,-1/)
    dex=real(ex,db)
    dey=real(ey,db)
    p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
    1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
    opp=(/0,3,4,1,2,7,8,5,6/)
    omega=1.0_db/tau
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:)=0 !EAST
    isfluid(nx,:)=0 !WEST
    isfluid(:,1)=0 !SOUTH 
    isfluid(:,ny)=0 !NORTH
    do ll=1,nlinks
        opposite(1:nx,1:ny,ll)=ll*isfluid(1:nx,1:ny) + opp(ll)*(1-isfluid(1:nx,1:ny))
    enddo
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    rho=0.0_db     !is to be intended as a delta rho
    !do ll=0,nlinks
    f0(1:nx,1:ny)=rho(1:nx,1:ny)*p(0)
    f1(1:nx,1:ny)=rho(1:nx,1:ny)*p(1)
    f2(1:nx,1:ny)=rho(1:nx,1:ny)*p(2)
    f3(1:nx,1:ny)=rho(1:nx,1:ny)*p(3)
    f4(1:nx,1:ny)=rho(1:nx,1:ny)*p(4)
    f5(1:nx,1:ny)=rho(1:nx,1:ny)*p(5)
    f6(1:nx,1:ny)=rho(1:nx,1:ny)*p(6)
    f7(1:nx,1:ny)=rho(1:nx,1:ny)*p(7)
    f8(1:nx,1:ny)=rho(1:nx,1:ny)*p(8)
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
    pause


    !$acc data copy(rho,u,v,ex,ey,dex,dey,p,opp, &
    !$acc & f0,f1,f2,f3,f4,f5,f6,f7,f8,isfluid)
    
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !!$acc update device(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !$acc kernels
        !**************************************** moments******************************
        do j=2,ny-1
            do i=2,nx-1         
                rho(i,j) = f0(i,j) &
                          +f1(i,j) &
                          +f2(i,j) &
                          +f3(i,j) &
                          +f4(i,j) &
                          +f5(i,j) &
                          +f6(i,j) &
                          +f7(i,j) &
                          +f8(i,j)

                u(i,j) = (f1(i,j) +f5(i,j) +f8(i,j) &
                         -f3(i,j) -f6(i,j) -f7(i,j)) !/rho(i,j)
                       
                v(i,j) = (f5(i,j) +f2(i,j) +f6(i,j) &
                         -f7(i,j) -f4(i,j) -f8(i,j)) !/rho(i,j)
                
            enddo
        enddo  
        !!$acc end kernels
        !
        ! if ( mod(step, stamp) .eq. 0 ) then
        !     ! !$acc update self(rho,u,v)
        !      write(6,*) 'timestep giusto ', step
        !      write(6,*) 'u=',u(nx/2,1+(ny-1)/2),'v=',v(nx/2,1+(ny-1)/2),'rho=',rho(nx/2,1+(ny-1)/2)
        !   endif 
        !
        !***********************************collision&forcing************************ 
        !!$acc update host(rho,u,v)
        !!$acc update device(rho,u,v)
        !!$acc kernels
        !$acc loop private(uu,temp,udotc,feq,dummy)
        do j=2,ny-1
           !$acc loop private(uu,temp,udotc,feq,dummy)
            do i=2,nx-1
                if(isfluid(i,j).eq.1)then
                    uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                    !oneminusuu= -uu !1.0_db - uu
                    
                    !0
                    feq=p(0)*(rho(i,j)-uu)
                    f0(i,j)=f0(i,j) + omega*(feq - f0(i,j)) 
                    
                    !1
                    udotc=u(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p(1)*(rho(i,j)+(temp + udotc))
                    f1(i,j)=f1(i,j) + omega*(feq - f1(i,j)) + fx*p(1)/cssq
                    
                    !3
                    feq=p(3)*(rho(i,j)+(temp - udotc))
                    f3(i,j)=f3(i,j) + omega*(feq - f3(i,j)) - fx*p(3)/cssq
                    
                    !2
                    udotc=v(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p(2)*(rho(i,j)+(temp + udotc))
                    f2(i,j)=f2(i,j) + omega*(feq - f2(i,j)) + fy*p(2)/cssq
                    
                    !4
                    feq=p(4)*(rho(i,j)+(temp - udotc))
                    f4(i,j)=f4(i,j) + omega*(feq - f4(i,j)) - fy*p(4)/cssq
                    
                    !5
                    udotc=(u(i,j)+v(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p(5)*(rho(i,j)+(temp + udotc))
                    f5(i,j)=f5(i,j) + omega*(feq - f5(i,j)) + fx*p(5)/cssq + fy*p(5)/cssq 
                    
                    !7
                    feq=p(7)*(rho(i,j)+(temp - udotc))
                    f7(i,j)=f7(i,j) + omega*(feq - f7(i,j)) - fx*p(7)/cssq - fy*p(7)/cssq
                    
                    !6
                    udotc=(-u(i,j)+v(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p(6)*(rho(i,j)+(temp + udotc))
                    f6(i,j)=f6(i,j) + omega*(feq - f6(i,j)) - fx*p(6)/cssq + fy*p(6)/cssq
                    
                    !8
                    feq=p(8)*(rho(i,j)+(temp - udotc))
                    f8(i,j)=f8(i,j) + omega*(feq - f8(i,j)) + fx*p(8)/cssq - fy*p(8)/cssq
                elseif(isfluid(i,j).eq.0)then
                    !1-3
                    dummy=f1(i,j)
					f1(i,j)=f3(i,j)
					f3(i,j)=dummy
                    !2-4
                    dummy=f2(i,j)
					f2(i,j)=f4(i,j)
					f4(i,j)=dummy
                    !5-7
                    dummy=f5(i,j)
					f5(i,j)=f7(i,j)
					f7(i,j)=dummy
                    !6-8
                    dummy=f6(i,j)
					f6(i,j)=f8(i,j)
					f8(i,j)=dummy
                endif
            enddo
        enddo
        !!$acc end kernels
        !******************************************call bcs************************
        !!$acc update host(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !periodic along y
        !!$acc kernels
        f0(2:nx-1,1)=f0(2:nx-1,ny-1)
        f1(2:nx-1,1)=f1(2:nx-1,ny-1)
        f2(2:nx-1,1)=f2(2:nx-1,ny-1)
        f3(2:nx-1,1)=f3(2:nx-1,ny-1)
        f4(2:nx-1,1)=f4(2:nx-1,ny-1)
        f5(2:nx-1,1)=f5(2:nx-1,ny-1)
        f6(2:nx-1,1)=f6(2:nx-1,ny-1)
        f7(2:nx-1,1)=f7(2:nx-1,ny-1)
        f8(2:nx-1,1)=f8(2:nx-1,ny-1)
        !
        f0(2:nx-1,ny)=f0(2:nx-1,2)
        f1(2:nx-1,ny)=f1(2:nx-1,2)
        f2(2:nx-1,ny)=f2(2:nx-1,2)
        f3(2:nx-1,ny)=f3(2:nx-1,2)
        f4(2:nx-1,ny)=f4(2:nx-1,2)
        f5(2:nx-1,ny)=f5(2:nx-1,2)
        f6(2:nx-1,ny)=f6(2:nx-1,2)
        f7(2:nx-1,ny)=f7(2:nx-1,2)
        f8(2:nx-1,ny)=f8(2:nx-1,2)
        !!$acc end kernels
        ! f(1:nx,1,:)=f(1:nx,ny-1,:)
        ! f(1:nx,ny,:)=f(1:nx,2,:)
        ! !******************************************streaming***************************  
        !
        !!$acc update host(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !!$acc update device(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !!$acc kernels
        ! right
        f1(2:nx,2:ny-1)=f1(1:nx,2:ny-1)
        f3(1:nx-1,2:ny-1)=f3(2:nx,2:ny-1)
        !up 
        f2(2:nx-1,2:ny)=f2(2:nx-1,1:ny-1)
        f4(2:nx-1,1:ny-1)=f4(2:nx-1,2:ny)
        ! up right
        f5(2:nx,2:ny)=f5(1:nx-1,1:ny-1)
        f7(1:nx-1,1:ny-1)=f7(2:nx,2:ny)
        ! up left
        f6(1:nx-1,2:ny)=f6(2:nx,1:ny-1)
        f8(2:nx,1:ny-1)=f8(1:nx,2:ny)
        !$acc end kernels
        
    enddo 
    call cpu_time(ts2)
    !$acc update host(rho,u,v)
    !$acc end data
    
    write(6,*) 'u=',u(nx/2,1+(ny-1)/2),'v=',v(nx/2,1+(ny-1)/2),'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'You''ve just wasted ', ts2-ts1, ' s of your life time' 

    
end program