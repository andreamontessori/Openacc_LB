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
    real(kind=db) :: fdum1,fdum2,fdum3,fdum4,fdum5,fdum6,fdum7,fdum8
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu,cssq,fx,temp,fdumm
    
    integer(kind=4), allocatable,  dimension(:)     :: ex,ey,opp
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    integer(kind=4), allocatable,  dimension(:,:,:) :: opposite
    
    real(kind=db), allocatable, dimension(:)     :: p,dex,dey,fdum
    real(kind=db), allocatable, dimension(:,:)   :: rho,u,v
    real(kind=db), allocatable, dimension(:,:,:) :: f,fbbk

       
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
    nx=65536
    ny=4096
    nsteps=1000
    stamp=1000
    fx=4.0_db*10.0**(-5)
    
    allocate(p(0:nlinks),ex(0:nlinks),ey(0:nlinks),dex(0:nlinks))
    allocate(dey(0:nlinks),opp(0:nlinks),fdum(0:nlinks))
    allocate(f(1:nx,1:ny,0:nlinks),fbbk(1:nx,1:ny,0:nlinks))
    allocate(u(1:nx,1:ny), v(1:nx,1:ny), rho(1:nx,1:ny))
    allocate(isfluid(1:nx,1:ny), opposite(1:nx,1:ny,1:nlinks)) !,omega_2d(1:nx,1:ny)) 
    
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
    f=0.0_db
    u=0.0_db
    v=0.0_db
    rho=0.0_db     !is to be meant as a delta rho
    do ll=0,nlinks
        f(1:nx,1:ny,ll)=rho(1:nx,1:ny)*p(ll)
    enddo
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

    !!$acc data copy(rho,u,v,ex,ey,dex,dey,f,p,opp) 
    !!!!!!!!!!!!!!!!!copyin(ex,ey,dex,dey,f,opposite,p)
    !$acc data copy(rho,u,v,ex,ey,dex,dey,f,p,opp,isfluid)
    !!$acc data copyin(rho,u,v,ex,ey,dex,dey,f,fbbk,opposite,p)
    !!$acc update device(rho,u,v,f)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !$acc kernels
        !**************************************** moments******************************
        do j=2,ny-1
            do i=2,nx-1         
                rho(i,j) = f(i,j,0) &
                          +f(i,j,1) &
                          +f(i,j,2) &
                          +f(i,j,3) &
                          +f(i,j,4) &
                          +f(i,j,5) &
                          +f(i,j,6) &
                          +f(i,j,7) &
                          +f(i,j,8)

                u(i,j) = (f(i,j,1) +f(i,j,5) +f(i,j,8) &
                         -f(i,j,3) -f(i,j,6) -f(i,j,7)) !/rho(i,j)
                       
                v(i,j) = (f(i,j,5) +f(i,j,2) +f(i,j,6) &
                         -f(i,j,7) -f(i,j,4) -f(i,j,8)) !/rho(i,j)
                
            enddo
        enddo  
        !
        ! if ( mod(step, stamp) .eq. 0 ) then
        !     ! !$acc update self(rho,u,v)
        !      write(6,*) 'timestep giusto ', step
        !      write(6,*) 'u=',u(nx/2,1+(ny-1)/2),'v=',v(nx/2,1+(ny-1)/2),'rho=',rho(nx/2,1+(ny-1)/2)
        !   endif 
        !

        !***********************************collision&forcing************************ 
        !$acc loop private(uu,temp,udotc,feq)
        do j=2,ny-1
            !$acc loop private(uu,temp,udotc,feq)
            do i=2,nx-1
                uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                !oneminusuu= -uu !1.0_db - uu
                
                !0
                feq=p(0)*(rho(i,j)-uu)
                f(i,j,0)=f(i,j,0) + omega*(feq - f(i,j,0)) 
                
                !1
                udotc=u(i,j)/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(1)*(rho(i,j)+(temp + udotc))
                f(i,j,1)=f(i,j,1) + omega*(feq - f(i,j,1)) + fx*p(1)/cssq
                
                !3
                feq=p(3)*(rho(i,j)+(temp - udotc))
                f(i,j,3)=f(i,j,3) + omega*(feq - f(i,j,3)) - fx*p(3)/cssq
                
                !2
                udotc=v(i,j)/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(2)*(rho(i,j)+(temp + udotc))
                f(i,j,2)=f(i,j,2) + omega*(feq - f(i,j,2))
                
                !4
                feq=p(4)*(rho(i,j)+(temp - udotc))
                f(i,j,4)=f(i,j,4) + omega*(feq - f(i,j,4))
                
                !5
                udotc=(u(i,j)+v(i,j))/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(5)*(rho(i,j)+(temp + udotc))
                f(i,j,5)=f(i,j,5) + omega*(feq - f(i,j,5)) + fx*p(5)/cssq
                
                !7
                feq=p(7)*(rho(i,j)+(temp - udotc))
                f(i,j,7)=f(i,j,7) + omega*(feq - f(i,j,7)) - fx*p(7)/cssq
                
                !6
                udotc=(-u(i,j)+v(i,j))/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(6)*(rho(i,j)+(temp + udotc))
                f(i,j,6)=f(i,j,6) + omega*(feq - f(i,j,6)) - fx*p(6)/cssq
                
                !8
                feq=p(8)*(rho(i,j)+(temp - udotc))
                f(i,j,8)=f(i,j,8) + omega*(feq - f(i,j,8)) + fx*p(8)/cssq
            enddo
        enddo
        !******************************************call bcs************************
        !bounceback! 
        !!$acc loop seq 
        do i=1,nx
            f(i,1,5)=f(i,1,7)
            f(i,1,6)=f(i,1,8)
            f(i,1,2)=f(i,1,4)
            f(i,ny,4)=f(i,ny,2)
            f(i,ny,7)=f(i,ny,5)
            f(i,ny,8)=f(i,ny,6)
        enddo  
        !!$acc loop seq 
        do j=1,ny
            f(1,j,1)=f(1,j,3)
            f(1,j,5)=f(1,j,7)
            f(1,j,8)=f(1,j,6)
            f(nx,j,2)=f(nx,j,1)
            f(nx,j,6)=f(nx,j,8)
            f(nx,j,7)=f(nx,j,5)
        enddo  
       
        ! Periodic
        !!$acc loop seq 
        do j=2,ny-1
            f(1,j,0)=f(nx-1,j,0)
            f(1,j,1)=f(nx-1,j,1)
            f(1,j,2)=f(nx-1,j,2)
            f(1,j,3)=f(nx-1,j,3)
            f(1,j,4)=f(nx-1,j,4)
            f(1,j,5)=f(nx-1,j,5)
            f(1,j,6)=f(nx-1,j,6)
            f(1,j,7)=f(nx-1,j,7)
            f(1,j,8)=f(nx-1,j,8)
            f(nx,j,0)=f(2,j,0)
            f(nx,j,1)=f(2,j,1)
            f(nx,j,2)=f(2,j,2)
            f(nx,j,3)=f(2,j,3)
            f(nx,j,4)=f(2,j,4)
            f(nx,j,5)=f(2,j,5)
            f(nx,j,6)=f(2,j,6)
            f(nx,j,7)=f(2,j,7)
            f(nx,j,8)=f(2,j,8)
        enddo
        
        ! f(1:nx,1,:)=f(1:nx,ny-1,:)
        ! f(1:nx,ny,:)=f(1:nx,2,:)
        ! !******************************************streaming***************************  
        !
         !$acc loop independent
        do j=2,ny-1        
             !$acc loop independent
            do i=nx,2,-1
                f(i,j,1)=f(i-1,j,1)
            enddo
             !$acc loop independent
            do i=1,nx-1
                f(i,j,3)=f(i+1,j,3)
            enddo
        enddo 
        ! 
         !$acc loop independent
        do j=ny,2,-1
             !$acc loop independent
            do i=2,nx-1
                f(i,j,2)=f(i,j-1,2) 
            enddo
             !$acc loop independent
            do i=1,nx-1
                f(i,j,6)=f(i+1,j-1,6)
            enddo
             !$acc loop independent
            do i=nx,2,-1
                f(i,j,5)=f(i-1,j-1,5)
            enddo
        enddo 
        !7
         !$acc loop independent
        do j=1,ny-1
            !$acc loop independent
            do i=2,nx-1
                f(i,j,4)=f(i,j+1,4)
            enddo
            !$acc loop independent
            do i=nx,2,-1
                f(i,j,8)=f(i-1,j+1,8)
            enddo
            !$acc loop independent
            do i=1,nx-1
                f(i,j,7)=f(i+1,j+1,7)
            enddo
        enddo  
        !$acc end kernels
    enddo 
    call cpu_time(ts2)
    !!$acc update host(rho,u,v)
    !$acc end data
    
    write(6,*) 'u=',u(nx/2,1+(ny-1)/2),'v=',v(nx/2,1+(ny-1)/2),'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'You''ve just wasted ', ts2-ts1, ' s of your life time' 

    
end program