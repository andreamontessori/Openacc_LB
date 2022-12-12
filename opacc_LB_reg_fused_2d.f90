program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    integer, parameter :: db=4 !kind(1.0)
    integer(kind=8) :: i,j,ll,l,dumm
    integer(kind=8) :: nx,ny,step,stamp,nlinks,nsteps,ngpus
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2 
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8
    real(kind=db) :: qxx1,qxx3,qxx5,qxx6,qxx7,qxx8
    real(kind=db) :: qyy2,qyy4,qyy5,qyy6,qyy7,qyy8
    real(kind=db) :: qxy5,qxy6,qxy7,qxy8,pi2cssq1,pi2cssq2
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,dummy
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p
    real(kind=db), allocatable, dimension(:,:) :: rho,u,v,pxx,pyy,pxy
    real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8

       
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
    nx=101
    ny=4024
    nsteps=40000
    stamp=1000
    fx=0.0_db*10.0**(-7)
    fy=1.0_db*10.0**(-6)
    allocate(p(0:nlinks))
    allocate(f0(0:nx+1,0:ny+1),f1(0:nx+1,0:ny+1),f2(0:nx+1,0:ny+1),f3(0:nx+1,0:ny+1),f4(0:nx+1,0:ny+1))
    allocate(f5(0:nx+1,0:ny+1),f6(0:nx+1,0:ny+1),f7(0:nx+1,0:ny+1),f8(0:nx+1,0:ny+1))
    allocate(rho(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny),pxx(1:nx,1:ny),pyy(1:nx,1:ny),pxy(1:nx,1:ny))
    allocate(isfluid(1:nx,1:ny)) !,omega_2d(1:nx,1:ny)) 
    
    !ex=(/0,1,0,-1,0,1,-1,-1,1/)
    !ey=(/0,0,1,0,-1,1,1,-1,-1/)

    p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
    1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
    omega=1.0_db/tau

    ! regularized: hermite 
    qxx1=1.0_db-cssq
    qxx3=1.0_db-cssq
    qxx5=1.0_db
    qxx6=-1.0_db
    qxx7=1.0_db
    qxx8=-1.0_db

    qyy2=1.0_db-cssq
    qyy4=1.0_db-cssq
    qyy5=1.0_db
    qyy6=-1.0_db
    qyy7=1.0_db
    qyy8=-1.0_db

    qxy5=1.0_db
    qxy6=-1.0_db
    qxy7=1.0_db
    qxy8=-1.0_db

    pi2cssq1=p(1)/(2.0_db*cssq**2)
    pi2cssq2=p(5)/(2.0_db*cssq**2)
    
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:)=0 !EAST
    isfluid(nx,:)=0 !WEST
    isfluid(:,1)=0 !SOUTH 
    isfluid(:,ny)=0 !NORTH
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    rho=0.0_db     !is to be intended as a delta rho
    !do ll=0,nlinks
    f0(1:nx,1:ny)=0.0_db
    f1(1:nx,1:ny)=0.0_db
    f2(1:nx,1:ny)=0.0_db
    f3(1:nx,1:ny)=0.0_db
    f4(1:nx,1:ny)=0.0_db
    f5(1:nx,1:ny)=0.0_db
    f6(1:nx,1:ny)=0.0_db
    f7(1:nx,1:ny)=0.0_db
    f8(1:nx,1:ny)=0.0_db
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


    !$acc data copy(p,rho,u,v,pxx,pxy,pyy,f0,f1,f2,f3,f4,f5,f6,f7,f8,isfluid)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moment + neq pressor*********
        !$acc kernels !present(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !$acc loop collapse(2) private(uu,temp,udotc,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8) 
        do j=1,ny
            do i=1,nx
                if(isfluid(i,j).eq.1.or.isfluid(i,j).eq.0)then
                    rho(i,j) = f0(i,j)+f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)

                    u(i,j) = (f1(i,j) +f5(i,j) +f8(i,j)-f3(i,j) -f6(i,j) -f7(i,j)) !/rho(i,j)
                        
                    v(i,j) = (f5(i,j) +f2(i,j) +f6(i,j)-f7(i,j) -f4(i,j) -f8(i,j))
                    ! non equilibrium pressor components
                    uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                    !1-3
                    udotc=u(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq1=f1(i,j)-p(1)*(rho(i,j)+(temp + udotc))
                    fneq3=f3(i,j)-p(3)*(rho(i,j)+(temp - udotc))
                    !2-4
                    udotc=v(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq2=f2(i,j)-p(2)*(rho(i,j)+(temp + udotc))
                    fneq4=f4(i,j)-p(4)*(rho(i,j)+(temp - udotc))
                    !5-7
                    udotc=(u(i,j)+v(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq5=f5(i,j)-p(5)*(rho(i,j)+(temp + udotc))
                    fneq7=f7(i,j)-p(7)*(rho(i,j)+(temp - udotc))
                    !6-8
                    udotc=(-u(i,j)+v(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq6=f6(i,j)-p(6)*(rho(i,j)+(temp + udotc))
                    fneq8=f8(i,j)-p(8)*(rho(i,j)+(temp - udotc))

                    pxx(i,j)= fneq1 + fneq3 + fneq5 + fneq6 + fneq7 + fneq8
                    pyy(i,j)= fneq2 + fneq4 + fneq5 + fneq6 + fneq7 + fneq8
                    pxy(i,j)= fneq5 - fneq6 + fneq7 - fneq8
                endif
                if(isfluid(i,j).eq.0)then
                     f1(i,j)=p(3)*rho(i,j) + pi2cssq1*qxx3*pxx(i,j)
                     f3(i,j)=p(1)*rho(i,j) + pi2cssq1*qxx1*pxx(i,j)
                     f2(i,j)=p(4)*rho(i,j) + pi2cssq1*qyy4*pyy(i,j)
                     f4(i,j)=p(2)*rho(i,j) + pi2cssq1*qyy2*pyy(i,j)
                     f5(i,j)=p(7)*rho(i,j) + pi2cssq2*(qxx7*pxx(i,j)+qyy7*pyy(i,j)+2.0_db*qxy7*pxy(i,j))
                     f7(i,j)=p(5)*rho(i,j) + pi2cssq2*(qxx5*pxx(i,j)+qyy5*pyy(i,j)+2.0_db*qxy5*pxy(i,j))
                     f6(i,j)=p(8)*rho(i,j) + pi2cssq2*(qxx8*pxx(i,j)+qyy8*pyy(i,j)+2.0_db*qxy8*pxy(i,j))
                     f8(i,j)=p(6)*rho(i,j) + pi2cssq2*(qxx6*pxx(i,j)+qyy6*pyy(i,j)+2.0_db*qxy6*pxy(i,j)) 
                endif
            enddo
        enddo
        !***********************************collision + no slip + forcing: fused implementation*********
        !$acc loop collapse(2) private(uu,temp,udotc,feq,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8) 
        do j=1,ny
           do i=1,nx    
                uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                !oneminusuu= -uu !1.0_db - uu
                !0
                feq=p(0)*(rho(i,j)-uu)
                f0(i,j)=f0(i,j) + omega*(feq - f0(i,j)) 
                !1
                udotc=u(i,j)/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(1)*(rho(i,j)+(temp + udotc))
                f1(i+1,j)= feq + (1.0_db-omega)*pi2cssq1*qxx1*pxx(i,j) + fx*p(1)/cssq !f1(i-1,j,nsp) + omega*(feq - f1(i-1,j,nsp)) + fx*p(1)/cssq
                !3
                feq=p(3)*(rho(i,j)+(temp - udotc))
                f3(i-1,j)=feq + (1.0_db-omega)*pi2cssq1*qxx3*pxx(i,j) - fx*p(3)/cssq !f3(i+1,j,nsp) + omega*(feq - f3(i+1,j,nsp)) - fx*p(3)/cssq
                !2
                udotc=v(i,j)/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(2)*(rho(i,j)+(temp + udotc))
                f2(i,j+1)= feq + (1.0_db-omega)*pi2cssq1*qyy2*pyy(i,j) + fy*p(2)/cssq !f2(i,j-1,nsp) + omega*(feq - f2(i,j-1,nsp)) + fy*p(2)/cssq
                !4
                feq=p(4)*(rho(i,j)+(temp - udotc))
                f4(i,j-1)=feq + (1.0_db-omega)*pi2cssq1*qyy4*pyy(i,j) - fy*p(4)/cssq !f4(i,j+1,nsp) + omega*(feq - f4(i,j+1,nsp)) - fy*p(4)/cssq
                !5
                udotc=(u(i,j)+v(i,j))/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(5)*(rho(i,j)+(temp + udotc))
                f5(i+1,j+1)= feq + (1.0_db-omega)*pi2cssq2*(qxx5*pxx(i,j)+qyy5*pyy(i,j)+2.0_db*qxy5*pxy(i,j)) + fx*p(5)/cssq + fy*p(5)/cssq!f5(i-1,j-1,nsp) + omega*(feq - f5(i-1,j-1,nsp)) + fx*p(5)/cssq + fy*p(5)/cssq 
                !7
                feq=p(7)*(rho(i,j)+(temp - udotc))
                f7(i-1,j-1)=feq + (1.0_db-omega)*pi2cssq2*(qxx7*pxx(i,j)+qyy7*pyy(i,j)+2.0_db*qxy7*pxy(i,j)) - fx*p(7)/cssq - fy*p(7)/cssq !f7(i+1,j+1,nsp) + omega*(feq - f7(i+1,j+1,nsp)) - fx*p(7)/cssq - fy*p(7)/cssq
                !6
                udotc=(-u(i,j)+v(i,j))/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(6)*(rho(i,j)+(temp + udotc))
                f6(i-1,j+1)= feq + (1.0_db-omega)*pi2cssq2*(qxx6*pxx(i,j)+qyy6*pyy(i,j)+2.0_db*qxy6*pxy(i,j)) - fx*p(6)/cssq + fy*p(6)/cssq !f6(i+1,j-1,nsp) + omega*(feq - f6(i+1,j-1,nsp)) - fx*p(6)/cssq + fy*p(6)/cssq
                !8
                feq=p(8)*(rho(i,j)+(temp - udotc))
                f8(i+1,j-1)=feq + (1.0_db-omega)*pi2cssq2*(qxx8*pxx(i,j)+qyy8*pyy(i,j)+2.0_db*qxy8*pxy(i,j)) + fx*p(8)/cssq - fy*p(8)/cssq !f8(i-1,j+1,nsp) + omega*(feq - f8(i-1,j+1,nsp)) + fx*p(8)/cssq - fy*p(8)/cssq
            enddo
        enddo
        !!$acc end kernels
        !******************************************call periodic bcs************************
        !periodic along y
        !!$acc kernels 
        f5(2:nx-1,2)=f5(2:nx-1,ny)
        f2(2:nx-1,2)=f2(2:nx-1,ny)
        f6(2:nx-1,2)=f6(2:nx-1,ny)
        f8(2:nx-1,ny-1)=f8(2:nx-1,1)
        f4(2:nx-1,ny-1)=f4(2:nx-1,1)
        f7(2:nx-1,ny-1)=f7(2:nx-1,1)
        !$acc end kernels 

    enddo 
    call cpu_time(ts2)
    !$acc update host(rho,u,v)

    !$acc end data


    !************************************************test points**********************************************!
    
    write(6,*) 'u=',u(nx/2,ny/2) ,'v=',v(nx/2,ny/2),'rho',rho(nx/2,ny/2) !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'u=',u(2,ny/2) ,'v=',v(2,ny/2),'rho',rho(2,ny/2)
    write(6,*) 'u=',u(1,ny/2) ,'v=',v(1,ny/2),'rho',rho(1,ny/2)
    write(6,*) 'You''ve just wasted ', ts2-ts1, ' s of your life time' 

    
end program