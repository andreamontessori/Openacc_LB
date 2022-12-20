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
    real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq1,pi2cssq2
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,dummy,fpc
    real(kind=db) :: addendum0,addendum1,addendum2,addendum3,addendum4,addendum5,addendum6,addendum7,addendum8
    real(kind=db) :: gaddendum0,gaddendum1,gaddendum2,gaddendum3,gaddendum4,gaddendum5,gaddendum6,gaddendum7,gaddendum8
    real(kind=db) :: psi_x,psi_y,mod_psi,mod_psi_sq,st_coeff,b0,b1,b2,beta,sigma
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p
    real(kind=db), allocatable, dimension(:,:) :: psi,rho,u,v,pxx,pyy,pxy
    real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8
    real(kind=db), allocatable, dimension(:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8

    !lattice pars  
    nlinks=8 !pari!
    tau=1.0_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
    omega=1.0_db/tau
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif

    !*******************************user parameters**************************
    nx=4096
    ny=4096
    nsteps=1000
    stamp=10000
    fx=0.0_db*10.0**(-7)
    fy=0.0_db*10.0**(-5)
    !**********************************allocation****************************
    allocate(p(0:nlinks))
    allocate(f0(0:nx+1,0:ny+1),f1(0:nx+1,0:ny+1),f2(0:nx+1,0:ny+1),f3(0:nx+1,0:ny+1),f4(0:nx+1,0:ny+1))
    allocate(f5(0:nx+1,0:ny+1),f6(0:nx+1,0:ny+1),f7(0:nx+1,0:ny+1),f8(0:nx+1,0:ny+1))
    allocate(g0(0:nx+1,0:ny+1),g1(0:nx+1,0:ny+1),g2(0:nx+1,0:ny+1),g3(0:nx+1,0:ny+1),g4(0:nx+1,0:ny+1))
    allocate(g5(0:nx+1,0:ny+1),g6(0:nx+1,0:ny+1),g7(0:nx+1,0:ny+1),g8(0:nx+1,0:ny+1))
    allocate(psi(0:nx+1,0:ny+1),rho(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny),pxx(1:nx,1:ny),pyy(1:nx,1:ny),pxy(1:nx,1:ny))
    allocate(isfluid(1:nx,1:ny)) !,omega_2d(1:nx,1:ny)) 
    p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
    1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
    !ex=(/0,1,0,-1,0,1,-1,-1,1/)
    !ey=(/0,0,1,0,-1,1,1,-1,-1/)

    ! regularized: hermite 
    qxx=1.0_db-cssq
    qyy=1.0_db-cssq
    qxy5_7=1.0_db
    qxy6_8=-1.0_db
    pi2cssq1=p(1)/(2.0_db*cssq**2)
    pi2cssq2=p(5)/(2.0_db*cssq**2)
    ! chromodynamic
    beta=1.5_db
    sigma=0.02_db
    st_coeff=(9.0_db/4.0_db)*sigma*omega
    b0=-4.0_db/27.0_db
    b1=2.0_db/27.0_db
    b2=5.0_db/108.0_db
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:)=0 !EAST
    isfluid(nx,:)=0 !WEST
    isfluid(:,1)=0 !SOUTH 
    isfluid(:,ny)=0 !NORTH
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    rho=1.0_db     !is to be intended as a delta rho
    psi=-1.0_db
    !do ll=0,nlinks
    f0(1:nx,1:ny)=p(0)*rho(:,:)
    f1(1:nx,1:ny)=p(1)*rho(:,:)
    f2(1:nx,1:ny)=p(2)*rho(:,:)
    f3(1:nx,1:ny)=p(3)*rho(:,:)
    f4(1:nx,1:ny)=p(4)*rho(:,:)
    f5(1:nx,1:ny)=p(5)*rho(:,:)
    f6(1:nx,1:ny)=p(6)*rho(:,:)
    f7(1:nx,1:ny)=p(7)*rho(:,:)
    f8(1:nx,1:ny)=p(8)*rho(:,:)
    !
    do i=nx/2-10,nx/2+10
        do j=ny/2-10,ny/2+10
            if ((i-nx/2)**2+(j-ny/2)**2<=10**2)then
                psi(i,j)=1.0_db
            endif
        enddo
    enddo
    g0(1:nx,1:ny)=p(0)*psi(1:nx,1:ny)
    g1(1:nx,1:ny)=p(1)*psi(1:nx,1:ny)
    g2(1:nx,1:ny)=p(2)*psi(1:nx,1:ny)
    g3(1:nx,1:ny)=p(3)*psi(1:nx,1:ny)
    g4(1:nx,1:ny)=p(4)*psi(1:nx,1:ny)
    g5(1:nx,1:ny)=p(5)*psi(1:nx,1:ny)
    g6(1:nx,1:ny)=p(6)*psi(1:nx,1:ny)
    g7(1:nx,1:ny)=p(7)*psi(1:nx,1:ny)
    g8(1:nx,1:ny)=p(8)*psi(1:nx,1:ny)
    !*************************************check data ************************ 
    write(6,*) '*******************LB data*****************'
    write(6,*) 'tau',tau
    write(6,*) 'omega',omega
    write(6,*) 'visc',visc_LB
    write(6,*) 'fx',fx
    write(6,*) 'cssq',cssq
    write(6,*) 'beta',beta
    write(6,*) 'sigma',sigma
    write(6,*) 'surface_tens_coeff',st_coeff
    write(6,*) '*******************INPUT data*****************'
    write(6,*) 'nx',nx
    write(6,*) 'ny',ny
    write(6,*) 'nsteps',nsteps
    write(6,*) 'stamp',stamp
    write(6,*) 'max fx',huge(fx)
    write(6,*) 'available gpus',ngpus
    write(6,*) '*******************************************'


    !$acc data copy(p,rho,u,v,pxx,pxy,pyy,f0,f1,f2,f3,f4,f5,f6,f7,f8,isfluid, &
    !$acc& g0,g1,g2,g3,g4,g5,g6,g7,g8,psi)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moment + neq pressor*********
        !$acc kernels 
        !$acc loop collapse(2) private(uu,temp,udotc,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8) 
        do j=1,ny
            do i=1,nx
                if(isfluid(i,j).eq.1.or.isfluid(i,j).eq.0)then
                    rho(i,j) = f0(i,j)+f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)
                    u(i,j) = (f1(i,j) +f5(i,j) +f8(i,j)-f3(i,j) -f6(i,j) -f7(i,j)) !/rho(i,j)
                    v(i,j) = (f5(i,j) +f2(i,j) +f6(i,j)-f7(i,j) -f4(i,j) -f8(i,j))
                    psi(i,j)= (g0(i,j)+g1(i,j)+g2(i,j)+g3(i,j)+g4(i,j)+g5(i,j)+g6(i,j)+g7(i,j)+g8(i,j))/rho(i,j)
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
                !no slip everywhere, always before fused: to be modified for generic pressure/velocity bcs
                if(isfluid(i,j).eq.0)then
                     f1(i,j)=p(3)*rho(i,j) + pi2cssq1*qxx*pxx(i,j)
                     f3(i,j)=p(1)*rho(i,j) + pi2cssq1*qxx*pxx(i,j)
                     f2(i,j)=p(4)*rho(i,j) + pi2cssq1*qyy*pyy(i,j)
                     f4(i,j)=p(2)*rho(i,j) + pi2cssq1*qyy*pyy(i,j)
                     f5(i,j)=p(7)*rho(i,j) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j))
                     f7(i,j)=p(5)*rho(i,j) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j))
                     f6(i,j)=p(8)*rho(i,j) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j))
                     f8(i,j)=p(6)*rho(i,j) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j))
                     psi(i,j)=-1.0_db

                endif
            enddo
        enddo
        !***********************************collision + no slip + forcing: fused implementation*********
        !$acc loop collapse(2)  
        do j=1,ny
           do i=1,nx   
                uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                !oneminusuu= -uu !1.0_db - uu
                !0
                !chromodynamic
                psi_x=(1.0_db/cssq)*(p(1)*(psi(i+1,j)-psi(i-1,j)) + p(5)*(psi(i+1,j+1) + psi(i+1,j-1)-psi(i-1,j+1)-psi(i-1,j-1)))
                psi_y=(1.0_db/cssq)*(p(1)*(psi(i,j+1)-psi(i,j-1)) + p(5)*(psi(i+1,j+1) - psi(i+1,j-1)+psi(i-1,j+1)-psi(i-1,j-1)))
                mod_psi=sqrt(psi_x**2+psi_y**2)
                mod_psi_sq=psi_x**2+psi_y**2 
                addendum0=0.0_db
                addendum1=0.0_db
                addendum2=0.0_db
                addendum3=0.0_db
                addendum4=0.0_db
                addendum5=0.0_db
                addendum6=0.0_db
                addendum7=0.0_db
                addendum8=0.0_db
                gaddendum0=0.0_db
                gaddendum1=0.0_db
                gaddendum2=0.0_db
                gaddendum3=0.0_db
                gaddendum4=0.0_db
                gaddendum5=0.0_db
                gaddendum6=0.0_db
                gaddendum7=0.0_db
                gaddendum8=0.0_db
                if(mod_psi>0.01)then
                addendum0=-st_coeff*mod_psi*b0
                addendum1=st_coeff*mod_psi*(p(1)*psi_x**2/mod_psi_sq - b1)
                addendum2=st_coeff*mod_psi*(p(2)*psi_y**2/mod_psi_sq - b1)
                addendum3=st_coeff*mod_psi*(p(3)*psi_x**2/mod_psi_sq - b1)
                addendum4=st_coeff*mod_psi*(p(4)*psi_y**2/mod_psi_sq - b1)
                addendum5=st_coeff*mod_psi*(p(5)*(psi_x+psi_y)**2/mod_psi_sq - b2)
                addendum6=st_coeff*mod_psi*(p(6)*(-psi_x+psi_y)**2/mod_psi_sq - b2)
                addendum7=st_coeff*mod_psi*(p(7)*(-psi_x-psi_y)**2/mod_psi_sq - b2)
                addendum8=st_coeff*mod_psi*(p(8)*(psi_x-psi_y)**2/mod_psi_sq - b2)
                !recolor
                gaddendum1=0.25*beta*rho(i,j)*p(1)*(1.0_db-psi(i,j)**2)*psi_x/mod_psi
                gaddendum2=0.25*beta*rho(i,j)*p(2)*(1.0_db-psi(i,j)**2)*psi_y/mod_psi
                gaddendum3=0.25*beta*rho(i,j)*p(3)*(1.0_db-psi(i,j)**2)*(-psi_x/mod_psi)
                gaddendum4=0.25*beta*rho(i,j)*p(4)*(1.0_db-psi(i,j)**2)*(-psi_y/mod_psi)
                gaddendum5=0.25*beta*rho(i,j)*p(5)*(1.0_db-psi(i,j)**2)*(psi_x/mod_psi + psi_y/mod_psi)
                gaddendum6=0.25*beta*rho(i,j)*p(6)*(1.0_db-psi(i,j)**2)*(-psi_x/mod_psi + psi_y/mod_psi)
                gaddendum7=0.25*beta*rho(i,j)*p(7)*(1.0_db-psi(i,j)**2)*(-psi_x/mod_psi - psi_y/mod_psi)
                gaddendum8=0.25*beta*rho(i,j)*p(8)*(1.0_db-psi(i,j)**2)*(psi_x/mod_psi - psi_y/mod_psi)
                endif
                !regularized collision + perturbation
                feq=p(0)*(rho(i,j)-uu)
                f0(i,j)=f0(i,j) + omega*(feq - f0(i,j)) + addendum0
                g0(i,j)=psi(i,j)*f0(i,j)
                !1
                udotc=u(i,j)/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(1)*(rho(i,j)+(temp + udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq1*qxx*pxx(i,j) + fx*p(1)/cssq + addendum1
                g1(i+1,j)=psi(i,j)*fpc + gaddendum1
                f1(i+1,j)= fpc !f1(i-1,j,nsp) + omega*(feq - f1(i-1,j,nsp)) + fx*p(1)/cssq
               
                !3
                feq=p(3)*(rho(i,j)+(temp - udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq1*qxx*pxx(i,j) - fx*p(3)/cssq + addendum3
                g3(i-1,j)=psi(i,j)*fpc + gaddendum3
                f3(i-1,j)= fpc !f3(i+1,j,nsp) + omega*(feq - f3(i+1,j,nsp)) - fx*p(3)/cssq
                !2
                udotc=v(i,j)/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(2)*(rho(i,j)+(temp + udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq1*qyy*pyy(i,j) + fy*p(2)/cssq + addendum2
                g2(i,j+1)=psi(i,j)*fpc + gaddendum2  
                f2(i,j+1)= fpc !f2(i,j-1,nsp) + omega*(feq - f2(i,j-1,nsp)) + fy*p(2)/cssq
                !4
                feq=p(4)*(rho(i,j)+(temp - udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq1*qyy*pyy(i,j) - fy*p(4)/cssq + addendum4 
                g4(i,j-1)=psi(i,j)*fpc +gaddendum4
                f4(i,j-1)=fpc!f4(i,j+1,nsp) + omega*(feq - f4(i,j+1,nsp)) - fy*p(4)/cssq
                !5
                udotc=(u(i,j)+v(i,j))/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(5)*(rho(i,j)+(temp + udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + fx*p(5)/cssq + fy*p(5)/cssq + addendum5
                g5(i+1,j+1)=psi(i,j)*fpc + gaddendum5
                f5(i+1,j+1)=fpc !f5(i-1,j-1,nsp) + omega*(feq - f5(i-1,j-1,nsp)) + fx*p(5)/cssq + fy*p(5)/cssq 
                !7
                feq=p(7)*(rho(i,j)+(temp - udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) - fx*p(7)/cssq - fy*p(7)/cssq + addendum7
                g7(i-1,j-1)=psi(i,j)*fpc + gaddendum7
                f7(i-1,j-1)=fpc !f7(i+1,j+1,nsp) + omega*(feq - f7(i+1,j+1,nsp)) - fx*p(7)/cssq - fy*p(7)/cssq
                !6
                udotc=(-u(i,j)+v(i,j))/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(6)*(rho(i,j)+(temp + udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) - fx*p(6)/cssq + fy*p(6)/cssq + addendum6
                g6(i-1,j+1)=psi(i,j)*fpc + gaddendum6
                f6(i-1,j+1)= fpc!f6(i+1,j-1,nsp) + omega*(feq - f6(i+1,j-1,nsp)) - fx*p(6)/cssq + fy*p(6)/cssq
                !8
                feq=p(8)*(rho(i,j)+(temp - udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) + fx*p(8)/cssq - fy*p(8)/cssq + addendum8
                g8(i+1,j-1)=psi(i,j)*fpc + gaddendum8
                f8(i+1,j-1)=fpc !f8(i-1,j+1,nsp) + omega*(feq - f8(i-1,j+1,nsp)) + fx*p(8)/cssq - fy*p(8)/cssq
            enddo
        enddo
        !
        !******************************************call periodic bcs: always after fused************************
        !periodic along y
        !!$acc kernels 
        ! f5(2:nx-1,2)=f5(2:nx-1,ny)
        ! f2(2:nx-1,2)=f2(2:nx-1,ny)
        ! f6(2:nx-1,2)=f6(2:nx-1,ny)
        ! f8(2:nx-1,ny-1)=f8(2:nx-1,1)
        ! f4(2:nx-1,ny-1)=f4(2:nx-1,1)
        ! f7(2:nx-1,ny-1)=f7(2:nx-1,1)
        !$acc end kernels 

    enddo 
    call cpu_time(ts2)
    !$acc update host(rho,u,v)

    !$acc end data


    !************************************************test points**********************************************!
    write(6,*) 'u=',u(nx/2,ny/2) ,'v=',v(nx/2,ny/2),'rho',rho(nx/2,ny/2),'psi',psi(nx/2,ny/2) !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'u=',u(2,ny/2) ,'v=',v(2,ny/2),'rho',rho(2,ny/2)
    write(6,*) 'u=',u(1,ny/2) ,'v=',v(1,ny/2),'rho',rho(1,ny/2)
    write(6,*) 'running time: ', ts2-ts1, 'seconds'     

    open(101, file = 'psi.out', status = 'replace')
    do j=1,ny
        do i=1,nx
            write(101,*) psi(i,j)        
        enddo
    enddo
    close(101)
end program