program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    integer, parameter :: db=4 !kind(1.0)
    integer(kind=8) :: i,j,ll,l,dumm
    integer(kind=8) :: nx,ny,step,stamp,nlinks,nsteps,ngpus,ncontact
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2,radius
    real(kind=db) :: visc_LB,uu,udotc,omega,feq,geq
    real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8
    real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq1,pi2cssq2,pi2cssq0
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,fpc
    real(kind=db) :: addendum0,addendum1,addendum2,addendum3,addendum4,addendum5,addendum6,addendum7,addendum8
    real(kind=db) :: gaddendum0,gaddendum1,gaddendum2,gaddendum3,gaddendum4,gaddendum5,gaddendum6,gaddendum7,gaddendum8
    real(kind=db) :: psi_x,psi_y,mod_psi,mod_psi_sq,st_coeff,b0,b1,b2,beta,sigma
    real(kind=db) :: one_ov_nu2,one_ov_nu1,nu_avg,rtot,rprod
    real(kind=db) :: press_excess,max_press_excess,psid1,psid2,psid3,psid4
    real(kind=db) :: rnd_n1,rnd_n2,drhok1,drhok2
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p
    real(kind=db), allocatable, dimension(:,:) :: psi,rhoA,rhoB,u,v,pxx,pyy,pxy
    real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8
    real(kind=db), allocatable, dimension(:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8

    !lattice pars  
    nlinks=8 !pari!
    cssq=1.0_db/3.0_db
    !fluid 1
    tau=1.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu1=1.0_db/visc_LB
    !fluid2
    tau=1.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu2=1.0_db/visc_LB
    omega=1.0_db/tau
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif

    !*******************************user parameters**************************
    nx=512
    ny=512
    nsteps=2000
    stamp=10000
    fx=0.0_db*10.0**(-7)
    fy=0.0_db*10.0**(-5)
    !**********************************allocation****************************
    allocate(p(0:nlinks))
    allocate(f0(0:nx+1,0:ny+1),f1(0:nx+1,0:ny+1),f2(0:nx+1,0:ny+1),f3(0:nx+1,0:ny+1),f4(0:nx+1,0:ny+1))
    allocate(f5(0:nx+1,0:ny+1),f6(0:nx+1,0:ny+1),f7(0:nx+1,0:ny+1),f8(0:nx+1,0:ny+1))
    allocate(g0(0:nx+1,0:ny+1),g1(0:nx+1,0:ny+1),g2(0:nx+1,0:ny+1),g3(0:nx+1,0:ny+1),g4(0:nx+1,0:ny+1))
    allocate(g5(0:nx+1,0:ny+1),g6(0:nx+1,0:ny+1),g7(0:nx+1,0:ny+1),g8(0:nx+1,0:ny+1))
    allocate(psi(0:nx+1,0:ny+1),rhoA(1:nx,1:ny),rhoB(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny),pxx(1:nx,1:ny),pyy(1:nx,1:ny),pxy(1:nx,1:ny))
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
    pi2cssq0=p(0)/(2.0_db*cssq**2)
    pi2cssq1=p(1)/(2.0_db*cssq**2)
    pi2cssq2=p(5)/(2.0_db*cssq**2)
    ! chromodynamic
    beta=0.95_db
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
    rhoA=0.0_db     !is to be intended as a delta rho
    rhoB=0.0_db     !is to be intended as a delta rho
    press_excess=0.0_db
    max_press_excess=0.01 !2
    ncontact=0
    psid1=0.0_db
    psid2=0.0_db
    psid3=0.0_db
    psid4=0.0_db
    psi=-1.0_db
    radius=55
    do i=(nx/2-radius-5)-radius,(nx/2-radius-5) +radius
        do j=ny/2-radius,ny/2+radius
            if ((i-(nx/2-radius-5))**2+(j-ny/2)**2<=radius**2)then
                psi(i,j)=1.0_db
                u(i,j)=0.05
            endif
        enddo
    enddo
    do i=(nx/2+radius+5)-radius,(nx/2+radius+5) +radius
        do j=ny/2-radius,ny/2+radius
            if ((i-(nx/2+radius+5))**2+(j-ny/2)**2<=radius**2)then
                psi(i,j)=1.0_db
                u(i,j)=-0.05
            endif
        enddo
    enddo
    
    !
    rhoB=0.5*(1.0_db-psi(1:nx,1:ny))
    rhoA=1.0_db-rhoB
    write(*,*) rhoB(nx/2,ny/2),rhoA(nx/2,ny/2)
    !do ll=0,nlinks
    f0(1:nx,1:ny)=p(0)*rhoA(:,:)*(1.0-u(:,:)*u(:,:)*0.5/cssq)
    f1(1:nx,1:ny)=p(1)*rhoA(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    f2(1:nx,1:ny)=p(2)*rhoA(:,:)*(1.0+v(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    f3(1:nx,1:ny)=p(3)*rhoA(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    f4(1:nx,1:ny)=p(4)*rhoA(:,:)*(1.0-v(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    f5(1:nx,1:ny)=p(5)*rhoA(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    f6(1:nx,1:ny)=p(6)*rhoA(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    f7(1:nx,1:ny)=p(7)*rhoA(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    f8(1:nx,1:ny)=p(8)*rhoA(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    !
    g0(1:nx,1:ny)=p(0)*rhoB(:,:)*(1.0-u(:,:)*u(:,:)*0.5/cssq)
    g1(1:nx,1:ny)=p(1)*rhoB(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    g2(1:nx,1:ny)=p(2)*rhoB(:,:)*(1.0+v(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    g3(1:nx,1:ny)=p(3)*rhoB(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    g4(1:nx,1:ny)=p(4)*rhoB(:,:)*(1.0-v(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    g5(1:nx,1:ny)=p(5)*rhoB(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    g6(1:nx,1:ny)=p(6)*rhoB(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    g7(1:nx,1:ny)=p(7)*rhoB(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    g8(1:nx,1:ny)=p(8)*rhoB(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
    !
    !pause
    !*************************************check data ************************ 
    write(6,*) '*******************LB data*****************'
    write(6,*) 'tau',tau
    write(6,*) 'omega',omega
    write(6,*) 'visc',one_ov_nu1,one_ov_nu2
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


    !$acc data copy(p,rhoA,rhoB,u,v,pxx,pxy,pyy,f0,f1,f2,f3,f4,f5,f6,f7,f8,isfluid, &
    !$acc& g0,g1,g2,g3,g4,g5,g6,g7,g8,psi)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moment + neq pressor*********
        !$acc kernels 
        !$acc loop collapse(2) !private(uu,temp,udotc,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8) 
        do j=1,ny
            do i=1,nx
                if(isfluid(i,j).eq.1.or.isfluid(i,j).eq.0)then
                    rhoA(i,j) = f0(i,j)+f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)
                    rhoB(i,j)=  g0(i,j)+g1(i,j)+g2(i,j)+g3(i,j)+g4(i,j)+g5(i,j)+g6(i,j)+g7(i,j)+g8(i,j)
                    psi(i,j)= (rhoA(i,j)-rhoB(i,j))/(rhoA(i,j)+rhoB(i,j))
                    u(i,j) = (f1(i,j)+f5(i,j) +f8(i,j)-f3(i,j) -f6(i,j) -f7(i,j) + &
                              g1(i,j)+g5(i,j) +g8(i,j)-g3(i,j) -g6(i,j) -g7(i,j)) !/rho(i,j)
                    v(i,j) = (f5(i,j) +f2(i,j) +f6(i,j)-f7(i,j) -f4(i,j) -f8(i,j) + &
                              g5(i,j) +g2(i,j) +g6(i,j)-g7(i,j) -g4(i,j) -g8(i,j) )
                    ! non equilibrium pressor components
                    uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                    !1-3
                    udotc=u(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    rtot=rhoA(i,j)+rhoB(i,j)
                    fneq1=f1(i,j)+g1(i,j)-p(1)*(rtot + (temp + udotc))
                    fneq3=f3(i,j)+g3(i,j)-p(3)*(rtot + (temp - udotc))
                    !2-4
                    udotc=v(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq2=f2(i,j)+g2(i,j)-p(2)*(rtot + (temp + udotc))
                    fneq4=f4(i,j)+g4(i,j)-p(4)*(rtot + (temp - udotc))
                    !5-7
                    udotc=(u(i,j)+v(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq5=f5(i,j)+g5(i,j)-p(5)*(rtot + (temp + udotc))
                    fneq7=f7(i,j)+g7(i,j)-p(7)*(rtot + (temp - udotc))
                    !6-8
                    udotc=(-u(i,j)+v(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq6=f6(i,j)+g6(i,j)-p(6)*(rtot + (temp + udotc))
                    fneq8=f8(i,j)+g8(i,j)-p(8)*(rtot + (temp - udotc))
                    !
                    pxx(i,j)= fneq1 + fneq3 + fneq5 + fneq6 + fneq7 + fneq8
                    pyy(i,j)= fneq2 + fneq4 + fneq5 + fneq6 + fneq7 + fneq8
                    pxy(i,j)= fneq5 - fneq6 + fneq7 - fneq8
                    if(isfluid(i,j).eq.0)then
                        !no slip everywhere, always before fused: to be modified for generic pressure/velocity bcs
                        rhoA(i,j)=0.0_db
                        rhoB(i,j)=1.0_db
                        f0(i,j)=((p(0))  + pi2cssq0*(-cssq*pxx(i,j)-cssq*pyy(i,j)))*rhoA(i,j)/rtot
                        f1(i,j)=((p(3))  + pi2cssq1*(qxx*pxx(i,j)-cssq*pyy(i,j)))*rhoA(i,j)/rtot
                        f3(i,j)=((p(1))  + pi2cssq1*(qxx*pxx(i,j)-cssq*pyy(i,j)))*rhoA(i,j)/rtot
                        f2(i,j)=((p(4))  + pi2cssq1*(qyy*pyy(i,j)-cssq*pxx(i,j)))*rhoA(i,j)/rtot
                        f4(i,j)=((p(2))  + pi2cssq1*(qyy*pyy(i,j)-cssq*pxx(i,j)))*rhoA(i,j)/rtot
                        f5(i,j)=((p(7)) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)))*rhoA(i,j)/rtot
                        f7(i,j)=((p(5)) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)))*rhoA(i,j)/rtot
                        f6(i,j)=((p(8)) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)))*rhoA(i,j)/rtot
                        f8(i,j)=((p(6)) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)))*rhoA(i,j)/rtot
                        !
                        g0(i,j)=((p(0))  + pi2cssq0*(-cssq*pxx(i,j)-cssq*pyy(i,j)))*rhob(i,j)/rtot
                        g1(i,j)=((p(3))  + pi2cssq1*(qxx*pxx(i,j)-cssq*pyy(i,j)))*rhob(i,j)/rtot
                        g3(i,j)=((p(1))  + pi2cssq1*(qxx*pxx(i,j)-cssq*pyy(i,j)))*rhob(i,j)/rtot
                        g2(i,j)=((p(4))  + pi2cssq1*(qyy*pyy(i,j)-cssq*pxx(i,j)))*rhob(i,j)/rtot
                        g4(i,j)=((p(2))  + pi2cssq1*(qyy*pyy(i,j)-cssq*pxx(i,j)))*rhob(i,j)/rtot
                        g5(i,j)=((p(7)) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)))*rhob(i,j)/rtot
                        g7(i,j)=((p(5)) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)))*rhob(i,j)/rtot
                        g6(i,j)=((p(8)) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)))*rhob(i,j)/rtot
                        g8(i,j)=((p(6)) + pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)))*rhob(i,j)/rtot
                        psi(i,j)=-1.0_db
                    endif
                endif
            enddo
        enddo
        !********************collision + no slip + forcing: fused implementation*********
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
                rtot=0.0_db
                rtot=rhoA(i,j)+rhoB(i,j)
                rprod=rhoA(i,j)*rhoB(i,j)
                nu_avg=1.0_db/(rhoA(i,j)*one_ov_nu1/rtot + rhoB(i,j)*one_ov_nu2/rtot)
                omega=2.0_db/(6.0_db*nu_avg + 1.0_db)
                st_coeff=(9.0_db/4.0_db)*sigma*omega
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
                if(mod_psi>0.0001)then ! i'm sitting on the interface 
                    addendum0=-st_coeff*mod_psi*b0
                    addendum1=st_coeff*mod_psi*(p(1)*psi_x**2/mod_psi_sq - b1)
                    addendum2=st_coeff*mod_psi*(p(2)*psi_y**2/mod_psi_sq - b1)
                    addendum3=st_coeff*mod_psi*(p(3)*psi_x**2/mod_psi_sq - b1)
                    addendum4=st_coeff*mod_psi*(p(4)*psi_y**2/mod_psi_sq - b1)
                    addendum5=st_coeff*mod_psi*(p(5)*(psi_x+psi_y)**2/mod_psi_sq - b2)
                    addendum6=st_coeff*mod_psi*(p(6)*(-psi_x+psi_y)**2/mod_psi_sq - b2)
                    addendum7=st_coeff*mod_psi*(p(7)*(-psi_x-psi_y)**2/mod_psi_sq - b2)
                    addendum8=st_coeff*mod_psi*(p(8)*(psi_x-psi_y)**2/mod_psi_sq - b2)
                    !recoloring
                    gaddendum1=p(1)*(rtot)*(rprod*beta*psi_x/mod_psi/rtot**2)
                    gaddendum2=p(2)*(rtot)*(rprod*beta*psi_y/mod_psi/rtot**2)
                    gaddendum3=p(3)*(rtot)*(rprod*beta*(-psi_x/mod_psi)/rtot**2)
                    gaddendum4=p(4)*(rtot)*(rprod*beta*(-psi_y/mod_psi)/rtot**2)
                    gaddendum5=p(5)*(rtot)*(rprod*beta*(psi_x/mod_psi + psi_y/mod_psi)/rtot**2)
                    gaddendum6=p(6)*(rtot)*(rprod*beta*(-psi_x/mod_psi + psi_y/mod_psi)/rtot**2)
                    gaddendum7=p(7)*(rtot)*(rprod*beta*(-psi_x/mod_psi - psi_y/mod_psi)/rtot**2)
                    gaddendum8=p(8)*(rtot)*(rprod*beta*(psi_x/mod_psi - psi_y/mod_psi)/rtot**2)
                    ! se psi interfaccia vicina a 3-4-5lu è piu' grande del mio valore allora applico nci
                endif
                !************************* pressure excess **********!
                press_Excess=0.0_db
                psid1=0.0_db
                psid2=0.0_db
                psid3=0.0_db
                psid4=0.0_db
                ncontact=0
                if(psi(i,j).lt.-0.9_db)then
                    
                    if (psi(i+3,j).gt.-0.85 .and. psi(i+3,j).lt.0.0_db .and. psi(i-3,j).gt.-0.85 .and. psi(i-3,j).lt.0.0_db)then !&
                        !.and. psi(i-3*nint(psi_x/mod_psi),j-3*nint(psi_y/mod_psi)).gt.-0.85)then
                        !press_excess=max_press_excess*(2.0_db+psi(i+2*nint(psi_x/mod_psi),j+2*nint(psi_y/mod_psi)) + psi(i-2*nint(psi_x/mod_psi),j-2*nint(psi_y/mod_psi)))
                        psid1=psi(i+2,j) + psi(i-2,j)
                        ncontact=ncontact+1
                        ! write(*,*) float(ncontact)
                        ! stop
                    else
                        psid1=-2.0_db
                    endif
                    if (psi(i,j+3).gt.-0.85 .and. psi(i,j+3).lt.0.0_db .and. psi(i,j-3).gt.-0.85 .and. psi(i,j-3).lt.0.0_db)then
                        psid2=psi(i,j+2) + psi(i,j-2)
                        ncontact=ncontact+1
                    else
                        psid2=-2.0_db
                    endif
                    if (psi(i+2,j+2).gt.-0.85 .and. psi(i+2,j+2).lt.0.0_db .and. psi(i-2,j-2).gt.-0.85 .and. psi(i-2,j-2).lt.0.0_db)then
                        psid3=psi(i+1,j+1) + psi(i-1,j-1)
                        ncontact=ncontact+1
                    else
                        psid3=-2.0_db
                    endif
                    if (psi(i-2,j+2).gt.-0.85 .and. psi(i-2,j+2).lt.0.0_db .and. psi(i+2,j-2).gt.-0.85 .and. psi(i+2,j-2).lt.0.0_db)then
                        psid4=psi(i+1,j-1) + psi(i-1,j+1)
                        ncontact=ncontact+1
                    else
                        psid4=-2.0_db
                    endif
                    if(ncontact.gt.0)then
                        press_excess=max_press_excess*(8.0_db+psid1+psid2+psid3+psid4)/float(ncontact)
                    endif
                endif
                !*************************pressure excess**********!
                !regularized collision + perturbation + recolouring
                feq=p(0)*(rtot+press_excess-uu)
                fpc=feq + (1.0_db-omega)*pi2cssq0*(- cssq*pyy(i,j)-cssq*pxx(i,j))  + addendum0
                f0(i,j)=fpc*rhoA(i,j)/rtot 
                g0(i,j)=fpc*rhoB(i,j)/rtot
                !1
                udotc=u(i,j)/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(1)*(rtot+press_excess+(temp + udotc))
                fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j) - cssq*pyy(i,j) ) + fx*p(1)/cssq + addendum1
                f1(i+1,j)= fpc*rhoA(i,j)/rtot + gaddendum1
                g1(i+1,j)= fpc*rhoB(i,j)/rtot - gaddendum1
                !3
                feq=p(3)*(rtot+press_excess+(temp - udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j) - cssq*pyy(i,j)) - fx*p(3)/cssq + addendum3
                f3(i-1,j)= fpc*rhoA(i,j)/rtot + gaddendum3 
                g3(i-1,j)= fpc*rhoB(i,j)/rtot - gaddendum3 
                !2
                udotc=v(i,j)/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(2)*(rtot+press_excess+(temp + udotc))
                fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j)) + fy*p(2)/cssq + addendum2 !
                f2(i,j+1)= fpc*rhoA(i,j)/rtot + gaddendum2  
                g2(i,j+1)= fpc*rhoB(i,j)/rtot - gaddendum2   
                !4
                feq=p(4)*(rtot+press_excess+(temp - udotc))
                fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j)) - fy*p(4)/cssq + addendum4 
                f4(i,j-1)=fpc*rhoA(i,j)/rtot + gaddendum4 
                g4(i,j-1)=fpc*rhoB(i,j)/rtot - gaddendum4 
                !5
                udotc=(u(i,j)+v(i,j))/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(5)*(rtot+press_excess+(temp + udotc))
                fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + fx*p(5)/cssq + fy*p(5)/cssq + addendum5
                f5(i+1,j+1)=fpc*rhoA(i,j)/rtot + gaddendum5
                g5(i+1,j+1)=fpc*rhoB(i,j)/rtot - gaddendum5
                !7
                feq=p(7)*(rtot+press_excess+(temp - udotc))
                fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) - fx*p(7)/cssq - fy*p(7)/cssq + addendum7
                f7(i-1,j-1)=fpc*rhoA(i,j)/rtot + gaddendum7
                g7(i-1,j-1)=fpc*rhoB(i,j)/rtot - gaddendum7
                !6
                udotc=(-u(i,j)+v(i,j))/cssq
                temp = -uu + 0.5_db*udotc*udotc
                feq=p(6)*(rtot+press_excess+(temp + udotc))
                fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) - fx*p(6)/cssq + fy*p(6)/cssq + addendum6
                f6(i-1,j+1)= fpc*rhoA(i,j)/rtot + gaddendum6
                g6(i-1,j+1)= fpc*rhoB(i,j)/rtot - gaddendum6
                !8
                feq=p(8)*(rtot+press_excess+(temp - udotc))
                fpc=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) + fx*p(8)/cssq - fy*p(8)/cssq + addendum8
                f8(i+1,j-1)=fpc*rhoA(i,j)/rtot + gaddendum8 
                g8(i+1,j-1)=fpc*rhoB(i,j)/rtot - gaddendum8 
            enddo
        enddo
        !
        !******************************************call periodic bcs: always after fused************************
        !periodic along y
        !!$acc kernels 
        ! f5(1:nx,1)=f5(1:nx,ny-1)
        ! f2(1:nx,1)=f2(1:nx,ny-1)
        ! f6(1:nx,1)=f6(1:nx,ny-1)
        ! f8(1:nx,ny)=f8(1:nx,2)
        ! f4(1:nx,ny)=f4(1:nx,2)
        ! f7(1:nx,ny)=f7(1:nx,2)

        ! f1(1,1:ny)=f1(nx-1,1:ny)
        ! f5(1,1:ny)=f5(nx-1,1:ny)
        ! f8(1,1:ny)=f8(nx-1,1:ny)

        ! f3(nx,1:ny)=f3(2,1:ny)
        ! f6(nx,1:ny)=f6(2,1:ny)
        ! f7(nx,1:ny)=f7(2,1:ny)
        !$acc end kernels 

    enddo 
    call cpu_time(ts2)
    !$acc update host(rhoA,rhoB,psi,u,v)

    !$acc end data


    !************************************************test points**********************************************!
    write(6,*) 'u=',u(nx/2,ny/2) ,'v=',v(nx/2,ny/2),'rho',rhoA(nx/2,ny/2),'psi',psi(nx/2,ny/2) !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'u=',u(2,ny/2) ,'v=',v(2,ny/2),'rho',rhoA(2,ny/2)
    write(6,*) 'u=',u(1,ny/2) ,'v=',v(1,ny/2),'rho',rhoA(1,ny/2)
    write(6,*) 'running time: ', ts2-ts1, 'seconds'     

    open(101, file = 'psi.out', status = 'replace')
    do i=1,nx
        do j=1,ny
            write(101,*) psi(i,j)!sqrt(u(i,j)**2+v(i,j)**2) !sqrt(u(i,j)**2+v(i,j)**2)!rhoB(i,j)!sqrt(u(i,j)**2+v(i,j)**2)!rhoB(i,j)+rhoA(i,j)!sqrt(u(i,j)**2+v(i,j)**2)!!!rhoB(i,j)+rhoA(i,j)!sqrt(u(i,j)**2+v(i,j)**2)  
        enddo
    enddo
    close(101)

    open(102, file = 'rhoB.out', status = 'replace')
    do i=1,nx
        do j=1,ny
            write(102,*) rhoB(i,j)!sqrt(u(i,j)**2+v(i,j)**2) !sqrt(u(i,j)**2+v(i,j)**2)!rhoB(i,j)!sqrt(u(i,j)**2+v(i,j)**2)!rhoB(i,j)+rhoA(i,j)!sqrt(u(i,j)**2+v(i,j)**2)!!!rhoB(i,j)+rhoA(i,j)!sqrt(u(i,j)**2+v(i,j)**2)  
        enddo
    enddo
    close(102)

    open(103, file = 'mu.out', status = 'replace')
    do i=1,nx
        do j=1,ny
            write(103,*) sqrt(u(i,j)**2+v(i,j)**2)!(i,j)!sqrt(u(i,j)**2+v(i,j)**2) !sqrt(u(i,j)**2+v(i,j)**2)!rhoB(i,j)!sqrt(u(i,j)**2+v(i,j)**2)!rhoB(i,j)+rhoA(i,j)!sqrt(u(i,j)**2+v(i,j)**2)!!!rhoB(i,j)+rhoA(i,j)!sqrt(u(i,j)**2+v(i,j)**2)  
        enddo
    enddo
    close(103)
end program