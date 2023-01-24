program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    integer, parameter :: db=4 !kind(1.0)
    integer(kind=8) :: i,j,k
    integer(kind=8) :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
    real(kind=db) :: visc_LB,uu,udotc,omega,feq,geq,fpc
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp

    real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq17
    real(kind=db) :: fneq9,fneq10,fneq11,fneq12,fneq13,fneq14,fneq15,fneq16,fneq18
    real(kind=db) :: ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft17
    real(kind=db) :: ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft18
    real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
    real(kind=db) :: pi2cssq1,pi2cssq2,pi2cssq0
    real(kind=db) :: addendum0,addendum1,addendum2,addendum3,addendum4,addendum5,addendum6,addendum7,addendum8
    real(kind=db) :: addendum9,addendum10,addendum11,addendum12,addendum13,addendum14,addendum15,addendum16,addendum17,addendum18
    real(kind=db) :: gaddendum0,gaddendum1,gaddendum2,gaddendum3,gaddendum4,gaddendum5,gaddendum6,gaddendum7,gaddendum8
    real(kind=db) :: gaddendum9,gaddendum10,gaddendum11,gaddendum12,gaddendum13,gaddendum14,gaddendum15,gaddendum16,gaddendum17,gaddendum18
    real(kind=db) :: psi_x,psi_y,pais_z,mod_psi,mod_psi_sq,st_coeff,b0,b1,b2,beta,sigma
    real(kind=db) :: one_ov_nu2,one_ov_nu1,nu_avg,rtot,rprod
    real(kind=db) :: press_excess,max_press_excess

    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    real(kind=db), allocatable, dimension(:,:,:) :: psi,rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
    real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18
    real(kind=db), allocatable, dimension(:,:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8,g9
    real(kind=db), allocatable, dimension(:,:,:) :: g10,g11,g12,g13,g14,g15,g16,g17,g18

       
    nlinks=18 !pari!
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
    nz=512
    nsteps=100
    stamp=1000
    fx=0.0_db*10.0**(-7)
    fy=0.0_db*10.0**(-5)
    fz=0.0_db*10.0**(-5)

    allocate(f0(0:nx+1,0:ny+1,0:nz+1),f1(0:nx+1,0:ny+1,0:nz+1),f2(0:nx+1,0:ny+1,0:nz+1),f3(0:nx+1,0:ny+1,0:nz+1))
    allocate(f4(0:nx+1,0:ny+1,0:nz+1),f5(0:nx+1,0:ny+1,0:nz+1),f6(0:nx+1,0:ny+1,0:nz+1),f7(0:nx+1,0:ny+1,0:nz+1))
    allocate(f8(0:nx+1,0:ny+1,0:nz+1),f9(0:nx+1,0:ny+1,0:nz+1),f10(0:nx+1,0:ny+1,0:nz+1),f11(0:nx+1,0:ny+1,0:nz+1))
    allocate(f12(0:nx+1,0:ny+1,0:nz+1),f13(0:nx+1,0:ny+1,0:nz+1),f14(0:nx+1,0:ny+1,0:nz+1),f15(0:nx+1,0:ny+1,0:nz+1))
    allocate(f16(0:nx+1,0:ny+1,0:nz+1),f17(0:nx+1,0:ny+1,0:nz+1),f18(0:nx+1,0:ny+1,0:nz+1))
    allocate(rho(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
    allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
    allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz),psi(0:nx+1,0:ny+1,0:nz+1))
    allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
    
    !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)

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
    !****************************************hermite projection vars**********
    pi2cssq0=p0/(2.0_db*cssq**2)
    pi2cssq1=p1/(2.0_db*cssq**2)
    pi2cssq2=p2/(2.0_db*cssq**2)

    qxx=1.0_db-cssq
    qyy=1.0_db-cssq
    qzz=1.0_db-cssq
    qxy_7_8=1.0_db
    qxy_9_10=-1.0_db
    qxz_15_16=1.0_db
    qxz_17_18=-1.0_db
    qyz_11_12=1.0_db
    qyz_13_14=-1.0_db
    ! chromodynamic
    beta=0.95_db
    sigma=0.07_db
    st_coeff=(9.0_db/4.0_db)*sigma*omega
    b0=-2.0_db/9.0_db
    b1=1.0_db/54.0_db
    b2=10_db/27.0_db
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    w=0.0_db
    rhoA=0.0_db  !total density
    rhoB=0.0_db  !total density
    press_excess=0.0_db
    max_press_excess=0.05
    psi=-1.0_db
    radius=14
    do i=35-radius,35+radius
        do j=ny/2-radius,ny/2+radius
            do k=nz/2-radius,nz/2+radius
                if ((i-35)**2+(j-ny/2)**2+(k-nz/2)**2<=radius**2)then
                    psi(i,j,k)=1.0_db
                endif
            enddo
        enddo
    enddo
    do i=66-radius,66+radius
        do j=ny/2-radius,ny/2+radius
            do k=nz/2-radius,nz/2+radius
                if ((i-66)**2+(j-ny/2)**2+(k-nz/2)**2<=radius**2)then
                    psi(i,j,k)=1.0_db
                endif
            enddo
        enddo
    enddo
    rhoB=0.5*(1.0_db-psi(1:nx,1:ny))
    rhoA=1.0_db-rhoB
    write(*,*) rhoB(nx/2,ny/2),rhoA(nx/2,ny/2)
    !!***************!!!
    f0(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p0
    f1(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p1
    f2(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p1
    f3(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p1
    f4(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p1
    f5(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p1
    f6(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p1
    f7(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f8(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f9(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f10(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f11(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f12(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f13(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f14(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f15(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f16(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f17(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    f18(1:nx,1:ny,1:nz)=rhoA(1:nx,1:ny,1:nz)*p2
    !
    g0(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p0
    g1(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1
    g2(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1
    g3(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1
    g4(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1
    g5(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1
    g6(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1
    g7(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g8(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g9(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g10(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g11(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g12(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g13(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g14(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g15(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g16(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g17(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
    g18(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2
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
    write(6,*) 'beta',beta
    write(6,*) 'sigma',sigma
    write(6,*) 'surface_tens_coeff',st_coeff
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
    !$acc data copy(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,isfluid,p0,p1,p2,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,psi, &
             !$acc& g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moments collision bbck + forcing************************ 
        
        !$acc kernels 
        !$acc loop collapse (3) 
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1.or.isfluid(i,j,k).eq.0)then
                        ft0=f0(i,j,k)+g0(i,j,k)
                        ft1=f1(i,j,k)+g1(i,j,k)
                        ft2=f2(i,j,k)+g2(i,j,k)
                        ft3=f3(i,j,k)+g3(i,j,k)
                        ft4=f4(i,j,k)+g4(i,j,k)
                        ft5=f5(i,j,k)+g5(i,j,k)
                        ft6=f6(i,j,k)+g6(i,j,k)
                        ft7=f7(i,j,k)+g7(i,j,k)
                        ft8=f8(i,j,k)+g8(i,j,k)
                        ft9=f9(i,j,k)+g9(i,j,k)
                        ft10=f10(i,j,k)+g10(i,j,k)
                        ft11=f11(i,j,k)+g11(i,j,k)
                        ft12=f12(i,j,k)+g12(i,j,k)
                        ft13=f13(i,j,k)+g13(i,j,k)
                        ft14=f14(i,j,k)+g14(i,j,k)
                        ft15=f15(i,j,k)+g15(i,j,k)
                        ft16=f16(i,j,k)+g16(i,j,k)
                        ft17=f17(i,j,k)+g17(i,j,k)
                        ft18=f18(i,j,k)+g18(i,j,k)
                        rhoA(i,j,k) = f0(i,j,k)+f1(i,j,k)+f2(i,j,k)+f3(i,j,k)+f4(i,j,k)+f5(i,j,k) &
                            +f6(i,j,k)+f7(i,j,k)+f8(i,j,k)+f9(i,j,k)+f10(i,j,k)+f11(i,j,k) &
                            +f12(i,j,k)+f13(i,j,k)+f14(i,j,k)+f15(i,j,k)+f16(i,j,k)+f17(i,j,k) &
                            +f18(i,j,k)

                        rhoB(i,j,k) = g0(i,j,k)+g1(i,j,k)+g2(i,j,k)+g3(i,j,k)+g4(i,j,k)+g5(i,j,k) &
                            +g6(i,j,k)+g7(i,j,k)+g8(i,j,k)+g9(i,j,k)+g10(i,j,k)+g11(i,j,k) &
                            +g12(i,j,k)+g13(i,j,k)+g14(i,j,k)+g15(i,j,k)+g16(i,j,k)+g17(i,j,k) &
                            +g18(i,j,k)

                        u(i,j,k) = (ft1+ft7+ft9+ft15+ft18)-(ft2+ft8+ft10+ft16+ft17) 
                        
                        v(i,j,k) = (ft3+ft7+ft10+ft11+ft13)-(ft4+ft8+ft9+ft12+ft14)

                        w(i,j,k) = (ft5+ft11+ft14+ft15+ft17)-(ft6+ft12+ft13+ft16+ft18)
                        
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !1-2
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        rtot=rhoA(i,j)+rhoB(i,j)
                        fneq1=ft1-p1*(rtot+(temp + udotc))
                        fneq2=ft2-p1*(rtot+(temp - udotc))
                        !3-4
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq3=ft3-p1*(rtot+(temp + udotc))
                        fneq4=ft4-p1*(rtot+(temp - udotc))
                        !5-6
                        udotc=w(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq5=ft5-p1*(rtot+(temp + udotc))
                        fneq6=ft6-p1*(rtot+(temp - udotc))
                        !7-8
                        udotc=(u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq7=ft7-p2*(rtot+(temp + udotc))
                        fneq8=ft8-p2*(rtot+(temp - udotc))
                        !10-9
                        udotc=(-u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq10=ft10-p2*(rtot+(temp + udotc))
                        fneq9=ft9-p2*(rtot+(temp - udotc))
                        !11-12
                        udotc=(v(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq11=ft11-p2*(rtot+(temp + udotc))
                        fneq12=ft12-p2*(rtot+(temp - udotc))
                        !13-14
                        udotc=(v(i,j,k)-w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq13=ft13 - p2*(rtot+(temp + udotc))
                        fneq14=ft14 - p2*(rtot+(temp - udotc))
                        !15-16
                        udotc=(u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq15=ft15-p2*(rtot+(temp + udotc))
                        fneq16=ft16-p2*(rtot+(temp - udotc))
                        !17-18
                        udotc=(-u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq17=ft17-p2*(rtot+(temp + udotc))
                        fneq18=ft18-p2*(rtot+(temp - udotc))
                        pxx(i,j,k)=fneq1+fneq2+fneq7+fneq8+fneq9+fneq10+fneq15+fneq16+fneq17+fneq18
                        pyy(i,j,k)=fneq3+fneq4+fneq7+fneq8+fneq9+fneq10+fneq11+fneq12+fneq13+fneq14
                        pzz(i,j,k)=fneq5+fneq6+fneq11+fneq12+fneq13+fneq14+fneq15+fneq16+fneq17+fneq18
                        pxy(i,j,k)= fneq7+fneq8-fneq9-fneq10
                        pxz(i,j,k)=fneq15+fneq16-fneq17-fneq18
                        pyz(i,j,k)=fneq11+fneq12-fneq13-fneq14
                    endif
                    !no slip everywhere, always before fused: to be modified for generic pressure/velocity bcs
                    if(isfluid(i,j,k).eq.0)then
                        f0(i,j,k)=(p0 + pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                        f1(i,j,k)=(p1 + pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                        f2(i,j,k)=(p1 + pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                        f3(i,j,k)=(p1 + pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                        f4(i,j,k)=(p1 + pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                        f5(i,j,k)=(p1 + pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))))*rhoA(i,j,k)/rtot
                        f6(i,j,k)=(p1 + pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))))*rhoA(i,j,k)/rtot

                        f7(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j)+2.0_db*qxy_7_8*pxy(i,j,k)))*rhoA(i,j,k)/rtot
                        f8(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j)+2.0_db*qxy_7_8*pxy(i,j,k)))*rhoA(i,j,k)/rtot
                        f9(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j)+2.0_db*qxy_9_10*pxy(i,j,k)))*rhoA(i,j,k)/rtot
                        f10(i,j,k)=(p2 +pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j)+2.0_db*qxy_9_10*pxy(i,j,k)))*rhoA(i,j,k)/rtot

                        f11(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j)+2.0_db*qyz_11_12*pyz(i,j,k)))*rhoA(i,j,k)/rtot
                        f12(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j)+2.0_db*qyz_11_12*pyz(i,j,k)))*rhoA(i,j,k)/rtot
                        f13(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j)+2.0_db*qyz_13_14*pyz(i,j,k)))*rhoA(i,j,k)/rtot
                        f14(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j)+2.0_db*qyz_13_14*pyz(i,j,k)))*rhoA(i,j,k)/rtot

                        f15(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j)+2.0_db*qxz_15_16*pxz(i,j,k)))*rhoA(i,j,k)/rtot
                        f16(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j)+2.0_db*qxz_15_16*pxz(i,j,k)))*rhoA(i,j,k)/rtot
                        f17(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j)+2.0_db*qxz_17_18*pxz(i,j,k)))*rhoA(i,j,k)/rtot
                        f18(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j)+2.0_db*qxz_17_18*pxz(i,j,k)))*rhoA(i,j,k)/rtot
                        !
                        g0(i,j,k)=(p0 + pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                        g1(i,j,k)=(p1 + pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                        g2(i,j,k)=(p1 + pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                        g3(i,j,k)=(p1 + pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                        g4(i,j,k)=(p1 + pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                        g5(i,j,k)=(p1 + pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))))*rhob(i,j,k)/rtot
                        g6(i,j,k)=(p1 + pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))))*rhob(i,j,k)/rtot

                        g7(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j)+2.0_db*qxy_7_8*pxy(i,j,k)))*rhob(i,j,k)/rtot
                        g8(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j)+2.0_db*qxy_7_8*pxy(i,j,k)))*rhob(i,j,k)/rtot
                        g9(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j)+2.0_db*qxy_9_10*pxy(i,j,k)))*rhob(i,j,k)/rtot
                        g10(i,j,k)=(p2 +pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j)+2.0_db*qxy_9_10*pxy(i,j,k)))*rhob(i,j,k)/rtot

                        g11(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j)+2.0_db*qyz_11_12*pyz(i,j,k)))*rhob(i,j,k)/rtot
                        g12(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j)+2.0_db*qyz_11_12*pyz(i,j,k)))*rhob(i,j,k)/rtot
                        g13(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j)+2.0_db*qyz_13_14*pyz(i,j,k)))*rhob(i,j,k)/rtot
                        g14(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j)+2.0_db*qyz_13_14*pyz(i,j,k)))*rhob(i,j,k)/rtot

                        g15(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j)+2.0_db*qxz_15_16*pxz(i,j,k)))*rhob(i,j,k)/rtot
                        g16(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j)+2.0_db*qxz_15_16*pxz(i,j,k)))*rhob(i,j,k)/rtot
                        g17(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j)+2.0_db*qxz_17_18*pxz(i,j,k)))*rhob(i,j,k)/rtot
                        g18(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j)+2.0_db*qxz_17_18*pxz(i,j,k)))*rhob(i,j,k)/rtot
                    endif
                enddo
            enddo
        enddo
        !!$acc end kernels
        !!$acc kernels
        !$acc loop collapse (3) 
        do k=1,nz
            do j=1,ny
                do i=1,nx
                        !
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq   
                        !!*******************************chromodynamics***************************!
                        psi_x=(1.0_db/cssq)*(p1*(psi(i+1,j,k)-psi(i-1,j,k)) + p2*(psi(i+1,j+1,k)+psi(i+1,j-1,k)+psi(i+1,j,k+1)+psi(i+1,j,k-1)) &
                                                -p2*(psi(i-1,j+1,k)+psi(i-1,j-1,k)+psi(i-1,j,k+1)+psi(i-1,j,k-1)))
                        psi_y=(1.0_db/cssq)*(p1*(psi(i,j+1,k)-psi(i,j-1,k)) + p2*(psi(i+1,j+1,k)+psi(i-1,j+1,k)+psi(i,j+1,k+1)+psi(i,j+1,k-1)) &
                                                -p2*(psi(i+1,j-1,k)+psi(i-1,j-1,k)+psi(i,j-1,k+1)+psi(i,j-1,k-1)))
                        psi_z=(1.0_db/cssq)*(p1*(psi(i,j+1,k)-psi(i,j-1,k)) + p2*(psi(i+1,j+1,k)+psi(i-1,j+1,k)+psi(i,j+1,k+1)+psi(i,j+1,k-1)) &
                                                -p2*(psi(i+1,j-1,k)+psi(i-1,j-1,k)+psi(i,j-1,k+1)+psi(i,j-1,k-1)))
                        !!************************************************************************!
                        !0
                        feq=p0*(rho(i,j,k)-uu)
                        f0(i,j,k)=feq + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k)))
                        
                        !1
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq1=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                        f1(i+1,j,k)=feq + fneq1 + fx*p1dcssq
                        
                        !2
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f2(i-1,j,k)=feq + fneq1 - fx*p1dcssq
                        
                        !3
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq3=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
                        f3(i,j+1,k)=feq+fneq3 + fy*p1dcssq
                        
                        !4
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f4(i,j-1,k)=feq+fneq3 - fy*p1dcssq
                        
                        !7
                        udotc=(u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq7=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k))
                        f7(i+1,j+1,k)=feq + fneq7 + (fx+fy)*p2dcssq 
                        
                        !8
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f8(i-1,j-1,k)=feq + fneq7 - (fx+fy)*p2dcssq
                        
                        !10
                        udotc=(-u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq10=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k))
                        f10(i-1,j+1,k)=feq+fneq10 +(fy-fx)*p2dcssq
                        
                        !9
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f9(i+1,j-1,k)=feq+fneq10 + (fx-fy)*p2dcssq

                        !5
                        udotc=w(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq5=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
                        f5(i,j,k+1)=feq+fneq5 + fz*p1dcssq
                        
                        !6
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f6(i,j,k-1)=feq+fneq5 - fz*p1dcssq

                        !15
                        udotc=(u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq15=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k))
                        f15(i+1,j,k+1)=feq+fneq15 + (fx+fz)*p2dcssq 
                        
                        !16
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f16(i-1,j,k-1)=feq+fneq15 - (fx+fz)*p2dcssq

                        !17
                        udotc=(-u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq17=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k))
                        f17(i-1,j,k+1)=feq+fneq17 +(fz-fx)*p2dcssq
                        
                        !18
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f18(i+1,j,k-1)=feq+fneq17 + (fx-fz)*p2dcssq

                        !11
                        udotc=(v(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq11=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k))
                        f11(i,j+1,k+1)=feq+fneq11+(fy+fz)*p2dcssq
                        
                        !12
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f12(i,j-1,k-1)=feq+fneq11 - (fy+fz)*p2dcssq

                        !13
                        udotc=(v(i,j,k)-w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq13=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k))
                        f13(i,j+1,k-1)=feq+fneq13 + (fy-fz)*p2dcssq
                        
                        !14
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f14(i,j-1,k+1)=feq+fneq13 + (fz-fy)*p2dcssq
                enddo
            enddo
        enddo
        !!$acc end kernels
        !******************************************call bcs************************
        !periodic along y
        !x=1     
        !!$acc kernels 
        f1(2,:,:)=f1(nx,:,:)
        f7(2,:,:)=f7(nx,:,:)
        f9(2,:,:)=f9(nx,:,:)
        f15(2,:,:)=f15(nx,:,:)
        f18(2,:,:)=f18(nx,:,:)
        !x=nx 
        f2(nx-1,:,:)=f2(1,:,:)
        f8(nx-1,:,:)=f8(1,:,:)
        f10(nx-1,:,:)=f10(1,:,:)
        f16(nx-1,:,:)=f16(1,:,:)
        f17(nx-1,:,:)=f17(1,:,:)

        !y=1
        f3(:,2,:)=f3(:,ny,:)
        f7(:,2,:)=f7(:,ny,:)
        f10(:,2,:)=f10(:,ny,:)
        f11(:,2,:)=f11(:,ny,:)
        f13(:,2,:)=f13(:,ny,:)
       
        !y=ny
        f4(:,ny-1,:)=f4(:,1,:)
        f8(:,ny-1,:)=f8(:,1,:)
        f9(:,ny-1,:)=f9(:,1,:)
        f12(:,ny-1,:)=f12(:,1,:)
        f14(:,ny-1,:)=f14(:,1,:)
        !$acc end kernels 
        
        !
    enddo 
    call cpu_time(ts2)
    !$acc update host(rho,u,v,w)
    !$acc end data
    

    write(6,*) 'u=',u(nx/2,ny/2,nz/2),'v=',v(nx/2,ny/2,nz/2),'w=',w(nx/2,ny/2,nz/2),'rho=',rho(nx/2,ny/2,nz/2)
    write(6,*) 'u=',u(nx/2,ny/2,1),'v=',v(nx/2,ny/2,1),'w=',w(nx/2,ny/2,1),'rho=',rho(nx/2,ny/2,1)
    write(6,*) 'You''ve just wasted ', ts2-ts1, ' s of your life time' 

    
end program