program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    !************************************block of vars*******************************************!
        integer, parameter :: db=4 !kind(1.0)
        integer(kind=8) :: i,j,k
        integer(kind=8) :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus
        integer,save :: iframe=0
        
        real(kind=db),parameter :: pi_greek=3.14159265359793234626433
        
        real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
        real(kind=db) :: visc_LB,uu,udotc,omega,feq,geq,fpc
        real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp
        real(kind=db) :: radius

        real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq17
        real(kind=db) :: fneq9,fneq10,fneq11,fneq12,fneq13,fneq14,fneq15,fneq16,fneq18
        real(kind=db) :: ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft17
        real(kind=db) :: ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft18
        real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
        real(kind=db) :: pi2cssq1,pi2cssq2,pi2cssq0
        real(kind=db) :: psid1,psid2,psid3,psid4,psid5,psid6,psid7,psid8,psid9,ncontact
        real(kind=db) :: addendum0,addendum1,addendum2,addendum3,addendum4,addendum5,addendum6,addendum7,addendum8
        real(kind=db) :: addendum9,addendum10,addendum11,addendum12,addendum13,addendum14,addendum15,addendum16,addendum17,addendum18
        real(kind=db) :: gaddendum1,gaddendum2,gaddendum3,gaddendum4,gaddendum5,gaddendum6,gaddendum7,gaddendum8
        real(kind=db) :: gaddendum9,gaddendum10,gaddendum11,gaddendum12,gaddendum13,gaddendum14,gaddendum15,gaddendum16,gaddendum17,gaddendum18
        real(kind=db) :: psi_x,psi_y,psi_z,mod_psi,mod_psi_sq,st_coeff,b0,b1,b2,beta,sigma
        real(kind=db) :: one_ov_nu2,one_ov_nu1,nu_avg,rtot,rprod
        real(kind=db) :: press_excess,max_press_excess

        integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
        real(kind=db), allocatable, dimension(:,:,:) :: psi,rhoA,rhoB,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
        real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
        real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18
        real(kind=db), allocatable, dimension(:,:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8,g9
        real(kind=db), allocatable, dimension(:,:,:) :: g10,g11,g12,g13,g14,g15,g16,g17,g18

       
   
    
    
    !*********************************** lattice parameters**************************************!
        nlinks=18 !pari!
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

    !**************************************user parameters**************************
        nx=400
        ny=400
        nz=400
        nsteps=10
        stamp=5
        fx=0.0_db*10.0**(-7)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
    !*****************************************allocation*******************************************************
        allocate(f0(0:nx+1,0:ny+1,0:nz+1),f1(0:nx+1,0:ny+1,0:nz+1),f2(0:nx+1,0:ny+1,0:nz+1),f3(0:nx+1,0:ny+1,0:nz+1))
        allocate(f4(0:nx+1,0:ny+1,0:nz+1),f5(0:nx+1,0:ny+1,0:nz+1),f6(0:nx+1,0:ny+1,0:nz+1),f7(0:nx+1,0:ny+1,0:nz+1))
        allocate(f8(0:nx+1,0:ny+1,0:nz+1),f9(0:nx+1,0:ny+1,0:nz+1),f10(0:nx+1,0:ny+1,0:nz+1),f11(0:nx+1,0:ny+1,0:nz+1))
        allocate(f12(0:nx+1,0:ny+1,0:nz+1),f13(0:nx+1,0:ny+1,0:nz+1),f14(0:nx+1,0:ny+1,0:nz+1),f15(0:nx+1,0:ny+1,0:nz+1))
        allocate(f16(0:nx+1,0:ny+1,0:nz+1),f17(0:nx+1,0:ny+1,0:nz+1),f18(0:nx+1,0:ny+1,0:nz+1))
        allocate(g0(0:nx+1,0:ny+1,0:nz+1),g1(0:nx+1,0:ny+1,0:nz+1),g2(0:nx+1,0:ny+1,0:nz+1),g3(0:nx+1,0:ny+1,0:nz+1))
        allocate(g4(0:nx+1,0:ny+1,0:nz+1),g5(0:nx+1,0:ny+1,0:nz+1),g6(0:nx+1,0:ny+1,0:nz+1),g7(0:nx+1,0:ny+1,0:nz+1))
        allocate(g8(0:nx+1,0:ny+1,0:nz+1),g9(0:nx+1,0:ny+1,0:nz+1),g10(0:nx+1,0:ny+1,0:nz+1),g11(0:nx+1,0:ny+1,0:nz+1))
        allocate(g12(0:nx+1,0:ny+1,0:nz+1),g13(0:nx+1,0:ny+1,0:nz+1),g14(0:nx+1,0:ny+1,0:nz+1),g15(0:nx+1,0:ny+1,0:nz+1))
        allocate(g16(0:nx+1,0:ny+1,0:nz+1),g17(0:nx+1,0:ny+1,0:nz+1),g18(0:nx+1,0:ny+1,0:nz+1))
        allocate(rhoA(1:nx,1:ny,1:nz),rhoB(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
        allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
        allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz),psi(0:nx+1,0:ny+1,0:nz+1))
        allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
    
    !***************************************lattice vars*************************************!
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)

        p0=(1.0_db/3.0_db)
        p1=(1.0_db/18.0_db)
        p2=(1.0_db/36.0_db)
        p1dcssq=p1/cssq
        p2dcssq=p2/cssq
    !****************************************geometry************************
        isfluid=1
        isfluid(1,:,:)=0 !left
        isfluid(nx,:,:)=0 !right
        isfluid(:,1,:)=0 !front 
        isfluid(:,ny,:)=0 !rear
        isfluid(:,:,1)=0 !bottom
        isfluid(:,:,nz)=0 !top
    !********************************hermite projection vars**********
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
    !*********************************chromodynamics vars*****************************
        ! 
        beta=0.95_db
        sigma=0.02_db
        st_coeff=(9.0_db/4.0_db)*sigma*omega
        b0=-1.0_db/3.0_db
        b1=1.0_db/18.0_db
        b2=1_db/36.0_db
    !********************************initialization of macrovars ************************    
        u=0.0_db
        v=0.0_db
        w=0.0_db
        rhoA(1:nx,1:ny,1:nz)=0.0_db  !total density
        rhoB(1:nx,1:ny,1:nz)=0.0_db  !total density
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
        rhoB=0.5*(1.0_db-psi(1:nx,1:ny,1:nz))
        rhoA=1.0_db-rhoB
        write(*,*) rhoB(nx/2,ny/2,nz/2),rhoA(nx/2,ny/2,nz/2)
    !*************************************set distros************************!
        f0(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p0
        f1(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1
        f2(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1
        f3(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1
        f4(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1
        f5(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1
        f6(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1
        f7(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f8(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f9(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f10(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f11(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f12(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f13(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f14(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f15(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f16(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f17(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
        f18(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2
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
    
    !***************************************check data ************************ 
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
    !*****************************************copy data on gpu*********************************************!
    !$acc data copy(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,isfluid,p0,p1,p2,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rhoA,rhoB,u,v,w,psi, &
             !$acc& g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18)
    
    
    !**********************************************************************!
    call cpu_time(ts1)
    !*************************************main loop*************************!
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

                        psi(i,j,k)= (rhoA(i,j,k)-rhoB(i,j,k))/(rhoA(i,j,k)+rhoB(i,j,k))

                        u(i,j,k) = (ft1+ft7+ft9+ft15+ft18)-(ft2+ft8+ft10+ft16+ft17) 
                        
                        v(i,j,k) = (ft3+ft7+ft10+ft11+ft13)-(ft4+ft8+ft9+ft12+ft14)

                        w(i,j,k) = (ft5+ft11+ft14+ft15+ft17)-(ft6+ft12+ft13+ft16+ft18)
                        
                        
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !1-2
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        rtot=rhoA(i,j,k)+rhoB(i,j,k)
                        fneq1=ft1-p1*(rtot+(temp + udotc))
                        fneq2=ft2-p1*(rtot+(temp - udotc))
                        !3-4
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rtot+(temp + udotc))
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
                        !no slip everywhere, always before fused: to be modified for generic pressure/velocity bcs
                        if(isfluid(i,j,k).eq.0)then
                            f0(i,j,k)=(p0 + pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                            f1(i,j,k)=(p1 + pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                            f2(i,j,k)=(p1 + pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                            f3(i,j,k)=(p1 + pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                            f4(i,j,k)=(p1 + pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))))*rhoA(i,j,k)/rtot
                            f5(i,j,k)=(p1 + pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))))*rhoA(i,j,k)/rtot
                            f6(i,j,k)=(p1 + pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))))*rhoA(i,j,k)/rtot

                            f7(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k)))*rhoA(i,j,k)/rtot
                            f8(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k)))*rhoA(i,j,k)/rtot
                            f9(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k)))*rhoA(i,j,k)/rtot
                            f10(i,j,k)=(p2 +pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k)))*rhoA(i,j,k)/rtot

                            f11(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k)))*rhoA(i,j,k)/rtot
                            f12(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k)))*rhoA(i,j,k)/rtot
                            f13(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k)))*rhoA(i,j,k)/rtot
                            f14(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k)))*rhoA(i,j,k)/rtot

                            f15(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k)))*rhoA(i,j,k)/rtot
                            f16(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k)))*rhoA(i,j,k)/rtot
                            f17(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k)))*rhoA(i,j,k)/rtot
                            f18(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k)))*rhoA(i,j,k)/rtot
                            !
                            g0(i,j,k)=(p0 + pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                            g1(i,j,k)=(p1 + pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                            g2(i,j,k)=(p1 + pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                            g3(i,j,k)=(p1 + pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                            g4(i,j,k)=(p1 + pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))))*rhob(i,j,k)/rtot
                            g5(i,j,k)=(p1 + pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))))*rhob(i,j,k)/rtot
                            g6(i,j,k)=(p1 + pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))))*rhob(i,j,k)/rtot

                            g7(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k)))*rhob(i,j,k)/rtot
                            g8(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k)))*rhob(i,j,k)/rtot
                            g9(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k)))*rhob(i,j,k)/rtot
                            g10(i,j,k)=(p2 +pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k)))*rhob(i,j,k)/rtot

                            g11(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k)))*rhob(i,j,k)/rtot
                            g12(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k)))*rhob(i,j,k)/rtot
                            g13(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k)))*rhob(i,j,k)/rtot
                            g14(i,j,k)=(p2 + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k)))*rhob(i,j,k)/rtot

                            g15(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k)))*rhob(i,j,k)/rtot
                            g16(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k)))*rhob(i,j,k)/rtot
                            g17(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k)))*rhob(i,j,k)/rtot
                            g18(i,j,k)=(p2 + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k)))*rhob(i,j,k)/rtot
                        endif
                    endif
                    
                enddo
            enddo
        enddo
        !$acc loop collapse (3) 
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq   
                    !!*******************************chromodynamics***************************!
                    psi_x=(1.0_db/cssq)*(p1*(psi(i+1,j,k)-psi(i-1,j,k)) + p2*(psi(i+1,j+1,k)+psi(i+1,j-1,k)+psi(i+1,j,k+1)+psi(i+1,j,k-1)) &
                                            -p2*(psi(i-1,j+1,k)+psi(i-1,j-1,k)+psi(i-1,j,k+1)+psi(i-1,j,k-1)))


                    psi_y=(1.0_db/cssq)*(p1*(psi(i,j+1,k)-psi(i,j-1,k)) + p2*(psi(i+1,j+1,k)+psi(i-1,j+1,k)+psi(i,j+1,k+1)+psi(i,j+1,k-1)) &
                                            -p2*(psi(i+1,j-1,k)+psi(i-1,j-1,k)+psi(i,j-1,k+1)+psi(i,j-1,k-1)))


                    psi_z=(1.0_db/cssq)*(p1*(psi(i,j,k+1)-psi(i,j,k-1)) + p2*(psi(i+1,j,k+1)+psi(i-1,j,k+1)+psi(i,j+1,k+1)+psi(i,j-1,k+1)) &
                                            -p2*(psi(i+1,j,k-1)+psi(i-1,j,k-1)+psi(i,j+1,k-1)+psi(i,j-1,k-1)))
                    mod_psi=sqrt(psi_x**2+psi_y**2+psi_z**2)

                    mod_psi_sq=psi_x**2+psi_y**2 +psi_z**2 
                    
                    rtot=0.0_db
                    
                    rtot=rhoA(i,j,k)+rhoB(i,j,k)
                    
                    rprod=rhoA(i,j,k)*rhoB(i,j,k)
                    
                    nu_avg=1.0_db/(rhoA(i,j,k)*one_ov_nu1/rtot + rhoB(i,j,k)*one_ov_nu2/rtot)
                    
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
                    addendum9=0.0_db
                    addendum10=0.0_db
                    addendum11=0.0_db
                    addendum12=0.0_db
                    addendum13=0.0_db
                    addendum14=0.0_db
                    addendum15=0.0_db
                    addendum16=0.0_db
                    addendum17=0.0_db
                    addendum18=0.0_db
                    gaddendum1=0.0_db
                    gaddendum2=0.0_db
                    gaddendum3=0.0_db
                    gaddendum4=0.0_db
                    gaddendum5=0.0_db
                    gaddendum6=0.0_db
                    gaddendum7=0.0_db
                    gaddendum8=0.0_db 
                    gaddendum9=0.0_db
                    gaddendum10=0.0_db
                    gaddendum11=0.0_db
                    gaddendum12=0.0_db
                    gaddendum13=0.0_db
                    gaddendum14=0.0_db
                    gaddendum15=0.0_db
                    gaddendum16=0.0_db
                    gaddendum17=0.0_db
                    gaddendum18=0.0_db
                    if (mod_psi>0.001)then
                        addendum0=-st_coeff*mod_psi*b0
                        addendum1=st_coeff*mod_psi*(p1*psi_x**2/mod_psi_sq - b1)
                        addendum2=addendum1
                        addendum3=st_coeff*mod_psi*(p1*psi_y**2/mod_psi_sq - b1)
                        addendum4=addendum3
                        addendum5=st_coeff*mod_psi*(p1*psi_z**2/mod_psi_sq - b1)
                        addendum6=addendum5
                        !
                        addendum7=st_coeff*mod_psi*(p2*(psi_x+psi_y)**2/mod_psi_sq - b2)
                        addendum8=addendum7
                        addendum9=st_coeff*mod_psi*(p2*(psi_x-psi_y)**2/mod_psi_sq - b2)
                        addendum10=addendum9
                        addendum11=st_coeff*mod_psi*(p2*(psi_y+psi_z)**2/mod_psi_sq - b2)
                        addendum12=addendum11
                        addendum13=st_coeff*mod_psi*(p2*(psi_y-psi_z)**2/mod_psi_sq - b2)
                        addendum14=addendum13
                        addendum15=st_coeff*mod_psi*(p2*(psi_x+psi_z)**2/mod_psi_sq - b2)
                        addendum16=addendum15
                        addendum17=st_coeff*mod_psi*(p2*(-psi_x+psi_z)**2/mod_psi_sq - b2)
                        addendum18=addendum17
                        !
                        gaddendum1=p1*(rtot)*(rprod*beta*psi_x/mod_psi/rtot**2)
                        gaddendum2=-gaddendum1
                        gaddendum3=p1*(rtot)*(rprod*beta*psi_y/mod_psi/rtot**2)
                        gaddendum4=-gaddendum3
                        gaddendum5=p1*(rtot)*(rprod*beta*psi_z/mod_psi/rtot**2)
                        gaddendum6=-gaddendum5
                        !
                        gaddendum7=p2*(rtot)*(rprod*beta*(psi_x+psi_y)/mod_psi/rtot**2)
                        gaddendum8=-gaddendum7
                        gaddendum9=p2*(rtot)*(rprod*beta*(psi_x-psi_y)/mod_psi/rtot**2)
                        gaddendum10=-gaddendum9
                        gaddendum11=p2*(rtot)*(rprod*beta*(psi_y+psi_z)/mod_psi/rtot**2)
                        gaddendum12=-gaddendum11
                        gaddendum13=p2*(rtot)*(rprod*beta*(psi_y-psi_z)/mod_psi/rtot**2)
                        gaddendum14=-gaddendum13
                        gaddendum15=p2*(rtot)*(rprod*beta*(psi_x+psi_z)/mod_psi/rtot**2)
                        gaddendum16=-gaddendum15
                        gaddendum17=p2*(rtot)*(rprod*beta*(-psi_x+psi_z)/mod_psi/rtot**2)
                        gaddendum18=-gaddendum17
                    endif
                    !!************************************************************************!
                    press_Excess=0.0_db
                    psid1=0.0_db
                    psid2=0.0_db
                    psid3=0.0_db
                    psid4=0.0_db
                    ncontact=0
                    if(psi(i,j,k).lt.-0.9_db)then
                        
                        if (psi(i+3,j,k).gt.-0.85 .and. psi(i+3,j,k).lt.0.0_db .and. psi(i-3,j,k).gt.-0.85 .and. psi(i-3,j,k).lt.0.0_db)then !&
                            psid1=psi(i+2,j,k) + psi(i-2,j,k)
                            ncontact=ncontact+1
                        else
                            psid1=-2.0_db
                        endif
                        if (psi(i,j+3,k).gt.-0.85 .and. psi(i,j+3,k).lt.0.0_db .and. psi(i,j-3,k).gt.-0.85 .and. psi(i,j-3,k).lt.0.0_db)then
                            psid2=psi(i,j+2,k) + psi(i,j-2,k)
                            ncontact=ncontact+1
                        else
                            psid2=-2.0_db
                        endif
                        if (psi(i,j,k+3).gt.-0.85 .and. psi(i,j,k+3).lt.0.0_db .and. psi(i,j,k-3).gt.-0.85 .and. psi(i,j,k-3).lt.0.0_db)then
                            psid3=psi(i,j,k+2) + psi(i,j,k-2)
                            ncontact=ncontact+1
                        else
                            psid3=-2.0_db
                        endif
                        if (psi(i+2,j+2,k).gt.-0.85 .and. psi(i+2,j+2,k).lt.0.0_db .and. psi(i-2,j-2,k).gt.-0.85 .and. psi(i-2,j-2,k).lt.0.0_db)then
                            psid4=psi(i+1,j+1,k) + psi(i-1,j-1,k)
                            ncontact=ncontact+1
                        else
                            psid4=-2.0_db
                        endif
                        if (psi(i+2,j,k+2).gt.-0.85 .and. psi(i+2,j,k+2).lt.0.0_db .and. psi(i-2,j,k-2).gt.-0.85 .and. psi(i-2,j,k-2).lt.0.0_db)then
                            psid5=psi(i+1,j,k+1) + psi(i-1,j,k-1)
                            ncontact=ncontact+1
                        else
                            psid5=-2.0_db
                        endif
                        if (psi(i,j+2,k+2).gt.-0.85 .and. psi(i,j+2,k+2).lt.0.0_db .and. psi(i,j-2,k-2).gt.-0.85 .and. psi(i,j-2,k-2).lt.0.0_db)then
                            psid6=psi(i,j+1,k+1) + psi(i,j-1,k-1)
                            ncontact=ncontact+1
                        else
                            psid6=-2.0_db
                        endif

                        if (psi(i-2,j+2,k).gt.-0.85 .and. psi(i-2,j+2,k).lt.0.0_db .and. psi(i+2,j-2,k).gt.-0.85 .and. psi(i+2,j-2,k).lt.0.0_db)then
                            psid7=psi(i+1,j-1,k) + psi(i-1,j+1,k)
                            ncontact=ncontact+1
                        else
                            psid7=-2.0_db
                        endif
                        if (psi(i-2,j,k+2).gt.-0.85 .and. psi(i-2,j,k+2).lt.0.0_db .and. psi(i+2,j,k-2).gt.-0.85 .and. psi(i+2,j,k-2).lt.0.0_db)then
                            psid8=psi(i+1,j,k-1) + psi(i-1,j,k+1)
                            ncontact=ncontact+1
                        else
                            psid8=-2.0_db
                        endif
                        if (psi(i,j-2,k+2).gt.-0.85 .and. psi(i,j-2,k+2).lt.0.0_db .and. psi(i,j+2,k-2).gt.-0.85 .and. psi(i,j+2,k-2).lt.0.0_db)then
                            psid9=psi(i,j+1,k-1) + psi(i,j-1,k+1)
                            ncontact=ncontact+1
                        else
                            psid9=-2.0_db
                        endif
                        if(ncontact.gt.0)then
                            press_excess=max_press_excess*(18.0_db+psid1+psid2+psid3+psid4+psid5+psid6+psid7+psid8+psid9)/real(ncontact,4)
                        endif
                    endif
                    !*************************************************************************!
                    !0
                    feq=p0*(rtot+press_excess-uu)
                    fpc=feq + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k))) + addendum0
                    f0(i,j,k)=fpc*rhoA(i,j,k)/rtot
                    g0(i,j,k)=fpc*rhoB(i,j,k)/rtot
                    
                    !1
                    udotc=u(i,j,k)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p1*(rtot+press_excess+(temp + udotc))
                    fneq1=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                    fpc=feq+fneq1+fx*p1dcssq + addendum1
                    f1(i+1,j,k)=fpc*rhoA(i,j,k)/rtot + gaddendum1
                    g1(i+1,j,k)=fpc*rhob(i,j,k)/rtot - gaddendum1
                    
                    !2
                    feq=p1*(rtot+press_excess+(temp - udotc))
                    fpc=feq + fneq1 - fx*p1dcssq + addendum2
                    f2(i-1,j,k)=fpc*rhoA(i,j,k)/rtot + gaddendum2
                    g2(i-1,j,k)=fpc*rhob(i,j,k)/rtot - gaddendum2
                    !3
                    udotc=v(i,j,k)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p1*(rtot+press_excess+(temp + udotc))
                    fneq3=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
                    fpc=feq+fneq3 + fy*p1dcssq+addendum3
                    f3(i,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum3
                    g3(i,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum3
                    
                    !4
                    feq=p1*(rtot+press_excess+(temp - udotc))
                    fpc=feq+fneq3- fy*p1dcssq+addendum4
                    f4(i,j-1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum4
                    g4(i,j-1,k)=fpc*rhob(i,j,k)/rtot - gaddendum4

                    !7
                    udotc=(u(i,j,k)+v(i,j,k))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p2*(rtot+press_excess+(temp + udotc))
                    fneq7=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k))
                    fpc=feq+fneq7+ (fx+fy)*p2dcssq + addendum7
                    f7(i+1,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum7
                    g7(i+1,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum7

                    !8
                    feq=p2*(rtot+press_excess+(temp - udotc))
                    fpc=feq + fneq7- (fx+fy)*p2dcssq + addendum8
                    f8(i-1,j-1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum8
                    g8(i-1,j-1,k)=fpc*rhob(i,j,k)/rtot - gaddendum8 

                    !10
                    udotc=(-u(i,j,k)+v(i,j,k))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p2*(rtot+press_excess+(temp + udotc))
                    fneq10=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k))
                    fpc=feq + fneq10 + (fy-fx)*p2dcssq + addendum10
                    f10(i-1,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum10
                    g10(i-1,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum10

                    !9
                    feq=p2*(rtot+press_excess+(temp - udotc))
                    fpc=feq + fneq10 + (fx-fy)*p2dcssq + addendum9
                    f9(i+1,j-1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum9
                    g9(i+1,j-1,k)=fpc*rhob(i,j,k)/rtot - gaddendum9

                    !5
                    udotc=w(i,j,k)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p1*(rtot+press_excess+(temp + udotc))
                    fneq5=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
                    fpc=feq+fneq5 + fz*p1dcssq + addendum5
                    f5(i,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum5
                    g5(i,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum5

                    !6
                    feq=p1*(rtot+press_excess+(temp - udotc))
                    fpc=feq+fneq5 - fz*p1dcssq + addendum6
                    f6(i,j,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum6
                    g6(i,j,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum6

                    !15
                    udotc=(u(i,j,k)+w(i,j,k))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p2*(rtot+press_excess+(temp + udotc))
                    fneq15=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k))
                    fpc=feq+fneq15 + (fx+fz)*p2dcssq + addendum15
                    f15(i+1,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum15
                    g15(i+1,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum15
                    
                    !16
                    feq=p2*(rtot+press_excess+(temp - udotc))
                    fpc=feq+fneq15 - (fx+fz)*p2dcssq + addendum16
                    f16(i-1,j,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum16
                    g16(i-1,j,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum16

                    !17
                    udotc=(-u(i,j,k)+w(i,j,k))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p2*(rtot+press_excess+(temp + udotc))
                    fneq17=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k))
                    fpc=feq+fneq17 +(fz-fx)*p2dcssq +addendum17
                    f17(i-1,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum17
                    g17(i-1,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum17
                    
                    !18
                    feq=p2*(rtot+press_excess+(temp - udotc))
                    fpc=feq+fneq17 + (fx-fz)*p2dcssq + addendum18
                    f18(i+1,j,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum18
                    g18(i+1,j,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum18

                    !11
                    udotc=(v(i,j,k)+w(i,j,k))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p2*(rtot+press_excess+(temp + udotc))
                    fneq11=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k))
                    fpc=feq+fneq11+(fy+fz)*p2dcssq + addendum11
                    f11(i,j+1,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum11
                    g11(i,j+1,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum11
                    
                    !12
                    feq=p2*(rtot+press_excess+(temp - udotc))
                    fpc=feq+fneq11 - (fy+fz)*p2dcssq + addendum12
                    f12(i,j-1,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum12
                    g12(i,j-1,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum12
                    
                    !13
                    udotc=(v(i,j,k)-w(i,j,k))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    feq=p2*(rtot+press_excess+(temp + udotc))
                    fneq13=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k))
                    fpc=feq+fneq13 + (fy-fz)*p2dcssq + addendum13
                    f13(i,j+1,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum13
                    g13(i,j+1,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum13
                    
                    !14
                    feq=p2*(rtot+press_excess+(temp - udotc))
                    fpc=feq+fneq13 + (fz-fy)*p2dcssq + addendum14
                    f14(i,j-1,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum14
                    g14(i,j-1,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum14
                enddo
            enddo
        enddo
        !******************************************call bcs************************
        !periodic along y
        !x=1     
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
        !****************************************writeonfile***************************************************!
            if(mod(step,stamp).eq.0)then
            !$acc update self(psi(1:nx,ny/2,1:nz),rhoB(1:nx,ny/2,1:nz),rhoA(1:nx,ny/2,1:nz),u(1:nx,ny/2,1:nz),w(1:nx,ny/2,1:nz))
            iframe=iframe+1
            !xz
            open(101, file = 'psi_xz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(102, file = 'rhoA_xz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(103, file = 'rhoB_xz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(104, file = 'u_xz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(105, file = 'w_xz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            do i=1,nx
                do j=1,nz
                    write(101,*) psi(i,ny/2,j)  
                enddo
            enddo
            close(101)
            do i=1,nx
                do j=1,nz
                    write(102,*) rhoA(i,ny/2,j) 
                enddo
            enddo
            close(102)
            do i=1,nx
                do j=1,nz
                    write(103,*) rhoB(i,ny/2,j)   
                enddo
            enddo
            close(103)
            do i=1,nx
                do j=1,nz
                    write(104,*) u(i,ny/2,j)   
                enddo
            enddo
            close(104)
            do i=1,nx
                do j=1,nz
                    write(105,*) w(i,ny/2,j)   
                enddo
            enddo
            close(105)
            ! xy
            !$acc update self(psi(1:nx,1:ny,nz/2),rhoB(1:nx,1:ny,nz/2),rhoA(1:nx,1:ny,nz/2),u(1:nx,1:ny,nz/2),v(1:nx,1:ny,nz/2))
            open(106, file = 'psi_xy'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(107, file = 'rhoA_xy'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(108, file = 'rhoB_xy'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(109, file = 'u_xy'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(110, file = 'v_xz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            do i=1,nx
                do j=1,ny
                    write(106,*) psi(i,j,nz/2)  
                enddo
            enddo
            close(106)
            do i=1,nx
                do j=1,ny
                    write(107,*) rhoA(i,j,nz/2)  
                enddo
            enddo
            close(107)
            do i=1,nx
                do j=1,ny
                    write(108,*) rhoB(i,j,nz/2)     
                enddo
            enddo
            close(108)
            do i=1,nx
                do j=1,ny
                    write(109,*) u(i,j,nz/2)     
                enddo
            enddo
            close(109)
            do i=1,nx
                do j=1,ny
                    write(110,*) v(i,j,nz/2)     
                enddo
            enddo
            close(110)
            ! yz
            !$acc update self(psi(nx/2,1:ny,1:nz),rhoB(nx/2,1:ny,1:nz),rhoA(nx/2,1:ny,1:nz),v(nx/2,1:ny,1:nz),w(nx/2,1:ny,1:nz))
            open(111, file = 'psi_yz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(112, file = 'rhoA_yz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(113, file = 'rhoB_yz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(114, file = 'v_yz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(115, file = 'w_xz'//write_fmtnumb(iframe)//'.out', status = 'replace')
            do i=1,nx
                do j=1,ny
                    write(111,*) psi(i,j,nz/2)  
                enddo
            enddo
            close(111)
            do i=1,nx
                do j=1,ny
                    write(112,*) rhoA(i,j,nz/2)  
                enddo
            enddo
            close(112)
            do i=1,nx
                do j=1,ny
                    write(113,*) rhoB(i,j,nz/2)     
                enddo
            enddo
            close(113)
            do i=1,nx
                do j=1,ny
                    write(114,*) v(i,j,nz/2)     
                enddo
            enddo
            close(114)
            do i=1,nx
                do j=1,ny
                    write(115,*) w(i,j,nz/2)     
                enddo
            enddo
            close(115)
            write(6,*) "files updated at t=", step
            endif
        
        
    enddo 
    call cpu_time(ts2)
    !$acc end data
   contains 

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