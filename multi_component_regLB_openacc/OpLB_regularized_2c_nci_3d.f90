program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    !************************************block of vars*******************************************!
        integer, parameter :: db=4 !kind(1.0)
        integer(kind=8) :: i,j,k,k_init,k_end,j_init,j_end,i_init,i_end
        integer(kind=8) :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus
        integer,save :: iframe=0
        
        real(kind=db),parameter :: pi_greek=3.14159265359793234626433
        
        real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
        real(kind=4)  :: p1cg,p2cg,p3cg
        real(kind=db) :: visc_LB,uu,udotc,omega,feq,geq,fpc,gpc,hpc
        real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp
        real(kind=db) :: radius

        real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq17
        real(kind=db) :: fneq9,fneq10,fneq11,fneq12,fneq13,fneq14,fneq15,fneq16,fneq18
        real(kind=db) :: ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft17
        real(kind=db) :: ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft18
        real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
        real(kind=db) :: pi2cssq1,pi2cssq2,pi2cssq0
        real(kind=db) :: addendum0,addendum1,addendum3,addendum5,addendum7
        real(kind=db) :: addendum10,addendum11,addendum13,addendum15,addendum17
        real(kind=db) :: gaddendum1,gaddendum2,gaddendum3,gaddendum4,gaddendum5,gaddendum6,gaddendum7,gaddendum8
        real(kind=db) :: gaddendum9,gaddendum10,gaddendum11,gaddendum12,gaddendum13,gaddendum14,gaddendum15,gaddendum16,gaddendum17,gaddendum18
        real(kind=db) :: psi_x,psi_y,psi_z,mod_psi,mod_psi_sq,st_coeff,b0,b1,b2,beta,sigma
        real(kind=db) :: one_ov_nu2,one_ov_nu1,nu_avg,rtot,rprod
        real(kind=db) :: max_press_excess,rr,ushifted,vshifted,wshifted,norm_x,norm_y,norm_z

        integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
        real(kind=db), allocatable, dimension(:,:,:) :: psi,rhoA,rhoB,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
        real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
        real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18
        real(kind=db), allocatable, dimension(:,:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8,g9
        real(kind=db), allocatable, dimension(:,:,:) :: g10,g11,g12,g13,g14,g15,g16,g17,g18
        real(kind=db), allocatable,dimension(:,:,:) :: random_field
        integer(kind=4), allocatable, dimension(:,:,:) :: nci_loc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!for custom geometry*********************************!
        integer:: ddrop,lc,jjd,jju
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!for custom geometry*********************************!
       
   
    
    
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
        nx=416
        ny=416
        nz=416
        nsteps=100
        stamp=20000
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
        allocate(random_field(25,25,25), nci_loc(1:nx,1:ny,1:nz))
    
    !***************************************lattice/vars*************************************!
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

    !***********************************define/read geometry if any**************************
            
            ! Lc=250+nz/2;
            ! do i=Lc,nz
            !     jjd=nint((i-Lc+1)*sind(30.0_db));
                
            !     if(jjd<ny/2-4)then
            !         isfluid(:,jjd:ny/2,i)=1; 
            !         jju=ny-jjd;
            !         isfluid(:,ny/2:jju,i)=1;
            !     endif
            ! enddo
            ! Ddrop=40;
            ! isfluid(:,2:ny-1,nz/2:Lc)=1;
            ! isfluid(:,ny/2-(Ddrop/(2)):ny/2+(Ddrop/(2)),nz/2:nz)=1;
            ! isfluid(:,:,1:nz/2)=isfluid(:,:,nz:nz/2+1:-1);
            ! isfluid(1,:,:)=0 !left
            ! isfluid(nx,:,:)=0 !right
            ! isfluid(:,1,:)=0 !front 
            ! isfluid(:,ny,:)=0 !rear
            ! isfluid(:,:,1)=0 !bottom
            ! isfluid(:,:,nz)=0 !top
            ! do k=1,nz
            !     do j=1,ny
            !         if(isfluid(nx/2,j,k).eq.3)then
            !             if(isfluid(nx/2,j+1,k).eq.1 .or. isfluid(nx/2,j-1,k).eq.1 .or. isfluid(nx/2,j,k+1).eq.1 .or. isfluid(nx/2,j,k-1).eq.1 &
            !                 .or. isfluid(nx/2,j+1,k+1).eq.1 .or. isfluid(nx/2,j+1,k-1).eq.1 .or. isfluid(nx/2,j-1,k+1).eq.1 .or. isfluid(nx/2,j-1,k-1).eq.1)then
            !                 isfluid(2:nx-1,j,k)=0
            !             endif
            !         endif
            !     enddo
            ! enddo
            ! open(231, file = 'isfluid.out', status = 'replace')
            !     do j=1,ny
            !         do k=1,nz
            !             write(231,*) isfluid(nx/2,j,k)  
            !         enddo
            !     enddo
            ! close(231)
        
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
        beta=0.95_db
        sigma=0.03_db
        st_coeff=(9.0_db/4.0_db)*sigma*omega
        b0=-2.0_db/9.0_db
        b1=1.0_db/54.0_db
        b2=1.0_db/27.0_db
        p1cg=(2.0_db/3.0_db)**2 * (1.0_db/6.0_db)
        p2cg=(1.0_db/6.0_db)**2 * (2.0_db/3.0_db)
        p3cg=(1.0_db/6.0_db)**3
        nci_loc=0
        max_press_excess=0.004
    !********************************initialization of macrovars ************************    
        u=0.0_db
        v=0.0_db
        w=0.0_db
        rhoA(1:nx,1:ny,1:nz)=0.0_db  !total density
        rhoB(1:nx,1:ny,1:nz)=0.0_db  !total density
        psi=-1.0_db
        radius=20
        !*****************************************Impacting droplets***************************!
            do i=(nx/2)-radius,(nx/2)+radius
                do j=ny/2-radius,ny/2+radius
                    do k=nz/2-radius,nz/2+radius
                        if ((i-(nx/2))**2+(j-ny/2)**2+(k-(nz/2))**2<=radius**2)then
                            psi(i,j,k)=1.0_db
                            u(i,j,k)=0.0_db
                        endif
                    enddo
                enddo
            enddo
            ! do i=(nx/2+radius+5)-radius,(nx/2+radius+5)+radius
            !     do j=ny/2-radius,ny/2+radius
            !         do k=nz/2-radius,nz/2+radius
            !             if ((i-(nx/2+radius+5))**2+(j-ny/2)**2+(k-(nz/2))**2<=radius**2)then
            !                 psi(i,j,k)=1.0_db
            !                 u(i,j,k)=-0.035/2.0_db
            !             endif
            !         enddo
            !     enddo
            ! enddo
        
        !*****************************************Spinodal decomposition***********************!
            ! do k=1,25
            !     do j=1,25
            !         do i=1,25
            !             call random_number(rr)
            !             random_field(i,j,k)=rr
            !         enddo
            !     enddo
            ! enddo
            ! !
            ! do k=1,25
            !     k_init=k_end+1
            !     k_end=k_end+10
            !     do j=1,25
            !         j_init=j_end+1
            !         j_end=j_end+10
            !         do i=1,25
            !             i_init=i_end+1
            !             i_end=i_end+10
            !             if(random_field(i,j,k).gt.0.5) psi(i_init:i_end,j_init:j_end,k_init:k_end)=1.0_db
            !         enddo
            !         i_init=0
            !         i_end=0
            !     enddo
            !     j_init=0
            !     j_end=0
            ! enddo
        !************************read initial conditions for macrovars*************************!
            ! open(231, file = 'psi.txt', status = 'old',action='read')
            ! do j=1,ny
            !     do k=1,nz
            !         read(231,*) psi(15,j,k)
            !     enddo
            ! enddo
            ! close(231)
            ! do i=4,nx-3
            !     psi(i,:,:)=psi(15,:,:)
            ! enddo
            
            ! open(231, file = 'psi.out', status = 'replace')
            ! do i=1,nx
            !     do j=1,ny
            !         do k=1,nz
            !         write(231,*) psi(i,j,k)  
            !         enddo
            !     enddo
            ! enddo
            ! close(231)
        !************************************single cylindrical droplets*******************!
            ! do i=5,nx-4
            !     do j=ny/2-radius,ny/2+radius
            !         do k=nz-15-radius-radius,nz-15-radius+radius
            !             if ((j-ny/2)**2+(k-(nz-15-radius))**2<=radius**2)then
            !                 psi(i,j,k)=1.0_db
            !             endif
            !         enddo
            !     enddo
            ! enddo
        !*****************************************dense emulsion in channel********************!
        !*****************************************turbulent emulsion********************!
        rhoB=0.5*(1.0_db-psi(1:nx,1:ny,1:nz))
        rhoA=1.0_db-rhoB
    !*************************************set distros************************!
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
        f0(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p0*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f1(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0+u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f2(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f3(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f4(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f5(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f6(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f7(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0+u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f8(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f9(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0+u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f10(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f11(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f12(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f13(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f14(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f15(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0 + u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2 - u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f16(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0 - u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2 - u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f17(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0 - u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2 - u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f18(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0 + u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2 - u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        !
        g0(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p0!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g1(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0+u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g2(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g3(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g4(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g5(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g6(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g7(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0+u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g8(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g9(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0+u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g10(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g11(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g12(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g13(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g14(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g15(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0 +u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2 - u(:,:,:)*u(:,:,:)*0.5/cssq)
        g16(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0 -u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2 - u(:,:,:)*u(:,:,:)*0.5/cssq)
        g17(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0 -u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2 - u(:,:,:)*u(:,:,:)*0.5/cssq)
        g18(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0 +u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2 - u(:,:,:)*u(:,:,:)*0.5/cssq)
    
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
        write(6,*) 'max press excess',max_press_excess
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
             !$acc& g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,nci_loc)
    
    
    !**********************************************************************!
    call cpu_time(ts1)
    !*************************************main loop*************************!
    do step=1,nsteps 
        !$acc kernels 
        !$acc loop collapse (3) 
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
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
                        nci_loc(i,j,k)=0
                        nci_loc(i,j,k)=0
                        rtot=rhoA(i,j,k)+rhoB(i,j,k)

                        psi(i,j,k)= (rhoA(i,j,k)-rhoB(i,j,k))/(rhoA(i,j,k)+rhoB(i,j,k))

                        u(i,j,k) = (ft1+ft7+ft9+ft15+ft18)-(ft2+ft8+ft10+ft16+ft17) 
                        
                        v(i,j,k) = (ft3+ft7+ft10+ft11+ft13)-(ft4+ft8+ft9+ft12+ft14)

                        w(i,j,k) = (ft5+ft11+ft14+ft15+ft17)-(ft6+ft12+ft13+ft16+ft18)
                        
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !1-2
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
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
                        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
                        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
                        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
                        pxx(i,j,k)=fneq1+fneq2+fneq7+fneq8+fneq9+fneq10+fneq15+fneq16+fneq17+fneq18
                        pyy(i,j,k)=fneq3+fneq4+fneq7+fneq8+fneq9+fneq10+fneq11+fneq12+fneq13+fneq14
                        pzz(i,j,k)=fneq5+fneq6+fneq11+fneq12+fneq13+fneq14+fneq15+fneq16+fneq17+fneq18
                        pxy(i,j,k)=fneq7+fneq8-fneq9-fneq10
                        pxz(i,j,k)=fneq15+fneq16-fneq17-fneq18
                        pyz(i,j,k)=fneq11+fneq12-fneq13-fneq14
                    endif
                enddo
            enddo
        enddo
        !$acc loop collapse (3) 
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                    !!*******************************near contact***************************!
                        if(isfluid(i,j,k).eq.1)then
                            
                            if(psi(i,j,k).lt.-0.9_db .and. i.gt.3 .and. i.lt.nx-2 .and. j.gt.3 .and. j.lt.ny-2 .and. k.gt.3 .and. k.lt.nz-2)then
                                
                                if (psi(i+3,j,k).gt.-0.85 .and. psi(i-3,j,k).gt.-0.85)then
                                    nci_loc(i+3,j,k)=1
                                    nci_loc(i-3,j,k)=1
                                endif
                                if (psi(i,j+3,k).gt.-0.85 .and. psi(i,j-3,k).gt.-0.85)then
                                    nci_loc(i,j+3,k)=1
                                    nci_loc(i,j-3,k)=1
                                endif
                                if (psi(i,j,k+3).gt.-0.85 .and. psi(i,j,k-3).gt.-0.85)then!(psi(i,j,k+3).gt.-0.85 .and. psi(i,j,k+3).lt.0.0_db .and. psi(i,j,k-3).gt.-0.85 .and. psi(i,j,k-3).lt.0.0_db)then
                                   nci_loc(i,j,k+3)=1
                                   nci_loc(i,j,k-3)=1
                                endif
                                if (psi(i+2,j+2,k).gt.-0.85 .and. psi(i-2,j-2,k).gt.-0.85)then!(psi(i+3,j+3,k).gt.-0.85 .and. psi(i+3,j+3,k).lt.0.0_db .and. psi(i-3,j-3,k).gt.-0.85 .and. psi(i-3,j-3,k).lt.0.0_db)then
                                    nci_loc(i+2,j+2,k)=1
                                    nci_loc(i-2,j-2,k)=1
                                endif
                                if (psi(i+2,j,k+2).gt.-0.85.and. psi(i-2,j,k-2).gt.-0.85)then!(psi(i+3,j,k+3).gt.-0.85 .and. psi(i+3,j,k+3).lt.0.0_db .and. psi(i-3,j,k-3).gt.-0.85 .and. psi(i-3,j,k-3).lt.0.0_db)then
                                    nci_loc(i+2,j,k+2)=1
                                    nci_loc(i-2,j,k-2)=1
                                endif
                                if (psi(i,j+2,k+2).gt.-0.85 .and. psi(i,j-2,k-2).gt.-0.85)then!(psi(i,j+3,k+3).gt.-0.85 .and. psi(i,j+3,k+3).lt.0.0_db .and. psi(i,j-3,k-3).gt.-0.85 .and. psi(i,j-3,k-3).lt.0.0_db)then
                                    nci_loc(i,j+2,k+2)=1
                                    nci_loc(i,j-2,k-2)=1
                                endif

                                if (psi(i-2,j+2,k).gt.-0.85.and. psi(i+2,j-2,k).gt.-0.85)then!(psi(i-3,j+3,k).gt.-0.85 .and. psi(i-3,j+3,k).lt.0.0_db .and. psi(i+3,j-3,k).gt.-0.85 .and. psi(i+3,j-3,k).lt.0.0_db)then
                                    nci_loc(i-2,j+2,k)=1
                                    nci_loc(i+2,j-2,k)=1
                                endif
                                if (psi(i-2,j,k+2).gt.-0.85.and. psi(i+2,j,k-2).gt.-0.85)then!(psi(i-3,j,k+3).gt.-0.85 .and. psi(i-3,j,k+3).lt.0.0_db .and. psi(i+3,j,k-3).gt.-0.85 .and. psi(i+3,j,k-3).lt.0.0_db)then
                                    nci_loc(i+2,j,k-2)=1
                                    nci_loc(i-2,j,k+2)=1
                                endif
                                if (psi(i,j-2,k+2).gt.-0.85.and. psi(i,j+2,k-2).gt.-0.85)then!(psi(i,j-3,k+3).gt.-0.85 .and. psi(i,j-3,k+3).lt.0.0_db .and. psi(i,j+3,k-3).gt.-0.85 .and. psi(i,j+3,k-3).lt.0.0_db)then
                                    nci_loc(i,j-2,k+2)=1
                                    nci_loc(i,j+2,k-2)=1
                                endif
                                
                            endif
                        endif
                    enddo
                enddo
            enddo
        !$acc loop collapse (3) 
        do k=1,nz
            do j=1,ny
                do i=1,nx  
                    if(isfluid(i,j,k).eq.1)then 
                    !!*******************************chromodynamics***************************!
                        psi_x=(1.0_db/cssq)*(p1cg*(psi(i+1,j,k)-psi(i-1,j,k)) + &
                                            p2cg*(psi(i+1,j+1,k)+psi(i+1,j-1,k)+psi(i+1,j,k+1)+psi(i+1,j,k-1)) &
                                        -p2cg*(psi(i-1,j+1,k)+psi(i-1,j-1,k)+psi(i-1,j,k+1)+psi(i-1,j,k-1)) &
                                        + p3cg*(psi(i+1,j+1,k+1)+psi(i+1,j+1,k-1)+psi(i+1,j-1,k+1)+psi(i+1,j-1,k-1)) &
                                        -p3cg*(psi(i-1,j+1,k+1)+psi(i-1,j+1,k-1)+psi(i-1,j-1,k+1)+psi(i-1,j-1,k-1))) 


                        psi_y=(1.0_db/cssq)*(p1cg*(psi(i,j+1,k)-psi(i,j-1,k)) &
                                            +p2cg*(psi(i+1,j+1,k)+psi(i-1,j+1,k)+psi(i,j+1,k+1)+psi(i,j+1,k-1)) &
                                        -p2cg*(psi(i+1,j-1,k)+psi(i-1,j-1,k)+psi(i,j-1,k+1)+psi(i,j-1,k-1)) &
                                        +p3cg*(psi(i+1,j+1,k+1)+psi(i-1,j+1,k+1)+psi(i+1,j+1,k-1)+psi(i-1,j+1,k-1)) &
                                        -p3cg*(psi(i+1,j-1,k+1)+psi(i-1,j-1,k+1)+psi(i+1,j-1,k-1)+psi(i-1,j-1,k-1)))


                        psi_z=(1.0_db/cssq)*(p1cg*(psi(i,j,k+1)-psi(i,j,k-1)) + &
                                            p2cg*(psi(i+1,j,k+1)+psi(i-1,j,k+1)+ psi(i,j+1,k+1)+psi(i,j-1,k+1)) &
                                        -p2cg*(psi(i+1,j,k-1)+psi(i-1,j,k-1)+psi(i,j+1,k-1)+psi(i,j-1,k-1)) &
                                        +p3cg*(psi(i+1,j+1,k+1)+psi(i-1,j+1,k+1)+psi(i+1,j-1,k+1)+psi(i-1,j-1,k+1)) &
                                        -p3cg*(psi(i+1,j+1,k-1)+psi(i-1,j+1,k-1)+psi(i+1,j-1,k-1)+psi(i-1,j-1,k-1)))
                        
                        mod_psi=sqrt(psi_x**2+psi_y**2+psi_z**2)

                        norm_x=0.0_db
                        norm_y=0.0_db
                        norm_z=0.0_db

                        mod_psi_sq=psi_x**2 + psi_y**2 +psi_z**2 
                        
                        rtot=0.0_db
                        
                        rtot=rhoA(i,j,k)+rhoB(i,j,k)
                        
                        rprod=rhoA(i,j,k)*rhoB(i,j,k)
                        
                        nu_avg=1.0_db/(rhoA(i,j,k)*one_ov_nu1/rtot + rhoB(i,j,k)*one_ov_nu2/rtot)
                        
                        omega=2.0_db/(6.0_db*nu_avg + 1.0_db)
                        
                        st_coeff=(9.0_db/4.0_db)*sigma*omega
                        addendum0=0.0_db
                        addendum1=0.0_db
                        addendum3=0.0_db
                        addendum5=0.0_db
                        addendum7=0.0_db
                        addendum10=0.0_db
                        addendum11=0.0_db
                        addendum13=0.0_db
                        addendum15=0.0_db
                        addendum17=0.0_db
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
                            addendum3=st_coeff*mod_psi*(p1*psi_y**2/mod_psi_sq - b1)
                            addendum5=st_coeff*mod_psi*(p1*psi_z**2/mod_psi_sq - b1)
                            !
                            addendum7=st_coeff*mod_psi*(p2*(psi_x+psi_y)**2/mod_psi_sq - b2)
                            addendum10=st_coeff*mod_psi*(p2*(-psi_x+psi_y)**2/mod_psi_sq - b2)
                            addendum11=st_coeff*mod_psi*(p2*(psi_y+psi_z)**2/mod_psi_sq - b2)
                            addendum13=st_coeff*mod_psi*(p2*(psi_y-psi_z)**2/mod_psi_sq - b2)
                            addendum15=st_coeff*mod_psi*(p2*(psi_x+psi_z)**2/mod_psi_sq - b2)
                            addendum17=st_coeff*mod_psi*(p2*(-psi_x+psi_z)**2/mod_psi_sq - b2)
                            !
                            gaddendum1=p1*(rtot)*(rprod*beta*psi_x/mod_psi/rtot**2)
                            gaddendum2=-gaddendum1 
                            gaddendum3=p1*(rtot)*(rprod*beta*psi_y/mod_psi/rtot**2)
                            gaddendum4=-gaddendum3 
                            gaddendum5=p1*(rtot)*(rprod*beta*psi_z/mod_psi/rtot**2)
                            gaddendum6=-gaddendum5 
                            !
                            gaddendum7=p2*(rtot)*(rprod*beta*(psi_x/mod_psi+psi_y/mod_psi)/rtot**2)
                            gaddendum8=-gaddendum7
                            gaddendum9=p2*(rtot)*(rprod*beta*(psi_x/mod_psi-psi_y/mod_psi)/rtot**2)
                            gaddendum10=-gaddendum9 
                            gaddendum11=p2*(rtot)*(rprod*beta*(psi_y/mod_psi+psi_z/mod_psi)/rtot**2)
                            gaddendum12=-gaddendum11 
                            gaddendum13=p2*(rtot)*(rprod*beta*(psi_y/mod_psi-psi_z/mod_psi)/rtot**2)
                            gaddendum14=-gaddendum13 
                            gaddendum15=p2*(rtot)*(rprod*beta*(psi_x/mod_psi+psi_z/mod_psi)/rtot**2)
                            gaddendum16=-gaddendum15 
                            gaddendum17=p2*(rtot)*(rprod*beta*(-psi_x/mod_psi+psi_z/mod_psi)/rtot**2)
                            gaddendum18=-gaddendum17
                            
                            norm_x=psi_x/mod_psi
                            norm_y=psi_y/mod_psi
                            norm_z=psi_z/mod_psi
                        endif
                    !******************************collision+stream*****************************!
                        !0
                        ushifted=u(i,j,k) + fx + float(nci_loc(i,j,k))*(norm_x)*max_press_excess*abs(rhoB(i,j,k))
                        vshifted=v(i,j,k) + fy + float(nci_loc(i,j,k))*(norm_y)*max_press_excess*abs(rhoB(i,j,k))
                        wshifted=w(i,j,k) + fz + float(nci_loc(i,j,k))*(norm_z)*max_press_excess*abs(rhoB(i,j,k))

                        uu=0.5_db*(ushifted*ushifted + vshifted*vshifted + wshifted*wshifted)/cssq 

                        feq=p0*(rtot-uu)
                        fpc=feq + addendum0 + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k))) 
                        f0(i,j,k)=fpc*rhoA(i,j,k)/rtot
                        g0(i,j,k)=fpc*rhoB(i,j,k)/rtot
                        !1
                        udotc=ushifted/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rtot+(temp + udotc))
                        fneq1=(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                        fpc=feq + addendum1+fneq1+fx*p1dcssq 
                        f1(i+1,j,k)=fpc*rhoA(i,j,k)/rtot + gaddendum1
                        g1(i+1,j,k)=fpc*rhob(i,j,k)/rtot - gaddendum1
                        
                        !2
                        feq=p1*(rtot+(temp - udotc))
                        fpc=feq + addendum1 + fneq1 - fx*p1dcssq 
                        f2(i-1,j,k)=fpc*rhoA(i,j,k)/rtot + gaddendum2
                        g2(i-1,j,k)=fpc*rhob(i,j,k)/rtot - gaddendum2
                        
                        !3
                        udotc=vshifted/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rtot+(temp + udotc))
                        fneq3=(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
                        fpc=feq +addendum3+fneq3 + fy*p1dcssq
                        f3(i,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum3
                        g3(i,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum3
                        
                        !4
                        feq=p1*(rtot+(temp - udotc))
                        fpc=feq +addendum3+fneq3- fy*p1dcssq
                        f4(i,j-1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum4
                        g4(i,j-1,k)=fpc*rhob(i,j,k)/rtot - gaddendum4

                        !7
                        udotc=(ushifted+vshifted)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rtot+(temp + udotc))
                        fneq7=(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*pxy(i,j,k))
                        fpc=feq + addendum7+fneq7+ (fx+fy)*p2dcssq 
                        f7(i+1,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum7
                        g7(i+1,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum7

                        !8
                        feq=p2*(rtot+(temp - udotc))
                        fpc=feq + addendum7+ fneq7- (fx+fy)*p2dcssq 
                        f8(i-1,j-1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum8
                        g8(i-1,j-1,k)=fpc*rhob(i,j,k)/rtot - gaddendum8 
                        
                        !10
                        udotc=(-ushifted+vshifted)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rtot+(temp + udotc))
                        fneq10=(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)-2.0_db*pxy(i,j,k))
                        fpc=feq + addendum10+ fneq10 + (fy-fx)*p2dcssq 
                        f10(i-1,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum10
                        g10(i-1,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum10

                        !9
                        feq=p2*(rtot+(temp - udotc))
                        fpc=feq + addendum10+ fneq10 + (fx-fy)*p2dcssq 
                        f9(i+1,j-1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum9
                        g9(i+1,j-1,k)=fpc*rhob(i,j,k)/rtot - gaddendum9

                        !5
                        udotc=wshifted/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rtot+(temp + udotc))
                        fneq5=(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
                        fpc=feq + addendum5+fneq5 + fz*p1dcssq 
                        f5(i,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum5
                        g5(i,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum5

                        !6
                        feq=p1*(rtot+(temp - udotc))
                        fpc=feq + addendum5 +fneq5 - fz*p1dcssq 
                        f6(i,j,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum6
                        g6(i,j,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum6

                        !15
                        udotc=(ushifted+wshifted)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rtot+(temp + udotc))
                        fneq15=(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*pxz(i,j,k))
                        fpc=feq + addendum15+fneq15 + (fx+fz)*p2dcssq 
                        f15(i+1,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum15
                        g15(i+1,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum15
                        
                        !16
                        feq=p2*(rtot+(temp - udotc))
                        fpc=feq + addendum15 +fneq15 - (fx+fz)*p2dcssq 
                        f16(i-1,j,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum16
                        g16(i-1,j,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum16

                        !17
                        udotc=(-ushifted+wshifted)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rtot+(temp + udotc))
                        fneq17=(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)-2.0_db*pxz(i,j,k))
                        fpc=feq +addendum17+fneq17 +(fz-fx)*p2dcssq 
                        f17(i-1,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum17
                        g17(i-1,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum17
                        
                        !18
                        feq=p2*(rtot+(temp - udotc))
                        fpc=feq + addendum17+fneq17 + (fx-fz)*p2dcssq 
                        f18(i+1,j,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum18
                        g18(i+1,j,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum18

                        !11
                        udotc=(vshifted+wshifted)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rtot+(temp + udotc))
                        fneq11=(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*pyz(i,j,k))
                        fpc=feq + addendum11+fneq11+(fy+fz)*p2dcssq 
                        f11(i,j+1,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum11
                        g11(i,j+1,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum11
                        
                        !12
                        feq=p2*(rtot+(temp - udotc))
                        fpc=feq + addendum11+fneq11 - (fy+fz)*p2dcssq 
                        f12(i,j-1,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum12
                        g12(i,j-1,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum12

                        !13
                        udotc=(vshifted-wshifted)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rtot+(temp + udotc))
                        fneq13=(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)-2.0_db*pyz(i,j,k))
                        fpc=feq + addendum13+fneq13 + (fy-fz)*p2dcssq 
                        f13(i,j+1,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum13
                        g13(i,j+1,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum13
                        
                        !14
                        feq=p2*(rtot+(temp - udotc))
                        fpc=feq + addendum13+fneq13 + (fz-fy)*p2dcssq 
                        f14(i,j-1,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum14
                        g14(i,j-1,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum14
                    endif
                enddo
            enddo
        enddo
        !********************************boundary conditions no slip everywhere********************************!
            !$acc loop independent 
            do k=1,nz
                !$acc loop independent 
                do j=1,ny
                    !$acc loop independent 
                    do i=1,nx
                        if(isfluid(i,j,k).eq.0)then
                            psi(i,j,k)=-1.0_db
                            f18(i+1,j,k-1)=f17(i,j,k) !gpc 
                            f17(i-1,j,k+1)=f18(i,j,k) !hpc

                            f16(i-1,j,k-1)=f15(i,j,k) !gpc 
                            f15(i+1,j,k+1)=f16(i,j,k) !hpc

                            f14(i,j-1,k+1)=f13(i,j,k)!gpc 
                            f13(i,j+1,k-1)=f14(i,j,k)!hpc
                            
                            f12(i,j-1,k-1)=f11(i,j,k)!gpc 
                            f11(i,j+1,k+1)=f12(i,j,k)!hpc

                            f10(i-1,j+1,k)=f9(i,j,k)!gpc 
                            f9(i+1,j-1,k)=f10(i,j,k)!hpc

                            f8(i-1,j-1,k)=f7(i,j,k)!gpc 
                            f7(i+1,j+1,k)=f8(i,j,k)!hpc

                            f6(i,j,k-1)=f5(i,j,k)!gpc 
                            f5(i,j,k+1)=f6(i,j,k)!hpc 


                            f4(i,j-1,k)=f3(i,j,k)!gpc 
                            f3(i,j+1,k)=f4(i,j,k)!hpc 

                            f2(i-1,j,k)=f1(i,j,k)!gpc 
                            f1(i+1,j,k)=f2(i,j,k)!hpc 
                            !************************!
                            g18(i+1,j,k-1)=g17(i,j,k) !gpc 
                            g17(i-1,j,k+1)=g18(i,j,k) !hpc

                            g16(i-1,j,k-1)=g15(i,j,k) !gpc 
                            g15(i+1,j,k+1)=g16(i,j,k) !hpc

                            g14(i,j-1,k+1)=g13(i,j,k)!gpc 
                            g13(i,j+1,k-1)=g14(i,j,k)!hpc
                            
                            g12(i,j-1,k-1)=g11(i,j,k)!gpc 
                            g11(i,j+1,k+1)=g12(i,j,k)!hpc

                            g10(i-1,j+1,k)=g9(i,j,k)!gpc 
                            g9(i+1,j-1,k)=g10(i,j,k)!hpc

                            g8(i-1,j-1,k)=g7(i,j,k)!gpc 
                            g7(i+1,j+1,k)=g8(i,j,k)!hpc

                            g6(i,j,k-1)=g5(i,j,k)!gpc 
                            g5(i,j,k+1)=g6(i,j,k)!hpc 


                            g4(i,j-1,k)=g3(i,j,k)!gpc 
                            g3(i,j+1,k)=g4(i,j,k)!hpc 

                            g2(i-1,j,k)=g1(i,j,k)!gpc 
                            g1(i+1,j,k)=g2(i,j,k)!hpc
                        endif
                    enddo
                enddo
            enddo
        !*********************************call bcs(other than no slip)************************
       
            ! !$acc loop independent 
            ! do j=1,ny
            !     !$acc loop independent 
            !     do i=1,nx
            !         psi(i,j,nz)=psi(i,j,2)
            !         psi(i,j,1)=psi(i,j,ny-1)
            !         !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
            !         !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
            !         !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
            !         f6(i,j,nz-1)=f6(i,j,1)
            !         f12(i,j,nz-1)=f12(i,j,1)
            !         f13(i,j,nz-1)=f13(i,j,1)
            !         f16(i,j,nz-1)=f16(i,j,1)
            !         f18(i,j,nz-1)=f18(i,j,1)

            !         g6(i,j,nz-1)=g6(i,j,1)
            !         g12(i,j,nz-1)=g12(i,j,1)
            !         g13(i,j,nz-1)=g13(i,j,1)
            !         g16(i,j,nz-1)=g16(i,j,1)
            !         g18(i,j,nz-1)=g18(i,j,1)

            !         f5(i,j,2)=f5(i,j,nz)
            !         f11(i,j,2)=f11(i,j,nz)
            !         f14(i,j,2)=f14(i,j,nz)
            !         f15(i,j,2)=f15(i,j,nz)
            !         f17(i,j,2)=f17(i,j,nz)

            !         g5(i,j,2)=g5(i,j,nz)
            !         g11(i,j,2)=g11(i,j,nz)
            !         g14(i,j,2)=g14(i,j,nz)
            !         g15(i,j,2)=g15(i,j,nz)
            !         g17(i,j,2)=g17(i,j,nz)
            !     enddo
            ! enddo
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
                open(110, file = 'v_xy'//write_fmtnumb(iframe)//'.out', status = 'replace')
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
                open(115, file = 'w_yz'//write_fmtnumb(iframe)//'.out', status = 'replace')
                do j=1,ny
                    do k=1,nz
                        write(111,*) psi(nx/2,j,k)  
                    enddo
                enddo
                close(111)
                do j=1,ny
                    do k=1,nz
                        write(112,*) rhoA(nx/2,j,k)   
                    enddo
                enddo
                close(112)
                do j=1,ny
                    do k=1,nz
                        write(113,*) rhoB(nx/2,j,k)    
                    enddo
                enddo
                close(113)
                do j=1,ny
                    do k=1,nz
                        write(114,*) v(nx/2,j,k)      
                    enddo
                enddo
                close(114)
                do j=1,ny
                    do k=1,nz
                        write(115,*) w(nx/2,j,k)     
                    enddo
                enddo
                close(115)
                write(6,*) "files updated at t=", step
                !
                ! open(117, file = 'psi3d'//write_fmtnumb(iframe)//'.out', status = 'replace')
                ! !$acc update self(psi)
                ! do i=1,nx,5
                !     do j=1,ny,5
                !         do k=1,nz,5
                !             write(117,*) psi(i,j,k)     
                !         enddo
                !     enddo
                ! enddo
                ! close(117)
            endif
        
    enddo 
    call cpu_time(ts2)
    !$acc end data
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)*real(1.d-9,kind=db)/(ts2-ts1)
   contains 
    !*************************************************functions************************************************!
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