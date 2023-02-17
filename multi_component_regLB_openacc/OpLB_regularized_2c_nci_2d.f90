program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    !*********************************variables
        integer, parameter :: db=4 !kind(1.0)
        integer(kind=8) :: i,j,ll,l,dumm
        integer(kind=8) :: nx,ny,step,stamp,nlinks,nsteps,ngpus,ncontact
        integer,save :: iframe=0
        
        real(kind=db),parameter :: pi_greek=3.14159265359793234626433
        
        real(kind=4)  :: ts1,ts2,radius
        real(kind=db) :: visc_LB,uu,udotc,omega,feq,geq
        real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8
        real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq1,pi2cssq2,pi2cssq0
        real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,fpc
        real(kind=db) :: addendum0,addendum1,addendum2,addendum3,addendum4,addendum5,addendum6,addendum7,addendum8
        real(kind=db) :: gaddendum0,gaddendum1,gaddendum2,gaddendum3,gaddendum4,gaddendum5,gaddendum6,gaddendum7,gaddendum8
        real(kind=db) :: psi_x,psi_y,mod_psi,mod_psi_sq,st_coeff,b0,b1,b2,beta,sigma,norm_x,norm_y
        real(kind=db) :: one_ov_nu2,one_ov_nu1,nu_avg,rtot,rprod
        real(kind=db) :: max_press_excess,ushifted,vshifted
        real(kind=db) :: rnd_n1,rnd_n2,drhok1,drhok2
        
        integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
        
        real(kind=db), allocatable, dimension(:)     :: p
        real(kind=db), allocatable, dimension(:,:) :: psi,rhoA,rhoB,u,v,pxx,pyy,pxy
        real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8
        real(kind=db), allocatable, dimension(:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8
        integer(kind=4), allocatable, dimension(:,:) :: nci_loc
        !*******************************for geometry
        integer:: ddrop,lc,jjd,jju,k

    !*********************************lattice pars  
        nlinks=8 !pari!
        cssq=1.0_db/3.0_db
        !fluid 1
        tau=0.650_db
        visc_LB=cssq*(tau-0.5_db)
        one_ov_nu1=1.0_db/visc_LB
        !fluid2
        tau=0.650_db
        visc_LB=cssq*(tau-0.5_db)
        one_ov_nu2=1.0_db/visc_LB
        omega=1.0_db/tau
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif
    !*******************************user parameters**************************
        nx=4016!500!500
        ny=4016 !500!600
        nsteps=1000
        stamp=100000
        fx=0.0_db*10.0**(-7)
        fy=-0.0_db*10.0**(-6)
    !**********************************allocation****************************
        allocate(p(0:nlinks))
        allocate(f0(0:nx+1,0:ny+1),f1(0:nx+1,0:ny+1),f2(0:nx+1,0:ny+1),f3(0:nx+1,0:ny+1),f4(0:nx+1,0:ny+1))
        allocate(f5(0:nx+1,0:ny+1),f6(0:nx+1,0:ny+1),f7(0:nx+1,0:ny+1),f8(0:nx+1,0:ny+1))
        allocate(g0(0:nx+1,0:ny+1),g1(0:nx+1,0:ny+1),g2(0:nx+1,0:ny+1),g3(0:nx+1,0:ny+1),g4(0:nx+1,0:ny+1))
        allocate(g5(0:nx+1,0:ny+1),g6(0:nx+1,0:ny+1),g7(0:nx+1,0:ny+1),g8(0:nx+1,0:ny+1))
        allocate(psi(0:nx+1,0:ny+1),rhoA(1:nx,1:ny),rhoB(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny),pxx(1:nx,1:ny),pyy(1:nx,1:ny),pxy(1:nx,1:ny))
        allocate(isfluid(1:nx,1:ny),nci_loc(1:nx,1:ny)) !,omega_2d(1:nx,1:ny)) 
        p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
        1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
        !ex=(/0,1,0,-1,0,1,-1,-1,1/)
        !ey=(/0,0,1,0,-1,1,1,-1,-1/)

    !***************************** regularized: hermite 
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
        radius=50
        isfluid=1
        isfluid(1,:)=0 !EAST
        isfluid(nx,:)=0 !WEST
        isfluid(:,1)=0 !SOUTH 
        isfluid(:,ny)=0 !NORTH
        !************************************constriction**********************!
            ! isfluid(nx/2-radius:nx/2+radius,2:ny-1)=1
            ! isfluid(2:nx-1,2:ny/2-2*radius)=1
            ! isfluid(2:nx-1,ny/2+2*radius:ny-1)=1
            ! do k=1,ny 
            !     do j=1,nx
            !         if(isfluid(j,k).eq.3)then
            !             if(isfluid(j+1,k).eq.1 .or. isfluid(j-1,k).eq.1 .or. isfluid(j,k+1).eq.1 .or. isfluid(j,k-1).eq.1 &
            !                 .or. isfluid(j+1,k+1).eq.1 .or. isfluid(j+1,k-1).eq.1 .or. isfluid(j-1,k+1).eq.1 .or. isfluid(j-1,k-1).eq.1)then
            !                 isfluid(j,k)=0
            !             endif
            !         endif
            !     enddo
            ! enddo
        !****************************************tapered**************************
            ! Lc=250+ny/2;
            ! do i=Lc,ny
            !     jjd=nint((i-Lc+1)*sind(30.0_db));
                
            !     if(jjd<nx/2-4)then
            !         isfluid(jjd:nx/2,i)=1; 
            !         jju=nx-jjd; 
            !         isfluid(nx/2:jju,i)=1;
            !     endif
            ! enddo 
            ! Ddrop=40;
            ! isfluid(2:nx-1,ny/2:Lc)=1;
            ! isfluid(nx/2-(Ddrop/(2)):nx/2+(Ddrop/(2)),ny/2:ny)=1;
            ! isfluid(:,1:ny/2)=isfluid(:,ny:ny/2+1:-1);
            ! isfluid(1,:)=0 !left
            ! isfluid(nx,:)=0 !right
            ! isfluid(:,1)=0 !front 
            ! isfluid(:,ny)=0 !rear z
            ! do k=1,ny 
            !     do j=1,nx
            !         if(isfluid(j,k).eq.3)then
            !             if(isfluid(j+1,k).eq.1 .or. isfluid(j-1,k).eq.1 .or. isfluid(j,k+1).eq.1 .or. isfluid(j,k-1).eq.1 &
            !                 .or. isfluid(j+1,k+1).eq.1 .or. isfluid(j+1,k-1).eq.1 .or. isfluid(j-1,k+1).eq.1 .or. isfluid(j-1,k-1).eq.1)then
            !                 isfluid(j,k)=0
            !             endif
            !         endif
            !     enddo
            ! enddo
        !****************************print geo**************************!
            !isfluid(2:nx-1,2:800)=1 !rear z
            open(231, file = 'isfluid.out', status = 'replace')
                do j=1,nx
                    do k=1,ny
                        write(231,*) isfluid(j,k)  
                    enddo
                enddo
            close(231)
    !*************************************initial conditions ************************    
        u=0.0_db
        v=0.0_db
        rhoA=0.0_db     !is to be intended as a delta rho
        rhoB=0.0_db     !is to be intended as a delta rho
        max_press_excess=0.002 !2
        nci_loc=0
        psi=-1.0_db
        !*******************************************impacting drops****************************!
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
        !***************************************constriction drops************************************!
            ! do i=(nx/2-radius-5)-radius,(nx/2-radius-5) +radius
            !     do j=2*radius+10-2*radius,2*radius+10+2*radius
            !         if ((i-(nx/2-radius-5))**2+(j-2*radius+10)**2<=radius**2)then
            !             psi(i,j)=1.0_db
            !         endif
            !     enddo
            ! enddo
            ! do i=(nx/2+radius+5)-radius,(nx/2+radius+5) +radius
            !     do j=2*radius+10-2*radius,2*radius+10+2*radius
            !         if ((i-(nx/2+radius+5))**2+(j-2*radius+10)**2<=radius**2)then
            !             psi(i,j)=1.0_db
            !         endif
            !     enddo
            ! enddo
        !*********************************************tapered single drops**********************************!
            ! do i=nx/2-radius,nx/2+radius
            !     do j=ny-20-radius-radius,ny-20-radius+radius
            !         if ((i-nx/2)**2+(j-(ny-20-radius))**2<=radius**2)then
            !             psi(i,j)=1.0_db
            !         endif
            !     enddo
            ! enddo
        !*********************************************tapered full system**********************!
            ! open(231, file = 'psi.txt', status = 'old',action='read')
            ! do j=1,nx
            !     do k=1,ny
            !         read(231,*) psi(j,k)
            !     enddo
            ! enddo
            ! close(231)
            ! !
            ! open(231, file = 'psi.out', status = 'replace')
            ! do i=1,nx
            !     do j=1,ny
            !         write(231,*) psi(i,j)  
            !     enddo
            ! enddo
            ! close(231)
    !*****************************************init distros*********************************!
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


    !*****************************************gpu copies**********************!
        !$acc data copy(p,rhoA,rhoB,u,v,pxx,pxy,pyy,f0,f1,f2,f3,f4,f5,f6,f7,f8,isfluid, &
        !$acc& g0,g1,g2,g3,g4,g5,g6,g7,g8,psi,nci_loc)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moment + neq pressor*********
            !$acc kernels 
            !$acc loop collapse(2) !private(uu,temp,udotc,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8) 
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j).eq.1)then
                        rhoA(i,j) = (f0(i,j)+f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j))
                        rhoB(i,j)=  (g0(i,j)+g1(i,j)+g2(i,j)+g3(i,j)+g4(i,j)+g5(i,j)+g6(i,j)+g7(i,j)+g8(i,j))
                        
                        u(i,j) = (f1(i,j)+f5(i,j) +f8(i,j)-f3(i,j) -f6(i,j) -f7(i,j) + &
                                g1(i,j)+g5(i,j) +g8(i,j)-g3(i,j) -g6(i,j) -g7(i,j))!/(rhoa(i,j)+rhoB(i,j))
                        v(i,j) = (f5(i,j) +f2(i,j) +f6(i,j)-f7(i,j) -f4(i,j) -f8(i,j) + &
                                g5(i,j) +g2(i,j) +g6(i,j)-g7(i,j) -g4(i,j) -g8(i,j) )!/(rhoa(i,j)+rhoB(i,j))

                        psi(i,j)= (rhoA(i,j)-rhoB(i,j))/(rhoA(i,j)+rhoB(i,j))  
                        nci_loc(i,j)=0    

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
                    endif
                enddo
            enddo
        !********************collision + no slip + forcing: fused implementation*********
            !$acc loop collapse(2)  
            do j=1,ny
                do i=1,nx 
                    if(isfluid(i,j).eq.1)then
                        !************************* pressure excess **********!
                        if(psi(i,j).lt.-0.9_db .and. i.gt.3 .and. i.lt.nx-2 .and. j.gt.3 .and. j.lt.ny-2)then
                            
                            if (psi(i+3,j).gt.-0.85 .and. psi(i-3,j).gt.-0.85)then !&
                                nci_loc(i+3,j)=1
                                nci_loc(i-3,j)=1
                            endif
                            if (psi(i,j+3).gt.-0.85 .and. psi(i,j-3).gt.-0.85)then
                               nci_loc(i,j+3)=1
                               nci_loc(i,j-3)=1
                            endif
                            if (psi(i+2,j+2).gt.-0.85 .and. psi(i-2,j-2).gt.-0.85 )then
                                nci_loc(i+2,j+2)=1
                                nci_loc(i-2,j-2)=1
                            endif
                            if (psi(i-2,j+2).gt.-0.85 .and. psi(i+2,j-2).gt.-0.85 )then
                                nci_loc(i-2,j+2)=1
                                nci_loc(i+2,j-2)=1
                            endif
                        endif
                        !*************************pressure excess**********!
                    endif
                enddo
            enddo
            !$acc loop collapse(2)  
            do j=1,ny
                do i=1,nx 
                    if(isfluid(i,j).eq.1)then 
                        !oneminusuu= -uu !1.0_db - uu
                        !0
                        !chromodynamic
                        psi_x=(1.0_db/cssq)*(p(1)*(psi(i+1,j)-psi(i-1,j)) + p(5)*(psi(i+1,j+1) + psi(i+1,j-1)-psi(i-1,j+1)-psi(i-1,j-1)))
                        psi_y=(1.0_db/cssq)*(p(1)*(psi(i,j+1)-psi(i,j-1)) + p(5)*(psi(i+1,j+1) - psi(i+1,j-1)+psi(i-1,j+1)-psi(i-1,j-1)))
                        mod_psi=sqrt(psi_x**2+psi_y**2)
                        mod_psi_sq=psi_x**2+psi_y**2 
                        norm_x=0.0_db
                        norm_y=0.0_db
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
                            norm_x=psi_x/mod_psi
                            norm_y=psi_y/mod_psi
                            ! se psi interfaccia vicina a 3-4-5lu è piu' grande del mio valore allora applico nci
                        endif
                        !regularized collision + perturbation + recolouring
                        ushifted=u(i,j) + fx + float(nci_loc(i,j))*(norm_x)*max_press_excess*abs(rhob(i,j))
                        vshifted=v(i,j) + fy + float(nci_loc(i,j))*(norm_y)*max_press_excess*abs(rhob(i,j))
                        uu=0.5_db*(ushifted*ushifted + vshifted*vshifted)/cssq
                        feq=p(0)*(rtot-uu)
                        fpc=feq + (1.0_db-omega)*pi2cssq0*(- cssq*pyy(i,j)-cssq*pxx(i,j))  + addendum0
                        f0(i,j)=fpc*(rhoA(i,j))/rtot 
                        g0(i,j)=fpc*(rhoB(i,j))/rtot
                        !1
                        udotc=ushifted/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(1)*(rtot+(temp + udotc))
                        fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j) - cssq*pyy(i,j) )+ addendum1 !+ (fx+float(nci_loc(i,j))*(norm_x)*max_press_excess*abs(rhoa(i,j)))*p(1)/cssq 
                        f1(i+1,j)= fpc*(rhoA(i,j))/rtot + gaddendum1
                        g1(i+1,j)= fpc*(rhoB(i,j))/rtot - gaddendum1
                        !3
                        feq=p(3)*(rtot+(temp - udotc))
                        fpc=feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j) - cssq*pyy(i,j)) + addendum3 !- (fx+float(nci_loc(i,j))*(norm_x)*max_press_excess*abs(rhoa(i,j)))*p(3)/cssq 
                        f3(i-1,j)= fpc*(rhoA(i,j))/rtot + gaddendum3 
                        g3(i-1,j)= fpc*(rhoB(i,j))/rtot - gaddendum3 
                        !2
                        udotc=vshifted/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(2)*(rtot+(temp + udotc))
                        fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j)) + addendum2 !+ (fy+float(nci_loc(i,j))*(norm_y)*max_press_excess*abs(rhoa(i,j)))*p(2)/cssq  !
                        f2(i,j+1)= fpc*(rhoA(i,j))/rtot + gaddendum2  
                        g2(i,j+1)= fpc*(rhoB(i,j))/rtot - gaddendum2   
                        !4
                        feq=p(4)*(rtot+(temp - udotc))
                        fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j)) + addendum4 !- (fy+float(nci_loc(i,j))*(norm_y)*max_press_excess*abs(rhoa(i,j)))*p(4)/cssq 
                        f4(i,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum4 
                        g4(i,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum4 
                        !5
                        udotc=(ushifted+vshifted)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(5)*(rtot+(temp + udotc))
                        fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + addendum5 !+ (fx + fy + float(nci_loc(i,j))*(norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(5)/cssq 
                        f5(i+1,j+1)=fpc*(rhoA(i,j))/rtot + gaddendum5
                        g5(i+1,j+1)=fpc*(rhoB(i,j))/rtot - gaddendum5
                        !7
                        feq=p(7)*(rtot+(temp - udotc))
                        fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + addendum7 !- (fx + fy + float(nci_loc(i,j))*(norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(7)/cssq 
                        f7(i-1,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum7
                        g7(i-1,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum7
                        !6
                        udotc=(-ushifted+vshifted)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(6)*(rtot+(temp + udotc))
                        fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) + addendum6 !+(-fx + fy + float(nci_loc(i,j))*(-norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(6)/cssq 
                        f6(i-1,j+1)= fpc*(rhoA(i,j))/rtot + gaddendum6
                        g6(i-1,j+1)= fpc*(rhoB(i,j))/rtot - gaddendum6
                        !8
                        feq=p(8)*(rtot+(temp - udotc))
                        fpc=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j))  + addendum8 !+( fx - fy + float(nci_loc(i,j))*(norm_x-norm_y)*max_press_excess*abs(rhoa(i,j)))*p(8)/cssq 
                        f8(i+1,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum8 
                        g8(i+1,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum8 
                    endif
                enddo
            enddo
        
        !********************************************bcs no slip*****************************************!
            !$acc loop independent 
            do j=1,ny
                !$acc loop independent 
                do i=1,nx
                    if(isfluid(i,j).eq.0)then
                        psi(i,j)=-1.0_db

                        f8(i+1,j-1)=f6(i,j)!gpc 
                        f7(i-1,j-1)=f5(i,j)!hpc

                        f6(i-1,j+1)=f8(i,j)!gpc 
                        f5(i+1,j+1)=f7(i,j)!hpc 


                        f4(i,j-1)=f2(i,j)!gpc 
                        f3(i-1,j)=f1(i,j)!hpc 

                        f2(i,j+1)=f4(i,j)!gpc 
                        f1(i+1,j)=f3(i,j)!hpc 
                        !************************!
                        g8(i+1,j-1)=g6(i,j)!gpc 
                        g7(i-1,j-1)=g5(i,j)!hpc

                        g6(i-1,j+1)=g8(i,j)!gpc 
                        g5(i+1,j+1)=g7(i,j)!hpc 


                        g4(i,j-1)=g2(i,j)!gpc 
                        g3(i-1,j)=g1(i,j)!hpc 

                        g2(i,j+1)=g4(i,j)!gpc 
                        g1(i+1,j)=g3(i,j)!hpc
                    endif
                enddo
            enddo
        
        !******************************************call periodic bcs************************
                  !$acc loop independent 
                   do i=1,nx  
                    psi(i,ny)=psi(i,2)
                    psi(i,1)=psi(i,ny-1)
                    !negative lungo y
                    f4(i,ny-1)=f4(i,1)
                    f7(i,ny-1)=f7(i,1)
                    f8(i,ny-1)=f8(i,1)

                    g4(i,ny-1)=g4(i,1)
                    g7(i,ny-1)=g7(i,1)
                    g8(i,ny-1)=g8(i,1)

                    f2(i,2)=f2(i,ny)
                    f5(i,2)=f5(i,ny)
                    f6(i,2)=f6(i,ny)

                    g2(i,2)=g2(i,ny)
                    g5(i,2)=g5(i,ny)
                    g6(i,2)=g6(i,ny)
                    
                !     psi(i,1)=psi(i,2)
                !     !pos lungo y
                !     f2(i,2)=0.0
                !     f5(i,2)=0.0
                !     f6(i,2)=0.0

                    ! g2(i,2)=p(2)*1.0_db 
                    ! g5(i,2)=p(5)*1.0_db 
                    ! g6(i,2)=p(6)*1.0_db  
                    
                    ! g4(i,ny-1)=p(2)*1.0_db 
                    ! g8(i,ny-1)=p(5)*1.0_db 
                    ! g7(i,ny-1)=p(6)*1.0_db 

                  enddo
            !$acc end kernels 
        !****************************************writeonfile***************************************************!
            if(mod(step,stamp).eq.0)then
            !$acc update self(psi(1:nx,1:ny),rhoB(1:nx,1:ny),rhoA(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny))
            iframe=iframe+1
            open(101, file = 'psi'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(102, file = 'rhoA'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(103, file = 'rhoB'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(104, file = 'u'//write_fmtnumb(iframe)//'.out', status = 'replace')
            open(105, file = 'v'//write_fmtnumb(iframe)//'.out', status = 'replace')
            do i=1,nx
                do j=1,ny
                    write(101,*) psi(i,j)  
                enddo
            enddo
            close(101)
            do i=1,nx
                do j=1,ny
                    write(102,*) rhoA(i,j)  
                enddo
            enddo
            close(102)
            do i=1,nx
                do j=1,ny
                    write(103,*) rhoB(i,j)  
                enddo
            enddo
            close(103)
            do i=1,nx
                do j=1,ny
                    write(104,*) u(i,j)  
                enddo
            enddo
            close(104)
            do i=1,nx
                do j=1,ny
                    write(105,*) v(i,j)  
                enddo
            enddo
            close(105)
            write(6,*) "files updated at t=", step
            endif
        !
    enddo 
    call cpu_time(ts2)

    !$acc end data

    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  nx*ny*nsteps/10.0_db**9/ts2-ts1

  contains 
  !*****************************************************functions********************************************************!
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
