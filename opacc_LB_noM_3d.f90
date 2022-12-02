program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    integer, parameter :: db=4 !kind(1.0)
    integer(kind=8) :: i,j,k
    integer(kind=8) :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2,p0,p1,p2
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp,dummy
    
    integer(kind=4), allocatable,  dimension(:)     :: ex,ey,ez,opp
    integer(kind=4), allocatable,  dimension(:,:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p,dex,dey,dez,fdum
    real(kind=db) :: rho,u,v,w
    real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
    real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18

       
    nlinks=18 !pari!
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
    nx=512
    ny=512
    nz=512
    nsteps=100
    stamp=1000
    fx=1.0_db*10.0**(-5)
    fy=0.0_db*10.0**(-5)
    fz=0.0_db*10.0**(-5)
    allocate(p(0:nlinks),ex(0:nlinks),ey(0:nlinks),dex(0:nlinks),ez(0:nlinks))
    allocate(dey(0:nlinks),dez(0:nlinks),opp(0:nlinks),fdum(0:nlinks))
    allocate(f0(1:nx,1:ny,1:nz),f1(1:nx,1:ny,1:nz),f2(1:nx,1:ny,1:nz),f3(1:nx,1:ny,1:nz))
    allocate(f4(1:nx,1:ny,1:nz),f5(1:nx,1:ny,1:nz),f6(1:nx,1:ny,1:nz),f7(1:nx,1:ny,1:nz))
    allocate(f8(1:nx,1:ny,1:nz),f9(1:nx,1:ny,1:nz),f10(1:nx,1:ny,1:nz),f11(1:nx,1:ny,1:nz))
    allocate(f12(1:nx,1:ny,1:nz),f13(1:nx,1:ny,1:nz),f14(1:nx,1:ny,1:nz),f15(1:nx,1:ny,1:nz))
    allocate(f16(1:nx,1:ny,1:nz),f17(1:nx,1:ny,1:nz),f18(1:nx,1:ny,1:nz))
    allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
    
    ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
	ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
	ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)

    dex=real(ex,db)
    dey=real(ey,db)
    dez=real(ez,db)

    p0=(1.0_db/3.0_db)
	p1=(1.0_db/18.0_db)
	p2=(1.0_db/36.0_db)
    p=(/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2/)
    opp=(/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
    omega=1.0_db/tau
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:,:)=0 !left
    isfluid(nx,:,:)=0 !right
    isfluid(:,1,:)=0 !front 
    isfluid(:,ny,:)=0 !rear
    isfluid(:,:,1)=0 !bottom
    isfluid(:,:,nz)=0 !top
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    w=0.0_db
    rho=0.0_db     !is to be intended as a delta rho
    !do ll=0,nlinks
    f0(1:nx,1:ny,1:nz)=0.0_db
    f1(1:nx,1:ny,1:nz)=0.0_db
    f2(1:nx,1:ny,1:nz)=0.0_db
    f3(1:nx,1:ny,1:nz)=0.0_db
    f4(1:nx,1:ny,1:nz)=0.0_db
    f5(1:nx,1:ny,1:nz)=0.0_db
    f6(1:nx,1:ny,1:nz)=0.0_db
    f7(1:nx,1:ny,1:nz)=0.0_db
    f8(1:nx,1:ny,1:nz)=0.0_db
    f9(1:nx,1:ny,1:nz)=0.0_db
    f10(1:nx,1:ny,1:nz)=0.0_db
    f11(1:nx,1:ny,1:nz)=0.0_db
    f12(1:nx,1:ny,1:nz)=0.0_db
    f13(1:nx,1:ny,1:nz)=0.0_db
    f14(1:nx,1:ny,1:nz)=0.0_db
    f15(1:nx,1:ny,1:nz)=0.0_db
    f16(1:nx,1:ny,1:nz)=0.0_db
    f17(1:nx,1:ny,1:nz)=0.0_db
    f18(1:nx,1:ny,1:nz)=0.0_db
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
    pause


    !$acc data copy(ex,ey,ez,dex,dey,dez,p, &
    !$acc & f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,isfluid)
    
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !**************************************** moments******************************
    !    !$acc kernels
    !     do k=2,nz-1
    !         do j=2,ny-1
    !             do i=2,nx-1         
    !                 rho(i,j,k) = f0(i,j,k) &
    !                         +f1(i,j,k) &
    !                         +f2(i,j,k) &
    !                         +f3(i,j,k) &
    !                         +f4(i,j,k) &
    !                         +f5(i,j,k) &
    !                         +f6(i,j,k) &
    !                         +f7(i,j,k) &
    !                         +f8(i,j,k) &
    !                         +f9(i,j,k) &
    !                         +f10(i,j,k) &
    !                         +f11(i,j,k) &
    !                         +f12(i,j,k) &
    !                         +f13(i,j,k) &
    !                         +f14(i,j,k) &
    !                         +f15(i,j,k) &
    !                         +f16(i,j,k) &
    !                         +f17(i,j,k) &
    !                         +f18(i,j,k)

    !                 u(i,j,k) = (f1(i,j,k)+f7(i,j,k)+f9(i,j,k)+f15(i,j,k)+f8(i,j,k)) &
    !                          -(f2(i,j,k)+f8(i,j,k)+f10(i,j,k)+f16(i,j,k)+f17(i,j,k)) 
                        
    !                 v(i,j,k) = (f3(i,j,k)+f7(i,j,k)+f10(i,j,k)+f11(i,j,k)+f13(i,j,k)) &
    !                         -(f4(i,j,k)+f8(i,j,k)+f9(i,j,k)+f12(i,j,k)+f14(i,j,k))

    !                 w(i,j,k) = (f5(i,j,k)+f11(i,j,k)+f14(i,j,k)+f15(i,j,k)+f17(i,j,k)) &
    !                         -(f6(i,j,k)+f12(i,j,k)+f13(i,j,k)+f16(i,j,k)+f18(i,j,k))
                    
    !             enddo
    !         enddo 
    !     enddo
    !     !$acc end kernels
        !***********************************collision&forcing************************ 
        !!$acc update host(rho,u,v)
        !!$acc update device(rho,u,v)
        !$acc kernels !!private(uu,temp,udotc,feq,dummy) 
        !$acc loop private(uu,temp,udotc,feq,dummy,u,v,w,rho)
        do k=2,nz-1
            do j=2,ny-1
                do i=2,nx-1
                    if(isfluid(i,j,k).eq.1)then
                        rho = f0(i,j,k) &
                            +f1(i,j,k) &
                            +f2(i,j,k) &
                            +f3(i,j,k) &
                            +f4(i,j,k) &
                            +f5(i,j,k) &
                            +f6(i,j,k) &
                            +f7(i,j,k) &
                            +f8(i,j,k) &
                            +f9(i,j,k) &
                            +f10(i,j,k) &
                            +f11(i,j,k) &
                            +f12(i,j,k) &
                            +f13(i,j,k) &
                            +f14(i,j,k) &
                            +f15(i,j,k) &
                            +f16(i,j,k) &
                            +f17(i,j,k) &
                            +f18(i,j,k)

                        u = (f1(i,j,k)+f7(i,j,k)+f9(i,j,k)+f15(i,j,k)+f8(i,j,k)) &
                             -(f2(i,j,k)+f8(i,j,k)+f10(i,j,k)+f16(i,j,k)+f17(i,j,k)) 
                        
                        v = (f3(i,j,k)+f7(i,j,k)+f10(i,j,k)+f11(i,j,k)+f13(i,j,k)) &
                            -(f4(i,j,k)+f8(i,j,k)+f9(i,j,k)+f12(i,j,k)+f14(i,j,k))

                        w = (f5(i,j,k)+f11(i,j,k)+f14(i,j,k)+f15(i,j,k)+f17(i,j,k)) &
                            -(f6(i,j,k)+f12(i,j,k)+f13(i,j,k)+f16(i,j,k)+f18(i,j,k))
                        !
                        uu=0.5_db*(u*u + v*v + w*w)/cssq
                        !0
                        feq=p(0)*(rho-uu)
                        f0(i,j,k)=f0(i,j,k) + omega*(feq - f0(i,j,k)) 
                        
                        !1
                        udotc=u/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(1)*(rho+(temp + udotc))
                        f1(i,j,k)=f1(i,j,k) + omega*(feq - f1(i,j,k)) + fx*p(1)/cssq
                        
                        !2
                        feq=p(2)*(rho+(temp - udotc))
                        f2(i,j,k)=f2(i,j,k) + omega*(feq - f2(i,j,k)) - fx*p(2)/cssq
                        
                        !3
                        udotc=v/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(3)*(rho+(temp + udotc))
                        f3(i,j,k)=f3(i,j,k) + omega*(feq - f3(i,j,k)) + fy*p(3)/cssq
                        
                        !4
                        feq=p(4)*(rho+(temp - udotc))
                        f4(i,j,k)=f4(i,j,k) + omega*(feq - f4(i,j,k)) - fy*p(4)/cssq
                        
                        !7
                        udotc=(u+v)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(7)*(rho+(temp + udotc))
                        f7(i,j,k)=f7(i,j,k) + omega*(feq - f7(i,j,k)) + (fx+fy)*p(7)/cssq  
                        
                        !8
                        feq=p(8)*(rho+(temp - udotc))
                        f8(i,j,k)=f8(i,j,k) + omega*(feq - f8(i,j,k)) - (fx+fy)*p(8)/cssq
                        
                        !10
                        udotc=(-u+v)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(10)*(rho+(temp + udotc))
                        f10(i,j,k)=f10(i,j,k) + omega*(feq - f10(i,j,k)) +(fy-fx)*p(10)/cssq
                        
                        !9
                        feq=p(9)*(rho+(temp - udotc))
                        f9(i,j,k)=f9(i,j,k) + omega*(feq - f9(i,j,k)) + (fx-fy)*p(9)/cssq

                        !5
                        udotc=w/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(5)*(rho+(temp + udotc))
                        f5(i,j,k)=f5(i,j,k) + omega*(feq - f5(i,j,k)) + fz*p(5)/cssq
                        
                        !6
                        feq=p(6)*(rho+(temp - udotc))
                        f6(i,j,k)=f6(i,j,k) + omega*(feq - f6(i,j,k)) - fz*p(6)/cssq

                        !15
                        udotc=(u+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(15)*(rho+(temp + udotc))
                        f15(i,j,k)=f15(i,j,k) + omega*(feq - f15(i,j,k)) + (fx+fz)*p(15)/cssq  
                        
                        !16
                        feq=p(16)*(rho+(temp - udotc))
                        f16(i,j,k)=f16(i,j,k) + omega*(feq - f16(i,j,k)) - (fx+fz)*p(16)/cssq

                        !17
                        udotc=(-u+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(17)*(rho+(temp + udotc))
                        f17(i,j,k)=f17(i,j,k) + omega*(feq - f17(i,j,k)) +(fz-fx)*p(17)/cssq
                        
                        !18
                        feq=p(18)*(rho+(temp - udotc))
                        f18(i,j,k)=f18(i,j,k) + omega*(feq - f18(i,j,k)) + (fx-fz)*p(18)/cssq

                        !11
                        udotc=(v+w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(11)*(rho+(temp + udotc))
                        f11(i,j,k)=f11(i,j,k) + omega*(feq - f11(i,j,k)) +(fy+fz)*p(11)/cssq
                        
                        !12
                        feq=p(12)*(rho+(temp - udotc))
                        f12(i,j,k)=f12(i,j,k) + omega*(feq - f12(i,j,k)) - (fy+fz)*p(12)/cssq

                        !13
                        udotc=(v-w)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(13)*(rho+(temp + udotc))
                        f13(i,j,k)=f13(i,j,k) + omega*(feq - f13(i,j,k)) + (fy-fz)*p(13)/cssq
                        
                        !14
                        feq=p(14)*(rho+(temp - udotc))
                        f14(i,j,k)=f14(i,j,k) + omega*(feq - f14(i,j,k)) + (fz-fy)*p(14)/cssq

                    elseif(isfluid(i,j,k).eq.0)then
                        !
                        dummy=f1(i,j,k)
                        f1(i,j,k)=f2(i,j,k)
                        f2(i,j,k)=dummy
                        !
                        dummy=f3(i,j,k)
                        f3(i,j,k)=f4(i,j,k)
                        f4(i,j,k)=dummy
                        !
                        dummy=f5(i,j,k)
                        f5(i,j,k)=f6(i,j,k)
                        f6(i,j,k)=dummy
                        !
                        dummy=f7(i,j,k)
                        f7(i,j,k)=f8(i,j,k)
                        f8(i,j,k)=dummy
                        !
                        dummy=f9(i,j,k)
                        f9(i,j,k)=f10(i,j,k)
                        f10(i,j,k)=dummy
                        !
                        dummy=f11(i,j,k)
                        f11(i,j,k)=f12(i,j,k)
                        f12(i,j,k)=dummy
                        !
                        dummy=f13(i,j,k)
                        f13(i,j,k)=f14(i,j,k)
                        f14(i,j,k)=dummy
                        !
                        dummy=f15(i,j,k)
                        f15(i,j,k)=f16(i,j,k)
                        f16(i,j,k)=dummy
                        !
                        dummy=f17(i,j,k)
                        f17(i,j,k)=f18(i,j,k)
                        f18(i,j,k)=dummy
                        
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        !******************************************call bcs************************
        !!$acc update host(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !periodic along y
        !$acc kernels 
        !x=1 
        f1(1,2:ny-1,2:nz-1)=f1(nx-1,2:ny-1,2:nz-1)
        f7(1,2:ny-1,2:nz-1)=f7(nx-1,2:ny-1,2:nz-1)
        f9(1,2:ny-1,2:nz-1)=f9(nx-1,2:ny-1,2:nz-1)
        f15(1,2:ny-1,2:nz-1)=f15(nx-1,2:ny-1,2:nz-1)
        f8(1,2:ny-1,2:nz-1)=f18(nx-1,2:ny-1,2:nz-1)
        !
        !x=nx 
        f2(nx,2:ny-1,2:nz-1)=f2(2,2:ny-1,2:nz-1)
        f8(nx,2:ny-1,2:nz-1)=f8(2,2:ny-1,2:nz-1)
        f10(1,2:ny-1,2:nz-1)=f10(2,2:ny-1,2:nz-1)
        f16(1,2:ny-1,2:nz-1)=f16(2,2:ny-1,2:nz-1)
        f17(1,2:ny-1,2:nz-1)=f17(2,2:ny-1,2:nz-1)
        !y=1
        f3(2:nx-1,1,1:nz)=f3(2:nx-1,ny-1,1:nz)
        f7(2:nx-1,1,1:nz)=f7(2:nx-1,ny-1,1:nz)
        f10(2:nx-1,1,1:nz)=f10(2:nx-1,ny-1,1:nz)
        f11(2:nx-1,1,1:nz)=f11(2:nx-1,ny-1,1:nz)
        f13(2:nx-1,1,1:nz)=f13(2:nx-1,ny-1,1:nz)
        !y=ny
        f4(2:nx-1,ny,1:nz)=f4(2:nx-1,2,1:nz)
        f8(2:nx-1,ny,1:nz)=f8(2:nx-1,2,1:nz)
        f9(2:nx-1,ny,1:nz)=f9(2:nx-1,2,1:nz)
        f12(2:nx-1,ny,1:nz)=f12(2:nx-1,2,1:nz)
        f14(2:nx-1,ny,1:nz)=f14(2:nx-1,2,1:nz)
        !$acc end kernels
        ! !******************************************streaming***************************  
        !
        !!$acc update host(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !!$acc update device(f0,f1,f2,f3,f4,f5,f6,f7,f8)
        !$acc kernels
        ! d2q9 xy
        f1(2:nx,2:ny-1,2:nz-1)=f1(1:nx-1,2:ny-1,2:nz-1)
        f2(1:nx-1,2:ny-1,2:nz-1)=f2(2:nx,2:ny-1,2:nz-1)
        ! 
        f3(2:nx-1,2:ny,2:nz-1)=f3(2:nx-1,1:ny-1,2:nz-1)
        f4(2:nx-1,1:ny-1,2:nz-1)=f4(2:nx-1,2:ny,2:nz-1)
        ! 
        f7(2:nx,2:ny,2:nz-1)=f7(1:nx-1,1:ny-1,2:nz-1)
        f8(1:nx-1,1:ny-1,2:nz-1)=f8(2:nx,2:ny,2:nz-1)
        ! 
        f9(1:nx-1,2:ny,2:nz-1)=f9(2:nx,1:ny-1,2:nz-1)
        f10(2:nx,1:ny-1,2:nz-1)=f10(1:nx,2:ny,2:nz-1)
        !
        ! d2q9 yz
        f5(2:nx-1,2:ny-1,2:nz)=f5(2:nx-1,2:ny-1,1:nz-1)
        f6(2:nx-1,2:ny-1,1:nz-1)=f6(2:nx-1,2:ny-1,2:nz)
        !
        f11(2:nx-1,2:ny,2:nz)=f11(2:nx-1,1:ny-1,1:nz-1)
        f12(2:nx-1,1:ny-1,1:nz-1)=f12(2:nx-1,2:ny,2:nz)
        !
        f13(2:nx-1,2:ny,1:nz-1)=f13(2:nx-1,1:ny-1,2:nz)
        f14(2:nx-1,1:ny-1,2:nz)=f14(2:nx-1,2:ny,1:nz-1)
        !
        !d2q9 xz
        f15(2:nx,2:ny-1,2:nz)=f15(1:nx-1,2:ny-1,1:nz-1)
        f16(1:nx-1,2:ny-1,1:nz-1)=f16(2:nx,2:ny-1,2:nz)
        !
        f18(2:nx,2:ny-1,1:nz-1)=f18(1:nx-1,2:ny-1,2:nz)
        f17(1:nx-1,2:ny-1,2:nz)=f17(2:nx,2:ny-1,1:nz-1)
        !$acc end kernels
        
    enddo 
    call cpu_time(ts2)
    !!$acc update host(rho,u,v)
    !$acc end data
    
    !write(6,*) 'u=',u(nx/2,ny/2,nz/2),'v=',v(nx/2,ny/2,nz/2),'w=',w(nx/2,ny/2,nz/2),'rho=',rho(nx/2,ny/2,nz/2)
    write(6,*) 'You''ve just wasted ', ts2-ts1, ' s of your life time' 

    
end program