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
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp

    real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq17
    real(kind=db) :: fneq9,fneq10,fneq11,fneq12,fneq13,fneq14,fneq15,fneq16,fneq18
    real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
    real(kind=db) :: pi2cssq1,pi2cssq2
    
    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    real(kind=db), allocatable, dimension(:,:,:) :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
    real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18

       
    nlinks=18 !pari!
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
    nx=512
    ny=512
    nz=512
    nsteps=100
    stamp=1000
    fx=1.0_db*10.0**(-7)
    fy=0.0_db*10.0**(-5)
    fz=0.0_db*10.0**(-5)

    allocate(f0(0:nx+1,0:ny+1,0:nz+1),f1(0:nx+1,0:ny+1,0:nz+1),f2(0:nx+1,0:ny+1,0:nz+1),f3(0:nx+1,0:ny+1,0:nz+1))
    allocate(f4(0:nx+1,0:ny+1,0:nz+1),f5(0:nx+1,0:ny+1,0:nz+1),f6(0:nx+1,0:ny+1,0:nz+1),f7(0:nx+1,0:ny+1,0:nz+1))
    allocate(f8(0:nx+1,0:ny+1,0:nz+1),f9(0:nx+1,0:ny+1,0:nz+1),f10(0:nx+1,0:ny+1,0:nz+1),f11(0:nx+1,0:ny+1,0:nz+1))
    allocate(f12(0:nx+1,0:ny+1,0:nz+1),f13(0:nx+1,0:ny+1,0:nz+1),f14(0:nx+1,0:ny+1,0:nz+1),f15(0:nx+1,0:ny+1,0:nz+1))
    allocate(f16(0:nx+1,0:ny+1,0:nz+1),f17(0:nx+1,0:ny+1,0:nz+1),f18(0:nx+1,0:ny+1,0:nz+1))
    allocate(rho(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
    allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
    allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz))
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

    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    w=0.0_db
    rho=0.0_db  !is to be intended as a delta rho
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
    !$acc data copy(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,isfluid,p0,p1,p2,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w)
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moments collision bbck + forcing************************ 
        
        !$acc kernels 
        !$acc loop collapse (3) !private(fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq9,fneq10,fneq11,&
        !!$acc& fneq12,fneq3,fneq14,fneq15,uu,temp,udotc)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1.or.isfluid(i,j,k).eq.0)then
                        rho(i,j,k) = f0(i,j,k)+f1(i,j,k)+f2(i,j,k)+f3(i,j,k)+f4(i,j,k)+f5(i,j,k) &
                            +f6(i,j,k)+f7(i,j,k)+f8(i,j,k)+f9(i,j,k)+f10(i,j,k)+f11(i,j,k) &
                            +f12(i,j,k)+f13(i,j,k)+f14(i,j,k)+f15(i,j,k)+f16(i,j,k)+f17(i,j,k) &
                            +f18(i,j,k)

                        u(i,j,k) = (f1(i,j,k)+f7(i,j,k)+f9(i,j,k)+f15(i,j,k)+f18(i,j,k)) &
                             -(f2(i,j,k)+f8(i,j,k)+f10(i,j,k)+f16(i,j,k)+f17(i,j,k)) 
                        
                        v(i,j,k) = (f3(i,j,k)+f7(i,j,k)+f10(i,j,k)+f11(i,j,k)+f13(i,j,k)) &
                            -(f4(i,j,k)+f8(i,j,k)+f9(i,j,k)+f12(i,j,k)+f14(i,j,k))

                        w(i,j,k) = (f5(i,j,k)+f11(i,j,k)+f14(i,j,k)+f15(i,j,k)+f17(i,j,k)) &
                            -(f6(i,j,k)+f12(i,j,k)+f13(i,j,k)+f16(i,j,k)+f18(i,j,k))
                        
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !1-2
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f1(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        fneq2=f2(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        !3-4
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq3=f3(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        fneq4=f4(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        !5-6
                        udotc=w(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq5=f5(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        fneq6=f6(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        !7-8
                        udotc=(u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq7=f7(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq8=f8(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        !10-9
                        udotc=(-u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq10=f10(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq9=f9(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        !11-12
                        udotc=(v(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq11=f11(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq12=f12(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        !13-14
                        udotc=(v(i,j,k)-w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq13=f13(i,j,k) - p2*(rho(i,j,k)+(temp + udotc))
                        fneq14=f14(i,j,k) - p2*(rho(i,j,k)+(temp - udotc))
                        !15-16
                        udotc=(u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq15=f15(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq16=f16(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        !17-18
                        udotc=(-u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq17=f17(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq18=f18(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        pxx(i,j,k)=fneq1+fneq2+fneq7+fneq8+fneq9+fneq10+fneq15+fneq16+fneq17+fneq18
                        pyy(i,j,k)=fneq3+fneq4+fneq7+fneq8+fneq9+fneq10+fneq11+fneq12+fneq13+fneq14
                        pzz(i,j,k)=fneq5+fneq6+fneq11+fneq12+fneq13+fneq14+fneq15+fneq16+fneq17+fneq18
                        pxy(i,j,k)= fneq7+fneq8-fneq9-fneq10
                        pxz(i,j,k)=fneq15+fneq16-fneq17-fneq18
                        pyz(i,j,k)=fneq11+fneq12-fneq13-fneq14
                    endif
                    !no slip everywhere, always before fused: to be modified for generic pressure/velocity bcs
                    if(isfluid(i,j,k).eq.0)then
                        f1(i,j,k)=p1*rho(i,j,k) + pi2cssq1*qxx*pxx(i,j,k)
                        f2(i,j,k)=p1*rho(i,j,k) + pi2cssq1*qxx*pxx(i,j,k)
                        f3(i,j,k)=p1*rho(i,j,k) + pi2cssq1*qyy*pyy(i,j,k)
                        f4(i,j,k)=p1*rho(i,j,k) + pi2cssq1*qyy*pyy(i,j,k)
                        f5(i,j,k)=p1*rho(i,j,k) + pi2cssq1*qzz*pzz(i,j,k)
                        f6(i,j,k)=p1*rho(i,j,k) + pi2cssq1*qzz*pzz(i,j,k)

                        f7(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k))
                        f8(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k))
                        f9(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k))
                        f10(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k))

                        f11(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k))
                        f12(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k))
                        f13(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k))
                        f14(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k))

                        f15(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k))
                        f16(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k))
                        f17(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k))
                        f18(i,j,k)=p2*rho(i,j,k) + pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k))
                    endif
                enddo
            enddo
        enddo
        !!$acc end kernels
        !!$acc kernels
        !$acc loop collapse (3) !private(feq,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq9,fneq10,fneq11,&
        !!$acc& fneq12,fneq3,fneq14,fneq15,uu,temp,udotc)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                        !
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !0
                        feq=p0*(rho(i,j,k)-uu)
                        f0(i,j,k)=f0(i,j,k) + omega*(feq - f0(i,j,k)) 
                        
                        !1
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq1=(1.0_db-omega)*pi2cssq1*qxx*pxx(i,j,k)
                        f1(i+1,j,k)=feq + fneq1 + fx*p1dcssq
                        
                        !2
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f2(i-1,j,k)=feq + fneq1 - fx*p1dcssq
                        
                        !3
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq3=(1.0_db-omega)*pi2cssq1*qyy*pyy(i,j,k)
                        f3(i,j+1,k)=feq+fneq3 + fy*p1dcssq
                        
                        !4
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f4(i,j-1,k)=feq+fneq3 - fy*p1dcssq
                        
                        !7
                        udotc=(u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq7=(1.0_db-omega)*(pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k)))
                        f7(i+1,j+1,k)=feq + fneq7 + (fx+fy)*p2dcssq 
                        
                        !8
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f8(i-1,j-1,k)=feq + fneq7 - (fx+fy)*p2dcssq
                        
                        !10
                        udotc=(-u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq10=(1.0_db-omega)*(pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k)))
                        f10(i-1,j+1,k)=feq+fneq10 +(fy-fx)*p2dcssq
                        
                        !9
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f9(i+1,j-1,k)=feq+fneq10 + (fx-fy)*p2dcssq

                        !5
                        udotc=w(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq5=(1.0_db-omega)*pi2cssq1*qzz*pzz(i,j,k)
                        f5(i,j,k+1)=feq+fneq5 + fz*p1dcssq
                        
                        !6
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f6(i,j,k-1)=feq+fneq5 - fz*p1dcssq

                        !15
                        udotc=(u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq15=(1.0_db-omega)*(pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k)))
                        f15(i+1,j,k+1)=feq+fneq15 + (fx+fz)*p2dcssq 
                        
                        !16
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f16(i-1,j,k-1)=feq+fneq15 - (fx+fz)*p2dcssq

                        !17
                        udotc=(-u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq17=(1.0_db-omega)*(pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k)))
                        f17(i-1,j,k+1)=feq+fneq17 +(fz-fx)*p2dcssq
                        
                        !18
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f18(i+1,j,k-1)=feq+fneq17 + (fx-fz)*p2dcssq

                        !11
                        udotc=(v(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq11=(1.0_db-omega)*(pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k)))
                        f11(i,j+1,k+1)=feq+fneq11+(fy+fz)*p2dcssq
                        
                        !12
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f12(i,j-1,k-1)=feq+fneq11 - (fy+fz)*p2dcssq

                        !13
                        udotc=(v(i,j,k)-w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq13=(1.0_db-omega)*(pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k)))
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