program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    implicit none
    
    include 'mpif.h'

    integer, parameter :: db=4 !kind(1.0)
    integer(kind=8) :: i,j,ll,l,dumm
    integer(kind=8) :: nx,ny,step,stamp,nlinks,nsteps,ngpus
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2 
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8
    real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq1,pi2cssq2,pi2cssq0
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,dummy
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p
    real(kind=db), allocatable, dimension(:,:) :: rho,u,v,pxx,pyy,pxy
    real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8

    ! mpi calls
    integer :: ierr, idrank, nprocss, bottom, top, ny_in, ny_end, nyrank
    real(kind=db), dimension(:,:) :: buf_send_to_up,buf_send_to_down,buf_rcv_from_up,buf_rcv_from_down

       
    nlinks=8 !pari!
    tau=1.0_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
    #ifdef _OPENACC
    ngpus=acc_get_num_devices(acc_device_nvidia)
    #endif
    omega=1.0_db/tau
    write(6,*) ngpus
    !*******************************user parameters**************************
    nx=4024
    ny=4024
    nsteps=1000
    stamp=1000
    fx=0.0_db*10.0**(-7)
    fy=1.0_db*10.0**(-8)
    !*******************************mpi calls and setup*******************************!
        call MPI_INIT( ierr ) 
		!save rank of this process
		call MPI_COMM_RANK( MPI_COMM_WORLD, idrank, ierr ) 
		! save number of processes
		call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocss, ierr )
		!finalize mpi
		if(idrank.eq.0)then
			write(6,*) 'rank 0 ','processes of', nprocss, 'processes'
		endif
        ! 1d decomposition
        if(idrank<mod(ny,nprocss))then
            nyrank=ny/nprocss+1
            ny_in=idrank*nyrank + 1
        else
            nyrank=ny/nprocss
            ny_in=ny-(nprocss-idrank)*nyrank+1 
        endif
        ny_end=ny_in+nyrank-1
       !extrema for message passing 
        top = mod((idrank + 1) , nprocss)
        bottom = idrank - 1;
        if (idrank.eq.nprocss-1)then
             bottom = idrank-1 !periodicity
             top=0
        endif
        if (idrank.eq.0)then
             bottom = nprocss-1 !periodicity
             top=1
        endif
		write(6,*) 'idrank',idrank
        write(6,*) 'ny_in ',ny_in
        write(6,*) 'ny_end',ny_end
        write(6,*) 'top   ',top
        write(6,*) 'bottom',bottom

        stop
    !********************************************allocation*********************
        allocate(p(0:nlinks))
        allocate(f0(0:nx+1,ny_in-1:ny_end+1),f1(0:nx+1,ny_in-1:ny_end+1),f2(0:nx+1,ny_in-1:ny_end+1),f3(0:nx+1,ny_in-1:ny_end+1),f4(0:nx+1,ny_in-1:ny_end+1))
        allocate(f5(0:nx+1,ny_in-1:ny_end+1),f6(0:nx+1,ny_in-1:ny_end+1),f7(0:nx+1,ny_in-1:ny_end+1),f8(0:nx+1,ny_in-1:ny_end+1))
        allocate(rho(1:nx,ny_in:ny_end),u(1:nx,ny_in:ny_end),v(1:nx,ny_in:ny_end),pxx(1:nx,ny_in:ny_end),pyy(1:nx,ny_in:ny_end),pxy(1:nx,ny_in:ny_end))
        allocate(isfluid(1:nx,ny_in:ny_end)) !,omega_2d(1:nx,1:ny)) 
        allocate(buf_send_to_up(3,1:nx-2),buf_send_to_down(3,1:nx-2),buf_rcv_from_up(3,1:nx-2),buf_rcv_from_down(3,1:nx-2))
        
        !ex=(/0,1,0,-1,0,1,-1,-1,1/)
        !ey=(/0,0,1,0,-1,1,1,-1,-1/)

        p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
        1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
    !*********************************** regularized: hermite ********************
        qxx=1.0_db-cssq
        qyy=1.0_db-cssq
        qxy5_7=1.0_db
        qxy6_8=-1.0_db
        pi2cssq0=p(0)/(2.0_db*cssq**2)
        pi2cssq1=p(1)/(2.0_db*cssq**2)
        pi2cssq2=p(5)/(2.0_db*cssq**2)
    
    !*****************************************geometry************************
        isfluid=1
        isfluid(1,ny_in:ny_end)=0 !EAST
        isfluid(nx,ny_in:ny_end)=0 !WEST
        do j=ny_in,ny_end
            if(j==1 .or. j==ny)then
            isfluid(1:nx,j)=0 !SOUTH 
            isfluid(1:nx,j)=0 !NORTH
            endif
        enddo
    !*************************************initial conditions ************************    
        u=0.0_db
        v=0.0_db
        rho=1.0_db     
        f0(1:nx,ny_in:ny_end)=p(0)*rho(1:nx,ny_in:ny_end)!0.0_db
        f1(1:nx,ny_in:ny_end)=p(1)*rho(1:nx,ny_in:ny_end)
        f2(1:nx,ny_in:ny_end)=p(2)*rho(1:nx,ny_in:ny_end)
        f3(1:nx,ny_in:ny_end)=p(3)*rho(1:nx,ny_in:ny_end)
        f4(1:nx,ny_in:ny_end)=p(4)*rho(1:nx,ny_in:ny_end)
        f5(1:nx,ny_in:ny_end)=p(5)*rho(1:nx,ny_in:ny_end)
        f6(1:nx,ny_in:ny_end)=p(6)*rho(1:nx,ny_in:ny_end)
        f7(1:nx,ny_in:ny_end)=p(7)*rho(1:nx,ny_in:ny_end)
        f8(1:nx,ny_in:ny_end)=p(8)*rho(1:nx,ny_in:ny_end)
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
        !$acc kernels 
        !$acc loop collapse(2) private(uu,temp,udotc,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8) 
        do j=ny_in,ny_end
            do i=1,nx
                if(isfluid(i,j).eq.1)then
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
            enddo 
        enddo
        !***********************************collision + no slip + forcing: fused implementation*********
            !$acc loop collapse(2) private(uu,temp,udotc,feq,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8) 
            do j=ny_in,ny_end
                do i=1,nx 
                    if(isfluid(i,j).eq.1)then   
                        uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                        !oneminusuu= -uu !1.0_db - uu
                        !0
                        feq=p(0)*(rho(i,j)-uu)
                        f0(i,j)=feq + (1.0_db-omega)*pi2cssq0*(-cssq*pxx(i,j)-cssq*pyy(i,j))
                        !1
                        udotc=u(i,j)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(1)*(rho(i,j)+(temp + udotc))
                        f1(i+1,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j)-cssq*pyy(i,j)) + fx*p(1)/cssq !f1(i-1,j,nsp) + omega*(feq - f1(i-1,j,nsp)) + fx*p(1)/cssq
                        !3
                        feq=p(3)*(rho(i,j)+(temp - udotc))
                        f3(i-1,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j)-cssq*pyy(i,j))  - fx*p(3)/cssq !f3(i+1,j,nsp) + omega*(feq - f3(i+1,j,nsp)) - fx*p(3)/cssq
                        !2
                        udotc=v(i,j)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(2)*(rho(i,j)+(temp + udotc))
                        f2(i,j+1)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j))  + fy*p(2)/cssq !f2(i,j-1,nsp) + omega*(feq - f2(i,j-1,nsp)) + fy*p(2)/cssq
                        !4
                        feq=p(4)*(rho(i,j)+(temp - udotc))
                        f4(i,j-1)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j))  - fy*p(4)/cssq !f4(i,j+1,nsp) + omega*(feq - f4(i,j+1,nsp)) - fy*p(4)/cssq
                        !5
                        udotc=(u(i,j)+v(i,j))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(5)*(rho(i,j)+(temp + udotc))
                        f5(i+1,j+1)= feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + fx*p(5)/cssq + fy*p(5)/cssq!f5(i-1,j-1,nsp) + omega*(feq - f5(i-1,j-1,nsp)) + fx*p(5)/cssq + fy*p(5)/cssq 
                        !7
                        feq=p(7)*(rho(i,j)+(temp - udotc))
                        f7(i-1,j-1)=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) - fx*p(7)/cssq - fy*p(7)/cssq !f7(i+1,j+1,nsp) + omega*(feq - f7(i+1,j+1,nsp)) - fx*p(7)/cssq - fy*p(7)/cssq
                        !6
                        udotc=(-u(i,j)+v(i,j))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p(6)*(rho(i,j)+(temp + udotc))
                        f6(i-1,j+1)= feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) - fx*p(6)/cssq + fy*p(6)/cssq !f6(i+1,j-1,nsp) + omega*(feq - f6(i+1,j-1,nsp)) - fx*p(6)/cssq + fy*p(6)/cssq
                        !8
                        feq=p(8)*(rho(i,j)+(temp - udotc))
                        f8(i+1,j-1)=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) + fx*p(8)/cssq - fy*p(8)/cssq !f8(i-1,j+1,nsp) + omega*(feq - f8(i-1,j+1,nsp)) + fx*p(8)/cssq - fy*p(8)/cssq
                    endif
                enddo
            enddo
        !********************************************bcs no slip*****************************************!
            !$acc loop independent 
            do j=1,ny
                !$acc loop independent 
                do i=1,nx
                    if(isfluid(i,j).eq.0)then

                        f8(i+1,j-1)=f6(i,j)!gpc 
                        f7(i-1,j-1)=f5(i,j)!hpc

                        f6(i-1,j+1)=f8(i,j)!gpc 
                        f5(i+1,j+1)=f7(i,j)!hpc 


                        f4(i,j-1)=f2(i,j)!gpc 
                        f3(i-1,j)=f1(i,j)!hpc 

                        f2(i,j+1)=f4(i,j)!gpc 
                        f1(i+1,j)=f3(i,j)!hpc 
                    endif
                enddo
            enddo
        !$acc end kernels
        !****************************************** mpi passing*****************************************
            ! update cpu with missing pops !acc update
            !$acc update host(f2(2:nx-1,ny_end+1),f5(2:nx-1,ny_end+1),f6(2:nx-1,ny_end+1))
            !fill the buffers
            buf_send_to_up(1,1:nx-2)=f2(2:nx-1,ny_end+1)
            buf_send_to_up(2,1:nx-2)=f5(2:nx-1,ny_end+1)
            buf_send_to_up(3,1:nx-2)=f6(2:nx-1,ny_end+1)
            ! update cpu with missing pops !acc update
            !$acc update host(f2(2:nx-1,ny_in-1),f5(2:nx-1,ny_in-1),f6(2:nx-1,ny_in-1))
            buf_send_to_down(1,1:nx-2)=f2(2:nx-1,ny_in-1)
            buf_send_to_down(2,1:nx-2)=f5(2:nx-1,ny_in-1)
            buf_send_to_down(3,1:nx-2)=f6(2:nx-1,ny_in-1)
            
            !send and receives
            call mpi_sendrecv(buf_send_to_up,3*(nx-2), MPI_double_precision, top, 334, buf_rcv_from_down, 3*(nx-2), MPI_double_precision, bottom, 334,MPI_COMM_WORLD, status, ierr)
            call mpi_sendrecv(buf_send_to_down,3*(nx-2), MPI_double_precision, bottom, 564, buf_rcv_from_up, 3*(nx-2), MPI_double_precision, top, 564,MPI_COMM_WORLD, status, ierr)

            if(idrank.gt.0)then
                f2(2:nx-1,ny_in)=buf_rcv_from_down(1,1:nx-2)
                f5(2:nx-1,ny_in)=buf_rcv_from_down(2,1:nx-2)
                f6(2:nx-1,ny_in)=buf_rcv_from_down(3,1:nx-2)
                !$acc update host(f2(2:nx-1,ny_in),f5(2:nx-1,ny_in),f6(2:nx-1,ny_in))
            elseif(idrank.eq.0)then !for periodic bcs: copy in the buffer and if needed this distro can be copied as in in ny_in+1
                f2(2:nx-1,ny_in-1)=buf_rcv_from_up(1,1:nx-2)
                f5(2:nx-1,ny_in-1)=buf_rcv_from_up(2,1:nx-2)
                f6(2:nx-1,ny_in-1)=buf_rcv_from_up(3,1:nx-2)
                !$acc update host(f2(2:nx-1,ny_in-1),f5(2:nx-1,ny_in-1),f6(2:nx-1,ny_in-1))
            endif
            if(idrank.lt.nprocss-1)then
                f4(2:nx-1,ny_end)=buf_rcv_from_up(1,1:nx-2)
                f7(2:nx-1,ny_end)=buf_rcv_from_up(2,1:nx-2)
                f8(2:nx-1,ny_end)=buf_rcv_from_up(3,1:nx-2)
                !$acc update host(f4(2:nx-1,ny_end),f7(2:nx-1,ny_end),f8(2:nx-1,ny_end))
            elseif(idrank.eq.0)then !for periodic bcs: copy in the buffer and if needed this distro can be copied as in in ny_end-1
                f4(2:nx-1,ny_end+1)=buf_rcv_from_up(1,1:nx-2)
                f7(2:nx-1,ny_end+1)=buf_rcv_from_up(2,1:nx-2)
                f8(2:nx-1,ny_end+1)=buf_rcv_from_up(3,1:nx-2)
                !$acc update host(f4(2:nx-1,ny_end+1),f7(2:nx-1,ny_end+1),f8(2:nx-1,ny_end+1))
            endif
        !******************************************call periodic bcs: always after fused************************
            !$acc kernels
            !periodic along y
            f5(2:nx-1,ny_in+1)=f5(2:nx-1,ny_in-1)
            f2(2:nx-1,ny_in+1)=f2(2:nx-1,ny_in-1)
            f6(2:nx-1,ny_in+1)=f6(2:nx-1,ny_in-1)

            f8(2:nx-1,ny_end-1)=f8(2:nx-1,ny_end+1)
            f4(2:nx-1,ny_end-1)=f4(2:nx-1,ny_end+1)
            f7(2:nx-1,ny_end-1)=f7(2:nx-1,ny_end+1)
            !$acc end kernels 

    enddo 
    call cpu_time(ts2)
    !$acc update host(rho,u,v)

    !$acc end data


    !************************************************test points**********************************************!
    write(6,*) 'u=',u(nx/2,ny/2) ,'v=',v(nx/2,ny/2),'rho',rho(nx/2,ny/2) !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'u=',u(2,ny/2) ,'v=',v(2,ny/2),'rho',rho(2,ny/2)
    write(6,*) 'u=',u(1,ny/2) ,'v=',v(1,ny/2),'rho',rho(1,ny/2)
    
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  nx*ny*nsteps/10.0_db**9/ts2-ts1

    open(101, file = 'v.out', status = 'replace')
    do j=1,ny
        do i=1,nx
            write(101,*) v(i,j) 
        enddo
    enddo
    close(101) 
end program