 
 module mysubs
   
    use cudafor
    
    implicit none 
    
    integer, parameter :: db=kind(1.e0)
    
    real(kind=db),parameter :: zero=real(0.d0,kind=db)
    real(kind=db),parameter :: one=real(1.d0,kind=db)
    real(kind=db),parameter :: half=real(0.5d0,kind=db)
    real(kind=db),parameter :: halfthree=real(1.5d0,kind=db)
    real(kind=db),parameter :: two=real(2.d0,kind=db)
    real(kind=db),parameter :: three=real(3.d0,kind=db)
    real(kind=db),parameter :: five=real(5.d0,kind=db)
    real(kind=db),parameter :: ten=real(10.d0,kind=db)
    real(kind=db),parameter :: eighteen=real(18.d0,kind=db)
    real(kind=db),parameter :: thirtysix=real(36.d0,kind=db)
    
    real(kind=db),parameter :: pi_greek=real(3.141592653589793238462643383279502884d0,kind=db)
    real(kind=db),parameter :: p0 = (one/three)
    real(kind=db),parameter :: p1 = (one/eighteen)
    real(kind=db),parameter :: p2 = (one/thirtysix)
    real(kind=db),parameter :: cssq = one/three
    real(kind=db),parameter :: onecssq = three
    real(kind=db),parameter :: halfonecssq = halfthree
    real(kind=db),parameter :: p1dcssq=p1/cssq
    real(kind=db),parameter :: p2dcssq=p2/cssq
    
    real(kind=db),parameter :: pi2cssq0=p0/(two*cssq**two)
    real(kind=db),parameter :: pi2cssq1=p1/(two*cssq**two)
    real(kind=db),parameter :: pi2cssq2=p2/(two*cssq**two)

    real(kind=db),parameter :: qxx=one-cssq
    real(kind=db),parameter :: qyy=one-cssq
    real(kind=db),parameter :: qzz=one-cssq
    real(kind=db),parameter :: qxy_7_8=one
    real(kind=db),parameter :: qxy_9_10=-one
    real(kind=db),parameter :: qxz_15_16=one
    real(kind=db),parameter :: qxz_17_18=-one
    real(kind=db),parameter :: qyz_11_12=one
    real(kind=db),parameter :: qyz_13_14=-one
    ! device arrays
      integer(kind=4), allocatable,  dimension(:,:,:), device   :: isfluid_d
      integer, constant :: nx_d,ny_d,nz_d,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
      
      real(kind=db), constant :: fx_d,fy_d,fz_d,omega_d,myrho_d,myu_d,myv_d,myw_d,oneminusomega_d
      real(kind=db), allocatable, dimension(:,:,:), device  :: rho_d,u_d,v_d,w_d,pxx_d,pyy_d,pzz_d,pxy_d,pxz_d,pyz_d
      real(kind=db), allocatable, dimension(:,:,:), device :: f0_d,f1_d,f2_d,f3_d,f4_d,f5_d,f6_d,f7_d,f8_d,f9_d
      real(kind=db), allocatable, dimension(:,:,:), device :: f10_d,f11_d,f12_d,f13_d,f14_d,f15_d,f16_d,f17_d,f18_d
      real(kind=db), allocatable, dimension(:,:,:), device :: rhoprint_d
      real(kind=db), allocatable, dimension(:,:,:,:), device :: velprint_d
      type (dim3) :: dimGrid,dimBlock,dimGridx,dimGridy,dimBlock2
   
      contains
      
      subroutine abortOnLastErrorAndSync(msg, step)
      implicit none
      integer, intent(in) :: step
      character(len=*), intent(in) :: msg
      integer :: istat0
    
      istat0 = cudaGetLastError()
    
      if (istat0/=0) then
        write(*,*) 'status after ',msg,':', cudaGetErrorString(istat0)
        write(*,*) 'Exiting at step:', step
        stop
      endif
    
      return
    
      end subroutine  abortOnLastErrorAndSync
  
      attributes(global) subroutine setup_pops()
      
            integer :: i,j,k
          
            
            i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
            j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
            k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
            
            !write(*,*)i,j,p_d(0)*myrho_d
            
            f0_d(i,j,k)=p0*myrho_d
            f1_d(i,j,k)=p1*myrho_d
            f2_d(i,j,k)=p1*myrho_d
            f3_d(i,j,k)=p1*myrho_d
            f4_d(i,j,k)=p1*myrho_d
            f5_d(i,j,k)=p1*myrho_d
            f6_d(i,j,k)=p1*myrho_d
            f7_d(i,j,k)=p2*myrho_d
            f8_d(i,j,k)=p2*myrho_d
            f9_d(i,j,k)=p2*myrho_d
            f10_d(i,j,k)=p2*myrho_d
            f11_d(i,j,k)=p2*myrho_d
            f12_d(i,j,k)=p2*myrho_d
            f13_d(i,j,k)=p2*myrho_d
            f14_d(i,j,k)=p2*myrho_d
            f15_d(i,j,k)=p2*myrho_d
            f16_d(i,j,k)=p2*myrho_d
            f17_d(i,j,k)=p2*myrho_d
            f18_d(i,j,k)=p2*myrho_d
     

      end subroutine setup_pops
  
      attributes(global) subroutine moments()
          
            integer :: i,j,k
            real(kind=db) :: uu,udotc,temp
            real(kind=db) ::fneq1
              
            i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
            j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
            k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
              
            if(isfluid_d(i,j,k).ne.1)return

            pxx_d(i,j,k)=zero
            pyy_d(i,j,k)=zero
            pzz_d(i,j,k)=zero
            pxy_d(i,j,k)=zero
            pxz_d(i,j,k)=zero
            pyz_d(i,j,k)=zero  

            rho_d(i,j,k) = f0_d(i,j,k)+f1_d(i,j,k)+f2_d(i,j,k)+f3_d(i,j,k)+f4_d(i,j,k)+f5_d(i,j,k) &
            +f6_d(i,j,k)+f7_d(i,j,k)+f8_d(i,j,k)+f9_d(i,j,k)+f10_d(i,j,k)+f11_d(i,j,k) &
            +f12_d(i,j,k)+f13_d(i,j,k)+f14_d(i,j,k)+f15_d(i,j,k)+f16_d(i,j,k)+f17_d(i,j,k) &
            +f18_d(i,j,k)

            u_d(i,j,k) = (f1_d(i,j,k)+f7_d(i,j,k)+f9_d(i,j,k)+f15_d(i,j,k)+f18_d(i,j,k)) &
              -(f2_d(i,j,k)+f8_d(i,j,k)+f10_d(i,j,k)+f16_d(i,j,k)+f17_d(i,j,k)) 
            
            v_d(i,j,k) = (f3_d(i,j,k)+f7_d(i,j,k)+f10_d(i,j,k)+f11_d(i,j,k)+f13_d(i,j,k)) &
              -(f4_d(i,j,k)+f8_d(i,j,k)+f9_d(i,j,k)+f12_d(i,j,k)+f14_d(i,j,k))

            w_d(i,j,k) = (f5_d(i,j,k)+f11_d(i,j,k)+f14_d(i,j,k)+f15_d(i,j,k)+f17_d(i,j,k)) &
              -(f6_d(i,j,k)+f12_d(i,j,k)+f13_d(i,j,k)+f16_d(i,j,k)+f18_d(i,j,k))
          
            uu=halfonecssq*(u_d(i,j,k)*u_d(i,j,k) + v_d(i,j,k)*v_d(i,j,k) + w_d(i,j,k)*w_d(i,j,k))
            !1-2
            udotc=u_d(i,j,k)*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pxx_d(i,j,k)=pxx_d(i,j,k)+f1_d(i,j,k)-p1*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pxx_d(i,j,k)=pxx_d(i,j,k)+f2_d(i,j,k)-p1*(rho_d(i,j,k)+(temp - udotc))
            !3-4
            udotc=v_d(i,j,k)*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pyy_d(i,j,k)=pyy_d(i,j,k)+f3_d(i,j,k)-p1*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pyy_d(i,j,k)=pyy_d(i,j,k)+f4_d(i,j,k)-p1*(rho_d(i,j,k)+(temp - udotc))
            !5-6
            udotc=w_d(i,j,k)*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f5_d(i,j,k)-p1*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f6_d(i,j,k)-p1*(rho_d(i,j,k)+(temp - udotc))
            !7-8
            udotc=(u_d(i,j,k)+v_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pxx_d(i,j,k)=pxx_d(i,j,k)+f7_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pyy_d(i,j,k)=pyy_d(i,j,k)+f7_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pxy_d(i,j,k)=pxy_d(i,j,k)+f7_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pxx_d(i,j,k)=pxx_d(i,j,k)+f8_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pyy_d(i,j,k)=pyy_d(i,j,k)+f8_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pxy_d(i,j,k)=pxy_d(i,j,k)+f8_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            !10-9
            udotc=(-u_d(i,j,k)+v_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pxx_d(i,j,k)=pxx_d(i,j,k)+f10_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pyy_d(i,j,k)=pyy_d(i,j,k)+f10_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pxy_d(i,j,k)=pxy_d(i,j,k)-f10_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pxx_d(i,j,k)=pxx_d(i,j,k)+f9_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pyy_d(i,j,k)=pyy_d(i,j,k)+f9_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pxy_d(i,j,k)=pxy_d(i,j,k)-f9_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            !11-12
            udotc=(v_d(i,j,k)+w_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f11_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pyy_d(i,j,k)=pyy_d(i,j,k)+f11_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pyz_d(i,j,k)=pyz_d(i,j,k)+f11_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f12_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pyy_d(i,j,k)=pyy_d(i,j,k)+f12_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pyz_d(i,j,k)=pyz_d(i,j,k)+f12_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            !13-14
            udotc=(v_d(i,j,k)-w_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f13_d(i,j,k) - p2*(rho_d(i,j,k)+(temp + udotc))
            pyy_d(i,j,k)=pyy_d(i,j,k)+f13_d(i,j,k) - p2*(rho_d(i,j,k)+(temp + udotc))
            pyz_d(i,j,k)=pyz_d(i,j,k)-f13_d(i,j,k) - p2*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f14_d(i,j,k) - p2*(rho_d(i,j,k)+(temp - udotc))
            pyy_d(i,j,k)=pyy_d(i,j,k)+f14_d(i,j,k) - p2*(rho_d(i,j,k)+(temp - udotc))
            pyz_d(i,j,k)=pyz_d(i,j,k)-f14_d(i,j,k) - p2*(rho_d(i,j,k)+(temp - udotc))
            !15-16
            udotc=(u_d(i,j,k)+w_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f15_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pxx_d(i,j,k)=pxx_d(i,j,k)+f15_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pxz_d(i,j,k)=pxz_d(i,j,k)+f15_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f16_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pxx_d(i,j,k)=pxx_d(i,j,k)+f16_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pxz_d(i,j,k)=pxz_d(i,j,k)+f16_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            !17-18
            udotc=(-u_d(i,j,k)+w_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f17_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pxx_d(i,j,k)=pxx_d(i,j,k)+f17_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            pxz_d(i,j,k)=pxz_d(i,j,k)-f17_d(i,j,k)-p2*(rho_d(i,j,k)+(temp + udotc))
            !fneq1=
            pzz_d(i,j,k)=pzz_d(i,j,k)+f18_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pxx_d(i,j,k)=pxx_d(i,j,k)+f18_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))
            pxz_d(i,j,k)=pxz_d(i,j,k)-f18_d(i,j,k)-p2*(rho_d(i,j,k)+(temp - udotc))

      end subroutine moments
  
      attributes(global) subroutine streamcoll()
          
            integer :: i,j,k
            real(kind=db) ::uu,udotc,temp,feq
            
              
            i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
            j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
            k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
              
            if(isfluid_d(i,j,k).ne.1)return
              
              
            uu=halfonecssq*(u_d(i,j,k)*u_d(i,j,k) + v_d(i,j,k)*v_d(i,j,k) + w_d(i,j,k)*w_d(i,j,k))
            !0
            feq=p0*(rho_d(i,j,k)-uu)
            f0_d(i,j,k)=feq + oneminusomega_d*pi2cssq0*(-cssq*(pyy_d(i,j,k)+pxx_d(i,j,k)+pzz_d(i,j,k)))
            
            !1
            udotc=u_d(i,j,k)*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p1*(rho_d(i,j,k)+(temp + udotc))
            f1_d(i+1,j,k)=feq + oneminusomega_d*pi2cssq1*(qxx*pxx_d(i,j,k)-cssq*(pyy_d(i,j,k)+pzz_d(i,j,k))) + fx_d*p1dcssq
            
            !2
            feq=p1*(rho_d(i,j,k)+(temp - udotc))
            f2_d(i-1,j,k)=feq + oneminusomega_d*pi2cssq1*(qxx*pxx_d(i,j,k)-cssq*(pyy_d(i,j,k)+pzz_d(i,j,k))) - fx_d*p1dcssq
            
            !3
            udotc=v_d(i,j,k)*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p1*(rho_d(i,j,k)+(temp + udotc))
            f3_d(i,j+1,k)=feq+oneminusomega_d*pi2cssq1*(qyy*pyy_d(i,j,k)-cssq*(pxx_d(i,j,k)+pzz_d(i,j,k))) + fy_d*p1dcssq
            
            !4
            feq=p1*(rho_d(i,j,k)+(temp - udotc))
            f4_d(i,j-1,k)=feq+oneminusomega_d*pi2cssq1*(qyy*pyy_d(i,j,k)-cssq*(pxx_d(i,j,k)+pzz_d(i,j,k))) - fy_d*p1dcssq
            
            !7
            udotc=(u_d(i,j,k)+v_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(rho_d(i,j,k)+(temp + udotc))
            f7_d(i+1,j+1,k)=feq + oneminusomega_d*pi2cssq2*(qxx*pxx_d(i,j,k)+qyy*pyy_d(i,j,k)-cssq*pzz_d(i,j,k)+two*qxy_7_8*pxy_d(i,j,k)) + (fx_d+fy_d)*p2dcssq 
            
            !8
            feq=p2*(rho_d(i,j,k)+(temp - udotc))
            f8_d(i-1,j-1,k)=feq + oneminusomega_d*pi2cssq2*(qxx*pxx_d(i,j,k)+qyy*pyy_d(i,j,k)-cssq*pzz_d(i,j,k)+two*qxy_7_8*pxy_d(i,j,k)) - (fx_d+fy_d)*p2dcssq
            
            !10
            udotc=(-u_d(i,j,k)+v_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(rho_d(i,j,k)+(temp + udotc))
            f10_d(i-1,j+1,k)=feq+oneminusomega_d*pi2cssq2*(qxx*pxx_d(i,j,k)+qyy*pyy_d(i,j,k)-cssq*pzz_d(i,j,k)+two*qxy_9_10*pxy_d(i,j,k)) +(fy_d-fx_d)*p2dcssq
            
            !9
            feq=p2*(rho_d(i,j,k)+(temp - udotc))
            f9_d(i+1,j-1,k)=feq+oneminusomega_d*pi2cssq2*(qxx*pxx_d(i,j,k)+qyy*pyy_d(i,j,k)-cssq*pzz_d(i,j,k)+two*qxy_9_10*pxy_d(i,j,k)) + (fx_d-fy_d)*p2dcssq

            !5
            udotc=w_d(i,j,k)*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p1*(rho_d(i,j,k)+(temp + udotc))
            f5_d(i,j,k+1)=feq+oneminusomega_d*pi2cssq1*(qzz*pzz_d(i,j,k)-cssq*(pxx_d(i,j,k)+pyy_d(i,j,k))) + fz_d*p1dcssq
            
            !6
            feq=p1*(rho_d(i,j,k)+(temp - udotc))
            f6_d(i,j,k-1)=feq+oneminusomega_d*pi2cssq1*(qzz*pzz_d(i,j,k)-cssq*(pxx_d(i,j,k)+pyy_d(i,j,k))) - fz_d*p1dcssq

            !15
            udotc=(u_d(i,j,k)+w_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(rho_d(i,j,k)+(temp + udotc))
            f15_d(i+1,j,k+1)=feq+oneminusomega_d*pi2cssq2*(qxx*pxx_d(i,j,k)+qzz*pzz_d(i,j,k)-cssq*pyy_d(i,j,k)+two*qxz_15_16*pxz_d(i,j,k)) + (fx_d+fz_d)*p2dcssq 
            
            !16
            feq=p2*(rho_d(i,j,k)+(temp - udotc))
            f16_d(i-1,j,k-1)=feq+ oneminusomega_d*pi2cssq2*(qxx*pxx_d(i,j,k)+qzz*pzz_d(i,j,k)-cssq*pyy_d(i,j,k)+two*qxz_15_16*pxz_d(i,j,k)) - (fx_d+fz_d)*p2dcssq

            !17
            udotc=(-u_d(i,j,k)+w_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(rho_d(i,j,k)+(temp + udotc))
            f17_d(i-1,j,k+1)=feq+oneminusomega_d*pi2cssq2*(qxx*pxx_d(i,j,k)+qzz*pzz_d(i,j,k)-cssq*pyy_d(i,j,k)+two*qxz_17_18*pxz_d(i,j,k)) +(fz_d-fx_d)*p2dcssq
            
            !18
            feq=p2*(rho_d(i,j,k)+(temp - udotc))
            f18_d(i+1,j,k-1)=feq+oneminusomega_d*pi2cssq2*(qxx*pxx_d(i,j,k)+qzz*pzz_d(i,j,k)-cssq*pyy_d(i,j,k)+two*qxz_17_18*pxz_d(i,j,k)) + (fx_d-fz_d)*p2dcssq

            !11
            udotc=(v_d(i,j,k)+w_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(rho_d(i,j,k)+(temp + udotc))
            f11_d(i,j+1,k+1)=feq+oneminusomega_d*pi2cssq2*(qyy*pyy_d(i,j,k)+qzz*pzz_d(i,j,k)-cssq*pxx_d(i,j,k)+two*qyz_11_12*pyz_d(i,j,k))+(fy_d+fz_d)*p2dcssq
            
            !12
            feq=p2*(rho_d(i,j,k)+(temp - udotc))
            f12_d(i,j-1,k-1)=feq+oneminusomega_d*pi2cssq2*(qyy*pyy_d(i,j,k)+qzz*pzz_d(i,j,k)-cssq*pxx_d(i,j,k)+two*qyz_11_12*pyz_d(i,j,k)) - (fy_d+fz_d)*p2dcssq

            !13
            udotc=(v_d(i,j,k)-w_d(i,j,k))*halfonecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(rho_d(i,j,k)+(temp + udotc))
            f13_d(i,j+1,k-1)=feq+oneminusomega_d*pi2cssq2*(qyy*pyy_d(i,j,k)+qzz*pzz_d(i,j,k)-cssq*pxx_d(i,j,k)+two*qyz_13_14*pyz_d(i,j,k)) + (fy_d-fz_d)*p2dcssq
            
            !14
            feq=p2*(rho_d(i,j,k)+(temp - udotc))
            f14_d(i,j-1,k+1)=feq+oneminusomega_d*pi2cssq2*(qyy*pyy_d(i,j,k)+qzz*pzz_d(i,j,k)-cssq*pxx_d(i,j,k)+two*qyz_13_14*pyz_d(i,j,k)) + (fz_d-fy_d)*p2dcssq

      end subroutine streamcoll
  
      attributes(global) subroutine bcs_no_slip()
          
          integer :: i,j,k
        
          
          i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
          j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
          k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
          
          !write(*,*)i,j,p_d(0)*myrho_d
          if(isfluid_d(i,j,k).ne.0)return
          

          f18_d(i+1,j,k-1)=f17_d(i,j,k) !gpc 
          f17_d(i-1,j,k+1)=f18_d(i,j,k) !hpc

          f16_d(i-1,j,k-1)=f15_d(i,j,k) !gpc 
          f15_d(i+1,j,k+1)=f16_d(i,j,k) !hpc

          f14_d(i,j-1,k+1)=f13_d(i,j,k)!gpc 
          f13_d(i,j+1,k-1)=f14_d(i,j,k)!hpc
        
          f12_d(i,j-1,k-1)=f11_d(i,j,k)!gpc 
          f11_d(i,j+1,k+1)=f12_d(i,j,k)!hpc

          f10_d(i-1,j+1,k)=f9_d(i,j,k)!gpc 
          f9_d(i+1,j-1,k)=f10_d(i,j,k)!hpc

          f8_d(i-1,j-1,k)=f7_d(i,j,k)!gpc 
          f7_d(i+1,j+1,k)=f8_d(i,j,k)!hpc

          f6_d(i,j,k-1)=f5_d(i,j,k)!gpc 
          f5_d(i,j,k+1)=f6_d(i,j,k)!hpc 

          f4_d(i,j-1,k)=f3_d(i,j,k)!gpc 
          f3_d(i,j+1,k)=f4_d(i,j,k)!hpc 

          f2_d(i-1,j,k)=f1_d(i,j,k)!gpc 
          f1_d(i+1,j,k)=f2_d(i,j,k)!hpc 

        

      end subroutine bcs_no_slip
  
      attributes(global) subroutine pbc_edge_x()
          
          integer :: j,k
        
          
          j = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
          k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
          
          if (j>ny_d .or. k>nz_d)return
          
          if(j>2 .and. j<ny_d-1 .and. k>2 .and. k<nz_d-1)then
            
            f1_d(2,j,k)=f1_d(nx_d,j,k)
            f7_d(2,j,k)=f7_d(nx_d,j,k)
            f9_d(2,j,k)=f9_d(nx_d,j,k)
            f15_d(2,j,k)=f15_d(nx_d,j,k)
            f18_d(2,j,k)=f18_d(nx_d,j,k)
            
            f2_d(nx_d-1,j,k)=f2_d(1,j,k)
            f8_d(nx_d-1,j,k)=f8_d(1,j,k)
            f10_d(nx_d-1,j,k)=f10_d(1,j,k)
            f16_d(nx_d-1,j,k)=f16_d(1,j,k)
            f17_d(nx_d-1,j,k)=f17_d(1,j,k)
          else
          
            if(j==2)then
              if(k==2)then
                f1_d(2,j,k)=f1_d(nx_d,j,k)
                f9_d(2,j,k)=f9_d(nx_d,j,k)
                f18_d(2,j,k)=f18_d(nx_d,j,k)
                
                f2_d(nx_d-1,j,k)=f2_d(1,j,k)
                f8_d(nx_d-1,j,k)=f8_d(1,j,k)
                f16_d(nx_d-1,j,k)=f16_d(1,j,k)
              elseif(k==nz_d-1)then
                f1_d(2,j,k)=f1_d(nx_d,j,k)
                f9_d(2,j,k)=f9_d(nx_d,j,k)
                f15_d(2,j,k)=f15_d(nx_d,j,k)
                
                f2_d(nx_d-1,j,k)=f2_d(1,j,k)
                f8_d(nx_d-1,j,k)=f8_d(1,j,k)
                f17_d(nx_d-1,j,k)=f17_d(1,j,k)
              else
                f1_d(2,j,k)=f1_d(nx_d,j,k)
                f9_d(2,j,k)=f9_d(nx_d,j,k)
                f15_d(2,j,k)=f15_d(nx_d,j,k)
                f18_d(2,j,k)=f18_d(nx_d,j,k)
                
                f2_d(nx_d-1,j,k)=f2_d(1,j,k)
                f8_d(nx_d-1,j,k)=f8_d(1,j,k)
                f16_d(nx_d-1,j,k)=f16_d(1,j,k)
                f17_d(nx_d-1,j,k)=f17_d(1,j,k)
              endif
            elseif(j==ny_d-1)then
              if(k==2)then
                f1_d(2,j,k)=f1_d(nx_d,j,k)
                f7_d(2,j,k)=f7_d(nx_d,j,k)
                f18_d(2,j,k)=f18_d(nx_d,j,k)
                
                f2_d(nx_d-1,j,k)=f2_d(1,j,k)
                f10_d(nx_d-1,j,k)=f10_d(1,j,k)
                f16_d(nx_d-1,j,k)=f16_d(1,j,k)
              elseif(k==nz_d-1)then
                f1_d(2,j,k)=f1_d(nx_d,j,k)
                f7_d(2,j,k)=f7_d(nx_d,j,k)
                f15_d(2,j,k)=f15_d(nx_d,j,k)
                
                f2_d(nx_d-1,j,k)=f2_d(1,j,k)
                f10_d(nx_d-1,j,k)=f10_d(1,j,k)
                f17_d(nx_d-1,j,k)=f17_d(1,j,k)
              else
                f1_d(2,j,k)=f1_d(nx_d,j,k)
                f7_d(2,j,k)=f7_d(nx_d,j,k)
                f15_d(2,j,k)=f15_d(nx_d,j,k)
                f18_d(2,j,k)=f18_d(nx_d,j,k)
                
                f2_d(nx_d-1,j,k)=f2_d(1,j,k)
                f10_d(nx_d-1,j,k)=f10_d(1,j,k)
                f16_d(nx_d-1,j,k)=f16_d(1,j,k)
                f17_d(nx_d-1,j,k)=f17_d(1,j,k)
              endif
            endif
          endif 
              
      end subroutine pbc_edge_x
  
      attributes(global) subroutine pbc_edge_y()
          
          integer :: i,k
        
          
          i = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
          k = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y

          if (i>nx_d .or. k>nz_d)return
          
          if(i>2 .and. i<nx_d-1 .and. k>2 .and. k<nz_d-1)then
          
            f3_d(i,2,k)=f3_d(i,ny_d,k)
            f7_d(i,2,k)=f7_d(i,ny_d,k)
            f10_d(i,2,k)=f10_d(i,ny_d,k)
            f11_d(i,2,k)=f11_d(i,ny_d,k)
            f13_d(i,2,k)=f13_d(i,ny_d,k)
            
            f4_d(i,ny_d-1,k)=f4_d(i,1,k)
            f8_d(i,ny_d-1,k)=f8_d(i,1,k)
            f9_d(i,ny_d-1,k)=f9_d(i,1,k)
            f12_d(i,ny_d-1,k)=f12_d(i,1,k)
            f14_d(i,ny_d-1,k)=f14_d(i,1,k)
          else
          
            if(i==2)then
              if(k==2)then
                f3_d(i,2,k)=f3_d(i,ny_d,k)
                f10_d(i,2,k)=f10_d(i,ny_d,k)
                f13_d(i,2,k)=f13_d(i,ny_d,k)
                
                f4_d(i,ny_d-1,k)=f4_d(i,1,k)
                f8_d(i,ny_d-1,k)=f8_d(i,1,k)
                f12_d(i,ny_d-1,k)=f12_d(i,1,k)
              elseif(k==nz_d-1)then
                f3_d(i,2,k)=f3_d(i,ny_d,k)
                f10_d(i,2,k)=f10_d(i,ny_d,k)
                f11_d(i,2,k)=f11_d(i,ny_d,k)

                f4_d(i,ny_d-1,k)=f4_d(i,1,k)
                f8_d(i,ny_d-1,k)=f8_d(i,1,k)
                f14_d(i,ny_d-1,k)=f14_d(i,1,k)
              else
                f3_d(i,2,k)=f3_d(i,ny_d,k)
                f10_d(i,2,k)=f10_d(i,ny_d,k)
                f11_d(i,2,k)=f11_d(i,ny_d,k)
                f13_d(i,2,k)=f13_d(i,ny_d,k)

                f4_d(i,ny_d-1,k)=f4_d(i,1,k)
                f8_d(i,ny_d-1,k)=f8_d(i,1,k)
                f12_d(i,ny_d-1,k)=f12_d(i,1,k)
                f14_d(i,ny_d-1,k)=f14_d(i,1,k)
              endif
            elseif(i==nx_d-1)then
              if(k==2)then
                f3_d(i,2,k)=f3_d(i,ny_d,k)
                f7_d(i,2,k)=f7_d(i,ny_d,k)
                f13_d(i,2,k)=f13_d(i,ny_d,k)

                f4_d(i,ny_d-1,k)=f4_d(i,1,k)
                f9_d(i,ny_d-1,k)=f9_d(i,1,k)
                f12_d(i,ny_d-1,k)=f12_d(i,1,k)
              elseif(k==nz_d-1)then 
                f3_d(i,2,k)=f3_d(i,ny_d,k)
                f7_d(i,2,k)=f7_d(i,ny_d,k)
                f11_d(i,2,k)=f11_d(i,ny_d,k)

                f4_d(i,ny_d-1,k)=f4_d(i,1,k)
                f9_d(i,ny_d-1,k)=f9_d(i,1,k)
                f14_d(i,ny_d-1,k)=f14_d(i,1,k)
              else
                f3_d(i,2,k)=f3_d(i,ny_d,k)
                f7_d(i,2,k)=f7_d(i,ny_d,k)
                f11_d(i,2,k)=f11_d(i,ny_d,k)
                f13_d(i,2,k)=f13_d(i,ny_d,k)

                f4_d(i,ny_d-1,k)=f4_d(i,1,k)
                f9_d(i,ny_d-1,k)=f9_d(i,1,k)
                f12_d(i,ny_d-1,k)=f12_d(i,1,k)
                f14_d(i,ny_d-1,k)=f14_d(i,1,k)
              endif
            endif
          endif
        

      end subroutine pbc_edge_y
  

  
  
      attributes(global) subroutine store_print()
          
          integer :: i,j,k
        
          
          i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
          j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
          k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
          
          
          !write(*,*)i,j,p_d(0)*myrho_d
          if(isfluid_d(i,j,k).eq.1)then
            rhoprint_d(i,j,k)=rho_d(i,j,k)
            velprint_d(1,i,j,k)=u_d(i,j,k)
            velprint_d(2,i,j,k)=v_d(i,j,k)
            velprint_d(3,i,j,k)=w_d(i,j,k)
            
          else
          
            rhoprint_d(i,j,k)=zero
            velprint_d(1,i,j,k)=zero
            velprint_d(2,i,j,k)=zero
            velprint_d(3,i,j,k)=zero
          
          endif
          
          return

      end subroutine store_print

 end module mysubs
 
 module prints
  
  use mysubs
  
  implicit none
  
	integer, parameter :: mxln=120
	character(len=8), allocatable, dimension(:) :: namevarvtk
	character(len=500), allocatable, dimension(:) :: headervtk
	character(len=30), allocatable, dimension(:) :: footervtk
	integer, allocatable, dimension(:) :: ndimvtk
	integer, allocatable, dimension(:) :: vtkoffset
	integer, allocatable, dimension(:) :: ndatavtk
	integer, allocatable, dimension(:) :: nheadervtk
	integer :: nfilevtk
	integer, allocatable, dimension(:) :: varlistvtk
	character :: delimiter
	character(len=*), parameter :: filenamevtk='out'
	
	real(kind=4), allocatable, dimension(:,:,:) :: rhoprint
	real(kind=4), allocatable, dimension(:,:,:,:) :: velprint
	logical :: lelittle
	character(len=mxln) :: dir_out
	character(len=mxln) :: extentvtk
	character(len=mxln) :: sevt1,sevt2
	character(len=1), allocatable, dimension(:) :: head1,head2
	
  
  contains

 subroutine header_vtk(nx,ny,nz,mystring500,namevar,extent,ncomps,iinisub,iend,myoffset, &
   new_myoffset,indent)
  
  implicit none
  
  integer, intent(in) :: nx,ny,nz
  character(len=8),intent(in) :: namevar
  character(len=120),intent(in) :: extent
  integer, intent(in) :: ncomps,iinisub,myoffset
  integer, intent(out) :: iend,new_myoffset
  integer, intent(inout) :: indent
  
  !namevar='density1'
  
  character(len=500), intent(out) :: mystring500
  ! End-character for binary-record finalize.
  character(1), parameter:: end_rec = char(10) 
  character(1) :: string1
  character(len=*),parameter :: topology='ImageData' 
  integer :: ioffset,nele,bytechar,byteint,byter4,byter8,iini
  
  iini=iinisub
  bytechar=kind(end_rec)
  byteint=kind(iini)
  byter4  = 4
  byter8  = 8
  
  mystring500=repeat(' ',500)
  
  iend=iini
  
  iini=iend+1
  nele=22
  iend=iend+nele
  mystring500(iini:iend)='<?xml version="1.0"?>'//end_rec
  
  new_myoffset=myoffset
  new_myoffset = new_myoffset + nele * bytechar
 
  
  iini=iend+1
  nele=67
  iend=iend+nele
  if(lelittle)then  
    mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
     '" version="0.1" byte_order="LittleEndian">'//end_rec
  else
    mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
     '" version="0.1" byte_order="BigEndian">   '//end_rec
  endif
  
  new_myoffset = new_myoffset + 67 * bytechar
 
  
  indent = indent + 2
  iini=iend+1
  nele=70
  iend=iend+nele
  mystring500(iini:iend) = repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="'//&
                 trim(extent)//'">'//end_rec
  

  new_myoffset = new_myoffset + 70 * bytechar
 
  
  indent = indent + 2
  iini=iend+1
  nele=63
  iend=iend+nele
  mystring500(iini:iend) = repeat(' ',indent)//'<Piece Extent="'//trim(extent)//'">'//end_rec
  
  new_myoffset = new_myoffset + 63 * bytechar
 
  
! initializing offset pointer
  ioffset = 0 
  
  indent = indent + 2
  iini=iend+1
  nele=18
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'<PointData>'//end_rec
  
  new_myoffset = new_myoffset + 18 * bytechar
  
  indent = indent + 2
  iini=iend+1
  nele=115
  iend=iend+nele
  
  if(ncomps/=1 .and. ncomps/=3)then
    write(6,'(a)')'ERROR in header_vtk'
    stop
  endif
  write(string1,'(i1)')ncomps
  mystring500(iini:iend)=repeat(' ',indent)//'<DataArray type="Float32" Name="'// &
   namevar//'" NumberOfComponents="'//string1// '" '//&
   'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec
  
  new_myoffset = new_myoffset + 115 * bytechar
 
  
  indent = indent - 2
  iini=iend+1
  nele=19
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'</PointData>'//end_rec
  
  new_myoffset = new_myoffset + 19 * bytechar
  
  
  indent = indent - 2
  iini=iend+1
  nele=13
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'</Piece>'//end_rec
  
  
  new_myoffset = new_myoffset + 13 * bytechar
 
  
  indent = indent - 2
  iini=iend+1
  nele=15
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec

  new_myoffset = new_myoffset + 15 * bytechar
 

  iini=iend+1
  nele=32
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec
  
  new_myoffset = new_myoffset + 32 * bytechar
  
  iini=iend+1
  nele=1
  iend=iend+nele
  mystring500(iini:iend)='_'
  
  new_myoffset = new_myoffset + 1 * bytechar
  
  return
  
 end subroutine header_vtk
 
 subroutine footer_vtk(nx,ny,nz,mystring30,iinisub,iend,myoffset, &
  new_myoffset,indent)
 
  implicit none
  
  integer, intent(in) :: nx,ny,nz
  integer, intent(in) :: iinisub,myoffset
  integer, intent(out) :: iend,new_myoffset
  integer, intent(inout) :: indent
  
  
  character(len=30), intent(out) :: mystring30
  ! End-character for binary-record finalize.
  character(1), parameter:: end_rec = char(10) 
  character(1) :: string1
  character(len=*),parameter :: topology='ImageData' 
  integer :: ioffset,nele,bytechar,byteint,byter4,byter8,iini
  
  iini=iinisub
  bytechar=kind(end_rec)
  byteint=kind(iini)
  byter4  = 4
  byter8  = 8
  
  mystring30=repeat(' ',30)
  
  iend=iini
  
  iini=iend+1
  nele=1
  iend=iend+nele
  mystring30(iini:iend)=end_rec
  
  new_myoffset = myoffset
  new_myoffset = new_myoffset + 1 * bytechar
 
  
  
  iini=iend+1
  nele=18
  iend=iend+nele
  mystring30(iini:iend)=repeat(' ',indent)//'</AppendedData>'//end_rec
  
  new_myoffset = new_myoffset + 18 * bytechar
  
  iini=iend+1
  nele=11
  iend=iend+nele
  mystring30(iini:iend)='</VTKFile>'//end_rec
  
  if(iend/=30)then
     write(6,'(a)')'ERROR in footer_vtk'
    stop
  endif
  
  return
  
 end subroutine footer_vtk
  
 subroutine test_little_endian(ltest)
 
!***********************************************************************
!     
!     LBsoft subroutine for checking if the computing architecture
!     is working in little-endian or big-endian
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none 
  integer, parameter :: ik1 = selected_int_kind(2) 
  integer, parameter :: ik4 = selected_int_kind(9) 
   
  logical, intent(out) :: ltest
   
  if(btest(transfer(int((/1,0,0,0/),ik1),1_ik4),0)) then 
    !it is little endian
    ltest=.true.
  else 
    !it is big endian
    ltest=.false.
  end if 
   
  return
   
 end subroutine test_little_endian 
 
 subroutine init_output(nx,ny,nz,ncomp,lvtk)
 
!***********************************************************************
!     
!     LBsoft subroutine for creating the folders containing the files
!     in image VTK legacy binary format in parallel IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  
  implicit none
  
  integer, intent(in) :: nx,ny,nz,ncomp
  logical, intent(in) :: lvtk
  character(len=255) :: path,makedirectory
  logical :: lexist
  
  integer :: i,j,k,nn,indent,myoffset,new_myoffset,iend
  integer, parameter :: byter4=4
  integer, parameter :: byteint=4
  integer, allocatable :: printlistvtk(:)
  integer, parameter :: ioxyz=54
  character(len=*), parameter :: filexyz='isfluid.xyz'
  character(len=120) :: mystring120
  
  call test_little_endian(lelittle)
  
  sevt1=repeat(' ',mxln)
  sevt2=repeat(' ',mxln)
  
  path = repeat(' ',255)
  call getcwd(path)
  
  !call get_environment_variable('DELIMITER',delimiter)
  path = trim(path)
  delimiter = path(1:1)
  if (delimiter==' ') delimiter='/'


  
  makedirectory=repeat(' ',255)
  makedirectory = 'output'//delimiter
  dir_out=trim(makedirectory)
#ifdef _INTEL
  inquire(directory=trim(makedirectory),exist=lexist)
#else
  inquire(file=trim(makedirectory),exist=lexist)
#endif
  
  if(.not. lexist)then
    makedirectory=repeat(' ',255)
    makedirectory = 'mkdir output'
    call system(makedirectory)
  endif
  mystring120=repeat(' ',120)
  
  
  makedirectory=repeat(' ',255)
  makedirectory=trim(path)//delimiter//'output'//delimiter
  
  extentvtk =  space_fmtnumb(1) // ' ' // space_fmtnumb(nx) // ' ' &
        // space_fmtnumb(1) // ' ' // space_fmtnumb(ny) // ' ' &
        // space_fmtnumb(1) // ' ' // space_fmtnumb(nz)
  
  if(ncomp==1)then
    nfilevtk=2
  elseif(ncomp==2)then
    nfilevtk=3
  endif
  
  allocate(printlistvtk(nfilevtk))
  do i=1,nfilevtk
    printlistvtk(i)=i
  enddo
  
  allocate(varlistvtk(nfilevtk))
  allocate(namevarvtk(nfilevtk))
  allocate(ndimvtk(nfilevtk))
  allocate(headervtk(nfilevtk))
  allocate(footervtk(nfilevtk))
  allocate(nheadervtk(nfilevtk))
  allocate(vtkoffset(nfilevtk))
  allocate(ndatavtk(nfilevtk))
  varlistvtk(1:nfilevtk)=printlistvtk(1:nfilevtk)
  
  if(ncomp==1)then
    do i=1,nfilevtk
      select case(printlistvtk(i))
      case(1)
        namevarvtk(i)='rho     '
        ndimvtk(i)=1
      case(2)
        namevarvtk(i)='vel     '
        ndimvtk(i)=3
      case default
        write(6,'(a)')'ERROR in init_output'
        stop
      end select
    enddo
  elseif(ncomp==2)then
    do i=1,nfilevtk
      select case(printlistvtk(i))
      case(1)
        namevarvtk(i)='rho1    '
        ndimvtk(i)=1
      case(2)
        namevarvtk(i)='rho2    '
        ndimvtk(i)=1
      case(3)
        namevarvtk(i)='vel     '
        ndimvtk(i)=3
      case default
        write(6,'(a)')'ERROR in init_output'
        stop
      end select
    enddo
  endif
  nn=nx*ny*nz
  
  do i=1,nfilevtk
    myoffset=0
    indent=0
    call header_vtk(nx,ny,nz,headervtk(i),namevarvtk(i),extentvtk,ndimvtk(i),0,iend,myoffset, &
    new_myoffset,indent)
    vtkoffset(i)=new_myoffset
    myoffset=new_myoffset+byteint+ndimvtk(i)*nn*byter4
    ndatavtk(i)=ndimvtk(i)*nn*byter4
    nheadervtk(i)=iend
    call footer_vtk(nx,ny,nz,footervtk(i),0,iend,myoffset, &
     new_myoffset,indent)
  enddo
  
  return

 end subroutine init_output
 
 subroutine string_char(mychar,nstring,mystring)
 
  implicit none
  
  integer :: i
  character(1), allocatable, dimension(:) :: mychar
  integer, intent(in) :: nstring
  character(len=*), intent(in) :: mystring
  
  allocate(mychar(nstring))
  
  do i=1,nstring
    mychar(i)=mystring(i:i)
  enddo
  
 end subroutine string_char
 
 function space_fmtnumb(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading spaces to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=6) :: space_fmtnumb
  integer :: numdigit,irest
  real(kind=8) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=6-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(space_fmtnumb,fmt=cnumberlabel)repeat(' ',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(space_fmtnumb,fmt=cnumberlabel)inum
  endif
  
  return

 end function space_fmtnumb
 
 function space_fmtnumb12(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading TWELVE spaces to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=12) :: space_fmtnumb12
  integer :: numdigit,irest
  real(kind=8) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=12-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(space_fmtnumb12,fmt=cnumberlabel)repeat(' ',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(space_fmtnumb12,fmt=cnumberlabel)inum
  endif
  
  return

 end function space_fmtnumb12
 
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
  real(kind=db) :: tmp

  i=1
  tmp=real(inum,kind=db)
  do
  if(tmp< ten )exit
	i=i+1
	tmp=tmp/ ten
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
    
 subroutine get_memory_gpu(fout,fout2)

!***********************************************************************
!     
!     LBsoft subroutine for register the memory usage
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************  
#ifdef _OPENACC
  use openacc
  use accel_lib
#elif defined _CUDA  
  use cudafor
#endif
  
  implicit none
  
  real(kind=db), intent(out) :: fout,fout2
  real(kind=db) :: myd(2),myd2(2)
  integer :: istat
#ifdef _OPENACC  
  integer :: myfree, total
#elif defined _CUDA  
  integer(kind=cuda_count_kind) :: myfree, total
#else
  integer :: myfree, total
#endif  
  
#ifdef _OPENACC
  myfree=acc_get_free_memory()
  total=acc_get_memory() 
#elif defined _CUDA
  istat = cudaMemGetInfo( myfree, total )
#else
  myfree=0
  total=0
#endif  
  fout = real(total-myfree,kind=4)/(1024.0**3.0)
  fout2 = real(total,kind=4)/(1024.0**3.0)
  
  return
  
 end subroutine get_memory_gpu
    
 subroutine print_memory_registration_gpu(iu,mybanner,mybanner2,&
  mymemory,totmem)
 
!***********************************************************************
!     
!     LBcuda subroutine for printing the memory registration
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification April 2022
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  character(len=*), intent(in) :: mybanner,mybanner2
  real(kind=db), intent(in) :: mymemory,totmem
  
  character(len=12) :: r_char,r_char2
  
  character(len=*),parameter :: of='(a)'
  
  
  
 
  write (r_char,'(f12.4)')mymemory
  write (r_char2,'(f12.4)')totmem
  write(iu,of)"                                                                               "
  write(iu,of)"******************************GPU MEMORY MONITOR*******************************"
  write(iu,of)"                                                                               "
  write(iu,'(4a)')trim(mybanner)," = ",trim(adjustl(r_char))," (GB)"
  write(iu,'(4a)')trim(mybanner2)," = ",trim(adjustl(r_char2))," (GB)"
  write(iu,of)"                                                                               "
  write(iu,of)"*******************************************************************************"
  write(iu,of)"                                                                               "
  
  return
  
 end subroutine print_memory_registration_gpu
 
 end module

program lb_openacc
    
    use cudafor
    use mysubs
    use prints
    
    implicit none
    
    
    integer :: i,j,k
    integer :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus
    
    
    real(kind=4)  :: ts1,ts2
    real(kind=db) :: visc_LB,omega,oneminusomega
    real(kind=db) :: tau,one_ov_nu,fx,fy,fz
    real(kind=db) :: time

    real(kind=db) :: myrho,myu,myv,myw
    
    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    real(kind=db), allocatable, dimension(:,:,:) :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    integer :: TILE_DIMx,TILE_DIMy,TILE_DIMz,TILE_DIM,istat,iframe
    logical :: lprint,lvtk,lpbc,lasync
    real(kind=db) :: mymemory,totmemory
    integer(kind=cuda_Stream_Kind) :: stream1,stream2
    type (cudaDeviceProp) :: prop
    type (cudaEvent) :: startEvent, stopEvent, dummyEvent, dummyEvent1, dummyEvent2
    
       
    nlinks=18 !pari!
    tau=one
    visc_LB=cssq*(tau-half)
    one_ov_nu=one/visc_LB
    
    istat = cudaGetDeviceCount(ngpus)
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif

    !*******************************user parameters and allocations**************************
        nx=256
        ny=256
        nz=256
        
        myu=zero
        myv=zero
        myw=zero
        myrho=one  !tot dens
        nsteps=1000
        stamp=100
        lprint=.false.
        lvtk=.false.
        lpbc=.true.
        lasync=.false.
        fx=one*10.0**(-5)
        fy=zero*10.0**(-5)
        fz=zero*10.0**(-5)
        
        TILE_DIMx=128
        TILE_DIMy=2
        TILE_DIMz=1
        TILE_DIM=16
        if (mod(nx, TILE_DIMx)/= 0) then
          write(*,*) 'nx must be a multiple of TILE_DIM'
          stop
        end if
        if (mod(ny, TILE_DIMy) /= 0) then
          write(*,*) 'ny must be a multiple of TILE_DIMy'
          stop
        end if
        if (mod(nz, TILE_DIMz) /= 0) then
          write(*,*) 'nz must be a multiple of TILE_DIMz'
          stop
        end if
        dimGrid  = dim3(nx/TILE_DIMx, ny/TILE_DIMy, nz/TILE_DIMz)
        dimBlock = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
        
        dimGridx  = dim3((ny+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
        dimGridy  = dim3((nx+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
        dimBlock2 = dim3(TILE_DIM, TILE_DIM, 1)

        !allocate(f0(0:nx+1,0:ny+1,0:nz+1),f1(0:nx+1,0:ny+1,0:nz+1),f2(0:nx+1,0:ny+1,0:nz+1),f3(0:nx+1,0:ny+1,0:nz+1))
        !allocate(f4(0:nx+1,0:ny+1,0:nz+1),f5(0:nx+1,0:ny+1,0:nz+1),f6(0:nx+1,0:ny+1,0:nz+1),f7(0:nx+1,0:ny+1,0:nz+1))
        !allocate(f8(0:nx+1,0:ny+1,0:nz+1),f9(0:nx+1,0:ny+1,0:nz+1),f10(0:nx+1,0:ny+1,0:nz+1),f11(0:nx+1,0:ny+1,0:nz+1))
        !allocate(f12(0:nx+1,0:ny+1,0:nz+1),f13(0:nx+1,0:ny+1,0:nz+1),f14(0:nx+1,0:ny+1,0:nz+1),f15(0:nx+1,0:ny+1,0:nz+1))
        !allocate(f16(0:nx+1,0:ny+1,0:nz+1),f17(0:nx+1,0:ny+1,0:nz+1),f18(0:nx+1,0:ny+1,0:nz+1))
        !allocate(rho(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
        !allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
        !allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz))
        allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
        
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
        
        omega=one/tau
        oneminusomega=one-omega
    !*****************************************geometry************************
        isfluid=1
        isfluid(1,:,:)=0 !left
        isfluid(nx,:,:)=0 !right
        isfluid(:,1,:)=0 !front 
        isfluid(:,ny,:)=0 !rear
        isfluid(:,:,1)=0 !bottom
        isfluid(:,:,nz)=0 !top

        
    !********************************copies to device*****************!
        nx_d=nx
        ny_d=ny
        nz_d=nz
        TILE_DIMx_d=TILE_DIMx
        TILE_DIMy_d=TILE_DIMy
        TILE_DIMz_d=TILE_DIMz
        TILE_DIM_d=TILE_DIM
        myrho_d=myrho
        myu_d=myu
        myv_d=myv
        myw_d=myw
        fx_d=fx
        fy_d=fy
        fz_d=fz
        omega_d=omega
        oneminusomega_d=oneminusomega
        
        allocate(isfluid_d(1:nx_d,1:ny_d,1:nz_d))
        istat = cudaDeviceSynchronize
        istat = cudaMemcpy(isfluid_d,isfluid,nx*ny*nz )
        istat = cudaDeviceSynchronize
        allocate(f0_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f1_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f2_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f3_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(f4_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f5_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f6_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f7_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(f8_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f9_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f10_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f11_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(f12_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f13_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f14_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f15_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(f16_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f17_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f18_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(rho_d(1:nx_d,1:ny_d,1:nz_d),u_d(1:nx_d,1:ny_d,1:nz_d),v_d(1:nx_d,1:ny_d,1:nz_d),w_d(1:nx_d,1:ny_d,1:nz_d), &
         pxx_d(1:nx_d,1:ny_d,1:nz_d),pyy_d(1:nx_d,1:ny_d,1:nz_d),pzz_d(1:nx_d,1:ny_d,1:nz_d), &
         pxy_d(1:nx_d,1:ny_d,1:nz_d),pxz_d(1:nx_d,1:ny_d,1:nz_d),pyz_d(1:nx_d,1:ny_d,1:nz_d))
        istat = cudaDeviceSynchronize
    !*************************************initial conditions ************************    
        
        call setup_pops<<<dimGrid,dimBlock>>>()
        
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz),velprint(3,1:nx,1:ny,1:nz))
          allocate(rhoprint_d(1:nx_d,1:ny_d,1:nz_d),velprint_d(3,1:nx_d,1:ny_d,1:nz_d))
        endif
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
        write(6,*) 'lpbc',lpbc
        write(6,*) 'lprint',lprint
        write(6,*) 'lvtk',lvtk
        write(6,*) 'lasync',lasync
        write(6,*) 'nsteps',nsteps
        write(6,*) 'stamp',stamp
        write(6,*) 'max fx',huge(fx)
        write(6,*) 'max fx',huge(fy)
        write(6,*) 'max fx',huge(fz)
        write(6,*) 'TILE_DIMx ',TILE_DIMx
        write(6,*) 'TILE_DIMy ',TILE_DIMy
        write(6,*) 'TILE_DIMz ',TILE_DIMz
        write(6,*) 'TILE_DIM ',TILE_DIM
        write(6,*) 'available gpus',ngpus
        write(6,*) '*******************************************'
        
        istat = cudaGetDeviceProperties(prop, 0)
    
        call printDeviceProperties(prop,6, 1)
        
        ! create events and streams
        istat = cudaStreamCreate(stream1)
        istat = cudaStreamCreate(stream2)
        istat = cudaforSetDefaultstream(stream1)
        istat = cudaDeviceSynchronize
    
  
        istat = cudaEventCreate(startEvent)
        istat = cudaEventCreate(stopEvent)  
        istat = cudaEventCreate(dummyEvent)
        istat = cudaEventCreate(dummyEvent1)
        istat = cudaEventCreate(dummyEvent2)
        
    if(lprint)then  
      call init_output(nx,ny,nz,1,lvtk)
      call string_char(head1,nheadervtk(1),headervtk(1))
      call string_char(head2,nheadervtk(2),headervtk(2))
    endif
    
    istat = cudaDeviceSynchronize
    iframe=0
    if(lprint)then
      call moments<<<dimGrid,dimBlock,0,stream1>>>()
      call store_print<<<dimGrid,dimBlock,0,stream1>>>()
      istat = cudaEventRecord(dummyEvent1, stream1)
      istat = cudaEventSynchronize(dummyEvent1)
      !write(6,*)'ciao 1',step,iframe
      if(lasync)then
        istat = cudaMemcpyAsync(rhoprint,rhoprint_d,nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
        istat = cudaMemcpyAsync(velprint,velprint_d,3*nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
      else
        istat = cudaMemcpy(rhoprint,rhoprint_d,nx*ny*nz )
        istat = cudaMemcpy(velprint,velprint_d,3*nx*ny*nz )
        istat = cudaEventRecord(dummyEvent, 0)
        istat = cudaEventSynchronize(dummyEvent)
        if(lvtk)then
          call print_vtk_sync
        else
          call print_raw_sync
        endif
      endif
    endif
    
    !*************************************time loop************************  
    call cpu_time(ts1)
    istat = cudaEventRecord(startEvent,0)
    do step=1,nsteps 
        !***********************************moments collision bbck + forcing************************ 

        !call moments<<<dimGrid,dimBlock,0,stream1>>>()
        
        !***********************************PRINT************************
        if(mod(step,stamp).eq.0)write(6,'(a,i8)')'step : ',step
        if(lprint)then
          if(mod(step,stamp).eq.0)then
            iframe=iframe+1
            !write(6,*)'ciao 1',step,iframe
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
            call store_print<<<dimGrid,dimBlock,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
            call abortOnLastErrorAndSync('after store_print', step)
            if(lasync)then
              call close_print_async
              istat = cudaMemcpyAsync(rhoprint,rhoprint_d,nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
              istat = cudaMemcpyAsync(velprint,velprint_d,3*nx*ny*nz,cudaMemcpyDeviceToHost,stream2)
            else
              istat = cudaMemcpy(rhoprint,rhoprint_d,nx*ny*nz )
              istat = cudaMemcpy(velprint,velprint_d,3*nx*ny*nz )
              istat = cudaEventRecord(dummyEvent, 0)
              istat = cudaEventSynchronize(dummyEvent)
              if(lvtk)then
                call print_vtk_sync
              else
                call print_raw_sync
              endif
            endif
          endif
          if(mod(step-stamp/4,stamp).eq.0 .and. lasync)then
            !write(6,*)'ciao 2',step,iframe
            istat = cudaEventRecord(dummyEvent2, stream2)
            istat = cudaEventSynchronize(dummyEvent2)
            if(lvtk)then
              call print_vtk_async
            else
              call print_raw_async
            endif
          endif
        endif
        
        !***********************************collision + no slip + forcing: fused implementation*********
        call  streamcoll<<<dimGrid,dimBlock,0,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after streamcoll', step)
       
       
        !********************************boundary conditions no slip everywhere********************************!
        call bcs_no_slip<<<dimGrid,dimBlock,0,stream1>>>()
        istat = cudaEventRecord(dummyEvent1, stream1)
        istat = cudaEventSynchronize(dummyEvent1)
        call abortOnLastErrorAndSync('after bcs_no_slip', step)
        
        
        !******************************************call other bcs************************
            !periodic along x
            !periodic along y

        
          if(lpbc)then
            call pbc_edge_x<<<dimGridx,dimBlock2,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
            call abortOnLastErrorAndSync('after pbc_edge_x', step)
            call pbc_edge_y<<<dimGridy,dimBlock2,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
            call abortOnLastErrorAndSync('after pbc_edge_y', step)
          endif

          if(lprint)then
            istat = cudaEventRecord(dummyEvent, stream1)
            istat = cudaEventSynchronize(dummyEvent)
          endif
        
    enddo 

    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time, startEvent, stopEvent)

    call cpu_time(ts2)

    if(lasync .and. lprint)then
      !write(6,*)'ciao 2',step,iframe
      istat = cudaEventRecord(dummyEvent2, stream2)
      istat = cudaEventSynchronize(dummyEvent2)
      if(lvtk)then
        call print_vtk_sync
      else
        call print_raw_sync
      endif
    endif

    istat = cudaDeviceSynchronize

    if(lprint)then
    call store_print<<<dimGrid,dimBlock>>>()
    istat = cudaDeviceSynchronize
    
      istat = cudaMemcpy(rhoprint,rhoprint_d,nx*ny*nz )
      istat = cudaMemcpy(velprint,velprint_d,3*nx*ny*nz )
      istat = cudaDeviceSynchronize
      write(6,*) 'u=',velprint(1,nx/2,ny/2,nz/2),'v=',velprint(2,nx/2,ny/2,nz/2),'w=',velprint(3,nx/2,ny/2,nz/2),'rho=',rhoprint(nx/2,ny/2,nz/2)
      write(6,*) 'u=',velprint(1,nx/2,ny/2,1),'v=',velprint(2,nx/2,ny/2,1),'w=',velprint(3,nx/2,ny/2,1),'rho=',rhoprint(nx/2,ny/2,1)

    endif
    
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time'
    write(6,*) 'cuda time elapsed: ', time, ' s of your life time'
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)*real(1.d-9,kind=db)/(ts2-ts1)
    write(6,*) 'glups cuda: ',  real(nx)*real(ny)*real(nz)*real(nsteps)*real(1.d-9,kind=db)*1000/time
    
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)

    contains
  
      subroutine printDeviceProperties(prop,iu,num)
      
          use cudafor
          type(cudadeviceprop) :: prop
          integer,intent(in) :: iu,num 
          
          write(iu,907)"                                                                               "
          write(iu,907)"*****************************GPU FEATURE MONITOR*******************************"
          write(iu,907)"                                                                               "
          
          write (iu,900) "Device Number: "      ,num
          write (iu,901) "Device Name: "        ,trim(prop%name)
          write (iu,903) "Total Global Memory: ",real(prop%totalGlobalMem)/1e9," Gbytes"
          write (iu,902) "sharedMemPerBlock: "  ,prop%sharedMemPerBlock," bytes"
          write (iu,900) "regsPerBlock: "       ,prop%regsPerBlock
          write (iu,900) "warpSize: "           ,prop%warpSize
          write (iu,900) "maxThreadsPerBlock: " ,prop%maxThreadsPerBlock
          write (iu,904) "maxThreadsDim: "      ,prop%maxThreadsDim
          write (iu,904) "maxGridSize: "        ,prop%maxGridSize
          write (iu,903) "ClockRate: "          ,real(prop%clockRate)/1e6," GHz"
          write (iu,902) "Total Const Memory: " ,prop%totalConstMem," bytes"
          write (iu,905) "Compute Capability Revision: ",prop%major,prop%minor
          write (iu,902) "TextureAlignment: "   ,prop%textureAlignment," bytes"
          write (iu,906) "deviceOverlap: "      ,prop%deviceOverlap
          write (iu,900) "multiProcessorCount: ",prop%multiProcessorCount
          write (iu,906) "integrated: "         ,prop%integrated
          write (iu,906) "canMapHostMemory: "   ,prop%canMapHostMemory
          write (iu,906) "ECCEnabled: "         ,prop%ECCEnabled
          write (iu,906) "UnifiedAddressing: "  ,prop%unifiedAddressing
          write (iu,900) "L2 Cache Size: "      ,prop%l2CacheSize
          write (iu,900) "maxThreadsPerSMP: "   ,prop%maxThreadsPerMultiProcessor
          
          write(iu,907)"                                                                               "
          write(iu,907)"*******************************************************************************"
          write(iu,907)"                                                                               "
          
          900 format (a,i0)
          901 format (a,a)
          902 format (a,i0,a)
          903 format (a,f16.8,a)
          904 format (a,2(i0,1x,'x',1x),i0)
          905 format (a,i0,'.',i0)
          906 format (a,l0)
          907 format (a)
          
          return
      
      end subroutine printDeviceProperties
  
      subroutine print_raw_sync
      
          implicit none
          
          
          sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
            '_'//trim(write_fmtnumb(iframe)) // '.raw'
          sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
            '_'//trim(write_fmtnumb(iframe)) // '.raw'
          open(unit=345,file=trim(sevt1), &
            status='replace',action='write',access='stream',form='unformatted')
          write(345)rhoprint
          close(345)
          open(unit=346,file=trim(sevt2), &
            status='replace',action='write',access='stream',form='unformatted')
          write(346)velprint
          close(346)
      
      end subroutine print_raw_sync
  
      subroutine print_vtk_sync
          implicit none
            
          sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
            '_'//trim(write_fmtnumb(iframe)) // '.vti'
          sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
            '_'//trim(write_fmtnumb(iframe)) // '.vti'
          open(unit=345,file=trim(sevt1), &
            status='replace',action='write',access='stream',form='unformatted')
          write(345)head1,ndatavtk(1),rhoprint,footervtk(1)
          close(345)
          open(unit=346,file=trim(sevt2), &
            status='replace',action='write',access='stream',form='unformatted')
          write(346)head2,ndatavtk(2),velprint,footervtk(2)
          close(346)
      
      end subroutine print_vtk_sync
  
      subroutine print_raw_async
      
          implicit none
          
          
          sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
            '_'//trim(write_fmtnumb(iframe)) // '.raw'
          sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
            '_'//trim(write_fmtnumb(iframe)) // '.raw'
          open(unit=345,file=trim(sevt1), &
            status='replace',action='write',access='stream',form='unformatted',&
            asynchronous='yes')
          write(345,asynchronous='yes')rhoprint
          
          open(unit=346,file=trim(sevt2), &
            status='replace',action='write',access='stream',form='unformatted',&
            asynchronous='yes')
          write(346,asynchronous='yes')velprint
      
      
      end subroutine print_raw_async
  
      subroutine print_vtk_async
          implicit none
            
          sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
            '_'//trim(write_fmtnumb(iframe)) // '.vti'
          sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
            '_'//trim(write_fmtnumb(iframe)) // '.vti'
            
          open(unit=345,file=trim(sevt1), &
            status='replace',action='write',access='stream',form='unformatted',&
            asynchronous='yes')
          write(345,asynchronous='yes')head1,ndatavtk(1),rhoprint
          
          
          open(unit=780,file=trim(sevt2), &
            status='replace',action='write',access='stream',form='unformatted',&
            asynchronous='yes')
          write(780,asynchronous='yes')head2,ndatavtk(2),velprint
      
      end subroutine print_vtk_async
  
      subroutine close_print_async
      
          implicit none
          
          wait(345)
          if(lvtk)write(345)footervtk(1)
          close(345)
          
          
          wait(780)
          if(lvtk)write(780)footervtk(2)
          close(780) 
      
      end subroutine close_print_async
    
end program
