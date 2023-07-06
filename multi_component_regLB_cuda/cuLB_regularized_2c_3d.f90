 
 module mysubs
   
   use cudafor
  
   implicit none
    integer, parameter :: db=4 !kind(1.0)
    ! device arrays and vars
    real(kind=db),parameter :: pi=3.14159265359793234626433
    integer(kind=4), allocatable,  dimension(:,:,:), device   :: isfluid_d
    integer, constant :: nx_d,ny_d,nz_d,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
    real(kind=db), constant :: p0_d,p1_d,p2_d,p1cg_d,p2cg_d,p3cg_d
    real(kind=db), constant :: fx_d,fy_d,fz_d,myrhoA_d,myrhoB_d,myu_d,myv_d,myw_d,cssq_d
    real(kind=db), constant :: p1dcssq_d,p2dcssq_d,pi2cssq0_d,pi2cssq1_d,pi2cssq2_d,qxx_d,qyy_d,qzz_d
    real(kind=db), constant :: one_ov_nu1_d, one_ov_nu2_d, sigma_d, b0_d, b1_d, b2_d, beta_d, max_press_excess_d
    real(kind=db), constant :: qxy_7_8_d,qxy_9_10_d,qxz_15_16_d,qxz_17_18_d,qyz_11_12_d,qyz_13_14_d
    real(kind=db), device :: omega_d
    !arrays
    integer*4, allocatable, dimension(:,:,:), device :: nci_loc_d
    real(kind=db), allocatable, dimension(:,:,:), device  :: rhoA_d,rhoB_d,psi_d,u_d,v_d,w_d,pxx_d,pyy_d,pzz_d,pxy_d,pxz_d,pyz_d
    real(kind=db), allocatable, dimension(:,:,:), device :: f0_d,f1_d,f2_d,f3_d,f4_d,f5_d,f6_d,f7_d,f8_d,f9_d
    real(kind=db), allocatable, dimension(:,:,:), device :: f10_d,f11_d,f12_d,f13_d,f14_d,f15_d,f16_d,f17_d,f18_d
    real(kind=db), allocatable, dimension(:,:,:), device :: g0_d,g1_d,g2_d,g3_d,g4_d,g5_d,g6_d,g7_d,g8_d,g9_d
    real(kind=db), allocatable, dimension(:,:,:), device :: g10_d,g11_d,g12_d,g13_d,g14_d,g15_d,g16_d,g17_d,g18_d
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

    attributes(device) function fcut(r, r1, r2)
  
        real(kind=db), intent(in) :: r, r1, r2
        real(kind=db) fcut

        if ( r <= r1 ) then
          fcut = 1.0_db
        elseif ( r > r2 ) then
          fcut = 0.0_db
        else
          fcut = 0.5_db * cos((r-r1)*Pi/(r2-r1)) + 0.5_db
        endif
    end function fcut
  
    attributes(global) subroutine setup_pops(myradius)
        
          real(kind=db), value :: myradius
          integer :: i,j,k
          real(kind=db) :: mydist,locpsi,locrhoA,locrhoB,tempr
      
        
          i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
          j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
          k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
        
          mydist=sqrt((real(i)-real(nx_d)*0.5_db)**2.0+(real(j)-real(ny_d)*0.5_db)**2.0 + (real(k)-real(nz_d)*0.5_db)**2.0)
              
          !tempr = fcut(mydist, myradius, myradius+4.0_db)
          
          locpsi=-1.0_db
          if (mydist<=20.0_db) locpsi=1.0_db
          
          !locrhoB = tempr
          !locrhoA = (1.0 - tempr)
          !locpsi = (locrhoA - locrhoB)/(locrhoA + locrhoB)
          
          locrhoB=0.5_db*(1.0_db-locpsi)
          locrhoA=1.0_db-locrhoB
          
          !if(i==nx_d/2 .and. j==ny_d/2)write(*,*)'cazzo',i,j
          
          u_d(i,j,k)=0.0_db
          v_d(i,j,k)=0.0_db
          w_d(i,j,k)=0.0_db

          !write(*,*)i,j,p_d(0)*myrho_d
          f0_d(i,j,k) = p0_d*locrhoA
          f1_d(i,j,k) = p1_d*locrhoA
          f2_d(i,j,k) = p1_d*locrhoA
          f3_d(i,j,k) = p1_d*locrhoA
          f4_d(i,j,k) = p1_d*locrhoA
          f5_d(i,j,k) = p1_d*locrhoA
          f6_d(i,j,k) = p1_d*locrhoA
          f7_d(i,j,k) = p2_d*locrhoA
          f8_d(i,j,k) = p2_d*locrhoA
          f9_d(i,j,k) = p2_d*locrhoA
          f10_d(i,j,k) = p2_d*locrhoA
          f11_d(i,j,k) = p2_d*locrhoA
          f12_d(i,j,k) = p2_d*locrhoA
          f13_d(i,j,k) = p2_d*locrhoA
          f14_d(i,j,k) = p2_d*locrhoA
          f15_d(i,j,k) = p2_d*locrhoA
          f16_d(i,j,k) = p2_d*locrhoA
          f17_d(i,j,k) = p2_d*locrhoA
          f18_d(i,j,k) = p2_d*locrhoA

          g0_d(i,j,k) = p0_d*locrhob
          g1_d(i,j,k) = p1_d*locrhob
          g2_d(i,j,k) = p1_d*locrhob
          g3_d(i,j,k) = p1_d*locrhob
          g4_d(i,j,k) = p1_d*locrhob
          g5_d(i,j,k) = p1_d*locrhob
          g6_d(i,j,k) = p1_d*locrhob
          g7_d(i,j,k) = p2_d*locrhob
          g8_d(i,j,k) = p2_d*locrhob
          g9_d(i,j,k) = p2_d*locrhob
          g10_d(i,j,k) = p2_d*locrhob
          g11_d(i,j,k) = p2_d*locrhob
          g12_d(i,j,k) = p2_d*locrhob
          g13_d(i,j,k) = p2_d*locrhob
          g14_d(i,j,k) = p2_d*locrhob
          g15_d(i,j,k) = p2_d*locrhob
          g16_d(i,j,k) = p2_d*locrhob
          g17_d(i,j,k) = p2_d*locrhob
          g18_d(i,j,k) = p2_d*locrhob
          
          rhoA_d(i,j,k)=locrhoA
          rhob_d(i,j,k)=locrhoB
          psi_d(i,j,k)=locpsi
        

    end subroutine setup_pops

    attributes(global) subroutine moments()
        
          integer :: i,j,k
          real(kind=db) :: uu,udotc,temp
          real(kind=db) :: fneq1,rtot
            
          i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
          j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
          k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
            
          if(isfluid_d(i,j,k).ne.1)return

          pxx_d(i,j,k)=0.0_db
          pyy_d(i,j,k)=0.0_db
          pzz_d(i,j,k)=0.0_db
          pxy_d(i,j,k)=0.0_db
          pxz_d(i,j,k)=0.0_db
          pyz_d(i,j,k)=0.0_db
            
          rhoA_d(i,j,k) = f0_d(i,j,k)+f1_d(i,j,k)+f2_d(i,j,k)+f3_d(i,j,k)+f4_d(i,j,k)+f5_d(i,j,k) &
          +f6_d(i,j,k)+f7_d(i,j,k)+f8_d(i,j,k)+f9_d(i,j,k)+f10_d(i,j,k)+f11_d(i,j,k) &
          +f12_d(i,j,k)+f13_d(i,j,k)+f14_d(i,j,k)+f15_d(i,j,k)+f16_d(i,j,k)+f17_d(i,j,k) &
          +f18_d(i,j,k)

          rhoB_d(i,j,k) = g0_d(i,j,k)+g1_d(i,j,k)+g2_d(i,j,k)+g3_d(i,j,k)+g4_d(i,j,k)+g5_d(i,j,k) &
          +g6_d(i,j,k)+g7_d(i,j,k)+g8_d(i,j,k)+g9_d(i,j,k)+g10_d(i,j,k)+g11_d(i,j,k) &
          +g12_d(i,j,k)+g13_d(i,j,k)+g14_d(i,j,k)+g15_d(i,j,k)+g16_d(i,j,k)+g17_d(i,j,k) &
          +g18_d(i,j,k)

          u_d(i,j,k) = (f1_d(i,j,k)+f7_d(i,j,k)+f9_d(i,j,k)+f15_d(i,j,k)+f18_d(i,j,k)) &
                      -(f2_d(i,j,k)+f8_d(i,j,k)+f10_d(i,j,k)+f16_d(i,j,k)+f17_d(i,j,k)) &
                      +(g1_d(i,j,k)+g7_d(i,j,k)+g9_d(i,j,k)+g15_d(i,j,k)+g18_d(i,j,k)) &
                      -(g2_d(i,j,k)+g8_d(i,j,k)+g10_d(i,j,k)+g16_d(i,j,k)+g17_d(i,j,k)) 
          
          v_d(i,j,k) = (f3_d(i,j,k)+f7_d(i,j,k)+f10_d(i,j,k)+f11_d(i,j,k)+f13_d(i,j,k)) &
                      -(f4_d(i,j,k)+f8_d(i,j,k)+f9_d(i,j,k)+f12_d(i,j,k)+f14_d(i,j,k)) &
                      +(g3_d(i,j,k)+g7_d(i,j,k)+g10_d(i,j,k)+g11_d(i,j,k)+g13_d(i,j,k)) &
                      -(g4_d(i,j,k)+g8_d(i,j,k)+g9_d(i,j,k)+g12_d(i,j,k)+g14_d(i,j,k))

          w_d(i,j,k) = (f5_d(i,j,k)+f11_d(i,j,k)+f14_d(i,j,k)+f15_d(i,j,k)+f17_d(i,j,k)) &
                      -(f6_d(i,j,k)+f12_d(i,j,k)+f13_d(i,j,k)+f16_d(i,j,k)+f18_d(i,j,k)) &
                      +(g5_d(i,j,k)+g11_d(i,j,k)+g14_d(i,j,k)+g15_d(i,j,k)+g17_d(i,j,k)) &
                      -(g6_d(i,j,k)+g12_d(i,j,k)+g13_d(i,j,k)+g16_d(i,j,k)+g18_d(i,j,k))
        
          uu=0.5_db*(u_d(i,j,k)*u_d(i,j,k) + v_d(i,j,k)*v_d(i,j,k) + w_d(i,j,k)*w_d(i,j,k))/cssq_d

          rtot=rhoA_d(i,j,k)+rhoB_d(i,j,k)

          psi_d(i,j,k)= (rhoA_d(i,j,k)-rhoB_d(i,j,k))/(rhoA_d(i,j,k)+rhoB_d(i,j,k))  
          nci_loc_d(i,j,k)=0 
          !1-2
          udotc=u_d(i,j,k)/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f1_d(i,j,k) + g1_d(i,j,k)-p1_d*(rtot+(temp + udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          fneq1=f2_d(i,j,k) + g2_d(i,j,k)-p1_d*(rtot+(temp - udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          !3-4
          udotc=v_d(i,j,k)/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f3_d(i,j,k)+g3_d(i,j,k)-p1_d*(rtot+(temp + udotc))
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          fneq1=f4_d(i,j,k)+g4_d(i,j,k)-p1_d*(rtot+(temp - udotc))
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          !5-6
          udotc=w_d(i,j,k)/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f5_d(i,j,k)+g5_d(i,j,k)-p1_d*(rtot+(temp + udotc))
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          fneq1=f6_d(i,j,k)+g6_d(i,j,k)-p1_d*(rtot+(temp - udotc))
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          !7-8
          udotc=(u_d(i,j,k)+v_d(i,j,k))/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f7_d(i,j,k)+g7_d(i,j,k)-p2_d*(rtot+(temp + udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          pxy_d(i,j,k)=pxy_d(i,j,k)+fneq1
          fneq1=f8_d(i,j,k)+g8_d(i,j,k)-p2_d*(rtot+(temp - udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          pxy_d(i,j,k)=pxy_d(i,j,k)+fneq1
          !10-9
          udotc=(-u_d(i,j,k)+v_d(i,j,k))/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f10_d(i,j,k)+g10_d(i,j,k)-p2_d*(rtot+(temp + udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          pxy_d(i,j,k)=pxy_d(i,j,k)-fneq1
          fneq1=f9_d(i,j,k)+g9_d(i,j,k)-p2_d*(rtot+(temp - udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          pxy_d(i,j,k)=pxy_d(i,j,k)-fneq1
          !11-12
          udotc=(v_d(i,j,k)+w_d(i,j,k))/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f11_d(i,j,k)+g11_d(i,j,k)-p2_d*(rtot+(temp + udotc))
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          pyz_d(i,j,k)=pyz_d(i,j,k)+fneq1
          fneq1=f12_d(i,j,k)+g12_d(i,j,k)-p2_d*(rtot+(temp - udotc))
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          pyz_d(i,j,k)=pyz_d(i,j,k)+fneq1
          !13-14
          udotc=(v_d(i,j,k)-w_d(i,j,k))/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f13_d(i,j,k)+g13_d(i,j,k) - p2_d*(rtot+(temp + udotc))
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          pyz_d(i,j,k)=pyz_d(i,j,k)-fneq1
          fneq1=f14_d(i,j,k)+g14_d(i,j,k) - p2_d*(rtot+(temp - udotc))
          pyy_d(i,j,k)=pyy_d(i,j,k)+fneq1
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          pyz_d(i,j,k)=pyz_d(i,j,k)-fneq1
          !15-16
          udotc=(u_d(i,j,k)+w_d(i,j,k))/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f15_d(i,j,k)+g15_d(i,j,k)-p2_d*(rtot+(temp + udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          pxz_d(i,j,k)=pyz_d(i,j,k)+fneq1
          fneq1=f16_d(i,j,k)+g16_d(i,j,k)-p2_d*(rtot+(temp - udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          pxz_d(i,j,k)=pyz_d(i,j,k)+fneq1
          !17-18
          udotc=(-u_d(i,j,k)+w_d(i,j,k))/cssq_d
          temp = -uu + 0.5_db*udotc*udotc
          fneq1=f17_d(i,j,k)+g17_d(i,j,k)-p2_d*(rtot+(temp + udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          pxz_d(i,j,k)=pyz_d(i,j,k)-fneq1
          fneq1=f18_d(i,j,k)+g18_d(i,j,k)-p2_d*(rtot+(temp - udotc))
          pxx_d(i,j,k)=pxx_d(i,j,k)+fneq1
          pzz_d(i,j,k)=pzz_d(i,j,k)+fneq1
          pxz_d(i,j,k)=pyz_d(i,j,k)-fneq1
  
    end subroutine moments
  
    attributes(global) subroutine streamcoll()
        
        integer :: i,j,k
        real(kind=db) :: uu, udotc, temp, feq, psi_x, psi_y,psi_z, rtot, st_coeff, mod_psi, mod_psi_sq, fpc
        real(kind=db) :: norm_x, norm_y, norm_z, rprod    
        real(kind=db) :: addendum0
        real(kind=db) :: gaddendum0
        real(kind=db) :: ushifted, vshifted, wshifted, nu_avg
        
          
        i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
        j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
        k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
          
        if(isfluid_d(i,j,k).ne.1)return

        psi_x=(1.0_db/cssq_d)*(p1cg_d*(psi_d(i+1,j,k)-psi_d(i-1,j,k)) + &
                            p2cg_d*(psi_d(i+1,j+1,k)+psi_d(i+1,j-1,k)+psi_d(i+1,j,k+1)+psi_d(i+1,j,k-1)) &
                        -p2cg_d*(psi_d(i-1,j+1,k)+psi_d(i-1,j-1,k)+psi_d(i-1,j,k+1)+psi_d(i-1,j,k-1)) &
                        + p3cg_d*(psi_d(i+1,j+1,k+1)+psi_d(i+1,j+1,k-1)+psi_d(i+1,j-1,k+1)+psi_d(i+1,j-1,k-1)) &
                        -p3cg_d*(psi_d(i-1,j+1,k+1)+psi_d(i-1,j+1,k-1)+psi_d(i-1,j-1,k+1)+psi_d(i-1,j-1,k-1))) 


        psi_y=(1.0_db/cssq_d)*(p1cg_d*(psi_d(i,j+1,k)-psi_d(i,j-1,k)) &
                            +p2cg_d*(psi_d(i+1,j+1,k)+psi_d(i-1,j+1,k)+psi_d(i,j+1,k+1)+psi_d(i,j+1,k-1)) &
                        -p2cg_d*(psi_d(i+1,j-1,k)+psi_d(i-1,j-1,k)+psi_d(i,j-1,k+1)+psi_d(i,j-1,k-1)) &
                        +p3cg_d*(psi_d(i+1,j+1,k+1)+psi_d(i-1,j+1,k+1)+psi_d(i+1,j+1,k-1)+psi_d(i-1,j+1,k-1)) &
                        -p3cg_d*(psi_d(i+1,j-1,k+1)+psi_d(i-1,j-1,k+1)+psi_d(i+1,j-1,k-1)+psi_d(i-1,j-1,k-1)))


        psi_z=(1.0_db/cssq_d)*(p1cg_d*(psi_d(i,j,k+1)-psi_d(i,j,k-1)) + &
                            p2cg_d*(psi_d(i+1,j,k+1)+psi_d(i-1,j,k+1)+ psi_d(i,j+1,k+1)+psi_d(i,j-1,k+1)) &
                        -p2cg_d*(psi_d(i+1,j,k-1)+psi_d(i-1,j,k-1)+psi_d(i,j+1,k-1)+psi_d(i,j-1,k-1)) &
                        +p3cg_d*(psi_d(i+1,j+1,k+1)+psi_d(i-1,j+1,k+1)+psi_d(i+1,j-1,k+1)+psi_d(i-1,j-1,k+1)) &
                        -p3cg_d*(psi_d(i+1,j+1,k-1)+psi_d(i-1,j+1,k-1)+psi_d(i+1,j-1,k-1)+psi_d(i-1,j-1,k-1)))   
        mod_psi=sqrt(psi_x**2+psi_y**2+psi_z**2)
        mod_psi_sq=psi_x**2+psi_y**2+psi_z**2
        norm_x=0.0_db
        norm_y=0.0_db
        norm_z=0.0_db
        rtot=0.0_db
        rtot=rhoA_d(i,j,k)+rhoB_d(i,j,k)
        rprod=rhoA_d(i,j,k)*rhoB_d(i,j,k)
        nu_avg=1.0_db/(rhoA_d(i,j,k)*one_ov_nu1_d/rtot + rhoB_d(i,j,k)*one_ov_nu2_d/rtot)
        omega_d=2.0_db/(6.0_db*nu_avg + 1.0_db)
        st_coeff=(9.0_db/4.0_db)*sigma_d*omega_d
        addendum0=0.0_db
        gaddendum0=0.0_db  

        if(mod_psi>0.0001)then
            norm_x=psi_x/mod_psi
            norm_y=psi_y/mod_psi
            norm_z=psi_z/mod_psi
                  
            ushifted=u_d(i,j,k) + fx_d + float(nci_loc_d(i,j,k))*(norm_x)*max_press_excess_d*abs(rhob_d(i,j,k))
            vshifted=v_d(i,j,k) + fy_d + float(nci_loc_d(i,j,k))*(norm_y)*max_press_excess_d*abs(rhob_d(i,j,k))
            wshifted=w_d(i,j,k) + fz_d + float(nci_loc_d(i,j,k))*(norm_z)*max_press_excess_d*abs(rhob_d(i,j,k))

            addendum0=-st_coeff*mod_psi*b0_d
            uu=0.5_db*(ushifted*ushifted + vshifted*vshifted + wshifted*wshifted)/cssq_d 
            !0
            feq=p0_d*(rtot-uu)
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq0_d*(-cssq_d*(pyy_d(i,j,k)+pxx_d(i,j,k)+pzz_d(i,j,k))) 
            f0_d(i,j,k)=fpc*rhoA_d(i,j,k)/rtot
            g0_d(i,j,k)=fpc*rhoB_d(i,j,k)/rtot
            !1
            udotc=ushifted/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p1_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p1_d*psi_x**2/mod_psi_sq - b1_d)
            gaddendum0=p1_d*(rtot)*(rprod*beta_d*psi_x/mod_psi/rtot**2)
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pxx_d(i,j,k)-cssq_d*(pyy_d(i,j,k)+pzz_d(i,j,k)))
            f1_d(i+1,j,k)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g1_d(i+1,j,k)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            !2
            feq=p1_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pxx_d(i,j,k)-cssq_d*(pyy_d(i,j,k)+pzz_d(i,j,k))) 
            f2_d(i-1,j,k)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g2_d(i-1,j,k)=fpc*rhob_d(i,j,k)/rtot + gaddendum0
            !3
            udotc=vshifted/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p1_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p1_d*psi_y**2/mod_psi_sq - b1_d)
            gaddendum0=p1_d*(rtot)*(rprod*beta_d*psi_y/mod_psi/rtot**2)
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*(pxx_d(i,j,k)+pzz_d(i,j,k))) 
            f3_d(i,j+1,k)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g3_d(i,j+1,k)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            !4
            feq=p1_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*(pxx_d(i,j,k)+pzz_d(i,j,k))) 
            f4_d(i,j-1,k)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g4_d(i,j-1,k)=fpc*rhob_d(i,j,k)/rtot + gaddendum0
            !7
            udotc=(ushifted+vshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p2_d*(psi_x+psi_y)**2/mod_psi_sq - b2_d)
            gaddendum0=p2_d*(rtot)*(rprod*beta_d*(psi_x/mod_psi+psi_y/mod_psi)/rtot**2)
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*pzz_d(i,j,k)+2.0_db*pxy_d(i,j,k))
            f7_d(i+1,j+1,k)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g7_d(i+1,j+1,k)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            !8
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*pzz_d(i,j,k)+2.0_db*pxy_d(i,j,k))
            f8_d(i-1,j-1,k)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g8_d(i-1,j-1,k)=fpc*rhob_d(i,j,k)/rtot + gaddendum0 
            !10
            udotc=(-ushifted+vshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p2_d*(-psi_x+psi_y)**2/mod_psi_sq - b2_d)
            gaddendum0=p2_d*(rtot)*(rprod*beta_d*(-psi_x/mod_psi+psi_y/mod_psi)/rtot**2)
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*pzz_d(i,j,k)-2.0_db*pxy_d(i,j,k)) 
            f10_d(i-1,j+1,k)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g10_d(i-1,j+1,k)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            !9
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*pzz_d(i,j,k)-2.0_db*pxy_d(i,j,k))  
            f9_d(i+1,j-1,k)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g9_d(i+1,j-1,k)=fpc*rhob_d(i,j,k)/rtot + gaddendum0
            !5
            udotc=wshifted/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p1_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p1_d*psi_z**2/mod_psi_sq - b1_d)
            gaddendum0=p1_d*(rtot)*(rprod*beta_d*psi_z/mod_psi/rtot**2)
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*(pxx_d(i,j,k)+pyy_d(i,j,k)))
            f5_d(i,j,k+1)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g5_d(i,j,k+1)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            !6
            feq=p1_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*(pxx_d(i,j,k)+pyy_d(i,j,k)))
            f6_d(i,j,k-1)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g6_d(i,j,k-1)=fpc*rhob_d(i,j,k)/rtot + gaddendum0
            !15
            udotc=(ushifted+wshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p2_d*(psi_x+psi_z)**2/mod_psi_sq - b2_d)
            gaddendum0=p2_d*(rtot)*(rprod*beta_d*(psi_x/mod_psi+psi_z/mod_psi)/rtot**2)
            fpc=feq + addendum0+(1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pyy_d(i,j,k)+2.0_db*pxz_d(i,j,k)) 
            f15_d(i+1,j,k+1)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g15_d(i+1,j,k+1)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            !16
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pyy_d(i,j,k)+2.0_db*pxz_d(i,j,k)) 
            f16_d(i-1,j,k-1)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g16_d(i-1,j,k-1)=fpc*rhob_d(i,j,k)/rtot + gaddendum0
            !17
            udotc=(-ushifted+wshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p2_d*(-psi_x+psi_z)**2/mod_psi_sq - b2_d)
            gaddendum0=p2_d*(rtot)*(rprod*beta_d*(-psi_x/mod_psi+psi_z/mod_psi)/rtot**2)
            fpc=feq +addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pyy_d(i,j,k)-2.0_db*pxz_d(i,j,k))
            f17_d(i-1,j,k+1)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g17_d(i-1,j,k+1)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            !18
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pyy_d(i,j,k)-2.0_db*pxz_d(i,j,k))
            f18_d(i+1,j,k-1)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g18_d(i+1,j,k-1)=fpc*rhob_d(i,j,k)/rtot + gaddendum0
            !11
            udotc=(vshifted+wshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p2_d*(psi_y+psi_z)**2/mod_psi_sq - b2_d)
            gaddendum0=p2_d*(rtot)*(rprod*beta_d*(psi_y/mod_psi+psi_z/mod_psi)/rtot**2)
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pyy_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pxx_d(i,j,k)+2.0_db*pyz_d(i,j,k))
            f11_d(i,j+1,k+1)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g11_d(i,j+1,k+1)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            
            !12
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pyy_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pxx_d(i,j,k)+2.0_db*pyz_d(i,j,k)) 
            f12_d(i,j-1,k-1)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g12_d(i,j-1,k-1)=fpc*rhob_d(i,j,k)/rtot + gaddendum0
            !13
            udotc=(vshifted-wshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            addendum0=st_coeff*mod_psi*(p2_d*(psi_y-psi_z)**2/mod_psi_sq - b2_d)
            gaddendum0=p2_d*(rtot)*(rprod*beta_d*(psi_y/mod_psi-psi_z/mod_psi)/rtot**2)
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pyy_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pxx_d(i,j,k)-2.0_db*pyz_d(i,j,k)) 
            f13_d(i,j+1,k-1)=fpc*rhoA_d(i,j,k)/rtot + gaddendum0
            g13_d(i,j+1,k-1)=fpc*rhob_d(i,j,k)/rtot - gaddendum0
            
            !14
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq + addendum0 + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pyy_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pxx_d(i,j,k)-2.0_db*pyz_d(i,j,k)) 
            f14_d(i,j-1,k+1)=fpc*rhoA_d(i,j,k)/rtot - gaddendum0
            g14_d(i,j-1,k+1)=fpc*rhob_d(i,j,k)/rtot + gaddendum0
        else
            ushifted=u_d(i,j,k) + fx_d 
            vshifted=v_d(i,j,k) + fy_d
            wshifted=w_d(i,j,k) + fz_d 
            uu=0.5_db*(ushifted*ushifted + vshifted*vshifted + wshifted*wshifted)/cssq_d
            !0
            feq=p0_d*(rtot-uu)
            fpc=feq  + (1.0_db-omega_d)*pi2cssq0_d*(-cssq_d*(pyy_d(i,j,k)+pxx_d(i,j,k)+pzz_d(i,j,k))) 
            f0_d(i,j,k)=fpc*rhoA_d(i,j,k)/rtot
            g0_d(i,j,k)=fpc*rhoB_d(i,j,k)/rtot
            !1
            udotc=ushifted/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p1_d*(rtot+(temp + udotc))
            fpc=feq +(1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pxx_d(i,j,k)-cssq_d*(pyy_d(i,j,k)+pzz_d(i,j,k)))
            f1_d(i+1,j,k)=fpc*rhoA_d(i,j,k)/rtot 
            g1_d(i+1,j,k)=fpc*rhob_d(i,j,k)/rtot 
            !2
            feq=p1_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pxx_d(i,j,k)-cssq_d*(pyy_d(i,j,k)+pzz_d(i,j,k))) 
            f2_d(i-1,j,k)=fpc*rhoA_d(i,j,k)/rtot 
            g2_d(i-1,j,k)=fpc*rhob_d(i,j,k)/rtot 
            !3
            udotc=vshifted/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p1_d*(rtot+(temp + udotc))
            fpc=feq +(1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*(pxx_d(i,j,k)+pzz_d(i,j,k))) 
            f3_d(i,j+1,k)=fpc*rhoA_d(i,j,k)/rtot 
            g3_d(i,j+1,k)=fpc*rhob_d(i,j,k)/rtot 
            !4
            feq=p1_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*(pxx_d(i,j,k)+pzz_d(i,j,k))) 
            f4_d(i,j-1,k)=fpc*rhoA_d(i,j,k)/rtot 
            g4_d(i,j-1,k)=fpc*rhob_d(i,j,k)/rtot 
            !7
            udotc=(ushifted+vshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*pzz_d(i,j,k)+2.0_db*pxy_d(i,j,k))
            f7_d(i+1,j+1,k)=fpc*rhoA_d(i,j,k)/rtot 
            g7_d(i+1,j+1,k)=fpc*rhob_d(i,j,k)/rtot 
            !8
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*pzz_d(i,j,k)+2.0_db*pxy_d(i,j,k))
            f8_d(i-1,j-1,k)=fpc*rhoA_d(i,j,k)/rtot 
            g8_d(i-1,j-1,k)=fpc*rhob_d(i,j,k)/rtot  
            !10
            udotc=(-ushifted+vshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*pzz_d(i,j,k)-2.0_db*pxy_d(i,j,k)) 
            f10_d(i-1,j+1,k)=fpc*rhoA_d(i,j,k)/rtot 
            g10_d(i-1,j+1,k)=fpc*rhob_d(i,j,k)/rtot 
            !9
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pyy_d(i,j,k)-cssq_d*pzz_d(i,j,k)-2.0_db*pxy_d(i,j,k))  
            f9_d(i+1,j-1,k)=fpc*rhoA_d(i,j,k)/rtot 
            g9_d(i+1,j-1,k)=fpc*rhob_d(i,j,k)/rtot 
            !5
            udotc=wshifted/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p1_d*(rtot+(temp + udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*(pxx_d(i,j,k)+pyy_d(i,j,k)))
            f5_d(i,j,k+1)=fpc*rhoA_d(i,j,k)/rtot 
            g5_d(i,j,k+1)=fpc*rhob_d(i,j,k)/rtot 
            !6
            feq=p1_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq1_d*((1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*(pxx_d(i,j,k)+pyy_d(i,j,k)))
            f6_d(i,j,k-1)=fpc*rhoA_d(i,j,k)/rtot 
            g6_d(i,j,k-1)=fpc*rhob_d(i,j,k)/rtot 
            !15
            udotc=(ushifted+wshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            fpc=feq +(1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pyy_d(i,j,k)+2.0_db*pxz_d(i,j,k)) 
            f15_d(i+1,j,k+1)=fpc*rhoA_d(i,j,k)/rtot 
            g15_d(i+1,j,k+1)=fpc*rhob_d(i,j,k)/rtot 
            !16
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pyy_d(i,j,k)+2.0_db*pxz_d(i,j,k)) 
            f16_d(i-1,j,k-1)=fpc*rhoA_d(i,j,k)/rtot 
            g16_d(i-1,j,k-1)=fpc*rhob_d(i,j,k)/rtot 
            !17
            udotc=(-ushifted+wshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pyy_d(i,j,k)-2.0_db*pxz_d(i,j,k))
            f17_d(i-1,j,k+1)=fpc*rhoA_d(i,j,k)/rtot 
            g17_d(i-1,j,k+1)=fpc*rhob_d(i,j,k)/rtot 
            !18
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pxx_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pyy_d(i,j,k)-2.0_db*pxz_d(i,j,k))
            f18_d(i+1,j,k-1)=fpc*rhoA_d(i,j,k)/rtot
            g18_d(i+1,j,k-1)=fpc*rhob_d(i,j,k)/rtot 
            !11
            udotc=(vshifted+wshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            fpc=feq +  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pyy_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pxx_d(i,j,k)+2.0_db*pyz_d(i,j,k))
            f11_d(i,j+1,k+1)=fpc*rhoA_d(i,j,k)/rtot 
            g11_d(i,j+1,k+1)=fpc*rhob_d(i,j,k)/rtot 
            
            !12
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pyy_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pxx_d(i,j,k)+2.0_db*pyz_d(i,j,k)) 
            f12_d(i,j-1,k-1)=fpc*rhoA_d(i,j,k)/rtot 
            g12_d(i,j-1,k-1)=fpc*rhob_d(i,j,k)/rtot 
            !13
            udotc=(vshifted-wshifted)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            feq=p2_d*(rtot+(temp + udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pyy_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pxx_d(i,j,k)-2.0_db*pyz_d(i,j,k)) 
            f13_d(i,j+1,k-1)=fpc*rhoA_d(i,j,k)/rtot 
            g13_d(i,j+1,k-1)=fpc*rhob_d(i,j,k)/rtot 
            
            !14
            feq=p2_d*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega_d)*pi2cssq2_d*((1.0_db-cssq_d)*pyy_d(i,j,k)+(1.0_db-cssq_d)*pzz_d(i,j,k)-cssq_d*pxx_d(i,j,k)-2.0_db*pyz_d(i,j,k)) 
            f14_d(i,j-1,k+1)=fpc*rhoA_d(i,j,k)/rtot 
            g14_d(i,j-1,k+1)=fpc*rhob_d(i,j,k)/rtot 
        
        endif
        
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
          !
          !
          !
          g18_d(i+1,j,k-1)=g17_d(i,j,k) !gpc 
          g17_d(i-1,j,k+1)=g18_d(i,j,k) !hpc

          g16_d(i-1,j,k-1)=g15_d(i,j,k) !gpc 
          g15_d(i+1,j,k+1)=g16_d(i,j,k) !hpc

          g14_d(i,j-1,k+1)=g13_d(i,j,k)!gpc 
          g13_d(i,j+1,k-1)=g14_d(i,j,k)!hpc
        
          g12_d(i,j-1,k-1)=g11_d(i,j,k)!gpc 
          g11_d(i,j+1,k+1)=g12_d(i,j,k)!hpc

          g10_d(i-1,j+1,k)=g9_d(i,j,k)!gpc 
          g9_d(i+1,j-1,k)=g10_d(i,j,k)!hpc

          g8_d(i-1,j-1,k)=g7_d(i,j,k)!gpc 
          g7_d(i+1,j+1,k)=g8_d(i,j,k)!hpc

          g6_d(i,j,k-1)=g5_d(i,j,k)!gpc 
          g5_d(i,j,k+1)=g6_d(i,j,k)!hpc 

          g4_d(i,j-1,k)=g3_d(i,j,k)!gpc 
          g3_d(i,j+1,k)=g4_d(i,j,k)!hpc 

          g2_d(i-1,j,k)=g1_d(i,j,k)!gpc 
          g1_d(i+1,j,k)=g2_d(i,j,k)!hpc


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
        rhoprint_d(i,j,k)=rhoA_d(i,j,k)
        velprint_d(1,i,j,k)=u_d(i,j,k)
        velprint_d(2,i,j,k)=v_d(i,j,k)
        velprint_d(3,i,j,k)=w_d(i,j,k)
        
      else
      
        rhoprint_d(i,j,k)=0.0_db
        velprint_d(1,i,j,k)=0.0_db
        velprint_d(2,i,j,k)=0.0_db
        velprint_d(3,i,j,k)=0.0_db
      
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
        if(tmp< 10.0_db )exit
          i=i+1
          tmp=tmp/ 10.0_db
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
  
    endsubroutine get_memory_gpu
    
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
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu1,one_ov_nu2,cssq,fx,fy,fz,temp
    real(kind=4)  :: p1cg,p2cg,p3cg
    real(kind=4) :: radius

    real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
    real(kind=db) :: pi2cssq1,pi2cssq2,pi2cssq0,myrhoA,myrhoB,myu,myv,myw
    real(kind=db) :: b0, b1, b2, beta, sigma
    real(kind=db) :: max_press_excess
    
    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    integer :: TILE_DIMx,TILE_DIMy,TILE_DIMz,TILE_DIM,istat,iframe
    logical :: lprint,lvtk,lpbc,lasync
    real(kind=db) :: mymemory,totmemory
    integer(kind=cuda_Stream_Kind) :: stream1,stream2
    type (cudaDeviceProp) :: prop
    type (cudaEvent) :: startEvent, stopEvent, dummyEvent, dummyEvent1, dummyEvent2
    
       
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
    istat = cudaGetDeviceCount(ngpus)
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif

    !*******************************user parameters and allocations**************************
        nx=64
        ny=64
        nz=64
        myu=0.0_db
        myv=0.0_db
        myw=0.0_db
        nsteps=1000
        stamp=100
        lprint=.false.
        lvtk=.false.
        lpbc=.false.
        lasync=.false.
        fx=0.0_db*10.0**(-4.0_db)
        fy=0.0_db*10.0**(-4.0_db)
        fz=0.0_db*10.0**(-5.0_db)
        TILE_DIMx=32
        TILE_DIMy=8
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

        allocate(isfluid(1:nx,1:ny,1:nz)) 

        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)

        p0=(1.0_db/3.0_db)
        p1=(1.0_db/18.0_db)
        p2=(1.0_db/36.0_db)
        p1dcssq=p1/cssq
        p2dcssq=p2/cssq
    !*****************************************geometry************************
        isfluid=1
        isfluid(1,:,:)=0 !left
        isfluid(nx,:,:)=0 !right
        isfluid(:,1,:)=0 !front 
        isfluid(:,ny,:)=0 !rear
        isfluid(:,:,1)=0 !bottom
        isfluid(:,:,nz)=0 !top
        myu = 0.0_db
        myv = 0.0_db
        myw = 0.0_db
        myrhoA = 1.0_db
        myrhoB = 0.0_db     
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
  !**************************************color gradient********************
        beta=0.95_db
        sigma=0.03_db
        b0=-2.0_db/9.0_db
        b1=1.0_db/54.0_db
        b2=1.0_db/27.0_db
        p1cg=(2.0_db/3.0_db)**2 * (1.0_db/6.0_db)
        p2cg=(1.0_db/6.0_db)**2 * (2.0_db/3.0_db)
        p3cg=(1.0_db/6.0_db)**3
        max_press_excess=0.0
  !************************************copies on gpu**************************
        nx_d=nx
        ny_d=ny
        nz_d=nz
        TILE_DIMx_d=TILE_DIMx
        TILE_DIMy_d=TILE_DIMy
        TILE_DIMz_d=TILE_DIMz
        TILE_DIM_d=TILE_DIM
        myrhoA_d=myrhoA
        myrhoB_d=myrhoB
        myu_d=myu
        myv_d=myv
        myw_d=myw
        fx_d=fx
        fy_d=fy
        fz_d=fz
        p0_d=p0
        p1_d=p1
        p2_d=p2
        one_ov_nu1_d=one_ov_nu1
        one_ov_nu2_d=one_ov_nu2
        omega_d=omega
        cssq_d=cssq
        p1dcssq_d=p1dcssq
        p2dcssq_d=p2dcssq
        pi2cssq0_d=pi2cssq0
        pi2cssq1_d=pi2cssq1
        pi2cssq2_d=pi2cssq2
        qxx_d=qxx
        qyy_d=qyy
        qzz_d=qzz
        qxy_7_8_d=qxy_7_8
        qxy_9_10_d=qxy_9_10
        qxz_15_16_d=qxz_15_16
        qxz_17_18_d=qxz_17_18
        qyz_11_12_d=qyz_11_12
        qyz_13_14_d=qyz_13_14
        sigma_d=sigma
        beta_d=beta
        b0_d=b0
        b1_d=b1
        b2_d=b2
        p1cg_d=p1cg
        p2cg_d=p1cg
        p3cg_d=p1cg
        max_press_excess_d=max_press_excess
        !*****************************************************************************
        allocate(isfluid_d(1:nx_d,1:ny_d,1:nz_d))
        istat = cudaDeviceSynchronize
        istat = cudaMemcpy(isfluid_d,isfluid,nx*ny*nz )
        istat = cudaDeviceSynchronize
        allocate(f0_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f1_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f2_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f3_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(f4_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f5_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f6_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f7_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(f8_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f9_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f10_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f11_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(f12_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f13_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f14_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f15_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(f16_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f17_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),f18_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(g0_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g1_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g2_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g3_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(g4_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g5_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g6_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g7_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(g8_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g9_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g10_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g11_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(g12_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g13_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g14_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g15_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(g16_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g17_d(0:nx_d+1,0:ny_d+1,0:nz_d+1),g18_d(0:nx_d+1,0:ny_d+1,0:nz_d+1))
        allocate(rhoA_d(1:nx_d,1:ny_d,1:nz_d),rhoB_d(1:nx_d,1:ny_d,1:nz_d),u_d(1:nx_d,1:ny_d,1:nz_d),v_d(1:nx_d,1:ny_d,1:nz_d),w_d(1:nx_d,1:ny_d,1:nz_d), &
         pxx_d(1:nx_d,1:ny_d,1:nz_d),pyy_d(1:nx_d,1:ny_d,1:nz_d),pzz_d(1:nx_d,1:ny_d,1:nz_d), &
         pxy_d(1:nx_d,1:ny_d,1:nz_d),pxz_d(1:nx_d,1:ny_d,1:nz_d),pyz_d(1:nx_d,1:ny_d,1:nz_d), psi_d(1:nx_d,1:ny_d,1:nz_d),nci_loc_d(1:nx_d,1:ny_d,1:nz_d))
        istat = cudaDeviceSynchronize
    !*************************************initial conditions ************************    
        call setup_pops<<<dimGrid,dimBlock>>>(radius)
        call abortOnLastErrorAndSync('after setup_pops', step)
        istat = cudaDeviceSynchronize
        
        allocate(rhoprint(1:nx,1:ny,1:nz),velprint(3,1:nx,1:ny,1:nz))
        allocate(rhoprint_d(1:nx_d,1:ny_d,1:nz_d), velprint_d(3,1:nx_d,1:ny_d,1:nz_d))
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
        write(6,*) 'max fy',huge(fy)
        write(6,*) 'max fz',huge(fz)
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
          istat = cudaEventRecord(dummyEvent1, stream1)
          istat = cudaEventSynchronize(dummyEvent1)
          call abortOnLastErrorAndSync('after moments', step)
          
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
        do step=1,nsteps 
        !***********************************moments collision bbck + forcing************************ 

            call moments<<<dimGrid,dimBlock,0,stream1>>>()
            istat = cudaEventRecord(dummyEvent1, stream1)
            istat = cudaEventSynchronize(dummyEvent1)
            call abortOnLastErrorAndSync('after moments', step)
            
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
            
            istat = cudaEventRecord(dummyEvent, stream1)
            istat = cudaEventSynchronize(dummyEvent)
        
        enddo 
    if(lasync)then
        istat = cudaEventRecord(dummyEvent2, stream2)
        istat = cudaEventSynchronize(dummyEvent2)
        if(lvtk)then
          call print_vtk_sync
        else
          call print_raw_sync
        endif
    endif
    istat = cudaDeviceSynchronize
    call cpu_time(ts2)
    call store_print<<<dimGrid,dimBlock>>>()
    istat = cudaDeviceSynchronize
    istat = cudaMemcpy(rhoprint,rhoprint_d,nx*ny*nz )
    istat = cudaMemcpy(velprint,velprint_d,3*nx*ny*nz )
    istat = cudaDeviceSynchronize
    

    write(6,*) 'u=',velprint(1,nx/2,ny/2,nz/2),'v=',velprint(2,nx/2,ny/2,nz/2),'w=',velprint(3,nx/2,ny/2,nz/2),'rho=',rhoprint(nx/2,ny/2,nz/2)
    write(6,*) 'u=',velprint(1,nx/2,ny/2,1),'v=',velprint(2,nx/2,ny/2,1),'w=',velprint(3,nx/2,ny/2,1),'rho=',rhoprint(nx/2,ny/2,1)
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)*real(1.d-9,kind=db)/(ts2-ts1)
    
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)

  contains
  
    subroutine printDeviceProperties(prop,iu,num)
  
      use cudafor
      type(cudadeviceprop) :: prop
      integer,intent(in) :: iu,num 
      !*****************************GPU FEATURE MONITOR*******************************
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
