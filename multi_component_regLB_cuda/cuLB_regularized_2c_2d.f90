#define HACK
module mysubs

   use cudafor

   integer, parameter :: db = 4 !kind(1.0)
   real(kind=db), parameter :: Pi = real(3.141592653589793238462643383279502884d0,kind=db)
   ! device arrays
   integer(kind=4), allocatable, dimension(:, :), device   :: isfluid_d
   integer, constant :: nx_d, ny_d, TILE_DIMx_d, TILE_DIMy_d, TILE_DIM_d
   integer, parameter :: nz = 1, nz_d = 1
   real(kind=db), dimension(0:8), constant :: p_d
   real(kind=db), constant :: fx_d, fy_d, qxx_d, qyy_d, qxy5_7_d, qxy6_8_d
   real(kind=db), device :: omega_d
   real(kind=db), constant :: pi2cssq0_d, pi2cssq1_d, pi2cssq2_d, myrhoA_d, myrhoB_d, myu_d, myv_d, cssq_d
   real(kind=db), constant :: one_ov_nu1_d, one_ov_nu2_d, sigma_d, b0_d, b1_d, b2_d, beta_d, max_press_excess_d
   integer*4, allocatable, dimension(:,:), device :: nci_loc_d
   real(kind=db), allocatable, dimension(:, :), device  :: rhoA_d, rhoB_d, psi_d, u_d, v_d, pxx_d, pyy_d, pxy_d
   real(kind=db), allocatable, dimension(:, :), device  :: f0_d, f1_d, f2_d, f3_d, f4_d, f5_d, f6_d, f7_d, f8_d
   real(kind=db), allocatable, dimension(:, :), device  :: g0_d, g1_d, g2_d, g3_d, g4_d, g5_d, g6_d, g7_d, g8_d
   real(kind=db), allocatable, dimension(:, :, :), device :: rhoprint_d
   real(kind=db), allocatable, dimension(:, :, :, :), device :: velprint_d
   

   contains
   
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
      
            implicit none
      
            real(kind=db), value :: myradius
            integer :: i, j
            real(kind=db) :: mydist,locpsi,locrhoA,locrhoB,tempr
            
            i = (blockIdx%x - 1)*TILE_DIMx_d + threadIdx%x
            j = (blockIdx%y - 1)*TILE_DIMy_d + threadIdx%y
            
            mydist=sqrt((real(i)-(real(nx_d)*0.5_db-myradius-5.0_db))**2.0+(real(j)-real(ny_d)*0.5_db)**2.0)
            
            tempr = fcut(mydist, myradius, myradius+4.0_db)
            
            locpsi=1.0_db
            if (mydist<=myradius)locpsi=-1.0_db
            
            !locrhoB = tempr
            !locrhoA = (1.0 - tempr)
            !locpsi = (locrhoA - locrhoB)/(locrhoA + locrhoB)
            
            locrhoB=0.5_db*(1.0_db-locpsi)
            locrhoA=1.0_db-locrhoB
            
            !if(i==nx_d/2 .and. j==ny_d/2)write(*,*)'cazzo',i,j
            
            u_d(i,j)=0.0_db
            v_d(i,j)=0.0_db

            !write(*,*)i,j,p_d(0)*myrho_d
            f0_d(i, j) = p_d(0)*locrhoA
            f1_d(i, j) = p_d(1)*locrhoA
            f2_d(i, j) = p_d(2)*locrhoA
            f3_d(i, j) = p_d(3)*locrhoA
            f4_d(i, j) = p_d(4)*locrhoA
            f5_d(i, j) = p_d(5)*locrhoA
            f6_d(i, j) = p_d(6)*locrhoA
            f7_d(i, j) = p_d(7)*locrhoA
            f8_d(i, j) = p_d(8)*locrhoA

            g0_d(i, j) = p_d(0)*locrhoB
            g1_d(i, j) = p_d(1)*locrhoB
            g2_d(i, j) = p_d(2)*locrhoB
            g3_d(i, j) = p_d(3)*locrhoB
            g4_d(i, j) = p_d(4)*locrhoB
            g5_d(i, j) = p_d(5)*locrhoB
            g6_d(i, j) = p_d(6)*locrhoB
            g7_d(i, j) = p_d(7)*locrhoB
            g8_d(i, j) = p_d(8)*locrhoB
            
            rhoA_d(i,j)=locrhoA
            rhob_d(i,j)=locrhoB
            psi_d(i,j)=locpsi
            
            return

      end subroutine setup_pops

      attributes(global) subroutine moments()
    
         implicit none
    
         integer :: i, j
         #ifdef HACK
            real(kind=db) ::uu, udotc, temp, fneq, locpxx, locpyy,locpxy,rtot
         #else
            real(kind=db) ::uu, udotc, temp, fneq1, fneq2, fneq3, fneq4, fneq5, fneq6, fneq7, fneq8, rtot
         #endif
         i = (blockIdx%x - 1)*TILE_DIMx_d + threadIdx%x
         j = (blockIdx%y - 1)*TILE_DIMy_d + threadIdx%y

         if (isfluid_d(i, j) .ne. 1) return
         
         rhoA_d(i, j) = f0_d(i, j) + f1_d(i, j) + f2_d(i, j) + f3_d(i, j) + f4_d(i, j) + f5_d(i, j) + f6_d(i, j) + f7_d(i, j) + f8_d(i, j)
         rhoB_d(i, j) = g0_d(i, j) + g1_d(i, j) + g2_d(i, j) + g3_d(i, j) + g4_d(i, j) + g5_d(i, j) + g6_d(i, j) + g7_d(i, j) + g8_d(i, j)
         u_d(i, j) = (f1_d(i, j) + f5_d(i, j) + f8_d(i, j) - f3_d(i, j) - f6_d(i, j) - f7_d(i, j)) + &
                     (g1_d(i, j) + g5_d(i, j) + g8_d(i, j) - g3_d(i, j) - g6_d(i, j) - g7_d(i, j))
         v_d(i, j) = (f5_d(i, j) + f2_d(i, j) + f6_d(i, j) - f7_d(i, j) - f4_d(i, j) - f8_d(i, j)) + &
                     (g5_d(i, j) + g2_d(i, j) + g6_d(i, j) - g7_d(i, j) - g4_d(i, j) - g8_d(i, j))

         psi_d(i, j) = (rhoA_d(i, j) - rhoB_d(i, j))/(rhoA_d(i, j) + rhoB_d(i, j))

         rtot = rhoA_d(i, j) + rhoB_d(i, j)

         ! non equilibrium pressor components
         uu = 0.5_db*(u_d(i, j)*u_d(i, j) + v_d(i, j)*v_d(i, j))/cssq_d
         #ifdef HACK
            !1-3
            udotc = u_d(i, j)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            fneq = (f1_d(i, j) + g1_d(i, j)) - p_d(1)*(rtot + (temp + udotc))
            locpxx = fneq
            fneq = (f3_d(i, j) + g3_d(i, j)) - p_d(3)*(rtot + (temp - udotc))
            locpxx = fneq + locpxx
            
            !2-4
            udotc = v_d(i, j)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            fneq = (f2_d(i, j) + g2_d(i, j)) - p_d(2)*(rtot + (temp + udotc))
            locpyy = fneq
            fneq = (f4_d(i, j) + g2_d(i, j)) - p_d(4)*(rtot + (temp - udotc))
            locpyy = fneq + locpyy
            !5-7
            udotc = (u_d(i, j) + v_d(i, j))/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            fneq = (f5_d(i, j) + g5_d(i, j)) - p_d(5)*(rtot + (temp + udotc))
            locpxx = fneq + locpxx
            locpyy = fneq + locpyy
            locpxy = fneq
            fneq = (f7_d(i, j) + g7_d(i, j)) - p_d(7)*(rtot + (temp - udotc))
            locpxx = fneq + locpxx
            locpyy = fneq + locpyy
            locpxy = fneq + locpxy
            !6-8
            udotc = (-u_d(i, j) + v_d(i, j))/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            fneq = (f6_d(i, j) + g6_d(i, j)) - p_d(6)*(rtot + (temp + udotc))
            locpxx = fneq + locpxx
            locpyy = fneq + locpyy
            locpxy = -fneq + locpxy
            fneq = (f8_d(i, j) + g8_d(i, j)) - p_d(8)*(rtot + (temp - udotc))
            locpxx = fneq + locpxx
            locpyy = fneq + locpyy
            locpxy = -fneq + locpxy
            
            pxx_d(i, j) = locpxx
            pyy_d(i, j) = locpyy
            pxy_d(i, j) = locpxy
         #else
            !1-3
            udotc = u_d(i, j)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            fneq1 = (f1_d(i, j) + g1_d(i, j)) - p_d(1)*(rtot + (temp + udotc))
            fneq3 = (f3_d(i, j) + g3_d(i, j)) - p_d(3)*(rtot + (temp - udotc))
            !2-4
            udotc = v_d(i, j)/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            fneq2 = (f2_d(i, j) + g2_d(i, j)) - p_d(2)*(rtot + (temp + udotc))
            fneq4 = (f4_d(i, j) + g2_d(i, j)) - p_d(4)*(rtot + (temp - udotc))
            !5-7
            udotc = (u_d(i, j) + v_d(i, j))/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            fneq5 = (f5_d(i, j) + g5_d(i, j)) - p_d(5)*(rtot + (temp + udotc))
            fneq7 = (f7_d(i, j) + g7_d(i, j)) - p_d(7)*(rtot + (temp - udotc))
            !6-8
            udotc = (-u_d(i, j) + v_d(i, j))/cssq_d
            temp = -uu + 0.5_db*udotc*udotc
            fneq6 = (f6_d(i, j) + g6_d(i, j)) - p_d(6)*(rtot + (temp + udotc))
            fneq8 = (f8_d(i, j) + g8_d(i, j)) - p_d(8)*(rtot + (temp - udotc))

            pxx_d(i, j) = fneq1 + fneq3 + fneq5 + fneq6 + fneq7 + fneq8
            pyy_d(i, j) = fneq2 + fneq4 + fneq5 + fneq6 + fneq7 + fneq8
            pxy_d(i, j) = fneq5 - fneq6 + fneq7 - fneq8
         #endif
      end subroutine moments

      attributes(global) subroutine streamcoll()
      
            implicit none
            
            integer :: i, j
            real(kind=db) :: uu, udotc, temp, feq, psi_x, psi_y, rtot, st_coeff, mod_psi, mod_psi_sq, fpc
            real(kind=db) :: norm_x, norm_y, rprod   
         #ifndef HACK   
            real(kind=db) :: addendum0, addendum1, addendum2, addendum3, addendum4, addendum5, addendum6, addendum7, addendum8
            real(kind=db) :: gaddendum0, gaddendum1, gaddendum2, gaddendum3, gaddendum4, gaddendum5, gaddendum6, gaddendum7, gaddendum8
         #endif
            real(kind=db) :: ushifted, vshifted,nu_avg

            i = (blockIdx%x - 1)*TILE_DIMx_d + threadIdx%x
            j = (blockIdx%y - 1)*TILE_DIMy_d + threadIdx%y
      
         if (isfluid_d(i, j) .ne. 1) return
         
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
            gaddendum0=0.0_db
            if(mod_psi>0.0001)then ! i'm sitting on the interface 
               norm_x=psi_x/mod_psi
               norm_y=psi_y/mod_psi
               ushifted=u(i,j) + fx + float(nci_loc(i,j))*(norm_x)*max_press_excess*abs(rhob(i,j))
               vshifted=v(i,j) + fy + float(nci_loc(i,j))*(norm_y)*max_press_excess*abs(rhob(i,j))

               addendum0=-st_coeff*mod_psi*b0
               uu=0.5_db*(ushifted*ushifted + vshifted*vshifted)/cssq
               feq=p(0)*(rtot-uu)
               fpc=feq + (1.0_db-omega)*pi2cssq0*(- cssq*pyy(i,j)-cssq*pxx(i,j))  + addendum0
               f0(i,j)=fpc*(rhoA(i,j))/rtot 
               g0(i,j)=fpc*(rhoB(i,j))/rtot

               addendum0=st_coeff*mod_psi*(p(1)*psi_x**2/mod_psi_sq - b1)
               gaddendum0=p(1)*(rtot)*(rprod*beta*psi_x/mod_psi/rtot**2)
               !1-3
               udotc=ushifted/cssq
               temp = -uu + 0.5_db*udotc*udotc
               feq=p(1)*(rtot+(temp + udotc))
               fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j) - cssq*pyy(i,j) )+ addendum0 !+ (fx+float(nci_loc(i,j))*(norm_x)*max_press_excess*abs(rhoa(i,j)))*p(1)/cssq 
               f1(i+1,j)= fpc*(rhoA(i,j))/rtot + gaddendum0
               g1(i+1,j)= fpc*(rhoB(i,j))/rtot - gaddendum0
               addendum0=st_coeff*mod_psi*(p(3)*psi_x**2/mod_psi_sq - b1)
               gaddendum0=p(3)*(rtot)*(rprod*beta*(-psi_x/mod_psi)/rtot**2)
               feq=p(3)*(rtot+(temp - udotc))
               fpc=feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j) - cssq*pyy(i,j)) + addendum0 !- (fx+float(nci_loc(i,j))*(norm_x)*max_press_excess*abs(rhoa(i,j)))*p(3)/cssq 
               f3(i-1,j)= fpc*(rhoA(i,j))/rtot + gaddendum0 
               g3(i-1,j)= fpc*(rhoB(i,j))/rtot - gaddendum0 

               !2-4
               addendum0=st_coeff*mod_psi*(p(2)*psi_y**2/mod_psi_sq - b1)
               gaddendum0=p(2)*(rtot)*(rprod*beta*psi_y/mod_psi/rtot**2)
               udotc=vshifted/cssq
               temp = -uu + 0.5_db*udotc*udotc
               feq=p(2)*(rtot+(temp + udotc))
               fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j)) + addendum0 !+ (fy+float(nci_loc(i,j))*(norm_y)*max_press_excess*abs(rhoa(i,j)))*p(2)/cssq  !
               f2(i,j+1)= fpc*(rhoA(i,j))/rtot + gaddendum0  
               g2(i,j+1)= fpc*(rhoB(i,j))/rtot - gaddendum0
               addendum0=st_coeff*mod_psi*(p(4)*psi_y**2/mod_psi_sq - b1)
               gaddendum0=p(4)*(rtot)*(rprod*beta*(-psi_y/mod_psi)/rtot**2)
               feq=p(4)*(rtot+(temp - udotc))
               fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j)) + addendum0 !- (fy+float(nci_loc(i,j))*(norm_y)*max_press_excess*abs(rhoa(i,j)))*p(4)/cssq 
               f4(i,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum0 
               g4(i,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum0 

               !5-7
               addendum0=st_coeff*mod_psi*(p(5)*(psi_x+psi_y)**2/mod_psi_sq - b2)
               gaddendum0=p(5)*(rtot)*(rprod*beta*(psi_x/mod_psi + psi_y/mod_psi)/rtot**2)
               udotc=(ushifted+vshifted)/cssq
               temp = -uu + 0.5_db*udotc*udotc
               feq=p(5)*(rtot+(temp + udotc))
               fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + addendum0 !+ (fx + fy + float(nci_loc(i,j))*(norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(5)/cssq 
               f5(i+1,j+1)=fpc*(rhoA(i,j))/rtot + gaddendum0
               g5(i+1,j+1)=fpc*(rhoB(i,j))/rtot - gaddendum0
               addendum0=st_coeff*mod_psi*(p(7)*(-psi_x-psi_y)**2/mod_psi_sq - b2)
               gaddendum0=p(7)*(rtot)*(rprod*beta*(-psi_x/mod_psi - psi_y/mod_psi)/rtot**2)
               feq=p(7)*(rtot+(temp - udotc))
               fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + addendum0 !- (fx + fy + float(nci_loc(i,j))*(norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(7)/cssq 
               f7(i-1,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum0
               g7(i-1,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum0

               !6-8
               addendum0=st_coeff*mod_psi*(p(6)*(-psi_x+psi_y)**2/mod_psi_sq - b2)
               gaddendum0=p(6)*(rtot)*(rprod*beta*(-psi_x/mod_psi + psi_y/mod_psi)/rtot**2)
               udotc=(-ushifted+vshifted)/cssq
               temp = -uu + 0.5_db*udotc*udotc
               feq=p(6)*(rtot+(temp + udotc))
               fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) + addendum0 !+(-fx + fy + float(nci_loc(i,j))*(-norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(6)/cssq 
               f6(i-1,j+1)= fpc*(rhoA(i,j))/rtot + gaddendum0
               g6(i-1,j+1)= fpc*(rhoB(i,j))/rtot - gaddendum0
               addendum0=st_coeff*mod_psi*(p(8)*(psi_x-psi_y)**2/mod_psi_sq - b2)
               gaddendum0=p(8)*(rtot)*(rprod*beta*(psi_x/mod_psi - psi_y/mod_psi)/rtot**2)
               feq=p(8)*(rtot+(temp - udotc))
               fpc=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j))  + addendum0 !+( fx - fy + float(nci_loc(i,j))*(norm_x-norm_y)*max_press_excess*abs(rhoa(i,j)))*p(8)/cssq 
               f8(i+1,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum0 
               g8(i+1,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum0 
               
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
            fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j) - cssq*pyy(i,j) )+ addendum0 !+ (fx+float(nci_loc(i,j))*(norm_x)*max_press_excess*abs(rhoa(i,j)))*p(1)/cssq 
            f1(i+1,j)= fpc*(rhoA(i,j))/rtot + gaddendum0
            g1(i+1,j)= fpc*(rhoB(i,j))/rtot - gaddendum0
            !3
            feq=p(3)*(rtot+(temp - udotc))
            fpc=feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j) - cssq*pyy(i,j)) + addendum0 !- (fx+float(nci_loc(i,j))*(norm_x)*max_press_excess*abs(rhoa(i,j)))*p(3)/cssq 
            f3(i-1,j)= fpc*(rhoA(i,j))/rtot + gaddendum0
            g3(i-1,j)= fpc*(rhoB(i,j))/rtot - gaddendum0 
            !2
            udotc=vshifted/cssq
            temp = -uu + 0.5_db*udotc*udotc
            feq=p(2)*(rtot+(temp + udotc))
            fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j)) + addendum0 !+ (fy+float(nci_loc(i,j))*(norm_y)*max_press_excess*abs(rhoa(i,j)))*p(2)/cssq  !
            f2(i,j+1)= fpc*(rhoA(i,j))/rtot + gaddendum0  
            g2(i,j+1)= fpc*(rhoB(i,j))/rtot - gaddendum0   
            !4
            feq=p(4)*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j)) + addendum0 !- (fy+float(nci_loc(i,j))*(norm_y)*max_press_excess*abs(rhoa(i,j)))*p(4)/cssq 
            f4(i,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum0 
            g4(i,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum0 
            !5
            udotc=(ushifted+vshifted)/cssq
            temp = -uu + 0.5_db*udotc*udotc
            feq=p(5)*(rtot+(temp + udotc))
            fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + addendum0 !+ (fx + fy + float(nci_loc(i,j))*(norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(5)/cssq 
            f5(i+1,j+1)=fpc*(rhoA(i,j))/rtot + gaddendum0
            g5(i+1,j+1)=fpc*(rhoB(i,j))/rtot - gaddendum0
            !7
            feq=p(7)*(rtot+(temp - udotc))
            fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + addendum0 !- (fx + fy + float(nci_loc(i,j))*(norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(7)/cssq 
            f7(i-1,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum0
            g7(i-1,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum0
            !6
            udotc=(-ushifted+vshifted)/cssq
            temp = -uu + 0.5_db*udotc*udotc
            feq=p(6)*(rtot+(temp + udotc))
            fpc=feq  + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) + addendum0 !+(-fx + fy + float(nci_loc(i,j))*(-norm_x+norm_y)*max_press_excess*abs(rhoa(i,j)))*p(6)/cssq 
            f6(i-1,j+1)= fpc*(rhoA(i,j))/rtot + gaddendum0
            g6(i-1,j+1)= fpc*(rhoB(i,j))/rtot - gaddendum0
            !8
            feq=p(8)*(rtot+(temp - udotc))
            fpc=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j))  + addendum0 !+( fx - fy + float(nci_loc(i,j))*(norm_x-norm_y)*max_press_excess*abs(rhoa(i,j)))*p(8)/cssq 
            f8(i+1,j-1)=fpc*(rhoA(i,j))/rtot + gaddendum0 
            g8(i+1,j-1)=fpc*(rhoB(i,j))/rtot - gaddendum0 
         endif
      end subroutine streamcoll

      attributes(global) subroutine bcs_no_slip()
      
            implicit none
            
            integer :: i, j

            i = (blockIdx%x - 1)*TILE_DIMx_d + threadIdx%x
            j = (blockIdx%y - 1)*TILE_DIMy_d + threadIdx%y

            !write(*,*)i,j,p_d(0)*myrho_d
            if (isfluid_d(i, j) .ne. 0) return
            
            psi_d(i,j)=-1.0_db
         

            f8_d(i + 1, j - 1) = f6_d(i, j)!gpc
            f7_d(i - 1, j - 1) = f5_d(i, j)!hpc

            f6_d(i - 1, j + 1) = f8_d(i, j)!gpc
            f5_d(i + 1, j + 1) = f7_d(i, j)!hpc

            f4_d(i, j - 1) = f2_d(i, j)!gpc
            f3_d(i - 1, j) = f1_d(i, j)!hpc

            f2_d(i, j + 1) = f4_d(i, j)!gpc
            f1_d(i + 1, j) = f3_d(i, j)!hpc
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            g8_d(i + 1, j - 1) = g6_d(i, j)!gpc
            g7_d(i - 1, j - 1) = g5_d(i, j)!hpc

            g6_d(i - 1, j + 1) = g8_d(i, j)!gpc
            g5_d(i + 1, j + 1) = g7_d(i, j)!hpc

            g4_d(i, j - 1) = g2_d(i, j)!gpc
            g3_d(i - 1, j) = g1_d(i, j)!hpc

            g2_d(i, j + 1) = g4_d(i, j)!gpc
            g1_d(i + 1, j) = g3_d(i, j)!hpc

      end subroutine bcs_no_slip

      attributes(global) subroutine pbc_edge_x()
            
            implicit none
            
            integer :: i, j

            j = (blockIdx%x - 1)*TILE_DIM_d + threadIdx%x
            
            if (j > ny_d) return

            if (j > 2 .and. j < ny_d - 1) then

               f1_d(2, j) = f1_d(nx_d, j)
               f5_d(2, j) = f5_d(nx_d, j)
               f8_d(2, j) = f8_d(nx_d, j)

               f3_d(nx_d - 1, j) = f3_d(1, j)
               f6_d(nx_d - 1, j) = f6_d(1, j)
               f7_d(nx_d - 1, j) = f7_d(1, j)

            else

               if (j == 2) then
                  f1_d(2, j) = f1_d(nx_d, j)
                  f8_d(2, j) = f8_d(nx_d, j)

                  f3_d(nx_d - 1, j) = f3_d(1, j)
                  f7_d(nx_d - 1, j) = f7_d(1, j)

               end if

               if (j == ny_d - 1) then
                  f1_d(2, j) = f1_d(nx_d, j)
                  f5_d(2, j) = f5_d(nx_d, j)

                  f3_d(nx_d - 1, j) = f3_d(1, j)
                  f6_d(nx_d - 1, j) = f6_d(1, j)
               end if
            end if
      end subroutine pbc_edge_x

      attributes(global) subroutine pbc_edge_y()

               integer :: i

               i = (blockIdx%x - 1)*TILE_DIM_d + threadIdx%x
               if (i > nx_d) return

               if (i > 2 .and. i < nx_d - 1) then

                  f2_d(i, 2) = f2_d(i, ny_d)
                  f5_d(i, 2) = f5_d(i, ny_d)
                  f6_d(i, 2) = f6_d(i, ny_d)

                  f4_d(i, ny_d - 1) = f4_d(i, 1)
                  f7_d(i, ny_d - 1) = f7_d(i, 1)
                  f8_d(i, ny_d - 1) = f8_d(i, 1)

               else

                  if (i == 2) then
                     f2_d(i, 2) = f2_d(i, ny_d)
                     f6_d(i, 2) = f6_d(i, ny_d)

                     f4_d(i, ny_d - 1) = f4_d(i, 1)
                     f7_d(i, ny_d - 1) = f7_d(i, 1)

                  end if

                  if (i == nx_d - 1) then
                     f2_d(i, 2) = f2_d(i, ny_d)
                     f5_d(i, 2) = f5_d(i, ny_d)

                     f4_d(i, ny_d - 1) = f4_d(i, 1)
                     f8_d(i, ny_d - 1) = f8_d(i, 1)

                  end if

               end if

      end subroutine pbc_edge_y

      attributes(global) subroutine store_print()

               integer :: i, j

               i = (blockIdx%x - 1)*TILE_DIMx_d + threadIdx%x
               j = (blockIdx%y - 1)*TILE_DIMy_d + threadIdx%y

               !write(*,*)i,j,p_d(0)*myrho_d
               if (isfluid_d(i, j) .eq. 1) then
                  rhoprint_d(i, j, 1) = rhoA_d(i, j)
                  velprint_d(1, i, j, 1) = u_d(i, j)
                  velprint_d(2, i, j, 1) = v_d(i, j)
                  velprint_d(3, i, j, 1) = 0

               else

                  rhoprint_d(i, j, 1) = 0.0_db
                  velprint_d(1, i, j, 1) = 0.0_db
                  velprint_d(2, i, j, 1) = 0.0_db
                  velprint_d(3, i, j, 1) = 0.0_db

               end if

               return

      end subroutine store_print

end module mysubs
!******************************prints*************************!
module prints

      use mysubs

      implicit none

      integer, parameter :: mxln = 120
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
      character(len=*), parameter :: filenamevtk = 'out'

      real(kind=4), allocatable, dimension(:, :, :) :: rhoprint
      real(kind=4), allocatable, dimension(:, :, :, :) :: velprint
      logical :: lelittle
      character(len=mxln) :: dir_out
      character(len=mxln) :: extentvtk
      character(len=mxln) :: sevt1, sevt2
      character(len=1), allocatable, dimension(:) :: head1, head2

      contains

      subroutine header_vtk(nx, ny, nz, mystring500, namevar, extent, ncomps, iinisub, iend, myoffset, &
                              new_myoffset, indent)

            implicit none

            integer, intent(in) :: nx, ny, nz
            character(len=8), intent(in) :: namevar
            character(len=120), intent(in) :: extent
            integer, intent(in) :: ncomps, iinisub, myoffset
            integer, intent(out) :: iend, new_myoffset
            integer, intent(inout) :: indent

            !namevar='density1'

            character(len=500), intent(out) :: mystring500
            ! End-character for binary-record finalize.
            character(1), parameter:: end_rec = char(10)
            character(1) :: string1
            character(len=*), parameter :: topology = 'ImageData'
            integer :: ioffset, nele, bytechar, byteint, byter4, byter8, iini

            iini = iinisub
            bytechar = kind(end_rec)
            byteint = kind(iini)
            byter4 = 4
            byter8 = 8

            mystring500 = repeat(' ', 500)

            iend = iini

            iini = iend + 1
            nele = 22
            iend = iend + nele
            mystring500(iini:iend) = '<?xml version="1.0"?>'//end_rec

            new_myoffset = myoffset
            new_myoffset = new_myoffset + nele*bytechar

            iini = iend + 1
            nele = 67
            iend = iend + nele
            if (lelittle) then
               mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
                                       '" version="0.1" byte_order="LittleEndian">'//end_rec
            else
               mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
                                       '" version="0.1" byte_order="BigEndian">   '//end_rec
            end if

            new_myoffset = new_myoffset + 67*bytechar

            indent = indent + 2
            iini = iend + 1
            nele = 70
            iend = iend + nele
            mystring500(iini:iend) = repeat(' ', indent)//'<'//trim(topology)//' WholeExtent="'// &
                                    trim(extent)//'">'//end_rec

            new_myoffset = new_myoffset + 70*bytechar

            indent = indent + 2
            iini = iend + 1
            nele = 63
            iend = iend + nele
            mystring500(iini:iend) = repeat(' ', indent)//'<Piece Extent="'//trim(extent)//'">'//end_rec

            new_myoffset = new_myoffset + 63*bytechar

            ! initializing offset pointer
            ioffset = 0

            indent = indent + 2
            iini = iend + 1
            nele = 18
            iend = iend + nele
            mystring500(iini:iend) = repeat(' ', indent)//'<PointData>'//end_rec

            new_myoffset = new_myoffset + 18*bytechar

            indent = indent + 2
            iini = iend + 1
            nele = 115
            iend = iend + nele

            if (ncomps /= 1 .and. ncomps /= 3) then
               write (6, '(a)') 'ERROR in header_vtk'
               stop
            end if
            write (string1, '(i1)') ncomps
            mystring500(iini:iend) = repeat(' ', indent)//'<DataArray type="Float32" Name="'// &
                                    namevar//'" NumberOfComponents="'//string1//'" '// &
                                    'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec

            new_myoffset = new_myoffset + 115*bytechar

            indent = indent - 2
            iini = iend + 1
            nele = 19
            iend = iend + nele
            mystring500(iini:iend) = repeat(' ', indent)//'</PointData>'//end_rec

            new_myoffset = new_myoffset + 19*bytechar

            indent = indent - 2
            iini = iend + 1
            nele = 13
            iend = iend + nele
            mystring500(iini:iend) = repeat(' ', indent)//'</Piece>'//end_rec

            new_myoffset = new_myoffset + 13*bytechar

            indent = indent - 2
            iini = iend + 1
            nele = 15
            iend = iend + nele
            mystring500(iini:iend) = repeat(' ', indent)//'</'//trim(topology)//'>'//end_rec

            new_myoffset = new_myoffset + 15*bytechar

            iini = iend + 1
            nele = 32
            iend = iend + nele
            mystring500(iini:iend) = repeat(' ', indent)//'<AppendedData encoding="raw">'//end_rec

            new_myoffset = new_myoffset + 32*bytechar

            iini = iend + 1
            nele = 1
            iend = iend + nele
            mystring500(iini:iend) = '_'

            new_myoffset = new_myoffset + 1*bytechar

            return

      end subroutine header_vtk

      subroutine footer_vtk(nx, ny, nz, mystring30, iinisub, iend, myoffset, &
                           new_myoffset, indent)

         implicit none

         integer, intent(in) :: nx, ny, nz
         integer, intent(in) :: iinisub, myoffset
         integer, intent(out) :: iend, new_myoffset
         integer, intent(inout) :: indent

         character(len=30), intent(out) :: mystring30
         ! End-character for binary-record finalize.
         character(1), parameter:: end_rec = char(10)
         character(1) :: string1
         character(len=*), parameter :: topology = 'ImageData'
         integer :: ioffset, nele, bytechar, byteint, byter4, byter8, iini

         iini = iinisub
         bytechar = kind(end_rec)
         byteint = kind(iini)
         byter4 = 4
         byter8 = 8

         mystring30 = repeat(' ', 30)

         iend = iini

         iini = iend + 1
         nele = 1
         iend = iend + nele
         mystring30(iini:iend) = end_rec

         new_myoffset = myoffset
         new_myoffset = new_myoffset + 1*bytechar

         iini = iend + 1
         nele = 18
         iend = iend + nele
         mystring30(iini:iend) = repeat(' ', indent)//'</AppendedData>'//end_rec

         new_myoffset = new_myoffset + 18*bytechar

         iini = iend + 1
         nele = 11
         iend = iend + nele
         mystring30(iini:iend) = '</VTKFile>'//end_rec

         if (iend /= 30) then
            write (6, '(a)') 'ERROR in footer_vtk'
            stop
         end if

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

            if (btest(transfer(int((/1, 0, 0, 0/), ik1), 1_ik4), 0)) then
               !it is little endian
               ltest = .true.
            else
               !it is big endian
               ltest = .false.
            end if

            return

      end subroutine test_little_endian

   !!!!!!!!!!!!!!!!!!!!!!init_output!!!!!!!!!!!!!!!!!!!
      subroutine init_output(nx, ny, nz, ncomp, lvtk)

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

            integer, intent(in) :: nx, ny, nz, ncomp
            logical, intent(in) :: lvtk
            character(len=255) :: path, makedirectory
            logical :: lexist

            integer :: i, j, k, nn, indent, myoffset, new_myoffset, iend
            integer, parameter :: byter4 = 4
            integer, parameter :: byteint = 4
            integer, allocatable :: printlistvtk(:)
            integer, parameter :: ioxyz = 54
            character(len=*), parameter :: filexyz = 'isfluid.xyz'
            character(len=120) :: mystring120

            call test_little_endian(lelittle)

            sevt1 = repeat(' ', mxln)
            sevt2 = repeat(' ', mxln)

            path = repeat(' ', 255)
            call getcwd(path)

            !call get_environment_variable('DELIMITER',delimiter)
            path = trim(path)
            delimiter = path(1:1)
            if (delimiter == ' ') delimiter = '/'

            makedirectory = repeat(' ', 255)
            makedirectory = 'output'//delimiter
            dir_out = trim(makedirectory)

            #ifdef _INTEL
               inquire (directory=trim(makedirectory), exist=lexist)
            #else
               inquire (file=trim(makedirectory), exist=lexist)
            #endif

            if (.not. lexist) then
                  makedirectory = repeat(' ', 255)
                  makedirectory = 'mkdir output'
                  call system(makedirectory)
            end if
            mystring120 = repeat(' ', 120)

            makedirectory = repeat(' ', 255)
            makedirectory = trim(path)//delimiter//'output'//delimiter

            extentvtk = space_fmtnumb(1)//' '//space_fmtnumb(nx)//' ' &
                        //space_fmtnumb(1)//' '//space_fmtnumb(ny)//' ' &
                        //space_fmtnumb(1)//' '//space_fmtnumb(nz)

            if (ncomp == 1) then
                  nfilevtk = 2
            elseif (ncomp == 2) then
                  nfilevtk = 3
            end if   
            allocate (printlistvtk(nfilevtk))
            do i = 1, nfilevtk
                  printlistvtk(i) = i
            end do

            allocate (varlistvtk(nfilevtk))
            allocate (namevarvtk(nfilevtk))
            allocate (ndimvtk(nfilevtk))
            allocate (headervtk(nfilevtk))
            allocate (footervtk(nfilevtk))
            allocate (nheadervtk(nfilevtk))
            allocate (vtkoffset(nfilevtk))
            allocate (ndatavtk(nfilevtk))
            varlistvtk(1:nfilevtk) = printlistvtk(1:nfilevtk)

            if (ncomp == 1) then
                  do i = 1, nfilevtk
                        select case (printlistvtk(i))
                        case (1)
                           namevarvtk(i) = 'rho     '
                           ndimvtk(i) = 1
                        case (2)
                           namevarvtk(i) = 'vel     '
                           ndimvtk(i) = 3
                        case default
                           write (6, '(a)') 'ERROR in init_output'
                           stop
                        end select
                  end do
            elseif (ncomp == 2) then
                  do i = 1, nfilevtk
                        select case (printlistvtk(i))
                        case (1)
                           namevarvtk(i) = 'rho1    '
                           ndimvtk(i) = 1
                        case (2)
                           namevarvtk(i) = 'rho2    '
                           ndimvtk(i) = 1
                        case (3)
                           namevarvtk(i) = 'vel     '
                           ndimvtk(i) = 3
                        case default
                           write (6, '(a)') 'ERROR in init_output'
                           stop
                        end select
                  end do
            end if
            nn = nx*ny*nz

            do i = 1, nfilevtk
               myoffset = 0
               indent = 0
               call header_vtk(nx, ny, nz, headervtk(i), namevarvtk(i), extentvtk, ndimvtk(i), 0, iend, myoffset, &
                              new_myoffset, indent)
               vtkoffset(i) = new_myoffset
               myoffset = new_myoffset + byteint + ndimvtk(i)*nn*byter4
               ndatavtk(i) = ndimvtk(i)*nn*byter4
               nheadervtk(i) = iend
               call footer_vtk(nx, ny, nz, footervtk(i), 0, iend, myoffset, &
                              new_myoffset, indent)
            end do

            return

      end subroutine init_output

      subroutine string_char(mychar, nstring, mystring)

         implicit none

         integer :: i
         character(1), allocatable, dimension(:) :: mychar
         integer, intent(in) :: nstring
         character(len=*), intent(in) :: mystring

         allocate (mychar(nstring))

         do i = 1, nstring
            mychar(i) = mystring(i:i)
         end do

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

         integer, intent(in) :: inum
         character(len=6) :: space_fmtnumb
         integer :: numdigit, irest
         real(kind=8) :: tmp
         character(len=22) :: cnumberlabel

         numdigit = dimenumb(inum)
         irest = 6 - numdigit
         if (irest > 0) then
            write (cnumberlabel, "(a,i8,a,i8,a)") "(a", irest, ",i", numdigit, ")"
            write (space_fmtnumb, fmt=cnumberlabel) repeat(' ', irest), inum
         else
            write (cnumberlabel, "(a,i8,a)") "(i", numdigit, ")"
            write (space_fmtnumb, fmt=cnumberlabel) inum
         end if

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

         integer, intent(in) :: inum
         character(len=12) :: space_fmtnumb12
         integer :: numdigit, irest
         real(kind=8) :: tmp
         character(len=22) :: cnumberlabel

         numdigit = dimenumb(inum)
         irest = 12 - numdigit
         if (irest > 0) then
            write (cnumberlabel, "(a,i8,a,i8,a)") "(a", irest, ",i", numdigit, ")"
            write (space_fmtnumb12, fmt=cnumberlabel) repeat(' ', irest), inum
         else
            write (cnumberlabel, "(a,i8,a)") "(i", numdigit, ")"
            write (space_fmtnumb12, fmt=cnumberlabel) inum
         end if

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

         integer, intent(in) :: inum
         integer :: dimenumb
         integer :: i
         real(kind=db) :: tmp

         i = 1
         tmp = real(inum, kind=db)
         do
            if (tmp < 10.0_db) exit
            i = i + 1
            tmp = tmp/10.0_db
         end do

         dimenumb = i

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

         integer, intent(in) :: inum
         character(len=6) :: write_fmtnumb
         integer :: numdigit, irest
         !real*8 :: tmp
         character(len=22) :: cnumberlabel

         numdigit = dimenumb(inum)
         irest = 6 - numdigit
         if (irest > 0) then
            write (cnumberlabel, "(a,i8,a,i8,a)") "(a", irest, ",i", numdigit, ")"
            write (write_fmtnumb, fmt=cnumberlabel) repeat('0', irest), inum
         else
            write (cnumberlabel, "(a,i8,a)") "(i", numdigit, ")"
            write (write_fmtnumb, fmt=cnumberlabel) inum
         end if

         return
      end function write_fmtnumb

      subroutine get_memory_gpu(fout, fout2)

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

               real(kind=db), intent(out) :: fout, fout2
               real(kind=db) :: myd(2), myd2(2)
               integer :: istat
         #ifdef _OPENACC
               integer :: myfree, total
         #elif defined _CUDA
               integer(kind=cuda_count_kind) :: myfree, total
         #else
               integer :: myfree, total
         #endif

         #ifdef _OPENACC
               myfree = acc_get_free_memory()
               total = acc_get_memory()
         #elif defined _CUDA
               istat = cudaMemGetInfo(myfree, total)
         #else
               myfree = 0
               total = 0
         #endif
         fout = real(total - myfree, kind=4)/(1024.0**3.0)
         fout2 = real(total, kind=4)/(1024.0**3.0)

         return

      end subroutine get_memory_gpu

      subroutine print_memory_registration_gpu(iu, mybanner, mybanner2, &
                                             mymemory, totmem)

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
         character(len=*), intent(in) :: mybanner, mybanner2
         real(kind=db), intent(in) :: mymemory, totmem

         character(len=12) :: r_char, r_char2

         character(len=*), parameter :: of = '(a)'

         write (r_char, '(f12.4)') mymemory
         write (r_char2, '(f12.4)') totmem
         write (iu, of) "                                                                               "
         write (iu, of) "******************************GPU MEMORY MONITOR*******************************"
         write (iu, of) "                                                                               "
         write (iu, '(4a)') trim(mybanner), " = ", trim(adjustl(r_char)), " (GB)"
         write (iu, '(4a)') trim(mybanner2), " = ", trim(adjustl(r_char2)), " (GB)"
         write (iu, of) "                                                                               "
         write (iu, of) "*******************************************************************************"
         write (iu, of) "                                                                               "

         return

      end subroutine print_memory_registration_gpu

end module

program lb_openacc

   use cudafor
   use mysubs
   use prints

   implicit none

   integer(kind=4) :: i, j, ll, l, dumm
   integer(kind=4) :: nx, ny, step, stamp, nlinks, nsteps, ngpus
   integer :: TILE_DIMx, TILE_DIMy, TILE_DIM, istat, iframe
   real(kind=db), parameter :: pi_greek = 3.141592653589793238462643383279502884_db
   logical :: lprint = .false.
   logical :: lvtk = .false.
   logical :: lpbc = .false.
   logical :: lasync = .false.
   real(kind=4)  :: ts1, ts2, time
   real(kind=db) :: visc_LB, uu, udotc, omega, feq
   !real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8
   real(kind=db) :: qxx, qyy, qxy5_7, qxy6_8, pi2cssq1, pi2cssq2, pi2cssq0
   real(kind=db) :: tau, one_ov_nu, cssq, fx, fy, temp, dummy, myrhoA,myrhoB, myu, myv
   real(kind=db) :: b0, b1, b2, beta, sigma
   real(kind=db) :: one_ov_nu2, one_ov_nu1
   real(kind=db) :: max_press_excess, ushifted, vshifted

   integer(kind=4), allocatable, dimension(:, :)   :: isfluid

   real(kind=db), allocatable, dimension(:)     :: p

   real(kind=db) :: mymemory, totmemory,radius
   integer(kind=cuda_Stream_Kind) :: stream1, stream2
   type(cudaDeviceProp) :: prop
   type(cudaEvent) :: startEvent, stopEvent, dummyEvent, dummyEvent1, dummyEvent2
   
   type(dim3) :: dimGrid, dimBlock

   nlinks = 8 !pari!
   tau = 1.0_db
   cssq = 1.0_db/3.0_db
   visc_LB = cssq*(tau - 0.5_db)
   one_ov_nu = 1.0_db/visc_LB

   istat = cudaGetDeviceCount(ngpus)
!#ifdef _OPENACC
!        ngpus=acc_get_num_devices(acc_device_nvidia)
!#else
!        ngpus=0
!#endif

   !*******************************user parameters**************************

   nx = 128
   ny = 64
   TILE_DIMx = 8
   TILE_DIMy = 1
   TILE_DIM = 8

   if (mod(nx, TILE_DIMx) /= 0) then
      write (*, *) 'nx must be a multiple of TILE_DIM'
      stop
   end if
   if (mod(ny, TILE_DIMy) /= 0) then
      write (*, *) 'ny must be a multiple of TILE_DIMy'
      stop
   end if
   dimGrid = dim3(nx/TILE_DIMx, ny/TILE_DIMy, 1)
   dimBlock = dim3(TILE_DIMx, TILE_DIMy, 1)
   
   radius=20
   nsteps = 10
   stamp = 1
   lprint = .true.
   lvtk = .false.
   lpbc = .false.
   lasync = .false.
   fx = 0.0_db*10.0_db**(-5.0_db)
   fy = 0.0_db*10.0_db**(-6.0_db)
   allocate (p(0:nlinks))
   !allocate(f0(0:nx+1,0:ny+1),f1(0:nx+1,0:ny+1),f2(0:nx+1,0:ny+1),f3(0:nx+1,0:ny+1),f4(0:nx+1,0:ny+1))
   !allocate(f5(0:nx+1,0:ny+1),f6(0:nx+1,0:ny+1),f7(0:nx+1,0:ny+1),f8(0:nx+1,0:ny+1))
   !allocate(rho(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny),pxx(1:nx,1:ny),pyy(1:nx,1:ny),pxy(1:nx,1:ny))
   allocate (isfluid(1:nx, 1:ny)) !,omega_2d(1:nx,1:ny))

   !ex=(/0,1,0,-1,0,1,-1,-1,1/)
   !ey=(/0,0,1,0,-1,1,1,-1,-1/)

   p = (/4.0_db/9.0_db, 1.0_db/9.0_db, 1.0_db/9.0_db, 1.0_db/9.0_db, 1.0_db/9.0_db, 1.0_db/36.0_db, &
         1.0_db/36.0_db, 1.0_db/36.0_db, 1.0_db/36.0_db/)

   omega = 1.0_db/tau

   ! regularized: hermite
   qxx = 1.0_db - cssq
   qyy = 1.0_db - cssq
   qxy5_7 = 1.0_db
   qxy6_8 = -1.0_db
   pi2cssq0 = p(0)/(2.0_db*cssq**2)
   pi2cssq1 = p(1)/(2.0_db*cssq**2)
   pi2cssq2 = p(5)/(2.0_db*cssq**2)
   !**************************************color gradient********************
   beta = 0.95_db
   sigma = 0.02_db
   b0 = -4.0_db/27.0_db
   b1 = 2.0_db/27.0_db
   b2 = 5.0_db/108.0_db
   max_press_excess=0.0_db
   !*****************************************geometry************************
   isfluid = 1
   isfluid(1, :) = 0 !EAST
   isfluid(nx, :) = 0 !WEST
   isfluid(:, 1) = 0 !SOUTH
   isfluid(:, ny) = 0 !NORTH
   !*************************************initial conditions ************************
   myu = 0.0_db
   myv = 0.0_db
   myrhoA = 1.0_db
   myrhoB = 0.0_db     !rho!
   !do ll=0,nlinks
!    f0(1:nx,1:ny)=p(0)*rho(:,:)!0.0_db
!    f1(1:nx,1:ny)=p(1)*rho(:,:)
!    f2(1:nx,1:ny)=p(2)*rho(:,:)
!    f3(1:nx,1:ny)=p(3)*rho(:,:)
!    f4(1:nx,1:ny)=p(4)*rho(:,:)
!    f5(1:nx,1:ny)=p(5)*rho(:,:)
!    f6(1:nx,1:ny)=p(6)*rho(:,:)
!    f7(1:nx,1:ny)=p(7)*rho(:,:)
!    f8(1:nx,1:ny)=p(8)*rho(:,:)
   !enddo
   !*************************************check data ************************
   write (6, '(a)') '*******************LB data*****************'
   write (6, *) 'tau', tau
   write (6, *) 'omega', omega
   write (6, *) 'visc', visc_LB
   write (6, *) 'fx', fx
   write (6, *) 'cssq', cssq
   write (6, '(a)') '*******************INPUT data*****************'
   write (6, *) 'nx', nx
   write (6, *) 'ny', ny
   write (6, *) 'lpbc', lpbc
   write (6, *) 'lprint', lprint
   write (6, *) 'lvtk', lvtk
   write (6, *) 'lasync', lasync
   write (6, *) 'nsteps', nsteps
   write (6, *) 'stamp', stamp
   write (6, *) 'max fx', huge(fx)
   write (6, *) 'TILE_DIMx ', TILE_DIMx
   write (6, *) 'TILE_DIMy ', TILE_DIMy
   write (6, *) 'TILE_DIM ', TILE_DIM
   write (6, *) 'available gpus', ngpus
   write (6, '(a)') '*******************************************'
   istat = cudaGetDeviceProperties(prop, 0)

   call printDeviceProperties(prop, 6, 1)
   ! create events and streams
   istat = cudaStreamCreate(stream1)
   if(lasync)istat = cudaStreamCreate(stream2)
   istat = cudaforSetDefaultstream(stream1)
   istat = cudaDeviceSynchronize

   istat = cudaEventCreate(startEvent)
   istat = cudaEventCreate(stopEvent)
   istat = cudaEventCreate(dummyEvent)
   istat = cudaEventCreate(dummyEvent1)
   istat = cudaEventCreate(dummyEvent2)

!!!!!!!!!!!!!from host to device!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nx_d = nx
   ny_d = ny
   TILE_DIMx_d = TILE_DIMx
   TILE_DIMy_d = TILE_DIMy
   TILE_DIM_d = TILE_DIM
   myrhoA_d = myrhoA
   myrhoB_d = myrhoB
   myu_d = myu
   myv_d = myv
   fx_d = fx
   fy_d = fy
   p_d = p
   omega_d = omega
   qxx_d = qxx
   qyy_d = qyy
   qxy5_7_d = qxy5_7
   qxy6_8_d = qxy6_8
   cssq_d = cssq
   pi2cssq0_d = pi2cssq0
   pi2cssq1_d = pi2cssq1
   pi2cssq2_d = pi2cssq2
   b0_d=b0
   b1_d=b1
   b2_d=b2
   sigma_d=sigma
   beta_d=beta
   max_press_excess_d=max_press_excess
   allocate (isfluid_d(1:nx_d, 1:ny_d))
   istat = cudaDeviceSynchronize
   istat = cudaMemcpy(isfluid_d, isfluid, nx*ny)
   istat = cudaDeviceSynchronize
   if (istat /= 0) write (*, *) 'status after copy isfluid:', istat
   allocate(rhoA_d(1:nx_d,1:ny_d),rhoB_d(1:nx_d,1:ny_d),u_d(1:nx_d,1:ny_d),v_d(1:nx_d,1:ny_d),pxx_d(1:nx_d,1:ny_d),pyy_d(1:nx_d,1:ny_d),pxy_d(1:nx_d,1:ny_d))
   allocate(f0_d(0:nx_d+1,0:ny_d+1),f1_d(0:nx_d+1,0:ny_d+1),f2_d(0:nx_d+1,0:ny_d+1),f3_d(0:nx_d+1,0:ny_d+1),f4_d(0:nx_d+1,0:ny_d+1))
   allocate (f5_d(0:nx_d + 1, 0:ny_d + 1), f6_d(0:nx_d + 1, 0:ny_d + 1), f7_d(0:nx_d + 1, 0:ny_d + 1), f8_d(0:nx_d + 1, 0:ny_d + 1))
   allocate(g0_d(0:nx_d+1,0:ny_d+1),g1_d(0:nx_d+1,0:ny_d+1),g2_d(0:nx_d+1,0:ny_d+1),g3_d(0:nx_d+1,0:ny_d+1),g4_d(0:nx_d+1,0:ny_d+1))
   allocate(g5_d(0:nx_d+1,0:ny_d+1),g6_d(0:nx_d+1,0:ny_d+1),g7_d(0:nx_d+1,0:ny_d+1),g8_d(0:nx_d+1,0:ny_d+1))
   allocate(psi_d(1:nx_d,1:ny_d))
   rhoA_d=0.0_db
   rhoB_d=0.0_db
   psi_d=0.0_db
   u_d=0.0_db
   v_d=0.0_db
   istat = cudaDeviceSynchronize
   
   call setup_pops <<< dimGrid, dimBlock >>> (radius)
   istat = cudaDeviceSynchronize
   
   allocate(rhoprint(1:nx, 1:ny, 1:nz), velprint(3, 1:nx, 1:ny, 1:nz))
   allocate(rhoprint_d(1:nx_d, 1:ny_d, 1:nz_d), velprint_d(3, 1:nx_d, 1:ny_d, 1:nz_d))
   if (lprint) then
         call init_output(nx, ny, nz, 1, lvtk)
         call string_char(head1, nheadervtk(1), headervtk(1))
         call string_char(head2, nheadervtk(2), headervtk(2))
   end if

   istat = cudaDeviceSynchronize
   iframe = 0
   step = 0
   if (lprint) then
         call moments <<< dimGrid, dimBlock, 0, stream1 >>> ()
         call store_print <<< dimGrid, dimBlock, 0, stream1 >>> ()
         istat = cudaEventRecord(dummyEvent1, stream1)
         istat = cudaEventSynchronize(dummyEvent1)
         !write(6,*)'ciao 1',step,iframe
         if (lasync) then
               istat = cudaMemcpyAsync(rhoprint, rhoprint_d, nx*ny*nz, cudaMemcpyDeviceToHost, stream2)
               istat = cudaMemcpyAsync(velprint, velprint_d, 3*nx*ny*nz, cudaMemcpyDeviceToHost, stream2)
         else
               istat = cudaMemcpy(rhoprint, rhoprint_d, nx*ny*nz)
               istat = cudaMemcpy(velprint, velprint_d, 3*nx*ny*nz)
               istat = cudaEventRecord(dummyEvent, 0)
               istat = cudaEventSynchronize(dummyEvent)
               if (lvtk) then
                  call print_vtk_sync
               else
                  call print_raw_sync
               end if
         end if
   end if

   !*************************************time loop************************
   call cpu_time(ts1)
   istat = cudaEventRecord(startEvent, 0)
   
   do step = 1, nsteps
      !***********************************moment + neq pressor*********

      call moments <<< dimGrid, dimBlock, 0, stream1 >>> ()

      !***********************************PRINT************************
         if (mod(step, stamp) .eq. 0) write (6, '(a,i8)') 'step : ', step
         if (lprint) then
            if (mod(step, stamp) .eq. 0) then
                  iframe = iframe + 1
                  !write(6,*)'ciao 1',step,iframe
                  istat = cudaEventRecord(dummyEvent1, stream1)
                  istat = cudaEventSynchronize(dummyEvent1)
                  call store_print <<< dimGrid, dimBlock, 0, stream1 >>> ()
                  istat = cudaEventRecord(dummyEvent1, stream1)
                  istat = cudaEventSynchronize(dummyEvent1)
                  if (lasync) then
                     call close_print_async
                     istat = cudaMemcpyAsync(rhoprint, rhoprint_d, nx*ny*nz, cudaMemcpyDeviceToHost, stream2)
                     istat = cudaMemcpyAsync(velprint, velprint_d, 3*nx*ny*nz, cudaMemcpyDeviceToHost, stream2)
                  else
                     istat = cudaMemcpy(rhoprint, rhoprint_d, nx*ny*nz, cudaMemcpyDeviceToHost)
                     istat = cudaMemcpy(velprint, velprint_d, 3*nx*ny*nz, cudaMemcpyDeviceToHost)
                     istat = cudaEventRecord(dummyEvent, 0)
                     istat = cudaEventSynchronize(dummyEvent)
                  if (lvtk) then
                     call print_vtk_sync
                  else
                     call print_raw_sync
                  end if
               end if
            end if
        
            if (mod(step - stamp/4, stamp) .eq. 0 .and. lasync) then
               !write(6,*)'ciao 2',step,iframe
               istat = cudaEventRecord(dummyEvent2, stream2)
               istat = cudaEventSynchronize(dummyEvent2)
               if (lvtk) then
                  call print_vtk_async
               else
                  call print_raw_async
               end if
            end if
         end if

         !***********************************collision + no slip + forcing: fused implementation*********
         call streamcoll <<< dimGrid, dimBlock, 0, stream1  >>> ()
         !********************************************bcs no slip*****************************************!
         call bcs_no_slip <<< dimGrid, dimBlock, 0, stream1 >>> ()
         !******************************************call periodic bcs: always after fused************************
         !periodic along y
         if (lpbc) then
            call pbc_edge_x <<< (ny + TILE_DIM - 1)/TILE_DIM, TILE_DIM, 0, stream1 >>> ()
            !call pbc_edge_y<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM,0,stream1>>>()
         end if

         istat = cudaEventRecord(dummyEvent, stream1)
         istat = cudaEventSynchronize(dummyEvent)

   end do

   if (lasync) then
      !write(6,*)'ciao 2',step,iframe
      istat = cudaEventRecord(dummyEvent2, stream2)
      istat = cudaEventSynchronize(dummyEvent2)
      if (lvtk) then
         call print_vtk_sync
      else
         call print_raw_sync
      end if
   end if
   
   istat = cudaDeviceSynchronize
   call cpu_time(ts2)
   istat = cudaEventRecord(stopEvent, 0)
   istat = cudaEventSynchronize(stopEvent)
   istat = cudaEventElapsedTime(time, startEvent, stopEvent)
   write (6, *) 'Time elapsed as measured with cuda    : ', time/1000.0, ' s of your life time'

   !************************************************test points**********************************************!
!    write(6,*) 'u=',u(nx/2,ny/2) ,'v=',v(nx/2,ny/2),'rho',rho(nx/2,ny/2) !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
!    write(6,*) 'u=',u(2,ny/2) ,'v=',v(2,ny/2),'rho',rho(2,ny/2)
!    write(6,*) 'u=',u(1,ny/2) ,'v=',v(1,ny/2),'rho',rho(1,ny/2)

   write (6, *) 'time elapsed as measured from cpu_time: ', ts2 - ts1, ' s of your life time'
   write (6, *) 'glups: ', real(nx)*real(ny)*real(nsteps)*real(1.d-9, kind=db)/(ts2 - ts1)

   ! istat = cudaDeviceSynchronize
   ! call store_print <<< dimGrid, dimBlock >>> ()
   ! istat = cudaDeviceSynchronize
   ! istat = cudaMemcpy(rhoprint, rhoprint_d, nx*ny*nz)
   ! istat = cudaMemcpy(velprint, velprint_d, 3*nx*ny*nz)
   ! istat = cudaDeviceSynchronize

   ! open (101, file='v.out', status='replace')
   ! do j = 1, ny
   !    !do i=1,nx
   !    i = nx/2
   !    write (101, *) velprint(2, i, j, 1)
   !    !enddo
   ! end do
   ! close (101)

   ! call get_memory_gpu(mymemory, totmemory)
   ! call print_memory_registration_gpu(6, 'DEVICE memory occupied at the end', &
   !                                    'total DEVICE memory', mymemory, totmemory)

contains

   subroutine printDeviceProperties(prop, iu, num)

         use cudafor
         type(cudadeviceprop) :: prop
         integer, intent(in) :: iu, num

         write (iu, 907) "                                                                               "
         write (iu, 907) "*****************************GPU FEATURE MONITOR*******************************"
         write (iu, 907) "                                                                               "

         write (iu, 900) "Device Number: ", num
         write (iu, 901) "Device Name: ", trim(prop%name)
         write (iu, 903) "Total Global Memory: ", real(prop%totalGlobalMem)/1e9, " Gbytes"
         write (iu, 902) "sharedMemPerBlock: ", prop%sharedMemPerBlock, " bytes"
         write (iu, 900) "regsPerBlock: ", prop%regsPerBlock
         write (iu, 900) "warpSize: ", prop%warpSize
         write (iu, 900) "maxThreadsPerBlock: ", prop%maxThreadsPerBlock
         write (iu, 904) "maxThreadsDim: ", prop%maxThreadsDim
         write (iu, 904) "maxGridSize: ", prop%maxGridSize
         write (iu, 903) "ClockRate: ", real(prop%clockRate)/1e6, " GHz"
         write (iu, 902) "Total Const Memory: ", prop%totalConstMem, " bytes"
         write (iu, 905) "Compute Capability Revision: ", prop%major, prop%minor
         write (iu, 902) "TextureAlignment: ", prop%textureAlignment, " bytes"
         write (iu, 906) "deviceOverlap: ", prop%deviceOverlap
         write (iu, 900) "multiProcessorCount: ", prop%multiProcessorCount
         write (iu, 906) "integrated: ", prop%integrated
         write (iu, 906) "canMapHostMemory: ", prop%canMapHostMemory
         write (iu, 906) "ECCEnabled: ", prop%ECCEnabled
         write (iu, 906) "UnifiedAddressing: ", prop%unifiedAddressing
         write (iu, 900) "L2 Cache Size: ", prop%l2CacheSize
         write (iu, 900) "maxThreadsPerSMP: ", prop%maxThreadsPerMultiProcessor

         write (iu, 907) "                                                                               "
         write (iu, 907) "*******************************************************************************"
         write (iu, 907) "                                                                               "

900   format(a, i0)
901   format(a, a)
902   format(a, i0, a)
903   format(a, f16.8, a)
904   format(a, 2(i0, 1x, 'x', 1x), i0)
905   format(a, i0, '.', i0)
906   format(a, l0)
907   format(a)
      return
   end subroutine printDeviceProperties

   subroutine print_raw_sync

      implicit none

      sevt1 = trim(dir_out)//trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
              '_'//trim(write_fmtnumb(iframe))//'.raw'
      sevt2 = trim(dir_out)//trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
              '_'//trim(write_fmtnumb(iframe))//'.raw'
      open (unit=345, file=trim(sevt1), &
            status='replace', action='write', access='stream', form='unformatted')
      write (345) rhoprint
      close (345)
      open (unit=346, file=trim(sevt2), &
            status='replace', action='write', access='stream', form='unformatted')
      write (346) velprint
      close (346)

   end subroutine print_raw_sync

   subroutine print_vtk_sync
      implicit none

      sevt1 = trim(dir_out)//trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
              '_'//trim(write_fmtnumb(iframe))//'.vti'
      sevt2 = trim(dir_out)//trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
              '_'//trim(write_fmtnumb(iframe))//'.vti'
      open (unit=345, file=trim(sevt1), &
            status='replace', action='write', access='stream', form='unformatted')
      write (345) head1, ndatavtk(1), rhoprint(1:nx, 1:ny, 1:nz), footervtk(1)
      close (345)
      open (unit=346, file=trim(sevt2), &
            status='replace', action='write', access='stream', form='unformatted')
      write (346) head2, ndatavtk(2), velprint(1:3, 1:nx, 1:ny, 1:nz), footervtk(2)
      close (346)

   end subroutine print_vtk_sync

   subroutine print_raw_async

      implicit none

      sevt1 = trim(dir_out)//trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
              '_'//trim(write_fmtnumb(iframe))//'.raw'
      sevt2 = trim(dir_out)//trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
              '_'//trim(write_fmtnumb(iframe))//'.raw'
      open (unit=345, file=trim(sevt1), &
            status='replace', action='write', access='stream', form='unformatted', &
            asynchronous='yes')
      write (345, asynchronous='yes') rhoprint

      open (unit=346, file=trim(sevt2), &
            status='replace', action='write', access='stream', form='unformatted', &
            asynchronous='yes')
      write (346, asynchronous='yes') velprint

   end subroutine print_raw_async

   subroutine print_vtk_async
      implicit none

      sevt1 = trim(dir_out)//trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
              '_'//trim(write_fmtnumb(iframe))//'.vti'
      sevt2 = trim(dir_out)//trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
              '_'//trim(write_fmtnumb(iframe))//'.vti'

      open (unit=345, file=trim(sevt1), &
            status='replace', action='write', access='stream', form='unformatted', &
            asynchronous='yes')
      write (345, asynchronous='yes') head1, ndatavtk(1), rhoprint

      open (unit=780, file=trim(sevt2), &
            status='replace', action='write', access='stream', form='unformatted', &
            asynchronous='yes')
      write (780, asynchronous='yes') head2, ndatavtk(2), velprint

   end subroutine print_vtk_async

   subroutine close_print_async

      implicit none

      wait(345)
      if (lvtk) write (345) footervtk(1)
      close (345)

      wait(780)
      if (lvtk) write (780) footervtk(2)
      close (780)

   end subroutine close_print_async

end program

