 
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
      integer(kind=1), allocatable,  dimension(:,:,:), device   :: isfluid_d
      integer, constant :: nx_d,ny_d,nz_d,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
      
      integer, constant :: nxblock_d
      integer, constant :: nxyblock_d
      integer, constant :: nblocks_d
      
      real(kind=db), constant :: fx,fy,fz,omega,myrho,myu,myv,myw,oneminusomega
!      real(kind=db), allocatable, dimension(:,:,:,:), device  :: rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
!      real(kind=db), allocatable, dimension(:,:,:,:), device :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
!      real(kind=db), allocatable, dimension(:,:,:,:), device :: f10,f11,f12,f13,f14,f15,f16,f17,f18
      
      real(kind=db), allocatable, dimension(:,:,:,:,:), device :: fpops,hfields
      
      real(kind=db), allocatable, dimension(:,:,:), device :: rhoprint_d
      real(kind=db), allocatable, dimension(:,:,:,:), device :: velprint_d
      type (dim3) :: dimGrid,dimBlock,dimGridx,dimGridy,dimBlock2, &
       dimGridhalo,dimBlockhalo
   
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
            integer :: gi,gj,gk,idblock
          
            
            gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
            gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
            gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	        i=threadIdx%x
	        j=threadIdx%y
	        k=threadIdx%z
	        
	        idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
            
            !write(*,*)i,j,p(0)*myrho
            
            fpops(i,j,k,0,idblock)=p0*myrho
            fpops(i,j,k,1,idblock)=p1*myrho
            fpops(i,j,k,2,idblock)=p1*myrho
            fpops(i,j,k,3,idblock)=p1*myrho
            fpops(i,j,k,4,idblock)=p1*myrho
            fpops(i,j,k,5,idblock)=p1*myrho
            fpops(i,j,k,6,idblock)=p1*myrho
            fpops(i,j,k,7,idblock)=p2*myrho
            fpops(i,j,k,8,idblock)=p2*myrho
            fpops(i,j,k,9,idblock)=p2*myrho
            fpops(i,j,k,10,idblock)=p2*myrho
            fpops(i,j,k,11,idblock)=p2*myrho
            fpops(i,j,k,12,idblock)=p2*myrho
            fpops(i,j,k,13,idblock)=p2*myrho
            fpops(i,j,k,14,idblock)=p2*myrho
            fpops(i,j,k,15,idblock)=p2*myrho
            fpops(i,j,k,16,idblock)=p2*myrho
            fpops(i,j,k,17,idblock)=p2*myrho
            fpops(i,j,k,18,idblock)=p2*myrho
     

      end subroutine setup_pops
      
      attributes(global) subroutine setup_pops_halo()
      
            integer :: i,j,k
            integer :: gi,gj,gk,idblock
          
            
            gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x -1
            gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y -1
            gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z -1
	
            if(gi>nx_d .or. gj>ny_d .or. gk>nz_d)return
	
            i=threadIdx%x
            j=threadIdx%y
            k=threadIdx%z
	        
            idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
            
            !write(*,*)i,j,p(0)*myrho
            
            fpops(i,j,k,0,idblock)=p0*myrho
            fpops(i,j,k,1,idblock)=p1*myrho
            fpops(i,j,k,2,idblock)=p1*myrho
            fpops(i,j,k,3,idblock)=p1*myrho
            fpops(i,j,k,4,idblock)=p1*myrho
            fpops(i,j,k,5,idblock)=p1*myrho
            fpops(i,j,k,6,idblock)=p1*myrho
            fpops(i,j,k,7,idblock)=p2*myrho
            fpops(i,j,k,8,idblock)=p2*myrho
            fpops(i,j,k,9,idblock)=p2*myrho
            fpops(i,j,k,10,idblock)=p2*myrho
            fpops(i,j,k,11,idblock)=p2*myrho
            fpops(i,j,k,12,idblock)=p2*myrho
            fpops(i,j,k,13,idblock)=p2*myrho
            fpops(i,j,k,14,idblock)=p2*myrho
            fpops(i,j,k,15,idblock)=p2*myrho
            fpops(i,j,k,16,idblock)=p2*myrho
            fpops(i,j,k,17,idblock)=p2*myrho
            fpops(i,j,k,18,idblock)=p2*myrho
     

      end subroutine setup_pops_halo
  
      attributes(global) subroutine moments()
          
            real(kind=db) :: uu,udotc,temp
            real(kind=db) ::fneq1,fneq2,fneq3,fneq4,fneq5,fneq6
            real(kind=db) ::fneq7,fneq8,fneq9,fneq10,fneq11,fneq12
            real(kind=db) ::fneq13,fneq14,fneq15,fneq16,fneq17,fneq18
            
            integer :: i,j,k
            integer :: gi,gj,gk,idblock
          
            
            gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	        gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	        gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	        i=threadIdx%x
            j=threadIdx%y
            k=threadIdx%z
	
            idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
            
            if(isfluid_d(gi,gj,gk).eq.1)then
                        hfields(i,j,k,1,idblock) = fpops(i,j,k,0,idblock)+fpops(i,j,k,1,idblock)+fpops(i,j,k,2,idblock)+fpops(i,j,k,3,idblock)+fpops(i,j,k,4,idblock)+fpops(i,j,k,5,idblock) &
                            +fpops(i,j,k,6,idblock)+fpops(i,j,k,7,idblock)+fpops(i,j,k,8,idblock)+fpops(i,j,k,9,idblock)+fpops(i,j,k,10,idblock)+fpops(i,j,k,11,idblock) &
                            +fpops(i,j,k,12,idblock)+fpops(i,j,k,13,idblock)+fpops(i,j,k,14,idblock)+fpops(i,j,k,15,idblock)+fpops(i,j,k,16,idblock)+fpops(i,j,k,17,idblock) &
                            +fpops(i,j,k,18,idblock)

                        hfields(i,j,k,2,idblock) = (fpops(i,j,k,1,idblock)+fpops(i,j,k,7,idblock)+fpops(i,j,k,9,idblock)+fpops(i,j,k,15,idblock)+fpops(i,j,k,18,idblock)) &
                             -(fpops(i,j,k,2,idblock)+fpops(i,j,k,8,idblock)+fpops(i,j,k,10,idblock)+fpops(i,j,k,16,idblock)+fpops(i,j,k,17,idblock)) 
                        
                        hfields(i,j,k,3,idblock) = (fpops(i,j,k,3,idblock)+fpops(i,j,k,7,idblock)+fpops(i,j,k,10,idblock)+fpops(i,j,k,11,idblock)+fpops(i,j,k,13,idblock)) &
                            -(fpops(i,j,k,4,idblock)+fpops(i,j,k,8,idblock)+fpops(i,j,k,9,idblock)+fpops(i,j,k,12,idblock)+fpops(i,j,k,14,idblock))

                        hfields(i,j,k,4,idblock) = (fpops(i,j,k,5,idblock)+fpops(i,j,k,11,idblock)+fpops(i,j,k,14,idblock)+fpops(i,j,k,15,idblock)+fpops(i,j,k,17,idblock)) &
                            -(fpops(i,j,k,6,idblock)+fpops(i,j,k,12,idblock)+fpops(i,j,k,13,idblock)+fpops(i,j,k,16,idblock)+fpops(i,j,k,18,idblock))
                        
                        uu=halfonecssq*(hfields(i,j,k,2,idblock)*hfields(i,j,k,2,idblock) + hfields(i,j,k,3,idblock)*hfields(i,j,k,3,idblock) + hfields(i,j,k,4,idblock)*hfields(i,j,k,4,idblock))
                        !1-2
                        udotc=hfields(i,j,k,2,idblock)*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=fpops(i,j,k,1,idblock)-p1*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq2=fpops(i,j,k,2,idblock)-p1*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        !3-4
                        udotc=hfields(i,j,k,3,idblock)*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq3=fpops(i,j,k,3,idblock)-p1*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq4=fpops(i,j,k,4,idblock)-p1*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        !5-6
                        udotc=hfields(i,j,k,4,idblock)*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq5=fpops(i,j,k,5,idblock)-p1*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq6=fpops(i,j,k,6,idblock)-p1*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        !7-8
                        udotc=(hfields(i,j,k,2,idblock)+hfields(i,j,k,3,idblock))*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq7=fpops(i,j,k,7,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq8=fpops(i,j,k,8,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        !10-9
                        udotc=(-hfields(i,j,k,2,idblock)+hfields(i,j,k,3,idblock))*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq10=fpops(i,j,k,10,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq9=fpops(i,j,k,9,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        !11-12
                        udotc=(hfields(i,j,k,3,idblock)+hfields(i,j,k,4,idblock))*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq11=fpops(i,j,k,11,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq12=fpops(i,j,k,12,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        !13-14
                        udotc=(hfields(i,j,k,3,idblock)-hfields(i,j,k,4,idblock))*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq13=fpops(i,j,k,13,idblock) - p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq14=fpops(i,j,k,14,idblock) - p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        !15-16
                        udotc=(hfields(i,j,k,2,idblock)+hfields(i,j,k,4,idblock))*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq15=fpops(i,j,k,15,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq16=fpops(i,j,k,16,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        !17-18
                        udotc=(-hfields(i,j,k,2,idblock)+hfields(i,j,k,4,idblock))*onecssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq17=fpops(i,j,k,17,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp + udotc))
                        fneq18=fpops(i,j,k,18,idblock)-p2*(hfields(i,j,k,1,idblock)+(temp - udotc))
                        hfields(i,j,k,5,idblock)=fneq1+fneq2+fneq7+fneq8+fneq9+fneq10+fneq15+fneq16+fneq17+fneq18
                        hfields(i,j,k,6,idblock)=fneq3+fneq4+fneq7+fneq8+fneq9+fneq10+fneq11+fneq12+fneq13+fneq14
                        hfields(i,j,k,7,idblock)=fneq5+fneq6+fneq11+fneq12+fneq13+fneq14+fneq15+fneq16+fneq17+fneq18
                        hfields(i,j,k,8,idblock)= fneq7+fneq8-fneq9-fneq10
                        hfields(i,j,k,9,idblock)=fneq15+fneq16-fneq17-fneq18
                        hfields(i,j,k,10,idblock)=fneq11+fneq12-fneq13-fneq14
                    endif

           
      end subroutine moments
  
      attributes(global) subroutine streamcoll()
          
            integer :: i,j,k
            real(kind=db) ::uu,udotc,temp,feq
            integer :: gi,gj,gk,myblock
            integer :: gii,gjj,gkk,xblock,yblock,zblock,iidblock,ii,jj,kk
            real(kind=db) ::fneq1,fneq2,fneq3,fneq4,fneq5,fneq6
            real(kind=db) ::fneq7,fneq8,fneq9,fneq10,fneq11,fneq12
            real(kind=db) ::fneq13,fneq14,fneq15,fneq16,fneq17,fneq18
          
            
    
            gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	        gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	        gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	        i=threadIdx%x
	        j=threadIdx%y
	        k=threadIdx%z
	
	        myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
            
              
            if(isfluid_d(gi,gj,gk).ne.1)return
              
              
            uu=halfonecssq*(hfields(i,j,k,2,myblock)*hfields(i,j,k,2,myblock) + hfields(i,j,k,3,myblock)*hfields(i,j,k,3,myblock) + hfields(i,j,k,4,myblock)*hfields(i,j,k,4,myblock))
            !0
            feq=p0*(hfields(i,j,k,1,myblock)-uu)
            fpops(i,j,k,0,myblock)=feq + oneminusomega*pi2cssq0*(-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock)))
            
            !1
            gii=gi+1
            gjj=gj
            gkk=gk
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=hfields(i,j,k,2,myblock)*onecssq
            temp = -uu + half*udotc*udotc
            feq=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,1,iidblock)=feq + oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock)-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) + fx*p1dcssq
            
            !2
            gii=gi-1
            gjj=gj
            gkk=gk
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,2,iidblock)=feq + oneminusomega*pi2cssq1*(qxx*hfields(i,j,k,5,myblock)-cssq*(hfields(i,j,k,6,myblock)+hfields(i,j,k,7,myblock))) - fx*p1dcssq
            
            !3
            gii=gi
            gjj=gj+1
            gkk=gk
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=hfields(i,j,k,3,myblock)*onecssq
            temp = -uu + half*udotc*udotc
            feq=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,3,iidblock)=feq+oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) + fy*p1dcssq
            
            !4
            gii=gi
            gjj=gj-1
            gkk=gk
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,4,iidblock)=feq+oneminusomega*pi2cssq1*(qyy*hfields(i,j,k,6,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,7,myblock))) - fy*p1dcssq
            
            !7
            gii=gi+1
            gjj=gj+1
            gkk=gk
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,7,iidblock)=feq + oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) + (fx+fy)*p2dcssq 
            
            !8
            gii=gi-1
            gjj=gj-1
            gkk=gk
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,8,iidblock)=feq + oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_7_8*hfields(i,j,k,8,myblock)) - (fx+fy)*p2dcssq
            
            !10
            gii=gi-1
            gjj=gj+1
            gkk=gk
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,3,myblock))*onecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,10,iidblock)=feq+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) +(fy-fx)*p2dcssq
            
            !9
            gii=gi+1
            gjj=gj-1
            gkk=gk
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,9,iidblock)=feq+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qyy*hfields(i,j,k,6,myblock)-cssq*hfields(i,j,k,7,myblock)+two*qxy_9_10*hfields(i,j,k,8,myblock)) + (fx-fy)*p2dcssq

            !5
            gii=gi
            gjj=gj
            gkk=gk+1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=hfields(i,j,k,4,myblock)*onecssq
            temp = -uu + half*udotc*udotc
            feq=p1*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,5,iidblock)=feq+oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) + fz*p1dcssq
            
            !6
            gii=gi
            gjj=gj
            gkk=gk-1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p1*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,6,iidblock)=feq+oneminusomega*pi2cssq1*(qzz*hfields(i,j,k,7,myblock)-cssq*(hfields(i,j,k,5,myblock)+hfields(i,j,k,6,myblock))) - fz*p1dcssq

            !15
            gii=gi+1
            gjj=gj
            gkk=gk+1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=(hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,15,iidblock)=feq+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) + (fx+fz)*p2dcssq 
            
            !16
            gii=gi-1
            gjj=gj
            gkk=gk-1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,16,iidblock)=feq+ oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_15_16*hfields(i,j,k,9,myblock)) - (fx+fz)*p2dcssq

            !17
            gii=gi-1
            gjj=gj
            gkk=gk+1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=(-hfields(i,j,k,2,myblock)+hfields(i,j,k,4,myblock))*onecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,17,iidblock)=feq+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) +(fz-fx)*p2dcssq
            
            !18
            gii=gi+1
            gjj=gj
            gkk=gk-1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,18,iidblock)=feq+oneminusomega*pi2cssq2*(qxx*hfields(i,j,k,5,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,6,myblock)+two*qxz_17_18*hfields(i,j,k,9,myblock)) + (fx-fz)*p2dcssq

            !11
            gii=gi
            gjj=gj+1
            gkk=gk+1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=(hfields(i,j,k,3,myblock)+hfields(i,j,k,4,myblock))*onecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,11,iidblock)=feq+oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock))+(fy+fz)*p2dcssq
            
            !12
            gii=gi
            gjj=gj-1
            gkk=gk-1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,12,iidblock)=feq+oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_11_12*hfields(i,j,k,10,myblock)) - (fy+fz)*p2dcssq

            !13
            gii=gi
            gjj=gj+1
            gkk=gk-1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            udotc=(hfields(i,j,k,3,myblock)-hfields(i,j,k,4,myblock))*onecssq
            temp = -uu + half*udotc*udotc
            feq=p2*(hfields(i,j,k,1,myblock)+(temp + udotc))
            fpops(ii,jj,kk,13,iidblock)=feq+oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) + (fy-fz)*p2dcssq
            
            !14
            gii=gi
            gjj=gj-1
            gkk=gk+1
            xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
            yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
            zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
            iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
            ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
            jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
            kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
            
            feq=p2*(hfields(i,j,k,1,myblock)+(temp - udotc))
            fpops(ii,jj,kk,14,iidblock)=feq+oneminusomega*pi2cssq2*(qyy*hfields(i,j,k,6,myblock)+qzz*hfields(i,j,k,7,myblock)-cssq*hfields(i,j,k,5,myblock)+two*qyz_13_14*hfields(i,j,k,10,myblock)) + (fz-fy)*p2dcssq

      end subroutine streamcoll
  
      attributes(global) subroutine bcs_no_slip()
          
          integer :: i,j,k
          integer :: gi,gj,gk,myblock
          integer :: gii,gjj,gkk,xblock,yblock,zblock,iidblock,ii,jj,kk
          
	      
	      gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	      gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	      gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	      i=threadIdx%x
	      j=threadIdx%y
	      k=threadIdx%z
	
	      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
        
          
          !write(*,*)i,j,p(0)*myrho_d
          if(isfluid_d(gi,gj,gk).ne.0)return
          
          gii=gi+1
          gjj=gj
          gkk=gk-1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,18,iidblock)=fpops(i,j,k,17,myblock) !gpc 
          
          gii=gi-1
          gjj=gj
          gkk=gk+1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,17,iidblock)=fpops(i,j,k,18,myblock) !hpc
          
          gii=gi-1
          gjj=gj
          gkk=gk-1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,16,iidblock)=fpops(i,j,k,15,myblock) !gpc 
          
          gii=gi+1
          gjj=gj
          gkk=gk+1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,15,iidblock)=fpops(i,j,k,16,myblock) !hpc
          
          gii=gi
          gjj=gj-1
          gkk=gk+1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,14,iidblock)=fpops(i,j,k,13,myblock)!gpc 
          
          gii=gi
          gjj=gj+1
          gkk=gk-1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,13,iidblock)=fpops(i,j,k,14,myblock)!hpc
          
          gii=gi
          gjj=gj-1
          gkk=gk-1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,12,iidblock)=fpops(i,j,k,11,myblock)!gpc 
          
          gii=gi
          gjj=gj+1
          gkk=gk+1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,11,iidblock)=fpops(i,j,k,12,myblock)!hpc
          
          gii=gi-1
          gjj=gj+1
          gkk=gk
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,10,iidblock)=fpops(i,j,k,9,myblock)!gpc 
          
          gii=gi+1
          gjj=gj-1
          gkk=gk
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,9,iidblock)=fpops(i,j,k,10,myblock)!hpc
          
          gii=gi-1
          gjj=gj-1
          gkk=gk
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,8,iidblock)=fpops(i,j,k,7,myblock)!gpc 
          
          gii=gi+1
          gjj=gj+1
          gkk=gk
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,7,iidblock)=fpops(i,j,k,8,myblock)!hpc
          
          gii=gi
          gjj=gj
          gkk=gk-1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,6,iidblock)=fpops(i,j,k,5,myblock)!gpc 
          
          gii=gi
          gjj=gj
          gkk=gk+1
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,5,iidblock)=fpops(i,j,k,6,myblock)!hpc 
          
          gii=gi
          gjj=gj-1
          gkk=gk
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,4,iidblock)=fpops(i,j,k,3,myblock)!gpc 
          
          gii=gi
          gjj=gj+1
          gkk=gk
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,3,iidblock)=fpops(i,j,k,4,myblock)!hpc 
          
          gii=gi-1
          gjj=gj
          gkk=gk
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d
          fpops(ii,jj,kk,2,iidblock)=fpops(i,j,k,1,myblock)!gpc
          
          gii=gi+1
          gjj=gj
          gkk=gk
          xblock=(gii+TILE_DIMx_d-1)/TILE_DIMx_d
          yblock=(gjj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gkk+TILE_DIMz_d-1)/TILE_DIMz_d
          iidblock=xblock+yblock*nxblock_d+zblock*nxyblock_d+1
          ii=gii-xblock*TILE_DIMx_d+TILE_DIMx_d
          jj=gjj-yblock*TILE_DIMy_d+TILE_DIMy_d
          kk=gkk-zblock*TILE_DIMz_d+TILE_DIMz_d 
          fpops(ii,jj,kk,1,iidblock)=fpops(i,j,k,2,myblock)!hpc 

        

      end subroutine bcs_no_slip
  
      attributes(global) subroutine pbc_edge_x()
          
          integer :: j,k
          integer :: gj,gk
          integer :: yblock,zblock
          integer :: gio,gie,xblocko,xblocke
          integer :: idblocko1,idblocke1,io1,ie1
          integer :: idblocko2,idblocke2,io2,ie2
          
          gj = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
	      gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
	      
          
          if (gj>ny_d .or. gk>nz_d)return
          
          yblock=(gj+TILE_DIMy_d-1)/TILE_DIMy_d
          zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
          j=gj-yblock*TILE_DIMy_d+TILE_DIMy_d
          k=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
          
          gio=nx_d
	      gie=2
	      xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
          idblocko1=xblocko+yblock*nxblock_d+zblock*nxyblock_d+1
            
          xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
	      idblocke1=xblocke+yblock*nxblock_d+zblock*nxyblock_d+1
	
	      io1=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
	      ie1=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
	      
	      gio=1
	      gie=nx_d-1
	      xblocko=(gio+TILE_DIMx_d-1)/TILE_DIMx_d
          idblocko2=xblocko+yblock*nxblock_d+zblock*nxyblock_d+1
            
          xblocke=(gie+TILE_DIMx_d-1)/TILE_DIMx_d
	      idblocke2=xblocke+yblock*nxblock_d+zblock*nxyblock_d+1
	
	      io2=gio-xblocko*TILE_DIMx_d+TILE_DIMx_d
	      ie2=gie-xblocke*TILE_DIMx_d+TILE_DIMx_d
          
          if(gj>2 .and. gj<ny_d-1 .and. gk>2 .and. gk<nz_d-1)then
            
            fpops(ie1,j,k,1,idblocke1)=fpops(io1,j,k,1,idblocko1)
            fpops(ie1,j,k,7,idblocke1)=fpops(io1,j,k,7,idblocko1)
            fpops(ie1,j,k,9,idblocke1)=fpops(io1,j,k,9,idblocko1)
            fpops(ie1,j,k,15,idblocke1)=fpops(io1,j,k,15,idblocko1)
            fpops(ie1,j,k,18,idblocke1)=fpops(io1,j,k,18,idblocko1)
            
            fpops(ie2,j,k,2,idblocke2)=fpops(io2,j,k,2,idblocko2)
            fpops(ie2,j,k,8,idblocke2)=fpops(io2,j,k,8,idblocko2)
            fpops(ie2,j,k,10,idblocke2)=fpops(io2,j,k,10,idblocko2)
            fpops(ie2,j,k,16,idblocke2)=fpops(io2,j,k,16,idblocko2)
            fpops(ie2,j,k,17,idblocke2)=fpops(io2,j,k,17,idblocko2)
          else
          
            if(gj==2)then
              if(gk==2)then
                
                fpops(ie1,j,k,1,idblocke1)=fpops(io1,j,k,1,idblocko1)
                fpops(ie1,j,k,9,idblocke1)=fpops(io1,j,k,9,idblocko1)
                fpops(ie1,j,k,18,idblocke1)=fpops(io1,j,k,18,idblocko1)
                
                fpops(ie2,j,k,2,idblocke2)=fpops(io2,j,k,2,idblocko2)
                fpops(ie2,j,k,8,idblocke2)=fpops(io2,j,k,8,idblocko2)
                fpops(ie2,j,k,16,idblocke2)=fpops(io2,j,k,16,idblocko2)
              elseif(gk==nz_d-1)then
                fpops(ie1,j,k,1,idblocke1)=fpops(io1,j,k,1,idblocko1)
                fpops(ie1,j,k,9,idblocke1)=fpops(io1,j,k,9,idblocko1)
                fpops(ie1,j,k,15,idblocke1)=fpops(io1,j,k,15,idblocko1)
                
                fpops(ie2,j,k,2,idblocke2)=fpops(io2,j,k,2,idblocko2)
                fpops(ie2,j,k,8,idblocke2)=fpops(io2,j,k,8,idblocko2)
                fpops(ie2,j,k,17,idblocke2)=fpops(io2,j,k,17,idblocko2)
              else
                fpops(ie1,j,k,1,idblocke1)=fpops(io1,j,k,1,idblocko1)
                fpops(ie1,j,k,9,idblocke1)=fpops(io1,j,k,9,idblocko1)
                fpops(ie1,j,k,15,idblocke1)=fpops(io1,j,k,15,idblocko1)
                fpops(ie1,j,k,18,idblocke1)=fpops(io1,j,k,18,idblocko1)
                
                fpops(ie2,j,k,2,idblocke2)=fpops(io2,j,k,2,idblocko2)
                fpops(ie2,j,k,8,idblocke2)=fpops(io2,j,k,8,idblocko2)
                fpops(ie2,j,k,16,idblocke2)=fpops(io2,j,k,16,idblocko2)
                fpops(ie2,j,k,17,idblocke2)=fpops(io2,j,k,17,idblocko2)
              endif
            elseif(gj==ny_d-1)then
              if(gk==2)then
                fpops(ie1,j,k,1,idblocke1)=fpops(io1,j,k,1,idblocko1)
                fpops(ie1,j,k,7,idblocke1)=fpops(io1,j,k,7,idblocko1)
                fpops(ie1,j,k,18,idblocke1)=fpops(io1,j,k,18,idblocko1)
                
                fpops(ie2,j,k,2,idblocke2)=fpops(io2,j,k,2,idblocko2)
                fpops(ie2,j,k,10,idblocke2)=fpops(io2,j,k,10,idblocko2)
                fpops(ie2,j,k,16,idblocke2)=fpops(io2,j,k,16,idblocko2)
              elseif(gk==nz_d-1)then
                fpops(ie1,j,k,1,idblocke1)=fpops(io1,j,k,1,idblocko1)
                fpops(ie1,j,k,7,idblocke1)=fpops(io1,j,k,7,idblocko1)
                fpops(ie1,j,k,15,idblocke1)=fpops(io1,j,k,15,idblocko1)
                
                fpops(ie2,j,k,2,idblocke2)=fpops(io2,j,k,2,idblocko2)
                fpops(ie2,j,k,10,idblocke2)=fpops(io2,j,k,10,idblocko2)
                fpops(ie2,j,k,17,idblocke2)=fpops(io2,j,k,17,idblocko2)
              else
                fpops(ie1,j,k,1,idblocke1)=fpops(io1,j,k,1,idblocko1)
                fpops(ie1,j,k,7,idblocke1)=fpops(io1,j,k,7,idblocko1)
                fpops(ie1,j,k,15,idblocke1)=fpops(io1,j,k,15,idblocko1)
                fpops(ie1,j,k,18,idblocke1)=fpops(io1,j,k,18,idblocko1)
                
                fpops(ie2,j,k,2,idblocke2)=fpops(io2,j,k,2,idblocko2)
                fpops(ie2,j,k,10,idblocke2)=fpops(io2,j,k,10,idblocko2)
                fpops(ie2,j,k,16,idblocke2)=fpops(io2,j,k,16,idblocko2)
                fpops(ie2,j,k,17,idblocke2)=fpops(io2,j,k,17,idblocko2)
              endif
            endif
          endif 
              
      end subroutine pbc_edge_x
  
      attributes(global) subroutine pbc_edge_y()
          
          integer :: i,k
          integer :: gi,gk
          integer :: xblock,zblock
          integer :: gjo,gje,yblocko,yblocke
          integer :: idblocko1,idblocke1,jo1,je1
          integer :: idblocko2,idblocke2,jo2,je2
        
          
          gi = (blockIdx%x-1) * TILE_DIM_d + threadIdx%x
          gk = (blockIdx%y-1) * TILE_DIM_d + threadIdx%y
          
          if (gi>nx_d .or. gk>nz_d)return
          
          xblock=(gi+TILE_DIMx_d-1)/TILE_DIMx_d
          zblock=(gk+TILE_DIMz_d-1)/TILE_DIMz_d
          i=gi-xblock*TILE_DIMx_d+TILE_DIMx_d
          k=gk-zblock*TILE_DIMz_d+TILE_DIMz_d
	
	      gjo=ny_d
	      gje=2
	
	      yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
          idblocko1=xblock+yblocko*nxblock_d+zblock*nxyblock_d+1
    
          yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
	      idblocke1=xblock+yblocke*nxblock_d+zblock*nxyblock_d+1
	      
	      jo1=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	      je1=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
	      
	      gjo=1
	      gje=ny_d-1
	
	      yblocko=(gjo+TILE_DIMy_d-1)/TILE_DIMy_d
          idblocko2=xblock+yblocko*nxblock_d+zblock*nxyblock_d+1
    
          yblocke=(gje+TILE_DIMy_d-1)/TILE_DIMy_d
	      idblocke2=xblock+yblocke*nxblock_d+zblock*nxyblock_d+1
	      
	      jo2=gjo-yblocko*TILE_DIMy_d+TILE_DIMy_d
	      je2=gje-yblocke*TILE_DIMy_d+TILE_DIMy_d
          
          if(gi>2 .and. gi<nx_d-1 .and. gk>2 .and. gk<nz_d-1)then
          
            fpops(i,je1,k,3,idblocke1)=fpops(i,jo1,k,3,idblocko1)
            fpops(i,je1,k,7,idblocke1)=fpops(i,jo1,k,7,idblocko1)
            fpops(i,je1,k,10,idblocke1)=fpops(i,jo1,k,10,idblocko1)
            fpops(i,je1,k,11,idblocke1)=fpops(i,jo1,k,11,idblocko1)
            fpops(i,je1,k,13,idblocke1)=fpops(i,jo1,k,13,idblocko1)
            
            fpops(i,je2,k,4,idblocke2)=fpops(i,jo2,k,4,idblocko2)
            fpops(i,je2,k,8,idblocke2)=fpops(i,jo2,k,8,idblocko2)
            fpops(i,je2,k,9,idblocke2)=fpops(i,jo2,k,9,idblocko2)
            fpops(i,je2,k,12,idblocke2)=fpops(i,jo2,k,12,idblocko2)
            fpops(i,je2,k,14,idblocke2)=fpops(i,jo2,k,14,idblocko2)
          else
          
            if(gi==2)then
              if(gk==2)then
                fpops(i,je1,k,3,idblocke1)=fpops(i,jo1,k,3,idblocko1)
                fpops(i,je1,k,10,idblocke1)=fpops(i,jo1,k,10,idblocko1)
                fpops(i,je1,k,13,idblocke1)=fpops(i,jo1,k,13,idblocko1)
                
                fpops(i,je2,k,4,idblocke2)=fpops(i,jo2,k,4,idblocko2)
                fpops(i,je2,k,8,idblocke2)=fpops(i,jo2,k,8,idblocko2)
                fpops(i,je2,k,12,idblocke2)=fpops(i,jo2,k,12,idblocko2)
              elseif(gk==nz_d-1)then
                fpops(i,je1,k,3,idblocke1)=fpops(i,jo1,k,3,idblocko1)
                fpops(i,je1,k,10,idblocke1)=fpops(i,jo1,k,10,idblocko1)
                fpops(i,je1,k,11,idblocke1)=fpops(i,jo1,k,11,idblocko1)

                fpops(i,je2,k,4,idblocke2)=fpops(i,jo2,k,4,idblocko2)
                fpops(i,je2,k,8,idblocke2)=fpops(i,jo2,k,8,idblocko2)
                fpops(i,je2,k,14,idblocke2)=fpops(i,jo2,k,14,idblocko2)
              else
                fpops(i,je1,k,3,idblocke1)=fpops(i,jo1,k,3,idblocko1)
                fpops(i,je1,k,10,idblocke1)=fpops(i,jo1,k,10,idblocko1)
                fpops(i,je1,k,11,idblocke1)=fpops(i,jo1,k,11,idblocko1)
                fpops(i,je1,k,13,idblocke1)=fpops(i,jo1,k,13,idblocko1)

                fpops(i,je2,k,4,idblocke2)=fpops(i,jo2,k,4,idblocko2)
                fpops(i,je2,k,8,idblocke2)=fpops(i,jo2,k,8,idblocko2)
                fpops(i,je2,k,12,idblocke2)=fpops(i,jo2,k,12,idblocko2)
                fpops(i,je2,k,14,idblocke2)=fpops(i,jo2,k,14,idblocko2)
              endif
            elseif(gi==nx_d-1)then
              if(gk==2)then
                fpops(i,je1,k,3,idblocke1)=fpops(i,jo1,k,3,idblocko1)
                fpops(i,je1,k,7,idblocke1)=fpops(i,jo1,k,7,idblocko1)
                fpops(i,je1,k,13,idblocke1)=fpops(i,jo1,k,13,idblocko1)

                fpops(i,je2,k,4,idblocke2)=fpops(i,jo2,k,4,idblocko2)
                fpops(i,je2,k,9,idblocke2)=fpops(i,jo2,k,9,idblocko2)
                fpops(i,je2,k,12,idblocke2)=fpops(i,jo2,k,12,idblocko2)
              elseif(gk==nz_d-1)then 
                fpops(i,je1,k,3,idblocke1)=fpops(i,jo1,k,3,idblocko1)
                fpops(i,je1,k,7,idblocke1)=fpops(i,jo1,k,7,idblocko1)
                fpops(i,je1,k,11,idblocke1)=fpops(i,jo1,k,11,idblocko1)

                fpops(i,je2,k,4,idblocke2)=fpops(i,jo2,k,4,idblocko2)
                fpops(i,je2,k,9,idblocke2)=fpops(i,jo2,k,9,idblocko2)
                fpops(i,je2,k,14,idblocke2)=fpops(i,jo2,k,14,idblocko2)
              else
                fpops(i,je1,k,3,idblocke1)=fpops(i,jo1,k,3,idblocko1)
                fpops(i,je1,k,7,idblocke1)=fpops(i,jo1,k,7,idblocko1)
                fpops(i,je1,k,11,idblocke1)=fpops(i,jo1,k,11,idblocko1)
                fpops(i,je1,k,13,idblocke1)=fpops(i,jo1,k,13,idblocko1)

                fpops(i,je2,k,4,idblocke2)=fpops(i,jo2,k,4,idblocko2)
                fpops(i,je2,k,9,idblocke2)=fpops(i,jo2,k,9,idblocko2)
                fpops(i,je2,k,12,idblocke2)=fpops(i,jo2,k,12,idblocko2)
                fpops(i,je2,k,14,idblocke2)=fpops(i,jo2,k,14,idblocko2)
              endif
            endif
          endif
        

      end subroutine pbc_edge_y
  

  
  
      attributes(global) subroutine store_print()
          
          integer :: i,j,k
          integer :: gi,gj,gk,idblock
          
          
	      gi = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
	      gj = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
	      gk = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
	
	      i=threadIdx%x
	      j=threadIdx%y
	      k=threadIdx%z
	
	      idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
          
          
          !write(*,*)i,j,p(0)*myrho_d
          if(isfluid_d(gi,gj,gk).eq.1)then
            rhoprint_d(gi,gj,gk)=hfields(i,j,k,1,idblock)
            velprint_d(1,gi,gj,gk)=hfields(i,j,k,2,idblock)
            velprint_d(2,gi,gj,gk)=hfields(i,j,k,3,idblock)
            velprint_d(3,gi,gj,gk)=hfields(i,j,k,4,idblock)
            
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
    real(kind=db) :: visc_LB,h_omega,h_oneminusomega
    real(kind=db) :: tau,one_ov_nu,h_fx,h_fy,h_fz
    real(kind=db) :: time

    real(kind=db) :: h_myrho,h_myu,h_myv,h_myw
    
    integer(kind=1), allocatable,dimension(:,:,:)   :: isfluid
    !real(kind=db), allocatable, dimension(:,:,:) :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    integer :: TILE_DIMx,TILE_DIMy,TILE_DIMz,TILE_DIM,istat,iframe
    integer :: nxblock,nyblock,nzblock,nblocks
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
        
        h_myu=zero
        h_myv=zero
        h_myw=zero
        h_myrho=one  !tot dens
        nsteps=1000
        stamp=100
        lprint=.false.
        lvtk=.false.
        lpbc=.true.
        lasync=.false.
        h_fx=one*10.0**(-5)
        h_fy=zero*10.0**(-5)
        h_fz=zero*10.0**(-5)
        
        TILE_DIMx=8
        TILE_DIMy=8
        TILE_DIMz=4
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
        
        dimGridhalo  = dim3((nx+2+TILE_DIMx-1)/TILE_DIMx,(ny+2+TILE_DIMy-1)/TILE_DIMy,(nz+2+TILE_DIMz-1)/TILE_DIMz)
        dimBlockhalo = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
        
        dimGridx  = dim3((ny+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
        dimGridy  = dim3((nx+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
        dimBlock2 = dim3(TILE_DIM, TILE_DIM, 1)
        
        !plus 2 for the halo forward and backward
        nxblock=nx/TILE_DIMx +2
        nyblock=ny/TILE_DIMy +2
        nzblock=nz/TILE_DIMz +2
        
        nblocks=nxblock*nyblock*nzblock
        
        write(6,*)'nx,ny,nz',nx,ny,nz
        write(6,*)'TILE_DIMx,TILE_DIMy,TILE_DIMz',TILE_DIMx,TILE_DIMy,TILE_DIMz
        write(6,*)'nxblock,nyblock,nzblock',nxblock,nyblock,nzblock
        write(6,*)'nblocks',nblocks

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
        
        h_omega=one/tau
        h_oneminusomega=one-h_omega
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
        myrho=h_myrho
        myu=h_myu
        myv=h_myv
        myw=h_myw
        fx=h_fx
        fy=h_fy
        fz=h_fz
        omega=h_omega
        oneminusomega=h_oneminusomega
        nxblock_d=nxblock
        nxyblock_d=nxblock*nyblock
        nblocks_d=nblocks
        
        allocate(isfluid_d(1:nx_d,1:ny_d,1:nz_d))
        istat = cudaDeviceSynchronize
        istat = cudaMemcpy(isfluid_d,isfluid,nx*ny*nz )
        istat = cudaDeviceSynchronize
        
        
        
!        allocate(f0(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f1(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f2(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f3(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(f4(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f5(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f6(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f7(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(f8(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f9(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f10(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f11(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(f12(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f13(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f14(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f15(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(f16(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f17(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),f18(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
!        allocate(rho(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),u(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),v(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),w(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks), &
!         pxx(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pyy(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pzz(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks), &
!         pxy(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pxz(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks),pyz(TILE_DIMx,TILE_DIMy,TILE_DIMz,nblocks))
        
        allocate(fpops(TILE_DIMx,TILE_DIMy,TILE_DIMz,0:18,nblocks))
        allocate(hfields(TILE_DIMx,TILE_DIMy,TILE_DIMz,10,nblocks))
         
        
         
        istat = cudaDeviceSynchronize
    !*************************************initial conditions ************************    
        
        !call setup_pops<<<dimGrid,dimBlock>>>()
        call setup_pops_halo<<<dimGridhalo,dimBlockhalo>>>()
        call abortOnLastErrorAndSync('after setup_pops', step)
        istat = cudaDeviceSynchronize
        
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz),velprint(3,1:nx,1:ny,1:nz))
          allocate(rhoprint_d(1:nx_d,1:ny_d,1:nz_d),velprint_d(3,1:nx_d,1:ny_d,1:nz_d))
        endif
    !*************************************check data ************************ 
        write(6,*) '*******************LB data*****************'
        write(6,*) 'tau',tau
        write(6,*) 'omega',h_omega
        write(6,*) 'visc',visc_LB
        write(6,*) 'fx',h_fx
        write(6,*) 'fy',h_fy
        write(6,*) 'fz',h_fz
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
        write(6,*) 'max fx',huge(h_fx)
        write(6,*) 'max fx',huge(h_fy)
        write(6,*) 'max fx',huge(h_fz)
        write(6,*) 'TILE_DIMx ',TILE_DIMx
        write(6,*) 'TILE_DIMy ',TILE_DIMy
        write(6,*) 'TILE_DIMz ',TILE_DIMz
        write(6,*) 'TILE_DIM ',TILE_DIM
        write(6,*) 'available gpus',ngpus
        write(6,*) '*******************************************'
	    write(6,*) 'TILE_DIMx',TILE_DIMx
        write(6,*) 'TILE_DIMy',TILE_DIMy
        write(6,*) 'TILE_DIMz',TILE_DIMz
        write(6,*) 'TILE_DIM ',TILE_DIM
        write(6,*)'nxblock,nyblock,nzblock',nxblock,nyblock,nzblock
        write(6,*)'nblocks',nblocks
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
    istat = cudaEventRecord(startEvent,0)
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
