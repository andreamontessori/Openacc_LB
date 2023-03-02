
 program main
 
  implicit none
  
   integer, parameter :: nx=8
   integer, parameter :: ny=8
   integer, parameter :: npops=8
                                                 !0 1 2  3  4 5  6  7  8
   integer, parameter, dimension(0:npops) :: ex=(/0,1,0,-1, 0,1,-1,-1, 1/)
   integer, parameter, dimension(0:npops) :: ey=(/0,0,1, 0,-1,1, 1,-1,-1/)
   integer, parameter, dimension(0:npops) ::opp=(/0,3,4, 1, 2,7, 8, 5, 6/)
   
   integer :: isfluid(0:nx+1,0:ny+1)
   

   
   integer :: i,j,k,ip,jp,kp,ii,jj,kk,io,jo,ko,l,lopp,idir
   character(len=1) :: sl,slopp
   character(len=8) :: sip,sio,sjp,sjo
   
   isfluid=3
   isfluid(2:nx-1,2:ny-1)=1
   do i=1,nx
      do j=1,ny
        !dead node with around a fluid node should be set equal to 0
        if(isfluid(i,j).eq.3)then
          do l=0,npops
            ii=i+ex(l)
            jj=j+ey(l)
            if(ii.gt.0 .and. ii.lt.nx+1 .and. jj.gt.0 .and. jj.lt.ny+1)then
              if(isfluid(ii,jj).eq.1)then
                isfluid(i,j)=0
              endif
            endif
           enddo
         endif
      enddo
    enddo
   
   idir=1
   select case(idir)
   case(1)
   
     io=nx
     ip=2
     sio='nx_d    '
     sip='2       '
     
     !io=1
     !ip=nx-1
     !sio='1       '
     !sip='nx_d-1  '
   
   
     j=ny-1
     sjo='j       '
     sjp=sjo
     i=io
   case(2)
     jo=ny
     jp=2
     sjo='ny_d    '
     sjp='2       '
     
     jo=1
     jp=ny-1
     sjo='1       '
     sjp='ny_d-1  '
   
   
     i=nx/2
     sio='i       '
     sip=sio
     j=jo
   end select
   
   write(6,*)'i j ',i,j
   do l=1,npops
     lopp=opp(l)
     write(sl,'(i1)')l
     write(slopp,'(i1)')lopp
     ii=i+ex(lopp)
     jj=j+ey(lopp)
     if(isfluid(ii,jj)==1)then
       write(6,'(13a)')'f',sl,'_d(',trim(sip),',',trim(sjp),')=f',sl,'_d(',trim(sio),',',trim(sjo),')'
     endif
   enddo
   
   
 end program
   
   
   
   
