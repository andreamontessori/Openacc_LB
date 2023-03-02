
 program main
 
  implicit none
  
   integer, parameter :: nx=8
   integer, parameter :: ny=8
   integer, parameter :: nz=8
   integer, parameter :: npops=18
                                                 !0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18
   integer, parameter, dimension(0:npops) :: ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
   integer, parameter, dimension(0:npops) :: ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
   integer, parameter, dimension(0:npops) :: ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
   integer, parameter, dimension(0:npops) ::opp=(/0, 2,  1, 4,  3,  6,  5,  8,   7, 10,   9, 12,  11, 14,  13, 16,  15,  18,  17/)
   
   integer :: isfluid(0:nx+1,0:ny+1,0:nz+1)
   

   
   integer :: i,j,k,ip,jp,kp,ii,jj,kk,io,jo,ko,l,lopp,idir
   character(len=2) :: sl,slopp
   character(len=8) :: sip,sio,sjp,sjo,skp,sko
   
   isfluid=3
   isfluid(2:nx-1,2:ny-1,2:nz-1)=1
   do k=1,nz
    do j=1,ny
      do i=1,nx
        !dead node with around a fluid node should be set equal to 0
        if(isfluid(i,j,k).eq.3)then
          do l=0,npops
            ii=i+ex(l)
            jj=j+ey(l)
            kk=k+ez(l)
            if(ii.gt.0 .and. ii.lt.nx+1 .and. jj.gt.0 .and. jj.lt.ny+1 .and. kk.gt.0 .and. kk.lt.nz+1)then
              if(isfluid(ii,jj,kk).eq.1)then
                isfluid(i,j,k)=0
              endif
            endif
           enddo
         endif
      enddo
    enddo
   enddo
   
   idir=1
   select case(idir)
   case(1)
     io=nx
     ip=2
     sio='nx_d    '
     sip='2       '
     
     io=1
     ip=nx-1
     sio='1       '
     sip='nx_d-1  '
   
     j=ny-1
     k=nz/2
     sjo='j       '
     sjp=sjo
     sko='k       '
     skp=sko
     i=io
   
   
   case(2)
   
     jo=ny
     jp=2
     sjo='ny_d    '
     sjp='2       '
   
     i=nx-1
     k=nz/2
     sio='i       '
     sip=sio
     sko='k       '
     skp=sko
     j=jo
   end select
   
   write(6,*)'i j k ',i,j,k
   do l=1,npops
     lopp=opp(l)
     write(sl,'(i2)')l
     write(slopp,'(i2)')lopp
     ii=i+ex(lopp)
     jj=j+ey(lopp)
     kk=k+ez(lopp)
     if(isfluid(ii,jj,kk)==1)then
       write(6,'(18a)')'f',trim(sl),'_d(',trim(sip),',',trim(sjp),',',trim(skp),')=f', &
        trim(sl),'_d(',trim(sio),',',trim(sjo),',',trim(sko),')'
     endif
   enddo
   
   
 end program
   
   
   
   
