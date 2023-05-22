 module prints
  
    implicit none
  
    integer, parameter :: db=4 !kind(1.0)
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
 
 endmodule

program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    use prints
    
    implicit none
    
    
    integer :: i,j,ll,l,dumm,ii,jj
    integer :: nx,ny,step,stamp,nlinks,nsteps,ngpus
    integer :: istat,iframe
    integer, parameter :: nz=1
    
    logical :: lprint,lvtk,lasync,lpbc
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2 
    real(kind=db) :: visc_LB,uu,udotc,omega,feq,uu0
    real(kind=db) :: fneq
    real(kind=db), parameter :: p0=(4.0_db/9.0_db)
    real(kind=db), parameter :: p1=(1.0_db/9.0_db)
    real(kind=db), parameter :: p2=(1.0_db/36.0_db)
    real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq1,pi2cssq2,pi2cssq0
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,dummy
    integer, parameter :: npops=8
                                                  !0 1 2  3  4 5  6  7  8
    integer, parameter, dimension(0:npops) :: ex=(/0,1,0,-1, 0,1,-1,-1, 1/)
    integer, parameter, dimension(0:npops) :: ey=(/0,0,1, 0,-1,1, 1,-1,-1/)
    integer, parameter, dimension(0:npops) ::opp=(/0,3,4, 1, 2,7, 8, 5, 6/)
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:,:) :: rho,u,v,pxx,pyy,pxy
    real(kind=db), allocatable, dimension(:,:) :: rhoh,uh,vh,pxxh,pyyh,pxyh
    real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8
    real(kind=db) :: mymemory,totmemory
    !$if _OPENACC
    integer :: devNum
    integer(acc_device_kind) :: devType
    devType = acc_get_device_type()
    devNum=acc_get_device_num(devType)
    !$endif
       
    nlinks=8 !pari!
    tau=1.0_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
#ifdef _OPENACC
    ngpus=acc_get_num_devices(acc_device_nvidia)
#else
    ngpus=0
#endif

    !*******************************user parameters**************************
    lpbc=.true.
    lprint=.true.
    lvtk=.true.
    lasync=.false.
    nx=256
    ny=64
    nsteps=100
    stamp=10
    fx=1.0_db*10.0**(-6)
    fy=0.0_db*10.0**(-8)
    
    allocate(f0(0:nx+1,0:ny+1),f1(0:nx+1,0:ny+1),f2(0:nx+1,0:ny+1),f3(0:nx+1,0:ny+1),f4(0:nx+1,0:ny+1))
    allocate(f5(0:nx+1,0:ny+1),f6(0:nx+1,0:ny+1),f7(0:nx+1,0:ny+1),f8(0:nx+1,0:ny+1))
    allocate(rho(0:nx+1,0:ny+1),u(0:nx+1,0:ny+1),v(0:nx+1,0:ny+1),pxx(0:nx+1,0:ny+1),pyy(0:nx+1,0:ny+1),pxy(0:nx+1,0:ny+1))
    allocate(rhoh(0:nx+1,0:ny+1),uh(0:nx+1,0:ny+1),vh(0:nx+1,0:ny+1),pxxh(0:nx+1,0:ny+1),pyyh(0:nx+1,0:ny+1),pxyh(0:nx+1,0:ny+1))
    allocate(isfluid(0:nx+1,0:ny+1)) !,omega_2d(1:nx,1:ny)) 
    if(lprint)then
      allocate(rhoprint(1:nx,1:ny,1:nz))
      allocate(velprint(1:3,1:nx,1:ny,1:nz))
      rhoprint(1:nx,1:ny,1:nz)=0.0
      velprint(1:3,1:nx,1:ny,1:nz)=0.0
    endif
    !ex=(/0,1,0,-1,0,1,-1,-1,1/)
    !ey=(/0,0,1,0,-1,1,1,-1,-1/)
    
    omega=1.0_db/tau

    ! regularized: hermite 
    qxx=1.0_db-cssq
    qyy=1.0_db-cssq
    qxy5_7=1.0_db
    qxy6_8=-1.0_db
    pi2cssq0=p0/(2.0_db*cssq**2)
    pi2cssq1=p1/(2.0_db*cssq**2)
    pi2cssq2=p2/(2.0_db*cssq**2)
    
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:)=0 !EAST
    isfluid(nx,:)=0 !WEST
    isfluid(:,1)=0 !SOUTH 
    isfluid(:,ny)=0 !NORTH
    if(lpbc)then
      isfluid=1
      isfluid(:,1)=0 !SOUTH 
      isfluid(:,ny)=0 !NORTH
    endif
	do j=1,ny
      do i=1,nx
	    if(isfluid(i,j).eq.1)then
		  do ll=1,npops
		    ii=i+ex(ll)
			jj=j+ey(ll)
		    if(ii.gt.0 .and. ii.lt.nx+1 .and. jj.gt.0 .and. jj.lt.ny+1)then
		      if(isfluid(ii,jj).eq.0)then
			    isfluid(i,j)=-1
			  endif
			endif
		  enddo
		endif
	  enddo
	enddo
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    rho=1.0_db     !rho!
    pxx=0.0_db
    pxy=0.0_db
    pyy=0.0_db
    
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
    write(6,*) 'lpbc',lpbc
    write(6,*) 'lprint',lprint
    write(6,*) 'lvtk',lvtk
    write(6,*) 'lasync',lasync
    write(6,*) 'nsteps',nsteps
    write(6,*) 'stamp',stamp
    write(6,*) 'max fx',huge(fx)
    write(6,*) 'available gpus',ngpus
    write(6,*) '*******************************************'
    step = 0

    !$acc data copy(rho,u,v,pxx,pxy,pyy,f0,f1,f2,f3,f4,f5,f6,f7,f8,isfluid,rhoprint,velprint) async(1)
    !$if _OPENACC        
    call printDeviceProperties(ngpus,devNum,devType,6)
    !$endif
    iframe=0
    write(6,'(a,i8,a,i8,3f16.4)')'start step : ',0,' frame ',iframe
    
    if(lprint)then  
      call init_output(nx,ny,nz,1,lvtk)
      call string_char(head1,nheadervtk(1),headervtk(1))
      call string_char(head2,nheadervtk(2),headervtk(2))
    endif
    
    !$acc wait(1)
    if(lprint)then
      !$acc kernels present(rhoprint,velprint,rho,u,v) async(1)
      !$acc loop independent collapse(2)  private(i,j)
      do j=1,ny
        do i=1,nx
          rhoprint(i,j,1)=real(rho(i,j),kind=4)
          velprint(1,i,j,1)=real(u(i,j),kind=4)
          velprint(2,i,j,1)=real(v(i,j),kind=4)
        enddo
      enddo
      !$acc end kernels 
      !$acc wait(1)
      if(lasync)then
        !$acc update host(rhoprint,velprint) async(2)
        continue
      else
        !$acc update host(rhoprint,velprint) async(2)
        !$acc wait(2)
        if(lvtk)then
          call print_vtk_sync(iframe)
        else
          call print_raw_sync(iframe)
        endif
      endif
    endif
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !******************************************call other bcs************************
        if(lpbc)then      
		  !periodic along x 
		  !$acc kernels async(1)
	      !$acc loop independent 
	       do j=2,ny-1
		     rho(1,j)=rho(nx-1,j)
		     u(1,j)=u(nx-1,j)
		     v(1,j)=v(nx-1,j)
		     pxx(1,j)=pxx(nx-1,j)
		     pyy(1,j)=pyy(nx-1,j)
		     pxy(1,j)=pxy(nx-1,j)
		     
		     rho(nx,j)=rho(2,j)
		     u(nx,j)=u(2,j)
		     v(nx,j)=v(2,j)
	         pxx(nx,j)=pxx(2,j)
	         pyy(nx,j)=pyy(2,j)
             pxy(nx,j)=pxy(2,j)
		  enddo
		!$acc end kernels
		endif
		
		!***********************************PRINT************************
        if(mod(step,stamp).eq.0) write(6,'(a,i8)')'step : ',step
        if(lprint)then
            if(mod(step,stamp).eq.0)then
                iframe=iframe+1
                !$acc wait(1)
                !$acc kernels present(rhoprint,velprint,rho,u,v) async(1)
                !$acc loop independent collapse(2)  private(i,j)
                do j=1,ny
                    do i=1,nx
                      rhoprint(i,j,1)=real(rho(i,j),kind=4)
                      velprint(1,i,j,1)=real(u(i,j),kind=4)
                      velprint(2,i,j,1)=real(v(i,j),kind=4)
                    enddo
                  enddo
                !$acc end kernels 
                !$acc wait(1)
                if(lasync)then
                    call close_print_async
                    !$acc update host(rhoprint,velprint) async(2)
                else
                    !$acc update host(rhoprint,velprint) async(2)
                    !$acc wait(2)
                    if(lvtk)then
                        call print_vtk_sync(iframe)
                    else
                        call print_raw_sync(iframe)
                    endif
                endif
            endif
            if(mod(step-stamp/4,stamp).eq.0 .and. lasync)then
                !write(6,*)'ciao 2',step,iframe
                !$acc wait(2)  
                if(lvtk)then
                    call print_vtk_async(iframe)
                else
                    call print_raw_async(iframe)
                endif
            endif
        endif
		
        
        !***********collision + no slip + forcing: fused implementation*********

        !$acc kernels async(1)
		!$acc loop collapse(2) private(uu,temp,udotc,feq,fneq,uu0) 
		do j=1,ny
		  do i=1,nx  
			  uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
			  !oneminusuu= -uu !1.0_db - uu
			  !0
			  feq=p0*(rho(i,j)-uu)
			  f0(i,j)=feq + (1.0_db-omega)*pi2cssq0*(-cssq*pxx(i,j)-cssq*pyy(i,j))
			  !1   -1  0
			  uu=0.5_db*(u(i-1,j)*u(i-1,j) + v(i-1,j)*v(i-1,j))/cssq
			  udotc=u(i-1,j)/cssq
			  temp = -uu + 0.5_db*udotc*udotc
			  feq=p1*(rho(i-1,j)+(temp + udotc))
			  f1(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i-1,j)-cssq*pyy(i-1,j)) + fx*p1/cssq 
			  !3   +1  0
			  uu=0.5_db*(u(i+1,j)*u(i+1,j) + v(i+1,j)*v(i+1,j))/cssq
			  udotc=u(i+1,j)/cssq
			  temp = -uu + 0.5_db*udotc*udotc
		  	  feq=p1*(rho(i+1,j)+(temp - udotc))
		  	  f3(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i+1,j)-cssq*pyy(i+1,j))  - fx*p1/cssq 
			  !2    0  -1
			  uu=0.5_db*(u(i,j-1)*u(i,j-1) + v(i,j-1)*v(i,j-1))/cssq
			  udotc=v(i,j-1)/cssq
			  temp = -uu + 0.5_db*udotc*udotc
			  feq=p1*(rho(i,j-1)+(temp + udotc))
			  f2(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j-1)-cssq*pxx(i,j-1))  + fy*p1/cssq 
			  !4    0   +1
			  uu=0.5_db*(u(i,j+1)*u(i,j+1) + v(i,j+1)*v(i,j+1))/cssq
			  udotc=v(i,j+1)/cssq
			  temp = -uu + 0.5_db*udotc*udotc
			  feq=p1*(rho(i,j+1)+(temp - udotc))
			  f4(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j+1)-cssq*pxx(i,j+1))  - fy*p1/cssq    
			  !5   -1   -1
			  uu=0.5_db*(u(i-1,j-1)*u(i-1,j-1) + v(i-1,j-1)*v(i-1,j-1))/cssq
			  udotc=(u(i-1,j-1)+v(i-1,j-1))/cssq
			  temp = -uu + 0.5_db*udotc*udotc
			  feq=p2*(rho(i-1,j-1)+(temp + udotc))
			  f5(i,j)= feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j-1)+qyy*pyy(i-1,j-1)+2.0_db*qxy5_7*pxy(i-1,j-1)) + &
			   fx*p2/cssq + fy*p2/cssq
			  !7   +1   +1
			  uu=0.5_db*(u(i+1,j+1)*u(i+1,j+1) + v(i+1,j+1)*v(i+1,j+1))/cssq
			  udotc=(u(i+1,j+1)+v(i+1,j+1))/cssq
			  temp = -uu + 0.5_db*udotc*udotc
			  feq=p2*(rho(i+1,j+1)+(temp - udotc))
			  f7(i,j)=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j+1)+qyy*pyy(i+1,j+1)+2.0_db*qxy5_7*pxy(i+1,j+1)) - &
			   fx*p2/cssq - fy*p2/cssq 
			  !6   +1   -1
			  uu=0.5_db*(u(i+1,j-1)*u(i+1,j-1) + v(i+1,j-1)*v(i+1,j-1))/cssq
			  udotc=(-u(i+1,j-1)+v(i+1,j-1))/cssq
			  temp = -uu + 0.5_db*udotc*udotc
			  feq=p2*(rho(i+1,j-1)+(temp + udotc))
			  f6(i,j)= feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j-1)+qyy*pyy(i+1,j-1)+2.0_db*qxy6_8*pxy(i+1,j-1)) - &
			   fx*p2/cssq + fy*p2/cssq 
			  !8   -1   +1
			  uu=0.5_db*(u(i-1,j+1)*u(i-1,j+1) + v(i-1,j+1)*v(i-1,j+1))/cssq
			  udotc=(-u(i-1,j+1)+v(i-1,j+1))/cssq
			  temp = -uu + 0.5_db*udotc*udotc
			  feq=p2*(rho(i-1,j+1)+(temp - udotc))
			  f8(i,j)=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j+1)+qyy*pyy(i-1,j+1)+2.0_db*qxy6_8*pxy(i-1,j+1)) + &
			   fx*p2/cssq - fy*p2/cssq
	     enddo
	   enddo
	  !$acc end kernels
      !$acc wait(1)
      
      !$acc kernels async(1)
	  !$acc loop collapse(2) private(uu,temp,udotc,feq,fneq,uu0) 
	   do j=1,ny
		  do i=1,nx  
			if(isfluid(i,j).eq.-1)then
			  uu0=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
			  !0
			  feq=p0*(rho(i,j)-uu0)
			  f0(i,j)=feq + (1.0_db-omega)*pi2cssq0*(-cssq*pxx(i,j)-cssq*pyy(i,j))
			  !1   -1  0
			  if(isfluid(i-1,j).eq.0)then
			    udotc=u(i,j)/cssq
			    temp = -uu0 + 0.5_db*udotc*udotc
			    feq=p1*(rho(i,j)+(temp - udotc))
			    f1(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j)-cssq*pyy(i,j)) - fx*p1/cssq 
			  else
			    uu=0.5_db*(u(i-1,j)*u(i-1,j) + v(i-1,j)*v(i-1,j))/cssq
			    udotc=u(i-1,j)/cssq
			    temp = -uu + 0.5_db*udotc*udotc
			    feq=p1*(rho(i-1,j)+(temp + udotc))
			    f1(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i-1,j)-cssq*pyy(i-1,j)) + fx*p1/cssq 
			  endif
			  !3   +1  0
			  if(isfluid(i+1,j).eq.0)then
			    udotc=u(i,j)/cssq
			    temp = -uu0 + 0.5_db*udotc*udotc
		  	    feq=p1*(rho(i,j)+(temp + udotc))
		  	    f3(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j)-cssq*pyy(i,j))  + fx*p1/cssq 
			  else
			    uu=0.5_db*(u(i+1,j)*u(i+1,j) + v(i+1,j)*v(i+1,j))/cssq
			    udotc=u(i+1,j)/cssq
			    temp = -uu + 0.5_db*udotc*udotc
		  	    feq=p1*(rho(i+1,j)+(temp - udotc))
		  	    f3(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i+1,j)-cssq*pyy(i+1,j))  - fx*p1/cssq 
		  	  endif
			  !2    0  -1
			  if(isfluid(i,j-1).eq.0)then
			    udotc=v(i,j)/cssq
			    temp = -uu0 + 0.5_db*udotc*udotc
			    feq=p1*(rho(i,j)+(temp - udotc))
			    f2(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j))  - fy*p1/cssq
			  else
			    uu=0.5_db*(u(i,j-1)*u(i,j-1) + v(i,j-1)*v(i,j-1))/cssq
			    udotc=v(i,j-1)/cssq
			    temp = -uu + 0.5_db*udotc*udotc
			    feq=p1*(rho(i,j-1)+(temp + udotc))
			    f2(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j-1)-cssq*pxx(i,j-1))  + fy*p1/cssq 
			  endif
			  !4    0   +1
			  if(isfluid(i,j+1).eq.0)then
			    udotc=v(i,j)/cssq
			    temp = -uu0 + 0.5_db*udotc*udotc
			    feq=p1*(rho(i,j)+(temp + udotc))
			    f4(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j)-cssq*pxx(i,j))  + fy*p1/cssq 
			  else
			    uu=0.5_db*(u(i,j+1)*u(i,j+1) + v(i,j+1)*v(i,j+1))/cssq
			    udotc=v(i,j+1)/cssq
			    temp = -uu + 0.5_db*udotc*udotc
			    feq=p1*(rho(i,j+1)+(temp - udotc))
			    f4(i,j)= feq + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j+1)-cssq*pxx(i,j+1))  - fy*p1/cssq 
			  endif
			  !5   -1   -1
			  if(isfluid(i-1,j-1).eq.0)then
			    udotc=(u(i,j)+v(i,j))/cssq
			    temp = -uu0 + 0.5_db*udotc*udotc
			    feq=p2*(rho(i,j)+(temp - udotc))
			    f5(i,j)= feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) - &
			     fx*p2/cssq - fy*p2/cssq
			  else
			    uu=0.5_db*(u(i-1,j-1)*u(i-1,j-1) + v(i-1,j-1)*v(i-1,j-1))/cssq
			    udotc=(u(i-1,j-1)+v(i-1,j-1))/cssq
			    temp = -uu + 0.5_db*udotc*udotc
			    feq=p2*(rho(i-1,j-1)+(temp + udotc))
			    f5(i,j)= feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j-1)+qyy*pyy(i-1,j-1)+2.0_db*qxy5_7*pxy(i-1,j-1)) + &
			     fx*p2/cssq + fy*p2/cssq
			  endif
			  !7   +1   +1
			  if(isfluid(i+1,j+1).eq.0)then
			    udotc=(u(i,j)+v(i,j))/cssq
			    temp = -uu0 + 0.5_db*udotc*udotc
			    feq=p2*(rho(i,j)+(temp + udotc))
			    f7(i,j)=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy5_7*pxy(i,j)) + &
			     fx*p2/cssq + fy*p2/cssq 
			  else
			    uu=0.5_db*(u(i+1,j+1)*u(i+1,j+1) + v(i+1,j+1)*v(i+1,j+1))/cssq
			    udotc=(u(i+1,j+1)+v(i+1,j+1))/cssq
			    temp = -uu + 0.5_db*udotc*udotc
			    feq=p2*(rho(i+1,j+1)+(temp - udotc))
			    f7(i,j)=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j+1)+qyy*pyy(i+1,j+1)+2.0_db*qxy5_7*pxy(i+1,j+1)) - &
			     fx*p2/cssq - fy*p2/cssq 
			  endif
			  !6   +1   -1
			  if(isfluid(i+1,j-1).eq.0)then
			    udotc=(-u(i,j)+v(i,j))/cssq
			    temp = -uu0 + 0.5_db*udotc*udotc
			    feq=p2*(rho(i,j)+(temp - udotc))
			    f6(i,j)= feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) + &
			     fx*p2/cssq - fy*p2/cssq 
			  else
			    uu=0.5_db*(u(i+1,j-1)*u(i+1,j-1) + v(i+1,j-1)*v(i+1,j-1))/cssq
			    udotc=(-u(i+1,j-1)+v(i+1,j-1))/cssq
			    temp = -uu + 0.5_db*udotc*udotc
			    feq=p2*(rho(i+1,j-1)+(temp + udotc))
			    f6(i,j)= feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j-1)+qyy*pyy(i+1,j-1)+2.0_db*qxy6_8*pxy(i+1,j-1)) - &
			     fx*p2/cssq + fy*p2/cssq 
			  endif
			  !8   -1   +1
			  if(isfluid(i-1,j+1).eq.0)then
			    udotc=(-u(i,j)+v(i,j))/cssq
			    temp = -uu + 0.5_db*udotc*udotc
			    feq=p2*(rho(i,j)+(temp + udotc))
			    f8(i,j)=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j)+qyy*pyy(i,j)+2.0_db*qxy6_8*pxy(i,j)) - &
			     fx*p2/cssq + fy*p2/cssq
			  else
			    uu=0.5_db*(u(i-1,j+1)*u(i-1,j+1) + v(i-1,j+1)*v(i-1,j+1))/cssq
			    udotc=(-u(i-1,j+1)+v(i-1,j+1))/cssq
			    temp = -uu + 0.5_db*udotc*udotc
			    feq=p2*(rho(i-1,j+1)+(temp - udotc))
			    f8(i,j)=feq + (1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j+1)+qyy*pyy(i-1,j+1)+2.0_db*qxy6_8*pxy(i-1,j+1)) + &
			     fx*p2/cssq - fy*p2/cssq
			  endif
	        endif
		  enddo
		enddo
		!$acc end kernels
        !$acc wait(1)

        
        !***********************************moment + neq pressor*********
        !$acc kernels async(1)
        !$acc loop collapse(2) private(uu,temp,udotc,fneq) 
        do j=1,ny
            do i=1,nx
                if(abs(isfluid(i,j)).eq.1)then
                    pxx(i,j)=0.0_db
                    pyy(i,j)=0.0_db
                    pxy(i,j)=0.0_db
                    rho(i,j) = f0(i,j)+f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)
                    u(i,j) = (f1(i,j) +f5(i,j) +f8(i,j)-f3(i,j) -f6(i,j) -f7(i,j)) !/rho(i,j)
                    v(i,j) = (f5(i,j) +f2(i,j) +f6(i,j)-f7(i,j) -f4(i,j) -f8(i,j))
                    ! non equilibrium pressor components
                    uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                    !1-3
                    udotc=u(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq=f1(i,j)-p1*(rho(i,j)+(temp + udotc))
                    pxx(i,j)=pxx(i,j)+fneq
                    fneq=f3(i,j)-p1*(rho(i,j)+(temp - udotc))
                    pxx(i,j)=pxx(i,j)+fneq
                    !2-4
                    udotc=v(i,j)/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq=f2(i,j)-p1*(rho(i,j)+(temp + udotc))
                    pyy(i,j)=pyy(i,j)+fneq
                    fneq=f4(i,j)-p1*(rho(i,j)+(temp - udotc))
                    pyy(i,j)=pyy(i,j)+fneq
                    !5-7
                    udotc=(u(i,j)+v(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq=f5(i,j)-p2*(rho(i,j)+(temp + udotc))
                    pxx(i,j)=pxx(i,j)+fneq
                    pyy(i,j)=pyy(i,j)+fneq
                    pxy(i,j)=pxy(i,j)+fneq
                    fneq=f7(i,j)-p2*(rho(i,j)+(temp - udotc))
                    pxx(i,j)=pxx(i,j)+fneq
                    pyy(i,j)=pyy(i,j)+fneq
                    pxy(i,j)=pxy(i,j)+fneq
                    !6-8
                    udotc=(-u(i,j)+v(i,j))/cssq
                    temp = -uu + 0.5_db*udotc*udotc
                    fneq=f6(i,j)-p2*(rho(i,j)+(temp + udotc))
                    pxx(i,j)=pxx(i,j)+fneq
                    pyy(i,j)=pyy(i,j)+fneq
                    pxy(i,j)=pxy(i,j)+fneq
                    fneq=f8(i,j)-p2*(rho(i,j)+(temp - udotc))
                    pxx(i,j)=pxx(i,j)+fneq
                    pyy(i,j)=pyy(i,j)+fneq
                    pxy(i,j)=pxy(i,j)+fneq

                endif
            enddo 
        enddo
        !$acc end kernels
        
        
    enddo 
    
    !$acc wait
    call cpu_time(ts2)
    !$acc update host(rho,u,v)
    !$acc end data


    !************************************************test points**********************************************!
    write(6,*) 'u=',u(nx/2,ny/2) ,'v=',v(nx/2,ny/2),'rho',rho(nx/2,ny/2) !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'u=',u(2,ny/2) ,'v=',v(2,ny/2),'rho',rho(2,ny/2)
    write(6,*) 'u=',u(1,ny/2) ,'v=',v(1,ny/2),'rho',rho(1,ny/2)
    
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nsteps)*real(1.d-9,kind=db)/(ts2-ts1)

    open(101, file = 'v.out', status = 'replace')
    do j=1,ny
        do i=1,nx
            write(101,*) v(i,j) 
        enddo
    enddo
    close(101) 
    
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)
    
    contains
  !$if _OPENACC  
  subroutine printDeviceProperties(ngpus,dev_Num,dev_Type,iu)
  
  
      use openacc
      
      integer :: ngpus,dev_Num
      integer(acc_device_kind) :: dev_Type
    
      integer,intent(in) :: iu 
      integer :: tot_mem,shared_mem
      character(len=255) :: myname,myvendor,mydriver
      
      call acc_get_property_string(dev_num,dev_type,acc_property_name,myname)
      tot_mem = acc_get_property(dev_num,dev_type,acc_property_memory)
      call acc_get_property_string(dev_num,dev_type,acc_property_vendor,myvendor)
      call acc_get_property_string(dev_num,dev_type,acc_property_driver,mydriver)
      
      write(iu,907)"                                                                               "
      write(iu,907)"*****************************GPU FEATURE MONITOR*******************************"
      write(iu,907)"                                                                               "
      
      write (iu,900) "Device Number: "      ,ngpus
      write (iu,901) "Device Name: "        ,trim(myname)
      write (iu,903) "Total Global Memory: ",real(tot_mem)/1e9," Gbytes"
      write (iu,901) "Vendor: "        ,trim(myvendor)
      write (iu,901) "Driver: "        ,trim(mydriver)
      
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
 !$endif  
  subroutine print_raw_sync(iframe)
  
      implicit none
      
      integer, intent(in) :: iframe
      
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
  
  subroutine print_vtk_sync(iframe)
      implicit none
      
      integer, intent(in) :: iframe
      
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
  
  subroutine print_raw_async(iframe)
  
      implicit none
      
      integer, intent(in) :: iframe
      
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
  
  subroutine print_vtk_async(iframe)
      implicit none
      
      integer, intent(in) :: iframe
      
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
