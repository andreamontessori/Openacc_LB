
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

    endsubroutine init_output
 
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
 
 end module


program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    use prints
    
    implicit none
    !*********************************variables
        
        integer :: i,j,ll,l,dumm
        integer :: nx,ny,step,stamp,nlinks,nsteps,ngpus,ncontact
        integer,save :: iframe=0
        integer, parameter :: nz=1
        
        logical :: lprint=.true.
        logical :: lvtk=.true.
        logical :: lasync=.false.
        
        real(kind=db),parameter :: pi_greek=3.14159265359793234626433
        
        real(kind=4)  :: ts1,ts2,radius
        real(kind=db) :: visc_LB,uu,udotc,omega,feq,geq
        real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8
        real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq1,pi2cssq2,pi2cssq0
        real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,fpc
        real(kind=db) :: addendum0,addendum1,addendum2,addendum3,addendum4,addendum5,addendum6,addendum7,addendum8
        real(kind=db) :: gaddendum0,gaddendum1,gaddendum2,gaddendum3,gaddendum4,gaddendum5,gaddendum6,gaddendum7,gaddendum8
        real(kind=db) :: psi_x,psi_y,mod_psi,mod_psi_sq,st_coeff,b0,b1,b2,beta,sigma,norm_x,norm_y
        real(kind=db) :: one_ov_nu2,one_ov_nu1,nu_avg,rtot,rprod
        real(kind=db) :: max_press_excess,ushifted,vshifted
        real(kind=db) :: rnd_n1,rnd_n2,drhok1,drhok2
        
        integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
        
        real(kind=db), allocatable, dimension(:)     :: p
        real(kind=db), allocatable, dimension(:,:) :: psi,rhoA,rhoB,u,v,pxx,pyy,pxy
        real(kind=db), allocatable, dimension(:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8
        real(kind=db), allocatable, dimension(:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8
        integer(kind=4), allocatable, dimension(:,:) :: nci_loc
        !*******************************for geometry
        integer:: ddrop,lc,jjd,jju,k

    !*********************************lattice pars  
        
        !$if _OPENACC
        integer :: devNum
        integer(acc_device_kind) :: devType
        devType = acc_get_device_type()
        devNum=acc_get_device_num(devType)
        !$endif
        
        nlinks=8 !pari!
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
#ifdef _OPENACC
        ngpus=acc_get_num_devices(acc_device_nvidia)
#else
        ngpus=0
#endif
    !*******************************user parameters**************************
        nx=256!500!500
        ny=256 !500!600
        nsteps=10
        stamp=10
        fx=0.0_db*10.0**(-7)
        fy=-0.0_db*10.0**(-6)
    !**********************************allocation****************************
        allocate(p(0:nlinks))
        allocate(f0(0:nx+1,0:ny+1),f1(0:nx+1,0:ny+1),f2(0:nx+1,0:ny+1),f3(0:nx+1,0:ny+1),f4(0:nx+1,0:ny+1))
        allocate(f5(0:nx+1,0:ny+1),f6(0:nx+1,0:ny+1),f7(0:nx+1,0:ny+1),f8(0:nx+1,0:ny+1))
        allocate(g0(0:nx+1,0:ny+1),g1(0:nx+1,0:ny+1),g2(0:nx+1,0:ny+1),g3(0:nx+1,0:ny+1),g4(0:nx+1,0:ny+1))
        allocate(g5(0:nx+1,0:ny+1),g6(0:nx+1,0:ny+1),g7(0:nx+1,0:ny+1),g8(0:nx+1,0:ny+1))
        allocate(psi(0:nx+1,0:ny+1),rhoA(1:nx,1:ny),rhoB(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny),pxx(1:nx,1:ny),pyy(1:nx,1:ny),pxy(1:nx,1:ny))
        allocate(isfluid(1:nx,1:ny),nci_loc(1:nx,1:ny)) !,omega_2d(1:nx,1:ny)) 
        p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
        1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
        
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz))
          allocate(velprint(1:3,1:nx,1:ny,1:nz))
          rhoprint(1:nx,1:ny,1:nz)=0.0
          velprint(1:3,1:nx,1:ny,1:nz)=0.0
        endif
        
        !ex=(/0,1,0,-1,0,1,-1,-1,1/)
        !ey=(/0,0,1,0,-1,1,1,-1,-1/)

    !***************************** regularized: hermite 
        qxx=1.0_db-cssq
        qyy=1.0_db-cssq
        qxy5_7=1.0_db
        qxy6_8=-1.0_db
        pi2cssq0=p(0)/(2.0_db*cssq**2)
        pi2cssq1=p(1)/(2.0_db*cssq**2)
        pi2cssq2=p(5)/(2.0_db*cssq**2)
        ! chromodynamic
        beta=0.95_db
        sigma=0.02_db
        st_coeff=(9.0_db/4.0_db)*sigma*omega
        b0=-4.0_db/27.0_db
        b1=2.0_db/27.0_db
        b2=5.0_db/108.0_db
    !*****************************************geometry************************
        radius=20
        isfluid=1
        isfluid(1,:)=0 !EAST
        isfluid(nx,:)=0 !WEST
        isfluid(:,1)=0 !SOUTH 
        isfluid(:,ny)=0 !NORTH
        !************************************constriction**********************!
            ! isfluid(nx/2-radius:nx/2+radius,2:ny-1)=1
            ! isfluid(2:nx-1,2:ny/2-2*radius)=1
            ! isfluid(2:nx-1,ny/2+2*radius:ny-1)=1
            ! do k=1,ny 
            !     do j=1,nx
            !         if(isfluid(j,k).eq.3)then
            !             if(isfluid(j+1,k).eq.1 .or. isfluid(j-1,k).eq.1 .or. isfluid(j,k+1).eq.1 .or. isfluid(j,k-1).eq.1 &
            !                 .or. isfluid(j+1,k+1).eq.1 .or. isfluid(j+1,k-1).eq.1 .or. isfluid(j-1,k+1).eq.1 .or. isfluid(j-1,k-1).eq.1)then
            !                 isfluid(j,k)=0
            !             endif
            !         endif
            !     enddo
            ! enddo
        !****************************************tapered**************************
            ! Lc=250+ny/2;
            ! do i=Lc,ny
            !     jjd=nint((i-Lc+1)*sind(30.0_db));
                
            !     if(jjd<nx/2-4)then
            !         isfluid(jjd:nx/2,i)=1; 
            !         jju=nx-jjd; 
            !         isfluid(nx/2:jju,i)=1;
            !     endif
            ! enddo 
            ! Ddrop=40;
            ! isfluid(2:nx-1,ny/2:Lc)=1;
            ! isfluid(nx/2-(Ddrop/(2)):nx/2+(Ddrop/(2)),ny/2:ny)=1;
            ! isfluid(:,1:ny/2)=isfluid(:,ny:ny/2+1:-1);
            ! isfluid(1,:)=0 !left
            ! isfluid(nx,:)=0 !right
            ! isfluid(:,1)=0 !front 
            ! isfluid(:,ny)=0 !rear z
            ! do k=1,ny 
            !     do j=1,nx
            !         if(isfluid(j,k).eq.3)then
            !             if(isfluid(j+1,k).eq.1 .or. isfluid(j-1,k).eq.1 .or. isfluid(j,k+1).eq.1 .or. isfluid(j,k-1).eq.1 &
            !                 .or. isfluid(j+1,k+1).eq.1 .or. isfluid(j+1,k-1).eq.1 .or. isfluid(j-1,k+1).eq.1 .or. isfluid(j-1,k-1).eq.1)then
            !                 isfluid(j,k)=0
            !             endif
            !         endif
            !     enddo
            ! enddo
        !****************************print geo**************************!
            !isfluid(2:nx-1,2:800)=1 !rear z
            open(231, file = 'isfluid.out', status = 'replace')
                do j=1,nx
                    do k=1,ny
                        write(231,*) isfluid(j,k)  
                    enddo
                enddo
            close(231)
    !*************************************initial conditions ************************    
        u=0.0_db
        v=0.0_db
        rhoA=0.0_db     !is to be intended as a delta rho
        rhoB=0.0_db     !is to be intended as a delta rho
        max_press_excess=0.002 !2
        nci_loc=0
        psi=-1.0_db
        !*******************************************impacting drops****************************!
            do i=(nx/2-radius-5)-radius,(nx/2-radius-5) +radius
                do j=ny/2-radius,ny/2+radius
                    if ((i-(nx/2-radius-5))**2+(j-ny/2)**2<=radius**2)then
                        psi(i,j)=1.0_db
                        u(i,j)=0.0
                    endif
                enddo
            enddo
            ! do i=(nx/2+radius+5)-radius,(nx/2+radius+5) +radius
            !     do j=ny/2-radius,ny/2+radius
            !         if ((i-(nx/2+radius+5))**2+(j-ny/2)**2<=radius**2)then
            !             psi(i,j)=1.0_db
            !             u(i,j)=-0.0
            !         endif
            !     enddo
            ! enddo
        !***************************************constriction drops************************************!
            ! do i=(nx/2-radius-5)-radius,(nx/2-radius-5) +radius
            !     do j=2*radius+10-2*radius,2*radius+10+2*radius
            !         if ((i-(nx/2-radius-5))**2+(j-2*radius+10)**2<=radius**2)then
            !             psi(i,j)=1.0_db
            !         endif
            !     enddo
            ! enddo
            ! do i=(nx/2+radius+5)-radius,(nx/2+radius+5) +radius
            !     do j=2*radius+10-2*radius,2*radius+10+2*radius
            !         if ((i-(nx/2+radius+5))**2+(j-2*radius+10)**2<=radius**2)then
            !             psi(i,j)=1.0_db
            !         endif
            !     enddo
            ! enddo
        !*********************************************tapered single drops**********************************!
            ! do i=nx/2-radius,nx/2+radius
            !     do j=ny-20-radius-radius,ny-20-radius+radius
            !         if ((i-nx/2)**2+(j-(ny-20-radius))**2<=radius**2)then
            !             psi(i,j)=1.0_db
            !         endif
            !     enddo
            ! enddo
        !*********************************************tapered full system**********************!
            ! open(231, file = 'psi.txt', status = 'old',action='read')
            ! do j=1,nx
            !     do k=1,ny
            !         read(231,*) psi(j,k)
            !     enddo
            ! enddo
            ! close(231)
            ! !
            ! open(231, file = 'psi.out', status = 'replace')
            ! do i=1,nx
            !     do j=1,ny
            !         write(231,*) psi(i,j)  
            !     enddo
            ! enddo
            ! close(231)
    !*****************************************init distros*********************************!
        rhoB=0.5*(1.0_db-psi(1:nx,1:ny))
        rhoA=1.0_db-rhoB
        write(*,*) rhoB(nx/2,ny/2),rhoA(nx/2,ny/2)
        !do ll=0,nlinks
        f0(1:nx,1:ny)=p(0)*rhoA(:,:)*(1.0-u(:,:)*u(:,:)*0.5/cssq)
        f1(1:nx,1:ny)=p(1)*rhoA(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        f2(1:nx,1:ny)=p(2)*rhoA(:,:)*(1.0+v(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        f3(1:nx,1:ny)=p(3)*rhoA(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        f4(1:nx,1:ny)=p(4)*rhoA(:,:)*(1.0-v(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        f5(1:nx,1:ny)=p(5)*rhoA(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        f6(1:nx,1:ny)=p(6)*rhoA(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        f7(1:nx,1:ny)=p(7)*rhoA(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        f8(1:nx,1:ny)=p(8)*rhoA(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        !
        g0(1:nx,1:ny)=p(0)*rhoB(:,:)*(1.0-u(:,:)*u(:,:)*0.5/cssq)
        g1(1:nx,1:ny)=p(1)*rhoB(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        g2(1:nx,1:ny)=p(2)*rhoB(:,:)*(1.0+v(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        g3(1:nx,1:ny)=p(3)*rhoB(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        g4(1:nx,1:ny)=p(4)*rhoB(:,:)*(1.0-v(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        g5(1:nx,1:ny)=p(5)*rhoB(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        g6(1:nx,1:ny)=p(6)*rhoB(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        g7(1:nx,1:ny)=p(7)*rhoB(:,:)*(1.0-u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)
        g8(1:nx,1:ny)=p(8)*rhoB(:,:)*(1.0+u(:,:)/cssq + 0.5*(u(:,:)/cssq)**2-u(:,:)*u(:,:)*0.5/cssq)

    !*************************************check data ************************ 
        write(6,*) '*******************LB data*****************'
        write(6,*) 'tau',tau
        write(6,*) 'omega',omega
        write(6,*) 'visc',one_ov_nu1,one_ov_nu2
        write(6,*) 'fx',fx
        write(6,*) 'cssq',cssq
        write(6,*) 'beta',beta
        write(6,*) 'sigma',sigma
        write(6,*) 'surface_tens_coeff',st_coeff
        write(6,*) '*******************INPUT data*****************'
        write(6,*) 'nx',nx
        write(6,*) 'ny',ny
        write(6,*) 'nsteps',nsteps
        write(6,*) 'stamp',stamp
        write(6,*) 'lprint',lprint
        write(6,*) 'lvtk',lvtk
        write(6,*) 'lasync',lasync
        write(6,*) 'max fx',huge(fx)
        write(6,*) 'available gpus',ngpus
        write(6,*) '*******************************************'
      

    !*****************************************gpu copies**********************!
        !$acc data copy(p,rhoA,rhoB,u,v,pxx,pxy,pyy,f0,f1,f2,f3,f4,f5,f6,f7,f8,isfluid, &
        !$acc& g0,g1,g2,g3,g4,g5,g6,g7,g8,psi,nci_loc,rhoprint,velprint)
    !*************************************time loop************************  
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
      !$acc kernels present(rhoprint,velprint,rhoa,u,v) async(1)
      !$acc loop independent collapse(2)  private(i,j)
      do j=1,ny
        do i=1,nx
          rhoprint(i,j,1)=real(rhoA(i,j),kind=4)
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
          call print_vtk_sync
        else
          call print_raw_sync
        endif
      endif
    endif
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moment + neq pressor*********
            !$acc kernels 
            !$acc loop collapse(2) !private(uu,temp,udotc,fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8) 
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j).eq.1)then
                        rhoA(i,j) = (f0(i,j)+f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j))
                        rhoB(i,j)=  (g0(i,j)+g1(i,j)+g2(i,j)+g3(i,j)+g4(i,j)+g5(i,j)+g6(i,j)+g7(i,j)+g8(i,j))
                        
                        u(i,j) = (f1(i,j)+f5(i,j) +f8(i,j)-f3(i,j) -f6(i,j) -f7(i,j) + &
                                g1(i,j)+g5(i,j) +g8(i,j)-g3(i,j) -g6(i,j) -g7(i,j))!/(rhoa(i,j)+rhoB(i,j))
                        v(i,j) = (f5(i,j) +f2(i,j) +f6(i,j)-f7(i,j) -f4(i,j) -f8(i,j) + &
                                g5(i,j) +g2(i,j) +g6(i,j)-g7(i,j) -g4(i,j) -g8(i,j) )!/(rhoa(i,j)+rhoB(i,j))

                        psi(i,j)= (rhoA(i,j)-rhoB(i,j))/(rhoA(i,j)+rhoB(i,j))  
                        nci_loc(i,j)=0    

                        ! non equilibrium pressor components
                        uu=0.5_db*(u(i,j)*u(i,j) + v(i,j)*v(i,j))/cssq
                        !1-3
                        udotc=u(i,j)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        rtot=rhoA(i,j)+rhoB(i,j)
                        fneq1=f1(i,j)+g1(i,j)-p(1)*(rtot + (temp + udotc))
                        fneq3=f3(i,j)+g3(i,j)-p(3)*(rtot + (temp - udotc))
                        !2-4
                        udotc=v(i,j)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq2=f2(i,j)+g2(i,j)-p(2)*(rtot + (temp + udotc))
                        fneq4=f4(i,j)+g4(i,j)-p(4)*(rtot + (temp - udotc))
                        !5-7
                        udotc=(u(i,j)+v(i,j))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq5=f5(i,j)+g5(i,j)-p(5)*(rtot + (temp + udotc))
                        fneq7=f7(i,j)+g7(i,j)-p(7)*(rtot + (temp - udotc))
                        !6-8
                        udotc=(-u(i,j)+v(i,j))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq6=f6(i,j)+g6(i,j)-p(6)*(rtot + (temp + udotc))
                        fneq8=f8(i,j)+g8(i,j)-p(8)*(rtot + (temp - udotc))
                        !
                        pxx(i,j)= fneq1 + fneq3 + fneq5 + fneq6 + fneq7 + fneq8
                        pyy(i,j)= fneq2 + fneq4 + fneq5 + fneq6 + fneq7 + fneq8
                        pxy(i,j)= fneq5 - fneq6 + fneq7 - fneq8
                    endif
                enddo
            enddo
            !$acc end kernels 
            
            if(mod(step,stamp).eq.0)write(6,'(a,i8)')'step : ',step
            if(lprint)then
              if(mod(step,stamp).eq.0)then
                iframe=iframe+1
                !$acc wait(1)
                !$acc kernels present(rhoprint,velprint,rhoa,u,v) async(1)
                !$acc loop independent collapse(2)  private(i,j)
                  do j=1,ny
                    do i=1,nx
                      rhoprint(i,j,1)=real(rhoA(i,j),kind=4)
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
                    call print_vtk_sync
                  else
                    call print_raw_sync
                  endif
                endif
              endif
              if(mod(step-stamp/4,stamp).eq.0 .and. lasync)then
                !write(6,*)'ciao 2',step,iframe
                !$acc wait(2)  
                if(lvtk)then
                  call print_vtk_async
                else
                  call print_raw_async
                endif
              endif
            endif
        !********************collision + no slip + forcing: fused implementation*********
            !$acc kernels 
            !$acc loop collapse(2)  
            do j=1,ny
                do i=1,nx 
                    if(isfluid(i,j).eq.1)then
                        !************************* pressure excess **********!
                        if(psi(i,j).lt.-0.9_db .and. i.gt.3 .and. i.lt.nx-2 .and. j.gt.3 .and. j.lt.ny-2)then
                            
                            if (psi(i+3,j).gt.-0.85 .and. psi(i-3,j).gt.-0.85)then !&
                                nci_loc(i+3,j)=1
                                nci_loc(i-3,j)=1
                            endif
                            if (psi(i,j+3).gt.-0.85 .and. psi(i,j-3).gt.-0.85)then
                               nci_loc(i,j+3)=1
                               nci_loc(i,j-3)=1
                            endif
                            if (psi(i+2,j+2).gt.-0.85 .and. psi(i-2,j-2).gt.-0.85 )then
                                nci_loc(i+2,j+2)=1
                                nci_loc(i-2,j-2)=1
                            endif
                            if (psi(i-2,j+2).gt.-0.85 .and. psi(i+2,j-2).gt.-0.85 )then
                                nci_loc(i-2,j+2)=1
                                nci_loc(i+2,j-2)=1
                            endif
                        endif
                        !*************************pressure excess**********!
                    endif
                enddo
            enddo
            !$acc loop collapse(2)  
            do j=1,ny
                do i=1,nx 
                    if(isfluid(i,j).eq.1)then 
                        !oneminusuu= -uu !1.0_db - uu
                        !0
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
                enddo
            enddo
        
        !********************************************bcs no slip*****************************************!
            !$acc loop independent 
            do j=1,ny
                !$acc loop independent 
                do i=1,nx
                    if(isfluid(i,j).eq.0)then
                        psi(i,j)=-1.0_db

                        f8(i+1,j-1)=f6(i,j)!gpc 
                        f7(i-1,j-1)=f5(i,j)!hpc

                        f6(i-1,j+1)=f8(i,j)!gpc 
                        f5(i+1,j+1)=f7(i,j)!hpc 


                        f4(i,j-1)=f2(i,j)!gpc 
                        f3(i-1,j)=f1(i,j)!hpc 

                        f2(i,j+1)=f4(i,j)!gpc 
                        f1(i+1,j)=f3(i,j)!hpc 
                        !************************!
                        g8(i+1,j-1)=g6(i,j)!gpc 
                        g7(i-1,j-1)=g5(i,j)!hpc

                        g6(i-1,j+1)=g8(i,j)!gpc 
                        g5(i+1,j+1)=g7(i,j)!hpc 


                        g4(i,j-1)=g2(i,j)!gpc 
                        g3(i-1,j)=g1(i,j)!hpc 

                        g2(i,j+1)=g4(i,j)!gpc 
                        g1(i+1,j)=g3(i,j)!hpc
                    endif
                enddo
            enddo
        
        !******************************************call periodic bcs************************
                !   !$acc loop independent 
                !    do i=1,nx  
                !     psi(i,ny)=psi(i,2)
                !     psi(i,1)=psi(i,ny-1)
                !     !negative lungo y
                !     f4(i,ny-1)=f4(i,1)
                !     f7(i,ny-1)=f7(i,1)
                !     f8(i,ny-1)=f8(i,1)

                !     g4(i,ny-1)=g4(i,1)
                !     g7(i,ny-1)=g7(i,1)
                !     g8(i,ny-1)=g8(i,1)

                !     f2(i,2)=f2(i,ny)
                !     f5(i,2)=f5(i,ny)
                !     f6(i,2)=f6(i,ny)

                !     g2(i,2)=g2(i,ny)
                !     g5(i,2)=g5(i,ny)
                !     g6(i,2)=g6(i,ny)
                    
                ! !     psi(i,1)=psi(i,2)
                ! !     !pos lungo y
                ! !     f2(i,2)=0.0
                ! !     f5(i,2)=0.0
                ! !     f6(i,2)=0.0

                !     ! g2(i,2)=p(2)*1.0_db 
                !     ! g5(i,2)=p(5)*1.0_db 
                !     ! g6(i,2)=p(6)*1.0_db  
                    
                !     ! g4(i,ny-1)=p(2)*1.0_db 
                !     ! g8(i,ny-1)=p(5)*1.0_db 
                !     ! g7(i,ny-1)=p(6)*1.0_db 

                !   enddo
            !$acc end kernels 
        !****************************************writeonfile***************************************************!
!            if(mod(step,stamp).eq.0)then
!            !$acc update self(psi(1:nx,1:ny),rhoB(1:nx,1:ny),rhoA(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny))
!           ! iframe=iframe+1
!            open(101, file = 'psi'//write_fmtnumb(iframe)//'.out', status = 'replace')
!            open(102, file = 'rhoA'//write_fmtnumb(iframe)//'.out', status = 'replace')
!            open(103, file = 'rhoB'//write_fmtnumb(iframe)//'.out', status = 'replace')
!            open(104, file = 'u'//write_fmtnumb(iframe)//'.out', status = 'replace')
!            open(105, file = 'v'//write_fmtnumb(iframe)//'.out', status = 'replace')
!            do i=1,nx
!                do j=1,ny
!                    write(101,*) psi(i,j)  
!                enddo
!            enddo
!            close(101)
!            do i=1,nx
!                do j=1,ny
!                    write(102,*) rhoA(i,j)  
!                enddo
!            enddo
!            close(102)
!            do i=1,nx
!                do j=1,ny
!                    write(103,*) rhoB(i,j)  
!                enddo
!            enddo
!            close(103)
!            do i=1,nx
!                do j=1,ny
!                    write(104,*) u(i,j)  
!                enddo
!            enddo
!            close(104)
!            do i=1,nx
!                do j=1,ny
!                    write(105,*) v(i,j)  
!                enddo
!            enddo
!            close(105)
!            write(6,*) "files updated at t=", step
!            endif
        !
    enddo 
    call cpu_time(ts2)
    if(lasync)then
      !$acc wait(2) 
      if(lvtk)then
        call print_vtk_sync
      else
        call print_raw_sync
      endif
    endif
    !$acc end data

    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nsteps)*real(1.d-9,kind=db)/(ts2-ts1)

  contains 
  !*****************************************************functions********************************************************!


    
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
   write(345)head1,ndatavtk(1),rhoprint(1:nx,1:ny,1:nz),footervtk(1)
   close(345)
   open(unit=346,file=trim(sevt2), &
    status='replace',action='write',access='stream',form='unformatted')
   write(346)head2,ndatavtk(2),velprint(1:3,1:nx,1:ny,1:nz),footervtk(2)
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
