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
 
 end module

program lb_openacc
    !$if _OPENACC
    use openacc
    !$endif
    
    use prints
    
    implicit none
    
    
    integer :: i,j,k,ii,jj,kk,ll
    integer :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus
    integer :: istat,iframe
    
    logical :: lprint,lvtk,lasync,lpbc
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
    real(kind=db) :: visc_LB,uu,udotc,omega,feq,uu0
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp

    real(kind=db) :: fneq
    real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
    real(kind=db) :: pi2cssq1,pi2cssq2,pi2cssq0
    integer, parameter :: npops=18
                                                  !0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18
    integer, parameter, dimension(0:npops) :: ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
    integer, parameter, dimension(0:npops) :: ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
    integer, parameter, dimension(0:npops) :: ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
    integer, parameter, dimension(0:npops) ::opp=(/0, 2,  1, 4,  3,  6,  5,  8,   7, 10,   9, 12,  11, 14,  13, 16,  15,  18,  17/)
    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    real(kind=db), allocatable, dimension(:,:,:) :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
    real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18
    real(kind=db) :: mymemory,totmemory
    !$if _OPENACC
    integer :: devNum
    integer(acc_device_kind) :: devType
    devType = acc_get_device_type()
    devNum=acc_get_device_num(devType)
    !$endif

       
    nlinks=18 !pari!
    tau=1.5_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
#ifdef _OPENACC
    ngpus=acc_get_num_devices(acc_device_nvidia)
#else
    ngpus=0
#endif

    !*******************************user parameters and allocations**************************m
        nx=32
        ny=32
        nz=32
        nsteps=10000
        stamp=100
        fx=1.0_db*10.0**(-5)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
        lpbc=.true.
        lprint=.true.
        lvtk=.true.
        lasync=.true.
        
        allocate(f0(0:nx+1,0:ny+1,0:nz+1),f1(0:nx+1,0:ny+1,0:nz+1),f2(0:nx+1,0:ny+1,0:nz+1),f3(0:nx+1,0:ny+1,0:nz+1))
        allocate(f4(0:nx+1,0:ny+1,0:nz+1),f5(0:nx+1,0:ny+1,0:nz+1),f6(0:nx+1,0:ny+1,0:nz+1),f7(0:nx+1,0:ny+1,0:nz+1))
        allocate(f8(0:nx+1,0:ny+1,0:nz+1),f9(0:nx+1,0:ny+1,0:nz+1),f10(0:nx+1,0:ny+1,0:nz+1),f11(0:nx+1,0:ny+1,0:nz+1))
        allocate(f12(0:nx+1,0:ny+1,0:nz+1),f13(0:nx+1,0:ny+1,0:nz+1),f14(0:nx+1,0:ny+1,0:nz+1),f15(0:nx+1,0:ny+1,0:nz+1))
        allocate(f16(0:nx+1,0:ny+1,0:nz+1),f17(0:nx+1,0:ny+1,0:nz+1),f18(0:nx+1,0:ny+1,0:nz+1))
        allocate(rho(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
        allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
        allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz))
        allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz))
          allocate(velprint(1:3,1:nx,1:ny,1:nz))
          rhoprint(1:nx,1:ny,1:nz)=0.0
          velprint(1:3,1:nx,1:ny,1:nz)=0.0
        endif
        

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
        if(lpbc)then
          isfluid=1
          isfluid(:,:,1)=0 !bottom
          isfluid(:,:,nz)=0 !top
        endif
        do k=1,nz
	      do j=1,ny
		    do i=1,nx
			  if(isfluid(i,j,k).eq.1)then
			    do ll=1,npops
				  ii=i+ex(ll)
			      jj=j+ey(ll)
				  kk=k+ez(ll) 
				  if(ii.gt.0 .and. ii.lt.nx+1 .and. jj.gt.0 .and. jj.lt.ny+1 .and. kk.gt.0 .and. kk.lt.nz+1)then
				    if(isfluid(ii,jj,kk).eq.0)then
					  isfluid(i,j,k)=-1
					endif
				  endif
				enddo
			  endif
		    enddo
		  enddo
	    enddo
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

    !*************************************initial conditions ************************    
        u=0.0_db
        v=0.0_db
        w=0.0_db
        rho=1.0_db  !tot dens
        pxx=0.0_db
        pxy=0.0_db
        pxz=0.0_db
        pyy=0.0_db
        pyz=0.0_db
        pzz=0.0_db
        !do ll=0,nlinks
!        f0(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p0
!        f1(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p1
!        f2(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p1
!        f3(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p1
!        f4(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p1
!        f5(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p1
!        f6(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p1
!        f7(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f8(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f9(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f10(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f11(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f12(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f13(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f14(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f15(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f16(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f17(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
!        f18(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)*p2
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
        write(6,*) 'lpbc',lpbc
        write(6,*) 'lprint',lprint
        write(6,*) 'lvtk',lvtk
        write(6,*) 'lasync',lasync
        write(6,*) 'nsteps',nsteps
        write(6,*) 'stamp',stamp
        write(6,*) 'max fx',huge(fx)
        write(6,*) 'max fx',huge(fy)
        write(6,*) 'max fx',huge(fz)
        write(6,*) '*******************************************'
    step = 0
    
    !$acc data copy(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,isfluid,p0,p1,p2,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,rhoprint,velprint) async(1)
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
      !$acc kernels present(rhoprint,velprint,rho,u,v,w) async(1)
      !$acc loop independent collapse(3)  private(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rhoprint(i,j,k)=real(rho(i,j,k),kind=4)
            velprint(1,i,j,k)=real(u(i,j,k),kind=4)
            velprint(2,i,j,k)=real(v(i,j,k),kind=4)
            velprint(3,i,j,k)=real(w(i,j,k),kind=4)
          enddo
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
            do k=2,nz-1
                !$acc loop independent 
                do j=2,ny-1
			      if(j>=2 .and. j<=ny-1 .and. k>=2 .and. k<=nz-1)then
			        rho(1,j,k)=rho(nx-1,j,k)
			        u(1,j,k)=u(nx-1,j,k)
			        v(1,j,k)=v(nx-1,j,k)
			        w(1,j,k)=w(nx-1,j,k)
			        pxx(1,j,k)=pxx(nx-1,j,k)
			        pyy(1,j,k)=pyy(nx-1,j,k)
			        pzz(1,j,k)=pzz(nx-1,j,k)
			        pxy(1,j,k)=pxy(nx-1,j,k)
			        pxz(1,j,k)=pxz(nx-1,j,k)
			        pyz(1,j,k)=pyz(nx-1,j,k)
			        
			        rho(nx,j,k)=rho(2,j,k)
			        u(nx,j,k)=u(2,j,k)
			        v(nx,j,k)=v(2,j,k)
			        w(nx,j,k)=w(2,j,k)
			        pxx(nx,j,k)=pxx(2,j,k)
			        pyy(nx,j,k)=pyy(2,j,k)
			        pzz(nx,j,k)=pzz(2,j,k)
			        pxy(nx,j,k)=pxy(2,j,k)
			        pxz(nx,j,k)=pxz(2,j,k)
			        pyz(nx,j,k)=pyz(2,j,k)
			      endif 
			    enddo
			enddo
            !$acc end kernels
            
            !periodic along y
            !$acc kernels async(1)
            !$acc loop independent 
            do k=2,nz-1
                !$acc loop independent 
                do i=2,nx-1
			      if(i>=2 .and. i<=nx-1 .and. k>=2 .and. k<=nz-1)then
			        rho(i,1,k)=rho(i,ny-1,k)
			        u(i,1,k)=u(i,ny-1,k)
			        v(i,1,k)=v(i,ny-1,k)
			        w(i,1,k)=w(i,ny-1,k)
			        pxx(i,1,k)=pxx(i,ny-1,k)
			        pyy(i,1,k)=pyy(i,ny-1,k)
			        pzz(i,1,k)=pzz(i,ny-1,k)
			        pxy(i,1,k)=pxy(i,ny-1,k)
			        pxz(i,1,k)=pxz(i,ny-1,k)
			        pyz(i,1,k)=pyz(i,ny-1,k)
			        
			        rho(i,ny,k)=rho(i,2,k)
			        u(i,ny,k)=u(i,2,k)
			        v(i,ny,k)=v(i,2,k)
			        w(i,ny,k)=w(i,2,k)
			        pxx(i,ny,k)=pxx(i,2,k)
			        pyy(i,ny,k)=pyy(i,2,k)
			        pzz(i,ny,k)=pzz(i,2,k)
			        pxy(i,ny,k)=pxy(i,2,k)
			        pxz(i,ny,k)=pxz(i,2,k)
			        pyz(i,ny,k)=pyz(i,2,k)
			      endif
			    enddo
			enddo
			!$acc end kernels
	    endif
        
        !***********************************PRINT************************
        if(mod(step,stamp).eq.0)write(6,'(a,i8)')'step : ',step
        if(lprint)then
          if(mod(step,stamp).eq.0)then
            iframe=iframe+1
            !$acc wait(1)
            !$acc kernels present(rhoprint,velprint,rho,u,v,w) async(1)
            !$acc loop independent collapse(3)  private(i,j,k)
            do k=1,nz
              do j=1,ny
                do i=1,nx
                  rhoprint(i,j,k)=real(rho(i,j,k),kind=4)
                  velprint(1,i,j,k)=real(u(i,j,k),kind=4)
                  velprint(2,i,j,k)=real(v(i,j,k),kind=4)
                  velprint(3,i,j,k)=real(w(i,j,k),kind=4)
                enddo
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
        
        !***********************************collision + no slip + forcing: fused implementation*********
        !$acc kernels async(1)
        !$acc loop collapse(3) private(feq,fneq,uu,uu0,temp,udotc)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
                        uu=0.5_db*(u(i,j,k)**2.0_db + v(i,j,k)**2.0_db + w(i,j,k)**2.0_db)/cssq
                        !0
                        feq=p0*(rho(i,j,k)-uu)
                        f0(i,j,k)=feq + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k)))
                        
                        !1
                        uu=0.5_db*(u(i-1,j,k)**2.0_db + v(i-1,j,k)**2.0_db + w(i-1,j,k)**2.0_db)/cssq
                        udotc=u(i-1,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i-1,j,k)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i-1,j,k)-cssq*(pyy(i-1,j,k)+pzz(i-1,j,k)))
                        f1(i,j,k)=feq + fneq + fx*p1dcssq
                        
                        !2
                        uu=0.5_db*(u(i+1,j,k)**2.0_db + v(i+1,j,k)**2.0_db + w(i+1,j,k)**2.0_db)/cssq
                        udotc=u(i+1,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i+1,j,k)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i+1,j,k)-cssq*(pyy(i+1,j,k)+pzz(i+1,j,k)))
                        f2(i,j,k)=feq + fneq - fx*p1dcssq
                        
                        !3
                        uu=0.5_db*(u(i,j-1,k)**2.0_db + v(i,j-1,k)**2.0_db + w(i,j-1,k)**2.0_db)/cssq
                        udotc=v(i,j-1,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j-1,k)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j-1,k)-cssq*(pxx(i,j-1,k)+pzz(i,j-1,k)))
                        f3(i,j,k)=feq+fneq + fy*p1dcssq
                        
                        !4
                        uu=0.5_db*(u(i,j+1,k)**2.0_db + v(i,j+1,k)**2.0_db + w(i,j+1,k)**2.0_db)/cssq
                        udotc=v(i,j+1,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j+1,k)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j+1,k)-cssq*(pxx(i,j+1,k)+pzz(i,j+1,k)))
                        f4(i,j,k)=feq+fneq - fy*p1dcssq
                        
                        !7
                        uu=0.5_db*(u(i-1,j-1,k)**2.0_db + v(i-1,j-1,k)**2.0_db + w(i-1,j-1,k)**2.0_db)/cssq
                        udotc=(u(i-1,j-1,k)+v(i-1,j-1,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i-1,j-1,k)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j-1,k)+qyy*pyy(i-1,j-1,k)-cssq*pzz(i-1,j-1,k)+2.0_db*qxy_7_8*pxy(i-1,j-1,k))
                        f7(i,j,k)=feq + fneq + (fx+fy)*p2dcssq 
                        
                        !8
                        uu=0.5_db*(u(i+1,j+1,k)**2.0_db + v(i+1,j+1,k)**2.0_db + w(i+1,j+1,k)**2.0_db)/cssq
                        udotc=(u(i+1,j+1,k)+v(i+1,j+1,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i+1,j+1,k)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j+1,k)+qyy*pyy(i+1,j+1,k)-cssq*pzz(i+1,j+1,k)+2.0_db*qxy_7_8*pxy(i+1,j+1,k))
                        f8(i,j,k)=feq + fneq - (fx+fy)*p2dcssq
                        
                        !10   +1 -1
                        uu=0.5_db*(u(i+1,j-1,k)**2.0_db + v(i+1,j-1,k)**2.0_db + w(i+1,j-1,k)**2.0_db)/cssq
                        udotc=(-u(i+1,j-1,k)+v(i+1,j-1,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i+1,j-1,k)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j-1,k)+qyy*pyy(i+1,j-1,k)-cssq*pzz(i+1,j-1,k)+2.0_db*qxy_9_10*pxy(i+1,j-1,k))
                        f10(i,j,k)=feq+fneq +(fy-fx)*p2dcssq
                        
                        !9  -1 +1
                        uu=0.5_db*(u(i-1,j+1,k)**2.0_db + v(i-1,j+1,k)**2.0_db + w(i-1,j+1,k)**2.0_db)/cssq
                        udotc=(-u(i-1,j+1,k)+v(i-1,j+1,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i-1,j+1,k)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j+1,k)+qyy*pyy(i-1,j+1,k)-cssq*pzz(i-1,j+1,k)+2.0_db*qxy_9_10*pxy(i-1,j+1,k))
                        f9(i,j,k)=feq+fneq + (fx-fy)*p2dcssq

                        !5 -1
                        uu=0.5_db*(u(i,j,k-1)**2.0_db + v(i,j,k-1)**2.0_db + w(i,j,k-1)**2.0_db)/cssq
                        udotc=w(i,j,k-1)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k-1)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k-1)-cssq*(pxx(i,j,k-1)+pyy(i,j,k-1)))
                        f5(i,j,k)=feq+fneq + fz*p1dcssq
                        
                        !6 +1
                        uu=0.5_db*(u(i,j,k+1)**2.0_db + v(i,j,k+1)**2.0_db + w(i,j,k+1)**2.0_db)/cssq
                        udotc=w(i,j,k+1)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k+1)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k+1)-cssq*(pxx(i,j,k+1)+pyy(i,j,k+1)))
                        f6(i,j,k)=feq+fneq - fz*p1dcssq

                        !15  -1  -1
                        uu=0.5_db*(u(i-1,j,k-1)**2.0_db + v(i-1,j,k-1)**2.0_db + w(i-1,j,k-1)**2.0_db)/cssq
                        udotc=(u(i-1,j,k-1)+w(i-1,j,k-1))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i-1,j,k-1)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j,k-1)+qzz*pzz(i-1,j,k-1)-cssq*pyy(i-1,j,k-1)+2.0_db*qxz_15_16*pxz(i-1,j,k-1))
                        f15(i,j,k)=feq+fneq + (fx+fz)*p2dcssq 
                        
                        !16  +1  +1
                        uu=0.5_db*(u(i+1,j,k+1)**2.0_db + v(i+1,j,k+1)**2.0_db + w(i+1,j,k+1)**2.0_db)/cssq
                        udotc=(u(i+1,j,k+1)+w(i+1,j,k+1))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i+1,j,k+1)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j,k+1)+qzz*pzz(i+1,j,k+1)-cssq*pyy(i+1,j,k+1)+2.0_db*qxz_15_16*pxz(i+1,j,k+1))
                        f16(i,j,k)=feq+fneq - (fx+fz)*p2dcssq

                        !17  +1   -1
                        uu=0.5_db*(u(i+1,j,k-1)**2.0_db + v(i+1,j,k-1)**2.0_db + w(i+1,j,k-1)**2.0_db)/cssq
                        udotc=(-u(i+1,j,k-1)+w(i+1,j,k-1))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i+1,j,k-1)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j,k-1)+qzz*pzz(i+1,j,k-1)-cssq*pyy(i+1,j,k-1)+2.0_db*qxz_17_18*pxz(i+1,j,k-1))
                        f17(i,j,k)=feq+fneq +(fz-fx)*p2dcssq
                        
                        !18   -1   +1
                        uu=0.5_db*(u(i-1,j,k+1)**2.0_db + v(i-1,j,k+1)**2.0_db + w(i-1,j,k+1)**2.0_db)/cssq
                        udotc=(-u(i-1,j,k+1)+w(i-1,j,k+1))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i-1,j,k+1)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j,k+1)+qzz*pzz(i-1,j,k+1)-cssq*pyy(i-1,j,k+1)+2.0_db*qxz_17_18*pxz(i-1,j,k+1))
                        f18(i,j,k)=feq+fneq + (fx-fz)*p2dcssq

                        !11  -1  -1
                        uu=0.5_db*(u(i,j-1,k-1)**2.0_db + v(i,j-1,k-1)**2.0_db + w(i,j-1,k-1)**2.0_db)/cssq
                        udotc=(v(i,j-1,k-1)+w(i,j-1,k-1))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j-1,k-1)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j-1,k-1)+qzz*pzz(i,j-1,k-1)-cssq*pxx(i,j-1,k-1)+2.0_db*qyz_11_12*pyz(i,j-1,k-1))
                        f11(i,j,k)=feq+fneq+(fy+fz)*p2dcssq
                        
                        !12   +1  +1
                        uu=0.5_db*(u(i,j+1,k+1)**2.0_db + v(i,j+1,k+1)**2.0_db + w(i,j+1,k+1)**2.0_db)/cssq
                        udotc=(v(i,j+1,k+1)+w(i,j+1,k+1))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j+1,k+1)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j+1,k+1)+qzz*pzz(i,j+1,k+1)-cssq*pxx(i,j+1,k+1)+2.0_db*qyz_11_12*pyz(i,j+1,k+1))
                        f12(i,j,k)=feq+fneq - (fy+fz)*p2dcssq

                        !13   -1   +1
                        uu=0.5_db*(u(i,j-1,k+1)**2.0_db + v(i,j-1,k+1)**2.0_db + w(i,j-1,k+1)**2.0_db)/cssq
                        udotc=(v(i,j-1,k+1)-w(i,j-1,k+1))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j-1,k+1)+(temp + udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j-1,k+1)+qzz*pzz(i,j-1,k+1)-cssq*pxx(i,j-1,k+1)+2.0_db*qyz_13_14*pyz(i,j-1,k+1))
                        f13(i,j,k)=feq+fneq + (fy-fz)*p2dcssq
                        
                        !14  +1 -1
                        uu=0.5_db*(u(i,j+1,k-1)**2.0_db + v(i,j+1,k-1)**2.0_db + w(i,j+1,k-1)**2.0_db)/cssq
                        udotc=(v(i,j+1,k-1)-w(i,j+1,k-1))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j+1,k-1)+(temp - udotc))
                        fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j+1,k-1)+qzz*pzz(i,j+1,k-1)-cssq*pxx(i,j+1,k-1)+2.0_db*qyz_13_14*pyz(i,j+1,k-1))
                        f14(i,j,k)=feq+fneq + (fz-fy)*p2dcssq
!********************************close to boundary conditions no slip everywhere********************************!
                    elseif(isfluid(i,j,k).eq.-1)then
                        uu0=0.5_db*(u(i,j,k)**2.0_db + v(i,j,k)**2.0_db + w(i,j,k)**2.0_db)/cssq
                        !0
                        feq=p0*(rho(i,j,k)-uu0)
                        f0(i,j,k)=feq + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k)))
                        
                        !1
                        if(isfluid(i-1,j,k).eq.0)then
                          udotc=u(i,j,k)/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                          f1(i,j,k)=feq + fneq - fx*p1dcssq
                        else
                          uu=0.5_db*(u(i-1,j,k)**2.0_db + v(i-1,j,k)**2.0_db + w(i-1,j,k)**2.0_db)/cssq
                          udotc=u(i-1,j,k)/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p1*(rho(i-1,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i-1,j,k)-cssq*(pyy(i-1,j,k)+pzz(i-1,j,k)))
                          f1(i,j,k)=feq + fneq + fx*p1dcssq
                        endif
                        
                        !2
                        if(isfluid(i+1,j,k).eq.0)then
                          udotc=u(i,j,k)/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                          f2(i,j,k)=feq + fneq + fx*p1dcssq
                        else
                          uu=0.5_db*(u(i+1,j,k)**2.0_db + v(i+1,j,k)**2.0_db + w(i+1,j,k)**2.0_db)/cssq
                          udotc=u(i+1,j,k)/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p1*(rho(i+1,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i+1,j,k)-cssq*(pyy(i+1,j,k)+pzz(i+1,j,k)))
                          f2(i,j,k)=feq + fneq - fx*p1dcssq
                        endif
                        
                        !3
                        if(isfluid(i,j-1,k).eq.0)then
                          udotc=v(i,j,k)/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
                          f3(i,j,k)=feq+fneq - fy*p1dcssq
                        else
                          uu=0.5_db*(u(i,j-1,k)**2.0_db + v(i,j-1,k)**2.0_db + w(i,j-1,k)**2.0_db)/cssq
                          udotc=v(i,j-1,k)/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j-1,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j-1,k)-cssq*(pxx(i,j-1,k)+pzz(i,j-1,k)))
                          f3(i,j,k)=feq+fneq + fy*p1dcssq
                        endif
                        
                        !4
                        if(isfluid(i,j+1,k).eq.0)then
                          udotc=v(i,j,k)/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
                          f4(i,j,k)=feq+fneq + fy*p1dcssq
                        else
                          uu=0.5_db*(u(i,j+1,k)**2.0_db + v(i,j+1,k)**2.0_db + w(i,j+1,k)**2.0_db)/cssq
                          udotc=v(i,j+1,k)/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j+1,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j+1,k)-cssq*(pxx(i,j+1,k)+pzz(i,j+1,k)))
                          f4(i,j,k)=feq+fneq - fy*p1dcssq
                        endif
                        
                        !7 -1 -1
                        if(isfluid(i-1,j-1,k).eq.0)then
                          udotc=(u(i,j,k)+v(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k))
                          f7(i,j,k)=feq + fneq - (fx+fy)*p2dcssq
                        else
                          uu=0.5_db*(u(i-1,j-1,k)**2.0_db + v(i-1,j-1,k)**2.0_db + w(i-1,j-1,k)**2.0_db)/cssq
                          udotc=(u(i-1,j-1,k)+v(i-1,j-1,k))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i-1,j-1,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j-1,k)+qyy*pyy(i-1,j-1,k)-cssq*pzz(i-1,j-1,k)+2.0_db*qxy_7_8*pxy(i-1,j-1,k))
                          f7(i,j,k)=feq + fneq + (fx+fy)*p2dcssq 
                        endif
                        
                        !8 +1 +1
                        if(isfluid(i+1,j+1,k).eq.0)then
                          udotc=(u(i,j,k)+v(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k))
                          f8(i,j,k)=feq + fneq + (fx+fy)*p2dcssq 
                        else
                          uu=0.5_db*(u(i+1,j+1,k)**2.0_db + v(i+1,j+1,k)**2.0_db + w(i+1,j+1,k)**2.0_db)/cssq
                          udotc=(u(i+1,j+1,k)+v(i+1,j+1,k))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i+1,j+1,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j+1,k)+qyy*pyy(i+1,j+1,k)-cssq*pzz(i+1,j+1,k)+2.0_db*qxy_7_8*pxy(i+1,j+1,k))
                          f8(i,j,k)=feq + fneq - (fx+fy)*p2dcssq
                        endif
                        
                        !10   +1 -1
                        if(isfluid(i+1,j-1,k).eq.0)then
                          udotc=(-u(i,j,k)+v(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k))
                          f10(i,j,k)=feq+fneq + (fx-fy)*p2dcssq
                        else
                          uu=0.5_db*(u(i+1,j-1,k)**2.0_db + v(i+1,j-1,k)**2.0_db + w(i+1,j-1,k)**2.0_db)/cssq
                          udotc=(-u(i+1,j-1,k)+v(i+1,j-1,k))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i+1,j-1,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j-1,k)+qyy*pyy(i+1,j-1,k)-cssq*pzz(i+1,j-1,k)+2.0_db*qxy_9_10*pxy(i+1,j-1,k))
                          f10(i,j,k)=feq+fneq +(fy-fx)*p2dcssq
                        endif
                        
                        !9  -1 +1
                        if(isfluid(i-1,j+1,k).eq.0)then
                          udotc=(-u(i,j,k)+v(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k))
                          f9(i,j,k)=feq+fneq +(fy-fx)*p2dcssq
                        else
                          uu=0.5_db*(u(i-1,j+1,k)**2.0_db + v(i-1,j+1,k)**2.0_db + w(i-1,j+1,k)**2.0_db)/cssq
                          udotc=(-u(i-1,j+1,k)+v(i-1,j+1,k))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i-1,j+1,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j+1,k)+qyy*pyy(i-1,j+1,k)-cssq*pzz(i-1,j+1,k)+2.0_db*qxy_9_10*pxy(i-1,j+1,k))
                          f9(i,j,k)=feq+fneq + (fx-fy)*p2dcssq
                        endif
                        
                        !5 -1
                        if(isfluid(i,j,k-1).eq.0)then
                          udotc=w(i,j,k)/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
                          f5(i,j,k)=feq+fneq - fz*p1dcssq
                        else
                          uu=0.5_db*(u(i,j,k-1)**2.0_db + v(i,j,k-1)**2.0_db + w(i,j,k-1)**2.0_db)/cssq
                          udotc=w(i,j,k-1)/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j,k-1)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k-1)-cssq*(pxx(i,j,k-1)+pyy(i,j,k-1)))
                          f5(i,j,k)=feq+fneq + fz*p1dcssq
                        endif
                        
                        !6 +1
                        if(isfluid(i,j,k+1).eq.0)then
                          udotc=w(i,j,k)/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
                          f6(i,j,k)=feq+fneq - fz*p1dcssq
                        else
                          uu=0.5_db*(u(i,j,k+1)**2.0_db + v(i,j,k+1)**2.0_db + w(i,j,k+1)**2.0_db)/cssq
                          udotc=w(i,j,k+1)/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p1*(rho(i,j,k+1)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k+1)-cssq*(pxx(i,j,k+1)+pyy(i,j,k+1)))
                          f6(i,j,k)=feq+fneq - fz*p1dcssq
                        endif

                        !15  -1  -1
                        if(isfluid(i-1,j,k-1).eq.0)then
                          udotc=(u(i,j,k)+w(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k))
                          f15(i,j,k)=feq+fneq - (fx+fz)*p2dcssq
                        else
                          uu=0.5_db*(u(i-1,j,k-1)**2.0_db + v(i-1,j,k-1)**2.0_db + w(i-1,j,k-1)**2.0_db)/cssq
                          udotc=(u(i-1,j,k-1)+w(i-1,j,k-1))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i-1,j,k-1)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j,k-1)+qzz*pzz(i-1,j,k-1)-cssq*pyy(i-1,j,k-1)+2.0_db*qxz_15_16*pxz(i-1,j,k-1))
                          f15(i,j,k)=feq+fneq + (fx+fz)*p2dcssq 
                        endif
                        
                        !16  +1  +1
                        if(isfluid(i+1,j,k+1).eq.0)then
                          udotc=(u(i,j,k)+w(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k))
                          f16(i,j,k)=feq+fneq + (fx+fz)*p2dcssq
                        else
                          uu=0.5_db*(u(i+1,j,k+1)**2.0_db + v(i+1,j,k+1)**2.0_db + w(i+1,j,k+1)**2.0_db)/cssq
                          udotc=(u(i+1,j,k+1)+w(i+1,j,k+1))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i+1,j,k+1)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j,k+1)+qzz*pzz(i+1,j,k+1)-cssq*pyy(i+1,j,k+1)+2.0_db*qxz_15_16*pxz(i+1,j,k+1))
                          f16(i,j,k)=feq+fneq - (fx+fz)*p2dcssq
                        endif

                        !17  +1   -1
                        if(isfluid(i+1,j,k-1).eq.0)then
                          udotc=(-u(i,j,k)+w(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k))
                          f17(i,j,k)=feq+fneq + (fx-fz)*p2dcssq
                        else
                          uu=0.5_db*(u(i+1,j,k-1)**2.0_db + v(i+1,j,k-1)**2.0_db + w(i+1,j,k-1)**2.0_db)/cssq
                          udotc=(-u(i+1,j,k-1)+w(i+1,j,k-1))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i+1,j,k-1)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i+1,j,k-1)+qzz*pzz(i+1,j,k-1)-cssq*pyy(i+1,j,k-1)+2.0_db*qxz_17_18*pxz(i+1,j,k-1))
                          f17(i,j,k)=feq+fneq +(fz-fx)*p2dcssq
                        endif
                        
                        !18   -1   +1
                        if(isfluid(i-1,j,k+1).eq.0)then
                          udotc=(-u(i,j,k)+w(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k))
                          f18(i,j,k)=feq+fneq +(fz-fx)*p2dcssq
                        else
                          uu=0.5_db*(u(i-1,j,k+1)**2.0_db + v(i-1,j,k+1)**2.0_db + w(i-1,j,k+1)**2.0_db)/cssq
                          udotc=(-u(i-1,j,k+1)+w(i-1,j,k+1))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i-1,j,k+1)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i-1,j,k+1)+qzz*pzz(i-1,j,k+1)-cssq*pyy(i-1,j,k+1)+2.0_db*qxz_17_18*pxz(i-1,j,k+1))
                          f18(i,j,k)=feq+fneq + (fx-fz)*p2dcssq
                        endif

                        !11  -1  -1
                        if(isfluid(i,j-1,k-1).eq.0)then
                          udotc=(v(i,j,k)+w(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k))
                          f11(i,j,k)=feq+fneq - (fy+fz)*p2dcssq
                        else
                          uu=0.5_db*(u(i,j-1,k-1)**2.0_db + v(i,j-1,k-1)**2.0_db + w(i,j-1,k-1)**2.0_db)/cssq
                          udotc=(v(i,j-1,k-1)+w(i,j-1,k-1))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j-1,k-1)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j-1,k-1)+qzz*pzz(i,j-1,k-1)-cssq*pxx(i,j-1,k-1)+2.0_db*qyz_11_12*pyz(i,j-1,k-1))
                          f11(i,j,k)=feq+fneq+(fy+fz)*p2dcssq
                        endif
                        
                        !12   +1  +1
                        if(isfluid(i,j+1,k+1).eq.0)then
                          udotc=(v(i,j,k)+w(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k))
                          f12(i,j,k)=feq+fneq+(fy+fz)*p2dcssq
                        else
                          uu=0.5_db*(u(i,j+1,k+1)**2.0_db + v(i,j+1,k+1)**2.0_db + w(i,j+1,k+1)**2.0_db)/cssq
                          udotc=(v(i,j+1,k+1)+w(i,j+1,k+1))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j+1,k+1)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j+1,k+1)+qzz*pzz(i,j+1,k+1)-cssq*pxx(i,j+1,k+1)+2.0_db*qyz_11_12*pyz(i,j+1,k+1))
                          f12(i,j,k)=feq+fneq - (fy+fz)*p2dcssq
                        endif

                        !13   -1   +1
                        if(isfluid(i,j-1,k+1).eq.0)then
                          udotc=(v(i,j,k)-w(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k))
                          f13(i,j,k)=feq+fneq + (fz-fy)*p2dcssq
                        else
                          uu=0.5_db*(u(i,j-1,k+1)**2.0_db + v(i,j-1,k+1)**2.0_db + w(i,j-1,k+1)**2.0_db)/cssq
                          udotc=(v(i,j-1,k+1)-w(i,j-1,k+1))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j-1,k+1)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j-1,k+1)+qzz*pzz(i,j-1,k+1)-cssq*pxx(i,j-1,k+1)+2.0_db*qyz_13_14*pyz(i,j-1,k+1))
                          f13(i,j,k)=feq+fneq + (fy-fz)*p2dcssq
                        endif
                        
                        !14  +1 -1
                        if(isfluid(i,j+1,k-1).eq.0)then
                          udotc=(v(i,j,k)-w(i,j,k))/cssq
                          temp = -uu0 + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j,k)+(temp + udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k))
                          f14(i,j,k)=feq+fneq + (fy-fz)*p2dcssq
                        else
                          uu=0.5_db*(u(i,j+1,k-1)**2.0_db + v(i,j+1,k-1)**2.0_db + w(i,j+1,k-1)**2.0_db)/cssq
                          udotc=(v(i,j+1,k-1)-w(i,j+1,k-1))/cssq
                          temp = -uu + 0.5_db*udotc*udotc
                          feq=p2*(rho(i,j+1,k-1)+(temp - udotc))
                          fneq=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j+1,k-1)+qzz*pzz(i,j+1,k-1)-cssq*pxx(i,j+1,k-1)+2.0_db*qyz_13_14*pyz(i,j+1,k-1))
                          f14(i,j,k)=feq+fneq + (fz-fy)*p2dcssq
                        endif
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        !$acc wait(1)
        !***********************************moment + neq pressor*********
        !$acc kernels async(1) 
        !$acc loop collapse (3) private(fneq,uu,temp,udotc)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(abs(isfluid(i,j,k)).eq.1)then
                        pxx(i,j,k)=0.0_db
                        pyy(i,j,k)=0.0_db
                        pzz(i,j,k)=0.0_db
                        pxy(i,j,k)=0.0_db
                        pxz(i,j,k)=0.0_db
                        pyz(i,j,k)=0.0_db

                        rho(i,j,k) = f0(i,j,k)+f1(i,j,k)+f2(i,j,k)+f3(i,j,k)+f4(i,j,k)+f5(i,j,k) &
                                     +f6(i,j,k)+f7(i,j,k)+f8(i,j,k)+f9(i,j,k)+f10(i,j,k)+f11(i,j,k) &
                                     +f12(i,j,k)+f13(i,j,k)+f14(i,j,k)+f15(i,j,k)+f16(i,j,k)+f17(i,j,k) &
                                     +f18(i,j,k)


                        u(i,j,k) = (f1(i,j,k)+f7(i,j,k)+f9(i,j,k)+&
                                    f15(i,j,k)+f18(i,j,k)) &
                                    -(f2(i,j,k)+f8(i,j,k)+f10(i,j,k) &
                                        +f16(i,j,k)+f17(i,j,k)) 
                        
                        v(i,j,k) = (f3(i,j,k)+f7(i,j,k)+f10(i,j,k)+ &
                                    f11(i,j,k)+f13(i,j,k)) &
                                    - (f4(i,j,k)+f8(i,j,k)+f9(i,j,k) &
                                        +f12(i,j,k)+f14(i,j,k))

                        w(i,j,k) = (f5(i,j,k)+f11(i,j,k)+f14(i,j,k)+ &
                                    f15(i,j,k)+f17(i,j,k))-&
                                        (f6(i,j,k)+f12(i,j,k)+f13(i,j,k) &
                                         +f16(i,j,k)+f18(i,j,k))
                        
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !1-2
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq=f1(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        fneq=f2(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq

                        !3-4
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq=f3(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        fneq=f4(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        !5-6
                        udotc=w(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq=f5(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        fneq=f6(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        !7-8
                        udotc=(u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq=f7(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        pxy(i,j,k)=pxy(i,j,k)+fneq
                        fneq=f8(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        pxy(i,j,k)=pxy(i,j,k)+fneq
                        !10-9
                        udotc=(-u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq=f10(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        pxy(i,j,k)=pxy(i,j,k)-fneq
                        fneq=f9(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        pxy(i,j,k)=pxy(i,j,k)-fneq
                        !11-12
                        udotc=(v(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq=f11(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        pyz(i,j,k)=pyz(i,j,k)+fneq
                        fneq=f12(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        pyz(i,j,k)=pyz(i,j,k)+fneq
                        !13-14
                        udotc=(v(i,j,k)-w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq=f13(i,j,k) - p2*(rho(i,j,k)+(temp + udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        pyz(i,j,k)=pyz(i,j,k)-fneq
                        fneq=f14(i,j,k) - p2*(rho(i,j,k)+(temp - udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        pyz(i,j,k)=pyz(i,j,k)-fneq
                        !15-16
                        udotc=(u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq=f15(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        pxz(i,j,k)=pxz(i,j,k)+fneq
                        fneq=f16(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        pxz(i,j,k)=pxz(i,j,k)+fneq
                        !17-18
                        udotc=(-u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq=f17(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        pxz(i,j,k)=pxz(i,j,k)-fneq
                        fneq=f18(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq
                        pzz(i,j,k)=pzz(i,j,k)+fneq
                        pxz(i,j,k)=pxz(i,j,k)-fneq
                        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
                        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
                        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels

        !$acc wait(1)
        
    enddo 
    !$acc wait
    if(lasync)then
      if(lvtk)then
        call print_vtk_sync(iframe)
      else
        call print_raw_sync(iframe)
      endif
    endif
    call cpu_time(ts2)
    !$acc update host(rho,u,v,w)
    !$acc end data
    

    write(6,*) 'u=',u(nx/2,ny/2,nz/2),'v=',v(nx/2,ny/2,nz/2),'w=',w(nx/2,ny/2,nz/2),'rho=',rho(nx/2,ny/2,nz/2)
    write(6,*) 'u=',u(nx/2,ny/2,1),'v=',v(nx/2,ny/2,1),'w=',w(nx/2,ny/2,1),'rho=',rho(nx/2,ny/2,1)
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)/1.0e9/(ts2-ts1)
    
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
