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
    
    
    integer :: i,j,k
    integer :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus
    integer :: istat,iframe
    
    logical :: lprint,lvtk,lasync,lpbc
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp,init_rho

    real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq17
    real(kind=db) :: fneq9,fneq10,fneq11,fneq12,fneq13,fneq14,fneq15,fneq16,fneq18
    real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
    real(kind=db) :: pi2cssq1,pi2cssq2,pi2cssq0
    
    integer :: TILE_DIMx,TILE_DIMy,TILE_DIMz
    integer :: nxblock,nyblock,nzblock,nblocks,nxyblock
    integer :: idblock,xblock,yblock,zblock
    integer :: ii,jj,kk
    
    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    real(kind=db), allocatable, dimension(:,:,:) :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
    real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18
    real(kind=db), allocatable, dimension(:,:,:,:,:) :: fpops,hfields
    
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
        nx=16
        ny=16
        nz=16
        nsteps=1000
        stamp=100
        fx=1.0_db*10.0**(-5)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
        lpbc=.true.
        lprint=.true.
        lvtk=.true.
        lasync=.false.
        
        init_rho=1.0_db
        
        TILE_DIMx=8
        TILE_DIMy=8
        TILE_DIMz=4
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
        
        nxblock=nx/TILE_DIMx
        nyblock=ny/TILE_DIMx
        nzblock=nz/TILE_DIMz
        
        nxyblock=nxblock*nyblock
        
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
        
        allocate(fpops(TILE_DIMx,TILE_DIMy,TILE_DIMz,0:18,nblocks))
        allocate(hfields(TILE_DIMx,TILE_DIMy,TILE_DIMz,10,nblocks))
        
        allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz))
          allocate(velprint(1:3,1:nx,1:ny,1:nz))
          rhoprint(1:nx,1:ny,1:nz)=0.0
          velprint(1:3,1:nx,1:ny,1:nz)=0.0
        endif
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)

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
        
        do zblock=1,nzblock
          do yblock=1,nyblock
            do xblock=1,nxblock
              idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,1,idblock)=0.0_db  !pxx
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,2,idblock)=0.0_db  !pyy
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,3,idblock)=0.0_db  !pzz
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,4,idblock)=0.0_db  !pxy
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,5,idblock)=0.0_db  !pxz
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,6,idblock)=0.0_db  !pyz
              
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,7,idblock)=0.0_db  !u
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,8,idblock)=0.0_db  !v
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,9,idblock)=0.0_db  !w
              
              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,10,idblock)=init_rho  !rho
              
              
              fpops(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,0,idblock)=init_rho*p0
              fpops(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,1:6,idblock)=init_rho*p1
              fpops(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,7:18,idblock)=init_rho*p2
              
            enddo
          enddo
        enddo
        
        write(6,*)'AHO',nblocks,nxblock,nyblock,nzblock
        do k=1,nz
          do j=1,ny
            do i=1,nx
              xblock=(i+TILE_DIMx-1)/TILE_DIMx
              yblock=(j+TILE_DIMy-1)/TILE_DIMy
              zblock=(k+TILE_DIMz-1)/TILE_DIMz
              idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
              ii=i-xblock*TILE_DIMx+TILE_DIMx
              jj=j-yblock*TILE_DIMy+TILE_DIMy
              kk=k-zblock*TILE_DIMz+TILE_DIMz
              if(i==nx/2 .and. j==ny/2 .and. k==nz/2)then
                init_rho=666.0_db
              else
                init_rho=1.0_db
              endif
              hfields(ii,jj,kk,1,idblock)=0.0_db  !pxx
              hfields(ii,jj,kk,2,idblock)=0.0_db  !pyy
              hfields(ii,jj,kk,3,idblock)=0.0_db  !pzz
              hfields(ii,jj,kk,4,idblock)=0.0_db  !pxy
              hfields(ii,jj,kk,5,idblock)=0.0_db  !pxz
              hfields(ii,jj,kk,6,idblock)=0.0_db  !pyz
              
              hfields(ii,jj,kk,7,idblock)=0.0_db  !u
              hfields(ii,jj,kk,8,idblock)=0.0_db  !v
              hfields(ii,jj,kk,9,idblock)=0.0_db  !w
              
              hfields(ii,jj,kk,10,idblock)=init_rho  !rho
              
              
              fpops(ii,jj,kk,0,idblock)=init_rho*p0
              fpops(ii,jj,kk,1:6,idblock)=init_rho*p1
              fpops(ii,jj,kk,7:18,idblock)=init_rho*p2
              
              write(6,'(10i8)')i,j,k,ii,jj,kk,xblock,yblock,zblock,idblock
              
            enddo
          enddo
        enddo
        
       stop
              
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
        
        
    !$acc data copy(fpops,isfluid,hfields,rhoprint,velprint,&
    !$acc& nx,ny,nz,fx,fy,fz,p0,p1,p2,tau,omega,visc_LB,cssq,&
    !$acc& p1dcssq,p2dcssq,pi2cssq0,pi2cssq1,pi2cssq2,qxx,qyy,qzz,qxy_7_8,&
    !$acc& qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14,nxblock,nyblock,&
    !$acc& nzblock,nxyblock,nblocks,TILE_DIMx,TILE_DIMy,TILE_DIMz) async(1)
    
    
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
      !$acc kernels present(rhoprint,velprint,hfields,TILE_DIMx,TILE_DIMy,TILE_DIMz) async(1)
      !$acc loop tile(TILE_DIMx,TILE_DIMy,TILE_DIMz) private(i,j,k,&
      !$acc& xblock,yblock,zblock,idblock,ii,jj,kk)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            xblock=(i+TILE_DIMx-1)/TILE_DIMx
            yblock=(j+TILE_DIMy-1)/TILE_DIMy
            zblock=(k+TILE_DIMz-1)/TILE_DIMz
            idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
            ii=i-xblock*TILE_DIMx+TILE_DIMx
            jj=j-yblock*TILE_DIMy+TILE_DIMy
            kk=k-zblock*TILE_DIMz+TILE_DIMz
            rhoprint(i,j,k)=real(hfields(ii,jj,kk,10,idblock),kind=4)
            velprint(1,i,j,k)=real(hfields(ii,jj,kk,7,idblock),kind=4)
            velprint(2,i,j,k)=real(hfields(ii,jj,kk,8,idblock),kind=4)
            velprint(3,i,j,k)=real(hfields(ii,jj,kk,9,idblock),kind=4)
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
    
    !$acc wait
    stop
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moments collision bbck + forcing************************ 
        !$acc kernels async(1)
        !$acc loop collapse(3) private(fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq9,fneq10,fneq11,&
        !$acc& fneq12,fneq3,fneq14,fneq15,uu,temp,udotc)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
                        rho(i,j,k) = f0(i,j,k)+f1(i,j,k)+f2(i,j,k)+f3(i,j,k)+f4(i,j,k)+f5(i,j,k) &
                            +f6(i,j,k)+f7(i,j,k)+f8(i,j,k)+f9(i,j,k)+f10(i,j,k)+f11(i,j,k) &
                            +f12(i,j,k)+f13(i,j,k)+f14(i,j,k)+f15(i,j,k)+f16(i,j,k)+f17(i,j,k) &
                            +f18(i,j,k)

                        u(i,j,k) = (f1(i,j,k)+f7(i,j,k)+f9(i,j,k)+f15(i,j,k)+f18(i,j,k)) &
                             -(f2(i,j,k)+f8(i,j,k)+f10(i,j,k)+f16(i,j,k)+f17(i,j,k)) 
                        
                        v(i,j,k) = (f3(i,j,k)+f7(i,j,k)+f10(i,j,k)+f11(i,j,k)+f13(i,j,k)) &
                            -(f4(i,j,k)+f8(i,j,k)+f9(i,j,k)+f12(i,j,k)+f14(i,j,k))

                        w(i,j,k) = (f5(i,j,k)+f11(i,j,k)+f14(i,j,k)+f15(i,j,k)+f17(i,j,k)) &
                            -(f6(i,j,k)+f12(i,j,k)+f13(i,j,k)+f16(i,j,k)+f18(i,j,k))
                        
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !1-2
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f1(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        fneq2=f2(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        !3-4
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq3=f3(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        fneq4=f4(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        !5-6
                        udotc=w(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq5=f5(i,j,k)-p1*(rho(i,j,k)+(temp + udotc))
                        fneq6=f6(i,j,k)-p1*(rho(i,j,k)+(temp - udotc))
                        !7-8
                        udotc=(u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq7=f7(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq8=f8(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        !10-9
                        udotc=(-u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq10=f10(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq9=f9(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        !11-12
                        udotc=(v(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq11=f11(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq12=f12(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        !13-14
                        udotc=(v(i,j,k)-w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq13=f13(i,j,k) - p2*(rho(i,j,k)+(temp + udotc))
                        fneq14=f14(i,j,k) - p2*(rho(i,j,k)+(temp - udotc))
                        !15-16
                        udotc=(u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq15=f15(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq16=f16(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        !17-18
                        udotc=(-u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq17=f17(i,j,k)-p2*(rho(i,j,k)+(temp + udotc))
                        fneq18=f18(i,j,k)-p2*(rho(i,j,k)+(temp - udotc))
                        pxx(i,j,k)=fneq1+fneq2+fneq7+fneq8+fneq9+fneq10+fneq15+fneq16+fneq17+fneq18
                        pyy(i,j,k)=fneq3+fneq4+fneq7+fneq8+fneq9+fneq10+fneq11+fneq12+fneq13+fneq14
                        pzz(i,j,k)=fneq5+fneq6+fneq11+fneq12+fneq13+fneq14+fneq15+fneq16+fneq17+fneq18
                        pxy(i,j,k)= fneq7+fneq8-fneq9-fneq10
                        pxz(i,j,k)=fneq15+fneq16-fneq17-fneq18
                        pyz(i,j,k)=fneq11+fneq12-fneq13-fneq14
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        
        
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
        !$acc loop collapse(3) private(feq,uu,temp,udotc)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !0
                        feq=p0*(rho(i,j,k)-uu)
                        f0(i,j,k)=feq + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k)))
                        
                        !1
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq1=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                        f1(i+1,j,k)=feq + fneq1 + fx*p1dcssq
                        
                        !2
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f2(i-1,j,k)=feq + fneq1 - fx*p1dcssq
                        
                        !3
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq3=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
                        f3(i,j+1,k)=feq+fneq3 + fy*p1dcssq
                        
                        !4
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f4(i,j-1,k)=feq+fneq3 - fy*p1dcssq
                        
                        !7
                        udotc=(u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq7=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k))
                        f7(i+1,j+1,k)=feq + fneq7 + (fx+fy)*p2dcssq 
                        
                        !8
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f8(i-1,j-1,k)=feq + fneq7 - (fx+fy)*p2dcssq
                        
                        !10
                        udotc=(-u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq10=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k))
                        f10(i-1,j+1,k)=feq+fneq10 +(fy-fx)*p2dcssq
                        
                        !9
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f9(i+1,j-1,k)=feq+fneq10 + (fx-fy)*p2dcssq

                        !5
                        udotc=w(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rho(i,j,k)+(temp + udotc))
                        fneq5=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
                        f5(i,j,k+1)=feq+fneq5 + fz*p1dcssq
                        
                        !6
                        feq=p1*(rho(i,j,k)+(temp - udotc))
                        f6(i,j,k-1)=feq+fneq5 - fz*p1dcssq

                        !15
                        udotc=(u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq15=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k))
                        f15(i+1,j,k+1)=feq+fneq15 + (fx+fz)*p2dcssq 
                        
                        !16
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f16(i-1,j,k-1)=feq+fneq15 - (fx+fz)*p2dcssq

                        !17
                        udotc=(-u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq17=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k))
                        f17(i-1,j,k+1)=feq+fneq17 +(fz-fx)*p2dcssq
                        
                        !18
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f18(i+1,j,k-1)=feq+fneq17 + (fx-fz)*p2dcssq

                        !11
                        udotc=(v(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq11=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k))
                        f11(i,j+1,k+1)=feq+fneq11+(fy+fz)*p2dcssq
                        
                        !12
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f12(i,j-1,k-1)=feq+fneq11 - (fy+fz)*p2dcssq

                        !13
                        udotc=(v(i,j,k)-w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(rho(i,j,k)+(temp + udotc))
                        fneq13=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k))
                        f13(i,j+1,k-1)=feq+fneq13 + (fy-fz)*p2dcssq
                        
                        !14
                        feq=p2*(rho(i,j,k)+(temp - udotc))
                        f14(i,j-1,k+1)=feq+fneq13 + (fz-fy)*p2dcssq
                    endif
                enddo
            enddo
        enddo
        !********************************boundary conditions no slip everywhere********************************!
            !$acc loop independent 
            do k=1,nz
                !$acc loop independent 
                do j=1,ny
                    !$acc loop independent 
                    do i=1,nx
                        if(isfluid(i,j,k).eq.0)then
                            f18(i+1,j,k-1)=f17(i,j,k) !gpc 
                            f17(i-1,j,k+1)=f18(i,j,k) !hpc

                            f16(i-1,j,k-1)=f15(i,j,k) !gpc 
                            f15(i+1,j,k+1)=f16(i,j,k) !hpc

                            f14(i,j-1,k+1)=f13(i,j,k)!gpc 
                            f13(i,j+1,k-1)=f14(i,j,k)!hpc
                            
                            f12(i,j-1,k-1)=f11(i,j,k)!gpc 
                            f11(i,j+1,k+1)=f12(i,j,k)!hpc

                            f10(i-1,j+1,k)=f9(i,j,k)!gpc 
                            f9(i+1,j-1,k)=f10(i,j,k)!hpc

                            f8(i-1,j-1,k)=f7(i,j,k)!gpc 
                            f7(i+1,j+1,k)=f8(i,j,k)!hpc

                            f6(i,j,k-1)=f5(i,j,k)!gpc 
                            f5(i,j,k+1)=f6(i,j,k)!hpc 


                            f4(i,j-1,k)=f3(i,j,k)!gpc 
                            f3(i,j+1,k)=f4(i,j,k)!hpc 

                            f2(i-1,j,k)=f1(i,j,k)!gpc 
                            f1(i+1,j,k)=f2(i,j,k)!hpc 
                        endif
                    enddo
                enddo
            enddo
        !$acc end kernels
        
        !******************************************call other bcs************************
         if(lpbc)then      
            !periodic along x 
            !$acc kernels async(1)
            !$acc loop independent 
            do k=2,nz-1
                !$acc loop independent 
                do j=2,ny-1
			      if(j>2 .and. j<ny-1 .and. k>2 .and. k<nz-1)then
			        f1(2,j,k)=f1(nx,j,k)
			        f7(2,j,k)=f7(nx,j,k)
			        f9(2,j,k)=f9(nx,j,k)
			        f15(2,j,k)=f15(nx,j,k)
			        f18(2,j,k)=f18(nx,j,k)
			        f2(nx-1,j,k)=f2(1,j,k)
			        f8(nx-1,j,k)=f8(1,j,k)
			        f10(nx-1,j,k)=f10(1,j,k)
			        f16(nx-1,j,k)=f16(1,j,k)
			        f17(nx-1,j,k)=f17(1,j,k)
			      else
			        if(j==2)then
			          if(k==2)then
			            f1(2,j,k)=f1(nx,j,k)
			            f9(2,j,k)=f9(nx,j,k)
			            f18(2,j,k)=f18(nx,j,k)
			            f2(nx-1,j,k)=f2(1,j,k)
			            f8(nx-1,j,k)=f8(1,j,k)
			            f16(nx-1,j,k)=f16(1,j,k)
			          elseif(k==nz-1)then
			            f1(2,j,k)=f1(nx,j,k)
			            f9(2,j,k)=f9(nx,j,k)
			            f15(2,j,k)=f15(nx,j,k)
			            f2(nx-1,j,k)=f2(1,j,k)
			            f8(nx-1,j,k)=f8(1,j,k)
			            f17(nx-1,j,k)=f17(1,j,k)
			          else
			            f1(2,j,k)=f1(nx,j,k)
			            f9(2,j,k)=f9(nx,j,k)
			            f15(2,j,k)=f15(nx,j,k)
			            f18(2,j,k)=f18(nx,j,k)
			            f2(nx-1,j,k)=f2(1,j,k)
			            f8(nx-1,j,k)=f8(1,j,k)
			            f16(nx-1,j,k)=f16(1,j,k)
			            f17(nx-1,j,k)=f17(1,j,k)
			          endif
			        elseif(j==ny-1)then
			          if(k==2)then
			            f1(2,j,k)=f1(nx,j,k)
			            f7(2,j,k)=f7(nx,j,k)
			            f18(2,j,k)=f18(nx,j,k)
			            f2(nx-1,j,k)=f2(1,j,k)
			            f10(nx-1,j,k)=f10(1,j,k)
			            f16(nx-1,j,k)=f16(1,j,k)
			          elseif(k==nz-1)then
			            f1(2,j,k)=f1(nx,j,k)
			            f7(2,j,k)=f7(nx,j,k)
			            f15(2,j,k)=f15(nx,j,k)
			            f2(nx-1,j,k)=f2(1,j,k)
			            f10(nx-1,j,k)=f10(1,j,k)
			            f17(nx-1,j,k)=f17(1,j,k)
			          else
			            f1(2,j,k)=f1(nx,j,k)
			            f7(2,j,k)=f7(nx,j,k)
			            f15(2,j,k)=f15(nx,j,k)
			            f18(2,j,k)=f18(nx,j,k)
			            f2(nx-1,j,k)=f2(1,j,k)
			            f10(nx-1,j,k)=f10(1,j,k)
			            f16(nx-1,j,k)=f16(1,j,k)
			            f17(nx-1,j,k)=f17(1,j,k)
			          endif
			        endif
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
			      if(i>2 .and. i<nx-1 .and. k>2 .and. k<nz-1)then
			        f3(i,2,k)=f3(i,ny,k)
			        f7(i,2,k)=f7(i,ny,k)
			        f10(i,2,k)=f10(i,ny,k)
			        f11(i,2,k)=f11(i,ny,k)
			        f13(i,2,k)=f13(i,ny,k)
			        f4(i,ny-1,k)=f4(i,1,k)
			        f8(i,ny-1,k)=f8(i,1,k)
			        f9(i,ny-1,k)=f9(i,1,k)
			        f12(i,ny-1,k)=f12(i,1,k)
			        f14(i,ny-1,k)=f14(i,1,k)
			      else
			        if(i==2)then
			          if(k==2)then
			            f3(i,2,k)=f3(i,ny,k)
			            f10(i,2,k)=f10(i,ny,k)
			            f13(i,2,k)=f13(i,ny,k)
			            f4(i,ny-1,k)=f4(i,1,k)
			            f8(i,ny-1,k)=f8(i,1,k)
			            f12(i,ny-1,k)=f12(i,1,k)
			          elseif(k==nz-1)then
			            f3(i,2,k)=f3(i,ny,k)
			            f10(i,2,k)=f10(i,ny,k)
			            f11(i,2,k)=f11(i,ny,k)
			            f4(i,ny-1,k)=f4(i,1,k)
			            f8(i,ny-1,k)=f8(i,1,k)
			            f14(i,ny-1,k)=f14(i,1,k)
			          else
			            f3(i,2,k)=f3(i,ny,k)
			            f10(i,2,k)=f10(i,ny,k)
			            f11(i,2,k)=f11(i,ny,k)
			            f13(i,2,k)=f13(i,ny,k)
			            f4(i,ny-1,k)=f4(i,1,k)
			            f8(i,ny-1,k)=f8(i,1,k)
			            f12(i,ny-1,k)=f12(i,1,k)
			            f14(i,ny-1,k)=f14(i,1,k)
			          endif
			        elseif(i==nx-1)then
			          if(k==2)then
			            f3(i,2,k)=f3(i,ny,k)
			            f7(i,2,k)=f7(i,ny,k)
			            f13(i,2,k)=f13(i,ny,k)
			            f4(i,ny-1,k)=f4(i,1,k)
			            f9(i,ny-1,k)=f9(i,1,k)
			            f12(i,ny-1,k)=f12(i,1,k)
			          elseif(k==nz-1)then 
			            f3(i,2,k)=f3(i,ny,k)
			            f7(i,2,k)=f7(i,ny,k)
			            f11(i,2,k)=f11(i,ny,k)
			            f4(i,ny-1,k)=f4(i,1,k)
			            f9(i,ny-1,k)=f9(i,1,k)
			            f14(i,ny-1,k)=f14(i,1,k)
			          else
			            f3(i,2,k)=f3(i,ny,k)
			            f7(i,2,k)=f7(i,ny,k)
			            f11(i,2,k)=f11(i,ny,k)
			            f13(i,2,k)=f13(i,ny,k)
			            f4(i,ny-1,k)=f4(i,1,k)
			            f9(i,ny-1,k)=f9(i,1,k)
			            f12(i,ny-1,k)=f12(i,1,k)
			            f14(i,ny-1,k)=f14(i,1,k)
			          endif
			        endif
			      endif
			    enddo
			enddo
			!$acc end kernels
			
          endif
!             f1(2,:,:)=f1(nx,:,:)
!             f7(2,:,:)=f7(nx,:,:)
!             f9(2,:,:)=f9(nx,:,:)
!             f15(2,:,:)=f15(nx,:,:)
!             f18(2,:,:)=f18(nx,:,:)
!             !x=nx 
!             f2(nx-1,:,:)=f2(1,:,:)
!             f8(nx-1,:,:)=f8(1,:,:)
!             f10(nx-1,:,:)=f10(1,:,:)
!             f16(nx-1,:,:)=f16(1,:,:)
!             f17(nx-1,:,:)=f17(1,:,:)

!             !y=1
!             f3(:,2,:)=f3(:,ny,:)
!             f7(:,2,:)=f7(:,ny,:)
!             f10(:,2,:)=f10(:,ny,:)
!             f11(:,2,:)=f11(:,ny,:)
!             f13(:,2,:)=f13(:,ny,:)
        
!             !y=ny
!             f4(:,ny-1,:)=f4(:,1,:)
!             f8(:,ny-1,:)=f8(:,1,:)
!             f9(:,ny-1,:)=f9(:,1,:)
!             f12(:,ny-1,:)=f12(:,1,:)
!             f14(:,ny-1,:)=f14(:,1,:)
       
        
    enddo 
    !$acc wait
    if(lasync)then
      if(lvtk)then
        call print_vtk_sync(iframe)
      else
        call print_raw_sync(iframe)
      endif
    endif
    !!$acc update host(rho,u,v,w)
    
    !$acc end data
    call cpu_time(ts2)
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
