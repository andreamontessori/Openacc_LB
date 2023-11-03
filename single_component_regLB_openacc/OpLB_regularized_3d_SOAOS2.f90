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
    
    integer :: itile,jtile,ktile
    
    logical :: lprint,lvtk,lasync,lpbc
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
    real(kind=db) :: visc_LB,uu,udotc,omega,feq
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp,init_rho

    real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq17
    real(kind=db) :: fneq9,fneq10,fneq11,fneq12,fneq13,fneq14,fneq15,fneq16,fneq18
    real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
    real(kind=db) :: pi2cssq1,pi2cssq2,pi2cssq0
    
    !integer :: TILE_DIMx,TILE_DIMy,TILE_DIMz
    integer :: nxblock,nyblock,nzblock,nblocks,nxyblock
    integer :: idblock,xblock,yblock,zblock
    integer :: ii,jj,kk
    integer :: oidblock,oxblock,oyblock,ozblock
    integer :: oii,ojj,okk
    integer :: oi,oj,ok
    
    integer(kind=4), allocatable,dimension(:,:,:,:) :: provace
    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    !real(kind=db), allocatable, dimension(:,:,:) :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    !real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
    !real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18
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

#define TILE_DIMx 4
#define TILE_DIMy 4
#define TILE_DIMz 4    

        if (mod(nx, TILE_DIMx)/= 0) then
          write(*,*) 'nx must be a multiple of TILEDIMx'
          stop
        end if
        if (mod(ny, TILE_DIMy) /= 0) then
          write(*,*) 'ny must be a multiple of TILEDIMy'
          stop
        end if
        if (mod(nz, TILE_DIMz) /= 0) then
          write(*,*) 'nz must be a multiple of TILEDIMz'
          stop
        end if
        
        nxblock=nx/TILE_DIMx
        nyblock=ny/TILE_DIMx
        nzblock=nz/TILE_DIMz
        
        nxyblock=nxblock*nyblock
        
        nblocks=nxblock*nyblock*nzblock
        
        write(6,*)'nx,ny,nz',nx,ny,nz
        write(6,*)'TILEDIMx,TILEDIMy,TILEDIMz',TILE_DIMx,TILE_DIMy,TILE_DIMz
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
        !isfluid(1,:,:)=0 !left
        !isfluid(nx,:,:)=0 !right
        !isfluid(:,1,:)=0 !front 
        !isfluid(:,ny,:)=0 !rear
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
        
!        do zblock=1,nzblock
!          do yblock=1,nyblock
!            do xblock=1,nxblock
!              idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,1,idblock)=0.0_db  !pxx
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,2,idblock)=0.0_db  !pyy
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,3,idblock)=0.0_db  !pzz
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,4,idblock)=0.0_db  !pxy
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,5,idblock)=0.0_db  !pxz
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,6,idblock)=0.0_db  !pyz
              
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,7,idblock)=0.0_db  !u
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,8,idblock)=0.0_db  !v
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,9,idblock)=0.0_db  !w
              
!              hfields(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,10,idblock)=init_rho  !rho
              
              
!              fpops(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,0,idblock)=init_rho*p0
!              fpops(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,1:6,idblock)=init_rho*p1
!              fpops(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz,7:18,idblock)=init_rho*p2
              
!            enddo
!          enddo
!        enddo
        
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
              if(i==nx/4 .and. j==ny/4 .and. k==nz/4)then
                init_rho=1.1_db
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
              fpops(ii,jj,kk,1,idblock)=init_rho*p1
              fpops(ii,jj,kk,2,idblock)=init_rho*p1
              fpops(ii,jj,kk,3,idblock)=init_rho*p1
              fpops(ii,jj,kk,4,idblock)=init_rho*p1
              fpops(ii,jj,kk,5,idblock)=init_rho*p1
              fpops(ii,jj,kk,6,idblock)=init_rho*p1
              fpops(ii,jj,kk,7,idblock)=init_rho*p2
              fpops(ii,jj,kk,8,idblock)=init_rho*p2
              fpops(ii,jj,kk,9,idblock)=init_rho*p2
              fpops(ii,jj,kk,10,idblock)=init_rho*p2
              fpops(ii,jj,kk,11,idblock)=init_rho*p2
              fpops(ii,jj,kk,12,idblock)=init_rho*p2
              fpops(ii,jj,kk,13,idblock)=init_rho*p2
              fpops(ii,jj,kk,14,idblock)=init_rho*p2
              fpops(ii,jj,kk,15,idblock)=init_rho*p2
              fpops(ii,jj,kk,16,idblock)=init_rho*p2
              fpops(ii,jj,kk,17,idblock)=init_rho*p2
              fpops(ii,jj,kk,18,idblock)=init_rho*p2
              
              !rhoprint(i,j,k)=init_rho
              !velprint(1,i,j,k)=0.0_db
              !velprint(2,i,j,k)=0.0_db
              !velprint(3,i,j,k)=0.0_db
              
              !write(6,'(10i8)')i,j,k,ii,jj,kk,xblock,yblock,zblock,idblock
              
            enddo
          enddo
        enddo
        
        allocate(provace(3,nx,ny,nz))
        
        !$acc loop collapse(6)
        !tile(8,16)
        do ktile = 1,nz,TILE_DIMz
          do jtile = 1,ny,TILE_DIMy
            do itile = 1,nx,TILE_DIMx
              do k = ktile,ktile+TILE_DIMz-1
                do j = jtile,jtile+TILE_DIMy-1
                  do i = itile,itile+TILE_DIMx-1
                   ! write(6,'(10i8)')i,j,k,itile,jtile,ktile
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
                   provace(1,i,j,k)=itile
                   provace(2,i,j,k)=jtile
                   provace(3,i,j,k)=ktile
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        
        
        !$acc loop tile(TILE_DIMx,TILE_DIMy,TILE_DIMz) 
        do k=1,nz
          do j=1,ny
            do i=1,nx
              itile=((i+TILE_DIMx-1)/TILE_DIMx-1)*TILE_DIMx+1
              jtile=((j+TILE_DIMy-1)/TILE_DIMy-1)*TILE_DIMy+1
              ktile=((k+TILE_DIMz-1)/TILE_DIMz-1)*TILE_DIMz+1
              
              !write(6,'(10i8)')i,j,k,itile,jtile,ktile
              if(provace(1,i,j,k).ne.itile .or. provace(2,i,j,k).ne.jtile .or. provace(3,i,j,k).ne.ktile)then
                write(6,'(a,10i8)')'CAZZO',i,j,k,itile,jtile,ktile
                stop
              endif
            enddo
          enddo
        enddo
        
!        iframe=0
!        call init_output(nx,ny,nz,1,lvtk)
!      call string_char(head1,nheadervtk(1),headervtk(1))
!      call string_char(head2,nheadervtk(2),headervtk(2))
!        call print_vtk_sync(iframe)
!      stop 
              
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
        write(6,*) 'TILEDIMx',TILE_DIMx
        write(6,*) 'TILEDIMy',TILE_DIMy
        write(6,*) 'TILEDIMz',TILE_DIMz
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
    !$acc& nzblock,nxyblock,nblocks) async(1)
    
    
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
      !$acc kernels present(hfields,rhoprint) async(1)
      !$acc loop independent tile(TILE_DIMx,TILE_DIMy,TILE_DIMz) private(i,j,k,&
      !$acc& xblock,yblock,zblock,idblock,ii,jj,kk)
      !tile(TILE_DIMx,TILE_DIMy,TILE_DIMz) 
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
            !write(*,*)i,j,k,rhoprint(i,j,k)
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
        !call sleep(1)
!        do k=1,nz
!         do j=1,ny
!          do i=1,nx
!            write(6,'(a,3i8,f16.8)')'CIAONE',i,j,k,rhoprint(i,j,k)
!          enddo
!        enddo
!      enddo
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
        !***********************************moments collision bbck + forcing************************ 
        !$acc kernels present(fpops,hfields,isfluid) async(1)
        !$acc loop tile(TILE_DIMx,TILE_DIMy,TILE_DIMz) private(fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq9,fneq10,fneq11,&
        !$acc& fneq12,fneq3,fneq14,fneq15,uu,temp,udotc,xblock,yblock,zblock,idblock,ii,jj,kk)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
                        xblock=(i+TILE_DIMx-1)/TILE_DIMx
                        yblock=(j+TILE_DIMy-1)/TILE_DIMy
                        zblock=(k+TILE_DIMz-1)/TILE_DIMz
                        idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                        ii=i-xblock*TILE_DIMx+TILE_DIMx
                        jj=j-yblock*TILE_DIMy+TILE_DIMy
                        kk=k-zblock*TILE_DIMz+TILE_DIMz
                        hfields(ii,jj,kk,10,idblock) = fpops(ii,jj,kk,0,idblock)+fpops(ii,jj,kk,1,idblock)+fpops(ii,jj,kk,2,idblock)+fpops(ii,jj,kk,3,idblock)+fpops(ii,jj,kk,4,idblock)+fpops(ii,jj,kk,5,idblock) &
                            +fpops(ii,jj,kk,6,idblock)+fpops(ii,jj,kk,7,idblock)+fpops(ii,jj,kk,8,idblock)+fpops(ii,jj,kk,9,idblock)+fpops(ii,jj,kk,10,idblock)+fpops(ii,jj,kk,11,idblock) &
                            +fpops(ii,jj,kk,12,idblock)+fpops(ii,jj,kk,13,idblock)+fpops(ii,jj,kk,14,idblock)+fpops(ii,jj,kk,15,idblock)+fpops(ii,jj,kk,16,idblock)+fpops(ii,jj,kk,17,idblock) &
                            +fpops(ii,jj,kk,18,idblock)

                        hfields(ii,jj,kk,7,idblock) = (fpops(ii,jj,kk,1,idblock)+fpops(ii,jj,kk,7,idblock)+fpops(ii,jj,kk,9,idblock)+fpops(ii,jj,kk,15,idblock)+fpops(ii,jj,kk,18,idblock)) &
                             -(fpops(ii,jj,kk,2,idblock)+fpops(ii,jj,kk,8,idblock)+fpops(ii,jj,kk,10,idblock)+fpops(ii,jj,kk,16,idblock)+fpops(ii,jj,kk,17,idblock)) 
                        
                        hfields(ii,jj,kk,8,idblock) = (fpops(ii,jj,kk,3,idblock)+fpops(ii,jj,kk,7,idblock)+fpops(ii,jj,kk,10,idblock)+fpops(ii,jj,kk,11,idblock)+fpops(ii,jj,kk,13,idblock)) &
                            -(fpops(ii,jj,kk,4,idblock)+fpops(ii,jj,kk,8,idblock)+fpops(ii,jj,kk,9,idblock)+fpops(ii,jj,kk,12,idblock)+fpops(ii,jj,kk,14,idblock))

                        hfields(ii,jj,kk,9,idblock) = (fpops(ii,jj,kk,5,idblock)+fpops(ii,jj,kk,11,idblock)+fpops(ii,jj,kk,14,idblock)+fpops(ii,jj,kk,15,idblock)+fpops(ii,jj,kk,17,idblock)) &
                            -(fpops(ii,jj,kk,6,idblock)+fpops(ii,jj,kk,12,idblock)+fpops(ii,jj,kk,13,idblock)+fpops(ii,jj,kk,16,idblock)+fpops(ii,jj,kk,18,idblock))
                        
                        uu=0.5_db*(hfields(ii,jj,kk,7,idblock)*hfields(ii,jj,kk,7,idblock) + &
                         hfields(ii,jj,kk,8,idblock)*hfields(ii,jj,kk,8,idblock) + &
                         hfields(ii,jj,kk,9,idblock)*hfields(ii,jj,kk,9,idblock))/cssq
                        !1-2
                        udotc=hfields(ii,jj,kk,7,idblock)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=fpops(ii,jj,kk,1,idblock)-p1*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq2=fpops(ii,jj,kk,2,idblock)-p1*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !3-4
                        udotc=hfields(ii,jj,kk,8,idblock)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq3=fpops(ii,jj,kk,3,idblock)-p1*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq4=fpops(ii,jj,kk,4,idblock)-p1*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !5-6
                        udotc=hfields(ii,jj,kk,9,idblock)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq5=fpops(ii,jj,kk,5,idblock)-p1*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq6=fpops(ii,jj,kk,6,idblock)-p1*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !7-8
                        udotc=(hfields(ii,jj,kk,7,idblock)+hfields(ii,jj,kk,8,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq7=fpops(ii,jj,kk,7,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq8=fpops(ii,jj,kk,8,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !10-9
                        udotc=(-hfields(ii,jj,kk,7,idblock)+hfields(ii,jj,kk,8,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq10=fpops(ii,jj,kk,10,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq9=fpops(ii,jj,kk,9,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !11-12
                        udotc=(hfields(ii,jj,kk,8,idblock)+hfields(ii,jj,kk,9,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq11=fpops(ii,jj,kk,11,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq12=fpops(ii,jj,kk,12,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !13-14
                        udotc=(hfields(ii,jj,kk,8,idblock)-hfields(ii,jj,kk,9,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq13=fpops(ii,jj,kk,13,idblock) - p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq14=fpops(ii,jj,kk,14,idblock) - p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !15-16
                        udotc=(hfields(ii,jj,kk,7,idblock)+hfields(ii,jj,kk,9,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq15=fpops(ii,jj,kk,15,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq16=fpops(ii,jj,kk,16,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !17-18
                        udotc=(-hfields(ii,jj,kk,7,idblock)+hfields(ii,jj,kk,9,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq17=fpops(ii,jj,kk,17,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq18=fpops(ii,jj,kk,18,idblock)-p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        hfields(ii,jj,kk,1,idblock)=fneq1+fneq2+fneq7+fneq8+fneq9+fneq10+fneq15+fneq16+fneq17+fneq18
                        hfields(ii,jj,kk,2,idblock)=fneq3+fneq4+fneq7+fneq8+fneq9+fneq10+fneq11+fneq12+fneq13+fneq14
                        hfields(ii,jj,kk,3,idblock)=fneq5+fneq6+fneq11+fneq12+fneq13+fneq14+fneq15+fneq16+fneq17+fneq18
                        hfields(ii,jj,kk,4,idblock)= fneq7+fneq8-fneq9-fneq10
                        hfields(ii,jj,kk,5,idblock)=fneq15+fneq16-fneq17-fneq18
                        hfields(ii,jj,kk,6,idblock)=fneq11+fneq12-fneq13-fneq14
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        !$acc wait(1)
        
        !***********************************PRINT************************
        if(mod(step,stamp).eq.0)write(6,'(a,i8)')'step : ',step
        if(lprint)then
          if(mod(step,stamp).eq.0)then
            iframe=iframe+1
            !$acc wait(1)
            !$acc kernels present(rhoprint,velprint,hfields) async(1)
            !$acc loop independent tile(TILE_DIMx,TILE_DIMy,TILE_DIMz) private(i,j,k,xblock,yblock,zblock,idblock,ii,jj,kk)
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
        !$acc kernels present(fpops,hfields,isfluid) async(1)
        !$acc loop independent tile(TILE_DIMx,TILE_DIMy,TILE_DIMz) private(feq,uu,temp,udotc,&
        !$acc& xblock,yblock,zblock,idblock,ii,jj,kk,oxblock,oyblock,ozblock,oidblock,oii,ojj,okk,&
        !$acc& oi,oj,ok)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
                        xblock=(i+TILE_DIMx-1)/TILE_DIMx
                        yblock=(j+TILE_DIMy-1)/TILE_DIMy
                        zblock=(k+TILE_DIMz-1)/TILE_DIMz
                        idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                        ii=i-xblock*TILE_DIMx+TILE_DIMx
                        jj=j-yblock*TILE_DIMy+TILE_DIMy
                        kk=k-zblock*TILE_DIMz+TILE_DIMz
                        uu=0.5_db*(hfields(ii,jj,kk,7,idblock)*hfields(ii,jj,kk,7,idblock) + &
                         hfields(ii,jj,kk,8,idblock)*hfields(ii,jj,kk,8,idblock) + &
                         hfields(ii,jj,kk,9,idblock)*hfields(ii,jj,kk,9,idblock))/cssq
                        !0
                        feq=p0*(hfields(ii,jj,kk,10,idblock)-uu)
                        fpops(ii,jj,kk,0,idblock)=feq + (1.0_db-omega)*pi2cssq0*&
                         (-cssq*(hfields(ii,jj,kk,2,idblock)+hfields(ii,jj,kk,1,idblock)&
                         +hfields(ii,jj,kk,3,idblock)))
                        
                        !1
                        oi=i+1
                        oi=mod(oi+nx-1,nx)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        oidblock=(oxblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        udotc=hfields(ii,jj,kk,7,idblock)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq1=(1.0_db-omega)*pi2cssq1*(qxx*hfields(ii,jj,kk,1,idblock)&
                         -cssq*(hfields(ii,jj,kk,2,idblock)+hfields(ii,jj,kk,3,idblock)))
                        !f1(i+1,j,k)
                        fpops(oii,jj,kk,1,oidblock)=feq + fneq1 + fx*p1dcssq
                        
                        !2
                        oi=i-1
                        oi=mod(oi+nx-1,nx)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        oidblock=(oxblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        feq=p1*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f2(i-1,j,k)
                        fpops(oii,jj,kk,2,oidblock)=feq + fneq1 - fx*p1dcssq
                        
                        !3
                        oj=j+1
                        oj=mod(oj+ny-1,ny)+1
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        oidblock=(xblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        udotc=hfields(ii,jj,kk,8,idblock)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq3=(1.0_db-omega)*pi2cssq1*(qyy*hfields(ii,jj,kk,2,idblock)&
                         -cssq*(hfields(ii,jj,kk,1,idblock)+hfields(ii,jj,kk,3,idblock)))
                        !f3(i,j+1,k)
                        fpops(ii,ojj,kk,3,oidblock)=feq+fneq3 + fy*p1dcssq
                        
                        !4
                        oj=j-1
                        oj=mod(oj+ny-1,ny)+1
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        oidblock=(xblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        feq=p1*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f4(i,j-1,k)
                        fpops(ii,ojj,kk,4,oidblock)=feq+fneq3 - fy*p1dcssq
                        
                        !7
                        oi=i+1
                        oj=j+1
                        oi=mod(oi+nx-1,nx)+1
                        oj=mod(oj+ny-1,ny)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        oidblock=(oxblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        udotc=(hfields(ii,jj,kk,7,idblock)+hfields(ii,jj,kk,8,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq7=(1.0_db-omega)*pi2cssq2*(qxx*hfields(ii,jj,kk,1,idblock)&
                         +qyy*hfields(ii,jj,kk,2,idblock)-cssq*hfields(ii,jj,kk,3,idblock)&
                         +2.0_db*qxy_7_8*hfields(ii,jj,kk,4,idblock))
                        !f7(i+1,j+1,k)
                        fpops(oii,ojj,kk,7,oidblock)=feq + fneq7 + (fx+fy)*p2dcssq 
                        
                        !8
                        oi=i-1
                        oj=j-1
                        oi=mod(oi+nx-1,nx)+1
                        oj=mod(oj+ny-1,ny)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        oidblock=(oxblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f8(i-1,j-1,k)
                        fpops(oii,ojj,kk,8,oidblock)=feq + fneq7 - (fx+fy)*p2dcssq
                        
                        !10
                        oi=i-1
                        oj=j+1
                        oi=mod(oi+nx-1,nx)+1
                        oj=mod(oj+ny-1,ny)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        oidblock=(oxblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        udotc=(-hfields(ii,jj,kk,7,idblock)+hfields(ii,jj,kk,8,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq10=(1.0_db-omega)*pi2cssq2*(qxx*hfields(ii,jj,kk,1,idblock)&
                         +qyy*hfields(ii,jj,kk,2,idblock)-cssq*hfields(ii,jj,kk,3,idblock)&
                         +2.0_db*qxy_9_10*hfields(ii,jj,kk,4,idblock))
                        !f10(i-1,j+1,k)
                        fpops(oii,ojj,kk,10,oidblock)=feq+fneq10 +(fy-fx)*p2dcssq
                        
                        !9
                        oi=i+1
                        oj=j-1
                        oi=mod(oi+nx-1,nx)+1
                        oj=mod(oj+ny-1,ny)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        oidblock=(oxblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f9(i+1,j-1,k)
                        fpops(oii,ojj,kk,9,oidblock)=feq+fneq10 + (fx-fy)*p2dcssq

                        !5
                        ok=k+1
                        ok=mod(ok+nz-1,nz)+1
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(xblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        udotc=hfields(ii,jj,kk,9,idblock)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq5=(1.0_db-omega)*pi2cssq1*(qzz*hfields(ii,jj,kk,3,idblock)&
                         -cssq*(hfields(ii,jj,kk,1,idblock)+hfields(ii,jj,kk,2,idblock)))
                        !f5(i,j,k+1)
                        fpops(ii,jj,okk,5,oidblock)=feq+fneq5 + fz*p1dcssq
                        
                        !6
                        ok=k-1
                        ok=mod(ok+nz-1,nz)+1
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(xblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        feq=p1*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f6(i,j,k-1)
                        fpops(ii,jj,okk,6,oidblock)=feq+fneq5 - fz*p1dcssq

                        !15
                        oi=i+1
                        ok=k+1
                        oi=mod(oi+nx-1,nx)+1
                        ok=mod(ok+nz-1,nz)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(oxblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        udotc=(hfields(ii,jj,kk,7,idblock)+hfields(ii,jj,kk,9,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq15=(1.0_db-omega)*pi2cssq2*(qxx*hfields(ii,jj,kk,1,idblock)&
                         +qzz*hfields(ii,jj,kk,3,idblock)-cssq*hfields(ii,jj,kk,2,idblock)&
                         +2.0_db*qxz_15_16*hfields(ii,jj,kk,5,idblock))
                        !f15(i+1,j,k+1)
                        fpops(oii,jj,okk,15,oidblock)=feq+fneq15 + (fx+fz)*p2dcssq 
                        
                        !16
                        oi=i-1
                        ok=k-1
                        oi=mod(oi+nx-1,nx)+1
                        ok=mod(ok+nz-1,nz)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(oxblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f16(i-1,j,k-1)
                        fpops(oii,jj,okk,16,oidblock)=feq+fneq15 - (fx+fz)*p2dcssq

                        !17
                        oi=i-1
                        ok=k+1
                        oi=mod(oi+nx-1,nx)+1
                        ok=mod(ok+nz-1,nz)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(oxblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        udotc=(-hfields(ii,jj,kk,7,idblock)+hfields(ii,jj,kk,9,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq17=(1.0_db-omega)*pi2cssq2*(qxx*hfields(ii,jj,kk,1,idblock)&
                         +qzz*hfields(ii,jj,kk,3,idblock)-cssq*hfields(ii,jj,kk,2,idblock)&
                         +2.0_db*qxz_17_18*hfields(ii,jj,kk,5,idblock))
                        !f17(i-1,j,k+1)
                        fpops(oii,jj,okk,17,oidblock)=feq+fneq17 +(fz-fx)*p2dcssq
                        
                        !18
                        oi=i+1
                        ok=k-1
                        oi=mod(oi+nx-1,nx)+1
                        ok=mod(ok+nz-1,nz)+1
                        oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(oxblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f18(i+1,j,k-1)
                        fpops(oii,jj,okk,18,oidblock)=feq+fneq17 + (fx-fz)*p2dcssq

                        !11
                        oj=j+1
                        ok=k+1
                        oj=mod(oj+ny-1,ny)+1
                        ok=mod(ok+nz-1,nz)+1
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(xblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        udotc=(hfields(ii,jj,kk,8,idblock)+hfields(ii,jj,kk,9,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq11=(1.0_db-omega)*pi2cssq2*(qyy*hfields(ii,jj,kk,2,idblock)&
                         +qzz*hfields(ii,jj,kk,3,idblock)-cssq*hfields(ii,jj,kk,1,idblock)&
                         +2.0_db*qyz_11_12*hfields(ii,jj,kk,6,idblock))
                        !f11(i,j+1,k+1)
                        fpops(ii,ojj,okk,11,oidblock)=feq+fneq11+(fy+fz)*p2dcssq
                        
                        !12
                        oj=j-1
                        ok=k-1
                        oj=mod(oj+ny-1,ny)+1
                        ok=mod(ok+nz-1,nz)+1
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(xblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f12(i,j-1,k-1)
                        fpops(ii,ojj,okk,12,oidblock)=feq+fneq11 - (fy+fz)*p2dcssq

                        !13
                        oj=j+1
                        ok=k-1
                        oj=mod(oj+ny-1,ny)+1
                        ok=mod(ok+nz-1,nz)+1
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(xblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        udotc=(hfields(ii,jj,kk,8,idblock)-hfields(ii,jj,kk,9,idblock))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp + udotc))
                        fneq13=(1.0_db-omega)*pi2cssq2*(qyy*hfields(ii,jj,kk,2,idblock)&
                         +qzz*hfields(ii,jj,kk,3,idblock)-cssq*hfields(ii,jj,kk,1,idblock)&
                         +2.0_db*qyz_13_14*hfields(ii,jj,kk,6,idblock))
                        !f13(i,j+1,k-1)
                        fpops(ii,ojj,okk,13,oidblock)=feq+fneq13 + (fy-fz)*p2dcssq
                        
                        !14
                        oj=j-1
                        ok=k+1
                        oj=mod(oj+ny-1,ny)+1
                        ok=mod(ok+nz-1,nz)+1
                        oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                        ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                        oidblock=(xblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                        ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                        okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                        feq=p2*(hfields(ii,jj,kk,10,idblock)+(temp - udotc))
                        !f14(i,j-1,k+1)
                        fpops(ii,ojj,okk,14,oidblock)=feq+fneq13 + (fz-fy)*p2dcssq
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        !$acc wait(1)
        !********************************boundary conditions no slip everywhere********************************!
        !$acc kernels present(fpops,hfields,isfluid) async(1)
        !$acc loop independent tile(TILE_DIMx,TILE_DIMy,TILE_DIMz) private(xblock,yblock,zblock,&
        !$acc& idblock,ii,jj,kk,oxblock,oyblock,ozblock,oidblock,oii,ojj,okk,&
        !$acc& oi,oj,ok)
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        if(isfluid(i,j,k).eq.0)then
                            xblock=(i+TILE_DIMx-1)/TILE_DIMx
                            yblock=(j+TILE_DIMy-1)/TILE_DIMy
                            zblock=(k+TILE_DIMz-1)/TILE_DIMz
                            idblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                            ii=i-xblock*TILE_DIMx+TILE_DIMx
                            jj=j-yblock*TILE_DIMy+TILE_DIMy
                            kk=k-zblock*TILE_DIMz+TILE_DIMz
                            
                            oi=i+1
                            ok=k-1
                            oi=mod(oi+nx-1,nx)+1
                            ok=mod(ok+nz-1,nz)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(oxblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(oii,jj,okk,18,oidblock)=fpops(ii,jj,kk,17,idblock) !gpc
                            
                            oi=i-1
                            ok=k+1
                            oi=mod(oi+nx-1,nx)+1
                            ok=mod(ok+nz-1,nz)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(oxblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz 
                            fpops(oii,jj,okk,17,oidblock)=fpops(ii,jj,kk,18,idblock) !hpc
                            
                            oi=i-1
                            ok=k-1
                            oi=mod(oi+nx-1,nx)+1
                            ok=mod(ok+nz-1,nz)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(oxblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(oii,jj,okk,16,oidblock)=fpops(ii,jj,kk,15,idblock) !gpc
                            
                            oi=i+1
                            ok=k+1
                            oi=mod(oi+nx-1,nx)+1
                            ok=mod(ok+nz-1,nz)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(oxblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(oii,jj,okk,15,oidblock)=fpops(ii,jj,kk,16,idblock) !hpc
                            
                            oj=j-1
                            ok=k+1
                            oj=mod(oj+ny-1,ny)+1
                            ok=mod(ok+nz-1,nz)+1
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(xblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(ii,ojj,okk,14,oidblock)=fpops(ii,jj,kk,13,idblock)!gpc 
                            
                            oj=j+1
                            ok=k-1
                            oj=mod(oj+ny-1,ny)+1
                            ok=mod(ok+nz-1,nz)+1
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(xblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(ii,ojj,okk,13,oidblock)=fpops(ii,jj,kk,14,idblock)!hpc
                            
                            oj=j-1
                            ok=k-1
                            oj=mod(oj+ny-1,ny)+1
                            ok=mod(ok+nz-1,nz)+1
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(xblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(ii,ojj,okk,12,oidblock)=fpops(ii,jj,kk,11,idblock)!gpc 
                            
                            oj=j+1
                            ok=k+1
                            oj=mod(oj+ny-1,ny)+1
                            ok=mod(ok+nz-1,nz)+1
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(xblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(ii,ojj,okk,11,oidblock)=fpops(ii,jj,kk,12,idblock)!hpc
                            
                            oi=i-1
                            oj=j+1
                            oi=mod(oi+nx-1,nx)+1
                            oj=mod(oj+ny-1,ny)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            oidblock=(oxblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            fpops(oii,ojj,kk,10,oidblock)=fpops(ii,jj,kk,9,idblock)!gpc 
                            
                            oi=i+1
                            oj=j-1
                            oi=mod(oi+nx-1,nx)+1
                            oj=mod(oj+ny-1,ny)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            oidblock=(oxblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            fpops(oii,ojj,kk,9,oidblock)=fpops(ii,jj,kk,10,idblock)!hpc
                            
                            oi=i-1
                            oj=j-1
                            oi=mod(oi+nx-1,nx)+1
                            oj=mod(oj+ny-1,ny)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            oidblock=(oxblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            fpops(oii,ojj,kk,8,oidblock)=fpops(ii,jj,kk,7,idblock)!gpc 
                            
                            oi=i+1
                            oj=j+1
                            oi=mod(oi+nx-1,nx)+1
                            oj=mod(oj+ny-1,ny)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            oidblock=(oxblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            fpops(oii,ojj,kk,7,oidblock)=fpops(ii,jj,kk,8,idblock)!hpc
                            
                            ok=k-1
                            ok=mod(ok+nz-1,nz)+1
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(xblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(ii,jj,okk,6,oidblock)=fpops(ii,jj,kk,5,idblock)!gpc 
                            
                            ok=k+1
                            ok=mod(ok+nz-1,nz)+1
                            ozblock=(ok+TILE_DIMz-1)/TILE_DIMz
                            oidblock=(xblock-1)+(yblock-1)*nxblock+(ozblock-1)*nxyblock+1
                            okk=ok-ozblock*TILE_DIMz+TILE_DIMz
                            fpops(ii,jj,okk,5,oidblock)=fpops(ii,jj,kk,6,idblock)!hpc 

                            oj=j-1
                            oj=mod(oj+ny-1,ny)+1
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            oidblock=(xblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            fpops(ii,ojj,kk,4,oidblock)=fpops(ii,jj,kk,3,idblock)!gpc 
                            
                            oj=j+1
                            oj=mod(oj+ny-1,ny)+1
                            oyblock=(oj+TILE_DIMy-1)/TILE_DIMy
                            oidblock=(xblock-1)+(oyblock-1)*nxblock+(zblock-1)*nxyblock+1
                            ojj=oj-oyblock*TILE_DIMy+TILE_DIMy
                            fpops(ii,ojj,kk,3,oidblock)=fpops(ii,jj,kk,4,idblock)!hpc 
                            
                            oi=i-1
                            oi=mod(oi+nx-1,nx)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            oidblock=(oxblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            fpops(oii,jj,kk,2,oidblock)=fpops(ii,jj,kk,1,idblock)!gpc 
                            
                            oi=i+1
                            oi=mod(oi+nx-1,nx)+1
                            oxblock=(oi+TILE_DIMx-1)/TILE_DIMx
                            oidblock=(oxblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                            oii=oi-oxblock*TILE_DIMx+TILE_DIMx
                            fpops(oii,jj,kk,1,oidblock)=fpops(ii,jj,kk,2,idblock)!hpc 
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
    !!$acc update host(rho,u,v,w)
    
    !$acc end data
    call cpu_time(ts2)
    write(6,*) 'u=',velprint(1,nx/2,ny/2,nz/2),'v=',velprint(2,nx/2,ny/2,nz/2),&
     'w=',velprint(3,nx/2,ny/2,nz/2),'rho=',rhoprint(nx/2,ny/2,nz/2)
    write(6,*) 'u=',velprint(1,nx/2,ny/2,1),'v=',velprint(2,nx/2,ny/2,1),&
     'w=',velprint(3,nx/2,ny/2,1),'rho=',rhoprint(nx/2,ny/2,1)
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
   integer :: i,j,k
   
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
