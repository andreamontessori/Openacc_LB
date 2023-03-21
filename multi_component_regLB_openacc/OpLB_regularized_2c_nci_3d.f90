 module prints
  
  implicit none
  
    integer, parameter, private :: db=4 !kind(1.0)
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
  
  allocate(rhoprint(1:nx,1:ny,1:nz))
  allocate(velprint(1:3,1:nx,1:ny,1:nz))
  rhoprint(1:nx,1:ny,1:nz)=0.0
  velprint(1:3,1:nx,1:ny,1:nz)=0.0
  
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
    !************************************block of vars*******************************************!
        integer, parameter :: db=4 !kind(1.0)
        integer :: i,j,k,k_init,k_end,j_init,j_end,i_init,i_end
        integer :: nx,ny,nz,step,stamp,nlinks,nsteps,ngpus
        integer,save :: iframe=0
        logical :: lprint,lvtk,lasync
        
        real(kind=db),parameter :: pi_greek=3.14159265359793234626433
        
        real(kind=4)  :: ts1,ts2,p0,p1,p2,p1dcssq,p2dcssq
        real(kind=4)  :: p1cg,p2cg,p3cg
        real(kind=db) :: visc_LB,uu,udotc,omega,feq,geq,fpc,gpc,hpc
        real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp
        real(kind=db) :: radius

        real(kind=db) :: fneq1,fneq2,fneq3,fneq4,fneq5,fneq6,fneq7,fneq8,fneq17
        real(kind=db) :: fneq9,fneq10,fneq11,fneq12,fneq13,fneq14,fneq15,fneq16,fneq18
        real(kind=db) :: ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft17
        real(kind=db) :: ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft18
        real(kind=db) :: qxx,qyy,qzz,qxy_7_8,qxy_9_10,qxz_15_16,qxz_17_18,qyz_11_12,qyz_13_14
        real(kind=db) :: pi2cssq1,pi2cssq2,pi2cssq0
        real(kind=db) :: addendum0,gaddendum0
        real(kind=db) :: psi_x,psi_y,psi_z,mod_psi,mod_psi_sq,st_coeff,b0,b1,b2,beta,sigma
        real(kind=db) :: one_ov_nu2,one_ov_nu1,nu_avg,rtot,rprod
        real(kind=db) :: max_press_excess,rr,ushifted,vshifted,wshifted,norm_x,norm_y,norm_z

        integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
        real(kind=db), allocatable, dimension(:,:,:) :: psi,rhoA,rhoB,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
        real(kind=db), allocatable, dimension(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
        real(kind=db), allocatable, dimension(:,:,:) :: f10,f11,f12,f13,f14,f15,f16,f17,f18
        real(kind=db), allocatable, dimension(:,:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8,g9
        real(kind=db), allocatable, dimension(:,:,:) :: g10,g11,g12,g13,g14,g15,g16,g17,g18
        real(kind=db), allocatable,dimension(:,:,:) :: random_field
        integer(kind=4), allocatable, dimension(:,:,:) :: nci_loc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!for custom geometry*********************************!
        integer:: ddrop,lc,jjd,jju
        !$if _OPENACC
        integer :: devNum
        integer(acc_device_kind) :: devType
        devType = acc_get_device_type()
        devNum=acc_get_device_num(devType)
        !$endif
        
   
    
    
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
        #ifdef _OPENACC
        ngpus=acc_get_num_devices(acc_device_nvidia)
        #else
        ngpus=0
        #endif

    !**************************************user parameters**************************
        nx=416
        ny=416
        nz=416
        nsteps=1000
        stamp=200000
        fx=0.0_db*10.0**(-7)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)

        lprint=.false.
        lvtk=.false.
        lasync=.false.
    !*****************************************allocation*******************************************************
        allocate(f0(0:nx+1,0:ny+1,0:nz+1),f1(0:nx+1,0:ny+1,0:nz+1),f2(0:nx+1,0:ny+1,0:nz+1),f3(0:nx+1,0:ny+1,0:nz+1))
        allocate(f4(0:nx+1,0:ny+1,0:nz+1),f5(0:nx+1,0:ny+1,0:nz+1),f6(0:nx+1,0:ny+1,0:nz+1),f7(0:nx+1,0:ny+1,0:nz+1))
        allocate(f8(0:nx+1,0:ny+1,0:nz+1),f9(0:nx+1,0:ny+1,0:nz+1),f10(0:nx+1,0:ny+1,0:nz+1),f11(0:nx+1,0:ny+1,0:nz+1))
        allocate(f12(0:nx+1,0:ny+1,0:nz+1),f13(0:nx+1,0:ny+1,0:nz+1),f14(0:nx+1,0:ny+1,0:nz+1),f15(0:nx+1,0:ny+1,0:nz+1))
        allocate(f16(0:nx+1,0:ny+1,0:nz+1),f17(0:nx+1,0:ny+1,0:nz+1),f18(0:nx+1,0:ny+1,0:nz+1))
        allocate(g0(0:nx+1,0:ny+1,0:nz+1),g1(0:nx+1,0:ny+1,0:nz+1),g2(0:nx+1,0:ny+1,0:nz+1),g3(0:nx+1,0:ny+1,0:nz+1))
        allocate(g4(0:nx+1,0:ny+1,0:nz+1),g5(0:nx+1,0:ny+1,0:nz+1),g6(0:nx+1,0:ny+1,0:nz+1),g7(0:nx+1,0:ny+1,0:nz+1))
        allocate(g8(0:nx+1,0:ny+1,0:nz+1),g9(0:nx+1,0:ny+1,0:nz+1),g10(0:nx+1,0:ny+1,0:nz+1),g11(0:nx+1,0:ny+1,0:nz+1))
        allocate(g12(0:nx+1,0:ny+1,0:nz+1),g13(0:nx+1,0:ny+1,0:nz+1),g14(0:nx+1,0:ny+1,0:nz+1),g15(0:nx+1,0:ny+1,0:nz+1))
        allocate(g16(0:nx+1,0:ny+1,0:nz+1),g17(0:nx+1,0:ny+1,0:nz+1),g18(0:nx+1,0:ny+1,0:nz+1))
        allocate(rhoA(1:nx,1:ny,1:nz),rhoB(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
        allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))  
        allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz),psi(0:nx+1,0:ny+1,0:nz+1))
        allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
        allocate(random_field(25,25,25), nci_loc(1:nx,1:ny,1:nz))
    
    !***************************************lattice/vars*************************************!
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
        p0=(1.0_db/3.0_db)
        p1=(1.0_db/18.0_db)
        p2=(1.0_db/36.0_db)
        p1dcssq=p1/cssq
        p2dcssq=p2/cssq
    !****************************************geometry************************
        isfluid=1
        isfluid(1,:,:)=0 !left
        isfluid(nx,:,:)=0 !right
        isfluid(:,1,:)=0 !front 
        isfluid(:,ny,:)=0 !rear
        isfluid(:,:,1)=0 !bottom
        isfluid(:,:,nz)=0 !top

    !***********************************define/read geometry if any**************************
            
            ! Lc=250+nz/2;
            ! do i=Lc,nz
            !     jjd=nint((i-Lc+1)*sind(30.0_db));
                
            !     if(jjd<ny/2-4)then
            !         isfluid(:,jjd:ny/2,i)=1; 
            !         jju=ny-jjd;
            !         isfluid(:,ny/2:jju,i)=1;
            !     endif
            ! enddo
            ! Ddrop=40;
            ! isfluid(:,2:ny-1,nz/2:Lc)=1;
            ! isfluid(:,ny/2-(Ddrop/(2)):ny/2+(Ddrop/(2)),nz/2:nz)=1;
            ! isfluid(:,:,1:nz/2)=isfluid(:,:,nz:nz/2+1:-1);
            ! isfluid(1,:,:)=0 !left
            ! isfluid(nx,:,:)=0 !right
            ! isfluid(:,1,:)=0 !front 
            ! isfluid(:,ny,:)=0 !rear
            ! isfluid(:,:,1)=0 !bottom
            ! isfluid(:,:,nz)=0 !top
            ! do k=1,nz
            !     do j=1,ny
            !         if(isfluid(nx/2,j,k).eq.3)then
            !             if(isfluid(nx/2,j+1,k).eq.1 .or. isfluid(nx/2,j-1,k).eq.1 .or. isfluid(nx/2,j,k+1).eq.1 .or. isfluid(nx/2,j,k-1).eq.1 &
            !                 .or. isfluid(nx/2,j+1,k+1).eq.1 .or. isfluid(nx/2,j+1,k-1).eq.1 .or. isfluid(nx/2,j-1,k+1).eq.1 .or. isfluid(nx/2,j-1,k-1).eq.1)then
            !                 isfluid(2:nx-1,j,k)=0
            !             endif
            !         endif
            !     enddo
            ! enddo
            ! open(231, file = 'isfluid.out', status = 'replace')
            !     do j=1,ny
            !         do k=1,nz
            !             write(231,*) isfluid(nx/2,j,k)  
            !         enddo
            !     enddo
            ! close(231)
        
    !********************************hermite projection vars**********
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
    !*********************************chromodynamics vars*****************************
        beta=0.95_db
        sigma=0.03_db
        st_coeff=(9.0_db/4.0_db)*sigma*omega
        b0=-2.0_db/9.0_db
        b1=1.0_db/54.0_db
        b2=1.0_db/27.0_db
        p1cg=(2.0_db/3.0_db)**2 * (1.0_db/6.0_db)
        p2cg=(1.0_db/6.0_db)**2 * (2.0_db/3.0_db)
        p3cg=(1.0_db/6.0_db)**3
        nci_loc=0
        max_press_excess=0.004
    !********************************initialization of macrovars ************************    
        u=0.0_db
        v=0.0_db
        w=0.0_db
        rhoA(1:nx,1:ny,1:nz)=0.0_db  !total density
        rhoB(1:nx,1:ny,1:nz)=0.0_db  !total density
        psi=-1.0_db
        radius=20
        !*****************************************Impacting droplets***************************!
            do i=(nx/2)-radius,(nx/2)+radius
                do j=ny/2-radius,ny/2+radius
                    do k=nz/2-radius,nz/2+radius
                        if ((i-(nx/2))**2+(j-ny/2)**2+(k-(nz/2))**2<=radius**2)then
                            psi(i,j,k)=1.0_db
                            u(i,j,k)=0.0_db
                        endif
                    enddo
                enddo
            enddo
            ! do i=(nx/2+radius+5)-radius,(nx/2+radius+5)+radius
            !     do j=ny/2-radius,ny/2+radius
            !         do k=nz/2-radius,nz/2+radius
            !             if ((i-(nx/2+radius+5))**2+(j-ny/2)**2+(k-(nz/2))**2<=radius**2)then
            !                 psi(i,j,k)=1.0_db
            !                 u(i,j,k)=-0.035/2.0_db
            !             endif
            !         enddo
            !     enddo
            ! enddo
        
        !*****************************************Spinodal decomposition***********************!
            ! do k=1,25
            !     do j=1,25
            !         do i=1,25
            !             call random_number(rr)
            !             random_field(i,j,k)=rr
            !         enddo
            !     enddo
            ! enddo
            ! !
            ! do k=1,25
            !     k_init=k_end+1
            !     k_end=k_end+10
            !     do j=1,25
            !         j_init=j_end+1
            !         j_end=j_end+10
            !         do i=1,25
            !             i_init=i_end+1
            !             i_end=i_end+10
            !             if(random_field(i,j,k).gt.0.5) psi(i_init:i_end,j_init:j_end,k_init:k_end)=1.0_db
            !         enddo
            !         i_init=0
            !         i_end=0
            !     enddo
            !     j_init=0
            !     j_end=0
            ! enddo
        !************************read initial conditions for macrovars*************************!
            ! open(231, file = 'psi.txt', status = 'old',action='read')
            ! do j=1,ny
            !     do k=1,nz
            !         read(231,*) psi(15,j,k)
            !     enddo
            ! enddo
            ! close(231)
            ! do i=4,nx-3
            !     psi(i,:,:)=psi(15,:,:)
            ! enddo
            
            ! open(231, file = 'psi.out', status = 'replace')
            ! do i=1,nx
            !     do j=1,ny
            !         do k=1,nz
            !         write(231,*) psi(i,j,k)  
            !         enddo
            !     enddo
            ! enddo
            ! close(231)
        !************************************single cylindrical droplets*******************!
            ! do i=5,nx-4
            !     do j=ny/2-radius,ny/2+radius
            !         do k=nz-15-radius-radius,nz-15-radius+radius
            !             if ((j-ny/2)**2+(k-(nz-15-radius))**2<=radius**2)then
            !                 psi(i,j,k)=1.0_db
            !             endif
            !         enddo
            !     enddo
            ! enddo
        !*****************************************dense emulsion in channel********************!
        !*****************************************turbulent emulsion********************!
        rhoB=0.5*(1.0_db-psi(1:nx,1:ny,1:nz))
        rhoA=1.0_db-rhoB
    !*************************************set distros************************!
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
        f0(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p0*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f1(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0+u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f2(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f3(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f4(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f5(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f6(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p1*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f7(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0+u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f8(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f9(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0+u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f10(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f11(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f12(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f13(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f14(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0-u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f15(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0 + u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2 - u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f16(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0 - u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2 - u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f17(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0 - u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2 - u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        f18(1:nx,1:ny,1:nz)=rhoa(1:nx,1:ny,1:nz)*p2*(1.0 + u(1:nx,1:ny,1:nz)/cssq + 0.5*(u(1:nx,1:ny,1:nz)/cssq)**2 - u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz)*0.5/cssq)
        !
        g0(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p0!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g1(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0+u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g2(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g3(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g4(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g5(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g6(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p1!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g7(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0+u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g8(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g9(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0+u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g10(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g11(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g12(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g13(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g14(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0-u(:,:,:)*u(:,:,:)*0.5/cssq)
        g15(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0 +u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2 - u(:,:,:)*u(:,:,:)*0.5/cssq)
        g16(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0 -u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2 - u(:,:,:)*u(:,:,:)*0.5/cssq)
        g17(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0 -u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2 - u(:,:,:)*u(:,:,:)*0.5/cssq)
        g18(1:nx,1:ny,1:nz)=rhoB(1:nx,1:ny,1:nz)*p2!*(1.0 +u(:,:,:)/cssq + 0.5*(u(:,:,:)/cssq)**2 - u(:,:,:)*u(:,:,:)*0.5/cssq)
    
    !***************************************check data ************************ 
        write(6,*) '*******************LB data*****************'
        write(6,*) 'tau',tau
        write(6,*) 'omega',omega
        write(6,*) 'visc',visc_LB
        write(6,*) 'fx',fx
        write(6,*) 'fy',fy
        write(6,*) 'fz',fz
        write(6,*) 'cssq',cssq
        write(6,*) 'beta',beta
        write(6,*) 'sigma',sigma
        write(6,*) 'surface_tens_coeff',st_coeff
        write(6,*) 'max press excess',max_press_excess
        write(6,*) '*******************INPUT data*****************'
        write(6,*) 'nx',nx
        write(6,*) 'ny',ny
        write(6,*) 'ny',nz
        write(6,*) 'nsteps',nsteps
        write(6,*) 'lprint',lprint
        write(6,*) 'lvtk',lvtk
        write(6,*) 'lasync',lasync
        write(6,*) 'stamp',stamp
        write(6,*) 'max fx',huge(fx)
        write(6,*) 'max fx',huge(fy)
        write(6,*) 'max fx',huge(fz)
        write(6,*) '*******************************************'
    !*****************************************copy data on gpu*********************************************!
        !$acc data copy(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,isfluid,p0,p1,p2,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rhoA,rhoB,u,v,w,psi, &
             !$acc& g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,nci_loc,rhoprint,velprint)
    
    
    !**********************************************************************!
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
		  !$acc kernels present(rhoprint,velprint,rhoa,u,v,w) async(1)
		  !$acc loop independent collapse(3)  private(i,j,k)
		  do k=1,nz
			do j=1,ny
			  do i=1,nx
				rhoprint(i,j,k)=real(rhoa(i,j,k),kind=4)
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
        
        call cpu_time(ts1)
    !*************************************main loop*************************!
    do step=1,nsteps 
        !$acc kernels async(1)
        !$acc loop collapse (3) 
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
                        pxx(i,j,k)=0.0_db
                        pyy(i,j,k)=0.0_db
                        pzz(i,j,k)=0.0_db
                        pxy(i,j,k)=0.0_db
                        pxz(i,j,k)=0.0_db
                        pyz(i,j,k)=0.0_db

                        rhoA(i,j,k) = f0(i,j,k)+f1(i,j,k)+f2(i,j,k)+f3(i,j,k)+f4(i,j,k)+f5(i,j,k) &
                                     +f6(i,j,k)+f7(i,j,k)+f8(i,j,k)+f9(i,j,k)+f10(i,j,k)+f11(i,j,k) &
                                     +f12(i,j,k)+f13(i,j,k)+f14(i,j,k)+f15(i,j,k)+f16(i,j,k)+f17(i,j,k) &
                                     +f18(i,j,k)

                        rhoB(i,j,k) = g0(i,j,k)+g1(i,j,k)+g2(i,j,k)+g3(i,j,k)+g4(i,j,k)+g5(i,j,k) &
                                     +g6(i,j,k)+g7(i,j,k)+g8(i,j,k)+g9(i,j,k)+g10(i,j,k)+g11(i,j,k) &
                                     +g12(i,j,k)+g13(i,j,k)+g14(i,j,k)+g15(i,j,k)+g16(i,j,k)+g17(i,j,k) &
                                     +g18(i,j,k)

                        nci_loc(i,j,k)=0
                        rtot=rhoA(i,j,k)+rhoB(i,j,k)

                        psi(i,j,k)= (rhoA(i,j,k)-rhoB(i,j,k))/(rhoA(i,j,k)+rhoB(i,j,k))

                        u(i,j,k) = (f1(i,j,k)+g1(i,j,k)+f7(i,j,k)+g7(i,j,k)+f9(i,j,k)+g9(i,j,k)+&
                                    f15(i,j,k)+g15(i,j,k)+f18(i,j,k)+g18(i,j,k)) &
                                    -(f2(i,j,k)+g2(i,j,k)+f8(i,j,k)+g8(i,j,k)+f10(i,j,k)+ &
                                        g10(i,j,k)+f16(i,j,k)+g16(i,j,k)+f17(i,j,k)+g17(i,j,k)) 
                        
                        v(i,j,k) = (f3(i,j,k)+g3(i,j,k)+f7(i,j,k)+g7(i,j,k)+f10(i,j,k)+ &
                                    g10(i,j,k)+f11(i,j,k)+g11(i,j,k)+f13(i,j,k)+g13(i,j,k)) &
                                    - (f4(i,j,k)+g4(i,j,k)+f8(i,j,k)+g8(i,j,k)+f9(i,j,k)+ &
                                        g9(i,j,k)+f12(i,j,k)+g12(i,j,k)+f14(i,j,k)+g14(i,j,k))

                        w(i,j,k) = (f5(i,j,k)+g5(i,j,k)+f11(i,j,k)+g11(i,j,k)+f14(i,j,k)+ &
                                    g14(i,j,k)+f15(i,j,k)+g15(i,j,k)+f17(i,j,k)+g17(i,j,k))-&
                                        (f6(i,j,k)+g6(i,j,k)+f12(i,j,k)+g12(i,j,k)+f13(i,j,k)+ &
                                         g13(i,j,k)+f16(i,j,k)+g16(i,j,k)+f18(i,j,k)+g18(i,j,k))
                        
                        uu=0.5_db*(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))/cssq
                        !1-2
                        udotc=u(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f1(i,j,k)+g1(i,j,k)-p1*(rtot+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        fneq1=f2(i,j,k)+g2(i,j,k)-p1*(rtot+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1

                        !3-4
                        udotc=v(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        feq=p1*(rtot+(temp + udotc))
                        fneq1=f3(i,j,k)+g3(i,j,k)-p1*(rtot+(temp + udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        fneq1=f4(i,j,k)+g4(i,j,k)-p1*(rtot+(temp - udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        !5-6
                        udotc=w(i,j,k)/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f5(i,j,k)+g5(i,j,k)-p1*(rtot+(temp + udotc))
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        fneq1=f6(i,j,k)+g6(i,j,k)-p1*(rtot+(temp - udotc))
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        !7-8
                        udotc=(u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f7(i,j,k)+g7(i,j,k)-p2*(rtot+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pxy(i,j,k)=pxy(i,j,k)+fneq1
                        fneq1=f8(i,j,k)+g8(i,j,k)-p2*(rtot+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pxy(i,j,k)=pxy(i,j,k)+fneq1
                        !10-9
                        udotc=(-u(i,j,k)+v(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f10(i,j,k)+g10(i,j,k)-p2*(rtot+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pxy(i,j,k)=pxy(i,j,k)-fneq1
                        fneq1=f9(i,j,k)+g9(i,j,k)-p2*(rtot+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pxy(i,j,k)=pxy(i,j,k)-fneq1
                        !11-12
                        udotc=(v(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f11(i,j,k)+g11(i,j,k)-p2*(rtot+(temp + udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pyz(i,j,k)=pyz(i,j,k)+fneq1
                        fneq1=f12(i,j,k)+g12(i,j,k)-p2*(rtot+(temp - udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pyz(i,j,k)=pyz(i,j,k)+fneq1
                        !13-14
                        udotc=(v(i,j,k)-w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f13(i,j,k)+g13(i,j,k) - p2*(rtot+(temp + udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pyz(i,j,k)=pyz(i,j,k)-fneq1
                        fneq1=f14(i,j,k)+g14(i,j,k) - p2*(rtot+(temp - udotc))
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pyz(i,j,k)=pyz(i,j,k)-fneq1
                        !15-16
                        udotc=(u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f15(i,j,k)+g15(i,j,k)-p2*(rtot+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pxz(i,j,k)=pxz(i,j,k)+fneq1
                        fneq1=f16(i,j,k)+g16(i,j,k)-p2*(rtot+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pxz(i,j,k)=pxz(i,j,k)+fneq1
                        !17-18
                        udotc=(-u(i,j,k)+w(i,j,k))/cssq
                        temp = -uu + 0.5_db*udotc*udotc
                        fneq1=f17(i,j,k)+g17(i,j,k)-p2*(rtot+(temp + udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pxz(i,j,k)=pxz(i,j,k)-fneq1
                        fneq1=f18(i,j,k)+g18(i,j,k)-p2*(rtot+(temp - udotc))
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pxz(i,j,k)=pxz(i,j,k)-fneq1
                        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
                        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
                        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        
        !****************************************writeonfile***************************************************!
        if(mod(step,stamp).eq.0)write(6,'(a,i8)')'step : ',step
        if(lprint)then
          if(mod(step,stamp).eq.0)then
            iframe=iframe+1
            !$acc wait(1)
            !$acc kernels present(rhoprint,velprint,rhoa,u,v,w) async(1)
            !$acc loop independent collapse(3)  private(i,j,k)
            do k=1,nz
              do j=1,ny
                do i=1,nx
                  rhoprint(i,j,k)=real(rhoa(i,j,k),kind=4)
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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !$acc kernels 
        !$acc loop collapse (3) 
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                    !!*******************************near contact***************************!
                        if(isfluid(i,j,k).eq.1)then
                            
                            if(psi(i,j,k).lt.-0.9_db .and. i.gt.3 .and. i.lt.nx-2 .and. j.gt.3 .and. j.lt.ny-2 .and. k.gt.3 .and. k.lt.nz-2)then
                                
                                if (psi(i+3,j,k).gt.-0.85 .and. psi(i-3,j,k).gt.-0.85)then
                                    nci_loc(i+3,j,k)=1
                                    nci_loc(i-3,j,k)=1
                                endif
                                if (psi(i,j+3,k).gt.-0.85 .and. psi(i,j-3,k).gt.-0.85)then
                                    nci_loc(i,j+3,k)=1
                                    nci_loc(i,j-3,k)=1
                                endif
                                if (psi(i,j,k+3).gt.-0.85 .and. psi(i,j,k-3).gt.-0.85)then!(psi(i,j,k+3).gt.-0.85 .and. psi(i,j,k+3).lt.0.0_db .and. psi(i,j,k-3).gt.-0.85 .and. psi(i,j,k-3).lt.0.0_db)then
                                   nci_loc(i,j,k+3)=1
                                   nci_loc(i,j,k-3)=1
                                endif
                                if (psi(i+2,j+2,k).gt.-0.85 .and. psi(i-2,j-2,k).gt.-0.85)then!(psi(i+3,j+3,k).gt.-0.85 .and. psi(i+3,j+3,k).lt.0.0_db .and. psi(i-3,j-3,k).gt.-0.85 .and. psi(i-3,j-3,k).lt.0.0_db)then
                                    nci_loc(i+2,j+2,k)=1
                                    nci_loc(i-2,j-2,k)=1
                                endif
                                if (psi(i+2,j,k+2).gt.-0.85.and. psi(i-2,j,k-2).gt.-0.85)then!(psi(i+3,j,k+3).gt.-0.85 .and. psi(i+3,j,k+3).lt.0.0_db .and. psi(i-3,j,k-3).gt.-0.85 .and. psi(i-3,j,k-3).lt.0.0_db)then
                                    nci_loc(i+2,j,k+2)=1
                                    nci_loc(i-2,j,k-2)=1
                                endif
                                if (psi(i,j+2,k+2).gt.-0.85 .and. psi(i,j-2,k-2).gt.-0.85)then!(psi(i,j+3,k+3).gt.-0.85 .and. psi(i,j+3,k+3).lt.0.0_db .and. psi(i,j-3,k-3).gt.-0.85 .and. psi(i,j-3,k-3).lt.0.0_db)then
                                    nci_loc(i,j+2,k+2)=1
                                    nci_loc(i,j-2,k-2)=1
                                endif

                                if (psi(i-2,j+2,k).gt.-0.85.and. psi(i+2,j-2,k).gt.-0.85)then!(psi(i-3,j+3,k).gt.-0.85 .and. psi(i-3,j+3,k).lt.0.0_db .and. psi(i+3,j-3,k).gt.-0.85 .and. psi(i+3,j-3,k).lt.0.0_db)then
                                    nci_loc(i-2,j+2,k)=1
                                    nci_loc(i+2,j-2,k)=1
                                endif
                                if (psi(i-2,j,k+2).gt.-0.85.and. psi(i+2,j,k-2).gt.-0.85)then!(psi(i-3,j,k+3).gt.-0.85 .and. psi(i-3,j,k+3).lt.0.0_db .and. psi(i+3,j,k-3).gt.-0.85 .and. psi(i+3,j,k-3).lt.0.0_db)then
                                    nci_loc(i+2,j,k-2)=1
                                    nci_loc(i-2,j,k+2)=1
                                endif
                                if (psi(i,j-2,k+2).gt.-0.85.and. psi(i,j+2,k-2).gt.-0.85)then!(psi(i,j-3,k+3).gt.-0.85 .and. psi(i,j-3,k+3).lt.0.0_db .and. psi(i,j+3,k-3).gt.-0.85 .and. psi(i,j+3,k-3).lt.0.0_db)then
                                    nci_loc(i,j-2,k+2)=1
                                    nci_loc(i,j+2,k-2)=1
                                endif
                                
                            endif
                        endif
                    enddo
                enddo
            enddo
        !$acc loop collapse (3) 
        do k=1,nz
            do j=1,ny
                do i=1,nx  
                    if(isfluid(i,j,k).eq.1)then 
                    !!*******************************chromodynamics***************************!
                        psi_x=(1.0_db/cssq)*(p1cg*(psi(i+1,j,k)-psi(i-1,j,k)) + &
                                            p2cg*(psi(i+1,j+1,k)+psi(i+1,j-1,k)+psi(i+1,j,k+1)+psi(i+1,j,k-1)) &
                                           -p2cg*(psi(i-1,j+1,k)+psi(i-1,j-1,k)+psi(i-1,j,k+1)+psi(i-1,j,k-1)) &
                                          + p3cg*(psi(i+1,j+1,k+1)+psi(i+1,j+1,k-1)+psi(i+1,j-1,k+1)+psi(i+1,j-1,k-1)) &
                                           -p3cg*(psi(i-1,j+1,k+1)+psi(i-1,j+1,k-1)+psi(i-1,j-1,k+1)+psi(i-1,j-1,k-1))) 


                        psi_y=(1.0_db/cssq)*(p1cg*(psi(i,j+1,k)-psi(i,j-1,k)) &
                                            +p2cg*(psi(i+1,j+1,k)+psi(i-1,j+1,k)+psi(i,j+1,k+1)+psi(i,j+1,k-1)) &
                                            -p2cg*(psi(i+1,j-1,k)+psi(i-1,j-1,k)+psi(i,j-1,k+1)+psi(i,j-1,k-1)) &
                                            +p3cg*(psi(i+1,j+1,k+1)+psi(i-1,j+1,k+1)+psi(i+1,j+1,k-1)+psi(i-1,j+1,k-1)) &
                                            -p3cg*(psi(i+1,j-1,k+1)+psi(i-1,j-1,k+1)+psi(i+1,j-1,k-1)+psi(i-1,j-1,k-1)))


                        psi_z=(1.0_db/cssq)*(p1cg*(psi(i,j,k+1)-psi(i,j,k-1)) + &
                                            p2cg*(psi(i+1,j,k+1)+psi(i-1,j,k+1)+ psi(i,j+1,k+1)+psi(i,j-1,k+1)) &
                                           -p2cg*(psi(i+1,j,k-1)+psi(i-1,j,k-1)+psi(i,j+1,k-1)+psi(i,j-1,k-1)) &
                                           +p3cg*(psi(i+1,j+1,k+1)+psi(i-1,j+1,k+1)+psi(i+1,j-1,k+1)+psi(i-1,j-1,k+1)) &
                                           -p3cg*(psi(i+1,j+1,k-1)+psi(i-1,j+1,k-1)+psi(i+1,j-1,k-1)+psi(i-1,j-1,k-1)))
                        
                        mod_psi=sqrt(psi_x**2+psi_y**2+psi_z**2)

                        norm_x=0.0_db
                        norm_y=0.0_db
                        norm_z=0.0_db

                        mod_psi_sq=psi_x**2 + psi_y**2 +psi_z**2 
                        
                        rtot=0.0_db
                        
                        rtot=rhoA(i,j,k)+rhoB(i,j,k)
                        
                        rprod=rhoA(i,j,k)*rhoB(i,j,k)
                        
                        nu_avg=1.0_db/(rhoA(i,j,k)*one_ov_nu1/rtot + rhoB(i,j,k)*one_ov_nu2/rtot)
                        
                        omega=2.0_db/(6.0_db*nu_avg + 1.0_db)
                        
                        st_coeff=(9.0_db/4.0_db)*sigma*omega
                        addendum0=0.0_db
                        gaddendum0=0.0_db  
                        !******************************collision+stream+perturbation*****************************!
                        if (mod_psi>0.001)then

                            norm_x=psi_x/mod_psi
                            norm_y=psi_y/mod_psi
                            norm_z=psi_z/mod_psi

                            ushifted=u(i,j,k) + fx + float(nci_loc(i,j,k))*(norm_x)*max_press_excess*abs(rhoB(i,j,k))
                            vshifted=v(i,j,k) + fy + float(nci_loc(i,j,k))*(norm_y)*max_press_excess*abs(rhoB(i,j,k))
                            wshifted=w(i,j,k) + fz + float(nci_loc(i,j,k))*(norm_z)*max_press_excess*abs(rhoB(i,j,k))
                            uu=0.5_db*(ushifted*ushifted + vshifted*vshifted + wshifted*wshifted)/cssq 
                           
                            !0
                            addendum0=-st_coeff*mod_psi*b0
                            feq=p0*(rtot-uu)
                            fpc=feq + addendum0 + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k))) 
                            f0(i,j,k)=fpc*rhoA(i,j,k)/rtot
                            g0(i,j,k)=fpc*rhoB(i,j,k)/rtot
                            !1
                            gaddendum0=p1*(rtot)*(rprod*beta*psi_x/mod_psi/rtot**2)    
                            addendum0=st_coeff*mod_psi*(p1*psi_x**2/mod_psi_sq - b1)
                            udotc=ushifted/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p1*(rtot+(temp + udotc))
                            fpc=feq + addendum0 + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                            f1(i+1,j,k)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g1(i+1,j,k)=fpc*rhob(i,j,k)/rtot - gaddendum0
                            !2
                            feq=p1*(rtot+(temp - udotc))
                            fpc=feq + addendum0 + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))) 
                            f2(i-1,j,k)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g2(i-1,j,k)=fpc*rhob(i,j,k)/rtot + gaddendum0
                            
                            !3
                            gaddendum0=p1*(rtot)*(rprod*beta*psi_y/mod_psi/rtot**2)
                            addendum0=st_coeff*mod_psi*(p1*psi_y**2/mod_psi_sq - b1)
                            udotc=vshifted/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p1*(rtot+(temp + udotc))
                            fpc=feq +addendum0+(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))) 
                            f3(i,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g3(i,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum0
                            !4
                            feq=p1*(rtot+(temp - udotc))
                            fpc=feq +addendum0+(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))) 
                            f4(i,j-1,k)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g4(i,j-1,k)=fpc*rhob(i,j,k)/rtot + gaddendum0
                            
                            !5
                            addendum0=st_coeff*mod_psi*(p1*psi_z**2/mod_psi_sq - b1)
                            gaddendum0=p1*(rtot)*(rprod*beta*psi_z/mod_psi/rtot**2)
                            udotc=wshifted/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p1*(rtot+(temp + udotc))
                            fpc=feq + addendum0+(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))) 
                            f5(i,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g5(i,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum0
                            !6
                            feq=p1*(rtot+(temp - udotc))
                            fpc=feq + addendum0 +(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))) 
                            f6(i,j,k-1)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g6(i,j,k-1)=fpc*rhob(i,j,k)/rtot + gaddendum0
                            
                            !7
                            addendum0=st_coeff*mod_psi*(p2*(psi_x+psi_y)**2/mod_psi_sq - b2)
                            gaddendum0=p2*(rtot)*(rprod*beta*(psi_x/mod_psi+psi_y/mod_psi)/rtot**2)
                            udotc=(ushifted+vshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq + addendum0+(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*pxy(i,j,k))
                            f7(i+1,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g7(i+1,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum0
                            !8
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq + addendum0+ (1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*pxy(i,j,k))
                            f8(i-1,j-1,k)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g8(i-1,j-1,k)=fpc*rhob(i,j,k)/rtot + gaddendum0
                            
                            !10
                            addendum0=st_coeff*mod_psi*(p2*(-psi_x+psi_y)**2/mod_psi_sq - b2)
                            gaddendum0=p2*(rtot)*(rprod*beta*(-psi_x/mod_psi+psi_y/mod_psi)/rtot**2)
                            udotc=(-ushifted+vshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq + addendum0 + (1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)-2.0_db*pxy(i,j,k)) 
                            f10(i-1,j+1,k)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g10(i-1,j+1,k)=fpc*rhob(i,j,k)/rtot - gaddendum0
                            !9
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq + addendum0 + (1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)-2.0_db*pxy(i,j,k))
                            f9(i+1,j-1,k)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g9(i+1,j-1,k)=fpc*rhob(i,j,k)/rtot + gaddendum0
                            
                            !15
                            addendum0=st_coeff*mod_psi*(p2*(psi_x+psi_z)**2/mod_psi_sq - b2)
                            gaddendum0=p2*(rtot)*(rprod*beta*(psi_x/mod_psi+psi_z/mod_psi)/rtot**2)
                            udotc=(ushifted+wshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq + addendum0+(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*pxz(i,j,k)) 
                            f15(i+1,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g15(i+1,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum0
                            !16
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq + addendum0 +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*pxz(i,j,k)) 
                            f16(i-1,j,k-1)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g16(i-1,j,k-1)=fpc*rhob(i,j,k)/rtot + gaddendum0
                            
                            !17
                            addendum0=st_coeff*mod_psi*(p2*(-psi_x+psi_z)**2/mod_psi_sq - b2)
                            gaddendum0=p2*(rtot)*(rprod*beta*(-psi_x/mod_psi+psi_z/mod_psi)/rtot**2)
                            udotc=(-ushifted+wshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq +addendum0+(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)-2.0_db*pxz(i,j,k)) 
                            f17(i-1,j,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g17(i-1,j,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum0   
                            !18
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq + addendum0+(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)-2.0_db*pxz(i,j,k)) 
                            f18(i+1,j,k-1)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g18(i+1,j,k-1)=fpc*rhob(i,j,k)/rtot + gaddendum0

                            !11
                            addendum0=st_coeff*mod_psi*(p2*(psi_y+psi_z)**2/mod_psi_sq - b2)
                            gaddendum0=p2*(rtot)*(rprod*beta*(psi_y/mod_psi+psi_z/mod_psi)/rtot**2)
                            udotc=(vshifted+wshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq + addendum0+(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*pyz(i,j,k))
                            f11(i,j+1,k+1)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g11(i,j+1,k+1)=fpc*rhob(i,j,k)/rtot - gaddendum0
                            !12
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq + addendum0+(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*pyz(i,j,k)) 
                            f12(i,j-1,k-1)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g12(i,j-1,k-1)=fpc*rhob(i,j,k)/rtot + gaddendum0
                            
                            !13
                            addendum0=st_coeff*mod_psi*(p2*(psi_y-psi_z)**2/mod_psi_sq - b2)
                            gaddendum0=p2*(rtot)*(rprod*beta*(psi_y/mod_psi-psi_z/mod_psi)/rtot**2)
                            udotc=(vshifted-wshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq + addendum0+(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)-2.0_db*pyz(i,j,k)) 
                            f13(i,j+1,k-1)=fpc*rhoA(i,j,k)/rtot + gaddendum0
                            g13(i,j+1,k-1)=fpc*rhob(i,j,k)/rtot - gaddendum0
                            !14
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq + addendum0+(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)-2.0_db*pyz(i,j,k)) 
                            f14(i,j-1,k+1)=fpc*rhoA(i,j,k)/rtot - gaddendum0
                            g14(i,j-1,k+1)=fpc*rhob(i,j,k)/rtot + gaddendum0

                        else   
                    !******************************collision+stream*****************************!
                            !0
                            ushifted=u(i,j,k) + fx 
                            vshifted=v(i,j,k) + fy 
                            wshifted=w(i,j,k) + fz 

                            uu=0.5_db*(ushifted*ushifted + vshifted*vshifted + wshifted*wshifted)/cssq 
                           
                            !0
                            feq=p0*(rtot-uu)
                            fpc=feq  + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k))) 
                            f0(i,j,k)=fpc*rhoA(i,j,k)/rtot
                            g0(i,j,k)=fpc*rhoB(i,j,k)/rtot
                            !1
                            udotc=ushifted/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p1*(rtot+(temp + udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                            f1(i+1,j,k)=fpc*rhoA(i,j,k)/rtot 
                            g1(i+1,j,k)=fpc*rhob(i,j,k)/rtot 
                            !2
                            feq=p1*(rtot+(temp - udotc))
                            fpc=feq  + (1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k))) 
                            f2(i-1,j,k)=fpc*rhoA(i,j,k)/rtot 
                            g2(i-1,j,k)=fpc*rhob(i,j,k)/rtot 
                            
                            !3
                            udotc=vshifted/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p1*(rtot+(temp + udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))) 
                            f3(i,j+1,k)=fpc*rhoA(i,j,k)/rtot 
                            g3(i,j+1,k)=fpc*rhob(i,j,k)/rtot 
                            !4
                            feq=p1*(rtot+(temp - udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k))) 
                            f4(i,j-1,k)=fpc*rhoA(i,j,k)/rtot 
                            g4(i,j-1,k)=fpc*rhob(i,j,k)/rtot 
                            
                            !5
                            udotc=wshifted/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p1*(rtot+(temp + udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))) 
                            f5(i,j,k+1)=fpc*rhoA(i,j,k)/rtot 
                            g5(i,j,k+1)=fpc*rhob(i,j,k)/rtot 
                            !6
                            feq=p1*(rtot+(temp - udotc))
                            fpc=feq  +(1.0_db-omega)*pi2cssq1*((1.0_db-cssq)*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k))) 
                            f6(i,j,k-1)=fpc*rhoA(i,j,k)/rtot 
                            g6(i,j,k-1)=fpc*rhob(i,j,k)/rtot 
                            
                            !7
                            udotc=(ushifted+vshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*pxy(i,j,k))
                            f7(i+1,j+1,k)=fpc*rhoA(i,j,k)/rtot 
                            g7(i+1,j+1,k)=fpc*rhob(i,j,k)/rtot 
                            !8
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq + (1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*pxy(i,j,k))
                            f8(i-1,j-1,k)=fpc*rhoA(i,j,k)/rtot 
                            g8(i-1,j-1,k)=fpc*rhob(i,j,k)/rtot 
                            
                            !10
                            udotc=(-ushifted+vshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq  + (1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)-2.0_db*pxy(i,j,k)) 
                            f10(i-1,j+1,k)=fpc*rhoA(i,j,k)/rtot 
                            g10(i-1,j+1,k)=fpc*rhob(i,j,k)/rtot 
                            !9
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq + + (1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pyy(i,j,k)-cssq*pzz(i,j,k)-2.0_db*pxy(i,j,k))
                            f9(i+1,j-1,k)=fpc*rhoA(i,j,k)/rtot 
                            g9(i+1,j-1,k)=fpc*rhob(i,j,k)/rtot 
                            
                            !15
                            udotc=(ushifted+wshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*pxz(i,j,k)) 
                            f15(i+1,j,k+1)=fpc*rhoA(i,j,k)/rtot 
                            g15(i+1,j,k+1)=fpc*rhob(i,j,k)/rtot 
                            !16
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq  +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*pxz(i,j,k)) 
                            f16(i-1,j,k-1)=fpc*rhoA(i,j,k)/rtot 
                            g16(i-1,j,k-1)=fpc*rhob(i,j,k)/rtot 
                            
                            !17
                            udotc=(-ushifted+wshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)-2.0_db*pxz(i,j,k)) 
                            f17(i-1,j,k+1)=fpc*rhoA(i,j,k)/rtot 
                            g17(i-1,j,k+1)=fpc*rhob(i,j,k)/rtot    
                            !18
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pxx(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pyy(i,j,k)-2.0_db*pxz(i,j,k)) 
                            f18(i+1,j,k-1)=fpc*rhoA(i,j,k)/rtot 
                            g18(i+1,j,k-1)=fpc*rhob(i,j,k)/rtot 

                            !11
                            udotc=(vshifted+wshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*pyz(i,j,k))
                            f11(i,j+1,k+1)=fpc*rhoA(i,j,k)/rtot 
                            g11(i,j+1,k+1)=fpc*rhob(i,j,k)/rtot 
                            !12
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*pyz(i,j,k)) 
                            f12(i,j-1,k-1)=fpc*rhoA(i,j,k)/rtot 
                            g12(i,j-1,k-1)=fpc*rhob(i,j,k)/rtot 
                            
                            !13
                            udotc=(vshifted-wshifted)/cssq
                            temp = -uu + 0.5_db*udotc*udotc
                            feq=p2*(rtot+(temp + udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)-2.0_db*pyz(i,j,k)) 
                            f13(i,j+1,k-1)=fpc*rhoA(i,j,k)/rtot 
                            g13(i,j+1,k-1)=fpc*rhob(i,j,k)/rtot 
                            !14
                            feq=p2*(rtot+(temp - udotc))
                            fpc=feq +(1.0_db-omega)*pi2cssq2*((1.0_db-cssq)*pyy(i,j,k)+(1.0_db-cssq)*pzz(i,j,k)-cssq*pxx(i,j,k)-2.0_db*pyz(i,j,k)) 
                            f14(i,j-1,k+1)=fpc*rhoA(i,j,k)/rtot 
                            g14(i,j-1,k+1)=fpc*rhob(i,j,k)/rtot 
                        endif
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
                            psi(i,j,k)=-1.0_db
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
                            !************************!
                            g18(i+1,j,k-1)=g17(i,j,k) !gpc 
                            g17(i-1,j,k+1)=g18(i,j,k) !hpc

                            g16(i-1,j,k-1)=g15(i,j,k) !gpc 
                            g15(i+1,j,k+1)=g16(i,j,k) !hpc

                            g14(i,j-1,k+1)=g13(i,j,k)!gpc 
                            g13(i,j+1,k-1)=g14(i,j,k)!hpc
                            
                            g12(i,j-1,k-1)=g11(i,j,k)!gpc 
                            g11(i,j+1,k+1)=g12(i,j,k)!hpc

                            g10(i-1,j+1,k)=g9(i,j,k)!gpc 
                            g9(i+1,j-1,k)=g10(i,j,k)!hpc

                            g8(i-1,j-1,k)=g7(i,j,k)!gpc 
                            g7(i+1,j+1,k)=g8(i,j,k)!hpc

                            g6(i,j,k-1)=g5(i,j,k)!gpc 
                            g5(i,j,k+1)=g6(i,j,k)!hpc 


                            g4(i,j-1,k)=g3(i,j,k)!gpc 
                            g3(i,j+1,k)=g4(i,j,k)!hpc 

                            g2(i-1,j,k)=g1(i,j,k)!gpc 
                            g1(i+1,j,k)=g2(i,j,k)!hpc
                        endif
                    enddo
                enddo
            enddo
        !*********************************call bcs(other than no slip)************************
       
            ! !$acc loop independent 
            ! do j=1,ny
            !     !$acc loop independent 
            !     do i=1,nx
            !         psi(i,j,nz)=psi(i,j,2)
            !         psi(i,j,1)=psi(i,j,ny-1)
            !         !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
            !         !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
            !         !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
            !         f6(i,j,nz-1)=f6(i,j,1)
            !         f12(i,j,nz-1)=f12(i,j,1)
            !         f13(i,j,nz-1)=f13(i,j,1)
            !         f16(i,j,nz-1)=f16(i,j,1)
            !         f18(i,j,nz-1)=f18(i,j,1)

            !         g6(i,j,nz-1)=g6(i,j,1)
            !         g12(i,j,nz-1)=g12(i,j,1)
            !         g13(i,j,nz-1)=g13(i,j,1)
            !         g16(i,j,nz-1)=g16(i,j,1)
            !         g18(i,j,nz-1)=g18(i,j,1)

            !         f5(i,j,2)=f5(i,j,nz)
            !         f11(i,j,2)=f11(i,j,nz)
            !         f14(i,j,2)=f14(i,j,nz)
            !         f15(i,j,2)=f15(i,j,nz)
            !         f17(i,j,2)=f17(i,j,nz)

            !         g5(i,j,2)=g5(i,j,nz)
            !         g11(i,j,2)=g11(i,j,nz)
            !         g14(i,j,2)=g14(i,j,nz)
            !         g15(i,j,2)=g15(i,j,nz)
            !         g17(i,j,2)=g17(i,j,nz)
            !     enddo
            ! enddo
        !$acc end kernels 


        
    enddo 
    call cpu_time(ts2)
    !$acc end data
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)*real(1.d-9,kind=db)/(ts2-ts1)
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
