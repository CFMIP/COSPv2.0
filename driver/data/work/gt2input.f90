program gt2input
  use netcdf
  implicit none
  character(64) :: dir
  integer :: tmax
  namelist /setup/ dir,tmax
  !---
  integer,parameter :: imax=256,jmax=128,kmax=40,bmax=41,hmax=9
  character(16) :: header(64)
  character(16) :: ikind,jkind,kkind,tkind
  integer :: igt,jgt,kgt,tgt
  real(4) :: var4s(imax,jmax),var4l(imax,jmax,kmax),var4b(imax,jmax,bmax)
  real(8) :: var8s(imax,jmax),var8l(imax,jmax,kmax),var8b(imax,jmax,bmax)
  !^^^
  integer,allocatable :: landmask(:,:,:)
  !---
  integer :: fid,vid,xid,yid,zid,bid,hid,tid
  !---
  integer :: i,j,k,t,lev

  open(unit=77,file='setup.nml',status='old',action='read')
  read(77,nml=setup)

  call nf90_check( nf90_create('../MIROC_inputs/'//trim(dir)//'_MIROCinput.nc', &
       &           NF90_NETCDF4 + NF90_CLOBBER, fid) )

  call nf90_check( nf90_def_dim(fid,'lon'  ,imax,  xid) )
  call nf90_check( nf90_def_dim(fid,'lat'  ,jmax,  yid) )
  call nf90_check( nf90_def_dim(fid,'level',kmax,  zid) )
  call nf90_check( nf90_def_dim(fid,'hydro',hmax,  hid) )
  call nf90_check( nf90_def_dim(fid,'time',NF90_UNLIMITED,tid) )

  ! ===== lendmask =====
  open(unit=10,file='../MIROC_inputs/'//trim(dir)//'/landmask.gt',&
       & form='unformatted',access='sequential',&
       & status='old',action='read',convert='big_endian')
  read(10) header
  ikind = header(29) ; read(header(31),'(I16)') igt
  jkind = header(32) ; read(header(34),'(I16)') jgt
  kkind = header(35) ; read(header(37),'(I16)') kgt
  !tkind = header(26) ; read(header(28),'(I16)') tgt
  if (igt /= imax .or. jgt /= jmax .or. kgt /= 1) then
     write(*,*) 'Ngrid error: ',igt,jgt,kgt,tgt ; stop
  end if
  rewind(10)

  allocate(landmask(imax,jmax,tmax))
  do t=1,tmax
     read(10) header
     
     if (trim(header(38)) == 'UR4') then
        read(10) var4s(1:imax,1:jmax)
        landmask(:,:,t) = nint(var4s)
     else if (trim(header(38)) == 'UR8') then
        read(10) var8s(1:imax,1:jmax)
        landmask(:,:,t) = nint(var8s)
     end if
  end do
  close(10)

  
  call nf90_check( nf90_def_var(fid,'landmask',NF90_INT,(/xid,yid,tid/),vid) )
  call nf90_check( nf90_put_var(fid,vid,landmask,start=(/1,1,1/), count=(/imax,jmax,tmax/) ) )
  
contains
  subroutine nf90_check(stat)
    implicit none
    integer,intent(in) :: stat
    
    if (stat /= nf90_noerr) then
       write(*,*) 'nf90 ERROR: ',trim(nf90_strerror(stat))
       stop
    end if
  end subroutine nf90_check
  
end program gt2input
