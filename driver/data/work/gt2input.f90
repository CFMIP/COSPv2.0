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
  real(8) :: clon(imax),clat(jmax)
  real(8),allocatable :: out2D(:,:,:),out3D(:,:,:,:),out3B(:,:,:,:)
  real(8),allocatable :: ReffIN(:,:,:,:,:)
  integer,parameter :: Nv2D=4
  character(16),parameter :: vname2D(2,Nv2D) = reshape( &
       & (/ 'landmask   ' , 'landmask   ' &
       & ,  'skt        ' , 'skt        ' &
       & ,  'surfelev   ' , 'orography  ' &
       & ,  'sunlit     ' , 'sunlit     ' &
       & /), shape=(/ 2,Nv2D /) )
  integer,parameter :: Nv3D=22
  character(16),parameter :: vname3D(2,Nv3D) = reshape( &
       & (/ 'gdp        ' , 'pfull      ' &
       & ,  'zlev       ' , 'height     ' &
       & ,  'zlevh      ' , 'height_half' &
       & ,  'temp       ' , 'T_abs      ' &
       & ,  'gdsh       ' , 'qv         ' &
       & ,  'gdrh       ' , 'rh         ' &
       & ,  'gwvel      ' , 'gwvel      ' &       
       & ,  'cf_ttl     ' , 'tca        ' &
       & ,  'cf_cnv     ' , 'cca        ' &
       & ,  'mr_lscliq  ' , 'mr_lsliq   ' &
       & ,  'mr_lscice  ' , 'mr_lsice   ' &
       & ,  'mr_cvcliq  ' , 'mr_ccliq   ' &
       & ,  'mr_cvcice  ' , 'mr_ccice   ' &
       & ,  'fl_lsrain  ' , 'fl_lsrain  ' &
       & ,  'fl_lssnow  ' , 'fl_lssnow  ' &
       & ,  'fl_cvrain  ' , 'fl_ccrain  ' &
       & ,  'fl_cvsnow  ' , 'fl_ccsnow  ' &
       & ,  'dtau_s     ' , 'dtau_s     ' &
       & ,  'dtau_c     ' , 'dtau_c     ' &
       & ,  'dem_s      ' , 'dem_s      ' &
       & ,  'dem_c      ' , 'dem_c      ' &
       & ,  'mr_ozone   ' , 'mr_ozone   ' &
       & /), shape=(/ 2,Nv3D /) )
  integer,parameter :: Nv3B=2
  character(16),parameter :: vname3B(2,Nv3B) = reshape( &
       & (/ 'gdpm       ' , 'phalf      ' &
       & ,  'gcumf      ' , 'gcumf      ' &
       & /), shape=(/ 2,Nv3B /) )
  !---
  integer :: fid,vid,xid,yid,zid,bid,hid,tid
  !---
  integer :: i,j,k,t,lev,n
  character(1) :: tp

  open(unit=77,file='setup.nml',status='old',action='read')
  read(77,nml=setup)
  allocate(out2D(imax,jmax,tmax))
  allocate(out3D(imax,jmax,kmax,tmax))
  allocate(out3B(imax,jmax,bmax,tmax))
  allocate(ReffIN(imax,jmax,kmax,9,tmax))

  ! @@@@@ output NetCDF @@@@@
  call nf90_check( nf90_create('../MIROC_inputs/'//trim(dir)//'_MIROCinput.nc', &
       &           IOR(NF90_NETCDF4, NF90_CLOBBER), fid) )
  call nf90_check( nf90_def_dim(fid,'lon'  ,imax,  xid) )
  call nf90_check( nf90_def_dim(fid,'lat'  ,jmax,  yid) )
  call nf90_check( nf90_def_dim(fid,'level',kmax,  zid) )
  !call nf90_check( nf90_def_dim(fid,'bound',kmax+1,bid) )
  call nf90_check( nf90_def_dim(fid,'hydro',hmax,  hid) )
  !call nf90_check( nf90_def_dim(fid,'time',NF90_UNLIMITED,tid) )
  write(*,*) 'output file create: '//trim(dir)//'_MIROCinput.nc'

  ! ===== lat / lon ====
  write(*,*) '+ lat/lon grid'
  open(unit=10,file='../MIROC_inputs/'//trim(dir)//'/mlon.gt',&
       & form='unformatted',access='sequential',&
       & status='old',action='read',convert='big_endian')
  read(10) header
  ikind = header(29) ; read(header(31),'(I16)') igt
  jkind = header(32) ; read(header(34),'(I16)') jgt
  kkind = header(35) ; read(header(37),'(I16)') kgt
  !tkind = header(26) ; read(header(28),'(I16)') tgt
  if (igt /= imax .or. jgt /= jmax .or. kgt /= 1) then
     write(*,*) 'Ngrid error: ',igt,jgt,kgt ; stop
  end if
  if (trim(header(38)) == 'UR4') then
     read(10) var4s(1:imax,1:jmax)
     clon(:) = dble(var4s(:,1))
  else if (trim(header(38)) == 'UR8') then
     read(10) var8s(1:imax,1:jmax)
     clon(:) = var8s(:,1)
  end if
  close(10)
  !-----
  open(unit=10,file='../MIROC_inputs/'//trim(dir)//'/mlat.gt',&
       & form='unformatted',access='sequential',&
       & status='old',action='read',convert='big_endian')
  read(10) header
  ikind = header(29) ; read(header(31),'(I16)') igt
  jkind = header(32) ; read(header(34),'(I16)') jgt
  kkind = header(35) ; read(header(37),'(I16)') kgt
  !tkind = header(26) ; read(header(28),'(I16)') tgt
  if (igt /= imax .or. jgt /= jmax .or. kgt /= 1) then
     write(*,*) 'Ngrid error: ',igt,jgt,kgt ; stop
  end if
  if (trim(header(38)) == 'UR4') then
     read(10) var4s(1:imax,1:jmax)
     clat(:) = dble(var4s(1,jmax:1:-1))
  else if (trim(header(38)) == 'UR8') then
     read(10) var8s(1:imax,1:jmax)
     clat(:) = var8s(1,jmax:1:-1)
  end if
  close(10)
  !-----
  call nf90_check( nf90_def_var(fid,'lon',NF90_DOUBLE,(/xid/),vid) )
  call nf90_check( nf90_put_var(fid,vid,clon,start=(/1/), count=(/imax/) ) )
  call nf90_check( nf90_def_var(fid,'lat',NF90_DOUBLE,(/yid/),vid) )
  call nf90_check( nf90_put_var(fid,vid,clat,start=(/1/), count=(/jmax/) ) )
  
  
  ! ===== variables 2D =====
  do n=1,Nv2D
     write(*,*) 'Now ... '//vname2D(1,n)//' -> '//vname2D(2,n)
     open(unit=10,file='../MIROC_inputs/'//trim(dir)//'/'//trim(vname2D(1,n))//'.gt',&
          & form='unformatted',access='sequential',&
          & status='old',action='read',convert='big_endian')
     read(10) header
     ikind = header(29) ; read(header(31),'(I16)') igt
     jkind = header(32) ; read(header(34),'(I16)') jgt
     kkind = header(35) ; read(header(37),'(I16)') kgt
     !tkind = header(26) ; read(header(28),'(I16)') tgt
     if (igt /= imax .or. jgt /= jmax .or. kgt /= 1) then
        write(*,*) 'Ngrid error: ',igt,jgt,kgt ; stop
     end if
     rewind(10)

     do t=1,tmax
        read(10) header

        if (trim(header(38)) == 'UR4') then
           read(10) var4s(1:imax,1:jmax)
           out2D(:,:,t) = dble(var4s(:,jmax:1:-1))
        else if (trim(header(38)) == 'UR8') then
           read(10) var8s(1:imax,1:jmax)
           out2D(:,:,t) = var8s(:,jmax:1:-1)
        end if
     end do
     close(10)

     call nf90_check( nf90_def_var(fid,trim(vname2D(2,n)),NF90_DOUBLE,(/xid,yid/),vid) )
     call nf90_check( nf90_put_var(fid,vid,out2D(:,:,tmax),start=(/1,1/), count=(/imax,jmax/) ) )
  end do
  
  ! ===== variables 3D =====
  do n=1,Nv3D
     write(*,*) 'Now ... '//vname3D(1,n)//' -> '//vname3D(2,n)
     open(unit=10,file='../MIROC_inputs/'//trim(dir)//'/'//trim(vname3D(1,n))//'.gt',&
          & form='unformatted',access='sequential',&
          & status='old',action='read',convert='big_endian')
     read(10) header
     ikind = header(29) ; read(header(31),'(I16)') igt
     jkind = header(32) ; read(header(34),'(I16)') jgt
     kkind = header(35) ; read(header(37),'(I16)') kgt
     !tkind = header(26) ; read(header(28),'(I16)') tgt
     if (igt /= imax .or. jgt /= jmax .or. kgt /= kmax) then
        write(*,*) 'Ngrid error: ',igt,jgt,kgt ; stop
     end if
     rewind(10)

     do t=1,tmax
        read(10) header

        if (trim(header(38)) == 'UR4') then
           read(10) var4l(1:imax,1:jmax,1:kmax)
           out3D(:,:,:,t) = dble(var4l(:,jmax:1:-1,:))
        else if (trim(header(38)) == 'UR8') then
           read(10) var8l(1:imax,1:jmax,1:kmax)
           out3D(:,:,:,t) = var8l(:,jmax:1:-1,:)
        end if
     end do
     close(10)
     
     !if (trim(vname3D(2,n)) == 'height' .or. trim(vname3D(2,n)) == 'height_half') then
     !   out3D = out3D/1000.  ! m -> km
     !end if
     
     call nf90_check( nf90_def_var(fid,trim(vname3D(2,n)),NF90_DOUBLE,(/xid,yid,zid/),vid) )
     call nf90_check( nf90_put_var(fid,vid,out3D(:,:,:,tmax),start=(/1,1,1/),count=(/imax,jmax,kmax/)) )
  end do
  
  ! ===== variables 3B =====
  do n=1,Nv3B
     write(*,*) 'Now ... '//vname3B(1,n)//' -> '//vname3B(2,n)
     open(unit=10,file='../MIROC_inputs/'//trim(dir)//'/'//trim(vname3B(1,n))//'.gt',&
          & form='unformatted',access='sequential',&
          & status='old',action='read',convert='big_endian')
     read(10) header
     ikind = header(29) ; read(header(31),'(I16)') igt
     jkind = header(32) ; read(header(34),'(I16)') jgt
     kkind = header(35) ; read(header(37),'(I16)') kgt
     !tkind = header(26) ; read(header(28),'(I16)') tgt
     if (igt /= imax .or. jgt /= jmax .or. kgt /= bmax) then
        write(*,*) 'Ngrid error: ',igt,jgt,kgt ; stop
     end if
     rewind(10)

     do t=1,tmax
        read(10) header

        if (trim(header(38)) == 'UR4') then
           read(10) var4b(1:imax,1:jmax,1:bmax)
           out3B(:,:,:,t) = dble(var4b(:,jmax:1:-1,:))
        else if (trim(header(38)) == 'UR8') then
           read(10) var8b(1:imax,1:jmax,1:bmax)
           out3B(:,:,:,t) = var8b(:,jmax:1:-1,:)
        end if
     end do
     close(10)

     !call nf90_check( nf90_def_var(fid,trim(vname3B(2,n)),NF90_DOUBLE,(/xid,yid,bid,tid/),vid) )
     !call nf90_check( nf90_put_var(fid,vid,out3B,start=(/1,1,1,1/), count=(/imax,jmax,bmax,tmax/) ) )

     call nf90_check( nf90_def_var(fid,trim(vname3B(2,n)),NF90_DOUBLE,(/xid,yid,zid/),vid) )
     if (trim(vname3B(2,n)) == 'phalf' ) then
        call nf90_check( nf90_put_var(fid,vid,out3B(:,:,1:kmax,tmax), &
             &           start=(/1,1,1/), count=(/imax,jmax,kmax/) ) )
     else if (trim(vname3B(2,n)) == 'gcumf' ) then
        call nf90_check( nf90_put_var(fid,vid,(out3B(:,:,1:kmax,tmax)+out3B(:,:,2:kmax+1,tmax))/2., &
             &           start=(/1,1,1/), count=(/imax,jmax,kmax/) ) )
     end if
  end do

  ! ===== Reff =====
  write(*,*) 'Now ... Reff 1-9'
  do n=1,hmax
     write(tp,'(I1)') n
     open(unit=10,file='../MIROC_inputs/'//trim(dir)//'/Reff'//tp//'.gt',&
          & form='unformatted',access='sequential',&
          & status='old',action='read',convert='big_endian')
     read(10) header
     ikind = header(29) ; read(header(31),'(I16)') igt
     jkind = header(32) ; read(header(34),'(I16)') jgt
     kkind = header(35) ; read(header(37),'(I16)') kgt
     !tkind = header(26) ; read(header(28),'(I16)') tgt
     if (igt /= imax .or. jgt /= jmax .or. kgt /= kmax) then
        write(*,*) 'Ngrid error: ',igt,jgt,kgt ; stop
     end if
     rewind(10)

     do t=1,tmax
        read(10) header

        if (trim(header(38)) == 'UR4') then
           read(10) var4l(1:imax,1:jmax,1:kmax)
           ReffIN(:,:,:,n,t) = dble(var4l(:,jmax:1:-1,:))
        else if (trim(header(38)) == 'UR8') then
           read(10) var8l(1:imax,1:jmax,1:kmax)
           ReffIN(:,:,:,n,t) = var8l(:,jmax:1:-1,:)
        end if
     end do
     close(10)
  end do
          
  call nf90_check( nf90_def_var(fid,'Reff',NF90_DOUBLE,(/xid,yid,zid,hid/),vid) )
  call nf90_check( nf90_put_var(fid,vid,ReffIN(:,:,:,:,tmax), &
       &           start=(/1,1,1,1/),count=(/imax,jmax,kmax,hmax/)) )

  
  ! END.
  call execute_command_line('chmod +x ../MIROC_inputs/'//trim(dir)//'_MIROCinput.nc')

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
