program test_isccp_cloud_types
  USE COSP_KINDS, ONLY: wp

  ! *****************************COPYRIGHT****************************
  ! (c) British Crown Copyright 2009, the Met Office.
  ! All rights reserved.
  ! $Revision: 23 $, $Date: 2011-03-31 07:41:37 -0600 (Thu, 31 Mar 2011) $
  ! $URL: https://cfmip-obs-sim.googlecode.com/svn/devel/branches/dustinswales/icarus-scops-4.1-bsd/test_isccp_cloud_types.f $
  ! 
  ! Redistribution and use in source and binary forms, with or without 
  ! modification, are permitted provided that the
  ! following conditions are met:
  ! 
  !     * Redistributions of source code must retain the above 
  !       copyright  notice, this list of conditions and the following 
  !       disclaimer.
  !     * Redistributions in binary form must reproduce the above 
  !       copyright notice, this list of conditions and the following 
  !       disclaimer in the documentation and/or other materials 
  !       provided with the distribution.
  !     * Neither the name of the Met Office nor the names of its 
  !       contributors may be used to endorse or promote products
  !       derived from this software without specific prior written 
  !       permission.
  ! 
  ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
  ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
  ! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
  ! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
  ! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
  ! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
  ! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
  ! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
  ! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
  ! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
  ! 
  ! *****************************COPYRIGHT*******************************
  ! *****************************COPYRIGHT*******************************
  ! *****************************COPYRIGHT*******************************
  implicit none 
  real(wp),parameter :: emsfc_lw = 0.99     ! LW emissivity of surface at 10.5 microns
  integer, parameter :: npoints  = 24,    & ! Number of model points in horizontal.
			nlev     = 6,     & ! Number of model levels.
                        ncol     = 40,    & ! Number of columns grid box is divided into.
                        ntau     = 7,     &
                        npres    = 7
  integer                                :: i,j,k,itau,top_height,top_height_direction, &
    					    overlap,ibox,debug,debugcol,ncolprint
  integer, dimension(npoints)            :: seed,sunlit
  real(wp)                                   :: lastboxptop
  real(wp),    dimension(npoints)            :: skt,totalcldarea,meantaucld,meanalbedocld, &
                                            meanptop,meantb,meantbclr
  real(wp),    dimension(npoints,nlev)       :: pfull,at,qv,cc,conv,dtau_s,dtau_c,dem_s,dem_c
  real(wp),    dimension(npoints,nlev+1)     :: phalf   
  real(wp),    dimension(npoints,ncol)       :: boxtau,boxptop,tb
  real(wp),    dimension(npoints,ntau,npres) :: fq_isccp
  real(wp),    dimension(npoints,ncol,nlev)  :: frac_out
  
  ! Get user supplied configuration
  write(6,*) 'Please enter top_height: '
  write(6,*) ' 1 = adjusted top IR+VIS '
  write(6,*) ' 2 = real top '
  write(6,*) ' 3 = adjusted top IR only'
  read(5,*) top_height
  write(6,'(i2)') top_height
  write(6,*) 'Please enter overlap type (1=max,2=rand,3=max/rand)'
  read(5,*) overlap
  write(6,'(i2)') overlap
  write(6,*) 'Please enter top_height_direction (1=original method,2=new method)'
  read(5,*) top_height_direction
  write(6,'(i2)') top_height_direction
  
  ncolprint=0
  do j=1,npoints 
     
     ! Open sample data files
     open(1,file='data/input/input.data',form='formatted')
     open(2,file='data/input/input.data.halved',form='formatted')
     
     if (j .lt. 13) then
        sunlit(j)=1
     else
        sunlit(j)=0
     endif
     
     k=mod(j-1,12) + 1  
     if (k.eq.2) then 
        read(1,*) skt(j)
        read(2,*) skt(j)
     else
        read(2,*) skt(j)
        read(1,*) skt(j)
     endif
     !        write(6,'(12X,3(f12.3,4X))') skt(j)
     
     if (k.eq.3) then 
        read(1,*) (pfull(j,i),i=1,nlev)
        read(2,*) (pfull(j,i),i=1,nlev)
     else
        read(2,*) (pfull(j,i),i=1,nlev)
        read(1,*) (pfull(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (pfull(j,i),i=1,nlev)
     
     if (k.eq.4) then 
        read(1,*) (at(j,i),i=1,nlev)
        read(2,*) (at(j,i),i=1,nlev)
     else
        read(2,*) (at(j,i),i=1,nlev)
        read(1,*) (at(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (at(j,i),i=1,nlev)
     
     if (k.eq.5) then 
        read(1,*) (qv(j,i),i=1,nlev)
        read(2,*) (qv(j,i),i=1,nlev)
     else
        read(2,*) (qv(j,i),i=1,nlev)
        read(1,*) (qv(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (qv(j,i),i=1,nlev)
     
     if (k.eq.6) then 
        read(1,*) (cc(j,i),i=1,nlev)
        read(2,*) (cc(j,i),i=1,nlev)
     else
        read(2,*) (cc(j,i),i=1,nlev)
        read(1,*) (cc(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (cc(j,i),i=1,nlev)
     
     if (k.eq.7) then 
        read(1,*) (conv(j,i),i=1,nlev)
        read(2,*) (conv(j,i),i=1,nlev)
     else
        read(2,*) (conv(j,i),i=1,nlev)
        read(1,*) (conv(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (conv(j,i),i=1,nlev)
     
     if (k.eq.8) then 
        read(1,*) (dtau_s(j,i),i=1,nlev)
        read(2,*) (dtau_s(j,i),i=1,nlev)
     else
        read(2,*) (dtau_s(j,i),i=1,nlev)
        read(1,*) (dtau_s(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (dtau_s(j,i),i=1,nlev)
     
     if (k.eq.9) then 
        read(1,*) (dtau_c(j,i),i=1,nlev)
        read(2,*) (dtau_c(j,i),i=1,nlev)
     else
        read(2,*) (dtau_c(j,i),i=1,nlev)
        read(1,*) (dtau_c(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (dtau_c(j,i),i=1,nlev)
     
     if (k.eq.10) then 
        read(1,*) (dem_s(j,i),i=1,nlev)
        read(2,*) (dem_s(j,i),i=1,nlev)
     else
        read(2,*) (dem_s(j,i),i=1,nlev)
        read(1,*) (dem_s(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (dem_s(j,i),i=1,nlev)
     
     if (k.eq.11) then 
        read(1,*) (dem_c(j,i),i=1,nlev)
        read(2,*) (dem_c(j,i),i=1,nlev)
     else
        read(2,*) (dem_c(j,i),i=1,nlev)
        read(1,*) (dem_c(j,i),i=1,nlev)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (dem_c(j,i),i=1,nlev)
     
     if (k.eq.12) then 
        read(1,*) (phalf(j,i),i=1,nlev+1)
        read(2,*) (phalf(j,i),i=1,nlev+1)
     else
        read(2,*) (phalf(j,i),i=1,nlev+1)
        read(1,*) (phalf(j,i),i=1,nlev+1)
     endif
     !        write(6,'(12X,3(f12.3,4X))') (phalf(j,i),i=1,nlev+1)
     
     ! Close files
     close(1)
     close(2)
     
     seed(j)=50  ! Note you want this to vary in a GCM - see the README file
  enddo

  debug=1     !  These should be zero when running for real
  debugcol=1  !  1 is required in the test suite, but you might want
              !  something like 1000 if debugging in a GCM

  call ISCCP_CLOUD_TYPES(debug, debugcol, npoints, sunlit, nlev, ncol, seed, pfull, &
                         phalf, qv, cc, conv, dtau_s, dtau_c, top_height,           &
                         top_height_direction, overlap, frac_out, skt, emsfc_lw,    &
                         at, dem_s, dem_c, fq_isccp, totalcldarea, meanptop,        &
                         meantaucld, meanalbedocld, meantb, meantbclr, boxtau, boxptop)
  write(6,*) 'START'
  do j=1,npoints 
     write (6,'(a5)') "j="
     write (6,'(I3)') j
     write(6,'(12X,3(f12.3,4X))') skt(j)
     write(6,'(12X,3(f12.3,4X))') (pfull(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (at(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (qv(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (cc(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (conv(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (dtau_s(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (dtau_c(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (dem_s(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (dem_c(j,i),i=1,nlev)
     write(6,'(12X,3(f12.3,4X))') (phalf(j,i),i=1,nlev+1)
     
     ! Print results
     if (top_height .eq. 1) then
        print *, '      ISCCP PCTAU DIAGRAM WITH IR-VIS CLOUD TOP HEIGHT ADJUSTMENT'
     else if (top_height .eq. 2) then
        print *, '      ISCCP PCTAU DIAGRAM WITHOUT CLOUD TOP HEIGHT ADJUSTMENT'
     else if (top_height .eq. 3) then
        print *, '      ISCCP PCTAU DIAGRAM WITH IR ONLY CLOUD TOP HEIGHT ADJUSTMENT'
     end if
     print *, ' '
     print *, '                            TAU '
     print *, ' '
     print *, '              taumin    1.3     3.6     9.4     23.     60.'
     print *, ' '
     if (totalcldarea(j) .ge. 0.) then
        print *, '     50 '
        write(6,'(12X,7(f12.3,4X))') (fq_isccp(j,itau,1),itau=1,7)
        print *, '    180 '
        write(6,'(12X,7(f12.3,4X))') (fq_isccp(j,itau,2),itau=1,7)
        print *, '    310 '
        write(6,'(12X,7(f12.3,4X))') (fq_isccp(j,itau,3),itau=1,7)
        print *, 'PC  440 '
        write(6,'(12X,7(f12.3,4X))') (fq_isccp(j,itau,4),itau=1,7)
        print *, '    560 '
        write(6,'(12X,7(f12.3,4X))') (fq_isccp(j,itau,5),itau=1,7)
        print *, '    680 '
        write(6,'(12X,7(f12.3,4X))') (fq_isccp(j,itau,6),itau=1,7)
        print *, '    800 '
        write(6,'(12X,7(f12.3,4X))') (fq_isccp(j,itau,7),itau=1,7)
        print *, '    sfc '
        print *, ' '
     else
        print *, '     50 '
        write(6,'(12X,7(e12.3,4X))') (fq_isccp(j,itau,1),itau=1,7)
        print *, '    180 '
        write(6,'(12X,7(e12.3,4X))') (fq_isccp(j,itau,2),itau=1,7)
        print *, '    310 '
        write(6,'(12X,7(e12.3,4X))') (fq_isccp(j,itau,3),itau=1,7)
        print *, 'PC  440 '
        write(6,'(12X,7(e12.3,4X))') (fq_isccp(j,itau,4),itau=1,7)
        print *, '    560 '
        write(6,'(12X,7(e12.3,4X))') (fq_isccp(j,itau,5),itau=1,7)
        print *, '    680 '
        write(6,'(12X,7(e12.3,4X))') (fq_isccp(j,itau,6),itau=1,7)
        print *, '    800 '
        write(6,'(12X,7(e12.3,4X))') (fq_isccp(j,itau,7),itau=1,7)
        print *, '    sfc '
        print *, ' '
     endif
     if (top_height .ne. 2. .and. top_height_direction .eq. 1) then
        print *, 'TOP HEIGHT IDENTIFIED AS LOWEST ALTITUDE OR HIGHEST'
        print *, 'PRESSURE WITH MATCHING CLOUD TOP TEMPERATURE'
     end if
     if (top_height .ne. 2. .and. top_height_direction .eq. 2) then
        print *, 'TOP HEIGHT IDENTIFIED AS HIGHEST ALTITUDE OR LOWEST'
        print *, 'PRESSURE WITH MATCHING CLOUD TOP TEMPERATURE'
     end if
     print *, ' '
     print *, ' '
     if (totalcldarea(j) .ge. 0.) then
        write (6,'(a,f10.3)') 'totalcldarea = ', totalcldarea(j)
     else
        write (6,'(a,e10.3)') 'totalcldarea = ', totalcldarea(j)
     end if
     if (totalcldarea(j) .gt. 0.) then
        write (6,'(a,f10.3)') 'meanptop     = ', meanptop(j)
        write (6,'(a,f10.3)') 'meantaucld   = ', meantaucld(j)
        write (6,'(a,f10.3)') 'meanalbedocld   = ', meanalbedocld(j)
     else
        write (6,'(a,e10.3)') 'meanptop     = ', meanptop(j)
        write (6,'(a,e10.3)') 'meantaucld   = ', meantaucld(j)
        write (6,'(a,e10.3)') 'meanalbedocld   = ', meanalbedocld(j)
     end if
     if (top_height .ne. 2.) then
        write (6,'(a,f10.3)') 'meantb = ', meantb(j)
        write (6,'(a,f10.3)') 'meantbclr = ', meantbclr(j)
     end if
     print *, '  box #         ptop           tau  '
     print *, ' -------       -------       -------'
     
     lastboxptop=-1
     do ibox = 1, ncol
        if (boxptop(j,ibox) .ne. lastboxptop) then
           if (boxptop(j,ibox) .gt. 0. .or. boxtau(j,ibox) .gt. 0.) then
              write(6,'(2X,i4,a,9X,f7.2,7X,f7.1)') ibox,' - ',boxptop(j,ibox),boxtau(j,ibox)
           else
              write(6,'(2X,i4,a,9X,e10.3,7X,e10.3)') ibox,' - ',boxptop(j,ibox),boxtau(j,ibox)
           end if
        end if
        lastboxptop=boxptop(j,ibox)
     enddo
     print *, ' '
     print *, ' '
     print *, ' '
  enddo

  ! END PROGRAM
end program test_isccp_cloud_types
