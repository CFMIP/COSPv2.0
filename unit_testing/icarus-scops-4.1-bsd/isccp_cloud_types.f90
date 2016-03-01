SUBROUTINE ISCCP_CLOUD_TYPES(debug,debugcol,npoints,sunlit,nlev,ncol,seed,pfull,phalf, &
                             qv,cc,conv,dtau_s,dtau_c,top_height,top_height_direction, &
                             overlap,frac_out,skt,emsfc_lw,at,dem_s,dem_c,fq_isccp,    &
                             totalcldarea,meanptop,meantaucld,meanalbedocld, meantb,   &
                             meantbclr,boxtau,boxptop)
  USE COSP_KINDS, ONLY: wp
  USE mod_icarus, ONLY: icarus
  USE mod_scops,  ONLY: scops
  USE mod_rng,    ONLY: rng_state, init_rng, get_rng

  !$Id: isccp_cloud_types.f,v 4.0 2009/03/06 11:05:11 hadmw Exp $
  
  ! *****************************COPYRIGHT****************************
  ! (c) British Crown Copyright 2009, the Met Office.
  ! All rights reserved.
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

      ! The maximum number of levels and columns 
      INTEGER, parameter :: ncolprint = 0

      ! INPUT
      INTEGER, intent(in)                               :: npoints,nlev,ncol,overlap,top_height,&
                                                           top_height_direction,debug,debugcol    
      INTEGER, intent(in),  dimension(npoints)          :: sunlit,seed   
      REAL(WP),intent(in)                               :: emsfc_lw                                      
      REAL(WP),intent(in), dimension(npoints)           :: skt
      REAL(WP),intent(in), dimension(npoints,nlev)      :: pfull,qv,cc,conv, dtau_s,dtau_c,dem_s,dem_c,at
      REAL(WP),intent(in), dimension(npoints,nlev+1)    :: phalf
      REAL(WP),intent(inout), dimension(npoints,ncol,nlev) :: frac_out  

      ! OUTPUT
      REAL(WP),intent(out), dimension(npoints,7,7)  :: fq_isccp 
      REAL(WP),intent(out), dimension(npoints)      ::  totalcldarea, meanptop, meantaucld,meanalbedocld,meantb,meantbclr 
      REAL(WP),intent(out), dimension(npoints,ncol) :: boxtau,boxptop  
                              
      ! Local variables
      INTEGER                            :: i,j,ilev,ibox
      INTEGER, dimension(nlev,ncol)      :: acc
      INTEGER, dimension(npoints,ncol)   :: levmatch
      character*10                       :: ftn09
      type(rng_state),dimension(npoints) :: rngs  ! Seeds for random number generator

      ! Parameters
      REAL(WP),parameter                 :: isccp_taumin=0.3
      character,parameter,dimension(6)   :: cchar=(/' ','-','1','+','I','+'/)
      character*1,parameter,dimension(6) :: cchar_realtops=(/ ' ',' ','1','1','I','I'/)

      call init_rng(rngs, seed)  
      CALL SCOPS(npoints,nlev,ncol,rngs,cc,conv,overlap,frac_out,ncolprint)
      
      CALL ICARUS(debug,debugcol,npoints,sunlit,nlev,ncol,pfull,phalf,qv,&
                  cc,conv,dtau_s,dtau_c,top_height,top_height_direction, &
                  frac_out,skt,emsfc_lw,at,dem_s,dem_c,fq_isccp,         &
                  totalcldarea,meanptop,meantaucld, meanalbedocld,meantb,&
                  meantbclr,boxtau,boxptop,levmatch)

      if ( debug.ne.0 ) then
         j=1
         write(6,'(a10)') 'j='
         write(6,'(8I10)') j
         write(6,'(a10)') 'debug='
         write(6,'(8I10)') debug
         write(6,'(a10)') 'debugcol='
         write(6,'(8I10)') debugcol
         write(6,'(a10)') 'npoints='
         write(6,'(8I10)') npoints
         write(6,'(a10)') 'nlev='
         write(6,'(8I10)') nlev
         write(6,'(a10)') 'ncol='
         write(6,'(8I10)') ncol
         write(6,'(a11)') 'top_height='
         write(6,'(8I10)') top_height
         write(6,'(a21)') 'top_height_direction='
         write(6,'(8I10)') top_height_direction
         write(6,'(a10)') 'overlap='
         write(6,'(8I10)') overlap
         write(6,'(a10)') 'emsfc_lw='
         write(6,'(8f10.2)') emsfc_lw
         do j=1,npoints,debug
            write(6,'(a10)') 'j='
            write(6,'(8I10)') j
            write(6,'(a10)') 'sunlit='
            write(6,'(8I10)') sunlit(j)
            write(6,'(a10)') 'pfull='
            write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
            write(6,'(a10)') 'phalf='
            write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
            write(6,'(a10)') 'qv='
            write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
            write(6,'(a10)') 'cc='
            write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
            write(6,'(a10)') 'conv='
            write(6,'(8f10.2)') (conv(j,i),i=1,nlev)
            write(6,'(a10)') 'dtau_s='
            write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
            write(6,'(a10)') 'dtau_c='
            write(6,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
            write(6,'(a10)') 'skt='
            write(6,'(8f10.2)') skt(j)
            write(6,'(a10)') 'at='
            write(6,'(8f10.2)') (at(j,i),i=1,nlev)
            write(6,'(a10)') 'dem_s='
            write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
            write(6,'(a10)') 'dem_c='
            write(6,'(8f10.3)') (dem_c(j,i),i=1,nlev)
         enddo
	     do j=1,npoints,debugcol
         
            ! Produce character output
            do ilev=1,nlev
               acc(ilev,1:ncol)=frac_out(j,1:ncol,ilev)*2
               where(levmatch(j,1:ncol) .eq. ilev) acc(ilev,1:ncol)=acc(ilev,1:ncol)+1
            enddo
         
            write(ftn09,11) j
11             format('ftn09.',i4.4)
            open(9, FILE=ftn09, FORM='FORMATTED')
         
            write(9,'(a1)') ' '
            write(9,'(10i5)') (ilev,ilev=5,nlev,5)
            write(9,'(a1)') ' '

            do ibox=1,ncol
               write(9,'(40(a1),1x,40(a1))') &
                    (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev),&
                    (cchar(acc(ilev,ibox)+1),ilev=1,nlev) 
            end do
            close(9)
         enddo
      endif
      
      return
    end SUBROUTINE ISCCP_CLOUD_TYPES

