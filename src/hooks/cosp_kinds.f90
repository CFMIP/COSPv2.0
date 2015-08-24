! ##########################################################################
!
! Copyright 2014, Regents of the University of Colorado. All right reserved.
! Use and duplication is permitted under the terms of the
! BSD 3-Clause License: https://opensource.org/licenses/BSD-3-Clause
!
! ##########################################################################

MODULE cosp_kinds
  implicit none
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND( 6, 37)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: wp = sp

END MODULE cosp_kinds
