module array_lib
  USE COSP_KINDS, ONLY: wp
  implicit none
contains

  ! ############################################################################
  !                               function INFIND
  ! ############################################################################
  function infind(list,val)
    implicit none
    ! ##########################################################################
    ! Purpose:
    !   Finds the index of an array that is closest to a value, plus the
    !   difference between the value found and the value specified
    !
    ! Inputs:
    !   [list]   an array of sequential values
    !   [val]    a value to locate
    ! Optional input:
    !   [sort]   set to 1 if [list] is in unknown/non-sequential order
    !
    ! Returns:
    !   index of [list] that is closest to [val]
    !
    ! Optional output:
    !   [dist]   set to variable containing [list([result])] - [val]
    !
    ! Requires:
    !   mrgrnk library
    !
    ! Created:
    !   10/16/03  John Haynes (haynes@atmos.colostate.edu)
    ! Modified:
    !   01/31/06  IDL to Fortran 90
    !   01/05/15  Removed optional arguements and corresponding dependencies Dustin Swales (dustin.swales@noaa.gov)
    !
    ! ##########################################################################

    ! INPUTS
    real(wp), dimension(:), intent(in) :: &
         list   ! An array of sequential values
    real(wp), intent(in) :: &
         val    ! A value to locate
    ! OUTPUTS
    integer :: &
         infind ! Index of [list] that is closest to [val]

    ! Internal Variables
    real(wp), dimension(size(list)) :: lists
    integer :: nlist, result, tmp(1), sort_list
    integer, dimension(size(list)) :: mask, idx
    
    sort_list = 0
    
    nlist = size(list)
    lists = list
    
    if (val >= lists(nlist)) then
       result = nlist
    else if (val <= lists(1)) then
       result = 1
    else
       mask(:) = 0
       where (lists < val) mask = 1
       tmp = minloc(mask,1)
       if (abs(lists(tmp(1)-1)-val) < abs(lists(tmp(1))-val)) then
          result = tmp(1) - 1
      else
         result = tmp(1)
      endif
   endif
   infind = result
 end function infind

end module array_lib
