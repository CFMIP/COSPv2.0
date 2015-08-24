module mod_cosp_error
  use cosp_kinds, ONLY: wp
contains
  ! ######################################################################################
  ! Subroutine errorMessage_print
  ! ######################################################################################
  subroutine errorMessage(message)
    ! Inputs
    character(len=*),intent(in) :: message

    print*,message
  end subroutine errorMessage
  
  ! ######################################################################################
  ! END MODULE
  ! ######################################################################################
end module mod_cosp_error
