module auxiliary_routines
  implicit none
  private

  public :: stop_program, get_current_time, write_time

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine stop_program(error_message, ioerror)
   use iso_fortran_env, only: error_unit

     character(len=*), intent(in), optional :: error_message
     integer, intent(in), optional :: ioerror

     if (present(error_message)) then
        write(unit=error_unit, fmt='(1x,a)') error_message
     end if

     if (present(ioerror)) then
        write(unit=error_unit, fmt='(1x,a,1x,g0)') &
                                        '[ERROR] Runtime Error Code :', ioerror
     end if

     write(unit=error_unit, fmt='(1x,a)') '[ERROR] Program stopped'
     stop

  end subroutine stop_program

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function get_current_time() result(current_time)
    real :: current_time

    integer, dimension(8) :: time_values

    call date_and_time(values=time_values)
    current_time = 86400.0 * time_values(3) + 3600.0 * time_values(5) + &
                   60.0 * time_values(6) + time_values(7) + &
                   0.001 * time_values(8)

  end function get_current_time

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine write_time(elapsed_time)
    use iso_fortran_env, only: output_unit
    real, intent(in) :: elapsed_time

    write(output_unit,'(1x,3(a,g0))') &
                            '[INFO] Elasped time (h:min:sec.msec) is ', &
                            int(elapsed_time / 3600.0), &
                            ':', int(mod(elapsed_time / 60.0, 60.0)), &
                            ':', mod(elapsed_time, 60.0)

  end subroutine write_time

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module auxiliary_routines