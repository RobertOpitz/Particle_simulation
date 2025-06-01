module auxiliary_routines
  use iso_fortran_env, only: real32, real64, output_unit, input_unit, error_unit

  implicit none
  private

  public :: stop_program
  public :: get_current_time
  public :: write_time
  public :: check_file
  public :: get_value

  interface get_value
     module procedure get_single_precision_real_value
     module procedure get_double_precision_real_value
     module procedure get_integer_value
     module procedure get_logical_value
  end interface

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine stop_program(error_message, ioerror)
  !subroutine stop_program(error_message, iomsg, ioerror)
  !
  !
  !
    character(len=*), intent(in), optional :: error_message
    !character(len=*), intent(in), optional :: iomsg
    integer, intent(in), optional :: ioerror
    
    if (present(error_message)) then
      write(unit=error_unit, fmt='(1x,a,1x,a)') '[ERROR]', error_message
    end if

    !if (present(iomsg)) then
    !  write(unit=error_unit, fmt='(1x,a,1x,a)') '[ERROR] Error message:', iomsg
    !end if
    
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
    real, intent(in) :: elapsed_time

    write(output_unit,'(1x,3(a,g0))') &
                            '[INFO] Elasped time (h:min:sec.msec) is ', &
                            int(elapsed_time / 3600.0), &
                            ':', int(mod(elapsed_time / 60.0, 60.0)), &
                            ':', mod(elapsed_time, 60.0)

  end subroutine write_time

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_file( file_name, overwrite_existing_file )
  !
  ! Purpose:
  !	Checks an Inputfile for it's Existence.
  !
  ! Record revisions:
  !	Date          Programmer      Desciption of Change
  !	====          ==========      ==================== 
  !	27.04.2008		Robert Opitz		Orignal Code
  !	03.09.2008		Robert Opitz		Extentension for "overwrite"
  !	06.12.2009		Robert Opitz		Structured programming, iso_fortran_env
  ! 24.05.2025    Robert Opitz    Minor adaptions
  !
    character( len=* ), intent( in out ) :: file_name
    logical, intent( in ), optional :: overwrite_existing_file
    
    logical :: overwrite
    logical :: file_exists
    logical :: yes_no
    
    intrinsic :: present
    
    ! set flag for overwritting
    if ( present( overwrite_existing_file ) .eqv. .true. ) then
      overwrite = overwrite_existing_file
    else
      overwrite = .true.
    end if
    
    !Check File
    inquire( file = file_name, exist = file_exists )
         
    if ( file_exists .eqv. .true. ) then
      if ( overwrite .eqv. .false. ) then
        ! do not overwrite existing file, if user want it so
        write( unit = output_unit, fmt='(1x,a,a,a)' ) 'File "', &
                                                      trim( file_name ), &
                                                      '" does already exist'
        call get_value( 'Overwrite File? (y/n)', yes_no )
        if ( yes_no .eqv. .false. ) call stop_program('Program stopped by user')
      end if
    else ! file does not exist
      if ( overwrite .eqv. .true. ) then
        ! trying to overwrite a non-existing file
        call stop_program( 'File "'//trim( file_name )//'" does not exist' )
      end if
    end if
    
  end subroutine check_file

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
  subroutine get_single_precision_real_value( lyrics, real_variable, &
    relational_operator, reference_value )
  !
  ! Purpose:
  !	Reads an real-Value. Optional: Compares the Value with a Reference.
  !	The variables "real_variable" and "reference_value" are both single precision.
  !
  ! Record revisions:
  !	Date			  Programmer    Desciption of Change
  !	====        ==========		====================
  !	21.05.2008  Robert Opitz  Orignal Code
  !	06.12.2009	Robert Opitz  iso_fortran_env
  !
    character( len=* ), intent( in ) :: lyrics
    real( kind=real32 ), intent( out ) :: real_variable
    character( len=* ), intent( in ), optional :: relational_operator
    real( kind=real32 ), intent( in ), optional :: reference_value

    real( kind=real32 ) :: input
    integer :: ioerror
    logical :: ok
    
    intrinsic :: present
    
    loop: do
      write( unit=output_unit, fmt='(1x,a)', advance='no' ) lyrics
      read( unit=input_unit, fmt=*, iostat=ioerror ) input
    
      if ( ioerror /= 0 ) then
        write( unit=output_unit, fmt='(1x,a)' ) 'Your Input was incorrect'
        cycle loop
      end if
        
      if ( present( relational_operator ) .and. present( reference_value ) ) then
        call check_value( ok )
        if ( .not. ok ) cycle loop
      else if ( present( relational_operator ) .or. present( reference_value ) ) then
        call stop_program( 'get_value error: relational_operator or reference_value missing' )
      end if
      
      real_variable = input
      exit loop
      
    end do loop
    
    contains
    
      subroutine check_value( ok )
      !
      !
      !
        logical, intent( out ) :: ok

        character(len=9), parameter :: fmt = '(1x,a,g0)'

        ok = .true.

        select case ( relational_operator )
          case( '<' )
            if ( input >= reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be smaller than : ', reference_value
              ok = .false.
            end if
          case( '<=' )
            if ( input > reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be smaller or equal than : ', reference_value
              ok = .false.
            end if
          case( '>' )
            if ( input <= reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be bigger than : ', reference_value
              ok = .false.
            end if
          case( '>=' )
            if ( input < reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be bigger or equal than : ', reference_value
              ok = .false.
            end if
          case( '/=' )
            if ( input == reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be unequal to : ', reference_value
              ok = .false.
            end if
          case default
            call stop_program( 'get_value error: relational_operator missing' )
        end select
        
      end subroutine check_value

  end subroutine get_single_precision_real_value
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_double_precision_real_value( lyrics, real_variable, relational_operator, reference_value )
  !
  ! Purpose:
  !	Reads an real-Value. Optional: Compares the Value with a Reference.
  !	The "real_variable" is double precision and the variable "reference_value" is single precision.
  !
  ! Record revisions: 
  !	Date        Programmer    Desciption of Change
  !	====        ==========    ====================
  !	02.08.2008  Robert Opitz  Orignal Code
  !	06.12.2009  Robert Opitz  iso_fortran_env
  !
    character( len=* ), intent( in ) :: lyrics
    real( kind=real64 ), intent( out ) :: real_variable
    character( len=* ), intent( in ), optional :: relational_operator
    real( kind=real64 ), intent( in ), optional :: reference_value
    
    real( kind=real64 ) :: value
    integer :: ioerror
    logical :: ok
  
    intrinsic :: present 

    loop: do
      write( unit=output_unit, fmt='(1x,a)', advance='no' ) lyrics
      read( unit=input_unit, fmt=*, iostat=ioerror ) value
      
      if ( ioerror /= 0 ) then
        write( unit=output_unit, fmt='(1x,a)' ) 'Your Input was incorrect'
        cycle loop
      end if
      
      if ( present( relational_operator ) .and. present( reference_value ) ) then
        call check_value( ok )
        if ( .not. ok ) cycle loop
      else if ( present( relational_operator ) .or. present( reference_value ) ) then
        call stop_program( 'get_value error: relational_operator or reference_value missing' )
      end if
      
      real_variable = value
      exit loop

    end do loop
    
    contains

      subroutine check_value( ok )
      !
      !
      !
        logical, intent( out ) :: ok
        
        character(len=9), parameter :: fmt = '(1x,a,g0)'

        ok = .true.
        
        select case ( relational_operator )
          case( '<' )
            if ( value >= reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be smaller than : ', reference_value
              ok = .false.
            end if
          case( '<=' )
            if ( value > reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be smaller or equal than : ', reference_value
              ok = .false.
            end if
          case( '>' )
            if ( value <= reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be bigger than : ', reference_value
              ok = .false.
            end if
          case( '>=' )
            if ( value < reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be bigger or equal than : ', reference_value
              ok = .false.
            end if
          case( '/=' )
            if ( value == reference_value ) then
              write( unit=output_unit, fmt=fmt ) 'The Number has to be unequal to : ', reference_value
              ok = .false.
            end if
          case default
            call stop_program( 'get_value error: relational_operator missing' )
        end select
        
      end subroutine check_value

  end subroutine get_double_precision_real_value
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_integer_value( lyrics, integer_variable, relational_operator, reference_value )
  !
  ! Purpose:
  !	Reads an integer-Value. Optional: Compares the Value with a Reference. 
  !
  ! Record revisions: 
  !	Date        Programmer    Desciption of Change
  !	====        ==========    ====================
  !	21.05.2008  Robert Opitz  Orignal Code
  !	06.12.2009  Robert Opitz  iso_fortran_env
  !
  character( len=* ), intent( in ) :: lyrics
  integer, intent( out ) :: integer_variable
  character( len=* ), intent( in ), optional :: relational_operator
  integer, intent( in ), optional :: reference_value
  
  integer :: value
  integer :: ioerror
  logical :: ok
  
  intrinsic :: present
  
  loop: do
    write( unit=output_unit, fmt='(1x,a)', advance='no' ) lyrics
    read( unit=input_unit, fmt=*, iostat=ioerror ) value
    
    if ( ioerror /= 0 ) then
      write( unit=output_unit, fmt='(1x,a)' ) 'Your Input was incorrect'
      cycle loop
    end if
    
    if ( present( relational_operator ) .and. present( reference_value ) ) then
      call check_value( ok )
      if ( .not. ok ) cycle loop
    else if ( present( relational_operator ) .or. present( reference_value ) ) then
      call stop_program( 'get_value error: relational_operator or reference_value missing' )
    end if
    
    integer_variable = value
    exit loop
    
  end do loop
  
  contains
  
    subroutine check_value( ok )
    !
    !
    !
      logical, intent( out ) :: ok
      
      character(len=9), parameter :: fmt = '(1x,a,g0)'
      
      ok = .true.
      
      select case ( relational_operator )
        case( '<' )
          if ( value >= reference_value ) then
            write( unit=output_unit, fmt=fmt ) 'The Number has to be smaller than : ', reference_value
            ok = .false.
          end if
        case( '<=' )
          if ( value > reference_value ) then
            write( unit=output_unit, fmt=fmt ) 'The Number has to be smaller or equal than : ', reference_value
            ok = .false.
          end if
        case( '>' )
          if ( value <= reference_value ) then
            write( unit=output_unit, fmt=fmt ) 'The Number has to be bigger than : ', reference_value
            ok = .false.
          end if
        case( '>=' )
          if ( value < reference_value ) then
            write( unit=output_unit, fmt=fmt ) 'The Number has to be bigger or equal than : ', reference_value
            ok = .false.
          end if
        case( '/=' )
          if ( value == reference_value ) then
            write( unit=output_unit, fmt=fmt ) 'The Number has to be unequal to : ', reference_value
            ok = .false.
          end if
        case default
          call stop_program( 'get_value error: relational_operator missing' )
      end select
      
    end subroutine check_value

  end subroutine get_integer_value

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_logical_value( lyrics, logical_variable, true_symbol, &
                                false_symbol )
  !
  ! Purpose:
  !	Reads an logical-Value.
  !
  ! Record revisions: 
  !	Date        Programmer    Desciption of Change
  !	====        ==========    ====================
  !	21.05.2008  Robert Opitz  Orignal Code
  !	01.08.2008  Robert Opitz  Extension to optional true- and false-symbols	
  !	06.12.2009  Robert Opitz  iso_fortran_env
  !
    character( len=* ), intent( in ) :: lyrics
    logical, intent( out ) :: logical_variable
    character( len=* ), intent( in ), optional :: true_symbol
    character( len=* ), intent( in ), optional :: false_symbol
  
    character( len=32 ) :: true_value
    character( len=32 ) :: false_value
    character( len=32 ) :: input
    integer :: ioerror
    
    intrinsic :: present, trim
        
    !---Begin: Initialisation---------------------------------------------------
    if ( present( true_symbol ) ) then
      true_value = true_symbol
    else
      true_value = 'y'
    end if
    
    if ( present( false_symbol) ) then
      false_value = false_symbol
    else
      false_value = 'n'
    end if
    !---End: Initialisation-----------------------------------------------------

    !---Begin: Evaluation-------------------------------------------------------
    loop: do
      write( unit=output_unit, fmt='(1x,a,1x)', advance='no' ) lyrics
      read( unit=input_unit, fmt=*, iostat=ioerror ) input
      
      if ( trim( input ) == trim( true_value ) ) then
        logical_variable = .true.
        exit loop
      else if ( trim( input ) == trim( false_value ) ) then
        logical_variable = .false.
         exit loop
      else
        write( unit=output_unit, fmt='(1x,a,a,1x,a)' ) 'Your Input has to be one of these two Signs : ', &
               trim( true_value ), trim( false_value )
      end if
    
    end do loop
    !---End: Evaluation---------------------------------------------------------

  end subroutine get_logical_value

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module auxiliary_routines