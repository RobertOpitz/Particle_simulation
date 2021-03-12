program particle_simulation
  implicit none

  real, parameter :: zero = 0.0
  real, parameter :: half = 0.5
  real, parameter :: one = 1.0
  real, parameter :: twelve = 12.0

  type particle_type
    integer :: id = -1
    real, dimension(3) :: location = zero
    real, dimension(3) :: displacement = zero
    real, dimension(3) :: velocity = zero
    real, dimension(3) :: force_old = zero
    real, dimension(3) :: force_new = zero
    !real, dimension(3) :: force_correction = zero
    integer, dimension(:), allocatable :: neighbor_list
    integer :: nb_of_neighbors = 0
  end type particle_type

  type(particle_type), dimension(:), allocatable :: particle_array

  !type particle_number
  !  integer :: i = -1 ! position of particle in particle_array
  !  type(particle_number), pointer :: next_particle => null()
  !end type particle_number

  !type :: cell_type
  !   type(particle_number), pointer :: first_particle => null()
  !end type cell_type

  !type(cell_type), dimension(:,:,:), allocatable :: linked_list

  real, dimension(3) :: length_of_area, length_of_cells
  real :: delta_t, time_mass_factor, epsilon_value, r_min, r_cut, F_rcut, &
          dr_max_tolerable, min_length_of_cell
  real :: start_time
  integer :: particle_file_unit, environment_file_unit, output_file_unit
  integer :: time_step, max_time_step
  integer :: i
  logical :: displacement_flag


  ! get all the files from the command link
  write(*,*) "[INFO] get command arguments"
  call get_command_arguments(environment_file_unit, &
                             particle_file_unit, &
                             output_file_unit)

  ! read fixed values from environment_file
  write(*,*) "[INFO] get parameters"
  call get_parameters(environment_file_unit, r_min, r_cut, &
                      length_of_area, length_of_cells, min_length_of_cell, &
                      epsilon_value, time_mass_factor, F_rcut, &
                      dr_max_tolerable, delta_t, max_time_step)

  ! read particle file and file particle_array with particles and their id,
  ! mass, location, and initial velocity
  write(*,*) "[INFO] initialize particle array"
  particle_array = initialize_particle_array(particle_file_unit)

  ! create linked_list space, and ill it with the particle_numbers
  !write(*,*) "[INFO] initialze cell space"
  !linked_list = initialze_linked_list(environment_file_unit, particle_array)

  ! update neighbor list for each particle in particle_array by using
  ! the linked_list space
  write(*,*) "[INFO] initialze neighbor lists"
  call update_neighbor_list(particle_array)!, linked_list)

  ! compute the initial force for all particles in particle list by using
  ! their neighborlist
  write(*,*) "[INFO] initialze forces"
  call compute_forces(particle_array)

  do i = 1, 3
    write(unit=output_file_unit, rec=i) length_of_area(i)
  end do

  write(*,*) "[INFO] start time integration"

  start_time = get_current_time()
  ! start loop for time integration
  time_integration: do time_step = 1, max_time_step

    ! compute locations for each particle in particle list
    call compute_locations(particle_array, displacement_flag)

    ! if the displacment is to large, then the neighbor_list and the
    ! linked_list for each particle needs to be updated
    if (displacement_flag .eqv. .TRUE.) then
    !  call update_linked_list(linked_list, particle_array)
      call update_neighbor_list(particle_array)!, linked_list)
    end if

    ! compute forces for each particle in the particle list using their
    ! neighbor list
    call compute_forces(particle_array)

    ! compute velocities for each particle
    call compute_velocities(particle_array)

    ! if it is time to write data, write data
    if (write_data_now(time_step)) then
      call write_data(output_file_unit, particle_array)!, &
                      !pack(linked_list, .true.))
    end if

  end do time_integration
  call write_time(get_current_time() - start_time)

  ! clean up
  call clean_up(particle_array)!, linked_list)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_command_arguments(environment_file_unit, &
                                   particle_file_unit, &
                                   output_file_unit)
     use iso_fortran_env, only: error_unit
     integer, intent(out) :: environment_file_unit
     integer, intent(out) :: particle_file_unit
     integer, intent(out) :: output_file_unit

     character(len=132) :: program_name, name_of_parameter_file, &
                           name_of_particle_file, name_of_output_file
     integer :: stat

     real, parameter :: dummy_value = zero
     integer :: length_of_record

     if (command_argument_count() /= 3) then
        if (command_argument_count() < 3) then
           write(error_unit, '(1x,a)') 'Not enough arguments!'
        else
           write(error_unit, '(1x,a)') 'To many arguments!'
        end if
        call get_command_argument(0, program_name)
        call stop_program('Usage is: '//trim(program_name)// &
                          ' <parameter-file> <particle-file> <output-file>')
     else
        call get_command_argument(1, name_of_parameter_file)
        ! Check file for existence
        !call check_file(name_of_parameter_file)
        open(newunit = environment_file_unit, file = name_of_parameter_file, &
             status = "old", action = "read", access = "sequential", &
             iostat = stat)
        if (stat /= 0) then
          call stop_program("[ERROR] Could not open parameter file.", stat)
        end if

        call get_command_argument(2, name_of_particle_file)
        ! Check file for existence
        !call check_file(name_of_particle_file)
        open(newunit = particle_file_unit, file = name_of_particle_file, &
             status = "old", action = "read", access = "sequential", &
             iostat = stat)
        if (stat /= 0) then
          call stop_program("[ERROR] Could not open particle file.", stat)
        end if

        call get_command_argument(3, name_of_output_file)
        ! Check file for existence
        ! if file exist, stop program -> do not overwrite old outputfile
        !call check_file(name_of_output_file, overwrite_existing_file = .false.)
        inquire (iolength = length_of_record) dummy_value
        open(newunit = output_file_unit, file = name_of_output_file, &
             status = "replace", access = "direct", action = "write", &
             form = "unformatted", recl = length_of_record, iostat=stat)
        if (stat /= 0) then
          call stop_program("[ERROR] Could not open output file.", stat)
        end if
     end if

  end subroutine get_command_arguments

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function write_data_now(i) result(write_now)
    integer, intent(in) :: i
    logical :: write_now
    integer, save :: write_step = 0

    write_now = .false.
    if (i == 1 .or. i == max_time_step .or. i == write_step) then
      write(*,*) "[INFO] write now, ", i
      write_step = write_step + 10
      write_now = .true.
    else
      write_now = .false.
    end if

  end function write_data_now

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_parameters(environment_file_unit, r_min, r_cut, &
                            length_of_area, length_of_cells, &
                            min_length_of_cell, epsilon_value, &
                            time_mass_factor, F_rcut, dr_max_tolerable, &
                            delta_t, max_time_step)
     integer, intent(in) :: environment_file_unit
     real, intent(out) :: r_min, r_cut, &
                          epsilon_value, time_mass_factor, delta_t, F_rcut, &
                          dr_max_tolerable, min_length_of_cell
     real, dimension(3), intent(out) :: length_of_area, length_of_cells
     integer, intent(out) :: max_time_step

     real :: particle_mass, particle_radius, r_cut_factor
     real :: s


     read(environment_file_unit, *) particle_radius, epsilon_value, particle_mass
     read(environment_file_unit, *) r_cut_factor
     read(environment_file_unit, *) length_of_area
     write(*,*) length_of_area
     read(environment_file_unit, *) length_of_cells
     read(environment_file_unit, *) delta_t
     read(environment_file_unit, *) max_time_step

     ! set the maximal acceptable displacement value (if displacement of one
     ! particle is larger than this value, the neighbor list for each particle
     ! needs to be updated)
     r_min = 2.0**(1.0/6.0) * particle_radius
     r_cut = r_cut_factor * particle_radius
     min_length_of_cell = minval(length_of_cells)
     dr_max_tolerable = min_length_of_cell - r_cut
     ! set fixed values for particle simulation
     s = (r_min / r_cut)**6
     F_rcut = (s / (r_cut * r_cut)) * (one - s)
     !
     time_mass_factor = half * delta_t / particle_mass

  end subroutine get_parameters

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function initialize_particle_array(particle_file_unit) result(particle_array)
    integer, intent(in) :: particle_file_unit
    type(particle_type), dimension(:), allocatable :: particle_array

    !integer :: nb_of_particles,
    integer :: id
    integer :: i, istat
    real, dimension(3) :: location, velocity

    ! write particle data from particle_file to particle_array
    i = 0
    do
       ! read particle data
       read(particle_file_unit, *, iostat = istat) id, location, velocity
       if (istat /= 0) exit ! end of file or something is wrong
       call increase_particle_array(particle_array)
       i = i + 1
       ! write particle data to particle array
       particle_array(i)%id = id
       particle_array(i)%location = location
       particle_array(i)%velocity = velocity
       allocate(particle_array(i)%neighbor_list(1))
    end do

    write(*,*) "[INFO] Nb of particles: ", i

    ! close particle file
    close(particle_file_unit)

  end function initialize_particle_array

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine increase_particle_array(p)
    type(particle_type), dimension(:), allocatable, intent(in out) :: p ! particle_array

    type(particle_type), dimension(:), allocatable :: temp
    integer :: n

    if (allocated(p)) then
      n = size(p)
      allocate(temp(n+1))
      temp(1:n) = p
      call move_alloc(from = temp, to = p)
    else
      allocate(p(1))
    end if

  end subroutine increase_particle_array

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! function initialze_linked_list(environment_file_unit, &
  !                                particle_array) result(cell_space)
  !   integer, intent(in) :: environment_file_unit
  !   type(particle_type), dimension(:), intent(in) :: particle_array
  !   type(cell_type), dimension(:,:,:), allocatable :: cell_space
  !
  !   integer, dimension(3) :: nc
  !   integer :: istat
  !   type(particle_number), pointer :: new_particle_number
  !   integer :: i
  !
  !   !
  !   nc = nint(length_of_area / length_of_cells)
  !   write(*,*) 'Number of cells ', nc
  !
  !   ! allocate linked_list
  !   allocate(cell_space(nc(1), nc(2), nc(3)), stat = istat)
  !   if (istat /= 0) then
  !     call stop_program('[ERROR] could not allocate linked_cells array', istat)
  !   end if
  !
  !   ! fill the linked_list_array with particle numbers
  !   do i = 1, size(particle_array)
  !     ! Create new particle
  !     allocate(new_particle_number, stat=istat)
  !     if (istat /= 0) then
  !       call stop_program('[ERROR] ALLOCATE failed for "new_particle_number"', &
  !                         istat)
  !     end if
  !
  !     new_particle_number%i = i
  !
  !     ! Compute position of target cell of particle
  !     nc = ceiling(particle_array(i)%location / length_of_cells)
  !
  !     ! connect particle with target cell in cell_space
  !     new_particle_number%next_particle => &
  !                                   cell_space(nc(1),nc(2),nc(3))%first_particle
  !     cell_space(nc(1),nc(2),nc(3))%first_particle => new_particle_number
  !   end do
  !
  !   ! close environment file
  !   close(environment_file_unit)
  !
  ! end function initialze_linked_list

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine update_neighbor_list(p)!, linked_list)
    type(particle_type), dimension(:), intent(in out) :: p ! particle_array
    !type(cell_type), dimension(:,:,:), intent(in) :: linked_list

    integer :: i, n, j

    n = size(p)
    do i = 1, n - 1
      particle_array(i)%displacement = zero
      particle_array(i)%nb_of_neighbors = 0
      do j = i + 1, n
        call add_neighbor_particle(p(i), p(j), j)
      end do
    end do

    particle_array(n)%displacement = zero
    particle_array(n)%nb_of_neighbors = 0

    !call collect_neighbors_within_cell(pack(linked_list, .TRUE.), &
    !                                   particle_array)

    !call collect_neighbors_from_other_cells(linked_list, particle_array)

  end subroutine update_neighbor_list

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! subroutine collect_neighbors_within_cell(cells, particle_array)
  !    type(cell_type), dimension(:), intent(in) :: cells
  !    type(particle_type), dimension(:), intent(in out) :: particle_array
  !
  !    type(particle_number), pointer :: p1
  !    type(particle_number), pointer :: p2
  !    integer :: ci
  !
  !    cell_loop: do ci = 1, size(cells)
  !
  !       p1 => cells(ci)%first_particle
  !       p1_loop: do while(associated(p1))
  !
  !         ! RESET DISPLACEMENT AND NUMBER OF NEIGHBORS
  !         particle_array(p1%i)%displacement = zero
  !         particle_array(p1%i)%nb_of_neighbors = 0
  !
  !         p2 => p1%next_particle
  !         p2_loop: do while(associated(p2))
  !           call add_neighbor_particle(particle_array(p1%i), &
  !                                      particle_array(p2%i), &
  !                                      p2%i)
  !           p2 => p2%next_particle
  !         end do p2_loop
  !
  !         p1 => p1%next_particle
  !     end do p1_loop
  !
  !   end do cell_loop
  !
  ! end subroutine collect_neighbors_within_cell

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! subroutine collect_neighbors_from_other_cells(cell_space, particle_array)
  ! !
  ! ! internal subroutine of subroutine compF_LC
  ! !
  ! ! Purpose
  ! !
  ! ! ----------------------------------------------------------------------------
  ! !      ---------------------------------
  ! !      !              !                !
  ! !      ! cx, cy+1, cz ! cx+1, cy+1, cz !
  ! !      !              !                !
  ! !      ---------------------------------
  ! !      !              !                !
  ! !      ! cx, cy, cz   ! cx+1, cy, cz   !
  ! !      !              !                !
  ! !      ---------------------------------
  ! !
  ! !      within cz Plane:
  ! !          (i)   cell(cx,cy,cz) with cell(cx,cy+1,cz)
  ! !          (ii)  cell(cx,cy,cz) with cell(cx+1,cy,cz)
  ! !          (iii) cell(cx,cy,cz) with cell(cx+1,cy+1,cz)
  ! !          (iv)  cell(cx+1,cy,cz) with cell(cx,cy+1,cz)
  ! !
  ! ! ----------------------------------------------------------------------------
  ! !  ---------------------------------   -------------------------------------
  ! !  !              !                !   !                !                  !
  ! !  ! cx, cy+1, cz ! cx+1, cy+1, cz !   ! cx, cy+1, cz+1 ! cx+1, cy+1, cz+1 !
  ! !  !              !                !   !                !                  !
  ! !  ---------------------------------   -------------------------------------
  ! !  !              !                !   !                !                  !
  ! !  ! cx, cy, cz   ! cx+1, cy, cz   !   ! cx, cy, cz+1   ! cx+1, cy, cz+1   !
  ! !  !              !                !   !                !                  !
  ! !  ---------------------------------   -------------------------------------
  ! !
  ! !      cz with cz+1 Plane:
  ! !            cz plane                cz+1 plane
  ! !      (1)   cell(cx,cy,cz)     with cell(cx,cy,cz+1)
  ! !      (2.a) cell(cx,cy,cz)     with cell(cx,cy+1,cz+1)
  ! !      (2.b) cell(cx,cy+1,cz)   with cell(cx,cy,cz+1)
  ! !      (3.a) cell(cx,cy,cz)     with cell cell(cx+1,cy,cz+1)
  ! !      (3.b) cell(cx+1,cy,cz)   with cell(cx,cy,cz+1)
  ! !      (4.a) cell(cx,cy,cz)     with cell(cx+1,cy+1,cz+1)
  ! !      (4.b) cell(cx+1,cy+1,cz) with cell(cx,cy,cz+1)
  ! !      (5.a) cell(cx+1,cy,cz)   with cell(cx,cy+1,cz+1)
  ! !      (5.b) cell(cx,cy+1,cz)   with cell(cx+1,cy,cz+1)
  ! !
  ! ! ----------------------------------------------------------------------------
  !    type(cell_type), dimension(:,:,:), intent(in) :: cell_space
  !    type(particle_type), dimension(:), intent(in out) :: particle_array
  !
  !    integer :: cx, cx_max
  !    integer :: cy, cy_max
  !    integer :: cz, cz_max
  !
  !    cx_max = size(cell_space, dim = 1)
  !    cy_max = size(cell_space, dim = 2)
  !    cz_max = size(cell_space, dim = 3)
  !
  !    !---cz Ebene---------------------------------------------------------------
  !    z_loop_a: do cz = 1, cz_max
  !
  !       ! (i) cell(cx,cy,cz) with cell(cx,cy+1,cz)
  !       do cy = 1, cy_max - 1
  !          do cx = 1, cx_max
  !             call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
  !                                                cell_space(cx, cy+1, cz), &
  !                                                particle_array)
  !          end do
  !       end do
  !
  !       ! (ii) cell(cx,cy,cz) with cell(cx+1,cy,cz)
  !       do cy = 1, cy_max
  !          do cx = 1, cx_max - 1
  !             call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
  !                                                cell_space(cx+1, cy, cz), &
  !                                                particle_array)
  !          end do
  !       end do
  !
  !       ! (iii) cell(cx,cy,cz)   with cell(cx+1,cy+1,cz)
  !       ! (iv)  cell(cx+1,cy,cz) with cell(cx,cy+1,cz)
  !       do cy = 1, cy_max - 1
  !          do cx = 1, cx_max - 1
  !             call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
  !                                                cell_space(cx+1, cy+1, cz), &
  !                                                particle_array)
  !             call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
  !                                                cell_space(cx, cy+1, cz), &
  !                                                particle_array)
  !          end do
  !       end do
  !
  !    end do z_loop_a
  !    !--------------------------------------------------------------------------
  !
  !    !---cz mit cz+1 Ebene------------------------------------------------------
  !    z_loop_b: do cz = 1, cz_max - 1
  !
  !       ! (1) cell(cx,cy,cz) with cell(cx,cy,cz+1)
  !       do cy = 1, cy_max
  !          do cx = 1, cx_max
  !             call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
  !                                                cell_space(cx, cy, cz+1), &
  !                                                particle_array)
  !          end do
  !       end do
  !
  !       ! (2)
  !       do cy = 1, cy_max - 1
  !          do cx = 1, cx_max
  !             ! (2.a) cell(cx,cy,cz) with cell(cx,cy+1,cz+1)
  !             call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
  !                                                cell_space(cx, cy+1, cz+1), &
  !                                                particle_array)
  !             ! (2.b) cell(cx,cy+1,cz) with cell(cx,cy,cz+1)
  !             call add_neighbors_from_other_cell(cell_space(cx, cy+1, cz), &
  !                                                cell_space(cx, cy, cz+1), &
  !                                                particle_array)
  !          end do
  !       end do
  !
  !       ! (3)
  !       do cy = 1, cy_max
  !          do cx = 1, cx_max - 1
  !             ! (3.a) cell(cx,cy,cz) with cell(cx+1,cy,cz+1)
  !             call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
  !                                                cell_space(cx+1, cy, cz+1), &
  !                                                particle_array)
  !             ! (3.b) cell(cx+1,cy,cz) with cell(cx,cy,cz+1)
  !             call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
  !                                                cell_space(cx, cy, cz+1), &
  !                                                particle_array)
  !          end do
  !       end do
  !
  !       ! (4) and (5)
  !       do cy = 1, cy_max - 1
  !          do cx = 1, cx_max - 1
  !             ! (4.a) cell(cx,cy,cz) with cell(cx+1,cy+1,cz+1)
  !             call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
  !                                                cell_space(cx+1, cy+1, cz+1), &
  !                                                particle_array)
  !             ! (4.b) cell(cx+1,cy+1,cz) with cell(cx,cy,cz+1)
  !             call add_neighbors_from_other_cell(cell_space(cx+1, cy+1, cz), &
  !                                                cell_space(cx, cy, cz+1), &
  !                                                particle_array)
  !             ! (5.a) cell(cx+1,cy,cz) with cell(cx,cy+1,cz+1)
  !             call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
  !                                                cell_space(cx, cy+1, cz+1), &
  !                                                particle_array)
  !             ! (5.b) cell(cx,cy+1,cz) with cell(cx+1,cy,cz)
  !             call add_neighbors_from_other_cell(cell_space(cx, cy+1, cz), &
  !                                                cell_space(cx+1, cy, cz+1), &
  !                                                particle_array)
  !          end do
  !       end do
  !
  !    end do z_loop_b
  !    !--------------------------------------------------------------------------
  !
  ! end subroutine collect_neighbors_from_other_cells

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! subroutine add_neighbors_from_other_cell(this_cell, other_cell, &
  !                                          particle_array)
  !    type(cell_type), intent(in) :: this_cell
  !    type(cell_type), intent(in) :: other_cell
  !    type(particle_type), dimension(:), intent(in out) :: particle_array
  !
  !    type(particle_number), pointer :: p1
  !    type(particle_number), pointer :: p2
  !
  !    if (.not. associated(other_cell%first_particle)) return
  !
  !    p1 => this_cell%first_particle
  !    p1_loop: do while(associated(p1))
  !
  !       p2 => other_cell%first_particle
  !       p2_loop: do while(associated(p2))
  !          call add_neighbor_particle(particle_array(p1%i), &
  !                                     particle_array(p2%i), &
  !                                     p2%i)
  !          p2 => p2%next_particle
  !       end do p2_loop
  !
  !       p1 => p1%next_particle
  !
  !    end do p1_loop
  !
  ! end subroutine add_neighbors_from_other_cell

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine add_neighbor_particle(p1, p2, p2_i)
  !
     type(particle_type), intent(in out) :: p1
     type(particle_type), intent(in) :: p2
     integer, intent(in) :: p2_i

     real :: r

     r = norm2(p1%location - p2%location)

     if (r <= min_length_of_cell) then

        p1%nb_of_neighbors = p1%nb_of_neighbors + 1

        !if (.not. allocated(p1%neighbor_list) .and. p1%nb_of_neighbors == 1) then
        !  allocate(p1%neighbor_list(1))
        if (p1%nb_of_neighbors > size(p1%neighbor_list)) then
           ! increase size of neighbor_list by one
           call increase_neighbor_list_size(p1%neighbor_list)
        end if

        p1%neighbor_list(p1%nb_of_neighbors) = p2_i

     end if

   end subroutine add_neighbor_particle

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   pure subroutine increase_neighbor_list_size(neighbor_list)
     integer, dimension(:), allocatable, intent(in out) :: neighbor_list

     integer, dimension(:), allocatable :: temp
     integer :: n

     n = size(neighbor_list)
     allocate(temp(n+1))
     temp(1:n) = neighbor_list
     call move_alloc(from = temp, to = neighbor_list)

   end subroutine increase_neighbor_list_size

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! subroutine update_linked_list(linked_cells, particle_array)
   !   type(cell_type), dimension(:,:,:), intent(in out) :: linked_cells
   !   type(particle_type), dimension(:), intent(in) :: particle_array
   !
   !   type(particle_number), pointer :: particle, previous_particle
   !   integer :: cx, cy, cz
   !   integer, dimension(3) :: new_cell
   !
   !   z_loop: do cz = 1, size(linked_cells, dim = 3)
   !     y_loop: do cy = 1, size(linked_cells, dim = 2)
   !       x_loop: do cx = 1, size(linked_cells, dim = 1)
   !
   !         ! check if no particles are in this cell
   !         if (.not. associated(linked_cells(cx, cy, cz)%first_particle)) then
   !           ! no particles, go for the next cell
   !           cycle x_loop
   !         end if
   !
   !         ! init
   !         particle => linked_cells(cx, cy, cz)%first_particle
   !         previous_particle => null()
   !
   !         particle_loop: do
   !           ! compute coordinates of new cell of particle
   !           new_cell = ceiling(particle_array(particle%i)%location / &
   !                              length_of_cells)
   !
   !           whereto: if (all(new_cell == (/cx, cy, cz/))) then
   !             ! Particle does not left its cell
   !             if (associated(particle%next_particle)) then
   !               ! go to next particle
   !               previous_particle => particle
   !               particle => particle%next_particle
   !             else
   !               ! end of particle list
   !               exit particle_loop
   !             end if
   !           else whereto
   !             ! Particle left its cell.
   !             ! Unhinge particle from its cell.
   !             ! Hinge particle to its new cell.
   !             ! Change the pointer settings.
   !             follower: if (associated(particle%next_particle)) then
   !               if (associated(previous_particle)) then
   !                 ! Particle has predecessor and successor
   !                 previous_particle%next_particle => particle%next_particle
   !                 !call entry_in_new_list(particle, new_cell, &
   !                !                        particle_buffer, linked_cells)
   !                 call entry_in_new_list(particle, new_cell, linked_cells)
   !                 particle => previous_particle%next_particle
   !               else
   !                 ! Particle is first particle in list
   !                 linked_cells(cx, cy, cz)%first_particle => particle%next_particle
   !                 call entry_in_new_list(particle, new_cell, linked_cells)
   !                 particle => linked_cells(cx, cy, cz)%first_particle
   !               end if
   !             else follower
   !               if (associated(previous_particle)) then
   !                 ! Particle is the last one in list
   !                 previous_particle%next_particle => null()
   !               else
   !                 ! Particle is the first and last one in list
   !                 linked_cells(cx, cy, cz)%first_particle => null()
   !               end if
   !               call entry_in_new_list(particle, new_cell, linked_cells)
   !               ! no more particles in list
   !               exit particle_loop
   !             end if follower
   !           end if whereto
   !         end do particle_loop
   !
   !       end do x_loop
   !     end do y_loop
   !   end do z_loop
   !
   ! end subroutine update_linked_list

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! subroutine entry_in_new_list(particle, new_cell, &
   !                              linked_cells)
   !   type(particle_number), pointer, intent(in out) :: particle
   !   integer, dimension(3), intent(in) :: new_cell
   !   type(cell_type), dimension(:,:,:), intent(in out) :: linked_cells
   !
   !   integer, dimension(3) :: max_cell
   !   integer :: x, y, z
   !
   !   max_cell = shape(linked_cells)
   !
   !   if (any(new_cell < (/1,1,1/)) .or. any(new_cell > max_cell)) then
   !     call delete_particle(particle)
   !   else
   !     x = new_cell(1)
   !     y = new_cell(2)
   !     z = new_cell(3)
   !     particle%next_particle => linked_cells(x,y,z)%first_particle
   !     linked_cells(x,y,z)%first_particle => particle
   !   end if
   !
   ! end subroutine entry_in_new_list

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! subroutine delete_particle(particle)
   !   type(particle_number), pointer, intent(in out) :: particle
   !
   !   type(particle_number), pointer :: tmp
   !   integer :: istat
   !
   !   tmp => particle
   !   deallocate(tmp, stat = istat)
   !   if (istat /= 0) then
   !     call stop_program("[ERROR] Deallocate failed for deleting "// &
   !                       "particle number.", istat)
   !   end if
   !
   ! end subroutine delete_particle

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine compute_forces(p)
     type(particle_type), dimension(:), intent(in out) :: p !particle_array

     integer :: i, j, n

     n = size(p)
     do i = 1, n - 1
       !if (i == 100) write(*,*) "nb of neighbors: ", p(1)%nb_of_neighbors, &
      !                          size(p(1)%neighbor_list)
       do j = i + 1, n!p(i)%nb_of_neighbors
         call compute_force_between_particles(p(i), p(j))
       end do
     end do

   end subroutine compute_forces

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   pure subroutine compute_force_between_particles(p1, p2)
     type(particle_type), intent(in out) :: p1, p2

     real, dimension(3) :: r_ij
     real, dimension(3) :: force
     real :: r
     real :: s
     real :: F_rij

     ! Direction Vector (Richtungsvektor)
     r_ij = p2%location - p1%location
     ! Distance between both Particles
     r = norm2(r_ij)

     if (r < r_cut) then
       ! compute epsilon for two particles
       !r_min = half * (p1%r_min, p2%r_min)
       !epsilon_value = sqrt(p1%epsilon * p2%epsilon)
       ! FORCE FOR TRUNCATED LENARD-JONES POTENTIAL
       s = (r_min / r)**6
       F_rij = (s / (r * r)) * (one - s)
       ! TRUNCATED FORCE
       force = r_ij * twelve * epsilon_value * (F_rij - F_rcut)
       ! add forces
       p1%force_new = p1%force_new + force
       p2%force_new = p2%force_new - force
     end if

   end subroutine compute_force_between_particles

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine compute_locations(p, displacement_flag)
     type(particle_type), dimension(:), intent(in out) :: p !particle_array
     logical, intent(out) :: displacement_flag

     real, dimension(3) :: dx
     real :: dr, dr_max, dr_max_first, dr_max_second
     integer :: i

     dr_max_first = zero
     dr_max_second = zero

     do i = 1, size(p)
       ! INIT FORCE
       ! the formerly new Force is now the old Force
       p(i)%force_old = p(i)%force_new
       ! Reset the new Force to dummy value
       p(i)%force_new = zero !particle%outer_force(:)
       ! reset correction value for the sum to zero
       !p%force_correction = zero
       ! COMPUTE LOCATION
       dx = delta_t * (p(i)%velocity + time_mass_factor * p(i)%force_old)
       p(i)%location = p(i)%location + dx
       p(i)%displacement = p(i)%displacement + dx
       ! ! GET MAX DISPLACEMENT OF PARTICLE
       dr = norm2(p(i)%displacement)
       if (dr > dr_max_first) then
           dr_max_second = dr_max_first
           dr_max_first = dr
       else if (dr > dr_max_second) then
           dr_max_second = dr
       end if
     end do

     dr_max = dr_max_first + dr_max_second

     ! if the maximal displacement becomes to large, update the neighbor list
     if (dr_max > dr_max_tolerable) then
       displacement_flag = .true.
     else
       displacement_flag = .false.
     end if

   end subroutine compute_locations

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine compute_velocities(p)
     type(particle_type), dimension(:), intent(in out) :: p ! particle_array

     integer :: i

     do i = 1, size(p)
        p(i)%velocity = p(i)%velocity + time_mass_factor * &
                        (p(i)%force_new + p(i)%force_old)
     end do

   end subroutine compute_velocities

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine write_data(ofu, p)!, linked_cells)
     integer, intent(in) :: ofu ! output_file_unit
     type(particle_type), dimension(:), intent(in) :: p !particle_array
     !type(cell_type), dimension(:), intent(in) :: linked_cells

     !type(particle_number), pointer :: p
     integer, save :: counter = 4, number_of_entries = 0
     integer :: d, j
     !integer :: ci

     !write(*,*) "Particles: ", &
    !            get_number_of_particles_in_area(linked_cells)

     ! number of entries
     number_of_entries = number_of_entries + 1
     write(unit=ofu, rec=4) number_of_entries

     counter = counter + 1
     d = counter
     do i = 1, size(p)
     !all_cells: do ci = 1, size(linked_cells)

       !p => linked_cells(ci)%first_particle
       !this_cell: do while(associated(p))
         ! write particle id
         d = d + 1
         write(unit=ofu, rec=d) p(i)%id

         ! write location of particle
         do j = 1, 3
           d = d + 1
           write(unit=ofu, rec=d) p(i)%location(j)
         end do

         ! write velocity of particle
         do j = 1, 3
           d = d + 1
           write(unit=ofu, rec=d) p(i)%velocity(j)
         end do

         !p => p%next_particle
       !end do this_cell
     !end do all_cells
     end do

     ! write number of data per time entry
     write(unit=ofu, rec=counter) d - counter

     counter = d

   end subroutine write_data

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine clean_up(particle_array)!, cell_space)
     type(particle_type), dimension(:), allocatable, intent(in out) :: particle_array
     !type(cell_type), dimension(:,:,:), allocatable, intent(in out) :: cell_space

     !type(particle_number), pointer :: particle
     !integer :: cz, cy, cx
     integer :: istat

     deallocate(particle_array, stat=istat)
     if (istat /= 0) then
        call stop_program("[ERROR] Failed deallocating particle_array.", istat)
     end if

     ! ! remove elements connected to the cell_space
     ! z_loop: do cz = 1, size(cell_space, dim=3)
     !    y_loop: do cy = 1, size(cell_space, dim=2)
     !       x_loop: do cx = 1, size(cell_space, dim=1)
     !
     !          particle_loop: do
     !             particle => cell_space(cx, cy, cz)%first_particle
     !             if (.not. associated(particle)) exit particle_loop
     !
     !             cell_space(cx, cy, cz)%first_particle => &
     !                         cell_space(cx, cy, cz)%first_particle%next_particle
     !
     !             deallocate(particle, stat=istat)
     !             if (istat /= 0) then
     !                call stop_program("[ERROR] Failed deallocating particle. ",&
     !                                  istat)
     !            end if
     !
     !          end do particle_loop
     !
     !       end do x_loop
     !    end do y_loop
     ! end do z_loop
     !
     ! ! deallocate the empty cell_space array
     ! deallocate(cell_space, stat=istat)
     ! if (istat /= 0) then
     !    call stop_program("[ERROR] Failed deallocating cell_space.", istat)
     ! end if

   end subroutine clean_up

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine stop_program(error_message, ioerror)
   !
   ! Purpose:
   !	   If an Error accured tell User and stop the Program.
   !
   ! Record revisions:
   !     Date			    Programmer		  Desciption of Change
   !	   ====			    ==========		  ====================
   !	   27.04.2008		Robert Opitz		Orignal Code
   !	   06.12.2009		Robert Opitz		iso_fortran_env
   !     28.09.2015   Robert Opitz    replaced i4 with 1x,g0
   !
    use iso_fortran_env, only: error_unit

      character(len=*), intent(in), optional :: error_message
      integer, intent(in), optional :: ioerror

      if (present(error_message)) then
         write(unit=error_unit, fmt='(1x,a)') error_message
      end if

      if (present(ioerror)) then
         write(unit=error_unit, fmt='(1x,a,1x,g0)') 'Runtime Error Code : ', &
                                                    ioerror
      end if

      write(unit=error_unit, fmt='(1x,a)') 'Program stopped'
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
      real :: elapsed_time

      write(*,'(1x,3(a,g0))') '[INFO] Elasped time (h:min:sec.msec) is ', &
                               int(elapsed_time / 3600.0), &
                               ':', int(mod(elapsed_time / 60.0, 60.0)), &
                               ':', mod(elapsed_time, 60.0)

   end subroutine write_time

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! function get_number_of_particles_in_area(cell_space) result(i)
   !   type(cell_type), dimension(:), intent(in) :: cell_space
   !   integer :: i
   !
   !   type(particle_number), pointer :: p
   !   integer :: ci
   !
   !   i = 0
   !   do ci = 1, size(cell_space)
   !     p => cell_space(ci)%first_particle
   !     do while(associated(p))
   !       i = i + 1
   !       p => p%next_particle
   !     end do
   !   end do
   ! end function get_number_of_particles_in_area

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end program particle_simulation
