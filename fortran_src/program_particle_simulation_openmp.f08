program particle_simulation
!
!
!
  use, intrinsic :: iso_fortran_env, only: rt => real32
  !$ use omp_lib
  use auxiliary_routines
  use stoermer_verlet_method, only: particle_type, &
                                    set_simulation_parameter, &
                                    compute_forces, compute_locations, &
                                    compute_velocities
  use particle_management, only: initialize_neighbor_lists, &
                                 update_neighbor_list_and_cell_space, &
                                 clean_up
  implicit none

  real(rt), parameter :: zero = 0.0_rt

  type(particle_type), dimension(:), allocatable :: particles

  real(rt), dimension(3) :: length_of_area, length_of_cell
  real :: start_time

  integer :: output_file_unit
  character(len = 132) :: environment_file, particle_file

  integer :: time_step, max_time_step, write_time_step
  integer :: i
  integer :: thread_id = 1, nb_of_threads
  integer, dimension(2) :: pindex
  logical, dimension(:), allocatable :: displacement_flag


  ! get all the files from the command link
  write(*,'(1x,a)') "[INFO] get command arguments"
  call get_command_arguments(environment_file, &
                             particle_file, &
                             output_file_unit)

  ! read fixed values from environment_file
  write(*,'(1x,a)') "[INFO] get and set simulation parameters"
  call get_and_set_parameters(environment_file, &
                              length_of_area, length_of_cell, &
                              write_time_step, max_time_step, nb_of_threads)

  ! read particle file and file particle_array with particles and their id,
  ! mass, location, and initial velocity
  write(*,'(1x,a)') "[INFO] initialze particle array"
  particles = initialize_particle_array(particle_file)

  allocate(displacement_flag(nb_of_threads))

  write(*,'(1x,a)') "[INFO] write size of simulation arae to output file"
  do i = 1, 3
    write(unit=output_file_unit, rec=i) length_of_area(i)
  end do

  ! COMBINE FOUR OR THREE PROCEDURES ABOVE IN ONE???
  call initialize_neighbor_lists(particles)

  write(*,'(1x,a)') "[INFO] start time integration"
  start_time = get_current_time()

  !$omp parallel num_threads(nb_of_threads) &
  !$omp& private(time_step, thread_id, pindex) &
  !$omp& firstprivate(max_time_step, write_time_step, nb_of_threads) &
  !$omp& shared(particles, displacement_flag)

  !$ thread_id = omp_get_thread_num() + 1

  pindex = get_particle_index(thread_id, nb_of_threads, size(particles))
  write(*,*) pindex

  ! compute the initial force for all particles in particle list by using
  ! their neighborlist
  write(*,'(1x,a,g0)') "[INFO] initialze forces, id = ", thread_id
  call compute_forces(particles, pindex)

  ! start loop for time integration
  time_integration: do time_step = 1, max_time_step

    !$omp barrier

    ! compute locations for each particle in particle list
    ! and check if particles are in need for updating their neighbor lists
    displacement_flag(thread_id) = compute_locations(particles, pindex)

    ! if the displacment is to large, then the neighbor_list and the
    ! cell_space for each particle needs to be updated
    !$omp barrier
    !$omp single
    if ( any(displacement_flag) .eqv. .true. ) then
      call update_neighbor_list_and_cell_space(particles)
      !if (no_active_particles(particles)) then
      !  write(*,'(1x,a)') "[INFO] no particles left in simulation area."
      !  !exit time_integration
      !end if
    end if
    !$omp end single

    ! compute forces for each particle in the particle list using their
    ! neighbor list
    call compute_forces(particles, pindex)

    ! compute velocities for each particle
    call compute_velocities(particles, pindex)

    ! if it is time to write data, write data
    if ( mod(time_step, write_time_step) == 0 ) then
      !$omp barrier
      !$omp single
      write(*,'(1x,a,1x,g0)') "[INFO] write for time step:", time_step
      call write_data(output_file_unit, particles)
      !$omp end single nowait
    end if

  end do time_integration
  !$omp end parallel

  call write_time(get_current_time() - start_time)

  write(*,'(1x,a)') "[INFO] clean up"
  call clean_up(particles)

  ! close output file
  close(output_file_unit)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_command_arguments(environment_file, &
                                   particle_file, &
                                   output_file_unit)
     use iso_fortran_env, only: error_unit
     character(len=*), intent(out) :: environment_file
     character(len=*), intent(out) :: particle_file
     integer, intent(out) :: output_file_unit

     character(len=132) :: program_name, &
                           name_of_output_file
     integer :: stat
     real(rt), parameter :: dummy_value = zero
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
        call get_command_argument(1, environment_file)
        ! Check file for existence
        call check_file(environment_file)

        call get_command_argument(2, particle_file)
        ! Check file for existence
        call check_file(particle_file)

        call get_command_argument(3, name_of_output_file)
        ! Check file for existence
        ! if file exist, stop program -> do not overwrite old outputfile
        inquire (iolength = length_of_record) dummy_value
        open(newunit = output_file_unit, file = name_of_output_file, &
             status = "replace", access = "direct", action = "write", &
             form = "unformatted", recl = length_of_record, iostat = stat)
        if (stat /= 0) then
          call stop_program("Could not open output file.", stat)
        end if
     end if

  end subroutine get_command_arguments

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_and_set_parameters(environment_file, &
                                    length_of_area, length_of_cell, &
                                    write_time_step, max_time_step, &
                                    nb_of_threads)
    character(len = 132), intent(in) :: environment_file
    real, dimension(3), intent(out) :: length_of_area, length_of_cell
    integer, intent(out) :: write_time_step, max_time_step
    integer, intent(out) :: nb_of_threads
    
    real(rt) :: epsilon_value, delta_t
    real(rt) :: particle_mass, particle_radius, r_cut_factor
    
    real(rt) :: sigma, epsilon, delta, mass, r_cut

    integer :: uef ! unit environment_file
    integer :: istat
    character(len=160) :: failure_msg

    namelist /ENVIRONMENTPARAMS/ nb_of_threads, sigma, epsilon, mass, &
                                 r_cut, length_of_area, length_of_cell, &
                                 delta_t, write_time_step, max_time_step

    open(newunit = uef, file = trim(environment_file), &
         status = 'old', action = 'read', iostat = istat, iomsg = failure_msg)
    if ( istat /= 0 ) then
      call stop_program( "OPEN failed for parameter-file '"//trim(environment_file)//"' with message: "//trim(failure_msg), istat )
    end if
    
    !read the namelist "paramterlist"
    read(unit = uef, nml = ENVIRONMENTPARAMS, iostat = istat, iomsg = failure_msg)
    if ( istat /= 0 ) then
      call stop_program( "READ failed for parameter-file '"//trim(environment_file)//"' with message: "//trim(failure_msg), istat )
    end if

    close(unit = uef, iostat = istat, iomsg = failure_msg)
    if ( istat /= 0 ) then
      call stop_program( "CLOSE failed for parameter-file '"//trim(environment_file)//"' with message: "//trim(failure_msg), istat )
    end if

    call set_simulation_parameter(delta_t, length_of_area, length_of_cell, &
                                  epsilon, sigma, mass, r_cut)

  end subroutine get_and_set_parameters

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function initialize_particle_array(particle_file) result(particles)
    character(len=*), intent(in) :: particle_file
    type(particle_type), dimension(:), allocatable :: particles !particle_array

    integer :: upf !unit_particle_file
    integer :: istat
    character(len=132) :: failure_msg

    integer :: id ! particle id
    integer :: n ! number of particles
    type(particle_type), dimension(:), allocatable :: temp
    real(rt), dimension(3) :: location, velocity

    ! open particle file unit
    open(newunit = upf, file = particle_file, &
         status = "old", action = "read", access = "sequential", &
         iostat = istat, iomsg = failure_msg)
    if (istat /= 0) then
      call stop_program("Could not open particle file '"//trim(particle_file)//"' with message: "//trim(failure_msg), istat)
    end if

    ! write particle data from particle_file to particle_array
    n = 0
    allocate(temp(0))
    do
      ! read particle data
      read(upf, *, iostat = istat, iomsg = failure_msg) id, location, velocity
      if (istat < 0) then
        exit ! end of file
      else if (istat > 0) then
        ! read error
        call stop_program("Could not read entry from particle file '"//trim(particle_file)//"' with message: "//trim(failure_msg), istat)
      end if
      ! add new particle
      n = n + 1
      if (n > size(temp)) call increase_particle_array(temp)
      ! write particle data to temp particle array
      temp(n)%id       = id
      temp(n)%location = location
      temp(n)%velocity = velocity
      allocate(temp(n)%neighbor_list(0))
    end do

    if (n == 0) then
      call stop_program("No particles in particle file")
    else
      write(*,'(1x,a,1x,g0)') "[INFO] Number of particles in file:", n
    end if

    ! create actual particle array
    allocate(particles(n))
    particles = temp(1:n)
    deallocate(temp)

    ! close particle file
    close(unit = upf, iostat = istat, iomsg = failure_msg)
    if ( istat /= 0 ) then
      call stop_program( "CLOSE failed for parameter-file '"//trim(environment_file)//"' with message: "//trim(failure_msg), istat )
    end if

  end function initialize_particle_array

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine increase_particle_array(p)
    type(particle_type), dimension(:), allocatable, intent(in out) :: p

    type(particle_type), dimension(:), allocatable :: temp
    integer :: n

    n = size(p)
    allocate(temp(n+1000))
    temp(1:n) = p
    call move_alloc(from = temp, to = p)

  end subroutine increase_particle_array

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  pure function get_particle_index(thread_id, num_threads, &
                                   nb_particles) result(pindex)
    integer, intent(in) :: thread_id
    integer, intent(in) :: num_threads
    integer, intent(in) :: nb_particles
    integer, dimension(2) :: pindex

    integer :: particles_per_thread

    particles_per_thread = nb_particles / num_threads

    pindex(1) = particles_per_thread * (thread_id - 1) + 1
    if (thread_id == num_threads) then
      pindex(2) = nb_particles
    else
      pindex(2) = particles_per_thread * (thread_id)
    end if

  end function get_particle_index

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine write_data(ofu, particles)
     integer, intent(in) :: ofu ! output_file_unit
     type(particle_type), dimension(:), intent(in) :: particles !particle_array

     integer, save :: counter = 4, number_of_entries = 0
     integer :: d
     integer :: j, i

     ! number of entries
     number_of_entries = number_of_entries + 1
     write(unit=ofu, rec=4) number_of_entries

     counter = counter + 1
     d = counter
     do i = 1, size(particles)
       if (.not. particles(i)%is_active) cycle

       d = d + 1
       write(unit=ofu, rec=d) particles(i)%id

       ! write location of particle
       do j = 1, 3
         d = d + 1
         write(unit=ofu, rec=d) particles(i)%location(j)
       end do

       ! write velocity of particle
       do j = 1, 3
         d = d + 1
         write(unit=ofu, rec=d) particles(i)%velocity(j)
       end do
     end do

     ! write number of data per time entry
     write(unit=ofu, rec=counter) d - counter

     counter = d

   end subroutine write_data

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function no_active_particles(particles) result(flag)
    type(particle_type), dimension(:), intent(in) :: particles
    logical :: flag

    flag = .true.
    loop: do i = 1, size(particles)
      if (particles(i)%is_active) then
        flag = .false.
        exit loop
      end if
    end do loop

  end function no_active_particles 

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end program particle_simulation
