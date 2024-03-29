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
        write(unit=error_unit, &
              fmt='(1x,a,1x,g0)') '[ERROR] Runtime Error Code : ', ioerror
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
    use iso_fortran_env, only: output_unit
    real :: elapsed_time

    write(output_unit,'(1x,3(a,g0))') &
                            '[INFO] Elasped time (h:min:sec.msec) is ', &
                            int(elapsed_time / 3600.0), &
                            ':', int(mod(elapsed_time / 60.0, 60.0)), &
                            ':', mod(elapsed_time, 60.0)

  end subroutine write_time

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module auxiliary_routines

!===============================================================================

module stoermer_verlet_method
  implicit none
  private

  real, parameter :: zero = 0.0
  real, parameter :: half = 0.5
  real, parameter :: one = 1.0
  real, parameter :: two = 2.0
  real, parameter :: six = 6.0
  real, parameter :: twelve = 12.0

  type particle_type
    integer :: id = -1
    real, dimension(3) :: location = zero
    real, dimension(3) :: displacement = zero
    real, dimension(3) :: velocity = zero
    real, dimension(3) :: force_old = zero
    real, dimension(3) :: force_new = zero
    !logical :: is_active = .true.
    integer, dimension(:), allocatable :: neighbor_list
    integer :: nb_of_neighbors = 0
  end type particle_type

  !type particle_array_type
  !  type(particle_type), dimension(:), allocatable :: particles
  !  integer, dimension(:), allocatable :: active_particles
  !end type particle_array_type

  ! constants to be set with the start of the simulation
  real :: time_mass_factor, delta_t
  real :: dr_max_tolerable
  real :: r_min, r_cut
  real :: F_rcut
  real :: epsilon_value_times_12
  real :: min_length_of_cell
  real, dimension(3) :: length_of_area, length_of_cell

  ! visible to the outside of module
  ! constante and variables
  public :: zero, length_of_area, min_length_of_cell, length_of_cell
  ! types
  public :: particle_type!, particles, active_particles
  ! procedures
  public :: set_simulation_parameter
  public :: compute_forces, compute_locations, compute_velocities

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_simulation_parameter(input_delta_t, input_length_of_area, &
                                      input_length_of_cell, &
                                      input_epsilon_value, particle_radius, &
                                      particle_mass, r_cut_factor)
    real, intent(in) :: input_delta_t
    real, dimension(3), intent(in) :: input_length_of_area, input_length_of_cell
    real, intent(in) :: input_epsilon_value
    real, intent(in) :: particle_radius
    real, intent(in) :: particle_mass
    real, intent(in) :: r_cut_factor

    real :: s

    ! set the maximal acceptable displacement value (if displacement of one
    ! particle is larger than this value, the neighbor list for each particle
    ! needs to be updated)
    delta_t = input_delta_t
    epsilon_value_times_12 = input_epsilon_value * twelve
    r_min = two**(one/six) * particle_radius
    r_cut = r_cut_factor * particle_radius
    length_of_area = input_length_of_area
    length_of_cell = input_length_of_cell
    min_length_of_cell = minval(input_length_of_cell)
    dr_max_tolerable = min_length_of_cell - r_cut
    ! set fixed values for particle simulation
    s = (r_min / r_cut)**6
    F_rcut = (s / (r_cut * r_cut)) * (one - s)
    !
    time_mass_factor = half * delta_t / particle_mass

  end subroutine set_simulation_parameter

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine compute_forces(particles, active_particles)
     type(particle_type), dimension(:), intent(in out) :: particles
     integer, dimension(:), intent(in) :: active_particles

     integer :: ii, jj, m, n

     ! loop over all active particles
     do m = 1, size(active_particles)
       ! get index for ith particle
       ii = active_particles(m)
       ! if the ith particle as neighbors, compute the force with them
       if (particles(ii)%nb_of_neighbors > 0) then
         ! loop over all neighbor particles of ith particle
         do n = 1, particles(ii)%nb_of_neighbors
           ! get index of jth neighbor particle
           jj = particles(ii)%neighbor_list(n)
           call compute_force_between_particles(particles(ii), &
                                                particles(jj))
         end do
       end if
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
       ! FORCE FOR TRUNCATED LENARD-JONES POTENTIAL
       s = (r_min / r)**6
       F_rij = (s / (r * r)) * (one - s)
       ! TRUNCATED FORCE
       !force = r_ij * twelve * epsilon_value * (F_rij - F_rcut)
       force = r_ij * epsilon_value_times_12 * (F_rij - F_rcut)
       ! add forces, using Newtons actio est reactio
       p1%force_new = p1%force_new + force
       !$omp atomic update
       p2%force_new = p2%force_new - force
     end if

   end subroutine compute_force_between_particles

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine compute_locations(particles, active_particles, displacement_flag)
     type(particle_type), dimension(:), intent(in out) :: particles
     integer, dimension(:), intent(in) :: active_particles
     logical, intent(out) :: displacement_flag

     real, dimension(3) :: dx
     real :: dr, dr_max, dr_max_first, dr_max_second
     integer :: ii, n

     dr_max_first = zero
     dr_max_second = zero

     !do i = 1, size(p)
     do n = 1, size(active_particles)
       ii = active_particles(n)

       ! INIT FORCE
       ! the formerly new Force is now the old Force
       particles(ii)%force_old = particles(ii)%force_new
       ! Reset the new Force to dummy value
       particles(ii)%force_new = zero !outer_force

       ! COMPUTE LOCATION
       dx = delta_t * (particles(ii)%velocity + &
                       time_mass_factor * particles(ii)%force_old)
       particles(ii)%location = particles(ii)%location + dx
       particles(ii)%displacement = particles(ii)%displacement + dx

       ! GET MAX DISPLACEMENT OF ith PARTICLE
       dr = norm2(particles(ii)%displacement)
       if (dr > dr_max_first) then
          dr_max_second = dr_max_first
          dr_max_first = dr
       else if (dr > dr_max_second) then
          dr_max_second = dr
       end if
     end do

     ! maximum displacement as first and second max. dispalcement found
     dr_max = dr_max_first + dr_max_second

     ! if the maximal displacement becomes to large, update the neighbor list
     if (dr_max > dr_max_tolerable) then
       displacement_flag = .true.
     else
       displacement_flag = .false.
     end if

   end subroutine compute_locations

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine compute_velocities(particles, active_particles)
     type(particle_type), dimension(:), intent(in out) :: particles
     integer, dimension(:), intent(in)  :: active_particles

     integer :: n, ii

     do n = 1, size(active_particles)
        ii = active_particles(n)
        particles(ii)%velocity = particles(ii)%velocity + time_mass_factor * &
                        (particles(ii)%force_new + particles(ii)%force_old)
     end do

   end subroutine compute_velocities

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module stoermer_verlet_method

!===============================================================================

module particle_management
!
! all the procedures and types for administrating the neighbor list
! of each particle. It's very bureaucratic stuff, actually.
!
  use iso_fortran_env, only: output_unit
  use auxiliary_routines, only: stop_program
  use stoermer_verlet_method, only : particle_type, &
                                     length_of_area, min_length_of_cell, &
                                     length_of_cell, zero
  implicit none
  private

  type particle_number_type
    integer :: i = -1 ! position of particle in particle_array
    type(particle_number_type), pointer :: next_number => null()
  end type particle_number_type

  type :: cell_type
     type(particle_number_type), pointer :: first_number => null()
  end type

  type(cell_type), dimension(:,:,:), allocatable :: cell_space

  ! public procedures
  public :: initialize_neighbor_lists, update_neighbor_list, clean_up

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine initialize_neighbor_lists(particles, active_particles)
    type(particle_type), dimension(:), intent(in out) :: particles
    integer, dimension(:), allocatable, intent(out) :: active_particles

    ! create cell_space space, and ill it with the particle_numbers
    write(output_unit,*) "[INFO] initialize cell space"
    cell_space = initialze_cell_space(particles)

    ! update neighbor list for each particle in particle_array by using
    ! the cell_space space
    write(output_unit,*) "[INFO] initialize neighbor lists"
    call update_actual_neighbor_list(particles, cell_space)
    !write(*,*) "fertig"

    write(output_unit,*) "[INFO] initialize active particles"
    active_particles = update_active_particles(pack(cell_space, .true.), &
                                               size(particles))

  end subroutine initialize_neighbor_lists

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function initialze_cell_space(particles) result(cell_space)
    type(particle_type), dimension(:), intent(in) :: particles
    type(cell_type), dimension(:,:,:), allocatable :: cell_space

    integer, dimension(3) :: nc
    integer :: istat
    type(particle_number_type), pointer :: new_particle_number
    integer :: ii

    !
    nc = nint(length_of_area / length_of_cell)
    write(output_unit,*) 'Number of cells ', nc
    write(output_unit,*) length_of_area
    write(output_unit,*) size(particles)

    ! allocate cell_space
    allocate(cell_space(nc(1), nc(2), nc(3)), stat = istat)
    if (istat /= 0) then
      call stop_program('[ERROR] could not allocate linked_cells array', istat)
    end if

    ! fill the cell_space_array with particle numbers
    do ii = 1, size(particles)
      ! Create new particle
      allocate(new_particle_number, stat=istat)
      if (istat /= 0) then
        call stop_program('[ERROR] ALLOCATE failed for "new_particle_number"', &
                          istat)
      end if

      if (any(particles(ii)%location > length_of_area) .or. &
          any(particles(ii)%location < (/0.0, 0.0, 0.0/))) then

          write(output_unit,*) "[INFO] particle outside of simulation area", &
                               particles(ii)%id
      else
          ! Compute position of target cell of particle
          nc = ceiling(particles(ii)%location / length_of_cell)

          !write(*,*) i, p(i)%location, nc
          ! set the particle number
          new_particle_number%i = ii

          ! connect particle with target cell in cell_space
          new_particle_number%next_number => &
                                    cell_space(nc(1),nc(2),nc(3))%first_number
          cell_space(nc(1),nc(2),nc(3))%first_number => new_particle_number
      end if
    end do

  end function initialze_cell_space

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine update_neighbor_list(particles, active_particles)
    type(particle_type), dimension(:), intent(in out) :: particles
    integer, dimension(:), intent(out) :: active_particles

    call update_cell_space(cell_space, particles)

    call update_actual_neighbor_list(particles, cell_space)

    active_particles = update_active_particles(pack(cell_space, .true.), &
                                               size(particles))
  end subroutine update_neighbor_list

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine update_actual_neighbor_list(particles, cell_space)
    type(particle_type), dimension(:), intent(in out) :: particles
    type(cell_type), dimension(:,:,:), intent(in) :: cell_space

    call collect_neighbors_within_cell(pack(cell_space, .TRUE.), particles)

    call collect_neighbors_from_other_cells(cell_space, particles)

  end subroutine update_actual_neighbor_list

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine collect_neighbors_within_cell(cells, particles)
     type(cell_type), dimension(:), intent(in) :: cells
     type(particle_type), dimension(:), intent(in out) :: particles

     type(particle_number_type), pointer :: p1
     type(particle_number_type), pointer :: p2
     integer :: ci

     cell_loop: do ci = 1, size(cells)

        p1 => cells(ci)%first_number
        p1_loop: do while(associated(p1))

          ! RESET DISPLACEMENT AND NUMBER OF NEIGHBORS
          particles(p1%i)%displacement = zero
          particles(p1%i)%nb_of_neighbors = 0

          p2 => p1%next_number
          p2_loop: do while(associated(p2))
            call add_neighbor_particle(particles(p1%i), &
                                       particles(p2%i), &
                                       p2%i)
            p2 => p2%next_number
          end do p2_loop

          p1 => p1%next_number
      end do p1_loop

    end do cell_loop

  end subroutine collect_neighbors_within_cell

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine collect_neighbors_from_other_cells(cell_space, particles)
  !
  ! internal subroutine of subroutine compF_LC
  !
  ! Purpose
  !
  ! ----------------------------------------------------------------------------
  !      ---------------------------------
  !      !              !                !
  !      ! cx, cy+1, cz ! cx+1, cy+1, cz !
  !      !              !                !
  !      ---------------------------------
  !      !              !                !
  !      ! cx, cy, cz   ! cx+1, cy, cz   !
  !      !              !                !
  !      ---------------------------------
  !
  !      within cz Plane:
  !          (i)   cell(cx,cy,cz) with cell(cx,cy+1,cz)
  !          (ii)  cell(cx,cy,cz) with cell(cx+1,cy,cz)
  !          (iii) cell(cx,cy,cz) with cell(cx+1,cy+1,cz)
  !          (iv)  cell(cx+1,cy,cz) with cell(cx,cy+1,cz)
  !
  ! ----------------------------------------------------------------------------
  !  ---------------------------------   -------------------------------------
  !  !              !                !   !                !                  !
  !  ! cx, cy+1, cz ! cx+1, cy+1, cz !   ! cx, cy+1, cz+1 ! cx+1, cy+1, cz+1 !
  !  !              !                !   !                !                  !
  !  ---------------------------------   -------------------------------------
  !  !              !                !   !                !                  !
  !  ! cx, cy, cz   ! cx+1, cy, cz   !   ! cx, cy, cz+1   ! cx+1, cy, cz+1   !
  !  !              !                !   !                !                  !
  !  ---------------------------------   -------------------------------------
  !
  !      cz with cz+1 Plane:
  !            cz plane                cz+1 plane
  !      (1)   cell(cx,cy,cz)     with cell(cx,cy,cz+1)
  !      (2.a) cell(cx,cy,cz)     with cell(cx,cy+1,cz+1)
  !      (2.b) cell(cx,cy+1,cz)   with cell(cx,cy,cz+1)
  !      (3.a) cell(cx,cy,cz)     with cell cell(cx+1,cy,cz+1)
  !      (3.b) cell(cx+1,cy,cz)   with cell(cx,cy,cz+1)
  !      (4.a) cell(cx,cy,cz)     with cell(cx+1,cy+1,cz+1)
  !      (4.b) cell(cx+1,cy+1,cz) with cell(cx,cy,cz+1)
  !      (5.a) cell(cx+1,cy,cz)   with cell(cx,cy+1,cz+1)
  !      (5.b) cell(cx,cy+1,cz)   with cell(cx+1,cy,cz+1)
  !
  ! ----------------------------------------------------------------------------
     type(cell_type), dimension(:,:,:), intent(in) :: cell_space
     type(particle_type), dimension(:), intent(in out) :: particles

     integer :: cx, cx_max
     integer :: cy, cy_max
     integer :: cz, cz_max

     cx_max = size(cell_space, dim = 1)
     cy_max = size(cell_space, dim = 2)
     cz_max = size(cell_space, dim = 3)

     !---cz Ebene---------------------------------------------------------------
     z_loop_a: do cz = 1, cz_max

        ! (i) cell(cx,cy,cz) with cell(cx,cy+1,cz)
        do cy = 1, cy_max - 1
           do cx = 1, cx_max
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx, cy+1, cz), &
                                                 particles)
           end do
        end do

        ! (ii) cell(cx,cy,cz) with cell(cx+1,cy,cz)
        do cy = 1, cy_max
           do cx = 1, cx_max - 1
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx+1, cy, cz), &
                                                 particles)
           end do
        end do

        ! (iii) cell(cx,cy,cz)   with cell(cx+1,cy+1,cz)
        ! (iv)  cell(cx+1,cy,cz) with cell(cx,cy+1,cz)
        do cy = 1, cy_max - 1
           do cx = 1, cx_max - 1
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx+1, cy+1, cz), &
                                                 particles)
              call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
                                                 cell_space(cx, cy+1, cz), &
                                                 particles)
           end do
        end do

     end do z_loop_a
     !--------------------------------------------------------------------------

     !---cz mit cz+1 Ebene------------------------------------------------------
     z_loop_b: do cz = 1, cz_max - 1

        ! (1) cell(cx,cy,cz) with cell(cx,cy,cz+1)
        do cy = 1, cy_max
           do cx = 1, cx_max
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx, cy, cz+1), &
                                                 particles)
           end do
        end do

        ! (2)
        do cy = 1, cy_max - 1
           do cx = 1, cx_max
              ! (2.a) cell(cx,cy,cz) with cell(cx,cy+1,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx, cy+1, cz+1), &
                                                 particles)
              ! (2.b) cell(cx,cy+1,cz) with cell(cx,cy,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx, cy+1, cz), &
                                                 cell_space(cx, cy, cz+1), &
                                                 particles)
           end do
        end do

        ! (3)
        do cy = 1, cy_max
           do cx = 1, cx_max - 1
              ! (3.a) cell(cx,cy,cz) with cell(cx+1,cy,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx+1, cy, cz+1), &
                                                 particles)
              ! (3.b) cell(cx+1,cy,cz) with cell(cx,cy,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
                                                 cell_space(cx, cy, cz+1), &
                                                 particles)
           end do
        end do

        ! (4) and (5)
        do cy = 1, cy_max - 1
           do cx = 1, cx_max - 1
              ! (4.a) cell(cx,cy,cz) with cell(cx+1,cy+1,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx+1, cy+1, cz+1), &
                                                 particles)
              ! (4.b) cell(cx+1,cy+1,cz) with cell(cx,cy,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx+1, cy+1, cz), &
                                                 cell_space(cx, cy, cz+1), &
                                                 particles)
              ! (5.a) cell(cx+1,cy,cz) with cell(cx,cy+1,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
                                                 cell_space(cx, cy+1, cz+1), &
                                                 particles)
              ! (5.b) cell(cx,cy+1,cz) with cell(cx+1,cy,cz)
              call add_neighbors_from_other_cell(cell_space(cx, cy+1, cz), &
                                                 cell_space(cx+1, cy, cz+1), &
                                                 particles)
           end do
        end do

     end do z_loop_b
     !--------------------------------------------------------------------------

  end subroutine collect_neighbors_from_other_cells

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine add_neighbors_from_other_cell(this_cell, other_cell, &
                                           particles)
     type(cell_type), intent(in) :: this_cell
     type(cell_type), intent(in) :: other_cell
     type(particle_type), dimension(:), intent(in out) :: particles

     type(particle_number_type), pointer :: p1
     type(particle_number_type), pointer :: p2

     if (.not. associated(other_cell%first_number)) return

     p1 => this_cell%first_number
     p1_loop: do while(associated(p1))

        p2 => other_cell%first_number
        p2_loop: do while(associated(p2))
           call add_neighbor_particle(particles(p1%i), &
                                      particles(p2%i), &
                                      p2%i)
           p2 => p2%next_number
        end do p2_loop

        p1 => p1%next_number
     end do p1_loop

  end subroutine add_neighbors_from_other_cell

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine add_neighbor_particle(p1, p2, index_p2)
  !
     type(particle_type), intent(in out) :: p1
     type(particle_type), intent(in) :: p2
     integer, intent(in) :: index_p2

     real :: r

     r = norm2(p1%location - p2%location)

     if (r <= min_length_of_cell) then
        ! increase the counter for the number of neighbors
        p1%nb_of_neighbors = p1%nb_of_neighbors + 1
        ! check if neighbor list needs to be elongated
        if (p1%nb_of_neighbors > size(p1%neighbor_list)) then
           ! increase size of neighbor_list by one
           call increase_neighbor_list_size(p1%neighbor_list)
        end if
        ! add index of neighbor particle to the list
        p1%neighbor_list(p1%nb_of_neighbors) = index_p2
     end if

   end subroutine add_neighbor_particle

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  pure subroutine increase_neighbor_list_size(neighbor_list)
     integer, dimension(:), allocatable, intent(in out) :: neighbor_list

     integer, dimension(:), allocatable :: temp
     integer :: n

     n = size(neighbor_list)
     allocate(temp(n+10))
     temp(1:n) = neighbor_list
     call move_alloc(from = temp, to = neighbor_list)

  end subroutine increase_neighbor_list_size

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine update_cell_space(linked_cells, particles)
     type(cell_type), dimension(:,:,:), intent(in out) :: linked_cells
     type(particle_type), dimension(:), intent(in) :: particles

     ! pointer the the current particle number, and the previous particle number
     type(particle_number_type), pointer :: p_number, prev_p_number
     integer :: cx, cy, cz
     integer, dimension(3) :: new_cell

     z_loop: do cz = 1, size(linked_cells, dim = 3)
       y_loop: do cy = 1, size(linked_cells, dim = 2)
         x_loop: do cx = 1, size(linked_cells, dim = 1)

           ! check if no particles are in this cell
           if (.not. associated(linked_cells(cx, cy, cz)%first_number)) then
             ! no particles, go for the next cell
             cycle x_loop
           end if

           ! init
           p_number => linked_cells(cx, cy, cz)%first_number
           prev_p_number => null()

           particle_loop: do
             ! compute coordinates of new cell of particle
             new_cell = ceiling(particles(p_number%i)%location / &
                                length_of_cell)

             whereto: if (all(new_cell == (/cx, cy, cz/))) then
               ! Particle does not left its cell
               if (associated(p_number%next_number)) then
                 ! go to next particle
                 prev_p_number => p_number
                 p_number => p_number%next_number
               else
                 ! end of particle list
                 exit particle_loop
               end if
             else whereto
               ! Particle left its cell.
               ! Unhinge particle from its cell.
               ! Hinge particle to its new cell.
               ! Change the pointer settings.
               follower: if (associated(p_number%next_number)) then
                 if (associated(prev_p_number)) then
                   ! Particle has predecessor and successor
                   prev_p_number%next_number => p_number%next_number
                   !call entry_in_new_list(particle, new_cell, &
                  !                        particle_buffer, linked_cells)
                   call entry_in_new_list(p_number, new_cell, linked_cells)
                   p_number => prev_p_number%next_number
                 else
                   ! Particle is first particle in list
                   linked_cells(cx, cy, cz)%first_number => p_number%next_number
                   call entry_in_new_list(p_number, new_cell, linked_cells)
                   p_number => linked_cells(cx, cy, cz)%first_number
                 end if
               else follower
                 if (associated(prev_p_number)) then
                   ! Particle is the last one in list
                   prev_p_number%next_number => null()
                 else
                   ! Particle is the first and last one in list
                   linked_cells(cx, cy, cz)%first_number => null()
                 end if
                 call entry_in_new_list(p_number, new_cell, linked_cells)
                 ! no more particles in list
                 exit particle_loop
               end if follower
             end if whereto
           end do particle_loop

         end do x_loop
       end do y_loop
     end do z_loop

  end subroutine update_cell_space

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine entry_in_new_list(p_number, new_cell, linked_cells)
     type(particle_number_type), pointer, intent(in out) :: p_number
     integer, dimension(3), intent(in) :: new_cell
     type(cell_type), dimension(:,:,:), intent(in out) :: linked_cells

     integer, dimension(3) :: max_cell
     integer :: x, y, z

     max_cell = shape(linked_cells)

     if (any(new_cell < (/1,1,1/)) .or. any(new_cell > max_cell)) then
       call delete_particle_number(p_number)
     else
       x = new_cell(1)
       y = new_cell(2)
       z = new_cell(3)
       p_number%next_number => linked_cells(x,y,z)%first_number
       linked_cells(x,y,z)%first_number => p_number
     end if

   end subroutine entry_in_new_list

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine delete_particle_number(p_number)
     type(particle_number_type), pointer, intent(in out) :: p_number

     type(particle_number_type), pointer :: tmp
     integer :: istat

     tmp => p_number
     deallocate(tmp, stat = istat)
     if (istat /= 0) then
       call stop_program("[ERROR] Deallocate failed for deleting "// &
                         "particle number.", istat)
     end if

   end subroutine delete_particle_number

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function update_active_particles(cell_space, max_nb_active_particles) &
                                  result(active_particles)
    type(cell_type), dimension(:), intent(in) :: cell_space
    integer, intent(in) :: max_nb_active_particles
    integer, dimension(:), allocatable :: active_particles

    integer, dimension(max_nb_active_particles) :: temp
    type(particle_number_type), pointer :: number
    integer :: ci, n

    n = 0
    do ci = 1, size(cell_space)
      number => cell_space(ci)%first_number
      do while(associated(number))
        n = n + 1
        temp(n) = number%i
        number => number%next_number
      end do
    end do

    allocate(active_particles(n))

    if (n > 0) active_particles = temp(1:n)

  end function update_active_particles

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine clean_up(particles, active_particles)
    type(particle_type), dimension(:), allocatable, intent(in out) :: particles
    integer, dimension(:), allocatable, intent(in out) :: active_particles

    type(particle_number_type), pointer :: p_number
    integer :: cz, cy, cx
    integer :: istat

    deallocate(particles, stat=istat)
    if (istat /= 0) then
       call stop_program("[ERROR] Failed deallocating particle_array.", istat)
    end if

    deallocate(active_particles, stat=istat)
    if (istat /= 0) then
       call stop_program("[ERROR] Failed deallocating active article array.", &
                         istat)
    end if

    ! remove elements connected to the cell_space
    z_loop: do cz = 1, size(cell_space, dim=3)
       y_loop: do cy = 1, size(cell_space, dim=2)
          x_loop: do cx = 1, size(cell_space, dim=1)

             particle_loop: do
                p_number => cell_space(cx, cy, cz)%first_number
                if (.not. associated(p_number)) exit particle_loop

                cell_space(cx, cy, cz)%first_number => &
                            cell_space(cx, cy, cz)%first_number%next_number

                deallocate(p_number, stat=istat)
                if (istat /= 0) then
                   call stop_program("[ERROR] Failed deallocating particle. ",&
                                     istat)
               end if

             end do particle_loop

          end do x_loop
       end do y_loop
    end do z_loop

    ! deallocate the empty cell_space array
    deallocate(cell_space, stat=istat)
    if (istat /= 0) then
       call stop_program("[ERROR] Failed deallocating cell_space.", istat)
    end if

  end subroutine clean_up

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module particle_management

!===============================================================================

program particle_simulation
  use iso_fortran_env, only: output_unit
  use auxiliary_routines
  use stoermer_verlet_method, only: particle_type, &
                                    set_simulation_parameter, &
                                    compute_forces, compute_locations, &
                                    compute_velocities
  use particle_management, only: initialize_neighbor_lists, &
                                 update_neighbor_list, clean_up
  implicit none

  real, parameter :: zero = 0.0

  ! contains all particles
  type(particle_type), dimension(:), allocatable :: particles
  ! contains the index for all active partciles,
  ! i.e. all particles that are inside the simulation area and can be moved
  integer, dimension(:), allocatable :: active_particles

  real, dimension(3) :: length_of_area, length_of_cell
  real :: start_time

  integer :: particle_file_unit, environment_file_unit, output_file_unit
  integer :: time_step, max_time_step, write_time_step
  integer :: i
  logical :: displacement_flag


  ! get all the files from the command link
  write(output_unit,*) "[INFO] get command arguments"
  call get_command_arguments(environment_file_unit, &
                             particle_file_unit, &
                             output_file_unit)

  ! read fixed values from environment_file
  write(output_unit,*) "[INFO] get and set simulation parameters"
  call get_and_set_parameters(environment_file_unit, &
                              length_of_area, length_of_cell, &
                              write_time_step, max_time_step)

  write(output_unit,*) "[INFO] write size of simulation area to output file"
  do i = 1, 3
    write(unit = output_file_unit, rec = i) length_of_area(i)
  end do

  ! read particle file and file particle_array with particles and their id,
  ! mass, location, and initial velocity
  write(output_unit,*) "[INFO] initialze particle array"
  particles = initialize_particle_array(particle_file_unit)

  ! set neighbors for each particle and get the index of the active particles
  ! in the array particles
  write(output_unit,*) "[INFO] initialize neighborlists, and set active particles"
  call initialize_neighbor_lists(particles, active_particles)

  ! compute the initial force for all active particles in particle list by using
  ! their neighborlist
  write(output_unit,*) "[INFO] initialze forces"
  call compute_forces(particles, active_particles)

  write(output_unit,*) "[INFO] start time integration"
  start_time = get_current_time()
  ! start loop for time integration
  time_integration: do time_step = 1, max_time_step

    ! compute locations for each particle in particle list
    call compute_locations(particles, active_particles, displacement_flag)

    ! if the displacment is to large, then the neighbor_list and the
    ! cell_space for each particle needs to be updated
    if (displacement_flag .eqv. .true.) then
      call update_neighbor_list(particles, active_particles)
      if (size(active_particles) == 0) then
        write(output_unit,*) "[INFO] no particles left in simulation area."
        exit time_integration
      end if
    end if

    ! compute forces for each particle in the particle list using their
    ! neighbor list
    call compute_forces(particles, active_particles)

    ! compute velocities for each particle
    call compute_velocities(particles, active_particles)

    ! if it is time to write data, write data
    if (write_data_now(time_step)) then
      call write_data(output_file_unit, particles, active_particles)
    end if

  end do time_integration
  call write_time(get_current_time() - start_time)

  write(output_unit,*) "[INFO] clean up"
  call clean_up(particles, active_particles)

  ! close output file
  close(output_file_unit)

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
           write(error_unit, '(1x,a)') '[ERROR] Not enough arguments!'
        else
           write(error_unit, '(1x,a)') '[ERROR] To many arguments!'
        end if
        call get_command_argument(0, program_name)
        call stop_program('[ERROR] Usage is: '//trim(program_name)// &
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

  subroutine get_and_set_parameters(efu, &
                                    length_of_area, length_of_cell, &
                                    write_time_step, max_time_step)
     integer, intent(in) :: efu ! environment_file_unit
     real, dimension(3), intent(out) :: length_of_area, length_of_cell
     !real, intent(out) :: min_length_of_cell
     integer, intent(out) :: write_time_step, max_time_step

     real :: epsilon_value, delta_t
     real :: particle_mass, particle_radius, r_cut_factor

     read(efu, *) particle_radius, epsilon_value, particle_mass
     read(efu, *) r_cut_factor
     read(efu, *) length_of_area
     write(output_unit,*) length_of_area
     read(efu, *) length_of_cell
     read(efu, *) delta_t
     read(efu, *) write_time_step, max_time_step


     !min_length_of_cell = minval(length_of_cells)

     call set_simulation_parameter(delta_t, length_of_area, length_of_cell, &
                                   epsilon_value, particle_radius, &
                                   particle_mass, r_cut_factor)

     ! close environment file
     close(efu)

  end subroutine get_and_set_parameters

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function initialize_particle_array(pfu) result(p)
    integer, intent(in) :: pfu !particle_file_unit
    type(particle_type), dimension(:), allocatable :: p !particle_array

    type(particle_type), dimension(:), allocatable :: temp
    integer :: id
    integer :: n, istat
    real, dimension(3) :: location, velocity

    ! write particle data from particle_file to particle_array
    n = 0
    allocate(temp(1000))
    do
      ! read particle data
      read(pfu, *, iostat = istat) id, location, velocity
      if (istat /= 0) exit ! end of file or somthing is wrong
      n = n + 1
      if (n > size(temp)) call increase_particle_array(temp)
      ! write particle data to temp particle array
      temp(n)%id = id
      temp(n)%location = location
      temp(n)%velocity = velocity
      allocate(temp(n)%neighbor_list(10))
    end do

    if (n == 0) then
      call stop_program("[ERROR] No particles in particle file.")
    else
      write(output_unit,*) "[INFO] Number of particles in file :", n
    end if

    ! create actual particle array
    allocate(p(n))
    p = temp(1:n)
    deallocate(temp)

    ! close particle file
    close(pfu)

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

  function write_data_now(time_step) result(write_now)
    integer, intent(in) :: time_step
    logical :: write_now

    integer, save :: write_step = 1

    if (time_step == write_step) then
      write(output_unit,*) "[INFO] write now for time step ", time_step
      write_step = write_step + write_time_step
      write_now = .true.
    else
      write_now = .false.
    end if

  end function write_data_now

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine write_data(ofu, particles, active_particles)
     integer, intent(in) :: ofu ! output_file_unit
     type(particle_type), dimension(:), intent(in) :: particles
     integer, dimension(:), intent(in) :: active_particles

     integer, save :: counter = 4, number_of_entries = 0
     integer :: d
     integer :: jj, ii, n, m

     m = size(active_particles)
     write(output_unit,*) "[INFO] Active particles: ", m

     ! number of entries
     number_of_entries = number_of_entries + 1
     write(unit=ofu, rec=4) number_of_entries

     counter = counter + 1
     d = counter
     do n = 1, m
       ii = active_particles(n)

       d = d + 1
       write(unit=ofu, rec=d) particles(ii)%id

       ! write location of particle
       do jj = 1, 3
         d = d + 1
         write(unit=ofu, rec=d) particles(ii)%location(jj)
       end do

       ! write velocity of particle
       do jj = 1, 3
         d = d + 1
         write(unit=ofu, rec=d) particles(ii)%velocity(jj)
       end do
     end do

     ! write number of data per time entry
     write(unit=ofu, rec=counter) d - counter

     counter = d

   end subroutine write_data

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end program particle_simulation
