module stoermer_verlet_method
!
!
!
  use, intrinsic :: iso_fortran_env, only: rt => real32
  implicit none
  private

  real(rt), parameter :: zero   = 0.0_rt
  real(rt), parameter :: half   = 0.5_rt
  real(rt), parameter :: one    = 1.0_rt
  real(rt), parameter :: two    = 2.0_rt
  real(rt), parameter :: six    = 6.0_rt
  real(rt), parameter :: twelve = 12.0_rt

  type particle_type
    integer :: id = -1
    !real(rt) :: epsilon, r_min 
    real(rt), dimension(3) :: location     = zero
    real(rt), dimension(3) :: displacement = zero
    real(rt), dimension(3) :: velocity     = zero
    real(rt), dimension(3) :: force_old    = zero
    real(rt), dimension(3) :: force_new    = zero
    logical :: is_active = .true. ! is particle inside of simulation area
    integer, dimension(:), allocatable :: neighbor_list
    integer :: nb_of_neighbors = 0
  end type particle_type

  ! constants to be set with the start of the simulation
  real(rt) :: time_mass_factor, delta_t
  real(rt) :: dr_max_tolerable
  real(rt) :: r_min, r_cut
  real(rt) :: F_rcut
  real(rt) :: epsilon_value
  real(rt) :: min_length_of_cell
  real(rt), dimension(3) :: length_of_area, length_of_cell

  ! visible to the outside of module
  ! constants and variables
  public :: length_of_area, min_length_of_cell, length_of_cell
  ! types
  public :: particle_type
  ! procedures
  public :: set_simulation_parameter
  public :: compute_forces, compute_locations, compute_velocities!, &
            !need_to_update_neighbor_list

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_simulation_parameter(input_delta_t, input_length_of_area, &
                                      input_length_of_cell, &
                                      input_epsilon_value, particle_radius, &
                                      particle_mass, r_cut_factor)
    real(rt), intent(in) :: input_delta_t
    real(rt), dimension(3), intent(in) :: input_length_of_area, input_length_of_cell
    real(rt), intent(in) :: input_epsilon_value
    real(rt), intent(in) :: particle_radius
    real(rt), intent(in) :: particle_mass
    real(rt), intent(in) :: r_cut_factor

    real(rt) :: s

    ! set the maximal acceptable displacement value (if displacement of one
    ! particle is larger than this value, the neighbor list for each particle
    ! needs to be updated)
    delta_t            = input_delta_t
    epsilon_value      = input_epsilon_value
    r_min              = two**(one/six) * particle_radius
    r_cut              = r_cut_factor * particle_radius
    length_of_area     = input_length_of_area
    length_of_cell     = input_length_of_cell
    min_length_of_cell = minval(input_length_of_cell)
    dr_max_tolerable   = min_length_of_cell - r_cut

    ! set fixed values for particle simulation
    s = (r_min / r_cut)**6
    F_rcut = (s / (r_cut * r_cut)) * (one - s)
    !
    time_mass_factor = half * delta_t / particle_mass

  end subroutine set_simulation_parameter

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  pure subroutine compute_forces(particles, pindex)
    type(particle_type), dimension(:), intent(in out):: particles
    integer, dimension(2), intent(in) :: pindex

    integer :: i, j, n

    particle_loop: do i = pindex(1), pindex(2)
      if (particles(i)%is_active .and. particles(i)%nb_of_neighbors > 0) then
        neighbor_loop: do n = 1, particles(i)%nb_of_neighbors
          j = particles(i)%neighbor_list(n)
          call compute_force_between_two_particles(particles(i), &
                                                   particles(j))
        end do neighbor_loop
      end if
    end do particle_loop

   end subroutine compute_forces

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  pure subroutine compute_force_between_two_particles(particle_i, particle_j)
     type(particle_type), intent(in out) :: particle_i, particle_j

     real(rt), dimension(3) :: r_ij
     real(rt), dimension(3) :: force
     real(rt) :: r
     real(rt) :: s
     real(rt) :: F_rij

     ! Direction Vector (Richtungsvektor)
     r_ij = particle_i%location - particle_j%location
     ! Distance between both Particles
     r = norm2(r_ij)

     if (r < r_cut) then
       ! FORCE FOR TRUNCATED LENARD-JONES POTENTIAL
       s = (r_min / r)**6
       F_rij = (s / (r * r)) * (one - s)
       ! TRUNCATED FORCE
       force = r_ij * twelve * epsilon_value * (F_rij - F_rcut)
       ! add forces
       particle_i%force_new = particle_i%force_new - force
       particle_j%force_new = particle_j%force_new + force
     end if

   end subroutine compute_force_between_particles

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function compute_locations(particles, pindex) result (displacement_flag)
     type(particle_type), dimension(:), intent(in out) :: particles
     integer, dimension(2), intent(in) :: pindex
     logical :: displacement_flag

     real(rt), dimension(3) :: dx
     real(rt) :: dr, dr_max, dr_max_first, dr_max_second
     integer :: i

     dr_max_first  = zero
     dr_max_second = zero

     do i = pindex(1), pindex(2)
       if (.not. particles(i)%is_active) cycle

       ! INIT FORCE
       ! the formerly new Force is now the old Force
       particles(i)%force_old = particles(i)%force_new
       ! Reset the new Force to dummy value
       particles(i)%force_new = zero !outer_force

       ! COMPUTE LOCATION
       dx = delta_t * (particles(i)%velocity + time_mass_factor * particles(i)%force_old)
       particles(i)%location     = particles(i)%location + dx
       particles(i)%displacement = particles(i)%displacement + dx

       ! GET MAX DISPLACEMENT OF PARTICLE
       dr = norm2(particles(i)%displacement)
       if (dr > dr_max_first) then
          dr_max_second = dr_max_first
          dr_max_first  = dr
       else if (dr > dr_max_second) then
          dr_max_second = dr
       end if
     end do

     dr_max = dr_max_first + dr_max_second

     ! if the maximal displacement becomes to large, 
     ! set flag for updating the neighbor list
     if (dr_max > dr_max_tolerable) then
       displacement_flag = .true.
     else
       displacement_flag = .false.
     end if

  !end subroutine compute_locations
  end function compute_locations

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  pure subroutine compute_velocities(particles, pindex)
    type(particle_type), dimension(:), intent(in out) :: particles
    integer, dimension(2), intent(in) :: pindex

    integer :: i

    do i = pindex(1), pindex(2)
      if (particles(i)%is_active) then
        particles(i)%velocity = particles(i)%velocity + time_mass_factor * &
                               (particles(i)%force_new + particles(i)%force_old)
      end if
    end do

  end subroutine compute_velocities
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module stoermer_verlet_method