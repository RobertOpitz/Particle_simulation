module particle_management
!
! all the procedures and types for administrating the neighbor list
! of each particle. It's very bureaucratic stuff, actually.
!
  use, intrinsic :: iso_fortran_env, only: rt => real32
  use auxiliary_routines, only: stop_program
  use stoermer_verlet_method, only : particle_type, &
                                     length_of_area, min_length_of_cell, &
                                     length_of_cell
  implicit none
  private

  type particle_number
    integer :: i = -1 ! position of particle in particle_array
    type(particle_number), pointer :: next_number => null()
  end type particle_number

  type :: cell_type
     type(particle_number), pointer :: first_number => null()
  end type

  type(cell_type), dimension(:,:,:), allocatable :: cell_space

  ! public procedures
  public :: initialize_neighbor_lists, update_neighbor_list_and_cell_space, &
            clean_up

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine initialize_neighbor_lists(particles)
    type(particle_type), dimension(:), intent(in out) :: particles

    ! create cell_space space, and ill it with the particle_numbers
    write(*,'(1x,a)') "[INFO] initialize cell space"
    cell_space = initialze_cell_space(particles)

    ! update neighbor list for each particle in particle_array by using
    ! the cell_space space
    write(*,'(1x,a)') "[INFO] initialize neighbor lists"
    call update_neighbor_list(particles, cell_space)

    write(*,'(1x,a)') "[INFO] initialize active particles"
    call update_active_particles(pack(cell_space, .true.), particles)

  end subroutine initialize_neighbor_lists

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function initialze_cell_space(particles) result(cell_space)
    type(particle_type), dimension(:), intent(in) :: particles !particle_array
    type(cell_type), dimension(:,:,:), allocatable :: cell_space

    integer, dimension(3) :: nc
    integer :: istat
    type(particle_number), pointer :: new_particle_number
    integer :: i

    !
    nc = nint(length_of_area / length_of_cell)
    write(*,*) 'Number of cells ', nc
    write(*,*) length_of_area
    write(*,*) size(particles)

    ! allocate cell_space
    allocate(cell_space(nc(1), nc(2), nc(3)), stat = istat)
    if (istat /= 0) then
      call stop_program('[ERROR] could not allocate linked_cells array', istat)
    end if

    ! fill the cell_space_array with particle numbers
    do i = 1, size(particles)
      ! Create new particle
      allocate(new_particle_number, stat=istat)
      if (istat /= 0) then
        call stop_program('[ERROR] ALLOCATE failed for "new_particle_number"', &
                          istat)
      end if

      if (any(particles(i)%location > length_of_area) .or. &
          any(particles(i)%location < (/0.0, 0.0, 0.0/))) then

          write(*,'(1x,a,1x,g0)') "[INFO] particle outside of simulation area:",&
                                  particles(i)%id

          ! delete this particles
          deallocate(new_particle_number, stat = istat)
          if (istat /= 0) then
            call stop_program("[ERROR] Failed deallocating particle. ", istat)
          end if
      else
          ! Compute position of target cell of particle
          nc = ceiling(particles(i)%location / length_of_cell)

          ! Initialize new_particle_number
          ! set the particle number
          ! connect particle with target cell in cell_space
          new_particle_number = particle_number(i, &
                                     cell_space(nc(1),nc(2),nc(3))%first_number)

          ! connected new new_particle_number to the linked list
          cell_space(nc(1),nc(2),nc(3))%first_number => new_particle_number
      end if
    end do

  end function initialze_cell_space

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine update_neighbor_list_and_cell_space(particles)
    !type(particle_array_type), intent(in out) :: particles
    type(particle_type), dimension(:), intent(in out) :: particles

    !write(*,*) length_of_area, min_length_of_cell, length_of_cell 

    call update_cell_space(cell_space, particles)

    call update_neighbor_list(particles, cell_space)

    call update_active_particles(pack(cell_space, .true.), particles)

  end subroutine update_neighbor_list_and_cell_space

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine update_neighbor_list(particle_array, cell_space)
    type(particle_type), dimension(:), intent(in out) :: particle_array
    type(cell_type), dimension(:,:,:), intent(in) :: cell_space

    call collect_neighbors_within_cell(pack(cell_space, .TRUE.), &
                                       particle_array)

    call collect_neighbors_from_other_cells(cell_space, particle_array)

  end subroutine update_neighbor_list

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine collect_neighbors_within_cell(cells, particle_array)
     type(cell_type), dimension(:), intent(in) :: cells
     type(particle_type), dimension(:), intent(in out) :: particle_array

     type(particle_number), pointer :: p1
     type(particle_number), pointer :: p2
     integer :: ci

     real(rt), parameter :: zero = 0.0_rt     

     cell_loop: do ci = 1, size(cells)

      p1 => cells(ci)%first_number
      p1_loop: do !while(associated(p1))
        if (.not. associated(p1)) exit p1_loop

        ! RESET DISPLACEMENT AND NUMBER OF NEIGHBORS
        particle_array(p1%i)%displacement    = zero
        particle_array(p1%i)%nb_of_neighbors = 0
        
        p2 => p1%next_number
        p2_loop: do !while(associated(p2))
          if (.not. associated(p2)) exit p2_loop
          call add_neighbor_particle(particle_array(p1%i), &
                                     particle_array(p2%i), &
                                     p2%i)
          p2 => p2%next_number
        end do p2_loop

        p1 => p1%next_number
      end do p1_loop

    end do cell_loop

  end subroutine collect_neighbors_within_cell

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine collect_neighbors_from_other_cells(cell_space, particle_array)
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
     type(particle_type), dimension(:), intent(in out) :: particle_array

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
                                                 particle_array)
           end do
        end do

        ! (ii) cell(cx,cy,cz) with cell(cx+1,cy,cz)
        do cy = 1, cy_max
           do cx = 1, cx_max - 1
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx+1, cy, cz), &
                                                 particle_array)
           end do
        end do

        ! (iii) cell(cx,cy,cz)   with cell(cx+1,cy+1,cz)
        ! (iv)  cell(cx+1,cy,cz) with cell(cx,cy+1,cz)
        do cy = 1, cy_max - 1
           do cx = 1, cx_max - 1
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx+1, cy+1, cz), &
                                                 particle_array)
              call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
                                                 cell_space(cx, cy+1, cz), &
                                                 particle_array)
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
                                                 particle_array)
           end do
        end do

        ! (2)
        do cy = 1, cy_max - 1
           do cx = 1, cx_max
              ! (2.a) cell(cx,cy,cz) with cell(cx,cy+1,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx, cy+1, cz+1), &
                                                 particle_array)
              ! (2.b) cell(cx,cy+1,cz) with cell(cx,cy,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx, cy+1, cz), &
                                                 cell_space(cx, cy, cz+1), &
                                                 particle_array)
           end do
        end do

        ! (3)
        do cy = 1, cy_max
           do cx = 1, cx_max - 1
              ! (3.a) cell(cx,cy,cz) with cell(cx+1,cy,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx+1, cy, cz+1), &
                                                 particle_array)
              ! (3.b) cell(cx+1,cy,cz) with cell(cx,cy,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
                                                 cell_space(cx, cy, cz+1), &
                                                 particle_array)
           end do
        end do

        ! (4) and (5)
        do cy = 1, cy_max - 1
           do cx = 1, cx_max - 1
              ! (4.a) cell(cx,cy,cz) with cell(cx+1,cy+1,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx, cy, cz), &
                                                 cell_space(cx+1, cy+1, cz+1), &
                                                 particle_array)
              ! (4.b) cell(cx+1,cy+1,cz) with cell(cx,cy,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx+1, cy+1, cz), &
                                                 cell_space(cx, cy, cz+1), &
                                                 particle_array)
              ! (5.a) cell(cx+1,cy,cz) with cell(cx,cy+1,cz+1)
              call add_neighbors_from_other_cell(cell_space(cx+1, cy, cz), &
                                                 cell_space(cx, cy+1, cz+1), &
                                                 particle_array)
              ! (5.b) cell(cx,cy+1,cz) with cell(cx+1,cy,cz)
              call add_neighbors_from_other_cell(cell_space(cx, cy+1, cz), &
                                                 cell_space(cx+1, cy, cz+1), &
                                                 particle_array)
           end do
        end do

     end do z_loop_b
     !--------------------------------------------------------------------------

  end subroutine collect_neighbors_from_other_cells

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine add_neighbors_from_other_cell(this_cell, other_cell, &
                                           particle_array)
     type(cell_type), intent(in) :: this_cell
     type(cell_type), intent(in) :: other_cell
     type(particle_type), dimension(:), intent(in out) :: particle_array

     type(particle_number), pointer :: p1
     type(particle_number), pointer :: p2

     if (.not. associated(other_cell%first_number)) return

     p1 => this_cell%first_number
     p1_loop: do !while(associated(p1))
        if (.not. associated(p1)) exit p1_loop

        p2 => other_cell%first_number
        p2_loop: do !while(associated(p2))
           if (.not. associated(p2)) exit p2_loop
           call add_neighbor_particle(particle_array(p1%i), &
                                      particle_array(p2%i), &
                                      p2%i)
           p2 => p2%next_number
        end do p2_loop

        p1 => p1%next_number
     end do p1_loop

  end subroutine add_neighbors_from_other_cell

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine add_neighbor_particle(p1, p2, p2_id)
  !
     type(particle_type), intent(in out) :: p1
     type(particle_type), intent(in) :: p2
     integer, intent(in) :: p2_id

     real(rt) :: r

     r = norm2(p1%location - p2%location)

     if (r <= min_length_of_cell) then

        p1%nb_of_neighbors = p1%nb_of_neighbors + 1

        if (p1%nb_of_neighbors > size(p1%neighbor_list)) then
           ! increase size of neighbor_list by one
           call increase_neighbor_list_size(p1%neighbor_list)
        end if

        p1%neighbor_list(p1%nb_of_neighbors) = p2_id

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

  subroutine update_cell_space(linked_cells, particle_array)
     type(cell_type), dimension(:,:,:), intent(in out) :: linked_cells
     type(particle_type), dimension(:), intent(in) :: particle_array

     type(particle_number), pointer :: particle, previous_particle
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
           particle => linked_cells(cx, cy, cz)%first_number
           previous_particle => null()

           particle_loop: do
             ! compute coordinates of new cell of particle
             new_cell = ceiling(particle_array(particle%i)%location / &
                                length_of_cell)

             whereto: if (all(new_cell == (/cx, cy, cz/))) then
               ! Particle does not left its cell
               if (associated(particle%next_number)) then
                 ! go to next particle
                 previous_particle => particle
                 particle => particle%next_number
               else
                 ! end of particle list
                 exit particle_loop
               end if
             else whereto
               ! Particle left its cell.
               ! Unhinge particle from its cell.
               ! Hinge particle to its new cell.
               ! Change the pointer settings.
               follower: if (associated(particle%next_number)) then
                 if (associated(previous_particle)) then
                   ! Particle has predecessor and successor
                   previous_particle%next_number => particle%next_number
                   call entry_in_new_list(particle, new_cell, linked_cells)
                   particle => previous_particle%next_number
                 else
                   ! Particle is first particle in list
                   linked_cells(cx, cy, cz)%first_number => particle%next_number
                   call entry_in_new_list(particle, new_cell, linked_cells)
                   particle => linked_cells(cx, cy, cz)%first_number
                 end if
               else follower
                 if (associated(previous_particle)) then
                   ! Particle is the last one in list
                   previous_particle%next_number => null()
                 else
                   ! Particle is the first and last one in list
                   linked_cells(cx, cy, cz)%first_number => null()
                 end if
                 call entry_in_new_list(particle, new_cell, linked_cells)
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

  subroutine entry_in_new_list(particle, new_cell, linked_cells)
     type(particle_number), pointer, intent(in out) :: particle
     integer, dimension(3), intent(in) :: new_cell
     type(cell_type), dimension(:,:,:), intent(in out) :: linked_cells

     integer, dimension(3) :: max_cell
     integer :: x, y, z

     max_cell = shape(linked_cells)

     if (any(new_cell < (/1,1,1/)) .or. any(new_cell > max_cell)) then
       call delete_particle_number(particle)
     else
       x = new_cell(1)
       y = new_cell(2)
       z = new_cell(3)
       particle%next_number => linked_cells(x,y,z)%first_number
       linked_cells(x,y,z)%first_number => particle
     end if

   end subroutine entry_in_new_list

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine delete_particle_number(particle)
     type(particle_number), pointer, intent(in out) :: particle

     type(particle_number), pointer :: tmp
     integer :: istat

     tmp => particle
     deallocate(tmp, stat = istat)
     if (istat /= 0) then
       call stop_program("[ERROR] Deallocate failed for deleting "// &
                         "particle number.", istat)
     end if

   end subroutine delete_particle_number

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine update_active_particles(cell_space, particles)
    type(cell_type), dimension(:), intent(in) :: cell_space
    type(particle_type), dimension(:), intent(in out) :: particles

    type(particle_number), pointer :: number
    integer :: ci

    particles(:)%is_active = .false.

    do ci = 1, size(cell_space)
      number => cell_space(ci)%first_number
      do !while(associated(number))
        if (.not. associated(number)) exit
        particles(number%i)%is_active = .true.
        number => number%next_number
      end do
    end do

  end subroutine update_active_particles

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine clean_up(particles)
    type(particle_type), dimension(:), allocatable, intent(in out) :: particles

    type(particle_number), pointer :: particle
    integer :: cz, cy, cx
    integer :: istat

    deallocate(particles, stat=istat)
    if (istat /= 0) then
       call stop_program("[ERROR] Failed deallocating particle_array.", istat)
    end if

    ! remove elements connected to the cell_space
    z_loop: do cz = 1, size(cell_space, dim=3)
       y_loop: do cy = 1, size(cell_space, dim=2)
          x_loop: do cx = 1, size(cell_space, dim=1)

             particle_loop: do
                particle => cell_space(cx, cy, cz)%first_number
                if (.not. associated(particle)) exit particle_loop

                cell_space(cx, cy, cz)%first_number => &
                            cell_space(cx, cy, cz)%first_number%next_number

                deallocate(particle, stat=istat)
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