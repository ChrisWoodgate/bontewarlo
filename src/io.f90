!----------------------------------------------------------------------!
! io.f90                                                               !
!                                                                      !
! Module containing input/output routines.                             !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2024 !
!----------------------------------------------------------------------!
module io

  use kinds
  use shared_data
  use command_line
  use display
  
  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Subroutine to write software version, date, and time               !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine write_info(point)

    character (len=1) point
    character (len=10) date,time

    ! Get date and time
    call date_and_time(date=date,time=time)

    write(6,'(/,72("="))')
    write(6,'(18x,"Bonte Warlo Version 0.2.1, 18.07.24")')
    write(6,'(72("-"))')

    if (point .eq. 's') then
      write(6,'(15x,"This run started at",1x,a," on",1x,a)')           &
               time(1:2)//":"//time(3:4)//":"//time(5:6),              &
               date(7:8)//"."//date(5:6)//"."//date(1:4)
    else if(point .eq. 'f') then
      write(6,'(15x,"This run finished at",1x,a," on",1x,a)')          &
               time(1:2)//":"//time(3:4)//":"//time(5:6),              &
               date(7:8)//"."//date(5:6)//"."//date(1:4)
    endif

    write(6,'(72("="),/)' )

  end subroutine write_info

  !--------------------------------------------------------------------!
  ! Subroutine to make directories for storing data                    !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine make_data_directories(my_rank)
    integer :: my_rank

    ! make a directory for the grid states, diagnostics, 
    ! and radial_densities for each thread
    if(my_rank == 0) call execute_command_line('mkdir -p grids')
    if(my_rank == 0) call execute_command_line('mkdir -p diagnostics')
    if(my_rank == 0) call execute_command_line('mkdir -p radial_densities')
  end subroutine make_data_directories

  !--------------------------------------------------------------------!
  ! Subroutine to parse control file                                   !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2024 !
  !--------------------------------------------------------------------!
  subroutine read_control_file(filename, parameters, my_rank)

    character(len=*), intent(in) :: filename
    logical, dimension(11) :: check
    type(run_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios, my_rank
    logical :: exists

    check = .false.

    ios=0; line=0

    ! Defaults if these are not specified.
    parameters%wc_range = 2
    parameters%delta_T = 1
    parameters%T_steps = 1
    parameters%lro = .false.
    parameters%nbr_swap = .false.
    parameters%sample_steps = 1000
    parameters%radial_sample_steps = 0
    parameters%burn_in = .false.
    parameters%dump_grids = .false.
    parameters%burn_in_steps = 1000

    inquire(file=trim(filename), exist=exists)

    if (.not. exists) then
      stop 'Could not find input file ' // trim(filename)
    end if

    if(my_rank == 0) then
      write(6,'(26("-"),x,"Parsing input file",x,26("-"),/)')
    end if

    open(15, file=filename, iostat=ios)

    do while (ios==0)

      read(15, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('mode')  
          read(buffer, *, iostat=ios) parameters%mode
          check(1) = .true.
        case ('n_1')  
          read(buffer, *, iostat=ios) parameters%n_1
          check(2) = .true.
        case ('n_2')  
          read(buffer, *, iostat=ios) parameters%n_2
          check(3) = .true.
        case ('n_3')  
          read(buffer, *, iostat=ios) parameters%n_3
          check(4) = .true.
        case ('burn_in')  
          read(buffer, *, iostat=ios) parameters%burn_in
        case ('burn_in_steps')  
          read(buffer, *, iostat=ios) parameters%burn_in_steps
        case ('n_mc_steps')  
          read(buffer, *, iostat=ios) parameters%mc_steps
          check(5) = .true.
        case ('sample_steps')  
          read(buffer, *, iostat=ios) parameters%sample_steps
        case ('radial_sample_steps')
          read(buffer, *, iostat=ios) parameters%radial_sample_steps
        case ('n_species')  
          read(buffer, *, iostat=ios) parameters%n_species
          check(6) = .true.
        case ('lattice')  
          read(buffer, *, iostat=ios) parameters%lattice
          check(7) = .true.
        case ('lattice_parameter')  
          read(buffer, *, iostat=ios) parameters%lattice_parameter
          check(8) = .true.
        case ('interaction_file')  
          read(buffer, *, iostat=ios) parameters%interaction_file
          check(9) = .true.
        case ('interaction_range')  
          read(buffer, *, iostat=ios) parameters%interaction_range
          check(10) = .true.
        case ('wc_range')  
          read(buffer, *, iostat=ios) parameters%wc_range
        case ('T')  
          read(buffer, *, iostat=ios) parameters%T
          check(11) = .true.
        case ('delta_T')  
          read(buffer, *, iostat=ios) parameters%delta_T
        case ('T_steps')  
          read(buffer, *, iostat=ios) parameters%T_steps
        case ('dump_grids')  
          read(buffer, *, iostat=ios) parameters%dump_grids
        case ('lro')  
          read(buffer, *, iostat=ios) parameters%lro
        case ('nbr_swap')  
          read(buffer, *, iostat=ios) parameters%nbr_swap
        case default
        end select
      end if
    end do


    allocate(parameters%species_names(parameters%n_species))
    allocate(parameters%species_concentrations(0:parameters%n_species))
    allocate(parameters%species_numbers(parameters%n_species))

    if (parameters%radial_sample_steps .lt. 1) then
       parameters%radial_sample_steps = parameters%sample_steps
    end if

    parameters%species_concentrations = 0.0_real64
    parameters%species_numbers = 0

    line=0

    close(15)

    open(15, file=filename, iostat=ios)


    do while (ios==0)

      read(15, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('species_names')  
          read(buffer, *, iostat=ios) parameters%species_names
          check(1) = .true.
        case ('species_concentrations')  
          read(buffer, *, iostat=ios) parameters%species_concentrations(1:)
        case ('species_numbers')  
          read(buffer, *, iostat=ios) parameters%species_numbers(:)
        case default
        end select
      end if
    end do

    close(15)

    if (.not. any(check)) then
      stop 'Missing parameter in input file'
    end if

  end subroutine read_control_file

  subroutine print_parse()

    print*, '###############################'
    print*, '#     Parsing input file      #'

  end subroutine print_parse

  subroutine echo_control_file(parameters)
    type(run_params) :: parameters
    integer :: i

    print*, ' Read mode = ', parameters%mode
    print*, ' Read n_1 = ', parameters%n_1
    print*, ' Read n_2 = ', parameters%n_2
    print*, ' Read n_3 = ', parameters%n_3
    print*, ' Read n_basis = ', parameters%n_basis
    print*, ' Read n_steps = ', parameters%mc_steps
    print*, ' Read n_species = ', parameters%n_species
    print*, ' Read lattice parameter = ', parameters%lattice_parameter
    print*, ' Read lattice = ', parameters%lattice
    print*, ' Read interaction_file = ', parameters%interaction_file
    print*, ' Read wc_range = ', parameters%wc_range
    print*, ' Read T = ', parameters%T
    print*, ' Read delta_T = ', parameters%delta_T
    print*, ' Read T_steps = ', parameters%T_steps
    print*, ' Read dump_grids = ', parameters%dump_grids
    print*, ' Read lro = ', parameters%lro
    print*, ' Read nbr_swap = ', parameters%nbr_swap

    ! Print specified concentrations/numbers of atoms
    if (abs(sum(parameters%species_concentrations)-1.0_real64) &
        .lt. 0.001) then
      do i=1, parameters%n_species
        print*, ' Read species ', i, ' = ', parameters%species_names(i), &
                ' at concentration ', parameters%species_concentrations(i)
      enddo
    else
      do i=1, parameters%n_species
        print*, ' Read species ', i, ' = ', parameters%species_names(i), &
                ' at ', parameters%species_numbers(i), ' atoms'
      enddo
    end if

  end subroutine echo_control_file

  !-----------------------------------------------------!
  ! Subroutine to read in exchange parameters from file !
  !-----------------------------------------------------!
  subroutine read_exchange(setup, my_rank)
    type(run_params) , intent(in) :: setup
    integer :: my_rank

    if(my_rank == 0) then
      write(6,'(15("-"),x,"Reading atom-atom interaction parameters",x,15("-"),/)')
    end if

    V_ex = 0.0_real64
    open(16, file=setup%interaction_file)
    read(16,*) V_ex   
    close(16)

    if (my_rank .eq. 0) then
      ! Print it out as a sanity check for user
      call pretty_print_exchange(setup)
    end if

    if(my_rank == 0) then
      write(6,'(72("-"),/)')
    end if

  end subroutine read_exchange

  !------------------------------------------------------!
  ! Subroutine to parse in exchange parameters from file !
  !------------------------------------------------------!
  subroutine parse_inputs(setup, my_rank)
    type(run_params) :: setup
    integer :: my_rank
    character(len=30) :: control

    ! Parse the command line arguments
    call parse_args()

    ! Parse the name of the input file
    if(my_rank == 0) then
      write(6,'(22("-"),x,"Parsing name of input file",x,22("-"),/)')
    end if
      if(.not. get_arg('control', control)) then
        if (my_rank == 0) then
          print*, 'Input file not specified with "control=<name>"'
          print*, 'Defaulting to searching for "input.txt"'
          print*, ' '
        end if
        control = 'input.txt'
      else
        if (my_rank == 0) then
          print*, 'Input file name is: ', control
          print*, ' '
        end if
      end if

    ! Read the input file
    call read_control_file(control, setup, my_rank)

    if(my_rank == 0) then
      call echo_control_file(setup)
      write(6,'(/,20("-"),x,"Parsed input file successfully",x,20("-"),/)')
    end if

  end subroutine parse_inputs

  
  !--------------------------------------------------------------------!
  ! Subroutine to read and parse nested sampling control file          !
  !                                                                    !
  ! L. B. Partay, Warwick                                         2024 !
  !--------------------------------------------------------------------!
  subroutine read_ns_file(filename, parameters)
    character(len=*), intent(in) :: filename
    logical, dimension(8) :: check
    type(ns_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios

    check = .false.

    ios=0; line=0

    open(25, file=filename, iostat=ios)

    if (ios .ne. 0) then
      stop 'Could not parse input file. Aborting...'
    end if

    write(*,'(a)', advance='no') new_line('a')
    print*, '###############################'
    print*, '#  Parsing NS input file      #'
    print*, '###############################'

    print*, '# NS input file name: ', filename

    do while (ios==0)

      read(25, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('n_walkers')
          read(buffer, *, iostat=ios) parameters%n_walkers
          print*, '# Read n_walkers = ', parameters%n_walkers
        case ('n_steps')
          read(buffer, *, iostat=ios) parameters%n_steps
          print*, '# Read n_steps = ', parameters%n_steps
        case ('n_iter')
          read(buffer, *, iostat=ios) parameters%n_iter
          print*, '# Read n_iter = ', parameters%n_iter
        case ('outfile_ener')
          read(buffer, *, iostat=ios) parameters%outfile_ener
          print*, '# Read outfile_ener = ', parameters%outfile_ener
        case ('outfile_traj')
          read(buffer, *, iostat=ios) parameters%outfile_traj
          print*, '# Read outfile_traj = ', parameters%outfile_traj
        case ('traj_freq')
          read(buffer, *, iostat=ios) parameters%traj_freq
          print*, '# Write configuration every n-th NS iteration = ', parameters%traj_freq
        case default
          print*, '# Skipping invalid label'
        end select
      end if
    end do

    print*, '# Finished parsing NS input file #'
    print*, '###############################', new_line('a')
    close(25)

  end subroutine read_ns_file

  !--------------------------------------------------------------------!
  ! Subroutine to read and parse nested tmmc control file              !
  !                                                                    !
  ! H. Naguszewski, Warwick                                       2024 !
  !--------------------------------------------------------------------!
  subroutine read_tmmc_file(filename, parameters, my_rank)
    integer :: my_rank
    character(len=*), intent(in) :: filename
    logical, dimension(8) :: check
    type(tmmc_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios

    check = .false.

    ios=0; line=0

    open(25, file=filename, iostat=ios)

    if (ios .ne. 0) then
      stop 'Could not parse tmmc input file. Aborting...'
    end if

    if (my_rank == 0) then
      write(*,'(a)', advance='no') new_line('a')
      print*, '###############################'
      print*, '#  Parsing tmmc input file    #'
      print*, '###############################'

      print*, '# tmmc input file name: ', filename
    end if

    do while (ios==0)

      read(25, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('burn_in')
          read(buffer, *, iostat=ios) parameters%burn_in
          if (my_rank == 0) then
            print*, '# Read burn_in = ', parameters%burn_in
          end if
          check(1) = .true.
        case ('burn_in_sweeps')
          read(buffer, *, iostat=ios) parameters%burn_in_sweeps
          if (my_rank == 0) then
            print*, '# Read burn_in_sweeps = ', parameters%burn_in_sweeps
          end if
          check(2) = .true.
        case ('mc_sweeps')
          read(buffer, *, iostat=ios) parameters%mc_sweeps
          if (my_rank == 0) then
            print*, '# Read mc_sweeps = ', parameters%mc_sweeps
          end if
          check(3) = .true.
        case ('bins')
          read(buffer, *, iostat=ios) parameters%bins
          if (my_rank == 0) then
            print*, '# Read bins = ', parameters%bins
          end if
          check(4) = .true.
        case ('weight_update')
          read(buffer, *, iostat=ios) parameters%weight_update
          if (my_rank == 0) then
            print*, '# Read weight_update = ', parameters%weight_update
          end if
          check(5) = .true.
        case ('energy_min')
          read(buffer, *, iostat=ios) parameters%energy_min
          if (my_rank == 0) then
            print*, '# Read energy_min = ', parameters%energy_min
          end if
          check(6) = .true.
        case ('energy_max')
          read(buffer, *, iostat=ios) parameters%energy_max
          if (my_rank == 0) then
            print*, '# Read energy_max = ', parameters%energy_max
          end if
          check(7) = .true.
        case ('bin_overlap')
          read(buffer, *, iostat=ios) parameters%bin_overlap
          if (my_rank == 0) then
            print*, '# Read bin_overlap = ', parameters%bin_overlap
          end if
          check(8) = .true.
        case default
          if (my_rank == 0) then
            print*, '# Skipping invalid label'
          end if
        end select
      end if
    end do

    if (my_rank == 0) then
      print*, '# Finished parsing tmmc input file #'
      print*, '####################################', new_line('a')
    end if
    close(25)

    if (.not. any(check)) then
      stop 'Missing parameter in tmmc input file'
    end if

  end subroutine read_tmmc_file

  !--------------------------------------------------------------------!
  ! Subroutine to read and parse nested wang landau control file       !
  !                                                                    !
  ! H. Naguszewski, Warwick                                       2024 !
  !--------------------------------------------------------------------!
  subroutine read_wl_file(filename, parameters, my_rank)
    integer :: my_rank
    character(len=*), intent(in) :: filename
    logical, dimension(11) :: check
    type(wl_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios

    check = .false.

    ios=0; line=0

    open(25, file=filename, iostat=ios)

    if (ios .ne. 0) then
      stop 'Could not parse wang landau input file. Aborting...'
    end if

    if (my_rank == 0) then
      write(*,'(a)', advance='no') new_line('a')
      print*, '##################################'
      print*, '# Parsing wang landau input file #'
      print*, '##################################'

      print*, '# wang landau input file name: ', filename
    end if

    do while (ios==0)

      read(25, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('burn_in')
          read(buffer, *, iostat=ios) parameters%burn_in
          if (my_rank == 0) then
            print*, '# Read burn_in = ', parameters%burn_in
          end if
          check(1) = .true.
        case ('burn_in_sweeps')
          read(buffer, *, iostat=ios) parameters%burn_in_sweeps
          if (my_rank == 0) then
            print*, '# Read burn_in_sweeps = ', parameters%burn_in_sweeps
          end if
          check(2) = .true.
        case ('mc_sweeps')
          read(buffer, *, iostat=ios) parameters%mc_sweeps
          if (my_rank == 0) then
            print*, '# Read mc_sweeps = ', parameters%mc_sweeps
          end if
          check(3) = .true.
        case ('bins')
          read(buffer, *, iostat=ios) parameters%bins
          if (my_rank == 0) then
            print*, '# Read bins = ', parameters%bins
          end if
          check(4) = .true.
        case ('num_windows')
          read(buffer, *, iostat=ios) parameters%num_windows
          if (my_rank == 0) then
            print*, '# Read num_windows = ', parameters%num_windows
          end if
          check(5) = .true.
        case ('bin_overlap')
          read(buffer, *, iostat=ios) parameters%bin_overlap
          if (my_rank == 0) then
            print*, '# Read bin_overlap = ', parameters%bin_overlap
          end if
          check(6) = .true.
        case ('tolerance')
          read(buffer, *, iostat=ios) parameters%tolerance
          if (my_rank == 0) then
            print*, '# Read tolerance = ', parameters%tolerance
          end if
          check(7) = .true.
        case ('flatness')
          read(buffer, *, iostat=ios) parameters%flatness
          if (my_rank == 0) then
            print*, '# Read flatness = ', parameters%flatness
          end if
          check(8) = .true.
        case ('wl_f')
          read(buffer, *, iostat=ios) parameters%wl_f
          if (my_rank == 0) then
            print*, '# Read wl_f = ', parameters%wl_f
          end if
          check(9) = .true.
        case ('energy_min')
          read(buffer, *, iostat=ios) parameters%energy_min
          if (my_rank == 0) then
            print*, '# Read energy_min = ', parameters%energy_min
          end if
          check(10) = .true.
        case ('energy_max')
          read(buffer, *, iostat=ios) parameters%energy_max
          if (my_rank == 0) then
            print*, '# Read energy_max = ', parameters%energy_max
          end if
          check(11) = .true.

        case default
          if (my_rank == 0) then
            print*, '# Skipping invalid label'
          end if
        end select
      end if
    end do

    if (my_rank == 0) then
      print*, '# Finished parsing wang landau input file #'
      print*, '###########################################', new_line('a')
    end if
    close(25)

    if (.not. any(check)) then
      stop 'Missing parameter in wang landau input file'
    end if

  end subroutine read_wl_file

  !--------------------------------------------------------------------!
  ! Subroutine to read and parse nested wang landau control file       !
  !                                                                    !
  ! H. Naguszewski, Warwick                                       2024 !
  !--------------------------------------------------------------------!
  subroutine read_es_file(filename, parameters, my_rank)
    integer :: my_rank
    character(len=*), intent(in) :: filename
    logical, dimension(2) :: check
    type(es_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios

    check = .false.

    ios=0; line=0

    open(25, file=filename, iostat=ios)

    if (ios .ne. 0) then
      stop 'Could not parse energy spectrum input file. Aborting...'
    end if

    if (my_rank == 0) then
      write(*,'(a)', advance='no') new_line('a')
      print*, '######################################'
      print*, '# Parsing energy spectrum input file #'
      print*, '######################################'

      print*, '# energy spectrum input file name: ', filename
    end if

    do while (ios==0)

      read(25, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('mc_sweeps')
          read(buffer, *, iostat=ios) parameters%mc_sweeps
          if (my_rank == 0) then
            print*, '# Read mc_sweeps = ', parameters%mc_sweeps
          end if
          check(1) = .true.
        case ('unique_energy_count')
          read(buffer, *, iostat=ios) parameters%unique_energy_count
          if (my_rank == 0) then
            print*, '# Read unique_energy_count = ', parameters%unique_energy_count
          end if
          check(2) = .true.
        case default
          if (my_rank == 0) then
            print*, '# Skipping invalid label'
          end if
        end select
      end if
    end do

    if (my_rank == 0) then
      print*, '# Finished parsing energy spectrum input file #'
      print*, '###############################################', new_line('a')
    end if
    close(25)

    if (.not. any(check)) then
      stop 'Missing parameter in energy spectrum input file'
    end if

  end subroutine read_es_file

end module io
