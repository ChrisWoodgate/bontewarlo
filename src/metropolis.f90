!----------------------------------------------------------------------!
! metropolis.f90                                                       !
!                                                                      !
! Module containing routines relating to the Metropolis algorithm      !
! using Kawasaki dynamics.                                             !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module metropolis

  use initialise
  use kinds
  use shared_data
  use io
  use comms
  use c_functions
  use energetics
  use random_site
  use analytics
  use write_netcdf
  use write_xyz
  use write_diagnostics
  
  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Runs Metropolis with Kawasaki dynamics and performs simulated      !
  ! annealing.                                                         !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine metropolis_simulated_annealing(setup, my_rank)

    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Arrays for storing data
    type(run_params) :: setup
  
    ! Integers used in calculations
    integer :: i,j, div_steps, accept, n_save !,ierr
    
    ! Temperature and temperature steps
    real(real64) :: beta, temp, sim_temp, current_energy, step_E,     &
                    step_Esq, C, acceptance
  
    ! Name of file for grid state and xyz file at this temperature
    !character(len=34) :: grid_file
    character(len=36) :: xyz_file

    ! Name of file for writing diagnostics at the end
    character(len=43) :: diagnostics_file
  
    ! Name of file for writing radial densities at the end
    character(len=37) :: radial_file

    ! Radial densities at each temperature step
    real(real64), allocatable, dimension(:,:,:) :: r_densities
  
    ! Set up the lattice
    call initial_setup(setup, config)

    call lattice_shells(setup, shells, config)

    n_save=floor(real(setup%mc_steps)/real(setup%sample_steps))
    div_steps = setup%mc_steps/1000
  
    ! Are we swapping neighbours or on the whole lattice?
    if (setup%nbr_swap) then
      setup%mc_step => monte_carlo_step_nbr
    else
      setup%mc_step => monte_carlo_step_lattice
    end if

    ! Allocate memory for radial densities
    allocate(r_densities(setup%n_species, setup%n_species, &
                         setup%wc_range))

    if(my_rank == 0) then
      write(6,'(/,72("-"),/)')
      write(6,'(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
    end if

    ! Loop over temperature steps
    do j=1, setup%T_steps
  
      step_E = 0.0_real64; step_Esq=0.0_real64
      acceptance = 0.0_real64; r_densities = 0.0_real64
    
      ! Work out the temperature and corresponding beta
      temp = setup%T + real(j-1, real64)*setup%delta_T
      sim_temp = temp*k_b_in_Ry
      beta = 1.0_real64/sim_temp
    
      ! Store this in an array
      temperature(j) = temp
    
      !---------!
      ! Burn in !
      !---------!
      if (setup%burn_in) then
        do i=1, setup%burn_in_steps
          ! Make one MC move
          accept = setup%mc_step(config, beta)
        end do 
      end if

      !-----------------------!
      ! Main Monte Carlo loop !
      !-----------------------!
      do i=1, setup%mc_steps
    
        ! Make one MC move
        accept = setup%mc_step(config, beta)
  
        acceptance = acceptance + accept

        ! Store data for averaging if requested
        if (mod(i, setup%sample_steps) .eq. 0) then

          ! Current energy
          current_energy = setup%full_energy(config)

          ! Add this to total for averaging
          step_E   = step_E + current_energy

          ! Add square to total for averaging
          step_Esq = step_Esq + current_energy**2

          ! Add radial densities for averaging
          r_densities = r_densities                     &
                      + radial_densities(setup, config, &
                                setup%wc_range, shells)
        end if
    
      end do

      ! Store the average radial densities at this temperature
      rho_of_T(:,:,:,j) = r_densities/n_save
  
      ! Store the average energy per atom at this temperature
      energies_of_T(j) = step_E/n_save/setup%n_atoms
  
      ! Heat capacity (per atom) at this temperature
      C = (step_Esq/n_save - (step_E/n_save)**2)/(sim_temp*temp)/setup%n_atoms

      ! Acceptance rate at this temperature
      acceptance_of_T(j) = acceptance/real(setup%mc_steps)
  
      ! Store the specific heat capacity at this temperature
      C_of_T(j) = C
    
      ! Dump grids if needed
      if (setup%dump_grids) then
         if (my_rank .le. 1) then
 !       write(grid_file, '(A11 I3.3 A11 I4.4 F2.1 A3)') 'grids/proc_', my_rank, '_grid_at_T_', &
 !                                            int(temp), temp-int(temp), '.nc'
          write(xyz_file, '(A11 I3.3 A12 I4.4 F2.1 A4)') 'grids/proc_', &
          my_rank, 'config_at_T_', int(temp), temp-int(temp),'.xyz'
  
          ! Write xyz file
          call xyz_writer(xyz_file, config, setup)
  
 !       ! Write grid to file
 !       call ncdf_grid_state_writer(grid_file , ierr, &
 !                              config, temp, setup)
        end if
      end if
  
      ! Compute the radial densities at the end of this temperature
      ! call radial_densities(setup, config, setup%wc_range, shells, rho_of_T, j)
  
      if (my_rank ==0) then
        ! Write that we have completed a particular temperature
        write(6,'(a,f7.2,a)',advance='yes') &
        " Temperature ", temp, " complete on process 0."
        write(6,'(a,f7.2,a)',advance='yes') &
        " Internal energy was ", 13.606_real64*1000*energies_of_T(j), " meV/atom"
        write(6,'(a,f9.4,a)',advance='yes') &
        " Heat capacity was ", C_of_T(j)/k_B_in_Ry, " kB/atom"
        write(6,'(a,f7.2,a,/)',advance='yes') &
        " Swap acceptance rate was ", 100.0*acceptance_of_T(j), "%"
      end if
    
    end do ! Loop over temperature

    write(radial_file, '(A22 I3.3 A12)') 'radial_densities/proc_', my_rank, '_rho_of_T.nc'
    write(diagnostics_file, '(A17 I4.4 A22)') 'diagnostics/proc_', my_rank, &
                                         'energy_diagnostics.dat'
  
    
    ! Write energy diagnostics
    call diagnostics_writer(diagnostics_file, temperature, &
                            energies_of_T, C_of_T, acceptance_of_T)

    ! Write the radial densities to file
    call ncdf_radial_density_writer(radial_file, rho_of_T, &
                                  shells, temperature, energies_of_T, setup)

    ! Average results across the simulation
    call comms_reduce_results(setup)

    if (my_rank .eq. 0) then
      ! Write energy diagnostics
      call diagnostics_writer('diagnostics/av_energy_diagnostics.dat', temperature, &
                              av_energies_of_T, av_C_of_T, av_acceptance_of_T)
      !Write the radial densities to file
      call ncdf_radial_density_writer('radial_densities/av_radial_density.nc', av_rho_of_T, &
                                    shells, temperature, av_energies_of_T, setup)
    end if

    if(my_rank == 0) then
      write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if
  
  end subroutine metropolis_simulated_annealing

  !--------------------------------------------------------------------!
  ! Runs Metropolis with Kawasaki dynamics. Simulates annealing down   !
  ! to a target temperature then draws samples N steps apart, where N  !
  ! is chosen by the user.                                             !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine metropolis_decorrelated_samples(setup, my_rank)

    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Arrays for storing data
    type(run_params) :: setup
  
    ! Integers used in calculations
    integer :: i,j, div_steps, accept, n_save
    
    ! Temperature and temperature steps
    real(real64) :: beta, temp, sim_temp, current_energy, acceptance
  
    ! Name of xyz file
    character(len=42) :: xyz_file

    ! Set up the lattice
    call initial_setup(setup, config)

    call lattice_shells(setup, shells, config)

    n_save=floor(real(setup%mc_steps)/real(setup%sample_steps))
    div_steps = setup%mc_steps/1000
  
    ! Are we swapping neighbours or on the whole lattice?
    if (setup%nbr_swap) then
      setup%mc_step => monte_carlo_step_nbr
    else
      setup%mc_step => monte_carlo_step_lattice
    end if

    if(my_rank == 0) then
      write(6,'(/,72("-"),/)')
      write(6,'(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
    end if

    !---------------------------------------------------!
    ! Burn-in at each temperature (simulated annealing) !
    !---------------------------------------------------!
    do j=1, setup%T_steps
  
      ! Work out the temperature and corresponding beta
      temp = setup%T + real(j-1, real64)*setup%delta_T
      sim_temp = temp*k_b_in_Ry
      beta = 1.0_real64/sim_temp
    
      ! Burn in
      if (setup%burn_in) then

        acceptance = 0.0_real64

        do i=1, setup%burn_in_steps
          ! Make one MC move
          accept = setup%mc_step(config, beta)
          acceptance = acceptance + accept
        end do

        if(my_rank == 0) then
          write(6,'(a,f7.2,a)',advance='yes') &
          " Burn-in complete at temperature ", temp, " on process 0."
          write(6,'(a,i7,a)',advance='yes') &
          " Accepted ", int(acceptance), " Monte Carlo moves at this temperature,"
          write(6,'(a,f7.2,a,/)',advance='yes') &
          " Corresponding to an acceptance rate of ", &
          100.0*acceptance/float(setup%burn_in_steps), " %"
        end if

      end if

    end do ! Loop over temperature

    !--------------------!
    ! Target Temperature !
    !--------------------!
 
    acceptance=0

    do i=1, setup%mc_steps
    
        ! Make one MC move
        accept = setup%mc_step(config, beta)
  
        acceptance = acceptance + accept

        ! Draw samples
        if (mod(i, setup%sample_steps) .eq. 0) then

          ! Get the energy of this configuration
          current_energy = setup%full_energy(config)

          write(xyz_file, '(A11 I3.3 A8 I4.4 A6 I4.4 F2.1 A4)') &
          'grids/proc_', my_rank, '_config_',                   &
          int(i/setup%sample_steps), '_at_T_', int(temp),       &
          temp-int(temp),'.xyz'

          ! Write xyz file
          call xyz_writer(xyz_file, config, setup)
  
          if (my_rank == 0) then
            write(6,'(a,i7,a,/)',advance='yes') &
            " Accepted an additional ", int(acceptance), " Monte Carlo moves before sample."
          end if

          acceptance=0
          
        end if
    
      end do

    if(my_rank == 0) then
      write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if
  
  end subroutine metropolis_decorrelated_samples

  !--------------------------------------------------------------------!
  ! Runs one MC step assuming pairs swapped across entire lattice      !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  function monte_carlo_step_lattice(setup, config, beta) result(accept)
    !integer(int16), allocatable, dimension(:,:,:,:) :: config
    integer(int16), dimension(:,:,:,:) :: config
    class(run_params), intent(in) :: setup
    integer :: accept
    integer, dimension(4) :: rdm1, rdm2
    real(real64) , intent(in) :: beta
    real(real64) :: e_unswapped, e_swapped, delta_e
    integer(int16) :: site1, site2

    ! Generate random numbers
    rdm1 = setup%rdm_site()
    rdm2 = setup%rdm_site()

    ! Get what is on those sites
    site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
    site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

    ! Don't do anything if same species
    if (site1 == site2) then
      accept = 0
      return
    end if

    e_unswapped = pair_energy(setup, config, rdm1, rdm2)

    call pair_swap(config, rdm1, rdm2)

    e_swapped = pair_energy(setup, config, rdm1, rdm2)   

    delta_e = e_swapped - e_unswapped

    ! MC swap
    if(delta_e .lt. 0.0_real64) then
      accept = 1
      return
    else if (genrand() .lt. exp(-beta*delta_E)) then
      accept = 1
      return
    else
      accept = 0
      call pair_swap(config, rdm1, rdm2)
    end if

  end function monte_carlo_step_lattice

  !--------------------------------------------------------------------!
  ! Runs one MC step assuming only neighbours can be swapped           !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  function monte_carlo_step_nbr(setup, config, beta) result(accept)
    !integer(int16), allocatable, dimension(:,:,:,:) :: config
    integer(int16), dimension(:,:,:,:) :: config
    class(run_params), intent(in) :: setup
    integer :: accept
    integer, dimension(4) :: rdm1, rdm2
    real(real64) , intent(in) :: beta
    real(real64) :: e_unswapped, e_swapped, delta_e
    integer(int16) :: site1, site2

    ! Generate random numbers
    rdm1 = setup%rdm_site()
    rdm2 = setup%rdm_nbr(rdm1)

    ! Get what is on those sites
    site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
    site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

    ! Don't do anything if same species
    if (site1 == site2) then
      accept=0
      return
    end if

    e_unswapped = pair_energy(setup, config, rdm1, rdm2)

    call pair_swap(config, rdm1, rdm2)

    e_swapped = pair_energy(setup, config, rdm1, rdm2)   

    delta_e = e_swapped - e_unswapped

    ! MC swap
    if(delta_e .lt. 0.0_real64) then
      accept = 1
      return
    else if (genrand() .lt. exp(-beta*delta_E)) then
      accept = 1
      return
    else
      accept = 0
      call pair_swap(config, rdm1, rdm2)
    end if

  end function monte_carlo_step_nbr

end module metropolis
