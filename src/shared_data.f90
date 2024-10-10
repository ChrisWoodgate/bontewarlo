!----------------------------------------------------------------------!
! shared_data.f90                                                      !
!                                                                      !
! Module containing important shared data.                             !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module shared_data

  use kinds
  
  implicit none

  save

  ! Internally the code works in Rydbergs (for historical reasons).
  ! These constants are hard-coded conversions used at various points.
  real(real64), parameter :: k_b_in_eV &
                             =8.167333262e-5_real64

  real(real64), parameter :: k_b_in_Ry &
                             =8.167333262e-5_real64/13.605693122_real64

  real(real64), parameter :: eV_to_Ry = 13.605693122_real64

  ! Random number seed
  integer(8) :: seed

  !--------------------------------------------------------------------!
  ! Type storing parameters defining simulation (used at runtime)      !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  type run_params
    ! Paramater defining the mode of operation
    integer :: mode
    ! Number of unit cells in each direction and
    ! number of atoms in the basis
    integer :: n_1, n_2, n_3, n_basis
    ! Number of chemical species
    integer :: n_species
    ! Number of atoms
    integer :: n_atoms
    ! Burn in if doing simulated annealing?
    logical :: burn_in
    ! Number of burn-in steps (at each temperature if annealing)
    integer :: burn_in_steps
    ! Number of monte carlo steps (at each temperature if annealing)
    integer :: mc_steps
    ! Number of monte carlo steps between drawing energies
    integer :: sample_steps
    ! Number of monte carlo steps between drawing radial densities
    ! (MUST be a multiple of sample_steps)
    integer :: radial_sample_steps
    ! Fixed or time-based random seed
    integer :: seedtime=1
    ! Lattice type - name, e.g. fcc, bcc, hcp, ...
    character(len=20) :: lattice
    ! Lattice parameter (for writing xyz file)
    real(real64) :: lattice_parameter
    ! Lattice vectors (for writing xyz file)
    real(real64), dimension(3,3) :: lattice_vectors
    ! Vector to second basis atom (for writing xyz file)
    real(real64), dimension(3) :: basis_vectors
    ! Names of the chemical species
    character(len=2), dimension(:), allocatable :: species_names
    ! Concentrations of the chemical species
    real(real64), dimension(:), allocatable :: species_concentrations
    ! Number of atoms of each chemical species
    integer(real64), dimension(:), allocatable :: species_numbers
    ! Atom-atom interchange interaction file name
    character(len=50) :: interaction_file
    ! Inverse temperature
    real(real64) :: beta
    ! Temperature of simulation (or start temperature if annealing)
    real(real64) :: T
    ! Temperature step size if annealing
    real(real64) :: delta_T
    ! Number of temperature steps (if annealing)
    integer :: T_steps
    ! Interaction range (number of coordination shells)
    integer :: interaction_range
    ! Number of coordination shells for Warren-Cowley parameters
    integer :: wc_range
    ! Hamiltonian to use
    procedure(hamiltonian), pointer, pass :: full_energy => null()
    ! neighbour energy function to use
    procedure(neighbour), pointer, pass :: nbr_energy => null()
    ! Random site function to use
    procedure(rand_site), pointer, pass :: rdm_site => null()
    ! Atom swap flag
    ! If True then only swap neighbouring pairs of atoms
    ! If False then permit swaps across the lattice
    logical :: nbr_swap
    ! Do we want to store long range order parameters
    logical :: lro
    ! Do we want to store grids
    logical :: dump_grids
    ! Random neighbour function to use
    procedure(rand_neighbour), pointer, pass :: rdm_nbr => null()
    ! Monte Carlo step to call. (Neighbour swap or whole lattice swap)
    procedure(monte_carlo), pointer, pass :: mc_step => null()
  end type run_params

  !--------------------------------------------------------------------!
  ! Type storing parameters defining NS simulation (used at runtime)   !
  !                                                                    !
  ! L. B. Partay,  Warwick                                        2024 !
  !--------------------------------------------------------------------!
  type ns_params
    ! Number of nested sampling walkers
    integer :: n_walkers
    ! Number of steps required for walk generating a new configuration
    integer :: n_steps
    ! Number of NS iterations before sampling is finished
    integer :: n_iter
    ! Frequency of writing configuration (every traj_freq-th NS iteration)
    integer :: traj_freq
    ! Output filenames
    character(len=100) :: outfile_ener, outfile_traj  
  end type ns_params

  !--------------------------------------------------------------------!
  ! Type storing parameters defining tmmc simulation (used at runtime) !
  !                                                                    !
  ! H. Naguszewski,  Warwick                                      2024 !
  !--------------------------------------------------------------------!
  type tmmc_params
    ! Burn in if doing simulated annealing?
    logical :: burn_in
    ! Number of burn-in sweeps (each sweep is n_atoms mc steps)
    integer :: burn_in_sweeps
    ! Number of sweeps (each sweep is n_atoms mc steps)
    integer :: mc_sweeps
    ! Number of bins across energy range
    integer :: bins
    ! Number of bias weight updates
    integer :: weight_update
    ! Energy range minimum
    real :: energy_min
    ! Energy range maximum
    real :: energy_max
    ! Number of bins in the overlap region
    integer :: bin_overlap
  end type tmmc_params

  !--------------------------------------------------------------------!
  ! Type storing parameters defining wang landau simulation            !
  !  (used at runtime)                                                 !
  !                                                                    !
  ! H. Naguszewski,  Warwick                                      2024 !
  !--------------------------------------------------------------------!
  type wl_params
    ! Burn in if doing simulated annealing?
    logical :: burn_in
    ! Number of burn-in sweeps (each sweep is n_atoms mc steps)
    integer :: burn_in_sweeps
    ! Number of  sweeps (each sweep is n_atoms mc steps)
    integer :: mc_sweeps
    ! Number of bins across energy range
    integer :: bins
    ! Number of energy windows
    integer :: num_windows
    ! Number of bins in the overlap region
    integer :: bin_overlap
    ! Tolerance for wang landau
    real :: tolerance
    ! Flatness for wang landau histogram
    real :: flatness
    ! Wang Landau density of states histogram tuning parameter
    real :: wl_f
    ! Energy range minimum
    real :: energy_min
    ! Energy range maximum
    real :: energy_max

  end type wl_params

  !--------------------------------------------------------------------!
  ! Type storing parameters defining energy spectrum determining       !
  !  (used at runtime)                                                 !
  !                                                                    !
  ! H. Naguszewski,  Warwick                                      2024 !
  !--------------------------------------------------------------------!
  type es_params
    ! Number of mc sweeps (each sweep is n_atoms mc steps)
    integer :: mc_sweeps
    ! Number of unique energies to be recorded
    integer :: unique_energy_count
  end type es_params

  !--------------------------------------------------------------------!
  ! Interface for neighbour implementation                             !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  interface

    ! Simulation energy
    function hamiltonian(setup, config)
      use kinds
      import :: run_params
      !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
      integer(int16), dimension(:,:,:,:), intent(in) :: config
      real(real64) :: hamiltonian
      class(run_params), intent(in) :: setup
    end function

    ! Neighbour
    function neighbour(setup, config, site_i, site_j, site_k)
      use kinds
      import :: run_params
      !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
      integer(int16), dimension(:,:,:,:), intent(in) :: config
      real(real64) :: neighbour
      class(run_params), intent(in) :: setup
      integer, intent(in) :: site_i, site_j, site_k
    end function

    ! Random site on the lattice
    function rand_site(setup)
      use kinds
      import :: run_params
      class(run_params), intent(in) :: setup
      integer, dimension(4) :: rand_site
    end function

    ! Random neighbour of that site
    function rand_neighbour(setup, site)
      use kinds
      import :: run_params
      class(run_params), intent(in) :: setup
      integer, dimension(4), intent(in) :: site
      integer, dimension(4) :: rand_neighbour
    end function

    ! Type of Monte Carlo step
    function monte_carlo(setup, config, beta) result(accept)
      use kinds
      import :: run_params
      !integer(int16), allocatable, dimension(:,:,:,:) :: config
      integer(int16), dimension(:,:,:,:) :: config
      class(run_params), intent(in) :: setup
      integer :: accept
      real(real64), intent(in) :: beta
    end function monte_carlo

  end interface

  ! Arrays that are only used by rank 1
  real(real64), allocatable, dimension(:) :: av_energies_of_T,        &
                                             av_C_of_T,               &
                                             av_acceptance_of_T
  real(real64), allocatable, dimension(:,:,:,:) :: av_rho_of_T

  ! Arrays used on all processors
  ! Indices run (basis, x, y, z)
  integer(int16), allocatable, dimension(:,:,:,:) :: config
  real(real64), allocatable, dimension(:) :: energies_of_T, C_of_T,   &
                                             acceptance_of_T,         &
                                             temperature
  real(real64), allocatable, dimension(:,:,:,:,:) :: order_of_T
  real(real64), allocatable, dimension(:,:,:,:) :: rho_of_T
  real(real64), dimension(:,:,:,:), allocatable :: order
  real(real64), dimension(:,:,:), allocatable :: V_ex
  real(real64), allocatable, dimension(:) :: shells
  
end module shared_data
