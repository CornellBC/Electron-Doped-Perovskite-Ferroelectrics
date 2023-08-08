module Coefficients
    !  This module constains all of the cofficients for the model.

    use Physical_Constants

    implicit none

    ! Model Details
    character (len=20) :: version = "v2.2 (07 Apr 23)" ! Version number for the model.

    ! Convergence Criteria
    double precision, parameter, public :: mixing_pol = 1.0d0 ! Simple mixing strength for polarization
    double precision, parameter, public :: mixing_pot = 1.0d0 ! Simple mixing strength for potential
    double precision, parameter, public :: mixing_rho = 1.0d0 ! Simple mixing strenght for the charge density
    double precision, parameter, public :: mixing_background = 1.0d0 ! Simple mixing strength for the background polarization
    ! If the mixing strength = 1.0d0 then there is no mixing done with that parameter.
    double precision, parameter, public :: eps_conv = 1.0d-6 ! Value for determining convergence
    logical, public :: done = .false. ! When this value becomes true, it indicates that convergence has been achieved.
    character(len=3), public :: conv_type = "pot" ! Two options - potential = "pot", or polarization = "pol" 

    ! Looping Constraint
    integer, parameter, public :: counter_max = 10 ! This is the number of times that the main module will loop before halting.
    ! We can build in an import function to take in the resulting data as the next initialize in order to continue with another
    ! pass through the data.
    integer, parameter, public :: archive_threshold = 1 ! This is the number of loops the main module will run before it archives
    ! the polarization, charge density, and potential.  This data can be used as input for a future loop (via the Initialize module)
    ! or to monitor performance over time.

    ! Import/Export Data
    logical, public :: import = .true. ! This will import data from a data file as part of the Initialize module.
    logical, public :: export = .true. ! This will export data after some archive_threshold (set above) number of runs through 
    ! the full loop to the archive.dat file.

    ! Geometry Coefficients
    double precision, parameter, public :: latt_dist = 3.95d-10 ! 1.0d-9 ! 3.95d-10  ! Width of the unit cell/lattice in m. 
    double precision, parameter, public :: cg_dist = 1.0d-9 ! 1.0d-9 ! Course graining width (CG) in m 
    !Taken from https://materialsproject.org/materials/mp-5229/
    integer, parameter, public :: dim_x = 28 ! Bill's model is 28
    integer, parameter, public :: dim_z_FEM = 46 ! On the secondary lattice, there is one fewer dimension.  Bill's model is 46.
    integer, parameter, public :: dim_z_PC = 5 ! Width in lattice cells of the polar cap.  The total size of the material in the z-
    ! direction is dim_z_FEM + dim_z_PC.  Bill's model is 5.
    integer, parameter, public :: dim_z_total = dim_z_FEM + dim_z_PC ! Total width of the lattice cell in the z-direction.
    integer, parameter, public :: dim_y = 28 ! Needed to determine L_y
    
    ! LGD Equation
    double precision, parameter, public :: a1 = 1.0d8 ! 1.0d8 ! 2.50d5 ! Coefficient for the square term in the LGD equation.
    double precision, parameter, public :: a3 = -8.0d7   !-8.0d7 ! Additional coefficient to match Bill's LGD equation.
    double precision, parameter, public :: a11 = 1.70d9 ! 1.70d9 ! 1.63d10  ! Coefficient for the quartic term - b
    double precision, parameter, public :: a12 = 1.47d9 ! 1.47d9 ! 1.37d9  ! Coefficient for the squared cross term - b'.  
    double precision, parameter, public :: g11 = 1.0d-10 ! 1.0d-10 ! 1.0d-10   ! 1d-10 - Coefficient for the derivatives - g11.  
    double precision, parameter, public :: g44 = 1.0d-10 !  1.0d-10 ! 1.0d-10   ! 1d-10 - Coefficient for the cross derivatives - g44.  
    double precision, parameter, public :: eta1 = 0.0d0 ! TBD ! Coefficient for the charge * (P^2 + P^2) terms
    double precision, parameter, public :: eta2 = 0.0d0 ! TBD ! Coefficient for the charge * P^2 terms
    double precision, parameter, public :: eta3 = 0.0d0 ! TBD ! Coefficient for the sum of the first derivative terms
    double precision, parameter, public :: eta4 = 0.0d0 ! TBD ! Coefficient for the first dreivative derivative terms
    double precision, parameter, public :: eta5 = 0.0d0 ! TBD ! Coefficient for the first derivative squared terms
    double precision, parameter, public :: eta6 = 0.0d0 ! TBD ! Coefficient for the first derivative squared terms
    double precision, parameter, public :: eta7 = 0.0d0 ! TBD ! Coefficient for the first derivative squared terms
    double precision, parameter, public :: eta8 = 0.0d0 ! TBD ! Coefficient for the first derivative squared terms
    double precision, parameter, public :: eta9 = 0.0d0 ! TBD ! Coefficient for the first derivative squared terms
    double precision, dimension(28,51), public :: P_x, P_z ! Polarization on the secondary lattice in the x and z-directions 
    double precision, dimension(28,51), public :: P_x_old, P_z_old ! Old polarization values.  Used for convergence.
    double precision, dimension(28,51), public :: Pbx, Pbz ! Background polarizations
    double precision, dimension(28,51), public :: Pbx_old, Pbz_old ! Old background polarization values
    double precision, public :: P_min = 1.0d-8 ! This is the smallest value of polarization possible.  Used to terminated the LGD equation
    ! if one of the values is going to zero.
    double precision, parameter, public :: delta = 0.5d0! Interval governing change in P from cycle-to-cycle
    double precision, parameter, public :: eps = 1.0d-8 ! Convergence value to terminate the GNGA
    
    ! Potential
    double precision, dimension(28,51), public :: potential ! Potential on the primary lattice.
    double precision, dimension(28,51), public :: potential_old ! Old version of potential.  Used for convergence.
    double precision, dimension(28,51), public :: potential_free ! Potential from the charge density
    double precision, dimension(28,51), public :: potential_back ! Potential from the background polarization
    double precision, dimension(28,51), public :: potential_ferro ! Potential arising from ferroelectricity
    double precision, parameter, public :: eps_p = 1.0d-8 ! Convergence value to terminate the Relaxation Method algorithm
    double precision, parameter, public :: potential_bottom = 0.0d0 ! Fixed potential next to the bottom capacitor plate
    double precision, parameter, public :: potential_top = 0.0d0! Fixed potential next to the top capacitor plate
    double precision, parameter, public :: chi_PC = 25.0d0 ! 25 ! Dielectric susceptibility of the polar cap ! LAO = 25 fm Bill's 
    !paper
    double precision, parameter, public :: chi_FEM = 4.5d0 ! 4.5 ! Dielectric susceptibility of the ferroelectric material ! 
    !STO = 4.5 fm Bill's paper

    ! Schrodinger 
    !double precision, parameter :: t = (hbar**2)/(2*elec_mass*(cg_dist**2)) ! Bill's hopping parameter
    double precision, parameter, public :: hop_stay_xy = 0 ! 4.0d0*t*cg_dist**2/latt_dist**2
    !-542.0d-3 / J_to_eV ! Energy associated with staying at the same site.
    double precision, parameter, public :: hop_stay_xz = 0 !-542.0d-3 / J_to_eV ! Energy associated with staying at the same site.
    double precision, parameter, public :: hop_stay_yz = 0 !-542.0d-3 / J_to_eV ! Energy associated with staying at the same site.  
    double precision, parameter, public :: hop_parallel = -236d-3/J_to_eV ! -t*(cg_dist**2)/(latt_dist**2) 
    ! -236.0d-3/J_to_eV ! 236 meV ! Energy associated with hopping to the same orbital type in the same plane.  
    double precision, parameter, public :: hop_perp = -35d-3/J_to_eV ! -t*(cg_dist**2)/(latt_dist**2) 
    ! 35 meV ! Energy associated with hopping to the same orbital type outside of the plane.
    double precision, dimension (1232,1232), public :: SE_matrix_xy, SE_matrix_xz, SE_matrix_yz ! Hamiltonian matrices
    double precision, dimension(1232), public :: W_xy, W_xz, W_yz ! Eigenenergies from the Hamiltonians
    double precision, parameter, public :: temp = 10 ! Temperature in Kelvin ! Where do we get this value?
    double precision, parameter, public :: eps_mu = 1.0d-12 ! Convergence value for the Bisection Method determining the chem pot.
    double precision, public :: chem_pot ! Chemical potential

    ! Electron Density
    double precision, parameter, public :: n2D = 0.1d0 ! 2D electron density per 2D unit cell.
    double precision, dimension(28,51), public :: rho_xy, rho_xz, rho_yz ! 2D free charge density for each orbital - xy, xz, yz.
    ! These are only defined in the FEM.
    double precision, dimension(28,51), public :: rho_xy_old, rho_xz_old, rho_yz_old ! Old versions of rho.  Used for convergence.

    ! Electric Field
    double precision, dimension(28,51), public :: E_x, E_z ! Electric field in the x and z-directions on the secondary lattice.
    ! Note that the z-direction has one fewer dimension than the primary lattice.  Defined over the entire lattice.

contains

    subroutine Convert_Primary_to_Secondary(primary,secondary)
        implicit none
        
        ! Variables
        double precision, dimension(28,51) :: primary  ! These dimensions are associated w/ the FEM only.
        double precision, dimension(28,51) :: secondary  ! These dimensions are associated w/ the FEM only.
        integer :: n, m

        ! Take a box of values on the primary lattice to determine the average value at the secondary lattice
        do n = 1, dim_x
            do m = 1, dim_z_FEM-1
                if (n==dim_x) then  ! Boundary condition where we need to wrap around the boundary in the x-direction
                    secondary(n,m) = 0.25*(primary(n,m) + primary(1,m) + primary(n,m+1) + primary(1,m+1))
                else
                    secondary(n,m) = 0.25*(primary(n,m) + primary(n+1,m) + primary(n,m+1) + primary(n+1,m+1))
                endif
            end do
        end do

    end subroutine
    
    subroutine Convert_Secondary_to_Primary(primary, secondary)  ! This doesn't work for values alongs the z-boundary
        implicit none
        
        ! Variables
        double precision, dimension(28,51) :: primary  ! These dimensions are associated w/ the FEM only.
        double precision, dimension(28,51) :: secondary ! These dimensions are associated w/ the FEM only.
        integer :: n, m

        ! Take a box of the secondary values to determine the average value on the primary lattice.  We can't define the values at
        ! the z-boundary (z = 1 or dim_z).
        do n = 1, dim_x
            do m = 1, dim_z_FEM-1 ! We can't generate a value for m = dim_z_total since we can't create a box for it.
                if (m==1) then
                    primary(n,m) = 0 ! This value is a placeholder, and shouldn't be used in future calculations.
                elseif ((n==1) .and. (m/=1)) then
                    primary(n,m) = 0.25*(secondary(dim_x,m-1) + secondary(n,m-1) + secondary(dim_x,m) + secondary(n,m))
                else
                    primary(n,m) = 0.25*(secondary(n-1,m-1) + secondary(n,m-1) + secondary(n-1,m) + secondary(n,m))
                end if
            end do
        end do
    
    end subroutine

end module Coefficients
