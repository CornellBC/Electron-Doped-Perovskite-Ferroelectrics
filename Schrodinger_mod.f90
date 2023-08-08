module Schrodinger
    ! Solves for the wavefunction using the matrix Hamiltonian for a grid (n,m) where 1<= n <= dim_x, and 1 <= m <= dim_z
    ! Using n2D, we calculate the chemical potential based on the energies from the Hamiltonians.
    ! Using the wavefunctions, energies and chemical potential, we calculate the electron density.
    
    use Physical_Constants
    use Coefficients
    use Initialize

    implicit none

    double precision, parameter :: beta = 1.0d0/(Boltzmann*temp) ! Beta coefficient for the Fermi function

    contains

        integer function pos(n,m,dim) result(W)
            ! Convert between (n,m) and a 1D position indicator
            implicit none
            integer, intent(in) :: n,m,dim
            W = n + (m-1)*dim
        end function

        double precision function fermi(energy1, energy2, mu) result(W)
            ! Calculate the Fermi function
            implicit none
            double precision, intent(in) :: energy1, energy2, mu
            W = 1/(1 + EXP(beta*(energy1 + energy2 - mu)))
        end function
        
        !double precision function fermi_int(energy) result(W)
            ! Calculate the Fermi integral.  Not used in this model.
        !    implicit none
        !    double precision, intent(in) :: energy
        !    double precision :: a, b, c, y0, y1, sqrtpi, sqrt2

        !    sqrtpi = SQRT(pi)
        !    sqrt2 = SQRT(2.0d0)

        !    a = SQRT(1.0d0 + 15.0d0/8.0d0 + 1.0d0/160.0d0)
        !    b = 1.8d0 - 0.61d0/2.0d0
        !    c = 2.0d0 + (2.0d0-sqrt2)*sqrt2
        !    y0 = 0.5d0*sqrt2 / SQRT(b + energy + (ABS(energy-b)**c + a**c)**(1.0d0/c))
        !    y1 = EXP(-energy)/sqrtpi
        !    W = 1.0d0/(y0+y1)/sqrtpi 

        !end function

        subroutine Solve_Schrodinger()
            
            implicit none
            
            ! Variables
            double precision :: L_x, L_y, n2D_temp, n2D_diff, n2D_real
            double precision :: chem_pot_max
            double precision :: chem_pot_min
            double precision :: sum1
            double precision, dimension (3) :: eigen
            integer :: n, m, en, k ! Counters for position (n,m), energy eigenvalues in x-z (en), and energy in y-direction (k)

            ! Variables for pre-calculation
            double precision :: Vol_cg
            double precision, dimension(:), allocatable :: cosk ! Values for COS given k = 1, dim_y

            ! Variables for DSYEVD
            character(1) :: JOBZ, UPLO
            integer :: LDA, dim_SE, LWORK_xy, LWORK_xz, LWORK_yz, INFO, LIWORK_xy, LIWORK_xz, LIWORK_yz
            integer, dimension(:), allocatable :: IWORK
            double precision, dimension(:), allocatable :: WORK
            
            ! Fermi Integral Variables
            double precision :: chem_pot_int, prefactor
            double precision, dimension(:), allocatable :: FI_xy, FI_xz, FI_yz ! Calculate the FI contribution for each energy

            ! Initial Values
            JOBZ = 'V' ! Eigenvalues and eigenvectors are created. 
            UPLO = 'U' ! Uses the upper triangle portion of the SE_matrix.
            dim_SE = dim_x*(dim_z_FEM-2) ! Dimension of the SE_matrix.  We're solving this over the primary lattice.
            LDA = dim_x*(dim_z_FEM-2) ! Number of rows of the SE_matrix.
            L_x = dim_x*cg_dist ! We coarse-grain the x-direction.
            L_y = dim_y*cg_dist ! We coarse-grain the y-direction.
            
            !print *, "Start Solve_Schrodinger."

            ! Allocate memory for the arrays
            allocate(FI_xy(dim_SE))
            allocate(FI_xz(dim_SE))
            allocate(FI_yz(dim_SE))
            allocate(cosk(dim_y))

            ! Calculate the prefactor for the COS correction
            prefactor = latt_dist**2/cg_dist**2

            ! Pre-calculate COS(k) for k = 1, dim_y
            do k = 1, dim_y
                cosk(k) = prefactor*COS(2*pi*k/dim_y)
            end do

            ! Set the old values for the electron density before they're updated
            rho_xy_old = rho_xy
            rho_xz_old = rho_xz
            rho_yz_old = rho_yz 
            
            ! Since we're solving for eigenvalues and eigenvectors, we don't need to define the wavefunction Z, or the energy.
            ! We will use the output from DSYEV to determine the different free charge densities.  
        
            ! Initialize the SE matrices with zeros.
            
            SE_matrix_xy = 0.0d0
            SE_matrix_xz = 0.0d0
            SE_matrix_yz = 0.0d0

            ! Create the SE matrices for H_xy, H_xz, H_yz
            do m = 2, dim_z_FEM-1
                do n = 1, dim_x                
                    !Case 1 - Bottom Left Hand Corner - The new boundary is at m=2.
                    if ((n==1) .and. (m==2)) then
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m)  + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xy(pos(n,m,dim_x),pos(n,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2) - Doesn't exist.
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m)  + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xz(pos(n,m,dim_x),pos(n,m-1,dim_x)) = hop_parallel - Doesn't exist for this case.
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m)  + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_yz(pos(n,m,dim_x),pos(n,m-1,dim_x)) = hop_parallel - Doesn't exist for this case.
                    !Case 2 - Left-hand Boundary
                    elseif ((n==1) .and. (m/=2) .and. (m/=dim_z_FEM-1)) then
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                    !Case 3 - Top Left-hand Corner
                    elseif ((n==1) .and. (m==dim_z_FEM-1)) then
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(dim_x,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                    !Case 4 - Top Boundary
                    elseif ((n/=1) .and. (n/=dim_x) .and. (m==dim_z_FEM-1)) then
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                    !Case 5 - Top Right-hand Corner
                    elseif ((n==dim_x) .and. (m==dim_z_FEM-1)) then
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(dim_x-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(dim_x-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(dim_x-1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                    !Case 6 - Right-hand Boundary
                    elseif ((n==dim_x) .and. (m/=2) .and. (m/=dim_z_FEM-1)) then
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(dim_x-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(dim_x-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(dim_x-1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                    !Case 7 - Bottom Right-hand Corner
                    elseif ((n==dim_x) .and. (m==2)) then
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(dim_x-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(dim_x-2,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(dim_x-1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                    !Case 8 - Bottom Boundary
                    elseif ((n/=1) .and. (n/=dim_x) .and. (m==2)) then
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            !SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                    !Case 9 - Bulk Value    
                    else
                        ! H_xy
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xy & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xy(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                        ! H_xz
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_xz & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_xz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                        ! H_yz
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-1,dim_x)) = hop_stay_yz & 
                            + 2.0d0*hop_perp*(1 - (latt_dist**2)/(cg_dist**2)) & 
                            + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2)) &
                            - funda_charge*potential(n,m) + 2.0d0*hop_parallel*(1 - (latt_dist**2)/(cg_dist**2))
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n+1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n-1,m-1,dim_x)) = hop_perp*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                            SE_matrix_yz(pos(n,m-1,dim_x),pos(n,m-2,dim_x)) = hop_parallel*(latt_dist**2)/(cg_dist**2)
                    end if
                end do
            end do     
            
            ! Calculate the eigenvalues and eigenvectors of the SE matrix.
            ! Use the DSYEVD package from LAPACK.

            ! Solve for the eigenvalue/vectors of the xy-orbital
            
            ! Allocate matrices
            allocate(WORK(1))
            allocate(IWORK(1))

            ! Optimize the size of the WORK and IWORK matrices for each orbital
            INFO = 0
            IWORK = 0
            
            call DSYEVD('V','U', dim_SE, SE_matrix_xy, dim_SE, W_xy, WORK, -1, IWORK, -1, INFO)
            LWORK_xy = INT(WORK(1))
            LIWORK_xy = MAX(IWORK(1),1)

            INFO = 0
            IWORK = 0
            
            call DSYEVD(JOBZ,UPLO, dim_SE, SE_matrix_xz,LDA,W_xz,WORK,-1,IWORK,-1,INFO)
            LWORK_xz = INT(WORK(1))
            LIWORK_xz = MAX(IWORK(1),1)

            INFO = 0
            IWORK = 0
            
            call DSYEVD(JOBZ,UPLO, dim_SE, SE_matrix_yz,LDA,W_xz,WORK,-1,IWORK,-1,INFO)
            LWORK_yz = INT(WORK(1))
            LIWORK_yz = MAX(IWORK(1),1)

            deallocate(WORK)
            deallocate(IWORK)
            
            ! Solve for xy-orbital
            
            ! Initialize the energies
            W_xy = 0.0d0  

            ! Allocate WORK matrix for xy
            allocate(WORK(LWORK_xy))
            allocate(IWORK(LIWORK_xy))

            INFO = 0
            ! Solve for the eigenvalues (output in W_xy) and eigenvectors (ouptut in SE_matrix_xy)
            call DSYEVD('V', 'U', dim_SE, SE_matrix_xy, dim_SE, W_xy, WORK, LWORK_xy, IWORK, LIWORK_xy, INFO)
                !if (INFO == 0) then
                !    print "(a38)", "Eigenvalues/vectors successful for xy."
                if (INFO > 0) then
                    print "(a43, i3)", "Unsuccessful calculations for xy.  Error = ", INFO
                !else
                !    print "(a25)", "Failed to converge in xy."
                end if

            ! Deallocate the WORK matrix for xy

            deallocate(WORK)
            deallocate(IWORK) 

            ! Solve for the xz-orbital
            
            ! Initialize the energies
            W_xz = 0.0d0

            ! Allocate WORK matrix for xz
            allocate(WORK(LWORK_xz))
            allocate(IWORK(LIWORK_xz))

            ! Solve for the eigenvalues (output in W_xz) and eigenvectors (ouptut in SE_matrix_xz)
           
            INFO = 0
            call DSYEVD(JOBZ,UPLO, dim_SE, SE_matrix_xz,LDA,W_xz,WORK,LWORK_xz,IWORK,LIWORK_xz,INFO)
                !if (INFO == 0) then
                !    print "(a38)", "Eigenvalues/vectors successful for xz."
                if (INFO > 0) then
                    print "(a43, i3)", "Unsuccessful calculations for xz.  Error = ", INFO
                !else
                !    print "(a25)", "Failed to converge in xz."
                end if

            ! Deallocate the work matrix
            deallocate(WORK)
            deallocate(IWORK)

            ! Solve for the yz-orbital
            
            ! Initialize the energies
            W_yz = 0.0d0

            ! Allocate WORK matrix for yz
            allocate(WORK(LWORK_yz))
            allocate(IWORK(LIWORK_yz))

            ! Solve for the eigenvalues (output in W_yz) and eigenvectors (ouptut in SE_matrix_yz)
            INFO = 0
            call DSYEVD(JOBZ,UPLO, dim_SE, SE_matrix_yz,LDA,W_yz,WORK,LWORK_yz,IWORK,LIWORK_yz,INFO)
                !if (INFO == 0) then
                !    print "(a38)", "Eigenvalues/vectors successful for yz."
                if (INFO > 0) then
                    print "(a43, i3)", "Unsuccessful calculations for yz.  Error = ", INFO
                !else
                !    print "(a25)", "Failed to converge in yz."
                end if

            ! Deallocate the work matrix
            deallocate(WORK)
            deallocate(IWORK)
            
            ! Solve for the chemical potential, chem_pot using the Bisection Method
            ! The solution using the "high" value must lead to a sum greater than n2D.  The solution for the "low" value must lead
            ! to a sum less than n2D.  We want to converge on n2D for some value of mu = chem_pot.
            
            ! Find values for the minimum and maximum chemical potential
            ! If the chemical potential is strongly negative, it makes the Fermi function zero.
            ! If the chemical potential is strongly positive, it makes the Fermi function one.
            ! The strength of the chemical potential is determined by its value relative to the energy.

            ! Find the largest value for the energies

            eigen(1) = MAXVAL(W_xy) 
            eigen(2) = MAXVAL(W_xz) 
            eigen(3) = MAXVAL(W_yz) 
            
            ! Initial guess for a maximal value for the chemical potential
            
            chem_pot_max = MAXVAL(eigen) 

            ! Find the smallest values for the energies

            eigen(1) = MINVAL(W_xy) 
            eigen(2) = MINVAL(W_xz) 
            eigen(3) = MINVAL(W_yz) 

            ! Initial guess for a minimal value for the chemical potential
            chem_pot_min = MINVAL(eigen) - 10.0d0/beta

            ! Check that these values correctly bracket the chemical potential.  We must have that the calculated solution for n2D
            ! using the minimum and maximum values are above and below the model value.
            
            ! Confirm that the min value for chemical potential leads to n2D_temp < n2D

            ! Generate the real n2D based on the lattice size.  We have n2D electrons per 2D unit cell (x and y-directions).
            ! Even though we're coarse grained in the x-direction, we assume that n2D is based on the physical unit cell.
            n2D_real = n2D/(latt_dist**2) 
            ! This is the 2D electron density (electrons per metres squared)                        

            ! Pre-calculate the volume of the coarse-grained unit cell.  Our model assumes that it's also coarse-grained 
            ! in the y-direction.
            Vol_cg = cg_dist**3 

            ! Confirm that the guesses for min and max chemical potential bracket the real value.

            ! Ensure looping works
            n2D_temp = ABS(n2D_real*10)

            n = 0  ! Prevent the loop from going forever
            do while ((n2D_temp > n2D_real) .and. (n<100000))
            
                ! Reset n2D_temp
                n2D_temp = 0.0d0

                ! Calculate the value for n2D_temp based on chem_pot_min
                do en = 1, dim_SE 
                    do k = 1, dim_y
                        n2D_temp = n2D_temp + 2.0d0/(L_x*L_y)*(fermi(W_xy(en), 2.0d0*hop_parallel*cosk(k),&
                         chem_pot_min) & 
                        + fermi(W_xz(en), 2.0d0*hop_perp*cosk(k), chem_pot_min) & 
                        + fermi(W_yz(en), 2.0d0*hop_parallel*cosk(k), chem_pot_min))                         
                    end do
                end do

                !if (n2D_temp <= n2D_real) then  
                !    print *, "Minimum chemical potential is good."
                if (n2D_temp >= n2D_real) then
                !    print *, "Minimum chemical potential is bad.  Set new value."
                    chem_pot_min = chem_pot_min -10.0d0/beta
                end if
                if (n==99999) then
                    print *, "Failed to find maximum chemical potential after 100000 cycles."
                endif 
                n = n+1

            end do

            ! Ensure looping works
            n2D_temp = -ABS(n2D_real*10)

            n = 0 ! Prevent the loop from going forever
            do while ((n2D_temp < n2D_real) .and. (n<100000))
            
                ! Reset n2D_temp
                n2D_temp = 0.0d0

                ! Calculate the value for n2D_temp based on chem_pot_max                
                do en = 1, dim_SE 
                    do k = 1, dim_y
                        n2D_temp = n2D_temp + 2.0d0/(L_x*L_y)*(fermi(W_xy(en), 2.0d0*hop_parallel*cosk(k),&
                         chem_pot_max) & 
                        + fermi(W_xz(en), 2.0d0*hop_perp*cosk(k), chem_pot_max) & 
                        + fermi(W_yz(en), 2.0d0*hop_parallel*cosk(k), chem_pot_max))                       
                    end do
                end do

                !if (n2D_temp >= n2D_real) then  
                !    print *, "Maximum chemical potential is good."
                if (n2D_temp <= n2D_real) then
                !    print *, "Maximum chemical potential is bad.  Set new value."
                    chem_pot_max = 2.0d0*ABS(chem_pot_max)
                endif
                if (n==99999) then
                    print *, "Failed to find maximum chemical potential after 100000 cycles."
                endif 
                n = n + 1
            end do
            
            ! Ensure looping works
            n2D_diff = ABS(eps_mu*n2D_real*10)
            chem_pot_int = 0.0d0

            n = 0 ! Prevent the loop from going forever
            do while (ABS(n2D_diff)>ABS(eps_mu*n2D_real) .and. (n<100000))
        
                ! Reset n2D_temporary value 
                n2D_temp = 0.0d0

                ! Mid-point value for chem_pot for evaluation
                chem_pot_int = (chem_pot_max + chem_pot_min)/2.0d0

                ! Calculate new n2D_temporary value                
                do en = 1, dim_SE 
                    do k = 1, dim_y
                        n2D_temp = n2D_temp + 2.0d0/(L_x*L_y)*(fermi(W_xy(en), 2.0d0*hop_parallel*cosk(k),&
                         chem_pot_int) & 
                        + fermi(W_xz(en), 2.0d0*hop_perp*cosk(k), chem_pot_int) & 
                        + fermi(W_yz(en), 2.0d0*hop_parallel*cosk(k), chem_pot_int))                     
                    end do
                end do

                n2D_diff = n2D_temp - n2D_real  

                if (n2D_diff<0.0d0) then
                    chem_pot_min = chem_pot_int
                else
                    chem_pot_max = chem_pot_int    
                end if

                if (n==99999) then
                    print *, "Failed to find chemical potential after 100000 cycles."
                endif 

                n = n+1
            end do

            ! Generate the electron density for the Fermi Integral          
            
            ! Reset charge densities
            rho_xy = 0.0d0
            rho_xz = 0.0d0
            rho_yz = 0.0d0

            ! Initialize Fermi functions (FI) for each orbital
            FI_xy = 0.0d0
            FI_xz = 0.0d0
            FI_yz = 0.0d0

            ! Pre-calculate the contribution from the Fermi Function
            do en = 1, dim_SE
                do k = 1, dim_y
                    FI_xy(en) = FI_xy(en) + 2.0*fermi(W_xy(en), 2.0d0*hop_parallel*cosk(k), &
                     chem_pot_int)
                    FI_xz(en) = FI_xz(en) + 2.0*fermi(W_xz(en), 2.0d0*hop_perp*cosk(k), &
                    chem_pot_int)
                    FI_yz(en) = FI_yz(en) + 2.0*fermi(W_yz(en), 2.0d0*hop_parallel*cosk(k), &
                    chem_pot_int)
                end do
            end do
            
            ! Calculate the charge density.  We divide by dim_y because this charge is spread evenly along the y-axis due to 
            ! translational invariance.
            do n = 1, dim_x
                do m = 2,dim_z_FEM-1
                    do en = 1, dim_SE
                            rho_xy(n,m) = rho_xy(n,m) -funda_charge/(Vol_cg*dim_y)* &
                            FI_xy(en)*(SE_Matrix_xy(pos(n,m-1,dim_x),en)**2)
                            rho_xz(n,m) = rho_xz(n,m) -funda_charge/(Vol_cg*dim_y)* &
                            FI_xz(en)*(SE_Matrix_xz(pos(n,m-1,dim_x),en)**2)
                            rho_yz(n,m) = rho_yz(n,m) -funda_charge/(Vol_cg*dim_y)* & 
                            FI_yz(en)*(SE_Matrix_yz(pos(n,m-1,dim_x),en)**2)               
                    end do
                end do
            end do

            !print *, "Total Electrons (from n2D) = ", n2D_real*L_x*cg_dist
            !print *, "Total Electrons (from rho) = ", sum1*Vol_cg/(-1.0d0*funda_charge)
            !pause
       
            ! Mix the charge densities.  We set mixing_rho = 1.0d0 to avoid this step.
            rho_xy = mixing_rho*rho_xy + (1.0d0 - mixing_rho)*rho_xy_old
            rho_xz = mixing_rho*rho_xz + (1.0d0 - mixing_rho)*rho_xz_old
            rho_yz = mixing_rho*rho_yz + (1.0d0 - mixing_rho)*rho_yz_old

            ! Set chem_pot
            chem_pot = chem_pot_int

            ! Deallocate memory
            deallocate(FI_xy)
            deallocate(FI_xz)
            deallocate(FI_yz)
            deallocate(cosk)
 
            !print *, "End Solve_Schrodinger."
        
        end subroutine Solve_Schrodinger
 
end module Schrodinger
