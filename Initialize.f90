module Initialize
    ! This module creates the initial values for the polarization, potential, and charge density.
    ! This module needs the "Coefficients Module" and the "Physical Constants Module."
    
    use Physical_Constants
    use Coefficients
    use WriteToFile

    implicit none

    integer, private :: n, m, nz
    
    contains

        subroutine initialize_model()
        
        implicit none

        !print *, "Start initialize."
     
        ! Condition the data
        ! To simplify writing to file, we have given all of the variables the same dimension.  We set all of the unused array 
        ! elements to zero.
        ! Polarization is only defined up to dim_z_FEM-1
        ! Charge density is only defined up to dim_z_FEM
        ! Potential is defined for all dim_z_total
        ! Variables are described in the Coefficients module.
        
        ! Set variables to zero
        P_x = 0.0d0
        P_x_old = 0.0d0
        P_z = 0.0d0
        P_z_old = 0.0d0
        Pbx = 0.0d0
        Pbz = 0.0d0
        Pbx_old = 0.0d0
        Pbz_old = 0.0d0
        potential = 0.0d0  
        potential_old = 0.0d0
        potential_ferro = 0.0d0
        potential_free = 0.0d0
        potential_back = 0.0d0
        rho_xy = 0.0d0
        rho_xy_old = 0.0d0
        rho_xz = 0.0d0
        rho_xz_old = 0.0d0
        rho_yz = 0.0d0
        rho_yz_old = 0.0d0
        E_x = 0.0d0
        E_z = 0.0d0

        ! Set the Initial Conditions

        if (import .eqv. .true.) then

            ! Import the initial conditions for polarization, potential and charge density from the import.dat file.
            ! Normally this will be a subset of the data contained in the archive.dat or data.dat files that will act as the new
            ! input for another run.

            call ImportFromFile()

            ! Ensure imported data satisfies our initial conditions.  If it is taken from a previous run, it should already satisfy
            ! these constraints.
            P_x(:,1) = 0.0d0
            P_x(:,dim_z_FEM:dim_z_total) = 0.0d0
            P_z(:,1) = 0.0d0
            P_z(:,dim_z_FEM:dim_z_total) = 0.0d0
            potential(:,1) = potential_bottom
            potential(:,dim_z_total) = potential_top
            potential_ferro(:,1) = potential_bottom
            potential_ferro(:,dim_z_total) = potential_top
            potential_free(:,1) = potential_bottom
            potential_free(:,dim_z_total) = potential_top
            potential_back(:,1) = potential_bottom
            potential_back(:,dim_z_total) = potential_top
            rho_xy(:,1) = 0.0d0
            rho_xz(:,1) = 0.0d0
            rho_yz(:,1) = 0.0d0
            rho_xy(:,dim_z_FEM:dim_z_total) = 0.0d0
            rho_xz(:,dim_z_FEM:dim_z_total) = 0.0d0
            rho_yz(:,dim_z_FEM:dim_z_total) = 0.0d0

        else           
            ! Initialize Polarization 

            ! Set polarization to zero in the arrays.  This initializes the polarization in the polar caps to be zero.
            P_x = 0.0d0
            P_z = 0.0d0

            ! Populate the polarization in the FEM
            do m = 2, dim_z_FEM-1 ! Secondary Lattice.  From the boundary conditions, P(m=1) = 0. 
                do n = 1, dim_x                
                    ! Kittel Domains
                        !if (n<6) then
                        !    P_z(n,m) = SQRT(-a3/(2.0d0*a11))
                        !elseif ((n==6) .and. (n==20)) then
                        !    P_z(n,m) = 0.0d0
                        !elseif ((n>6) .and. (n<20)) then
                        !    P_z(n,m) = -SQRT(-a3/(2.0d0*a11))
                        !else
                        !    P_z(n,m) = SQRT(-a3/(2.0d0*a11))
                        !end if
                    ! Sinusoid
                        P_z(n,m) = 0.1d0*SQRT(-a3/(2.0d0*a11))*SIN(6.28*(n-1)/dim_x+1.5d0)*SIN(3.14*(m-1)/dim_z_FEM)
                    ! Random
                        !call RANDOM_NUMBER(P_x(n,m))
                        !call RANDOM_NUMBER(P_z(n,m))
                    ! Uniform
                        !P_x(n,m) = 0.2d0
                        !P_z(n,m) = 0.2d0
                end do
            end do

            ! Initial Charge Density

            ! Only one of polarization or charge density needs to be initialized, although you can do both.  If you initialize 
            ! the charge density, leave the polarization as zero.
        
            !do m = 2, dim_z_FEM-1
            !    do n = 1, dim_x                
                    ! Plane of Charge
                    !    if (m==dim_z_FEM/2) then ! Plane of charge along the centre. 
                        !    rho_xy(n,m) = -1.0d0*n2D*funda_charge/((latt_dist**2)*cg_dist)/3.0d0
                        !    rho_xz(n,m) = -1.0d0*n2D*funda_charge/((latt_dist**2)*cg_dist)/3.0d0
                        !    rho_yz(n,m) = -1.0d0*n2D*funda_charge/((latt_dist**2)*cg_dist)/3.0d0
                        ! The division by 3 shares the electrons evenly between each t2g orbital. 
                        !else
                        !    rho_xy(n,m) = 0.0d0
                        !    rho_xz(n,m) = 0.0d0
                        !    rho_yz(n,m) = 0.0d0
                        !endif
                    ! Uniform Density                  
                        !rho_xy(n,m) = -1.0d0*funda_charge*n2D/((latt_dist**2)*cg_dist*(dim_z_FEM-2))/3.0d0
                        !rho_xz(n,m) = -1.0d0*funda_charge*n2D/((latt_dist**2)*cg_dist*(dim_z_FEM-2))/3.0d0
                        !rho_yz(n,m) = -1.0d0*funda_charge*n2D/((latt_dist**2)*cg_dist*(dim_z_FEM-2))/3.0d0
                !end do
            !end do
            
            ! Initialize Potential - This is generally left as zero, and produced via the Solve_Potential subroutine.  It could also
            ! be initialized using other means, but wasn't in this case.

            potential_back = 0.0d0
            potential_ferro = 0.0d0
            potential_free = 0.0d0
            potential = potential_back + potential_ferro + potential_free

            ! Calculate the Electric Field - This was stolen from Bill's code and slightly adjusted.
            nz = dim_z_total
            E_z(:,1:nz-1) = -( potential(:,2:nz) - potential(:,1:nz-1) &
                + cshift(potential(:,2:nz), shift=1, dim=1) &
                - cshift(potential(:,1:nz-1), shift=1, dim=1) )/(2*cg_dist)
            E_x(:,1:nz-1) = -( cshift(potential(:,1:nz-1),shift=1,dim=1) - potential(:,1:nz-1) &
               + cshift(potential(:,2:nz), shift=1, dim=1) &
                - (potential(:,2:nz-1) ))/(2*cg_dist)
            E_z(:,nz) = 0.0d0  ! There are no m = dim_z_total points on the secondary lattice.
            E_x(:,nz) = 0.0d0  ! There are no m = dim_z_total points on the secondary lattice.

        end if
        
        ! Write the Coefficients and Initial Conditions to the data.dat file

        call WriteInitialize()

        ! Set old values for polarization, potential, and charge density
        P_x_old = P_x
        P_z_old = P_z
        potential_old = potential  
        rho_xy_old = rho_xy
        rho_xz_old = rho_xz
        rho_yz_old = rho_yz
        Pbx_old = 0.0d0
        Pbz_old = 0.0d0

        end subroutine initialize_model

end module Initialize
