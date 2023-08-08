module Convergence

    ! This module determines whether the model has converged or not.
    ! The convergence criteria are modular, and can be done for polarization, potential or charge density.

    use Coefficients
    use Anderson_Mix

    implicit none

    contains

    subroutine Test_Convergence()
        ! We can check convergence using potential or polarization.

        implicit none

        ! Variables
        double precision :: max_pot_diff, max_pol_diff, max_P ! Maximum value for potential and polarization
        double precision :: max_rho_xy_diff, max_rho_xz_diff, max_rho_yz_diff ! Maximum value for charge density 
        double precision, dimension(:,:), allocatable :: potential_new ! Interim value for potential.  Used for convergence. 
        double precision, dimension(:,:), allocatable :: P_x_new, P_z_new ! Interim value for polarization.
        double precision, dimension(:,:), allocatable :: rho_xy_new, rho_xz_new, rho_yz_new ! Interim value for potential.  
        ! Used for convergence.  
        integer :: n, m, nz
        integer :: counter=0

        ! Anderson Mixing
        double precision, dimension(:), allocatable :: X, Y, P_Anderson

        ! Allocate memory

        ! Potential convergence
        allocate(potential_new(dim_x,dim_z_total))

        ! Polarization convergence
        allocate(P_x_new(dim_x,dim_z_total))
        allocate(P_z_new(dim_x,dim_z_total))

        ! Charge density convergence
        allocate(rho_xy_new(dim_x,dim_z_total))
        allocate(rho_xz_new(dim_x,dim_z_total))
        allocate(rho_yz_new(dim_x,dim_z_total))

        ! Anderson Mixing
        allocate(P_Anderson(dim_x*dim_z_total))
        allocate(X(dim_x*dim_z_total))
        allocate(Y(dim_x*dim_z_total))

        ! Select different forms of convergence - potential or polarization

        if (conv_type=="pot") then
        
            ! Determine the difference between the new and old potentials.

            ! For the first cycle, we set potential_old = potential derived from the original polarizations.  This allows for
            ! subsequent mixes to be close to the potential arising from the initial polarizations.

            max_pot_diff = MAXVAL(ABS(potential - potential_old))
            
            ! Write results to screen to monitor the results while the program is running.
            if (MODULO(counter,20)==0) then
                print "(a10,i6,a11,e21.15)", "Counter = ", counter, ", % Diff = ", max_pot_diff/MAXVAL(ABS(potential))
            end if

            ! Test for convergence

            if (max_pot_diff <= eps_conv*MAXVAL(ABS(potential))) then 
                done = .true.
                print *, "Model has converged."
            elseif (MAXVAL(ABS(potential))<= 1.0d-8) then
                ! This keeps the model from chasing a potential that is going to zero.
                done = .true.
                print *, "Potential is tiny."
            else ! Mix the potential            

                ! Two mixing choices.  You need to comment out the mixing you don't want. 
                ! Simple Mixing
                !    potential_new = mixing_pot*potential + (1.0d0 - mixing_pot)*potential_old ! Mixed value.
                !    potential = potential_new ! Set new potential for testing.

                ! Anderson Mixing for P_x
                    X = 0.0d0
                    Y = 0.0d0
                    P_Anderson = 0.0d0
                                    
                    X = reshape(potential_old, (/ dim_x*dim_z_total /))
                    Y = reshape(potential, (/ dim_x*dim_z_total /))
                    call AM(X, Y, P_Anderson)
                    potential = reshape(P_Anderson, (/ dim_x, dim_z_total /))

                ! Calculate the new electric field (Bill's version)
                E_x = 0.0d0
                E_z = 0.0d0
                
                nz = dim_z_total
                E_z(:,1:nz-1) = -( potential(:,2:nz) - potential(:,1:nz-1) &
                    + cshift(potential(:,2:nz), shift=1, dim=1) &
                    - cshift(potential(:,1:nz-1), shift=1, dim=1) )/(2*cg_dist)
                E_x(:,1:nz-1) = -( cshift(potential(:,1:nz-1),shift=1,dim=1) - potential(:,1:nz-1) &
                    + cshift(potential(:,2:nz), shift=1, dim=1) &
                    - (potential(:,2:nz-1) ))/(2*cg_dist)
                E_z(:,nz) = 0d0
                E_x(:,nz) = 0d0              

            end if   

        elseif (conv_type=="pol") then

            ! Determine the difference between the new and old potentials.
            
            max_pol_diff = MAXVAL(SQRT((P_x-P_x_old)**2 + (P_z - P_z_old)**2))

            max_P = MAXVAL(SQRT(P_x**2 + P_z**2))

            ! Write results to screen to monitor the results while the program is running.
            if (MODULO(counter,20)==0) then
                print "(a10,i6,a11,e21.15)", "Counter = ", counter, ", % Diff = ", max_pol_diff/max_P
            end if
            
            ! Test for convergence

            if (max_pol_diff <= eps_conv*max_P) then
                done = .true.
                print *, "Model has converged."
            elseif (max_P <= P_min) then
                done = .true.
                print *, "P is tiny."
            else ! Mix the new terms                       
                ! Simple Mixing
                    P_x_new = mixing_pol*P_x + (1.0d0 - mixing_pol)*P_x_old 
                    P_x = P_x_new
            
                    P_z_new = mixing_pol*P_z + (1.0d0 - mixing_pol)*P_z_old     
                    P_z = P_z_new   
                
                ! Anderson Mixing for P_x
                    !X = 0
                    !Y = 0
                    !P_x_Anderson = 0
                                    
                    !X = reshape(P_x, (/ dim_x*dim_z_total /))
                    !Y = reshape(P_x_old, (/ dim_x*dim_z_total /))
                    !call AM(X, Y, P_Anderson)
                    !P_x = reshape(P_Anderson, (/ dim_x, dim_z_total /))

                ! Anderson Mixing for P_z
                    !X = 0
                    !Y = 0
                    !P_Anderson = 0

                    !X = reshape(P_z, (/ dim_x*dim_z_total /))
                    !Y = reshape(P_z_old, (/ dim_x*dim_z_total /))
                    !call AM(X, Y, P_Anderson)
                    !P_x = reshape(P_Anderson, (/ dim_x, dim_z_total /))

            end if   

        endif 

        counter = counter +1
        
        ! Deallocate memory
        deallocate(potential_new)
        
        deallocate(P_x_new)
        deallocate(P_z_new)
        
        deallocate(rho_xy_new)
        deallocate(rho_xz_new)
        deallocate(rho_yz_new)

        ! Anderson Mixing
        deallocate(P_Anderson)
        deallocate(X)
        deallocate(Y)

        !print *, "End test convergence."

    end subroutine Test_Convergence

end module Convergence
