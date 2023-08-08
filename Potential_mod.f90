module Potential_Mod
    ! Using the charge density (rho) and polarization (P), we calculate the potential.
    ! Each location is determined by a grid of values - (n,m) - where 1<= n <= dim_x, and 1<= m <= dim_z.
    ! dim_x and dim_z are the number of lattice units in the x and z-direction respectively.
    ! Physical_Constants and Coefficients have the parameters required for this module.
   
    use Physical_Constants
    use Coefficients

    implicit none

    ! Variables
    double precision, dimension (:,:), allocatable, private :: rho_total, chi_mat  ! Total charge density, 
    ! susceptibility of material
    double precision, dimension (:,:), allocatable, private :: dPxx, dPzz, dQxx, dQzz ! Derivatives of polarization
    double precision, dimension (:,:), allocatable, private :: Px, Pz ! Net polarization, i.e. Px = P_x + Pbx
    integer, private :: n, m ! Looping parameters
 
contains
    integer function pos(n,m,dim) result(W)
        ! Convert between (n,m) and a 1D coordinate
        implicit none
        integer, intent(in) :: n,m,dim
        W = n + (m-1)*dim
    end function pos

    double precision function temp_pot(n, m, potential, rho, switch) result(W) 
        ! Calculate a new potential at (n,m) based on surrounding values
        implicit none
        integer, intent(in) :: n,m, switch
        double precision, intent(in) :: rho
        double precision, dimension(:,:), intent(in) :: potential
        
        if (switch==1) then ! This is the case at the left-hand edge (n = 1)
            W = 0.25d0*( potential(n+1,m) + potential(dim_x,m) + potential(n,m+1) + potential(n,m-1) + &
            (cg_dist**2)/perm_of_FS*rho )
        elseif (switch==2) then ! This is the case at the right-hand edge (n = dim_x)
            W = 0.25d0*( potential(1,m) + potential(n-1,m) + potential(n,m+1) + potential(n,m-1) + &
            (cg_dist**2)/perm_of_FS*rho )
        else ! This is the case when we're not at either edge in the x-direction (n <> 1 or dim_x)
            W = 0.25d0*( potential(n+1,m) + potential(n-1,m) + potential(n,m+1) + potential(n,m-1) + &
            (cg_dist**2)/perm_of_FS*rho )
        endif
    
    end function temp_pot

    subroutine background_charge(Px, Pz, dPxx, dPzz)
        ! Calculate the bound charge from the partial derivatives dPx/dx (dPxx) and dPz/dz (dPzz)
        implicit none

        ! Variables
        double precision, dimension (:,:), intent(in) :: Px, Pz
        double precision, dimension (:,:), intent(out) :: dPxx, dPzz 

        do m = 1, dim_z_total
            do n = 1, dim_x            
                ! Bottom left(n,m) = (1,1)
                if ((n==1) .and. (m==1)) then !  We assume that P and E below the FEM are both zero, so the m-1 terms disappear.
                    dPxx(n,m) = 0.5d0*(Px(1,m) - Px(dim_x,m))/cg_dist
                    dPzz(n,m) = 0.5d0*(Pz(1,m) + Pz(dim_x,m))/cg_dist
                ! Left-hand edge (n,m) = (1,m) & dim_z_PC/=0
                elseif ((n==1) .and. (m/=1) .and. (m<dim_z_total)) then
                    dPxx(n,m) = 0.5d0*(Px(1,m) - Px(dim_x,m) + Px(1,m-1) - Px(dim_x,m-1))/cg_dist
                    dPzz(n,m) = 0.5d0*(Pz(dim_x,m) - Pz(dim_x,m-1) + Pz(1,m) - Pz(1,m-1))/cg_dist               
                ! Top left-hand boundary between the PC and capacitor (n=1) and (m=dim_z_total)
                elseif ((n==1) .and. (m==dim_z_total)) then
                    dPxx(n,m) = 0.5d0*(Px(1,m-1) - Px(dim_x,m-1))/cg_dist
                    dPzz(n,m) = -0.5d0*(Pz(1,m-1) + Pz(dim_x,m-1))/cg_dist
                ! Top boundary between between the PC and capacitor
                elseif ((n/=1) .and. (n/=dim_x) .and. (m==dim_z_total)) then
                    dPxx(n,m) = 0.5d0*(Px(n,m-1) - Px(n-1,m-1))/cg_dist
                    dPzz(n,m) = -0.5d0*(Pz(n,m-1) + Pz(n-1,m-1))/cg_dist
                ! Top right-hand boundary between FEM & PC (n,m) = (dim_x,dim_z_total)
                elseif ((n==dim_x) .and. (m==dim_z_total)) then
                    dPxx(n,m) = 0.5d0*(Px(n,m-1) - Px(n-1,m-1))/cg_dist
                    dPzz(n,m) = -0.5d0*(Pz(n,m-1) + Pz(n-1,m-1))/cg_dist
                ! Right-hand edge (n,m) = (dim_x,m)
                elseif ((n==dim_x) .and. (m<dim_z_total) .and. (m/=1)) then
                    dPxx(n,m) = 0.5d0*(Px(n,m) - Px(n-1,m) + Px(n,m-1) - Px(n-1,m-1))/cg_dist
                    dPzz(n,m) = 0.5d0*(Pz(n,m) - Pz(n,m-1) + Pz(n-1,m) - Pz(n-1,m-1))/cg_dist
                ! Bottom right-hand edge (n,m) = (dim_x,1)
                elseif ((n==dim_x) .and. (m==1)) then
                    dPxx(n,m) = 0.5d0*(Px(n,m) - Px(n-1,m))/cg_dist
                    dPzz(n,m) = 0.5d0*(Pz(n,m) + Pz(n-1,m))/cg_dist
                ! Bottom edge (n,m) = (n,1)
                elseif ((n/=1) .and. (n/=dim_x) .and. (m==1)) then
                    dPxx(n,m) = 0.5d0*(Px(n,m) - Px(n-1,m))/cg_dist
                    dPzz(n,m) = 0.5d0*(Pz(n,m) + Pz(n-1,m))/cg_dist
                ! Bulk value n /= 1 or dim_x, m /= 1, dim_z_total
                else
                    dPxx(n,m) = 0.5d0*(Px(n,m) - Px(n-1,m) + Px(n,m-1) - Px(n-1,m-1))/cg_dist
                    dPzz(n,m) = 0.5d0*(Pz(n,m) - Pz(n,m-1) + Pz(n-1,m) - Pz(n-1,m-1))/cg_dist
                end if
                !print "(a4,i3,a6,i3,a9,e21.15,a9,e21.15)", "n = ", n, ", m = ", m, ", dPxx = ", dPxx(n,m), ", dPzz = ", dPzz(n,m)
                !pause
            end do
        end do

    end subroutine background_charge

    subroutine solve_phi(potential, rho)
        ! Implement the relaxation algorithm
        implicit none

        ! Variables
        double precision, dimension(:,:) :: potential, rho
        double precision, dimension(:,:), allocatable :: temp_potential
        double precision :: max_conv_p

        allocate (temp_potential(dim_x,dim_z_total))

        ! Reset temp potentials
        temp_potential = 0.0d0

        ! The Relaxation Method

        ! Initial value to start loop
        max_conv_p = eps_p*1.0d100

        do while (max_conv_p > eps_p)

            ! Calculate new value for the potential using a temporary variable.
            do m = 1, dim_z_total
                do n = 1, dim_x                
                    ! Bottom Left Corner - (n,m) = (1,1)
                    if ((n==1) .and. (m==1)) then
                        temp_potential(n,m) = potential_bottom
                    ! Left-hand edge - (n,m) = (1,m)
                    elseif ((n==1) .and. (m /= 1) .and. (m /= dim_z_total)) then
                        temp_potential(n,m) = temp_pot(n, m, potential, rho(n,m),1)
                    ! Top left corner (n,m) = (1,dim_z)
                    elseif ((n==1) .and. (m==dim_z_total)) then
                        temp_potential(n,m) = potential_top
                    ! Top boundary (n,m) = (n, dim_z)
                    elseif ((n /= 1) .and. (n /= dim_x) .and. (m==dim_z_total)) then
                        temp_potential(n,m) = potential_top
                    ! Right-hand Corner - (n,m) = (dim_x,dim_z)
                    elseif ((n==dim_x) .and. (m==dim_z_total)) then
                        temp_potential(n,m) = potential_top
                    ! Right-hand edge - (n,m) = (dim_x,m)
                    elseif ((n==dim_x) .and. (m /= 1) .and. (m /= dim_z_total)) then
                        temp_potential(n,m) = temp_pot(n, m, potential, rho(n,m),2)
                    ! Bottom right corner - (n,m) = (dim_x,1)
                    elseif ((n==dim_x) .and. (m==1)) then
                        temp_potential(n,m) = potential_bottom
                    ! Bottom edge - (n,m) = (n,1)
                    elseif ((n/=1) .and. (n/=dim_x) .and. (m==1)) then
                        temp_potential(n,m) = potential_bottom
                    ! Bulk Case - (n,m) where n \= 1, dim_x, and m \= 1, dim_z
                    else
                        temp_potential(n,m) = temp_pot(n, m, potential, rho(n,m),0)
                    end if
                    !print "(a16,i3,a2,i3,a4,e21.15)", "Temp Potential (", n, ", ", m, ") = ", temp_potential(n,m)
                    !pause
                end do
            end do

            ! Determine the maximum relative difference in the change in potential
            max_conv_p = MAXVAL(ABS(temp_potential - potential))/MAXVAL(ABS(temp_potential))

            !print *, "Max Conv P = ", max_conv_p

            ! Set the potential to the temporary potential.  If they have converged, then the loop will stop. 
            potential = temp_potential
      
        end do

        deallocate (temp_potential)

    end subroutine solve_phi

    subroutine Solve_Potential()
        implicit none

        !print *, "Start Solve_Potential."

        ! Allocate memory

        allocate (rho_total(dim_x,dim_z_total))
        allocate (chi_mat(dim_x,dim_z_total)) ! Secondary lattice.  Reduce z-dimension by 1.
        allocate (dPxx(dim_x,dim_z_total)) ! Primary lattice.  
        allocate (dPzz(dim_x,dim_z_total)) ! Primary lattice. 
        allocate (Px(dim_x,dim_z_total))
        allocate (Pz(dim_x,dim_z_total)) 
        allocate (dQxx(dim_x,dim_z_total))   
        allocate (dQzz(dim_x,dim_z_total))
        
        ! Set the old potential to the current potential before the latter is updated
        potential_old = potential
        
        ! Initialize the susceptibility.  We assume that the dielectric constant is defined on the secondary lattice with E and P.  

        do m = 1, dim_z_total ! Secondary lattice.  Values above dim_z_total-1 won't be used.
            do n = 1, dim_x            
                if (m<dim_z_FEM) then
                    chi_mat(n,m) = chi_FEM
                else
                    chi_mat(n,m) = chi_PC
                end if
            end do
        end do

        ! Our polarizations P_x and P_z only account for the polarization associated with the ferroelectric effect.  We need to
        ! generate the background polarization as well.  We can combine these into new variables, Px and Pz.  

        ! NOTE - P_x and P_z must only be non-zero in the FEM.  They must be zero whenever m>=dim_z_FEM

        ! Let's build a quick error check for this requirement.
        do m = dim_z_FEM, dim_z_total 
            do n = 1, dim_x
                if ((P_x(n,m)/=0.0d0) .or. (P_z(n,m)/=0.0d0)) then
                    print *, "Polarizations are non-zero in the PC."
                    !pause
                end if
            end do
        end do

        ! Initialize Px and Pz
        Px = 0.0d0
        Pz = 0.0d0
        
        ! Calculate the background polarization - Pbx, Pbz
        Pbx = perm_of_FS*chi_mat*E_x
        Pbz = perm_of_FS*chi_mat*E_z

        ! Mix the new and old background polarizations.  This is effectively the same as mixing the new and old electric fields.
        ! We set mixing_background = 1.0d0 to skip this simple mixing.
        Pbx = mixing_background*Pbx + (1.0d0 - mixing_background)*Pbx_old
        Pbz = mixing_background*Pbz + (1.0d0 - mixing_background)*Pbz_old
        
        ! Set old background polarization
        Pbx_old = Pbx
        Pbz_old = Pbz        

        ! Apply boundary conditions on the different potentials
        potential(:,1) = potential_bottom
        potential(:,dim_z_total) = potential_top
        potential_free(:,1) = potential_bottom
        potential_free(:,dim_z_total) = potential_top
        potential_back(:,1) = potential_bottom
        potential_back(:,dim_z_total) = potential_top
        potential_ferro(:,1) = potential_bottom
        potential_ferro(:,dim_z_total) = potential_top

        ! Solve for different components of the potential

        ! Potential due to free charge
        call solve_phi (potential_free, rho_xy + rho_xz + rho_yz)
        
        ! Potential due to background polarization
        dPxx = 0.0d0
        dPzz = 0.0d0
        
        call background_charge(Pbx, Pbz, dPxx, dPzz)
        call solve_phi (potential_back, -dPxx - dPzz)
        
        ! Potential due to ferroelectric polarization
        dPxx = 0.0d0
        dPzz = 0.0d0

        call background_charge(P_x, P_z, dPxx, dPzz)
        call solve_phi (potential_ferro, -dPxx - dPzz)
        
        ! The true potential is the sum of all of these components
        potential = potential_back + potential_ferro + potential_free

        ! Alternately, we can solve for the total potential directly by adding together the free and bound charge densities.

        ! Calculate Px and Pz.  These are the total polarizations.
        !Px = P_x + Pbx
        !Pz = P_z + Pbz

        ! Initialize dPxx and dPzz
        !dPxx = 0.0d0
        !dPzz = 0.0d0

        !call background_charge(Px,Pz, dPxx, dPzz)
        !call solve_phi(potential,rho_xy + rho_xz + rho_yz - dPxx - dPzz)

        ! Enforce boundary conditions on the potential
        potential(:,1) = potential_bottom
        potential(:,dim_z_total) = potential_top

        ! Re-initialize E_x and E_z
        E_x = 0.0d0
        E_z = 0.0d0

        ! Calculate the electric field

        do m = 1, dim_z_total
            do n = 1, dim_x
                if (m==dim_z_total) then
                    E_x(n,m) = 0.0d0
                    E_z(n,m) = 0.0d0
                elseif ((n==dim_x) .and. (m<dim_z_total)) then
                    E_x(n,m) = -0.50d0/cg_dist*(potential(1,m+1) - potential(n,m+1) + potential(1,m) - potential(n,m))
                    E_z(n,m) = -0.50d0/cg_dist*(potential(1,m+1) - potential(1,m) + potential(n,m+1) - potential(n,m))                
                else
                    E_x(n,m) = -0.50d0/cg_dist*(potential(n+1,m+1) - potential(n,m+1) + potential(n+1,m) - potential(n,m))
                    E_z(n,m) = -0.50d0/cg_dist*(potential(n+1,m+1) - potential(n+1,m) + potential(n,m+1) - potential(n,m))
                end if
            end do
        end do

        !print *, "End Solve_Potential."

        ! Deallocate memory
        deallocate (rho_total)
        deallocate (chi_mat)
        deallocate (dPxx)
        deallocate (dPzz) 
        deallocate (Px)
        deallocate (Pz)
        deallocate(dQxx)
        deallocate(dQzz)

    end subroutine Solve_Potential

end module Potential_Mod
