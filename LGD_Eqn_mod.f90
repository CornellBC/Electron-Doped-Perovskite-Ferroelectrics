module LGDEqn
    ! This module relies on the "Coefficients" module for many of the coefficients.
    ! This version of the model does not assume that we wrap the model in the z-direction.  We only impose that the derivative
    ! goes to zero along that boundary. 
    ! The Free Energy Equation has been changed to match Bill's.   
    
    use Physical_Constants
    use Coefficients
    use Anderson_Mix

    implicit none

    !Variables - Defined for the entire module
    logical, private,save :: conv_lgd ! Converged when flag is true.
    double precision, dimension(:,:), allocatable, private, save :: g_x, g_z ! Defined on the secondary lattice
    double precision, private :: conv_max ! Maximum value of the convergence array
    double precision, dimension(:,:), allocatable, private, save :: d2Pxx, d2Pxz, d2Pzz, d2Pzx ! Defined on the secondary lattice
    double precision, dimension(:,:), allocatable, private, save :: rho_xy_sec, rho_xz_sec, rho_yz_sec 
    ! Charge density on the secondary lattice
    double precision, dimension(:,:), allocatable, private, save :: drho_xy_x_sec, drho_xz_x_sec, drho_yz_x_sec ! Derivative of the 
    ! charge density on the secondary lattice WRT x 
    double precision, dimension(:,:), allocatable, private, save :: drho_xy_z_sec, drho_xz_z_sec, drho_yz_z_sec ! Derivative of the  
    ! charge density on the secondary lattice WRT z
    double precision, dimension(:,:), allocatable, private, save :: Hess_x, Hess_z 
    double precision, dimension(:,:), allocatable, private, save :: chi_x, chi_z 
    ! In this model, each unit cell has an x and z coordinate.  We will change this to one dimension using
    ! (n,m) = n + (m-1)*dim_x.  We will need to do this for the Hessian.
    ! dPxx = derivative of P_x WRT x, dPxz = derivative of P_x WRT z. 
    ! d2Pxx = second derivative of P_x WRT x, d2Pxz = second derivative of P_x WRT z
    integer, private :: n, m  ! Looping variables
    
    ! Electric Potential
    double precision, dimension(:,:), allocatable :: D_x, D_z ! Electric displacement defined on the secondary lattice.

    ! Anderson Mixing
    double precision, dimension(:), allocatable :: X, Y, P_Anderson  ! For Anderson Mixing of polarzation.  Not used.

    !Constraints/Boundary Conditions

    !1.  We must have that dP_x/dz = 0 at z = 0, and Lz.  In FORTAN, our first unit cell is at i = 1, and our last one is at 
    !i = dim_x.  So, we must have that dP_x/dz(n,1) = dP_x/dz(n,dim_z) = 0 for all 1 <= n <= dim_x.  We assume that P(n,0) = P(n,1)
    ! and P(n,dim_z) = P(n, dim_z+1) for all 1 <= n <= dim_x when approximating the derivatives.
    !2.  Similarly, we must have that dP_z/dz(n,1) = dP_z/dz(n, dim_z) for all 1 <= n <= dim_x.  
    !3.  Since the material repeats in the x-direction, we can wrap around the values for the polarization to help determine 
    ! the derivatives.

contains
    double precision function dLGD_x_E(P1, P2, d2P11, d2P12, rho1, rho2, rho3, drho11, drho21, drho31, DF) result(W)
        ! First variational derivative of the free energy WRT P_x
        implicit none
        double precision, intent(in) :: P1, P2, d2P11, d2P12 ! Polarizations and their second derivatives
        double precision, intent(in) :: rho1, rho2, rho3, drho11, drho21, drho31 ! Charge densities and their first derivatives
        double precision, intent(in) :: DF ! Electric potential
            W = 2.0d0*a1*P1 + 4.0d0*a11*(P1**3) + 2.0d0*a12*P1*(P2**2) + eta1*P1*(rho1 + rho2) + eta2*P1*rho3 &
            -1.0d0*eta3*(drho11 + drho21) - 1.0d0*eta4*drho31 &
            -1.0d0*d2P11*(g11 + eta5*(rho1 + rho2) + eta9*rho3) &
            -1.0d0*d2P12*(g44 + eta6*rho2 + eta7*rho1 + eta8*rho3) &
            -1.0d0/perm_of_FS*(DF - P1)
         
    end function dLGD_x_E

    double precision function dLGD_z_E(P1, P2, d2P11, d2P12, rho1, rho2, rho3, drho11, drho21, drho31, DF) result(W)
        ! First variational derivative of the free energy WRT P_z
        implicit none
        double precision, intent(in) :: P1, P2, d2P11, d2P12 ! Polarizations and their second derivatives
        double precision, intent(in) :: rho1, rho2, rho3, drho11, drho21, drho31 ! Charge densities and their first derivatives
        double precision, intent(in) :: DF ! Electric potential
            W = 2.0d0*a3*P1 + 4.0d0*a11*(P1**3) + 2.0d0*a12*P1*(P2**2) + eta1*P1*(rho1 + rho2) + eta2*P1*rho3 &
            -1.0d0*eta3*(drho11 + drho21) - 1.0d0*eta4*drho31 &
            -1.0d0*d2P11*(g11 + eta5*(rho1 + rho2) + eta9*rho3) &
            -1.0d0*d2P12*(g44 + eta6*rho1 + eta7*rho2 + eta8*rho3) &
            -1.0d0/perm_of_FS*(DF - P1)  
            ! Note:  The eta6 and eta7 terms are reversed from dLGD_x.                               
        end function dLGD_z_E

    double precision function d2LGD_x_E(P1, P2, rho1, rho2, rho3, switch) result(W)
        ! Second variational derivative WRT P_x
        implicit none
        double precision, intent(in) :: P1, P2 ! Polarizations
        double precision, intent(in) :: rho1, rho2, rho3 ! Charge densities and their first derivatives
        integer, intent(in) :: switch
        if (switch==0) then ! Bulk Case
            W = 2.0d0*a1 + 12.0d0*a11*(P1**2) + 2.0d0*a12*(P2**2) + eta1*(rho1 + rho2) + eta2*rho3 &
            + 2.0d0/(cg_dist**2)*(g11 + eta5*(rho1 + rho2) + eta9*rho3) &
            + 2.0d0/(cg_dist**2)*(g44 + eta6*rho2 + eta7*rho1 + eta8*rho3) &
            + 1.0d0/perm_of_FS  
        elseif (switch==1) then ! (m==dim_z_FEM-1) or (m==1)
            W = 2.0d0*a1 + 12.0d0*a11*(P1**2) + 2.0d0*a12*(P2**2) + eta1*(rho1 + rho2) + eta2*rho3 &
            + 2.0d0/(cg_dist**2)*(g11 + eta5*(rho1 + rho2) + eta9*rho3) &
            + 1.0d0/(cg_dist**2)*(g44 + eta6*rho2 + eta7*rho1 + eta8*rho3) & 
            + 1.0d0/perm_of_FS 
        end if 
    end function d2LGD_x_E

    double precision function d2LGD_z_E(P1, P2, rho1, rho2, rho3, switch) result(W)
        ! Second variational derivative WRT P_z
        implicit none
        double precision, intent(in) :: P1, P2 ! Polarizations
        double precision, intent(in) :: rho1, rho2, rho3 ! Charge densities and their first derivatives
        integer, intent(in) :: switch
        if (switch==0) then ! Bulk Case
            W = 2.0d0*a3 + 12.0d0*a11*(P1**2) + 2.0d0*a12*(P2**2) + eta1*(rho1 + rho2) + eta2*rho3 &
            + 2.0d0/(cg_dist**2)*(g11 + eta5*(rho1 + rho2) + eta9*rho3) &
            + 2.0d0/(cg_dist**2)*(g44 + eta6*rho1 + eta7*rho2 + eta8*rho3) &
            + 1.0d0/perm_of_FS
        elseif (switch==1) then ! (m==dim_z_FEM-1) or (m==1)
            W = 2.0d0*a3 + 12.0d0*a11*(P1**2) + 2.0d0*a12*(P2**2) + eta1*(rho1 + rho2) + eta2*rho3 &
            + 1.0d0/(cg_dist**2)*(g11 + eta5*(rho1 + rho2) + eta9*rho3) &
            + 2.0d0/(cg_dist**2)*(g44 + eta6*rho1 + eta7*rho2 + eta8*rho3) &
            + 1.0d0/perm_of_FS 
        end if 
    end function d2LGD_z_E
    
    integer function pos(n,m,dim) result(W)
        ! Convert between (n,m) and the single coordinate i
        implicit none
        integer, intent(in) :: n,m,dim
        W = n + (m-1)*dim
    end function
   
    subroutine Solve_LGD()
        implicit none

        ! Allocate memory for the arrays

        allocate (g_x(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (g_z(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (d2Pxx(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (d2Pxz(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (d2Pzz(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (d2Pzx(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (rho_xy_sec(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (rho_xz_sec(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (rho_yz_sec(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (drho_xy_x_sec(dim_x, dim_z_FEM-1)) ! Secondary lattice
        allocate (drho_xz_x_sec(dim_x, dim_z_FEM-1)) ! Secondary lattice
        allocate (drho_yz_x_sec(dim_x, dim_z_FEM-1)) ! Secondary lattice
        allocate (drho_xy_z_sec(dim_x, dim_z_FEM-1)) ! Secondary lattice
        allocate (drho_xz_z_sec(dim_x, dim_z_FEM-1)) ! Secondary lattice
        allocate (drho_yz_z_sec(dim_x, dim_z_FEM-1)) ! Secondary lattice
        allocate (Hess_x(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (Hess_z(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (chi_x(dim_x,dim_z_FEM-1)) ! Secondary lattice
        allocate (chi_z(dim_x,dim_z_FEM-1)) ! Secondary lattice

        ! Electric Potential
        allocate (D_x(dim_x,dim_z_total)) ! Secondary lattice
        allocate (D_z(dim_x,dim_z_total)) ! Secondary lattice

        !Check this module is being called.

        !print *,"Starting LGD solve."

        ! Apply the boundary condition.  Set P_x and P_z = 0.0d0 along the bottom
        P_x(:,1) = 0.0d0
        P_z(:,1) = 0.0d0
        
        ! Set the old polarizations before P_x and P_z are updated
        P_x_old = P_x
        P_z_old = P_z
        
        ! Convert the charge densities (rho_xy, rho_xz, rho_yz) to the secondary lattice.

        call Convert_Primary_to_Secondary(rho_xy, rho_xy_sec)
        call Convert_Primary_to_Secondary(rho_xz, rho_xz_sec)
        call Convert_Primary_to_Secondary(rho_yz, rho_yz_sec)
        
        ! Calculate the derivatives of the charge densities on the secondary lattice
        ! We do these outside of the loop because they don't change once they're calculated.

        ! drho WRT x

        do m = 1, dim_z_FEM-1 ! Secondary lattice
            do n = 1, dim_x
                if (n==dim_x) then
                    drho_xy_x_sec(n,m) = 0.5d0*(rho_xy(1,m) - rho_xy(dim_x,m) + rho_xy(1,m+1) - rho_xy(dim_x,m+1))/cg_dist
                    drho_xz_x_sec(n,m) = 0.5d0*(rho_xz(1,m) - rho_xz(dim_x,m) + rho_xz(1,m+1) - rho_xz(dim_x,m+1))/cg_dist
                    drho_yz_x_sec(n,m) = 0.5d0*(rho_yz(1,m) - rho_yz(dim_x,m) + rho_yz(1,m+1) - rho_yz(dim_x,m+1))/cg_dist
                else
                    drho_xy_x_sec(n,m) = 0.5d0*(rho_xy(n+1,m) - rho_xy(n,m) + rho_xy(n+1,m+1) - rho_xy(n,m+1))/cg_dist
                    drho_xz_x_sec(n,m) = 0.5d0*(rho_xz(n+1,m) - rho_xz(n,m) + rho_xz(n+1,m+1) - rho_xz(n,m+1))/cg_dist
                    drho_yz_x_sec(n,m) = 0.5d0*(rho_yz(n+1,m) - rho_yz(n,m) + rho_yz(n+1,m+1) - rho_yz(n,m+1))/cg_dist
                end if
            end do
        end do

        ! drho WRT z

        do m = 1, dim_z_FEM-1 ! Secondary lattice
            do n = 1, dim_x
                drho_xy_z_sec(n,m) = 0.5d0*(rho_xy(n,m+1) - rho_xy(n,m) + rho_xy(n+1,m+1) - rho_xy(n+1,m))/cg_dist
                drho_xz_z_sec(n,m) = 0.5d0*(rho_xz(n,m+1) - rho_xz(n,m) + rho_xz(n+1,m+1) - rho_xz(n+1,m))/cg_dist
                drho_yz_z_sec(n,m) = 0.5d0*(rho_yz(n,m+1) - rho_yz(n,m) + rho_yz(n+1,m+1) - rho_yz(n+1,m))/cg_dist
            end do
        end do
                    
        ! Calculate the electric potential
        
        ! Reset values
        D_x = 0.0d0
        D_z = 0.0d0
        
        ! Bill's Method - No eps_inf
        D_x = perm_of_FS*E_x + P_x
        D_z = perm_of_FS*E_z + P_z

        !Gradient Newton Galerkin Algorithm (GNGA)
        
        !Step 0 - Establish the looping of the algorithm

        conv_lgd = .false.

        do while (conv_lgd .eqv. .false.)
          
            !Step 0.5 - Initialize the 2nd order derivatives of P_x and P_z WRT x and z
            
            ! Zeroize the 2nd order derivatives
            d2Pxx = 0.0d0
            d2Pxz = 0.0d0
            d2Pzz = 0.0d0
            d2Pzx = 0.0d0
            
            !Create 2nd derivative - P_x WRT x.
            ! The first derivative takes us from secondary -> primary.  The second derivative takes us from primary -> secondary.
            
            do m = 1,dim_z_FEM-1 !  This derivative goes from secondary lattice to secondary lattice.
                do n = 1,dim_x                
                    if (n==1) then
                        d2Pxx(n,m) = (P_x(2,m) + P_x(dim_x,m)- 2.0d0*P_x(1,m))/(cg_dist**2)                
                    elseif (n==dim_x) then
                        d2Pxx(n,m) = (P_x(1,m) + P_x(dim_x-1,m) - 2.0d0*P_x(dim_x,m))/(cg_dist**2)
                    else
                        d2Pxx(n,m) = (P_x(n+1,m) + P_x(n-1,m) - 2.0d0*P_x(n,m))/(cg_dist**2)
                    end if
                    !print "(a4,i3,a6,i3,a10,e21.15)", "n = ", n, ", m = ", m, ", d2Pxx = ", d2Pxx(n,m)
                    !pause
                end do
            end do

            !Create 2nd derivative - P_z WRT z
            
            do m = 1,dim_z_FEM-1 ! Secondary lattice.
                do n = 1,dim_x                
                    if (m==1) then  ! We assume that the derivative across the boundary is 0, so P_z(n,1) = P_z(n,0)
                        d2Pzz(n,m) = (P_z(n,2) - 1.0d0*P_z(n,1))/(cg_dist**2)                
                    elseif (m==dim_z_FEM-1) then !    
                        d2Pzz(n,m) = (P_z(n,dim_z_FEM-2) - 1.0d0*P_z(n,dim_z_FEM-1))/(cg_dist**2) ! We must have the boundary
                        ! condition - dP_x/d_z = 0 at the boundary => P_z(n,dim_z_FEM-1) = P_z(n,dim_x_FEM)
                    else
                        d2Pzz(n,m) = (P_z(n,m+1) + P_z(n,m-1) - 2.0d0*P_z(n,m))/(cg_dist**2)
                    end if
                    !print "(a4,i3,a6,i3,a10,e21.15)", "n = ", n, ", m = ", m, ", d2Pzz = ", d2Pzz(n,m)
                    !pause
                end do
            end do

            !Create the cross-derivatives - P_x WRT z

            do m = 1,dim_z_FEM-1 ! Secondary lattice.
                do n = 1,dim_x
                    if (m==1) then
                        d2Pxz(n,m) = (P_x(n,m+1) - 1.0d0*P_x(n,m))/(cg_dist**2)                
                    elseif (m==dim_z_FEM-1) then ! The polarization in the PC must be zero.
                        d2Pxz(n,m) = (P_x(n,m-1) - 1.0d0*P_x(n,m))/(cg_dist**2) ! We must have the boundary
                        ! condition - dP_x/d_z = 0.  Check that this is true, and not that P_x = 0 in the FEM.
                    else
                        d2Pxz(n,m) = (P_x(n,m+1) + P_x(n,m-1) - 2.0d0*P_x(n,m))/(cg_dist**2)
                    end if
                    !print "(a4,i3,a6,i3,a10,e21.15)", "n = ", n, ", m = ", m, ", d2Pxz = ", d2Pxz(n,m)
                    !pause
                end do
            end do

            !Create the cross-derivatives - P_z WRT x

            do m = 1,dim_z_FEM-1 ! Secondary lattice.
                do n = 1,dim_x                
                    if (n==1) then
                        d2Pzx(n,m) = (P_z(2,m) + P_z(dim_x,m)- 2.0d0*P_z(1,m))/(cg_dist**2)                
                    elseif (n==dim_x) then
                        d2Pzx(n,m) = (P_z(1,m) + P_z(dim_x-1,m) - 2.0d0*P_z(dim_x,m))/(cg_dist**2)
                    else
                        d2Pzx(n,m) = (P_z(n+1,m) + P_z(n-1,m) - 2.0d0*P_z(n,m))/(cg_dist**2)
                    end if
                    !print "(a4,i3,a6,i3,a10,e21.15)", "n = ", n, ", m = ", m, ", d2Pzx = ", d2Pzx(n,m)
                    !pause
                end do
            end do            

            ! Take the variational derivative of the free energy with WRT P_x and P_z

            ! Zeroize the derivatives
            g_x = 0.0d0
            g_z = 0.0d0

            do m = 1, dim_z_FEM-1 ! Secondary lattice
                do n = 1, dim_x                
                    if (m==1) then
                        g_x(n,m) = 0.0d0
                        g_z(n,m) = 0.0d0
                    else
                        g_x(n,m) = dLGD_x_E(P_x(n,m),P_z(n,m),d2Pxx(n,m),d2Pxz(n,m),rho_xy_sec(n,m), rho_xz_sec(n,m), & 
                        rho_yz_sec(n,m), drho_xy_x_sec(n,m), drho_xz_x_sec(n,m), drho_yz_x_sec(n,m), D_x(n,m))
                        g_z(n,m) = dLGD_z_E(P_z(n,m),P_x(n,m),d2Pzz(n,m),d2Pzx(n,m),rho_xz_sec(n,m), rho_yz_sec(n,m), & 
                        rho_xy_sec(n,m), drho_xz_z_sec(n,m), drho_yz_z_sec(n,m), drho_xy_z_sec(n,m), D_z(n,m))                        
                    end if
                    !print "(a4,i3,a6,i3,a8,e21.15,a8,e21.15)", "n = ", n, ", m = ", m, ", g_x = ", g_x(n,m), ", g_z = ", g_z(n,m)
                    !pause
                end do    
            end do

            !Zeroize the Hessian.  We will populate it with non-zero values later.
            
            Hess_x = 0.0d0
            Hess_z = 0.0d0

            !Step 2 - Calculate the Hessian for each cell.
            
            ! We calculate only the diagonal terms in the Hessian, and we apply them to the corresponding (n,m) location. 
            do m = 1, dim_z_FEM-1 ! Secondary Lattice
                do n = 1, dim_x
                    if ((m==1) .or. (m==dim_z_FEM-1)) then
                        Hess_x(n,m) = d2LGD_x_E(P_x(n,m),P_z(n,m), rho_xy_sec(n,m), rho_xz_sec(n,m), & 
                        rho_yz_sec(n,m), 1)
                        Hess_z(n,m) = d2LGD_z_E(P_z(n,m),P_x(n,m), rho_xz_sec(n,m), rho_yz_sec(n,m), &
                        rho_xy_sec(n,m), 1)
                    else
                        Hess_x(n,m) = d2LGD_x_E(P_x(n,m),P_z(n,m), rho_xy_sec(n,m), rho_xz_sec(n,m), & 
                        rho_yz_sec(n,m), 0)
                        Hess_z(n,m) = d2LGD_z_E(P_z(n,m),P_x(n,m), rho_xz_sec(n,m), rho_yz_sec(n,m), &
                        rho_xy_sec(n,m), 0)
                    end if
                    !print "(a4,i3,a6,i3,a11,e21.15,a11,e21.15)", "n = ", n, ", m = ", m, ", Hess_x = ", &
                    !Hess_x(n,m), ", Hess_z = ", Hess_z(n,m)
                    !pause
                end do
            end do
            
            ! Define chi using the 1st and 2nd derivatives
            
            ! Zeroize inital values
            chi_x = 0.0d0
            chi_z = 0.0d0
            
            ! Calculate chi 
            chi_x = g_x/Hess_x
            chi_z = g_z/Hess_z
            
            ! Step 4 - Update values for the polarization
                        
            do m = 1, dim_z_total
                do n = 1, dim_x
                    if (m <= dim_z_FEM-1) then
                        P_x(n,m) = P_x(n,m) - delta*chi_x(n,m)
                        P_z(n,m) = P_z(n,m) - delta*chi_z(n,m)
                    else
                        P_x(n,m) = 0.0d0
                        P_z(n,m) = 0.0d0
                    end if
                    !print "(a4,i3,a6,i3,a8,e21.15,a8,e21.15)", "n = ", n, ", m = ", m, ", P_x = ", P_x(n,m), ", P_z = ", P_z(n,m) 
                    !pause
                end do
            end do

            ! Step 5(a) - Determine convergence.  Consider the change in the magnitude of P = sqrt(P_x^2 + P_z^2).  
            
            conv_max = MAXVAL(SQRT((delta*chi_x)**2 + (delta*chi_z)**2))

            ! Step 6 - Test for convergence. Use the change in the magnitude.

            if (conv_max < eps*MAXVAL(SQRT(P_x**2 + P_z**2))) then  
                ! We multiply eps*MAXVAL because we want to use relative magnitudes.
                conv_lgd = .true.       
            elseif (MAXVAL(SQRT(P_x**2 + P_z**2))<= P_min) then
                ! This is used to keep the polarization from converging to zero, which can take a very long time.
                conv_lgd = .true.
                print *, "P is tiny."
            end if

        end do

        ! Simple Mixing - We set mixing_pol = 1.0d0 in the Coefficients to disable this.
        P_x = mixing_pol*P_x + (1.0d0 - mixing_pol)*P_x_old
        P_z = mixing_pol*P_z + (1.0d0 - mixing_pol)*P_z_old

        ! Anderson Mixing - This code allows for the Anderson Mixing of the polarization, but wasn't used.
        !allocate (X(dim_x*dim_z_total))
        !allocate (Y(dim_x*dim_z_total))
        !allocate (P_Anderson(dim_x*dim_z_total))
        
        !X = 0.0d0
        !Y = 0.0d0
        !P_Anderson = 0.0d0
                                
        !X = reshape(P_x_old, (/ dim_x*dim_z_total /))
        !Y = reshape(P_x, (/ dim_x*dim_z_total /))
        !call AM(X, Y, P_Anderson)
        !P_x = reshape(P_Anderson, (/ dim_x, dim_z_total /))

        !X = 0.0d0
        !Y = 0.0d0
        !P_Anderson = 0.0d0
                                
        !X = reshape(P_z_old, (/ dim_x*dim_z_total /))
        !Y = reshape(P_z, (/ dim_x*dim_z_total /))
        !call AM(X, Y, P_Anderson)
        !P_z = reshape(P_Anderson, (/ dim_x, dim_z_total /))

        !deallocate (X)
        !deallocate (Y)
        !deallocate (P_Anderson)       

        !print *, "End LGD solve."

        ! Deallocate memory for these arrays 
        deallocate (g_x)
        deallocate (g_z)
        deallocate (d2Pxx)
        deallocate (d2Pxz)
        deallocate (d2Pzz)
        deallocate (d2Pzx)
        deallocate (rho_xy_sec) 
        deallocate (rho_xz_sec)
        deallocate (rho_yz_sec) 
        deallocate (drho_xy_x_sec) 
        deallocate (drho_xz_x_sec) 
        deallocate (drho_yz_x_sec)
        deallocate (drho_xy_z_sec) 
        deallocate (drho_xz_z_sec) 
        deallocate (drho_yz_z_sec) 
        deallocate (Hess_x)
        deallocate (Hess_z)
        deallocate (chi_x)
        deallocate (chi_z)

        ! Electric Potential
        deallocate (D_x)
        deallocate (D_z)

    end subroutine Solve_LGD

end module LGDEqn