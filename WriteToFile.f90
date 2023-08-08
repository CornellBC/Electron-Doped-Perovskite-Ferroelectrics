Module WriteToFile

    ! This module contains does input and output to external files - data.dat, import.dat and archive.dat.

    use Coefficients  ! We need the dimensions for the FEM and PC

    implicit none

    contains

        subroutine WriteInitialize()
            ! Write parameters from Coefficients module and initial values from Initialize module to data.dat.
            implicit none

            ! Error Checking for Write-to-File
            integer :: err_status
            character (len=256) :: err_iomsg

            ! Looping
            integer :: n, m, UN
            
            ! Print the intial conditions to the data file
                
            print *, "Start printing initial conditions and coefficients to file."
            
            open(UN, file='data.dat', status='old', position='append', iostat = err_status, iomsg = err_iomsg)
            
            ! Write the Configuration/Coefficients to the data file

            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a10, a20)") "Version = ", version
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a9,l1)") "Import = ", import
            write(UN, "(a10,l1)") "Archive = ", export
            write(UN, "(a18,i5)") "Maximum Counter = ", counter_max
            write(UN, "(a20,i4)") "Archive Threshold = ", archive_threshold 
            write(UN, "(a26,a3)") "Convergence/Mixing Type = ", conv_type
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a18)") "Start coefficients"
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a19)") "Mixing Coefficients"
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a24, e21.15)") "Mixing - Polarization = ", mixing_pol
            write(UN, "(a21, e21.15)") "Mixing - Potential = ", mixing_pol
            write(UN, "(a26, e21.15)") "Mixing - Charge Density = ", mixing_rho
            write(UN, "(a35, e21.15)") "Mixing - Background Polarization = ", mixing_background
            write(UN, "(a24, e21.15)") "Epsilon - Convergence = ", eps_conv
            write(UN, "(a35)") "-----------------------------------"          
            write(UN, "(a21)") "Geometry Coefficients"
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a19, e21.15)") "Lattice distance = ", latt_dist
            write(UN, "(a27, e21.15)") "Coarse Graining Distance = ", cg_dist
            write(UN, "(a8, i4)") "Dim_x = ", dim_x
            write(UN, "(a8, i4)") "Dim_z - FEM = ", dim_z_FEM
            write(UN, "(a20, i4)") "Dim_z - Polar Cap = ", dim_z_PC
            write(UN, "(a16, i4)") "Dim_z - Total = ", dim_z_total
            write(UN, "(a8, i4)") "Dim_y = ", dim_y
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a16)") "LGD Coefficients"
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a5, e21.15)") "a1 = ", a1
            write(UN, "(a5, e21.15)") "a3 = ", a3
            write(UN, "(a6, e21.15)") "a11 = ", a11 
            write(UN, "(a6, e21.15)") "a12 = ", a12
            write(UN, "(a6, e21.15)") "g11 = ", g11
            write(UN, "(a6, e21.15)") "g44 = ", g44
            write(UN, "(a7, e21.15)") "eta1 = ", eta1
            write(UN, "(a7, e21.15)") "eta2 = ", eta2
            write(UN, "(a7, e21.15)") "eta3 = ", eta3
            write(UN, "(a7, e21.15)") "eta4 = ", eta4
            write(UN, "(a7, e21.15)") "eta5 = ", eta5
            write(UN, "(a7, e21.15)") "eta6 = ", eta6
            write(UN, "(a7, e21.15)") "eta7 = ", eta7
            write(UN, "(a7, e21.15)") "eta8 = ", eta8
            write(UN, "(a7, e21.15)") "eta9 = ", eta9
            write(UN, "(a8, e21.15)") "delta = ", delta
            write(UN, "(a10, e21.15)") "epsilon = ", eps
            write(UN, "(a8, e21.15)") "P_min = ", P_min
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a22)") "Potential Coefficients"
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a20, e21.15)") "epsilon_potential = ", eps_p
            write(UN, "(a21, e21.15)") "Potential - Bottom = ", potential_bottom
            write(UN, "(a18, e21.15)") "Potential - Top = ", potential_top
            write(UN, "(a18, e21.15)") "Chi - Polar Cap = ", chi_PC
            write(UN, "(a12, e21.15)") "Chi - FEM = ", chi_FEM
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a24)") "Schrodinger Coefficients"
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a16, e21.15)") "Hop - Stay_xy = ", hop_stay_xy
            write(UN, "(a16, e21.15)") "Hop - Stay_xz = ", hop_stay_xz
            write(UN, "(a16, e21.15)") "Hop - Stay_yz = ", hop_stay_yz
            write(UN, "(a17, e21.15)") "Hop - Parallel = ", hop_parallel
            write(UN, "(a22, e21.15)") "Hop - Perpendicular = ", hop_perp
            write(UN, "(a14, e21.15)") "Temperature = ", temp
            write(UN, "(a19, e21.15)") "epsilon_chem pot = ", eps_mu
            write(UN, "(a27)") "Charge Density Coefficients"
            write(UN, "(a22, e21.15)") "2D electron density = ", n2D
            write(UN, "(a35)") "-----------------------------------"
            write(UN, "(a16)") "End coefficients"
            write(UN, "(a35)") "-----------------------------------"
            
            ! Write the initial conditions to the data file
            write(UN, "(a24)") "Start initial conditions"

            do n = 1, dim_x
                do m = 1, dim_z_total
                    write(UN, "(i3, a1, i3, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, & 
                    a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15)") n, ",", m, ",", &
                    P_x(n,m), ",", P_z(n,m), ",", &
                    potential(n,m), ",", potential_ferro(n,m), ",", potential_back(n,m), ",", potential_free(n,m), ",",&
                    rho_xy(n,m), ",", rho_xz(n,m), ",", rho_yz(n,m), ",", & 
                    E_x(n,m), ",", E_z(n,m)
                end do
            end do 

            write(UN, "(a22)") "End initial conditions"

            close(UN)

        end subroutine WriteInitialize
    
        subroutine Write2File(counter)
            ! Write final results to data.dat.
            implicit none

            ! Loop Countain
            integer :: counter

            ! Error Checking for Write-to-File
            integer :: err_status
            character (len=256) :: err_iomsg

            ! Looping
            integer :: n, m, UN, en
                        
            print *, "Start print to file."
    
            open(newunit=UN, file='data.dat', status='old', position='append', iostat = err_status, iomsg = err_iomsg)

            write(UN, "(a26,i6)") "Start data set.  Counter =", counter

            do n = 1, dim_x
                do m = 1, dim_z_total
                    write(UN, "(i3, a1, i3, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, & 
                    a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15)") n, ",", m, ",", &
                    P_x(n,m), ",", P_z(n,m), ",", &
                    potential(n,m), ",", potential_ferro(n,m), ",", potential_back(n,m), ",", potential_free(n,m), ",",&
                    rho_xy(n,m), ",", rho_xz(n,m), ",", rho_yz(n,m), ",", & 
                    E_x(n,m), ",", E_z(n,m)
                end do
            end do

            write(UN, "(a12)") "End data set"

            print *, "End print to file."

        end subroutine Write2File

        subroutine Archive2File(counter)
            ! Write data from intermediate steps to archive.dat.  This data is always appended to preserve a history.
            ! The archive_threshold parameters determines how often the data will be written to the file.
            implicit none

            ! Loop Countain
            integer :: counter

            ! Error Checking for Write-to-File
            integer :: err_status
            character (len=256) :: err_iomsg

            ! Looping
            integer :: n, m
            integer :: UN
                        
            print *, "Start archive to file for counter =", counter

            open(newunit=UN, file='archive.dat', status='old', position='append', iostat = err_status, iomsg = err_iomsg)

            write(UN, "(a29, i6)") "Start data set for counter = ", counter

            do n = 1, dim_x
                do m = 1, dim_z_total
                    write(UN, "(i3, a1, i3, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15, & 
                    a1, e21.15, a1, e21.15, a1, e21.15, a1, e21.15)") n, ",", m, ",", &
                    P_x(n,m), ",", P_z(n,m), ",", &
                    potential(n,m), ",", potential_ferro(n,m), ",", potential_back(n,m), ",", potential_free(n,m), ",",&
                    rho_xy(n,m), ",", rho_xz(n,m), ",", rho_yz(n,m), ",", & 
                    E_x(n,m), ",", E_z(n,m)
                end do
            end do

            write(UN, "(a12)") "End data set"

            close(UN)

            print *, "End archive to file for counter = ", counter

        end subroutine Archive2File

        subroutine ImportFromFile()
            ! Used by the Initialize module to import data from the import.dat file.
            implicit none

            ! Error Checking for Write-to-File
            integer :: err_status
            character (len=256) :: err_iomsg

            ! Looping
            integer :: n, m, UN
            integer :: x = 0 ! x position that we import from a data file
            integer :: z = 0 ! z position that we import from a data file
         
            ! Import the initial conditions for polarization, potential and charge density from a data file
            ! Normally this will be a subset of the data contained in the archive.dat or data.dat files that will act as the new
            ! input for another run.

            print *, "Start import from file." 

            open (newunit=UN, file = 'import.dat', status = 'old', iostat = err_status, iomsg = err_iomsg)

            do n = 1, dim_x
                do m = 1, dim_z_total
                    read(UN, *) x,z, P_x(x,z), P_z(x,z), potential(x,z), potential_ferro(x,z), potential_back(x,z), & 
                    potential_free(x,z), rho_xy(x,z), rho_xz(x,z), rho_yz(x,z), E_x(x,z), E_z(x,z)                    
                end do
            end do
        
            close (UN)

            print *, "End import from file."

        end subroutine ImportFromFile

end module WriteToFile
