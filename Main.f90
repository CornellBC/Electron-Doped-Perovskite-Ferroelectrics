program Model
    ! The program is broken down as follows:
    !  a.  Main program (called Model) which orchestrates all of the steps.
    !  b.  Physical Constants Module (called Physical_Constants) contains common geometric and physical constants.  
    !  c.  Coefficients Module (called Coefficients) contains all of the coefficients and parameters related to the model.
    !  d.  Initialize Module (called Initialize) initializes the polarization, charge density, and potential values.
    !  d.  LGD Equation Module (called LGDEqn) solves for new values of the polarization given the charge density.
    !  e.  Potential Equation Module (called Potential) solves for new values of the potential given the polarization, 
    !      energies, and wavefunctions.
    !  f.  Schrodinger Equation Module (called Schrodinger) solves for new energy/eigenfunction pairs given the potential.
    !  g.  Convergence Module (called Convergence) tests for convergence.
    !  h.  WriteToFile Module (called WriteToFile) handles import and export of data with external files.
    !  i.  Anderson Mixing Module (called Anderson_Mix) implements the Anderson mixing algorithm.
    !
    ! The program has the following additional dependencies/modules:
    !  a.  BLAS Library (called BLAS.lib)
    !  b.  LAPACK Library (called LAPACK.lib)
    !  c.  (Optional) FlexiBLAS (combined BLAS and LAPACK library used by SHARCNet)
       
    use Physical_Constants
    use Coefficients
    use Initialize
    use LGDEqn
    use Potential_Mod
    use Schrodinger
    use Convergence
    use WriteToFile
    use Anderson_Mix

    implicit none

    ! Looping
    integer :: counter=0

    ! Step 1 - Initialize the Model
    ! Since the variable are public, they don't need to be passed individually to this subroutine.
    ! We only need to initialize the polarizations and charge densities.  Potential will be solved using these as input 
    ! in the next step.

    call initialize_model()

    ! We only need to solve for the potential if it's not imported in the Initialize module.
    if (import .eqv. .false.) then
        call Solve_Potential()
    end if

    if ((MODULO(counter, archive_threshold)==0) .and. (export .eqv. .true.)) then
        call Archive2File(counter)
    end if

    do while ((done .eqv. .false.) .and. (counter<counter_max))      

        counter = counter+1

        ! Step 2 - Solve the LGD Equation for a given polarization P = (P_x,P_z).  

        call Solve_LGD()
        
        ! Step 3 - Solve the Schrodinger Equation to get energies and free charge densities.
        
        call Solve_Schrodinger()

        ! Step 4 - Solve for the Potential for a given polarization P and charge density rho
        ! Since the polarization values are public, they don't need to be passed individually to this subroutine.

        call Solve_Potential()

        ! Output mid-calculation results to a data file.  This can be used as a starting point for future iterations, or for 
        ! investigating bugs. 
        
        ! Take the counter modulo the archive_threshold = 0 as the trigger.

        if ((MODULO(counter, archive_threshold)==0) .and. (export .eqv. .true.)) then
            call Archive2File(counter)
        end if

        ! Step 5 - Check convergence.  We use potential to check for convergence, and for mixing.  

        call Test_Convergence()

        ! Loop repeats if the done flag isn't true.
   
    end do

    ! Step 6 - Write Polarization and Potential values to the data file

    print *, "Write final data."

    call Write2File(counter)

    print *, "End program."

end program Model
