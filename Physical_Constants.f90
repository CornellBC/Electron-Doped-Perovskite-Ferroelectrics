module Physical_Constants

    implicit none

    ! Physical Constants
    double precision, parameter :: funda_charge = 1.60217663d-19 ! Fundamental Charge, Units - Coulombs (C)
    double precision, parameter :: elec_mass = 9.1093837015d-31 ! Electron Mass, Units - kg
    double precision, parameter :: prot_mass = 1.67262191d-27 ! Proton Mass, Units - kg
    double precision, parameter :: neut_mass = 1.67492749804d-27 ! Neutron Mass, Units - kg
    double precision, parameter :: perm_of_FS = 8.85418782d-12 ! Permittivity of Free Space, Units - m^-3 kg^-1 s^4 A^2
    double precision, parameter :: h = 6.62607015d-34 ! Planck's Constant, Units - m^2 kg s^-1
    double precision, parameter :: J_to_eV = 6.24150907446d18 ! eV-per-Joule, Units = eV
    double precision, parameter :: h_ev = h*(J_to_eV)! Planck's Constant in eV, Units - 4.1357d-15 eV s^-1
    double precision, parameter :: hbar = 1.054571817d-34 ! Planck's Constant divided by 2*Pi, Units - J s or kg m^2 s^-1
    double precision, parameter :: hbar_ev = hbar*(J_to_eV) ! Planck's Constant divided by 2*Pi in eV, Units - eV s^-1
    double precision, parameter :: Boltzmann = 1.380649d-23 ! Boltzmann's Constant, Units - J / K 

    ! Geometric Constants
    double precision, parameter :: pi = 4*ATAN(1.d0)
    
end module Physical_Constants

