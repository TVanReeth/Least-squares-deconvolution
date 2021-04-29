! -*- f90 -*-
MODULE PhysicalConstants
!==============================================================================
! This file contains the values of several widely used astrophysical constants,
! in cgs units where possible, unless mentioned otherwise. As a safety precaution, 
! the used units are given for each constant if they have any.


    IMPLICIT NONE
    
    real(selected_real_kind(14)),    parameter :: pi = 3.14159265358979312 
    
    real(selected_real_kind(8,5)),   parameter :: c = 2.99792458d5         !light velocity (km/s)
    real(selected_real_kind(5,27)),  parameter :: h = 6.62620d-27          !Planck cst (erg s)
    real(selected_real_kind(6,16)),  parameter :: k_b = 1.380662d-16       !Bolzmann cst (erg K^-1)
    real(selected_real_kind(5,5)),   parameter :: sigma = 5.66956d-5       !Stefan-Bolzmann cst (erg cm^-2 s^-1 K^-4)
    real(selected_real_kind(3,8)),   parameter :: G = 6.673d-8             !gravitation cst (dyn cm^2 g^-2)
    real(selected_real_kind(6,10)),  parameter :: qe = 4.803203d-10        !electron charge (esu)
    real(selected_real_kind(6,28)),  parameter :: me = 9.109382d-28        !electron mass (g)
    real(selected_real_kind(6,24)),  parameter :: mp = 1.672622d-24        !proton mass (g)
    real(selected_real_kind(6,24)),  parameter :: mn = 1.674927d-24        !neutron mass (g)
    real(selected_real_kind(7,23)),  parameter :: N0 = 6.0221367d23        !Avogadro number (mol^-1)
    real(selected_real_kind(5,7)),   parameter :: R = 8.31441d7            !Universal gas constant (erg mol^-1 K^-1)
    real(selected_real_kind(5,12)),  parameter :: k_e = 8.98755d12         !Coulomb cst (dyn cm^2 g^-2)
    real(selected_real_kind(8,5)),   parameter :: RH = 1.09677583d5        !Rydberg cst for 1H (cm^-1)
    real(selected_real_kind(14,20)), parameter :: mu0 = 1.3981972969d-20   !relative permeability of vacuum (dyn s^2 esu^-2)
    real(selected_real_kind(8,2)),   parameter :: eps0 = 7.9577472d-2      !relative permittivity of vacuum (s^4 esu^2 s^-2 cm^-3 g^-1)
    
    real(selected_real_kind(4,33)),  parameter :: M_sol = 1.9891d33        !Solar mass (g)
    real(selected_real_kind(3,6)),   parameter :: F_sol = 1.365d6          !Solar irradiance (erg s^-1 cm^-2)
    real(selected_real_kind(3,33)),  parameter :: L_sol = 3.839d33         !Solar luminosity (erg s^-1)
    real(selected_real_kind(5,10)),  parameter :: R_sol = 6.95508d10       !Solar radius (cm)
    real(selected_real_kind(3,3)),   parameter :: T_sol = 5.777d3          !Solar effective temperature (K)
    real(selected_real_kind(2,0)),   parameter :: M_bol_abs = 4.74d0       !Solar absolute bolometric magnitude (magn)
    real(selected_real_kind(3,1)),   parameter :: m_bol_app = -2.683d1     !Solar apparent bolometric magnitude (magn)
    real(selected_real_kind(3,1)),   parameter :: m_U_app   = -2.591d1     !Solar apparent ultraviolet magnitude (magn)
    real(selected_real_kind(3,1)),   parameter :: m_B_app   = -2.610d1     !Solar apparent blue magnitude (magn)
    real(selected_real_kind(3,1)),   parameter :: m_V_app   = -2.675d1     !Solar apparent visual magnitude (magn)
    real(selected_real_kind(2,0)),   parameter :: BC_sol = -0.08d0         !Solar bolometric coorection (magn)
    real(selected_real_kind(4,27)),  parameter :: M_earth = 5.9736d27      !Earth mass (g)
    real(selected_real_kind(6,8)),   parameter :: R_earth = 6.378136d8     !Earth radius (equatorial) (cm)
    real(selected_real_kind(10,13)), parameter :: AU = 1.4959787066d13     !Astronomical unit (cm)
    real(selected_real_kind(9,17)),  parameter :: ly = 9.460730472d17      !light year (cm)
    real(selected_real_kind(7,18)),  parameter :: pc = 3.0856776d18        !parsec (cm)
    real(selected_real_kind(11,4)),  parameter :: d_sid = 8.61640905309d4  !Sidereal day (s)
    real(selected_real_kind(4,4)),   parameter :: d_sol = 8.6400d4         !Solar day (s)
    real(selected_real_kind(8,7)),   parameter :: y_sid = 3.15581450d7     !Sidereal year (s)
    real(selected_real_kind(9,7)),   parameter :: y_trop = 3.155692519d7   !Tropical year (s)
    real(selected_real_kind(7,7)),   parameter :: y_jul = 3.1557600d7      !Julian year (s)
    real(selected_real_kind(7,7)),   parameter :: y_greg = 3.1556952d7     !Gregorian year (s)
    
END MODULE PhysicalConstants