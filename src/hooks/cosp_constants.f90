MODULE cosp_math_constants
  USE cosp_kinds, only: wp
  IMPLICIT NONE
  REAL(wp), PARAMETER :: pi =  3.14159265358979323846264338327950288419717_wp

END MODULE cosp_math_constants

MODULE cosp_phys_constants
  USE cosp_kinds, only: wp
  IMPLICIT NONE
  
  REAL(wp), PARAMETER :: &
       tmelt  = 273.15_wp,      & ! Melting temperature of ice/snow [K]
       rhoice = 917._wp,        & ! Density of ice [kg/m3]
       rholiq = 1000._wp          ! Density of liquid water [kg/m3]

  ! Molecular weights
  REAL(wp), PARAMETER :: &
       amw   = 18.01534_wp,     & ! Water   [g/mol]
       amd   = 28.9644_wp,      & ! Dry air [g/mol]
       amO3  = 47.9983_wp,      & ! Ozone   [g/mol]
       amCO2 = 44.0096_wp,      & ! CO2     [g/mol]
       amCH4 = 16.0426_wp,      & ! Methane [g/mol]
       amN2O = 44.0129_wp,      & ! N2O     [g/mol]
       amCO  = 28.0102_wp         ! CO      [g/mol]

  ! WMO/SI value
  REAL(wp), PARAMETER :: &
       avo   = 6.023E23_wp,     & ! Avogadro constant used by ISCCP simulator [1/mol]
       grav  = 9.806650_wp        ! Av. gravitational acceleration used by ISCCP simulator [m/s2]

  ! Thermodynamic constants for the dry and moist atmosphere
  REAL(wp), PARAMETER :: &
       rd  = 287.04_wp,         & ! Gas constant for dry air [J/K/Kg]
       cpd = 1004.64_wp,        & ! Specific heat at constant pressure for dry air [J/K/Kg]
       rv  = 461.51_wp,         & ! Gas constant for water vapor [J/K/Kg]
       cpv = 1869.46_wp,        & ! Specific heat at constant pressure for water vapor [J/K/Kg]
       km  = 1.38e-23_wp          ! Boltzmann constant [J/K]

END MODULE cosp_phys_constants


