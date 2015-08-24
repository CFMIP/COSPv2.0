! $Revision: 23 $, $Date: 2011-03-31 07:41:37 -0600 (Thu, 31 Mar 2011) $
! $URL: https://cfmip-obs-sim.googlecode.com/svn/devel/branches/dustinswales/MODIS_simulator/test_modis_simulator.F90 $
program test_modis_simulator
  USE COSP_KINDS, ONLY: wp
  use mod_modis_sim
  implicit none
  ! 
  ! Tests cases for MODIS L2 (pixel) and L3 (grid-scale) simulators
  !
  
  !
  ! Inputs
  !
  integer, parameter :: nLayers = 6, nSubCols = 10
  real(wp), dimension(1, nLayers)    :: temp, pressureLayers
  real(wp), dimension(1, nLayers+1)  :: pressureLevels
  real(wp), dimension(1, nSubCols, nLayers) :: &
                                    opticalThickness, cloudWater, cloudIce, waterSize, iceSize, &
                                    liquid_opticalThickness, ice_opticalThickness
  !
  ! Pixel-scale retreivals
  ! 
  integer, dimension(1, nSubCols) :: retrievedPhase
  real(wp),    dimension(1, nSubCols) :: isccpTau, isccpCloudTopPressure, &
                                  retrievedCloudTopPressure, retrievedTau, retrievedSize
  integer :: i, k
  character(len = 6) :: phase
  
  !
  ! Grid mean properties 
  !
  real(wp), dimension(1) :: &
    Cloud_Fraction_Total_Mean,       Cloud_Fraction_Water_Mean,       Cloud_Fraction_Ice_Mean,       &
    Optical_Thickness_Total_Mean,    Optical_Thickness_Water_Mean,    Optical_Thickness_Ice_Mean,    &
    Optical_Thickness_Total_LogMean, Optical_Thickness_Water_LogMean, Optical_Thickness_Ice_LogMean, &
                                     Cloud_Particle_Size_Water_Mean,  Cloud_Particle_Size_Ice_Mean,  &
    Cloud_Top_Pressure_Total_Mean,                                                                   &
                                     Liquid_Water_Path_Mean,          Ice_Water_Path_Mean
  real(wp), dimension(1, numTauHistogramBins, numPressureHistogramBins) :: Optical_Thickness_vs_Cloud_Top_Pressure

  ! ---------------------------------------------------------------------------------------------
  pressureLevels(1, :) = (/   0., 200., 400., 600., 800., 1000., 1200. /) 
  pressureLayers(1, :) = (pressureLevels(1, 2:) + pressureLevels(1, 1:nLayers)) / 2
  temp(1, :) =           (/ 200., 220., 240., 260., 280., 300.0 /) 
  isccpCloudTopPressure(:, :) = 999. 
  

  !
  ! Tests:  these 9 from Steve Platnick
  !
  ice_opticalThickness = 0.; liquid_opticalThickness = 0.
  iceSize              = 0.; waterSize = 0. 
  
  ice_opticalThickness   (1, 1, 3) = 0.2
  iceSize                (1, 1, 3) = 30
  
  ice_opticalThickness   (1, 2, 3) = 1.
  iceSize                (1, 2, 3) = 30
  
  liquid_opticalThickness(1, 3, 5) = 0.2
  waterSize              (1, 3, 5) = 10 

  liquid_opticalThickness(1, 4, 5) = 1.
  waterSize              (1, 4, 5) = 10 
  
  ice_opticalThickness   (1, 5, 3) = 0.2
  iceSize                (1, 5, 3) = 30
  liquid_opticalThickness(1, 5, 5) = 5.
  waterSize              (1, 5, 5) = 10 

  ice_opticalThickness   (1, 6, 3) = 1. 
  iceSize                (1, 6, 3) = 30
  liquid_opticalThickness(1, 6, 5) = 5.
  waterSize              (1, 6, 5) = 10
  
  ice_opticalThickness   (1, 7, 3) = 10. 
  iceSize                (1, 7, 3) = 30
  liquid_opticalThickness(1, 7, 5) = 5.
  waterSize              (1, 7, 5) = 10
  
  ice_opticalThickness   (1, 8, 2) = 1. 
  iceSize                (1, 8, 2) = 15.
  ice_opticalThickness   (1, 8, 3) = 1. 
  iceSize                (1, 8, 3) = 30.
  ice_opticalThickness   (1, 8, 4) = 1. 
  iceSize                (1, 8, 4) = 60.
  
  liquid_opticalThickness(1, 9, 4) = 1. 
  waterSize              (1, 9, 4) = 15
  liquid_opticalThickness(1, 9, 5) = 3. 
  waterSize              (1, 9, 5) = 10.
  liquid_opticalThickness(1, 9, 6) = 1. 
  waterSize              (1, 9, 6) = 50.
  
  call modis_l2_simulator(temp(1,:), pressureLayers(1,:), pressureLevels(1,:),                       &
                         liquid_opticalThickness(1,:, :), ice_opticalThickness(1,:, :), &
                         waterSize(1,:, :), iceSize(1,:, :), &
                         isccpTau(1, :), isccpCloudTopPressure(1, :),                            &
                         retrievedPhase(1, :), retrievedCloudTopPressure(1, :), retrievedTau(1, :), retrievedSize(1, :))
  do i = 1, nSubcols
    print * 
    print *, "Profile ", i
    if(any(ice_opticalThickness(1, i, :) > 0.) .or. any(liquid_opticalThickness(1, i, :) > 0.)) & 
      print *, "center p  ice: tau   size    water: tau    size" 
    do k = 1, nLayers
      if(ice_opticalThickness(1, i, k) > 0. .or. liquid_opticalThickness(1, i, k) > 0.) &
        write(*, '(3x, f5.0, 7x, 2(f4.1, 3x, f4.1, 10x))'), &
          pressureLayers(1, k), ice_opticalThickness(1, i, k), iceSize(1, i, k), liquid_opticalThickness(1, i, k), waterSize(1, i, k)
    end do
    select case(retrievedPhase(1, i)) 
      case(phaseIsNone) 
        phase = "None" 
      case(phaseIsIce) 
        phase = "Ice" 
      case(phaseIsLiquid)
        phase = "Liquid" 
      case(phaseIsUndetermined)
        phase = "Undet."
      case default
        phase = "???" 
    end select
    print *, "Retrievals" 
    write(*, '("  Phase: ", a6, 2x, "CTP: ", f5.0, ", tau = ", f4.1, ", size = ", f4.1)') &
          phase, retrievedCloudTopPressure(1, i), retrievedTau(1, i), retrievedSize(1, i)
  end do 
  
  call modis_L3_simulator(retrievedPhase, retrievedCloudTopPressure, retrievedTau, retrievedSize,            &
       Cloud_Fraction_Total_Mean,       Cloud_Fraction_Water_Mean,       Cloud_Fraction_Ice_Mean,       &
       Optical_Thickness_Total_Mean,    Optical_Thickness_Water_Mean,    Optical_Thickness_Ice_Mean,    &
       Optical_Thickness_Total_LogMean, Optical_Thickness_Water_LogMean, Optical_Thickness_Ice_LogMean, &
                                        Cloud_Particle_Size_Water_Mean,  Cloud_Particle_Size_Ice_Mean,  &
       Cloud_Top_Pressure_Total_Mean,                                                                   &
                                        Liquid_Water_Path_Mean,          Ice_Water_Path_Mean,           &    
       Optical_Thickness_vs_Cloud_Top_Pressure)

  print *; print *; print *, "Grid means" 
  do i = 1, nSubcols
    write(*, '("  Phase: ", i2, 2x, "CTP: ", f5.0, ", tau = ", f4.1, ", size = ", f4.1)') &
          retrievedPhase(1, i), retrievedCloudTopPressure(1, i), retrievedTau(1, i), retrievedSize(1, i)
  end do
  print *, "         Total    Liquid     Ice" 
  write(*, '("Fraction: ", 3(f4.3, 6x))') Cloud_Fraction_Total_Mean,       Cloud_Fraction_Water_Mean,       Cloud_Fraction_Ice_Mean
  write(*, '("Tau:      ", 3(f4.1, 6x))')  Optical_Thickness_Total_Mean,    Optical_Thickness_Water_Mean,    Optical_Thickness_Ice_Mean
  write(*, '("TauLog:   ", 3(f4.1, 6x))')  Optical_Thickness_Total_LogMean, Optical_Thickness_Water_LogMean, Optical_Thickness_Ice_LogMean
  write(*, '("Size:     ", 10x, 2(f4.1, 6x))')  Cloud_Particle_Size_Water_Mean,  Cloud_Particle_Size_Ice_Mean
  write(*, '("WaterPath:",  9x, 2(f5.1, 5x))')  Liquid_Water_Path_Mean,          Ice_Water_Path_Mean
  write(*, '("CTP:      ", f5.1)') Cloud_Top_Pressure_Total_Mean
  print *; print *, "Histogram" 
  write(*, '(8x, 7(f4.1, 3x))')  tauHistogramBoundaries(:)
  do k = 1, numPressureHistogramBins
     write(*, '(f5.0, 3x, 7(f4.2, 3x))') pressureHistogramBoundaries(k), Optical_Thickness_vs_Cloud_Top_Pressure(1, :, k)
  end do 
  print *, "Total fraction from histogram:", sum(Optical_Thickness_vs_Cloud_Top_Pressure(1, :, :))
end program   test_modis_simulator

