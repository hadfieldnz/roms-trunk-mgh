!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Fennel et al. (2006) model:                          !
!                                                                      !
!   AttSW    Light attenuation due to sea water [1/m].                 !
!   AttChl   Light attenuation by Chlorophyll [1/(mg_Chl m2)].         !
!   BioIter  Maximum number of iterations to achieve convergence       !
!              of the nonlinear solution.                              !
!   Chl2C_m  Maximum chlorophyll to carbon ratio [mg_Chl/mg_C].        !
!   ChlMin   Chlorophill minimum threshold value [mg_Chl/m3].          !
!   CoagR    Coagulation rate: agregation rate of SDeN + Phyt ==> LDeN !
!              [1/day].                                                !
!   D_p5NH4  Half-saturation radiation for nitrification inhibition    !
!              [Watts/m2].                                             !
!   I_thNH4  Radiation threshold for nitrification inhibition          !
!              [Watts/m2].                                             !
!   K_NH4    Inverse half-saturation for Phytoplankton NH4 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_NO3    Inverse half-saturation for Phytoplankton NO3 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_Phy    Zooplankton half-saturation, squared constant for         !
!              ingestion [mmol_N/m3]^2.                                !
!   LDeRR    Large Detrital re-mineralization rate [1/day].            !
!   NitriR   Nitrification rate: oxidation of NH4 to NO3 [1/day].      !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   PhyCN    Phytoplankton Carbon:Nitrogen ratio [mol_C/mol_N].        !
!   PhyIP    Phytoplankton NH4 inhibition parameter [1/(mmol_N)].      !
!   PhyIS    Phytoplankton, initial slope of the P-I curve             !
!              [mg_C/(mg_Chl W m-2 day)].                              !
!   ZooMin   Phytoplankton minimum threshold value [mmol_N/m3].        !
!   PhyMR    Phytoplankton mortality rate [1/day] to small detritus.   !
!   SDeAR    Small detritus aggregation rate into Large detritus       !
!              [1/day].                                                !
!   SDeBR    Small Detrital breakdown to NH4 rate [1/day].             !
!   SDeRR    Large Detrital re-mineralization rate [1/day].            !
!   Vp0      Eppley temperature-limited and light-limited growth       !
!              tuning parameter [nondimensional].                      !
!   wLDet    Vertical sinking velocities for Large Detritus            !
!              fraction [m/day].                                       !
!   wPhy     Vertical sinking velocity for Phytoplankton               !
!              fraction [m/day].                                       !
!   wSDet    Vertical sinking velocities for Small Detritus            !
!              fraction [m/day].                                       !
!   ZooAE_N  Zooplankton nitrogen assimilation efficiency fraction     !
!              [nondimensional].                                       !
!   ZooBM    Zooplankton basal metabolism [1/day].                     !
!   ZooCN    Zooplankton Carbon:Nitrogen ratio [mol_C/mol_N].          !
!   ZooER    Zooplankton specific excretion rate [1/day].              !
!   ZooGR    Zooplankton maximum growth rate [1/day].                  !
!   ZooMin   Zooplankton minimum threshold value [mmol_N/m3].          !
!   ZooMR    Zooplankton mortality to Detritus [1/day].                !
!   pCO2air  CO2 partial pressure in the air [ppmv].                   !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none
!
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: AttChl        ! 1/(mg_Chl m2)
        real(r8), dimension(Ngrids) :: Chl2C_m       ! mg_Chl/mg_C
        real(r8), dimension(Ngrids) :: ChlMin        ! mg_Chl/m3
        real(r8), dimension(Ngrids) :: CoagR         ! 1/day
        real(r8), dimension(Ngrids) :: D_p5NH4       ! Watts/m2
        real(r8), dimension(Ngrids) :: I_thNH4       ! Watts/m2
        real(r8), dimension(Ngrids) :: K_NH4         ! m3/mmol_N
        real(r8), dimension(Ngrids) :: K_NO3         ! m3/mmol_N
        real(r8), dimension(Ngrids) :: K_Phy         ! (mmol_N/m3)^2
        real(r8), dimension(Ngrids) :: LDeRRN        ! 1/day
        real(r8), dimension(Ngrids) :: LDeRRC        ! 1/day
        real(r8), dimension(Ngrids) :: NitriR        ! 1/day
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
        real(r8), dimension(Ngrids) :: PhyCN         ! mol_C/mol_N
        real(r8), dimension(Ngrids) :: PhyIP         ! 1/mmol_N
        real(r8), dimension(Ngrids) :: PhyIS         ! 1/(Watts m-2 day)
        real(r8), dimension(Ngrids) :: PhyMin        ! mmol_N/m3
        real(r8), dimension(Ngrids) :: PhyMR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeAR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeBR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeRRN        ! 1/day
        real(r8), dimension(Ngrids) :: SDeRRC        ! 1/day
        real(r8), dimension(Ngrids) :: Vp0           ! nondimensional
        real(r8), dimension(Ngrids) :: wLDet         ! m/day
        real(r8), dimension(Ngrids) :: wPhy          ! m/day
        real(r8), dimension(Ngrids) :: wSDet         ! m/day
        real(r8), dimension(Ngrids) :: ZooAE_N       ! nondimensional
        real(r8), dimension(Ngrids) :: ZooBM         ! 1/day
        real(r8), dimension(Ngrids) :: ZooCN         ! mol_C/mol_N
        real(r8), dimension(Ngrids) :: ZooER         ! 1/day
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMin        ! mmol_N/m3
        real(r8), dimension(Ngrids) :: ZooMR         ! 1/day
        real(r8), dimension(Ngrids) :: pCO2air       ! ppmv
