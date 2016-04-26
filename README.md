# sury

    """
    Purpose: Derive bulk parameters from urban-canopy parameters according
             the Semi-Empirical URban canopY parametrization (SURY).  
    Reference:
             Wouters, H., Demuzere, M., Blahak, U., Fortuniak, K., Maiheu., B., 
             Camps, J., Tielemans, and N. P. M. van Lipzig, 2016.  Efficient
             urban canopy parametrization for atmospheric modelling:
             description and application with the COSMO-CLM model (version
             5.0_clm6) for a Belgian Summer, Geosci. Model Dev.  Discuss., in
             review, 2016.
    Version: 1.0
    Author: Hendrik Wouters <hendrik.wouters@kuleuven.be> 
                                                                              
    License: GNU GENERAL PUBLIC LICENSE Version 3 (license is provided along)
    input (default values are listed above):
        alb:      substrate albedo      
        emi:      substrate emissivity
        lam_s:    substrate heat conductivity
        Cv_s:     substrate heat capacity
        H:        building height
        HtW:      canyon height-to-width ratio; validity range: [0..2]
        roof_f:   roof fraction
        depths:   depths of the ground layers in the vertical column of the bulk land-surface module
        lam_soil: heat conductivity of the soil below
        Cv_soil:  heat capacity of the soil below
        alb_snow: albedo of snow
        emi_snow: emissivity of snow 
        snow_f:   snow fraction
        ustar:    friction velocity
    output:
        alb_bulk: bulk albedo
        emi_bulk: bulk emissivity
        lam_bulk: vertical profile bulk heat conductivity
        Cv_bulk: vertical profile substrate heat capacity
        z0 : aerodynamic roughness length
        kBm1: kB^{-1} = log(z0/z0H), where z0H is the thermal roughness length
    Usage:
        import sury
        
        # default urban-canopy parameters:
        sury.sury()
        # user-specified urban-canopy parameters, eg.,:
        sury.sury(H = 20., HtW = 10., roof_f = 0.5)
                                                                              
    """
