
from numpy import exp

def sury(alb=0.101,\
         emi=0.86,\
         lam_s = 0.777,\
         Cv_s=1.25e6,\
         H=15.,\
         HtW = 1.5,\
         roof_f = 0.667,\
         depths = [0.,0.01,0.035,0.08,0.17,0.35,0.71,1.43,2.87,5.75,11.51],\
         lam_soil = 0.28,\
         Cv_soil = 1.35e6,\
         alb_snow = 0.70,\
         emi_snow = 0.997,\
         snow_f = 0.0,\
         ustar=0.25):

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

    if HtW > 2.0:
        print('Warning, specified canyon height-to-width-ratio falls outside of validity range. In order to extend the formulation, please contact the author.')
        
    # surface area index of a perfect canyon: parallel canyons and flat roofs
    SAI = (1. + 2. * HtW)*(1. - roof_f) + roof_f

    # surface-level bulk thermal properties
    Cv_bulk_s  = SAI * Cv_s
    lam_bulk_s = SAI * lam_s

    # vertical profiles for bulk heat conductivity and heat capacity
    Cv_bulk = []
    lam_bulk = []
    for d in depths:
        weight = min(d,H) 
        Cv_bulk.append( (1.-d/H) * Cv_bulk_s + d/H * Cv_soil )
        lam_bulk.append( (1.-d/H) * lam_bulk_s + d/H * lam_soil )

    # bulk radiative properties
    psi_canyon = exp(-0.6 * HtW)
    psi_bulk = roof_f + (1.-roof_f)*psi_canyon
    alb_bulk = ((1.-snow_f) * alb + snow_f * alb_snow) * psi_bulk
    emi_bulk = 1. - psi_bulk*(1. - ((1. - snow_f) * emi + snow_f * emi_snow))

    # surface-layer turbulence properties
    z_0 = 0.075 * H
    nu = 1.461e-5
    Re = ustar * z_0/nu
    kBm1 = 1.29 * Re**0.25 - 2.
    
    return {'alb_bulk':alb_bulk,\
            'emi_bulk':emi_bulk,\
            'lam_bulk':lam_bulk,\
            'Cv_bulk' :Cv_bulk,\
            'z_0'     :z_0,\
            'kBm1'    :kBm1\
           }






