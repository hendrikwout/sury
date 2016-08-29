from numpy import exp,amin,zeros_like

def sury(alb=0.101,\
         emi=0.86,\
         lam_s = 0.767,\
         Cv_s=1.25e6,\
         h=15.,\
         htw = 1.5,\
         roof_f = 0.667,\
         z = 0.,\
         lam_soil = 0.28,\
         Cv_soil = 1.35e6,\
         alb_snow = 0.70,\
         emi_snow = 0.997,\
         snow_f = 0.0,\
         ustar=0.25,\
         **kwargs):

    """
    Purpose: SURY converts urban canopy parameters – containing the three-dimensional 
	     information such as from WUDAPT – into bulk parameters. It can be used 
	     with existing bulk urban land-surface schemes for providing canopy-
	     dependent urban physics to atmospheric modelling with a low computation cost.
             
    Reference:
             Wouters, H., Demuzere, M., Blahak, U., Fortuniak, K., Maiheu., B., 
             Camps, J., Tielemans, and N. P. M. van Lipzig, 2016.  The efficient
             urban-canopy dependency parametrization SURY (v1.0) for atmospheric modelling:
             description and application with the COSMO-CLM model (v5.0_clm6) for a 
             Belgian Summer, Geosci. Model Dev., 2016.

    Version: 1.0

    Author: Hendrik Wouters <hendrik.wouters@kuleuven.be> 
                                                                              
    License: GNU GENERAL PUBLIC LICENSE Version 3 (license is provided along)


    Input (default values are listed above):
        alb:      albedo [ 1 ] 
        emi:      emissivity [ 1 ]
        lam_s:    surface heat conductivity [W m-1 K-1] 
        Cv_s:     surface heat capacity [J m-3 K-1]
        H:        building height [m]
        htw:      canyon height-to-width ratio [1],  validity range: 0..2
        roof_f:   roof fraction [ 1 ]
        z:        depths of the ground layers in the vertical column of 
                  the bulk land-surface module. Please note that this can be specified
                  as a list or array, eg., [0.01,0.035,0.08,0.17,0.35,0.71,1.43,2.87,5.75,11.51]
        lam_soil: heat conductivity of the soil below [W m-1 K-1] 
        Cv_soil:  heat capacity of the soil below [J m-3 K-1]
        alb_snow: albedo of snow [ 1 ]
        emi_snow: emissivity of snow [ 1 ]
        snow_f:   snow fraction [ 1 ]
        ustar:    friction velocity [ m s-1 ]

        --- When provided, Cv_s is obtained from the weighted average of the
            values $C_v,i} and $\lambda_i$) for wall, roof and road according to their
            surface fractions in the urban canopy (see Eqns. 10 of Wouters et al., 2016)
        Cv_roof: heat capacity of the roofs [J m-3 K-1]
        Cv_wall: heat capacity of the walls [J m-3 K-1]
        Cv_road: heat capacity of the roads [J m-3 K-1]

        --- When provided, lam_s is obtained from the weighted average of the
            values $C_v,i} and $\lambda_i$) for wall, roof and road according to their
            surface fractions in the urban canopy, see Eq.  11 of Wouters et al., 2016)
        lam_roof: heat conductivity of the roofs [W m-1 K-1]
        lam_wall: heat conductivity of the walls [W m-1 K-1]
        lam_road: heat conductivity of the roads [W m-1 K-1]

        --- When provided, alb is obtained from the
            values  for wall, roof and road according to  Eq. 16 of Wouters et al. (2016)
        alb_roof: heat conductivity of the roofs [W m-1 K-1]
        alb_wall: heat conductivity of the walls [W m-1 K-1]
        alb_road: heat conductivity of the roads [W m-1 K-1]

        --- When provided, emi is obtained from the
            values  for wall, roof and road analogous to  Eq. 16 of Wouters et al. (2016)
        emi_roof: heat conductivity of the roofs [W m-1 K-1]
        emi_wall: heat conductivity of the walls [W m-1 K-1]
        emi_road: heat conductivity of the roads [W m-1 K-1]



    Output:
        alb_bulk: bulk albedo [ 1 ]
        emi_bulk: bulk emissivity [ 1 ]
        lam_bulk: vertical profile bulk heat conductivity [W m-1 K-1]
        Cv_bulk: vertical profile substrate heat capacity [J m-3 K-1]
        z0 : aerodynamic roughness length [ m ]
        kBm1: kB^{-1} = log(z0/z0H), where z0H is the thermal roughness length [ 1 ]


    Usage:
        import sury
        
        # Use the default urban-canopy parameters:
        sury.sury()
        # Some user-specified urban-canopy parameters (As always in Python, please explicitly 
        # specify the parameter names when deviating from the default argument order):
        sury.sury(h = 20., htw = 10., roof_f = 0.5)

    """

    if htw > 2.0:
        print('Warning, specified canyon height-to-width-ratio falls outside of validity range. In order to extend the formulation, please contact the author.')
        
    # surface area index of a perfect canyon: parallel canyons and flat roofs
    SAI = (1. + 2. * htw)*(1. - roof_f) + roof_f


    if ('Cv_roof' in kwargs.keys()) or ('Cv_wall' in kwargs.keys()) or ('Cv_road' in kwargs.keys()):
        Cv_s =(kwargs['Cv_roof'] * roof_f + (kwargs['Cv_wall']*2.*htw+kwargs['Cv_road'])*(1.-roof_f)  )/(roof_f+ (2.*htw+1.)*(1.-roof_f))
        print('Warning, Cv_s is calculated from the values for heat capacity from roads, walls and roofs.'+str(Cv_s))
    if ('lam_roof' in kwargs.keys()) or ('lam_wall' in kwargs.keys()) or ('lam_road' in kwargs.keys()):
        lam_s =(kwargs['lam_roof'] * roof_f + (kwargs['lam_wall']*2.*htw+kwargs['lam_road'])*(1.-roof_f)  )/(roof_f+ (2.*htw+1.)*(1.-roof_f))
        print('Warning, lam_s is calculated from the values for heat capacity from roads, walls and roofs: '+str(lam_s))


    # surface-level bulk thermal properties
    Cv_bulk_s  = SAI * Cv_s
    lam_bulk_s = SAI * lam_s

    # vertical profiles for bulk heat conductivity and heat capacity

    weight = amin([z,zeros_like(z)+h],axis=0)
    Cv_bulk=  (1.-weight/h) * Cv_bulk_s + weight/h * Cv_soil 
    lam_bulk = (1.-weight/h) * lam_bulk_s + weight/h * lam_soil 

    # bulk radiative properties
    psi_canyon = exp(-0.6 * htw)
    psi_bulk = roof_f + (1.-roof_f)*psi_canyon

    if ('alb_roof' in kwargs.keys()) or ('alb_wall' in kwargs.keys()) or ('alb_road' in kwargs.keys()):
        alb_roof_snow   = kwargs['alb_roof']  *(1. - snow_f) +  alb_snow* snow_f
        alb_road_snow = kwargs['alb_road']*(1. - snow_f) +  alb_snow* snow_f
        alb_wall_snow   = kwargs['alb_wall']  *(1. - snow_f) +  alb_snow* snow_f
        alb_bulk = (alb_road_snow  + 2. *htw * alb_wall_snow )/(1.+2. * htw) * psi_canyon * (1. - roof_f) + alb_roof_snow * roof_f 

    else:
        alb_bulk = ((1.-snow_f) * alb + snow_f * alb_snow) * psi_bulk


    if ('emi_roof' in kwargs.keys()) or ('emi_wall' in kwargs.keys()) or ('emi_road' in kwargs.keys()):
        albth_roof_snow   = (1.-kwargs['emi_roof']) * (1. - snow_f) + (1. - alb_snow) * snow_f
        albth_road_snow   = (1.-kwargs['emi_road']) * (1. - snow_f) + (1. - emi_snow)* snow_f
        albth_wall_snow   = (1.-kwargs['emi_wall'])  *(1. - snow_f) + (1. - emi_snow) * snow_f
        emi_bulk = 1. - ((albth_road_snow  + 2. *htw * albth_wall_snow )/(1.+2. * htw) * psi_canyon * (1. - roof_f) + albth_roof_snow * roof_f)

    else:
         emi_bulk = 1. - psi_bulk*(1. - ((1. - snow_f) * emi + snow_f * emi_snow))

    # surface-layer turbulence properties
    z_0 = 0.075 * h
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

