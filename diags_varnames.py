"""
List of variables names to be used
==============================

"""

Var_list = { # paXXX
              'm01s16i222': 'pressure_at_sea_level',
              'm01s03i236': 'temperature',
              'm01s03i250': 'dew_point_temperature',
              'm01s03i248': 'fog_fraction',
              'm01s03i245': 'rh',
              'm01s03i237': 'q',
              'm01s00i409': 'sfc_pressure',
              'm01s02i207': 'sfc_downwelling_LW',
              'm01s01i235': 'sfc_downwelling_SW',
              'm01s02i201': 'sfc_net_LW',                         ## swath and full nest
              'm01s01i201': 'sfc_net_SW',                         ## swath and full nest
              'm01s00i024': 'sfc_temperature',
              'm01s01i207': 'toa_incoming_shortwave_flux',
              'm01s02i205': 'toa_outgoing_longwave_flux',
              'm01s01i208': 'toa_outgoing_shortwave_flux',
              'm01s03i225': 'u_10m',
              'm01s03i226': 'v_10m',
              'm01s01i504': 'seaice_albedo_cats',                         ## swath and full nest
              'm01s01i505': 'seaice_albedo_agg',                         ## swath and full nest
               # pbXXX
              'm01s03i304': 'turbulent_mixing_height_after_bl',
              'm01s03i360': 'h_decoupled_layer_base',
              'm01s03i361': 'h_sc_cloud_base',
              'm01s03i476': 'bl_type',
              'm01s09i216': 'cloud_area_fraction_assuming_random_overlap',
              'm01s09i217': 'cloud_area_fraction_assuming_maximum_random_overlap',
              'm01s09i221': 'wet_bulb_freezing_level_altitude',
              'm01s30i461': 'total_column_q',                                   # th 1-70
              'm01s03i241': 'water_evaporation_amount',
              'm01s02i392': 'LWP',
              'm01s02i391': 'IWP',
              'm01s16i222': 'air_pressure_at_sea_level',                        ## swath and full nest
              'm01s03i236': 'air_temperature_at_1.5m',                          ## swath and full nest
              'm01s03i025': 'bl_depth',
              'm01s03i250': 'dew_point_temperature_at_1.5m',
              'm01s03i248': 'fog_fraction',
              'm01s09i205': 'high_cloud',
              'm01s09i203': 'low_cloud',
              'm01s09i204': 'medium_cloud',
              'm01s03i245': 'rh_1.5m',                        ## swath and full nest
              'm01s03i237': 'q_1.5m',                        ## swath and full nest
              'm01s04i203': 'rainfall_flux',
              'm01s04i204': 'snowfall_flux',
              'm01s00i409': 'sfc_pressure',                             ## swath and full nest
              'm01s00i024': 'sfc_temperature',                              ## swath and full nest
              'm01s03i234': 'latent_heat_flux',
              'm01s03i217': 'sensible_heat_flux',
              'm01s03i247': 'visibility',
              'm01s03i225': 'u_10m',                             ## swath and full nest
              'm01s03i226': 'v_10m',                            ## swath and full nest
               # pcXXX  -- CLOUDNET
              'm01s04i118': 'radr_refl',                          # th 1-70
              'm01s00i266': 'cloud_fraction',                  # th 1-70 - pc
              'm01s00i267': 'liquid_cloud_fraction',
              'm01s00i268': 'ice_cloud_fraction',
              'm01s00i004': 'theta',
              'm01s00i408': 'pressure',                                     # th 1-70 - pc
              'm01s16i004': 'temperature',                                  # th 1-70 - pc
              'm01s00i012': 'qice',                # th 1-70 - pc
              'm01s00i254': 'qliq',       # th 1-70 - pc
              'm01s00i010': 'q',                                # th 1-70 - pc
              'm01s00i150': 'w',                              # th 1-70 - pc
              'm01s00i002': 'u',                                    # th 1-70 - pc
              'm01s00i003': 'v',                                   # th 1-70 - pc
              # pdXXX -- BOUNDARY LAYER
              'm01s03i362': 'entrainment_rate_SML',
              'm01s03i363': 'entrainment_rate_BL',
              'm01s03i464': 'obukhov_length',
              'm01s03i465': 'explicit_friction_velocity',
              'm01s04i298': 'diagnosed_turbulent_dissipation_rate',
              'm01s03i208': 'bulk_richardson_number',
              'm01s03i219': 'atmosphere_downward_eastward_stress',              # ro 1-70
              'm01s03i220': 'atmosphere_downward_northward_stress',             # ro 1-70
              'm01s03i473': 'tke',                         # ro 1-70
              'm01s00i004': 'theta',
              'm01s00i031': 'sea_ice_fraction',
              'm01s00i010': 'q',
              'm01s03i460': 'surface_downward_eastward_stress',
              'm01s03i461': 'surface_downward_northward_stress',
              'm01s00i026': 'surface_roughness_length',
              'm01s03i223': 'surface_upward_water_flux',
              'm01s03i501': 'mixing_length_for_momentum',
              'm01s03i471': 'BL_momentum_diffusion',
              'm01s03i469': 'vertical_buoyancy_gradient',
              'm01s03i135': 'production_rate_of_tke_by_shear',
              'm01s03i136': 'production_rate_of_tke_by_buoyancy',
              'm01s03i137': 'dissipation_rate_of_tke',
              # peXXX -- CASIM MICROPHYSICS
              'm01s00i075': 'qnliq',
              'm01s00i078': 'qnice',
              'm01s00i083': 'qnsolaeroliq',
              'm01s00i084': 'qnsolaerorain',
              'm01s00i088': 'qnsolaero',
              'm01s00i271': 'qicecrystals'
              }


def returnWantedVarnames():

    return Var_list.keys()


def findfieldName(stash):

    if stash in Var_list.keys():
        return Var_list[stash]
    else:
        return None
