# Other libraries
import numpy as np
import lmfit as lsq

# Models
from .models import Molecule
from .models import Data_Sample_Lipid_Augment

# Other imports
from .symprobabilities import *
from scipy import signal as sig

# Asymmetrical model
# The math only - used by LMFit to run each iteration
def asym_model(
    q,          # independant
    Vci,        # chain_volume
    Vhi,        # headgroup_volume
    Vti,        # terminal_methyl_volume
    Ali,        # area_per_lipid
    Dhi,        # headgroup_thickness
    Vco,        # chain_volume
    Vho,        # headgroup_volume
    Vto,        # terminal_methyl_volume
    Alo,        # area_per_lipid
    Dho,        # headgroup_thickness
    Vw,         # water_volume
    sig,        # smearing factor
    bci,        # chain_b
    bhi,        # headgroup_b
    bti,        # terminal_methyl_b
    bco,        # chain_b
    bho,        # headgroup_b
    bto,        # terminal_methyl_b
    bw,         # water_b
    scale,      # scale
    bg          # bg
):
    return (
        q**(-2) * scale * (
            (
                ( 2 * ( np.exp( - ((q * sig)**2) / 2)) ) 
                * ( 
                    (
                        (
                            (
                                ( ( bhi - Vhi * (bw/Vw) ) * ( np.cos(-q*Dhi-(q*Vci/Ali))-np.cos(q*Vci/Ali) ) )
                                / (q*Ali*Dhi) 
                            ) + (
                                ( (bho-Vho*bw/Vw) * (np.cos(q*Vco/Alo) - np.cos(q*Dho+(q*Vco/Alo))) ) 
                                / (q*Alo*Dho)
                            ) + (
                                ( ((bci-2*bti) / (Vci-2*Vti)-bw/Vw) * (np.cos(q*Vci/Ali) - np.cos(2*q*Vti/Ali)) ) 
                                /q
                            ) + (
                                ( ((bco-2*bto) / (Vco-2*Vto)-bw/Vw) * (np.cos(2*q*Vto/Alo) - np.cos(q*Vco/Alo)) )
                                /q
                            ) + (
                                ( ((bti/Vti) - bw/Vw) * (np.cos(2*q*Vti/Ali) - 1) )
                                /q
                            ) + (
                                ( ((bto/Vto) - bw/Vw) * (1 - np.cos(2*q*Vto/Alo)) )
                                /q
                            )
                        )**2
                    ) + (
                        (
                            (
                                ( (bhi-Vhi*(bw/Vw)) * (-np.sin(-q*Dhi-(q*Vci/Ali)) - np.sin(q*Vci/Ali)) )
                                /(q*Ali*Dhi)
                            ) + (
                                ( (bho-Vho*(bw/Vw)) * (-np.sin(q*Vco/Alo) + np.sin(q*Dho+(q*Vco/Alo))) )
                                /(q*Alo*Dho)
                            ) + (
                                ( (( (bci-2*bti) / (Vci-2*Vti) ) - (bw/Vw)) * (np.sin(q*Vci/Ali) - np.sin(2*q*Vti/Ali)) )
                                /q
                            ) + (
                                ( (( (bco-2*bto) / (Vco-2*Vto) ) - (bw/Vw)) * (np.sin(q*Vco/Alo) - np.sin(2*q*Vto/Alo)) )
                                /q
                            ) + (
                                ( ((bti/Vti) - bw/Vw) * (np.sin(2*q*Vti/Ali)) )
                                /q
                            ) + (
                                ( ((bto/Vto) - bw/Vw) * (np.sin(2*q*Vto/Alo)) )
                                /q
                            )
                        )
                    )**2
                )**(1/2)
            )**2
        ) + bg
    )

# Asymmetrical model and separated form factor
# The same math, but including the separated form factor equation
def asym_model_separated(
    q,          # independant
    Vci,        # chain_volume
    Vhi,        # headgroup_volume
    Vti,        # terminal_methyl_volume
    Ali,        # area_per_lipid
    Dhi,        # headgroup_thickness
    Vco,        # chain_volume
    Vho,        # headgroup_volume
    Vto,        # terminal_methyl_volume
    Alo,        # area_per_lipid
    Dho,        # headgroup_thickness
    Vw,         # water_volume
    sig,        # smearing factor
    r,          # Average vesicle radius
    rs,         # Relative size polydispersity
    bci,        # chain_b
    bhi,        # headgroup_b
    bti,        # terminal_methyl_b
    bco,        # chain_b
    bho,        # headgroup_b
    bto,        # terminal_methyl_b
    bw,         # water_b
    scale,      # scale
    bg          # bg
):

    return (

        q**(-2) * scale * (
                (4 * 10**-8) * (
                    (
                        (8*(np.pi**2)*(r**2)*(rs**4)) * (1 + rs**-2) * (2 + rs**-2)
                    ) * (
                        1 - (
                            (
                                ((1 + 4*(q**2)*(r**2)*(rs**4))**(-1 / (2*rs**2))) * (np.cos((2 + rs**-2) * np.arctan(2 * q * r * (rs**2))))
                            ) / (1 + 4*(q**2)*(r**2)*(rs**4))
                        )
                    )
            ) * (
                ( 2 * ( np.exp( - ((q * sig)**2) / 2)) ) 
                * ( 
                    (
                        (
                            (
                                ( ( bhi - Vhi * (bw/Vw) ) * ( np.cos(-q*Dhi-(q*Vci/Ali))-np.cos(q*Vci/Ali) ) )
                                / (q*Ali*Dhi) 
                            ) + (
                                ( (bho-Vho*bw/Vw) * (np.cos(q*Vco/Alo) - np.cos(q*Dho+(q*Vco/Alo))) ) 
                                / (q*Alo*Dho)
                            ) + (
                                ( ((bci-2*bti) / (Vci-2*Vti)-bw/Vw) * (np.cos(q*Vci/Ali) - np.cos(2*q*Vti/Ali)) ) 
                                /q
                            ) + (
                                ( ((bco-2*bto) / (Vco-2*Vto)-bw/Vw) * (np.cos(2*q*Vto/Alo) - np.cos(q*Vco/Alo)) )
                                /q
                            ) + (
                                ( ((bti/Vti) - bw/Vw) * (np.cos(2*q*Vti/Ali) - 1) )
                                /q
                            ) + (
                                ( ((bto/Vto) - bw/Vw) * (1 - np.cos(2*q*Vto/Alo)) )
                                /q
                            )
                        )**2
                    ) + (
                        (
                            (
                                ( (bhi-Vhi*(bw/Vw)) * (-np.sin(-q*Dhi-(q*Vci/Ali)) - np.sin(q*Vci/Ali)) )
                                /(q*Ali*Dhi)
                            ) + (
                                ( (bho-Vho*(bw/Vw)) * (-np.sin(q*Vco/Alo) + np.sin(q*Dho+(q*Vco/Alo))) )
                                /(q*Alo*Dho)
                            ) + (
                                ( (( (bci-2*bti) / (Vci-2*Vti) ) - (bw/Vw)) * (np.sin(q*Vci/Ali) - np.sin(2*q*Vti/Ali)) )
                                /q
                            ) + (
                                ( (( (bco-2*bto) / (Vco-2*Vto) ) - (bw/Vw)) * (np.sin(q*Vco/Alo) - np.sin(2*q*Vto/Alo)) )
                                /q
                            ) + (
                                ( ((bti/Vti) - bw/Vw) * (np.sin(2*q*Vti/Ali)) )
                                /q
                            ) + (
                                ( ((bto/Vto) - bw/Vw) * (np.sin(2*q*Vto/Alo)) )
                                /q
                            )
                        )
                    )**2
                )**(1/2)
            )**2
        ) + bg
    )

# Calculate result from model for an individual dataset
# Connetc the model, the dataset and the parameters
def calc_asym_model(fit_parameters, q, data, sff):
    # Convert array
    q_array = np.array(q)

    ### Unpack parameters
    ## Inner
    Vci = fit_parameters['in_chain_volume'].value
    Vhi = fit_parameters['in_headgroup_volume'].value
    Vti = fit_parameters['in_terminal_methyl_volume'].value
    # Unknown
    Ali = fit_parameters['in_area_per_lipid'].value
    Dhi = fit_parameters['in_headgroup_thickness'].value
    ## Outer
    Vco = fit_parameters['out_chain_volume'].value
    Vho = fit_parameters['out_headgroup_volume'].value
    Vto = fit_parameters['out_terminal_methyl_volume'].value
    # Unknown
    Alo = fit_parameters['out_area_per_lipid'].value
    Dho = fit_parameters['out_headgroup_thickness'].value
    ## Shared
    # Smearing factor
    sig = fit_parameters['sigma'].value
    # Separated form factor
    r = fit_parameters['average_vesicle_radius'].value
    rs = fit_parameters['relative_size'].value
    ### Per dataset
    ## Inner
    bci = fit_parameters['in_chain_b_%i' % data.id].value
    bhi = fit_parameters['in_headgroup_b_%i' % data.id].value
    bti = fit_parameters['in_terminal_methyl_b_%i' % data.id].value
    ## Inner
    bco = fit_parameters['out_chain_b_%i' % data.id].value
    bho = fit_parameters['out_headgroup_b_%i' % data.id].value
    bto = fit_parameters['out_terminal_methyl_b_%i' % data.id].value
    ## Shared
    bw = fit_parameters['water_b_%i' % data.id].value
    Vw = fit_parameters['combined_water_volume_%i' % data.id].value
    # Tweaks
    scale = fit_parameters['scale_%i' % data.id].value
    bg = fit_parameters['background_%i' % data.id].value

    # Return the calculated model
    if sff:
        calc_result = asym_model_separated(q_array, Vci, Vhi, Vti, Ali, Dhi, Vco, Vho, Vto, Alo, Dho, Vw, sig, r, rs, bci, bhi, bti, bco, bho, bto, bw, scale, bg)
    else:
        calc_result = asym_model(q_array, Vci, Vhi, Vti, Ali, Dhi, Vco, Vho, Vto, Alo, Dho, Vw, sig, bci, bhi, bti, bco, bho, bto, bw, scale, bg)

    return calc_result


# Objective function
# Create a residual for each dataset, then flatten into a single array for minimize()
def asymmetrical_objective_function(fit_parameters, x, datas, sff):
    # Delcare
    current_residual = []
    combined_residuals = []
    scaled_water = []

    # Check if the water prob is negative - if it is, impose a penalty

    ### Unpack parameters
    ## Inner
    Vci = fit_parameters['in_chain_volume'].value
    Vhi = fit_parameters['in_headgroup_volume'].value
    Vti = fit_parameters['in_terminal_methyl_volume'].value
    # Unknown
    Ali = fit_parameters['in_area_per_lipid'].value
    Dhi = fit_parameters['in_headgroup_thickness'].value
    ## Outer
    Vco = fit_parameters['out_chain_volume'].value
    Vho = fit_parameters['out_headgroup_volume'].value
    Vto = fit_parameters['out_terminal_methyl_volume'].value
    # Unknown
    Alo = fit_parameters['out_area_per_lipid'].value
    Dho = fit_parameters['out_headgroup_thickness'].value
    ## Shared
    # Smearing factor
    sig = fit_parameters['sigma'].value

    in_x_values = np.arange(-40, 0.2, 0.2)
    out_x_values = np.arange(-0.2, 40, 0.2)
    in_water_prob = water(Vci, Vhi, Ali, Dhi, sig, in_x_values)
    out_water_prob = water(Vco, Vho, Alo, Dho, sig, out_x_values)

    # If the probability is above 0 - do nothing. If it's less than 0, scale the value (penalty)
    for value in in_water_prob:
        if value >= 0:
            scaled_water.append(0)
        else:
            scaled_water.append((value**2)*(10**7))
    for value in out_water_prob:
        if value >= 0:
            scaled_water.append(0)
        else:
            scaled_water.append((value**2)*(10**7))

    # Make an array of residuals
    for data in datas:
        # Get error
        current_error = []
        # Check for 0's
        for value in data.error_value[data.min_index:data.max_index]:
            if value == 0:
                value = 1
            
            current_error.append(value)

        # Do math
        current_residual = data.intensity_value[data.min_index:data.max_index] - calc_asym_model(fit_parameters, data.q_value[data.min_index:data.max_index], data, sff)
        # Weight for error
        weighted_residual = np.power(current_residual, 2) / np.power(current_error, 2)
        # Append
        combined_residuals.extend(weighted_residual)

    combined_residuals.extend(scaled_water)

    return combined_residuals

# Augments per data-set
# Calculate augmentations for each dataset as requried
def adjust_b_values(data, in_sample_lipids, out_sample_lipids, water, d_water, temp):
    # Query
    augments = Data_Sample_Lipid_Augment.objects.filter(data_set_title=data)

    in_augment_dict = {}
    out_augment_dict = {}

    for sample_lipid in in_sample_lipids:
        for augment in augments:
            if augment.sample_lipid_name == sample_lipid:
                in_augment_dict[sample_lipid] = augment

    for sample_lipid in out_sample_lipids:
        for augment in augments:
            if augment.sample_lipid_name == sample_lipid:
                out_augment_dict[sample_lipid] = augment

    # Temp
    x = temp

    # Declare
    in_terminal_methyl_b = 0
    in_chain_b = 0
    in_headgroup_b = 0

    out_terminal_methyl_b = 0
    out_chain_b = 0
    out_headgroup_b = 0

    water_b = 0
    calculated_water_volume = 0

    # Calculate water volume
    calculated_water_volume = (
        (
            eval(d_water.total_volume_equation) * data.d2o_mol_fraction
        ) + (
            eval(water.total_volume_equation) * (1 - data.d2o_mol_fraction)
        )
    )

    # If it's x-ary data, work with the electron values
    if data.data_type == 'XR':
        # bw
        water_b = water.electrons

        ## Inner
        # Check each lipid in the sample
        for sample_lipid in in_sample_lipids:
            # If it's a standard lipid, take values from there [b_value = electrons * mol fraction]
            if sample_lipid.sample_lipid_name.project_lipid_name:
                # bt
                in_terminal_methyl_b = in_terminal_methyl_b + (sample_lipid.sample_lipid_name.project_lipid_name.tm_electrons * sample_lipid.lipid_mol_fraction)
                # bc
                in_chain_b = in_chain_b + (sample_lipid.sample_lipid_name.project_lipid_name.tg_electrons * sample_lipid.lipid_mol_fraction)
                # bh
                in_headgroup_b = in_headgroup_b + (sample_lipid.sample_lipid_name.project_lipid_name.hg_electrons * sample_lipid.lipid_mol_fraction)
            # If it's a custom lipid check for user entered values [b_value = electrons * mol fraction]
            else:
                # bt
                in_terminal_methyl_b = in_terminal_methyl_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tm_electrons * sample_lipid.lipid_mol_fraction)
                # bc
                in_chain_b = in_chain_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tg_electrons * sample_lipid.lipid_mol_fraction)
                # bh
                in_headgroup_b = in_headgroup_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.hg_electrons * sample_lipid.lipid_mol_fraction)
        ## Outer
        # Check each lipid in the sample
        for sample_lipid in out_sample_lipids:
            # If it's a standard lipid, take values from there [b_value = electrons * mol fraction]
            if sample_lipid.sample_lipid_name.project_lipid_name:
                # bt
                out_terminal_methyl_b = out_terminal_methyl_b + (sample_lipid.sample_lipid_name.project_lipid_name.tm_electrons * sample_lipid.lipid_mol_fraction)
                # bc
                out_chain_b = out_chain_b + (sample_lipid.sample_lipid_name.project_lipid_name.tg_electrons * sample_lipid.lipid_mol_fraction)
                # bh
                out_headgroup_b = out_headgroup_b + (sample_lipid.sample_lipid_name.project_lipid_name.hg_electrons * sample_lipid.lipid_mol_fraction)
            # If it's a custom lipid check for user entered values [b_value = electrons * mol fraction]
            else:
                # bt
                out_terminal_methyl_b = out_terminal_methyl_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tm_electrons * sample_lipid.lipid_mol_fraction)
                # bc
                out_chain_b = out_chain_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tg_electrons * sample_lipid.lipid_mol_fraction)
                # bh
                out_headgroup_b = out_headgroup_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.hg_electrons * sample_lipid.lipid_mol_fraction)
    
    # If it's neturon data, work with scattering values
    else:
        # bw
        water_b = (d_water.scattering_length * data.d2o_mol_fraction) + (water.scattering_length * (1 - data.d2o_mol_fraction))

        # Inner
        # Check each lipid in the sample
        for sample_lipid in in_sample_lipids:
            
            # If it's a standard lipid, take values from there
            if sample_lipid.sample_lipid_name.project_lipid_name:
                
                # If there is an augmentation (i.e. the user has indicated it's deuterated), calculate the scaled values [b_value = (scattering length + net change) * mol fraction]
                if sample_lipid in in_augment_dict.keys():
                    
                    # For standard augmentations
                    if in_augment_dict[sample_lipid].sample_lipid_augment:
                        # bt
                        in_terminal_methyl_b = in_terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering + in_augment_dict[sample_lipid].sample_lipid_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bc
                        in_chain_b = in_chain_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering + in_augment_dict[sample_lipid].sample_lipid_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bh
                        in_headgroup_b = in_headgroup_b + ((sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering + in_augment_dict[sample_lipid].sample_lipid_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    
                    # For custom augmentations
                    else:
                        # bt
                        in_terminal_methyl_b = in_terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering + in_augment_dict[sample_lipid].sample_lipid_custom_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bc
                        in_chain_b = in_chain_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering + in_augment_dict[sample_lipid].sample_lipid_custom_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bh
                        in_headgroup_b = in_headgroup_b + ((sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering + in_augment_dict[sample_lipid].sample_lipid_custom_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                
                # If there isn't any augmentation, take the scatting length without scaling it [b_value = scattering length * mol fraction]
                else:
                    # bt
                    in_terminal_methyl_b = in_terminal_methyl_b + (sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering * sample_lipid.lipid_mol_fraction)
                    # bc
                    in_chain_b = in_chain_b + (sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering * sample_lipid.lipid_mol_fraction)
                    # bh
                    in_headgroup_b = in_headgroup_b + (sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering * sample_lipid.lipid_mol_fraction)
            
            # If it's a custom lipid check for user entered values
            else:
                
                # If there is an augmentation (i.e. the user has indicated it's deuterated), calculate the scaled values [b_value = (scattering length + net change) * mol fraction]
                if sample_lipid in in_augment_dict.keys():
                    # bt
                    in_terminal_methyl_b = in_terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.tm_scattering + in_augment_dict[sample_lipid].sample_lipid_custom_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    # bc
                    in_chain_b = in_chain_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.tg_scattering + in_augment_dict[sample_lipid].sample_lipid_custom_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    # bh
                    in_headgroup_b = in_headgroup_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.hg_scattering + in_augment_dict[sample_lipid].sample_lipid_custom_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                
                # If there isn't any augmentation, take the scatting length without scaling it [b_value = scattering length * mol fraction]
                else:
                    # bt
                    in_terminal_methyl_b = in_terminal_methyl_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tm_scattering * sample_lipid.lipid_mol_fraction)
                    # bc
                    in_chain_b = in_chain_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tg_scattering * sample_lipid.lipid_mol_fraction)
                    # bh
                    in_headgroup_b = in_headgroup_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.hg_scattering * sample_lipid.lipid_mol_fraction)

        # Outer
        # Check each lipid in the sample
        for sample_lipid in out_sample_lipids:
            
            # If it's a standard lipid, take values from there
            if sample_lipid.sample_lipid_name.project_lipid_name:
                
                # If there is an augmentation (i.e. the user has indicated it's deuterated), calculate the scaled values [b_value = (scattering length + net change) * mol fraction]
                if sample_lipid in out_augment_dict.keys():
                    
                    # For standard augmentations
                    if out_augment_dict[sample_lipid].sample_lipid_augment:
                        # bt
                        out_terminal_methyl_b = out_terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering + out_augment_dict[sample_lipid].sample_lipid_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bc
                        out_chain_b = out_chain_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering + out_augment_dict[sample_lipid].sample_lipid_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bh
                        out_headgroup_b = out_headgroup_b + ((sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering + out_augment_dict[sample_lipid].sample_lipid_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    
                    # For custom augmentations
                    else:
                        # bt
                        out_terminal_methyl_b = out_terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering + out_augment_dict[sample_lipid].sample_lipid_custom_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bc
                        out_chain_b = out_chain_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering + out_augment_dict[sample_lipid].sample_lipid_custom_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bh
                        out_headgroup_b = out_headgroup_b + ((sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering + out_augment_dict[sample_lipid].sample_lipid_custom_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                
                # If there isn't any augmentation, take the scatting length without scaling it [b_value = scattering length * mol fraction]
                else:
                    # bt
                    out_terminal_methyl_b = out_terminal_methyl_b + (sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering * sample_lipid.lipid_mol_fraction)
                    # bc
                    out_chain_b = out_chain_b + (sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering * sample_lipid.lipid_mol_fraction)
                    # bh
                    out_headgroup_b = out_headgroup_b + (sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering * sample_lipid.lipid_mol_fraction)
            
            # If it's a custom lipid check for user entered values
            else:
                
                # If there is an augmentation (i.e. the user has indicated it's deuterated), calculate the scaled values [b_value = (scattering length + net change) * mol fraction]
                if  sample_lipid in out_augment_dict.keys():
                    # bt
                    out_terminal_methyl_b = out_terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.tm_scattering + out_augment_dict[sample_lipid].sample_lipid_custom_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    # bc
                    out_chain_b = out_chain_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.tg_scattering + out_augment_dict[sample_lipid].sample_lipid_custom_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    # bh
                    out_headgroup_b = out_headgroup_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.hg_scattering + out_augment_dict[sample_lipid].sample_lipid_custom_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                
                # If there isn't any augmentation, take the scatting length without scaling it [b_value = scattering length * mol fraction]
                else:
                    # bt
                    out_terminal_methyl_b = out_terminal_methyl_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tm_scattering * sample_lipid.lipid_mol_fraction)
                    # bc
                    out_chain_b = out_chain_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tg_scattering * sample_lipid.lipid_mol_fraction)
                    # bh
                    out_headgroup_b = out_headgroup_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.hg_scattering * sample_lipid.lipid_mol_fraction)

    # Return all the b values in an array to be referenced later on
    b_values = [in_chain_b, in_headgroup_b, in_terminal_methyl_b, out_chain_b, out_headgroup_b, out_terminal_methyl_b, water_b, calculated_water_volume]
    
    return(b_values)

# Parameters
def asymmetrical_paramitize(parameter, in_sample_lipids, out_sample_lipids, datas, temp):
    ## DELCARE
    # Other molecules
    water = Molecule.objects.get(compound_name='water')
    d_water = Molecule.objects.get(compound_name='deuterated_water')
    
    # Parameters
    # Check each value in the database and prepare a parameter (lmfit) object that includes each of them
    fit_parameters = lsq.Parameters()
    fit_parameters.add_many(
        # Inner
        ( # Vc
            'in_chain_volume',
            parameter.in_chain_volume,
            not(parameter.in_chain_volume_lock),
        ),
        ( # Vh
            'in_headgroup_volume',
            parameter.in_headgroup_volume,
            not(parameter.in_headgroup_volume_lock),
        ),
        ( # Vt
            'in_terminal_methyl_volume',
            parameter.in_terminal_methyl_volume,
            not(parameter.in_terminal_methyl_volume_lock),
            parameter.in_terminal_methyl_volume_lowerbound,
            parameter.in_terminal_methyl_volume_upperbound
        ),
        # Unknown
        ( # Al
            'in_area_per_lipid',
            parameter.in_lipid_area,
            not(parameter.in_lipid_area_lock),
            parameter.in_lipid_area_lowerbound,
            parameter.in_lipid_area_upperbound
        ),
        ( # Dh
            'in_headgroup_thickness',
            parameter.in_headgroup_thickness,
            not(parameter.in_headgroup_thickness_lock),
            parameter.in_headgroup_thickness_lowerbound,
            parameter.in_headgroup_thickness_upperbound
        ),
        # Outer
        ( # Vc
            'out_chain_volume',
            parameter.out_chain_volume,
            not(parameter.out_chain_volume_lock),
        ),
        ( # Vh
            'out_headgroup_volume',
            parameter.out_headgroup_volume,
            not(parameter.out_headgroup_volume_lock),
        ),
        ( # Vt
            'out_terminal_methyl_volume',
            parameter.out_terminal_methyl_volume,
            not(parameter.out_terminal_methyl_volume_lock),
            parameter.out_terminal_methyl_volume_lowerbound,
            parameter.out_terminal_methyl_volume_upperbound
        ),
        # Unknown
        ( # Al
            'out_area_per_lipid',
            parameter.out_lipid_area,
            not(parameter.out_lipid_area_lock),
            parameter.out_lipid_area_lowerbound,
            parameter.out_lipid_area_upperbound
        ),
        ( # Dh
            'out_headgroup_thickness',
            parameter.out_headgroup_thickness,
            not(parameter.out_headgroup_thickness_lock),
            parameter.out_headgroup_thickness_lowerbound,
            parameter.out_headgroup_thickness_upperbound
        ),
        # Smearing factor
        ( # Sigma
            'sigma',
            parameter.sigma,
            not(parameter.sigma_lock),
            parameter.sigma_lowerbound,
            parameter.sigma_upperbound
        ),
        ## Separated form factor
        ( # Average vesicle radius
            'average_vesicle_radius',
            parameter.average_vesicle_radius,
            not(parameter.average_vesicle_radius_lock),
            parameter.average_vesicle_radius_upperbound,
            parameter.average_vesicle_radius_lowerbound
        ),
        ( # Relative size polydispersity
            'relative_size',
            parameter.relative_size,
            not(parameter.relative_size_lock),
            parameter.relative_size_upperbound,
            parameter.relative_size_lowerbound
        )
    )

    # Multiple datasets
    # Create b value parameters for each dataset using adjuest_b_values() - add these values to the parameter object and label them with the id number for their respective dataset
    try:
        for data in datas:
            # Get values for this data
            b_values = adjust_b_values(data, in_sample_lipids, out_sample_lipids, water, d_water, temp)

            fit_parameters.add_many(
                # Inner
                ( # bc
                    'in_chain_b_%i' % data.id,
                    b_values[0],
                    False
                ),
                ( # bh
                    'in_headgroup_b_%i' % data.id,
                    b_values[1],
                    False
                ),
                ( # bt
                    'in_terminal_methyl_b_%i' % data.id,
                    b_values[2],
                    False
                ),
                # Outer
                ( # bc
                    'out_chain_b_%i' % data.id,
                    b_values[3],
                    False
                ),
                ( # bh
                    'out_headgroup_b_%i' % data.id,
                    b_values[4],
                    False
                ),
                ( # bt
                    'out_terminal_methyl_b_%i' % data.id,
                    b_values[5],
                    False
                ),
                # Shared
                ( # bw
                    'water_b_%i' % data.id,
                    b_values[6],
                    False
                ),
                ( # Vw
                    'combined_water_volume_%i' % data.id,
                    b_values[7],
                    False
                ),
                # Tweaks
                ( # Scale
                    'scale_%i' % data.id,
                    data.scale,
                    not(data.scale_lock),
                    data.scale_lowerbound,
                    data.scale_upperbound
                ),
                ( # BG
                    'background_%i' % data.id,
                    data.background,
                    not(data.background_lock),
                    data.background_lowerbound,
                    data.background_upperbound
                )
            )
    # Single dataset
    # If there is only one dataset, the above loop will break and you will cry - therefore do the same thing but without a loop
    except TypeError:
        # Get values for this data
        b_values = adjust_b_values(datas, in_sample_lipids, out_sample_lipids, water, d_water, temp)

        fit_parameters.add_many(
            # Inner
            ( # bc
                'in_chain_b_%i' % datas.id,
                b_values[0],
                False
            ),
            ( # bh
                'in_headgroup_b_%i' % datas.id,
                b_values[1],
                False
            ),
            ( # bt
                'in_terminal_methyl_b_%i' % datas.id,
                b_values[2],
                False
            ),
            # Outer
            ( # bc
                'out_chain_b_%i' % datas.id,
                b_values[3],
                False
            ),
            ( # bh
                'out_headgroup_b_%i' % datas.id,
                b_values[4],
                False
            ),
            ( # bt
                'out_terminal_methyl_b_%i' % datas.id,
                b_values[5],
                False
            ),
            # Shared
            ( # bw
                'water_b_%i' % datas.id,
                b_values[6],
                False
            ),
            ( # Vw
                'combined_water_volume_%i' % datas.id,
                b_values[7],
                False
            ),
            # Tweaks
            ( # Scale
                'scale_%i' % datas.id,
                datas.scale,
                not(datas.scale_lock),
                datas.scale_lowerbound,
                datas.scale_upperbound
            ),
            ( # BG
                'background_%i' % datas.id,
                datas.background,
                not(datas.background_lock),
                datas.background_lowerbound,
                datas.background_upperbound
            )
        )

    return fit_parameters

# Graphs / fit / probabilities / etc
def asymmetrical_graph(parameter, in_sample_lipids, out_sample_lipids, data, temp):
    # Get parameters
    fit_parameters = asymmetrical_paramitize(parameter, in_sample_lipids, out_sample_lipids, data, temp)

    # Get result
    model_result = calc_asym_model(
        fit_parameters,
        data.q_value[data.min_index:data.max_index],
        data,
        parameter.separated
    )

    return model_result

def asymmetrical_fit(parameter, in_sample_lipids, out_sample_lipids, datas, temp):
    # Get parameters
    fit_parameters = asymmetrical_paramitize(parameter, in_sample_lipids, out_sample_lipids, datas, temp)

    # Get result
    x = None

    fit_result = lsq.minimize(
        asymmetrical_objective_function,
        fit_parameters,
        args=(x, datas, parameter.separated)
    )

    return fit_result

def asymmetrical_sdp(parameter, in_head_prob, in_methyl_prob, in_tm_prob, in_water_prob, out_head_prob, out_methyl_prob, out_tm_prob, out_water_prob, in_sample_lipids, out_sample_lipids, data, temp):
    # Declare
    in_sdp_final = []
    out_sdp_final = []

    # Get parameters
    fit_parameters = asymmetrical_paramitize(parameter, in_sample_lipids, out_sample_lipids, data, temp)

    ## Unpack parameters
    # Inner
    Vci = fit_parameters['in_chain_volume'].value
    Vhi = fit_parameters['in_headgroup_volume'].value
    Vti = fit_parameters['in_terminal_methyl_volume'].value
    # Outer
    Vco = fit_parameters['out_chain_volume'].value
    Vho = fit_parameters['out_headgroup_volume'].value
    Vto = fit_parameters['out_terminal_methyl_volume'].value
    # Water
    Vw = fit_parameters['combined_water_volume_%i' % data.id].value
    ## b values
    # Inner
    bci = fit_parameters['in_chain_b_%i' % data.id].value
    bhi = fit_parameters['in_headgroup_b_%i' % data.id].value
    bti = fit_parameters['in_terminal_methyl_b_%i' % data.id].value
    # Outter
    bco = fit_parameters['out_chain_b_%i' % data.id].value
    bho = fit_parameters['out_headgroup_b_%i' % data.id].value
    bto = fit_parameters['out_terminal_methyl_b_%i' % data.id].value
    # Water
    bw = fit_parameters['water_b_%i' % data.id].value

    if Vhi == 0 or Vci == 0 or Vti == 0 or Vho == 0 or Vco == 0 or Vto == 0 or Vw == 0:
        combined_sdp = 0
        return(combined_sdp)

    ## Scale probabilities
    # Inner
    in_scaled_head_prob = (np.asarray(in_head_prob)*(bhi/Vhi))
    in_scaled_methyl_prob = (np.asarray(in_methyl_prob)*(bci/Vci))
    in_scaled_tm_prob = (np.asarray(in_tm_prob)*(bti/Vti))
    # Outter
    out_scaled_head_prob = (np.asarray(out_head_prob)*(bho/Vho))
    out_scaled_methyl_prob = (np.asarray(out_methyl_prob)*(bco/Vco))
    out_scaled_tm_prob = (np.asarray(out_tm_prob)*(bto/Vto))
    # Water
    in_scaled_water_prob = (np.asarray(in_water_prob)*(bw/Vw))
    out_scaled_water_prob = (np.asarray(out_water_prob)*(bw/Vw))

    # Add them together
    for in_h_value, in_m_value, in_tm_value, in_w_value in zip(in_scaled_head_prob, in_scaled_methyl_prob, in_scaled_tm_prob, in_scaled_water_prob):
        in_sdp_final.append(in_h_value+in_m_value+in_tm_value+in_w_value)

    for out_h_value, out_m_value, out_tm_value, out_w_value in zip(out_scaled_head_prob, out_scaled_methyl_prob, out_scaled_tm_prob, out_scaled_water_prob):
        out_sdp_final.append(out_h_value+out_m_value+out_tm_value+out_w_value)

    combined_sdp = [in_sdp_final, out_sdp_final, in_scaled_head_prob, in_scaled_methyl_prob, in_scaled_tm_prob, out_scaled_head_prob, out_scaled_methyl_prob, out_scaled_tm_prob, in_scaled_water_prob, out_water_prob]

    return(combined_sdp)

def asym_additional_parameters(parameter, in_sample_lipids, out_sample_lipids, data, temp, in_head_prob, out_head_prob, in_x_values, out_x_values):
    # Declare
    additional_parameters = []

    # Get parameters
    fit_parameters = asymmetrical_paramitize(parameter, in_sample_lipids, out_sample_lipids, data, temp)

    ## Inner
    # Calculated
    Ali = fit_parameters['in_area_per_lipid'].value
    # Shared
    Vhi = fit_parameters['in_headgroup_volume'].value
    Vci = fit_parameters['in_chain_volume'].value
    Vti = fit_parameters['in_terminal_methyl_volume'].value
    ## Outter
    # Calculated
    Alo = fit_parameters['out_area_per_lipid'].value
    # Shared
    Vho = fit_parameters['out_headgroup_volume'].value
    Vco = fit_parameters['out_chain_volume'].value
    Vto = fit_parameters['out_terminal_methyl_volume'].value

    # If Al for the inner or the outer layer is 0 don't bother calculating
    if Ali == 0 or Alo == 0:
        Dbi = 0
        Dci = 0
        Dbo = 0
        Dco = 0
    # Calculate normal paramteres with regular old math
    else:
        Dbi = (2 * (Vci + Vhi)) / Ali
        Dci = ((2 * Vci) / Ali) / 2
        Dbo = (2 * (Vco + Vho)) / Alo
        Dco = ((2 * Vco) / Alo) / 2

    # Find peaks
    # Take the head probability values for each half and use scipy 'signal' to check for the peak in each
    in_peak_index = sig.argrelextrema(in_head_prob, np.greater)[0]
    out_peak_index = sig.argrelextrema(out_head_prob, np.greater)[0]

    # Calculate distance
    # Take the x value for the peaks and do your regular old distance calculation for each side
    # inner peak -> 0 && 0 <- outer peak
    in_distance = (np.sqrt( (in_x_values[in_peak_index] - 0)**2 + (in_head_prob[in_peak_index] - in_head_prob[-1])**2 ))
    out_distance = (np.sqrt( (0 - out_x_values[out_peak_index])**2 + (out_head_prob[0] - out_head_prob[out_peak_index])**2 ))

    additional_parameters.append(round(Dbi, 2))
    additional_parameters.append(round(Dci, 2))
    additional_parameters.append(round(Dbo, 2))
    additional_parameters.append(round(Dco, 2))
    additional_parameters.append(round(abs(in_distance[0]), 2))
    additional_parameters.append(round(abs(out_distance[0]), 2))

    return(additional_parameters)
