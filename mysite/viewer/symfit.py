# Other libraries
import numpy as np
import lmfit as lsq

# Models
from .models import Molecule
from .models import Data_Sample_Lipid_Augment

# Other imports
from .probabilities import *
from scipy import signal as sig

# Symmetrical model
# The math only - used by LMFit to run each iteration
def sym_model(
    q,      # independant
    Vc,     # chain_volume
    Vh,     # headgroup_volume
    Vt,     # terminal_methyl_volume
    Vw,     # water_volume
    Al,     # area_per_lipid
    Dh,     # headgroup_thickness
    sig,    # smearing factor
    bc,     # chain_b
    bh,     # headgroup_b
    bt,     # terminal_methyl_b
    bw,     # water_b
    scale,  # scale
    bg      # bg
):
    return (
        (q**(-2) * scale*(
            (
                ((2*(np.exp(-((q*sig)**2)/2)))/(q*Dh*Al*Vt*Vw*(Vc-2*Vt)))
                * np.abs(
                    Vt*(bw*(Al*Dh-Vh)*(Vc-2*Vt)+Vw*bh*(Vc-2*Vt)-Vw*Al*Dh*(bc-2*bt))
                    * np.sin(q*Vc/Al) 
                    + Vt*(Vc-2*Vt)*(bw*Vh-bh*Vw)*np.sin(q*Dh+q*Vc/Al)
                    + Vw*Al*Dh*(bc*Vt-bt*Vc)*np.sin(2*q*Vt/Al)
                )
            ) **2
        ) + bg)
    )

# Symmetrical model with structure factor Sn(q)
def sym_model_structure_factor(
    q,      # independant
    Vc,     # chain_volume
    Vh,     # headgroup_volume
    Vt,     # terminal_methyl_volume
    Vw,     # water_volume
    Al,     # area_per_lipid
    Dh,     # headgroup_thickness
    sig,    # smearing factor
    r,      # Average vesicle radius
    rs,     # Relative size polydispersity
    bc,     # chain_b
    bh,     # headgroup_b
    bt,     # terminal_methyl_b
    bw,     # water_b
    scale,  # scale
    bg,     # bg
    
    # structure factor input
    N,      # number of bilayers 
):

    # I(q) = q^-2 * P(q)                # original
    # I(q) = q^-2 * P(q) * scale + bg   # modified for the plot
    # I(q) = q^-2 * P(q) * S(q)         # extra data for structure factor
    print("sf")
    return (
        q**(-2) * scale * (
            # Pv(r,rs,q)
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
            # Fb(q)^2
            ) * (
                ((2*(np.exp(-((q*sig)**2)/2)))/(q*Dh*Al*Vt*Vw*(Vc-2*Vt)))
                * np.abs(
                    Vt*(bw*(Al*Dh-Vh)*(Vc-2*Vt)+Vw*bh*(Vc-2*Vt)-Vw*Al*Dh*(bc-2*bt))
                    * np.sin(q*Vc/Al) 
                    + Vt*(Vc-2*Vt)*(bw*Vh-bh*Vw)*np.sin(q*Dh+q*Vc/Al)
                    + Vw*Al*Dh*(bc*Vt-bt*Vc)*np.sin(2*q*Vt/Al)
                )
            ) **2
            # Sn(q)
            * (2)
            # * (N + 2 * np.sum([(N - k) * np.cos(k*q*d) * np.exp(-(d * q / 2 * np.pi)) for k in range(1, N-1)]))
        ) + bg
    )

# Symmetrical model and separated form factor
# The same math, but including the separated form factor equation
def sym_model_separated(
    q,      # independant
    Vc,     # chain_volume
    Vh,     # headgroup_volume
    Vt,     # terminal_methyl_volume
    Vw,     # water_volume
    Al,     # area_per_lipid
    Dh,     # headgroup_thickness
    sig,    # smearing factor
    r,      # Average vesicle radius
    rs,     # Relative size polydispersity
    bc,     # chain_b
    bh,     # headgroup_b
    bt,     # terminal_methyl_b
    bw,     # water_b
    scale,  # scale
    bg      # bg
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
                ((2*(np.exp(-((q*sig)**2)/2)))/(q*Dh*Al*Vt*Vw*(Vc-2*Vt)))
                * np.abs(
                    Vt*(bw*(Al*Dh-Vh)*(Vc-2*Vt)+Vw*bh*(Vc-2*Vt)-Vw*Al*Dh*(bc-2*bt))
                    * np.sin(q*Vc/Al) 
                    + Vt*(Vc-2*Vt)*(bw*Vh-bh*Vw)*np.sin(q*Dh+q*Vc/Al)
                    + Vw*Al*Dh*(bc*Vt-bt*Vc)*np.sin(2*q*Vt/Al)
                )
            ) **2
        ) + bg
    )

# Calculate result from model for an individual dataset
# Connetc the model, the dataset and the parameters
def calc_sym_model(fit_parameters, q, data, sff, use_structure_factor):
    print("calculating sym model...")
    # Convert array
    q_array = np.array(q)

    ## Unpack parameters
    # Shared
    Vc = fit_parameters['chain_volume'].value
    Vh = fit_parameters['headgroup_volume'].value
    Vt = fit_parameters['terminal_methyl_volume'].value
    # Unknown
    Al = fit_parameters['area_per_lipid'].value
    Dh = fit_parameters['headgroup_thickness'].value
    # Smearing factor
    sig = fit_parameters['sigma'].value
    # Separated form factor
    r = fit_parameters['average_vesicle_radius'].value
    rs = fit_parameters['relative_size'].value
    # Per dataset
    bc = fit_parameters['chain_b_%i' % data.id].value
    bh = fit_parameters['headgroup_b_%i' % data.id].value
    bt = fit_parameters['terminal_methyl_b_%i' % data.id].value
    bw = fit_parameters['water_b_%i' % data.id].value
    Vw = fit_parameters['combined_water_volume_%i' % data.id].value
    # Tweaks
    scale = fit_parameters['scale_%i' % data.id].value
    bg = fit_parameters['background_%i' % data.id].value

    # Return the calculated model
    if sff:
        calc_result = sym_model_separated(q_array, Vc, Vh, Vt, Vw, Al, Dh, sig, r, rs, bc, bh, bt, bw, scale, bg)
    elif use_structure_factor:
        calc_result = sym_model_structure_factor(q_array, Vc, Vh, Vt, Vw, Al, Dh, sig, bc, bh, bt, bw, scale, bg)
    else:
        calc_result = sym_model(q_array, Vc, Vh, Vt, Vw, Al, Dh, sig, bc, bh, bt, bw, scale, bg)

    # Make whole result negative if water prob is neg
    return calc_result


# Objective function 
# Create a residual for each dataset, then flatten into a single array for minimize()
def symmetrical_objective_function(fit_parameters, x, datas, sff, use_structure_factor):
    # Delcare
    current_residual = []
    combined_residuals = []
    scaled_water = []

    # Check if the water prob is negative - if it is, impose a penalty

    ## Unpack parameters
    # Shared
    Vc = fit_parameters['chain_volume'].value
    Vh = fit_parameters['headgroup_volume'].value
    Vt = fit_parameters['terminal_methyl_volume'].value
    # Unknown
    Al = fit_parameters['area_per_lipid'].value
    Dh = fit_parameters['headgroup_thickness'].value
    # Smearing factor
    sig = fit_parameters['sigma'].value

    x_values = np.arange(-40, 40, 0.2)
    water_prob = water(Vc, Vh, Al, Dh, sig, x_values)

    # If the probability is above 0 - do nothing. If it's less than 0, scale the value (penalty)
    for value in water_prob:
        if value >= 0:
            scaled_water.append(0)
        else:
            scaled_water.append((value**2)*(10**7))

    # Make an array of residuals
    for data in datas:
        # Get error
        current_error = []
        # Check for 0's - make them 1 so they have no effect
        for value in data.error_value[data.min_index:data.max_index]:
            if value == 0:
                value = 1
            
            current_error.append(value)

        # Do math
        current_residual = data.intensity_value[data.min_index:data.max_index] - calc_sym_model(fit_parameters, data.q_value[data.min_index:data.max_index], data, sff, use_structure_factor)
        # Weight for error
        weighted_residual = np.power(current_residual, 2) / np.power(current_error, 2)
        # Append
        combined_residuals.extend(weighted_residual)

    combined_residuals.extend(scaled_water)

    return combined_residuals

# Augments per data-set
# Calculate augmentations for each dataset as requried
def adjust_b_values(data, sample_lipids, water, d_water, temp):
    # Query
    augments = Data_Sample_Lipid_Augment.objects.filter(data_set_title=data)

    augment_dict = {}

    for sample_lipid in sample_lipids:
        for augment in augments:
            if augment.sample_lipid_name == sample_lipid:
                augment_dict[sample_lipid] = augment

    # Temp
    x = temp

    # Declare
    terminal_methyl_b = 0
    chain_b = 0
    headgroup_b = 0
    water_b = 0
    calculated_water_volume = 0

    # Calculate water volume based on the D2O percentage
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

        # Check each lipid in the sample 
        for sample_lipid in sample_lipids:
            # If it's a standard lipid, take values from there [b_value = electrons * mol fraction]
            if sample_lipid.sample_lipid_name.project_lipid_name:
                # bt
                terminal_methyl_b = terminal_methyl_b + (sample_lipid.sample_lipid_name.project_lipid_name.tm_electrons * sample_lipid.lipid_mol_fraction)
                # bc
                chain_b = chain_b + (sample_lipid.sample_lipid_name.project_lipid_name.tg_electrons * sample_lipid.lipid_mol_fraction)
                # bh
                headgroup_b = headgroup_b + (sample_lipid.sample_lipid_name.project_lipid_name.hg_electrons * sample_lipid.lipid_mol_fraction)
            # If it's a custom lipid check for user entered values [b_value = electrons * mol fraction]
            else:
                # bt
                terminal_methyl_b = terminal_methyl_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tm_electrons * sample_lipid.lipid_mol_fraction)
                # bc
                chain_b = chain_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tg_electrons * sample_lipid.lipid_mol_fraction)
                # bh
                headgroup_b = headgroup_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.hg_electrons * sample_lipid.lipid_mol_fraction)
    
    # If it's neturon data, work with scattering values
    else:
        # bw
        water_b = (d_water.scattering_length * data.d2o_mol_fraction) + (water.scattering_length * (1 - data.d2o_mol_fraction))

        # Check each lipid in the sample
        for sample_lipid in sample_lipids:

            # If it's a standard lipid, take values from there
            if sample_lipid.sample_lipid_name.project_lipid_name:

                # If there is an augmentation (i.e. the user has indicated it's deuterated), calculate the scaled values [b_value = (scattering length + net change) * mol fraction]
                if sample_lipid in augment_dict.keys():

                    # For standard augmentations
                    if augment_dict[sample_lipid].sample_lipid_augment:
                        # bt
                        terminal_methyl_b = terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering + augment_dict[sample_lipid].sample_lipid_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bc
                        chain_b = chain_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering + augment_dict[sample_lipid].sample_lipid_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bh
                        headgroup_b = headgroup_b + ((sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering + augment_dict[sample_lipid].sample_lipid_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    
                    # For custom augmentations
                    else:
                        # bt
                        terminal_methyl_b = terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering + augment_dict[sample_lipid].sample_lipid_custom_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bc
                        chain_b = chain_b + ((sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering + augment_dict[sample_lipid].sample_lipid_custom_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                        # bh
                        headgroup_b = headgroup_b + ((sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering + augment_dict[sample_lipid].sample_lipid_custom_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                
                # If there isn't any augmentation, take the scatting length without scaling it [b_value = scattering length * mol fraction]
                else:
                    # bt
                    terminal_methyl_b = terminal_methyl_b + (sample_lipid.sample_lipid_name.project_lipid_name.tm_scattering * sample_lipid.lipid_mol_fraction)
                    # bc
                    chain_b = chain_b + (sample_lipid.sample_lipid_name.project_lipid_name.tg_scattering * sample_lipid.lipid_mol_fraction)
                    # bh
                    headgroup_b = headgroup_b + (sample_lipid.sample_lipid_name.project_lipid_name.hg_scattering * sample_lipid.lipid_mol_fraction)
            
            
            # If it's a custom lipid check for user entered values
            else:

                # If there is an augmentation (i.e. the user has indicated it's deuterated), calculate the scaled values [b_value = (scattering length + net change) * mol fraction]
                if sample_lipid in augment_dict.keys():
                    # bt
                    terminal_methyl_b = terminal_methyl_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.tm_scattering + augment_dict[sample_lipid].sample_lipid_custom_augment.tmg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    # bc
                    chain_b = chain_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.tg_scattering + augment_dict[sample_lipid].sample_lipid_custom_augment.tg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                    # bh
                    headgroup_b = headgroup_b + ((sample_lipid.sample_lipid_name.project_user_lipid_name.hg_scattering + augment_dict[sample_lipid].sample_lipid_custom_augment.hg_scattering_net_change) * sample_lipid.lipid_mol_fraction)
                
                # If there isn't any augmentation, take the scatting length without scaling it [b_value = scattering length * mol fraction]
                else:
                    # bt
                    terminal_methyl_b = terminal_methyl_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tm_scattering * sample_lipid.lipid_mol_fraction)
                    # bc
                    chain_b = chain_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.tg_scattering * sample_lipid.lipid_mol_fraction)
                    # bh
                    headgroup_b = headgroup_b + (sample_lipid.sample_lipid_name.project_user_lipid_name.hg_scattering * sample_lipid.lipid_mol_fraction)

    # Return all the b values in an array to be referenced later on
    b_values = [chain_b, headgroup_b, terminal_methyl_b, water_b, calculated_water_volume]
    
    return(b_values)

# Parameters
def symmetrical_paramitize(parameter, sample_lipids, datas, temp, advanced):
    ## DELCARE
    # Other molecules
    water = Molecule.objects.get(compound_name='water')
    d_water = Molecule.objects.get(compound_name='deuterated_water')
    
    # If advanced options are not on, automatically lock volumes
    # Vc
    if advanced:
        vc_lock = not(parameter.chain_volume_lock)
    else:
        vc_lock = False
    # Vh
    if advanced:
        vh_lock = not(parameter.headgroup_volume_lock)
    else:
        vh_lock = False

    # Parameters
    # Check each value in the database and prepare a parameter (lmfit) object that includes each of them
    fit_parameters = lsq.Parameters()
    fit_parameters.add_many(
        # Shared
        ( # Vc
            'chain_volume',
            parameter.chain_volume,
            vc_lock,
        ),
        ( # Vh
            'headgroup_volume',
            parameter.headgroup_volume,
            vh_lock,
        ),
        ( # Vt
            'terminal_methyl_volume',
            parameter.terminal_methyl_volume,
            not(parameter.terminal_methyl_volume_lock),
            parameter.terminal_methyl_volume_lowerbound,
            parameter.terminal_methyl_volume_upperbound
        ),
        # Unknown
        ( # Al
            'area_per_lipid',
            parameter.lipid_area,
            not(parameter.lipid_area_lock),
            parameter.lipid_area_lowerbound,
            parameter.lipid_area_upperbound
        ),
        ( # Dh
            'headgroup_thickness',
            parameter.headgroup_thickness,
            not(parameter.headgroup_thickness_lock),
            parameter.headgroup_thickness_lowerbound,
            parameter.headgroup_thickness_upperbound
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

    #Multiple datasets
    # Create b value parameters for each dataset using adjuest_b_values() - add these values to the parameter object and label them with the id number for their respective dataset
    try:
        for data in datas:
            # Get values for this data
            b_values = adjust_b_values(data, sample_lipids, water, d_water, temp)

            fit_parameters.add_many(
                ( # bc
                    'chain_b_%i' % data.id,
                    b_values[0],
                    False
                ),
                ( # bh
                    'headgroup_b_%i' % data.id,
                    b_values[1],
                    False
                ),
                ( # bt
                    'terminal_methyl_b_%i' % data.id,
                    b_values[2],
                    False
                ),
                ( # bw
                    'water_b_%i' % data.id,
                    b_values[3],
                    False
                ),
                ( # Vw
                    'combined_water_volume_%i' % data.id,
                    b_values[4],
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
        b_values = adjust_b_values(datas, sample_lipids, water, d_water, temp)

        fit_parameters.add_many(
            ( # bc
                'chain_b_%i' % datas.id,
                b_values[0],
                False
            ),
            ( # bh
                'headgroup_b_%i' % datas.id,
                b_values[1],
                False
            ),
            ( # bt
                'terminal_methyl_b_%i' % datas.id,
                b_values[2],
                False
            ),
            ( # bw
                'water_b_%i' % datas.id,
                b_values[3],
                False
            ),
            ( # Vw
                'combined_water_volume_%i' % datas.id,
                b_values[4],
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

# Fit / graphs / etc
def symmetrical_graph(parameter, sample_lipids, data, temp, advanced):
    # Get parameters
    fit_parameters = symmetrical_paramitize(parameter, sample_lipids, data, temp, advanced)

    # Get result
    model_result = calc_sym_model(
        fit_parameters,
        data.q_value[data.min_index:data.max_index],
        data,
        parameter.separated
    )

    return model_result

def symmetrical_fit(parameter, sample_lipids, datas, temp, advanced):
    # Get parameters
    fit_parameters = symmetrical_paramitize(parameter, sample_lipids, datas, temp, advanced)

    # Get result
    x = None

    fit_result = lsq.minimize(
        symmetrical_objective_function,
        fit_parameters,
        args=(x, datas, parameter.separated, parameter.use_structure_factor)
    )

    return fit_result

def symmetrical_sdp(parameter, head_prob, methyl_prob, tm_prob, water_prob, sample_lipids, data, temp, advanced):
    # Declare
    sdp_final = []

    # Get parameters
    fit_parameters = symmetrical_paramitize(parameter, sample_lipids, data, temp, advanced)

    # Shared
    Vh = fit_parameters['headgroup_volume'].value
    Vc = fit_parameters['chain_volume'].value
    Vt = fit_parameters['terminal_methyl_volume'].value
    Vw = fit_parameters['combined_water_volume_%i' % data.id].value

    # Unpack b values
    bh = fit_parameters['headgroup_b_%i' % data.id].value
    bc = fit_parameters['chain_b_%i' % data.id].value
    bt = fit_parameters['terminal_methyl_b_%i' % data.id].value
    bw = fit_parameters['water_b_%i' % data.id].value

    if Vh == 0 or Vc == 0 or Vt == 0 or Vw == 0:
        combined_sdp = 0
        return(combined_sdp)

    # Scale probabilities
    scaled_head_prob = (np.asarray(head_prob)*(bh/Vh))
    scaled_methyl_prob = (np.asarray(methyl_prob)*(bc/Vc))
    scaled_tm_prob = (np.asarray(tm_prob)*(bt/Vt))
    scaled_water_prob = (np.asarray(water_prob)*(bw/Vw))

    # Add them together
    for h_value, m_value, tm_value, w_value in zip(scaled_head_prob, scaled_methyl_prob, scaled_tm_prob, scaled_water_prob):
        sdp_final.append(h_value+m_value+tm_value+w_value)

    combined_sdp = [sdp_final, scaled_head_prob, scaled_methyl_prob, scaled_tm_prob, scaled_water_prob]

    return(combined_sdp)

def sym_additional_parameters(parameter, sample_lipids, data, temp, x_values, head_prob, advanced):
    # Declare
    additional_parameters = []

    # Get parameters
    fit_parameters = symmetrical_paramitize(parameter, sample_lipids, data, temp, advanced)

    # Calculated
    Al = fit_parameters['area_per_lipid'].value
    # Shared
    Vh = fit_parameters['headgroup_volume'].value
    Vc = fit_parameters['chain_volume'].value
    Vt = fit_parameters['terminal_methyl_volume'].value

    # If Al is 0 don't bother calculating
    if Al == 0:
        Db = 0
        Dc = 0
    # Calculate normal paramteres with regular old math
    else:
        Db = (2 * (Vc + Vh)) / Al
        Dc = (2 * Vc) / Al

    # Dhh is calculated based on the peak to peak distance for the headgroup SDP profile. More on why this decision was made in the paper

    # Find peaks
    # Take the head probability values and use scipy 'signal' to check for peaks
    peak_indexes = sig.argrelextrema(head_prob, np.greater)[0]

    # Calculate distance
    # Take the x value for the peaks and do your regular old distance calculation... if no prob has been calculated, skip it
    if peak_indexes.size == 0:
        Dhh = 0
    else:
        Dhh = (np.sqrt( (x_values[peak_indexes[1]] - x_values[peak_indexes[0]])**2 + (head_prob[peak_indexes[1]] - head_prob[peak_indexes[0]])**2 ))

    additional_parameters.append(round(Db, 2))
    additional_parameters.append(round(Dc, 2))
    additional_parameters.append(round(Dhh, 2))

    return(additional_parameters)
