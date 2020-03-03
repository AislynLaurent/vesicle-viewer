# Other libraries
import numpy as np
import lmfit as lsq

# Models
from .models import Data_Lipid
from .models import Data_Lipid_Atom
from .models import Molecule

# Asymmetrical model
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
        return asym_model_separated(q_array, Vci, Vhi, Vti, Ali, Dhi, Vco, Vho, Vto, Alo, Dho, Vw, sig, r, rs, bci, bhi, bti, bco, bho, bto, bw, scale, bg)
    else:
        return asym_model(q_array, Vci, Vhi, Vti, Ali, Dhi, Vco, Vho, Vto, Alo, Dho, Vw, sig, bci, bhi, bti, bco, bho, bto, bw, scale, bg)


# Objective function create a residual for each, then flatten
def asymmetrical_objective_function(fit_parameters, x, datas, sff):
    # Delcare
    residuals = []
    combined_residuals = []

    # Make an array of residuals
    for data in datas:
        # Do math
        residuals = data.intensity_value[data.min_index:data.max_index] - calc_asym_model(fit_parameters, data.q_value[data.min_index:data.max_index], data, sff)

        # Append
        combined_residuals.extend(residuals)

    return combined_residuals

# Augments per data-set
def adjust_b_values(data, in_project_lipids, out_project_lipids, water, d_water, temp):
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
        (eval(d_water.total_volume_equation) * data.d2o_mol_fraction) 
        + (eval(water.total_volume_equation) * (1 - data.d2o_mol_fraction))
    )

    if data.data_type == 'XR':
        # bw
        water_b = water.electrons

        ## Inner
        for project_lipid in in_project_lipids:
            # bt
            in_terminal_methyl_b = in_terminal_methyl_b + (project_lipid.project_lipid_name.tm_electrons * project_lipid.lipid_mol_fraction)
            # bc
            in_chain_b = in_chain_b + (project_lipid.project_lipid_name.tg_electrons * project_lipid.lipid_mol_fraction)
            # bh
            in_headgroup_b = in_headgroup_b + (project_lipid.project_lipid_name.hg_electrons * project_lipid.lipid_mol_fraction)
        ## Outer
        for project_lipid in out_project_lipids:
            ## Outer
            # bt
            out_terminal_methyl_b = out_terminal_methyl_b + (project_lipid.project_lipid_name.tm_electrons * project_lipid.lipid_mol_fraction)
            # bc
            out_chain_b = out_chain_b + (project_lipid.project_lipid_name.tg_electrons * project_lipid.lipid_mol_fraction)
            # bh
            out_headgroup_b = out_headgroup_b + (project_lipid.project_lipid_name.hg_electrons * project_lipid.lipid_mol_fraction)
    else:
        # bw
        water_b = (d_water.scattering_length * data.d2o_mol_fraction) + (water.scattering_length * (1 - data.d2o_mol_fraction))

        # Inner
        for project_lipid in in_project_lipids:
            # bt
            in_terminal_methyl_b = in_terminal_methyl_b + (project_lipid.project_lipid_name.tm_scattering * project_lipid.lipid_mol_fraction)

            # Get data lipid
            try:
                data_lipid = Data_Lipid.objects.get(data_lipid_name=project_lipid)
            except Data_Lipid.DoesNotExist:
                data_lipid = None

            if not data_lipid:
                # bc
                in_chain_b = in_chain_b + (project_lipid.project_lipid_name.tg_scattering * project_lipid.lipid_mol_fraction)
                # bh
                in_headgroup_b = in_headgroup_b + (project_lipid.project_lipid_name.hg_scattering * project_lipid.lipid_mol_fraction)
            else:
                # Delcare
                total_hg_adjustment = 0
                total_tg_adjustment = 0

                # Get atoms
                atoms = Data_Lipid_Atom.objects.filter(data_lipid_name=data_lipid).select_related('data_lipid_atom_name')
                
                if atoms:
                    # Find adjustments
                    for atom in atoms:
                        if atom.atom_location == 'HG':
                            total_hg_adjustment = total_hg_adjustment + (atom.data_lipid_atom_name.scattering_length_adj * atom.data_lipid_atom_ammount)

                        if atom.atom_location == 'TG':
                            total_tg_adjustment = total_tg_adjustment + (atom.data_lipid_atom_name.scattering_length_adj * atom.data_lipid_atom_ammount)

                # bc
                in_chain_b = in_chain_b +  ( (in_chain_b - total_tg_adjustment) * project_lipid.lipid_mol_fraction)
                # bh
                in_headgroup_b = in_headgroup_b + ( (in_headgroup_b - total_hg_adjustment) * project_lipid.lipid_mol_fraction)

        # Outer
        for project_lipid in in_project_lipids:
            # bt
            out_terminal_methyl_b = out_terminal_methyl_b + (project_lipid.project_lipid_name.tm_scattering * project_lipid.lipid_mol_fraction)

            # Get data lipid
            try:
                data_lipid = Data_Lipid.objects.get(data_lipid_name=project_lipid)
            except Data_Lipid.DoesNotExist:
                data_lipid = None

            if not data_lipid:
                # bc
                out_chain_b = out_chain_b + (project_lipid.project_lipid_name.tg_scattering * project_lipid.lipid_mol_fraction)
                # bh
                out_headgroup_b = out_headgroup_b + (project_lipid.project_lipid_name.hg_scattering * project_lipid.lipid_mol_fraction)
            else:
                # Delcare
                total_hg_adjustment = 0
                total_tg_adjustment = 0

                # Get atoms
                atoms = Data_Lipid_Atom.objects.filter(data_lipid_name=data_lipid).select_related('data_lipid_atom_name')
                
                if atoms:
                    # Find adjustments
                    for atom in atoms:
                        if atom.atom_location == 'HG':
                            total_hg_adjustment = total_hg_adjustment + (atom.data_lipid_atom_name.scattering_length_adj * atom.data_lipid_atom_ammount)

                        if atom.atom_location == 'TG':
                            total_tg_adjustment = total_tg_adjustment + (atom.data_lipid_atom_name.scattering_length_adj * atom.data_lipid_atom_ammount)

                # bc
                out_chain_b = out_chain_b +  ( (out_chain_b - total_tg_adjustment) * project_lipid.lipid_mol_fraction)
                # bh
                out_headgroup_b = out_headgroup_b + ( (out_headgroup_b - total_hg_adjustment) * project_lipid.lipid_mol_fraction)

    b_values = [in_chain_b, in_headgroup_b, in_terminal_methyl_b, out_chain_b, out_headgroup_b, out_terminal_methyl_b, water_b, calculated_water_volume]
    
    return(b_values)

# Fit function
def asymmetrical_paramitize(parameter, in_project_lipids, out_project_lipids, datas, temp):
    ## DELCARE
    # Other molecules
    water = Molecule.objects.get(compound_name='water')
    d_water = Molecule.objects.get(compound_name='deuterated_water')
    
    # Parameters
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
    try:
        for data in datas:
            # Get values for this data
            b_values = adjust_b_values(data, in_project_lipids, out_project_lipids, water, d_water, temp)

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
    except TypeError:
        # Get values for this data
        b_values = adjust_b_values(datas, in_project_lipids, out_project_lipids, water, d_water, temp)

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

def asymmetrical_graph(parameter, in_project_lipids, out_project_lipids, data, temp):
    # Get parameters
    fit_parameters = asymmetrical_paramitize(parameter, in_project_lipids, out_project_lipids, data, temp)

    # Get result
    model_result = calc_asym_model(
        fit_parameters,
        data.q_value[data.min_index:data.max_index],
        data,
        parameter.separated
    )

    return model_result

def asymmetrical_fit(parameter, in_project_lipids, out_project_lipids, datas, temp):
    # Get parameters
    fit_parameters = asymmetrical_paramitize(parameter, in_project_lipids, out_project_lipids, datas, temp)

    # Get result
    x = None

    fit_result = lsq.minimize(
        asymmetrical_objective_function,
        fit_parameters,
        args=(x, datas, parameter.separated)
    )

    return fit_result
