# Other libraries
import numpy as np
import lmfit as lsq

# Models
from .models import Data_Lipid
from .models import Data_Lipid_Atom
from .models import Molecule

# Symmetrical model
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
        q**(-2) * scale*(
            (
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
def calc_sym_model(fit_parameters, q, data):
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
    return sym_model(q_array, Vc, Vh, Vt, Vw, Al, Dh, sig, bc, bh, bt, bw, scale, bg)


# Objective function create a residual for each, then flatten
def symmetrical_objective_function(fit_parameters, x, datas):
    # Delcare
    residuals = []
    combined_residuals = []

    # Make an array of residuals
    for data in datas:
        # Do math
        residuals = data.intensity_value[data.min_index:data.max_index] - calc_sym_model(fit_parameters, data.q_value[data.min_index:data.max_index], data)

        # Append
        combined_residuals.extend(residuals)

    return combined_residuals

# Augments per data-set
def adjust_b_values(data, project_lipids, water, d_water, temp):
    # Temp
    x = temp

    # Declare
    terminal_methyl_b = 0
    chain_b = 0
    headgroup_b = 0
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

        for project_lipid in project_lipids:
            # bt
            terminal_methyl_b = terminal_methyl_b + (project_lipid.project_lipid_name.tm_electrons * project_lipid.lipid_mol_fraction)
            # bc
            chain_b = chain_b +  (project_lipid.project_lipid_name.tg_electrons * project_lipid.lipid_mol_fraction)
            # bh
            headgroup_b = headgroup_b + (project_lipid.project_lipid_name.hg_electrons * project_lipid.lipid_mol_fraction)
    else:
        # bw
        water_b = (d_water.scattering_length * data.d2o_mol_fraction) + (water.scattering_length * (1 - data.d2o_mol_fraction))

        for project_lipid in project_lipids:
            # bt
            terminal_methyl_b = terminal_methyl_b + (project_lipid.project_lipid_name.tm_scattering * project_lipid.lipid_mol_fraction)

            # Get data lipid
            try:
                data_lipid = Data_Lipid.objects.get(data_lipid_name=project_lipid)
            except Data_Lipid.DoesNotExist:
                data_lipid = None

            if not data_lipid:
                # bc
                chain_b = chain_b +  (project_lipid.project_lipid_name.tg_scattering * project_lipid.lipid_mol_fraction)
                # bh
                headgroup_b = headgroup_b + (project_lipid.project_lipid_name.hg_scattering * project_lipid.lipid_mol_fraction)
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
                chain_b = chain_b +  ( (chain_b - total_tg_adjustment) * project_lipid.lipid_mol_fraction)
                # bh
                headgroup_b = headgroup_b + ( (headgroup_b - total_hg_adjustment) * project_lipid.lipid_mol_fraction)

    b_values = [chain_b, headgroup_b, terminal_methyl_b, water_b, calculated_water_volume]
    
    return(b_values)

# Fit function
def symmetrical_paramitize(parameter, project_lipids, datas, temp):
    ## DELCARE
    # Other molecules
    water = Molecule.objects.get(compound_name='water')
    d_water = Molecule.objects.get(compound_name='deuterated_water')

    # Parameters
    fit_parameters = lsq.Parameters()
    fit_parameters.add_many(
        # Shared
        ( # Vc
            'chain_volume',
            parameter.chain_volume,
            not(parameter.chain_volume_lock),
            parameter.chain_volume_lowerbound,
            parameter.chain_volume_upperbound
        ),
        ( # Vh
            'headgroup_volume',
            parameter.headgroup_volume,
            not(parameter.headgroup_volume_lock),
            parameter.headgroup_volume_lowerbound,
            parameter.headgroup_volume_upperbound
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
        )
    )

    # Per dataset
    try:
        for data in datas:
            # Get values for this data
            b_values = adjust_b_values(data, project_lipids, water, d_water, temp)

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
    except TypeError:
        # Get values for this data
        b_values = adjust_b_values(datas, project_lipids, water, d_water, temp)

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

def symmetrical_graph(parameter, project_lipids, data, temp):
    # Get parameters
    fit_parameters = symmetrical_paramitize(parameter, project_lipids, data, temp)

    # Get result
    model_result = calc_sym_model(
        fit_parameters,
        data.q_value[data.min_index:data.max_index],
        data
    )

    return model_result

def symmetrical_fit(parameter, project_lipids, datas, temp):
    # Get parameters
    fit_parameters = symmetrical_paramitize(parameter, project_lipids, datas, temp)

    # Get result
    x = None

    fit_result = lsq.minimize(
        symmetrical_objective_function,
        fit_parameters,
        args=(x, datas)
    )

    return fit_result
