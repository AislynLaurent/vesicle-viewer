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
    bc,     # chain_b
    bh,     # headgroup_b
    bt,     # terminal_methyl_b
    bw,     # water_b
    scale,  # scale
    bg      # bg
):
    return (
        ((
            scale*(
                2*(
                    Vt*(
                        bw*(
                            ((Al*Dh) - Vh) * (Vc - 2*Vt)
                        )
                        + ((Vw*bh) * (Vc - 2*Vt)) - ((Vw*Al*Dh) * (bc - 2*bt))
                    )
                    * np.sin(
                        (q*Vc) / Al
                    )
                    + Vt*(
                        (Vc - 2*Vt) * (bw*Vh - bh*Vw)
                    )
                    * np.sin(
                        q*Dh + ((q*Vc) / Al)
                    )
                    + (Vw*Al*Dh)*(bc*Vt - bt*Vc)
                    * np.sin(
                        ((2*q) * Vc) / Al
                    )
                ) 
                / (
                    (q*Dh*Al*Vt*Vw) * (Vc - 2*Vt)
                )
            )
        ) ** 2 )
        + bg
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
    Vw = fit_parameters['water_volume'].value
    # Unknown
    Al = fit_parameters['area_per_lipid'].value
    Dh = fit_parameters['headgroup_thickness'].value
    # Per dataset
    bc = fit_parameters['chain_b_%i' % data.id].value
    bh = fit_parameters['headgroup_b_%i' % data.id].value
    bt = fit_parameters['terminal_methyl_b_%i' % data.id].value
    bw = fit_parameters['water_b_%i' % data.id].value
    # Adjust
    scale = fit_parameters['scale'].value
    bg = fit_parameters['background'].value

    # Return the calculated model
    return sym_model(q_array, Vc, Vh, Vt, Vw, Al, Dh, bc, bh, bt, bw, scale, bg)


# Objective function create a residual for each, then flatten
def symmetrical_objective_function(fit_parameters, x, datas):
    # Delcare
    residuals = []
    combined_residuals = []

    # Make an array of residuals
    for data in datas:
        # Do math
        residuals = data.intensity_value - calc_sym_model(fit_parameters, data.q_value, data)

        # Append
        combined_residuals.extend(residuals)

    return combined_residuals

# Augments per data-set
def adjust_b_values(data, project_lipids, water, d_water):
    # Declare
    terminal_methyl_b = 0
    chain_b = 0
    headgroup_b = 0
    water_b = 0

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

    b_values = [chain_b, headgroup_b, terminal_methyl_b, water_b]
    
    return(b_values)

# Fit function
def symmetrical_paramitize(parameter, project_lipids, datas):
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
        ( # Vw
            'water_volume',
            water.total_volume,
            False
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
        # Tweaks
        ( # Scale
            'scale',
            parameter.scale,
            not(parameter.scale_lock),
            parameter.scale_lowerbound,
            parameter.scale_upperbound
        ),
        ( # BG
            'background',
            parameter.background,
            not(parameter.background_lock),
            parameter.background_lowerbound,
            parameter.background_upperbound
        )
    )

    # Per dataset
    try:
        for data in datas:
            # Get values for this data
            b_values = adjust_b_values(data, project_lipids, water, d_water)

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
                )
            )
    except TypeError:
        # Get values for this data
        b_values = adjust_b_values(datas, project_lipids, water, d_water)

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
            )
        )

    return fit_parameters

def symmetrical_graph(parameter, project_lipids, data):
    # Get parameters
    fit_parameters = symmetrical_paramitize(parameter, project_lipids, data)

    # Get result
    model_result = calc_sym_model(
        fit_parameters,
        data.q_value,
        data
    )

    return model_result

def symmetrical_fit(parameter, project_lipids, datas):
    # Get parameters
    fit_parameters = symmetrical_paramitize(parameter, project_lipids, datas)

    # Get result
    x = None

    fit_result = lsq.minimize(
        symmetrical_objective_function,
        fit_parameters,
        args=(x, datas)
    )

    return fit_result
