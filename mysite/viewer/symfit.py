# Other libraries
import numpy as np
import lmfit as lsq

# Models
from .models import Project
from .models import Project_Lipid
from .models import Symmetrical_Parameters
from .models import Data_Set
from .models import Data_Lipid
from .models import Data_Lipid_Atom
from .models import Molecule
from .models import Lipid

# Symmetrical model
def sym_model(
    q,  # independant
    Vc, # chain_volume
    Vh, # headgroup_volume
    Vt, # terminal_methyl_volume
    Vw, # water_volume
    Al, # area_per_lipid
    Dh, # headgroup_thickness
    bc, # chain_b
    bh, # headgroup_b
    bt, # terminal_methyl_b
    bw  # water_b
):
    return (
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

    # Return the calculated model
    return sym_model(q_array, Vc, Vh, Vt, Vw, Al, Dh, bc, bh, bt, bw)


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
def symmetrical_fit(project_id, parameter_id):
    # DELCARE
    # Project
    parameter = Symmetrical_Parameters.objects.get(id=parameter_id)

    # Project lipids
    project_lipids = Project_Lipid.objects.filter(project_title_id=project_id).select_related('project_lipid_name')
    
    # Other molecules
    water = Molecule.objects.get(compound_name='water')
    d_water = Molecule.objects.get(compound_name='deuterated_water')

    # Data
    datas = Data_Set.objects.filter(project_title_id=project_id)

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
        )
    )

    # Per dataset
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

    # Fit!
    x = None

    fit_result = lsq.minimize(
        symmetrical_objective_function,
        fit_parameters,
        args=(x, datas)
    )
    
    # print('report')
    # print(lsq.fit_report(fit_result))
    # print('\n')

    # print('var_names')
    # print(fit_result.var_names)
    # print('covar')
    # print(fit_result.covar)
    # print('init_vals')
    # print(fit_result.init_vals)
    # print('init_values')
    # print(fit_result.init_values)
    # print('ier')
    # print(fit_result.ier)
    # print('residual')
    # print(fit_result.residual)
    # print('chisqr')
    # print(fit_result.chisqr)
    # print('redchi')
    # print(fit_result.redchi)
    # print('aic')
    # print(fit_result.aic)
    # print('bic')
    # print(fit_result.bic)
    # print('flatchain')
    # print(fit_result.flatchain)

    return fit_result
