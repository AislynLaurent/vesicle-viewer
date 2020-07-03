# Other libraries
import math as mt

# Headgroup probability function
def head(parameters, data):
    ## Delcare
    headgroup_probabilities = []

    ## Unpack parameters
    # Shared
    Vc = parameters.chain_volume
    Vh = parameters.headgroup_volume
    # Unknown
    Al = parameters.lipid_area
    Dh = parameters.headgroup_thickness
    # Smearing factor
    sig = parameters.sigma

    for z in data:
        calculated_value = (
            (Vh/(2*Al*Dh)) * (
                mt.erf(
                    (Vc+Al*(Dh-z)) / (mt.sqrt(2)*Al*sig)
                ) + mt.erf(
                    (Vc+Al*(Dh+z)) / (mt.sqrt(2)*Al*sig)
                ) - mt.erf(
                    (Vc-Al*z) / (mt.sqrt(2)*Al*sig)
                ) - mt.erf(
                    (Vc+Al*z) / (mt.sqrt(2)*Al*sig)
                )
            )
        )

        headgroup_probabilities.append(calculated_value)

    return headgroup_probabilities

# Total hydrocarbon chain probability function
def chain(parameters, data):
    ## Delcare
    chain_probabilities = []

    ## Unpack parameters
    # Shared
    Vc = parameters.chain_volume
    # Unknown
    Al = parameters.lipid_area
    # Smearing factor
    sig = parameters.sigma

    for z in data.q_value[data.min_index:data.max_index]:
        calculated_value = (
            (1/2) * (
                mt.erf(
                    (Vc-Al*z) / (mt.sqrt(2)*Al*sig)
                ) + mt.erf(
                    (Vc+Al*z) / (mt.sqrt(2)*Al*sig)
                )
            )
        )

        chain_probabilities.append(calculated_value)

    return chain_probabilities

# Terminal methyl probability function
def terminal(parameters, data):
    ## Delcare
    tm_probabilities = []

    ## Unpack parameters
    # Shared
    Vt = parameters.terminal_methyl_volume
    # Unknown
    Al = parameters.lipid_area
    # Smearing factor
    sig = parameters.sigma

    for z in data.q_value[data.min_index:data.max_index]:
        calculated_value = (
            (1/2) * (
                mt.erf(
                    (2*Vt-Al*z) / (mt.sqrt(2)*Al*sig)
                ) + mt.erf(
                    (2*Vt+Al*z) / (mt.sqrt(2)*Al*sig)
                )
            )
        )

        tm_probabilities.append(calculated_value)

    return tm_probabilities

# Methylene probability function
def methylene(parameters, data):
    ## Delcare
    methylene_probabilities = []

    ## Unpack parameters
    # Shared
    Vc = parameters.chain_volume
    Vt = parameters.terminal_methyl_volume
    # Unknown
    Al = parameters.lipid_area
    # Smearing factor
    sig = parameters.sigma

    for z in data.q_value[data.min_index:data.max_index]:
        calculated_value = (
            (1/2) * (
                mt.erf(
                    (Vc-Al*z) / (mt.sqrt(2)*Al*sig)
                ) +
                mt.erf(
                    (Vc+Al*z) / (mt.sqrt(2)*Al*sig)
                ) +
                mt.erfc(
                    (2*Vt-Al*z) / (mt.sqrt(2)*Al*sig)
                ) +
                mt.erfc(
                    (2*Vt+Al*z) / (mt.sqrt(2)*Al*sig)
                ) -2
            )
        )

        methylene_probabilities.append(calculated_value)

    return methylene_probabilities

# Water probability function
def water(parameters, data):
    ## Delcare
    water_probabilities = []

    ## Unpack parameters
    # Shared
    Vc = parameters.chain_volume
    Vh = parameters.headgroup_volume
    # Unknown
    Al = parameters.lipid_area
    Dh = parameters.headgroup_thickness
    # Smearing factor
    sig = parameters.sigma

    for z in data.q_value[data.min_index:data.max_index]:
        calculated_value = (
            (Vh / (2*Al*Dh)) * (
                mt.erfc(
                    (Vc+Al*(Dh-z)) / (mt.sqrt(2)*Al*sig)
                ) + mt.erfc(
                    (Vc+Al*(Dh+z)) / (mt.sqrt(2)*Al*sig)
                )
            ) + ((1/2) - Vh / (2*Al*Dh)) * (
                mt.erfc(
                    (Vc-Al*z) / (mt.sqrt(2)*Al*sig)
                ) + mt.erfc(
                    (Vc+Al*z)/(mt.sqrt(2)*Al*sig)
                )
            )
        )

        water_probabilities.append(calculated_value)

    return water_probabilities