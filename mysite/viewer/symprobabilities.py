# Other libraries
import math as mt

# Headgroup probability function
def head(
    Vc,
    Vh,
    Al,
    Dh,
    sig,
    data
):
    ## Delcare
    headgroup_probabilities = []

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
def chain(
    Vc,
    Al,
    sig,
    data
):
    ## Delcare
    chain_probabilities = []

    for z in data:
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
def terminal(
    Vt,
    Al,
    sig,
    data
):
    ## Delcare
    tm_probabilities = []

    for z in data:
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
def methylene(
    Vc,
    Vt,
    Al,
    sig,
    data
):
    ## Delcare
    methylene_probabilities = []

    for z in data:
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
def water(
    Vc,
    Vh,
    Al,
    Dh,
    sig,
    data
):
    ## Delcare
    water_probabilities = []

    for z in data:
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