# create_cigale_cat.py: Script to create a file with the fluxes and uncertainties, in the format needed for CIGALE.
# This script also performs a MW extinction correction on the fluxes, and calculates the TIR luminosity.
# 27-05-2019
# Updated on 13-07-2022, to work on aenerys
# Marjorie Decleir

# Import the necessary packages
import numpy as np
import pandas as pd

from astropy.constants import c
from astropy.table import Table


def ext_corr(flux_table, filters, ext_factor, galaxy):
    """
    Function to correct the fluxes and uncertainties for the MW foreground extinction
    """
    EBV = flux_table.loc[galaxy, "E(B-V)"]
    for i, filter in enumerate(filters):
        # Correct the flux and uncertainty for the MW extinction
        A_lambda = ext_factor[i] * EBV * 3.1
        flux_table.loc[galaxy, filter] *= 10 ** (0.4 * A_lambda)
        flux_table.loc[galaxy, filter + "_err"] *= 10 ** (0.4 * A_lambda)

        # Convert the units from Jy to mJy
        flux_table.loc[galaxy, filter] *= 10 ** 3.0
        flux_table.loc[galaxy, filter + "_err"] *= 10 ** 3.0


def optical(flux_table, galaxy):
    """
    Function to exclude the BVRI fluxes if SDSS fluxes are available
    """
    BVRI = ["B", "V", "R", "I", "B_err", "V_err", "R_err", "I_err"]
    if not np.isnan(flux_table.loc[galaxy, "u"]):
        for f in BVRI:
            flux_table.loc[galaxy, f] = "NaN"


def TIR(flux_table, filters, galaxy):
    """
    Function to calculate the TIR luminosity based on the Galametz+2013 calibrations (table 3, set 2)
    """
    # Define the filters to be used in the calculation, the wavelengths of these filters (from Dale+2017), and the coefficients used in the calculation
    filters = ["MIPS24", "PACS70", "PACS100", "PACS160", "SPIRE250"]
    wavelengths = [23.8e-6, 71.8e-6, 103.0e-6, 167.0e-6, 252.0e-6]
    coeffs = np.array([2.023, 0.523, 0.390, 0.577, 0.721])
    coeffs_err = np.array([0.082, 0.028, 0.021, 0.036, 0.070])
    L = np.zeros(5)
    L_err = np.zeros(5)

    # For NGC0584, use the formulation with only 4 bands (exclude the upper limit for SPIRE 250)
    if galaxy == "NGC0584":
        filters = ["MIPS24", "PACS70", "PACS100", "PACS160"]
        wavelengths = [23.8e-6, 71.8e-6, 103.0e-6, 167.0e-6]
        coeffs = np.array([2.064, 0.539, 0.277, 0.938])
        coeffs_err = np.array([0.091, 0.030, 0.017, 0.012])
        L = np.zeros(4)
        L_err = np.zeros(4)

    # Convert the IR flux densities and their uncertainties (in Jy) to solar luminosities (L_sol)
    # 1 Jy = 1e-26 * 4 * pi * (D[m])^2 * c / lambda / 3.828e26 L_sol
    for i, filter in enumerate(filters):
        L[i] = (
            flux_table.loc[galaxy, filter]
            * 1e-26
            * 4
            * np.pi
            * (flux_table.loc[galaxy, "D(Mpc)"] * 1e6 * 3.08567758128e16) ** 2
            * c.value
            / wavelengths[i]
            / 3.828e26
        )
        L_err[i] = (
            flux_table.loc[galaxy, filter + "_err"]
            * 1e-26
            * 4
            * np.pi
            * (flux_table.loc[galaxy, "D(Mpc)"] * 1e6 * 3.08567758128e16) ** 2
            * c.value
            / wavelengths[i]
            / 3.828e26
        )

    # Calculate the TIR luminosity (in L_sol)
    LTIR = np.sum(coeffs * L)
    LTIR_err = np.sqrt(np.sum(L ** 2 * coeffs_err ** 2 + coeffs ** 2 * L_err ** 2))

    # For galaxies that have an upper limit for LTIR (i.e. the uncertainty is nan), take -LTIR as the uncertainty so that CIGALE will treat it as an upper limit
    if np.isnan(LTIR_err):
        LTIR_err = -LTIR

    # Convert the units of the TIR luminosity and its uncertainty from L_sol to W, and save the luminosity and uncertainty in the table
    flux_table.loc[galaxy, "TIR(W)"] = LTIR * 3.828e26
    flux_table.loc[galaxy, "TIR_err(W)"] = LTIR_err * 3.828e26


def write_line(frame, flux_table, galaxy):
    """
    Function to write all fluxes, uncertainties and other input parameters of a certain galaxy into the data frame
    """
    frame.loc[galaxy] = [
        galaxy,
        0.0,
        flux_table.loc[galaxy, "FUV"],
        flux_table.loc[galaxy, "UVW2"],
        flux_table.loc[galaxy, "UVM2"],
        flux_table.loc[galaxy, "NUV"],
        flux_table.loc[galaxy, "UVW1"],
        flux_table.loc[galaxy, "u"],
        flux_table.loc[galaxy, "B"],
        flux_table.loc[galaxy, "g"],
        flux_table.loc[galaxy, "V"],
        flux_table.loc[galaxy, "r"],
        flux_table.loc[galaxy, "R"],
        flux_table.loc[galaxy, "i"],
        flux_table.loc[galaxy, "I"],
        flux_table.loc[galaxy, "z"],
        flux_table.loc[galaxy, "J"],
        flux_table.loc[galaxy, "H"],
        flux_table.loc[galaxy, "K"],
        flux_table.loc[galaxy, "I1"],
        flux_table.loc[galaxy, "I2"],
        flux_table.loc[galaxy, "FUV_err"],
        flux_table.loc[galaxy, "UVW2_err"],
        flux_table.loc[galaxy, "UVM2_err"],
        flux_table.loc[galaxy, "NUV_err"],
        flux_table.loc[galaxy, "UVW1_err"],
        flux_table.loc[galaxy, "u_err"],
        flux_table.loc[galaxy, "B_err"],
        flux_table.loc[galaxy, "g_err"],
        flux_table.loc[galaxy, "V_err"],
        flux_table.loc[galaxy, "r_err"],
        flux_table.loc[galaxy, "R_err"],
        flux_table.loc[galaxy, "i_err"],
        flux_table.loc[galaxy, "I_err"],
        flux_table.loc[galaxy, "z_err"],
        flux_table.loc[galaxy, "J_err"],
        flux_table.loc[galaxy, "H_err"],
        flux_table.loc[galaxy, "K_err"],
        flux_table.loc[galaxy, "I1_err"],
        flux_table.loc[galaxy, "I2_err"],
        flux_table.loc[galaxy, "TIR(W)"],
        flux_table.loc[galaxy, "TIR_err(W)"],
        flux_table.loc[galaxy, "D(Mpc)"],
    ]


def write_line_without_SWIFT(frame, flux_table, galaxy):
    """
    Function to write all fluxes, uncertainties and other input parameters of a certain galaxy into the data frame, but without Swift values
    """
    frame.loc[galaxy] = [
        galaxy,
        0.0,
        flux_table.loc[galaxy, "FUV"],
        "NaN",
        "NaN",
        flux_table.loc[galaxy, "NUV"],
        "NaN",
        flux_table.loc[galaxy, "u"],
        flux_table.loc[galaxy, "B"],
        flux_table.loc[galaxy, "g"],
        flux_table.loc[galaxy, "V"],
        flux_table.loc[galaxy, "r"],
        flux_table.loc[galaxy, "R"],
        flux_table.loc[galaxy, "i"],
        flux_table.loc[galaxy, "I"],
        flux_table.loc[galaxy, "z"],
        flux_table.loc[galaxy, "J"],
        flux_table.loc[galaxy, "H"],
        flux_table.loc[galaxy, "K"],
        flux_table.loc[galaxy, "I1"],
        flux_table.loc[galaxy, "I2"],
        flux_table.loc[galaxy, "FUV_err"],
        "NaN",
        "NaN",
        flux_table.loc[galaxy, "NUV_err"],
        "NaN",
        flux_table.loc[galaxy, "u_err"],
        flux_table.loc[galaxy, "B_err"],
        flux_table.loc[galaxy, "g_err"],
        flux_table.loc[galaxy, "V_err"],
        flux_table.loc[galaxy, "r_err"],
        flux_table.loc[galaxy, "R_err"],
        flux_table.loc[galaxy, "i_err"],
        flux_table.loc[galaxy, "I_err"],
        flux_table.loc[galaxy, "z_err"],
        flux_table.loc[galaxy, "J_err"],
        flux_table.loc[galaxy, "H_err"],
        flux_table.loc[galaxy, "K_err"],
        flux_table.loc[galaxy, "I1_err"],
        flux_table.loc[galaxy, "I2_err"],
        flux_table.loc[galaxy, "TIR(W)"],
        flux_table.loc[galaxy, "TIR_err(W)"],
        flux_table.loc[galaxy, "D(Mpc)"],
    ]


def main():
    # Define the data path and the output path
    data_path = "/Users/mdecleir/Documents/DustKING/Data/"
    out_path = "/Users/mdecleir/Documents/DustKING/Cigale_fitting/"

    # Read the Excel file with the fluxes, uncertainties, and other input parameters
    flux_table = pd.read_excel(
        data_path + "KINGFISH_fluxes.xlsx",
        sheet_name="fluxes_new",
        index_col=0,
    )

    # List the filters
    filters = [
        "FUV",
        "UVW2",
        "UVM2",
        "NUV",
        "UVW1",
        "u",
        "B",
        "g",
        "V",
        "r",
        "R",
        "i",
        "I",
        "z",
        "J",
        "H",
        "K",
        "I1",
        "I2",
    ]

    # List the MW extinction correction factors A_lambda/A_V for these filters
    ext_factor = [
        2.5119958270240605,
        2.6456038211331627,
        2.7802372893255503,
        2.5732355432191647,
        2.1981597308149308,
        1.5670190297253372,
        1.3552756971936588,
        1.1730679870240717,
        1.023940272798787,
        0.8505686494605011,
        0.8099883501761922,
        0.6353203721851299,
        0.5934366673386918,
        0.4703845665725045,
        0.2858477781819144,
        0.1806693322493994,
        0.11673558516750641,
        0.05206923563335898,
        0.035659195499728694,
    ]

    # Create a new data frame (one with Swift and one without Swift)
    frame_with = pd.DataFrame(
        columns=[
            "#id",
            "redshift",
            "FUV",
            "SWIFT_UVW2",
            "SWIFT_UVM2",
            "NUV",
            "SWIFT_UVW1",
            "u_prime",
            "B_Harris",
            "g_prime",
            "V_Harris",
            "r_prime",
            "R_Harris",
            "i_prime",
            "I_Harris",
            "z_prime",
            "J_2mass",
            "H_2mass",
            "Ks_2mass",
            "IRAC1",
            "IRAC2",
            "FUV_err",
            "SWIFT_UVW2_err",
            "SWIFT_UVM2_err",
            "NUV_err",
            "SWIFT_UVW1_err",
            "u_prime_err",
            "B_Harris_err",
            "g_prime_err",
            "V_Harris_err",
            "r_prime_err",
            "R_Harris_err",
            "i_prime_err",
            "I_Harris_err",
            "z_prime_err",
            "J_2mass_err",
            "H_2mass_err",
            "Ks_2mass_err",
            "IRAC1_err",
            "IRAC2_err",
            "dust.luminosity",
            "dust.luminosity_err",
            "distance",
        ]
    )
    frame_without = pd.DataFrame(
        columns=[
            "#id",
            "redshift",
            "FUV",
            "SWIFT_UVW2",
            "SWIFT_UVM2",
            "NUV",
            "SWIFT_UVW1",
            "u_prime",
            "B_Harris",
            "g_prime",
            "V_Harris",
            "r_prime",
            "R_Harris",
            "i_prime",
            "I_Harris",
            "z_prime",
            "J_2mass",
            "H_2mass",
            "Ks_2mass",
            "IRAC1",
            "IRAC2",
            "FUV_err",
            "SWIFT_UVW2_err",
            "SWIFT_UVM2_err",
            "NUV_err",
            "SWIFT_UVW1_err",
            "u_prime_err",
            "B_Harris_err",
            "g_prime_err",
            "V_Harris_err",
            "r_prime_err",
            "R_Harris_err",
            "i_prime_err",
            "I_Harris_err",
            "z_prime_err",
            "J_2mass_err",
            "H_2mass_err",
            "Ks_2mass_err",
            "IRAC1_err",
            "IRAC2_err",
            "dust.luminosity",
            "dust.luminosity_err",
            "distance",
        ]
    )

    for galaxy in flux_table.index:
        # Apply the MW extinction correction
        ext_corr(flux_table, filters, ext_factor, galaxy)

        # Avoid overweight in the optical bands
        optical(flux_table, galaxy)
        # Calculate the TIR luminosity
        TIR(flux_table, filters, galaxy)

        # Add the information to the table
        write_line(frame_with, flux_table, galaxy)
        write_line_without_SWIFT(frame_without, flux_table, galaxy)

    Table.from_pandas(frame_with).write(
        out_path + "DustKING.cat", format="ascii", overwrite=True
    )
    Table.from_pandas(frame_without).write(
        out_path + "DustKING_noSWIFT.cat", format="ascii", overwrite=True
    )


if __name__ == "__main__":
    main()
