# compare_runs.py: Script to compare the Cigale results of a run with and without Swift data

# Import the necessary packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table


def plot_comp(quantity, ax, xvalue, xerr, yvalue, yerr):
    """
    Function to plot the comparison between two runs
    """
    # Plot the data points
    ax.scatter(xvalue, yvalue, s=15)

    # Add the 1-1 line
    line = np.linspace(
        np.min(xvalue) - np.ptp(xvalue) * 0.04,
        np.max(xvalue) + np.ptp(xvalue) * 0.04,
        1000,
    )
    ax.plot(line, line, c="g")

    # Add a median error bar
    ax.errorbar(
        np.min(xvalue) + np.median(xerr),
        np.max(yvalue) - np.median(yerr),
        xerr=np.median(xerr),
        yerr=np.median(yerr),
        ls="none",
        c="k",
        elinewidth=0.5,
        capsize=2,
        capthick=0.5,
    )

    # Finalize the plot
    ax.set_aspect("equal")
    ax.set_xlabel(quantity + " with SWIFT")
    ax.set_ylabel(quantity + " without SWIFT")


def print_comp(quantity, xvalue, xerr, yvalue, yerr):
    """
    Function to print some useful information about the results"
    """
    print(
        "Median " + quantity + " uncertainty with and without Swift:",
        np.median(xerr),
        np.median(yerr),
    )
    print(
        "Median of relative " + quantity + " uncertainty with and without Swift:",
        np.median(xerr / xvalue),
        np.median(yerr / yvalue),
    )


def compare(path, results_with, results_without, catalog):
    """
    Function to create a plot with a comparison between the results with and without Swift, and print useful results
    """
    fig, ax = plt.subplots(2, 1, figsize=(5, 9.5))

    # Check which galaxies do not have SWIFT at all, and exclude them from the comparison
    noSWIFT_mask = (
        catalog["SWIFT_UVM2"].mask
        * catalog["SWIFT_UVW2"].mask
        * catalog["SWIFT_UVW1"].mask
    )
    bumps_with = results_with["bayes.attenuation.uv_bump_amplitude"][~noSWIFT_mask]
    bumps_err_with = results_with["bayes.attenuation.uv_bump_amplitude_err"][
        ~noSWIFT_mask
    ]
    bumps_without = results_without["bayes.attenuation.uv_bump_amplitude"][
        ~noSWIFT_mask
    ]
    bumps_err_without = results_without["bayes.attenuation.uv_bump_amplitude_err"][
        ~noSWIFT_mask
    ]
    slopes_with = results_with["bayes.attenuation.powerlaw_slope"][~noSWIFT_mask]
    slopes_err_with = results_with["bayes.attenuation.powerlaw_slope_err"][
        ~noSWIFT_mask
    ]
    slopes_without = results_without["bayes.attenuation.powerlaw_slope"][~noSWIFT_mask]
    slopes_err_without = results_without["bayes.attenuation.powerlaw_slope_err"][
        ~noSWIFT_mask
    ]
    for i, gal in enumerate(results_with["id"][~noSWIFT_mask]):
        print(i, gal)

    # Compare the bumps
    plot_comp(
        "$B$",
        ax[0],
        bumps_with,
        bumps_err_with,
        bumps_without,
        bumps_err_without,
    )
    print_comp(
        "bump",
        bumps_with,
        bumps_err_with,
        bumps_without,
        bumps_err_without,
    )

    # Compare the slopes
    plot_comp(
        "$\delta$",
        ax[1],
        slopes_with,
        slopes_err_with,
        slopes_without,
        slopes_err_without,
    )

    print_comp(
        "slope",
        slopes_with,
        slopes_err_with,
        slopes_without,
        slopes_err_without,
    )

    # Finalize and save the figure
    plt.tight_layout()
    plt.savefig(path + "comp_SWIFT.pdf", bbox_inches="tight")


def main():
    # Define the path
    path = "/Users/mdecleir/Documents/DustKING/Cigale_fitting/Aug2022/"

    # Read the files with the results
    results_with = Table.read(
        path + "out_err0/results.txt",
        format="ascii",
    )
    results_without = Table.read(
        path + "out_noSWIFT/results.txt",
        format="ascii",
    )

    # Read the catalog
    catalog = Table.read(path + "DustKING.cat", format="ascii")

    # Plotting parameters.
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = "True"
    plt.rcParams["ytick.right"] = "True"
    plt.rcParams["xtick.labelsize"] = 12
    plt.rcParams["ytick.labelsize"] = 12
    plt.rcParams["axes.labelsize"] = 16

    # Compare the results
    compare(path, results_with, results_without, catalog)


if __name__ == "__main__":
    main()
