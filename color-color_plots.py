# color-color_plots.py: Script to create UV color-color plots.
# 14-05-2020
# Marjorie Decleir
# Updated on 19-08-2022 to be more readable

# Import the necessary packages
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from matplotlib.lines import Line2D


def plot_data(ax, x1, x2, y1, y2, x1_err, x2_err, y1_err, y2_err):
    """
    Function to plot the data points
    """
    # Calculate the uncertainties on the colors
    err_x = 2.5 / np.log(10) * np.sqrt((x1_err / x1) ** 2 + (x2_err / x2) ** 2)
    err_y = 2.5 / np.log(10) * np.sqrt((y1_err / y1) ** 2 + (y2_err / y2) ** 2)

    # Plot the colors
    ax.errorbar(
        -2.5 * np.log10(x1 / x2),
        -2.5 * np.log10(y1 / y2),
        xerr=err_x,
        yerr=err_y,
        fmt=".",
        alpha=0.7,
        c="maroon",
        ms=3,
        elinewidth=0.2,
    )


def plot_curves(ax, curves, xcol, ycol):
    """
    Function to plot the theoretical curves. xcol and ycol give the names of the columns in the .dat files that need to be plotted.
    """
    colors = ["dodgerblue", "lime", "magenta", "black"]
    linestyles = ["--", "-", ":"]
    lines = []
    AVs = [0, 0.5, 1, 1.5, 2]
    markers = ["o", "^", "D", "s", "*"]

    # Plot the theoretical curves
    for i, curve in enumerate(curves):
        lines.append(
            ax.plot(
                curve[xcol], curve[ycol], c=colors[i % 4], ls=linestyles[i // 4], lw=1
            )[0]
        )
        # Put markers on the curves indicating the A(V)s
        for j, AV in enumerate(AVs):
            row_mask = curve["A_V"] == AV
            ax.scatter(
                curve[xcol][row_mask],
                curve[ycol][row_mask],
                c=colors[i % 4],
                marker=markers[j % 5],
                edgecolor="k",
                lw=0.15,
                s=15,
                zorder=3,
            )

    return lines


def main():
    # Define the paths
    inpath = "/Users/mdecleir/Documents/DustKING/Cigale_fitting/Aug2022/"
    outpath = "/Users/mdecleir/Documents/DustKING/Cigale_fitting/Aug2022/"

    # Define plot properties
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = "True"
    plt.rcParams["ytick.right"] = "True"
    size = 15

    # Create the plot
    fig, ax = plt.subplots(3, 3, figsize=(12, 9), sharex="col", sharey="row")

    # Obtain the catalog with data
    catalog = Table.read(inpath + "DustKING.cat", format="ascii")

    # Exclude galaxies
    row_mask = np.full_like(catalog, False, dtype=bool)
    exclude = ["IC0342"]
    for galaxy in exclude:
        row_mask += catalog["id"] == galaxy
    catalog = catalog[~row_mask]

    # Obtain the data
    uw2 = catalog["SWIFT_UVW2"]
    um2 = catalog["SWIFT_UVM2"]
    uw1 = catalog["SWIFT_UVW1"]
    FUV = catalog["FUV"]
    NUV = catalog["NUV"]
    uw2_err = catalog["SWIFT_UVW2_err"]
    um2_err = catalog["SWIFT_UVM2_err"]
    uw1_err = catalog["SWIFT_UVW1_err"]
    FUV_err = catalog["FUV_err"]
    NUV_err = catalog["NUV_err"]

    # Get the theoretical curves
    curve_path = "/Users/mdecleir/Documents/DustKING/Data/SSPcurves/"
    cal_10 = Table.read(curve_path + "SSP_colors_Cal_10.dat", format="ascii")
    cal_100 = Table.read(curve_path + "SSP_colors_Cal_102.dat", format="ascii")
    cal_300 = Table.read(curve_path + "SSP_colors_Cal_286.dat", format="ascii")
    cal = Table.read(curve_path + "colors_avSED_Cal.dat", format="ascii")
    mw = Table.read(curve_path + "colors_avSED_MW.dat", format="ascii")
    mw_10 = Table.read(curve_path + "SSP_colors_MW_10.dat", format="ascii")
    mw_100 = Table.read(curve_path + "SSP_colors_MW_102.dat", format="ascii")
    mw_300 = Table.read(curve_path + "SSP_colors_MW_286.dat", format="ascii")
    smc_10 = Table.read(curve_path + "SSP_colors_SMC_10.dat", format="ascii")
    smc_100 = Table.read(curve_path + "SSP_colors_SMC_102.dat", format="ascii")
    smc_300 = Table.read(curve_path + "SSP_colors_SMC_286.dat", format="ascii")
    smc = Table.read(curve_path + "colors_avSED_SMC.dat", format="ascii")
    curve_tables = [
        cal_10,
        cal_100,
        cal_300,
        cal,
        mw_10,
        mw_100,
        mw_300,
        mw,
        smc_10,
        smc_100,
        smc_300,
        smc,
    ]

    # Plot W2-M2 vs W2-W1
    plot_data(ax[0, 0], uw2, uw1, uw2, um2, uw2_err, uw1_err, uw2_err, um2_err)
    lines = plot_curves(ax[0, 0], curve_tables, "UVW2-UVW1", "UVW2-UVM2")
    # Plot M2-W1 vs W2-W1
    plot_data(ax[1, 0], uw2, uw1, um2, uw1, uw2_err, uw1_err, um2_err, uw1_err)
    plot_curves(ax[1, 0], curve_tables, "UVW2-UVW1", "UVM2-UVW1")
    # Plot M2-W1 vs W2-M2
    plot_data(ax[1, 1], uw2, um2, um2, uw1, uw2_err, um2_err, um2_err, uw1_err)
    plot_curves(ax[1, 1], curve_tables, "UVW2-UVM2", "UVM2-UVW1")
    # Plot FUV-NUV vs W2-W1
    plot_data(ax[2, 0], uw2, uw1, FUV, NUV, uw2_err, uw1_err, FUV_err, NUV_err)
    plot_curves(ax[2, 0], curve_tables, "UVW2-UVW1", "FUV-NUV")
    # Plot FUV-NUV vs W2-M2
    plot_data(ax[2, 1], uw2, um2, FUV, NUV, uw2_err, um2_err, FUV_err, NUV_err)
    plot_curves(ax[2, 1], curve_tables, "UVW2-UVM2", "FUV-NUV")
    # Plot FUV-NUV vs M2-W1
    plot_data(ax[2, 2], um2, uw1, FUV, NUV, um2_err, uw1_err, FUV_err, NUV_err)
    plot_curves(ax[2, 2], curve_tables, "UVM2-UVW1", "FUV-NUV")

    # Finish the figure
    ax[0, 0].set_ylabel("UVW2 - UVM2 \n (1991 $\AA$ - 2221 $\AA$)", fontsize=size)
    ax[1, 0].set_ylabel("UVM2 - UVW1\n (2221 $\AA$ - 2486 $\AA$)", fontsize=size)
    ax[2, 0].set_xlabel("UVW2 - UVW1\n (1991 $\AA$ - 2486 $\AA$)", fontsize=size)
    ax[2, 0].set_ylabel("FUV - NUV\n (1528 $\AA$ - 2271 $\AA$)", fontsize=size)
    ax[2, 1].set_xlabel("UVW2 - UVM2\n (1991 $\AA$ - 2221 $\AA$)", fontsize=size)
    ax[2, 2].set_xlabel("UVM2 - UVW1\n (2221 $\AA$ - 2486 $\AA$)", fontsize=size)
    ax[0, 1].axis("off")
    ax[0, 2].axis("off")
    ax[1, 2].axis("off")
    ax[0, 0].set_xlim(-0.2, 1.6)
    ax[0, 0].set_ylim(-0.8, 0.4)
    ax[1, 0].set_ylim(-0.2, 1.8)
    ax[2, 0].set_ylim(-0.4, 2.4)
    ax[2, 1].set_xlim(-0.8, 0.25)
    ax[2, 2].set_xlim(-0.2, 1.75)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    # Create a legend
    legend_curves = plt.figlegend(
        lines,
        (
            "Cal 10 Myr",
            "Cal 100 Myr",
            "Cal 300 Myr",
            "Cal cst SFR",
            "MW 10 Myr",
            "MW 100 Myr",
            "MW 300 Myr",
            "MW cst SFR",
            "SMC 10 Myr",
            "SMC 100 Myr",
            "SMC 300 Myr",
            "SMC cst SFR",
        ),
        fontsize=size,
        bbox_to_anchor=(0.06, -0.05, 0.9, 1.0),
        loc=1,
        ncol=3,
        borderaxespad=0.0,
        fancybox=True,
        shadow=True,
        columnspacing=0.8,
        labelspacing=0.5,
    )

    marker_list = ["o", "^", "D", "s", "*"]
    markers = [
        Line2D(
            range(1),
            range(1),
            c="w",
            marker=item,
            ms=8,
            mfc="grey",
            mec="k",
            mew=1,
        )
        for item in marker_list
    ]
    legend_markers = plt.figlegend(
        markers,
        ("$A_V=0.0$", "$A_V=0.5$", "$A_V=1.0$", "$A_V=1.5$", "$A_V=2.0$"),
        fontsize=size,
        numpoints=1,
        bbox_to_anchor=(0.9, 0.54, -0.02, 0.0),
        loc=7,
        fancybox=True,
        shadow=True,
    )

    # Save the figure
    plt.savefig(outpath + "KF_color_plots.pdf")


if __name__ == "__main__":
    main()
