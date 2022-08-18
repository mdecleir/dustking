# Script to plot probability distribution functions
# Based on the __init__.py script in CIGALE

# Import the necessary packages
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path


def plot_pdf(outdir, ax, obj_name, var_name, x_label, ls, color):
    """
    Function to plot the pdf of one variable
    """
    fnames = outdir.glob(f"{obj_name}_{var_name}_chi2-block-*.npy")

    likelihood = []
    model_variable = []
    for fname in fnames:
        data = np.memmap(fname, dtype=np.float64)
        data = np.memmap(fname, dtype=np.float64, shape=(2, data.size // 2))

        likelihood.append(np.exp(-data[0, :] / 2.0))
        model_variable.append(data[1, :])

    if len(likelihood) > 0:
        likelihood = np.concatenate(likelihood)
        model_variable = np.concatenate(model_variable)
        w = np.where(np.isfinite(likelihood) & np.isfinite(model_variable))
        likelihood = likelihood[w]
        model_variable = model_variable[w]

        Npdf = 100
        min_hist = np.min(model_variable)
        max_hist = np.max(model_variable)
        Nhist = min(Npdf, len(np.unique(model_variable)))

        if min_hist == max_hist:
            pdf_grid = np.array([min_hist, max_hist])
            pdf_prob = np.array([1.0, 1.0])
        else:
            pdf_prob, pdf_grid = np.histogram(
                model_variable,
                Nhist,
                (min_hist, max_hist),
                weights=likelihood,
                density=True,
            )
            pdf_x = (pdf_grid[1:] + pdf_grid[:-1]) / 2.0

            pdf_grid = np.linspace(min_hist, max_hist, Npdf)
            pdf_prob = np.interp(pdf_grid, pdf_x, pdf_prob)

        ax.plot(pdf_grid, pdf_prob, ls=ls, color=color)
        ax.set_xlabel(x_label, fontsize=15)
        ax.set_ylabel("Probability density", fontsize=15)
        ax.minorticks_on()


def plot_4pdfs(path, object):
    """
    Function to plot the pdfs of 4 parameters in one figure
    """
    # Create the figure
    figure, axes = plt.subplots(2, 2, figsize=(8, 8))

    # Plot the pdfs for 4 parameters
    params = [
        "attenuation.uv_bump_amplitude",
        "attenuation.powerlaw_slope",
        "attenuation.FUV",
        "sfh.sfr",
    ]
    labels = ["$B$", "$\delta$", "$A_{FUV}$", "SFR ($M_{\odot}$/yr)"]
    for i, ax in enumerate(axes.flatten()):
        # For the run without Swift
        plot_pdf(
            path / "out_noSWIFT",
            ax,
            object,
            params[i],
            labels[i],
            ls="--",
            color="orange",
        )

        # For the run with Swift
        plot_pdf(
            path / "out_err0",
            ax,
            object,
            params[i],
            labels[i],
            ls="-",
            color="k",
        )

    # Finalize and save the figure
    axes[0, 0].set_xlim([-1, 6])
    axes[0, 1].set_xlim([-1.2, 0.1])
    axes[1, 0].set_xlim([0, 5.8])
    axes[1, 1].set_xlim([0, 12.5])
    plt.tight_layout()
    figure.savefig(path / "KFpdfs22.pdf")


def main():
    path = "/Users/mdecleir/Documents/DustKING/Cigale_fitting/Jul2022/"
    plot_4pdfs(Path(path), "NGC4254")


if __name__ == "__main__":
    main()
