
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.ticker import MaxNLocator


def plot_hrd(data, members):

    x, y, c = data["TEFF"], data["LOGG"], data["FEH"]
    x_err, y_err = data["E_TEFF"], data["E_LOGG"]

    member_scatter_kwds = {
        "edgecolor": "#000000",
        "linewidths": 2,
        "s": 50,
        "zorder": 2
    }

    uves = np.array(["U580" in _ for _ in data["SETUP"]])
    giraffe = ~uves

    fig, ax = plt.subplots()
    scat = ax.scatter(x[members * uves], y[members * uves], c=c[members * uves],
        marker="s", label="UVES", **member_scatter_kwds)
    
    scat = ax.scatter(x[members * giraffe], y[members * giraffe],
        c=c[members * giraffe], marker="o", label="GIRAFFE",
        **member_scatter_kwds)
    
    ax.errorbar(x[members], y[members],
        xerr=x_err[members], yerr=y_err[members],
        fmt=None, ecolor="#000000", zorder=-1, elinewidth=1.5)

    #ax.legend(loc="upper left", frameon=False)

    cbar = plt.colorbar(scat)
    cbar.set_label(r"$[{\rm Fe/H}]$")

    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylim(ax.get_ylim()[::-1])

    ax.set_xlabel(r"$T_{\rm eff}$ $(K)$")
    ax.set_ylabel(r"$\log{g}$")

    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))

    ax.set(adjustable='box-forced',
        aspect=np.ptp(ax.get_xlim())/np.ptp(ax.get_ylim()))

    fig.tight_layout()

    # Load isochrone?
    """
    isochrone = Table.read("basti.isochrone", format="ascii",
        names=(
            "(M/Mo)in", "(M/Mo)", "log(L/Lo)", "logTe", "Mv", 
            "(U-B)", "(B-V)", "(V-I)", "(V-R)", "(V-J)", 
            "(V-K)", "(V-L)", "(H-K)"))
    """

    return fig




def plot_selection(data, figsize=None):

    candidates = data["GES_FLD"] == "NGC2808"
    members = (data["FEH"] < -0.7) * (data["VRAD"] > 80) * (data["VRAD"] < 130)\
            * candidates


    fig, axes = plt.subplots(2, figsize=figsize or (5.3, 8.3))

    x, y = data["VRAD"], data["FEH"]
    x_err, y_err = data["E_VRAD"], data["E_FEH"]

    candidate_scatter_kwds = {
        "facecolor": "#FFFFFF",
        "edgecolor": "#666666",
        "s": 30,
        "zorder": -1
    }
    member_scatter_kwds = {
        "facecolor": "#FFFFFF",
        "edgecolor": "#000000",
        "linewidths": 2,
        "s": 30,
        "zorder": 2
    }

    # Show candidates.
    axes[0].errorbar(x[candidates], y[candidates], fmt=None,
        xerr=x_err[candidates], yerr=y_err[candidates],
        ecolor="#666666", zorder=-2)
    axes[0].scatter(x[candidates], y[candidates],
        label=r"Candidates", **candidate_scatter_kwds)

    # Show members.
    axes[0].errorbar(x[members], y[members], fmt=None,
        xerr=x_err[members], yerr=y_err[members],
        ecolor="#000000", zorder=1, elinewidth=2)
    axes[0].scatter(x[members], y[members],
        label=r"Members", **member_scatter_kwds)

    # Show Caretta (2015) region?
    # http://adsabs.harvard.edu/abs/2015ApJ...810..148C
    caretta_2015_feh = (-1.129, (0.005**2 + 0.034**2 + 0.03**2)**0.5)
    axes[0].axhspan(
        caretta_2015_feh[0] - caretta_2015_feh[1],
        caretta_2015_feh[0] + caretta_2015_feh[1],
        facecolor="blue", alpha=0.5, zorder=-10,
        edgecolor="None", label=r"Caretta (2015)")

    # Show Harris line?
    # http://www.physics.mcmaster.ca/~harris/mwgc.dat
    harris_feh = -1.14
    axes[0].axhline(harris_feh, c="b", lw=2, alpha=0.5, zorder=-5,
        label=r"Harris (1996)")

    # Axes, ranges, etc.
    axes[0].legend(loc="upper right", fontsize=10, frameon=False)
    axes[0].set_xlabel(r"$V_{\rm rad}$ $({\rm km s}^{-1})$")
    axes[0].set_ylabel(r"$[{\rm Fe/H}]$ $({\rm dex})$")

    axes[0].xaxis.set_major_locator(MaxNLocator(4))
    axes[0].yaxis.set_major_locator(MaxNLocator(4))



    # Show their on-sky position.
    x, y = data["RA"], data["DEC"]
    axes[1].scatter(x[candidates], y[candidates], **candidate_scatter_kwds)
    axes[1].scatter(x[members], y[members], **member_scatter_kwds)

    # Axes, ranges, etc.
    axes[1].set_xlabel(r"$\alpha$ $(^\circ)$")
    axes[1].set_ylabel(r"$\delta$ $(^\circ)$")

    axes[1].xaxis.set_major_locator(MaxNLocator(4))
    axes[1].yaxis.set_major_locator(MaxNLocator(4))

    axes[1].set_xticklabels(
        [r"${0:.1f}$".format(_) for _ in axes[1].get_xticks()])
    axes[1].set_yticklabels(
        [r"${0:.1f}$".format(_) for _ in axes[1].get_yticks()])


    # Set equal limits.
    ptp = np.max([
        np.ptp(axes[1].get_xlim()),
        np.ptp(axes[1].get_ylim())
    ])
    axes[1].set_xlim(
        np.mean(axes[1].get_xlim()) - 0.5 * ptp,
        np.mean(axes[1].get_xlim()) + 0.5 * ptp)
    axes[1].set_ylim(
        np.mean(axes[1].get_ylim()) - 0.5 * ptp,
        np.mean(axes[1].get_ylim()) + 0.5 * ptp)


    # Set equal aspect.
    axes[0].set(adjustable='box-forced', aspect=np.ptp(axes[0].get_xlim())/np.ptp(axes[0].get_ylim()))
    axes[1].set(adjustable='box-forced', aspect=np.ptp(axes[1].get_xlim())/np.ptp(axes[1].get_ylim()))

    fig.tight_layout()


    return (fig, candidates, members)

if __name__ == "__main__":

    image = fits.open("GES_iDR4_WG15_Recommended_Abundances_20042016.fits")

    fig_selection, candidates, members = plot_selection(image[1].data)
    fig_selection.savefig("members.pdf", dpi=300)
    

    fig_hrd = plot_hrd(image[1].data, members)
    fig_hrd.savefig("hrd.pdf", dpi=300)

