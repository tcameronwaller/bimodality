"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard
import os
import pickle

# Relevant
import numpy
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams.update({'figure.max_open_warning': 0})
import matplotlib.pyplot

# Custom

import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def define_font_properties():
    """
    Defines font properties.

    arguments:

    raises:

    returns:
        (dict<object>): references to definitions of font properties

    """

    # Define font values.
    values_one = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 1000,
        "weight": 1000,
        "size": 30
    }
    values_two = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 1000,
        "size": 25
    }
    values_three = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 1000,
        "size": 20
    }
    values_four = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 500,
        "size": 15
    }
    # Define font properties.
    properties_one = matplotlib.font_manager.FontProperties(
        family=values_one["family"],
        style=values_one["style"],
        variant=values_one["variant"],
        stretch=values_one["stretch"],
        weight=values_one["weight"],
        size=values_one["size"]
    )
    properties_two = matplotlib.font_manager.FontProperties(
        family=values_two["family"],
        style=values_two["style"],
        variant=values_two["variant"],
        stretch=values_two["stretch"],
        weight=values_two["weight"],
        size=values_two["size"]
    )
    properties_three = matplotlib.font_manager.FontProperties(
        family=values_three["family"],
        style=values_three["style"],
        variant=values_three["variant"],
        stretch=values_three["stretch"],
        weight=values_three["weight"],
        size=values_three["size"]
    )
    properties_four = matplotlib.font_manager.FontProperties(
        family=values_four["family"],
        style=values_four["style"],
        variant=values_four["variant"],
        stretch=values_four["stretch"],
        weight=values_four["weight"],
        size=values_four["size"]
    )
    # Compile and return references.
    return {
        "values": {
            "one": values_one,
            "two": values_two,
            "three": values_three,
            "four": values_four
        },
        "properties": {
            "one": properties_one,
            "two": properties_two,
            "three": properties_three,
            "four": properties_four
        }
    }


def define_color_properties():
    """
    Defines color properties.

    arguments:

    raises:

    returns:
        (dict<tuple>): references to definitions of color properties

    """

    # Black.
    black = (0.0, 0.0, 0.0, 1.0)
    # White.
    white = (1.0, 1.0, 1.0, 1.0)
    white_faint = (1.0, 1.0, 1.0, 0.75)
    # Blue.
    blue = (0.0, 0.2, 0.5, 1.0)
    blue_faint = (0.0, 0.2, 0.5, 0.75)
    # Orange.
    orange = (1.0, 0.6, 0.2, 1.0)
    orange_faint = (1.0, 0.6, 0.2, 0.75)
    # Compile and return references.
    return {
        "black": black,
        "white": white,
        "white_faint": white_faint,
        "blue": blue,
        "blue_faint": blue_faint,
        "orange": orange,
        "orange_faint": orange_faint
    }


def plot_distribution_histogram(
    series=None,
    name=None,
    bin_method=None,
    bin_count=None,
    label_bins=None,
    label_counts=None,
    fonts=None,
    colors=None,
    line=None,
    position=None,
    text=None,
):
    """
    Creates a figure of a chart of type histogram to represent the frequency
    distribution of a single series of values.

    arguments:
        series (list<float>): series of counts
        name (str): name of series
        bin_method (str): method to define bins
        bin_count (int): count of bins to define and populate
        label_bins (str): label for bins
        label_counts (str): label for counts
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        line (bool): whether to draw a vertical line
        position (float): position for vertical line
        text (str): text to include on figure

    raises:

    returns:
        (object): figure object

    """

    # Define and populate bins.
    # Bin method "auto" is useful.
    #values, edges = numpy.histogram(series, bins=count_bins)
    if bin_method == "count":
        bin_edges = numpy.histogram_bin_edges(series, bins=bin_count)
    else:
        bin_edges = numpy.histogram_bin_edges(series, bins=bin_method)

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = matplotlib.pyplot.axes()
    values, bins, patches = axes.hist(
        series,
        bins=bin_edges,
        histtype="bar",
        align="left",
        orientation="vertical",
        rwidth=0.35,
        log=False,
        color=colors["blue"],
        label=name,
        stacked=False
    )
    if False:
        axes.legend(
            loc="upper right",
            markerscale=2.5,
            markerfirst=True,
            prop=fonts["properties"]["one"],
            edgecolor=colors["black"]
        )
    axes.set_xlabel(
        xlabel=label_bins,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=label_counts,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["two"]["size"],
        labelcolor=colors["black"]
    )
    if line:
        axes.axvline(
            x=position,
            ymin=0,
            ymax=1,
            alpha=1.0,
            color=colors["orange"],
            linestyle="--",
            linewidth=3.0,
        )
    if len(text) > 0:
        matplotlib.pyplot.text(
            1,
            1,
            text,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
        )
    return figure


# Procedures


def read_source_analysis(dock=None):
    """
    Reads and organizes source information from file

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_analysis = os.path.join(dock, "analysis")
    #path_genes_report = os.path.join(
    #    path_analysis, "genes_report.txt"
    #)
    path_reports = os.path.join(
        path_analysis, "reports.pickle"
    )
    # Read information from file.
    #genes = utility.read_file_text_list(path_genes_report)
    with open(path_reports, "rb") as file_source:
        reports = pickle.load(file_source)
    # Compile and return information.
    return {
    #    "genes": genes,
        "reports": reports,
    }


def plot_charts_analysis(
    dock=None
):
    """
    Plots charts from the analysis process.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.confirm_path_directory(path_plot)
    path_analysis = os.path.join(path_plot, "analysis")
    # Remove previous files since they change from run to run.
    utility.remove_directory(path=path_analysis)
    utility.confirm_path_directory(path_analysis)

    # Read source information from file.
    source = read_source_analysis(dock=dock)
    #source["genes"]
    #source["reports"]

    # Extract genes' identifiers.
    genes = list(source["reports"].keys())
    # Iterate on genes.
    for gene in genes:
        # Access information.
        name = source["reports"][gene]["name"]
        # Specify directories and files.
        path_report = os.path.join(
            path_plot, name
        )
        utility.confirm_path_directory(path_report)

        # Create figures.
        figure_real = plot_distribution_histogram(
            series=source["reports"][gene]["abundances_real"],
            name="",
            bin_method="count",
            bin_count=50,
            label_bins="Bins",
            label_counts="Counts",
            fonts=fonts,
            colors=colors,
            line=False,
            position=0,
            text="",
        )
        # Specify directories and files.
        file = ("abundance_distribution_real.svg")
        path_file = os.path.join(path_report, file)
        # Write figure.
        write_figure(
            path=path_file,
            figure=figure_real
        )

        figure_shuffle = plot_distribution_histogram(
            series=source["reports"][gene]["abundances_shuffle"],
            name="",
            bin_method="count",
            bin_count=50,
            label_bins="Bins",
            label_counts="Counts",
            fonts=fonts,
            colors=colors,
            line=False,
            position=0,
            text="",
        )
        # Specify directories and files.
        file = ("abundance_distribution_shuffle.svg")
        path_file = os.path.join(path_report, file)
        # Write figure.
        write_figure(
            path=path_file,
            figure=figure_shuffle
        )

        for type in ["coefficient", "dip", "mixture", "combination"]:
            print(source["reports"][gene]["scores_distributions"][type]["score"])
            # Create figure.
            figure_score = plot_distribution_histogram(
                series=source["reports"][gene]["scores_distributions"][type]["distribution"],
                name="",
                bin_method="count",
                bin_count=50,
                label_bins="Bins",
                label_counts="Counts",
                fonts=fonts,
                colors=colors,
                line=True,
                position=source["reports"][gene]["scores_distributions"][type]["score"],
                text="",
            )
            # Specify directories and files.
            file = ("score_distribution_" + type + ".svg")
            path_file = os.path.join(path_report, file)
            # Write figure.
            write_figure(
                path=path_file,
                figure=figure_score
            )

    pass





def write_figure(path=None, figure=None):
    """
    Writes figure to file.

    arguments:
        path (str): path to directory and file
        figure (object): figure object

    raises:

    returns:

    """

    # Write information to file.
    figure.savefig(
        path,
        format="svg",
        #dpi=600,
        facecolor="w",
        edgecolor="w",
        transparent=False
    )
    pass





###############################################################################
# Procedure


def execute_procedure(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """


    plot_charts_analysis(dock=dock)

    pass


if (__name__ == "__main__"):
    execute_procedure()
