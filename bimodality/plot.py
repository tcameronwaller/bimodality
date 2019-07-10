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
import itertools
import statistics

# Relevant
import pandas
import numpy
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams.update({'figure.max_open_warning': 0})
import matplotlib.pyplot
import matplotlib.lines
import seaborn

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


def plot_bar_stack(
    data=None,
    label_vertical=None,
    label_horizontal=None,
    fonts=None,
    colors=None,
    rotation=None,
    legend=None,
):
    """
    Creates a figure of a chart of type histogram to represent the frequency
    distribution of a single series of values.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        label_vertical (str): label for vertical axis
        label_horizontal (str): label for horizontal axis
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        rotation (str): value for rotation of labels
        legend (bool): whether to include a legend for series on the chart

    raises:

    returns:
        (object): figure object

    """

    # Define groups and series.
    groups = list(data.loc[:, "groups"])
    data.set_index(
        ["groups"],
        append=False,
        drop=True,
        inplace=True
    )
    series_names = list(data.columns)

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()

    # Initialize bases.
    bases = list(itertools.repeat(0, len(groups)))
    # Initialize bars.
    bars = []
    # Iterate on series.
    for series_name in series_names:
        # Access series values.
        values = list(data.loc[:, series_name])
        print(series_name)
        print(values)
        # Create charts on axes.
        series_bars = axes.bar(
            range(len(groups)),
            values,
            width=0.8,
            bottom=bases,
            align="center",
        )
        bars.append(series_bars[0])
        # Restore bases.
        bases = list(map(sum, zip(bases, values)))
        print(bases)
    axes.set_ylabel(
        ylabel=label_vertical,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes.set_xlabel(
        xlabel=label_horizontal,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"],
        rotation="horizontal",
    )
    axes.set_xticks(range(len(groups)), minor=False)
    axes.set_xticklabels(
        groups,
        minor=False,
        rotation=rotation,
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
    if legend:
        axes.legend(
            bars,
            series_names,
            loc="upper left",
            markerscale=2.5,
            markerfirst=True,
            prop=fonts["properties"]["one"],
            edgecolor=colors["black"]
        )
    return figure


def plot_scatter_cluster(
    data=None,
    abscissa=None,
    ordinate=None,
    label_horizontal=None,
    label_vertical=None,
    factor=None,
    fonts=None,
    colors=None,
    legend=None,
):
    """
    Creates a figure of a chart of type histogram to represent the frequency
    distribution of a single series of values.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with independent variable
        ordinate (str): name of data column with dependent variable
        label_horizontal (str): label for horizontal axis
        label_vertical (str): label for vertical axis
        factor (str): name of data column with groups or factors of samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        legend (bool): whether to include a legend for series on the chart

    raises:

    returns:
        (object): figure object

    """

    ##########
    # Organize data.
    # Separate points by groups.
    # Define groups.
    data = data.copy(deep=True)
    data.set_index(
        factor,
        append=False,
        drop=True,
        inplace=True
    )
    print(data)
    groups = data.groupby(level=[factor])
    print("Count of groups by factor: " + str(len(groups)))
    colors_series = list(seaborn.color_palette("hls", n_colors=len(groups)))

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlabel(
        xlabel=label_horizontal,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=label_vertical,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
        labelcolor=colors["black"]
    )
    # Plot points for values from each group.
    index = 0
    for name, group in groups:
        values_x = group[abscissa].to_list()
        values_y = group[ordinate].to_list()
        handle = axes.plot(
            values_x,
            values_y,
            linestyle="",
            marker="o",
            markersize=2.5,
            markeredgecolor=colors_series[index],
            markerfacecolor=colors_series[index]
        )
        index += 1
        pass
    # Plot labels for each group.
    labels = []
    index = 0
    for name, group in groups:
        values_x = group[abscissa].to_list()
        mean_x = statistics.median(values_x)
        values_y = group[ordinate].to_list()
        mean_y = statistics.median(values_y)
        axes.text(
            mean_x,
            mean_y,
            str(index+1),
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["three"],
            horizontalalignment="center",
            verticalalignment="center"
        )
        label = str(index+1) + ": " + name
        labels.append(label)
        index += 1
        pass
    # Create legend.
    # Create custome elements for the legend.
    elements = create_legend_elements(
        colors=colors_series,
        labels=labels
    )
    axes.legend(
        handles=elements,
        loc="upper right",
        prop=fonts["properties"]["four"],
    )
    return figure


def create_legend_elements(
    colors=None,
    labels=None,
):
    """
    Creates custom elements for legend.

    arguments:
        colors (list<dict>): colors
        labels (str): name of data column with independent variable

    raises:

    returns:
        (list<object>): elements for legend

    """

    elements = []
    for index in range(len(labels)):
        element = matplotlib.lines.Line2D(
            [0],
            [0],
            marker="o",
            color=colors[index],
            label=labels[index],
            markerfacecolor=colors[index],
            markersize=15,
        )
        elements.append(element)
    return elements


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
        #format="png",
        format="svg",
        #dpi=600,
        facecolor="w",
        edgecolor="w",
        transparent=False
    )
    pass


# Procedures


# Analysis

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

# Sample

def read_source_sample(dock=None):
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
    path_assembly = os.path.join(dock, "assembly")
    #path_genes_report = os.path.join(
    #    path_analysis, "genes_report.txt"
    #)
    path_samples = os.path.join(
        path_assembly, "data_samples_tissues_patients.pickle"
    )
    # Read information from file.
    #genes = utility.read_file_text_list(path_genes_report)
    with open(path_samples, "rb") as file_source:
        data_samples = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_samples": data_samples,
    }


def plot_chart_sample_patients(
    data=None,
    fonts=None,
    colors=None,
    path=None,
):
    """
    Plots charts from the sample process.

    arguments:
        data (object): Pandas data frame with information about samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        path (str): path to directory

    raises:

    returns:

    """

    data_patients = data.copy(deep=True)
    data_patients.reset_index(level=["sample"], inplace=True)
    # Remove relevant redundancy.
    print(
        "Data dimensions before removal of redundancy: " +
        str(data_patients.shape)
    )
    data_patients.drop_duplicates(
        subset=["patient"],
        keep="first",
        inplace=True,
    )
    print(
        "Data dimensions after removal of redundancy: " +
        str(data_patients.shape)
    )
    # Remove irrelevant columns.
    data_patients.drop(
        labels=["sample", "tissue_major", "tissue_minor"],
        axis="columns",
        inplace=True,
    )
    data_patients.set_index(
        ["sex", "age"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data_patients)
    data_patients_counts = data_patients.groupby(
        level=["sex", "age"],
        sort=True,
        as_index=False
    ).size().to_frame(
        name="count"
    )
    data_patients_counts.reset_index(
        level=["sex", "age"], inplace=True
    )
    print(data_patients_counts)
    # Organize data for chart.
    data_pivot = data_patients_counts.pivot_table(
        values="count",
        index="sex",
        columns="age",
        aggfunc="sum",
        fill_value=0,
    )
    print(data_pivot)
    data_pivot.rename_axis(
        index="groups",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_pivot.reset_index(
        level=["groups"], inplace=True
    )
    print(data_pivot)

    # Create figures.
    figure = plot_bar_stack(
        data=data_pivot,
        label_vertical="Patients",
        label_horizontal="Sex",
        fonts=fonts,
        colors=colors,
        rotation="horizontal",
        legend=True,
    )
    # Specify directories and files.
    file = ("sample_patients.svg")
    path_file = os.path.join(path, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def plot_chart_sample_tissues_samples(
    data=None,
    fonts=None,
    colors=None,
    path=None,
):
    """
    Plots charts from the sample process.

    arguments:
        data (object): Pandas data frame with information about samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        path (str): path to directory

    raises:

    returns:

    """

    data_tissues = data.copy(deep=True)
    data_tissues.reset_index(level=["sample"], inplace=True)
    # Remove relevant redundancy.
    print(
        "Data dimensions before removal of redundancy: " +
        str(data_tissues.shape)
    )
    # Some patients have samples for multiple minor types of the same major
    # type of tissue.
    # Preserve this representation of coverage at the sample level.
    data_tissues.drop_duplicates(
        subset=["patient", "tissue_major", "tissue_minor"],
        keep="first",
        inplace=True,
    )
    print(
        "Data dimensions after removal of redundancy: " +
        str(data_tissues.shape)
    )
    # Remove irrelevant columns.
    data_tissues.drop(
        labels=["sample", "age", "sex"],
        axis="columns",
        inplace=True,
    )
    data_tissues.set_index(
        ["tissue_major", "tissue_minor"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data_tissues)
    data_tissues_counts = data_tissues.groupby(
        level=["tissue_major", "tissue_minor"],
        sort=True,
        as_index=False
    ).size().to_frame(
        name="count"
    )
    data_tissues_counts.reset_index(
        level=["tissue_major", "tissue_minor"], inplace=True
    )
    print(data_tissues_counts)
    # Organize data for chart.
    data_pivot = data_tissues_counts.pivot_table(
        values="count",
        index="tissue_major",
        columns="tissue_minor",
        aggfunc="sum",
        fill_value=0,
    )
    print(data_pivot)
    data_pivot.rename_axis(
        index="groups",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_pivot.reset_index(
        level=["groups"], inplace=True
    )
    print(data_pivot)
    tissues_filter = [
        "breast",
        "cervix",
        "fallopius",
        "ovary",
        "prostate",
        "testis",
        "uterus",
        "vagina",
    ]
    data_filter = data_pivot.loc[~data_pivot["groups"].isin(tissues_filter)]

    # Create figures.
    figure = plot_bar_stack(
        data=data_filter,
        label_vertical="Samples per tissue type",
        label_horizontal="Major tissue type",
        fonts=fonts,
        colors=colors,
        rotation="vertical",
        legend=False,
    )
    # Specify directories and files.
    file = ("sample_tissues_samples.svg")
    path_file = os.path.join(path, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def plot_chart_sample_tissues_patients(
    data=None,
    fonts=None,
    colors=None,
    path=None,
):
    """
    Plots charts from the sample process.

    arguments:
        data (object): Pandas data frame with information about samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        path (str): path to directory

    raises:

    returns:

    """

    data_copy = data.copy(deep=True)
    data_copy.reset_index(level=["sample"], inplace=True)
    # Remove relevant redundancy.
    print(
        "Data dimensions before removal of redundancy: " +
        str(data_copy.shape)
    )
    # Some patients have samples for multiple minor types of the same major
    # type of tissue.
    # Remove this redundancy.
    data_copy.drop_duplicates(
        subset=["patient", "tissue_major"],
        keep="first",
        inplace=True,
    )
    print(
        "Data dimensions after removal of redundancy: " +
        str(data_copy.shape)
    )
    # Remove irrelevant columns.
    data_copy.drop(
        labels=["sample", "age", "sex", "tissue_minor"],
        axis="columns",
        inplace=True,
    )
    data_copy.set_index(
        ["tissue_major"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data_copy)
    data_count = data_copy.groupby(
        level=["tissue_major"],
        sort=True,
        as_index=False
    ).size().to_frame(
        name="count"
    )
    print(data_count)
    # Organize data for chart.
    data_count.rename_axis(
        index="groups",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_count.reset_index(
        level=["groups"], inplace=True
    )
    print(data_count)
    tissues_filter = [
        "breast",
        "cervix",
        "fallopius",
        "ovary",
        "prostate",
        "testis",
        "uterus",
        "vagina",
    ]
    data_filter = data_count.loc[~data_count["groups"].isin(tissues_filter)]

    # Create figures.
    figure = plot_bar_stack(
        data=data_filter,
        label_vertical="Patients per tissue type",
        label_horizontal="Major tissue type",
        fonts=fonts,
        colors=colors,
        rotation="vertical",
        legend=False,
    )
    # Specify directories and files.
    file = ("sample_tissues_patients.svg")
    path_file = os.path.join(path, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def plot_chart_sample_patients_tissues(
    data=None,
    fonts=None,
    colors=None,
    path=None,
):
    """
    Plots charts from the sample process.

    arguments:
        data (object): Pandas data frame with information about samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        path (str): path to directory

    raises:

    returns:

    """

    data_copy = data.copy(deep=True)
    data_copy.reset_index(level=["sample"], inplace=True)
    # Remove relevant redundancy.
    print(
        "Data dimensions before removal of redundancy: " +
        str(data_copy.shape)
    )
    # Some patients have samples for multiple minor types of the same major
    # type of tissue.
    # Remove this redundancy.
    data_copy.drop_duplicates(
        subset=["patient", "tissue_major"],
        keep="first",
        inplace=True,
    )
    print(
        "Data dimensions after removal of redundancy: " +
        str(data_copy.shape)
    )
    # Remove irrelevant columns.
    data_copy.drop(
        labels=["sample", "age", "sex", "tissue_minor"],
        axis="columns",
        inplace=True,
    )
    data_copy.set_index(
        ["patient"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data_copy)
    data_count = data_copy.groupby(
        level=["patient"],
        sort=True,
        as_index=False
    ).size().to_frame(
        name="count"
    )
    print(data_count)
    # Organize data for chart.
    counts = list(data_count.loc[:, "count"])
    # Create figures.
    figure = plot_distribution_histogram(
        series=counts,
        name="",
        bin_method="count",
        bin_count=21,
        label_bins="Bins by tissues sampled for each patient",
        label_counts="Counts of patients in each bin",
        fonts=fonts,
        colors=colors,
        line=True,
        position=7,
        text="",
    )
    # Specify directories and files.
    file = ("sample_patients_tissues.svg")
    path_file = os.path.join(path, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def plot_charts_sample(
    dock=None
):
    """
    Plots charts from the sample process.

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
    path_sample = os.path.join(path_plot, "sample")
    utility.confirm_path_directory(path_sample)
    # Read source information from file.
    source = read_source_sample(dock=dock)
    #source["data_samples"]
    print(source["data_samples"])
    # Patients' sex and age.
    plot_chart_sample_patients(
        data=source["data_samples"],
        fonts=fonts,
        colors=colors,
        path=path_sample
    )
    # Samples per major and minor tissue.
    plot_chart_sample_tissues_samples(
        data=source["data_samples"],
        fonts=fonts,
        colors=colors,
        path=path_sample
    )
    # Patients per major tissue.
    plot_chart_sample_tissues_patients(
        data=source["data_samples"],
        fonts=fonts,
        colors=colors,
        path=path_sample
    )
    # Tissues per patient.
    plot_chart_sample_patients_tissues(
        data=source["data_samples"],
        fonts=fonts,
        colors=colors,
        path=path_sample
    )

    pass


# Tissue


def read_source_tissue(dock=None):
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
    path_tissue = os.path.join(dock, "tissue")
    path_sample_component = os.path.join(
        path_tissue, "data_sample_component.pickle"
    )
    # Read information from file.
    data_sample_component = pandas.read_pickle(path_sample_component)
    # Compile and return information.
    return {
        "data_sample_component": data_sample_component,
    }


def plot_chart_tissue_component(
    data=None,
    fonts=None,
    colors=None,
    path=None,
):
    """
    Plots charts from the tissue process.

    arguments:
        data (object): Pandas data frame with information about samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        path (str): path to directory

    raises:

    returns:

    """

    # Organize data.
    # Drop the sample designations.
    data.drop(
        labels=["tissue_minor"],
        axis="columns",
        inplace=True
    )
    # Create charts for combinations of principal components.
    ##########
    # Component 1 and Component 2.
    # Create figures.
    figure = plot_scatter_cluster(
        data=data,
        abscissa="component_1",
        ordinate="component_2",
        label_horizontal="Principal Component 1",
        label_vertical="Principal Component 2",
        factor="tissue_major",
        fonts=fonts,
        colors=colors,
        legend=True,
    )
    # Specify directories and files.
    file = ("tissue_components_1_2.svg")
    path_file = os.path.join(path, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )
    ##########
    # Component 1 and Component 3.
    # Create figures.
    figure = plot_scatter_cluster(
        data=data,
        abscissa="component_1",
        ordinate="component_3",
        label_horizontal="Principal Component 1",
        label_vertical="Principal Component 3",
        factor="tissue_major",
        fonts=fonts,
        colors=colors,
        legend=True,
    )
    # Specify directories and files.
    file = ("tissue_components_1_3.svg")
    path_file = os.path.join(path, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )
    ##########
    # Component 2 and Component 3.
    # Create figures.
    figure = plot_scatter_cluster(
        data=data,
        abscissa="component_2",
        ordinate="component_3",
        label_horizontal="Principal Component 2",
        label_vertical="Principal Component 3",
        factor="tissue_major",
        fonts=fonts,
        colors=colors,
        legend=True,
    )
    # Specify directories and files.
    file = ("tissue_components_2_3.svg")
    path_file = os.path.join(path, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def plot_charts_tissue(
    dock=None
):
    """
    Plots charts from the tissue process.

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
    path_tissue = os.path.join(path_plot, "tissue")
    utility.confirm_path_directory(path_tissue)
    # Read source information from file.
    source = read_source_tissue(dock=dock)
    # Principal components.
    data = source["data_sample_component"].reset_index(
        level=["sample", "person", "tissue_major", "tissue_minor"],
        inplace=False
    )
    # Create charts for combinations of principal components.
    plot_chart_tissue_component(
        data=data,
        fonts=fonts,
        colors=colors,
        path=path_tissue
    )

    pass





def read_source_restriction(dock=None):
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
    path_collection = os.path.join(dock, "collection")
    path_genes_patients = os.path.join(
        path_collection, "genes_patients.pickle"
    )
    # Read information from file.
    with open(path_genes_patients, "rb") as file_source:
        genes_patients = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes_patients": genes_patients,
    }


def plot_charts_restriction(
    dock=None
):
    """
    Plots charts from the tissue process.

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
    path_restriction = os.path.join(path_plot, "restriction")
    utility.confirm_path_directory(path_restriction)
    # Read source information from file.
    source = read_source_restriction(dock=dock)
    # Principal components.
    print(source["genes_patients"])
    print(str(len(source["genes_patients"])))
    # Create charts for combinations of principal components.
    plot_chart_restriction_genes_patients(
        data=data,
        fonts=fonts,
        colors=colors,
        path=path_tissue
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


    #plot_charts_analysis(dock=dock)
    #plot_charts_sample(dock=dock)
    #plot_charts_tissue(dock=dock)
    plot_charts_restriction(dock=dock)

    pass


if (__name__ == "__main__"):
    execute_procedure()
