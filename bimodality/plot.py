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
import matplotlib_venn
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


def plot_heatmap(
    data=None,
    label_rows=None,
    label_columns=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type heatmap.

    arguments:
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap
        label_rows (bool): whether to include explicit labels on heatmap's rows
        label_columns (bool): whether to include explicit labels on heatmap's
            columns
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    # Map data's columns to heatmap's rows.
    # Map data's rows to heatmap's columns.
    labels_rows = data.columns.to_list()
    labels_columns = data.index.to_list()
    matrix = numpy.transpose(data.values)

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = matplotlib.pyplot.axes()

    # Represent data as a color grid.
    image = axes.imshow(
        matrix,
        #cmap= ,
        aspect="auto",
    )

    # Create legend for color map.
    label_bar = "legend"
    bar = axes.figure.colorbar(image, ax=axes)
    bar.ax.set_ylabel(label_bar, rotation=-90, va="bottom")

    # Create ticks and labels for each grid.
    # Let the horizontal axes labeling appear on top.
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["two"]["size"],
        labelcolor=colors["black"],
        top=True,
        bottom=False,
        labeltop=True,
        labelbottom=False
    )
    # Rotate the tick labels and set their alignment.
    matplotlib.pyplot.setp(
        axes.get_xticklabels(),
        rotation=-30,
        ha="right",
        rotation_mode="anchor"
    )

    if label_columns:
        axes.set_xticks(numpy.arange(matrix.shape[1]))
        axes.set_xticklabels(labels_columns)
    if label_rows:
        axes.set_yticks(numpy.arange(matrix.shape[0]))
        axes.set_yticklabels(
            labels_rows,
            minor=False,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )

    return figure

# TODO: This function is in progress...
# TODO: I'm trying to include a bar or bars to designate groups of persons
def plot_heatmap_groups(
    data=None,
    label_rows=None,
    label_columns=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type heatmap.

    arguments:
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap
        label_rows (bool): whether to include explicit labels on heatmap's rows
        label_columns (bool): whether to include explicit labels on heatmap's
            columns
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    # Map data's columns to heatmap's rows.
    # Map data's rows to heatmap's columns.
    labels_rows = data.columns.to_list()
    labels_columns = data.index.to_list()
    matrix = numpy.transpose(data.values)

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = figure.subplots(
        nrows=2,
        ncols=1,
        sharex=True,
        squeeze=True,
        gridspec_kw=dict(height_ratios=[1, 5])
    )

    axes[0].axis("off")

    # Represent data as a color grid.
    image = axes[1].imshow(
        matrix,
        #cmap= ,
        aspect="auto",
    )

    # Create legend for color map.
    label_bar = "legend"
    bar = axes[1].figure.colorbar(image, ax=axes[1])
    bar.ax.set_ylabel(label_bar, rotation=-90, va="bottom")

    # Create ticks and labels for each grid.
    # Let the horizontal axes labeling appear on top.
    axes[1].tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["two"]["size"],
        labelcolor=colors["black"],
        top=False,
        bottom=True,
        labeltop=False,
        labelbottom=True
    )
    # Rotate the tick labels and set their alignment.
    matplotlib.pyplot.setp(
        axes[1].get_xticklabels(),
        rotation=-30,
        ha="right",
        rotation_mode="anchor"
    )

    if label_columns:
        axes[1].set_xticks(numpy.arange(matrix.shape[1]))
        axes[1].set_xticklabels(labels_columns)
    if label_rows:
        axes[1].set_yticks(numpy.arange(matrix.shape[0]))
        axes[1].set_yticklabels(
            labels_rows,
            minor=False,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )

    return figure


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
    axes.set_title(
        name,
        #fontdict=fonts["properties"]["one"],
        loc="center",
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
    color_count=None,
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
        color_count (int): count of colors for subcategories (stacks)
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

    if color_count > 1:
        colors_series = list(
            seaborn.color_palette("hls", n_colors=color_count)
        )
    else:
        colors_series = [colors["blue"]]

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
    index = 0
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
            color=colors_series[index],
        )
        bars.append(series_bars[0])
        # Restore bases.
        bases = list(map(sum, zip(bases, values)))
        print(bases)
        index += 1
        pass

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
        mean_x = statistics.mean(values_x)
        values_y = group[ordinate].to_list()
        mean_y = statistics.mean(values_y)
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


def plot_overlap_sets(
    sets=None,
    count=None,
    fonts=None,
    colors=None,
):
    """
    Creates a Venn Diagram to represent overlap between multiple sets.

    arguments:
        sets (dict<list<str>>): values in sets
        count (int): count of sets, either 2 or 3
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize information.
    names = list(sets.keys())

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    if count == 2:
        venn = matplotlib_venn.venn2(
            subsets=[
                set(sets[names[0]]),
                set(sets[names[1]]),
            ],
            set_labels=names,
        )
    elif count == 3:
        venn = matplotlib_venn.venn3(
            subsets=[
                set(sets[names[0]]),
                set(sets[names[1]]),
                set(sets[names[2]]),
            ],
            set_labels=names,
        )

    return figure


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


##########
# Sex and age of persons in GTEx
# Status: working


def read_source_person_sex_age(dock=None):
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
    path_selection = os.path.join(dock, "selection")
    path_data_persons_sex_age_counts = os.path.join(
        path_selection, "data_persons_sex_age_counts.pickle"
    )
    # Read information from file.
    with open(path_data_persons_sex_age_counts, "rb") as file_source:
        data_persons_sex_age_counts = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_persons_sex_age_counts": data_persons_sex_age_counts,
    }


def prepare_chart_person_sex_age(
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

    # Read source information from file.
    source = read_source_person_sex_age(dock=dock)
    print(source["data_persons_sex_age_counts"])

    # Copy data.
    data_persons_sex_age_counts = (
        source["data_persons_sex_age_counts"].copy(deep=True)
    )
    # Organize data.
    data_pivot = data_persons_sex_age_counts.pivot_table(
        values="count",
        index="sex",
        columns="age_decade",
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
    print("these are the data passed to the plot function")
    print(data_pivot)

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    path_persons = os.path.join(path_plot, "persons")
    utility.create_directory(path_persons)

    # Create figures.
    figure = plot_bar_stack(
        data=data_pivot,
        label_vertical="Persons",
        label_horizontal="Sex",
        fonts=fonts,
        colors=colors,
        color_count=6,
        rotation="horizontal",
        legend=True,
    )
    # Specify directories and files.
    file = ("persons_sex_age.svg")
    path_file = os.path.join(path_persons, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


##########
# Counts of persons per major tissue type
# Status: working

def read_source_persons_per_tissue(dock=None):
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
    path_selection = os.path.join(dock, "selection")
    path_data_persons_per_tissue = os.path.join(
        path_selection, "data_persons_per_tissue.pickle"
    )
    # Read information from file.
    with open(path_data_persons_per_tissue, "rb") as file_source:
        data_persons_per_tissue = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_persons_per_tissue": data_persons_per_tissue,
    }


def prepare_chart_persons_per_tissue(dock=None):
    """
    Plots charts from the sample process.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_persons_per_tissue(dock=dock)
    print(source["data_persons_per_tissue"])

    # Copy data.
    data_persons_per_tissue = (
        source["data_persons_per_tissue"].copy(deep=True)
    )
    # Organize data.
    data_persons_per_tissue.set_index(
        "tissue",
        drop=True,
        inplace=True,
    )
    data_persons_per_tissue.rename_axis(
        index="groups",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_persons_per_tissue.reset_index(
        level=["groups"], inplace=True
    )
    print(data_persons_per_tissue)

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    path_persons = os.path.join(path_plot, "persons")
    utility.create_directory(path_persons)

    # Create figures.
    figure = plot_bar_stack(
        data=data_persons_per_tissue,
        label_vertical="Persons per tissue type",
        label_horizontal="Major tissue type",
        fonts=fonts,
        colors=colors,
        color_count=1,
        rotation="vertical",
        legend=False,
    )
    # Specify directories and files.
    file = ("persons_per_tissue.svg")
    path_file = os.path.join(path_persons, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


##########
# Counts of tissues per person, histogram with bins by counts
# Status: working


def read_source_tissues_per_person(dock=None):
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
    path_selection = os.path.join(dock, "selection")
    path_data_tissues_per_person = os.path.join(
        path_selection, "data_tissues_per_person.pickle"
    )
    path_tissues_per_person = os.path.join(
        path_selection, "tissues_per_person.pickle"
    )
    # Read information from file.
    with open(path_data_tissues_per_person, "rb") as file_source:
        data_tissues_per_person = pickle.load(file_source)
    with open(path_tissues_per_person, "rb") as file_source:
        tissues_per_person = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_tissues_per_person": data_tissues_per_person,
        "tissues_per_person": tissues_per_person,
    }


def prepare_chart_tissues_per_person(dock=None):
    """
    Plots charts from the sample process.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_tissues_per_person(dock=dock)
    print(source["data_tissues_per_person"])
    print(source["tissues_per_person"])

    # Organize data.

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    path_persons = os.path.join(path_plot, "persons")
    utility.create_directory(path_persons)

    # Create figures.
    figure = plot_distribution_histogram(
        series=source["tissues_per_person"],
        name="",
        bin_method="count",
        bin_count=17,
        label_bins="Bins by count of tissues per person",
        label_counts="Counts of persons in each bin",
        fonts=fonts,
        colors=colors,
        line=True,
        position=11,
        text="",
    )
    # Specify directories and files.
    file = ("tissues_per_person.svg")
    path_file = os.path.join(path_persons, file)
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


##########
# Distributions of each bimodality measure's scores across genes
# Status: working

def read_source_modality_gene_distribution(
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_distribution = os.path.join(dock, "distribution")
    path_genes_scores = os.path.join(
        path_distribution, "collection", "genes_scores.pickle"
    )
    path_scores = os.path.join(
        path_distribution, "collection", "scores.pickle"
    )

    path_candidacy = os.path.join(dock, "candidacy")
    path_measures_thresholds = os.path.join(
        path_candidacy, "measures_thresholds.pickle"
    )

    # Read information from file.
    with open(path_genes_scores, "rb") as file_source:
        genes_scores = pickle.load(file_source)
    with open(path_scores, "rb") as file_source:
        scores = pickle.load(file_source)
    with open(path_measures_thresholds, "rb") as file_source:
        measures_thresholds = pickle.load(file_source)

    # Compile and return information.
    return {
        "genes_scores": genes_scores,
        "scores": scores,
        "measures_thresholds": measures_thresholds,
    }


def plot_chart_modality_gene_distribution(
    values=None,
    threshold=None,
    path=None
):
    """
    Plots charts from the analysis process.

    arguments:
        values (list<float>): values
        threshold (float): value of threshold for which to draw line
        path (str): path to directory and file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_distribution_histogram(
        series=values,
        name="",
        bin_method="count",
        bin_count=50,
        label_bins="Bins",
        label_counts="Counts",
        fonts=fonts,
        colors=colors,
        line=True,
        position=threshold,
        text="",
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_charts_modality_gene_distribution(
    dock=None
):
    """
    Plots charts.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    print("going to plot distribution of bimodality metrics across genes")

    # Read source information from file.
    source = read_source_modality_gene_distribution(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_candidacy = os.path.join(path_plot, "candidacy")
    path_measure = os.path.join(path_candidacy, "measure")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_measure)
    utility.create_directories(path=path_measure)

    # Specify bimodality measures.
    measures = ["dip", "mixture", "coefficient"]
    # Iterate on bimodality measures.
    for measure in measures:

        utility.print_terminal_partition(level=3)
        print(measure)

        # Calculate selection threshold.
        threshold = source["measures_thresholds"][measure]

        # Prepare charts.
        # Specify directories and files.
        path_measure_figure = os.path.join(
            path_measure, str(measure + ".svg")
        )

        plot_chart_modality_gene_distribution(
            values=source["scores"][measure]["values"],
            threshold=threshold,
            path=path_measure_figure
        )

        pass

    pass


##########
# Overlap between sets in selection of genes by bimodality
# Status: working


def read_source_gene_sets_candidacy(
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_candidacy = os.path.join(dock, "candidacy")
    path_sets_genes_measures = os.path.join(
        path_candidacy, "sets_genes_measures.pickle"
    )
    # Read information from file.
    with open(path_sets_genes_measures, "rb") as file_source:
        sets_genes_measures = pickle.load(file_source)
    # Compile and return information.
    return {
        "sets_genes_measures": sets_genes_measures,
    }


def plot_chart_gene_sets_candidacy(
    sets=None,
    path=None
):
    """
    Plots charts from the analysis process.

    arguments:
        sets (dict<list<str>>): values in sets
        path (str): path to directory and file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_overlap_sets(
        sets=sets,
        count=3,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_charts_gene_sets_candidacy(
    dock=None
):
    """
    Plots charts.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    print("going to plot overlap between sets of genes by bimodality")

    # Read source information from file.
    source = read_source_gene_sets_candidacy(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_candidacy = os.path.join(path_plot, "candidacy")
    path_set = os.path.join(path_candidacy, "set")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_set)
    utility.create_directories(path=path_set)
    path_set_figure = os.path.join(path_set, "sets.svg")

    plot_chart_gene_sets_candidacy(
        sets=source["sets_genes_measures"],
        path=path_set_figure,
    )

    pass


##########
# Overlap between sets in selection of genes by permutation probability of
# bimodality
# Status: working


def read_source_gene_sets_probability(
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_probability = os.path.join(dock, "probability")
    path_sets_genes_measures = os.path.join(
        path_probability, "sets_genes_measures.pickle"
    )
    # Read information from file.
    with open(path_sets_genes_measures, "rb") as file_source:
        sets_genes_measures = pickle.load(file_source)
    # Compile and return information.
    return {
        "sets_genes_measures": sets_genes_measures,
    }


def plot_chart_gene_sets_probability(
    sets=None,
    path=None
):
    """
    Plots charts from the analysis process.

    arguments:
        sets (dict<list<str>>): values in sets
        path (str): path to directory and file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_overlap_sets(
        sets=sets,
        count=3,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_charts_gene_sets_probability(
    dock=None
):
    """
    Plots charts.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    print("going to plot overlap between sets of genes by probability")

    # Read source information from file.
    source = read_source_gene_sets_probability(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_probability = os.path.join(path_plot, "probability")
    path_set = os.path.join(path_probability, "set")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_set)
    utility.create_directories(path=path_set)
    path_set_figure = os.path.join(path_set, "sets.svg")

    plot_chart_gene_sets_probability(
        sets=source["sets_genes_measures"],
        path=path_set_figure,
    )

    pass


##########
# Overlap between sets in selection of genes by heritability
# Status: working


def read_source_gene_heritability(
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_split = os.path.join(dock, "split")
    path_genes_split = os.path.join(path_split, "genes.txt")

    path_heritability = os.path.join(dock, "heritability")
    path_collection = os.path.join(path_heritability, "collection")
    path_sets_genes_models = os.path.join(
        path_collection, "sets_genes_models.pickle"
    )
    # Read information from file.
    genes_split = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_genes_split,
    )
    with open(path_sets_genes_models, "rb") as file_source:
        sets_genes_models = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes_split": genes_split,
        "sets_genes_models": sets_genes_models,
    }


def plot_chart_gene_heritability(
    sets=None,
    path=None
):
    """
    Plots charts from the analysis process.

    arguments:
        sets (dict<list<str>>): values in sets
        path (str): path to directory and file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_overlap_sets(
        sets=sets,
        count=3,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_charts_gene_heritability(
    dock=None
):
    """
    Plots charts.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    print("going to plot overlap between sets of genes by heritability")

    # Read source information from file.
    source = read_source_gene_heritability(dock=dock)

    # Organize information.
    sets = dict()
    sets["total"] = source["genes_split"]
    sets["simple"] = source["sets_genes_models"]["simple"]
    sets["complex"] = source["sets_genes_models"]["complex"]

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_heritability = os.path.join(path_plot, "heritability")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_heritability)
    utility.create_directory(path=path_heritability)
    path_sets = os.path.join(path_heritability, "sets.svg")

    plot_chart_gene_heritability(
        sets=sets,
        path=path_sets,
    )

    pass


##########
# Overlap between sets in selection of genes by integration
# Status: working


def read_source_gene_sets_integration(
    dock=None
):
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
    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_candidacy = os.path.join(
        path_candidacy, "genes_candidacy.pickle"
    )

    path_probability = os.path.join(dock, "probability")
    path_genes_probability = os.path.join(
        path_probability, "genes_probability.pickle"
    )

    path_heritability = os.path.join(dock, "heritability")
    path_collection = os.path.join(path_heritability, "collection")
    path_genes_heritability = os.path.join(
        path_collection, "genes_heritability.pickle"
    )
    # Read information from file.
    with open(path_genes_candidacy, "rb") as file_source:
        genes_candidacy = pickle.load(file_source)
    with open(path_genes_probability, "rb") as file_source:
        genes_probability = pickle.load(file_source)
    with open(path_genes_heritability, "rb") as file_source:
        genes_heritability = pickle.load(file_source)

    # Compile and return information.
    return {
        "genes_candidacy": genes_candidacy,
        "genes_probability": genes_probability,
        "genes_heritability": genes_heritability,
    }


def plot_chart_gene_sets_integration(
    sets=None,
    path=None
):
    """
    Plots charts from the analysis process.

    arguments:
        sets (dict<list<str>>): values in sets
        path (str): path to directory and file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_overlap_sets(
        sets=sets,
        count=3,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_charts_gene_sets_integration(
    dock=None
):
    """
    Plots charts.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    print("going to plot overlap between sets of genes by probability")

    # Read source information from file.
    source = read_source_gene_sets_integration(dock=dock)

    # Organize information.
    sets = dict()
    sets["candidacy"] = source["genes_candidacy"]
    sets["probability"] = source["genes_probability"]
    sets["heritability"] = source["genes_heritability"]

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_integration = os.path.join(path_plot, "integration")
    path_set = os.path.join(path_integration, "set")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_set)
    utility.create_directories(path=path_set)
    path_set_figure = os.path.join(path_set, "sets.svg")

    plot_chart_gene_sets_probability(
        sets=sets,
        path=path_set_figure,
    )

    pass


##########
# Distributions of genes' pan-tissue aggregate signals across persons
# Status: working


def read_source_genes_persons_signals_initial(
    dock=None
):
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
    path_integration = os.path.join(dock, "integration")
    path_genes_integration = os.path.join(
        path_integration, "genes_integration.pickle"
    )

    # Read information from file.
    with open(path_genes_integration, "rb") as file_source:
        genes_integration = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes_integration": genes_integration,
    }


def read_source_genes_persons_signals(
    gene=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_distribution = os.path.join(dock, "distribution")
    path_gene = os.path.join(path_distribution, "genes", gene)
    path_data_gene_persons_signals = os.path.join(
        path_gene, "data_gene_persons_signals.pickle"
    )

    # Read information from file.
    data_gene_persons_signals = (
        pandas.read_pickle(path_data_gene_persons_signals)
    )
    # Compile and return information.
    return {
        "data_gene_persons_signals": data_gene_persons_signals,
    }


def plot_chart_genes_persons_signals(
    gene=None,
    values=None,
    path=None
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a single gene
        values (list<float>): values of a gene's pan-tissue aggregate signals
            across persons
        path (str): path for file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_distribution_histogram(
        series=values,
        name=gene,
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
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_charts_genes_persons_signals(
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

    # Read source information from file.
    source_initial = read_source_genes_persons_signals_initial(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_distribution = os.path.join(path_plot, "distribution")
    path_genes_persons = os.path.join(path_distribution, "genes_persons")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_genes_persons)
    utility.create_directories(path=path_genes_persons)

    # Iterate on genes.
    for gene in source_initial["genes_integration"]:

        # Read source information from file.
        source_gene = read_source_genes_persons_signals(
            gene=gene,
            dock=dock
        )
        # Access information.
        data_gene_persons_signals = source_gene["data_gene_persons_signals"]
        # Distribution of gene's signals across persons.
        values = data_gene_persons_signals["value"].to_list()

        # Create charts for the gene.
        path_gene_figure = os.path.join(path_genes_persons, str(gene + ".svg"))
        plot_chart_genes_persons_signals(
            gene=gene,
            values=values,
            path=path_gene_figure
        )

        pass

    pass


##########
# Genes' signals across tissues and persons
# Status: working


def read_source_genes_signals_tissues_persons_initial(
    dock=None
):
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
    path_integration = os.path.join(dock, "integration")
    path_genes_integration = os.path.join(
        path_integration, "genes_integration.pickle"
    )

    # Read information from file.
    with open(path_genes_integration, "rb") as file_source:
        genes_integration = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes_integration": genes_integration,
    }


def read_source_genes_signals_tissues_persons(
    gene=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_distribution = os.path.join(dock, "distribution")
    path_gene = os.path.join(path_distribution, "genes", gene)
    path_data_gene_persons_signals = os.path.join(
        path_gene, "data_gene_persons_signals.pickle"
    )
    path_data_gene_signals_tissues_persons = os.path.join(
        path_gene, "data_gene_signals_tissues_persons_normal_standard.pickle"
    )

    # Read information from file.
    data_gene_persons_signals = (
        pandas.read_pickle(path_data_gene_persons_signals)
    )
    data_gene_signals_tissues_persons = (
        pandas.read_pickle(path_data_gene_signals_tissues_persons)
    )
    # Compile and return information.
    return {
        "data_gene_persons_signals": data_gene_persons_signals,
        "data_gene_signals_tissues_persons": data_gene_signals_tissues_persons,
    }


def organize_genes_signals_tissues_persons(
    data_gene_persons_signals=None,
    data_gene_signals_tissues_persons=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons
        data_gene_signals_tissues_persons (object): Pandas data frame of a
            gene's signals across tissues and persons

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across tissues and
            persons

    """

    # Rename the aggregate, pantissue signals.
    data_gene_persons_signals.rename(
        columns={
            "value": "signal",
        },
        inplace=True,
    )
    # Introduce aggregate, pantissue signals to tissue-person matrix.
    # These aggregate signals will be useful to sort the data.
    # Join
    data_hybrid = data_gene_signals_tissues_persons.join(
        data_gene_persons_signals,
        how="left",
        on="person"
    )
    # Sort data by the aggregate, pantissue signals.
    data_hybrid.sort_values(
        by=["signal"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Remove the column for aggregate, pantissue signals.
    data_hybrid.drop(
        labels="signal",
        axis="columns",
        inplace=True
    )
    # Return information
    return data_hybrid


def plot_chart_genes_signals_tissues_persons(
    gene=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a single gene
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(gene + ".svg")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_heatmap(
        data=data,
        label_columns=False,
        label_rows=True,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_genes_signals_tissues_persons(
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

    # Read source information from file.
    source_initial = (
        read_source_genes_signals_tissues_persons_initial(dock=dock)
    )

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_distribution = os.path.join(path_plot, "distribution")
    path_genes_tissues_persons = os.path.join(
        path_distribution, "genes_tissues_persons"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_genes_tissues_persons)
    utility.create_directories(path=path_genes_tissues_persons)

    # Iterate on genes.
    for gene in source_initial["genes_integration"]:

        # Read source information from file.
        source_gene = read_source_genes_signals_tissues_persons(
            gene=gene,
            dock=dock
        )
        # Organize data.
        # Sort persons by their pantissue aggregate signals for the gene.
        # The same order is important to compare the heatmap to the histogram.
        data_gene_signals_tissues_persons = (
            organize_genes_signals_tissues_persons(
                data_gene_persons_signals=source_gene["data_gene_persons_signals"],
                data_gene_signals_tissues_persons=(
                    source_gene["data_gene_signals_tissues_persons"]
                ),
            )
        )
        # Create charts for the gene.
        plot_chart_genes_signals_tissues_persons(
            gene=gene,
            data=data_gene_signals_tissues_persons,
            path_directory=path_genes_tissues_persons
        )
        pass
    pass


##########
# Genes' signals across groups of persons
# Status: working


def read_source_signals_genes_persons_groups(
    dock=None
):
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
    path_selection = os.path.join(dock, "selection")
    path_persons_properties = os.path.join(
        path_selection, "data_persons_properties.pickle"
    )
    path_distribution = os.path.join(dock, "distribution")
    path_collection = os.path.join(path_distribution, "collection")
    path_data_signals_genes_persons = os.path.join(
        path_collection, "data_signals_genes_persons.pickle"
    )
    path_integration = os.path.join(dock, "integration")
    path_genes_integration = os.path.join(
        path_integration, "genes_integration.pickle"
    )

    # Read information from file.
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    with open(path_genes_integration, "rb") as file_source:
        genes_integration = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_persons_properties": data_persons_properties,
        "data_signals_genes_persons": data_signals_genes_persons,
        "genes_integration": genes_integration,
    }


def organize_signals_genes_persons(
    genes=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        genes (list<str>): identifiers of genes
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        data_signals_genes_persons (object): Pandas data frame of genes'
            pantissue signals across persons

    raises:

    returns:
        (object): Pandas data frame of genes' pantissue signals across persons

    """

    # Copy data.
    data_persons = data_persons_properties.copy(deep=True)
    data_signals = data_signals_genes_persons.copy(deep=True)
    # Organize data.
    data_persons = data_persons.loc[
        :, data_persons.columns.isin([
            "age", "body", "delay",
            "hardiness", "season", "sex",
        ])
    ]
    data_signals = data_signals.loc[
        :, data_signals.columns.isin(genes)
    ]

    # Introduce aggregate, pantissue signals for each gene to person data.
    # Join
    data_hybrid = data_signals.join(
        data_persons,
        how="left",
        on="person"
    )

    # Sort data by the aggregate, pantissue signals.
    data_hybrid.sort_values(
        by=["sex", "age"], # ["sex", "age"]
        axis="index",
        ascending=True,
        inplace=True,
    )

    if False:
        utility.print_terminal_partition(level=2)
        print("check sort order...")
        print(data_hybrid)

    # Remove unnecessary columns.
    data_hybrid.drop(
        labels=[
            "age", "body", "delay",
            "hardiness", "season", "sex",
        ],
        axis="columns",
        inplace=True
    )

    # Remove persons without signal for any genes.
    data_hybrid.dropna(
        axis="index",
        how="all",
        thresh=1,
        inplace=True,
    )
    # Return information
    return data_hybrid


# TODO: I eventually want to use function "plot_heatmap_groups"
def plot_chart_signals_genes_persons_groups(
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str("genes_signals_persons.svg")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    # TODO: I eventually want to use function "plot_heatmap_groups"
    figure = plot_heatmap(
        data=data,
        label_columns=False,
        label_rows=True,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_signals_genes_persons_groups(
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

    # Read source information from file.
    source = read_source_signals_genes_persons_groups(dock=dock)

    # Organize data.
    # Combine all genes' signals in single data frame.
    # Sort persons by sex and age.

    ######
    # TODO: Might be good to include a designation of sort groups in the arguments to this function...
    ######
    data = organize_signals_genes_persons(
        genes=source["genes_integration"],
        data_persons_properties=source["data_persons_properties"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
    )

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_distribution = os.path.join(path_plot, "distribution")
    path_genes_persons = os.path.join(
        path_distribution, "genes_persons_groups"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_genes_persons)
    utility.create_directories(path=path_genes_persons)
    # Create chart.
    plot_chart_signals_genes_persons_groups(
        data=data,
        path_directory=path_genes_persons
    )

    pass


##########
# Correlations in signals between pairs of genes
# Status: in progres...

# TODO: I need to move this functionality from integration procedure to here


################# Need to Update ###############################

# Sample

# TODO: this obsolete function plots the minor tissue categories within each major tissue category...
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

    path_variance_all = os.path.join(
        path_tissue, "data_component_variance_all.pickle"
    )
    path_component_all = os.path.join(
        path_tissue, "data_factor_component_all.pickle"
    )
    path_variance_few = os.path.join(
        path_tissue, "data_component_variance_few.pickle"
    )
    path_component_few = os.path.join(
        path_tissue, "data_factor_component_few.pickle"
    )

    # Read information from file.
    data_variance_all = pandas.read_pickle(path_variance_all)
    data_component_all = pandas.read_pickle(path_component_all)
    data_variance_few = pandas.read_pickle(path_variance_few)
    data_component_few = pandas.read_pickle(path_component_few)
    # Compile and return information.
    return {
        "data_variance_all": data_variance_all,
        "data_component_all": data_component_all,
        "data_variance_few": data_variance_few,
        "data_component_few": data_component_few,
    }


def plot_chart_tissue_component(
    data_component=None,
    data_variance=None,
    file_suffix=None,
    fonts=None,
    colors=None,
    path=None,
):
    """
    Plots charts from the tissue process.

    arguments:
        data_component (object): Pandas data frame with information about
            components
        data_variance (object): Pandas data frame with information about
            variance of each component
        file_suffix (str): suffix to append to the end of file names
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        path (str): path to directory

    raises:

    returns:

    """

    # Organize data.
    # Extract the variance of each component.
    variance_one = data_variance.loc[
        data_variance["component"] == 1,
        "variance"
    ]
    variance_one = round(float(variance_one * 100), 1)

    variance_two = data_variance.loc[
        data_variance["component"] == 2,
        "variance"
    ]
    variance_two = round(float(variance_two * 100), 1)

    variance_three = data_variance.loc[
        data_variance["component"] == 3,
        "variance"
    ]
    variance_three = round(float(variance_three * 100), 1)

    # Principal components for all tissues.
    data_component = data_component.reset_index(
        level=["sample", "person", "tissue_major", "tissue_minor"],
        inplace=False
    )
    # Drop the sample designations.
    data_component.drop(
        labels=["tissue_minor"],
        axis="columns",
        inplace=True
    )
    # Create charts for combinations of principal components.
    ##########
    # Component 1 and Component 2.
    # Create figures.
    figure = plot_scatter_cluster(
        data=data_component,
        abscissa="component_1",
        ordinate="component_2",
        label_horizontal=("Component 1 (" + str(variance_one) + "%)"),
        label_vertical=("Component 2 (" + str(variance_two) + "%)"),
        factor="tissue_major",
        fonts=fonts,
        colors=colors,
        legend=True,
    )
    # Specify directories and files.
    file = ("tissue_components_1_2" + file_suffix + ".svg")
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
        data=data_component,
        abscissa="component_1",
        ordinate="component_3",
        label_horizontal=("Component 1 (" + str(variance_one) + "%)"),
        label_vertical=("Component 3 (" + str(variance_three) + "%)"),
        factor="tissue_major",
        fonts=fonts,
        colors=colors,
        legend=True,
    )
    # Specify directories and files.
    file = ("tissue_components_1_3" + file_suffix + ".svg")
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
        data=data_component,
        abscissa="component_2",
        ordinate="component_3",
        label_horizontal=("Component 2 (" + str(variance_two) + "%)"),
        label_vertical=("Component 3 (" + str(variance_three) + "%)"),
        factor="tissue_major",
        fonts=fonts,
        colors=colors,
        legend=True,
    )
    # Specify directories and files.
    file = ("tissue_components_2_3" + file_suffix + ".svg")
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
    utility.remove_directory(path=path_tissue)
    utility.create_directory(path_tissue)
    # Read source information from file.
    source = read_source_tissue(dock=dock)

    # All tissues.
    # Create charts for combinations of principal components.
    path_study = os.path.join(path_tissue, "all")
    utility.create_directory(path_study)
    plot_chart_tissue_component(
        data_component=source["data_component_all"],
        data_variance=source["data_variance_all"],
        file_suffix="",
        fonts=fonts,
        colors=colors,
        path=path_study
    )

    # Few tissues.
    # Create charts for combinations of principal components.
    path_study = os.path.join(path_tissue, "few")
    utility.create_directory(path_study)
    plot_chart_tissue_component(
        data_component=source["data_component_few"],
        data_variance=source["data_variance_few"],
        file_suffix="",
        fonts=fonts,
        colors=colors,
        path=path_study
    )

    pass


def read_source_tissues_persons(dock=None):
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
    path_manual = os.path.join(dock, "manual")
    path_tissues_persons = os.path.join(
        path_manual, "tissues_persons.tsv"
    )

    # Read information from file.
    data_tissues_persons = pandas.read_csv(
        path_tissues_persons,
        sep="\t",
        header=0,
    )
    # Compile and return information.
    return {
        "data_tissues_persons": data_tissues_persons,
    }


def plot_chart_tissues_persons(
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
    utility.create_directory(path_plot)
    path_tissues_persons = os.path.join(path_plot, "tissues_persons")
    # Remove previous files since they change from run to run.
    utility.remove_directory(path=path_tissues_persons)
    utility.create_directory(path_tissues_persons)
    # Read source information from file.
    source = read_source_tissues_persons(dock=dock)

    # TODO: create scatter plot...


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
    utility.create_directory(path_restriction)
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

    # Remove previous files to avoid version or batch confusion.
    path_plot = os.path.join(dock, "plot")
    utility.remove_directory(path=path_plot)
    utility.create_directory(path_plot)

    if False:
        # Plot chart for sex and age of persons in GTEx consortium.
        prepare_chart_person_sex_age(dock=dock)

        # Plot chart for counts of persons per major tissue type.
        prepare_chart_persons_per_tissue(dock=dock)

        # Plot chart for counts of major tissue type per person, a histogram.
        prepare_chart_tissues_per_person(dock=dock)

        # Plot charts of distributions of each bimodality measure's scores across
        # genes.
        #- distribution of a single bimodality measure's scores across all genes
        #-- histogram
        #-- bins by bimodality measure's score (z-score first?)
        #-- count of genes in each bin
        prepare_charts_modality_gene_distribution(dock=dock)

    # Plot charts of distributions of genes' pan-tissue aggregate signals
    # across persons.
    prepare_charts_genes_persons_signals(dock=dock)

    # Plot charts, heatmaps, for each gene's signals across persons (columns)
    # and tissues (rows).
    prepare_charts_genes_signals_tissues_persons(dock=dock)

    # TODO: plot aggregate, pantissue signals across persons for each hit gene
    # TODO: plot these as a heatmap with genes across rows and persons across columns
    # TODO: figure out how to sort people columns by 1. sex, 2. age
    # Plot charts, heatmaps, for each gene's aggregate, pantissue signals
    # across persons.
    # Genes will be across rows, and persons will be across columns.
    # This chart will illustrate whether the same persons are in low and high
    # groups for multiple genes.
    prepare_charts_signals_genes_persons_groups(dock=dock)

    if False:
        # Plot charts of overlap between sets in selection of genes by bimodality.
        prepare_charts_gene_sets_candidacy(dock=dock)

        # Plot charts of overlap between sets in selection of genes by permutation
        # probability of bimodality.
        prepare_charts_gene_sets_probability(dock=dock)

        # Plot charts of overlap between sets in selection of genes by
        # heritability.
        prepare_charts_gene_heritability(dock=dock)

        # Plot charts of overlap between sets in selection of genes by
        # integration.
        prepare_charts_gene_sets_integration(dock=dock)



    #plot_charts_analysis(dock=dock)
    #plot_charts_tissue(dock=dock)
    #plot_charts_restriction(dock=dock)
    #plot_chart_tissues_persons(dock=dock)

    # Plot charts of distributions of genes' bimodality measurement scores and
    # permutations.
    #plot_charts_distribution(dock=dock)

    pass


if (__name__ == "__main__"):
    execute_procedure()
