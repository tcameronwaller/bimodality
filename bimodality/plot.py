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
import copy

# Relevant
import pandas
import numpy
import scipy
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams.update({'figure.max_open_warning': 0})
import matplotlib.pyplot
import matplotlib.lines
import matplotlib_venn
import seaborn
import sklearn

# Custom

import assembly
import selection
import integration
import prediction
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
    values_five = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 300,
        "weight": 300,
        "size": 10
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
    properties_five = matplotlib.font_manager.FontProperties(
        family=values_five["family"],
        style=values_five["style"],
        variant=values_five["variant"],
        stretch=values_five["stretch"],
        weight=values_five["weight"],
        size=values_five["size"]
    )
    # Compile and return references.
    return {
        "values": {
            "one": values_one,
            "two": values_two,
            "three": values_three,
            "four": values_four,
            "five": values_five,
        },
        "properties": {
            "one": properties_one,
            "two": properties_two,
            "three": properties_three,
            "four": properties_four,
            "five": properties_five,
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
    # Gray.
    gray = (0.7, 0.7, 0.7, 1.0)
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
        "gray": gray,
        "white": white,
        "white_faint": white_faint,
        "blue": blue,
        "blue_faint": blue_faint,
        "orange": orange,
        "orange_faint": orange_faint
    }


def plot_heatmap_symmetric_diverge(
    data=None,
    minimum=None,
    maximum=None,
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
        minimum (float): minimal value
        maximum (float): maximal value
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
    matrix = numpy.transpose(data.to_numpy())

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = matplotlib.pyplot.axes()

    # Represent data as a color grid.
    image = axes.imshow(
        matrix,
        cmap="PuOr", # "PuOr", "RdBu" diverging color map
        vmin=minimum,
        vmax=maximum,
        aspect="equal",
        origin="upper",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )

    # Create legend for color map.
    label_bar = "correlation"
    bar = axes.figure.colorbar(
        image,
        orientation="vertical",
        ax=axes,
    )
    bar.ax.set_ylabel(
        label_bar,
        rotation=-90,
        va="bottom",
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["three"],
    )
    bar.ax.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=2.0,
        color=colors["black"],
        pad=3,
        labelsize=fonts["values"]["three"]["size"],
        labelcolor=colors["black"],
    )

    # Create ticks and labels for each grid.
    # Let the horizontal axes labeling appear on top.
    axes.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=10,
        labelcolor=colors["black"],
        top=True,
        bottom=False, # False
        left=True,
        right=False,
        labeltop=True,
        labelbottom=False,
        labelleft=True,
        labelright=False,
    )

    # Create ticks and labels.
    if (label_columns and (data.shape[1] <= 50)):
        axes.set_xticks(numpy.arange(matrix.shape[1]))
        axes.set_xticklabels(
            labels_columns,
            #minor=False,
            rotation=-30,
            rotation_mode="anchor",
            ha="right", # horizontal alignment
            va="bottom", # vertical alignment
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )
    if (label_rows and (data.shape[0] <= 50)):
        axes.set_yticks(numpy.arange(matrix.shape[0]))
        axes.set_yticklabels(
            labels_rows,
            #minor=False,
            ha="right", # horizontal alignment
            va="center", # vertical alignment
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["three"]
        )

    # Return figure.
    return figure


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
            fontproperties=fonts["properties"]["five"]
        )

    return figure


def plot_heatmap_asymmetric_groups(
    title=None,
    property=None,
    type=None,
    properties=None,
    data=None,
    fonts=None,
    colors=None,
    title_signal=None,
):
    """
    Creates a figure of a chart of type heatmap.

    Data must have observations across rows.
    Data's observations must already be in sort order.
    Data must have features across columns.
    Data must also have a single column of name "group".
    Values in column "group" must be integers.

    The chart will organize features across rows and observations across
    columns.
    The chart will represent integer values of the group of each observation in
    a separate chart across columns.


    arguments:
        title (str): title for chart
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        properties (list<str>): unique values of property
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        title_signal (str): title for color bar on main chart

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    # Copy data.
    #data_selection = data.iloc[:, :]
    data_property = data.copy(deep=True)
    data_value = data.copy(deep=True)
    # Groups.
    property_values = data_property["property"].to_list()
    # Map data's columns to heatmap's rows.
    # Map data's rows to heatmap's columns.
    data_value.drop(
        labels=["property"],
        axis="columns",
        inplace=True
    )
    labels_rows = data_value.columns.to_list()
    labels_columns = data_value.index.to_list()
    matrix = numpy.transpose(data_value.to_numpy())
    # Determine maximal and minimal values.
    minimum_property = round(min(property_values), 2)
    maximum_property = round(max(property_values), 2)
    #utility.print_terminal_partition(level=2)
    #print(minimum_property)
    #print(maximum_property)
    array = matrix.flatten()
    minimum_value = round(numpy.nanmin(array), 2)
    maximum_value = round(numpy.nanmax(array), 2)
    #utility.print_terminal_partition(level=3)
    #print(minimum_value)
    #print(maximum_value)

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811), # 40 cm, 30 cm
        #tight_layout=True,
    )
    figure.suptitle(
        title,
        x=0.01,
        y=0.99,
        ha="left",
        va="top",
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes = figure.subplots(
        nrows=2,
        ncols=2,
        sharex=False,
        sharey=False,
        #squeeze=True,
        gridspec_kw=dict(
            hspace=0.05,
            wspace=0.05,
            height_ratios=[1, 10],
            width_ratios=[50, 1],
            left=0.12,
            right=0.94,
            top=0.94,
            bottom=0.05,
        ),
    )

    ##########
    # Top axes.
    if (type == "category"):
        #color_map = matplotlib.colors.ListedColormap(
        #    [colors["blue"], colors["orange"]], 2
        #)
        color_map = matplotlib.pyplot.get_cmap("GnBu", len(properties))
        axis_property_labels = properties
    elif type == "binary":
        color_map = matplotlib.pyplot.get_cmap("GnBu", 2)
        axis_property_labels = [minimum_property, maximum_property]
    elif type == "continuity":
        color_map = "GnBu"
        axis_property_labels = [minimum_property, maximum_property]
    image_property = axes[0, 0].imshow(
        [property_values],
        cmap=color_map,
        vmin=minimum_property,
        vmax=maximum_property,
        aspect="auto",
        origin="upper",
        # Extent: (left, right, bottom, top)
        extent=(-0.5, (len(property_values) - 0.5), (1 + 0.5), -0.5),
    )
    axes[0, 0].set_yticks(numpy.arange(1))
    axes[0, 0].set_yticklabels(
        [str(property)],
        #minor=False,
        ha="right", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    axes[0, 0].set_xlabel(
        "persons' property",
        rotation=0,
        ha="right",
        va="bottom",
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["three"],
    )
    #bar_group.ax.xaxis.set_label_position("top")
    axes[0, 0].xaxis.set_label_coords(1.0, 1.2)
    axes[0, 0].tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=3.0,
        pad=7,
        top=False,
        bottom=False,
        left=True,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=True,
        labelright=False,
        #color=colors["black"],
        #labelsize=fonts["values"]["four"]["size"],
        #labelcolor=colors["black"],
    )
    # Create legend for color map.
    bar_property = figure.colorbar(
        image_property,
        cax=axes[0, 1],
        ticks=[minimum_property, maximum_property],
        orientation="vertical",
        use_gridspec=True,
    )
    bar_property.ax.set_yticklabels(
        axis_property_labels,
        #minor=False,
        ha="left", # horizontal alignment
        va="bottom", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    bar_property.ax.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=3.0,
        pad=7,
        top=False,
        bottom=False,
        left=False,
        right=True,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=True,
        #color=colors["black"],
        #labelsize=fonts["values"]["four"]["size"],
        #labelcolor=colors["black"],
    )


    ##########
    # Bottom axes.

    # Represent data as a color grid.
    two_slope_normalization = matplotlib.colors.TwoSlopeNorm(
        vmin=minimum_value,
        vcenter=0.0,
        vmax=maximum_value,
    )
    image_value = axes[1, 0].imshow(
        matrix,
        cmap="PuOr", # sequential: "RdPu", diverging: "PuOr"
        #vmin=minimum_value,
        #vmax=maximum_value,
        aspect="auto", # "auto", "equal"
        origin="upper", # "lower" or "upper"
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
        alpha=1.0,
        norm=two_slope_normalization,
        filternorm=True,
        resample=True,
    )
    #axes[2, 0].set_xlim(-10, (matrix.shape[1] + 10))
    #axes[2, 0].set_ylim(-3, (matrix.shape[0] + 3))

    # Create ticks and labels for each grid.
    # Let the horizontal axes labeling appear on top.
    if matrix.shape[0] > 25:
        size_count = "five"
    elif matrix.shape[0] <= 25:
        size_count = "four"
    axes[1, 0].tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        pad=7,
        top=False,
        bottom=True,
        left=True,
        right=False,
        labeltop=False,
        labelbottom=True,
        labelleft=True,
        labelright=False,
        color=colors["black"],
        labelsize=fonts["values"][size_count]["size"],
        labelcolor=colors["black"],
    )
    axes[1, 0].set_yticks(numpy.arange(matrix.shape[0]))
    axes[1, 0].set_yticklabels(
        labels_rows,
        #minor=False,
        ha="right", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"][size_count]
    )
    # Create legend for color map.
    bar_value = figure.colorbar(
        image_value,
        cax=axes[1, 1],
        ticks=[minimum_value, maximum_value],
        orientation="vertical",
        use_gridspec=True,
    )
    bar_value.ax.set_ylabel(
        title_signal,
        rotation=-90,
        ha="center",
        va="bottom",
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["three"],
    )
    bar_value.ax.yaxis.set_label_coords(2, 0.5)
    bar_value.ax.set_yticklabels(
        [minimum_value, maximum_value],
        #minor=False,
        ha="left", # horizontal alignment
        va="top", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    bar_value.ax.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=3.0,
        pad=7,
        top=False,
        bottom=False,
        left=False,
        right=True,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=True,
        #color=colors["black"],
        #labelsize=fonts["values"]["four"]["size"],
        #labelcolor=colors["black"],
    )
    # Return figure.
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
    if False:
        axes.set_title(
            name,
            #fontdict=fonts["properties"]["one"],
            loc="center",
        )
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
            linewidth=5.0, # 3.0, 7.5
        )
    if len(text) > 0:
        matplotlib.pyplot.text(
            1,
            1,
            text,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["one"]
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


def plot_scatter(
    data=None,
    abscissa=None,
    ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    fonts=None,
    colors=None,
    size=None,
):
    """
    Creates a figure of a chart of type scatter.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with independent variable
        ordinate (str): name of data column with dependent variable
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        factor (str): name of data column with groups or factors of samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        size (int): size of marker

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    data_copy = data.copy(deep=True)
    data_selection = data_copy.loc[:, [abscissa, ordinate]]
    data_selection.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    values_abscissa = data_selection[abscissa].values
    values_ordinate = data_selection[ordinate].values

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlabel(
        xlabel=title_abscissa,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=title_ordinate,
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
    handle = axes.plot(
        values_abscissa,
        values_ordinate,
        linestyle="",
        marker="o",
        markersize=size, # 5, 15
        markeredgecolor=colors["blue"],
        markerfacecolor=colors["blue"]
    )

    # Return figure.
    return figure


def plot_scatter_threshold(
    data=None,
    abscissa=None,
    ordinate=None,
    threshold_abscissa=None,
    selection_abscissa=None,
    threshold_ordinate=None,
    selection_ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type scatter with thresholds on each
        dimension.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with independent variable
        ordinate (str): name of data column with dependent variable
        threshold_abscissa (float): threshold for abscissa
        selection_abscissa (str): selection criterion for abscissa's values
            against threshold
        threshold_ordinate (float): threshold for ordinate
        selection_ordinate (str): selection criterion for ordinate's values
            against threshold
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        factor (str): name of data column with groups or factors of samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    data = data.copy(deep=True)
    data = data.loc[:, [abscissa, ordinate]]
    data.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    # Divide values by whether they pass thresholds on both dimensions.
    collection = utility.segregate_data_two_thresholds(
        data=data,
        abscissa=abscissa,
        ordinate=ordinate,
        threshold_abscissa=threshold_abscissa,
        selection_abscissa=selection_abscissa,
        threshold_ordinate=threshold_ordinate,
        selection_ordinate=selection_ordinate,
    )

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlabel(
        xlabel=title_abscissa,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=title_ordinate,
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
    handle = axes.plot(
        collection["fail"][abscissa].values,
        collection["fail"][ordinate].values,
        linestyle="",
        marker="o",
        markersize=2.5,
        markeredgecolor=colors["gray"],
        markerfacecolor=colors["gray"]
    )
    handle = axes.plot(
        collection["pass"][abscissa].values,
        collection["pass"][ordinate].values,
        linestyle="",
        marker="o",
        markersize=5,
        markeredgecolor=colors["blue"],
        markerfacecolor=colors["blue"]
    )

    # Plot lines for each threshold value...
    # Create lines for thresholds.
    axes.axvline(
        x=threshold_abscissa,
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["orange"],
        linestyle="--",
        linewidth=3.0,
    )
    axes.axhline(
        y=threshold_ordinate,
        xmin=0,
        xmax=1,
        alpha=1.0,
        color=colors["orange"],
        linestyle="--",
        linewidth=3.0,
    )

    # Return figure.
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


def write_figure_png(path=None, figure=None):
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
        format="png",
        dpi=300,
        facecolor="w",
        edgecolor="w",
        transparent=False
    )
    pass


##################################################
# Procedures


##########
# Selection of health history variables for persons respiration, inflammation,
# infection, and steroid.
# heatmaps
# Status: working


def read_source_persons_health_variables(
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
    path_collections_health = os.path.join(
        dock, "assembly", "sample", "collections_health_variables.pickle"
    )
    path_persons_selection = os.path.join(
        dock, "selection", "tight", "samples_genes_signals", "persons.pickle"
    )
    path_persons_sets = os.path.join(
        dock, "selection", "tight", "persons_properties", "charts",
        "persons_sets.pickle"
    )
    path_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", "regression",
        "data_persons_properties.pickle"
    )

    # Read information from file.
    with open(path_collections_health, "rb") as file_source:
        collections_health_variables = pickle.load(file_source)
    with open(path_persons_selection, "rb") as file_source:
        persons_selection = pickle.load(file_source)
    with open(path_persons_sets, "rb") as file_source:
        persons_sets = pickle.load(file_source)
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    # Compile and return information.
    return {
        "collections_health_variables": collections_health_variables,
        "persons_selection": persons_selection,
        "persons_sets": persons_sets,
        "data_persons_properties": data_persons_properties,
    }


def organize_persons_health_variables(
    property=None,
    type=None,
    sequence=None,
    persons=None,
    data_persons_properties=None,
    health_collection=None,
    data_health_collection=None,
):
    """
    Organize information for chart.

    Notice that the data have features (health variables) across columns and
    instances (persons) across rows.

    Sequence of genes across rows depends on hierarchical cluster by their
    similarities across persons.
    Sequence of persons across columns depends either on sort by values of
    property or on hierarchical cluster by their similarities across genes.

    arguments:
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        sequence (str): method for sequence of persons, sort by property's
            values or cluster by similarities across genes
        persons (list<str>): identifiers of persons
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        health_collection (str): name of collection of health variables
        data_health_collection (object): Pandas data frame of health variables
            across persons

    raises:

    returns:
        (dict): information for chart

    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    data_health_collection = data_health_collection.copy(deep=True)
    # Organize health variables.
    # Select data for persons.
    data_health_selection = data_health_collection.loc[
        data_health_collection.index.isin(persons), :
    ]
    # Replace missing values with zero.
    data_health_selection.fillna(
        value=0.0,
        #axis="columns",
        inplace=True,
    )
    # Cluster data.
    # Cluster features.
    data_health_cluster = utility.cluster_data_columns(
        data=data_health_selection,
    )
    # Cluster instances.
    data_health_sequence = utility.cluster_data_rows(
        data=data_health_cluster,
    )
    # Organize properties.
    if type == "category":
        data_persons_properties["property"], properties = pandas.factorize(
            data_persons_properties[property],
            sort=True
        )
        properties = list(properties)
    elif type == "continuity":
        data_persons_properties["property"] = data_persons_properties[property]
        properties = None
    data_persons_properties = data_persons_properties.loc[
        :, data_persons_properties.columns.isin(["property", property])
    ]
    data_hybrid = data_health_sequence.join(
        data_persons_properties,
        how="left",
        on="person"
    )
    # Determine whether to sort by persons' values of property.
    if sequence == "sort":
        data_hybrid.sort_values(
            by=["property"],
            axis="index",
            ascending=True,
            inplace=True,
        )
    # Remove the column for the named property.
    data_hybrid.drop(
        labels=[property],
        axis="columns",
        inplace=True
    )
    # Compile information.
    bin = dict()
    bin["properties"] = properties
    bin["data"] = data_hybrid
    # Return information
    return bin


def plot_chart_persons_health_variables(
    name=None,
    property=None,
    type=None,
    properties=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        name (str): name of file and chart
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        properties (list<str>): unique values of property
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(name + ".png")
    )
    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Create figure.
    figure = plot_heatmap_asymmetric_groups(
        title=name,
        property=property,
        type=type,
        properties=properties,
        data=data,
        fonts=fonts,
        colors=colors,
        title_signal="persons' health history variables",
    )
    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_persons_health_variables_variable(
    property=None,
    type=None,
    persons=None,
    data_persons_properties=None,
    health_collection=None,
    data_health_collection=None,
    path_sort=None,
    path_cluster=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        persons (list<str>): identifiers of persons
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        health_collection (str): name of collection of health variables
        data_health_collection (object): Pandas data frame of health variables
            across persons
        path_sort (str): path to directory
        path_cluster (str): path to directory

    raises:

    returns:

    """

    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_persons_health_variables(
        property=property,
        type=type,
        sequence="sort",
        persons=persons,
        data_persons_properties=data_persons_properties,
        health_collection=health_collection,
        data_health_collection=data_health_collection,
    )
    # Create charts for the gene.
    plot_chart_persons_health_variables(
        name=str(health_collection + "_" + property),
        property=property,
        type=type,
        properties=bin["properties"],
        data=bin["data"],
        path_directory=path_sort,
    )

    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_persons_health_variables(
        property=property,
        type=type,
        sequence="cluster",
        persons=persons,
        data_persons_properties=data_persons_properties,
        health_collection=health_collection,
        data_health_collection=data_health_collection,
    )
    # Create charts for the gene.
    plot_chart_persons_health_variables(
        name=str(health_collection + "_" + property),
        property=property,
        type=type,
        properties=bin["properties"],
        data=bin["data"],
        path_directory=path_cluster,
    )
    pass


def prepare_charts_persons_health_variables(
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
    source = read_source_persons_health_variables(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_assembly = os.path.join(path_plot, "assembly")
    path_persons_health = os.path.join(
        path_assembly, "persons_health_variables"
    )
    path_sort = os.path.join(
        path_persons_health, "sort_persons"
    )
    path_cluster = os.path.join(
        path_persons_health, "cluster_persons"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_persons_health)
    utility.create_directories(path=path_sort)
    utility.create_directories(path=path_cluster)

    # Iterate on categorical and ordinal groups of properties.
    properties = list()
    properties.append(dict(name="sex_text", type="category"))
    properties.append(dict(name="age", type="continuity"))
    properties.append(dict(name="body", type="continuity"))
    properties.append(dict(name="hardiness", type="continuity"))
    properties.append(dict(name="climate", type="continuity"))
    properties.append(dict(name="ventilation", type="category"))
    properties.append(dict(name="refrigeration", type="category"))
    # Iterate on collections of health variables.
    for collection in source["collections_health_variables"]:
        # Iterate on categorical variables.
        for property in properties:
            # Prepare charts for genes.
            prepare_charts_persons_health_variables_variable(
                property=property["name"],
                type=property["type"],
                persons=source["persons_sets"]["selection"],
                data_persons_properties=source["data_persons_properties"],
                health_collection=collection,
                data_health_collection=(
                    source["collections_health_variables"][collection]
                ),
                path_sort=path_sort,
                path_cluster=path_cluster,
            )
            pass
    pass


##########
# Selection of health history variables for persons respiration, inflammation,
# infection, and steroid.
# heatmaps
# Status: working


def read_source_persons_properties_adjacency(
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
    path_persons_selection = os.path.join(
        dock, "selection", "tight", "samples_genes_signals", "persons.pickle"
    )
    path_persons_sets = os.path.join(
        dock, "selection", "tight", "persons_properties", "charts",
        "persons_sets.pickle"
    )
    path_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", "regression",
        "data_persons_properties.pickle"
    )

    # Read information from file.
    with open(path_persons_selection, "rb") as file_source:
        persons_selection = pickle.load(file_source)
    with open(path_persons_sets, "rb") as file_source:
        persons_sets = pickle.load(file_source)
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    # Compile and return information.
    return {
        "persons_selection": persons_selection,
        "persons_sets": persons_sets,
        "data_persons_properties": data_persons_properties,
    }


def organize_persons_properties_adjacency(
    master=None,
    properties=None,
    sequence=None,
    persons=None,
    data_persons_properties=None,
):
    """
    Organize information for chart.

    Notice that the data have features (health variables) across columns and
    instances (persons) across rows.

    Sequence of genes across rows depends on hierarchical cluster by their
    similarities across persons.
    Sequence of persons across columns depends either on sort by values of
    property or on hierarchical cluster by their similarities across genes.

    arguments:
        master (str): name of feature from persons' properties to use for
            groups
        properties (dict<str>): names and types of persons' properties,
            category or continuity
        sequence (str): method for sequence of persons, sort by property's
            values or cluster by similarities across genes
        persons (list<str>): identifiers of persons
        data_persons_properties (object): Pandas data frame of persons and
            their properties

    raises:

    returns:
        (dict): information for chart

    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    # Organize health variables.
    # Select data for persons.
    data_persons = data_persons_properties.loc[
        data_persons_properties.index.isin(persons), :
    ]
    data_properties = data_persons.loc[
        :, data_persons.columns.isin(properties.keys())
    ]
    # Scale feature values.
    data_scale = data_properties.apply(
        lambda column:
            sklearn.preprocessing.minmax_scale(
                column.to_numpy(),
                feature_range=(0, 1),
                axis=0,
                copy=True,
            ) if (properties[column.name]["type"] == "continuity") else column,
        axis="index",
    )
    # Replace missing values with zero.
    data_scale.fillna(
        value=0.0,
        #axis="columns",
        inplace=True,
    )
    # Cluster data.
    # Cluster features.
    data_cluster = utility.cluster_data_columns(
        data=data_scale,
    )
    # Cluster instances.
    data_sequence = utility.cluster_data_rows(
        data=data_cluster,
    )
    # Organize properties.
    if properties[master]["type"] == "category":
        data_sequence["property"], values_category = pandas.factorize(
            data_sequence[master],
            sort=True
        )
        values_category = list(values_category)
    else:
        data_sequence["property"] = data_sequence[master]
        values_category = None
    # Determine whether to sort by persons' values of property.
    if sequence == "sort":
        data_sequence.sort_values(
            by=["property"],
            axis="index",
            ascending=True,
            inplace=True,
        )
    # Remove the column for the named property.
    data_sequence.drop(
        labels=[master],
        axis="columns",
        inplace=True
    )
    # Prepare translations of properties' names.
    translations = dict()
    for property in properties.keys():
        translations[property] = properties[property]["name"]
        pass
    # Rename genes in data.
    data_sequence.rename(
        columns=translations,
        inplace=True,
    )

    # Compile information.
    bin = dict()
    bin["values_category"] = values_category
    bin["data"] = data_sequence
    # Return information
    return bin


def plot_chart_persons_properties_adjacency(
    name=None,
    property=None,
    type=None,
    properties=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        name (str): name of file and chart
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        properties (list<str>): unique values of property
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(name + ".png")
    )
    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Create figure.
    figure = plot_heatmap_asymmetric_groups(
        title=name,
        property=property,
        type=type,
        properties=properties,
        data=data,
        fonts=fonts,
        colors=colors,
        title_signal="persons' properties",
    )
    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_persons_properties_adjacency_property(
    master=None,
    properties=None,
    persons=None,
    data_persons_properties=None,
    path_sort=None,
    path_cluster=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        master (str): name of feature from persons' properties to use for
            groups
        properties (dict<str>): names and types of persons' properties,
            category or continuity
        persons (list<str>): identifiers of persons
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        path_sort (str): path to directory
        path_cluster (str): path to directory

    raises:

    returns:

    """

    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_persons_properties_adjacency(
        master=master,
        properties=properties,
        sequence="sort",
        persons=persons,
        data_persons_properties=data_persons_properties,
    )
    # Create charts for the gene.
    plot_chart_persons_properties_adjacency(
        name=properties[master]["name"],
        property=properties[master]["name"],
        type=properties[master]["type"],
        properties=bin["values_category"],
        data=bin["data"],
        path_directory=path_sort,
    )

    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_persons_properties_adjacency(
        master=master,
        properties=properties,
        sequence="cluster",
        persons=persons,
        data_persons_properties=data_persons_properties,
    )
    # Create charts for the gene.
    plot_chart_persons_properties_adjacency(
        name=properties[master]["name"],
        property=properties[master]["name"],
        type=properties[master]["type"],
        properties=bin["values_category"],
        data=bin["data"],
        path_directory=path_cluster,
    )
    pass


def prepare_charts_persons_properties_adjacency(
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
    source = read_source_persons_properties_adjacency(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_assembly = os.path.join(path_plot, "assembly")
    path_persons_properties = os.path.join(
        path_assembly, "persons_properties"
    )
    path_sort = os.path.join(
        path_persons_properties, "sort_persons"
    )
    path_cluster = os.path.join(
        path_persons_properties, "cluster_persons"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_persons_properties)
    utility.create_directories(path=path_sort)
    utility.create_directories(path=path_cluster)

    # Define persons of interest.
    persons = source["persons_sets"]["selection"]
    #persons = source["persons_sets"]["non_ventilation"]

    # Iterate on categorical and ordinal groups of properties.
    properties = dict()
    properties["sex_y"] = dict(name="sex_y", type="binary")
    properties["age"] = dict(name="age", type="continuity")
    properties["body"] = dict(name="bmi", type="continuity")
    properties["hardiness"] = dict(name="hardy", type="continuity")
    properties["climate"] = dict(name="climate", type="continuity")
    properties["respiration_binary"] = dict(name="respiration", type="binary")
    properties["inflammation_binary"] = dict(
        name="inflammation", type="binary"
    )
    properties["infection_binary"] = dict(name="infection", type="binary")
    properties["mononucleosis_binary"] = dict(
        name="mononucleosis", type="binary"
    )
    properties["steroid_binary"] = dict(name="steroid", type="binary")
    properties["ventilation_binary"] = dict(name="ventilation", type="binary")
    #properties["ventilation_duration"] = dict(
    #    name="vent_dur", type="continuity"
    #)
    properties["delay"] = dict(name="delay", type="continuity")
    properties["refrigeration_binary"] = dict(
        name="refrigeration", type="binary"
    )
    #properties["refrigeration_duration"] = dict(
    #    name="refrig_dur", type="continuity"
    #)
    # Iterate on categorical variables.
    for property in properties.keys():
        # Prepare charts for genes.
        prepare_charts_persons_properties_adjacency_property(
            master=property,
            properties=properties,
            persons=persons,
            data_persons_properties=source["data_persons_properties"],
            path_sort=path_sort,
            path_cluster=path_cluster,
        )
        pass
    pass





##########
# Sex, age, and hardiness of persons in GTEx cohort
# Status: in progress


def read_source_persons_sex_age_hardiness(dock=None):
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
    path_selection = os.path.join(dock, "selection", "tight")
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


def prepare_charts_persons_sex_age_hardiness(
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
    source = read_source_persons_sex_age_hardiness(dock=dock)

    if False:

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
# Counts of batches per person, histogram with bins by counts
# Status: in progress


def read_source_batches_per_person(dock=None):
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_data_persons_properties_raw = os.path.join(
        path_selection, "data_persons_properties_raw.pickle"
    )
    # Read information from file.
    with open(path_data_persons_properties_raw, "rb") as file_source:
        data_persons_properties_raw = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_persons_properties_raw": data_persons_properties_raw,
    }


def prepare_chart_batches_per_person(dock=None):
    """
    Plots charts from the sample process.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    def split_batches(characters):
        values_split = characters.split(",")
        count = len(values_split)
        return count

    # Read source information from file.
    source = read_source_batches_per_person(dock=dock)
    print(source["data_persons_properties_raw"])

    # Organize data.
    data = source["data_persons_properties_raw"].copy(deep=True)
    data["count_batches_isolation"] = data["batches_isolation"].map(split_batches)
    counts = data["count_batches_isolation"].to_list()

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
        series=counts,
        name="",
        bin_method="count",
        bin_count=20,
        label_bins="Bins by count of batches per person",
        label_counts="Counts of persons in each bin",
        fonts=fonts,
        colors=colors,
        line=False,
        position=11,
        text="",
    )
    # Specify directories and files.
    file = ("batches_per_person.svg")
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
    path_selection = os.path.join(dock, "selection", "tight")
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
    path_selection = os.path.join(dock, "selection", "tight")
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
        bin_count=10,
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
# Overlap in sets from selection of genes and samples
# Status: working


def read_source_sets_selection_genes_samples(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_genes_gtex = os.path.join(
        path_selection, "genes_gtex.pickle"
    )
    path_genes_selection = os.path.join(
        path_selection, "genes_selection.pickle"
    )
    path_samples_gtex = os.path.join(
        path_selection, "samples_gtex.pickle"
    )
    path_samples_selection = os.path.join(
        path_selection, "samples_selection.pickle"
    )

    # Read information from file.
    with open(path_genes_gtex, "rb") as file_source:
        genes_gtex = pickle.load(file_source)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_samples_gtex, "rb") as file_source:
        samples_gtex = pickle.load(file_source)
    with open(path_samples_selection, "rb") as file_source:
        samples_selection = pickle.load(file_source)

    # Compile and return information.
    return {
        "genes_gtex": genes_gtex,
        "genes_selection": genes_selection,
        "samples_gtex": samples_gtex,
        "samples_selection": samples_selection,
    }


def plot_chart_sets_selection_genes_samples(
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
        count=2,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_charts_sets_selection_genes_samples(
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

    # Read source information from file.
    source = read_source_sets_selection_genes_samples(dock=dock)

    # Organize information.
    sets_genes = dict()
    sets_genes["gtex"] = source["genes_gtex"]
    sets_genes["selection"] = source["genes_selection"]

    sets_samples = dict()
    sets_samples["gtex"] = source["samples_gtex"]
    sets_samples["selection"] = source["samples_selection"]

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_selection = os.path.join(path_plot, "selection")
    path_set = os.path.join(path_selection, "set")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_set)
    utility.create_directories(path=path_set)
    path_figure_genes = os.path.join(path_set, "genes.svg")
    path_figure_samples = os.path.join(path_set, "samples.svg")

    plot_chart_sets_selection_genes_samples(
        sets=sets_genes,
        path=path_figure_genes,
    )
    plot_chart_sets_selection_genes_samples(
        sets=sets_samples,
        path=path_figure_samples,
    )

    pass


##########
# Overlap in sets from selection of persons and their properties
# Status: in progress


def read_source_sets_selection_persons_properties(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_persons_properties = os.path.join(
        path_selection, "persons_properties"
    )
    path_charts = os.path.join(
        path_persons_properties, "charts"
    )

    path_persons_sets = os.path.join(
        path_charts, "persons_sets.pickle"
    )

    # Read information from file.
    with open(path_persons_sets, "rb") as file_source:
        persons_sets = pickle.load(file_source)

    # Compile and return information.
    return {
        "persons_sets": persons_sets,
    }


def plot_chart_sets_selection_persons_properties(
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
        count=len(sets.keys()),
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_charts_sets_selection_persons_properties(
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

    # Read source information from file.
    source = read_source_sets_selection_persons_properties(dock=dock)

    # Organize information.
    sets_one = dict()
    sets_one["gtex"] = source["persons_sets"]["original"]
    sets_one["selection"] = source["persons_sets"]["selection"]
    sets_one["genotype"] = source["persons_sets"]["genotype"]

    sets_two = dict()
    sets_two["ventilation"] = source["persons_sets"]["ventilation"]
    sets_two["respiration"] = source["persons_sets"]["respiration"]
    sets_two["inflammation"] = source["persons_sets"]["inflammation"]

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_selection = os.path.join(path_plot, "selection")
    path_directory = os.path.join(path_selection, "persons_sets")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    path_figure_one = os.path.join(path_directory, "one.svg")
    path_figure_two = os.path.join(path_directory, "two.svg")

    plot_chart_sets_selection_persons_properties(
        sets=sets_one,
        path=path_figure_one,
    )
    plot_chart_sets_selection_persons_properties(
        sets=sets_two,
        path=path_figure_two,
    )

    pass


##########
# Selection of principal components for regression
# Status: working


def read_source_regression_principal_components(
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
    path_assembly = os.path.join(dock, "assembly")
    path_data_genotype_variance_plink = os.path.join(
        path_assembly, "data_genotype_variance_plink.pickle"
    )

    path_selection = os.path.join(dock, "selection", "tight")
    path_data_tissues_variance = os.path.join(
        path_selection, "data_tissues_variance.pickle"
    )
    path_data_facilities_variance = os.path.join(
        path_selection, "data_facilities_variance.pickle"
    )
    path_data_batches_isolation_variance = os.path.join(
        path_selection, "data_batches_isolation_variance.pickle"
    )
    path_data_batches_sequence_variance = os.path.join(
        path_selection, "data_batches_sequence_variance.pickle"
    )

    # Read information from file.
    data_genotype_variance_plink = pandas.read_pickle(
        path_data_genotype_variance_plink
    )
    data_tissues_variance = pandas.read_pickle(
        path_data_tissues_variance
    )
    data_facilities_variance = pandas.read_pickle(
        path_data_facilities_variance
    )
    data_batches_isolation_variance = pandas.read_pickle(
        path_data_batches_isolation_variance
    )
    data_batches_sequence_variance = pandas.read_pickle(
        path_data_batches_sequence_variance
    )

    # Compile and return information.
    return {
        "data_genotype_variance_plink": data_genotype_variance_plink,
        "data_tissues_variance": data_tissues_variance,
        "data_facilities_variance": data_facilities_variance,
        "data_batches_isolation_variance": data_batches_isolation_variance,
        "data_batches_sequence_variance": data_batches_sequence_variance,
    }


def plot_chart_regression_principal_components(
    data=None,
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        data (object): Pandas data frame of variables
        path_file (str): path for file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_scatter(
        data=data,
        abscissa="component",
        ordinate="variance",
        title_abscissa="Principal components",
        title_ordinate="Variances",
        fonts=fonts,
        colors=colors,
        size=15,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_regression_principal_components(
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
    source = read_source_regression_principal_components(dock=dock)

    # Organize data.
    types = dict()
    types["component"] = "int32"
    types["variance"] = "float32"

    data_genotype = source["data_genotype_variance_plink"].iloc[0:25,:]
    data_genotype = data_genotype.astype(
        types,
    )

    data_tissues = source["data_tissues_variance"]
    data_tissues["variance"] = data_tissues["variance"].apply(
        lambda value: (value * 100)
    )
    data_tissues = data_tissues.astype(
        types,
    )

    data_facilities = source["data_facilities_variance"]
    data_facilities["variance"] = data_facilities["variance"].apply(
        lambda value: (value * 100)
    )
    data_facilities = data_facilities.astype(
        types,
    )

    data_isolation = source["data_batches_isolation_variance"]
    data_isolation["variance"] = data_isolation["variance"].apply(
        lambda value: (value * 100)
    )
    data_isolation = data_isolation.astype(
        types,
    )

    data_sequence = source["data_batches_sequence_variance"]
    data_sequence["variance"] = data_sequence["variance"].apply(
        lambda value: (value * 100)
    )
    data_sequence = data_sequence.astype(
        types,
    )

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_selection = os.path.join(path_plot, "selection")
    path_component = os.path.join(path_selection, "component")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_component)
    utility.create_directories(path=path_component)

    plot_chart_regression_principal_components(
        data=data_genotype,
        path_file=os.path.join(path_component, "genotype.svg"),
    )
    plot_chart_regression_principal_components(
        data=data_tissues,
        path_file=os.path.join(path_component, "tissues.svg"),
    )
    plot_chart_regression_principal_components(
        data=data_facilities,
        path_file=os.path.join(path_component, "facilities.svg"),
    )
    plot_chart_regression_principal_components(
        data=data_isolation,
        path_file=os.path.join(path_component, "batches_isolation.svg"),
    )
    plot_chart_regression_principal_components(
        data=data_sequence,
        path_file=os.path.join(path_component, "batches_sequence.svg"),
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

    # Set measures of modality.
    measures = list(source["scores"].keys())
    # Iterate on bimodality measures.
    for measure in measures:

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
# Thresholds on heritability for selection of genes
# Status: working


def read_source_gene_heritability(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )

    path_heritability = os.path.join(dock, "heritability")
    path_collection = os.path.join(path_heritability, "collection")

    path_data_genes_heritabilities_simple = os.path.join(
        path_collection, "data_genes_heritabilities_simple.pickle"
    )
    path_data_genes_heritabilities_complex = os.path.join(
        path_collection, "data_genes_heritabilities_complex.pickle"
    )
    path_data_genes_unimodal_heritabilities_simple = os.path.join(
        path_collection, "data_genes_unimodal_heritabilities_simple.pickle"
    )
    path_data_genes_unimodal_heritabilities_complex = os.path.join(
        path_collection, "data_genes_unimodal_heritabilities_complex.pickle"
    )
    path_data_genes_multimodal_heritabilities_simple = os.path.join(
        path_collection, "data_genes_multimodal_heritabilities_simple.pickle"
    )
    path_data_genes_multimodal_heritabilities_complex = os.path.join(
        path_collection, "data_genes_multimodal_heritabilities_complex.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_genes_heritabilities_simple = pandas.read_pickle(
        path_data_genes_heritabilities_simple
    )
    data_genes_heritabilities_complex = pandas.read_pickle(
        path_data_genes_heritabilities_complex
    )
    data_genes_unimodal_heritabilities_simple = pandas.read_pickle(
        path_data_genes_unimodal_heritabilities_simple
    )
    data_genes_unimodal_heritabilities_complex = pandas.read_pickle(
        path_data_genes_unimodal_heritabilities_complex
    )
    data_genes_multimodal_heritabilities_simple = pandas.read_pickle(
        path_data_genes_multimodal_heritabilities_simple
    )
    data_genes_multimodal_heritabilities_complex = pandas.read_pickle(
        path_data_genes_multimodal_heritabilities_complex
    )

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_genes_heritabilities_simple": data_genes_heritabilities_simple,
        "data_genes_heritabilities_complex": data_genes_heritabilities_complex,
        "data_genes_unimodal_heritabilities_simple": (
            data_genes_unimodal_heritabilities_simple
        ),
        "data_genes_unimodal_heritabilities_complex": (
            data_genes_unimodal_heritabilities_complex
        ),
        "data_genes_multimodal_heritabilities_simple": (
            data_genes_multimodal_heritabilities_simple
        ),
        "data_genes_multimodal_heritabilities_complex": (
            data_genes_multimodal_heritabilities_complex
        ),
    }


def organize_gene_heritability(
    data_gene_annotation=None,
    data_genes_heritabilities=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities

    raises:

    returns:
        (dict): information about genes' heritabilities

    """

    # Organize data.
    data = data_genes_heritabilities.copy(deep=True)
    # Calculate percentages.
    data["percentage"] = data["proportion"].apply(
        lambda value: (value * 100)
    )
    data = data.loc[
        :, data.columns.isin([
            "percentage", "discovery_log",
        ])
    ]
    data.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    # Prepare translation of genes' identifiers to names.
    identifiers = data.index.to_list()
    translations = dict()
    for identifier in identifiers:
        translations[identifier] = assembly.access_gene_name(
            identifier=identifier,
            data_gene_annotation=data_gene_annotation,
        )
        pass
    # Rename genes in data.
    data.rename(
        index=translations,
        inplace=True,
    )
    # Return information.
    return data


def plot_chart_gene_heritability(
    data=None,
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        data (object): Pandas data frame of variables
        path_file (str): path for file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_scatter_threshold(
        data=data,
        abscissa="discovery_log",
        ordinate="percentage",
        title_abscissa="-1 * log10(FDR)",
        title_ordinate="% heritability",
        threshold_abscissa=1.30,
        selection_abscissa=">=",
        threshold_ordinate=10,
        selection_ordinate=">=",
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_gene_heritability(
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
    source = read_source_gene_heritability(dock=dock)

    # Define parameters for each chart.
    charts = list()
    record = {
        "title": "genes_selection_simple",
        "data": source["data_genes_heritabilities_simple"],
    }
    charts.append(record)
    record = {
        "title": "genes_selection_complex",
        "data": source["data_genes_heritabilities_complex"],
    }
    charts.append(record)
    record = {
        "title": "genes_unimodal_simple",
        "data": source["data_genes_unimodal_heritabilities_simple"],
    }
    charts.append(record)
    record = {
        "title": "genes_unimodal_complex",
        "data": source["data_genes_unimodal_heritabilities_complex"],
    }
    charts.append(record)
    record = {
        "title": "genes_multimodal_simple",
        "data": source["data_genes_multimodal_heritabilities_simple"],
    }
    charts.append(record)
    record = {
        "title": "genes_multimodal_complex",
        "data": source["data_genes_multimodal_heritabilities_complex"],
    }
    charts.append(record)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_heritability = os.path.join(path_plot, "heritability")
    path_directory = os.path.join(
        path_heritability, "thresholds"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Iterate on charts.
    for chart in charts:
        # Organize data.
        data = organize_gene_heritability(
            data_gene_annotation=source["data_gene_annotation"],
            data_genes_heritabilities=chart["data"],
        )
        # Create chart.
        path_file = os.path.join(
            path_directory, str(chart["title"] + ".svg")
        )
        plot_chart_gene_heritability(
            data=data,
            path_file=path_file
        )

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
    path_sets_unimodal = os.path.join(
        path_candidacy, "sets_unimodal.pickle"
    )
    path_sets_multimodal = os.path.join(
        path_candidacy, "sets_multimodal.pickle"
    )

    # Read information from file.
    with open(path_sets_unimodal, "rb") as file_source:
        sets_unimodal = pickle.load(file_source)
    with open(path_sets_multimodal, "rb") as file_source:
        sets_multimodal = pickle.load(file_source)

    # Compile and return information.
    return {
        "sets_unimodal": sets_unimodal,
        "sets_multimodal": sets_multimodal,
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
    path_figure_unimodal = os.path.join(path_set, "sets_unimodal.svg")
    path_figure_multimodal = os.path.join(path_set, "sets_multimodal.svg")

    plot_chart_gene_sets_candidacy(
        sets=source["sets_unimodal"],
        path=path_figure_unimodal,
    )
    plot_chart_gene_sets_candidacy(
        sets=source["sets_multimodal"],
        path=path_figure_multimodal,
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


def read_source_sets_gene_heritability(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_genes_selection = os.path.join(
        path_selection, "genes_selection.pickle"
    )

    path_heritability = os.path.join(dock, "heritability")
    path_collection = os.path.join(path_heritability, "collection")
    path_sets_genes_selection = os.path.join(
        path_collection, "sets_genes_selection.pickle"
    )
    path_sets_genes_unimodal = os.path.join(
        path_collection, "sets_genes_unimodal.pickle"
    )
    path_sets_genes_multimodal = os.path.join(
        path_collection, "sets_genes_multimodal.pickle"
    )
    # Read information from file.
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_sets_genes_selection, "rb") as file_source:
        sets_genes_selection = pickle.load(file_source)
    with open(path_sets_genes_unimodal, "rb") as file_source:
        sets_genes_unimodal = pickle.load(file_source)
    with open(path_sets_genes_multimodal, "rb") as file_source:
        sets_genes_multimodal = pickle.load(file_source)

    # Compile and return information.
    return {
        "genes_selection": genes_selection,
        "sets_genes_selection": sets_genes_selection,
        "sets_genes_unimodal": sets_genes_unimodal,
        "sets_genes_multimodal": sets_genes_multimodal,
    }


def plot_chart_sets_gene_heritability(
    sets=None,
    count=None,
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        sets (dict<list<str>>): values in sets
        count (int): count of sets
        path_file (str): path to directory and file

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
        count=count,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_sets_gene_heritability(
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

    # Read source information from file.
    source = read_source_sets_gene_heritability(dock=dock)

    # Organize information.
    # Define parameters for each chart.
    charts = list()
    sets_selection = source["sets_genes_selection"]
    sets_selection["total"] = source["genes_selection"]
    record = {
        "title": "genes_selection",
        "count": 3,
        "sets": sets_selection,
    }
    charts.append(record)
    record = {
        "title": "genes_unimodal",
        "count": 2,
        "sets": source["sets_genes_unimodal"],
    }
    charts.append(record)
    record = {
        "title": "genes_multimodal",
        "count": 2,
        "sets": source["sets_genes_multimodal"],
    }
    charts.append(record)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_heritability = os.path.join(path_plot, "heritability")
    path_directory = os.path.join(
        path_heritability, "sets"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Iterate on charts.
    for chart in charts:
        # Create chart.
        path_file = os.path.join(
            path_directory, str(chart["title"] + ".svg")
        )
        plot_chart_sets_gene_heritability(
            sets=chart["sets"],
            count=chart["count"],
            path_file=path_file,
        )

    pass


##########
# Overlap between sets in selection of genes by permutation
# Status: working


def read_source_sets_gene_permutation(
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
    path_probability = os.path.join(dock, "probability")
    path_sets_genes_measures = os.path.join(
        path_probability, "sets_genes_measures.pickle"
    )
    path_sets_genes_selection = os.path.join(
        path_probability, "sets_genes_selection.pickle"
    )
    path_sets_genes_unimodal = os.path.join(
        path_probability, "sets_genes_unimodal.pickle"
    )
    path_sets_genes_multimodal = os.path.join(
        path_probability, "sets_genes_multimodal.pickle"
    )

    # Read information from file.
    with open(path_sets_genes_measures, "rb") as file_source:
        sets_genes_measures = pickle.load(file_source)
    with open(path_sets_genes_selection, "rb") as file_source:
        sets_genes_selection = pickle.load(file_source)
    with open(path_sets_genes_unimodal, "rb") as file_source:
        sets_genes_unimodal = pickle.load(file_source)
    with open(path_sets_genes_multimodal, "rb") as file_source:
        sets_genes_multimodal = pickle.load(file_source)

    # Compile and return information.
    return {
        "sets_genes_measures": sets_genes_measures,
        "sets_genes_selection": sets_genes_selection,
        "sets_genes_unimodal": sets_genes_unimodal,
        "sets_genes_multimodal": sets_genes_multimodal,
    }


def plot_chart_sets_gene_permutation(
    sets=None,
    count=None,
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        sets (dict<list<str>>): values in sets
        count (int): count of sets
        path_file (str): path to directory and file

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
        count=count,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_sets_gene_permutation(
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

    # Read source information from file.
    source = read_source_sets_gene_permutation(dock=dock)

    # Organize information.
    # Define parameters for each chart.
    charts = list()
    record = {
        "title": "sets_genes_measures",
        "count": 3,
        "sets": source["sets_genes_measures"],
    }
    charts.append(record)
    record = {
        "title": "sets_genes_selection",
        "count": 2,
        "sets": source["sets_genes_selection"],
    }
    charts.append(record)
    record = {
        "title": "sets_genes_unimodal",
        "count": 2,
        "sets": source["sets_genes_unimodal"],
    }
    charts.append(record)
    record = {
        "title": "sets_genes_multimodal",
        "count": 2,
        "sets": source["sets_genes_multimodal"],
    }
    charts.append(record)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_permutation = os.path.join(path_plot, "permutation")
    path_directory = os.path.join(
        path_permutation, "sets"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Iterate on charts.
    for chart in charts:
        # Create chart.
        path_file = os.path.join(
            path_directory, str(chart["title"] + ".svg")
        )
        plot_chart_sets_gene_permutation(
            sets=chart["sets"],
            count=chart["count"],
            path_file=path_file,
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
# Histograms
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )
    path_genes_selection = os.path.join(
        path_selection, "genes_selection.pickle"
    )
    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_unimodal = os.path.join(
        path_candidacy, "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)

    genes_query = integration.read_source_annotation_query_genes_all(
        dock=dock,
    )

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_selection": genes_selection,
        "genes_unimodal": genes_unimodal,
        "genes_multimodal": genes_multimodal,
        "genes_query": genes_query,
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
    name=None,
    values=None,
    path=None
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a gene
        name (str): name of a gene
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
        text=name,
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
    path_candidacy = os.path.join(path_plot, "candidacy")
    path_distribution = os.path.join(path_candidacy, "distribution")
    path_unimodal = os.path.join(path_distribution, "unimodal")
    path_multimodal = os.path.join(path_distribution, "multimodal")
    path_query = os.path.join(path_distribution, "query")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_distribution)
    utility.create_directories(path=path_unimodal)
    utility.create_directories(path=path_multimodal)
    utility.create_directories(path=path_query)

    # Iterate on genes.
    for gene in source_initial["genes_unimodal"]:
        # Access gene's name.
        name = assembly.access_gene_name(
            identifier=gene,
            data_gene_annotation=source_initial["data_gene_annotation"],
        )
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
        path_gene_figure = os.path.join(path_unimodal, str(name + ".svg"))
        plot_chart_genes_persons_signals(
            gene=gene,
            name=name,
            values=values,
            path=path_gene_figure
        )
        pass

    for gene in source_initial["genes_multimodal"]:
        # Access gene's name.
        name = assembly.access_gene_name(
            identifier=gene,
            data_gene_annotation=source_initial["data_gene_annotation"],
        )
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
        path_gene_figure = os.path.join(path_multimodal, str(name + ".svg"))
        plot_chart_genes_persons_signals(
            gene=gene,
            name=name,
            values=values,
            path=path_gene_figure
        )
        pass

    # Plot charts for genes of interest.
    for gene in source_initial["genes_query"]:
        if gene in source_initial["genes_selection"]:
            # Access gene's name.
            name = assembly.access_gene_name(
                identifier=gene,
                data_gene_annotation=source_initial["data_gene_annotation"],
            )
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
            path_gene_figure = os.path.join(path_query, str(name + ".svg"))
            plot_chart_genes_persons_signals(
                gene=gene,
                name=name,
                values=values,
                path=path_gene_figure
            )
            pass

    pass


##########
# Distributions of a gene's pan-tissue aggregate signals across persons after
# permutation
# Histograms
# Status: in progress


def read_source_gene_signals_persons_permutation(
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
    path_probability = os.path.join(dock, "probability")
    path_gene_signals_permutation = os.path.join(
        path_probability, "gene_signals_permutation.pickle"
    )

    # Read information from file.
    with open(path_gene_signals_permutation, "rb") as file_source:
        gene_signals_permutation = pickle.load(file_source)

    # Compile and return information.
    return {
        "gene_signals_permutation": gene_signals_permutation,
    }


def plot_chart_gene_signals_persons_permutation(
    gene=None,
    name=None,
    values=None,
    path=None
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a gene
        name (str): name of a gene
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
        text=name,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )

    pass


def prepare_chart_gene_signals_persons_permutation(
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
    source = read_source_gene_signals_persons_permutation(dock=dock)

    # Define paths.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_permutation = os.path.join(path_plot, "permutation")
    path_directory = os.path.join(path_permutation, "distribution")
    path_file = os.path.join(
        path_directory, "gene_permutation.svg"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    plot_chart_gene_signals_persons_permutation(
        gene="",
        name="",
        values=source["gene_signals_permutation"],
        path=path_file
    )
    pass


##########
# Genes' signals across tissues (rows) and persons (columns), with a chart
# for those persons' properties
# Sort order of persons is by pan-tissue signal
# heatmaps
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )
    path_persons_properties = os.path.join(
        path_selection, "data_persons_properties.pickle"
    )
    path_genes_selection = os.path.join(
        path_selection, "genes_selection.pickle"
    )
    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)

    genes_query = integration.read_source_annotation_query_genes_all(
        dock=dock,
    )

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "genes_selection": genes_selection,
        "genes_multimodal": genes_multimodal,
        "genes_query": genes_query,
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
    property=None,
    type=None,
    data_persons_properties=None,
    data_gene_persons_signals=None,
    data_gene_signals_tissues_persons=None,
):
    """
    Organize information for chart.

    arguments:
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        data_persons_properties (object): Pandas data frame of persons'
            properties
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons
        data_gene_signals_tissues_persons (object): Pandas data frame of a
            gene's signals across tissues and persons

    raises:

    returns:
        (dict): information for chart

    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    data_gene_persons_signals = data_gene_persons_signals.copy(deep=True)
    data_gene_signals_tissues_persons = data_gene_signals_tissues_persons.copy(
        deep=True
    )
    # Rename the aggregate, pantissue signals.
    data_gene_persons_signals.rename(
        columns={
            "value": "signal",
        },
        inplace=True,
    )
    # Introduce aggregate, pantissue signals to tissue-person matrix.
    # These aggregate signals will be useful to sort the data.
    data_sort = data_gene_signals_tissues_persons.join(
        data_gene_persons_signals,
        how="left",
        on="person"
    )
    # Introduce groups to data.
    if type == "category":
        data_persons_properties["property"], properties = pandas.factorize(
            data_persons_properties[property],
            sort=True
        )
        properties = list(properties)
    elif type == "continuity":
        data_persons_properties["property"] = data_persons_properties[property]
        properties = None
    data_persons_properties = data_persons_properties.loc[
        :, data_persons_properties.columns.isin(["property", property])
    ]
    data_hybrid = data_sort.join(
        data_persons_properties,
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
        labels=["signal", property],
        axis="columns",
        inplace=True
    )
    # Compile information.
    bin = dict()
    bin["properties"] = properties
    bin["data"] = data_hybrid
    # Return information
    return bin


def plot_chart_genes_signals_tissues_persons(
    gene=None,
    name=None,
    property=None,
    type=None,
    properties=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a single gene
        name (str): name of gene
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        properties (list<str>): unique values of property
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(name + ".png")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_heatmap_asymmetric_groups(
        title=name,
        property=property,
        type=type,
        properties=properties,
        data=data,
        fonts=fonts,
        colors=colors,
        title_signal="genes' pan-tissue signals across persons",
    )
    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_genes_signals_tissues_persons_gene(
    property=None,
    type=None,
    gene=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    path_directory=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        genes_multimodal (list<str>): identifiers of genes
        genes_query (list<str>): identifiers of genes
        genes_selection (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        path_directory (str): path to directory
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Name.
    name = assembly.access_gene_name(
        identifier=gene,
        data_gene_annotation=data_gene_annotation,
    )
    # Read source information from file.
    source_gene = read_source_genes_signals_tissues_persons(
        gene=gene,
        dock=dock
    )
    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = (
        organize_genes_signals_tissues_persons(
            property=property,
            type=type,
            data_persons_properties=data_persons_properties,
            data_gene_persons_signals=(
                source_gene["data_gene_persons_signals"]
            ),
            data_gene_signals_tissues_persons=(
                source_gene["data_gene_signals_tissues_persons"]
            ),
        )
    )
    # Create charts for the gene.
    plot_chart_genes_signals_tissues_persons(
        gene=gene,
        name=name,
        property=property,
        type=type,
        properties=bin["properties"],
        data=bin["data"],
        path_directory=path_directory
    )
    pass


def prepare_charts_genes_signals_tissues_persons_variable(
    property=None,
    type=None,
    genes_multimodal=None,
    genes_query=None,
    genes_selection=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    path_tissues_persons=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        genes_multimodal (list<str>): identifiers of genes
        genes_query (list<str>): identifiers of genes
        genes_selection (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        path_tissues_persons (str): path to parent directory
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Define paths.
    path_variable = os.path.join(
        path_tissues_persons, property,
    )
    path_multimodal = os.path.join(path_variable, "multimodal")
    path_query = os.path.join(path_variable, "query")
    utility.create_directories(path=path_multimodal)
    utility.create_directories(path=path_query)

    # Iterate on multimodal genes.
    for gene in genes_multimodal:
        prepare_charts_genes_signals_tissues_persons_gene(
            property=property,
            type=type,
            gene=gene,
            data_gene_annotation=data_gene_annotation,
            data_persons_properties=data_persons_properties,
            path_directory=path_multimodal,
            dock=dock,
        )
        pass

    # Iterate on interest genes.
    for gene in genes_query:
        if gene in genes_selection:
            prepare_charts_genes_signals_tissues_persons_gene(
                property=property,
                type=type,
                gene=gene,
                data_gene_annotation=data_gene_annotation,
                data_persons_properties=data_persons_properties,
                path_directory=path_query,
                dock=dock,
            )
            pass
        pass
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
    path_candidacy = os.path.join(path_plot, "candidacy")
    path_tissues_persons = os.path.join(
        path_candidacy, "tissues_persons"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_tissues_persons)

    # Iterate on categorical and ordinal groups of variables.
    variables = list()
    #variables.append(dict(name="sex", type="category"))
    variables.append(dict(name="age", type="continuity"))
    #variables.append(dict(name="body", type="continuity"))
    #variables.append(dict(name="hardiness", type="continuity"))
    #variables.append(dict(name="season_sequence", type="continuity"))
    # Iterate on categorical variables.
    for variable in variables:
        # Prepare charts for genes.
        prepare_charts_genes_signals_tissues_persons_variable(
            property=variable["name"],
            type=variable["type"],
            genes_multimodal=source_initial["genes_multimodal"],
            genes_query=source_initial["genes_query"],
            genes_selection=source_initial["genes_selection"],
            data_gene_annotation=source_initial["data_gene_annotation"],
            data_persons_properties=source_initial["data_persons_properties"],
            path_tissues_persons=path_tissues_persons,
            dock=dock,
        )
        pass
    pass


##########
# Multiple genes' (rows) signals across persons (columns), with a chart
# for those persons' properties
# Sort order of persons is by property and my hierarchical clustering.
# heatmaps
# Status: working


def read_source_genes_signals_persons_properties(
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
    path_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_genes_selection = os.path.join(
        dock, "selection", "tight", "samples_genes_signals",
        "genes.pickle"
    )
    path_persons_sets = os.path.join(
        dock, "selection", "tight", "persons_properties", "charts",
        "persons_sets.pickle"
    )
    path_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", "regression",
        "data_persons_properties.pickle"
    )

    path_distribution = os.path.join(dock, "distribution")
    path_collection = os.path.join(path_distribution, "collection")
    path_data_signals_genes_persons = os.path.join(
        path_collection, "data_signals_genes_persons.pickle"
    )

    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )

    path_prediction_genes_sets = os.path.join(
        dock, "prediction", "sets.pickle"
    )
    path_data_regression_genes = os.path.join(
        dock, "prediction", "data_regression_genes.pickle"
    )


    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)
    with open(path_persons_sets, "rb") as file_source:
        persons_sets = pickle.load(file_source)


    with open(path_prediction_genes_sets, "rb") as file_source:
        prediction_genes_sets = pickle.load(file_source)
    data_regression_genes = pandas.read_pickle(
        path_data_regression_genes
    )

    genes_query_population = (
        integration.read_source_annotation_query_genes_set(
            set="population",
            dock=dock,
    ))

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "data_signals_genes_persons": data_signals_genes_persons,
        "genes_selection": genes_selection,
        "persons_sets": persons_sets,
        "genes_multimodal": genes_multimodal,
        "genes_query_population": genes_query_population,
        "prediction_genes_sets": prediction_genes_sets,
        "data_regression_genes": data_regression_genes,
    }


def organize_genes_signals_persons_properties(
    property=None,
    type=None,
    persons=None,
    sequence=None,
    genes_query=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
):
    """
    Organize information for chart.

    Notice that the data have features (genes) across columns and instances
    (persons) across rows.

    Sequence of genes across rows depends on hierarchical cluster by their
    similarities across persons.
    Sequence of persons across columns depends either on sort by values of
    property or on hierarchical cluster by their similarities across genes.

    arguments:
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        persons (list<str>): identifiers of persons
        sequence (str): method for sequence of persons, sort by property's
            values or cluster by similarities across genes
        genes_query (list<str>): identifiers of genes
        genes_selection (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons

    raises:

    returns:
        (dict): information for chart

    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    data_signals_genes_persons = data_signals_genes_persons.copy(deep=True)

    # Organize signals.
    # Select data for query genes.
    data_signals_genes = data_signals_genes_persons.loc[
        :, data_signals_genes_persons.columns.isin(genes_query)
    ]
    # Select data for persons.
    data_signals_selection = data_signals_genes.loc[
        data_signals_genes.index.isin(persons), :
    ]
    data_persons_properties = data_persons_properties.loc[
        data_persons_properties.index.isin(persons), :
    ]
    # Translate identifiers of genes.
    identifiers = data_signals_selection.columns.to_list()
    translations = dict()
    for identifier in identifiers:
        translations[identifier] = assembly.access_gene_name(
            identifier=identifier,
            data_gene_annotation=data_gene_annotation,
        )
        pass
    data_signals_selection.rename(
        columns=translations,
        inplace=True,
    )
    # Cluster data.
    data_signals_cluster = utility.cluster_data_columns(
        data=data_signals_selection,
    )
    data_signals_sequence = utility.cluster_data_rows(
        data=data_signals_cluster,
    )
    # Organize properties.
    if type == "category":
        data_persons_properties["property"], properties = pandas.factorize(
            data_persons_properties[property],
            sort=True
        )
        properties = list(properties)
    elif type == "continuity":
        data_persons_properties["property"] = data_persons_properties[property]
        properties = None
    data_persons_properties = data_persons_properties.loc[
        :, data_persons_properties.columns.isin(["property", property])
    ]
    data_hybrid = data_signals_sequence.join(
        data_persons_properties,
        how="left",
        on="person"
    )
    # Determine whether to sort by persons' values of property.
    if sequence == "sort":
        data_hybrid.sort_values(
            by=["property"],
            axis="index",
            ascending=True,
            inplace=True,
        )
    # Remove the column for aggregate, pantissue signals.
    data_hybrid.drop(
        labels=[property],
        axis="columns",
        inplace=True
    )
    # Compile information.
    bin = dict()
    bin["properties"] = properties
    bin["data"] = data_hybrid
    # Return information
    return bin


def plot_chart_genes_signals_persons_properties(
    name=None,
    property=None,
    type=None,
    properties=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        name (str): name of file and chart
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        properties (list<str>): unique values of property
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(name + ".png")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_heatmap_asymmetric_groups(
        title=name,
        property=property,
        type=type,
        properties=properties,
        data=data,
        fonts=fonts,
        colors=colors,
        title_signal="genes' pan-tissue signals across persons",
    )
    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_genes_signals_persons_properties_variable(
    property=None,
    type=None,
    persons=None,
    genes_query=None,
    genes_selection=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
    path_sort=None,
    path_cluster=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        property (str): name of feature from persons' properties to use for
            groups
        type (str): type of property, category or continuity
        persons (list<str>): identifiers of persons
        genes_query (list<str>): identifiers of genes
        genes_selection (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        path_sort (str): path to directory
        path_cluster (str): path to directory

    raises:

    returns:

    """

    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_genes_signals_persons_properties(
        property=property,
        type=type,
        persons=persons,
        sequence="sort",
        genes_query=genes_query,
        data_gene_annotation=data_gene_annotation,
        data_persons_properties=data_persons_properties,
        data_signals_genes_persons=data_signals_genes_persons,
    )
    # Create charts for the gene.
    plot_chart_genes_signals_tissues_persons(
        name=property,
        property=property,
        type=type,
        properties=bin["properties"],
        data=bin["data"],
        path_directory=path_sort,
    )

    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_genes_signals_persons_properties(
        property=property,
        type=type,
        persons=persons,
        sequence="cluster",
        genes_query=genes_query,
        data_gene_annotation=data_gene_annotation,
        data_persons_properties=data_persons_properties,
        data_signals_genes_persons=data_signals_genes_persons,
    )
    # Create charts for the gene.
    plot_chart_genes_signals_tissues_persons(
        name=property,
        property=property,
        type=type,
        properties=bin["properties"],
        data=bin["data"],
        path_directory=path_cluster,
    )
    pass


def prepare_charts_genes_signals_persons_properties(
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
    source = read_source_genes_signals_persons_properties(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_integration = os.path.join(path_plot, "integration")
    path_persons_properties = os.path.join(
        path_integration, "persons_properties"
    )
    path_sort = os.path.join(
        path_persons_properties, "sort_persons"
    )
    path_cluster = os.path.join(
        path_persons_properties, "cluster_persons"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_persons_properties)
    utility.create_directories(path=path_sort)
    utility.create_directories(path=path_cluster)

    variables_regression = selection.define_variables()
    sets = prediction.select_genes_by_group_association(
        variables=variables_regression["hypothesis"],
        genes_selection=source["genes_selection"],
        genes_unimodal=list(),
        genes_multimodal=list(),
        threshold_r_square=0.15,
        count_selection=25,
        data_regression_genes=source["data_regression_genes"],
    )

    # Define genes of interest.
    #genes_interest = source["genes_query_population"]
    #genes_interest = sets["selection"]["sex_y_scale"]
    #genes_interest = sets["selection"]["age_scale"]
    #genes_interest = sets["selection"]["respiration_binary_scale"]
    genes_interest = sets["selection"]["ventilation_duration_scale"]
    utility.print_terminal_partition(level=2)
    print(len(genes_interest))

    # Iterate on categorical and ordinal groups of variables.
    variables = list()
    #variables.append(dict(name="sex_text", type="category"))
    #variables.append(dict(name="age", type="continuity"))
    #variables.append(dict(name="respiration", type="category"))
    variables.append(dict(name="ventilation_duration", type="continuity"))
    # Iterate on categorical variables.

    # TODO: I need to specify the persons of interest here...

    for variable in variables:
        # Prepare charts for genes.
        prepare_charts_genes_signals_persons_properties_variable(
            property=variable["name"],
            type=variable["type"],
            persons=source["persons_sets"]["ventilation"],
            genes_query=genes_interest,
            genes_selection=source["genes_selection"],
            data_gene_annotation=source["data_gene_annotation"],
            data_persons_properties=source["data_persons_properties"],
            data_signals_genes_persons=source["data_signals_genes_persons"],
            path_sort=path_sort,
            path_cluster=path_cluster,
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_persons_properties = os.path.join(
        path_selection, "data_persons_properties.pickle"
    )
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
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
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    with open(path_genes_integration, "rb") as file_source:
        genes_integration = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "data_signals_genes_persons": data_signals_genes_persons,
        "genes_integration": genes_integration,
    }


def organize_signals_genes_persons_groups(
    genes=None,
    sorts=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        genes (list<str>): identifiers of genes
        sorts (list<str>): names of columns by which to sort data
        data_gene_annotation (object): Pandas data frame of genes' annotations
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
        by=sorts,
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

    # Prepare translation of genes' identifiers to names.
    translations = dict()
    for identifier in genes:
        translations[identifier] = assembly.access_gene_name(
            identifier=identifier,
            data_gene_annotation=data_gene_annotation,
        )
        pass

    # Rename genes in data.
    data_hybrid.rename(
        columns=translations,
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
    data = organize_signals_genes_persons_groups(
        genes=source["genes_integration"],
        sorts=["sex", "age"],
        data_gene_annotation=source["data_gene_annotation"],
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
# Principal components groups of persons
# Scatter plot
# Status: working


def read_source_persons_genes_components(
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
    path_data_persons_genes_components = os.path.join(
        path_integration, "data_persons_genes_components.pickle"
    )
    path_data_persons_genes_variances = os.path.join(
        path_integration, "data_persons_genes_variances.pickle"
    )

    # Read information from file.
    data_persons_genes_components = pandas.read_pickle(
        path_data_persons_genes_components
    )
    data_persons_genes_variances = pandas.read_pickle(
        path_data_persons_genes_variances
    )
    # Compile and return information.
    return {
        "data_persons_genes_components": data_persons_genes_components,
        "data_persons_genes_variances": data_persons_genes_variances,
    }

# TODO: enable this function to determine variances automatically

def plot_chart_persons_genes_components(
    data_components=None,
    data_variances=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        data_components (object): Pandas data frame of principal components
            across instances
        data_variances (object): Pandas data frame of variances for principal
            components
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Organize data.
    data_components = data_components.copy(deep=True)
    data_components.reset_index(
        level=None,
        inplace=True
    )
    data_variances.set_index(
        "component",
        drop=True,
        inplace=True,
    )
    #print(data_variances)
    title_1 = (
        "Component 1 (" +
        str(round((data_variances.at[1, "variance"] * 100), 1)) +
        "%)"
    )
    title_2 = (
        "Component 2 (" +
        str(round((data_variances.at[2, "variance"] * 100), 1)) +
        "%)"
    )
    title_3 = (
        "Component 3 (" +
        str(round((data_variances.at[3, "variance"] * 100), 1)) +
        "%)"
    )

    # Define file name.
    title = str("components_1_2")
    path_file = os.path.join(
        path_directory, str(title + ".svg")
    )
    # Create figure.
    figure = plot_scatter(
        data=data_components,
        abscissa="component_1",
        ordinate="component_2",
        title_abscissa=title_1,
        title_ordinate=title_2,
        fonts=fonts,
        colors=colors,
        size=5,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    # Define file name.
    title = str("components_1_3")
    path_file = os.path.join(
        path_directory, str(title + ".svg")
    )
    # Create figure.
    figure = plot_scatter(
        data=data_components,
        abscissa="component_1",
        ordinate="component_3",
        title_abscissa=title_1,
        title_ordinate=title_3,
        fonts=fonts,
        colors=colors,
        size=5,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    # Define file name.
    title = str("components_2_3")
    path_file = os.path.join(
        path_directory, str(title + ".svg")
    )
    # Create figure.
    figure = plot_scatter(
        data=data_components,
        abscissa="component_2",
        ordinate="component_3",
        title_abscissa=title_2,
        title_ordinate=title_3,
        fonts=fonts,
        colors=colors,
        size=5,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_persons_genes_components(
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
    source = read_source_persons_genes_components(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_integration = os.path.join(path_plot, "integration")
    path_directory = os.path.join(
        path_integration, "persons_components"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Plot.
    plot_chart_persons_genes_components(
        data_components=source["data_persons_genes_components"],
        data_variances=source["data_persons_genes_variances"],
        path_directory=path_directory,
    )

    pass



##########
# Correlations in signals between pairs of genes
# Heatmap
# Status: working


def read_source_signals_genes_correlations(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )
    path_integration = os.path.join(dock, "integration")
    path_data_correlation_genes_unimodal = os.path.join(
        path_integration, "data_correlation_genes_unimodal.pickle"
    )
    path_data_correlation_genes_multimodal = os.path.join(
        path_integration, "data_correlation_genes_multimodal.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_correlation_genes_unimodal = pandas.read_pickle(
        path_data_correlation_genes_unimodal
    )
    data_correlation_genes_multimodal = pandas.read_pickle(
        path_data_correlation_genes_multimodal
    )

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_correlation_genes_unimodal": data_correlation_genes_unimodal,
        "data_correlation_genes_multimodal": data_correlation_genes_multimodal,
    }


def organize_signals_genes_correlations(
    data_gene_annotation=None,
    data_genes_correlations=None,
):
    """
    Organizes data for plot.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_genes_correlations (object): Pandas data frame of correlations
            between genes

    raises:

    returns:
        (object): Pandas data frame of correlations between genes

    """

    # Prepare translation of genes' identifiers to names.
    identifiers = data_genes_correlations.index.to_list()
    translations = dict()
    for identifier in identifiers:
        translations[identifier] = assembly.access_gene_name(
            identifier=identifier,
            data_gene_annotation=data_gene_annotation,
        )
        #print(identifier + " ... " + translations[identifier])
        #print(identifier)
        pass

    # Rename genes in data.
    data_genes_correlations.rename(
        index=translations,
        inplace=True,
    )
    data_genes_correlations.rename(
        columns=translations,
        inplace=True,
    )

    return data_genes_correlations


def plot_chart_signals_genes_correlations(
    data=None,
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_file (str): path for file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_heatmap_symmetric_diverge(
        data=data,
        minimum=-1.0,
        maximum=1.0,
        label_columns=True,
        label_rows=True,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_signals_genes_correlations(
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
    source = read_source_signals_genes_correlations(dock=dock)

    # Organize data.
    data_unimodal = organize_signals_genes_correlations(
        data_gene_annotation=source["data_gene_annotation"],
        data_genes_correlations=source["data_correlation_genes_unimodal"],
    )
    utility.print_terminal_partition(level=2)
    data_multimodal = organize_signals_genes_correlations(
        data_gene_annotation=source["data_gene_annotation"],
        data_genes_correlations=source["data_correlation_genes_multimodal"],
    )

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_integration = os.path.join(path_plot, "integration")
    path_directory = os.path.join(
        path_integration, "correlations_genes"
    )
    # Define file name.
    path_unimodal = os.path.join(
        path_directory, str("genes_unimodal.svg")
    )
    path_multimodal = os.path.join(
        path_directory, str("genes_multimodal.svg")
    )

    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Create chart.
    plot_chart_signals_genes_correlations(
        data=data_unimodal,
        path_file=path_unimodal
    )
    plot_chart_signals_genes_correlations(
        data=data_multimodal,
        path_file=path_multimodal
    )

    pass


##########
# Correlations in signals between pairs of genes
# Heatmap
# Status: working


def read_source_signals_genes_correlations_prediction(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )

    path_integration = os.path.join(dock, "integration")
    path_data_correlation_multimodal_hardiness = os.path.join(
        path_integration, "data_correlation_multimodal_hardiness.pickle"
    )
    path_data_correlation_multimodal_sex_age_body = os.path.join(
        path_integration, "data_correlation_multimodal_sex_age_body.pickle"
    )
    path_data_correlation_multimodal_union = os.path.join(
        path_integration, "data_correlation_multimodal_union.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_correlation_multimodal_hardiness = pandas.read_pickle(
        path_data_correlation_multimodal_hardiness
    )
    data_correlation_multimodal_sex_age_body = pandas.read_pickle(
        path_data_correlation_multimodal_sex_age_body
    )
    data_correlation_multimodal_union = pandas.read_pickle(
        path_data_correlation_multimodal_union
    )

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_correlation_multimodal_hardiness": (
            data_correlation_multimodal_hardiness
        ),
        "data_correlation_multimodal_sex_age_body": (
            data_correlation_multimodal_sex_age_body
        ),
        "data_correlation_multimodal_union": data_correlation_multimodal_union,
    }


def prepare_charts_signals_genes_correlations_prediction(
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
    source = read_source_signals_genes_correlations_prediction(dock=dock)

    # Organize data.
    data_hardiness = organize_signals_genes_correlations(
        data_gene_annotation=source["data_gene_annotation"],
        data_genes_correlations=(
            source["data_correlation_multimodal_hardiness"]
        ),
    )
    data_sex_age_body = organize_signals_genes_correlations(
        data_gene_annotation=source["data_gene_annotation"],
        data_genes_correlations=(
            source["data_correlation_multimodal_sex_age_body"]
        ),
    )
    data_union = organize_signals_genes_correlations(
        data_gene_annotation=source["data_gene_annotation"],
        data_genes_correlations=source["data_correlation_multimodal_union"],
    )

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_prediction = os.path.join(path_plot, "prediction")
    path_directory = os.path.join(
        path_prediction, "correlations_genes"
    )
    # Define file name.
    path_hardiness = os.path.join(
        path_directory, str("genes_multimodal_hardiness.svg")
    )
    path_sex_age_body = os.path.join(
        path_directory, str("genes_multimodal_sex_age_body.svg")
    )
    path_union = os.path.join(
        path_directory, str("genes_multimodal_union.svg")
    )

    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Create chart.
    plot_chart_signals_genes_correlations(
        data=data_hardiness,
        path_file=path_hardiness
    )
    plot_chart_signals_genes_correlations(
        data=data_sex_age_body,
        path_file=path_sex_age_body
    )
    plot_chart_signals_genes_correlations(
        data=data_union,
        path_file=path_union
    )

    pass


##########
# Correlations in signals between pairs of genes from Gene Ontology enrichment
# Heatmap
# Status: in progress


def read_source_signals_genes_correlations_query(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
    }


def prepare_charts_signals_genes_correlations_query(
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
    source = read_source_signals_genes_correlations_query(dock=dock)

    # Specify directory for source files.
    path_integration = os.path.join(dock, "integration")
    path_query_correlations = os.path.join(
        path_integration, "query_correlations"
    )
    files_query_correlations = os.listdir(path=path_query_correlations)

    # Specify directory for product files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_plot_integration = os.path.join(path_plot, "integration")
    path_directory = os.path.join(
        path_plot_integration, "correlations_query_genes"
    )

    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Iterate on source files.
    for file in files_query_correlations:
        # Extract name of set.
        file_first = file.replace("data_", "")
        file_last = file_first.replace(".pickle", "")
        name = file_last
        # Read information.
        path_data = os.path.join(
            path_query_correlations, file
        )
        data = pandas.read_pickle(path_data)
        # Organize data.
        data_organization = organize_signals_genes_correlations(
            data_gene_annotation=source["data_gene_annotation"],
            data_genes_correlations=data,
        )
        # Define file name.
        path_file = os.path.join(
            path_directory, str(name + ".svg")
        )
        # Create chart.
        plot_chart_signals_genes_correlations(
            data=data_organization,
            path_file=path_file
        )

    pass


##########
# Correlations in pantissue signals across persons between pairs of specific
# genes of interest
# scatter plots
# Status: working


def read_source_signals_persons_gene_pairs(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )
    path_distribution = os.path.join(dock, "distribution")
    path_collection = os.path.join(path_distribution, "collection")
    path_data_signals_genes_persons = os.path.join(
        path_collection, "data_signals_genes_persons.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_signals_genes_persons": data_signals_genes_persons,
    }


def organize_signals_persons_gene_pairs(
    pairs_genes=None,
    data_gene_annotation=None,
    data_signals_genes_persons=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        pairs_genes (list<tuple>): pairs of genes' identifiers
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_signals_genes_persons (object): Pandas data frame of genes'
            pantissue signals across persons

    raises:

    returns:
        (list<dict>): information about each pair of genes

    """

    # Collect information about pairs of genes.
    records = list()
    # Iterate on pairs.
    for genes in pairs_genes:
        # Determine genes' names.
        names = dict()
        names[genes[0]] = assembly.access_gene_name(
            identifier=genes[0],
            data_gene_annotation=data_gene_annotation,
        )
        names[genes[1]] = assembly.access_gene_name(
            identifier=genes[1],
            data_gene_annotation=data_gene_annotation,
        )
        # Organize genes' signals.
        # Remove persons with null signals for both genes.
        data_pair = data_signals_genes_persons.loc[:, list(genes)]
        data_pair.dropna(
            axis="index",
            how="any",
            inplace=True,
        )
        signals = dict()
        signals[genes[0]] = data_pair[genes[0]].values
        signals[genes[1]] = data_pair[genes[1]].values
        # Collect information.
        record = dict()
        record["genes"] = genes
        record["names"] = names
        record["signals"] = signals
        record["data"] = data_pair
        records.append(record)

        pass

    # Return information.
    return records


def plot_chart_signals_persons_gene_pair(
    genes=None,
    names=None,
    signals=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        genes (tuple<str>): identifiers of genes
        names (dict<str>): names of genes
        signals (dict<array>): values of genes' pantissue signals across
            persons
        data (object): Pandas data frame of variables
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    title = str(names[genes[0]] + "_" + names[genes[1]])
    path_file = os.path.join(
        path_directory, str(title + ".svg")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_scatter(
        data=data,
        abscissa=genes[0],
        ordinate=genes[1],
        title_abscissa=names[genes[0]],
        title_ordinate=names[genes[1]],
        fonts=fonts,
        colors=colors,
        size=5,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_signals_persons_gene_pairs(
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

    # Specify pairs of genes of interest.
    # PWP2: ENSG00000241945
    # GATD3A: ENSG00000160221
    # GATD3B: ENSG00000280071
    # LOC102724159: ENSG00000275464
    # SLC7A11: ENSG00000151012
    # ADM2: ENSG00000128165
    # TBC1D3D: ENSG00000274419
    # TBC1D3L: ENSG00000274512
    # NPIPB2: ENSG00000234719
    # NPIPB15: ENSG00000196436
    pairs_genes = [
        # GATD3A, GATD3B... strong correlation
        ("ENSG00000160221", "ENSG00000280071"),
        # IL1R2, LILRA2... strong correlation
        ("ENSG00000115590", "ENSG00000239998"),
        # XKR4, PTHLH... strong correlation
        ("ENSG00000206579", "ENSG00000087494"),
        # XKR4, ALOX15B... strong correlation
        ("ENSG00000206579", "ENSG00000179593"),
        # CCL20, IL1R2... strong correlation
        ("ENSG00000115009", "ENSG00000115590"),
        # COL27A1, CEP57L1... strong correlation
        ("ENSG00000196739", "ENSG00000183137"),
        # PHC2, SORT1... strong correlation
        ("ENSG00000134686", "ENSG00000134243"),

        # CD163, LILRA2... strong correlation
        ("ENSG00000177575", "ENSG00000239998"),
        # PHC2, CKAP4... strong correlation
        ("ENSG00000134686", "ENSG00000136026"),
        # CKAP4, CXCL5... strong correlation
        ("ENSG00000136026", "ENSG00000163735"),
        # GCKR, RETN... strong correlation, both relate to diabetes risk
        ("ENSG00000084734", "ENSG00000104918"),
        # CD163, HPR... strong correlation, both genes are relevant to blood hemoglobin and haptoglobin
        ("ENSG00000177575", "ENSG00000261701"),
        # CERKL, CNGA4... strong correlation, both relate to neurons
        ("ENSG00000188452", "ENSG00000132259"),
        # PHC2, CLP1... strong correlation, both relate to transcriptional regulation
        ("ENSG00000134686", "ENSG00000172409"),

        # SPRED3, GCKR... strong correlation
        ("ENSG00000188766", "ENSG00000084734"),
        # SPRED3, RETN... strong correlation
        ("ENSG00000188766", "ENSG00000104918"),

        # NPIPA5, NPIPB2... weak correlation
        #("ENSG00000183793", "ENSG00000234719"),
        # NPIPA5, NPIPB15... weak correlation
        #("ENSG00000183793", "ENSG00000196436"),
        # NPIPB2, NPIPB15... weak correlation
        #("ENSG00000234719", "ENSG00000196436"),

        # LIPG, LILRA2... weak correlation
        #("ENSG00000101670", "ENSG00000239998"),
        # GTF2A1L, ZBTB16... weak correlation
        #("ENSG00000242441", "ENSG00000109906"),
        # GTF2A1L, ALOX15B... weak correlation
        #("ENSG00000242441", "ENSG00000179593"),
        # CLP1, ZBTB16... weak correlation
        #("ENSG00000172409", "ENSG00000109906"),
        # GTF2A1L, PHC2... weak correlation
        #("ENSG00000242441", "ENSG00000134686"),
        # CD163, LIPG... weak correlation
        #("ENSG00000177575", "ENSG00000101670"),
        # SPRED3, PTHLH... weak correlation
        #("ENSG00000188766", "ENSG00000087494"),
        # XKR4, XKR9... weak correlation
        #("ENSG00000206579", "ENSG00000221947"),
        # FAM20C, FKBP5... weak correlation
        #("ENSG00000177706", "ENSG00000096060"),

    ]

    # Read source information from file.
    source = read_source_signals_persons_gene_pairs(dock=dock)

    # Organize data.
    pairs = organize_signals_persons_gene_pairs(
        pairs_genes=pairs_genes,
        data_gene_annotation=source["data_gene_annotation"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
    )

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_candidacy = os.path.join(path_plot, "candidacy")
    path_scatter = os.path.join(
        path_candidacy, "scatter"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_scatter)
    utility.create_directories(path=path_scatter)

    # Iterate on pairs of genes.
    for pair in pairs:
        # Create chart.
        plot_chart_signals_persons_gene_pair(
            genes=pair["genes"],
            names=pair["names"],
            signals=pair["signals"],
            data=pair["data"],
            path_directory=path_scatter
        )

        pass

    pass


##########
# Regressions on genes' pan-tissue signals across persons
# Status: in progress


def read_source_regressions_genes(
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
    path_prediction = os.path.join(dock, "prediction")
    path_data_regression_genes_selection_scale = os.path.join(
        path_prediction, "data_regression_genes_selection_scale.pickle"
    )
    path_data_regression_genes_unimodal_scale = os.path.join(
        path_prediction, "data_regression_genes_unimodal_scale.pickle"
    )
    path_data_regression_genes_multimodal_scale = os.path.join(
        path_prediction, "data_regression_genes_multimodal_scale.pickle"
    )

    # Read information from file.
    data_regression_genes_selection_scale = pandas.read_pickle(
        path_data_regression_genes_selection_scale
    )
    data_regression_genes_unimodal_scale = pandas.read_pickle(
        path_data_regression_genes_unimodal_scale
    )
    data_regression_genes_multimodal_scale = pandas.read_pickle(
        path_data_regression_genes_multimodal_scale
    )

    # Compile and return information.
    return {
        "data_regression_genes_selection_scale": (
            data_regression_genes_selection_scale
        ),
        "data_regression_genes_unimodal_scale": (
            data_regression_genes_unimodal_scale
        ),
        "data_regression_genes_multimodal_scale": (
            data_regression_genes_multimodal_scale
        ),
    }


def plot_charts_regressions_genes(
    path=None,
    data=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        path (str): path to directory and file
        data (object): Pandas data frame

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Charts for each set of genes...
    # 1. distribution of R-square adjust
    # 2. swarm plots for each biological covariate
    # 3. violin plots?

    # Specify directories and files.
    path_figure_swarm = os.path.join(
        path, "swarm.svg"
    )

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


def prepare_charts_regressions_genes(
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

    print("going to plot charts for regressions on genes")

    # Read source information from file.
    source = read_source_regressions_genes(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_prediction = os.path.join(path_plot, "prediction")
    path_selection = os.path.join(path_prediction, "selection")
    path_unimodal = os.path.join(path_prediction, "unimodal")
    path_multimodal = os.path.join(path_prediction, "multimodal")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_prediction)
    utility.create_directories(path=path_selection)
    utility.create_directories(path=path_unimodal)
    utility.create_directories(path=path_multimodal)

    # Prepare charts.
    plot_charts_regressions_genes(
        path_directory=path_selection,
        data=source["data_regression_genes_selection_scale"],
    )
    plot_charts_regressions_genes(
        path_directory=path_unimodal,
        data=source["data_regression_genes_unimodal_scale"],
    )
    plot_charts_regressions_genes(
        path_directory=path_multimodal,
        data=source["data_regression_genes_multimodal_scale"],
    )
    pass


##########
# Residuals from regression on genes' pan-tissue signals
# scatter plots
# Status: in progress


def read_source_genes_regression_residuals(
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )
    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )

    path_prediction = os.path.join(dock, "prediction")
    path_residuals_genes = os.path.join(
        path_prediction, "residuals_genes.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)
    with open(path_residuals_genes, "rb") as file_source:
        residuals_genes = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_multimodal": genes_multimodal,
        "residuals_genes": residuals_genes,
    }


def organize_gene_regression_residuals(
    gene=None,
    residuals_genes=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a gene
        residuals_genes (dict<array>): collection of residuals from regressions
            for each gene

    raises:

    returns:
        (object): Pandas data frame of a gene's residuals across observations

    """

    # Access gene's residuals.
    residuals = copy.deepcopy(residuals_genes[gene])

    # Organize data.
    collection = dict()
    collection["residual"] = residuals
    collection["observation"] = numpy.arange(
        1,
        (len(residuals) + 1),
        1,
    )
    data = pandas.DataFrame(
        data=collection,
    )
    # Return information
    return data


def plot_chart_gene_regression_residuals(
    gene=None,
    name=None,
    data=None,
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a single gene
        name (str): name of gene
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_file (str): path for file

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_scatter(
        data=data,
        abscissa="observation",
        ordinate="residual",
        title_abscissa="Person",
        title_ordinate="Residuals",
        fonts=fonts,
        colors=colors,
        size=7,
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_genes_regression_residuals(
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
    source = read_source_genes_regression_residuals(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_prediction = os.path.join(path_plot, "prediction")
    path_directory = os.path.join(
        path_prediction, "residuals"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Iterate on genes.
    for gene in source["genes_multimodal"]:
        # Name.
        name = assembly.access_gene_name(
            identifier=gene,
            data_gene_annotation=source["data_gene_annotation"],
        )
        # Define file name.
        path_file = os.path.join(
            path_directory, str(name + ".svg")
        )
        # Organize data.
        data_gene = organize_gene_regression_residuals(
            gene=gene,
            residuals_genes=source["residuals_genes"],
        )
        # Create charts for the gene.
        plot_chart_gene_regression_residuals(
            gene=gene,
            name=name,
            data=data_gene,
            path_file=path_file
        )
        pass
    pass





###################################################################
################# Need to Update ###############################



def report_metrics_from_modality(
    name=None,
    series=None,
    dock=None,
):
    """
    Calculates and prints reports on multiple metrics.

    arguments:
        name (str): name of series
        series (list<float>): series of values of type float
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:

    """

    # Calculate skewness.
    skewness = scipy.stats.skew(series, axis=None)
    # Calculate kurtosis.
    kurtosis = scipy.stats.kurtosis(series, axis=None, fisher=True)
    # Calculate bimodality coefficient.
    coefficient = calculate_bimodality_coefficient(series)
    # Calculate dip statistic.
    dip = calculate_dip_statistic(series=series)
    # Calculate mixture model score.
    mixture = calculate_mixture_model_score(series=series)

    # Prepare metric report text.
    text = prepare_metric_report_text(
        skewness=skewness,
        kurtosis=kurtosis,
        coefficient=coefficient,
        dip=dip,
        mixture=mixture,
    )
    print(text)

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()

    # Create figure.
    figure = plot.plot_distribution_histogram(
        series=series,
        name="",
        bin_method="count",
        bin_count=50,
        label_bins="Bins",
        label_counts="Counts",
        fonts=fonts,
        colors=colors,
        line=False,
        position=0,
        text=text,
    )
    # Specify directories and files.
    path_metric = os.path.join(dock, "metric")
    path_figure = os.path.join(path_metric, "figure")
    utility.create_directory(path_figure)
    file = (name + "_distribution.svg")
    path_file = os.path.join(path_figure, file)
    # Write figure.
    plot.write_figure(
        path=path_file,
        figure=figure
    )

    pass



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

    ##########
    ##########
    ##########
    # Assembly procedure

    # Plot charts, heatmaps, for values of GTEx clinical health variables
    # across persons.
    # Clinical health variables will be across rows, and persons will be across
    # columns.
    # Sort order across rows depends on hierarchical clustering.
    # Sort order across columns depends on hierarchical clustering.
    #prepare_charts_persons_health_variables(dock=dock)

    # Plot charts, heatmaps, for adjacent summaries of properties across
    # persons.
    #prepare_charts_persons_properties_adjacency(dock=dock)

    ##########
    ##########
    ##########
    # Selection procedure

    # Plot charts for sex, age, and hardiness of persons in GTEx cohort.
    #prepare_charts_persons_sex_age_hardiness(dock=dock)

    # Plot charts of overlap between sets of persons by clinical categories.
    #prepare_charts_sets_selection_persons_properties(dock=dock)



    if False:

        # Plot chart for counts of persons per major tissue type.
        prepare_chart_persons_per_tissue(dock=dock)

        # Plot chart for counts of major tissue type per person, a histogram.
        prepare_chart_tissues_per_person(dock=dock)


        # Plot charts of overlap between sets in selection of genes and samples.
        prepare_charts_sets_selection_genes_samples(dock=dock)

        # Plot charts for selection of principal components on categorical
        # variables for regression.
        prepare_charts_regression_principal_components(dock=dock)

        ##########
        ##########
        ##########
        # Distribution procedure

        # Plot charts of distributions of genes' pan-tissue aggregate signals
        # across persons.
        prepare_charts_genes_persons_signals(dock=dock)

        # Plot charts of distributions of each modality measure's scores across
        # genes.
        prepare_charts_modality_gene_distribution(dock=dock)


        ##########
        ##########
        ##########
        # Candidacy procedure

    if False:


        # Plot charts of overlap between sets in selection of genes by bimodality.
        prepare_charts_gene_sets_candidacy(dock=dock)

        # Plot charts for correlations in signals between pairs of genes.
        # Charts will be scatter plots.
        # Each point will represent a person with pantissue signals from each gene.
        prepare_charts_signals_persons_gene_pairs(dock=dock)

        # Plot charts, heatmaps, for each gene's signals across persons (columns)
        # and tissues (rows).
        # Include a chart to depict persons' properties.
        # Sort order of persons is by pan-tissue signal.
        prepare_charts_genes_signals_tissues_persons(dock=dock)

    if False:

        ##########
        ##########
        ##########
        # Heritability procedure

        # Plot charts for heritability of genes' pantissue signals.
        # Charts will be scatter plots.
        prepare_charts_gene_heritability(dock=dock)

        # Plot charts of overlap between sets in selection of genes by
        # heritability.
        prepare_charts_sets_gene_heritability(dock=dock)

        ##########
        ##########
        ##########
        # Probability procedure

        # Plot charts of overlap between sets in significance of genes by
        # permutation.
        prepare_charts_sets_gene_permutation(dock=dock)

        # Plot chart of the distribution of a gene's pan-tissue aggregate signals
        # across persons after permutation.
        prepare_chart_gene_signals_persons_permutation(dock=dock)

        ##########
        ##########
        ##########
        # Prediction procedure

        # Plot charts, scatter plots, for residuals from regressions on each gene's
        # pan-tissue signals across persons.
        prepare_charts_genes_regression_residuals(dock=dock)


        ##########
        ##########
        ##########
        # Integration procedure

    # Plot charts, scatter plots, for components by genes' pan-tissue
    # signals across groups of persons.
    #prepare_charts_persons_genes_components(dock=dock)

    # Plot charts for correlations between pairs of genes of interest from
    # Gene Ontology enrichment.
    # Chart is adjacency matrix heatmap.
    #prepare_charts_signals_genes_correlations_query(dock=dock)

    # Plot charts, heatmaps, for multiple genes' pan-tissue signals across
    # persons along with those persons' properties.
    # Genes will be across rows, and persons will be across columns.
    # Sort order across rows depends on hierarchical clustering.
    # In some charts, sort order across columns depends on persons' properties
    # (sex, age, body mass index, hardiness).
    # In other charts, sort order across columns depends on hierarchical
    # clustering.
    prepare_charts_genes_signals_persons_properties(dock=dock)

    if False:

        # Plot charts for correlations between pairs of all genes of interest.
        # Chart is adjacency matrix heatmap.
        prepare_charts_signals_genes_correlations(dock=dock)

        # Plot charts for correlations between pairs of all genes of interest.
        # Specific to multimodal genes that associate with hypothetical variables
        # in prediction procedure.
        # Chart is adjacency matrix heatmap.
        prepare_charts_signals_genes_correlations_prediction(dock=dock)

        # TODO: I think this chart might be obsolete?
        prepare_charts_signals_genes_persons_groups(dock=dock)


    if False:

        # Plot chart for batches per person.
        prepare_chart_batches_per_person(dock=dock)

        # Plot charts for regressions on genes.
        prepare_charts_regressions_genes(dock=dock)

        # Plot charts of distributions of genes' bimodality measurement scores and
        # permutations.
        plot_charts_distribution(dock=dock)

        # Plot charts of overlap between sets in selection of genes by permutation
        # probability of bimodality.
        prepare_charts_gene_sets_probability(dock=dock)

        # Plot charts of overlap between sets in selection of genes by
        # integration.
        prepare_charts_gene_sets_integration(dock=dock)

    #plot_charts_analysis(dock=dock)
    #plot_charts_tissue(dock=dock)
    #plot_charts_restriction(dock=dock)
    #plot_chart_tissues_persons(dock=dock)

    pass


if (__name__ == "__main__"):
    execute_procedure()
