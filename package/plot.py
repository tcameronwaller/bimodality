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
        "stretch": 900,
        "weight": 1000,
        "size": 25
    }
    values_three = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 750,
        "weight": 1000,
        "size": 20
    }
    values_four = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 500,
        "size": 17
    }
    values_five = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 400,
        "weight": 400,
        "size": 15
    }
    values_six = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 300,
        "weight": 300,
        "size": 13
    }
    values_seven = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 300,
        "weight": 300,
        "size": 10
    }
    values_eight = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 200,
        "weight": 200,
        "size": 7
    }
    values_nine = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 150,
        "weight": 150,
        "size": 5
    }
    values_ten = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 100,
        "weight": 100,
        "size": 3
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
    properties_six = matplotlib.font_manager.FontProperties(
        family=values_six["family"],
        style=values_six["style"],
        variant=values_six["variant"],
        stretch=values_six["stretch"],
        weight=values_six["weight"],
        size=values_six["size"]
    )
    properties_seven = matplotlib.font_manager.FontProperties(
        family=values_seven["family"],
        style=values_seven["style"],
        variant=values_seven["variant"],
        stretch=values_seven["stretch"],
        weight=values_seven["weight"],
        size=values_seven["size"]
    )
    properties_eight = matplotlib.font_manager.FontProperties(
        family=values_eight["family"],
        style=values_eight["style"],
        variant=values_eight["variant"],
        stretch=values_eight["stretch"],
        weight=values_eight["weight"],
        size=values_eight["size"]
    )
    properties_nine = matplotlib.font_manager.FontProperties(
        family=values_nine["family"],
        style=values_nine["style"],
        variant=values_nine["variant"],
        stretch=values_nine["stretch"],
        weight=values_nine["weight"],
        size=values_nine["size"]
    )
    properties_ten = matplotlib.font_manager.FontProperties(
        family=values_ten["family"],
        style=values_ten["style"],
        variant=values_ten["variant"],
        stretch=values_ten["stretch"],
        weight=values_ten["weight"],
        size=values_ten["size"]
    )
    # Compile and return references.
    return {
        "values": {
            "one": values_one,
            "two": values_two,
            "three": values_three,
            "four": values_four,
            "five": values_five,
            "six": values_six,
            "seven": values_seven,
            "eight": values_eight,
            "nine": values_nine,
            "ten": values_ten,
        },
        "properties": {
            "one": properties_one,
            "two": properties_two,
            "three": properties_three,
            "four": properties_four,
            "five": properties_five,
            "six": properties_six,
            "seven": properties_seven,
            "eight": properties_eight,
            "nine": properties_nine,
            "ten": properties_ten,
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
    # Clear.
    clear = (1.0, 1.0, 1.0, 0.0)
    clear_faint = (1.0, 1.0, 1.0, 0.25)
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
        "clear": clear,
        "clear_faint": clear_faint,
        "blue": blue,
        "blue_faint": blue_faint,
        "orange": orange,
        "orange_faint": orange_faint
    }


# TODO: pass variable for label on scale bar...

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
    label_bar = "Spearman correlation (FDR <= 0.05)"
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
        pad=2,
        labelsize=fonts["values"]["five"]["size"],
        labelcolor=colors["black"],
    )

    # Create ticks and labels for each grid.
    # Let the horizontal axes labeling appear on top.
    axes.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=1.0,
        color=colors["black"],
        pad=5,
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
    if (label_columns and data.shape[1] <= 100):
        if (100 >= data.shape[1] and data.shape[1] >= 90):
            size_column = "ten"
        elif (90 > data.shape[1] and data.shape[1] >= 75):
            size_column = "nine"
        elif (75 > data.shape[1] and data.shape[1] >= 50):
            size_column = "eight"
        elif (50 > data.shape[1] and data.shape[1] >= 25):
            size_column = "seven"
        elif (25 > data.shape[1]):
            size_column = "six"
        axes.set_xticks(numpy.arange(matrix.shape[1]))
        axes.set_xticklabels(
            labels_columns,
            #minor=False,
            rotation=-60,
            rotation_mode="anchor",
            ha="right", # horizontal alignment
            va="bottom", # vertical alignment
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_column]
        )
    if (label_rows and data.shape[0] <= 100):
        if (100 >= data.shape[0] and data.shape[0] >= 90):
            size_row = "nine"
        elif (90 > data.shape[0] and data.shape[0] >= 75):
            size_row = "eight"
        elif (75 > data.shape[0] and data.shape[0] >= 50):
            size_row = "six"
        elif (50 > data.shape[0] and data.shape[0] >= 25):
            size_row = "five"
        elif (25 > data.shape[0]):
            size_row = "three"
        axes.set_yticks(numpy.arange(matrix.shape[0]))
        axes.set_yticklabels(
            labels_rows,
            #minor=False,
            ha="right", # horizontal alignment
            va="center", # vertical alignment
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_row]
        )
    if False:
        # Rotate the tick labels and set their alignment.
        matplotlib.pyplot.setp(
            axes.get_xticklabels(),
            rotation=-30,
            ha="right",
            rotation_mode="anchor"
        )
        pass
    # Return figure.
    return figure


# TODO: accommodate different variable types... categorical, continuous divergent, continuous, ordinal, etc...


def organize_heatmap_asymmetric_data(
    data=None,
):
    """
    Organizes data for chart.

    arguments:
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap

    raises:

    returns:
        (dict): collection of information for chart

    """

    # Copy data.
    data = data.copy(deep=True)
    # Map data's columns to heatmap's rows.
    # Map data's rows to heatmap's columns.
    labels_rows = data.columns.to_list()
    labels_columns = data.index.to_list()
    matrix = numpy.transpose(data.to_numpy())
    # Determine maximal and minimal values.
    array = matrix.flatten()
    minimum = round(numpy.nanmin(array), 2)
    maximum = round(numpy.nanmax(array), 2)
    # Collect information.
    bin = dict()
    bin["data"] = data
    bin["matrix"] = matrix
    bin["minimum"] = minimum
    bin["maximum"] = maximum
    bin["labels_rows"] = labels_rows
    bin["labels_columns"] = labels_columns
    # Return information.
    return bin


def plot_heatmap_asymmetric(
    data=None,
    title=None,
    label_scale=None,
    type=None,
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
        title (str): title for top left of figure
        label_scale (str): label for heatmap scale bar
        type (str): type of property, category, binary, ordinal, continuous, or
            continuous_divergent
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
    bucket = organize_heatmap_asymmetric_data(
        data=data,
    )
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811), # 40 cm, 30 cm
        tight_layout=True,
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
    axes = matplotlib.pyplot.axes()
     # Main chart in bottom panel.
    organize_heatmap_asymmetric_master_main_bottom(
        label_scale=label_scale,
        type=type,
        matrix=bucket["matrix"],
        minimum=bucket["minimum"],
        maximum=bucket["maximum"],
        labels_rows=bucket["labels_rows"],
        labels_columns=bucket["labels_columns"],
        fonts=fonts,
        colors=colors,
        axis_main=axes,
        axis_scale=axes,
        figure=figure,
     )

    return figure


# Asymmetric heatmap with master and main panels
# Sort by master variable and cluster by similarities


def organize_data_master_main_sort_cluster(
    type_master=None,
    sequence=None,
    data=None,
    index=None,
):
    """
    Organize sort and cluster of data.

    arguments:
        type_master (str): type of master variable, category, binary, ordinal,
            or continuous
        data (object): Pandas data frame of main and master feature variables
            across instances
        index (str): name of index that is common between data_master and
            data_main
        sequence (str): method to determine sequence of instances, sort by
            master variable's values or cluster by similarities across features


    raises:

    returns:
        (object): Pandas data frame of pan-tissue signals across genes and
            persons with property

    """

    if (sequence == "sort"):
        data.sort_values(
            by=["master"],
            axis="index",
            ascending=True,
            inplace=True,
        )
        data_cluster_columns = utility.cluster_data_columns(
            data=data,
        )
        data_cluster_columns.reset_index(
            level=None,
            inplace=True
        )
        data_cluster_columns.set_index(
            [index],
            append=False,
            drop=True,
            inplace=True
        )
        if (type_master in ["category", "binary", "ordinal"]):
            data_cluster_columns.reset_index(
                level=None,
                inplace=True
            )
            data_cluster_columns.set_index(
                [index, "master"],
                append=False,
                drop=True,
                inplace=True
            )
            data_cluster_rows = utility.cluster_data_rows_by_group(
                group="master",
                index=index,
                data=data_cluster_columns,
            )
            # Organize data.
            data_cluster_rows.reset_index(
                level=None,
                inplace=True
            )
            data_cluster_rows.set_index(
                [index],
                append=False,
                drop=True,
                inplace=True
            )
            data_sequence = data_cluster_rows.copy(deep=True)
        else:
            data_sequence = data_cluster_columns.copy(deep=True)
    elif (sequence == "cluster"):
        data_cluster_columns = utility.cluster_data_columns(
            data=data,
        )
        data_cluster_rows = utility.cluster_data_rows(
            data=data_cluster_columns,
        )
        data_sequence = data_cluster_rows.copy(deep=True)
    # Return information.
    return data_sequence


def organize_data_master_main(
    data_master=None,
    master=None,
    type_master=None,
    data_main=None,
    type_main=None,
    scale_unit_main=None,
    columns_main_scale_unit=None,
    fill_missing=None,
    index=None,
    sequence=None,
):
    """
    Organize information to sort multiple main data variables by a single
    master variable.

    Data have features across columns and instances across rows.

    Sequence of features depends on hierarchical cluster by their
    similarities across instances.
    Sequence of instances depends either on sort by values of master variable
    or on hierarchical cluster by their similarities across features.

    arguments:
        data_master (object): Pandas data frame including master variable
            across instances that match those of data_main
        master (str): name of feature in data_master that is master variable
        type_master (str): type of master variable, category, binary, ordinal,
            or continuous
        data_main (object): Pandas data frame of feature variables across
            instances that match those of data_master
        type_main (str): type of main variables, category, binary, ordinal,
            continuous, or continuous_divergent
        scale_unit_main (bool): whether to scale columns in data_main
        columns_main_scale_unit (list<str>): names of columns in data_main to
            scale to unit distribution between 0 and 1
        fill_missing (bool): whether to fill missing values with zero
        index (str): name of index that is common between data_master and
            data_main
        sequence (str): method to determine sequence of instances, sort by
            master variable's values or cluster by similarities across features

    raises:

    returns:
        (dict): information for chart

    """

    # Copy data.
    data_master = data_master.copy(deep=True)
    data_main = data_main.copy(deep=True)

    # Organize main variables.
    # Drop any rows or columns in main data with only missing values.
    data_main.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    data_main.dropna(
        axis="columns",
        how="all",
        inplace=True,
    )
    # Scale main variable values.
    if scale_unit_main:
        data_main_scale = data_main.apply(
            lambda column:
                sklearn.preprocessing.minmax_scale(
                    column.to_numpy(),
                    feature_range=(0, 1),
                    axis=0,
                    copy=True,
                ) if (column.name in columns_main_scale_unit) else column,
            axis="index",
        )
    else:
        data_main_scale = data_main
    # Replace missing values in main data with zero.
    if fill_missing:
        # Set any infinite values to missing.
        data_main_scale[data_main_scale == numpy.inf] = numpy.nan
        data_main_scale.fillna(
            value=0.0,
            #axis="columns",
            inplace=True,
        )
    else:
        data_main_scale.fillna(
            value=-1.0,
            #axis="columns",
            inplace=True,
        )
    # Organize master variable.
    if (type_master in ["category", "binary", "ordinal"]):
        data_master["master"], labels_categories_master = pandas.factorize(
            data_master[master],
            sort=True
        )
        labels_categories_master = list(labels_categories_master)
    elif type_master == "continuous":
        data_master["master"] = data_master[master]
        labels_categories_master = list()
    data_master = data_master.loc[
        :, data_master.columns.isin(["master", master])
    ]
    # Join master and main data.
    data_hybrid = data_main_scale.join(
        data_master,
        how="left",
        on=index
    )
    data_hybrid.drop(
        labels=[master],
        axis="columns",
        inplace=True
    )
    # Drop any rows in hybrid data with missing values in master.
    data_hybrid.dropna(
        axis="index",
        how="all",
        subset=["master"],
        inplace=True,
    )
    # Sort and cluster data.
    if sequence != "none":
        data_hybrid_sequence = organize_data_master_main_sort_cluster(
            type_master=type_master,
            sequence=sequence,
            data=data_hybrid,
            index=index,
        )
    else:
        data_hybrid_sequence = data_hybrid
    # Compile information.
    bin = dict()
    bin["labels_categories_master"] = labels_categories_master
    bin["data"] = data_hybrid_sequence
    # Return information
    return bin


def organize_heatmap_asymmetric_master_main_data(
    data=None,
):
    """
    Organizes data for chart.

    arguments:
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap

    raises:

    returns:
        (dict): collection of information for chart

    """

    # Copy data.
    data_master = data.copy(deep=True)
    data_main = data.copy(deep=True)
    # Groups.
    values_master = data_master["master"].to_list()
    values_master_unique = utility.collect_unique_elements(
        elements_original=values_master,
    )
    # Map data's columns to heatmap's rows.
    # Map data's rows to heatmap's columns.
    data_main.drop(
        labels=["master"],
        axis="columns",
        inplace=True
    )
    labels_rows_main = data_main.columns.to_list()
    labels_columns_main = data_main.index.to_list()
    matrix_main = numpy.transpose(data_main.to_numpy())
    # Determine maximal and minimal values.
    minimum_master = round(min(values_master), 2)
    maximum_master = round(max(values_master), 2)
    array_main = matrix_main.flatten()
    minimum_main = round(numpy.nanmin(array_main), 2)
    maximum_main = round(numpy.nanmax(array_main), 2)

    # Collect information.
    bin = dict()
    bin["data_master"] = data_master
    bin["values_master"] = values_master
    bin["values_master_unique"] = values_master_unique
    bin["minimum_master"] = minimum_master
    bin["maximum_master"] = maximum_master
    bin["data_main"] = data_main
    bin["matrix_main"] = matrix_main
    bin["minimum_main"] = minimum_main
    bin["maximum_main"] = maximum_main
    bin["labels_rows_main"] = labels_rows_main
    bin["labels_columns_main"] = labels_columns_main
    # Return information.
    return bin


def initialize_heatmap_asymmetric_master_main_figure_axes(
    title=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type heatmap.

    arguments:
        title (str): title for chart
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object, object): instances of figure and axes objects

    """

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
            left=0.11,
            right=0.94,
            top=0.94,
            bottom=0.05,
        ),
    )
    # Return information.
    return figure, axes


def organize_heatmap_asymmetric_master_main_top(
    title_subordinate=None,
    label=None,
    type=None,
    values=None,
    values_unique=None,
    labels_categories=None,
    minimum=None,
    maximum=None,
    fonts=None,
    colors=None,
    axes=None,
    figure=None,
):
    """
    Organizes top panel of figure.

    arguments:
        title_subordinate (str): subordinate title for top right of figure
        label (str): label for master heatmap
        type (str): type of property, category or continuity
        values (list): values of property
        values_unique (list): unique values of property
        labels_categories (list<str>): labels for categorical ticks
        minimum (float): minimal value of property
        maximum (float): maximal value of property
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        axes (object): instance axis object
        figure (object): instance figure object

    raises:

    returns:
        (object): figure object

    """

    # Define color map and labels.
    if ((type == "category") or (type == "binary") or (type == "ordinal")):
        #color_map = matplotlib.colors.ListedColormap(
        #    [colors["blue"], colors["orange"]], 2
        #)
        color_map = matplotlib.pyplot.get_cmap("GnBu", len(values_unique))
        ticks = values_unique
        labels_ticks = labels_categories
    elif (type == "continuous"):
        color_map = "GnBu"
        ticks = [minimum, maximum]
        labels_ticks = [minimum, maximum]
    # Initialize chart.
    image = axes[0, 0].imshow(
        [values],
        cmap=color_map,
        vmin=minimum,
        vmax=maximum,
        aspect="auto",
        origin="upper",
        # Extent: (left, right, bottom, top)
        extent=(-0.5, (len(values) - 0.5), (1 + 0.5), -0.5),
    )
    axes[0, 0].set_yticks(numpy.arange(1))
    axes[0, 0].set_yticklabels(
        [str(label)],
        #minor=False,
        ha="right", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    axes[0, 0].set_xlabel(
        str(title_subordinate),
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
    bar = figure.colorbar(
        image,
        cax=axes[0, 1],
        ticks=ticks,
        orientation="vertical",
        use_gridspec=True,
    )
    bar.ax.set_yticklabels(
        labels_ticks,
        #minor=False,
        ha="left", # horizontal alignment
        va="bottom", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    bar.ax.tick_params(
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
    pass


def organize_heatmap_asymmetric_master_main_bottom(
    label_scale=None,
    type=None,
    matrix=None,
    minimum=None,
    maximum=None,
    labels_rows=None,
    labels_columns=None,
    fonts=None,
    colors=None,
    axis_main=None,
    axis_scale=None,
    figure=None,
):
    """
    Organizes top panel of figure.

    arguments:
        label_scale (str): label for heatmap scale bar
        type (str): type of property, category, binary, ordinal, continuous, or
            continuous_divergent
        matrix (array<array>): values of properties
        minimum (float): minimal value of property
        maximum (float): maximal value of property
        labels_rows (list<str>): labels for matrix rows
        labels_columns (list<str>): labels for matrix columns
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        axis_main (object): instance of axis for main plot
        axis_scale (object): instance of axis for scale bar plot
        figure (object): instance figure object

    raises:

    returns:
        (object): figure object

    """

    # Define color map and labels.
    if (
        (type == "category") or
        (type == "binary") or
        (type == "ordinal") or
        (type == "continuous")
    ):
        color_map = "RdPu"
        scale = matplotlib.colors.Normalize(
            vmin=minimum,
            vmax=maximum,
        )
        ticks = [minimum, maximum]
        labels_ticks = [minimum, maximum]
    elif (type == "continuous_divergent"):
        color_map = "PuOr"
        scale = matplotlib.colors.TwoSlopeNorm(
            vmin=minimum,
            vcenter=0.0,
            vmax=maximum,
        )
        ticks = [minimum, 0.0, maximum]
        labels_ticks = [minimum, 0.0, maximum]

    # Initialize chart.
    image = axis_main.imshow(
        matrix,
        cmap=color_map, # sequential: "RdPu", diverging: "PuOr"
        aspect="auto", # "auto", "equal"
        origin="upper", # "lower" or "upper"
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
        alpha=1.0,
        norm=scale,
        filternorm=True,
        resample=True,
    )
    #axes[2, 0].set_xlim(-10, (matrix.shape[1] + 10))
    #axes[2, 0].set_ylim(-3, (matrix.shape[0] + 3))

    # Create ticks and labels for each grid.
    # Let the horizontal axes labeling appear on top.
    if matrix.shape[0] <= 100:
        if (100 >= matrix.shape[0] and matrix.shape[0] >= 90):
            size_count = "ten"
        elif (90 > matrix.shape[0] and matrix.shape[0] >= 75):
            size_count = "nine"
        elif (75 > matrix.shape[0] and matrix.shape[0] >= 50):
            size_count = "eight"
        elif (50 > matrix.shape[0] and matrix.shape[0] >= 25):
            size_count = "seven"
        elif (25 > matrix.shape[0]):
            size_count = "three"
        axis_main.tick_params(
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
        axis_main.set_yticks(numpy.arange(matrix.shape[0]))
        axis_main.set_yticklabels(
            labels_rows,
            #minor=False,
            ha="right", # horizontal alignment
            va="center", # vertical alignment
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_count]
        )
        pass
    # Create legend for color map.
    if axis_scale is axis_main:
        bar = axis_main.figure.colorbar(
            image,
            orientation="vertical",
            ax=axis_main,
            ticks=ticks,
        )
    else:
        bar = figure.colorbar(
            image,
            cax=axis_scale,
            ticks=ticks,
            orientation="vertical",
            use_gridspec=True,
        )
        pass
    bar.ax.set_ylabel(
        label_scale,
        rotation=-90,
        ha="center",
        va="bottom",
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["three"],
    )
    bar.ax.yaxis.set_label_coords(2.5, 0.5)
    bar.ax.set_yticklabels(
        labels_ticks,
        #minor=False,
        ha="left", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    bar.ax.tick_params(
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
    pass


def plot_heatmap_asymmetric_master_main_top_bottom(
    title=None,
    title_subordinate=None,
    label_master=None,
    labels_categories_master=None,
    label_main=None,
    type_master=None,
    type_main=None,
    data=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type heatmap.

    Data must have observations across rows.
    Data's observations must already be in sort order.
    Data must have features across columns.
    Data must also have a single column of name "master".
    Values in column "master" must be integers.

    The chart will organize features across rows and observations across
    columns.
    The chart will represent integer values of the group of each observation in
    a separate chart across columns.


    arguments:
        title (str): title for top left of figure
        title_subordinate (str): subordinate title for top right of figure
        label_master (str): label for left of master heatmap row
        labels_categories_master (list<str>): labels for scale ticks of
            categorical master variable
        label_main (str): label for main heatmap scale
        type_master (str): type of master variable, category, binary, ordinal,
            or continuous
        type_main (str): type of main variables, category, binary, ordinal,
            continuous, or continuous_divergent
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    bin_data = organize_heatmap_asymmetric_master_main_data(
        data=data,
    )
    # Initialize figure.
    figure, axes = initialize_heatmap_asymmetric_master_main_figure_axes(
        title=title,
        fonts=fonts,
        colors=colors,
    )

    # Master chart in top panel.
    organize_heatmap_asymmetric_master_main_top(
        title_subordinate=title_subordinate,
        label=label_master,
        type=type_master,
        values=bin_data["values_master"],
        values_unique=bin_data["values_master_unique"],
        labels_categories=labels_categories_master,
        minimum=bin_data["minimum_master"],
        maximum=bin_data["maximum_master"],
        fonts=fonts,
        colors=colors,
        axes=axes,
        figure=figure,
     )

     # Main chart in bottom panel.
    organize_heatmap_asymmetric_master_main_bottom(
        label_scale=label_main,
        type=type_main,
        matrix=bin_data["matrix_main"],
        minimum=bin_data["minimum_main"],
        maximum=bin_data["maximum_main"],
        labels_rows=bin_data["labels_rows_main"],
        labels_columns=bin_data["labels_columns_main"],
        fonts=fonts,
        colors=colors,
        axis_main=axes[1, 0],
        axis_scale=axes[1, 1],
        figure=figure,
     )

    # Return figure.
    return figure


###


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

    if color_count == 1:
        colors_series = [colors["blue"]]
    elif color_count == 2:
        colors_series = [colors["blue"], colors["orange"]]
    elif color_count > 2:
        colors_series = list(
            seaborn.color_palette("hls", n_colors=color_count)
        )
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
        rotation=-30,
        rotation_mode="anchor",
        ha="left", # horizontal alignment
        va="top", # vertical alignment
        minor=False,
        #rotation=rotation,
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


def plot_boxes(
    arrays=None,
    labels_groups=None,
    label_vertical=None,
    label_horizontal=None,
    label_top_center=None,
    label_top_left=None,
    label_top_right=None,
    orientation=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type histogram to represent the frequency
    distribution of a single series of values.

    arguments:
        arrays (list<array>): NumPy arrays of values for each group
        labels_groups (list<str>): labels for each group
        label_vertical (str): label for vertical axis
        label_horizontal (str): label for horizontal axis
        label_top_center (str): label for top center of plot area
        label_top_left (str): label for top left of plot area
        label_top_right (str): label for top right of plot area
        orientation (str): orientation of figure, either "portrait" or
            "landscape"
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Define colors.
    color_count = len(arrays)
    if color_count == 1:
        colors_series = [colors["blue"]]
    elif color_count == 2:
        colors_series = [colors["gray"], colors["blue"]]
    elif color_count == 3:
        colors_series = [colors["gray"], colors["blue"], colors["orange"]]
    elif color_count > 3:
        colors_series = list(
            seaborn.color_palette("hls", n_colors=color_count)
        )
    # Create figure.
    if orientation == "portrait":
        figure = matplotlib.pyplot.figure(
            figsize=(11.811, 15.748),
            tight_layout=True
        )
    elif orientation == "landscape":
        figure = matplotlib.pyplot.figure(
            figsize=(15.748, 11.811),
            tight_layout=True
        )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    # Create boxes.
    handle = axes.boxplot(
        arrays,
        notch=False,
        vert=True,
        widths=0.7,
        patch_artist=True,
        labels=labels_groups,
        manage_ticks=True,
    )
    # Fill boxes with colors.
    for box_patch, color_box in zip(handle["boxes"], colors_series):
        box_patch.set_facecolor(color_box)
        pass
    # Label axes.
    axes.set_ylabel(
        ylabel=label_vertical,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    if len(label_horizontal) > 0:
        axes.set_xlabel(
            xlabel=label_horizontal,
            labelpad=20,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["one"],
            rotation="horizontal",
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
    if len(label_top_center) > 0:
        axes.text(
            0.5,
            0.9,
            label_top_center,
            horizontalalignment="center",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["one"]
        )
    if len(label_top_left) > 0:
        axes.text(
            0.25,
            0.9,
            label_top_left,
            horizontalalignment="center",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["two"]
        )
    if len(label_top_right) > 0:
        axes.text(
            0.75,
            0.9,
            label_top_right,
            horizontalalignment="center",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["two"]
        )
    return figure


def plot_scatter_factor_groups(
    data=None,
    abscissa=None,
    ordinate=None,
    label_horizontal=None,
    label_vertical=None,
    factor=None,
    label_factor=None,
    labels_factors=None,
    fonts=None,
    colors=None,
    point_size=None,
    plot_factor_labels=None,
    legend=None,
):
    """
    Creates a figure of a chart of type scatter to represent the association
    of two variables.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with variable for horizontal (x)
            axis
        ordinate (str): name of data column with variable for vertical (y) axis
        label_horizontal (str): label for horizontal axis
        label_vertical (str): label for vertical axis
        factor (str): name of data column with categorical variable to
            distinguish groups of instances
        label_factor (str): label to describe the factor
        labels_factors (list<str>): labels to describe the factor's values
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        point_size (float): size for scatter points
        plot_factor_labels (bool): whether to plot factor labels on chart
        legend (bool): whether to include a legend for series on the chart

    raises:

    returns:
        (object): figure object

    """

    ##########
    # Organize data.
    data = data.copy(deep=True)
    data.reset_index(
        level=None,
        inplace=True
    )
    data.set_index(
        factor,
        append=False,
        drop=True,
        inplace=True
    )
    groups = data.groupby(level=[factor])
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
            markersize=point_size,
            markeredgecolor=colors_series[index],
            markerfacecolor=colors_series[index]
        )
        index += 1
        pass
    # Plot labels for each group.
    labels = []
    index = 0
    for name, group in groups:
        if plot_factor_labels:
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
        label = str(str(index+1) + ": " + str(labels_factors[index]))
        labels.append(label)
        index += 1
        pass
    # Create legend.
    # Create custome elements for the legend.
    if legend:
        elements = create_legend_elements(
            colors=colors_series,
            labels=labels
        )
        axes.legend(
            handles=elements,
            loc="upper right",
            prop=fonts["properties"]["four"],
            title=label_factor,
            title_fontsize=fonts["values"]["three"]["size"]
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
        abscissa (str): name of data column with variable for horizontal (x)
            axis
        ordinate (str): name of data column with variable for vertical (y) axis
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

# TODO: probably obsolete?
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


def plot_scatter_label_emphasis_points(
    emphasis_keys=None,
    label_keys=None,
    column_key=None,
    column_label=None,
    line_abscissa=None,
    line_ordinate=None,
    line_ordinate_origin=None,
    data=None,
    abscissa=None,
    ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type scatter with thresholds on each
        dimension.

    arguments:
        emphasis_keys (list<str>): keys for rows of points for emphasis
        label_keys (list<str>): keys for rows of points for labels
        column_key (str): name of column in data with keys
        column_label (str): name of column in data with labels
        line_abscissa (float): point on abscissa for horizontal line
        line_ordinate (float): point on ordinate for vertical line
        line_ordinate_origin (bool): whether to draw vertical origin line
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with independent variable
        ordinate (str): name of data column with dependent variable
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    data = data.copy(deep=True)
    ordinate_minimum = min(data[ordinate].to_list())
    ordinate_maximum = max(data[ordinate].to_list())
    data_bore = data.loc[~data[column_key].isin(emphasis_keys), :]
    data_emphasis = data.loc[data[column_key].isin(emphasis_keys), :]
    data_label = data.loc[data[column_key].isin(label_keys), :]
    data_bore.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    data_emphasis.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    data_label.dropna(
        axis="index",
        how="any",
        inplace=True,
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
        data_bore[abscissa].to_numpy(),
        data_bore[ordinate].to_numpy(),
        linestyle="",
        marker="o",
        markersize=5,
        markeredgecolor=colors["gray"],
        markerfacecolor=colors["gray"]
    )
    handle = axes.plot(
        data_emphasis[abscissa].to_numpy(),
        data_emphasis[ordinate].to_numpy(),
        linestyle="",
        marker="o",
        markersize=7,
        markeredgecolor=colors["blue"],
        markerfacecolor=colors["blue"]
    )
    handle = axes.plot(
        data_label[abscissa].to_numpy(),
        data_label[ordinate].to_numpy(),
        linestyle="",
        marker="o",
        markersize=9,
        markeredgecolor=colors["orange"],
        markerfacecolor=colors["orange"]
    )

    # Plot lines for each threshold value...
    # Create lines for thresholds.
    if line_ordinate_origin:
        axes.axvline(
            x=0.0,
            ymin=0,
            ymax=1,
            alpha=1.0,
            color=colors["black"],
            linestyle="--",
            linewidth=3.0,
        )
    if line_abscissa is not None:
        axes.axvline(
            x=line_abscissa,
            ymin=0,
            ymax=1,
            alpha=1.0,
            color=colors["orange"],
            linestyle="--",
            linewidth=5.0,
        )
    if line_ordinate is not None:
        axes.axhline(
            y=line_ordinate,
            xmin=0,
            xmax=1,
            alpha=1.0,
            color=colors["orange"],
            linestyle="--",
            linewidth=5.0,
        )

    # Plot labels.
    # Place bottom of label above point by 5% of maximal y value.
    for label_key in label_keys:
        data_label = data.loc[data[column_key].isin([label_key]), :]
        if (data_label.shape[0] > 0):
            for index_point in data_label.index.to_list():
                axes.text(
                    (data_label.at[index_point, abscissa]),
                    (
                        data_label.at[index_point, ordinate] +
                        (0.02 * ordinate_maximum)
                    ),
                    data_label[column_label].to_list()[0],
                    backgroundcolor=colors["white_faint"],
                    color=colors["black"],
                    fontproperties=fonts["properties"]["three"],
                    horizontalalignment="center",
                    verticalalignment="bottom"
                )
                pass
            pass
        pass
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
    fonts=None,
    colors=None,
):
    """
    Creates a Venn Diagram to represent overlap between multiple sets.

    arguments:
        sets (dict<list<str>>): values in sets
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
    if len(names) == 2:
        venn = matplotlib_venn.venn2(
            subsets=[
                set(sets[names[0]]),
                set(sets[names[1]]),
            ],
            set_labels=names,
        )
    elif len(names) == 3:
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


def read_source_persons_health_collection_variables(
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
        dock, "selection", "tight", "persons_properties",
        "persons_sets.pickle"
    )
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", "selection",
        "data_persons_properties.pickle"
    )

    # Read information from file.
    with open(path_collections_health, "rb") as file_source:
        collections_health_variables = pickle.load(file_source)
    with open(path_persons_selection, "rb") as file_source:
        persons_selection = pickle.load(file_source)
    with open(path_persons_sets, "rb") as file_source:
        persons_sets = pickle.load(file_source)
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    # Compile and return information.
    return {
        "collections_health_variables": collections_health_variables,
        "persons_selection": persons_selection,
        "persons_sets": persons_sets,
        "data_persons_properties": data_persons_properties,
    }


def plot_chart_persons_health_collection_variables(
    title=None,
    label=None,
    master=None,
    type_master=None,
    type_main=None,
    labels_categories_master=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        type_main (str): type of main variables
        labels_categories_master (list<str>): labels for scale ticks of
            categorical master variable
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(title + ".png")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_heatmap_asymmetric_master_main_top_bottom(
        title=title,
        title_subordinate="",
        label_master=label,
        labels_categories_master=labels_categories_master,
        label_main="clinical annotations across persons",
        type_master=type_master,
        type_main=type_main,
        data=data,
        fonts=fonts,
        colors=colors,
    )

    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_persons_health_collection_variables_property(
    title=None,
    label=None,
    master=None,
    type_master=None,
    type_main=None,
    health_collection=None,
    data_health_collection=None,
    persons=None,
    data_persons_properties=None,
    path_directory=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        type_main (str): type of main variables
        health_collection (str): name of collection of health variables
        data_health_collection (object): Pandas data frame of health variables
            across persons
        persons (list<str>): identifiers of persons
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        path_directory (str): path to directory

    raises:

    returns:

    """

    # Copy data.
    data_health_collection = data_health_collection.copy(deep=True)
    # Select data for persons.
    data_health_persons = data_health_collection.loc[
        data_health_collection.index.isin(persons), :
    ]
    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_data_master_main(
        data_master=data_persons_properties,
        master=master,
        type_master=type_master,
        data_main=data_health_persons,
        type_main=type_main,
        scale_unit_main=True,
        columns_main_scale_unit=data_health_persons.columns.tolist(),
        fill_missing=True,
        index="person",
        sequence="sort",
    )
    # Create charts for the gene.
    plot_chart_persons_health_collection_variables(
        title=title,
        label=label,
        master=master,
        type_master=type_master,
        type_main=type_main,
        labels_categories_master=bin["labels_categories_master"],
        data=bin["data"],
        path_directory=path_directory,
    )
    pass


def define_parameters_persons_health_collection_variables():
    """
    Defines parameters for plots of persons' properties.

    arguments:

    raises:

    returns:
        (dict): collection of parameters

    """

    parameters = dict()
    if False:
        parameters["sex_text"] = dict(
            title="sex", label="sex", property="sex_text", type="category",
            master=True, main=False,
        )
        parameters["age_grade"] = dict(
            title="age", label="age_grade", property="age", type="ordinal",
            master=True, main=True,
        )
        parameters["ventilation_duration_grade"] = dict(
            title="ventilation", label="ventilation",
            property="ventilation_duration_grade", type="ordinal",
            master=True, main=True,
        )
    parameters["ventilation_binary"] = dict(
        title="ventilation_binary", label="ventilation",
        property="ventilation_binary", type="binary",
        master=True, main=False,
    )
    # Return information.
    return parameters


def prepare_charts_persons_health_collection_variables(
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
    source = read_source_persons_health_collection_variables(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_assembly = os.path.join(path_plot, "assembly")
    path_directory = os.path.join(
        path_assembly, "persons_health_variables"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Define properties for master parameters.
    parameters = define_parameters_persons_health_collection_variables()
    # Iterate on collections of health variables.
    for collection in source["collections_health_variables"].keys():
        # Iterate on categorical variables.
        for parameter in parameters.keys():
            # Prepare charts for properties.
            print(collection)
            print(parameters[parameter]["title"])
            title = str(collection + "_" + parameters[parameter]["title"])
            prepare_charts_persons_health_collection_variables_property(
                title=title,
                label=parameters[parameter]["label"],
                master=parameters[parameter]["property"],
                type_master=parameters[parameter]["type"],
                type_main="binary",
                health_collection=collection,
                data_health_collection=(
                    source["collections_health_variables"][collection]
                ),
                persons=source["persons_sets"]["selection"],
                data_persons_properties=source["data_persons_properties"],
                path_directory=path_directory,
            )
            pass
    pass


##########
# Selection of health history variables for persons respiration, inflammation,
# infection, and steroid.
# heatmaps
# Status: working


def read_source_persons_properties_adjacency(
    group=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        group (str): group of persons, either selection or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_persons_sets = os.path.join(
        dock, "selection", "tight", "persons_properties",
        "persons_sets.pickle"
    )
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", group,
        "data_persons_properties.pickle"
    )

    # Read information from file.
    with open(path_persons_sets, "rb") as file_source:
        persons_sets = pickle.load(file_source)
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    # Compile and return information.
    return {
        "persons_sets": persons_sets,
        "data_persons_properties": data_persons_properties,
    }


def plot_chart_persons_properties_adjacency(
    title=None,
    label=None,
    master=None,
    type_master=None,
    labels_categories_master=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        labels_categories_master (list<str>): labels for scale ticks of
            categorical master variable
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(title + ".png")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_heatmap_asymmetric_master_main_top_bottom(
        title=title,
        title_subordinate="",
        label_master=label,
        labels_categories_master=labels_categories_master,
        label_main="properties across persons",
        type_master=type_master,
        type_main="continuous",
        data=data,
        fonts=fonts,
        colors=colors,
    )

    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def select_translate_persons_properties_data_main(
    persons=None,
    master=None,
    mains=None,
    parameters=None,
    data_persons_properties=None,
):
    """
    Organize information for chart.

    arguments:
        persons (list<str>): identifiers of persons
        master (str): name of feature from persons' properties to use for
            master variable
        mains (list<str>): names of properties for main data
        parameters (dict): collection of parameters
        data_persons_properties (object): Pandas data frame of persons and
            their properties

    raises:

    returns:
        (object): Pandas data frame of persons and their properties
    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    # Select data for persons and properties of interest.
    data_persons = data_persons_properties.loc[
        data_persons_properties.index.isin(persons), :
    ]
    columns_main = list(filter(lambda value: value != master, mains))
    data_properties = data_persons.loc[
        :, data_persons.columns.isin(columns_main)
    ]
    # Organize signals.
    # Translate identifiers of genes.
    #identifiers = data_signals_genes_persons.columns.to_list()
    translations = dict()
    for property in mains:
        translations[property] = parameters[property]["label"]
        pass
    data_properties.rename(
        columns=translations,
        inplace=True,
    )
    # Return information.
    return data_properties


def prepare_charts_persons_properties_adjacency_property(
    title=None,
    label=None,
    master=None,
    type_master=None,
    parameters=None,
    mains=None,
    mains_labels=None,
    persons=None,
    data_persons_properties=None,
    path_sort=None,
    path_cluster=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        parameters (dict): collection of parameters
        mains (list<str>): names of properties for main data
        mains_labels (list<str>): labels of properties for main data
        persons (list<str>): identifiers of persons
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        path_sort (str): path to directory
        path_cluster (str): path to directory

    raises:

    returns:

    """

    # Select and translate data for genes.
    data_main = select_translate_persons_properties_data_main(
        persons=persons,
        master=master,
        mains=mains,
        parameters=parameters,
        data_persons_properties=data_persons_properties,
    )

    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_data_master_main(
        data_master=data_persons_properties,
        master=master,
        type_master=type_master,
        data_main=data_main,
        type_main="countinuous",
        scale_unit_main=True,
        columns_main_scale_unit=mains_labels,
        fill_missing=True,
        index="person",
        sequence="sort",
    )
    # Create charts for the gene.
    plot_chart_persons_properties_adjacency(
        title=title,
        label=label,
        master=master,
        type_master=type_master,
        labels_categories_master=bin["labels_categories_master"],
        data=bin["data"],
        path_directory=path_sort,
    )
    pass


def define_parameters_persons_properties_adjacency():
    """
    Defines parameters for plots of persons' properties.

    arguments:

    raises:

    returns:
        (dict): collection of parameters

    """

    parameters = dict()
    parameters["sex_y"] = dict(
        title="sex", label="sex_y", property="sex_y", type="binary",
        master=False, main=True,
    )
    parameters["sex_text"] = dict(
        title="sex", label="sex", property="sex_text", type="category",
        master=True, main=False,
    )
    parameters["age"] = dict(
        title="age", label="age", property="age", type="continuous",
        master=False, main=True,
    )
    parameters["age_grade"] = dict(
        title="age", label="age_grade", property="age", type="ordinal",
        master=True, main=True,
    )
    parameters["body"] = dict(
        title="BMI", label="BMI", property="body", type="continuous",
        master=False, main=True,
    )
    parameters["climate"] = dict(
        title="season", label="season", property="climate", type="ordinal",
        master=False, main=True,
    )
    parameters["hardiness"] = dict(
        title="hardiness", label="hardiness", property="hardiness",
        type="ordinal",
        master=False, main=True,
    )
    parameters["respiration_binary"] = dict(
        title="respiration", label="respiration",
        property="respiration_binary", type="binary",
        master=False, main=True,
    )
    parameters["smoke"] = dict(
        title="smoke", label="smoke", property="smoke", type="ordinal",
        master=True, main=True,
    )
    parameters["inflammation_binary"] = dict(
        title="inflammation", label="inflammation",
        property="inflammation_binary", type="binary",
        master=False, main=True,
    )
    parameters["leukocyte_binary"] = dict(
        title="leukocyte", label="leukocyte",
        property="leukocyte_binary", type="binary",
        master=False, main=True,
    )
    parameters["infection_binary"] = dict(
        title="infection", label="infection",
        property="infection_binary", type="binary",
        master=False, main=True,
    )
    parameters["mononucleosis_binary"] = dict(
        title="CMV-EBV", label="CMV-EBV",
        property="mononucleosis_binary", type="binary",
        master=False, main=True,
    )
    parameters["steroid_binary"] = dict(
        title="steroid", label="steroid",
        property="steroid_binary", type="binary",
        master=False, main=True,
    )
    parameters["heart_binary"] = dict(
        title="circulation", label="circulation",
        property="heart_binary", type="binary",
        master=False, main=True,
    )
    parameters["diabetes_binary"] = dict(
        title="diabetes", label="diabetes",
        property="diabetes_binary", type="binary",
        master=False, main=True,
    )
    parameters["ventilation_duration_grade"] = dict(
        title="ventilation", label="ventilation",
        property="ventilation_duration_grade", type="ordinal",
        master=True, main=True,
    )
    parameters["ventilation_binary"] = dict(
        title="ventilation_binary", label="ventilation",
        property="ventilation_binary", type="binary",
        master=True, main=False,
    )
    parameters["delay"] = dict(
        title="delay", label="delay",
        property="delay", type="continuous",
        master=False, main=True,
    )
    parameters["refrigeration_duration"] = dict(
        title="refrig_dur", label="refrig_dur",
        property="refrigeration_duration", type="continuous",
        master=False, main=True,
    )
    # Return information.
    return parameters


def organize_parameters_persons_properties_adjacency():
    """
    Organizes parameters for plots of persons' properties.

    arguments:

    raises:

    returns:
        (dict): collection of parameters

    """

    # Define parameters.
    parameters = define_parameters_persons_properties_adjacency()
    # Collect parameters.
    masters = list()
    mains = list()
    mains_labels = list()
    for property in parameters:
        if parameters[property]["master"]:
            masters.append(property)
        if parameters[property]["main"]:
            mains.append(property)
            mains_labels.append(parameters[property]["label"])
    # Compile information.
    bin = dict()
    bin["parameters"] = parameters
    bin["masters"] = masters
    bin["mains"] = mains
    bin["mains_labels"] = mains_labels
    # Return information.
    return bin


def prepare_charts_persons_properties_adjacency(
    dock=None,
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
    source = read_source_persons_properties_adjacency(
        group="selection",
        dock=dock,
    )

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

    # Iterate on categorical and ordinal groups of properties.
    bin = organize_parameters_persons_properties_adjacency()
    # Iterate on categorical variables.
    for master in bin["masters"]:
        # Prepare charts for properties.
        prepare_charts_persons_properties_adjacency_property(
            title=bin["parameters"][master]["title"],
            label=bin["parameters"][master]["label"],
            master=master,
            type_master=bin["parameters"][master]["type"],
            parameters=bin["parameters"],
            mains=bin["mains"],
            mains_labels=bin["mains_labels"],
            persons=source["persons_sets"]["selection"],
            data_persons_properties=source["data_persons_properties"],
            path_sort=path_sort,
            path_cluster=path_cluster,
        )
        pass
    pass



##########
# Distributions of each bimodality measure's scores across genes
# Status: working


def read_source_ventilation_duration_distribution(
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
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", "selection",
        "data_persons_properties.pickle"
    )

    # Read information from file.
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)

    # Compile and return information.
    return {
        "data_persons_properties": data_persons_properties,
    }


def plot_chart_ventilation_duration_distribution(
    values=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        values (list<float>): values
        threshold (float): value of threshold for which to draw line
        path_directory (str): path to directory

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Specify directories and files.
    path_file = os.path.join(
        path_directory, "ventilation_duration.svg"
    )

    # Create figure.
    figure = plot_distribution_histogram(
        series=values,
        name="",
        bin_method="count",
        bin_count=70,
        label_bins="duration on ventilation (hours)",
        label_counts="counts of persons per bin",
        fonts=fonts,
        colors=colors,
        line=False,
        position=1,
        text="",
    )
    # Write figure.
    write_figure(
        path=path_file,
        figure=figure
    )

    pass


def prepare_chart_ventilation_duration_distribution(
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
    source = read_source_ventilation_duration_distribution(dock=dock)

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_assembly = os.path.join(path_plot, "assembly")
    path_directory = os.path.join(
        path_assembly, "ventilation_duration"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Prepare charts.
    plot_chart_ventilation_duration_distribution(
        values=source["data_persons_properties"]["ventilation_duration"].to_list(),
        path_directory=path_directory
    )

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


def read_source_selection_dimension_principal_components(
    cohort=None,
    stringency=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        stringency (str): category, loose or tight, of selection criteria
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_genotype_variance = os.path.join(
        dock, "assembly", "genotype", cohort,
        "data_genotype_variance.pickle"
    )
    path_variances_data = os.path.join(
        dock, "selection", stringency, "persons_properties", cohort,
        "dimension", "variances_data.pickle"
    )
    # Read information from file.
    data_genotype_variance = pandas.read_pickle(
        path_data_genotype_variance
    )
    variances_data = pandas.read_pickle(
        path_variances_data
    )
    # Compile and return information.
    return {
        "data_genotype_variance": data_genotype_variance,
        "variances_data": variances_data,
    }


def plot_chart_selection_dimension_principal_components(
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
    write_figure_png(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_selection_dimension_principal_components_cohort(
    cohort=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_selection_dimension_principal_components(
        cohort=cohort,
        stringency="tight",
        dock=dock
    )
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    path_dimension = os.path.join(path_plot, "selection", "dimension")
    path_directory = os.path.join(path_dimension, cohort)
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Organize data.
    types = dict()
    types["component"] = "int32"
    types["variance"] = "float32"
    data_genotype = source["data_genotype_variance"].iloc[0:10,:]
    data_genotype = data_genotype.astype(
        types,
    )
    plot_chart_selection_dimension_principal_components(
        data=data_genotype,
        path_file=os.path.join(path_directory, "genotype.png"),
    )
    # Iterate on categorical variables.
    for property in source["variances_data"].keys():
        # Organize data.
        data = source["variances_data"][property]
        data = data.astype(
            types,
        )
        data["variance"] = data["variance"].apply(
            lambda value: (value * 100)
        )
        plot_chart_selection_dimension_principal_components(
            data=data,
            path_file=os.path.join(path_directory, str(property + ".png")),
        )
    pass


def prepare_charts_selection_dimension_principal_components(
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    prepare_charts_selection_dimension_principal_components_cohort(
        cohort="selection",
        dock=dock,
    )
    prepare_charts_selection_dimension_principal_components_cohort(
        cohort="respiration",
        dock=dock,
    )
    prepare_charts_selection_dimension_principal_components_cohort(
        cohort="ventilation",
        dock=dock,
    )
    pass


##########
# Distributions of each bimodality measure's scores across genes
# Status: working


def read_source_modality_measure_genes_distribution(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_genes_scores = os.path.join(
        dock, "distribution", cohort, "collection", "genes_scores.pickle"
    )
    path_scores = os.path.join(
        dock, "distribution", cohort, "collection", "scores.pickle"
    )
    path_measures_thresholds = os.path.join(
        dock, "candidacy", cohort, "threshold", "measures_thresholds.pickle"
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


def plot_chart_modality_measure_genes_distribution(
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


def prepare_charts_modality_measure_genes_distribution_cohort(
    cohort=None,
    dock=None
):
    """
    Plots charts.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_modality_measure_genes_distribution(
        cohort=cohort,
        dock=dock,
    )
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_candidacy = os.path.join(path_plot, "candidacy")
    path_measure = os.path.join(path_candidacy, "measure")
    path_directory = os.path.join(path_measure, cohort)
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    # Set measures of modality.
    measures = list(source["scores"].keys())
    # Iterate on bimodality measures.
    for measure in measures:
        # Prepare charts.
        # Specify directories and files.
        path_figure = os.path.join(
            path_directory, str(measure + ".svg")
        )
        plot_chart_modality_measure_genes_distribution(
            values=source["scores"][measure]["values"],
            threshold=source["measures_thresholds"][measure],
            path=path_figure
        )
        pass
    pass


def prepare_charts_modality_measure_genes_distribution(
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

    # Organize data and regress cases.
    cohorts = [
        "selection",
        "respiration",
        "ventilation",
    ]
    for cohort in cohorts:
        prepare_charts_modality_measure_genes_distribution_cohort(
            cohort=cohort,
            dock=dock,
        )
    pass



##########
# Thresholds on heritability for selection of genes
# Status: working


def read_source_gene_heritability(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Read genes sets.
    genes_sets = integration.read_source_genes_sets_combination_integration(
        cohort="selection",
        dock=dock,
    )
    # Specify directories and files.
    path_data_genes_heritability = os.path.join(
        dock, "heritability", cohort, "collection",
        "data_genes_heritability.pickle"
    )
    path_threshold_proportion = os.path.join(
        dock, "heritability", cohort, "collection",
        "threshold_proportion.pickle"
    )
    path_threshold_probability_log = os.path.join(
        dock, "heritability", cohort, "collection",
        "threshold_probability_log.pickle"
    )

    # Read information from file.
    data_genes_heritability = pandas.read_pickle(
        path_data_genes_heritability
    )
    with open(path_threshold_proportion, "rb") as file_source:
        threshold_proportion = pickle.load(file_source)
    with open(path_threshold_probability_log, "rb") as file_source:
        threshold_probability_log = pickle.load(file_source)

    # Compile and return information.
    return {
        "genes_sets": genes_sets,
        "data_genes_heritability": data_genes_heritability,
        "threshold_proportion": threshold_proportion,
        "threshold_probability_log": threshold_probability_log,
    }


def organize_gene_heritability(
    data_genes_heritability=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        data_genes_heritability (object): Pandas data frame of genes'
            heritabilities

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    # Organize data.
    data = data_genes_heritability.copy(deep=True)
    data.reset_index(
        level=None,
        inplace=True
    )
    # Calculate percentages.
    data["percentage"] = data["proportion"].apply(
        lambda value: (value * 100)
    )
    data = data.loc[
        :, data.columns.isin([
            "identifier", "name", "proportion", "percentage",
            "probability_log",
        ])
    ]
    data.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    # Return information.
    return data


def plot_chart_gene_heritability(
    genes_multimodal=None,
    genes_feature=None,
    threshold_proportion=None,
    threshold_probability_log=None,
    data=None,
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        genes_multimodal (list<str>): identifiers of multimodal genes
        genes_feature (list<str>): identifiers of feature genes
        threshold_proportion (float): threshold by proportion of phenotypic
            variance attributable to genotype
        threshold_probability_log (float): threshold by probability of
            heritability estimate
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
    figure = plot_scatter_label_emphasis_points(
        emphasis_keys=genes_multimodal,
        label_keys=genes_feature,
        column_key="identifier",
        column_label="name",
        line_abscissa=threshold_probability_log,
        line_ordinate=threshold_proportion,
        line_ordinate_origin=False,
        data=data,
        abscissa="probability_log",
        ordinate="proportion",
        title_abscissa="-1 * log10(p)",
        title_ordinate="V(g) / V(p)",
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure_png(
        path=path_file,
        figure=figure
    )

    pass

# BLARG
# TODO:
# read in the threshold information from heritability.py
# read in heritable genes along with multimodal and feature genes from a new function in integration.py


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
    source = read_source_gene_heritability(
        cohort="selection",
        dock=dock,
    )

    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_directory = os.path.join(path_plot, "heritability")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    path_file = os.path.join(
        path_directory, str("genes_heritability.png")
    )
    # Organize data for figure.
    data = organize_gene_heritability(
        data_genes_heritability=source["data_genes_heritability"],
    )
    # Plot chart and create figure.
    plot_chart_gene_heritability(
        genes_multimodal=source["genes_sets"]["candidacy"]["multimodal"],
        genes_feature=source["genes_sets"]["query"]["heritability"],
        threshold_proportion=source["threshold_proportion"],
        threshold_probability_log=source["threshold_probability_log"],
        data=data,
        path_file=path_file
    )
    pass


##########
# Collection of COVID-19 genes
# Status: working


def read_source_collection_comparisons_folds(
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

    # Read genes sets.
    genes_sets = integration.read_source_genes_sets_combination_integration(
        cohort="selection",
        dock=dock,
    )
    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_data_covid19 = os.path.join(
        dock, "collection", "covid19", "data_genes_comparisons_studies.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_covid19 = pandas.read_pickle(path_data_covid19)

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_covid19": data_covid19,
        "genes_sets": genes_sets,
    }


def plot_chart_collection_comparisons_folds(
    emphasis_keys=None,
    label_keys=None,
    column_key=None,
    column_label=None,
    data=None,
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        emphasis_keys (list<str>): keys for rows of points for emphasis
        label_keys (list<str>): keys for rows of points for labels
        column_key (str): name of column in data with keys
        column_label (str): name of column in data with labels
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
    figure = plot_scatter_label_emphasis_points(
        emphasis_keys=emphasis_keys,
        label_keys=label_keys,
        column_key=column_key,
        column_label=column_label,
        line_abscissa=None,
        line_ordinate=None,
        line_ordinate_origin=True,
        data=data,
        abscissa="log2_fold_direction",
        ordinate="comparisons",
        title_abscissa="log2 fold change",
        title_ordinate="count DE comparisons",
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure_png(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_collection_comparisons_folds(
    dock=None,
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
    source = read_source_collection_comparisons_folds(dock=dock)
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_directory = os.path.join(path_plot, "collection")
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    path_file = os.path.join(
        path_directory, "covid19_compilation.png"
    )

    # Copy data.
    data = source["data_covid19"].copy(deep=True)
    # Replicate records for positive and negative logarithmic fold changes of
    # genes in COVID-19 patients.
    data.reset_index(
        level=None,
        inplace=True
    )
    data_up = data.loc[
        :, data.columns.isin([
            "identifier", "name", "accumulations", "log2_fold_positive",
        ])
    ]
    data_down = data.loc[
        :, data.columns.isin([
            "identifier", "name", "depletions", "log2_fold_negative",
        ])
    ]
    data_up["direction"] = "up"
    data_down["direction"] = "down"
    data_up.rename(
        columns={
            "log2_fold_positive": "log2_fold_direction",
            "accumulations": "comparisons",
        },
        inplace=True,
    )
    data_down.rename(
        columns={
            "log2_fold_negative": "log2_fold_direction",
            "depletions": "comparisons",
        },
        inplace=True,
    )
    data_direction = data_up.append(
        data_down,
        ignore_index=True,
    )
    data_direction.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    # Plot chart and create figure.
    plot_chart_collection_comparisons_folds(
        emphasis_keys=source["genes_sets"]["candidacy"]["multimodal"],
        label_keys=source["genes_sets"]["query"]["feature"],
        column_key="identifier",
        column_label="name",
        data=data_direction,
        path_file=path_file
    )
    pass


##########
# Box plots for genes' signals between groups of persons
# Status: in progress


def read_source_genes_signals_persons_groups(
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
    path_comparisons_two = os.path.join(
        dock, "integration", "groups", "genes_comparisons_two_groups.pickle"
    )
    path_comparisons_four = os.path.join(
        dock, "integration", "groups", "genes_comparisons_four_groups.pickle"
    )
    # Read information from file.
    with open(path_comparisons_two, "rb") as file_source:
        comparisons_two = pickle.load(file_source)
    with open(path_comparisons_four, "rb") as file_source:
        comparisons_four = pickle.load(file_source)
    # Compile and return information.
    return {
        "comparisons_two": comparisons_two,
        "comparisons_four": comparisons_four,
    }


def plot_chart_genes_signals_persons_two_groups_comparison(
    comparison=None,
    path_directory=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        comparison (dict): information for chart
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Specify path to directory and file.
    path_file = os.path.join(
        path_directory, str(comparison["title"] + ".png")
    )
    # Organize group labels.
    label_1 = str(
        comparison["group_1_label"] +
        " (" + str(comparison["group_1_valids"]) + ")"
    )
    label_2 = str(
        comparison["group_2_label"] +
        " (" + str(comparison["group_2_valids"]) + ")"
    )
    labels_groups = [label_1, label_2]

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_boxes(
        arrays=[comparison["group_1_values"], comparison["group_2_values"]],
        labels_groups=labels_groups,
        label_vertical=str(comparison["gene_name"] + " pan-tissue signal"),
        label_horizontal="",
        label_top_center=str("p: " + numpy.format_float_scientific(
            comparison["probability"], precision=3
        )),
        label_top_left="",
        label_top_right="",
        orientation="portrait",
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure_png(
        path=path_file,
        figure=figure
    )

    pass


def plot_chart_genes_signals_persons_four_groups_comparison(
    comparison=None,
    path_directory=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        comparison (dict): information for chart
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Specify path to directory and file.
    path_file = os.path.join(
        path_directory, str(comparison["title"] + ".png")
    )
    # Organize group labels.
    label_1 = str(
        comparison["group_1_label"] +
        " (" + str(comparison["group_1_valids"]) + ")"
    )
    label_2 = str(
        comparison["group_2_label"] +
        " (" + str(comparison["group_2_valids"]) + ")"
    )
    label_3 = str(
        comparison["group_3_label"] +
        " (" + str(comparison["group_3_valids"]) + ")"
    )
    label_4 = str(
        comparison["group_4_label"] +
        " (" + str(comparison["group_4_valids"]) + ")"
    )
    labels_groups = [label_1, label_2, label_3, label_4]

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_boxes(
        arrays=[
            comparison["group_1_values"], comparison["group_2_values"],
            comparison["group_3_values"], comparison["group_4_values"],
        ],
        labels_groups=labels_groups,
        label_vertical=str(comparison["gene_name"] + " pan-tissue signal"),
        label_horizontal="",
        label_top_center="",
        label_top_left=str("p: " + numpy.format_float_scientific(
            comparison["probability_1_2"], precision=3
        )),
        label_top_right=str("p: " + numpy.format_float_scientific(
            comparison["probability_3_4"], precision=3
        )),
        orientation="landscape",
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure_png(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_genes_signals_persons_groups(
    dock=None,
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
    source = read_source_genes_signals_persons_groups(dock=dock)
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_integration = os.path.join(path_plot, "integration")
    path_directory_two = os.path.join(path_integration, "groups_two")
    utility.remove_directory(path=path_directory_two)
    utility.create_directories(path=path_directory_two)
    path_directory_four = os.path.join(path_integration, "groups_four")
    utility.remove_directory(path=path_directory_four)
    utility.create_directories(path=path_directory_four)
    # Create figures.
    for comparison in source["comparisons_two"]:
        plot_chart_genes_signals_persons_two_groups_comparison(
            comparison=comparison,
            path_directory=path_directory_two,
        )
        pass
    for comparison in source["comparisons_four"]:
        plot_chart_genes_signals_persons_four_groups_comparison(
            comparison=comparison,
            path_directory=path_directory_four,
        )
        pass
    pass



##########
# Overlap between sets in selection of genes by bimodality
# Status: working


def read_source_candidacy_gene_sets_overlap_measures(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_sets_genes_multimodal = os.path.join(
        dock, "candidacy", cohort, "distribution", "multimodal",
        "sets_multimodal.pickle"
    )
    path_sets_genes_unimodal = os.path.join(
        dock, "candidacy", cohort, "distribution", "unimodal",
        "sets_unimodal.pickle"
    )
    # Read information from file.
    with open(path_sets_genes_unimodal, "rb") as file_source:
        sets_genes_unimodal = pickle.load(file_source)
    with open(path_sets_genes_multimodal, "rb") as file_source:
        sets_genes_multimodal = pickle.load(file_source)

    # Compile and return information.
    return {
        "sets_genes_unimodal": sets_genes_unimodal,
        "sets_genes_multimodal": sets_genes_multimodal,
    }


def read_source_candidacy_gene_sets_overlap_multimodal_cohorts(
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
    path_selection = os.path.join(
        dock, "candidacy", "selection", "distribution", "multimodal",
        "genes.pickle"
    )
    path_respiration = os.path.join(
        dock, "candidacy", "respiration", "distribution", "multimodal",
        "genes.pickle"
    )
    path_ventilation = os.path.join(
        dock, "candidacy", "ventilation", "distribution", "multimodal",
        "genes.pickle"
    )
    # Read information from file.
    with open(path_selection, "rb") as file_source:
        genes_multimodal_selection = pickle.load(file_source)
    with open(path_respiration, "rb") as file_source:
        genes_multimodal_respiration = pickle.load(file_source)
    with open(path_ventilation, "rb") as file_source:
        genes_multimodal_ventilation = pickle.load(file_source)
    # Compile and return information.
    return {
        "selection": genes_multimodal_selection,
        "respiration": genes_multimodal_respiration,
        "ventilation": genes_multimodal_ventilation,
    }


def plot_chart_candidacy_gene_sets_overlap(
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
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )
    pass


def prepare_charts_candidacy_gene_sets_overlap_measures(
    cohort=None,
    dock=None
):
    """
    Plots charts.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_candidacy_gene_sets_overlap_measures(
        cohort=cohort,
        dock=dock,
    )
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_candidacy = os.path.join(path_plot, "candidacy")
    path_sets = os.path.join(path_candidacy, "sets")
    path_directory = os.path.join(path_sets, "measures")
    # Remove previous files to avoid version or batch confusion.
    #utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    path_multimodal = os.path.join(
        path_directory, str(cohort + "_multimodal_genes_measures.svg")
    )
    plot_chart_candidacy_gene_sets_overlap(
        sets=source["sets_genes_multimodal"],
        path=path_multimodal,
    )
    pass


def prepare_chart_candidacy_gene_sets_overlap_cohorts(
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
    source = read_source_candidacy_gene_sets_overlap_multimodal_cohorts(
        dock=dock,
    )
    sets_multimodal_genes_cohorts = dict()
    sets_multimodal_genes_cohorts["selection"] = source["selection"]
    sets_multimodal_genes_cohorts["respiration"] = source["respiration"]
    sets_multimodal_genes_cohorts["ventilation"] = source["ventilation"]
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_candidacy = os.path.join(path_plot, "candidacy")
    path_sets = os.path.join(path_candidacy, "sets")
    path_directory = os.path.join(path_sets, "cohorts")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    path_figure = os.path.join(
        path_directory, "multimodal_genes_cohorts.svg"
    )
    plot_chart_candidacy_gene_sets_overlap(
        sets=sets_multimodal_genes_cohorts,
        path=path_figure,
    )
    pass


def prepare_charts_candidacy_gene_sets_overlap(
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

    # Prepare charts for overlap of genes by multiple modality measures.
    cohorts = [
        "selection",
        "respiration",
        "ventilation",
    ]
    for cohort in cohorts:
        prepare_charts_candidacy_gene_sets_overlap_measures(
            cohort=cohort,
            dock=dock,
        )
    # Prepare chart for overlap of multimodal genes for each cohort.
    prepare_chart_candidacy_gene_sets_overlap_cohorts(
        dock=dock,
    )

    pass


##########
# Overlap between sets of genes by association to predictor variables in
# regression.
# Status: working


# TODO: obsolete for now...
def read_source_prediction_genes_set(
    cohort=None,
    distribution=None,
    variable=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        distribution (str): group of genes by their distribution of pan-tissue
            signals
        variable (str): name of predictor variable from regression that
            corresponds to set of genes
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_sets_genes = os.path.join(
        dock, "prediction", cohort, "genes",
        "sets_genes.pickle"
    )
    # Read information from file.
    with open(path_sets_genes, "rb") as file_source:
        sets_genes = pickle.load(file_source)
    # Compile and return information.
    return sets_genes[distribution][variable]


def define_parameters_prediction_genes_sets_overlap():
    """
    Defines parameters for plots of persons' properties.

    arguments:

    raises:

    returns:
        (dict): collection of parameters

    """

    # Use "source["sets_genes"]["multimodal"][parameter["set"]]" to iterate
    # through multiple sets of genes.
    parameters = list()
    if True:
        parameters.append(dict(
            title="sex_ventilation",
            cohort_one="selection",
            cohort_two="selection",
            cohort_three="selection",
            distribution_one="any",
            distribution_two="any",
            distribution_three="sex_ventilation",
            variable_one="sex_risk*ventilation_binary_scale",
            variable_two="sex_y_scale",
            variable_three="ventilation_binary_scale",
            label_one="sex*ventilation",
            label_two="sex",
            label_three="ventilation",
        ))
    if True:
        parameters.append(dict(
            title="age_ventilation",
            cohort_one="selection",
            cohort_two="selection",
            cohort_three="selection",
            distribution_one="any",
            distribution_two="any",
            distribution_three="age_ventilation",
            variable_one="age*ventilation_binary_scale",
            variable_two="age_scale",
            variable_three="ventilation_binary_scale",
            label_one="age*ventilation",
            label_two="age",
            label_three="ventilation",
        ))
    # Return information.
    return parameters


def plot_chart_prediction_genes_sets_overlap(
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
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure(
        path=path,
        figure=figure
    )
    pass


def prepare_chart_prediction_gene_sets_overlap(
    title=None,
    cohort_one=None,
    cohort_two=None,
    cohort_three=None,
    distribution_one=None,
    distribution_two=None,
    distribution_three=None,
    variable_one=None,
    variable_two=None,
    variable_three=None,
    label_one=None,
    label_two=None,
    label_three=None,
    path_directory=None,
    dock=None,
):
    """
    Plots charts.

    arguments:
        title (str): name for chart's file
        cohort_one (str): cohort of persons for set one
        cohort_two (str): cohort of persons for set two
        cohort_three (str): cohort of persons for set three
        distribution_one (str): group of genes by distribution for set one
        distribution_two (str): group of genes by distribution for set two
        distribution_three (str): group of genes by distribution for set three
        variable_one (str): name of predictor variable from regression that
            corresponds to set one
        variable_two (str): name of predictor variable from regression that
            corresponds to set two
        variable_three (str): name of predictor variable from regression that
            corresponds to set three
        label_one (str): label on chart for set one
        label_two (str): label on chart for set two
        label_three (str): label on chart for set three
        path_directory (str): path to directory
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # TODO: "distribution_three" is a temporary stand-in for "model"

    # Read lists of genes' identifiers from file.
    source = integration.read_source_genes_sets_prediction_interaction(
        cohort="selection",
        dock=dock,
    )
    genes_sets = dict()
    genes_sets[label_one] = source[distribution_three]["multimodal"][variable_one]
    genes_sets[label_two] = source["main"]["multimodal"][variable_two]
    genes_sets[label_three] = source["main"]["multimodal"][variable_three]
    # Define path to file.
    path_file = os.path.join(
        path_directory, str(title + ".svg")
    )
    # Plot figure.
    plot_chart_prediction_genes_sets_overlap(
        sets=genes_sets,
        path=path_file,
    )
    pass


def prepare_charts_prediction_gene_sets_overlap(
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

    # Specify directories and files.
    path_directory = os.path.join(
        dock, "plot", "prediction", "gene_sets"
    )
    utility.create_directories(path=path_directory)
    # Specify combinations of parameters for charts.
    parameters = define_parameters_prediction_genes_sets_overlap()
    for parameter in parameters:
        # Prepare charts for genes.
        prepare_chart_prediction_gene_sets_overlap(
            title=parameter["title"],
            cohort_one=parameter["cohort_one"],
            cohort_two=parameter["cohort_two"],
            cohort_three=parameter["cohort_three"],
            distribution_one=parameter["distribution_one"],
            distribution_two=parameter["distribution_two"],
            distribution_three=parameter["distribution_three"],
            variable_one=parameter["variable_one"],
            variable_two=parameter["variable_two"],
            variable_three=parameter["variable_three"],
            label_one=parameter["label_one"],
            label_two=parameter["label_two"],
            label_three=parameter["label_three"],
            path_directory=path_directory,
            dock=dock,
        )
        pass
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
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        sets (dict<list<str>>): values in sets
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
        "sets": sets_selection,
    }
    charts.append(record)
    record = {
        "title": "genes_unimodal",
        "sets": source["sets_genes_unimodal"],
    }
    charts.append(record)
    record = {
        "title": "genes_multimodal",
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
    path_file=None
):
    """
    Plots charts from the analysis process.

    arguments:
        sets (dict<list<str>>): values in sets
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
        "sets": source["sets_genes_measures"],
    }
    charts.append(record)
    record = {
        "title": "sets_genes_selection",
        "sets": source["sets_genes_selection"],
    }
    charts.append(record)
    record = {
        "title": "sets_genes_unimodal",
        "sets": source["sets_genes_unimodal"],
    }
    charts.append(record)
    record = {
        "title": "sets_genes_multimodal",
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


def read_source_genes_persons_signals(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", cohort,
        "data_persons_properties.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )
    # Read information from file.
    query_gene_sets = integration.read_source_annotation_query_genes_sets(
        dock=dock
    )
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    # Compile and return information.
    return {
        "query_gene_sets": query_gene_sets,
        "data_gene_annotation": data_gene_annotation,
        "data_signals_genes_persons": data_signals_genes_persons,
    }


def read_source_genes_persons_signals_gene(
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
    write_figure_png(
        path=path,
        figure=figure
    )

    pass


def prepare_chart_gene_persons_signals(
    gene=None,
    data_gene_annotation=None,
    data_signals_genes_persons=None,
    path_directory=None,
):
    """
    Prepares and plots data.

    arguments:
        gene (str): identifier of a gene
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        path_directory (str): path to directory

    raises:

    returns:

    """

    if gene in data_signals_genes_persons.columns.tolist():
        # Access gene's name.
        name = assembly.access_gene_name(
            identifier=gene,
            data_gene_annotation=data_gene_annotation,
        )
        # Access information.
        signals = data_signals_genes_persons[gene].dropna().to_list()
        # Create charts for the gene.
        path_figure = os.path.join(path_directory, str(name + ".png"))
        plot_chart_genes_persons_signals(
            gene=gene,
            name=name,
            values=signals,
            path=path_figure
        )
    pass


def prepare_charts_genes_persons_signals_cohort(
    cohort=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_genes_persons_signals(
        cohort=cohort,
        dock=dock,
    )
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_distribution = os.path.join(path_plot, "distribution")
    path_directory = os.path.join(path_distribution, "genes_persons", cohort)
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Plot charts for genes of interest.
    genes_interest = source["query_gene_sets"]["distribution"]
    for gene in genes_interest:
        prepare_chart_gene_persons_signals(
            gene=gene,
            data_gene_annotation=source["data_gene_annotation"],
            data_signals_genes_persons=source["data_signals_genes_persons"],
            path_directory=path_directory,
        )
        pass
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

    cohorts = [
        "selection",
        #"respiration",
        #"ventilation",
    ]
    for cohort in cohorts:
        prepare_charts_genes_persons_signals_cohort(
            cohort=cohort,
            dock=dock,
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
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", cohort,
        "data_persons_properties.pickle"
    )
    path_genes_distribution = os.path.join(
        dock, "distribution", cohort, "genes"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    genes_distribution = utility.extract_subdirectory_names(
        path=path_genes_distribution
    )
    query_gene_sets = integration.read_source_annotation_query_genes_sets(
        dock=dock
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "genes_distribution": genes_distribution,
        "query_gene_sets": query_gene_sets,
    }


def read_source_genes_signals_tissues_persons(
    gene=None,
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_gene_persons_signals = os.path.join(
        dock, "distribution", cohort, "genes", gene,
        "data_gene_persons_signals.pickle"
    )
    path_data_gene_signals_tissues_persons = os.path.join(
        dock, "distribution", cohort, "genes", gene,
        "data_gene_signals_tissues_persons_normal_standard.pickle"
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


def sort_gene_signals_tissues_persons_by_pantissue(
    data_gene_persons_signals=None,
    data_gene_signals_tissues_persons=None,
):
    """
    Organize information for chart.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons
        data_gene_signals_tissues_persons (object): Pandas data frame of a
            gene's signals across tissues and persons

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across tissues and
            persons after sort on pan-tissue signals

    """

    # Copy data.
    data_gene_persons_signals = data_gene_persons_signals.copy(deep=True)
    data_gene_signals_tissues_persons = data_gene_signals_tissues_persons.copy(
        deep=True
    )
    # Rename the aggregate, pantissue signals.
    data_gene_persons_signals.rename(
        columns={"value": "signal",},
        inplace=True,
    )
    # Introduce aggregate, pantissue signals to tissue-person matrix.
    # These aggregate signals will be useful to sort the data.
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
        labels=["signal"],
        axis="columns",
        inplace=True
    )
    # Return information
    return data_hybrid


def define_parameters_genes_signals_tissues_persons():
    """
    Defines parameters for plots of persons' properties.

    arguments:

    raises:

    returns:
        (dict): collection of parameters

    """

    # Use "source["sets_genes"]["multimodal"][parameter["set"]]" to iterate
    # through multiple sets of genes.
    parameters = list()
    parameters.append(dict(
        title="sex", label="sex",
        set="sex_y_scale",
        property="sex_y", type="binary",
    ))
    parameters.append(dict(
        title="ventilation_binary", label="ventilation",
        set="ventilation_binary_scale",
        property="ventilation_binary", type="binary",
    ))
    parameters.append(dict(
        title="leukocyte_binary", label="leukocyte",
        set="leukocyte_binary_scale",
        property="leukocyte_binary", type="binary",
    ))
    # Return information.
    return parameters


def plot_chart_genes_signals_tissues_persons_by_property(
    gene=None,
    gene_name=None,
    title=None,
    label=None,
    master=None,
    type_master=None,
    labels_categories_master=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a gene
        gene_name (str): name of gene
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        labels_categories_master (list<str>): labels for scale ticks of
            categorical master variable
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(gene_name + ".png")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_heatmap_asymmetric_master_main_top_bottom(
        title=gene_name,
        title_subordinate="",
        label_master=label,
        labels_categories_master=labels_categories_master,
        label_main="genes' signals across tissues and persons",
        type_master=type_master,
        type_main="continuous_divergent",
        data=data,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_genes_signals_tissues_persons_by_property_gene(
    title=None,
    label=None,
    master=None,
    type_master=None,
    gene=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    cohort=None,
    path_directory=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        gene (str): identifier of a gene
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        cohort (str): cohort of persons--selection, respiration, or ventilation
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
    source = read_source_genes_signals_tissues_persons(
        gene=gene,
        cohort=cohort,
        dock=dock
    )
    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    data_sort = sort_gene_signals_tissues_persons_by_pantissue(
        data_gene_persons_signals=source["data_gene_persons_signals"],
        data_gene_signals_tissues_persons=source[
            "data_gene_signals_tissues_persons"
        ],
    )
    # Organize master and main data.
    bin = organize_data_master_main(
        data_master=data_persons_properties,
        master=master,
        type_master=type_master,
        data_main=data_sort,
        type_main="countinuous_divergent",
        scale_unit_main=False,
        columns_main_scale_unit=list(),
        fill_missing=True,
        index="person",
        sequence="sort", # "none", "sort", or "cluster"
    )
    # Create charts for the gene.
    plot_chart_genes_signals_tissues_persons_by_property(
        gene=gene,
        gene_name=name,
        title=title,
        label=label,
        master=master,
        type_master=type_master,
        labels_categories_master=bin["labels_categories_master"],
        data=bin["data"],
        path_directory=path_directory,
    )
    pass


def prepare_charts_genes_signals_tissues_persons_by_property(
    title=None,
    label=None,
    master=None,
    type_master=None,
    genes_query=None,
    genes_distribution=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    cohort=None,
    path_parent=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        genes_query (list<str>): identifiers of genes for which to plot signals
        genes_distribution (list<str>): identifiers of genes with signals
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        cohort (str): cohort of persons--selection, respiration, or ventilation
        path_parent (str): path to directory
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Filter query genes to ensure that they have signals.
    genes_query_valid = utility.filter_common_elements(
        list_one=genes_query,
        list_two=genes_distribution,
    )
    # Specify directories and files.
    path_directory = os.path.join(path_parent, title)
    utility.create_directories(path=path_directory)
    # Iterate on genes.
    for gene in genes_query_valid:
        prepare_charts_genes_signals_tissues_persons_by_property_gene(
            title=title,
            label=label,
            master=master,
            type_master=type_master,
            gene=gene,
            data_gene_annotation=data_gene_annotation,
            data_persons_properties=data_persons_properties,
            cohort=cohort,
            path_directory=path_directory,
            dock=dock,
        )
        pass
    pass


def prepare_charts_genes_signals_tissues_persons_by_properties(
    genes_query=None,
    genes_distribution=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    cohort=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        genes_query (list<str>): identifiers of genes for which to plot signals
        genes_distribution (list<str>): identifiers of genes with signals
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Specify directories and files.
    path_directory = os.path.join(
        dock, "plot", "distribution", "tissues_persons_by_property", cohort,
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    # Specify combinations of parameters for charts.
    parameters = define_parameters_genes_signals_tissues_persons()
    # Iterate on parameters.
    for parameter in parameters:
        # Prepare charts for genes.
        prepare_charts_genes_signals_tissues_persons_by_property(
            title=parameter["title"],
            label=parameter["label"],
            master=parameter["property"],
            type_master=parameter["type"],
            genes_query=genes_query,
            genes_distribution=genes_distribution,
            data_gene_annotation=data_gene_annotation,
            data_persons_properties=data_persons_properties,
            cohort=cohort,
            path_parent=path_directory,
            dock=dock,
        )
        pass
    pass


def plot_chart_genes_signals_tissues_persons(
    gene=None,
    gene_name=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a gene
        gene_name (str): name of gene
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(gene_name + ".png")
    )
    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    label_scale = str(gene_name + "   signals across tissues and persons")

    # Create figure.
    figure = plot_heatmap_asymmetric(
        title="",
        label_scale=label_scale,
        type="continuous_divergent",
        label_rows=True,
        label_columns=True,
        data=data,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


# BLARG
# TODO: organize data in a new general and versatile function...
# TODO: drop rows columns with all missing, fill missing, scale columns, cluster


def prepare_charts_genes_signals_tissues_persons_gene(
    gene=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    cohort=None,
    path_directory=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        gene (str): identifier of a gene
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        cohort (str): cohort of persons--selection, respiration, or ventilation
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
    source = read_source_genes_signals_tissues_persons(
        gene=gene,
        cohort=cohort,
        dock=dock
    )
    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    data_sort = sort_gene_signals_tissues_persons_by_pantissue(
        data_gene_persons_signals=source["data_gene_persons_signals"],
        data_gene_signals_tissues_persons=source[
            "data_gene_signals_tissues_persons"
        ],
    )
    # Drop any rows or columns in main data with only missing values.
    data_sort.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    data_sort.dropna(
        axis="columns",
        how="all",
        inplace=True,
    )
    # Replace missing values in main data with zero.
    # Set any infinite values to missing.
    data_sort[data_sort == numpy.inf] = numpy.nan
    data_sort.fillna(
        value=0.0,
        #axis="columns",
        inplace=True,
    )
    # Cluster columns of data.
    data_cluster_columns = utility.cluster_data_columns(
        data=data_sort,
    )
    data_cluster_columns.reset_index(
        level=None,
        inplace=True
    )
    data_cluster_columns.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    data_cluster_rows = utility.cluster_data_rows(
        data=data_cluster_columns,
    )

    # Scale data?

    # TODO: write a function similar to "organize_data_master_main" but for the simpler scenario
    # --> yeah... probably worth it
    # enable cluster by columns flag
    # enable cluster by rows flag
    # enable scale adjustment for specific columns as in the master_main function

    # Create charts for the gene.
    plot_chart_genes_signals_tissues_persons(
        gene=gene,
        gene_name=name,
        data=data_cluster_rows,
        path_directory=path_directory,
    )
    pass


def prepare_charts_genes_signals_tissues_persons_simple(
    genes_query=None,
    genes_distribution=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    cohort=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        genes_query (list<str>): identifiers of genes for which to plot signals
        genes_distribution (list<str>): identifiers of genes with signals
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Specify directories and files.
    path_directory = os.path.join(
        dock, "plot", "distribution", "tissues_persons", cohort,
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    # Filter query genes to ensure that they have signals.
    genes_query_valid = utility.filter_common_elements(
        list_one=genes_query,
        list_two=genes_distribution,
    )
    # Iterate on genes.
    for gene in genes_query_valid:
        prepare_charts_genes_signals_tissues_persons_gene(
            gene=gene,
            data_gene_annotation=data_gene_annotation,
            data_persons_properties=data_persons_properties,
            cohort=cohort,
            path_directory=path_directory,
            dock=dock,
        )
        pass
    pass


def prepare_charts_genes_signals_tissues_persons_cohort(
    cohort=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_genes_signals_tissues_persons_initial(
        cohort=cohort,
        dock=dock,
    )
    # Prepare plots with persons sorted and clustered by persons' properties.
    prepare_charts_genes_signals_tissues_persons_by_properties(
        genes_query=source["query_gene_sets"]["distribution"],
        genes_distribution=source["genes_distribution"],
        data_gene_annotation=source["data_gene_annotation"],
        data_persons_properties=source["data_persons_properties"],
        cohort=cohort,
        dock=dock,
    )
    # Prepare plots with persons sorted by pan-tissue signal.
    prepare_charts_genes_signals_tissues_persons_simple(
        genes_query=source["query_gene_sets"]["distribution"],
        genes_distribution=source["genes_distribution"],
        data_gene_annotation=source["data_gene_annotation"],
        data_persons_properties=source["data_persons_properties"],
        cohort=cohort,
        dock=dock,
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

    cohorts = [
        "selection",
        #"respiration",
        #"ventilation",
    ]
    for cohort in cohorts:
        prepare_charts_genes_signals_tissues_persons_cohort(
            cohort=cohort,
            dock=dock,
        )
        pass
    pass


##########
# Genes' signals across persons
# Sort order of persons is by pan-tissue signal
# heatmaps
# Status: working


def read_source_genes_signals_persons(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Read genes sets.
    genes_sets = integration.read_source_genes_sets_combination_integration(
            cohort=cohort,
            dock=dock,
    )
    genes_correlation = integration.read_source_genes_sets_correlation(
            cohort=cohort,
            dock=dock,
    )
    path_genes_distribution = os.path.join(
        dock, "distribution", cohort, "genes"
    )
    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    genes_distribution = utility.extract_subdirectory_names(
        path=path_genes_distribution
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_signals_genes_persons": data_signals_genes_persons,
        "genes_sets": genes_sets,
        "genes_correlation": genes_correlation,
        "genes_distribution": genes_distribution,
    }


def plot_chart_genes_signals_persons(
    set_name=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        set_name (str): name of set of genes
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(set_name + ".png")
    )
    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    label_scale = str("genes' signals across tissues and persons")

    # Create figure.
    figure = plot_heatmap_asymmetric(
        title="",
        label_scale=label_scale,
        type="continuous_divergent",
        label_rows=True,
        label_columns=True,
        data=data,
        fonts=fonts,
        colors=colors,
    )
    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_chart_genes_signals_persons(
    set_name=None,
    genes=None,
    data_gene_annotation=None,
    data_signals_genes_persons=None,
    path_directory=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        set_name (str): name of set of genes
        genes (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        path_directory (str): path to directory


    raises:

    returns:

    """

    # Organize data.
    # Select and translate data for genes.
    data_signals = (
        integration.select_translate_gene_identifiers_data_columns(
            genes_query=genes,
            data_gene_annotation=data_gene_annotation,
            data_signals_genes_persons=data_signals_genes_persons,
    ))
    # Drop any rows or columns in main data with only missing values.
    data_signals.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    data_signals.dropna(
        axis="columns",
        how="all",
        inplace=True,
    )
    # Replace missing values in main data with zero.
    # Set any infinite values to missing.
    data_signals[data_signals == numpy.inf] = numpy.nan
    data_signals.fillna(
        value=0.0,
        #axis="columns",
        inplace=True,
    )
    # Cluster columns of data.
    data_cluster_columns = utility.cluster_data_columns(
        data=data_signals,
    )
    data_cluster_columns.reset_index(
        level=None,
        inplace=True
    )
    data_cluster_columns.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    data_cluster_rows = utility.cluster_data_rows(
        data=data_cluster_columns,
    )

    # Scale data?

    # TODO: write a function similar to "organize_data_master_main" but for the simpler scenario
    # --> yeah... probably worth it
    # enable cluster by columns flag
    # enable cluster by rows flag
    # enable scale adjustment for specific columns as in the master_main function

    # Create charts for the gene.
    plot_chart_genes_signals_persons(
        set_name=set_name,
        data=data_cluster_rows,
        path_directory=path_directory,
    )
    pass


def prepare_charts_genes_signals_persons(
    dock=None,
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
    source = read_source_genes_signals_persons(
        cohort="selection",
        dock=dock,
    )
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    path_integration = os.path.join(path_plot, "integration")
    path_directory = os.path.join(
        path_integration, "genes_across_persons"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)
    # Filter query genes to ensure that they have signals.
    genes_query = source["genes_correlation"]
    genes_query_valid = utility.filter_common_elements(
        list_one=genes_query,
        list_two=source["genes_distribution"],
    )
    # Prepare and plot data for sets of genes.
    prepare_chart_genes_signals_persons(
        set_name="correlation",
        genes=genes_query_valid,
        data_gene_annotation=source["data_gene_annotation"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
        path_directory=path_directory,
    )


    pass



##########
# Multiple genes' (rows) signals across persons (columns), with a chart
# for those persons' properties
# Sort order of persons is by property and my hierarchical clustering.
# heatmaps
# Status: working


def read_source_prediction_genes_signals_persons_properties(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Read genes sets.
    genes_sets = integration.read_source_genes_sets_combination_integration(
        cohort="selection",
        dock=dock,
    )
    genes_correlation = integration.read_source_genes_sets_correlation(
            cohort="selection",
            dock=dock,
    )

    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", cohort,
        "data_persons_properties.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "data_signals_genes_persons": data_signals_genes_persons,
        "genes_sets": genes_sets,
        "genes_correlation": genes_correlation,
    }


def define_parameters_prediction_genes_signals_persons_properties():
    """
    Defines parameters for plots of persons' properties.

    arguments:

    raises:

    returns:
        (dict): collection of parameters

    """

    # Use "source["sets_genes"]["multimodal"][parameter["set"]]" to iterate
    # through multiple sets of genes.
    parameters = list()
    if False:
        parameters.append(dict(
            name="smoke", set="smoke_scale",
            property="smoke", type="ordinal",
        ))
        parameters.append(dict(
            name="climate", set="climate_scale",
            property="climate", type="category",
        ))
        parameters.append(dict(
            title="leukocyte", label="leukocyte",
            set="leukocyte_binary_scale",
            property="leukocyte_binary", type="binary",
        ))
        parameters.append(dict(
            title="inflammation", label="inflammation",
            set="inflammation_binary_scale",
            property="inflammation_binary", type="binary",
        ))
        parameters.append(dict(
            title="cmv_ebv", label="CMV-EBV",
            set="mononucleosis_binary_scale",
            property="mononucleosis_binary", type="binary",
        ))
        parameters.append(dict(
            title="ventilation_duration", label="ventilation",
            set="ventilation_binary_scale",
            property="ventilation_duration_scale", type="continuous",
        ))
        parameters.append(dict(
            title="ventilation_grade", label="ventilation",
            set="ventilation_binary_scale",
            property="ventilation_duration_grade", type="ordinal",
        ))
    parameters.append(dict(
        title="sex", label="sex", set="sex_y_scale", property="sex_text",
        type="category",
    ))
    parameters.append(dict(
        title="age", label="age", set="age_scale",
        property="age_grade", type="ordinal",
    ))
    parameters.append(dict(
        title="race", label="race", set="race_white_scale", property="race_white",
        type="category",
    ))
    parameters.append(dict(
        title="ventilation_binary", label="ventilation",
        set="ventilation_binary_scale",
        property="ventilation_binary", type="binary",
    ))
    # Return information.
    return parameters


def plot_chart_prediction_genes_signals_persons_properties(
    title=None,
    label=None,
    master=None,
    type_master=None,
    labels_categories_master=None,
    data=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        labels_categories_master (list<str>): labels for scale ticks of
            categorical master variable
        data (object): Pandas data frame of a gene's aggregate, pantissue
            signals across tissues and persons
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define file name.
    path_figure = os.path.join(
        path_directory, str(title + ".png")
    )

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Create figure.
    figure = plot_heatmap_asymmetric_master_main_top_bottom(
        title=title,
        title_subordinate="",
        label_master=label,
        labels_categories_master=labels_categories_master,
        label_main="genes' pan-tissue signals across persons",
        type_master=type_master,
        type_main="continuous_divergent",
        data=data,
        fonts=fonts,
        colors=colors,
    )

    # Write figure.
    write_figure_png(
        path=path_figure,
        figure=figure
    )

    pass


def prepare_charts_query_genes_signals_persons_properties_variable(
    title=None,
    label=None,
    master=None,
    type_master=None,
    genes_query=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
    path_sort=None,
    path_cluster=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        title (str): title for chart
        label (str): label for property feature
        master (str): name of feature from persons' properties to use for
            master variable
        type_master (str): type of master variable
        genes_query (list<str>): identifiers of genes
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

    # Select and translate data for genes.
    data_signals_genes_persons = (
        integration.select_translate_gene_identifiers_data_columns(
            genes_query=genes_query,
            data_gene_annotation=data_gene_annotation,
            data_signals_genes_persons=data_signals_genes_persons,
    ))

    # Organize data.
    # Sort persons by their pantissue aggregate signals for the gene.
    # The same order is important to compare the heatmap to the histogram.
    bin = organize_data_master_main(
        data_master=data_persons_properties,
        master=master,
        type_master=type_master,
        data_main=data_signals_genes_persons,
        type_main="countinuous_divergent",
        scale_unit_main=False,
        columns_main_scale_unit=list(),
        fill_missing=True,
        index="person",
        sequence="sort", # "sort" or "cluster"
    )
    # Create charts for the gene.
    plot_chart_prediction_genes_signals_persons_properties(
        title=title,
        label=label,
        master=master,
        type_master=type_master,
        labels_categories_master=bin["labels_categories_master"],
        data=bin["data"],
        path_directory=path_sort,
    )

    if False:
        # Organize data.
        # Sort persons by their pantissue aggregate signals for the gene.
        # The same order is important to compare the heatmap to the histogram.
        bin_cluster = organize_prediction_genes_signals_persons_properties(
            property=property,
            type=type,
            sequence="cluster",
            genes_query=genes_query,
            data_gene_annotation=data_gene_annotation,
            data_persons_properties=data_persons_properties,
            data_signals_genes_persons=data_signals_genes_persons,
        )
        # Create charts for the gene.
        plot_chart_prediction_genes_signals_persons_properties(
            title=title,
            label=label,
            property=property,
            type=type,
            properties=bin_cluster["properties"],
            data=bin_cluster["data"],
            path_directory=path_cluster,
        )
    pass


def prepare_charts_query_genes_signals_persons_properties_genes(
    set_name=None,
    genes=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
    path_parent=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        set_name (str): name of set of genes
        genes (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        path_parent (str): path to directory

    raises:

    returns:

    """

    # Filter query genes to ensure that they have signals.
    genes_valid = utility.filter_common_elements(
        list_one=genes,
        list_two=data_signals_genes_persons.columns.tolist(),
    )
    # Specify directories and files.
    path_sort = os.path.join(
        path_parent, set_name, "sort_persons"
    )
    path_cluster = os.path.join(
        path_parent, set_name, "cluster_persons"
    )
    utility.create_directories(path=path_sort)
    utility.create_directories(path=path_cluster)
    # Specify combinations of parameters for charts.
    parameters = (
        define_parameters_prediction_genes_signals_persons_properties()
    )
    for parameter in parameters:
        # Report.
        #utility.print_terminal_partition(level=3)
        #print(str(len(source["sets_genes"]["multimodal"][parameter["set"]])))
        # Prepare charts for genes.
        prepare_charts_query_genes_signals_persons_properties_variable(
            title=parameter["title"],
            label=parameter["label"],
            master=parameter["property"],
            type_master=parameter["type"],
            genes_query=genes_valid,
            data_gene_annotation=data_gene_annotation,
            data_persons_properties=data_persons_properties,
            data_signals_genes_persons=data_signals_genes_persons,
            path_sort=path_sort,
            path_cluster=path_cluster,
        )
        pass
    pass


def prepare_charts_query_genes_signals_persons_properties_cohort(
    cohort=None,
    dock=None
):
    """
    Plots charts from the analysis process.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_prediction_genes_signals_persons_properties(
        cohort=cohort,
        dock=dock,
    )
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    path_prediction = os.path.join(path_plot, "prediction")
    path_persons_properties = os.path.join(
        path_prediction, "persons_properties"
    )
    path_cohort = os.path.join(
        path_persons_properties, cohort
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_cohort)
    utility.create_directories(path=path_cohort)
    # Select relevant persons.
    # Prepare and plot data for sets of genes.
    prepare_charts_query_genes_signals_persons_properties_genes(
        set_name="covid19_multimodal",
        genes=source["genes_sets"]["combination"]["covid19_multimodal"],
        data_gene_annotation=source["data_gene_annotation"],
        data_persons_properties=source["data_persons_properties"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
        path_parent=path_cohort,
    )
    prepare_charts_query_genes_signals_persons_properties_genes(
        set_name="covid19_multimodal_prediction",
        genes=(
            source["genes_sets"]["prediction"]["priority"]
        ),
        data_gene_annotation=source["data_gene_annotation"],
        data_persons_properties=source["data_persons_properties"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
        path_parent=path_cohort,
    )
    prepare_charts_query_genes_signals_persons_properties_genes(
        set_name="correlation",
        genes=source["genes_correlation"],
        data_gene_annotation=source["data_gene_annotation"],
        data_persons_properties=source["data_persons_properties"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
        path_parent=path_cohort,
    )
    pass


def prepare_charts_query_genes_signals_persons_properties(
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

    cohorts = [
        "selection",
        #"respiration",
        #"ventilation",
    ]
    for cohort in cohorts:
        prepare_charts_query_genes_signals_persons_properties_cohort(
            cohort=cohort,
            dock=dock,
        )
        pass
    pass


##########
# Principal components groups of persons
# Scatter plot
# Status: working


def read_source_persons_genes_components(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", cohort,
        "data_persons_properties.pickle"
    )
    path_data_persons_genes_components = os.path.join(
        dock, "stratification", cohort, "component",
        "data_persons_genes_components.pickle"
    )
    path_data_persons_genes_variances = os.path.join(
        dock, "stratification", cohort, "component",
        "data_persons_genes_variances.pickle"
    )
    # Read information from file.
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    data_persons_genes_components = pandas.read_pickle(
        path_data_persons_genes_components
    )
    data_persons_genes_variances = pandas.read_pickle(
        path_data_persons_genes_variances
    )
    # Compile and return information.
    return {
        "data_persons_properties": data_persons_properties,
        "data_persons_genes_components": data_persons_genes_components,
        "data_persons_genes_variances": data_persons_genes_variances,
    }


def organize_data_persons_genes_components(
    factor=None,
    data_persons_properties=None,
    data_components=None,
    data_variances=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        factor (str): name of categorical property to distinguish groups
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        data_components (object): Pandas data frame of principal components
            across instances
        data_variances (object): Pandas data frame of variances for principal
            components

    raises:

    returns:
        (dict): collection of information for charts

    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    data_components = data_components.copy(deep=True)
    data_variances = data_variances.copy(deep=True)
    # Organize data.
    data_variances.set_index(
        ["component"],
        append=False,
        drop=True,
        inplace=True,
    )
    data_components.rename_axis(
        index="person",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_persons_properties["factor"], labels_categories = pandas.factorize(
        data_persons_properties[factor],
        sort=True,
    )
    labels_categories = list(labels_categories)
    data_factor = data_persons_properties.loc[
        :, data_persons_properties.columns.isin(["factor", factor])
    ]
    # Join master and main data.
    data_hybrid = data_components.join(
        data_factor,
        how="left",
        on="person"
    )
    data_hybrid.drop(
        labels=[factor],
        axis="columns",
        inplace=True
    )
    # Organize labels for axes.
    variance_1 = round((data_variances.at[1, "variance"] * 100), 1)
    variance_2 = round((data_variances.at[2, "variance"] * 100), 1)
    variance_3 = round((data_variances.at[3, "variance"] * 100), 1)
    label_1 = ("Component 1 (" + str(variance_1) + "%)")
    label_2 = ("Component 2 (" + str(variance_2) + "%)")
    label_3 = ("Component 3 (" + str(variance_3) + "%)")
    # Compile information.
    bin = dict()
    bin["data_factor_components"] = data_hybrid
    bin["labels_factors"] = labels_categories
    bin["label_1"] = label_1
    bin["label_2"] = label_2
    bin["label_3"] = label_3
    # Return information.
    return bin


def plot_charts_persons_genes_components(
    data_factor_components=None,
    factor=None,
    label_factor=None,
    labels_factors=None,
    label_1=None,
    label_2=None,
    label_3=None,
    path_directory=None
):
    """
    Plots charts from the analysis process.

    arguments:
        data_factor_components (object): Pandas data frame of factor and
            components across observations
        factor (str): name of categorical property to distinguish groups
        label_factor (str): name of factor
        labels_factors (list<str>): names of categorical factor
        label_1 (str): label depicting variance of component 1
        label_2 (str): label depicting variance of component 2
        label_3 (str): label depicting variance of component 3
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Define path to file.
    file = str("components_1_2.png")
    path_file = os.path.join(path_directory, file)
    # Create figure.
    figure = plot_scatter_factor_groups(
        data=data_factor_components,
        abscissa="component_1", # variable for horizontal (x) axis
        ordinate="component_2", # variable for vertical (y) axis
        label_horizontal=label_1,
        label_vertical=label_2,
        factor=factor,
        label_factor=label_factor,
        labels_factors=labels_factors,
        fonts=fonts,
        colors=colors,
        point_size=7.5,
        plot_factor_labels=True,
        legend=True,
    )
    # Write figure.
    write_figure_png(
        path=path_file,
        figure=figure
    )

    # Define path to file.
    file = str("components_1_3.png")
    path_file = os.path.join(path_directory, file)
    # Create figure.
    figure = plot_scatter_factor_groups(
        data=data_factor_components,
        abscissa="component_1", # variable for horizontal (x) axis
        ordinate="component_3", # variable for vertical (y) axis
        label_horizontal=label_1,
        label_vertical=label_3,
        factor=factor,
        label_factor=label_factor,
        labels_factors=labels_factors,
        fonts=fonts,
        colors=colors,
        point_size=7.5,
        plot_factor_labels=True,
        legend=True,
    )
    # Write figure.
    write_figure_png(
        path=path_file,
        figure=figure
    )

    # Define path to file.
    file = str("components_2_3.png")
    path_file = os.path.join(path_directory, file)
    # Create figure.
    figure = plot_scatter_factor_groups(
        data=data_factor_components,
        abscissa="component_2", # variable for horizontal (x) axis
        ordinate="component_3", # variable for vertical (y) axis
        label_horizontal=label_2,
        label_vertical=label_3,
        factor=factor,
        label_factor=label_factor,
        labels_factors=labels_factors,
        fonts=fonts,
        colors=colors,
        point_size=7.5,
        plot_factor_labels=True,
        legend=True,
    )
    # Write figure.
    write_figure_png(
        path=path_file,
        figure=figure
    )

    pass


def prepare_charts_persons_genes_components_cohort(
    cohort=None,
    dock=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_persons_genes_components(
        cohort=cohort,
        dock=dock,
    )
    # Specify group variables.
    if cohort == "selection":
        factor = "ventilation"
    elif cohort == "respiration":
        factor = "hardiness" # hardiness
    elif cohort == "ventilation":
        #factor = "ventilation_duration_grade" # "sex_text", "ventilation", "age_grade"
        factor = "sex_text"
        #factor = "age_grade"
        pass
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_integration = os.path.join(path_plot, "integration")
    path_directory = os.path.join(
        path_integration, "persons_components", cohort, factor
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Organize information for charts.
    bin = organize_data_persons_genes_components(
        factor=factor,
        data_persons_properties=source["data_persons_properties"],
        data_components=source["data_persons_genes_components"],
        data_variances=source["data_persons_genes_variances"],
    )
    # Plot.
    plot_charts_persons_genes_components(
        data_factor_components=bin["data_factor_components"],
        factor="factor",
        label_factor=factor,
        labels_factors=bin["labels_factors"],
        label_1=bin["label_1"],
        label_2=bin["label_2"],
        label_3=bin["label_3"],
        path_directory=path_directory,
    )
    pass


def prepare_charts_persons_genes_components(
    dock=None
):
    """prepare_charts_persons_genes_components
    Plots charts from the analysis process.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    prepare_charts_persons_genes_components_cohort(
        cohort="selection",
        dock=dock,
    )
    prepare_charts_persons_genes_components_cohort(
        cohort="respiration",
        dock=dock,
    )
    prepare_charts_persons_genes_components_cohort(
        cohort="ventilation",
        dock=dock,
    )
    pass



##########
# Correlations in signals between pairs of genes
# Heatmap
# Status: working


def read_source_signals_genes_correlations(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_data_genes_correlations = os.path.join(
        dock, "integration", cohort, "correlation",
        "data_genes_correlations.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_genes_correlations = pandas.read_pickle(
        path_data_genes_correlations
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_genes_correlations": data_genes_correlations,
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
    # Return information.
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
    write_figure_png(
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
    source = read_source_signals_genes_correlations(
        cohort="selection",
        dock=dock,
    )
    # Organize data.
    data = organize_signals_genes_correlations(
        data_gene_annotation=source["data_gene_annotation"],
        data_genes_correlations=source["data_genes_correlations"],
    )
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    utility.create_directory(path_plot)
    path_integration = os.path.join(path_plot, "integration")
    path_directory = os.path.join(
        path_integration, "correlations_genes"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=path_directory)
    utility.create_directories(path=path_directory)

    # Define file name.
    path_file = os.path.join(
        path_directory, str("genes_correlations.png")
    )
    # Create chart.
    plot_chart_signals_genes_correlations(
        data=data,
        path_file=path_file,
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


##########
# Major gene ontology functional categories of multimodal genes
# Highlight genes of interest from regression
# Status: in progress


# TODO: read in all functional sets from integration procedure...

# ventilation_binary_scale
# mononucleosis_binary_scale
# leukocyte_binary_scale
# inflammation_binary_scale
def read_source_multimodal_genes_ontology_sets(dock=None):
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
    path_data = os.path.join(
        dock, "integration", "selection", "set", "cardinality",
        "ventilation_binary_scale.pickle"
    )
    # Read information from file.
    data = pandas.read_pickle(
        path_data
    )
    # Compile and return information.
    return {
        "data": data,
    }


def prepare_charts_multimodal_genes_ontology_sets(
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
    source = read_source_multimodal_genes_ontology_sets(dock=dock)

    # Organize data.
    data = source["data"].copy(deep=True)
    data.drop(
        labels=["total"],
        axis="columns",
        inplace=True
    )
    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Specify directories and files.
    path_plot = os.path.join(dock, "plot")
    path_directory = os.path.join(path_plot, "function")
    utility.create_directory(path_directory)

    # Create figures.
    figure = plot_bar_stack(
        data=data,
        label_vertical="genes in each category",
        label_horizontal="biological process parent sets",
        fonts=fonts,
        colors=colors,
        color_count=2,
        rotation="horizontal",
        legend=True,
    )
    # Specify directories and files.
    file = ("test.png")
    path_file = os.path.join(path_directory, file)
    # Write figure.
    write_figure_png(
        path=path_file,
        figure=figure
    )

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

# Tissue        print(group)



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
    # This chart depicts groups of individual variables that relate, such as
    # all variables for inflammation.
    #prepare_charts_persons_health_collection_variables(dock=dock)

    # Plot charts, heatmaps, for adjacent summaries of properties across
    # persons.
    #prepare_charts_persons_properties_adjacency(dock=dock)

    # Plot chart, histogram, for distribution of values of ventilation
    # duration.
    #prepare_chart_ventilation_duration_distribution(dock=dock)

    ##########
    ##########
    ##########
    # Selection procedure

    # TODO: what is this? probably obsolete...
    # Plot charts for sex, age, and hardiness of persons in GTEx cohort.
    #prepare_charts_persons_sex_age_hardiness(dock=dock)

    # Plot charts of overlap between sets of persons by clinical categories.
    #prepare_charts_sets_selection_persons_properties(dock=dock)


    #TODO: need to save the necessary file in selection procedure for this plot to work again.
    # Plot chart for counts of major tissue type per person, a histogram.
    #prepare_chart_tissues_per_person(dock=dock)


    if False:

        # Plot chart for counts of persons per major tissue type.
        prepare_chart_persons_per_tissue(dock=dock)

        # Plot charts of overlap between sets in selection of genes and samples.
        prepare_charts_sets_selection_genes_samples(dock=dock)

        # Plot charts for selection of principal components on categorical
        # variables for regression.
        prepare_charts_selection_dimension_principal_components(dock=dock)

    ##########
    ##########
    ##########
    # Collection procedure

    # Plot charts to summarize positive and negative logarithmic fold changes
    # of genes in multiple studies on COVID-19 patients.
    if False:
        prepare_charts_collection_comparisons_folds(dock=dock)

    ##########
    ##########
    ##########
    # Distribution procedure

    # Plot charts of distributions of genes' pan-tissue aggregate signals
    # across persons.
    if False:
        prepare_charts_genes_persons_signals(dock=dock)

    # Plot charts, heatmaps, for each gene's signals across persons (columns)
    # and tissues (rows).
    # Include a chart to depict persons' properties.
    # Sort order of persons is by pan-tissue signal.
    if False:
        prepare_charts_genes_signals_tissues_persons(dock=dock)

    if False:

        ##########
        ##########
        ##########
        # Candidacy procedure

        # Plot charts of overlap between sets in selection of genes by bimodality.
        prepare_charts_candidacy_gene_sets_overlap(dock=dock)

        # Plot charts of distributions of each modality measure's scores across
        # genes.
        prepare_charts_modality_measure_genes_distribution(dock=dock)


    if False:

        # Plot charts for correlations in signals between pairs of genes.
        # Charts will be scatter plots.
        # Each point will represent a person with pantissue signals from each gene.
        prepare_charts_signals_persons_gene_pairs(dock=dock)

    ##########
    ##########
    ##########
    # Heritability procedure

    # Plot charts for heritability of genes' pantissue signals.
    # Charts will be scatter plots.
    prepare_charts_gene_heritability(dock=dock)

    if False:
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

    if False:

        # Plot charts, scatter plots, for residuals from regressions on each gene's
        # pan-tissue signals across persons.
        #prepare_charts_genes_regression_residuals(dock=dock)

        # Plot charts, set overlap charts, for genes in multiple sets from
        # regression and association to predictor variables.
        prepare_charts_prediction_gene_sets_overlap(dock=dock)



    ##########
    ##########
    ##########
    # Function procedure

    # Major gene ontology functional categories of multimodal genes
    # Highlight genes of interest from regression
    #prepare_charts_multimodal_genes_ontology_sets(dock=dock)

    ##########
    ##########
    ##########
    # Integration procedure

    # Plot charts, box plots, for comparisons of genes' pan-tissue signals
    # between groups of persons.
    if False:
        prepare_charts_genes_signals_persons_groups(dock=dock)

    # Plot charts, heatmaps for multiple genes' pan-tissue signals across
    # persons.
    if True:
        prepare_charts_genes_signals_persons(dock=dock)

    # Plot charts, heatmaps, for multiple genes' pan-tissue signals across
    # persons along with those persons' properties.
    # Genes will be across rows, and persons will be across columns.
    # Sort order across rows depends on hierarchical clustering.
    # In some charts, sort order across columns depends on persons' properties
    # (sex, age, body mass index, hardiness).
    # In other charts, sort order across columns depends on hierarchical
    # clustering.
    if True:
        prepare_charts_query_genes_signals_persons_properties(dock=dock)

    # Plot charts for correlations between pairs of genes of interest.
    # Chart is adjacency matrix heatmap.
    # This function reads a single data table from the integration procedure.
    if True:
        prepare_charts_signals_genes_correlations(dock=dock)


    if False:

        # Plot charts, scatter plots, for components by genes' pan-tissue
        # signals across groups of persons.
        prepare_charts_persons_genes_components(dock=dock)

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
