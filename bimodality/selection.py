"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os
import math
import statistics
import pickle
import itertools

# Relevant

import numpy
import pandas

# Custom

import assembly
import measurement
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source(dock=None):
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
    path_gene_annotation_gtex = os.path.join(
        path_assembly, "data_gene_annotation_gtex.pickle"
    )
    path_gene_annotation_gencode = os.path.join(
        path_assembly, "data_gene_annotation_gencode.pickle"
    )
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.feather"
    )

    # Read information from file.
    data_gene_annotation_gtex = pandas.read_pickle(
        path_gene_annotation_gtex
    )
    data_gene_annotation_gencode = pandas.read_pickle(
        path_gene_annotation_gencode
    )
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_signal = pandas.read_feather(
        path=path_gene_signal,
    )
    # Compile and return information.
    return {
        "data_gene_annotation_gtex": data_gene_annotation_gtex,
        "data_gene_annotation_gencode": data_gene_annotation_gencode,
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_signal": data_gene_signal,
    }


def organize_data_axes_indices(data=None):
    """
    Organizes data with names and indices.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples.

    """

    # The Pandas function to rename axis copies data deeply by default and
    # requires a lot of memory.
    utility.print_terminal_partition(level=2)

    # Organize data.

    print(data.iloc[0:10, 0:10])
    print("Organize data with names and indices.")
    data.set_index(
        "gene",
        drop=True,
        inplace=True,
    )
    data.rename_axis(
        index="gene",
        axis="index",
        copy=False,
        inplace=True,
    )
    data.rename_axis(
        columns="sample",
        axis="columns",
        copy=False,
        inplace=True
    )
    print(data.iloc[0:10, 0:10])
    return data


# Summary.


def summarize_samples_genes(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Summarize selection of samples and genes.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_signal (object): Pandas data frame of genes' counts or
            signals across samples

    raises:

    returns:

    """

    # Copy data.
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    data_gene_signal = data_gene_signal.copy(deep=True)

    utility.print_terminal_partition(level=1)
    print("Summary of selection of samples and genes.")
    # Counts of samples and genes.
    print("Genes' signals...")
    print("Count of samples: " + str(data_gene_signal.shape[1]))
    print("Count of genes: " + str(data_gene_signal.shape[0]))

    # Major tissue categories.
    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_signal.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_sample=data_transposition,
    )

    groups = data_factor.groupby(level=["person"])
    print("Count of persons: " + str(len(groups)))

    groups = data_factor.groupby(level=["tissue_minor"])
    print("Count of minor tissues: " + str(len(groups)))

    groups = data_factor.groupby(level=["tissue_major"])
    print("Count of major tissues: " + str(len(groups)))

    print("Count of persons with samples for each major tissue.")
    records = []
    for name, group in groups:
        data = group.reset_index(
            level=None,
            inplace=False
        )
        data.drop_duplicates(
            subset=["person", "tissue_major"],
            keep="first",
            inplace=True,
        )
        # Compile information.
        record = dict()
        record["tissue_major"] = name
        record["persons"] = data.shape[0]
        records.append(record)

    data_summary = utility.convert_records_to_dataframe(
        records=records
    )
    print(data_summary)

    data_factor.reset_index(
        level=["tissue_major"],
        inplace=True
    )
    data_adipose = data_factor.loc[data_factor["tissue_major"] == "adipose", :]
    print(data_adipose)

    pass


##########
# Select genes and samples.

# Genes.


def select_gene_annotation(
    data_gene_annotation=None,
):
    """
    Selects genes that encode proteins.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' annotations

    """

    utility.print_terminal_partition(level=2)
    print(
        "Select annotations for genes that encode proteins and are on " +
        "autosomes."
    )

    # Select entries for genes.
    print(data_gene_annotation.shape)
    data_gene_annotation = (
        data_gene_annotation.loc[data_gene_annotation["feature"] == "gene", :]
    )
    data_gene_annotation.drop(
        labels="feature",
        axis="columns",
        inplace=True
    )
    print(data_gene_annotation.shape)

    # Select entries for genes on autosomes (non-sex chromosomes).
    # There are 2422 genes on the X chromosome, of which 850 encode proteins.
    # There are 567 genes on the Y chromosome, of which 66 encode proteins.
    # There are 37 genes on the mitochondrial chromosome, of which 13 encode
    # proteins.
    print(
        "count of genes in reference annotations: " +
        str(data_gene_annotation.shape[0])
    )
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["seqname"] != "chrX", :
        ]
    )
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["seqname"] != "chrY", :
        ]
    )
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["seqname"] != "chrM", :
        ]
    )
    print(
        "count of reference genes on nuclear (not mito) autosomes: " +
        str(data_gene_annotation.shape[0])
    )
    # Select entries for genes that encode proteins.
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["gene_type"] == "protein_coding", :
        ]
    )
    print(
        "count of reference genes that encode proteins: " +
        str(data_gene_annotation.shape[0])
    )

    print(data_gene_annotation.iloc[0:10, 0:7])

    return data_gene_annotation


def select_genes(
    data_gene_annotation=None,
    data_gene_signal=None
):
    """
    Selects genes that encode proteins.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and patients

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

    # Copy data.
    data_gene_annotation = data_gene_annotation.copy(deep=True)
    data_gene_signal = data_gene_signal.copy(deep=True)

    # Select protein-coding genes from signal data.
    # Signal matrix has genes on the index dimension.
    utility.print_terminal_partition(level=2)
    print(
        "Selection of signal genes that encode proteins."
    )
    # Describe original count of genes.
    print("signal total genes: " + str(data_gene_signal.shape[0]))
    genes_signal = data_gene_signal.index.to_list()
    print("signal total genes: " + str(len(genes_signal)))

    # Extract identifiers of protein-coding genes.
    genes_protein = data_gene_annotation.index.to_list()
    print(
        "count of reference genes of type 'protein_coding': " +
        str(len(genes_protein))
    )

    # Filter gene signals.
    # Remove indices.
    #data_gene_signal.reset_index(
    #    level=None,
    #    inplace=True
    #)
    #data_gene_signal_protein = data_gene_signal.loc[
    #    data_gene_signal["gene"].isin(genes_protein), :
    #]
    data_gene_signal_protein = data_gene_signal.loc[
        genes_protein, :
    ]
    genes = data_gene_signal_protein.index.to_list()

    utility.print_terminal_partition(level=2)
    print(
        "signal genes that encode proteins and on autosomes: " +
        str(data_gene_signal_protein.shape[0])
    )
    print(
        "signal genes that encode proteins and on autosomes: " +
        str(len(genes))
    )
    print(
        "unique signal genes that encode proteins and on autosomes: " +
        str(len(utility.collect_unique_elements(elements_original=genes)))
    )

    utility.print_terminal_partition(level=2)

    # Return information.
    return data_gene_signal_protein


# Samples.


def exclude_samples_by_tissues(
    tissues_minor=None,
    tissues_major=None,
    data_samples_tissues_persons=None,
):
    """
    Exclude samples for specific major and minor tissues.

    arguments:
        tissues_minor (list<str>): names of minor tissues to exclude
        tissues_major (list<str>): names of major tissues to exclude
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_minor = data_samples.loc[
        ~data_samples["tissue_minor"].isin(tissues_minor),
        :
    ]
    data_major = data_minor.loc[
        ~data_minor["tissue_major"].isin(tissues_major),
        :
    ]
    # Return information.
    return data_major


def include_samples_by_tissues(
    tissues_major=None,
    data_samples_tissues_persons=None,
):
    """
    Include samples for specific major and minor tissues.

    arguments:
        tissues_major (list<str>): names of major tissues to exclude
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_major = data_samples.loc[
        data_samples["tissue_major"].isin(tissues_major), :
    ]
    # Return information.
    return data_major


# Signal.


def filter_rows_columns_by_threshold_proportion(
    data=None,
    dimension=None,
    threshold=None,
    proportion=None,
):
    """
    Filters either rows or columns.

    Persistence of a row or column requires at least a specific count of values
    beyond a specific threshold.

    Filter rows by consideration of values across columns in each row.
    Filter columns by consideration of values across rows in each column.

    arguments:
        data (object): Pandas data frame of genes' signals across samples
        dimension (str): dimension to filter, either "row" or "column"
        threshold (float): minimal signal
        proportion (float): minimal proportion of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Determine count from proportion.
    if dimension == "row":
        # Filter rows by consideration of columns for each row.
        columns = data.shape[1]
        count = round(proportion * columns)
    elif dimension == "column":
        # Filter columns by consideration of rows for each columns.
        rows = data.shape[0]
        count = round(proportion * rows)

    # Determine whether values exceed threshold.
    data_threshold = (data >= threshold)
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    # This aggregation operation produces a series.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )

    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count, : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count]
    return data_pass


def select_samples_genes_signals(
    samples=None,
    data_gene_signal=None
):
    """
    Selects genes' signals by samples.

    Data format should have genes across rows and samples across columns.

    arguments:
        samples (list<str): identifiers of samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    data_gene_signal = data_gene_signal.copy(deep=True)
    utility.print_terminal_partition(level=2)
    print("shape of original data...")
    print(data_gene_signal.shape)
    # Select samples on the basis of tissues.
    data_gene_signal_sample = data_gene_signal.loc[
        :, data_gene_signal.columns.isin(samples)
    ]
    print("shape of data after selection of samples...")
    print(data_gene_signal_sample.shape)
    # Return information.
    return data_gene_signal_sample


def select_samples_genes_signals_coverage(
    threshold=None,
    proportion_gene=None,
    proportion_sample=None,
    data_gene_signal=None
):
    """
    Selects samples and genes by signals beyond threshold.

    Data format should have genes across rows and samples across columns.

    arguments:
        threshold (float): minimal gene's signal for a sample
        proportion_gene (float): minimal proportion of samples that must have
            signals beyond threshold for a gene to pass
        proportion_sample (float): minimal proportion of genes that must have
            signals beyond threshold for a sample to pass
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    # Filter genes by their signals across samples.
    # Filter to keep only genes with signals beyond threshold in proportion of
    # samples.
    data_row = filter_rows_columns_by_threshold_proportion(
        data=data_gene_signal,
        dimension="row",
        threshold=threshold,
        proportion=proportion_gene
    )
    print("shape of data after signal filter on genes...")
    print(data_row.shape)
    # Filter samples by their signals across genes.
    # Filter to keep only samples with signals beyond threshold in proportion
    # of genes.
    data_column = filter_rows_columns_by_threshold_proportion(
        data=data_row,
        dimension="column",
        threshold=threshold,
        proportion=proportion_sample
    )
    print("shape of data after signal filter on samples...")
    print(data_column.shape)

    # Return information.
    return data_column


# Filter information about samples, tissues, and persons.


def select_samples_tissues_persons(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Selects samples, tissues, persons, and their properties from the selection
    of samples on the basis of genes' signals.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples
    """

    # Extract selection samples.
    samples = utility.collect_unique_elements(
        elements_original=data_gene_signal.columns.tolist()
    )
    # Select information about samples.
    data_samples_copy = data_samples_tissues_persons.copy(deep=True)
    data_samples_selection = data_samples_copy.loc[
        data_samples_copy.index.isin(samples), :
    ]

    # Return information.
    return data_samples_selection


##########
# Organize information for further analyses.


# Extract information about samples, tissues, and persons.


def extract_persons_tissues_samples(
    data_samples_tissues_persons=None,
):
    """
    Extracts information about persons, tissues, and samples.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (dict<list<str>>): identifiers of persons, tissues, and samples
    """

    # Extract samples.
    samples = utility.collect_unique_elements(
        elements_original=data_samples_tissues_persons.index.to_list()
    )
    # Extract selection persons.
    persons = utility.collect_unique_elements(
        elements_original=data_samples_tissues_persons["person"].to_list()
    )
    # Extract selection tissues.
    tissues = utility.collect_unique_elements(
        elements_original=(
            data_samples_tissues_persons["tissue_major"].to_list()
        )
    )

    # Summary.
    utility.print_terminal_partition(level=2)
    print("count of unique samples: " + str(len(samples)))
    print("count of unique persons: " + str(len(persons)))
    print("count of unique tissues: " + str(len(tissues)))

    # Compile information.
    information = {
        "persons": persons,
        "tissues": tissues,
        "samples": samples,
    }

    # Return information.
    return information


# Count tissues per person and persons per tissue.


def count_factor_group_elements(
    factor=None,
    name_factor=None,
    name_elements=None,
    data=None,
):
    """
    Counts elements in groups by factor.

    arguments:
        factor (str): name of factor column in data
        name_factor (str): name for factor column in count data
        name_elements (str): name for element column in count data
        data (object): Pandas data frame of elements in groups by factors

    raises:

    returns:
        (object): Pandas data frame of counts of elements in each factor group
    """

    data_copy = data.copy(deep=True)

    data_copy.set_index(
        [factor],
        append=False,
        drop=True,
        inplace=True
    )
    groups = data_copy.groupby(level=[factor])

    records = []
    for name, group in groups:
        data_group = group.reset_index(
            level=None,
            inplace=False
        )
        data_group.drop_duplicates(
            subset=None,
            keep="first",
            inplace=True,
        )
        # Compile information.
        record = dict()
        record[name_factor] = name
        record[name_elements] = data_group.shape[0]
        records.append(record)

    data_count = utility.convert_records_to_dataframe(
        records=records
    )

    return data_count


def count_tissues_persons_groups(
    data_samples_tissues_persons=None,
):
    """
    Counts tissues per person and persons per tissue.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (dict): information about persons, tissues, and samples
    """

    utility.print_terminal_partition(level=1)

    # Copy data.
    data_samples_copy = data_samples_tissues_persons.copy(deep=True)
    # Organize data.
    data_samples_copy.reset_index(
        level=None,
        inplace=True
    )
    data_samples = data_samples_copy.loc[
        :, data_samples_copy.columns.isin([
            "sample", "person", "tissue_major"
        ])
    ]
    data_samples.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_samples.reindex()
    # At this point, data include multiple samples for minor tissues.
    data_samples.drop(
        labels=["sample"],
        axis="columns",
        inplace=True
    )
    data_samples.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    # At this point, data represent unique pairs of persons and major tissues
    # for which they have samples.

    # Count tissues per person and persons per tissue.
    data_tissues_per_person = count_factor_group_elements(
        factor="person",
        name_factor="person",
        name_elements="tissues",
        data=data_samples,
    )
    data_persons_per_tissue = count_factor_group_elements(
        factor="tissue_major",
        name_factor="tissue",
        name_elements="persons",
        data=data_samples,
    )

    # Compile information.
    information = {
        "data_tissues_per_person": data_tissues_per_person,
        "data_persons_per_tissue": data_persons_per_tissue,
        "tissues_per_person": data_tissues_per_person["tissues"].to_list(),
        "persons_per_tissue": data_persons_per_tissue["persons"].to_list(),
    }

    # Return information.
    return information


# Extract persons' properties.


def collect_person_unique_sample_values(
    values=None,
    delimiter=None,
):
    """
    Collects unique values from a person's samples.

    arguments:
        values (list<str>): values from a person's samples
        delimiter (str): character for separation of values

    raises:

    returns:
        (str): values

    """

    #values_valid = list(filter(lambda value: not math.isnan(value), values))
    #values_type = list(map(lambda value: str(value), values_valid))
    values_combination = delimiter.join(values)
    values_split = values_combination.split(delimiter)
    values_strip = list(map(lambda value: value.strip(), values_split))
    values_unique = utility.collect_unique_elements(
        elements_original=values_strip
    )
    values_string = delimiter.join(values_unique)
    return values_string


def extract_persons_properties(
    data_samples_tissues_persons=None,
):
    """
    Extracts information about persons.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (object): Pandas data frame of information about persons

    """

    # Organize data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_samples.rename_axis(
        "",
        axis="columns",
        inplace=True,
    )
    data_samples.reset_index(
        level=None,
        inplace=True
    )
    data_samples.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )

    # Split data by person.
    groups = data_samples.groupby(
        level=["person"],
    )
    data_collection = pandas.DataFrame()
    for name, data_original in groups:
        data_original.reset_index(
            level=None,
            inplace=True
        )
        data_novel = data_original.copy(deep=True)
        data_novel.drop(
            labels=[
                "sample",
                "tissue_major",
                "tissue_minor",
                "facilities",
                "batch_isolation",
                "batches_analysis",
            ],
            axis="columns",
            inplace=True
        )
        data_novel.drop_duplicates(
            subset=None,
            keep="first",
            inplace=True,
        )
        data_novel["facilities"] = collect_person_unique_sample_values(
            values=data_original["facilities"].dropna().to_list(),
            delimiter=",",
        )
        data_novel["batches_isolation"] = collect_person_unique_sample_values(
            values=data_original["batch_isolation"].dropna().to_list(),
            delimiter=",",
        )
        data_novel["batches_analysis"] = collect_person_unique_sample_values(
            values=data_original["batches_analysis"].dropna().to_list(),
            delimiter=",",
        )
        data_collection = data_collection.append(
            data_novel,
            ignore_index=True,
        )

    # Organize data.
    data_collection.reindex()
    data_collection.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    # Return information.
    return data_collection


# Organize persons' properties as covariates.


def collect_persons_unique_values(
    values=None,
    delimiter=None,
):
    """
    Collects unique values from a person's samples.

    arguments:
        values (list<str>): values from a person's samples
        delimiter (str): character for separation of values

    raises:

    returns:
        (str): values

    """

    #values_valid = list(filter(lambda value: not math.isnan(value), values))
    #values_type = list(map(lambda value: str(value), values_valid))
    values_combination = delimiter.join(values)
    values_split = values_combination.split(delimiter)
    values_strip = list(map(lambda value: value.strip(), values_split))
    values_unique = utility.collect_unique_elements(
        elements_original=values_strip
    )
    return values_unique


def expand_persons_categories_matrix(
    category=None,
    values=None,
    delimiter=None,
    data=None,
):
    """
    Expands a person's categorical properties to a binary matrix.

    arguments:
        category (str): name of a categorical property
        values (list<str>): all unique values of the categorical property
            across persons
        delimiter (str): character for separation of values
        data (object): Pandas data frame of persons' properties

    raises:

    returns:
        (object): Pandas data frame of a binary matrix

    """

    # Iterate on persons.
    persons = data.index.to_list()
    data_collection = pandas.DataFrame()
    for person in persons:
        # Create an empty matrix with values of zero for all facilities or
        # batches.
        data_empty = pandas.DataFrame(
            numpy.int(0),
            index=[0],
            columns=values
        )
        # Adjust matrix to match person's values of the categorical property.
        data_person = data_empty.copy(deep=True)
        data_person["person"] = person
        # Iterate on the person's values of the categorical property.
        values_person_raw = data.at[person, category]
        values_person = values_person_raw.split(delimiter)
        for value in values_person:
            data_person[value] = 1

        # Collect person's values.
        data_collection = data_collection.append(
            data_person,
            ignore_index=True,
        )

    # Organize data.
    data_collection.reset_index(
        level=None,
        inplace=True
    )
    data_collection.drop(
        labels="index",
        axis="columns",
        inplace=True
    )
    data_collection.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    # Return information.
    return data_collection


def organize_persons_properties_categories(
    data_persons_properties_raw=None,
):
    """
    Organizes information about persons' categorical properties.

    arguments:
        data_persons_properties_raw (object): Pandas data frame of persons'
            properties

    raises:

    returns:
        (dict): information about persons' categorical properties

    """

    # Organize data.
    data_copy = data_persons_properties_raw.copy(deep=True)

    # Extract and collect all unique values of facilities and batches across
    # all persons.
    facilities = collect_persons_unique_values(
        values=data_copy["facilities"].dropna().to_list(),
        delimiter=",",
    )
    batches_isolation = collect_persons_unique_values(
        values=data_copy["batches_isolation"].dropna().to_list(),
        delimiter=",",
    )
    batches_analysis = collect_persons_unique_values(
        values=data_copy["batches_analysis"].dropna().to_list(),
        delimiter=",",
    )

    # Summarize facilities and batches.
    utility.print_terminal_partition(level=2)
    print("Count of unique facilities: " + str(len(facilities)))
    print("Count of unique isolation batches: " + str(len(batches_isolation)))
    print("Count of unique analysis batches: " + str(len(batches_analysis)))
    utility.print_terminal_partition(level=2)

    # Create two-dimensional matrices to represent persons' categories.
    # Persons have values of zero (False) or one (True) to indicate whether
    # they had a sample in each facility or batch.
    data_facilities = expand_persons_categories_matrix(
        category="facilities",
        values=facilities,
        delimiter=",",
        data=data_copy,
    )
    data_batches_isolation = expand_persons_categories_matrix(
        category="batches_isolation",
        values=batches_isolation,
        delimiter=",",
        data=data_copy,
    )
    data_batches_analysis = expand_persons_categories_matrix(
        category="batches_analysis",
        values=batches_analysis,
        delimiter=",",
        data=data_copy,
    )
    # Merge batches to a single matrix.
    data_batches = data_batches_isolation.join(
        data_batches_analysis,
        how="left",
        on="person"
    )

    # Reduce dimensionality.
    report_facilities = utility.calculate_principal_components(
        data=data_facilities,
        components=3,
        report=True,
    )
    report_batches_isolation = utility.calculate_principal_components(
        data=data_batches_isolation,
        components=10,
        report=True,
    )
    report_batches_analysis = utility.calculate_principal_components(
        data=data_batches_analysis,
        components=5,
        report=True,
    )
    report_batches = utility.calculate_principal_components(
        data=data_batches,
        components=10,
        report=True,
    )

    # Compile information.
    information = dict()
    information["facilities"] = report_facilities
    information["batches_isolation"] = report_batches_isolation
    information["batches_analysis"] = report_batches_analysis
    information["batches"] = report_batches
    # Return information.
    return information


def insert_components(
    prefix=None,
    data_components=None,
    data_collection=None,
):
    """
    Inserts principal components of a property into a data frame.

    arguments:
        prefix (str): name to use as a prefix for components
        data_components (object): Pandas data frame of a property's principal
            components across observations
        data_collection (object): Pandas data frame in which to insert the
            property's principal components

    raises:

    returns:
        (object): Pandas data frame including property's principal components

    """

    # Copy data.
    data_components = data_components.copy(deep=True)
    data_collection = data_collection.copy(deep=True)

    # Change names of columns for property's principal components.
    components = data_components.shape[1]
    translations = dict()
    for count in range(components):
        original = ("component_" + str(count + 1))
        novel = (prefix + "_" + str(count + 1))
        translations[original] = novel
    data_components.rename(
        columns=translations,
        inplace=True,
    )

    # Insert components.
    data_combination = data_collection.join(
        data_components,
        how="left",
    )

    return data_combination


def organize_persons_properties(
    data_persons_properties_raw=None,
):
    """
    Organizes information about persons.

    arguments:
        data_persons_properties_raw (object): Pandas data frame of persons'
            properties

    raises:

    returns:
        (object): Pandas data frame of information about persons

    """

    # Organize data.
    data_persons_properties = data_persons_properties_raw.copy(deep=True)
    # Transform and reduce dimensionality on persons' categorical properties.
    bin = organize_persons_properties_categories(
        data_persons_properties_raw=data_persons_properties_raw,
    )
    data_facilities = bin["facilities"]["data_observations_components"]
    data_batches_isolation = bin["batches_isolation"][
        "data_observations_components"
    ]
    data_batches_analysis = bin["batches_analysis"][
        "data_observations_components"
    ]
    # Replace persons' raw properties.
    data_persons_properties = insert_components(
        prefix="facilities",
        data_components=data_facilities,
        data_collection=data_persons_properties,
    )
    data_persons_properties = insert_components(
        prefix="batches_isolation",
        data_components=data_batches_isolation,
        data_collection=data_persons_properties,
    )
    data_persons_properties = insert_components(
        prefix="batches_analysis",
        data_components=data_batches_analysis,
        data_collection=data_persons_properties,
    )

    # Organize variables for regression.
    data_persons_properties["female"] = data_persons_properties.apply(
        lambda row: 1 if (row["sex"] == "female") else 0,
        axis="columns",
    )
    data_persons_properties["season_sequence"] = data_persons_properties.apply(
        lambda row:
            1 if (row["season"] == "spring") else
            (2 if (row["season"] == "summer") else
            (3 if (row["season"] == "fall") else
            (4 if (row["season"] == "winter") else
            float("nan")))),
        axis="columns",
    )
    data_persons_properties["spring"] = data_persons_properties.apply(
        lambda row: 1 if (row["season"] == "spring") else 0,
        axis="columns",
    )
    data_persons_properties["summer"] = data_persons_properties.apply(
        lambda row: 1 if (row["season"] == "summer") else 0,
        axis="columns",
    )
    data_persons_properties["fall"] = data_persons_properties.apply(
        lambda row: 1 if (row["season"] == "fall") else 0,
        axis="columns",
    )
    data_persons_properties["winter"] = data_persons_properties.apply(
        lambda row: 1 if (row["season"] == "winter") else 0,
        axis="columns",
    )

    # Return information.
    return data_persons_properties


# Count persons' by categories of sex and age.


def count_persons_sex_age(
    data_persons_properties=None,
):
    """
    Counts persons in each group by sex and age.

    arguments:
        data_persons_properties (object): Pandas data frame of persons and
            their properties

    raises:

    returns:
        (object): Pandas data frame of information about persons
    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    # Organize data.
    data_persons_properties.reset_index(
        level=None,
        inplace=True
    )
    data_persons_properties = data_persons_properties.loc[
        :, data_persons_properties.columns.isin([
            "person", "sex", "decade",
        ])
    ]
    data_persons_properties.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_persons_properties.reindex()
    data_persons_properties.set_index(
        ["sex", "decade"],
        append=False,
        drop=True,
        inplace=True
    )

    data_persons_sex_age_counts = data_persons_properties.groupby(
        level=["sex", "decade"],
        sort=True,
        as_index=False
    ).size().to_frame(
        name="count"
    )
    data_persons_sex_age_counts.reset_index(
        level=["sex", "decade"], inplace=True
    )

    # Return information.
    return data_persons_sex_age_counts


# Organize covariates for heritability analysis.


def organize_heritability_covariates(
    data_persons_properties=None,
):
    """
    Organizes information about families and persons for heritability analysis.

    As the data from GTEx do not include families, generate distinct families
    for each person for use in heritability analysis in GCTA.

    arguments:
        data_persons_properties (object): Pandas data frame of persons their
            properties

    raises:

    returns:
        (dict): information for heritability analysis
    """

    # Introduce empty identifiers for families for compatibility with Plink and
    # GCTA.
    # Organize data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    data_persons_properties["family"] = 0
    data_persons_properties.reset_index(
        level=None,
        inplace=True
    )
    # Extract families and persons.
    data_copy = data_persons_properties.copy(deep=True)
    data_families_persons = data_copy.loc[
        :, data_copy.columns.isin(["person", "family"])
    ]
    data_families_persons = data_families_persons[[
        "family",
        "person",
    ]]
    # Extract categorical covariates.
    data_copy = data_persons_properties.copy(deep=True)
    data_persons_categories = data_copy.loc[
        :, data_copy.columns.isin([
            "person", "family", "sex", "season"
        ])
    ]
    data_persons_categories = data_persons_categories[[
        "family",
        "person",
        "sex",
        "season",
    ]]
    # Extract quantitative covariates.
    data_copy = data_persons_properties.copy(deep=True)
    data_persons_quantities = data_copy.loc[
        :, data_copy.columns.isin([
            "person", "family", "age", "body", "hardiness", "delay"
        ])
    ]
    data_persons_quantities = data_persons_quantities[[
        "family",
        "person",
        "age",
        "body",
        "hardiness",
        "delay",
    ]]

    # Compile information.
    information = {
        "data_families_persons": data_families_persons,
        "data_persons_categories": data_persons_categories,
        "data_persons_quantities": data_persons_quantities,
    }

    # Return information.
    return information


def extract_organize_information(
    data_samples_tissues_persons=None,
    data_gene_signal=None,
):
    """
    Extracts and organizes information about samples and genes for further
    analysis.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (dict): information
    """

    # Copy data.
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    data_gene_signal = data_gene_signal.copy(deep=True)

    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_signal.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_gene_signal_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_sample=data_transposition,
    )

    # Extract information about samples, tissues, and persons.
    collection = extract_persons_tissues_samples(
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    # Count tissues per person and persons per tissue.
    counts = count_tissues_persons_groups(
        data_samples_tissues_persons=data_samples_tissues_persons,
    )

    # Extract information about persons' properties.
    data_persons_properties_raw = extract_persons_properties(
        data_samples_tissues_persons=data_samples_tissues_persons,
    )

    # Expand covariates.
    # Prepare covariates for regression.
    data_persons_properties = organize_persons_properties(
        data_persons_properties_raw=data_persons_properties_raw,
    )

    # Count persons in groups by sex and age.
    data_persons_sex_age_counts = count_persons_sex_age(
        data_persons_properties=data_persons_properties_raw,
    )
    # Organize information for heritability analysis.
    heritability = organize_heritability_covariates(
        data_persons_properties=data_persons_properties_raw,
    )

    # Compile information.
    information = {
        "data_gene_signal_factor": data_gene_signal_factor,

        "persons": collection["persons"],
        "tissues": collection["tissues"],
        "samples": collection["samples"],

        "data_tissues_per_person": counts["data_tissues_per_person"],
        "tissues_per_person": counts["tissues_per_person"],
        "data_persons_per_tissue": counts["data_persons_per_tissue"],
        "persons_per_tissue": counts["persons_per_tissue"],

        "data_persons_sex_age_counts": data_persons_sex_age_counts,

        "data_persons_properties_raw": data_persons_properties_raw,
        "data_persons_properties": data_persons_properties,
        "data_families_persons": heritability["data_families_persons"],
        "data_persons_categories": heritability["data_persons_categories"],
        "data_persons_quantities": heritability["data_persons_quantities"],
    }
    # Return information.
    return information


##########
# Loose and tight routines.


def select_organize_samples_genes_signals(
    data_gene_annotation_gtex=None,
    data_gene_annotation_gencode=None,
    data_samples_tissues_persons=None,
    data_gene_signal=None,
    stringency=None,
):
    """
    Selects samples and genes by signals beyond threshold.

    Data format should have genes across rows and samples across columns.

    arguments:
        data_gene_annotation_gtex (object): Pandas data frame of genes'
            annotations from GTEx
        data_gene_annotation_gencode (object): Pandas data frame of genes'
            annotations from GENCODE
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples
        stringency (str): category, loose or tight, of selection criteria

    raises:

    returns:
        (dict): information about samples, genes, and signals
    """

    utility.print_terminal_partition(level=1)
    print("Sample selection criteria: " + stringency)
    utility.print_terminal_partition(level=1)

    # Copy data.
    data_gene_annotation_gtex = data_gene_annotation_gtex.copy(deep=True)
    data_gene_annotation_gencode = data_gene_annotation_gencode.copy(deep=True)
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    data_gene_signal = data_gene_signal.copy(deep=True)

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_signal=data_gene_signal,
    )

    # Select genes that encode proteins and are on autosomes.
    # Count of genes: 18353
    utility.print_terminal_partition(level=1)
    print("Selection of genes.")
    utility.print_terminal_partition(level=2)
    print("GTEx reference gene annotations.")
    # GTEx reference: 18380 genes on autosomes that encode proteins.
    data_gene_annotation_gtex = select_gene_annotation(
        data_gene_annotation=data_gene_annotation_gtex,
    )
    utility.print_terminal_partition(level=2)
    print("GENCODE reference gene annotations.")
    # GENCODE reference: 19035 genes on autosomes that encode proteins.
    data_gene_annotation_gencode = select_gene_annotation(
        data_gene_annotation=data_gene_annotation_gencode,
    )
    # Use genes' annotations from GTEx.
    utility.print_terminal_partition(level=2)
    print("Select signal genes by GTEx reference.")
    # GTEx signal genes on autosomes that encode proteins: 18380
    data_gene_signal_gene_gtex = select_genes(
        data_gene_annotation=data_gene_annotation_gtex,
        data_gene_signal=data_gene_signal,
    )

    # Use genes' annotations from GENCODE.
    utility.print_terminal_partition(level=2)
    print("Select signal genes by GENCODE reference.")
    # GTEx signal genes on autosomes that encode proteins: 19035
    data_gene_signal_gene_gencode = select_genes(
        data_gene_annotation=data_gene_annotation_gencode,
        data_gene_signal=data_gene_signal,
    )

    # Continue analyses with selection of genes by GENCODE.
    data_gene_signal_gene = data_gene_signal_gene_gencode.copy(deep=True)

    # Select samples by tissues of interest.
    # Exclude samples for cell lines that are not whole tissues.
    # Exclude samples for tissues with coverage by <100 samples.
    data_samples_exclusion = exclude_samples_by_tissues(
        tissues_minor=["fibroblast", "lymphocyte"],
        tissues_major=["bladder", "cervix", "fallopius", "kidney"],
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    if stringency == "loose":
        # Select all tissues with adequate coverage.
        tissues_selection = [
            "adipose", "adrenal", "artery", "blood", "brain", "breast",
            "colon", "esophagus", "heart", "intestine", "liver", "lung",
            "muscle", "nerve", "ovary", "pancreas", "pituitary", "prostate",
            "salivary", "skin", "spleen", "stomach", "testis", "thyroid",
            "uterus", "vagina",
        ]
    elif stringency == "tight":
        # Select sex-neutral tissues.
        tissues_selection = [
            "adipose", "adrenal", "artery", "blood", "brain", "colon",
            "esophagus", "heart", "intestine", "liver", "lung", "muscle",
            "nerve", "pancreas", "pituitary", "salivary", "skin", "spleen",
            "stomach", "thyroid",
        ]
    data_samples_inclusion = include_samples_by_tissues(
        tissues_major=tissues_selection,
        data_samples_tissues_persons=data_samples_exclusion,
    )
    samples_preliminary = data_samples_inclusion.index.to_list()
    data_gene_signal_sample = select_samples_genes_signals(
        samples=samples_preliminary,
        data_gene_signal=data_gene_signal_gene,
    )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=data_samples_inclusion,
        data_gene_signal=data_gene_signal_sample,
    )

    utility.print_terminal_partition(level=1)
    print("Selection of genes and samples by signal coverage.")
    utility.print_terminal_partition(level=1)

    # Select genes and samples by signals and coverage.
    if stringency == "loose":
        # Each gene must have signal greater than or equal to 0.1 in at least
        # 1% of samples.
        # Each sample must have signal greater than or equal to 0.1 in at least
        # 10% of genes.
        data_gene_signal_selection = select_samples_genes_signals_coverage(
            threshold=0.1,
            proportion_gene=0.01, # 0.01
            proportion_sample=0.1, # 0.1
            data_gene_signal=data_gene_signal_sample,
        )
    elif stringency == "tight":
        # Each gene must have signal greater than or equal to 0.1 TPM in at
        # least 10% of samples.
        # Each sample must have signal greater than or equal to 0.1 TPM in at
        # least 50% of genes.
        data_gene_signal_selection = select_samples_genes_signals_coverage(
            threshold=0.1,
            proportion_gene=0.1, # 0.1
            proportion_sample=0.5, # 0.5
            data_gene_signal=data_gene_signal_sample,
        )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=data_samples_inclusion,
        data_gene_signal=data_gene_signal_selection,
    )

    # Select samples, tissues, and persons from filtered genes' signals.
    data_samples_tissues_persons_selection = select_samples_tissues_persons(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_signal=data_gene_signal_selection,
    )
    utility.print_terminal_partition(level=2)
    print("selection of samples, tissues, and persons")
    print(data_samples_tissues_persons_selection)

    # Extract and organize information for subsequent analyses.
    information = extract_organize_information(
        data_samples_tissues_persons=data_samples_tissues_persons_selection,
        data_gene_signal=data_gene_signal_selection,
    )
    # Compile information.
    information["data_gene_annotation_gtex"] = data_gene_annotation_gtex
    information["data_gene_annotation_gencode"] = data_gene_annotation_gencode
    information["data_samples_tissues_persons"] = (
        data_samples_tissues_persons_selection
    )
    information["data_gene_signal"] = data_gene_signal_selection
    # Return information.
    return information


##########
# Product.


def write_product(
    dock=None,
    stringency=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    utility.create_directories(path_selection)

    path_gene_annotation_gtex = os.path.join(
        path_selection, "data_gene_annotation_gtex.pickle"
    )
    path_gene_annotation_gtex_text = os.path.join(
        path_selection, "data_gene_annotation_gtex.tsv"
    )
    path_gene_annotation_gencode = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )
    path_gene_annotation_gencode_text = os.path.join(
        path_selection, "data_gene_annotation_gencode.tsv"
    )

    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    path_gene_signal_factor = os.path.join(
        path_selection, "data_gene_signal_factor.pickle"
    )

    path_samples_tissues_persons = os.path.join(
        path_selection, "data_samples_tissues_persons.pickle"
    )
    path_samples_tissues_persons_text = os.path.join(
        path_selection, "data_samples_tissues_persons.tsv"
    )
    path_persons = os.path.join(
        path_selection, "persons.txt"
    )
    path_tissues = os.path.join(
        path_selection, "tissues.txt"
    )
    path_samples = os.path.join(
        path_selection, "samples.txt"
    )

    path_data_tissues_per_person = os.path.join(
        path_selection, "data_tissues_per_person.pickle"
    )
    path_tissues_per_person = os.path.join(
        path_selection, "tissues_per_person.pickle"
    )
    path_data_persons_per_tissue = os.path.join(
        path_selection, "data_persons_per_tissue.pickle"
    )
    path_persons_per_tissue = os.path.join(
        path_selection, "persons_per_tissue.pickle"
    )

    path_data_persons_sex_age_counts = os.path.join(
        path_selection, "data_persons_sex_age_counts.pickle"
    )

    path_persons_properties_raw = os.path.join(
        path_selection, "data_persons_properties_raw.pickle"
    )
    path_persons_properties_raw_text = os.path.join(
        path_selection, "data_persons_properties_raw.tsv"
    )
    path_persons_properties = os.path.join(
        path_selection, "data_persons_properties.pickle"
    )
    path_persons_properties_text = os.path.join(
        path_selection, "data_persons_properties.tsv"
    )

    path_families_persons = os.path.join(
        path_selection, "data_families_persons.pickle"
    )
    path_families_persons_text = os.path.join(
        path_selection, "families_persons.tsv"
    )
    path_persons_categories = os.path.join(
        path_selection, "persons_categories.tsv"
    )
    path_persons_quantities = os.path.join(
        path_selection, "persons_quantities.tsv"
    )

    # Write information to file.

    pandas.to_pickle(
        information["data_gene_annotation_gtex"],
        path_gene_annotation_gtex
    )
    pandas.to_pickle(
        information["data_gene_annotation_gencode"],
        path_gene_annotation_gencode
    )
    information["data_gene_annotation_gtex"].to_csv(
        path_or_buf=path_gene_annotation_gtex_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_gene_annotation_gencode"].to_csv(
        path_or_buf=path_gene_annotation_gencode_text,
        sep="\t",
        header=True,
        index=True,
    )

    pandas.to_pickle(
        information["data_gene_signal"],
        path_gene_signal
    )
    pandas.to_pickle(
        information["data_gene_signal_factor"],
        path_gene_signal_factor
    )

    pandas.to_pickle(
        information["data_samples_tissues_persons"],
        path_samples_tissues_persons
    )
    information["data_samples_tissues_persons"].to_csv(
        path_or_buf=path_samples_tissues_persons_text,
        sep="\t",
        header=True,
        index=True,
    )
    utility.write_file_text_list(
        elements=information["persons"],
        delimiter="\n",
        path_file=path_persons
    )
    utility.write_file_text_list(
        elements=information["tissues"],
        delimiter="\n",
        path_file=path_tissues
    )
    utility.write_file_text_list(
        elements=information["samples"],
        delimiter="\n",
        path_file=path_samples
    )

    pandas.to_pickle(
        information["data_tissues_per_person"],
        path_data_tissues_per_person
    )
    pandas.to_pickle(
        information["tissues_per_person"],
        path_tissues_per_person
    )
    pandas.to_pickle(
        information["data_persons_per_tissue"],
        path_data_persons_per_tissue
    )
    pandas.to_pickle(
        information["persons_per_tissue"],
        path_persons_per_tissue
    )

    pandas.to_pickle(
        information["data_persons_sex_age_counts"],
        path_data_persons_sex_age_counts
    )

    pandas.to_pickle(
        information["data_persons_properties_raw"],
        path_persons_properties_raw
    )
    information["data_persons_properties_raw"].to_csv(
        path_or_buf=path_persons_properties_raw_text,
        sep="\t",
        header=True,
        index=True,
    )
    pandas.to_pickle(
        information["data_persons_properties"],
        path_persons_properties
    )
    information["data_persons_properties"].to_csv(
        path_or_buf=path_persons_properties_text,
        sep="\t",
        header=True,
        index=True,
    )

    pandas.to_pickle(
        information["data_families_persons"],
        path_families_persons
    )
    information["data_families_persons"].to_csv(
        path_or_buf=path_families_persons_text,
        sep="\t",
        na_rep="NA",
        header=False,
        index=False,
    )
    information["data_persons_categories"].to_csv(
        path_or_buf=path_persons_categories,
        sep="\t",
        na_rep="NA",
        header=False,
        index=False,
    )
    information["data_persons_quantities"].to_csv(
        path_or_buf=path_persons_quantities,
        sep="\t",
        na_rep="NA",
        header=False,
        index=False,
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
    path_selection = os.path.join(dock, "selection")
    utility.remove_directory(path=path_selection)

    # Read source information from file.
    source = read_source(dock=dock)
    # Organize data.
    data_gene_signal = organize_data_axes_indices(
        data=source["data_gene_signal"]
    )

    # Selection 1: loose.
    # Select data for general analysis.
    # 1. select genes on autosomes that encode proteins
    # 2. select samples from whole tissues, not cell lines
    # 3. select samples from tissues with coverage of >=100 samples
    # 4. select genes with signal in >=1% of samples
    # 5. select samples with signal in >=10% of genes
    if False:
        loose = select_organize_samples_genes_signals(
            data_gene_annotation_gtex=source["data_gene_annotation_gtex"],
            data_gene_annotation_gencode=source["data_gene_annotation_gencode"],
            data_samples_tissues_persons=source["data_samples_tissues_persons"],
            data_gene_signal=data_gene_signal,
            stringency="loose",
        )
    # Write product information to file.
    if False:
        write_product(
            dock=dock,
            stringency="loose",
            information=loose,
        )

    # Selection 2: tight.
    # Select data for pantissue aggregate analysis.
    # Criteria:
    # 1. select genes on autosomes that encode proteins
    # 2. select samples from whole tissues, not cell lines
    # 3. select samples from tissues with coverage of >=100 samples
    # 4. select samples from sex-neutral tissues
    # 5. select genes with signal in >=10% of samples
    # 6. select samples with signal in >=50% of genes
    tight = select_organize_samples_genes_signals(
        data_gene_annotation_gtex=source["data_gene_annotation_gtex"],
        data_gene_annotation_gencode=source["data_gene_annotation_gencode"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal,
        stringency="tight",
    )
    # Write product information to file.
    write_product(
        dock=dock,
        stringency="tight",
        information=tight,
    )

    pass


if (__name__ == "__main__"):
    execute_procedure()
