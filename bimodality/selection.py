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
import organization
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
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_gene_count = os.path.join(
        path_assembly, "data_gene_count.feather"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.feather"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(
        path_gene_annotation
    )
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_count = pandas.read_feather(
        path=path_gene_count
    )
    data_gene_signal = pandas.read_feather(
        path=path_gene_signal,
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_count": data_gene_count,
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


# Genes.


def select_genes(
    data_gene_annotation=None,
    data_gene_signal=None
):
    """
    Selects genes that encode proteins.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and patients.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

    # Due to previous filter in assembly procedure, annotation genes are only
    # protein-coding.

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
        "count of GENCODE genes of type 'protein_coding': " +
        str(len(genes_protein))
    )

    # Filter gene signals.
    # Remove indices.
    data_gene_signal.reset_index(
        level=None,
        inplace=True
    )
    data_gene_signal_protein = data_gene_signal.loc[
        data_gene_signal["gene"].isin(genes_protein), :
    ]

    genes = data_gene_signal_protein.index.to_list()
    print(
        "signal genes that encode proteins: " +
        str(data_gene_signal_protein.shape[0])
    )
    print(
        "signal genes that encode proteins: " +
        str(len(genes))
    )
    print(
        "unique signal genes that encode proteins: " +
        str(len(utility.collect_unique_elements(elements_original=genes)))
    )

    utility.print_terminal_partition(level=2)

    # Organize data.
    data_gene_signal_protein = organize_data_axes_indices(
        data=data_gene_signal_protein
    )

    return data_gene_signal_protein


# Samples.


def select_samples_by_exclusion(
    data_samples_tissues_persons=None,
):
    """
    Selects samples by exclusion of specific tissues or persons.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (list<str>): identifiers of samples

    """

    # Define tissues for exclusion.
    # These tissues are extra-corporeal tissues, such as cell lines.
    # skin, fibroblast
    # blood, lymphocyte
    tissues_extracorp = []
    tissues_extracorp.append("fibroblast")
    tissues_extracorp.append("lymphocyte")
    data_samples_corp = data_samples_tissues_persons.loc[
        ~data_samples_tissues_persons["tissue_minor"].isin(tissues_extracorp),
    ]
    # These tissues have fewer than 100 total samples.
    tissues_coverage = []
    tissues_coverage.append("bladder")
    tissues_coverage.append("cervix")
    tissues_coverage.append("fallopius")
    tissues_coverage.append("kidney")
    data_samples_coverage = data_samples_corp.loc[
        ~data_samples_corp["tissue_major"].isin(tissues_coverage),
    ]
    # Extract identifiers of samples.
    samples = data_samples_coverage.index.to_list()
    # Return information.
    return samples


def select_samples_by_inclusion(
    data_samples_tissues_persons=None,
):
    """
    Selects samples by inclusion of specific tissues or persons.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (list<str>): identifiers of samples

    """

    # Define tissues for inclusion.
    # These tissues are not sex-specific.
    tissues_sex = [
        "adipose", # 797
        "adrenal", # 258
        "artery", # 786
        "blood", # 767
        "brain", # 382
        #"breast", # 459
        "colon", # 555
        "esophagus", # 710
        "heart", # 561
        "intestine", # 187
        "liver", # 226
        "lung", # 578
        "muscle", # 803
        "nerve", # 619
        #"ovary", # 180
        "pancreas", # 328
        "pituitary", # 283
        #"prostate", # 245
        "salivary", # 162
        "skin", # 912
        "spleen", # 241
        "stomach", # 359
        #"testis", # 361
        "thyroid", # 653
        #"uterus", # 142
        #"vagina", # 156
    ]

    samples = select_samples_by_tissue(
        tissues=tissues_sex,
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    # Return information.
    return samples


def select_samples_by_tissue(
    tissues=None,
    data_samples_tissues_persons=None,
):
    """
    Selects samples by specific tissues.

    arguments:
        tissues (list<str>): identifiers of tissues
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (list<str>): identifiers of samples

    """

    # Select samples by major tissues.
    data_samples = data_samples_tissues_persons.loc[
        data_samples_tissues_persons["tissue_major"].isin(tissues), :
    ]
    # Extract identifiers of samples.
    samples = data_samples.index.to_list()
    # Return information.
    return samples


def select_samples(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Selects genes' signals for samples of interest.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    # Select samples by exclusion.
    # Exclude samples for extra-corporeal tissues, such as cell lines.
    # Exclude samples for tissues with fewer than 50 total samples.
    samples_exclusion = select_samples_by_exclusion(
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    #print(len(samples_exclusion))
    # Select samples by inclusion.
    samples_inclusion = select_samples_by_inclusion(
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    #print(len(samples_inclusion))
    # Select samples that are in both exclusion and inclusion lists.
    samples = utility.filter_common_elements(
        list_one=samples_exclusion,
        list_two=samples_inclusion,
    )
    utility.print_terminal_partition(level=2)
    print("Count of samples to select: " + str(len(samples)))

    # Select genes' signals for samples of interest.
    utility.print_terminal_partition(level=1)
    data_selection = data_gene_signal.loc[:, samples]

    # Return information.
    return data_selection


def select_samples_genes_by_few_tissues(
    tissues_major=None,
    data_gene_sample=None
):
    """
    Selects samples and genes by tissues of interest.

    arguments:
        tissues_major (list<str>): major tissues of interest
        data_gene_sample (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons

    """

    data_gene_sample.reset_index(
        level=["tissue_major", "tissue_minor"],
        inplace=True
    )
    data_selection = data_gene_sample.loc[
        data_gene_sample["tissue_major"].isin(tissues_major),
    ]
    data_selection.set_index(
        ["tissue_major", "tissue_minor"],
        append=False,
        drop=True,
        inplace=True
    )
    return data_selection


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
    threshold=None,
    proportion=None,
    data_gene_signal=None
):
    """
    Selects samples and genes by signals beyond threshold.

    Data format should have genes across rows and samples across columns.

    arguments:
        threshold (float): minimal gene's signal for a sample
        proportion (float): minimal proportion of samples or genes that must
            pass threshold
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    utility.print_terminal_partition(level=2)
    print("shape of original data...")
    print(data_gene_signal.shape)

    # Filter genes by their signals across samples.
    # Filter to keep only genes with signals beyond threshold in at least one
    # sample.
    data_row = filter_rows_columns_by_threshold_proportion(
        data=data_gene_signal,
        dimension="row",
        threshold=threshold,
        proportion=proportion
    )
    print("shape of data after signal filter on genes...")
    print(data_row.shape)
    # Filter samples by their signals across genes.
    # Filter to keep only samples with signals beyond threshold in at least one
    # gene.
    data_column = filter_rows_columns_by_threshold_proportion(
        data=data_row,
        dimension="column",
        threshold=threshold,
        proportion=proportion
    )
    print("shape of data after signal filter on samples...")
    print(data_column.shape)

    # Return information.
    return data_column


# Persons, tissues, samples.


def extract_gene_signal_persons_tissues_samples(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Extracts persons, tissues, and samples from selection of genes' signals.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (dict<list<str>>): collection of persons, tissues, and samples
    """

    utility.print_terminal_partition(level=2)
    print("extract samples, persons, and tissues")

    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_signal.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_gene_signal_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_sample=data_transposition,
    )
    print(data_gene_signal_factor)
    # Remove indices.
    data_gene_signal_factor.reset_index(
        level=None,
        inplace=True
    )

    # Extract selection samples.
    samples = utility.collect_unique_elements(
        elements_original=data_gene_signal_factor["sample"].to_list()
    )
    # Extract selection persons.
    persons = utility.collect_unique_elements(
        elements_original=data_gene_signal_factor["person"].to_list()
    )
    # Extract selection tissues.
    tissues_major = utility.collect_unique_elements(
        elements_original=data_gene_signal_factor["tissue_major"].to_list()
    )

    # Summary.
    utility.print_terminal_partition(level=2)
    print("count of unique samples: " + str(len(samples)))
    print("count of unique persons: " + str(len(persons)))
    print("count of unique tissues_major: " + str(len(tissues_major)))


    # Compile information.
    information = {
        "samples": samples,
        "persons": persons,
        "tissues_major": tissues_major,
    }

    # Return information.
    return information


##########
# Product.


def write_product(dock=None, information=None):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files.
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection")
    utility.create_directory(path_selection)
    path_gene_count = os.path.join(
        path_selection, "data_gene_count.pickle"
    )
    path_gene_count_factor = os.path.join(
        path_selection, "data_gene_count_factor.pickle"
    )
    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    path_gene_signal_factor = os.path.join(
        path_selection, "data_gene_signal_factor.pickle"
    )

    path_samples = os.path.join(
        path_selection, "samples.txt"
    )
    path_persons = os.path.join(
        path_selection, "persons.txt"
    )
    path_tissues = os.path.join(
        path_selection, "tissues.txt"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_gene_count"],
        path_gene_count
    )
    pandas.to_pickle(
        information["data_gene_count_factor"],
        path_gene_count_factor
    )
    pandas.to_pickle(
        information["data_gene_signal"],
        path_gene_signal
    )
    pandas.to_pickle(
        information["data_gene_signal_factor"],
        path_gene_signal_factor
    )

    utility.write_file_text_list(
        elements=information["samples"],
        delimiter="\t",
        path_file=path_samples
    )
    utility.write_file_text_list(
        elements=information["persons"],
        delimiter="\t",
        path_file=path_persons
    )
    utility.write_file_text_list(
        elements=information["tissues_major"],
        delimiter="\t",
        path_file=path_tissues
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
    data_gene_count = organize_data_axes_indices(
        data=source["data_gene_count"]
    )
    data_gene_signal = organize_data_axes_indices(
        data=source["data_gene_signal"]
    )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal,
    )
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_count,
    )

    # Select genes that encode proteins.
    data_gene_signal_gene = select_genes(
        data_gene_annotation=source["data_gene_annotation"],
        data_gene_signal=data_gene_signal,
    )
    data_gene_count_gene = select_genes(
        data_gene_annotation=source["data_gene_annotation"],
        data_gene_signal=data_gene_count,
    )

    # Select samples by tissues.
    data_gene_signal_sample = select_samples(
        data_gene_signal=data_gene_signal_gene,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
    )
    data_gene_count_sample = select_samples(
        data_gene_signal=data_gene_count_gene,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
    )

    # Select genes and samples by signals.
    data_gene_signal_selection = select_samples_genes_signals(
        threshold=1.0,
        proportion=0.5,
        data_gene_signal=data_gene_signal_sample,
    )
    data_gene_count_selection = select_samples_genes_signals(
        threshold=1.0,
        proportion=0.1,
        data_gene_signal=data_gene_count_sample,
    )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_count_selection,
    )
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal_selection,
    )

    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_count_selection.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_gene_count_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_sample=data_transposition,
    )
    print(data_gene_count_factor)

    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_signal_selection.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_gene_signal_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_sample=data_transposition,
    )
    print(data_gene_signal_factor)

    # Extract gene signal persons, tissues, samples.
    collection = extract_gene_signal_persons_tissues_samples(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal_selection,
    )

    # Compile information.
    information = {
        "data_gene_count": data_gene_count_selection,
        "data_gene_count_factor": data_gene_count_factor,
        "data_gene_signal": data_gene_signal_selection,
        "data_gene_signal_factor": data_gene_signal_factor,
        "samples": collection["samples"],
        "persons": collection["persons"],
        "tissues_major": collection["tissues_major"],
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
