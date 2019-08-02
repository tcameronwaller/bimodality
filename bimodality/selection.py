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
        path_assembly, "data_gene_count.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(
        path_gene_annotation
    )
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_count = pandas.read_pickle(path_gene_count)
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_count": data_gene_count,
        "data_gene_signal": data_gene_signal,
    }


def summarize_samples_genes(
    data_samples_tissues_persons=None,
    data_gene_count=None,
    data_gene_signal=None
):
    """
    Summarize selection of samples and genes.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_count (object): Pandas data frame of genes' counts across
            samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:

    """

    # Copy data.
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    data_gene_count = data_gene_count.copy(deep=True)
    data_gene_signal = data_gene_signal.copy(deep=True)

    utility.print_terminal_partition(level=1)
    print("Summary of selection of samples and genes.")
    # Counts of samples and genes.
    print("Genes' counts...")
    print("Count of samples: " + str(data_gene_count.shape[1]))
    print("Count of genes: " + str(data_gene_count.shape[0]))

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
    groups = data_factor.groupby(level=["tissue_major"])
    print("Count of major tissues: " + str(len(groups)))
    for name, group in groups:
        print(name)
    data_factor.reset_index(
        level=["tissue_major"],
        inplace=True
    )
    data_adipose = data_factor.loc[data_factor["tissue_major"] == "adipose", :]
    print(data_adipose)

    pass


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

    def check_gene_type(identifier=None):
        type = data_gene_annotation.at[identifier, "gene_type"]
        return type == "protein_coding"
    #data = data.loc[lambda identifier: check_gene_type(identifier)]

    utility.print_terminal_partition(level=2)
    print(
        "Selection of genes that encode proteins."
    )
    # Describe original count of genes.
    # Signal matrix has genes on the index dimension.


    print("signal genes, original: " + str(data_gene_signal.shape[0]))
    genes_signal = data_gene_signal.index.to_list()
    print("signal genes, original: " + str(len(genes_signal)))
    # Filter genes by their annotations.
    #print(data_gene_annotation.loc["ENSG00000223972", "gene_type"])
    genes_protein = data_gene_annotation.index.to_list()
    print(
        "count of GENCODE genes of type 'protein_coding': " +
        str(len(genes_protein))
    )

    # Filter gene signals.
    genes_signal_protein = utility.filter_common_elements(
        list_one=genes_protein, list_two=genes_signal
    )
    print(
        "signal genes that encode proteins: " +
        str(len(genes_signal_protein))
    )
    data_gene_signal = data_gene_signal.loc[genes_signal_protein, :]
    print(
        "signal genes that encode proteins: " + str(data_gene_signal.shape[0])
    )
    return data_gene_signal


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
    # These tissues have fewer than 50 total samples.
    tissues_coverage = []
    tissues_coverage.append("kidney")
    tissues_coverage.append("bladder")
    tissues_coverage.append("cervix")
    tissues_coverage.append("fallopius")
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
        "adipose", # 552
        "adrenal", # 190
        "artery", # 551
        "blood", # 407
        "brain", # 254
        "colon", # 371
        "esophagus", # 513
        "heart", # 399
        "liver", # 175
        "lung", # 427
        "muscle", # 564
        "nerve", # 414
        "pancreas", # 248
        "pituitary", # 183
        "skin", # 583
        "intestine", # 137
        #"salivary", # 97
        "spleen", # 162
        "stomach", # 262
        "thyroid", # 446
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
        data_samples_tissues_persons["tissue_major"].isin(tissues),
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
    print(len(samples_exclusion))
    # Select samples by inclusion.
    samples_inclusion = select_samples_by_inclusion(
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    print(len(samples_inclusion))
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


def filter_rows_columns_by_threshold_count(
    data=None,
    dimension=None,
    threshold=None,
    count=None,
):
    """
    Filters either rows or columns.

    Persistence of a row or column requires at least a specific count of values
    beyond a specific threshold.

    arguments:
        data (object): Pandas data frame of genes' signals across samples
        dimension (str): dimension to filter, either "row" or "column"
        threshold (float): minimal signal
        count (int): minimal count

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = values[values]
        return (len(values_true) >= count)

    # Determine whether values exceed threshold.
    data_threshold = (data >= threshold)
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )
    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count.any(axis="columns"), : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count.any(axis="index")]
    return data_pass


def select_samples_genes_signals(
    threshold=None,
    count=None,
    data_gene_signal=None
):
    """
    Selects samples and genes by signals beyond threshold.

    arguments:
        threshold (float): minimal signal
        count (int): minimal count of samples or genes with minimal signals
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    # Filter genes by signal.
    # Filter to keep only genes with signals beyond threshold in at least one
    # sample.
    data_row = filter_rows_columns_by_threshold_count(
        data=data_gene_signal,
        dimension="row",
        threshold=threshold,
        count=count
    )
    # Filter samples by signal.
    # Filter to keep only samples with signals beyond threshold in at least one
    # gene.
    data_column = filter_rows_columns_by_threshold_count(
        data=data_row,
        dimension="column",
        threshold=threshold,
        count=count
    )

    # Return information.
    return data_column


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
    utility.confirm_path_directory(path_selection)
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

    # Read source information from file.
    source = read_source(dock=dock)

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_count=source["data_gene_count"],
        data_gene_signal=source["data_gene_signal"],
    )

    # Select genes that encode proteins.
    data_gene_signal_gene = select_genes(
        data_gene_annotation=source["data_gene_annotation"],
        data_gene_signal=source["data_gene_signal"]
    )
    data_gene_count_gene = select_genes(
        data_gene_annotation=source["data_gene_annotation"],
        data_gene_signal=source["data_gene_count"]
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
        count=10,
        data_gene_signal=data_gene_signal_sample,
    )
    data_gene_count_selection = select_samples_genes_signals(
        threshold=1.0,
        count=10,
        data_gene_signal=data_gene_count_sample,
    )


    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_count=data_gene_count_selection,
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

    # Compile information.
    information = {
        "data_gene_count": data_gene_count_selection,
        "data_gene_count_factor": data_gene_count_factor,
        "data_gene_signal": data_gene_signal_selection,
        "data_gene_signal_factor": data_gene_signal_factor,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
