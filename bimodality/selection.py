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
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_measurement = os.path.join(dock, "measurement")
    path_gene_count = os.path.join(
        path_measurement, "data_gene_count.pickle"
    )
    path_gene_signal = os.path.join(
        path_measurement, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_count = pandas.read_pickle(path_gene_count)
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
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
        "adipose",
        "adrenal",
        "artery",
        "blood",
        "brain",
        "colon",
        "esophagus",
        "heart",
        "liver",
        "lung",
        "muscle",
        "nerve",
        "pancreas",
        "pituitary",
        "skin",
        "intestine",
        "salivary",
        "spleen",
        "stomach",
        "thyroid",
    ]
    data_samples_sex = data_samples_tissues_persons.loc[
        data_samples_tissues_persons["tissue_major"].isin(tissues_sex),
    ]
    # Extract identifiers of samples.
    samples = data_samples_sex.index.to_list()
    # Return information.
    return samples


def select_samples_genes(
    samples=None,
    data_gene_signal=None
):
    """
    Selects genes' signals for samples of interest.

    arguments:
        samples (list<str>): identifiers of samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    # Select genes' signals for samples of interest.
    utility.print_terminal_partition(level=1)
    data_sample = data_gene_signal.loc[:, samples]

    # Filter genes by signal.
    # Filter to keep only genes with signals beyond threshold in at least one
    # sample.
    data_detection_gene = measurement.filter_genes_by_signal_threshold(
        data=data_sample,
        threshold=1,
    )

    # Filter samples by signal.
    # Filter to keep only samples with signals beyond threshold in at least one
    # gene.
    data_detection_sample = measurement.filter_samples_by_signal_threshold(
        data=data_detection_gene,
        threshold=1,
    )

    # Return information.
    return data_detection_sample


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

    # Select samples by exclusion.
    # Exclude samples for extra-corporeal tissues, such as cell lines.
    # Exclude samples for tissues with fewer than 50 total samples.
    samples_exclusion = select_samples_by_exclusion(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
    )
    print(len(samples_exclusion))
    # Select samples by inclusion.
    samples_inclusion = select_samples_by_inclusion(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
    )
    print(len(samples_inclusion))
    # Select samples that are in both exclusion and inclusion lists.
    samples = utility.filter_common_elements(
        list_one=samples_exclusion,
        list_two=samples_inclusion,
    )
    utility.print_terminal_partition(level=2)
    print("Count of samples to select: " + str(len(samples)))

    # Select genes' signals for specific samples.
    data_gene_count_selection = select_samples_genes(
        samples=samples,
        data_gene_signal=source["data_gene_count"],
    )
    print(data_gene_count_selection)

    # Select genes' signals for specific samples.
    data_gene_signal_selection = select_samples_genes(
        samples=samples,
        data_gene_signal=source["data_gene_signal"],
    )
    print(data_gene_signal_selection)

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
