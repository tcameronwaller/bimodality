"""
...

"""

###############################################################################
# Notes

# TODO: move this procedure to the organization procedure... ?

###############################################################################
# Installation and importation

# Standard

import os
import math
import statistics
import pickle
import json

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import assembly
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
    path_data_samples_tissues_persons = os.path.join(
        dock, "selection", "tight", "samples_genes_signals",
        "data_samples_tissues_persons.pickle"
    )
    path_data_gene_signal = os.path.join(
        dock, "selection", "tight", "samples_genes_signals",
        "data_gene_signal.pickle"
    )
    path_persons_sets = os.path.join(
        dock, "selection", "tight", "persons_properties",
        "persons_sets.pickle",
    )
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_data_samples_tissues_persons
    )
    data_gene_signal = pandas.read_pickle(path_data_gene_signal)
    with open(path_persons_sets, "rb") as file_source:
        persons_sets = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_signal": data_gene_signal,
        "persons_sets": persons_sets,
    }


def select_samples_signals_persons(
    persons=None,
    data_samples_tissues_persons=None,
    data_gene_signal=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        persons (list<str>): identifiers of persons
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues across samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_signals = data_gene_signal.copy(deep=True)
    # Select samples.
    data_samples_persons = data_samples.loc[
        data_samples["person"].isin(persons), :
    ]
    samples = data_samples_persons.index.to_list()
    # Select signals.
    data_signal_persons = data_signals.loc[
        :, data_signals.columns.isin(samples)
    ]
    # Collect information.
    bin = dict()
    bin["data_samples_tissues_persons"] = data_samples_persons
    bin["data_gene_signal"] = data_signal_persons
    # Return information.
    return bin


def associate_factors(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Splits genes' signals across samples by gene.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' signals
            across samples

    """

    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_signal_transposition = data_gene_signal.transpose(copy=True)
    data_signal_transposition.reset_index(
        level=None,
        inplace=True
    )
    data_signal_transposition.set_index(
        ["sample"],
        append=False,
        drop=True,
        inplace=True
    )
    # Associate samples to persons and tissues.
    data_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_sample=data_signal_transposition,
    )
    # Return information.
    return data_factor


def split_genes_signals(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Splits genes' signals across samples by gene.

    Splitting the data by gene is more efficient for parallel processing in
    batches by gene. Each process only needs to read in data for its gene.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' signals
            across samples

    """

    # Copy data.
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    # Associate samples to factors.
    data_factor = associate_factors(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_signal=data_gene_signal,
    )
    # Organize data.
    data_factor.reset_index(
        level=None,
        inplace=True
    )
    data_factor.drop(
        labels=["sample", "tissue_minor"],
        axis="columns",
        inplace=True
    )
    data_factor = data_factor.astype(
        {
            "person": "str",
            "tissue_major": "str",
        },
        copy=True,
    )
    data_factor.set_index(
        ["person", "tissue_major"],
        append=False,
        drop=True,
        inplace=True
    )
    # Optimize data types.
    data_factor = assembly.convert_data_types(
        data=data_factor,
        type="float32"
    )

    # Split data by gene.
    # An alternative option is to split data by groups.
    #data_long = data_factor.stack("gene").to_frame(name="signal")
    #groups = data_type.groupby("gene")
    #for name, group in groups:
    #    pass
    #This option is inefficient and overloads memory.
    # Instead, slice the data by columns for genes.

    # Collect data frames by gene.
    genes_samples_signals = dict()
    for gene in data_factor.columns.to_list():
        data_gene = data_factor[gene].to_frame(name="signal")
        genes_samples_signals[gene] = data_gene
    # Return information.
    return genes_samples_signals


def summarize_genes_samples_signals(
    genes_samples_signals=None,
):
    """
    Summarize information about a gene's samples and signals.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (dict): counts of persons and tissues

    """

    # Report.
    utility.print_terminal_partition(level=2)
    print("Count of data by genes: " + str(len(genes_samples_signals.keys())))
    print("Access data for a single gene.")
    utility.print_terminal_partition(level=2)
    data = genes_samples_signals["ENSG00000231925"]
    print(data)
    utility.print_terminal_partition(level=2)
    print("Determine counts of persons and tissues.")
    print("Split gene's signals by person.")
    groups = data.groupby("person")
    persons = len(groups)
    print("Count of groups by person: " + str(persons))
    print("Split gene's signals by major tissue category.")
    groups = data.groupby("tissue_major")
    tissues = len(groups)
    print("Count of groups by tissue: " + str(tissues))
    pass


def split_report_write_genes_signals(
    cohort=None,
    persons=None,
    data_samples_tissues_persons=None,
    data_gene_signal=None,
    path_directory=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        persons (list<str>): identifiers of persons
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues across samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples
        path_directory (str): path to directory for product directories and
            files
        report (bool): whether to print reports about the selection

    raises:

    returns:

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=1)
        print("... Split procedure for: " + str(cohort) + " persons...")
        print("Count persons: " + str(len(persons)))
        utility.print_terminal_partition(level=2)

    # Copy data.
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    data_gene_signal = data_gene_signal.copy(deep=True)
    # Select samples for relevant persons.
    bin = select_samples_signals_persons(
        persons=persons,
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_signal=data_gene_signal,
    )
    # Split genes' signals across tissues and patients by gene.
    genes_samples_signals = split_genes_signals(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_signal=bin["data_gene_signal"],
    )

    # Organize genes' identifiers.
    # Format of genes' identifiers needs to be readable by Bash as an array.
    genes = utility.collect_unique_elements(
        elements_original=list(genes_samples_signals.keys())
    )

    # Summarize information for a single gene.
    # Access data for single gene for demonstration.
    if report:
        summarize_genes_samples_signals(
            genes_samples_signals=genes_samples_signals,
        )
    # Write the entire collection of all genes' signals to a single file.
    # Also write each gene's signals to a separate file.
    # Conserve memory in parallel pipeline by reading data for each gene
    # separately.
    # Compile information.
    information = {
        "genes": genes,
        "genes_samples_signals": genes_samples_signals,
    }
    # Write product information to file.
    write_product(
        information=information,
        path_directory=path_directory,
    )
    pass


def write_product(
    information=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        information (dict): information to write to file
        path_directory (str): path to directory for product directories and
            files

    raises:

    returns:

    """

    # Specify directories and files.
    utility.create_directories(path_directory)
    path_genes = os.path.join(
        path_directory, "genes.pickle"
    )
    path_genes_samples_signals = os.path.join(
        path_directory, "genes_samples_signals.pickle"
    )
    # Write information to file.
    with open(path_genes, "wb") as file_product:
        pickle.dump(information["genes"], file_product)
    with open(path_genes_samples_signals, "wb") as file_product:
        pickle.dump(information["genes_samples_signals"], file_product)

    # Write separate files for genes.
    # Specify directories and files.
    path_collection = os.path.join(path_directory, "collection")
    utility.create_directory(path_collection)
    # Iterate on genes.
    for gene in information["genes"]:
        # Specify directories and files.
        path_gene = os.path.join(path_collection, (gene + ".pickle"))
        # Write information to file.
        pandas.to_pickle(
            information["genes_samples_signals"][gene],
            path_gene
        )
        pass
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
    path_split = os.path.join(dock, "split")
    utility.remove_directory(path=path_split)

    # Read source information from file.
    source = read_source(dock=dock)

    split_report_write_genes_signals(
        cohort="selection",
        persons=source["persons_sets"]["selection"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=source["data_gene_signal"],
        path_directory=os.path.join(dock, "split", "selection"),
        report=True,
    )
    split_report_write_genes_signals(
        cohort="respiration",
        persons=source["persons_sets"]["respiration"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=source["data_gene_signal"],
        path_directory=os.path.join(dock, "split", "respiration"),
        report=True,
    )
    split_report_write_genes_signals(
        cohort="ventilation",
        persons=source["persons_sets"]["ventilation"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=source["data_gene_signal"],
        path_directory=os.path.join(dock, "split", "ventilation"),
        report=True,
    )

    pass


if (__name__ == "__main__"):
    execute_procedure()
