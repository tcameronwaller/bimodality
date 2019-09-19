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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import assembly
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
    path_selection = os.path.join(dock, "selection")
    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_signal": data_gene_signal,
    }


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
    data_transposition = data_gene_signal.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_sample=data_transposition,
    )
    return data_factor


def split_genes_signals(
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

    # Split genes' signals across tissues and patients by gene.
    utility.print_terminal_partition(level=2)
    print(
        "Split of genes' signals across samples by gene."
    )

    # Optimize data types.
    data_gene_signal = assembly.convert_data_types(
        data=data_gene_signal,
        type="float32"
    )

    # Associate samples to factors.
    utility.print_terminal_partition(level=2)
    print(
        "Association of samples to factors."
    )
    data_factor = associate_factors(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_signal=data_gene_signal,
    )
    print(data_factor)
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

    # TODO: Wait... I don't need to stack genes or groupby genes...
    # I just need to split by gene columns...
    # for column in columns...



    # Stack by gene.
    utility.print_terminal_partition(level=2)
    print("Stack of genes to factors.")
    #data_long = data_factor.stack("gene").to_frame(name="signal")
    signals_long = data_factor.stack("gene")
    data_long = signals_long.to_frame(name="signal")
    # Organize data.
    data_long.reset_index(
        level=None,
        inplace=True
    )
    # Convert data types.
    data_type = data_long.astype(
        {
            "gene": "str",
            "person": "str",
            "tissue_major": "str",
            "signal": "float32",
        },
        copy=True,
    )
    data_type.set_index(
        ["gene", "person", "tissue_major"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data_type)

    # Split by gene.
    utility.print_terminal_partition(level=2)
    print("Split of genes' signals across samples by gene.")
    groups = data_type.groupby("gene")
    print("Count of groups by genes: " + str(len(groups)))
    # Collect data frames by gene.
    genes_samples_signals = dict()
    for name, group in groups:
        data = group.copy(deep=True)
        genes_samples_signals[name] = data
        if name == "ENSG00000231925":
            print(name)
            print(group)
    # Access data by gene.
    print("Access data for a single gene.")
    print(genes_samples_signals["ENSG00000231925"])

    # Return information.
    return genes_samples_signals


def summarize_gene_samples_signals(
    data_gene_samples_signals=None,
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

    print(data_gene_samples_signals)
    utility.print_terminal_partition(level=2)
    print("Determine counts of persons and tissues.")
    print("Split gene's signals by person.")
    groups = data_gene_samples_signals.groupby("person")
    persons = len(groups)
    print("Count of groups by person: " + str(persons))
    # Count of persons = 710
    print("Split gene's signals by major tissue category.")
    groups = data_gene_samples_signals.groupby("tissue_major")
    tissues = len(groups)
    print("Count of groups by tissue: " + str(tissues))
    # Count of tissues = 20
    # Compile information.
    information = {
        "persons": persons,
        "tissues": tissues,
    }
    # Return information.
    return information


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
    path_split = os.path.join(dock, "split")
    utility.create_directory(path_split)
    path_genes = os.path.join(
        path_split, "genes.txt"
    )
    path_persons = os.path.join(
        path_split, "persons.pickle"
    )
    path_tissues = os.path.join(
        path_split, "tissues.pickle"
    )
    path_signal = os.path.join(
        path_split, "genes_samples_signals.pickle"
    )
    path_demo = os.path.join(
        path_split, "demonstration_gene_samples_signals.pickle"
    )
    # Write information to file.

    # List of genes needs to be easy to read in Bash.
    utility.write_file_text_list(
        elements=information["genes"],
        delimiter="\n",
        path_file=path_genes
    )
    with open(path_persons, "wb") as file_product:
        pickle.dump(information["persons"], file_product)
    with open(path_tissues, "wb") as file_product:
        pickle.dump(information["tissues"], file_product)
    with open(path_signal, "wb") as file_product:
        pickle.dump(information["genes_samples_signals"], file_product)
    pandas.to_pickle(
        information["data_gene_samples_signals"],
        path_demo
    )

    # Write separate files for genes.
    # Specify directories and files.
    path_collection = os.path.join(path_split, "collection")
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

    # Read source information from file.
    source = read_source(dock=dock)

    # Summary.
    print(source["data_gene_signal"].iloc[0:10, 0:10])
    print(source["data_gene_signal"].shape)

    # Split genes' signals across tissues and patients by gene.
    genes_samples_signals = split_genes_signals(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=source["data_gene_signal"]
    )

    # Organize genes' identifiers.
    # Format of genes' identifiers needs to be readable by Bash as an array.
    genes = source["data_gene_signal"].index.to_list()
    print("count of genes: " + str(len(genes)))

    # Summarize information for a single gene.
    # Access data for single gene for demonstration.
    data_gene_samples_signals = genes_samples_signals["ENSG00000231925"]
    persons_tissues = summarize_gene_samples_signals(
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Write the entire collection of all genes' signals to a single file.
    # Also write each gene's signals to a separate file.
    # Conserve memory in parallel pipeline by reading data for each gene
    # separately.

    # Compile information.
    information = {
        "genes": genes,
        "persons": persons_tissues["persons"],
        "tissues": persons_tissues["tissues"],
        "genes_samples_signals": genes_samples_signals,
        "data_gene_samples_signals": data_gene_samples_signals,
    }

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
