"""
...

"""

###############################################################################
# Notes

# If the count of bootstrap shuffles is too small, then it is likely not to
# have any values equal to or greater than the real value, in which case the
# probability is zero.

###############################################################################
# Installation and importation

# Standard

import os
import math
import statistics
import pickle
import copy
import random
import gc

# Relevant

import numpy
import pandas
import scipy.stats
import statsmodels

# Custom

import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


# Source


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
    path_split = os.path.join(dock, "split")
    path_genes = os.path.join(path_split, "genes.txt")
    path_heritability = os.path.join(dock, "heritability")

    # Read information from file.
    genes_split = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_genes,
    )
    genes_heritability = utility.extract_subdirectory_names(
        path=path_heritability
    )

    # Compile and return information.
    return {
        "genes_split": genes_split,
        "genes_heritability": genes_heritability,
    }


# Collection


def collect_successful_genes(
    genes=None,
    path_heritability=None,
    method=None,
):
    """
    Checks genes from split and heritability batch procedures.

    arguments:
        genes (list<str>): identifiers of genes
        path_heritability (str): path to heritability directory
        method (str): method from which to collect, "simple" or "complex"

    raises:

    returns:
        (list<str>): identifiers of genes for which heritability analysis was
            successful

    """

    # Collect information about genes.
    successes = list()
    # Iterate on genes.
    for gene in genes:
        # Specify directories and files.
        path_heritability_gene = os.path.join(path_heritability, gene)
        path_method = os.path.join(path_heritability_gene, method)
        path_report = os.path.join(
            path_method, "out.hsq"
        )
        # Determine whether heritability produced a successful report.
        if os.path.exists(path_report):
            successes.append(gene)
    # Return information.
    return successes


def check_genes(
    genes_split=None,
    genes_heritability=None,
    path_heritability=None,
):
    """
    Checks genes from split and heritability batch procedures.

    arguments:
        genes_split (list<str>): identifiers of genes
        genes_heritability (list<str>): identifiers of genes
        path_heritability (str): path to heritability directory

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print(
        "Compare lists of genes from split and heritability " +
        "procedures."
    )

    print("count of genes from split (master): " + str(len(genes_split)))
    print("count of genes from heritability: " + str(len(genes_heritability)))

    print(
        "Check whether all split list genes are in heritability list."
    )
    # Returns True if all elements in list_two are in list_one.
    match = utility.compare_lists_by_inclusion(
        list_one=genes_split,
        list_two=genes_heritability,
    )
    print("split genes match heritability genes: " + str(match))

    utility.print_terminal_partition(level=2)

    # Determine counts of genes for which heritability analysis converged
    # successfully.
    genes_simple = collect_successful_genes(
        genes=genes_heritability,
        path_heritability=path_heritability,
        method="simple",
    )
    print(
        "count of successful heritability genes (simple): " +
        str(len(genes_simple))
    )
    genes_complex = collect_successful_genes(
        genes=genes_heritability,
        path_heritability=path_heritability,
        method="complex",
    )
    print(
        "count of successful heritability genes (complex): " +
        str(len(genes_complex))
    )

    utility.print_terminal_partition(level=2)

    pass


def collect_organize_genes_heritabilities(
    genes=None,
    path_heritability=None,
    method=None,
):
    """
    Collects and organizes information about genes.

    Heritability analysis in GCTA does not converge successfully for all genes.

    Data structure.

    - genes_heritabilities (dict)
    -- gene (dict)
    --- proportion (float)
    --- probability (float)

    arguments:
        genes (list<str>): identifiers of genes
        path_heritability (str): path to heritability directory
        method (str): method from which to collect, "simple" or "complex"

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    # Collect genes' heritabilities.
    utility.print_terminal_partition(level=3)
    print("Collect genes' heritabilities.")
    # Collect information about genes.
    duds = list()
    records = list()
    # Iterate on genes.
    for gene in genes:
        # Specify directories and files.
        path_heritability_gene = os.path.join(path_heritability, gene)
        path_method = os.path.join(path_heritability_gene, method)
        path_report = os.path.join(
            path_method, "out.hsq"
        )
        # Determine whether heritability produced a successful report.
        if os.path.exists(path_report):
            # Read information from file.
            data_report = pandas.read_csv(
                path_report,
                sep="\t",
                header=0,
                low_memory=False,
            )
            # Organize data.
            data_report.set_index(
                "Source",
                drop=True,
                inplace=True,
            )
            # Access values.
            proportion = data_report.at["V(G)/Vp", "Variance"]
            probability = data_report.at["Pval", "Variance"]
        else:
            duds.append(gene)
            proportion = float("nan")
            probability = float("nan")
        # Compile information.
        record = dict()
        record["gene"] = gene
        record["proportion"] = proportion
        record["probability"] = probability
        records.append(record)

    utility.print_terminal_partition(level=3)
    print("collection complete")
    print("Duds: ")
    print(duds)
    # Organize data.
    data = utility.convert_records_to_dataframe(
        records=records,
    )
    # Return information.
    return data


def calculate_false_discovery_rate(
    data_genes_heritabilities=None,
):
    """
    Calculates false discover rate.

    arguments:
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    # Copy data.
    data_copy = data_genes_heritabilities.copy(deep=True)
    # False discovery rate method cannot accommodate missing values.
    # Remove null values.
    data_null_boolean = pandas.isna(data_copy["probability"])
    data_null = data_copy.loc[data_null_boolean]
    data_valid = data_copy.dropna(
        axis="index",
        how="any",
        inplace=False,
    )
    # Calculate false discovery rates from probabilities.
    probabilities = data_valid["probability"].to_numpy()
    report = statsmodels.stats.multitest.multipletests(
        probabilities,
        alpha=0.05,
        method="fdr_bh",
        is_sorted=False,
    )
    discoveries = report[1]
    data_valid["discovery"] = discoveries
    data_null["discovery"] = float("nan")
    # Combine null and valid portions of data.
    data_discoveries = data_valid.append(
        data_null,
        ignore_index=True,
    )

    return data_discoveries


def select_heritable_genes_by_threshold(
    data_genes_heritabilities=None,
    property=None,
    threshold=None,
):
    """
    Collects and organizes information about genes.

    arguments:
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities
        property (str): name of column by which to set threshold
        threshold (float): maximal value of property for selection

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    # Copy genes' heritabilities.
    data_copy = data_genes_heritabilities.copy(deep=True)

    # Set threshold.
    data_threshold = data_copy.loc[data_copy[property] < threshold]

    # Extract identifiers of genes.
    genes = data_threshold["gene"].to_list()

    # Return information.
    return genes


def convert_genes_heritabilities(
    data_genes_heritabilities=None,
):
    """
    Collects and organizes information about genes.

    arguments:
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities

    raises:

    returns:
        (dict<dict<float>>): information about genes' heritabilities

    """

    data_genes_heritabilities.reset_index(
        level=None,
        drop=True,
        inplace=True
    )
    data_genes_heritabilities.set_index(
        "gene",
        drop=True,
        inplace=True,
    )
    genes_heritabilities = data_genes_heritabilities.to_dict(
        "index"
    )
    return genes_heritabilities





# Product


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
    path_heritability = os.path.join(dock, "heritability")
    path_collection = os.path.join(path_heritability, "collection")
    utility.create_directory(path_collection)
    path_sets_genes_models = os.path.join(
        path_collection, "sets_genes_models.pickle"
    )
    path_simple = os.path.join(
        path_collection, "genes_heritabilities_simple.pickle"
    )
    path_complex = os.path.join(
        path_collection, "genes_heritabilities_complex.pickle"
    )
    path_data_simple = os.path.join(
        path_collection, "data_genes_heritabilities_simple.pickle"
    )
    path_data_complex = os.path.join(
        path_collection, "data_genes_heritabilities_complex.pickle"
    )
    path_genes_heritability = os.path.join(
        path_collection, "genes_heritability.pickle"
    )
    # Write information to file.
    with open(path_sets_genes_models, "wb") as file_product:
        pickle.dump(information["sets_genes_models"], file_product)
    with open(path_simple, "wb") as file_product:
        pickle.dump(
            information["genes_heritabilities_simple"], file_product
        )
    with open(path_complex, "wb") as file_product:
        pickle.dump(
            information["genes_heritabilities_complex"], file_product
        )
    with open(path_genes_heritability, "wb") as file_product:
        pickle.dump(
            information["genes_heritability"], file_product
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
    path_collection = os.path.join(dock, "heritability", "collection")
    utility.remove_directory(path=path_collection)

    # Read source information from file.
    source = read_source(dock=dock)

    # Verify that heritabilities are available for all genes.
    # simple method unavailable: ENSG00000022567
    # complex method unavailable: ENSG00000214860, ENSG00000180354
    path_heritability = os.path.join(dock, "heritability")
    check_genes(
        genes_split=source["genes_split"],
        genes_heritability=source["genes_heritability"],
        path_heritability=path_heritability,
    )

    # Collect and organize information about genes' heritabilities.
    data_genes_heritabilities_simple = collect_organize_genes_heritabilities(
        genes=source["genes_heritability"],
        path_heritability=path_heritability,
        method="simple",
    )
    data_genes_heritabilities_complex = collect_organize_genes_heritabilities(
        genes=source["genes_heritability"],
        path_heritability=path_heritability,
        method="complex",
    )

    # Correct probabilities for false discovery rate.
    data_genes_heritabilities_simple = calculate_false_discovery_rate(
        data_genes_heritabilities=data_genes_heritabilities_simple
    )
    data_genes_heritabilities_complex = calculate_false_discovery_rate(
        data_genes_heritabilities=data_genes_heritabilities_complex
    )

    # Summary.
    utility.print_terminal_partition(level=2)
    print("simple method:")
    print(data_genes_heritabilities_simple)
    utility.print_terminal_partition(level=2)
    print("complex method:")
    print(data_genes_heritabilities_complex)

    # Select genes with significant heritabilities.
    utility.print_terminal_partition(level=2)
    print("Counts of heritable genes...")
    genes_heritable_simple = select_heritable_genes_by_threshold(
        data_genes_heritabilities=data_genes_heritabilities_simple,
        property="discovery",
        threshold=0.05,
    )
    utility.print_terminal_partition(level=2)
    print("simple method: " + str(len(genes_heritable_simple)))
    genes_heritable_complex = select_heritable_genes_by_threshold(
        data_genes_heritabilities=data_genes_heritabilities_complex,
        property="discovery",
        threshold=0.05,
    )
    utility.print_terminal_partition(level=2)
    print("complex method: " + str(len(genes_heritable_complex)))
    sets_genes_models = dict()
    sets_genes_models["simple"] = genes_heritable_simple
    sets_genes_models["complex"] = genes_heritable_complex
    #print(sets_genes_models)

    # Convert data to dictionaries.
    # Organize data.
    genes_heritabilities_simple = convert_genes_heritabilities(
        data_genes_heritabilities=data_genes_heritabilities_simple,
    )
    genes_heritabilities_complex = convert_genes_heritabilities(
        data_genes_heritabilities=data_genes_heritabilities_complex,
    )

    utility.print_terminal_partition(level=2)
    print(genes_heritabilities_simple["ENSG00000161634"])
    print(genes_heritabilities_complex["ENSG00000161634"])

    # Select genes that pass filters by multiple measures of bimodality.
    genes_heritability = utility.select_elements_by_sets(
        names=["simple", "complex"],
        sets=sets_genes_models,
        count=1,
    )
    utility.print_terminal_partition(level=2)
    print(
        "selection of genes by multiple measurements: " +
        str(len(genes_heritability))
    )
    utility.print_terminal_partition(level=2)


    # Compile information.
    information = {
        "genes_heritabilities_simple": genes_heritabilities_simple,
        "genes_heritabilities_complex": genes_heritabilities_complex,
        "sets_genes_models": sets_genes_models,
        "genes_heritability": genes_heritability,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
