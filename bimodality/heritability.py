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

# Custom

import assembly
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality



##########
# Initialization


def initialize_directories(dock=None):
    """
    Initialize directories for procedure's product files.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = dock
    paths["heritability"] = os.path.join(paths["dock"], "heritability")
    # Remove previous files to avoid version or batch confusion.

    # Define paths for groups of persons.
    groups = list()
    groups.append("selection")
    groups.append("ventilation")
    for group in groups:
        paths[group] = dict()
        paths[group]["genes"] = os.path.join(
            paths["heritability"], group, "genes"
        )
        paths[group]["collection"] = os.path.join(
            paths["heritability"], group, "collection"
        )
        # Initialize directories.
        #utility.remove_directory(path=paths[group]["collection"])
        utility.create_directories(path=paths[group]["collection"])
    # Return information.
    return paths


# Source


def read_source_initial(
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
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_genes_selection = os.path.join(
        dock, "selection", "tight", "samples_genes_signals",
        "genes.pickle"
    )

    path_heritability_genes = os.path.join(
        dock, "heritability", group, "genes"
    )

    path_genes_unimodal = os.path.join(
        dock, "candidacy", group, "unimodal", "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        dock, "candidacy", group, "multimodal", "genes_multimodal.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    genes_heritability = utility.extract_subdirectory_names(
        path=path_heritability_genes
    )
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_selection": genes_selection,
        "genes_heritability": genes_heritability,
        "genes_unimodal": genes_unimodal,
        "genes_multimodal": genes_multimodal,
    }


# Collection


def collect_successful_genes(
    genes=None,
    path_genes=None,
):
    """
    Checks genes from split and heritability batch procedures.

    arguments:
        genes (list<str>): identifiers of genes
        path_genes (str): path to heritability genes directory

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
        path_gene = os.path.join(path_genes, gene)
        path_report = os.path.join(
            path_gene, "report.hsq"
        )
        # Determine whether heritability produced a successful report.
        if os.path.exists(path_report):
            successes.append(gene)
    # Return information.
    return successes


def check_genes(
    genes_selection=None,
    genes_heritability=None,
    path_genes=None,
):
    """
    Checks genes from split and heritability batch procedures.

    arguments:
        genes_selection (list<str>): identifiers of genes
        genes_heritability (list<str>): identifiers of genes
        path_genes (str): path to heritability genes directory

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print(
        "Compare lists of genes from selection and heritability " +
        "procedures."
    )

    print(
        "count of genes from selection (master): " + str(len(genes_selection))
    )
    print("count of genes from heritability: " + str(len(genes_heritability)))

    print(
        "Check whether all selection list genes are in heritability list."
    )
    # Returns True if all elements in list_two are in list_one.
    match = utility.compare_lists_by_inclusion(
        list_one=genes_selection,
        list_two=genes_heritability,
    )
    print("selection genes match heritability genes: " + str(match))
    utility.print_terminal_partition(level=2)

    # Determine counts of genes for which heritability analysis converged
    # successfully.
    genes_valid = collect_successful_genes(
        genes=genes_heritability,
        path_genes=path_genes,
    )
    print(
        "count of successful heritability genes: " +
        str(len(genes_valid))
    )
    utility.print_terminal_partition(level=2)

    pass


def read_collect_organize_genes_heritabilities(
    genes=None,
    data_gene_annotation=None,
    path_genes=None,
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
        data_gene_annotation (object): Pandas data frame of genes' annotations
        path_genes (str): path to heritability genes directory

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
        name = assembly.access_gene_name(
            identifier=gene,
            data_gene_annotation=data_gene_annotation,
        )
        # Specify directories and files.
        path_gene = os.path.join(path_genes, gene)
        path_report = os.path.join(path_gene, "report.hsq")
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
            genotype = data_report.at["V(G)", "Variance"]
            residual = data_report.at["V(e)", "Variance"]
            phenotype = data_report.at["Vp", "Variance"]
            proportion = data_report.at["V(G)/Vp", "Variance"]
            probability = data_report.at["Pval", "Variance"]
        else:
            duds.append(gene)
            genotype = float("nan")
            residual = float("nan")
            phenotype = float("nan")
            proportion = float("nan")
            probability = float("nan")
        # Compile information.
        record = dict()
        record["gene"] = gene
        record["name"] = name
        record["genotype"] = genotype
        record["residual"] = residual
        record["phenotype"] = phenotype
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
    data.set_index(
        "gene",
        drop=True,
        inplace=True,
    )
    # Return information.
    return data


def calculate_organize_discoveries(
    threshold=None,
    data_genes_heritabilities=None,
):
    """
    Calculates false discover rate.

    arguments:
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    def calculate_negative_logarithm(value=None):
        if math.isnan(value) or (value == 0):
            logarithm = float("nan")
        else:
            logarithm = (-1 * math.log(value, 10))
        return logarithm
    # Copy data.
    data_copy = data_genes_heritabilities.copy(deep=True)
    # Calculate false discovery rates.
    data_discoveries = utility.calculate_false_discovery_rate(
        threshold=threshold,
        probability="probability",
        discovery="discovery",
        significance="significance",
        data_probabilities=data_copy,
    )
    # Calculate logarithm of false discovery rate.
    data_discoveries["discovery_log"] = data_discoveries["discovery"].apply(
        lambda value: calculate_negative_logarithm(value=value)
    )
    return data_discoveries



# Segregation


def select_heritabilities_by_genes(
    genes=None,
    data_genes_heritabilities=None,
):
    """
    Selects information on the heritabilities of specific genes.

    arguments:
        genes (list<str>): identifiers of genes to select
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities
    """

    # Select information about samples.
    data_genes_heritabilities = data_genes_heritabilities.copy(deep=True)
    data_selection = data_genes_heritabilities.loc[
        data_genes_heritabilities.index.isin(genes), :
    ]
    # Return information.
    return data_selection


def select_heritable_genes_by_threshold(
    data_genes_heritabilities=None,
    threshold_proportion=None,
):
    """
    Collects and organizes information about genes.

    arguments:
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities
        threshold_proportion (float): threshold by proportion of phenotypic
            variance attributable to genotype

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    # Copy genes' heritabilities.
    data_copy = data_genes_heritabilities.copy(deep=True)
    # Set threshold.
    data_proportion = data_copy.loc[
        data_copy["proportion"] >= threshold_proportion
    ]
    data_discovery = data_proportion.loc[data_proportion["significance"], :]
    # Extract identifiers of genes.
    genes = data_discovery.index.to_list()
    # Return information.
    return genes


def select_genes_heritabilities_sets(
    genes=None,
    threshold_proportion=None,
    threshold_discovery=None,
    data_genes_heritabilities=None,
):
    """
    Selects information on the heritabilities of specific genes.

    arguments:
        genes (list<str>): identifiers of genes to select
        threshold_proportion (float): threshold by proportion of phenotypic
            variance attributable to genotype
        threshold_discovery (float): threshold by false discovery rate
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities

    raises:

    returns:
        (dict): information about genes' heritabilities
    """

    # Copy data.
    data_genes_heritabilities = data_genes_heritabilities.copy(deep=True)
    # Select data for genes of interest.
    data_genes = select_heritabilities_by_genes(
        genes=genes,
        data_genes_heritabilities=data_genes_heritabilities,
    )
    # Correct probabilities for false discovery rate across query genes.
    data_discovery = calculate_organize_discoveries(
        threshold=threshold_discovery,
        data_genes_heritabilities=data_genes,
    )
    # Determine genes with heritable distributions.
    utility.print_terminal_partition(level=2)
    print("Counts of heritable genes...")
    genes_heritable = select_heritable_genes_by_threshold(
        data_genes_heritabilities=data_discovery,
        threshold_proportion=threshold_proportion,
    )
    utility.print_terminal_partition(level=2)
    print("heritable genes: " + str(len(genes_heritable)))
    # Compile information.
    information = dict()
    information["data"] = data_discovery
    information["genes_heritable"] = genes_heritable
    # Return information.
    return information


# Product


def write_product(
    group=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        group (str): group of persons, either selection or ventilation
        information (object): information to write to file.
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_data_genes_heritability = os.path.join(
        paths[group]["collection"], "data_genes_heritability.pickle"
    )
    path_data_genes_heritability_text = os.path.join(
        paths[group]["collection"], "data_genes_heritability.tsv"
    )

    path_genes_selection = os.path.join(
        paths[group]["collection"], "genes_selection.pickle"
    )
    path_genes_selection_text = os.path.join(
        paths[group]["collection"], "genes_selection.txt"
    )
    path_genes_unimodal = os.path.join(
        paths[group]["collection"], "genes_unimodal.pickle"
    )
    path_genes_unimodal_text = os.path.join(
        paths[group]["collection"], "genes_unimodal.txt"
    )
    path_genes_multimodal = os.path.join(
        paths[group]["collection"], "genes_multimodal.pickle"
    )
    path_genes_multimodal_text = os.path.join(
        paths[group]["collection"], "genes_multimodal.txt"
    )

    # Write information to file.
    information["data_genes_heritability"].to_pickle(
        path=path_data_genes_heritability
    )
    information["data_genes_heritability"].to_csv(
        path_or_buf=path_data_genes_heritability_text,
        sep="\t",
        header=True,
        index=True,
    )

    with open(path_genes_selection, "wb") as file_product:
        pickle.dump(information["genes_selection"], file_product)
    utility.write_file_text_list(
        elements=information["genes_selection"],
        delimiter="\n",
        path_file=path_genes_selection_text
    )
    with open(path_genes_unimodal, "wb") as file_product:
        pickle.dump(information["genes_unimodal"], file_product)
    utility.write_file_text_list(
        elements=information["genes_unimodal"],
        delimiter="\n",
        path_file=path_genes_unimodal_text
    )
    with open(path_genes_multimodal, "wb") as file_product:
        pickle.dump(information["genes_multimodal"], file_product)
    utility.write_file_text_list(
        elements=information["genes_multimodal"],
        delimiter="\n",
        path_file=path_genes_multimodal_text
    )

    pass


###############################################################################
# Procedure


def collect_select_report_write_heritability_genes(
    group=None,
    paths=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        group (str): group of persons, either selection or ventilation
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Report.
    utility.print_terminal_partition(level=1)
    print("... Heritability procedure for: " + str(group) + " persons...")
    utility.print_terminal_partition(level=2)

    # Read source information from file.
    source = read_source_initial(
        group="selection",
        dock=paths["dock"],
    )

    # Verify that heritabilities are available for all genes.
    # simple method unavailable: ENSG00000022567
    # complex method unavailable: ENSG00000214860, ENSG00000180354
    check_genes(
        genes_selection=source["genes_selection"],
        genes_heritability=source["genes_heritability"],
        path_genes=paths[group]["genes"],
    )

    # Collect and organize information about genes' heritabilities.
    data_genes_heritabilities = read_collect_organize_genes_heritabilities(
        genes=source["genes_heritability"],
        data_gene_annotation=source["data_gene_annotation"],
        path_genes=paths[group]["genes"],
    )
    # Summary.
    utility.print_terminal_partition(level=2)
    print("data_genes_heritabilities")
    print(data_genes_heritabilities)

    # Segregate heritability reports for selection genes, unimodal genes, and
    # multimodal genes.
    utility.print_terminal_partition(level=2)
    print("Count of selection genes: " + str(len(source["genes_selection"])))
    print("Count of unimodal genes: " + str(len(source["genes_unimodal"])))
    print("Count of multimodal genes: " + str(len(source["genes_multimodal"])))
    utility.print_terminal_partition(level=2)

    # Select sets of heritable genes.
    # The false discover rate correction here is relative to each query set of
    # genes.
    bin_selection = select_genes_heritabilities_sets(
        genes=source["genes_selection"],
        threshold_proportion=0.1,
        threshold_discovery=0.05,
        data_genes_heritabilities=data_genes_heritabilities,
    )
    bin_unimodal = select_genes_heritabilities_sets(
        genes=source["genes_unimodal"],
        threshold_proportion=0.1,
        threshold_discovery=0.05,
        data_genes_heritabilities=data_genes_heritabilities,
    )
    utility.print_terminal_partition(level=3)
    print("group: " + group)
    print("multimodal genes")
    bin_multimodal = select_genes_heritabilities_sets(
        genes=source["genes_multimodal"],
        threshold_proportion=0.1,
        threshold_discovery=0.05,
        data_genes_heritabilities=data_genes_heritabilities,
    )

    # Compile information.
    information = {
        "data_genes_heritability": bin_selection["data"],
        "genes_selection": bin_selection["genes_heritable"],
        "genes_unimodal": bin_unimodal["genes_heritable"],
        "genes_multimodal": bin_multimodal["genes_heritable"],
    }
    #Write product information to file.
    write_product(
        group=group,
        information=information,
        paths=paths,
    )

    pass


def execute_procedure(
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Initialize directories.
    paths = initialize_directories(dock=dock)
    # Organize data and regress cases.
    collect_select_report_write_heritability_genes(
        group="selection",
        paths=paths,
    )
    collect_select_report_write_heritability_genes(
        group="ventilation",
        paths=paths,
    )
    pass



if (__name__ == "__main__"):
    execute_procedure()
