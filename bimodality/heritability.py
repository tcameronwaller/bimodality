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


def read_source_initial(dock=None):
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
    path_genes_selection = os.path.join(
        path_selection, "genes_selection.pickle"
    )

    path_heritability = os.path.join(dock, "heritability")
    path_heritability_genes = os.path.join(path_heritability, "genes")

    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_unimodal = os.path.join(
        path_candidacy, "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )


    # Read information from file.
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
        "genes_selection": genes_selection,
        "genes_heritability": genes_heritability,
        "genes_unimodal": genes_unimodal,
        "genes_multimodal": genes_multimodal,
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
        path_heritability_gene = os.path.join(path_heritability, "genes", gene)
        path_method = os.path.join(path_heritability_gene, method)
        path_report = os.path.join(
            path_method, "report.hsq"
        )
        # Determine whether heritability produced a successful report.
        if os.path.exists(path_report):
            successes.append(gene)
    # Return information.
    return successes


def check_genes(
    genes_selection=None,
    genes_heritability=None,
    path_heritability=None,
):
    """
    Checks genes from split and heritability batch procedures.

    arguments:
        genes_selection (list<str>): identifiers of genes
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


def read_collect_organize_genes_heritabilities(
    genes=None,
    method=None,
    dock=None,
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
        method (str): method from which to collect, "simple" or "complex"
        dock (str): path to root or dock directory for source and product
            directories and files


    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    path_heritability = os.path.join(dock, "heritability")

    # Collect genes' heritabilities.
    utility.print_terminal_partition(level=3)
    print("Collect genes' heritabilities.")
    # Collect information about genes.
    duds = list()
    records = list()
    # Iterate on genes.
    for gene in genes:
        # Specify directories and files.
        path_heritability_gene = os.path.join(path_heritability, "genes", gene)
        path_method = os.path.join(path_heritability_gene, method)
        path_report = os.path.join(
            path_method, "report.hsq"
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


def calculate_false_discovery_rate(
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
    data_copy.reset_index(
        level=None,
        inplace=True
    )

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
        alpha=threshold,
        method="fdr_bh",
        is_sorted=False,
    )
    significances = report[0]
    discoveries = report[1]
    data_valid["significance"] = significances
    data_valid["discovery"] = discoveries
    data_null["significance"] = False
    data_null["discovery"] = float("nan")

    # Combine null and valid portions of data.
    data_discoveries = data_valid.append(
        data_null,
        ignore_index=True,
    )

    # Organize data.
    data_discoveries.set_index(
        "gene",
        drop=True,
        inplace=True,
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
    threshold_discovery=None,
):
    """
    Collects and organizes information about genes.

    arguments:
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities
        threshold_proportion (float): threshold by proportion of phenotypic
            variance attributable to genotype
        threshold_discovery (float): threshold by false discovery rate

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
    data_discovery = data_proportion.loc[
        data_proportion["discovery"] <= threshold_discovery
    ]

    # Extract identifiers of genes.
    genes = data_discovery["gene"].to_list()

    # Return information.
    return genes




def select_genes_heritabilities_sets(
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

    pass



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
    path_data_genes_heritabilities_simple = os.path.join(
        path_collection, "data_genes_heritabilities_simple.pickle"
    )
    path_data_genes_heritabilities_complex = os.path.join(
        path_collection, "data_genes_heritabilities_complex.pickle"
    )
    path_data_genes_heritabilities_complex_text = os.path.join(
        path_collection, "data_genes_heritabilities_complex.tsv"
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
    information["data_genes_heritabilities_simple"].to_pickle(
        path=path_data_genes_heritabilities_simple
    )
    information["data_genes_heritabilities_complex"].to_pickle(
        path=path_data_genes_heritabilities_complex
    )
    information["data_genes_heritabilities_complex"].to_csv(
        path_or_buf=path_data_genes_heritabilities_complex_text,
        sep="\t",
        header=True,
        index=True,
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
    source = read_source_initial(dock=dock)

    # Verify that heritabilities are available for all genes.
    # simple method unavailable: ENSG00000022567
    # complex method unavailable: ENSG00000214860, ENSG00000180354
    path_heritability = os.path.join(dock, "heritability")
    check_genes(
        genes_selection=source["genes_selection"],
        genes_heritability=source["genes_heritability"],
        path_heritability=path_heritability,
    )

    # Collect and organize information about genes' heritabilities.
    data_genes_heritabilities_simple = (
        read_collect_organize_genes_heritabilities(
            genes=source["genes_heritability"],
            method="simple",
            dock=dock,
    ))
    data_genes_heritabilities_complex = (
        read_collect_organize_genes_heritabilities(
            genes=source["genes_heritability"],
            method="complex",
            dock=dock,
    ))

    # Correct probabilities for false discovery rate.
    utility.print_terminal_partition(level=2)
    data_genes_heritabilities_simple = calculate_false_discovery_rate(
        threshold=0.05,
        data_genes_heritabilities=data_genes_heritabilities_simple
    )
    data_genes_heritabilities_complex = calculate_false_discovery_rate(
        threshold=0.05,
        data_genes_heritabilities=data_genes_heritabilities_complex
    )

    # Summary.
    utility.print_terminal_partition(level=2)
    print("simple method:")
    print(data_genes_heritabilities_simple)
    utility.print_terminal_partition(level=2)
    print("complex method:")
    print(data_genes_heritabilities_complex)

    # Select sets of heritable genes.
    bin_genes_selection = select_genes_heritabilities_sets(
        genes=source["genes_selection"],
        data_genes_heritabilities_simple=data_genes_heritabilities_simple,
        data_genes_heritabilities_complex=data_genes_heritabilities_complex,
    )
    bin_genes_unimodal = select_genes_heritabilities_sets(
        genes=source["genes_unimodal"],
        data_genes_heritabilities_simple=data_genes_heritabilities_simple,
        data_genes_heritabilities_complex=data_genes_heritabilities_complex,
    )
    bin_genes_multimodal = select_genes_heritabilities_sets(
        genes=source["genes_multimodal"],
        data_genes_heritabilities_simple=data_genes_heritabilities_simple,
        data_genes_heritabilities_complex=data_genes_heritabilities_complex,
    )

    if False:
        # TODO: make this cleaner by generating a driver function
        # pass the function
        # 1. genes of interest
        # 2. data_genes_heritabilities_simple
        # 3. data_genes_heritabilities_complex
        # return from function a dictionary...
        # 1. data_genes_selection_heritabilities_simple
        # 2. identifiers of the heritable genes organized in the sets format for charting...

        # Segregate heritability reports for selection genes, unimodal genes, and
        # multimodal genes.
        print("Count of unimodal genes: " + str(len(source["genes_unimodal"])))
        print("Count of multimodal genes: " + str(len(source["genes_multimodal"])))
        data_genes_unimodal_heritabilities_simple = select_heritabilities_by_genes(
            genes=source["genes_unimodal"],
            data_genes_heritabilities=data_genes_heritabilities_simple,
        )
        data_genes_unimodal_heritabilities_complex = select_heritabilities_by_genes(
            genes=source["genes_unimodal"],
            data_genes_heritabilities=data_genes_heritabilities_complex,
        )
        data_genes_multimodal_heritabilities_simple = select_heritabilities_by_genes(
            genes=source["genes_multimodal"],
            data_genes_heritabilities=data_genes_heritabilities_simple,
        )
        data_genes_multimodal_heritabilities_complex = select_heritabilities_by_genes(
            genes=source["genes_multimodal"],
            data_genes_heritabilities=data_genes_heritabilities_complex,
        )

        # TODO: now select "heritable" genes by thresholds on FDR and proportion of variance explained (PVE)


        # TODO: I need sets of genes that are heritable for...
        # 1. unimodal simple
        # 2. unimodal complex
        # 3. multimodal simple
        # 4. multimodal complex
        # How do unimodal genes differ from multimodal genes on the basis of their heritability in simple and complex models?


        # TODO: create set overlap charts (simple versus complex) for selection genes, unimodal, and multimodal


    if False:
        # Select genes with significant heritabilities.
        utility.print_terminal_partition(level=2)
        print("Counts of heritable genes...")
        genes_heritable_simple = select_heritable_genes_by_threshold(
            data_genes_heritabilities=data_genes_heritabilities_simple,
            threshold_proportion=0.5,
            threshold_discovery=0.05,
        )
        utility.print_terminal_partition(level=2)
        print("simple method: " + str(len(genes_heritable_simple)))
        genes_heritable_complex = select_heritable_genes_by_threshold(
            data_genes_heritabilities=data_genes_heritabilities_complex,
            threshold_proportion=0.5,
            threshold_discovery=0.05,
        )
        utility.print_terminal_partition(level=2)
        print("complex method: " + str(len(genes_heritable_complex)))
        sets_genes_models = dict()
        sets_genes_models["simple"] = genes_heritable_simple
        sets_genes_models["complex"] = genes_heritable_complex
        #print(sets_genes_models)

        pass


    # Compile information.
    information = {
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

        #"genes_heritabilities_simple": genes_heritabilities_simple,
        #"genes_heritabilities_complex": genes_heritabilities_complex,
        #"sets_genes_models": sets_genes_models,
        #"genes_heritability": genes_heritable_complex, # can also use "genes_heritability" which is union
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
