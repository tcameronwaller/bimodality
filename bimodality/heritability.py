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
    cohorts = list()
    cohorts.append("selection")
    cohorts.append("respiration")
    cohorts.append("ventilation")
    for cohort in cohorts:
        paths[cohort] = dict()
        paths[cohort]["genes"] = os.path.join(
            paths["heritability"], cohort, "genes"
        )
        paths[cohort]["collection"] = os.path.join(
            paths["heritability"], cohort, "collection"
        )
        # Initialize directories.
        utility.remove_directory(path=paths[cohort]["collection"])
        utility.create_directories(path=paths[cohort]["collection"])
    # Return information.
    return paths


# Source


def read_source_initial(
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
    path_genes_selection = os.path.join(
        dock, "selection", "tight", "samples_genes_signals",
        "genes.pickle"
    )
    path_genes_distribution = os.path.join(
        dock, "distribution", cohort, "genes"
    )
    path_genes_heritability_complete = os.path.join(
        dock, "heritability", cohort, "genes"
    )

    path_genes_unimodal = os.path.join(
        dock, "candidacy", cohort, "unimodal", "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        dock, "candidacy", cohort, "multimodal", "genes_multimodal.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    genes_distribution = utility.extract_subdirectory_names(
        path=path_genes_distribution
    )
    genes_heritability_complete = utility.extract_subdirectory_names(
        path=path_genes_heritability_complete
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_selection": genes_selection,
        "genes_distribution": genes_distribution,
        "genes_heritability_complete": genes_heritability_complete,
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
    genes_distribution=None,
    genes_heritability_complete=None,
    path_genes=None,
):
    """
    Checks genes from split and heritability batch procedures.

    arguments:
        genes_selection (list<str>): identifiers of genes from selection
            procedure
        genes_distribution (list<str>): identifiers of genes with valid
            pan-tissue signal distributions
        genes_heritability_complete (list<str>): identifiers of genes for which
            the heritability procedure completed
        path_genes (str): path to heritability genes directory

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print(
        "Compare lists of genes from distribution and heritability " +
        "procedures."
    )

    print(
        "count of genes from selection: " + str(len(genes_selection))
    )
    print(
        "count of genes from distribution (master): " +
        str(len(genes_distribution))
    )
    print(
        "count of genes from heritability: " +
        str(len(genes_heritability_complete))
    )

    print(
        "Check whether all selection list genes are in heritability list."
    )
    # Returns True if all elements in list_two are in list_one.
    match = utility.compare_lists_by_inclusion(
        list_one=genes_distribution,
        list_two=genes_heritability_complete,
    )
    print("distribution genes match heritability genes: " + str(match))
    utility.print_terminal_partition(level=2)

    # Determine counts of genes for which heritability analysis converged
    # successfully.
    genes_valid = collect_successful_genes(
        genes=genes_heritability_complete,
        path_genes=path_genes,
    )
    print(
        "count of successful heritability genes: " +
        str(len(genes_valid))
    )
    utility.print_terminal_partition(level=2)

    pass


def find_intersection_heritability_genes(
    genes_selection=None,
    genes_distribution=None,
    genes_heritability_complete=None,
    path_genes=None,
    report=None,
):
    """
    Reads and organizes source information from file

    arguments:
        genes_selection (list<str>): identifiers of genes from selection
            procedure
        genes_distribution (list<str>): identifiers of genes with valid
            pan-tissue signal distributions
        genes_heritability_complete (list<str>): identifiers of genes for which
            the heritability procedure completed
        path_genes (str): path to heritability genes directory
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of genes of interest from selection that also
            have valid heritability measurements

    """

    # Determine genes for which heritability analysis converged successfully.
    genes_heritability_valid = collect_successful_genes(
        genes=genes_heritability_complete,
        path_genes=path_genes,
    )
    # Determine intersection genes of interest.
    genes_interest = utility.filter_common_elements(
        list_one=genes_selection,
        list_two=genes_heritability_valid,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "genes of interest with valid heritabilities: " +
            str(len(genes_interest))
        )
        utility.print_terminal_partition(level=2)
    # Return information.
    return genes_interest


def read_collect_organize_genes_heritabilities(
    genes=None,
    data_gene_annotation=None,
    path_genes=None,
    report=None,
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
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

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
            error = data_report.at["V(G)/Vp", "SE"]
            confidence = (1.96 * error)
            confidence_low = proportion - confidence
            confidence_high = proportion + confidence
            probability = data_report.at["Pval", "Variance"]
            count = data_report.at["n", "Variance"]
        else:
            duds.append(gene)
            genotype = float("nan")
            residual = float("nan")
            phenotype = float("nan")
            proportion = float("nan")
            error = float("nan")
            confidence = float("nan")
            confidence_low = float("nan")
            confidence_high = float("nan")
            probability = float("nan")
            count = float("nan")
        # Compile information.
        record = dict()
        record["gene"] = gene
        record["name"] = name
        record["genotype"] = genotype
        record["residual"] = residual
        record["phenotype"] = phenotype
        record["proportion"] = proportion
        record["error"] = error
        record["confidence_95_interval"] = confidence
        record["confidence_95_low"] = confidence_low
        record["confidence_95_high"] = confidence_high
        record["probability"] = probability
        record["count"] = count
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

    # Report.
    if report:
        utility.print_terminal_partition(level=3)
        print("Collect genes' heritabilities.")
        utility.print_terminal_partition(level=2)
        print("data_genes_heritabilities")
        print(data)
    # Return information.
    return data


# Filter


def determine_confidence_threshold_pass(
    value=None,
    threshold=None,
):
    """
    Determine whether confidence interval is within threshold.

    arguments:
        value (float): value of half the range of the 95% confidence interval
        threshold (float): maximal value for half the range of the 95%
            confidence interval

    raises:

    returns:
        (bool): whether confidence interval

    """

    if not math.isnan(value):
        if (value <= threshold):
            return True
        else:
            return False
    else:
        return False


def filter_heritabilities_confidence(
    data_genes_heritability=None,
    threshold=None,
    report=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        data_genes_heritability (object): Pandas data frame of genes'
            heritabilities
        threshold (float): maximal confidence interval
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    # Remove all columns from persons properties except the covariates
    # Copy data.
    data = data_genes_heritability.copy(deep=True)
    # Organize data.
    data["threshold"] = data["confidence_95_interval"].apply(
        lambda value:
            determine_confidence_threshold_pass(
                value=value,
                threshold=threshold,
            )
    )
    data_confidence = data.loc[
        data["threshold"] == True, :
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("data after filter by confidence interval")
        print(data_confidence)
        utility.print_terminal_partition(level=2)
        print("count of candidate genes': " + str(data_confidence.shape[0]))

    # Return information.
    return data_confidence


def calculate_organize_discoveries(
    threshold=None,
    data_genes_heritability=None,
):
    """
    Calculates false discover rate.

    arguments:
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries
        data_genes_heritability (object): Pandas data frame of genes'
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
    data_copy = data_genes_heritability.copy(deep=True)
    # Calculate false discovery rates.
    data_discovery = utility.calculate_false_discovery_rate(
        threshold=threshold,
        probability="probability",
        discovery="discovery",
        significance="significance",
        data_probabilities=data_copy,
    )
    # Calculate logarithm of false discovery rate.
    data_discovery["probability_log"] = data_discovery["probability"].apply(
        lambda value: calculate_negative_logarithm(value=value)
    )
    data_discovery["discovery_log"] = data_discovery["discovery"].apply(
        lambda value: calculate_negative_logarithm(value=value)
    )
    return data_discovery


def organize_genes_heritability_data(
    data_genes_heritability=None,
    report=None,
):
    """
    Organize data summarizing genes' heritabilities.

    arguments:
        data_genes_heritability (object): Pandas data frame of genes'
            heritabilities
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of genes' heritabilities

    """

    # Copy data.
    data = data_genes_heritability.copy(deep=True)
    columns = list()
    columns.append("name")
    columns.append("proportion")
    columns.append("count")
    columns.append("probability")
    columns.append("probability_log")
    columns.append("discovery")
    columns.append("discovery_log")
    columns.append("significance")
    columns.append("error")
    columns.append("confidence_95_interval")
    columns.append("confidence_95_low")
    columns.append("confidence_95_high")
    columns.append("residual")
    columns.append("genotype")
    columns.append("phenotype")
    data = data[[*columns]]
    data.sort_values(
        by=["probability"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("data after organization of columns")
        print(data)
    return data


def select_heritable_genes(
    data_genes_heritability=None,
    threshold_proportion=None,
    threshold_probability=None,
    report=None,
):
    """
    Collects and organizes information about genes.

    arguments:
        data_genes_heritability (object): Pandas data frame of genes'
            heritabilities
        threshold_proportion (float): threshold by proportion of phenotypic
            variance attributable to genotype
        threshold_probability (float): threshold by probability of heritability
            estimate
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of heritable genes

    """

    # Copy genes' heritabilities.
    data_copy = data_genes_heritability.copy(deep=True)
    # Set threshold.
    data_proportion = data_copy.loc[
        data_copy["proportion"] >= threshold_proportion
    ]
    data_probability = data_proportion.loc[
        data_proportion["probability"] <= threshold_probability
    ]
    # Extract identifiers of genes.
    genes = data_probability.index.to_list()
    # Report.
    if report:
        percentage = round((len(genes) / data_copy.shape[0]) * 100, 2)
        utility.print_terminal_partition(level=2)
        print(
            "count of 'heritable' genes': " +
            str(len(genes)) + " (" + str(percentage) + " %)"
        )
    # Return information.
    return genes


# Product


def write_product(
    cohort=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        information (object): information to write to file.
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_data_genes_heritability = os.path.join(
        paths[cohort]["collection"], "data_genes_heritability.pickle"
    )
    path_data_genes_heritability_text = os.path.join(
        paths[cohort]["collection"], "data_genes_heritability.tsv"
    )
    path_genes_heritability = os.path.join(
        paths[cohort]["collection"], "genes.pickle"
    )
    path_genes_heritability_text = os.path.join(
        paths[cohort]["collection"], "genes.txt"
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

    with open(path_genes_heritability, "wb") as file_product:
        pickle.dump(information["genes_heritability"], file_product)
    utility.write_file_text_list(
        elements=information["genes_heritability"],
        delimiter="\n",
        path_file=path_genes_heritability_text
    )

    pass


###############################################################################
# Procedure


def collect_select_report_write_heritability_genes(
    cohort=None,
    paths=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Report.
    utility.print_terminal_partition(level=1)
    print("... Heritability procedure for: " + str(cohort) + " persons...")
    utility.print_terminal_partition(level=2)

    # Read source information from file.
    source = read_source_initial(
        cohort=cohort,
        dock=paths["dock"],
    )

    # Verify that heritabilities are available for all genes.
    # simple method unavailable: ENSG00000022567
    # complex method unavailable: ENSG00000214860, ENSG00000180354
    check_genes(
        genes_selection=source["genes_selection"],
        genes_distribution=source["genes_distribution"],
        genes_heritability_complete=source["genes_heritability_complete"],
        path_genes=paths[cohort]["genes"],
    )

    # Calculate intersection of genes from selection, distribution, and
    # heritability procedures.
    genes_interest = find_intersection_heritability_genes(
        genes_selection=source["genes_selection"],
        genes_distribution=source["genes_distribution"],
        genes_heritability_complete=source["genes_heritability_complete"],
        path_genes=paths[cohort]["genes"],
        report=False,
    )
    # Collect and organize information about genes' heritabilities.
    data_genes_heritability = read_collect_organize_genes_heritabilities(
        genes=genes_interest,
        data_gene_annotation=source["data_gene_annotation"],
        path_genes=paths[cohort]["genes"],
        report=False,
    )
    # Filter heritabilities by their confidence intervals.
    # Allow 95% confidence intervals that span up to 33% of the reasonable
    # range (interval is 0.33, threshold is 0.33 / 2 = 0.167).
    # (interval is 0.50, threshold is 0.5 / 2 = 0.25)
    data_confidence = filter_heritabilities_confidence(
        data_genes_heritability=data_genes_heritability,
        threshold=0.25, # half the range of the 95% confidence interval
        report=True,
    )
    # Correct probabilities for false discovery rate across genes.
    data_discovery = calculate_organize_discoveries(
        threshold=0.05,
        data_genes_heritability=data_confidence,
    )
    # Organize data for genes' heritabilities.
    data_organization = organize_genes_heritability_data(
        data_genes_heritability=data_discovery,
        report=False,
    )

    # Select genes with most heritability.
    genes_heritability = select_heritable_genes(
        data_genes_heritability=data_organization,
        threshold_proportion=0.05,
        threshold_probability=0.1,
        report=True,
    )
    # Compile information.
    information = {
        "data_genes_heritability": data_organization,
        "genes_heritability": genes_heritability,
    }
    #Write product information to file.
    write_product(
        cohort=cohort,
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
        cohort="selection",
        paths=paths,
    )
    pass



if (__name__ == "__main__"):
    execute_procedure()
