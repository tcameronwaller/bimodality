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
import functools
import multiprocessing
import datetime
import gc

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import permutation
import organization
import restriction
import aggregation
import metric
import candidacy
import category
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source_local_initial(
    source_genes=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        source_genes (str): name of directory from which to obtain genes list
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    if source_genes == "split":
        path_source = os.path.join(dock, "split")
    elif source_genes == "candidacy":
        path_source = os.path.join(dock, "candidacy")
    elif source_genes == "combination":
        path_source = os.path.join(dock, "combination")
    path_genes = os.path.join(
        path_source, "genes.txt"
    )
    # Read information from file.
    genes = utility.read_file_text_list(path_genes)
    # Compile and return information.
    return {
        "genes": genes,
    }


def read_source(
    gene=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_permutation = os.path.join(dock, "permutation")
    path_permutations = os.path.join(
        path_permutation, "permutations.pickle"
    )
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    with open(path_permutations, "rb") as file_source:
        permutations = pickle.load(file_source)
    data_gene_samples_signals = pandas.read_pickle(path_gene)
    # Compile and return information.
    return {
        "data_gene_samples_signals": data_gene_samples_signals,
        "permutations": permutations,
    }


##########
# Distribution


def calculate_bimodality_metrics(
    values=None
):
    """
    Calculates bimodality metrics for a distribution.

    arguments:
        values (list): values of a distribution

    raises:

    returns:
        (dict<float>): values of metrics of bimodality

    """

    # Calculate metrics of bimodality.
    coefficient = metric.calculate_bimodality_coefficient(series=values)
    dip = metric.calculate_dip_statistic(series=values)
    mixture = metric.calculate_mixture_model_score(series=values)
    # Compile information.
    information = {
        "coefficient": coefficient,
        "dip": dip,
        "mixture": mixture,
    }

    # Return information.
    return information


def generate_null_metrics():
    """
    Generates null metrics.

    arguments:

    raises:

    returns:
        (dict<float>): values of metrics of bimodality

    """

    # Compile information.
    information = {
        "coefficient": float("nan"),
        "dip": float("nan"),
        "mixture": float("nan"),
    }

    # Return information.
    return information


def determine_gene_distributions(
    gene=None,
    modality=None,
    data_gene_samples_signals=None,
):
    """
    Prepares and describes the distributions of a gene's signals across
    persons.

    arguments:
        gene (str): identifier of gene
        modality (bool): whether to calculate metrics for the modality of
            gene's distribution
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (dict): information about the distribution of a gene's signals across
            persons

    """

    # Organize data for analysis.
    collection_organization = organization.execute_procedure(
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Determine distributions of gene's signals by multiple methods.
    imputation = determine_gene_distribution(
        gene=gene,
        modality=modality,
        method="imputation",
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        ),
    )
    availability = determine_gene_distribution(
        gene=gene,
        modality=modality,
        method="availability",
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        ),
    )

    # Compile information.
    information = {
        "organization": collection_organization,
        "imputation": imputation,
        "availability": availability,
    }

    # Return information.
    return information


# This function includes important parameters.
# Tissues for selection in restriction procedure.
# Count of tissues to require coverage in restriction procedure.
def determine_gene_distribution(
    gene=None,
    modality=None,
    method=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Prepares and describes the distributions of a gene's signals across
    persons.

    arguments:
        gene (str): identifier of gene
        modality (bool): whether to calculate metrics for the modality of
            gene's distribution
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information about the distribution of a gene's signals across
            persons

    """

    # Specify tissues for the imputation method of the restriction procedure.
    # This selection includes all non-sexual tissues with coverage of samples
    # from at least 200 persons.
    tissues = [
        "adipose", # 552
        #"adrenal", # 190
        "artery", # 551
        "blood", # 407
        "brain", # 254
        "colon", # 371
        "esophagus", # 513
        "heart", # 399
        #"liver", # 175
        "lung", # 427
        "muscle", # 564
        "nerve", # 414
        "pancreas", # 248
        #"pituitary", # 183
        "skin", # 583
        #"intestine", # 137
        #"spleen", # 162
        "stomach", # 262
        "thyroid", # 446
    ]

    # Restriction
    collection_restriction = restriction.execute_procedure(
        method=method,
        count=10,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            data_gene_persons_tissues_signals
        ),
    )

    # Aggregation
    # Use the final "data_gene_persons_tissues_signals" from restriction
    # procedure after imputation.
    collection_aggregation = aggregation.execute_procedure(
        data_gene_persons_tissues_signals=(
            collection_restriction["data_gene_persons_tissues_signals"]
        ),
    )

    # Determine distribution modality metrics.
    # Determine whether to calculate metrics for gene's distribution.
    if modality:
        population = candidacy.evaluate_gene_population(
            threshold=100,
            data_gene_persons_signals=collection_aggregation["data"],
        )
        signal = candidacy.evaluate_gene_signal(
            threshold=0.0,
            data_gene_persons_signals=collection_aggregation["data"],
        )
        if signal and population:
            # Calculate metrics.
            # Calculate metrics of bimodality.
            scores = calculate_bimodality_metrics(
                values=collection_aggregation["values"]
            )
            pass
        else:
            # Generate missing values.
            scores = generate_null_metrics()
            pass
    else:
        # Generate missing values.
        scores = generate_null_metrics()
        pass

    # Compile information.
    information = {
        "report_restriction": collection_restriction["report_gene"],
        "values": collection_aggregation["values"],
        "scores": scores,
        "report_aggregation": collection_aggregation["report_gene"],
    }

    # Return information.
    return information


# Use this function as a reference to extract information for a gene's
# distribution.
def extract_gene_distribution_information(
    method=None,
    observation=None,
):
    """
    Extracts information about a gene's distribution of pan-tissue signals
    across persons.

    arguments:
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        observation (dict): information about a gene's actual distribution of
            signals across persons and tissues

    raises:

    returns:

    """

    # Access information.
    report_organization = observation["organization"]["report_gene"]
    report_restriction = observation[method]["report_restriction"]
    report_aggregation = observation[method]["report_aggregation"]



    # Gene's signals across original, complete set of persons and tissues.
    data_gene_persons_tissues_signals_complete = (
        observation["organization"]["data_gene_persons_tissues_signals"]
    )
    # Mean of gene's signals across major tissue categories.
    data_gene_tissue_mean = report_organization["data_gene_tissue_mean"]
    # Variance of gene's signals across major tissue categories.
    data_gene_tissue_variance = (
        report_organization["data_gene_tissue_variance"]
    )
    # Standard deviation of gene's signals across major tissue categories.
    data_gene_tissue_deviation = (
        report_organization["data_gene_tissue_deviation"]
    )
    # Coefficient of variation of gene's signals across major tissue
    # categories.
    data_gene_tissue_variation = (
        report_organization["data_gene_tissue_variation"]
    )
    # Dispersion in gene's signals across major tissue categories.
    data_gene_tissue_dispersion = (
        report_organization["data_gene_tissue_dispersion"]
    )



    # Gene's signals across tissues and persons, before any imputation.
    data_gene_persons_tissues_signals = (
        report_restriction["data_gene_persons_tissues_signals"]
    )
    # Boolean availability of valid gene's signals across tissues and persons.
    data_gene_persons_tissues = (
        report_restriction["data_gene_persons_tissues"]
    )
    # Count of tissues for which gene has valid signals across persons.
    data_gene_persons_tissues_count = (
        report_restriction["data_gene_persons_tissues_count"]
    )
    # Count of persons for which gene has valid signals across adequate
    # tissues.
    restriction_persons = report_restriction["persons"]
    # Mean count of tissues per person.
    restriction_tissues_mean = report_restriction["tissues_mean"]
    # Median count of tissues per person.
    restriction_tissues_median = report_restriction["tissues_median"]



    # Gene's pan-tissue aggregate signals across persons.
    data_gene_persons_signals = (
        report_aggregation["data_gene_persons_signals"]
    )
    # Count of persons for which gene has valid signals across adequate
    # tissues.
    aggregation_persons = report_aggregation["persons"]

    pass



##########
# Shuffle


def shuffle_gene_distributions(
    gene=None,
    modality=None,
    permutations=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Prepares and describes the distribution of a gene's signals across persons.

    arguments:
        gene (str): identifier of gene
        modality (bool): whether to calculate metrics for the modality of
            gene's distribution
        permutations (list<list<list<int>>>): matrices of permutation indices
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information about the distribution of a gene's signals across
            persons

    """


    # Determine distributions of gene's signals by multiple methods.
    imputation = shuffle_gene_distribution(
        gene=gene,
        modality=modality,
        method="imputation",
        permutations=permutations,
        data_gene_persons_tissues_signals=data_gene_persons_tissues_signals,
    )
    availability = shuffle_gene_distribution(
        gene=gene,
        modality=modality,
        method="availability",
        permutations=permutations,
        data_gene_persons_tissues_signals=data_gene_persons_tissues_signals,
    )

    # Compile information.
    information = {
        "imputation": imputation,
        "availability": availability,
    }

    # Return information.
    return information


def shuffle_gene_distribution(
    gene=None,
    modality=None,
    method=None,
    permutations=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Prepares and describes the distribution of a gene's signals across persons.

    arguments:
        gene (str): identifier of gene
        modality (bool): whether to calculate metrics for the modality of
            gene's distribution
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        permutations (list<list<list<int>>>): matrices of permutation indices
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information about the distribution of a gene's signals across
            persons

    """

    # Initialize collections of metrics.
    coefficients = list()
    dips = list()
    mixtures = list()
    # Iterate on permutations.
    for shuffle_matrix in permutations:
        # Shuffle gene's signals.
        data_shuffle = permutation.shuffle_gene_signals(
            data_gene_signals=data_gene_persons_tissues_signals,
            shuffle=shuffle_matrix
        )
        # Determine distributions of gene's signals
        collection = determine_gene_distribution(
            gene=gene,
            modality=modality,
            method=method,
            data_gene_persons_tissues_signals=data_shuffle,
        )

        # Collect metrics.
        coefficients.append(
            collection["scores"]["coefficient"]
        )
        dips.append(
            collection["scores"]["dip"]
        )
        mixtures.append(
            collection["scores"]["mixture"]
        )
        pass
    # Compile information.
    information = {
        "coefficient": coefficients,
        "dip": dips,
        "mixture": mixtures,
    }

    # Return information.
    return information


##########
# Product


def write_product(gene=None, dock=None, information=None):
    """
    Writes product information to file.

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files.
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_distribution = os.path.join(dock, "distribution")
    utility.create_directory(path_distribution)
    path_gene = os.path.join(path_distribution, gene)
    utility.create_directory(path_gene)
    path_imputation = os.path.join(path_gene, "imputation")
    utility.create_directory(path_imputation)
    path_availability = os.path.join(path_gene, "availability")
    utility.create_directory(path_availability)

    path_report_organization = os.path.join(
        path_gene, "report_organization.pickle"
    )
    # Imputation method
    path_imputation_report_restriction = os.path.join(
        path_imputation, "report_restriction.pickle"
    )
    path_imputation_report_aggregation = os.path.join(
        path_imputation, "report_aggregation.pickle"
    )
    path_imputation_scores = os.path.join(
        path_imputation, "scores.pickle"
    )
    path_imputation_permutations = os.path.join(
        path_imputation, "permutations.pickle"
    )
    # Availability method
    path_availability_report_restriction = os.path.join(
        path_availability, "report_restriction.pickle"
    )
    path_availability_report_aggregation = os.path.join(
        path_availability, "report_aggregation.pickle"
    )
    path_availability_scores = os.path.join(
        path_availability, "scores.pickle"
    )
    path_availability_permutations = os.path.join(
        path_availability, "permutations.pickle"
    )

    # Write information to file.
    with open(path_report_organization, "wb") as file_product:
        pickle.dump(information["report_organization"], file_product)

    # Imputation method
    with open(path_imputation_report_restriction, "wb") as file_product:
        pickle.dump(information["imputation"]["report_restriction"], file_product)
    with open(path_imputation_report_aggregation, "wb") as file_product:
        pickle.dump(
            information["imputation"]["report_aggregation"], file_product
        )
    with open(path_imputation_scores, "wb") as file_product:
        pickle.dump(information["imputation"]["scores"], file_product)
    with open(path_imputation_permutations, "wb") as file_product:
        pickle.dump(information["imputation"]["permutations"], file_product)

    # Imputation method
    with open(path_availability_report_restriction, "wb") as file_product:
        pickle.dump(
            information["availability"]["report_restriction"], file_product
        )
    with open(path_availability_report_aggregation, "wb") as file_product:
        pickle.dump(
            information["availability"]["report_aggregation"], file_product
        )
    with open(path_availability_scores, "wb") as file_product:
        pickle.dump(information["availability"]["scores"], file_product)
    with open(path_availability_permutations, "wb") as file_product:
        pickle.dump(information["availability"]["permutations"], file_product)

    pass


###############################################################################
# Procedure

# Change procedure to run and collect both by "availability" and "imputation" methods...


def execute_procedure(
    gene=None,
    permutations=None,
    data_gene_samples_signals=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of gene
        permutations (list<list<list<int>>>): matrices of permutation indices
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Organization



    # Restriction
    # Support either availability or imputation methods.
    # (I wish I could do this without removing data for persons who don't qualify... just assign missing values?)


    # Aggregation
    # TODO: apply the candidacy quality checks (population, signal, tissue)


    # Modality


    # TODO: prepare a table with all original persons with missing values if they were'nt used in final distribution


    # Compilation


    ##################Old Stuff...###########################################3

    # Prepare and describe distribution of real gene's signals.

    # Determine gene's distributions of aggregate tissue signals across
    # persons.
    observation = determine_gene_distributions(
        gene=gene,
        modality=True,
        data_gene_samples_signals=data_gene_samples_signals,
    )

    # Shuffle gene's distributions.
    shuffle = shuffle_gene_distributions(
        gene=gene,
        modality=True,
        permutations=permutations,
        data_gene_persons_tissues_signals=(
            observation["organization"]["data_gene_persons_tissues_signals"]
        ),
    )

    # Compile information.
    information_imputation = {
        "report_restriction": (
            observation["imputation"]["report_restriction"]
        ),
        "report_aggregation": (
            observation["imputation"]["report_aggregation"]
        ),
        "scores": observation["imputation"]["scores"],
        "permutations": shuffle["imputation"],
    }
    information_availability = {
        "report_restriction": (
            observation["availability"]["report_restriction"]
        ),
        "report_aggregation": (
            observation["availability"]["report_aggregation"]
        ),
        "scores": observation["availability"]["scores"],
        "permutations": shuffle["availability"],
    }
    information = {
        "report_organization": observation["organization"]["report_gene"],
        "imputation": information_imputation,
        "availability": information_availability,
    }
    #Write product information to file.
    write_product(gene=gene, dock=dock, information=information)

    pass


def execute_procedure_local(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_local_initial(source_genes="split", dock=dock)

    print("count of genes: " + str(len(source["genes"])))
    #print("count of permutations: " + str(len(source["permutations"])))

    # Report date and time.
    start = datetime.datetime.now()
    print(start)

    # Remove previous files to avoid version or batch confusion.
    path_distribution = os.path.join(dock, "distribution")
    utility.remove_directory(path=path_distribution)

    # Set up partial function for iterative execution.
    # Each iteration uses the same values for "genes_signals", "permutations", and
    # "dock" variables.
    execute_procedure_gene = functools.partial(
        execute_procedure_local_sub,
        dock=dock
    )

    # Initialize multiprocessing pool.
    # Iterative shuffle procedure.
    # 4 concurrent processes require approximately 60% of 32 Gigabytes memory.
    # 5 concurrent processes require approximately 75% of 32 Gigabytes memory.
    # 6 concurrent processes require approximately 90% of 32 Gigabytes memory.
    # Index shuffle procedure.
    # 6 concurrent processes require approximately 50% of 32 Gigabytes memory.
    # 7 concurrent processes require approximately 65% of 32 Gigabytes memory.
    #pool = multiprocessing.Pool(processes=os.cpu_count())
    pool = multiprocessing.Pool(processes=100)

    # Iterate on genes.
    check_genes=[
        "ENSG00000198965",
    ]
    #report = pool.map(execute_procedure_gene, check_genes)
    #report = pool.map(execute_procedure_gene, source["genes"][0:32])
    report = pool.map(execute_procedure_gene, source["genes"])

    # Report.
    #print("Process complete for the following genes...")
    #print(str(len(report)))

    # Report date and time.
    end = datetime.datetime.now()
    print(end)
    print("duration: " + str(end - start))

    pass


def execute_procedure_local_sub(gene=None, dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        gene (str): identifier of single gene for which to execute the process.

    raises:

    returns:

    """

    # Enable automatic garbage collection to clear memory.
    #gc.enable()

    # Report gene.
    #print("gene: " + gene)

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Collect garbage to clear memory.
    #gc.collect()

    # Execute procedure.
    execute_procedure(
        gene=gene,
        permutations=source["permutations"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        dock=dock
    )

    # Report contents of directory.
    path_distribution = os.path.join(dock, "distribution")
    directories = os.listdir(path_distribution)
    count = len(directories)
    if (count % 10 == 0):
        print("complete genes: " + str(len(directories)))

    pass


def execute_procedure_remote(dock=None, gene=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        gene (str): identifier of single gene for which to execute the process.

    raises:

    returns:

    """

    # Execute procedure.
    execute_procedure_remote_sub(
        gene=gene,
        dock=dock
    )

    pass


def execute_procedure_remote_sub(gene=None, dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        gene (str): identifier of single gene for which to execute the process.

    raises:

    returns:

    """

    # Enable automatic garbage collection to clear memory.
    gc.enable()

    # Report gene.
    #print("gene: " + gene)

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Collect garbage to clear memory.
    gc.collect()

    # Execute procedure.
    execute_procedure(
        gene=gene,
        permutations=source["permutations"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        dock=dock
    )

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
