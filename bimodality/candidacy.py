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
import time

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import organization
import pipe
import regression
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Source


def read_source_initial(
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
    path_assembly = os.path.join(dock, "assembly")
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_samples_signals = pandas.read_pickle(path_gene)
    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_samples_signals": data_gene_samples_signals,
    }


##########
# Candidacy


def determine_gene_signal_distribution(
    gene=None,
    method=None,
    count=None,
    data_gene_samples_signals=None,
):
    """
    Determine gene's distribution of aggregate tissue scores across persons.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (bool): whether the gene is a candidate

    """

    # Organize data for analysis.
    collection_organization = organization.execute_procedure(
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Prepare and describe distribution of real gene's signals.
    # Method for selection is either "availability" or "imputation".
    # Specific tissues are only relevant to "imputation" method.

    if False:
        tissues = [
            "adipose", "blood", "colon", "esophagus", "heart", "muscle",
            "lung", "nerve", "skin", "thyroid",
        ] # 416
    # This strategy sets a simple threshold at 200 persons with samples for
    # each of these tissues.
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
        #"salivary", # 97 # now excluded at "selection" procedure... threshold 100
        #"spleen", # 162
        "stomach", # 262
        "thyroid", # 446
    ] # 300
    distribution = pipe.describe_gene_signal_distribution(
        method=method,
        count=count,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        )
    )

    # Compile information.
    information = {
        "data_gene_persons_tissues": (
            distribution["report_restriction"]["data_gene_persons_tissues"]
        ),
        "data_gene_persons_signals": (
            distribution["report_distribution"]["data_gene_persons_signals"]
        )
    }
    # Return information.
    return information


def evaluate_gene_signals(
    data_gene_persons_signals=None
):
    """
    Evaluates the suitability of a gene's aggregate tissue signals across
    persons.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate signals across persons

    raises:

    returns:
        (bool): whether the gene's signals are adequate for analysis

    """

    # Extract values of distribution.
    # Remove missing values.
    data = data_gene_persons_signals.copy(deep=True)

    # Count non-missing values.
    data_valid = data.dropna(
        axis="index",
        how="any",
        inplace=False,
    )
    count_valid = data_valid.size

    # Count non-zero values.
    data_nonzero = (data != 0)
    data_signal = data.loc[data_nonzero.any(axis="columns"), : ]
    count_signal = data_signal.size

    # Calculate standard deviation of values.
    deviation = numpy.std(data["value"].values, axis=0)

    # Evaluate distribution's suitability for analysis.
    if (
        (count_valid > 10) and
        (deviation > 0)
    ):
        return True
    else:
        return False


def evaluate_gene_tissues(
    data_samples_tissues_persons=None,
    data_gene_persons_tissues=None,
    data_gene_persons_signals=None,
):
    """
    Evaluates the selection of tissues for calculation of gene's aggregate
    scores.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_persons_tissues (object): Pandas data frame of a gene's
            selection of tissues across persons
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate signals across persons

    raises:

    returns:
        (bool): whether the gene's selection of tissues is acceptable

    """

    # Organize information about gene's scores across persons and tissues.
    data_persons_tissues_signals = regression.organize_gene_persons_signals(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_tissues=data_gene_persons_tissues,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    # Determine whether tissue selection confounds gene's aggregate scores
    # across persons.
    #print(data_persons_tissues_signals)
    tissues = data_persons_tissues_signals.columns.to_list()
    tissues.remove("value")
    tissues.remove("female")
    tissues.remove("age")
    report = regression.evaluate_variance_by_binary_groups(
        groups=tissues,
        data=data_persons_tissues_signals,
    )

    # Evaluate whether gene's scores differ by tissue selection.
    if (
        (math.isnan(report["probability"])) or
        (report["probability"] > 0.07)
    ):
        return True
    else:
        return False


def evaluate_gene_candidacy(
    gene=None,
    method=None,
    count=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (bool): whether the gene is a candidate

    """

    # Determine gene's distribution of aggregate tissue scores across persons.
    collection = determine_gene_signal_distribution(
        gene=gene,
        method=method,
        count=count,
        data_gene_samples_signals=data_gene_samples_signals,
    )

    # Determine whether the gene's distribution is adequate.
    # Determine whether gene's aggregate scores across tissues have adequate
    # variance to apply metrics to the distribution.
    signal = evaluate_gene_signals(
        data_gene_persons_signals=collection["data_gene_persons_signals"],
    )

    # Determine whether the gene's distribution is dependent on tissue
    # selection.
    tissue = evaluate_gene_tissues(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_tissues=collection["data_gene_persons_tissues"],
        data_gene_persons_signals=collection["data_gene_persons_signals"],
    )

    # Evaluate gene's candidacy.
    if (
        signal and
        tissue
    ):
        return True
    else:
        if (not signal):
            #print("Gene: " + gene + " failed signal candidacy!")
            pass
        if (not tissue):
            #print("Gene: " + gene + " failed tissue candidacy!")
            pass
        return False


##########
# Product


def write_product(
    gene=None,
    method=None,
    dock=None,
):
    """
    Writes product information to file.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        method (str): list to which to append the gene, "all", "imputation", or
            "availability"
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:

    """

    # Specify directories and files.
    path_candidacy = os.path.join(dock, "candidacy")
    utility.confirm_path_directory(path_candidacy)
    path_method = os.path.join(path_candidacy, method)
    utility.confirm_path_directory(path_method)
    path_gene = os.path.join(path_method, (gene))

    # Write information to file.
    utility.write_file_text_list(
        information=["null"], path_file=path_gene
    )

    pass


def read_collect_write_genes(
    dock=None,
):
    """
    Collects information about genes.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:
        (list<str>): information about genes

    """

    # Specify directories and files.
    path_candidacy = os.path.join(dock, "candidacy")
    path_all = os.path.join(path_candidacy, "all")
    path_imputation = os.path.join(path_candidacy, "imputation")
    path_availability = os.path.join(path_candidacy, "availability")
    path_genes = os.path.join(path_candidacy, "genes.txt")
    path_genes_imputation = os.path.join(
        path_candidacy, "genes_imputation.txt"
    )
    path_genes_availability = os.path.join(
        path_candidacy, "genes_availability.txt"
    )
    # Read genes.
    genes_imputation = os.listdir(path_imputation)
    genes_availability = os.listdir(path_availability)
    utility.remove_directory(path=path_all)
    utility.remove_directory(path=path_imputation)
    utility.remove_directory(path=path_availability)
    print("Count of imputation genes: " + str(len(genes_imputation)))
    print("Count of availability genes: " + str(len(genes_availability)))
    # Combine genes.
    genes_combination = utility.filter_unique_union_elements(
        list_one=genes_imputation,
        list_two=genes_availability,
    )
    print("Count of consensus genes: " + str(len(genes_combination)))

    # Write information to file.
    utility.write_file_text_list(
        information=genes_imputation, path_file=path_genes_imputation
    )
    utility.write_file_text_list(
        information=genes_availability, path_file=path_genes_availability
    )
    utility.write_file_text_list(
        information=genes_combination, path_file=path_genes
    )

    pass


def report_gene_tissue_persons(
    gene=None,
    count=None,
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of a gene
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )
    # Determine gene's distribution of aggregate tissue scores across persons.
    imputation = determine_gene_signal_distribution(
        gene=gene,
        method="imputation",
        count=count,
        data_gene_samples_signals=source["data_gene_samples_signals"],
    )
    availability = determine_gene_signal_distribution(
        gene=gene,
        method="availability",
        count=count,
        data_gene_samples_signals=source["data_gene_samples_signals"],
    )
    utility.print_terminal_partition(level=2)
    print("Persons by imputation and availability methods...")
    print(imputation["data_gene_persons_signals"].size)
    print(availability["data_gene_persons_signals"].size)

    pass



###############################################################################
# Procedure

def execute_procedure(
    gene=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of a gene
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Evaluate gene's candidacy.
    # Imputation
    candidacy_imputation = evaluate_gene_candidacy(
        gene=gene,
        method="imputation",
        count=10,
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_samples_signals=data_gene_samples_signals,
    )
    if candidacy_imputation:
        #Write product information to file.
        write_product(
            gene=gene,
            method="imputation",
            dock=dock,
        )
    # Availability
    candidacy_availability = evaluate_gene_candidacy(
        gene=gene,
        method="availability",
        count=10,
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_samples_signals=data_gene_samples_signals,
    )
    if candidacy_availability:
        #Write product information to file.
        write_product(
            gene=gene,
            method="availability",
            dock=dock,
        )

    #Write product information to file.
    write_product(
        gene=gene,
        method="all",
        dock=dock,
    )

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
    source = read_source_initial(
        source_genes="split",
        dock=dock
    )
    print("count of genes: " + str(len(source["genes"])))

    # Report date and time.
    start = datetime.datetime.now()
    print(start)

    # Remove previous files to avoid version or batch confusion.
    path_candidacy = os.path.join(dock, "candidacy")
    utility.remove_directory(path=path_candidacy)

    # Set up partial function for iterative execution.
    # Each iteration uses the same values for "genes_signals", "shuffles", and
    # "dock" variables.
    execute_procedure_gene = functools.partial(
        execute_procedure_local_sub,
        dock=dock
    )

    # Initialize multiprocessing pool.
    #pool = multiprocessing.Pool(processes=os.cpu_count())
    pool = multiprocessing.Pool(processes=8)

    # Iterate on genes.
    check_genes=[
        "ENSG00000198965",
    ]
    #report = pool.map(execute_procedure_gene, check_genes)
    #report = pool.map(execute_procedure_gene, source["genes"][0:10])
    report = pool.map(execute_procedure_gene, source["genes"])

    # Pause procedure.
    time.sleep(5.0)

    # Report.
    #print("Process complete for the following genes...")
    #print(str(len(report)))

    report_gene_tissue_persons(
        gene="ENSG00000231925", # TAPBP
        count=10,
        dock=dock,
    )

    # Collect genes.
    read_collect_write_genes(dock=dock)

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

    # Report gene.
    #print("gene: " + gene)

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Execute procedure.
    execute_procedure(
        gene=gene,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        dock=dock
    )

    # Report contents of directory.
    path_candidacy = os.path.join(dock, "candidacy")
    path_all = os.path.join(path_candidacy, "all")
    files = os.listdir(path_all)
    count = len(files)
    if (count % 10 == 0):
        print("complete genes: " + str(len(files)))

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
