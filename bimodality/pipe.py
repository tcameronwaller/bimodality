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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import metric
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
    path_split = os.path.join(dock, "split")
    path_genes = os.path.join(
        path_split, "genes.txt"
    )
    path_signal = os.path.join(
        path_split, "genes_signals_patients_tissues.pickle"
    )
    path_shuffle = os.path.join(dock, "shuffle")
    path_shuffles = os.path.join(
        path_shuffle, "shuffles.pickle"
    )
    # Read information from file.
    genes = utility.read_file_text_list(path_genes)
    with open(path_signal, "rb") as file_source:
        genes_signals_patients_tissues = pickle.load(file_source)
    with open(path_shuffles, "rb") as file_source:
        shuffles = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes": genes,
        "genes_signals_patients_tissues": genes_signals_patients_tissues,
        "shuffles": shuffles,
    }


def calculate_sum_gene_signal_tissue(data_gene_signals=None):
    """
    Calculates the sum of genes' signals across tissues for each patient.

    arguments:
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues.

    raises:

    returns:
        (object): Pandas data frame of sum of standard score signals for all
            genes across tissues for each patient.

    """

    # Drop the tissue designations.
    if False:
        data_gene_signals.reset_index(level=["tissue"], inplace=True)
        data_gene_signals.drop(
            labels="tissue",
            axis="columns",
            inplace=True
        )
    # Apply calculation to groups.
    data_sum = data_gene_signals.sum(axis="columns").to_frame(name="value")
    return data_sum


def calculate_mean_gene_signal_patient(data_gene_signals=None):
    """
    Calculates the sum of genes' signals across patients for each tissue.

    arguments:
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues.

    raises:

    returns:
        (object): Pandas data frame of sum of standard score signals for all
            genes across tissues for each patient.

    """

    data = data_gene_signals.aggregate(statistics.mean, axis="index")
    return data


def calculate_bimodality_metrics(values=None):
    """
    Calculates metrics of bimodality.

    arguments:
        values (list<float>): Values for which to evaluate distribution.

    raises:

    returns:
        (dict<float>): Values of metrics of bimodality.

    """

    # Calculate metrics of bimodality.
    coefficient = metric.calculate_bimodality_coefficient(series=values)
    dip = metric.calculate_dip_statistic(series=values)

    # Compile information.
    information = {
        "coefficient": coefficient,
        "dip": dip
    }

    # Return information.
    return information


def analyze_gene_signal_distribution(data_gene_signals=None):
    """
    Analyzes the distribution of a single gene's signals across tissues and
    patients.

    arguments:
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues.

    raises:

    returns:
        (object): Pandas data frame of sum of standard score signals for all
            genes across tissues for each patient.

    """

    # Check that the mean of gene's signals across patients are consistent with
    # shuffles.
    # As the signals are in standard score space, the means center around zero.
    # Consistency in the actual values suggests that the shuffles were correct.
    if False:
        data_check = calculate_mean_gene_signal_patient(
            data_gene_signals=data_gene_signals
        )
        print(
            "Check that sum of all patients signals for each tissue does not " +
            "change across shuffles."
        )
        print(data_check)
        print(data_gene_signals.shape)

    # Aggregate signals by tissues.
    data_gene_signal_aggregation = calculate_sum_gene_signal_tissue(
        data_gene_signals=data_gene_signals
    )
    #print(data_gene_signal_aggregation)
    values = data_gene_signal_aggregation["value"].to_list()
    #print(values)

    # Calculate metrics of bimodality.
    scores = calculate_bimodality_metrics(values=values)
    #print(scores)

    # Return information.
    return scores


def shuffle_gene_signals(
    data_gene_signals=None,
    shuffle=None
):
    """
    Shuffles the association between tissues and patients of gene's signals.

    arguments:
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues.
        shuffle (list<list<int>>): Matrix of indices.

    raises:

    returns:
        (object): Pandas data frame of a single gene's signals across specific
            patients and tissues.

    """

    # Determine identities and counts of tissues and patients.
    tissues = data_gene_signals.columns.to_list()
    patients = data_gene_signals.index.to_list()
    count_tissues = len(tissues)
    count_patients = len(patients)

    # Collect gene's signals.
    gene_signals_shuffle = list()
    for index_patient in range(count_patients):
        patient = patients[index_patient]
        record = dict()
        record["patient"] = patient
        for index_tissue in range(count_tissues):
            tissue = tissues[index_tissue]
            index_patient_shuffle = shuffle[index_tissue][index_patient]
            signal = (
                data_gene_signals.iloc[index_patient_shuffle, index_tissue]
            )
            record[tissue] = signal
        gene_signals_shuffle.append(record)

    # Convert records to data frame.
    data = utility.convert_records_to_dataframe(
        records=gene_signals_shuffle
    )
    data.rename_axis(
        columns="tissue",
        axis="columns",
        copy=False,
        inplace=True
    )
    data.set_index(
        "patient",
        drop=True,
        inplace=True,
    )
    return data


def collect_shuffle_distributions(
    data_gene_signals=None,
    shuffles=None
):
    """
    Analyzes the distribution of a single gene's signals across tissues and
    patients.

    arguments:
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues.
        shuffles (list<list<list<int>>>): Matrices of indices.

    raises:

    returns:
        (dict<list<float>>): Values of metrics of bimodality.

    """

    # Initialize collection of metrics.
    distributions = dict()
    distributions["coefficient"] = list()
    distributions["dip"] = list()
    # Iterate on shuffles.
    for shuffle in shuffles:

        # Shuffle gene's signals.
        data_gene_signals_shuffle = shuffle_gene_signals(
            data_gene_signals=data_gene_signals,
            shuffle=shuffle
        )
        #print(data_gene_signals_shuffle.iloc[0:10, :])

        # Analyze gene's signals.
        scores = analyze_gene_signal_distribution(
            data_gene_signals=data_gene_signals_shuffle
        )

        # Collect metrics.
        distributions["coefficient"].append(scores["coefficient"])
        distributions["dip"].append(scores["dip"])

    # Return information.
    return distributions


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
    path_pipe = os.path.join(dock, "pipe")
    utility.confirm_path_directory(path_pipe)
    path_gene = os.path.join(path_pipe, gene)
    utility.confirm_path_directory(path_gene)
    path_scores = os.path.join(
        path_gene, "scores.pickle"
    )
    path_distributions = os.path.join(
        path_gene, "distributions.pickle"
    )
    # Write information to file.
    with open(path_scores, "wb") as file_product:
        pickle.dump(information["scores"], file_product)
    with open(path_distributions, "wb") as file_product:
        pickle.dump(information["distributions"], file_product)

    pass


def execute_process(gene=None, genes_signals=None, shuffles=None, dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        genes_signals (dict<object>): Collection of matrices.
        shuffles (list<list<list<int>>>): Matrices of indices.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Access gene's signals.
    # "ENSG00000029363"
    data_gene_signals = (
        genes_signals[gene].copy(deep=True)
    )
    #print(gene)
    #print(data_gene_signals.iloc[0:10, :])

    # Analyze the gene's real signals.
    scores = analyze_gene_signal_distribution(
        data_gene_signals=data_gene_signals
    )
    #print(scores)

    # Collect random distributions of scores.
    distributions = collect_shuffle_distributions(
        data_gene_signals=data_gene_signals,
        shuffles=shuffles
    )
    #print(distributions["coefficient"][0:100])
    #print(str(len(distributions["coefficient"])))

    # Compile information.
    information = {
        "scores": scores,
        "distributions": distributions
    }
    #Write product information to file.
    write_product(gene=gene, dock=dock, information=information)

    # Return gene identifier as verification of process completion.
    #return gene

    pass


###############################################################################
# Procedure


def execute_procedure_local(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # 4. Collect shuffle metrics for gene
    # 4.1. Apply shuffles to the gene's matrix... iterate on the previously-computed shuffle indices
    # 4.2. For each shuffle, aggregate and calculate metrics
    # 4.3. Collect metrics for the gene... probably within a list or dictionary?
    # 5. Compile a report for the gene, probably consisting of multiple files.
    # 6. Create a directory on file for the gene.
    # 7. Save to file the report for the individual gene: real metrics, shuffle distribution, etc.

    # Read source information from file.
    source = read_source(dock=dock)

    # Report date and time.
    print(datetime.datetime.now())

    # Running 8 processes concurrently requires approximately 22.5 Gigabytes of
    # memory.
    # Each process, then, requires approximately 3 Gigabytes of memory.
    # Running each batch of processes equal to the count of cores requires
    # approximately 12 minutes.
    # (18,000 genes) / (8 cores) = 2,250 batches
    # (12 minutes / batch) * (2,250 batches) = (27,000 minutes) = (450 hours)

    # Set up partial function for iterative execution.
    # Each iteration uses the same values for "genes_signals", "shuffles", and
    # "dock" variables.
    execute_process_gene = functools.partial(
        execute_process,
        genes_signals=source["genes_signals_patients_tissues"],
        shuffles=source["shuffles"],
        dock=dock
    )

    # Initialize multiprocessing pool.
    pool = multiprocessing.Pool(processes=os.cpu_count())

    # Iterate on genes.
    #report = pool.map(execute_process_gene, source["genes"][0:5])
    report = pool.map(execute_process_gene, source["genes"][0:100])

    # Report.
    #print("Process complete for the following genes...")
    #print(str(len(report)))

    # Report date and time.
    print(datetime.datetime.now())

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

    # 4. Collect shuffle metrics for gene
    # 4.1. Apply shuffles to the gene's matrix... iterate on the previously-computed shuffle indices
    # 4.2. For each shuffle, aggregate and calculate metrics
    # 4.3. Collect metrics for the gene... probably within a list or dictionary?
    # 5. Compile a report for the gene, probably consisting of multiple files.
    # 6. Create a directory on file for the gene.
    # 7. Save to file the report for the individual gene: real metrics, shuffle distribution, etc.

    # Read source information from file.
    source = read_source(dock=dock)

    # Execute procedure.
    execute_process(
        gene=gene,
        genes_signals=source["genes_signals_patients_tissues"],
        shuffles=source["shuffles"],
        dock=dock
    )

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
