"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os
import statistics
import pickle

# Relevant

import numpy
import pandas
#import cmapPy.pandasGEXpress.parse
#import cmapPy.pandasGEXpress.gct2gctx

# Custom

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
    path_organization = os.path.join(dock, "organization")
    path_gene_signal = os.path.join(
        path_organization, "data_gene_persons_tissues_signals.pickle"
    )
    # Read information from file.
    data_gene_persons_tissues_signals = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_gene_persons_tissues_signals": data_gene_persons_tissues_signals,
    }



##########
# Report.


##################################################
# Scrap

def calculate_gene_signal_mean_by_persons_tissues_old(
    persons_tissues_samples=None,
    data_gene_signal=None
):
    """
    Calculates the mean values of genes' signals in each tissue of each
    person.

    Redundant samples exist for many persons and tissues.
    Calculate the mean signal for each gene from all of these redundant
    samples.

    Product data format.
    Indices by person and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # person tissue
    # bob     kidney  5.0    10.0   7.0
    # bob     skin    5.0    10.0   7.0
    # bob     liver   5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # april   liver   2.0    3.0    4.0
    # april   skin    2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        data_gene_signal (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of mean values of genes' signals for all
            tissues of all persons.

    """

    # Create index for convenient access.
    data_gene_signal_index = data_gene_signal.set_index("Name")
    # Determine list of all genes with signals.
    genes = data_gene_signal["Name"].to_list()
    # Collect mean signals of all genes for each person and tissue.
    summary = list()
    # Iterate on persons.
    for person in persons_tissues_samples:
        # Iterate on tissues for which person has samples.
        for tissue in persons_tissues_samples[person]:
            # Collect signals for all genes in each sample.
            signals_genes = dict()
            for sample in persons_tissues_samples[person][tissue]:
                for gene in genes:
                    signal = data_gene_signal_index.at[gene, sample]
                    # Determine whether an entry already exists for the gene.
                    if gene in signals_genes:
                        signals_genes[gene].append(signal)
                    else:
                        signals_genes[gene] = list([signal])
            # Compile record.
            record = dict()
            record["person"] = person
            record["tissue"] = tissue
            # Calculate mean signal for each gene in each person and tissue.
            for gene in genes:
                signal_mean = statistics.mean(signals_genes[gene])
                record[gene] = signal_mean
            # Include record.
            summary.append(record)
    # Organize data.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    #data_summary_index = data_summary.set_index(
    #    ["person", "tissue"], append=False, drop=False
    #)
    return data_summary

def calculate_gene_signal_median_by_tissues_old(
    tissues=None,
    persons=None,
    tissues_persons_samples=None,
    data_gene_signal_mean=None
):
    """
    Calculates the median values of signals for all genes across specific
    persons and tissues.

    Redundant samples exist for many persons and tissues.
    Calculate the mean signal for each gene from all of these redundant
    samples.

    Product data format.
    Index by gene allows convenient access.
    ##################################################
    # gene   brain colon liver heart kidney
    # abc1   7     3     2     6     3
    # mpc1   5     3     2     6     3
    # rip1   3     3     2     6     3
    # rnt1   2     3     2     6     3
    # rnx3   9     3     2     6     3
    # sst5   3     3     2     6     3
    # sph6   2     3     2     6     3
    # xtr4   8     3     2     6     3
    ##################################################

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons with signal for tissues of interest.
        tissues_persons_samples (dict<dict<list<str>>): Samples for each
            person of each tissue.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of median values of genes' signals for all
            tissues.

    """

    # Create index for convenient access.
    data_gene_signal_mean_index = data_gene_signal_mean.set_index(
        ["person", "tissue"], append=False, drop=True
    )
    # Extract genes from names of columns.
    genes = data_gene_signal_mean_index.columns.to_list()
    # Collect median signals of all genes for each tissue.
    summary = list()
    # Iterate on genes.
    for gene in genes:
        # Compile record.
        record = dict()
        record["gene"] = gene
        # Iterate on tissues.
        for tissue in tissues:
            # Collect gene's signals in the tissue across all persons.
            signals_tissue = list()
            # Iterate on persons.
            for person in tissues_persons_samples[tissue]:
                signal_tissue_person = data_gene_signal_mean_index.loc[
                    (person, tissue), gene
                ]
                signals_tissue.append(signal_tissue_person)
            # Calculate median value of gene's signal in tissue across all
            # persons.
            signal_tissue = statistics.median(signals_tissue)
            record[tissue] = signal_tissue
        # Include record.
        summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    return data_summary

def collect_gene_signal_by_persons_tissues_old(
    tissues=None,
    persons=None,
    persons_tissues_samples=None,
    data_gene_signal_median_tissue=None,
    data_gene_signal_mean=None
):
    """
    Collects the signals for all genes in each tissue of each person.

    Product data format.
    Indices by person and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # person tissue
    # bob     kidney  5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons with signal for tissues of interest.
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        data_gene_signal_median_tissue (object): Pandas data frame of median
            values of genes' signals for all tissues.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific persons and tissues.

    """

    # Create indices for convenient access.
    data_gene_signal_mean_index = data_gene_signal_mean.set_index(
        ["person", "tissue"], append=False, drop=True
    )
    data_gene_signal_median_tissue_index = (
        data_gene_signal_median_tissue.set_index("gene")
    )
    # Extract genes from names of columns.
    genes = data_gene_signal_mean_index.columns.to_list()
    # Collect signals of all genes for each person and tissue.
    summary = list()
    # Iterate on persons of interest.
    for person in persons:
        # Iterate on tissues of interest.
        for tissue in tissues:
            # Compile record.
            record = dict()
            record["person"] = person
            record["tissue"] = tissue
            # Determine whether the person has samples for the tissue.
            if tissue in persons_tissues_samples[person]:
                # Collect mean values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_mean_index.loc[
                        (person, tissue), gene
                    ]
                    record[gene] = signal
            else:
                # Collect imputation values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_median_tissue_index.loc[
                        gene, tissue
                    ]
                    record[gene] = signal
            # Include record.
            summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    #data_summary_index = data_summary.set_index(
    #    ["person", "tissue"], append=False, drop=False
    #)
    return data_summary

##################################################


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
    path_restriction = os.path.join(dock, "restriction")
    utility.create_directory(path_restriction)
    path_gene_report = os.path.join(
        path_restriction, "gene_report.pickle"
    )
    path_gene_signal = os.path.join(
        path_restriction, "data_gene_persons_tissues_signals.pickle"
    )
    # Write information to file.
    with open(path_gene_report, "wb") as file_product:
        pickle.dump(information["report_gene"], file_product)
    pandas.to_pickle(
        information["data_gene_persons_tissues_signals"], path_gene_signal
    )

    pass


###############################################################################
# Procedure


def test(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    utility.print_terminal_partition(level=1)
    print("Welcome to the test, demonstration of the restriction procedure!")

    # Read source information from file.
    source = read_source(dock=dock)
    print(source["data_gene_persons_tissues_signals"])

    # Call procedure.
    # Method for selection is either "availability" or "imputation".
    # Specific tissues are only relevant to "imputation" method.
    tissues = [
        "adipose", "blood", "colon", "esophagus", "heart", "muscle",
        "lung", "nerve", "skin", "thyroid",
    ]
    information = execute_procedure(
        method="availability",
        count=10,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            source["data_gene_persons_tissues_signals"]
        ),
    )
    print(information["report_gene"]["data_persons_tissues"])
    print(information["report_gene"]["data_persons_tissues_count"])
    print(information["report_gene"]["persons_tissues_mean"])
    print(information["report_gene"]["persons_tissues_median"])

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


def execute_procedure(
    method=None,
    count=None,
    tissues=None,
    data_gene_persons_tissues_signals=None
):
    """
    Function to execute module's main behavior.

    arguments:
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        count (int): minimal count of tissues with signal coverage for
            selection by "availability" or "imputation" methods
        tissues (list<str>): specific tissues to select by "imputation" method
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): report about gene, Pandas data frame of gene's signals across
            persons and tissues

    """



if (__name__ == "__main__"):
    execute_procedure()
