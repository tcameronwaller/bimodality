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
    path_access = os.path.join(dock, "access")
    path_attribute_sample = os.path.join(path_access, "attribute_sample.txt")
    path_attribute_patient = os.path.join(path_access, "attribute_patient.txt")
    path_selection = os.path.join(dock, "selection_small")
    path_samples = os.path.join(
        path_selection, "samples_tissues_patients.pickle"
    )
    path_gene_signal = os.path.join(path_selection, "gene_signal.pickle")
    path_tissues = os.path.join(path_selection, "tissues.pickle")
    path_patients = os.path.join(path_selection, "patients.pickle")
    path_tissues_samples = os.path.join(
        path_selection, "tissues_samples.pickle"
    )
    path_patients_samples = os.path.join(
        path_selection, "patients_samples.pickle"
    )
    # Read information from file.
    data_samples_tissues_patients = pandas.read_pickle(path_samples)
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    data_attribute_patient = pandas.read_csv(
        path_attribute_patient,
        sep="\t",
        header=0,
        #nrows=1000,
    )
    data_attribute_sample = pandas.read_csv(
        path_attribute_sample,
        sep="\t",
        header=0,
        #nrows=1000,
    )
    with open(path_tissues, "rb") as file_source:
        tissues = pickle.load(file_source)
    with open(path_patients, "rb") as file_source:
        patients = pickle.load(file_source)
    with open(path_tissues_samples, "rb") as file_source:
        tissues_patients_samples = pickle.load(file_source)
    with open(path_patients_samples, "rb") as file_source:
        patients_tissues_samples = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_samples_tissues_patients": data_samples_tissues_patients,
        "data_gene_signal": data_gene_signal,
        "data_attribute_patient": data_attribute_patient,
        "data_attribute_sample": data_attribute_sample,
        "tissues": tissues,
        "patients": patients,
        "tissues_patients_samples": tissues_patients_samples,
        "patients_tissues_samples": patients_tissues_samples
    }

##################################################
# Scrap

def calculate_gene_signal_mean_by_patients_tissues_old(
    patients_tissues_samples=None,
    data_gene_signal=None
):
    """
    Calculates the mean values of genes' signals in each tissue of each
    patient.

    Redundant samples exist for many patients and tissues.
    Calculate the mean signal for each gene from all of these redundant
    samples.

    Product data format.
    Indices by patient and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # patient tissue
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
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patient.
        data_gene_signal (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of mean values of genes' signals for all
            tissues of all patients.

    """

    # Create index for convenient access.
    data_gene_signal_index = data_gene_signal.set_index("Name")
    # Determine list of all genes with signals.
    genes = data_gene_signal["Name"].to_list()
    # Collect mean signals of all genes for each patient and tissue.
    summary = list()
    # Iterate on patients.
    for patient in patients_tissues_samples:
        # Iterate on tissues for which patient has samples.
        for tissue in patients_tissues_samples[patient]:
            # Collect signals for all genes in each sample.
            signals_genes = dict()
            for sample in patients_tissues_samples[patient][tissue]:
                for gene in genes:
                    signal = data_gene_signal_index.at[gene, sample]
                    # Determine whether an entry already exists for the gene.
                    if gene in signals_genes:
                        signals_genes[gene].append(signal)
                    else:
                        signals_genes[gene] = list([signal])
            # Compile record.
            record = dict()
            record["patient"] = patient
            record["tissue"] = tissue
            # Calculate mean signal for each gene in each patient and tissue.
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
    #    ["patient", "tissue"], append=False, drop=False
    #)
    return data_summary


def calculate_gene_signal_median_by_tissues_old(
    tissues=None,
    patients=None,
    tissues_patients_samples=None,
    data_gene_signal_mean=None
):
    """
    Calculates the median values of signals for all genes across specific
    patients and tissues.

    Redundant samples exist for many patients and tissues.
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
        patients (list<str>): Patients with signal for tissues of interest.
        tissues_patients_samples (dict<dict<list<str>>): Samples for each
            patient of each tissue.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all patients.

    raises:

    returns:
        (object): Pandas data frame of median values of genes' signals for all
            tissues.

    """

    # Create index for convenient access.
    data_gene_signal_mean_index = data_gene_signal_mean.set_index(
        ["patient", "tissue"], append=False, drop=True
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
            # Collect gene's signals in the tissue across all patients.
            signals_tissue = list()
            # Iterate on patients.
            for patient in tissues_patients_samples[tissue]:
                signal_tissue_patient = data_gene_signal_mean_index.loc[
                    (patient, tissue), gene
                ]
                signals_tissue.append(signal_tissue_patient)
            # Calculate median value of gene's signal in tissue across all
            # patients.
            signal_tissue = statistics.median(signals_tissue)
            record[tissue] = signal_tissue
        # Include record.
        summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    return data_summary


def collect_gene_signal_by_patients_tissues_old(
    tissues=None,
    patients=None,
    patients_tissues_samples=None,
    data_gene_signal_median_tissue=None,
    data_gene_signal_mean=None
):
    """
    Collects the signals for all genes in each tissue of each patient.

    Product data format.
    Indices by patient and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # patient tissue
    # bob     kidney  5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        tissues (list<str>): Tissues of interest.
        patients (list<str>): Patients with signal for tissues of interest.
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patient.
        data_gene_signal_median_tissue (object): Pandas data frame of median
            values of genes' signals for all tissues.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all patients.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific patients and tissues.

    """

    # Create indices for convenient access.
    data_gene_signal_mean_index = data_gene_signal_mean.set_index(
        ["patient", "tissue"], append=False, drop=True
    )
    data_gene_signal_median_tissue_index = (
        data_gene_signal_median_tissue.set_index("gene")
    )
    # Extract genes from names of columns.
    genes = data_gene_signal_mean_index.columns.to_list()
    # Collect signals of all genes for each patient and tissue.
    summary = list()
    # Iterate on patients of interest.
    for patient in patients:
        # Iterate on tissues of interest.
        for tissue in tissues:
            # Compile record.
            record = dict()
            record["patient"] = patient
            record["tissue"] = tissue
            # Determine whether the patient has samples for the tissue.
            if tissue in patients_tissues_samples[patient]:
                # Collect mean values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_mean_index.loc[
                        (patient, tissue), gene
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
    #    ["patient", "tissue"], append=False, drop=False
    #)
    return data_summary

##################################################


def collect_samples_tissues_patients_reference(
    data_samples_tissues_patients=None,
):
    """
    Collects tissues, patients, and genes' signals for each sample.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<str>>): Tissue and patient for each sample.

    """

    samples_tissues_patients = utility.convert_dataframe_to_records(
        data=data_samples_tissues_patients
    )
    reference = dict()
    for record in samples_tissues_patients:
        sample = record["sample"]
        tissue = record["tissue"]
        patient = record["patient"]
        if sample not in reference:
            reference[sample] = dict()
            reference[sample]["tissue"] = tissue
            reference[sample]["patient"] = patient
    return reference


def collect_genes_patients_tissues_samples(
    data_samples_tissues_patients=None,
    data_gene_signal=None
):
    """
    Collects tissues, patients, and genes' signals for each sample.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

    # Create reference for tissues and patients of each sample.
    reference = collect_samples_tissues_patients_reference(
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    def match_tissue(sample):
        return reference[sample]["tissue"]
    def match_patient(sample):
        return reference[sample]["patient"]
    # Designate patient and tissue of each sample.
    data = data_gene_signal.transpose(copy=True)
    data.reset_index(level="sample", inplace=True)
    data["tissue"] = (
        data["sample"].apply(match_tissue)
    )
    data["patient"] = (
        data["sample"].apply(match_patient)
    )
    data.set_index(
        ["patient", "tissue", "sample"],
        append=False,
        drop=True,
        inplace=True
    )
    return data


def calculate_gene_signal_mean_by_patients_tissues(
    data=None
):
    """
    Calculates the mean values of genes' signals across all samples from each
    tissue of each patient.

    Redundant samples exist for many patients and tissues.
    Calculate the mean signal for each gene from all of these redundant
    samples.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples,
            tissues, and patients.

    raises:

    returns:
        (object): Pandas data frame of mean values of genes' signals for all
            tissues of all patients.

    """

    # Define groups.
    groups = data.groupby(level=["patient", "tissue"])
    # Apply calculation to groups.
    #data_mean = groups.aggregate(lambda x: x.mean())
    data_mean = groups.aggregate(statistics.mean)
    data_mean = data_mean.reset_index(level=["sample"])
    # TODO: do I need to drop the sample column?
    return data_mean


def calculate_gene_signal_median_by_tissues(
    data=None
):
    """
    Calculates the median values of signals for all genes across specific
    patients and tissues.

    Redundant samples exist for many patients and tissues.
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
        data (object): Pandas data frame of mean values of
            genes' signals for all tissues of all patients.

    raises:

    returns:
        (object): Pandas data frame of median values of genes' signals for all
            tissues.

    """

    # Define groups.
    groups = data.groupby(level=["tissue"])
    # Apply calculation to groups.
    #data_median = groups.aggregate(lambda x: x.median())
    data_median = groups.aggregate(statistics.median())
    data_median = data_median.reset_index(level=["patient"])
    # TODO: do I need to drop the patient column?
    return data_median


def collect_gene_signal_imputation(
    tissues=None,
    patients=None,
    patients_tissues_samples=None,
    data_gene_signal_median=None,
    data_gene_signal_mean=None
):
    """
    Collects the signals for all genes in each tissue of each patient.

    Product data format.
    Indices by patient and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # patient tissue
    # bob     kidney  5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        tissues (list<str>): Tissues of interest.
        patients (list<str>): Patients with signal for tissues of interest.
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patient.
        data_gene_signal_median_tissue (object): Pandas data frame of median
            values of genes' signals for all tissues.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all patients.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific patients and tissues.

    """

    # Extract genes from names of columns.
    genes = data_gene_signal_mean.columns.to_list()
    # Collect signals of all genes for each patient and tissue.
    summary = list()
    # Iterate on patients of interest.
    for patient in patients:
        # Iterate on tissues of interest.
        for tissue in tissues:
            # Compile record.
            record = dict()
            record["patient"] = patient
            record["tissue"] = tissue
            # Determine whether the patient has samples for the tissue.
            if tissue in patients_tissues_samples[patient]:
                # Collect mean values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_mean.loc[
                        (patient, tissue), gene
                    ]
                    record[gene] = signal
            else:
                # Collect imputation values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_median.loc[
                        tissue, gene
                    ]
                    record[gene] = signal
            # Include record.
            summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    # Create indices for convenient access.
    # TODO: also name the series of genes...
    data_gene_signal_imputation = (
        data_summary.set_index(
            ["patient", "tissue"], append=False, drop=True
        )
    )
    return data_summary


########################################################
# TODO: go ahead and calculate z-scores...
########################################################


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
    path_organization = os.path.join(dock, "organization")
    utility.confirm_path_directory(path_organization)
    path_median_text = os.path.join(
        path_organization, "data_gene_signal_median.tsv"
    )
    path_imputation_text = os.path.join(
        path_organization, "data_gene_signal_imputation.tsv"
    )
    path_imputation = os.path.join(
        path_organization, "data_gene_signal_imputation.pickle"
    )
    # Write information to file.
    information["data_gene_signal_median_tissue"].to_csv(
        path_median_text, sep="\t"
    )
    information["data_gene_signal_imputation"].to_csv(
        path_imputation_text, sep="\t"
    )
    pandas.to_pickle(
        information["data_gene_signal_imputation"], path_imputation
    )


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

    # Organization of signals.
    utility.print_terminal_partition(level=1)
    print("Organization of signals.")
    print(
        "Imputation of genes' signals for some patients and tissues will " +
        "enhance coverage."
    )

    # Include identifiers for patient and tissue for each sample in data for
    # genes' signals.
    #data = source["data_gene_signal"].transpose(copy=True)
    #print(data.iloc[0:10, 0:10])
    data = collect_genes_patients_tissues_samples(
        data_samples_tissues_patients=source["data_samples_tissues_patients"],
        data_gene_signal=source["data_gene_signal"]
    )
    print(data.iloc[0:10, 0:7])

    # Calculate mean signal for each gene in all samples for each patient and
    # tissue.
    utility.print_terminal_partition(level=2)
    print(
        "Calculation of mean values of signals for each gene across " +
        "redundant samples for each patient and tissue."
    )
    print(
        "Calculate mean values before imputation to avoid bias for patients " +
        "and tissues with extra samples."
    )
    data_gene_signal_mean = calculate_gene_signal_mean_by_patients_tissues(data=data)
    print(data_gene_signal_mean.iloc[0:10, 0:7])

    # Calculation of imputation values.
    # Is it more correct to impute values before or after transformation to
    # logarithmic space?
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
    utility.print_terminal_partition(level=2)
    print("Calculation of imputation values.")
    print(
        "For imputation, use the median value of each gene's signal " +
        "within each tissue across all patients."
    )
    print(
        "Calculate imputation values before selection of specific patients " +
        "and tissues so as to derive these values from more measurements."
    )
    # The map of associations between patients, tissues, and samples does not
    # need filtration.
    # The lists of tissues and patients of interest control the selection of
    # samples.
    data_gene_signal_median = calculate_gene_signal_median_by_tissues(
        data=data_gene_signal_mean
    )
    print(data_gene_signal_median.iloc[0:25, 0:10])

    # Collection of genes' aggregate signals for all patients and tissues.
    utility.print_terminal_partition(level=2)
    print(
        "Collection of genes' signals for all patients and tissues with " +
        "imputation of misses."
    )
    data_gene_signal_imputation = collect_gene_signal_imputation(
        tissues=source["tissues"],
        patients=source["patients"],
        patients_tissues_samples=source["patients_tissues_samples"],
        data_gene_signal_median=data_gene_signal_median,
        data_gene_signal_mean=data_gene_signal_mean
    )
    print(data_gene_signal_imputation.iloc[0:25, 0:5])






    if False:
        print(data_gene_signal_mean.shape)
        data_gene_signal_mean_index = data_gene_signal_mean.set_index(
            ["patient", "tissue"], append=False, drop=True
        )
        print(data_gene_signal_mean_index.iloc[0:25, 0:5])
        # Demonstration of access to single values with multiple levels of indices.
        print("Access to values at degrees of reference.")
        print(data_gene_signal_mean_index.loc[
            ("GTEX-111CU", "Skin"), "ENSG00000078808.12"
        ])
        print(data_gene_signal_mean_index.loc[
            ("GTEX-111CU", "Muscle"), "ENSG00000078808.12"
        ])
        print(data_gene_signal_mean_index.loc[
            ("GTEX-111CU", "Nerve"), "ENSG00000078808.12"
        ])
        print(data_gene_signal_mean_index.loc["GTEX-111CU", ])
        # Extract genes from names of columns.
        #headers = data_gene_signal_mean_index.columns.to_list()
        #print(headers)
        #print("count of genes: " + str(len(headers)))

        # Reverse an index to a column.
        #dataframe["new_column"] = dataframe.index
        #dataframe.reset_index(level=["index_name"])

        # Compile information.
        information = {
            "data_gene_signal_median_tissue": data_gene_signal_median_tissue,
            "data_gene_signal_imputation": data_gene_signal_imputation,
        }
        #Write product information to file.
        write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
