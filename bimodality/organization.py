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
    path_signal_gene = os.path.join(path_access, "signal_gene.gct")
    path_selection = os.path.join(dock, "selection")
    path_tissues = os.path.join(path_selection, "tissues.pickle")
    path_patients = os.path.join(path_selection, "patients.pickle")
    path_tissues_samples = os.path.join(
        path_selection, "tissues_samples.pickle"
    )
    path_patients_samples = os.path.join(
        path_selection, "patients_samples.pickle"
    )
    # Read information from file.
    data_signal_gene = pandas.read_csv(
        path_signal_gene,
        sep="\t",
        header=2,
        nrows=100,
    )
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
        "data_attribute_patient": data_attribute_patient,
        "data_attribute_sample": data_attribute_sample,
        "data_signal_gene": data_signal_gene,
        "tissues": tissues,
        "patients": patients,
        "tissues_patients_samples": tissues_patients_samples,
        "patients_tissues_samples": patients_tissues_samples
    }


def calculate_gene_signal_mean_by_patients_tissues(
    patients_tissues_samples=None,
    data_signal_gene=None
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
        data_signal_gene (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of mean values of genes' signals for all
            tissues of all patients.

    """

    # Create index for convenient access.
    data_signal_gene_index = data_signal_gene.set_index("Name")
    # Determine list of all genes with signals.
    genes = data_signal_gene["Name"].to_list()
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
                    signal = data_signal_gene_index.at[gene, sample]
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


def calculate_gene_signal_median_by_tissues(
    tissues=None,
    patients=None,
    tissues_patients_samples=None,
    data_signal_gene_mean=None
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
        data_signal_gene_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all patients.

    raises:

    returns:
        (object): Pandas data frame of median values of genes' signals for all
            tissues.

    """

    # Create index for convenient access.
    data_signal_gene_mean_index = data_signal_gene_mean.set_index(
        ["patient", "tissue"], append=False, drop=True
    )
    # Extract genes from names of columns.
    genes = data_signal_gene_mean_index.columns.to_list()
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
                signal_tissue_patient = data_signal_gene_mean_index.loc[
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


def collect_gene_signal_by_patients_tissues(
    tissues=None,
    patients=None,
    patients_tissues_samples=None,
    data_signal_gene_tissue_median=None,
    data_signal_gene_mean=None
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
        data_signal_gene_tissue_median (object): Pandas data frame of median
            values of genes' signals for all tissues.
        data_signal_gene_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all patients.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific patients and tissues.

    """

    # Create indices for convenient access.
    data_signal_gene_mean_index = data_signal_gene_mean.set_index(
        ["patient", "tissue"], append=False, drop=True
    )
    data_signal_gene_tissue_median_index = (
        data_signal_gene_tissue_median.set_index("gene")
    )
    # Extract genes from names of columns.
    genes = data_signal_gene_mean_index.columns.to_list()
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
                    signal = data_signal_gene_mean_index.loc[
                        (patient, tissue), gene
                    ]
                    record[gene] = signal
            else:
                # Collect imputation values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_signal_gene_tissue_median_index.loc[
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


def calculate_logarithm_gene_signal(data_signal_gene=None):
    """
    Calculates the base-2 logarithm of genes' signals in each sample.

    Original gene signals are in transcript counts per million (TPM).

    To accommodate gene signals of 0.0, add a pseudo count of 1.0 to all counts
    before calculation of the base-2 logarithm.

    arguments:
        data_signal_gene (object): Pandas data frame of signals for all genes
            across specific patients and tissues.

    raises:

    returns:
        (object): Pandas data frame of base-2 logarithmic signals for all genes
            across specific patients and tissues.

    """

    # lambda x: math.log((x + 1), 2)

    # An alternative approach would be to set the label columns to indices.
    if False:
        data_signal_gene_index = data_signal_gene.set_index(
            ["Name", "Description"], append=True, drop=True
        )
        data_signal_gene_log = data_signal_gene_index.apply(
            lambda value: 2 * value
        )
        data_log = data_signal_index.copy()
        data_log.iloc[0:, 2:] = data_log.iloc[0:, 2:].applymap(
            lambda value: math.log((value + 1.0), 2)
        )

    data_signal_gene_index = data_signal_gene.set_index(
        ["patient", "tissue"], append=False, drop=True
    )
    data_log = data_signal_gene_index.applymap(
        lambda value: math.log((value + 1.0), 2)
    )
    return data_log


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

    # Samples beginning with code patient "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    ##################################################
    ##################################################
    ##################################################

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of genes' signals in all samples.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' signals in all samples.")
    print("Genes' signals are in Transcripts Per Million (TPM).")
    print(source["data_signal_gene"].iloc[0:10, 0:10])
    #print(source["data_signal_gene"].loc[:,"GTEX-1117F-0226-SM-5GZZ7"])

    # Summarize table of patients' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of patients' attributes.")
    print(source["data_attribute_patient"].iloc[0:10, 0:10])

    # Summarize table of samples' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of samples' attributes.")
    print(source["data_attribute_sample"].iloc[0:10, 0:7])

    # Determine whether data have any missing values.
    print("Determine whether data have any missing values.")
    print(
        "original data frame shape: " + str(source["data_signal_gene"].shape)
    )
    print(
        "shape without missing axis 0: " +
        str(source["data_signal_gene"].dropna(axis=0).shape)
    )
    print(
        "shape without missing axis 1: " +
        str(source["data_signal_gene"].dropna(axis=1).shape)
    )

    # How to reference information about a specific sample.
    if False:
        print(
            source["data_attribute_sample"].loc[
                (
                    source["data_attribute_sample"]["SAMPID"] ==
                    "GTEX-1117F-0003-SM-58Q7G"
                )
            ]
        )

    ##################################################
    ##################################################
    ##################################################

    # Organization of signals.
    utility.print_terminal_partition(level=1)
    print("Organization of signals.")
    print(
        "Imputation of genes' signals for some patients and tissues will " +
        "enhance coverage."
    )

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
    data_signal_gene_mean = calculate_gene_signal_mean_by_patients_tissues(
        patients_tissues_samples=source["patients_tissues_samples"],
        data_signal_gene=source["data_signal_gene"]
    )
    print(data_signal_gene_mean.shape)
    data_signal_gene_mean_index = data_signal_gene_mean.set_index(
        ["patient", "tissue"], append=False, drop=True
    )
    print(data_signal_gene_mean_index.iloc[0:25, 0:5])
    # Demonstration of access to single values with multiple levels of indices.
    print("Access to values at degrees of reference.")
    print(data_signal_gene_mean_index.loc[
        ("GTEX-111CU", "Skin"), "ENSG00000078808.12"
    ])
    print(data_signal_gene_mean_index.loc[
        ("GTEX-111CU", "Muscle"), "ENSG00000078808.12"
    ])
    print(data_signal_gene_mean_index.loc[
        ("GTEX-111CU", "Nerve"), "ENSG00000078808.12"
    ])
    print(data_signal_gene_mean_index.loc["GTEX-111CU", ])
    # Extract genes from names of columns.
    headers = data_signal_gene_mean_index.columns.to_list()
    print(headers)
    print("count of genes: " + str(len(headers)))

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
    data_signal_gene_tissue_median = calculate_gene_signal_median_by_tissues(
        tissues=source["tissues"],
        patients=source["patients"],
        tissues_patients_samples=source["tissues_patients_samples"],
        data_signal_gene_mean=data_signal_gene_mean
    )
    print(data_signal_gene_tissue_median.iloc[0:25, :])

    # Collection of genes' aggregate signals for all patients and tissues.
    utility.print_terminal_partition(level=2)
    print(
        "Collection of genes' aggregate signals for all patients and " +
        "tissues across redundant samples and with imputation of misses."
    )
    data_signal_gene_aggregation = collect_gene_signal_by_patients_tissues(
        tissues=source["tissues"],
        patients=source["patients"],
        patients_tissues_samples=source["patients_tissues_samples"],
        data_signal_gene_tissue_median=data_signal_gene_tissue_median,
        data_signal_gene_mean=data_signal_gene_mean
    )
    data_signal_gene_aggregation_index = (
        data_signal_gene_aggregation.set_index(
            ["patient", "tissue"], append=False, drop=True
        )
    )
    print(data_signal_gene_aggregation_index.iloc[0:25, 0:5])

    # Transform genes' signals to base-2 logarithmic space.
    # Transform before calculation of any median or mean values. <- maybe false
    # Logarithms are not distributive.
    # To accommodate gene signals of 0.0, add a pseudo count of 1.0 to all
    # counts before calculation of the base-2 logarithm.
    # log-2 signal = log2(TPM + pseudo-count)
    # pseudo-count = 1.0
    utility.print_terminal_partition(level=2)
    print("Transformation of genes' signals to base-2 logarithmic space.")
    print(
        "To accommodate gene signals of 0.0, add a pseudo count of 1.0 to " +
        "all counts before calculation of the base-2 logarithm."
    )
    data_signal_gene_log = calculate_logarithm_gene_signal(
        data_signal_gene=data_signal_gene_aggregation
    )
    print(data_signal_gene_log.iloc[0:10, 0:10])

    # Reverse an index to a column.
    #dataframe["new_column"] = dataframe.index
    #dataframe.reset_index(level=["index_name"])


    ##################################################
    ##################################################
    ##################################################

if (__name__ == "__main__"):
    execute_procedure()
