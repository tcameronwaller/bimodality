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
    path_assembly = os.path.join(dock, "assembly")
    path_samples_tissues_patients = os.path.join(
        path_assembly, "data_samples_tissues_patients.pickle"
    )
    path_patients_tissues_samples = os.path.join(
        path_assembly, "patients_tissues_samples.pickle"
    )
    path_tissues_patients_samples = os.path.join(
        path_assembly, "tissues_patients_samples.pickle"
    )
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_samples_tissues_patients = pandas.read_pickle(
        path_samples_tissues_patients
    )
    patients_tissues_samples = pandas.read_pickle(
        path_patients_tissues_samples
    )
    tissues_patients_samples = pandas.read_pickle(
        path_tissues_patients_samples
    )
    data_gene_annotation = pandas.read_pickle(
        path_gene_annotation
    )
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_samples_tissues_patients": data_samples_tissues_patients,
        "patients_tissues_samples": patients_tissues_samples,
        "tissues_patients_samples": tissues_patients_samples,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_signal": data_gene_signal,
    }


##########
# Tissues and patients.

def select_tissues_patients(
    data_samples_tissues_patients=None,
    patients_tissues_samples=None,
    tissues_patients_samples=None
):
    """
    Selects tissues and patients of interest for further analyses.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patient.
        tissues_patients_samples (dict<dict<list<str>>): Samples for each
            patient of each tissue.

    raises:

    returns:
        (dict<list<str>>): Information about tissues and patients.

    """

    utility.print_terminal_partition(level=1)
    print("Selection of tissues and patients of interest.")

    #

    # Count patients with samples for each tissue.
    # How many patients have samples for each tissue?
    ##################################################
    # tissue patients
    # brain   50
    # liver   10
    # heart   5
    ##################################################
    utility.print_terminal_partition(level=2)
    print("Counts of patients with samples for each tissue.")
    data_tissues_patients_counts = count_tissues_patients(
        filter=False,
        tissues=[],
        patients=[],
        tissues_patients_samples=tissues_patients_samples
    )
    print(data_tissues_patients_counts)

    # Count tissues for which each patient has samples.
    ##################################################
    # patient tissues
    # alice   3
    # john    5
    # paul    10
    # mary    13
    ##################################################
    utility.print_terminal_partition(level=2)
    print("Counts of tissues for which each patient has samples.")
    data_patients_tissues_counts = count_patients_tissues(
        filter=False,
        tissues=[],
        patients_tissues_samples=patients_tissues_samples
    )
    print(data_patients_tissues_counts)

    # How many patients have samples for each count of tissues?
    # Populate bins of patients with each count of tissues.
    # Plot histogram.
    ##################################################
    # tissues patients
    # 3       50
    # 5       10
    # 7       5
    ##################################################
    utility.print_terminal_partition(level=2)
    print("Counts of patients with samples for each count of tissues.")
    data_patients_tissues_counts_bins = count_bins_patients_tissues(
        data_patients_tissues_counts=data_patients_tissues_counts
    )
    print(data_patients_tissues_counts_bins)

    ##################################################
    ##################################################
    ##################################################

    # Define tissues of interest.
    utility.print_terminal_partition(level=2)
    print("Selection of tissues of interest.")
    print("Call these tissues the selection tissues.")
    #tissues_sets = define_test_tissue_sets()
    tissues = [
        "Skin", "Muscle", "Adipose Tissue", "Esophagus", "Thyroid", "Lung",
        "Blood", "Nerve", "Heart", "Colon",
    ]
    tissues_count = len(tissues)
    print(tissues)
    print("count of selection tissues: " + str(tissues_count))

    # Collect patients with samples for all of the selection tissues.
    utility.print_terminal_partition(level=2)
    print(
        "Counts of selection tissues for which each patient has samples."
    )
    data_patients_tissues_coverage = count_patients_tissues(
        filter=True,
        tissues=tissues,
        patients_tissues_samples=patients_tissues_samples
    )
    print(data_patients_tissues_coverage)
    #print(data_gene_signal_mean.loc["GTEX-12WSD", ])
    # Filter for those patients which miss fewer than 4 selection tissues.
    print("Allow each patient to miss fewer than 4 selection tissues.")
    print(
        "Hence the extent of imputation will be less than 40% of tissues " +
        "for each patient."
    )
    misses = 3
    matches = tissues_count - misses
    summary_patients_coverage = (
        data_patients_tissues_coverage.loc[
            data_patients_tissues_coverage["count"] >= matches
        ]
    )
    print(summary_patients_coverage)
    print(
        "count of selection patients: " +
        str(summary_patients_coverage.shape[0])
    )

    # Filter for patients with adequate coverage of tissues of interest.
    patients = filter_patients_tissues_coverage(
        threshold=matches,
        data_patients_tissues_counts=data_patients_tissues_coverage
    )
    print(patients)
    print("count of selection patients: " + str(len(patients)))
    print("Call these patients the selection patients.")

    # Collect selection patients that have samples for each tissue.
    # Summarize the extent of imputation necessary for each tissue.
    ##################################################
    # tissue  patients
    # brain   125
    # liver   175
    # kidney  150
    # lung    255
    # skin    300
    ##################################################
    utility.print_terminal_partition(level=2)
    print(
        "Collection of selection patients with samples for each selection " +
        "tissue."
    )
    print(
        "Summary of extent of imputation necessary for each tissue."
    )
    data_tissues_patients_coverage = count_tissues_patients(
        filter=True,
        tissues=tissues,
        patients=patients,
        tissues_patients_samples=tissues_patients_samples
    )
    print(data_tissues_patients_coverage)

    # Compile information.
    information = {
        "tissues": tissues,
        "patients": patients
    }

    # Return information.
    return information


def filter_samples_tissues_patients(
    tissues=None,
    patients=None,
    data_samples_tissues_patients=None
):
    """
    Filters samples by tissues and patients.

    arguments:
        tissues (list<str>): Tissues of interest.
        patients (list<str>): Patients with signal for tissues of interest.
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:
        (object): Pandas data frame of patients and tissues for all samples.

    """

    samples_tissues_patients = utility.convert_dataframe_to_records(
        data=data_samples_tissues_patients
    )
    # Determine whether each sample passes filters by tissues and patients.
    samples_tissues_patients_filter = list()
    for record in samples_tissues_patients:
        sample = record["sample"]
        tissue = record["tissue"]
        patient = record["patient"]
        if (patient in patients) and (tissue in tissues):
            samples_tissues_patients_filter.append(record)
    return utility.convert_records_to_dataframe(
        records=samples_tissues_patients_filter
    )


def count_tissues_patients(
    filter=None, tissues=None, patients=None, tissues_patients_samples=None
):
    """
    Counts patients with samples for each tissue.

    arguments:
        filter (bool): Whether to filter the tissues and patients for counts.
        tissues (list<str>): Tissues of interest.
        patients (list<str>): Patients of interest.
        tissues_patients_samples (dict<dict<list<str>>): Samples for each
            patient of each tissue.

    raises:

    returns:
        (object): Pandas data frame of counts of patients with samples for each
            tissue.

    """

    # Count unique patients with samples for each tissue.
    coverage = list()
    # Determine whether to filter tissues.
    if filter:
        tissues_relevant = tissues
    else:
        tissues_relevant = list(tissues_patients_samples.keys())
    # Iterate on tissues.
    for tissue in tissues_relevant:
        # Collect unique patients with samples for tissue.
        patients_tissue = list(tissues_patients_samples[tissue].keys())
        # Determine whether to filter patients.
        if filter:
            patients_relevant = utility.filter_common_elements(
                list_one=patients, list_two=patients_tissue
            )
        else:
            patients_relevant = patients_tissue
        # Count patients.
        patients_count = len(patients_relevant)
        # Compile record.
        record = dict()
        record["tissue"] = tissue
        record["patients"] = patients_relevant
        record["count"] = patients_count
        # Include record.
        coverage.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(records=coverage)
    data_sort = data.sort_values("count", axis=0, ascending=False)
    return data_sort


def count_patients_tissues(
    filter=None, tissues=None, patients_tissues_samples=None
):
    """
    Counts tissues for which each patient has samples.

    arguments:
        filter (bool): Whether to filter the tissues and patients for counts.
        tissues (list<str>): Tissues of interest.
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patient.

    raises:

    returns:
        (object): Pandas data frame of counts of tissues for which each patient
            has samples.

    """

    # Count unique tissues for which each patient has samples.
    coverage = list()
    # Iterate on patients.
    for patient in patients_tissues_samples:
        # Collect unique tissues for which patient has samples.
        tissues_patient = list(patients_tissues_samples[patient].keys())
        # Determine whether to filter tissues.
        if filter:
            tissues_relevant = utility.filter_common_elements(
                list_one=tissues, list_two=tissues_patient
            )
        else:
            tissues_relevant = tissues_patient
        # Count tissues.
        tissues_count = len(tissues_relevant)
        # Compile record.
        record = dict()
        record["patient"] = patient
        record["tissues"] = tissues_relevant
        record["count"] = tissues_count
        # Include record.
        coverage.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(records=coverage)
    data_sort = data.sort_values("count", axis=0, ascending=False)
    return data_sort


def count_bins_patients_tissues(data_patients_tissues_counts=None):
    """
    Count patients in each bin by counts of tissues.

    arguments:
        data_patients_tissues_counts (object): Pandas data frame of counts of
            tissues for which each patient has samples.

    raises:

    returns:
        (object): Pandas data frame of counts of patients in each bin by counts
            of tissues.

    """

    patients_tissues_counts = utility.convert_dataframe_to_records(
        data=data_patients_tissues_counts
    )
    # Compile summary.
    bins = dict()
    # Iterate on patients.
    for record in patients_tissues_counts:
        patient = record["patient"]
        count = record["count"]
        # Determine whether an entry already exists for the bin.
        if count in bins:
            bins[count] += 1
        else:
            bins[count] = 1
    return bins


def define_test_tissue_sets():
    """
    Defines sets of tissue for testing.

    arguments:

    raises:

    returns:
        (dict<list<str>>): Lists of tissues.

    """

    # tissue combination: counts of patients which miss fewer than 3 tissues
    # 1: 11
    # 2: 48
    # 3: 73
    # 4: 77
    # 5: 95
    # 6: 100
    # 7: 110
    # 8: 135
    # 9: 137
    # 10: 141
    # 11: 194
    # 12: 204
    # 13: 228
    # 14: 230
    # 15: 245
    # 16: 291
    # 17: 294
    # 18: 296
    # 19: 299
    # 20: 312
    # 21: 346
    # 22: 391
    tissues = dict()
    tissues[1] = [
        "Skin", "Muscle", "Adipose Tissue", "Blood Vessel", "Esophagus",
        "Thyroid", "Lung", "Blood", "Nerve", "Heart", "Colon", "Stomach",
        "Brain", "Pancreas", "Adrenal Gland", "Pituitary", "Liver"
    ]
    tissues[2] = [
        "Skin", "Muscle", "Adipose Tissue", "Esophagus", "Thyroid", "Lung",
        "Blood", "Nerve", "Heart", "Colon", "Stomach", "Brain", "Pancreas",
        "Liver"
    ]
    tissues[3] = [
        "Skin", "Muscle", "Adipose Tissue", "Esophagus", "Thyroid", "Lung",
        "Blood", "Nerve", "Heart", "Colon", "Stomach", "Brain", "Liver"
    ]
    tissues[4] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Blood",
        "Nerve", "Heart", "Colon", "Stomach", "Brain", "Liver"
    ]
    tissues[5] = [
        "Skin", "Muscle", "Adipose Tissue", "Blood Vessel", "Esophagus",
        "Thyroid", "Lung", "Blood", "Nerve", "Heart", "Colon", "Stomach",
        "Brain", "Pancreas"
    ]
    tissues[6] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Blood",
        "Nerve", "Heart", "Colon", "Stomach", "Brain", "Pancreas"
    ]
    tissues[7] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Nerve",
        "Heart", "Colon", "Stomach", "Brain", "Liver"
    ]
    tissues[8] = [
        "Skin", "Muscle", "Adipose Tissue", "Blood Vessel", "Esophagus",
        "Thyroid", "Lung", "Blood", "Nerve", "Heart", "Colon", "Stomach",
        "Brain",
    ]
    tissues[9] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Nerve",
        "Heart", "Colon", "Stomach", "Brain", "Pancreas",
    ]
    tissues[10] = [
        "Skin", "Muscle", "Adipose Tissue", "Esophagus", "Thyroid", "Lung",
        "Blood", "Nerve", "Heart", "Colon", "Stomach", "Brain",
    ]

    tissues[11] = [
        "Skin", "Muscle", "Adipose Tissue", "Blood Vessel", "Esophagus",
        "Thyroid", "Lung", "Blood", "Nerve", "Heart", "Colon", "Brain"
    ]
    # 200 - 250
    tissues[12] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Nerve",
        "Heart", "Colon", "Stomach", "Brain",
    ]
    tissues[13] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Nerve",
        "Heart", "Colon", "Stomach", "Pancreas",
    ]
    tissues[14] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Blood",
        "Nerve", "Heart", "Colon", "Brain"
    ]
    tissues[15] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung",
        "Blood", "Nerve", "Heart", "Colon", "Stomach"
    ]

    # 250 - 300
    tissues[16] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Nerve",
        "Heart", "Colon", "Brain"
    ]
    tissues[17] = [
        "Skin", "Muscle", "Adipose Tissue", "Blood Vessel", "Esophagus",
        "Thyroid", "Lung", "Blood", "Nerve", "Heart", "Colon",
    ]
    tissues[18] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Nerve",
        "Heart", "Colon", "Pancreas"
    ]
    tissues[19] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Nerve",
        "Heart", "Colon", "Stomach"
    ]
    # 300 - 400
    # Favorite.
    tissues[20] = [
        "Skin", "Muscle", "Adipose Tissue", "Esophagus", "Thyroid", "Lung",
        "Blood", "Nerve", "Heart", "Colon",
    ]
    tissues[21] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Blood",
        "Nerve", "Heart", "Colon",
    ]
    tissues[22] = [
        "Skin", "Muscle", "Adipose Tissue", "Thyroid", "Lung", "Nerve",
        "Heart", "Colon",
    ]

    return tissues


def filter_patients_tissues_coverage(
    threshold=None, data_patients_tissues_counts=None
):
    """
    Filters patients by coverage of selection tissues.

    Includes counts of tissues of interest for which each patient misses
    signal.

    arguments:
        threshold (int): Count of selection tissues which a patient must have.
        data_patients_tissues_counts (object): Pandas data frame of counts of
            tissues for which each patient has samples.

    raises:

    returns:
        (list<str>): Patients with signal for selection tissues.

    """

    patients_tissues_counts = utility.convert_dataframe_to_records(
        data=data_patients_tissues_counts
    )
    patients = list()
    for record in patients_tissues_counts:
        patient = record["patient"]
        count = record["count"]
        # Determine whether the patient has few enough misses.
        if count >= threshold:
            patients.append(patient)
    return patients


##########
# Genes.


def select_samples(
    tissues=None,
    patients=None,
    data_gene_signal=None
):
    """
    Selects samples of interest for further analyses.

    arguments:
        tissues (list<str>): Tissues of interest.
        patients (list<str>): Patients of interest.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and patients.

    raises:

    returns:
        (dict): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

    # Select samples from patients and tissues of interest.
    utility.print_terminal_partition(level=2)
    print("Selection of samples from patients and tissues of interest.")
    print("count of samples, original: " + str(data_gene_signal.shape[0]))
    data_gene_signal.reset_index(
        level=["patient", "tissue", "sample"], inplace=True
    )
    data_gene_signal.set_index(
        ["patient"],
        append=False,
        drop=True,
        inplace=True
    )
    data_gene_signal = data_gene_signal.loc[patients, : ]
    print(
        "count of samples from patients of interest: " +
        str(data_gene_signal.shape[0])
    )
    data_gene_signal.reset_index(level=["patient"], inplace=True)
    data_gene_signal.set_index(
        ["tissue"],
        append=False,
        drop=True,
        inplace=True
    )
    data_gene_signal = data_gene_signal.loc[tissues, : ]
    print(
        "count of samples from tissues of interest: " +
        str(data_gene_signal.shape[0])
    )
    data_gene_signal.reset_index(level=["tissue"], inplace=True)
    data_gene_signal.set_index(
        ["patient", "tissue", "sample"],
        append=False,
        drop=True,
        inplace=True
    )

    return data_gene_signal


def select_genes_detection(data_gene_signal=None):
    """
    Selects detectable genes with nonzero signals.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and patients.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

    utility.print_terminal_partition(level=2)
    print(
        "Selection of detectable genes with nonzero signals in patients " +
        "and tissues of interest."
    )
    print("genes, original: " + str(data_gene_signal.shape[1]))
    data_nonzero = (data_gene_signal != 0)
    data_signal = data_gene_signal.loc[ : , data_nonzero.any(axis="index")]
    print("genes, detection: " + str(data_signal.shape[1]))
    return data_signal


def select_samples_genes(
    patients=None,
    tissues=None,
    data_gene_annotation=None,
    data_gene_signal=None
):
    """
    Selects samples and genes of interest for further analyses.

    arguments:
        patients (list<str>): Patients of interest.
        tissues (list<str>): Tissues of interest.
        data_gene_annotation (object): Pandas data frame of genes' annotations.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and patients.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

    utility.print_terminal_partition(level=1)
    print("Selection of samples and genes of interest.")

    # Select samples from patients and tissues of interest.
    data_gene_signal = select_samples(
        patients=patients,
        tissues=tissues,
        data_gene_signal=data_gene_signal,
    )

    # Select genes with detectable, non-zero signal in tissues and patients of
    # interest.
    data_gene_signal = select_genes_detection(
        data_gene_signal=data_gene_signal
    )

    # Select genes that encode proteins.
    data_gene_signal = select_genes_protein(
        data_gene_annotation=data_gene_annotation,
        data_gene_signal=data_gene_signal
    )

    # Return information.
    return data_gene_signal


##########
# Product.


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
    path_selection = os.path.join(dock, "selection")
    utility.confirm_path_directory(path_selection)
    path_tissues = os.path.join(path_selection, "tissues.pickle")
    path_patients = os.path.join(path_selection, "patients.pickle")
    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    # Write information to file.
    with open(path_tissues, "wb") as file_product:
        pickle.dump(information["tissues"], file_product)
    with open(path_patients, "wb") as file_product:
        pickle.dump(information["patients"], file_product)
    pandas.to_pickle(
        information["data_gene_signal"],
        path_gene_signal
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

    # Read source information from file.
    source = read_source(dock=dock)

    ##################################################
    ##################################################
    ##################################################

    # Selection of tissues and patients of interest.
    tissues_patients = select_tissues_patients(
        data_samples_tissues_patients=source["data_samples_tissues_patients"],
        patients_tissues_samples=source["patients_tissues_samples"],
        tissues_patients_samples=source["tissues_patients_samples"]
    )

    ##################################################
    ##################################################
    ##################################################

    # TODO:
    # TODO:
    # TODO:
    # 1. Divide tissues into 2 groups
    # - 1.1. major tissues without minor categories
    # -- select minor categories of tissues that are not cell lines specifically for skin and blood
    # - 1.2. major tissues with minor categories
    # -- select whether to combine the minor categories...

    # TODO:
    # TODO:
    # TODO:



    # Selection of samples and genes of interest.
    data_gene_signal = select_samples_genes(
        tissues=tissues_patients["tissues"],
        patients=tissues_patients["patients"],
        data_gene_annotation=source["data_gene_annotation"],
        data_gene_signal=source["data_gene_signal"]
    )

    ##################################################
    ##################################################
    ##################################################

    utility.print_terminal_partition(level=2)
    print(
        "Summary of signals."
    )
    print(data_gene_signal.iloc[0:10, 0:10])
    print(data_gene_signal.shape)

    # Compile information.
    information = {
        "tissues": tissues_patients["tissues"],
        "patients": tissues_patients["patients"],
        "data_gene_signal": data_gene_signal,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)


    # Scrap.

    if False:
        # Summarize tissues and patients.
        # Extract unique values of tissue.
        tissues_test = list(source["data_gene_signal_imputation"]["tissue"].unique())
        patients_test = list(source["data_gene_signal_imputation"]["patient"].unique())


    pass


if (__name__ == "__main__"):
    execute_procedure()
