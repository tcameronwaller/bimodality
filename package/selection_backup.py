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
import gc

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


# TODO: when ready, read in all data...

def read_source(dock=None):
    """
    Reads and organizes source information from file.

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
    # Read information from file.
    data_gene_signal = pandas.read_csv(
        path_signal_gene,
        sep="\t",
        header=2,
        #nrows=100,
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
    # Compile and return information.
    return {
        "data_attribute_patient": data_attribute_patient,
        "data_attribute_sample": data_attribute_sample,
        "data_gene_signal": data_gene_signal,
    }


def collect_samples_tissues_patients(
    data_gene_signal=None, data_attribute_sample=None
):
    """
    Collects matches of samples, tissues, and patients.

    Product data format.
    Indices by patient and tissue allow convenient access.
    ##################################################
    # sample   patient   tissue
    # sm_1     bob       brain
    # sm_2     bob       heart
    # sm_3     bob       liver
    # sm_4     julie     brain
    # sm_5     julie     kidney
    # sm_6     betty     pancreas
    ##################################################

    arguments:
        data_gene_signal (object): Pandas data frame of gene's signals for all
            samples.
        data_attribute_sample (object): Pandas data frame of attributes for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of patients and tissues for all samples.

    """

    # Determine samples with signals for genes.
    # Extract samples' identifiers from data of gene's signals.
    # Extract names of columns.
    # Pandas series
    #headers = data_gene_signal.columns
    # NumPy array
    #headers = data_gene_signal.columns.values
    # List
    headers = data_gene_signal.columns.to_list()
    # Exclude name and description.
    identifiers_sample = headers[2:]

    # Create index for samples' attributes by samples' identifiers.
    data_attribute_sample_index = data_attribute_sample.set_index("SAMPID")

    # Collect tissues and patients for each sample.
    samples_tissues_patients = list()
    for identifier_sample in identifiers_sample:
        #tissue = data_attribute_sample.loc[
        #    data_attribute_sample["SAMPID"] == identifier_sample,
        #].at[0, "SMTS"]
        tissue = data_attribute_sample_index.at[identifier_sample, "SMTS"]
        identifier_patient = utility.extract_gtex_sample_patient_identifier(
            identifier_sample
        )
        record = {
            "sample": identifier_sample,
            "tissue": tissue,
            "patient": identifier_patient,
        }
        samples_tissues_patients.append(record)
    return utility.convert_records_to_dataframe(
        records=samples_tissues_patients
    )


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


def collect_patients_tissues_samples(data_samples_tissues_patients=None):
    """
    Collect hierarchical structure of patients, tissues, and samples.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<list<str>>): Samples for each tissue of each patient.

    """

    samples_tissues_patients = utility.convert_dataframe_to_records(
        data=data_samples_tissues_patients
    )
    # Collect unique tissues and samples for each patient.
    patients_tissues_samples = dict()
    for record in samples_tissues_patients:
        sample = record["sample"]
        tissue = record["tissue"]
        patient = record["patient"]
        # Determine whether an entry already exists for the patient.
        if patient in patients_tissues_samples:
            # Determine whether an entry already exists for the tissue.
            if tissue in patients_tissues_samples[patient]:
                patients_tissues_samples[patient][tissue].append(sample)
            else:
                patients_tissues_samples[patient][tissue] = list([sample])
        else:
            patients_tissues_samples[patient] = dict()
            patients_tissues_samples[patient][tissue] = list([sample])
    return patients_tissues_samples


def collect_tissues_patients_samples(data_samples_tissues_patients=None):
    """
    Collect hierarchical structure of tissues, patients, and samples.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<list<str>>): Samples for each patient of each tissue.

    """

    samples_tissues_patients = utility.convert_dataframe_to_records(
        data=data_samples_tissues_patients
    )
    # Collect unique patients and samples for each tissue.
    tissues_patients_samples = dict()
    for record in samples_tissues_patients:
        sample = record["sample"]
        tissue = record["tissue"]
        patient = record["patient"]
        # Determine whether an entry already exists for the tissue.
        if tissue in tissues_patients_samples:
            # Determine whether an entry already exists for the patient.
            if patient in tissues_patients_samples[tissue]:
                tissues_patients_samples[tissue][patient].append(sample)
            else:
                tissues_patients_samples[tissue][patient] = list([sample])
        else:
            tissues_patients_samples[tissue] = dict()
            tissues_patients_samples[tissue][patient] = list([sample])
    return tissues_patients_samples


def expand_print_patients_tissues_samples(patients_tissues_samples=None):
    """
    Collects tissues and samples for each patient.

    arguments:
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patients.

    raises:

    returns:


    """

    print(list(patients_tissues_samples.keys()))
    for patient in patients_tissues_samples:
        print("patient: " + patient)
        for tissue in patients_tissues_samples[patient]:
            print("tissue: " + tissue)
            for sample in patients_tissues_samples[patient][tissue]:
                print("sample: " + sample)
    pass


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
    utility.create_directory(path_selection)
    path_tissues = os.path.join(path_selection, "tissues.pickle")
    path_patients = os.path.join(path_selection, "patients.pickle")
    path_tissues_samples = os.path.join(
        path_selection, "tissues_samples.pickle"
    )
    path_patients_samples = os.path.join(
        path_selection, "patients_samples.pickle"
    )
    # Write information to file.
    with open(path_tissues, "wb") as file_product:
        pickle.dump(information["tissues"], file_product)
    with open(path_patients, "wb") as file_product:
        pickle.dump(information["patients"], file_product)
    with open(path_tissues_samples, "wb") as file_product:
        pickle.dump(information["tissues_patients_samples"], file_product)
    with open(path_patients_samples, "wb") as file_product:
        pickle.dump(information["patients_tissues_samples"], file_product)



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

    # Enable automatic garbage collection to clear memory.
    gc.enable()

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
    print(source["data_gene_signal"].iloc[0:10, 0:10])
    #print(source["data_gene_signal"].loc[:,"GTEX-1117F-0226-SM-5GZZ7"])

    # Summarize table of patients' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of patients' attributes.")
    print(source["data_attribute_patient"].iloc[0:10, 0:10])

    # Summarize table of samples' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of samples' attributes.")
    print(source["data_attribute_sample"].iloc[0:10, 0:7])

    ##################################################
    ##################################################
    ##################################################

    # TODO: filter out samples with values of 0 for all genes...
    # TODO: do so before parsing patients -> tissues -> samples
    # TODO: do so before selecting patients and tissues.

    #print(source["data_gene_signal"].loc[:,"GTEX-1117F-0226-SM-5GZZ7"])

    # Check for data quality.
    utility.print_terminal_partition(level=1)
    print("Check for quality of genes' signals.")

    # Temporarily drop column for genes' descriptions.
    utility.print_terminal_partition(level=2)
    print("Drop column for genes' description.")
    data_gene_signal = source["data_gene_signal"].drop(
        labels="Description",
        axis="columns",
    )
    print(data_gene_signal.iloc[0:10, 0:10])

    # Organize data with names and indices.
    # The Pandas function to rename axis copies data deeply by default and
    # requires a lot of memory.
    utility.print_terminal_partition(level=2)
    print("Organize data with names and indices.")
    data_gene_signal = data_gene_signal.set_index("Name")
    data_gene_signal_axis_temporary = data_gene_signal_index.rename_axis(
        index="gene",
        axis="index",
        copy=False,
    )
    data_gene_signal_axis = data_gene_signal_axis_temporary.rename_axis(
        columns="sample",
        axis="columns",
        copy=False,
    )
    print(data_gene_signal_axis.iloc[0:10, 0:10])

    # Collect garbage to clear memory.
    gc.collect()

    # Check for redundant genes.
    utility.print_terminal_partition(level=2)
    print("Check for redundant genes in genes' signals.")
    print("Consider names of genes.")
    # Reset indices to consider names of genes.
    data_gene_signal.reset_index()
    print(data_gene_signal.iloc[0:10, 0:10])
    data_redundancy = (
        data_gene_signal.duplicated(subset=None, keep="first")
    )
    data_redundancy_list = data_redundancy.to_list()
    if any(data_redundancy_list):
        print("Redundancy in genes: Yes")
    else:
        print("Redundancy in genes: No")

    # Check for missing values of genes' signals.
    utility.print_terminal_partition(level=2)
    print("Check for missing values in genes' signals.")
    print(
        "shape of original data frame: " + str(data_gene_signal_axis.shape)
    )
    print(
        "shape without missing axis 0: " +
        str(data_gene_signal_axis.dropna(axis=0).shape)
    )
    print(
        "shape without missing axis 1: " +
        str(data_gene_signal_axis.dropna(axis=1).shape)
    )

    # Collect garbage to clear memory.
    gc.collect()

    # Check for samples with values of 0 for all genes' signals.
    utility.print_terminal_partition(level=2)
    print("Check for samples with values of 0 for all genes' signals.")
    print(
        "shape of original data frame: " +
        str(data_gene_signal_axis.shape)
    )
    data_nonzero = (
        data_gene_signal_axis
        .loc[ : , (data_gene_signal_axis != 0).any(axis="index")]
    )
    print(
        "shape of data frame without zero samples: " + str(data_nonzero.shape)
    )

    # Collect garbage to clear memory.
    gc.collect()

    # Check for genes with values of 0 for signals across all samples.
    utility.print_terminal_partition(level=2)
    print("Check for genes with values of 0 for signals across all samples.")
    print("These genes are undetectable.")
    print(
        "shape of original data frame: " +
        str(data_gene_signal_axis.shape)
    )
    data_nonzero = (
        data_gene_signal_axis
        .loc[(data_gene_signal_axis != 0).any(axis="columns"), : ]
    )
    print(
        "shape of data frame without zero samples: " + str(data_nonzero.shape)
    )

    # Collect garbage to clear memory.
    gc.collect()



    ##################################################
    ##################################################
    ##################################################

    # Summarize organization of data.
    utility.print_terminal_partition(level=1)
    print("Hierarchical organization of samples.")

    # Organize samples by hierarchy of patients and tissues.
    # Collect tissues and patients for each sample.
    utility.print_terminal_partition(level=2)
    print("Collection of tissues and patients for each sample.")
    print("Extract patient identifiers from sample identifiers.")
    #GTEX-1117F-0226-SM-5GZZ7
    identifier = (
        utility.extract_gtex_sample_patient_identifier(
            "GTEX-1117F-0226-SM-5GZZ7"
        )
    )
    print("GTEX-1117F-0226-SM-5GZZ7: " + identifier)
    data_samples_tissues_patients = collect_samples_tissues_patients(
        data_gene_signal=source["data_gene_signal"],
        data_attribute_sample=source["data_attribute_sample"]
    )
    print(data_samples_tissues_patients.iloc[0:10, :])
    print(data_samples_tissues_patients.shape)

    ##################################################
    ##################################################
    ##################################################


    # Organization of samples by patients and tissues.
    # Collect unique patients, unique tissues for each patient, and unique
    # samples for each tissue of each patient.
    utility.print_terminal_partition(level=2)
    print("Organization of samples by hierarchy of patients and tissues.")
    print("Collection of hierarchical groups by patient, tissue, and sample.")
    utility.print_terminal_partition(level=4)
    print(
        "data hierarchy is " +
        "patients (714) -> tissues (30) -> samples (11688) -> genes (~20,000)"
    )
    patients_tissues_samples = collect_patients_tissues_samples(
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    if False:
        expand_print_patients_tissues_samples(
            patients_tissues_samples=patients_tissues_samples
        )
    print("Printing first 10 unique patients...")
    print(list(patients_tissues_samples.keys())[:9])
    print("Printing groups for patient 'GTEX-14LZ3'... ")
    print(patients_tissues_samples["GTEX-14LZ3"])

    # Organization of samples by tissues and patients.
    # Collect unique tissues, unique patients for each tissue, and unique
    # samples for each patient of each tissue.
    utility.print_terminal_partition(level=2)
    print("Organization of samples by hierarchy of tissues and patients.")
    print("Collection of hierarchical groups by tissue, patient, and sample.")
    utility.print_terminal_partition(level=4)
    print(
        "data hierarchy is " +
        "tissue (30) -> patients (714) -> samples (11688) -> genes (~20,000)"
    )
    tissues_patients_samples = collect_tissues_patients_samples(
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    print("Printing first 10 unique tissues...")
    print(list(tissues_patients_samples.keys())[:9])
    print("Printing groups for tissue 'Bladder'... ")
    print(tissues_patients_samples["Bladder"])

    ##################################################
    ##################################################
    ##################################################

    # Selection of tissues and patients of interest.
    utility.print_terminal_partition(level=1)
    print("Selection of tissues and patients of interest.")

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

    ##################################################
    ##################################################
    ##################################################

    # Compile information.
    information = {
        "patients_tissues_samples": patients_tissues_samples,
        "tissues_patients_samples": tissues_patients_samples,
        "tissues": tissues,
        "patients": patients,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
