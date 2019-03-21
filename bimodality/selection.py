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
    # Compile and return information.
    return {
        "data_attribute_patient": data_attribute_patient,
        "data_attribute_sample": data_attribute_sample,
        "data_signal_gene": data_signal_gene,
    }


def collect_samples_tissues_patients(
    data_signal_gene=None, data_attribute_sample=None
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
        data_signal_gene (object): Pandas data frame of gene's signals for all
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
    #headers = data_signal_gene.columns
    # NumPy array
    #headers = data_signal_gene.columns.values
    # List
    headers = data_signal_gene.columns.to_list()
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
    Collects tissues and samples for each patient.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<list<str>>): Samples for each tissue of each patients.

    """

    samples_tissues_patients = utility.convert_dataframe_to_records(
        data=data_samples_tissues_patients
    )
    # Collect tissues and samples for each patient.
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


def collect_tissues_patients(data_samples_tissues_patients=None):
    """
    Collects patients with samples for each tissue.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<int>>): Counts of samples for each tissue from each patient.

    """

    samples_tissues_patients = utility.convert_dataframe_to_records(
        data=data_samples_tissues_patients
    )
    # Collect unique patients with samples for each tissue.
    tissues_patients = dict()
    for record in samples_tissues_patients:
        sample = record["sample"]
        tissue = record["tissue"]
        patient = record["patient"]
        # Determine whether an entry already exists for the tissue.
        if tissue in tissues_patients:
            # Determine whether an entry already exists for the patient.
            if patient in tissues_patients[tissue]:
                tissues_patients[tissue][patient] += 1
            else:
                tissues_patients[tissue][patient] = 1
        else:
            tissues_patients[tissue] = dict()
            tissues_patients[tissue][patient] = 1
    return tissues_patients


def count_tissues_patients(tissues_patients=None):
    """
    Counts patients with samples for each tissue.

    arguments:
        tissues_patients (dict<dict<int>>): Counts of samples for each tissue
            from each patient.

    raises:

    returns:
        (object): Pandas data frame of counts of patients with samples for each
            tissue.

    """

    # Count unique patients with samples for each tissue.
    # Compile summary.
    summary = list()
    # Iterate on tissues.
    for tissue in tissues_patients:
        # Count patients.
        patients = len(list(tissues_patients[tissue].keys()))
        # Compile entry.
        entry = dict()
        entry["tissue"] = tissue
        entry["patients"] = patients
        # Include entry.
        summary.append(entry)
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    data_summary_sort = data_summary.sort_values(
        "patients", axis=0, ascending=False
    )
    return data_summary_sort


def count_patients_tissues(patients_tissues_samples=None):
    """
    Counts tissues for which each patient has samples.

    arguments:
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patient.

    raises:

    returns:
        (object): Pandas data frame of counts of tissues for which each patient
            has samples.

    """

    # Count unique tissues for which each patient has samples.
    # Compile summary.
    summary = list()
    # Iterate on patients.
    for patient in patients_tissues_samples:
        # Count tissues.
        tissues = len(list(patients_tissues_samples[patient].keys()))
        # Compile entry.
        entry = dict()
        entry["patient"] = patient
        entry["tissues"] = tissues
        # Include entry.
        summary.append(entry)
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    data_summary_sort = data_summary.sort_values(
        "tissues", axis=0, ascending=False
    )
    return data_summary_sort


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
    summary = dict()
    # Iterate on patients.
    for record in patients_tissues_counts:
        patient = record["patient"]
        tissues = record["tissues"]
        # Determine whether an entry already exists for the bin.
        if tissues in summary:
            summary[tissues] += 1
        else:
            summary[tissues] = 1
    return summary


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


def collect_coverage_patients(
    tissues=None, patients_tissues_samples=None
):
    """
    Collects patients with coverage of signals for specific tissues.

    Includes counts of tissues of interest for which each patient misses
    signal.

    arguments:
        tissues (list<str>): Tissues of interest.
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patient.

    raises:

    returns:
        (object): Pandas data frame of counts of tissues for which each patient
            misses signal.

    """

    # Collect tissues and samples for each patient.
    coverage_patients = list()
    # Iterate on patients.
    for patient in patients_tissues_samples:
        # Count tissues.
        tissues_patient = list(patients_tissues_samples[patient].keys())
        # Determine count of tissues for which patient misses signal.
        misses = list()
        for tissue in tissues:
            if tissue not in tissues_patient:
                misses.append(tissue)
        misses_count = len(misses)
        record = {
            "patient": patient,
            "misses": misses_count,
        }
        coverage_patients.append(record)
    # Create a data frame and sort by counts.
    data = utility.convert_records_to_dataframe(
        records=coverage_patients
    )
    data_sort = data.sort_values(
        "misses", axis=0, ascending=True
    )
    return data_sort


def collect_coverage_tissues(
    tissues=None, patients=None, data_samples_tissues_patients=None
):
    """
    Collects patients with coverage of signals for specific tissues.

    Includes counts of tissues of interest for which each patient misses
    signal.

    arguments:
        tissues (list<str>): Tissues of interest.
        patients (list<str>): Patients of interest.
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:
        (object): Pandas data frame of counts of patients with samples for each
            tissue.

    """

    # Collect unique patients with samples for each tissue.
    tissues_patients = collect_tissues_patients(
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    # Count unique patients with samples for each tissue.
    coverage_tissues = list()
    # Iterate on tissues.
    for tissue in tissues:
        # Filter patients.
        tissue_patients = list(tissues_patients[tissue].keys())
        tissue_patients_filter = utility.filter_common_elements(
            list_one=patients, list_two=tissue_patients
        )
        count_tissue_patients = len(tissue_patients_filter)
        # Compile record.
        record = dict()
        record["tissue"] = tissue
        record["patients"] = tissue_patients_filter
        record["count"] = count_tissue_patients
        # Include record.
        coverage_tissues.append(record)
    data_coverage_tissues = utility.convert_records_to_dataframe(
        records=coverage_tissues
    )
    data_coverage_tissues_sort = data_coverage_tissues.sort_values(
        "count", axis=0, ascending=False
    )
    return data_coverage_tissues_sort


def filter_tissues_coverage_patients(
    threshold=None, data_coverage_patients=None
):
    """
    Filters patients by coverage of specific tissues.

    Includes counts of tissues of interest for which each patient misses
    signal.

    arguments:
        threshold (int): Count of permissible tissues which a patient can miss.
        data_coverage_patients (object): Pandas data frame of counts of
            tissues for which each patient misses signal.

    raises:

    returns:
        (list<str>): Patients with signal for tissues of interest.

    """

    tissues_coverage_patients = utility.convert_dataframe_to_records(
        data=data_coverage_patients
    )
    patients = list()
    for record in tissues_coverage_patients:
        patient = record["patient"]
        misses = record["misses"]
        # Determine whether the patient has few enough misses.
        if misses <= threshold:
            patients.append(patient)
    return patients


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

    ##################################################
    ##################################################
    ##################################################

    # Summarize organization of data.
    utility.print_terminal_partition(level=1)
    print("Selection of tissues and patients of interest.")

    # Collect samples by their patients and tissues.
    # Collect matches of samples, tissues, and patients.
    utility.print_terminal_partition(level=2)
    print("Collection of patients and tissues that correspond to each sample.")
    print(
        "data hierarchy is " +
        "patients (714) -> tissues (?) -> samples (?) -> genes (~20,000)"
    )
    print("Extract patient identifiers from sample identifiers.")
    #GTEX-1117F-0226-SM-5GZZ7
    identifier = (
        utility.extract_gtex_sample_patient_identifier(
            "GTEX-1117F-0226-SM-5GZZ7"
        )
    )
    print("GTEX-1117F-0226-SM-5GZZ7: " + identifier)
    data_samples_tissues_patients = collect_samples_tissues_patients(
        data_signal_gene=source["data_signal_gene"],
        data_attribute_sample=source["data_attribute_sample"]
    )
    print(data_samples_tissues_patients.iloc[0:10, :])

    # Collect unique patients, unique tissues for each patient, and unique
    # samples for each tissue of each patient.
    utility.print_terminal_partition(level=2)
    print("Collection of hierarchical groups by patient, tissue, and sample.")
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

    # Summarize counts of patients for each tissue.
    # How many patients have samples for each tissue?
    ##################################################
    # tissue patients
    # brain   50
    # liver   10
    # heart   5
    ##################################################
    utility.print_terminal_partition(level=2)
    print("Counts of patients with samples for each tissue.")
    tissues_patients = collect_tissues_patients(
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    #print(tissues_patients)
    #print(tissues_patients["Blood"])
    tissues_patients_counts = count_tissues_patients(
        tissues_patients=tissues_patients
    )
    print(tissues_patients_counts)

    # Summarize counts of patients for each count of tissues.
    # How many patients have samples for each count of tissues?
    # Count tissues for each patient.
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
        patients_tissues_samples=patients_tissues_samples
    )
    print(data_patients_tissues_counts)

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
    patients_tissues_counts_bins = count_bins_patients_tissues(
        data_patients_tissues_counts=data_patients_tissues_counts
    )
    print(patients_tissues_counts_bins)

    # Define tissues of interest.
    utility.print_terminal_partition(level=2)
    print("Define tissues of interest.")
    print("Call these tissues the selection tissues.")
    #tissues_sets = define_test_tissue_sets()

    tissues = [
        "Skin", "Muscle", "Adipose Tissue", "Esophagus", "Thyroid", "Lung",
        "Blood", "Nerve", "Heart", "Colon",
    ]
    print(tissues)

    # Collect patients with samples for all of a combination of tissues.
    # Collect information about how many tissues each patient is missing.
    ##################################################
    # patient missing_tissues
    # alice   0
    # john    0
    # paul    0
    # mary    1
    # jim     1
    # susan   1
    # phil    2
    # rachel  3
    ##################################################
    utility.print_terminal_partition(level=2)
    print(
        "Collection of patients with adequate coverage of samples for " +
        "tissues of interest."
    )
    print("Allow for fewer than 4 missing tissues for each patient.")
    print(
        "Hence the extent of imputation will be less than 40% of tissues " +
        "for each patient."
    )
    print("Also consider extent of imputation for each tissue.")
    data_coverage_patients = collect_coverage_patients(
        tissues=tissues,
        patients_tissues_samples=patients_tissues_samples
    )
    #print(data_coverage_patients)
    #print(data_signal_gene_mean.loc["GTEX-12WSD", ])
    # Filter for those patients which miss less than 4 tissues of interest.
    summary_coverage_patients = (
        data_coverage_patients.loc[data_coverage_patients["misses"] < 4]
    )
    print(summary_coverage_patients)
    print(
        "count of patients: " + str(summary_coverage_patients.shape[0])
    )

    # Filter for patients with adequate coverage of tissues of interest.

    patients = filter_tissues_coverage_patients(
        threshold=3,
        data_coverage_patients=data_coverage_patients
    )
    print(patients)
    print("count of patients: " + str(len(patients)))
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
    data_coverage_tissues = collect_coverage_tissues(
        tissues=tissues,
        patients=patients,
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    print("**************************")
    print(data_coverage_tissues)






    # Filter patients, tissues, and samples.
    # This step is unnecessary.
    if False:
        data_samples_tissues_patients_filter = filter_samples_tissues_patients(
            patients=patients,
            tissues=tissues,
            data_samples_tissues_patients=data_samples_tissues_patients
        )
        patients_tissues_samples_filter = collect_patients_tissues_samples(
            data_samples_tissues_patients=data_samples_tissues_patients_filter
        )
        if False:
            expand_print_patients_tissues_samples(
                patients_tissues_samples=patients_tissues_samples
            )
        print("Printing first 10 unique patients...")
        print(list(patients_tissues_samples_filter.keys())[:9])
        print("Printing groups for patient 'GTEX-131XE'... ")
        print(patients_tissues_samples_filter["GTEX-131XE"])

    ##################################################
    ##################################################
    ##################################################

    pass


if (__name__ == "__main__"):
    execute_procedure()
