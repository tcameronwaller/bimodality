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


def calculate_logarithm_gene_signal(data_signal_gene=None):
    """
    Calculates the base-2 logarithm of genes' signals in each sample.

    Original gene signals are in transcript counts per million (TPM).

    To accommodate gene signals of 0.0, add a pseudo count of 1.0 to all counts
    before calculation of the base-2 logarithm.

    arguments:
        data_signal_gene (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of base-2 logarithmic genes' signals for
            all samples.

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
    data_log = data_signal_gene.copy()
    data_log.iloc[0:, 2:] = data_log.iloc[0:, 2:].applymap(
        lambda value: math.log((value + 1.0), 2)
    )
    return data_log


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
    # Count unique patients with samples for each tissue.
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


def calculate_gene_signal_median_by_patients_tissues(
    tissues=None,
    patients=None,
    patients_tissues_samples=None,
    data_signal_gene=None
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
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patients.
        data_signal_gene (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of median signals for all genes across
            specific patients and tissues.

    """

    # Create index for convenient access.
    data_signal_gene_index = data_signal_gene.set_index("Name")
    # Determine list of all genes with signals.
    genes = data_signal_gene["Name"].to_list()
    # Compile summary.
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
            for patient in patients:
                # Determine whether the patient has samples for the tissue.
                if tissue in patients_tissues_samples[patient]:
                    # Iterate on samples for the patient and tissue.
                    signals_tissue_patient = list()
                    for sample in patients_tissues_samples[patient][tissue]:
                        # Collect signals for each sample.
                        # Calculate mean signal.
                        signal = data_signal_gene_index.at[gene, sample]
                        signals_tissue_patient.append(signal)
                    # Calculate gene's mean signal across all samples.
                    signal_tissue_patient = statistics.mean(
                        signals_tissue_patient
                    )
                    signals_tissue.append(signal_tissue_patient)
            # Calculate gene's median signal in tissue across all patients.
            signal_tissue = statistics.median(signals_tissue)
            record[tissue] = signal_tissue
        # Include record.
        summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    return data_summary


# TODO: this function might be obsolete...

def calculate_gene_signal_mean_by_patients_tissues(
    patients_tissues_samples=None,
    data_signal_gene=None
):
    """
    Calculates the mean values of signals for all genes in each tissue of each
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
    # april   brain   2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patients.
        data_signal_gene (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of mean signals for all genes in all
            tissues of all patients.

    """

    # Create index for convenient access.
    data_signal_gene_index = data_signal_gene.set_index("Name")
    # Determine list of all genes with signals.
    genes = data_signal_gene["Name"].to_list()
    # Compile summary.
    summary = list()
    # Iterate on patients.
    for patient in patients_tissues_samples:
        # Iterate on tissues for which patient has samples.
        for tissue in patients_tissues_samples[patient]:
            # Collect signals for all genes in each sample.
            samples_genes_signals = dict()
            for sample in patients_tissues_samples[patient][tissue]:
                for gene in genes:
                    signal = data_signal_gene_index.at[gene, sample]
                    # Determine whether an entry already exists for the gene.
                    if gene in samples_genes_signals:
                        samples_genes_signals[gene].append(signal)
                    else:
                        samples_genes_signals[gene] = list([signal])
            # Compile entry.
            entry = dict()
            entry["patient"] = patient
            entry["tissue"] = tissue
            # Calculate mean signal for each gene in tissue.
            for gene in genes:
                mean = statistics.mean(samples_genes_signals[gene])
                entry[gene] = mean
            # Include entry.
            summary.append(entry)
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    data_summary_index = data_summary.set_index(
        ["patient", "tissue"], append=False, drop=False
    )
    return data_summary_index


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


def collect_tissues_coverage_patients(
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
    tissues_coverage_patients = list()
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
        tissues_coverage_patients.append(record)
    # Create a data frame and sort by counts.
    data = utility.convert_records_to_dataframe(
        records=tissues_coverage_patients
    )
    data_sort = data.sort_values(
        "misses", axis=0, ascending=True
    )
    return data_sort


def filter_tissues_coverage_patients(
    threshold=None, data_tissues_coverage_patients=None
):
    """
    Filters patients by coverage of specific tissues.

    Includes counts of tissues of interest for which each patient misses
    signal.

    arguments:
        threshold (int): Count of permissible tissues which a patient can miss.
        data_tissues_coverage_patients (object): Pandas data frame of counts of
            tissues for which each patient misses signal.

    raises:

    returns:
        (list<str>): Patients with signal for tissues of interest.

    """

    tissues_coverage_patients = utility.convert_dataframe_to_records(
        data=data_tissues_coverage_patients
    )
    patients = list()
    for record in tissues_coverage_patients:
        patient = record["patient"]
        misses = record["misses"]
        # Determine whether the patient has few enough misses.
        if misses <= threshold:
            patients.append(patient)
    return patients


def collect_gene_signal_by_patients_tissues(
    tissues=None,
    patients=None,
    data_signal_gene_tissue_median=None,
    patients_tissues_samples=None,
    data_signal_gene=None
):
    """
    Calculates the mean values of signals for all genes in each tissue of each
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
    # april   brain   2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        patients_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each patients.
        data_signal_gene (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of mean signals for all genes in all
            tissues of all patients.

    """

    # Create indices for convenient access.
    data_signal_gene_index = data_signal_gene.set_index("Name")
    data_signal_gene_tissue_median_index = (
        data_signal_gene_tissue_median.set_index("gene")
    )
    # Determine list of all genes with signals.
    genes = data_signal_gene["Name"].to_list()
    # Compile summary.
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
                # Collect signals for all genes across samples.
                signals_samples = dict()
                # Iterate on samples.
                for sample in patients_tissues_samples[patient][tissue]:
                    # Iterate on genes.
                    for gene in genes:
                        signal = data_signal_gene_index.at[gene, sample]
                        # Determine whether an entry already exists for the
                        # gene.
                        if gene in signals_samples:
                            signals_samples[gene].append(signal)
                        else:
                            signals_samples[gene] = list([signal])
                # Collect mean signals for each gene across samples.
                # Iterate on genes.
                for gene in genes:
                    signal_mean = statistics.mean(signals_samples[gene])
                    # Compile record.
                    record[gene] = signal_mean
            else:
                # Collect imputation values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = (
                        data_signal_gene_tissue_median_index.at[gene, tissue]
                    )
                    # Compile record.
                    record[gene] = signal
            # Include record.
            summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    data_summary_index = data_summary.set_index(
        ["patient", "tissue"], append=False, drop=False
    )
    return data_summary_index


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
    #tissues_sets = define_test_tissue_sets()

    tissues = [
        "Skin", "Muscle", "Adipose Tissue", "Esophagus", "Thyroid", "Lung",
        "Blood", "Nerve", "Heart", "Colon",
    ]
    print(tissues)

    # Collect patients with gene signals for all of a combination of tissues.
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
        "Collection of patients with gene signals for tissues of interest."
    )
    print("Allow for fewer than 3 missing tissues for each patient.")
    print("Hence the extent of imputation will be less than 33%.")
    data_tissues_coverage_patients = collect_tissues_coverage_patients(
        tissues=tissues,
        patients_tissues_samples=patients_tissues_samples
    )
    #print(data_tissues_coverage_patients)
    #print(data_signal_gene_mean.loc["GTEX-12WSD", ])
    # Filter for those patients which miss less than 3 tissues of interest.
    summary = (
        data_tissues_coverage_patients.loc[
            data_tissues_coverage_patients["misses"] < 3
        ]
    )
    print(summary)
    print(summary.shape)
    # Filter patients.
    patients = filter_tissues_coverage_patients(
        threshold=2,
        data_tissues_coverage_patients=data_tissues_coverage_patients
    )
    print(patients)
    print("count of patients: " + str(len(patients)))

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

    # Organization of signals for tissues and patients of interest.
    utility.print_terminal_partition(level=1)
    print("Organization of signals.")
    print(
        "To enhance coverage, it might be necessary to impute values of " +
        "genes' signals for some patients and tissues."
    )

    print("**********")
    print(
        "I think it is most correct to transform genes' signals to " +
        "logarithmic space before calculation of median or mean values."
    )
    print("**********")

    # Transform genes' signals to base-2 logarithmic space.
    # Transform before calculation of any median or mean values.
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
        data_signal_gene=source["data_signal_gene"]
    )
    print(data_signal_gene_log.iloc[0:10, 0:10])

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
    print("**********")
    print(
        "I am unsure whether it is more accurate to impute values before " +
        "or after transformation to logarithmic space."
    )
    print("**********")
    # The map of associations between patients, tissues, and samples does not
    # need filtration.
    # The lists of tissues and patients of interest control the selection of
    # samples.
    data_signal_gene_tissue_median = (
        calculate_gene_signal_median_by_patients_tissues(
            tissues=tissues,
            patients=patients,
            patients_tissues_samples=patients_tissues_samples,
            data_signal_gene=data_signal_gene_log
        )
    )
    print(data_signal_gene_tissue_median.iloc[0:25, :])

    # Collection of genes' aggregate signals for all patients and tissues.
    utility.print_terminal_partition(level=2)
    print(
        "Collection of genes' aggregate signals for all patients and " +
        "tissues across redundant samples and with imputation of misses."
    )
    data_signal_gene_aggregation = collect_gene_signal_by_patients_tissues(
        tissues=tissues,
        patients=patients,
        data_signal_gene_tissue_median=data_signal_gene_tissue_median,
        patients_tissues_samples=patients_tissues_samples,
        data_signal_gene=data_signal_gene_log
    )
    print(data_signal_gene_aggregation.iloc[0:25, 0:5])






    ##################################################
    ##################################################
    ##################################################


    utility.print_terminal_partition(level=2)






    # Calculate mean signal for each gene in all samples for each patient and
    # tissue.
    utility.print_terminal_partition(level=2)
    print(
        "Calculation of mean values of signals for each gene across " +
        "redundant samples for each patient and tissue."
    )
    data_signal_gene_mean = calculate_gene_signal_mean_by_patients_tissues(
        patients_tissues_samples=patients_tissues_samples,
        data_signal_gene=data_signal_gene_log
    )
    print(data_signal_gene_mean.iloc[0:25, 0:5])
    # Demonstration of access to single values with multiple levels of indices.
    print("Access to values at degrees of reference.")
    print(data_signal_gene_mean.loc[
        ("GTEX-111CU", "Skin"), "ENSG00000078808.12"
    ])
    print(data_signal_gene_mean.loc[
        ("GTEX-111CU", "Muscle"), "ENSG00000078808.12"
    ])
    print(data_signal_gene_mean.loc[
        ("GTEX-111CU", "Nerve"), "ENSG00000078808.12"
    ])
    print(data_signal_gene_mean.loc["GTEX-111CU", ])


    # TODO: I need to figure out how to impute values for missing tissue measurements...







    #####################################################################
    ####################################################################


    utility.print_terminal_partition(level=1)

    print("current new stuff...")

    # Now consider coverage of patients and tissues for the sake of filtering...



    ######################################################################



    # TODO: Now collect counts of patients in bins for each count of tissues (1-21).


    ####################################################
    ######################################################
    ####################################################3


    # Determine whether data have any missing values.
    print("original data frame shape: ")
    print(source["data_signal_gene"].shape)
    print("shape without missing axis 0: ")
    print(source["data_signal_gene"].dropna(axis=0).shape)
    print("shape without missing axis 1: ")
    print(source["data_signal_gene"].dropna(axis=1).shape)
    # Explore data a little...

    # Experiment with dataframe.iterrows() and dataframe.itertuples()
    # data.itertuples(index=True)

    # How to reference information about a specific sample... use pandas.loc()?

    print("Select row")
    print(source["data_attribute_sample"].loc[source["data_attribute_sample"]["SAMPID"] == "GTEX-1117F-0003-SM-58Q7G"])
    print("select column from row...")
    print(source["data_attribute_sample"].loc[source["data_attribute_sample"]["SAMPID"] == "GTEX-1117F-0003-SM-58Q7G"]["SMTS"][0])
    print(source["data_attribute_sample"].loc[source["data_attribute_sample"]["SAMPID"] == "GTEX-1117F-0003-SM-58Q7G"].at[0, "SMTS"])


    ####################################################
    ######################################################
    ####################################################3

    ##################################################




    # TODO: Scrap below here...


    if False:

        counter = 0
        for row in data.iterrows():
            # apply function here to determine if the sample passes criteria
            # I assume that I can access all of the gene values this way...
            # if not, then consider using iterrows
            if counter < 3:
                print(row.Description)
            counter += 1

    pass


if (__name__ == "__main__"):
    execute_procedure()
