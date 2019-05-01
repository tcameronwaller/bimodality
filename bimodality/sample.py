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

##########
# Source


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
    path_gene_count = os.path.join(
        path_assembly, "data_gene_count.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_samples_tissues_patients = pandas.read_pickle(
        path_samples_tissues_patients
    )
    data_gene_count = pandas.read_pickle(path_gene_count)
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_samples_tissues_patients": data_samples_tissues_patients,
        "data_gene_count": data_gene_count,
        "data_gene_signal": data_gene_signal,
    }


##########
# Patients

# patients female: 257
# patients male: 494


# Consider the distribution of patients by sex... and age?



##########
# Tissues


def collect_categories_tissues_patients_samples(
    data_samples_tissues_patients=None
):
    """
    Collect hierarchical structure of tissues, patients, and samples.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<list<str>>): Samples for each patient of each tissue.

    """

    # Organize data.
    data = data_samples_tissues_patients.copy(deep=True)
    data.reset_index(
        level=["sample"], inplace=True
    )
    samples_tissues_patients = utility.convert_dataframe_to_records(
        data=data
    )
    # Collect unique patients and samples for each tissue.
    tissues_patients_samples = dict()
    for record in samples_tissues_patients:
        sample = record["sample"]
        tissue_major = record["tissue_major"]
        tissue_minor = record["tissue_minor"]
        patient = record["patient"]
        # Determine whether an entry already exists for the major tissue.
        if tissue_major in tissues_patients_samples:
            # Determine whether an entry already exists for the minor tissue.
            if tissue_minor in tissues_patients_samples[tissue_major]:
                # Determine whether an entry already exists for the patient.
                if patient in tissues_patients_samples[tissue_major][tissue_minor]:
                    tissues_patients_samples[tissue_major][tissue_minor][patient].append(sample)
                else:
                    tissues_patients_samples[tissue_major][tissue_minor][patient] = list([sample])
            else:
                tissues_patients_samples[tissue_major][tissue_minor] = dict()
                tissues_patients_samples[tissue_major][tissue_minor][patient] = list([sample])
        else:
            tissues_patients_samples[tissue_major] = dict()
            tissues_patients_samples[tissue_major][tissue_minor] = dict()
            tissues_patients_samples[tissue_major][tissue_minor][patient] = list([sample])
    return tissues_patients_samples


def count_categories_tissues_patients(
    tissues_patients_samples=None
):
    """
    Counts patients with samples for each tissue.

    arguments:
        tissues_patients_samples (dict<dict<list<str>>): Samples for each
            patient of each tissue.

    raises:

    returns:
        (object): Pandas data frame of counts of patients with samples for each
            tissue.

    """

    # Count unique patients with samples for each tissue.
    coverage = list()
    # Access major tissues.
    tissues_major = list(tissues_patients_samples.keys())
    # Iterate on major tissues.
    for tissue_major in tissues_major:
        # Access minor tissues.
        tissues_minor = list(tissues_patients_samples[tissue_major].keys())
        # Iterate on minor tissues.
        for tissue_minor in tissues_minor:
            # Count unique patients with samples for the tissue.
            patients_tissue = list(
                tissues_patients_samples[tissue_major][tissue_minor].keys()
            )
            count_patients = len(patients_tissue)
            # Compile record.
            record = dict()
            record["tissue_major"] = tissue_major
            record["tissue_minor"] = tissue_minor
            record["patients"] = count_patients
            # Include record.
            coverage.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(records=coverage)
    data_sort = data.sort_values(
        ["tissue_major", "tissue_minor"],
        axis=0,
        ascending=True
    )
    return data_sort


def describe_samples_tissues(
    data_samples_tissues_patients=None
):
    """
    Describe the distribution of samples across tissues.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:

    """

    # Print terminal partition.
    utility.print_terminal_partition(level=2)
    # Report.
    print("Description of tissues of samples.")

    # Collect categories of samples.
    tissues_patients_samples = collect_categories_tissues_patients_samples(
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    # Report.
    if False:
        utility.print_terminal_partition(level=3)
        print("Printing unique major tissues...")
        print(list(tissues_patients_samples.keys())[:])
        utility.print_terminal_partition(level=3)
        print("Printing unique minor tissues for major tissue 'Brain'... ")
        print(list(tissues_patients_samples["Brain"].keys())[:])

    # Count categories of samples.
    data_tissues_patients_counts = count_categories_tissues_patients(
        tissues_patients_samples=tissues_patients_samples
    )
    # Report.
    utility.print_terminal_partition(level=3)
    print("Printing counts of samples in unique categories.")
    print(data_tissues_patients_counts)

    pass


def describe_samples_categories(
    data_samples_tissues_patients=None
):
    """
    Describe the distribution of samples across tissues.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples.

    raises:

    returns:

    """

    # Explore samples by patient attributes.
    # Describe distributions of patients by sex and age.

    # Explore samples by tissue attributes.
    # Describe categories of all samples.
    describe_samples_tissues(
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    # Describe categories of specific samples.
    # Print terminal partition.
    utility.print_terminal_partition(level=1)
    # Report.
    print("Description of categories of samples for male patients.")
    data_samples_male = data_samples_tissues_patients.loc[
        data_samples_tissues_patients["sex"] == "male", :
    ]
    #print(data_samples_male)
    describe_samples_tissues(
        data_samples_tissues_patients=data_samples_male
    )
    # Print terminal partition.
    utility.print_terminal_partition(level=1)
    # Report.
    print("Description of categories of samples for female patients.")
    data_samples_female = data_samples_tissues_patients.loc[
        data_samples_tissues_patients["sex"] == "female", :
    ]
    #print(data_samples_female)
    describe_samples_tissues(
        data_samples_tissues_patients=data_samples_female
    )

    # Describe similarity between minor categories of tissues.
    # Print terminal partition.
    utility.print_terminal_partition(level=1)
    # Report.
    print("Description of similarity between minor tissues.")
    print("Relevant major tissues, asexual with > 135 patients:")
    print(
        "Thyroid, Stomach, Spleen, Intestine, Skin, Pituitary, Pancreas, " +
        "Nerve, Muscle, Lung, Liver, Heart, Esophagus, Colon, Brain, " +
        "Vessel, Blood, Adrenal, Adipose"
    )
    utility.print_terminal_partition(level=3)
    print("Remove data for cell lines, relevant to Skin and Blood.")
    utility.print_terminal_partition(level=3)
    print("Relevant major tissues with minor categories:")
    print("skin, heart, esophagus, colon, brain, artery, adipose")

    pass



####################################################
# Move to Selection procedure... ?


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

########################################################


def define_tissue_comparisons():
    """
    Defines minor categories to compare for each major category of tissue.

    arguments:

    raises:

    returns:
        (dict<list<str>>): Minor categories to compare for major categories of
            tissue

    """

    comparisons = dict()
    comparisons["adipose"] = ["subcutane", "viscera"]
    comparisons["artery"] = ["aorta", "coronary", "tibial"]
    comparisons["brain"] = [
        "accumbens", "amygdala", "caudate", "cerebellum", "cingulate",
        "cortex", "frontal", "hemisphere", "hippocampus", "hypothalamus",
        "nigra", "putamen", "spine",
    ]
    comparisons["colon"] = ["sigmoid", "transverse"]
    comparisons["esophagus"] = ["junction", "mucosa", "muscularis"]
    comparisons["heart"] = ["atrium", "ventricle"]
    comparisons["skin"] = ["dark", "light"]
    return comparisons


def organize_differential_expression_data_sets(
    data_samples_tissues_patients=None,
    data_gene_count=None,
):
    """
    Collect hierarchical structure of tissues, patients, and samples.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples
        data_gene_count (object): Pandas data frame of genes' counts for all
            samples

    raises:

    returns:
        (list<dict>): Collections of data sets for differential expression
            analyses

    """

    # Collect data sets.
    sets = list()
    comparisons = define_tissue_comparisons()
    comparisons_list = list(comparisons.keys())
    for comparison in comparisons_list:
        data = organize_differential_expression_data_set(
            tissue_major=comparison,
            tissues_minor=comparisons[comparison],
            data_samples_tissues_patients=data_samples_tissues_patients,
            data_gene_count=data_gene_count,
        )
        # Collect the data set.
        sets.append(data)
    return sets


def organize_differential_expression_data_set(
    data_samples_tissues_patients=None,
    data_gene_count=None,
):
    """
    Collect hierarchical structure of tissues, patients, and samples.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples
        data_gene_count (object): Pandas data frame of genes' counts for all
            samples

    raises:

    returns:
        (list<dict>): Collections of data sets for differential expression
            analyses

    """

    return 0




##########
# Product


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
    path_sample = os.path.join(dock, "sample")
    utility.confirm_path_directory(path_sample)
    path_patients_tissues_samples = os.path.join(
        path_sample, "patients_tissues_samples.pickle"
    )
    path_tissues_patients_samples = os.path.join(
        path_sample, "tissues_patients_samples.pickle"
    )

    # Write information to file.
    with open(path_patients_tissues_samples, "wb") as file_product:
        pickle.dump(information["patients_tissues_samples"], file_product)
    with open(path_tissues_patients_samples, "wb") as file_product:
        pickle.dump(information["tissues_patients_samples"], file_product)
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

    print(source["data_samples_tissues_patients"].iloc[0:10, :])
    print(source["data_gene_count"].iloc[0:10, 0:10])

    # Describe coverage of tissues by samples.
    describe_samples_categories(
        data_samples_tissues_patients=source["data_samples_tissues_patients"]
    )

    # Organize data sets for comparison by differential expression of genes.
    # 1. Define samples relevant for comparison of minor categories of tissues.
    # 2. Organize a matrix to designate groups for each sample.
    # - Filter "data_samples_tissues_patients" and organize.
    ###########################################################################
    # sample                   patient    tissue
    # GTEX-1117F-0226-SM-5GZZ7 GTEX-1117F subcutane
    # GTEX-111CU-1826-SM-5GZYN GTEX-111CU viscera
    # GTEX-111FC-0226-SM-5N9B8 GTEX-111FC subcutane
    # GTEX-111YS-2426-SM-5GZZQ GTEX-111YS viscera
    # GTEX-1122O-2026-SM-5NQ91 GTEX-1122O subcutane
    ###########################################################################
    # 3. Filter data of gene's signals to include only relevant samples.
    # - Filter "data_gene_count".
    ###########################################################################
    # sample          GTEX-1117F-...  GTEX-111CU-...  GTEX-111FC-...
    # gene
    # ENSG00000186092             53             125             534
    # ENSG00000187634             53             125             534
    # ENSG00000188976             53             125             534
    # ENSG00000187961             53             125             534
    # ENSG00000187583             53             125             534
    ###########################################################################
    # 4. Sort sequences of sample matrix and gene matrix to match.

    # New function
    # parameters...
    # major tissue
    # minor tissues
    # data_samples_tissues_patients
    # data_gene_count
    # returns...
    # sample matrix and gene matrix in proper format and sort sequence for DESeq2

    data_sets = organize_differential_expression_data_sets(
        data_samples_tissues_patients=source["data_samples_tissues_patients"],
        data_gene_count=source["data_gene_count"],
    )








    ########################################################################


    ################################################

    # Organization of samples by patients and tissues.
    # Collect unique patients, unique tissues for each patient, and unique
    # samples for each tissue of each patient.
    utility.print_terminal_partition(level=2)
    print("Organization of samples by hierarchy of patients and tissues.")
    print("Collection of hierarchical groups by patient, tissue, and sample.")
    utility.print_terminal_partition(level=4)
    print(
        "data hierarchy is " +
        "patients (714) -> tissues (30) -> samples (11688) -> genes (~50,000)"
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
        "tissue (30) -> patients (714) -> samples (11688) -> genes (~50,000)"
    )
    tissues_patients_samples = collect_tissues_patients_samples(
        data_samples_tissues_patients=data_samples_tissues_patients
    )
    print("Printing first 10 unique tissues...")
    print(list(tissues_patients_samples.keys())[:9])
    print("Printing groups for tissue 'Bladder'... ")
    print(tissues_patients_samples["Bladder"])






    # Compile information.
    information = {
        "patients_tissues_samples": patients_tissues_samples,
        "tissues_patients_samples": tissues_patients_samples,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
