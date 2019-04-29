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
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_samples_tissues_patients = pandas.read_pickle(
        path_samples_tissues_patients
    )
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_samples_tissues_patients": data_samples_tissues_patients,
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
        ascending=False
    )
    return data_sort


def describe_samples_categories(
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

    # Print terminal partition.
    utility.print_terminal_partition(level=2)
    # Report.
    print("Description of categories of samples.")

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




####################################################
# Move to Selection procedure... ?


def define_tissue_translations():
    """
    Defines and translates names of tissues.

    arguments:

    raises:

    returns:
        (dict): definitions of tissue names

    """

    # Define names of tissues.
    tissues_major = dict()
    tissues_major["Adipose Tissue"] = "adipose"
    tissues_major["Adrenal Gland"] = "adrenal"
    tissues_major["Bladder"] = "bladder"
    tissues_major["Blood"] = "blood"
    tissues_major["Blood Vessel"] = "vessel"
    tissues_major["Bone Marrow"] = "marrow"
    tissues_major["Brain"] = "brain"
    tissues_major["Breast"] = "breast"
    tissues_major["Cervix Uteri"] = "cervix"
    tissues_major["Colon"] = "colon"
    tissues_major["Esophagus"] = "esophagus"
    tissues_major["Fallopian Tube"] = "fallopian"
    tissues_major["Heart"] = "Heart"
    tissues_minor = dict()
    tissues_minor["Adipose - Subcutaneous"] = "subcutane"
    # Compile information.
    information = dict()
    information["major"] = tissues_major
    information["minor"] = tissues_minor
    # Return information.
    return information


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

    # Explore samples by patient attributes.
    # Describe distributions of patients by sex and age.

    # Explore samples by tissue attributes.
    # Describe categories of all samples.
    describe_samples_categories(
        data_samples_tissues_patients=source["data_samples_tissues_patients"]
    )
    # Describe categories of specific samples.
    # Print terminal partition.
    utility.print_terminal_partition(level=1)
    # Report.
    print("Description of categories of samples for male patients.")
    data_samples_male = source["data_samples_tissues_patients"].loc[
        source["data_samples_tissues_patients"]["sex"] == "male", :
    ]
    #print(data_samples_male)
    describe_samples_categories(
        data_samples_tissues_patients=data_samples_male
    )
    # Print terminal partition.
    utility.print_terminal_partition(level=1)
    # Report.
    print("Description of categories of samples for female patients.")
    data_samples_female = source["data_samples_tissues_patients"].loc[
        source["data_samples_tissues_patients"]["sex"] == "female", :
    ]
    #print(data_samples_female)
    describe_samples_categories(
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
    print("Skin, Heart, Esophagus, Colon, Brain, Vessel, Adipose")


    # Organize data sets for comparison by differential expression of genes.





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
