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


def define_tissue_comparisons():
    """
    Defines minor categories to compare for each major category of tissue.

    arguments:

    raises:

    returns:
        (dict<list<str>>): Minor categories to compare for each major category
            of tissue

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


def check_index_column_sequence_match(
    data_index=None,
    data_column=None,
):
    """
    Check that counts and sequences of indexes and columns match.

    arguments:
        data_index (object): Pandas data frame with relevant values in index
        data_column (object): Pandas data frame with relevant values in columns

    raises:

    returns:
        (bool): Whether indices and columns match

    """

    # Extract indices and columns.
    indices = data_index.index.to_list()
    columns = data_column.columns.to_list()
    # Check counts.
    length = (len(indices) == len(columns))
    # Check sequences.
    matches = list()
    for index in range(len(indices)):
        value_index = indices[index]
        value_column = columns[index]
        matches.append(value_index == value_column)
    sequence = all(matches)
    return (length and sequence)


def organize_differential_expression_data_sets(
    comparisons=None,
    data_samples_tissues_patients=None,
    data_gene_count=None,
):
    """
    Collect hierarchical structure of tissues, patients, and samples.

    arguments:
        comparisons (dict<list<str>>): Minor categories to compare for each
            major category of tissue
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples
        data_gene_count (object): Pandas data frame of genes' counts for all
            samples

    raises:

    returns:
        (list<dict>): Collections of data sets for differential expression
            analyses

    """

    # Print terminal partition.
    utility.print_terminal_partition(level=2)
    # Report.
    print(
        "Organization of data sets for differential gene expression " +
        "comparison of minor categories of tissues."
    )

    # Collect data sets.
    sets = list()
    tissues_major = list(comparisons.keys())
    for tissue_major in tissues_major:
        set = organize_differential_expression_data_set(
            tissue_major=tissue_major,
            tissues_minor=comparisons[tissue_major],
            data_samples_tissues_patients=data_samples_tissues_patients,
            data_gene_count=data_gene_count,
        )
        # Collect the data set.
        sets.append(set)

    # Print terminal partition.
    utility.print_terminal_partition(level=2)
    # Report.
    print(
        "Data sets by major tissues:"
    )
    for set in sets:
        print(set["tissue"])
    return sets


def organize_differential_expression_data_set(
    tissue_major=None,
    tissues_minor=None,
    data_samples_tissues_patients=None,
    data_gene_count=None,
):
    """
    Collect hierarchical structure of tissues, patients, and samples.

    arguments:
        tissue_major (str): Name of a major category of tissue
        tissues_minor (list<str>): Names of minor categories of tissue
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples
        data_gene_count (object): Pandas data frame of genes' counts for all
            samples

    raises:

    returns:
        (dict): Information about samples and genes for a single
            comparison by differential expression

    """

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

    # Copy data.
    data_sample = data_samples_tissues_patients.copy(deep=True)
    data_gene = data_gene_count.copy(deep=True)
    # Filter samples by categories of tissues.
    # Filter samples by groups.
    data_sample = data_sample.loc[data_sample["tissue_major"] == tissue_major]
    data_sample = (
        data_sample.loc[data_sample["tissue_minor"].isin(tissues_minor), :]
    )
    #print(data_sample.iloc[0:10, 0:10])
    samples = data_sample.index.to_list()
    # Filter measurements of genes by samples.
    data_gene = data_gene.loc[:, data_gene.columns.isin(samples)]
    #print(data_gene.iloc[0:10, 0:10])

    # Sort sequence of samples.
    # Sort indices of data frames.
    #data = data.reindex(sorted(data.columns), axis="columns")
    data_sample.sort_index(axis="index", ascending=True, inplace=True)
    data_gene.sort_index(axis="columns", ascending=True, inplace=True)

    # Confirm that sequences of samples match.
    check = check_index_column_sequence_match(
        data_index=data_sample,
        data_column=data_gene,
    )
    print("Counts and sequences match: " + str(check))

    # Compile information.
    information = dict()
    information["tissue"] = tissue_major
    information["sample"] = data_sample
    information["gene"] = data_gene
    # Return information.
    return information


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

    write_product_sets(dock=dock, information=information)

    # Specify directories and files.
    path_tissue = os.path.join(dock, "tissue")
    utility.confirm_path_directory(path_tissue)
    path_tissues = os.path.join(
        path_tissue, "tissues.txt"
    )
    # Write information to file.
    utility.write_file_text_list(
        information=information["tissues"],
        path_file=path_tissues
    )
    pass


def write_product_sets(dock=None, information=None):
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
    path_tissue = os.path.join(dock, "tissue")
    utility.confirm_path_directory(path_tissue)
    path_sets = os.path.join(path_tissue, "sets")
    # Remove previous files since they change from run to run.
    utility.remove_directory(path=path_sets)
    utility.confirm_path_directory(path_sets)
    # Iterate on sets.
    for set in information["sets"]:
        # Access information.
        tissue = set["tissue"]
        data_sample = set["sample"]
        data_gene = set["gene"]
        # Specify directories and files.
        path_sample = os.path.join(
            path_sets, (tissue + "_samples.tsv")
        )
        path_gene = os.path.join(
            path_sets, (tissue + "_genes.tsv")
        )
        # Write information to file.
        data_sample.to_csv(
            path_or_buf=path_sample,
            sep="\t",
            header=True,
            index=True,
        )
        data_gene.to_csv(
            path_or_buf=path_gene,
            sep="\t",
            header=True,
            index=True,
        )
        pass
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

    # Organize data sets for differential gene expression comparisons of minor
    # categories of tissues.
    comparisons = define_tissue_comparisons()
    tissues_major = list(comparisons.keys())
    sets = organize_differential_expression_data_sets(
        comparisons=comparisons,
        data_samples_tissues_patients=source["data_samples_tissues_patients"],
        data_gene_count=source["data_gene_count"],
    )

    # Compile information.
    information = {
        "tissues": tissues_major,
        "sets": sets,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
