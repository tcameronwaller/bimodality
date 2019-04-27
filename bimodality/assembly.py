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
import gtfparse

# Custom

import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Source.


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
    path_gene_annotation = os.path.join(
        path_access, "annotation_gene_gencode.gtf"
    )
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")
    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_patient_attribute = pandas.read_csv(
        path_attribute_patient,
        sep="\t",
        header=0,
        #nrows=1000,
    )
    data_sample_attribute = pandas.read_csv(
        path_attribute_sample,
        sep="\t",
        header=0,
        #nrows=1000,
    )
    data_gene_annotation = gtfparse.read_gtf(path_gene_annotation)
    data_gene_signal = pandas.read_csv(
        path_gene_signal,
        sep="\t",
        header=2,
        #nrows=1000,
    )
    # Compile and return information.
    return {
        "data_patient_attribute": data_patient_attribute,
        "data_sample_attribute": data_sample_attribute,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_signal": data_gene_signal,
    }


##########
# Summary.


def summarize_raw_data(
    data_patient_attribute=None,
    data_sample_attribute=None,
    data_gene_annotation=None,
    data_gene_signal=None
):
    """
    Optimizes data types.

    arguments:
        data_patient_attribute (object): Pandas data frame of attributes for all
            samples.
        data_sample_attribute (object): Pandas data frame of attributes for all
            samples.
        data_gene_annotation (object): Pandas data frame of annotation of
            genes.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples.

    raises:

    returns:

    """

    # Samples beginning with code patient "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of patients' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of patients' attributes.")
    print(data_patient_attribute)
    print(data_patient_attribute.iloc[0:10, 0:10])

    # Summarize table of samples' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of samples' attributes.")
    print(data_sample_attribute)
    print(data_sample_attribute.iloc[0:10, 0:7])

    # Summarize table of genes' annotations.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' annotations from GENCODE version 29.")
    print(data_gene_annotation)
    print(data_gene_annotation.iloc[0:10, 0:10])

    # Summarize table of genes' signals in all samples.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' signals in all samples.")
    print("Genes' signals are in Transcripts Per Million (TPM).")
    print(data_gene_signal)
    print(data_gene_signal.iloc[0:10, 0:10])
    #print(source["data_gene_signal"].loc[:,"GTEX-1117F-0226-SM-5GZZ7"])

    # How to reference information about a specific sample.
    if False:
        print(data_sample_attribute.loc[
            data_sample_attribute["SAMPID"] == "GTEX-1117F-0003-SM-58Q7G"
        ]
        )
    pass


##########
# Organization of samples' attributes, specifically associations to patients
# and tissues.


def extract_gtex_sample_patient_identifier(sample=None):
    """
    Extracts the patient's identifier from a sample's identifier.

    arguments:
        sample (str): identifier of a sample

    raises:

    returns:
        (str): identifier of a patient

    """

    split_strings = sample.split("-")
    patient = "-".join(split_strings[0:2])
    return patient


def translate_sex(value=None):
    """
    Translates annotations of sex.

    arguments:
        value (int): Annotation of sex

    raises:

    returns:
        (str): Name of sex

    """

    if value == 1:
        sex = "male"
    elif value == 2:
        sex = "female"
    return sex


def collect_samples_tissues_patients(
    samples=None,
    data_patient_attribute=None,
    data_sample_attribute=None
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
        samples (list<str>): identifiers of samples
        data_patient_attribute (object): Pandas data frame of attributes for
            all patients
        data_sample_attribute (object): Pandas data frame of attributes for all
            samples

    raises:

    returns:
        (list<dict<str>>): information about patients and tissues for samples

    """

    # Organize data.
    data_sample_attribute.set_index(
        ["SAMPID"],
        append=False,
        drop=True,
        inplace=True
    )
    data_patient_attribute.set_index(
        ["SUBJID"],
        append=False,
        drop=True,
        inplace=True
    )
    # Collect tissues and patients for each sample.
    samples_tissues_patients = list()
    for sample in samples:
        #tissue = data_sample_attribute.loc[
        #    data_sample_attribute["SAMPID"] == identifier_sample,
        #].at[0, "SMTS"]
        # Access tissue attributes.
        tissue_major = data_sample_attribute.at[sample, "SMTS"]
        tissue_minor = data_sample_attribute.at[sample, "SMTSD"]
        # Access patient attributes.
        patient = extract_gtex_sample_patient_identifier(sample=sample)
        sex_raw = data_patient_attribute.at[patient, "SEX"]
        sex = translate_sex(value=sex_raw)
        age = data_patient_attribute.at[patient, "AGE"]
        record = {
            "sample": sample,
            "tissue_major": tissue_major,
            "tissue_minor": tissue_minor,
            "patient": patient,
            "sex": sex,
            "age": age,
        }
        samples_tissues_patients.append(record)
    # Return information.
    return samples_tissues_patients


def organize_samples_tissues_patients(
    data_patient_attribute=None,
    data_sample_attribute=None,
    data_gene_signal=None
):
    """
    Optimizes data types.

    arguments:
        data_patient_attribute (object): Pandas data frame of attributes for
            all patients
        data_sample_attribute (object): Pandas data frame of attributes for all
            samples
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples

    raises:

    returns:
        (dict): information about samples, tissues, and patients

    """

    # Extract association of samples to patients and tissues.
    utility.print_terminal_partition(level=1)
    print("Association of samples to patients and tissues.")

    # Organize samples by hierarchy of patients and tissues.
    # Collect tissues and patients for each sample.
    utility.print_terminal_partition(level=2)
    print("Collection of tissues and patients for each sample.")
    print("Extract patient identifiers from sample identifiers.")
    # Extract identifiers of samples with measurements for genes.
    # Extract names of columns.
    # Pandas series
    #headers = data_gene_signal.columns
    # NumPy array
    #headers = data_gene_signal.columns.values
    # List
    headers = data_gene_signal.columns.to_list()
    # Exclude name and description.
    samples = headers[2:]
    # Collect information about samples.
    samples_tissues_patients = collect_samples_tissues_patients(
        samples=samples,
        data_patient_attribute=data_patient_attribute,
        data_sample_attribute=data_sample_attribute
    )
    data_samples_tissues_patients = utility.convert_records_to_dataframe(
        records=samples_tissues_patients
    )
    data_samples_tissues_patients.set_index(
        ["sample"],
        append=False,
        drop=True,
        inplace=True
    )
    data_samples_tissues_patients.rename_axis(
        index="sample",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_samples_tissues_patients.rename_axis(
        columns="properties",
        axis="columns",
        copy=False,
        inplace=True
    )
    print(data_samples_tissues_patients.iloc[0:10, :])
    print(data_samples_tissues_patients.shape)

    # Return information.
    return data_samples_tissues_patients


##########
# Organization of genes' annotations.


def extract_gene_identifier(string):
    """
    Removes designation of version from Ensembl gene identifier.

    arguments:
        string (str): Ensembl identifier of gene.

    raises:

    returns:
        (str): Ensemble identifier of gene without designation of version.

    """

    string_split = string.split(".")
    identifier_gene = string_split[0]
    return identifier_gene


def define_redundant_records():
    """
    Defines genes for which records are redundant.

    arguments:

    raises:

    returns:
        (dict<str>): genes with redundant records

    """

    correction = dict()
    correction["ENSG00000002586"] = "CD99"
    correction["ENSG00000196433"] = "ASMT"
    correction["ENSG00000197976"] = "AKAP17A"
    correction["ENSG00000182162"] = "P2RY8"
    correction["ENSG00000167393"] = "PPP2R3B"
    correction["ENSG00000185291"] = "IL3RA"
    correction["ENSG00000205755"] = "CRLF2"
    correction["ENSG00000182378"] = "PLCXD1"
    correction["ENSG00000198223"] = "CSF2RA"
    correction["ENSG00000214717"] = "ZBED1"
    correction["ENSG00000185960"] = "SHOX"
    correction["ENSG00000169084"] = "DHRSX"
    correction["ENSG00000168939"] = "SPRY3"
    correction["ENSG00000169093"] = "ASMTL"
    correction["ENSG00000178605"] = "GTPBP6"
    correction["ENSG00000124333"] = "VAMP7"
    correction["ENSG00000169100"] = "SLC25A6"
    return correction


def organize_genes_annotations(
    data=None
):
    """
    Organizes genes' annotations.

    arguments:
        data (object): Pandas data frame of genes' annotations from GENCODE.

    raises:

    returns:
        (object): Pandas data frame of genes' annotations.

    """

    # Organize annotations of genes.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' annotations.")

    # Define and select relevant columns.
    print(data.shape)
    columns = ["feature", "gene_id", "gene_type", "gene_name"]
    data_interest = data.loc[ :, columns]
    print(data_interest.shape)
    print(data_interest.iloc[0:10, 0:15])
    # Select entries for genes.
    data_gene_annotation = (
        data_interest.loc[data_interest["feature"] == "gene"]
    )
    data_gene_annotation.drop(
        labels="feature",
        axis="columns",
        inplace=True
    )
    # Remove version designations from genes' identifiers.
    data_gene_annotation["identifier"] = (
        data_gene_annotation["gene_id"].apply(extract_gene_identifier)
    )
    data_gene_annotation.drop(
        labels="gene_id",
        axis="columns",
        inplace=True
    )
    # Remove redundant records.
    # Some genes' records are redundant.
    redundancy = define_redundant_records()
    data_gene_annotation.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    # Organize axes.
    data_gene_annotation.set_index(
        ["identifier"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data_gene_annotation.iloc[0:10, 0:15])

    utility.print_terminal_partition(level=1)

    # Some genes' records are redundant.
    print("check for redundancy in genes' annotation records...")
    genes_redundant = list(redundancy.keys())
    for gene in genes_redundant:
        print(data_gene_annotation.loc[gene])
        print(data_gene_annotation.loc[gene, "gene_name"])
    # Return information.
    return data_gene_annotation


##########
# Organization of genes' signals.


def organize_data_axes_indices(data=None):
    """
    Organizes data with names and indices.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples.

    """

    # The Pandas function to rename axis copies data deeply by default and
    # requires a lot of memory.
    utility.print_terminal_partition(level=2)

    # Remove version designations from genes' identifiers.
    print("Remove version designations from genes' identifiers.")
    data["gene"] = data["Name"].apply(extract_gene_identifier)
    data.drop(
        labels="Name",
        axis="columns",
        inplace=True
    )
    # Drop column for genes' descriptions.
    print("Drop column for genes' description.")
    data.drop(
        labels="Description",
        axis="columns",
        inplace=True,
    )
    print(data.iloc[0:10, 0:10])
    print("Organize data with names and indices.")
    data.set_index(
        "gene",
        drop=True,
        inplace=True,
    )
    data.rename_axis(
        index="gene",
        axis="index",
        copy=False,
        inplace=True,
    )
    data.rename_axis(
        columns="sample",
        axis="columns",
        copy=False,
        inplace=True
    )
    print(data.iloc[0:10, 0:10])
    return data


def optimize_data_types(data=None):
    """
    Optimizes data types.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples.

    """

    utility.print_terminal_partition(level=2)
    print("Optimize efficiency of storage of data.")

    # Summarize data type and size.
    print("Print technical information about data.")
    print(data.info())

    # All values are of the same type.

    if False:
        # data.select_dtypes(include=["float"])
        data_type = data.apply(pandas.to_numeric, downcast="float")

    # Optimize data types.
    # Store genes' signals as numeric values of type float in 32 bits.
    data_type = data.astype("float32")

    # Summarize data type and size.
    print("Print technical information about data.")
    print(data_type.info())

    return data_type


def check_missing_values(data=None):
    """
    Checks data for missing values and prints reports.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print("Check for missing values in genes' signals.")
    print("shape of original data frame: " + str(data.shape))
    print("shape without missing axis 0: " + str(data.dropna(axis=0).shape))
    print("shape without missing axis 1: " + str(data.dropna(axis=1).shape))
    pass


def check_redundancy_genes(data=None):
    """
    Checks data for redundancy in genes.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print("Check for redundant genes in genes' signals.")
    print("Consider names of genes.")
    # Reset indices to consider names of genes.
    data = data.reset_index()
    print(data.iloc[0:10, 0:10])
    data_redundancy = data.duplicated(subset=None, keep="first")
    data_redundancy_list = data_redundancy.to_list()
    if any(data_redundancy_list):
        print("Redundancy in genes: Yes")
    else:
        print("Redundancy in genes: No")
    pass


def check_zero_samples(data=None):
    """
    Checks data for samples with values of 0 for all genes' signals.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print("Check for samples with values of 0 for all genes' signals.")
    print("shape of original data frame: " + str(data.shape))
    data_nonzero = (data != 0)
    print(
        "shape of data frame without zero samples: " +
        str(data.loc[ : , data_nonzero.any(axis="index")].shape)
    )
    pass


def check_zero_genes(data=None):
    """
    Checks data for genes with values of 0 for signals across all samples.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print("Check for genes with values of 0 for signals across all samples.")
    print("These genes are undetectable.")
    print("shape of original data frame: " + str(data.shape))
    data_nonzero = (data != 0)
    print(
        "shape of data frame without zero samples: " +
        str(data.loc[data_nonzero.any(axis="columns"), : ].shape)
    )
    print("Now printing a summary of data for genes with all zero signals.")
    data_zero = (data == 0)
    data_signal_zero = data.loc[data_zero.all(axis="columns"), : ]
    print(data_signal_zero.iloc[0:10, 0:10])
    #groups = data_signal_zero.groupby(level="gene")
    #print(groups.describe())
    pass


def drop_undetectable_genes(data=None):
    """
    Drops genes with values of 0 for signals across all samples.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples.

    """

    utility.print_terminal_partition(level=2)
    print("Drop genes that are undetectable.")
    data_nonzero = (data != 0)
    data_signal = data.loc[data_nonzero.any(axis="columns"), : ]
    print("Data without undetectable genes.")
    print(data_signal.iloc[0:10, 0:10])
    print("data dimensions: " + str(data_signal.shape))
    return data_signal


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
    reference=None,
    data_gene_signal=None
):
    """
    Collects tissues, patients, and genes' signals for each sample.

    arguments:
        reference (dict<dict<str>>): Tissue and patient for each sample.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

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


def organize_genes_signals(
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

    # Organize genes' signals.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' signals.")

    # Organize data with names and indices.
    data_gene_signal = organize_data_axes_indices(
        data=data_gene_signal
    )

    # Optimize data types.
    data_gene_signal = optimize_data_types(data=data_gene_signal)

    # Check for data quality.
    utility.print_terminal_partition(level=2)
    print("Check for quality of genes' signals.")

    # Check for missing values of genes' signals.
    check_missing_values(data=data_gene_signal)

    # Check for redundant genes.
    check_redundancy_genes(data=data_gene_signal)

    # Check for samples with values of 0 for all genes' signals.
    check_zero_samples(data=data_gene_signal)

    # Check for genes with values of 0 for signals across all samples.
    check_zero_genes(data=data_gene_signal)

    # Remove irrelevant signals for genes.
    utility.print_terminal_partition(level=2)
    print("Removal of signals for undetectable genes.")
    # Drop undetectable genes.
    data_gene_signal = drop_undetectable_genes(data=data_gene_signal)

    print("Original count of genes was 56202.")
    print("Count of genes with nonzero signal is 55863.")

    # Associate genes' signals to patients and tissues.
    utility.print_terminal_partition(level=2)
    print(
        "Association of samples to patients and tissues along with genes' " +
        "signals."
    )

    # Create reference for tissues and patients of each sample.
    reference = collect_samples_tissues_patients_reference(
        data_samples_tissues_patients=data_samples_tissues_patients
    )

    # Include identifiers for patient and tissue for each sample in data for
    # genes' signals.
    data_gene_signal = collect_genes_patients_tissues_samples(
        reference=reference,
        data_gene_signal=data_gene_signal
    )
    #print(data_gene_signal)
    print(data_gene_signal.iloc[0:10, 0:7])

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
    path_assembly = os.path.join(dock, "assembly")
    utility.confirm_path_directory(path_assembly)
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

    # Write information to file.
    pandas.to_pickle(
        information["data_samples_tissues_patients"],
        path_samples_tissues_patients
    )
    with open(path_patients_tissues_samples, "wb") as file_product:
        pickle.dump(information["patients_tissues_samples"], file_product)
    pandas.to_pickle(
        information["tissues_patients_samples"],
        path_tissues_patients_samples
    )
    pandas.to_pickle(
        information["data_gene_annotation"],
        path_gene_annotation
    )
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

    # Data for genes' signals is extensive.
    # Conserve memory.
    # Avoid unnecessary copies of the data.
    # Containerize portions of script within separate functions.
    # Collect garbage frequently.
    # Optimize data types of genes' signals.

    # Enable automatic garbage collection to clear memory.
    gc.enable()

    # Read source information from file.
    source = read_source(dock=dock)

    ##################################################
    ##################################################
    ##################################################

    # Summarize the structures of the raw data.
    if False:
        summarize_raw_data(
            data_patient_attribute=source["data_patient_attribute"],
            data_sample_attribute=source["data_sample_attribute"],
            data_gene_annotation=source["data_gene_annotation"],
            data_gene_signal=source["data_gene_signal"]
        )

    ##################################################
    ##################################################
    ##################################################

    # Organize associations of samples to patients and tissues.
    data_samples_tissues_patients = organize_samples_tissues_patients(
        data_patient_attribute=source["data_patient_attribute"],
        data_sample_attribute=source["data_sample_attribute"],
        data_gene_signal=source["data_gene_signal"]
    )

    ##################################################
    ##################################################
    ##################################################

    # Organize genes' annotations.
    data_gene_annotation = organize_genes_annotations(
        data=source["data_gene_annotation"]
    )

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # TODO: Problem... at this point "organize_genes_signals" includes identifiers for patients and tissues for each sample... I'm not ready for that yet...
    # TODO: do that later... like before it's time to split the gene signals or something...

    # Organize genes' signals.
    data_gene_signal = organize_genes_signals(
        data_samples_tissues_patients=data_samples_tissues_patients,
        data_gene_signal=source["data_gene_signal"]
    )

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # Compile information.
    information = {
        "data_samples_tissues_patients": data_samples_tissues_patients,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_signal": data_gene_signal
    }

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
