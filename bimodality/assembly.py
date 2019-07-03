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
    path_attribute_person = os.path.join(path_access, "attribute_person.txt")
    path_gene_annotation = os.path.join(
        path_access, "annotation_gene_gencode.gtf"
    )
    path_gene_count = os.path.join(path_access, "count_gene.gct")
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")
    path_customization = os.path.join(dock, "customization")
    path_tissues_major = os.path.join(
        path_customization, "translation_tissues_major.tsv"
    )
    path_tissues_minor = os.path.join(
        path_customization, "translation_tissues_minor.tsv"
    )
    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_person_attribute = pandas.read_csv(
        path_attribute_person,
        sep="\t",
        header=0,
    )
    data_sample_attribute = pandas.read_csv(
        path_attribute_sample,
        sep="\t",
        header=0,
    )
    data_gene_annotation = gtfparse.read_gtf(path_gene_annotation)
    data_gene_count = pandas.read_csv(
        path_gene_count,
        sep="\t",
        header=2,
        #nrows=1000,
    )
    data_gene_signal = pandas.read_csv(
        path_gene_signal,
        sep="\t",
        header=2,
        #nrows=1000,
    )
    data_tissues_major = pandas.read_csv(
        path_tissues_major,
        sep="\t",
        header=0,
    )
    data_tissues_minor = pandas.read_csv(
        path_tissues_minor,
        sep="\t",
        header=0,
    )
    # Compile and return information.
    return {
        "data_person_attribute": data_person_attribute,
        "data_sample_attribute": data_sample_attribute,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_count": data_gene_count,
        "data_gene_signal": data_gene_signal,
        "data_tissues_major": data_tissues_major,
        "data_tissues_minor": data_tissues_minor,
    }


##########
# Organization of samples' attributes, specifically associations to persons
# and tissues.


def read_source_sample(dock=None):
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
    path_attribute_person = os.path.join(path_access, "attribute_person.txt")
    path_customization = os.path.join(dock, "customization")
    path_tissues_major = os.path.join(
        path_customization, "translation_tissues_major.tsv"
    )
    path_tissues_minor = os.path.join(
        path_customization, "translation_tissues_minor.tsv"
    )
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")
    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_person_attribute = pandas.read_csv(
        path_attribute_person,
        sep="\t",
        header=0,
    )
    data_sample_attribute = pandas.read_csv(
        path_attribute_sample,
        sep="\t",
        header=0,
    )
    data_tissues_major = pandas.read_csv(
        path_tissues_major,
        sep="\t",
        header=0,
    )
    data_tissues_minor = pandas.read_csv(
        path_tissues_minor,
        sep="\t",
        header=0,
    )
    data_gene_signal = pandas.read_csv(
        path_gene_signal,
        sep="\t",
        header=2,
        #nrows=1000,
    )
    # Compile and return information.
    return {
        "data_person_attribute": data_person_attribute,
        "data_sample_attribute": data_sample_attribute,
        "data_tissues_major": data_tissues_major,
        "data_tissues_minor": data_tissues_minor,
        "data_gene_signal": data_gene_signal,
    }


def summarize_raw_data_sample(
    data_person_attribute=None,
    data_sample_attribute=None,
):
    """
    Optimizes data types.

    arguments:
        data_person_attribute (object): Pandas data frame of attributes for all
            samples.
        data_sample_attribute (object): Pandas data frame of attributes for all
            samples.

    raises:

    returns:

    """

    # Samples beginning with code person "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of persons' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of persons' attributes.")
    print(data_person_attribute)
    print(data_person_attribute.iloc[0:10, 0:10])

    # Summarize table of samples' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of samples' attributes.")
    print(data_sample_attribute)
    print(data_sample_attribute.iloc[0:10, 0:7])

    pass


def extract_gtex_sample_person_identifier(sample=None):
    """
    Extracts the person's identifier from a sample's identifier.

    arguments:
        sample (str): identifier of a sample

    raises:

    returns:
        (str): identifier of a person

    """

    split_strings = sample.split("-")
    person = "-".join(split_strings[0:2])
    return person


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


def collect_samples_tissues_persons(
    samples=None,
    data_person_attribute=None,
    data_sample_attribute=None,
    data_tissues_major=None,
    data_tissues_minor=None,
):
    """
    Collects matches of samples, tissues, and persons.

    Product data format.
    Indices by person and tissue allow convenient access.
    ##################################################
    # sample   person   tissue
    # sm_1     bob       brain
    # sm_2     bob       heart
    # sm_3     bob       liver
    # sm_4     julie     brain
    # sm_5     julie     kidney
    # sm_6     betty     pancreas
    ##################################################

    arguments:
        samples (list<str>): identifiers of samples
        data_person_attribute (object): Pandas data frame of attributes for
            all persons
        data_sample_attribute (object): Pandas data frame of attributes for all
            samples
        data_tissues_major (object): Pandas data frame of translations for
            names of major tissues
        data_tissues_minor (object): Pandas data frame of translations for
            names of minor tissues

    raises:

    returns:
        (list<dict<str>>): information about persons and tissues for samples

    """

    # Organize data.
    data_sample_attribute.set_index(
        ["SAMPID"],
        append=False,
        drop=True,
        inplace=True
    )
    data_person_attribute.set_index(
        ["SUBJID"],
        append=False,
        drop=True,
        inplace=True
    )
    data_tissues_major.set_index(
        ["source"],
        append=False,
        drop=True,
        inplace=True
    )
    data_tissues_minor.set_index(
        ["source"],
        append=False,
        drop=True,
        inplace=True
    )
    # Collect tissues and persons for each sample.
    samples_tissues_persons = list()
    for sample in samples:
        #tissue = data_sample_attribute.loc[
        #    data_sample_attribute["SAMPID"] == identifier_sample,
        #].at[0, "SMTS"]
        # Access tissue attributes.
        major = data_sample_attribute.at[sample, "SMTS"]
        minor = data_sample_attribute.at[sample, "SMTSD"]
        tissue_major = data_tissues_major.at[major, "product"]
        tissue_minor = data_tissues_minor.at[minor, "product"]
        # Access person attributes.
        person = extract_gtex_sample_person_identifier(sample=sample)
        sex_raw = data_person_attribute.at[person, "SEX"]
        sex = translate_sex(value=sex_raw)
        age = data_person_attribute.at[person, "AGE"]
        record = {
            "sample": sample,
            "tissue_major": tissue_major,
            "tissue_minor": tissue_minor,
            "person": person,
            "sex": sex,
            "age": age,
        }
        samples_tissues_persons.append(record)
    # Return information.
    return samples_tissues_persons


def organize_samples_tissues_persons(
    dock=None,
):
    """
    Organize information about samples.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    # Read source information from file.
    source = read_source_sample(dock=dock)

    # Summarize the structures of the raw data.
    summarize_raw_data_sample(
        data_person_attribute=source["data_person_attribute"],
        data_sample_attribute=source["data_sample_attribute"],
    )

    # Extract association of samples to persons and tissues.
    utility.print_terminal_partition(level=1)
    print("Association of samples to persons and tissues.")

    # Organize samples by hierarchy of persons and tissues.
    # Collect tissues and persons for each sample.
    utility.print_terminal_partition(level=2)
    print("Collection of tissues and persons for each sample.")
    print(
        "Extract sample identifiers from the original data for genes' " +
        "signals."
    )
    print("These samples are comprehensive before any filters or selection.")
    print("Extract person identifiers from sample identifiers.")
    # Extract identifiers of samples with measurements for genes.
    # Extract names of columns.
    # Pandas series
    #headers = data_gene_signal.columns
    # NumPy array
    #headers = data_gene_signal.columns.values
    # List
    headers = source["data_gene_signal"].columns.to_list()
    # Exclude name and description.
    samples = headers[2:]
    # Collect information about samples.
    samples_tissues_persons = collect_samples_tissues_persons(
        samples=samples,
        data_person_attribute=source["data_person_attribute"],
        data_sample_attribute=source["data_sample_attribute"],
        data_tissues_major=source["data_tissues_major"],
        data_tissues_minor=source["data_tissues_minor"],
    )
    data_samples_tissues_persons = utility.convert_records_to_dataframe(
        records=samples_tissues_persons
    )
    data_samples_tissues_persons.set_index(
        ["sample"],
        append=False,
        drop=True,
        inplace=True
    )
    data_samples_tissues_persons.rename_axis(
        index="sample",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_samples_tissues_persons.rename_axis(
        columns="properties",
        axis="columns",
        copy=False,
        inplace=True
    )
    print(data_samples_tissues_persons.iloc[0:10, :])
    print(data_samples_tissues_persons.shape)

    # Return information.
    return data_samples_tissues_persons


##########
# Organization of genes' annotations.


def read_source_gene_annotation(dock=None):
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
    path_gene_annotation = os.path.join(
        path_access, "annotation_gene_gencode.gtf"
    )
    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_gene_annotation = gtfparse.read_gtf(path_gene_annotation)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
    }


def summarize_raw_data_gene_annotation(
    data_gene_annotation=None,
):
    """
    Optimizes data types.

    arguments:
        data_gene_annotation (object): Pandas data frame of annotation of
            genes.

    raises:

    returns:

    """

    # Samples beginning with code person "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of genes' annotations.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' annotations from GENCODE version 29.")
    print(data_gene_annotation)
    print(data_gene_annotation.iloc[0:10, 0:10])

    pass


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
    dock=None
):
    """
    Organizes genes' annotations.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of genes' annotations.

    """

    # Read source information from file.
    source = read_source_gene_annotation(dock=dock)

    # Summarize the structures of the raw data.
    summarize_raw_data_gene_annotation(
        data_gene_annotation=source["data_gene_annotation"],
    )

    # Organize annotations of genes.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' annotations.")

    # Define and select relevant columns.
    print(source["data_gene_annotation"].shape)
    columns = ["feature", "gene_id", "gene_type", "gene_name"]
    data_interest = source["data_gene_annotation"].loc[ :, columns]
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
    # Select entries for genes that encode proteins.
    print(
        "count of genes in GENCODE annotations: " +
        str(data_gene_annotation.shape[0])
    )
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["gene_type"] == "protein_coding"
        ]
    )
    print(data_gene_annotation.iloc[0:10, 0:7])
    print(
        "count of GENCODE genes of type 'protein_coding': " +
        str(data_gene_annotation.shape[0])
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
# Organization of genes' counts.


def read_source_gene_count(dock=None):
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
    path_gene_count = os.path.join(path_access, "count_gene.gct")
    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_gene_count = pandas.read_csv(
        path_gene_count,
        sep="\t",
        header=2,
        #nrows=1000,
    )
    # Compile and return information.
    return {
        "data_gene_count": data_gene_count,
    }


def summarize_raw_data_gene_count(
    data_gene_count=None,
):
    """
    Optimizes data types.

    arguments:
        data_gene_count (object): Pandas data frame of genes' counts across
            samples.

    raises:

    returns:

    """

    # Samples beginning with code person "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of genes' signals in all samples.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' counts in all samples.")
    print("Genes' counts represent reads mapped to genes.")
    print(data_gene_count)
    print(data_gene_count.iloc[0:10, 0:10])
    print("Count of genes: " + str(data_gene_count.shape[0]))
    print("Count of samples: " + str(data_gene_count.shape[1]))

    pass


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


def convert_data_types(data=None, type=None):
    """
    Optimizes data types.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples
        type (str): variable type to which to convert values in data

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
    data_type = data.astype(type)

    # Summarize data type and size.
    print("Print technical information about data.")
    print(data_type.info())

    return data_type


def organize_genes_counts(
    dock=None,
):
    """
    Organizes counts of genes.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of genes' counts across samples

    """

    # Read source information from file.
    source = read_source_gene_count(dock=dock)

    # Summarize the structures of the raw data.
    summarize_raw_data_gene_count(
        data_gene_count=source["data_gene_count"],
    )

    # Organize genes' signals.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' counts.")

    # Organize data with names and indices.
    data_gene_count = organize_data_axes_indices(
        data=source["data_gene_count"]
    )

    # Optimize data types.
    data_gene_count = convert_data_types(
        data=data_gene_count,
        type="int32"
    )

    # Return information.
    return data_gene_count



##########
# Organization of genes' signals.


def read_source_gene_signal(dock=None):
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
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")
    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_gene_signal = pandas.read_csv(
        path_gene_signal,
        sep="\t",
        header=2,
        #nrows=1000,
    )
    # Compile and return information.
    return {
        "data_gene_signal": data_gene_signal,
    }


def summarize_raw_data_gene_signal(
    data_gene_signal=None,
):
    """
    Optimizes data types.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples.

    raises:

    returns:

    """

    # Samples beginning with code person "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of genes' signals in all samples.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' signals in all samples.")
    print("Genes' signals are in Transcripts Per Million (TPM).")
    print(data_gene_signal)
    print(data_gene_signal.iloc[0:10, 0:10])
    print("Count of genes: " + str(data_gene_signal.shape[0]))
    print("Count of samples: " + str(data_gene_signal.shape[1]))
    #print(source["data_gene_signal"].loc[:,"GTEX-1117F-0226-SM-5GZZ7"])

    pass


def organize_genes_signals(
    dock=None,
):
    """
    Collects tissues, persons, and genes' signals for each sample.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Read source information from file.
    source = read_source_gene_signal(dock=dock)

    # Summarize the structures of the raw data.
    summarize_raw_data_gene_signal(
        data_gene_signal=source["data_gene_signal"],
    )

    # Organize genes' signals.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' signals.")

    # Organize data with names and indices.
    data_gene_signal = organize_data_axes_indices(
        data=source["data_gene_signal"]
    )

    # Optimize data types.
    data_gene_signal = convert_data_types(
        data=data_gene_signal,
        type="float32"
    )

    # Return information.
    return data_gene_signal


##########
# Association of samples, tissues, and persons


def collect_samples_persons_tissues_reference(
    data_samples_tissues_persons=None,
):
    """
    Collects person and tissue for each sample.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<str>>): person and tissue for each sample

    """

    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_samples.reset_index(level="sample", inplace=True)
    samples_tissues_persons = utility.convert_dataframe_to_records(
        data=data_samples
    )
    reference = dict()
    for record in samples_tissues_persons:
        sample = record["sample"]
        person = record["person"]
        tissue_major = record["tissue_major"]
        tissue_minor = record["tissue_minor"]
        if sample not in reference:
            reference[sample] = dict()
            reference[sample]["person"] = person
            reference[sample]["tissue_major"] = tissue_major
            reference[sample]["tissue_minor"] = tissue_minor
    return reference


def collect_genes_samples_persons_tissues(
    reference=None,
    data_gene_sample=None
):
    """
    Collects samples, persons, and tissues for each gene's signal.

    arguments:
        reference (dict<dict<str>>): Tissue and person for each sample.
        data_gene_sample (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons.

    """

    def match_person(sample):
        return reference[sample]["person"]
    def match_tissue_major(sample):
        return reference[sample]["tissue_major"]
    def match_tissue_minor(sample):
        return reference[sample]["tissue_minor"]
    # Designate person and tissue of each sample.
    data = data_gene_sample.copy(deep=True)
    data.reset_index(level="sample", inplace=True)
    data["person"] = (
        data["sample"].apply(match_person)
    )
    data["tissue_major"] = (
        data["sample"].apply(match_tissue_major)
    )
    data["tissue_minor"] = (
        data["sample"].apply(match_tissue_minor)
    )
    data.set_index(
        ["sample", "person", "tissue_major", "tissue_minor"],
        append=False,
        drop=True,
        inplace=True
    )
    return data


def associate_samples_persons_tissues(
    data_samples_tissues_persons=None,
    data_gene_sample=None
):
    """
    Associates samples, persons, and tissues for each gene's signal.

    The table data_gene_sample needs to have samples across rows and genes
    across columns.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_sample (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Associate genes' signals to persons and tissues.
    print("Associate samples to factors.")
    #print(data_samples_tissues_persons)
    #print(data_gene_sample)
    # Create reference for tissues and persons of each sample.
    reference = collect_samples_persons_tissues_reference(
        data_samples_tissues_persons=data_samples_tissues_persons
    )
    # Include identifiers for person and tissue for each sample in data for
    # genes' signals.
    data_gene_sample_person_tissue = collect_genes_samples_persons_tissues(
        reference=reference,
        data_gene_sample=data_gene_sample
    )
    #print(data_gene_signal)
    print(data_gene_sample_person_tissue.iloc[0:10, 0:7])
    # Return information.
    return data_gene_sample_person_tissue


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
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_samples_tissues_persons_text = os.path.join(
        path_assembly, "data_samples_tissues_persons.txt"
    )
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_gene_count = os.path.join(
        path_assembly, "data_gene_count.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.pickle"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_samples_tissues_persons"],
        path_samples_tissues_persons
    )
    information["data_samples_tissues_persons"].to_csv(
        path_or_buf=path_samples_tissues_persons_text,
        sep="\t",
        header=True,
        index=False,
    )
    pandas.to_pickle(
        information["data_gene_annotation"],
        path_gene_annotation
    )
    pandas.to_pickle(
        information["data_gene_count"],
        path_gene_count
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

    # Memory conservation
    # Data for genes' signals is extensive.
    # Conserve memory.
    # Avoid unnecessary copies of the data.
    # Containerize portions of script within separate functions.
    # Collect garbage frequently.
    # Optimize data types of genes' signals.
    # Organize independent portions of data separately to avoid storing data
    # within memory unnecessarily.

    # Enable automatic garbage collection to clear memory.
    gc.enable()

    ##################################################
    ##################################################
    ##################################################

    # Organize associations of samples to persons and tissues.
    data_samples_tissues_persons = organize_samples_tissues_persons(
        dock=dock,
    )

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # Organize genes' annotations.
    data_gene_annotation = organize_genes_annotations(
        dock=dock,
    )

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # Organize genes' counts.
    data_gene_count = organize_genes_counts(
        dock=dock,
    )

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # Organize genes' signals.
    data_gene_signal = organize_genes_signals(
        dock=dock,
    )

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # Compile information.
    information = {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_count": data_gene_count,
        "data_gene_signal": data_gene_signal
    }

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
