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
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_annotation = pandas.read_pickle(
        path_gene_annotation
    )
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_signal": data_gene_signal,
    }


##########
# Scrap from previous Assembly procedure...

####################################################
# Move to Selection procedure... ?


def collect_persons_tissues_samples(data_samples_tissues_persons=None):
    """
    Collect hierarchical structure of persons, tissues, and samples.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<list<str>>): Samples for each tissue of each person.

    """

    samples_tissues_persons = utility.convert_dataframe_to_records(
        data=data_samples_tissues_persons
    )
    # Collect unique tissues and samples for each person.
    persons_tissues_samples = dict()
    for record in samples_tissues_persons:
        sample = record["sample"]
        tissue = record["tissue"]
        person = record["person"]
        # Determine whether an entry already exists for the person.
        if person in persons_tissues_samples:
            # Determine whether an entry already exists for the tissue.
            if tissue in persons_tissues_samples[person]:
                persons_tissues_samples[person][tissue].append(sample)
            else:
                persons_tissues_samples[person][tissue] = list([sample])
        else:
            persons_tissues_samples[person] = dict()
            persons_tissues_samples[person][tissue] = list([sample])
    return persons_tissues_samples


def collect_tissues_persons_samples(data_samples_tissues_persons=None):
    """
    Collect hierarchical structure of tissues, persons, and samples.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<list<str>>): Samples for each person of each tissue.

    """

    samples_tissues_persons = utility.convert_dataframe_to_records(
        data=data_samples_tissues_persons
    )
    # Collect unique persons and samples for each tissue.
    tissues_persons_samples = dict()
    for record in samples_tissues_persons:
        sample = record["sample"]
        tissue = record["tissue"]
        person = record["person"]
        # Determine whether an entry already exists for the tissue.
        if tissue in tissues_persons_samples:
            # Determine whether an entry already exists for the person.
            if person in tissues_persons_samples[tissue]:
                tissues_persons_samples[tissue][person].append(sample)
            else:
                tissues_persons_samples[tissue][person] = list([sample])
        else:
            tissues_persons_samples[tissue] = dict()
            tissues_persons_samples[tissue][person] = list([sample])
    return tissues_persons_samples


def expand_print_persons_tissues_samples(persons_tissues_samples=None):
    """
    Collects tissues and samples for each person.

    arguments:
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each persons.

    raises:

    returns:


    """

    print(list(persons_tissues_samples.keys()))
    for person in persons_tissues_samples:
        print("person: " + person)
        for tissue in persons_tissues_samples[person]:
            print("tissue: " + tissue)
            for sample in persons_tissues_samples[person][tissue]:
                print("sample: " + sample)
    pass

########################################################

##########


##########
# Tissues and persons.

def select_tissues_persons(
    data_samples_tissues_persons=None,
    persons_tissues_samples=None,
    tissues_persons_samples=None
):
    """
    Selects tissues and persons of interest for further analyses.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        tissues_persons_samples (dict<dict<list<str>>): Samples for each
            person of each tissue.

    raises:

    returns:
        (dict<list<str>>): Information about tissues and persons.

    """

    utility.print_terminal_partition(level=1)
    print("Selection of tissues and persons of interest.")

    #

    # Count persons with samples for each tissue.
    # How many persons have samples for each tissue?
    ##################################################
    # tissue persons
    # brain   50
    # liver   10
    # heart   5
    ##################################################
    utility.print_terminal_partition(level=2)
    print("Counts of persons with samples for each tissue.")
    data_tissues_persons_counts = count_tissues_persons(
        filter=False,
        tissues=[],
        persons=[],
        tissues_persons_samples=tissues_persons_samples
    )
    print(data_tissues_persons_counts)

    # Count tissues for which each person has samples.
    ##################################################
    # person tissues
    # alice   3
    # john    5
    # paul    10
    # mary    13
    ##################################################
    utility.print_terminal_partition(level=2)
    print("Counts of tissues for which each person has samples.")
    data_persons_tissues_counts = count_persons_tissues(
        filter=False,
        tissues=[],
        persons_tissues_samples=persons_tissues_samples
    )
    print(data_persons_tissues_counts)

    # How many persons have samples for each count of tissues?
    # Populate bins of persons with each count of tissues.
    # Plot histogram.
    ##################################################
    # tissues persons
    # 3       50
    # 5       10
    # 7       5
    ##################################################
    utility.print_terminal_partition(level=2)
    print("Counts of persons with samples for each count of tissues.")
    data_persons_tissues_counts_bins = count_bins_persons_tissues(
        data_persons_tissues_counts=data_persons_tissues_counts
    )
    print(data_persons_tissues_counts_bins)

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

    # Collect persons with samples for all of the selection tissues.
    utility.print_terminal_partition(level=2)
    print(
        "Counts of selection tissues for which each person has samples."
    )
    data_persons_tissues_coverage = count_persons_tissues(
        filter=True,
        tissues=tissues,
        persons_tissues_samples=persons_tissues_samples
    )
    print(data_persons_tissues_coverage)
    #print(data_gene_signal_mean.loc["GTEX-12WSD", ])
    # Filter for those persons which miss fewer than 4 selection tissues.
    print("Allow each person to miss fewer than 4 selection tissues.")
    print(
        "Hence the extent of imputation will be less than 40% of tissues " +
        "for each person."
    )
    misses = 3
    matches = tissues_count - misses
    summary_persons_coverage = (
        data_persons_tissues_coverage.loc[
            data_persons_tissues_coverage["count"] >= matches
        ]
    )
    print(summary_persons_coverage)
    print(
        "count of selection persons: " +
        str(summary_persons_coverage.shape[0])
    )

    # Filter for persons with adequate coverage of tissues of interest.
    persons = filter_persons_tissues_coverage(
        threshold=matches,
        data_persons_tissues_counts=data_persons_tissues_coverage
    )
    print(persons)
    print("count of selection persons: " + str(len(persons)))
    print("Call these persons the selection persons.")

    # Collect selection persons that have samples for each tissue.
    # Summarize the extent of imputation necessary for each tissue.
    ##################################################
    # tissue  persons
    # brain   125
    # liver   175
    # kidney  150
    # lung    255
    # skin    300
    ##################################################
    utility.print_terminal_partition(level=2)
    print(
        "Collection of selection persons with samples for each selection " +
        "tissue."
    )
    print(
        "Summary of extent of imputation necessary for each tissue."
    )
    data_tissues_persons_coverage = count_tissues_persons(
        filter=True,
        tissues=tissues,
        persons=persons,
        tissues_persons_samples=tissues_persons_samples
    )
    print(data_tissues_persons_coverage)

    # Compile information.
    information = {
        "tissues": tissues,
        "persons": persons
    }

    # Return information.
    return information


def filter_samples_tissues_persons(
    tissues=None,
    persons=None,
    data_samples_tissues_persons=None
):
    """
    Filters samples by tissues and persons.

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons with signal for tissues of interest.
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples.

    """

    samples_tissues_persons = utility.convert_dataframe_to_records(
        data=data_samples_tissues_persons
    )
    # Determine whether each sample passes filters by tissues and persons.
    samples_tissues_persons_filter = list()
    for record in samples_tissues_persons:
        sample = record["sample"]
        tissue = record["tissue"]
        person = record["person"]
        if (person in persons) and (tissue in tissues):
            samples_tissues_persons_filter.append(record)
    return utility.convert_records_to_dataframe(
        records=samples_tissues_persons_filter
    )


def count_tissues_persons(
    filter=None, tissues=None, persons=None, tissues_persons_samples=None
):
    """
    Counts persons with samples for each tissue.

    arguments:
        filter (bool): Whether to filter the tissues and persons for counts.
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons of interest.
        tissues_persons_samples (dict<dict<list<str>>): Samples for each
            person of each tissue.

    raises:

    returns:
        (object): Pandas data frame of counts of persons with samples for each
            tissue.

    """

    # Count unique persons with samples for each tissue.
    coverage = list()
    # Determine whether to filter tissues.
    if filter:
        tissues_relevant = tissues
    else:
        tissues_relevant = list(tissues_persons_samples.keys())
    # Iterate on tissues.
    for tissue in tissues_relevant:
        # Collect unique persons with samples for tissue.
        persons_tissue = list(tissues_persons_samples[tissue].keys())
        # Determine whether to filter persons.
        if filter:
            persons_relevant = utility.filter_common_elements(
                list_one=persons, list_two=persons_tissue
            )
        else:
            persons_relevant = persons_tissue
        # Count persons.
        persons_count = len(persons_relevant)
        # Compile record.
        record = dict()
        record["tissue"] = tissue
        record["persons"] = persons_relevant
        record["count"] = persons_count
        # Include record.
        coverage.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(records=coverage)
    data_sort = data.sort_values("count", axis=0, ascending=False)
    return data_sort


def count_persons_tissues(
    filter=None, tissues=None, persons_tissues_samples=None
):
    """
    Counts tissues for which each person has samples.

    arguments:
        filter (bool): Whether to filter the tissues and persons for counts.
        tissues (list<str>): Tissues of interest.
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.

    raises:

    returns:
        (object): Pandas data frame of counts of tissues for which each person
            has samples.

    """

    # Count unique tissues for which each person has samples.
    coverage = list()
    # Iterate on persons.
    for person in persons_tissues_samples:
        # Collect unique tissues for which person has samples.
        tissues_person = list(persons_tissues_samples[person].keys())
        # Determine whether to filter tissues.
        if filter:
            tissues_relevant = utility.filter_common_elements(
                list_one=tissues, list_two=tissues_person
            )
        else:
            tissues_relevant = tissues_person
        # Count tissues.
        tissues_count = len(tissues_relevant)
        # Compile record.
        record = dict()
        record["person"] = person
        record["tissues"] = tissues_relevant
        record["count"] = tissues_count
        # Include record.
        coverage.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(records=coverage)
    data_sort = data.sort_values("count", axis=0, ascending=False)
    return data_sort


def count_bins_persons_tissues(data_persons_tissues_counts=None):
    """
    Count persons in each bin by counts of tissues.

    arguments:
        data_persons_tissues_counts (object): Pandas data frame of counts of
            tissues for which each person has samples.

    raises:

    returns:
        (object): Pandas data frame of counts of persons in each bin by counts
            of tissues.

    """

    persons_tissues_counts = utility.convert_dataframe_to_records(
        data=data_persons_tissues_counts
    )
    # Compile summary.
    bins = dict()
    # Iterate on persons.
    for record in persons_tissues_counts:
        person = record["person"]
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

    # tissue combination: counts of persons which miss fewer than 3 tissues
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


def filter_persons_tissues_coverage(
    threshold=None, data_persons_tissues_counts=None
):
    """
    Filters persons by coverage of selection tissues.

    Includes counts of tissues of interest for which each person misses
    signal.

    arguments:
        threshold (int): Count of selection tissues which a person must have.
        data_persons_tissues_counts (object): Pandas data frame of counts of
            tissues for which each person has samples.

    raises:

    returns:
        (list<str>): persons with signal for selection tissues.

    """

    persons_tissues_counts = utility.convert_dataframe_to_records(
        data=data_persons_tissues_counts
    )
    persons = list()
    for record in persons_tissues_counts:
        person = record["person"]
        count = record["count"]
        # Determine whether the person has few enough misses.
        if count >= threshold:
            persons.append(person)
    return persons


##########
# Genes. ... I think this is scrap...


def select_samples(
    tissues=None,
    persons=None,
    data_gene_signal=None
):
    """
    Selects samples of interest for further analyses.

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons of interest.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons.

    raises:

    returns:
        (dict): Pandas data frame of genes' signals for all samples, tissues,
            and persons.

    """

    # Select samples from persons and tissues of interest.
    utility.print_terminal_partition(level=2)
    print("Selection of samples from persons and tissues of interest.")
    print("count of samples, original: " + str(data_gene_signal.shape[0]))
    data_gene_signal.reset_index(
        level=["person", "tissue", "sample"], inplace=True
    )
    data_gene_signal.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    data_gene_signal = data_gene_signal.loc[persons, : ]
    print(
        "count of samples from persons of interest: " +
        str(data_gene_signal.shape[0])
    )
    data_gene_signal.reset_index(level=["person"], inplace=True)
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
        ["person", "tissue", "sample"],
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
            samples, tissues, and persons.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons.

    """

    utility.print_terminal_partition(level=2)
    print(
        "Selection of detectable genes with nonzero signals in persons " +
        "and tissues of interest."
    )
    print("genes, original: " + str(data_gene_signal.shape[1]))
    data_nonzero = (data_gene_signal != 0)
    data_signal = data_gene_signal.loc[ : , data_nonzero.any(axis="index")]
    print("genes, detection: " + str(data_signal.shape[1]))
    return data_signal


def select_samples_genes(
    persons=None,
    tissues=None,
    data_gene_annotation=None,
    data_gene_signal=None
):
    """
    Selects samples and genes of interest for further analyses.

    arguments:
        persons (list<str>): persons of interest.
        tissues (list<str>): Tissues of interest.
        data_gene_annotation (object): Pandas data frame of genes' annotations.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons.

    """

    utility.print_terminal_partition(level=1)
    print("Selection of samples and genes of interest.")

    # Select samples from persons and tissues of interest.
    data_gene_signal = select_samples(
        persons=persons,
        tissues=tissues,
        data_gene_signal=data_gene_signal,
    )

    # Select genes with detectable, non-zero signal in tissues and persons of
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
    path_persons = os.path.join(path_selection, "persons.pickle")
    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    # Write information to file.
    with open(path_tissues, "wb") as file_product:
        pickle.dump(information["tissues"], file_product)
    with open(path_persons, "wb") as file_product:
        pickle.dump(information["persons"], file_product)
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

    ########################################################################
    # Move this stuff to the "selection" procedure?
    ################################################
    # Organization of samples by persons and tissues.
    # Collect unique persons, unique tissues for each person, and unique
    # samples for each tissue of each person.
    utility.print_terminal_partition(level=2)
    print("Organization of samples by hierarchy of persons and tissues.")
    print("Collection of hierarchical groups by person, tissue, and sample.")
    utility.print_terminal_partition(level=4)
    print(
        "data hierarchy is " +
        "persons (714) -> tissues (30) -> samples (11688) -> genes (~50,000)"
    )
    persons_tissues_samples = collect_persons_tissues_samples(
        data_samples_tissues_persons=data_samples_tissues_persons
    )
    if False:
        expand_print_persons_tissues_samples(
            persons_tissues_samples=persons_tissues_samples
        )
    print("Printing first 10 unique persons...")
    print(list(persons_tissues_samples.keys())[:9])
    print("Printing groups for person 'GTEX-14LZ3'... ")
    print(persons_tissues_samples["GTEX-14LZ3"])

    # Organization of samples by tissues and persons.
    # Collect unique tissues, unique persons for each tissue, and unique
    # samples for each person of each tissue.
    utility.print_terminal_partition(level=2)
    print("Organization of samples by hierarchy of tissues and persons.")
    print("Collection of hierarchical groups by tissue, person, and sample.")
    utility.print_terminal_partition(level=4)
    print(
        "data hierarchy is " +
        "tissue (30) -> persons (714) -> samples (11688) -> genes (~50,000)"
    )
    tissues_persons_samples = collect_tissues_persons_samples(
        data_samples_tissues_persons=data_samples_tissues_persons
    )
    print("Printing first 10 unique tissues...")
    print(list(tissues_persons_samples.keys())[:9])
    print("Printing groups for tissue 'Bladder'... ")
    print(tissues_persons_samples["Bladder"])
    ########################################################################


    ##################################################
    ##################################################
    ##################################################

    # Selection of tissues and persons of interest.
    tissues_persons = select_tissues_persons(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        persons_tissues_samples=source["persons_tissues_samples"],
        tissues_persons_samples=source["tissues_persons_samples"]
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
        tissues=tissues_persons["tissues"],
        persons=tissues_persons["persons"],
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
        "tissues": tissues_persons["tissues"],
        "persons": tissues_persons["persons"],
        "data_gene_signal": data_gene_signal,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)


    # Scrap.

    if False:
        # Summarize tissues and persons.
        # Extract unique values of tissue.
        tissues_test = list(source["data_gene_signal_imputation"]["tissue"].unique())
        persons_test = list(source["data_gene_signal_imputation"]["person"].unique())


    pass


if (__name__ == "__main__"):
    execute_procedure()
