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
import itertools
import copy
import gc
import time

# Relevant

import numpy
import pandas
import scipy

# Custom

import assembly
import measurement
import prediction
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
##########
##########
# Genes


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
    path_assembly_annotation = os.path.join(dock, "assembly", "annotation")
    path_assembly_signal = os.path.join(dock, "assembly", "signal")

    path_gene_annotation_gtex = os.path.join(
        path_assembly_annotation, "data_gene_annotation_gtex.pickle"
    )
    path_gene_annotation_gencode = os.path.join(
        path_assembly_annotation, "data_gene_annotation_gencode.pickle"
    )
    path_genes_gtex = os.path.join(
        path_assembly_signal, "genes_gtex.pickle"
    )

    # Read information from file.
    data_gene_annotation_gtex = pandas.read_pickle(
        path_gene_annotation_gtex
    )
    data_gene_annotation_gencode = pandas.read_pickle(
        path_gene_annotation_gencode
    )
    with open(path_genes_gtex, "rb") as file_source:
        genes_gtex = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_gene_annotation_gtex": data_gene_annotation_gtex,
        "data_gene_annotation_gencode": data_gene_annotation_gencode,
        "genes_gtex": genes_gtex,
    }


def select_gene_annotation(
    genes_gtex=None,
    data_gene_annotation=None,
):
    """
    Selects annotations for genes that are on nuclear autosomes, encode
    proteins, and have signals in GTEx data.

    arguments:
        genes_gtex (list<str>): identifiers of genes with signals in GTEx data
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' annotations

    """

    utility.print_terminal_partition(level=2)
    print("Original shape of data: " + str(data_gene_annotation.shape))

    # Select entries for genes.
    data_gene_annotation = (
        data_gene_annotation.loc[data_gene_annotation["feature"] == "gene", :]
    )
    data_gene_annotation.drop(
        labels="feature",
        axis="columns",
        inplace=True
    )
    utility.print_terminal_partition(level=2)
    print("count of reference genes: " + str(data_gene_annotation.shape[0]))

    # Select entries for genes on nuclear autosomes (non-sex chromosomes).
    # There are 2422 genes on the X chromosome, of which 850 encode proteins.
    # There are 567 genes on the Y chromosome, of which 66 encode proteins.
    # There are 37 genes on the mitochondrial chromosome, of which 13 encode
    # proteins.
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["seqname"] != "chrX", :
        ]
    )
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["seqname"] != "chrY", :
        ]
    )
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["seqname"] != "chrM", :
        ]
    )
    utility.print_terminal_partition(level=2)
    print(
        "count of reference genes on nuclear (not mitochondrial) autosomes: " +
        str(data_gene_annotation.shape[0])
    )

    # Select entries for genes that encode proteins.
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["gene_type"] == "protein_coding", :
        ]
    )
    utility.print_terminal_partition(level=2)
    print(
        "count of reference genes that encode proteins: " +
        str(data_gene_annotation.shape[0])
    )

    # Select entries for genes that have signal in GTEx data.
    data_gene_annotation = data_gene_annotation.loc[
        data_gene_annotation.index.isin(genes_gtex), :
    ]
    genes_gtex_annotation = data_gene_annotation.index.to_list()

    # Report.
    utility.print_terminal_partition(level=2)
    print(
        "count of reference genes that have signal in GTEx data: " +
        str(data_gene_annotation.shape[0])
    )
    utility.print_terminal_partition(level=2)
    print(data_gene_annotation.iloc[0:10, 0:7])

    # Compile information.
    bin = dict()
    bin["data"] = data_gene_annotation
    bin["genes_gtex_annotation"] = genes_gtex_annotation
    # Return information.
    return bin


def select_organize_genes_annotations(
    stringency=None,
    dock=None,
):
    """
    Selects and organizes information about genes.

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
    """

    # Remove previous files to avoid version or batch confusion.
    path_selection = os.path.join(dock, "selection", str(stringency))
    path_gene_annotation = os.path.join(path_selection, "gene_annotation")
    utility.remove_directory(path=path_gene_annotation)

    # Read source information from file.
    source = read_source_gene_annotation(dock=dock)

    ##########
    ##########
    ##########
    # Select protein-coding genes on nuclear autosomes

    # Select genes that encode proteins and are on autosomes.
    # Count of genes: 18353
    utility.print_terminal_partition(level=1)
    print("Selection of genes.")

    # GTEx reference.
    utility.print_terminal_partition(level=2)
    print("GTEx reference gene annotations.")
    # GTEx reference: 18380 genes on autosomes that encode proteins.
    # GTEx reference: 18380 genes with signal.
    bin_gtex = select_gene_annotation(
        genes_gtex=source["genes_gtex"],
        data_gene_annotation=source["data_gene_annotation_gtex"],
    )

    # GENCODE reference.
    utility.print_terminal_partition(level=2)
    print("GENCODE reference gene annotations.")
    # GENCODE reference: 19027 genes on autosomes that encode proteins.
    # GENCODE reference: 18324 genes with signal.
    bin_gencode = select_gene_annotation(
        genes_gtex=source["genes_gtex"],
        data_gene_annotation=source["data_gene_annotation_gencode"],
    )

    # Continue analyses with selection of genes by GENCODE.
    utility.print_terminal_partition(level=2)
    print(
        "count of protein-coding genes with signals: " +
        str(len(bin_gencode["genes_gtex_annotation"]))
    )
    utility.print_terminal_partition(level=3)

    # Compile information.
    information = dict()
    information["data_gene_annotation_gtex"] = bin_gtex["data"]
    information["data_gene_annotation_gencode"] = bin_gencode["data"]
    information["genes_gtex_annotation"] = bin_gencode["genes_gtex_annotation"]

    # Write product information to file.
    write_product_gene_annotation(
        stringency=stringency,
        dock=dock,
        information=information
    )
    pass


##########
##########
##########
# Samples, genes, signals


def read_source_selection(
    stringency=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_assembly_sample = os.path.join(dock, "assembly", "sample")
    path_assembly_signal = os.path.join(dock, "assembly", "signal")
    path_samples_tissues_persons = os.path.join(
        path_assembly_sample, "data_samples_tissues_persons.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly_signal, "data_gene_signal.feather"
    )
    path_samples_gtex = os.path.join(
        path_assembly_signal, "samples_gtex.pickle"
    )

    path_selection = os.path.join(dock, "selection", str(stringency))
    path_gene_annotation = os.path.join(path_selection, "gene_annotation")
    path_genes_gtex_annotation = os.path.join(
        path_gene_annotation, "genes_gtex_annotation.pickle"
    )

    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_signal_feather = pandas.read_feather(
        path=path_gene_signal,
    )
    data_gene_signal = organize_data_axes_indices(
        data=data_gene_signal_feather
    )
    with open(path_samples_gtex, "rb") as file_source:
        samples_gtex = pickle.load(file_source)
    with open(path_genes_gtex_annotation, "rb") as file_source:
        genes_gtex_annotation = pickle.load(file_source)


    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_signal": data_gene_signal,
        "samples_gtex": samples_gtex,
        "genes_gtex_annotation": genes_gtex_annotation,
    }


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
    # Organize data.
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
    return data


def summarize_samples_genes(
    samples_gtex=None,
    genes_gtex_annotation=None,
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Summarize selection of samples and genes.

    arguments:
        samples_gtex (list<str>): identifiers of samples with signals in GTEx
            data
        genes_gtex_annotation (list<str>): identifiers of protein-coding genes
            on nuclear autosomes with signals in GTEx data
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_signal (object): Pandas data frame of genes' counts or
            signals across samples

    raises:

    returns:

    """

    # Copy data.
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    data_gene_signal = data_gene_signal.copy(deep=True)

    utility.print_terminal_partition(level=1)
    print("Summary of selection of samples and genes.")

    utility.print_terminal_partition(level=2)
    print("Count of samples from GTEx: " + str(len(samples_gtex)))
    utility.print_terminal_partition(level=2)
    print(
        "Count of protein-coding genes on nuclear autosomes with signals in " +
        "GTEx: " + str(len(genes_gtex_annotation))
    )

    # Counts of samples and genes.
    utility.print_terminal_partition(level=2)
    print("Genes' signals...")
    print("Count of samples: " + str(data_gene_signal.shape[1]))
    print("Count of genes: " + str(data_gene_signal.shape[0]))

    # Major tissue categories.
    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_signal.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_sample=data_transposition,
    )

    groups = data_factor.groupby(level=["person"])
    print("Count of persons: " + str(len(groups)))

    groups = data_factor.groupby(level=["tissue_minor"])
    print("Count of minor tissues: " + str(len(groups)))

    groups = data_factor.groupby(level=["tissue_major"])
    print("Count of major tissues: " + str(len(groups)))

    print("Count of persons with samples for each major tissue.")
    records = []
    for name, group in groups:
        data = group.reset_index(
            level=None,
            inplace=False
        )
        data.drop_duplicates(
            subset=["person", "tissue_major"],
            keep="first",
            inplace=True,
        )
        # Compile information.
        record = dict()
        record["tissue_major"] = name
        record["persons"] = data.shape[0]
        records.append(record)

    data_summary = utility.convert_records_to_dataframe(
        records=records
    )
    print(data_summary)

    data_factor.reset_index(
        level=["tissue_major"],
        inplace=True
    )
    data_adipose = data_factor.loc[data_factor["tissue_major"] == "adipose", :]
    #print(data_adipose)

    pass


def select_signals_by_genes(
    genes=None,
    data_gene_signal=None,
):
    """
    Selects genes that encode proteins.

    arguments:
        genes (list<str>): identifiers of genes
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

    # Copy data.
    data_gene_signal = data_gene_signal.copy(deep=True)
    # Filter gene signals.
    data_gene_signal_annotation = data_gene_signal.loc[
        genes, :
    ]
    # Return information.
    return data_gene_signal_annotation


def exclude_samples_by_tissues(
    tissues_minor=None,
    tissues_major=None,
    data_samples_tissues_persons=None,
):
    """
    Exclude samples for specific major and minor tissues.

    arguments:
        tissues_minor (list<str>): names of minor tissues to exclude
        tissues_major (list<str>): names of major tissues to exclude
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_minor = data_samples.loc[
        ~data_samples["tissue_minor"].isin(tissues_minor),
        :
    ]
    data_major = data_minor.loc[
        ~data_minor["tissue_major"].isin(tissues_major),
        :
    ]
    # Return information.
    return data_major


def include_samples_by_tissues(
    tissues_major=None,
    data_samples_tissues_persons=None,
):
    """
    Include samples for specific major and minor tissues.

    arguments:
        tissues_major (list<str>): names of major tissues to exclude
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_major = data_samples.loc[
        data_samples["tissue_major"].isin(tissues_major), :
    ]
    # Return information.
    return data_major


def select_samples_persons_by_tissues_count(
    count=None,
    data_samples_tissues_persons=None,
):
    """
    Select persons for whom samples are available for a minimal count of
        tissues.

    arguments:
        count (int): Minimal count of major tissues for which each person must
            have samples
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across samples
    """

    # Select persons by tissues.
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    # Count unique major tissues per person.
    data_tissues_per_person = utility.count_data_factors_groups_elements(
        factors=["person"],
        element="tissue_major",
        count="counts",
        data=data_samples_tissues_persons,
    )
    # Select persons by counts of tissues for which they have samples.
    data_persons_selection = data_tissues_per_person.loc[
        (data_tissues_per_person["counts"] >= count), :
    ]
    # Extract identifiers of persons.
    persons = utility.collect_unique_elements(
        elements_original=data_persons_selection["person"].to_list()
    )
    # Select samples by persons.
    data_samples_selection = data_samples_tissues_persons.loc[
        data_samples_tissues_persons["person"].isin(persons), :
    ]
    # Return information.
    return data_samples_selection


def report_selection_samples_persons_by_tissues(
    samples_gtex=None,
    genes_gtex_annotation=None,
    samples_tissues_persons=None,
    data_samples_original=None,
    data_samples_tissues=None,
    data_samples_persons=None,
    data_gene_signal=None,
):
    """
    Select samples and persons by tissues.

    arguments:
        samples_gtex (list<str>): identifiers of samples with signals in GTEx
            data
        genes_gtex_annotation (list<str>): identifiers of protein-coding genes
            on nuclear autosomes with signals in GTEx data
        samples_tissues_persons (list<str>): identifiers of samples from
            persons with adequate coverage of tissues of interest
        data_samples_original (object): Pandas data frame of persons and
            tissues across samples
        data_samples_tissues (object): Pandas data frame of persons and
            tissues across samples
        data_samples_persons (object): Pandas data frame of persons and
            tissues across samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
    """

    # Copy data.
    data_samples_persons = data_samples_persons.copy(deep=True)

    # Summarize original counts of samples and genes.
    data_gene_signal_temporary = select_signals_by_samples(
        samples=data_samples_tissues.index.to_list(),
        data_gene_signal=data_gene_signal,
    )
    utility.print_terminal_partition(level=1)
    print(
        "Summary after selection of samples by exclusion and inclusion " +
        "tissues."
    )
    summarize_samples_genes(
        samples_gtex=samples_gtex,
        genes_gtex_annotation=genes_gtex_annotation,
        data_samples_tissues_persons=data_samples_tissues,
        data_gene_signal=data_gene_signal_temporary,
    )
    # Select final samples on basis of persons' eligibility.
    utility.print_terminal_partition(level=1)
    print("Selection of samples by persons' eligibility.")
    print("Persons must have samples for adequate count of tissues.")
    # Count tissues per person.
    data_tissues_per_person_initial = (
        utility.count_data_factors_groups_elements(
            factors=["person"],
            element="tissue_major",
            count="counts",
            data=data_samples_original,
    ))
    utility.print_terminal_partition(level=2)
    print("Tissues per person in the original samples before any selection.")
    mean_tissues = data_tissues_per_person_initial["counts"].mean()
    print("Mean tissues per person (initial): " + str(mean_tissues))
    count_persons = len(data_tissues_per_person_initial["counts"].to_list())
    print("Count of persons (initial): " + str(count_persons))

    # Count tissues per person.
    data_tissues_per_person_final = (
        utility.count_data_factors_groups_elements(
            factors=["person"],
            element="tissue_major",
            count="counts",
            data=data_samples_persons,
    ))
    utility.print_terminal_partition(level=2)
    print(
        "Tissues per person in the samples after selection by tissue and by " +
        "count of tissues per person."
    )
    mean_tissues = data_tissues_per_person_final["counts"].mean()
    print("Mean tissues per person (final): " + str(mean_tissues))
    count_persons = len(data_tissues_per_person_final["counts"].to_list())
    print("Count of persons (final): " + str(count_persons))

    pass


def select_samples_persons_by_tissues(
    stringency=None,
    report=None,
    samples_gtex=None,
    genes_gtex_annotation=None,
    data_samples_tissues_persons=None,
    data_gene_signal=None,
):
    """
    Select samples and persons by tissues.

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        report (bool): whether to print reports about the selection
        samples_gtex (list<str>): identifiers of samples with signals in GTEx
            data
        genes_gtex_annotation (list<str>): identifiers of protein-coding genes
            on nuclear autosomes with signals in GTEx data
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (dict): information about selection of samples by tissues
    """

    # Select persons by tissues.
    data_samples_original = data_samples_tissues_persons.copy(deep=True)
    data_gene_signal = data_gene_signal.copy(deep=True)

    # Select samples from persons with at least 10 out of 20 sex-neutral
    # tissues.
    # Select samples by tissues of interest.
    # Exclude samples for cell lines that are not whole tissues.
    # Exclude samples for tissues with coverage by <100 samples.
    data_samples_exclusion = exclude_samples_by_tissues(
        tissues_minor=["fibroblast", "lymphocyte"],
        tissues_major=["bladder", "cervix", "fallopius", "kidney"],
        data_samples_tissues_persons=data_samples_original,
    )
    if stringency == "loose":
        # Select all tissues with adequate coverage.
        tissues_selection = [
            "adipose", "adrenal", "artery", "blood", "brain", "breast",
            "colon", "esophagus", "heart", "intestine", "liver", "lung",
            "muscle", "nerve", "ovary", "pancreas", "pituitary", "prostate",
            "salivary", "skin", "spleen", "stomach", "testis", "thyroid",
            "uterus", "vagina",
        ]
    elif stringency == "tight":
        # Select sex-neutral tissues.
        tissues_selection = [
            "adipose", "adrenal", "artery", "blood", "brain", "colon",
            "esophagus", "heart", "intestine", "liver", "lung", "muscle",
            "nerve", "pancreas", "pituitary", "salivary", "skin", "spleen",
            "stomach", "thyroid",
        ]
    data_samples_tissues = include_samples_by_tissues(
        tissues_major=tissues_selection,
        data_samples_tissues_persons=data_samples_exclusion,
    )

    # Select persons with adequate sample coverage of multiple tissues.
    data_samples_persons = select_samples_persons_by_tissues_count(
        count=10,
        data_samples_tissues_persons=data_samples_tissues,
    )
    samples_tissues_persons = utility.collect_unique_elements(
        elements_original=data_samples_persons.index.to_list()
    )

    # Report.
    if report:
        report_selection_samples_persons_by_tissues(
            samples_gtex=samples_gtex,
            genes_gtex_annotation=genes_gtex_annotation,
            samples_tissues_persons=samples_tissues_persons,
            data_samples_original=data_samples_original,
            data_samples_tissues=data_samples_tissues,
            data_samples_persons=data_samples_persons,
            data_gene_signal=data_gene_signal,
        )
    # Compile information.
    bin = dict()
    bin["samples"] = samples_tissues_persons
    bin["data_samples_tissues_persons"] = data_samples_persons
    # Return information.
    return bin


def select_signals_by_samples(
    samples=None,
    data_gene_signal=None
):
    """
    Selects genes' signals by samples.

    Data format should have genes across rows and samples across columns.

    arguments:
        samples (list<str): identifiers of samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    # Copy data.
    data_gene_signal = data_gene_signal.copy(deep=True)
    # Select samples on the basis of tissues.
    data_gene_signal_sample = data_gene_signal.loc[
        :, data_gene_signal.columns.isin(samples)
    ]
    # Return information.
    return data_gene_signal_sample


def normalize_collect_report_gene_signals(
    data_gene_signal=None,
    threshold=None,
    report=None,
):
    """
    Normalizes genes' signals by logarithmic transformation and collects
    values.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples
        threshold (float): value for threshold on genes' signals
        report (bool): whether to print reports about the selection

    raises:

    returns:
        (array): values of genes' signals
    """

    data = data_gene_signal.copy(deep=True)
    # Transform genes' signals to logarithmic space.
    # Transform values of genes' signals to base-2 logarithmic space.
    data_log = utility.calculate_pseudo_logarithm_signals(
        pseudo_count=1.0,
        base=2.0,
        data=data,
    )
    # Flatten signals to one dimensional array.
    signals = data_log.to_numpy().flatten()

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Normalize and collect distribution of genes' signals.")
        print("Threshold: " + str(threshold) + " = log2(1.0 pseudo + 0.1 TPM)")
        array_threshold = signals[signals >= threshold]
        totals = signals.size
        passes = array_threshold.size
        print(
            "Percent signals >= threshold: " +
            str(round((passes / totals) * 100)) + "%"
        )
    # Return information.
    return signals


def select_samples_genes_by_signals_coverage(
    threshold=None,
    proportion_gene=None,
    proportion_sample=None,
    data_gene_signal=None,
    report=None,
):
    """
    Selects samples and genes by signals beyond threshold.

    Data format should have genes across rows and samples across columns.

    arguments:
        threshold (float): minimal gene's signal for a sample
        proportion_gene (float): minimal proportion of samples that must have
            signals beyond threshold for a gene to pass
        proportion_sample (float): minimal proportion of genes that must have
            signals beyond threshold for a sample to pass
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples
        report (bool): whether to print reports about the selection

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    # Filter genes by their signals across samples.
    # Filter to keep only genes with signals beyond threshold in proportion of
    # samples.
    data_row = utility.filter_rows_columns_by_threshold_proportion(
        data=data_gene_signal,
        dimension="row",
        threshold=threshold,
        proportion=proportion_gene
    )

    # Filter samples by their signals across genes.
    # Filter to keep only samples with signals beyond threshold in proportion
    # of genes.
    data_column = utility.filter_rows_columns_by_threshold_proportion(
        data=data_row,
        dimension="column",
        threshold=threshold,
        proportion=proportion_sample
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Selection of samples and genes by proportion of signals beyond " +
            "threshold."
        )
        utility.print_terminal_partition(level=3)
        print("shape of data after signal filter on genes...")
        print(data_row.shape)
        utility.print_terminal_partition(level=3)
        print("shape of data after signal filter on samples...")
        print(data_column.shape)

    # Return information.
    return data_column


def report_selection_samples_genes_by_signals(
    samples_gtex=None,
    genes_gtex_annotation=None,
    samples_tissues_persons=None,
    data_samples_original=None,
    data_samples_tissues=None,
    data_samples_persons=None,
    data_gene_signal=None,
):
    """
    Select samples and persons by tissues.

    arguments:
        samples_gtex (list<str>): identifiers of samples with signals in GTEx
            data
        genes_gtex_annotation (list<str>): identifiers of protein-coding genes
            on nuclear autosomes with signals in GTEx data
        samples_tissues_persons (list<str>): identifiers of samples from
            persons with adequate coverage of tissues of interest
        data_samples_original (object): Pandas data frame of persons and
            tissues across samples
        data_samples_tissues (object): Pandas data frame of persons and
            tissues across samples
        data_samples_persons (object): Pandas data frame of persons and
            tissues across samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
    """

    # Copy data.
    data_samples_persons = data_samples_persons.copy(deep=True)

    # Summarize original counts of samples and genes.
    data_gene_signal_temporary = select_signals_by_samples(
        samples=data_samples_tissues.index.to_list(),
        data_gene_signal=data_gene_signal,
    )
    utility.print_terminal_partition(level=1)
    print(
        "Summary after selection of samples by exclusion and inclusion " +
        "tissues."
    )
    summarize_samples_genes(
        samples_gtex=samples_gtex,
        genes_gtex_annotation=genes_gtex_annotation,
        data_samples_tissues_persons=data_samples_tissues,
        data_gene_signal=data_gene_signal_temporary,
    )
    # Select final samples on basis of persons' eligibility.
    utility.print_terminal_partition(level=1)
    print("Selection of samples by persons' eligibility.")
    print("Persons must have samples for adequate count of tissues.")
    # Count tissues per person.
    data_tissues_per_person_initial = (
        utility.count_data_factors_groups_elements(
            factors=["person"],
            element="tissue_major",
            count="counts",
            data=data_samples_original,
    ))
    utility.print_terminal_partition(level=2)
    print("Tissues per person in the original samples before any selection.")
    mean_tissues = data_tissues_per_person_initial["counts"].mean()
    print("Mean tissues per person (initial): " + str(mean_tissues))
    count_persons = len(data_tissues_per_person_initial["counts"].to_list())
    print("Count of persons (initial): " + str(count_persons))

    # Count tissues per person.
    data_tissues_per_person_final = (
        utility.count_data_factors_groups_elements(
            factors=["person"],
            element="tissue_major",
            count="counts",
            data=data_samples_persons,
    ))
    utility.print_terminal_partition(level=2)
    print(
        "Tissues per person in the samples after selection by tissue and by " +
        "count of tissues per person."
    )
    mean_tissues = data_tissues_per_person_final["counts"].mean()
    print("Mean tissues per person (final): " + str(mean_tissues))
    count_persons = len(data_tissues_per_person_final["counts"].to_list())
    print("Count of persons (final): " + str(count_persons))

    pass


def select_samples_genes_by_signals(
    stringency=None,
    report=None,
    data_gene_signal=None,
):
    """
    Select samples and persons by tissues.

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        report (bool): whether to print reports about the selection
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (dict): information about selection of samples and genes by signals
    """

    # Select persons by tissues.
    data_gene_signal = data_gene_signal.copy(deep=True)

    # Select genes and samples by their signal coverage.
    utility.print_terminal_partition(level=1)
    print("Selection of genes and samples by signal coverage.")
    utility.print_terminal_partition(level=1)
    # Fill missing signals with values of zero.
    data_gene_signal_fill = data_gene_signal.fillna(
        value=0.0,
        inplace=False,
    )
    # Collect initial distribution of genes' signals for selection of signal
    # threshold.
    signals_initial = normalize_collect_report_gene_signals(
        data_gene_signal=data_gene_signal_fill,
        threshold=math.log((0.1 + 1.0), 2), # pseudo count 1.0, 0.1 TPM
        report=report,
    )
    # Select genes and samples by signals and coverage.
    if stringency == "loose":
        # Each gene must have signal greater than or equal to 0.1 in at least
        # 1% of samples.
        # Each sample must have signal greater than or equal to 0.1 in at least
        # 10% of genes.
        data_gene_signal_selection = select_samples_genes_by_signals_coverage(
            threshold=0.1,
            proportion_gene=0.01, # 0.01
            proportion_sample=0.1, # 0.1
            data_gene_signal=data_gene_signal_fill,
            report=report,
        )
    elif stringency == "tight":
        # Each gene must have signal greater than or equal to 0.1 TPM in at
        # least 50% of samples.
        # This threshold removes genes that have strong specificity to a few
        # tissues, which might introduce error in subsequent analyses.
        # This threshold also reduces the count of genes for subsequent
        # analyses, hence increasing statistical power.
        # Each sample must have signal greater than or equal to 0.1 TPM in at
        # least 50% of genes.
        # TCW 8 February 2020
        # threshold ... proportion_gene ... proportion_sample ... genes ... samples
        # 0.0           0.5                 0.5
        # 0.1           0.1                 0.1
        # 0.1           0.5                 0.1
        # 0.1           0.5                 0.5                   15130     11775 <-- Priority!
        # 1.0           0.1                 0.1
        # 1.0           0.5                 0.1
        # 1.0           0.5                 0.5
        data_gene_signal_selection = select_samples_genes_by_signals_coverage(
            threshold=0.1,
            proportion_gene=0.5, # 0.5
            proportion_sample=0.5, # 0.5
            data_gene_signal=data_gene_signal_fill,
            report=report,
        )
    # Collect initial distribution of genes' signals for selection of signal
    # threshold.
    signals_final = normalize_collect_report_gene_signals(
        data_gene_signal=data_gene_signal_selection,
        threshold=math.log((0.1 + 1.0), 2), # pseudo count 1.0, 0.1 TPM
        report=report,
    )

    # Compile information.
    bin = dict()
    bin["signals_initial"] = signals_initial
    bin["signals_final"] = signals_final
    bin["data_gene_signal"] = data_gene_signal_selection
    # Return information.
    return bin


def extract_final_genes_samples_tissues_persons(
    report=None,
    data_gene_signal=None,
    data_samples_tissues_persons=None,
):
    """
    Extracts final selections of genes, samples, tissues, and persons from
    data signals across samples and genes.

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        report (bool): whether to print reports
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (dict): information about selection of samples and genes
    """

    # Select persons by tissues.
    data_gene_signal = data_gene_signal.copy(deep=True)
    data_samples = data_samples_tissues_persons.copy(deep=True)

    # Extract identifiers of genes and samples.
    genes = utility.collect_unique_elements(
        elements_original=data_gene_signal.index.to_list()
    )
    samples = utility.collect_unique_elements(
        elements_original=data_gene_signal.columns.to_list()
    )
    # Select samples, tissues, and persons.
    data_samples_selection = data_samples_tissues_persons.loc[
        data_samples_tissues_persons.index.isin(samples), :
    ]
    # Extract identifiers of persons and tissues.
    persons = utility.collect_unique_elements(
        elements_original=(
            data_samples_selection["person"].to_list()
        )
    )
    tissues = utility.collect_unique_elements(
        elements_original=(
            data_samples_selection["tissue_major"].to_list()
        )
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Final selections of genes, samples, persons, and tissues.")
        utility.print_terminal_partition(level=3)
        print("count of genes: " + str(len(genes)))
        print("count of samples: " + str(len(samples)))
        print("count of persons: " + str(len(persons)))
        print("count of tissues: " + str(len(tissues)))
        utility.print_terminal_partition(level=3)

    # Compile information.
    bin = dict()
    bin["genes"] = genes
    bin["samples"] = samples
    bin["persons"] = persons
    bin["tissues"] = tissues
    bin["data_samples_tissues_persons"] = data_samples_selection
    # Return information.
    return bin


def select_organize_samples_genes_signals(
    stringency=None,
    dock=None,
):
    """
    Selects samples and genes by signals beyond threshold.

    Data format should have genes across rows and samples across columns.

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
    """

    # Remove previous files to avoid version or batch confusion.
    path_selection = os.path.join(dock, "selection", str(stringency))
    path_samples_genes_signals = os.path.join(
        path_selection, "samples_genes_signals"
    )
    utility.remove_directory(path=path_samples_genes_signals)

    # Read source information from file.
    source = read_source_selection(
        stringency=stringency,
        dock=dock,
    )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        samples_gtex=source["samples_gtex"],
        genes_gtex_annotation=source["genes_gtex_annotation"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=source["data_gene_signal"],
    )

    # Select signals for protein-coding genes on nuclear autosomes.
    data_gene_signal_annotation = select_signals_by_genes(
        genes=source["genes_gtex_annotation"],
        data_gene_signal=source["data_gene_signal"],
    )

    # Select samples and persons by tissues.
    bin_tissues = select_samples_persons_by_tissues(
        stringency=stringency,
        report=True,
        samples_gtex=source["samples_gtex"],
        genes_gtex_annotation=source["genes_gtex_annotation"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal_annotation,
    )

    # Select samples in signal data.
    data_gene_signal_sample = select_signals_by_samples(
        samples=bin_tissues["samples"],
        data_gene_signal=data_gene_signal_annotation,
    )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        samples_gtex=source["samples_gtex"],
        genes_gtex_annotation=source["genes_gtex_annotation"],
        data_samples_tissues_persons=(
            bin_tissues["data_samples_tissues_persons"]
        ),
        data_gene_signal=data_gene_signal_sample,
    )

    # Select samples and genes by their signals.
    bin_signal = select_samples_genes_by_signals(
        stringency=stringency,
        report=True,
        data_gene_signal=data_gene_signal_sample,
    )

    # Extract final selections of samples, genes, tissues, and persons.
    bin_selection = extract_final_genes_samples_tissues_persons(
        data_gene_signal=bin_signal["data_gene_signal"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        report=True,
    )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        samples_gtex=source["samples_gtex"],
        genes_gtex_annotation=source["genes_gtex_annotation"],
        data_samples_tissues_persons=(
            bin_selection["data_samples_tissues_persons"]
        ),
        data_gene_signal=bin_signal["data_gene_signal"],
    )

    # Compile information.
    information = dict()
    information["data_gene_signal"] = bin_signal["data_gene_signal"]
    information["signals_initial"] = bin_signal["signals_initial"]
    information["signals_final"] = bin_signal["signals_final"]

    information["data_samples_tissues_persons"] = (
        bin_selection["data_samples_tissues_persons"]
    )
    information["genes"] = bin_selection["genes"]
    information["samples"] = bin_selection["samples"]
    information["persons"] = bin_selection["persons"]
    information["tissues"] = bin_selection["tissues"]
    # Write product information to file.
    write_product_samples_genes_signals(
        stringency=stringency,
        dock=dock,
        information=information
    )
    pass


##########
##########
##########
# Persons' properties


def read_source_persons_properties(
    stringency=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_assembly_sample = os.path.join(dock, "assembly", "sample")
    path_samples_tissues_persons = os.path.join(
        path_assembly_sample, "data_samples_tissues_persons.pickle"
    )

    path_selection = os.path.join(dock, "selection", str(stringency))
    path_samples_genes_signals = os.path.join(
        path_selection, "samples_genes_signals"
    )
    path_samples_tissues_persons_selection = os.path.join(
        path_samples_genes_signals, "data_samples_tissues_persons.pickle"
    )

    path_persons = os.path.join(
        dock, "selection", str(stringency), "samples_genes_signals",
        "persons.pickle"
    )

    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_samples_tissues_persons_selection = pandas.read_pickle(
        path_samples_tissues_persons_selection
    )
    with open(path_persons, "rb") as file_source:
        persons = pickle.load(file_source)

    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_samples_tissues_persons_selection": (
            data_samples_tissues_persons_selection
        ),
        "persons_selection": persons,
    }


def extract_genotypes(
    data_samples_tissues_persons=None,
):
    """
    Extracts persons' genotypes.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues across samples

    raises:

    returns:
        (object): Pandas data frame of persons' genotypes

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    # Organize data.
    data_samples.reset_index(
        level=None,
        inplace=True
    )
    genotypes = [
        "genotype_1",
        "genotype_2",
        "genotype_3",
        "genotype_4",
        "genotype_5",
        "genotype_6",
        "genotype_7",
        "genotype_8",
        "genotype_9",
        "genotype_10",
        "genotype_11",
        "genotype_12",
        "genotype_13",
        "genotype_14",
        "genotype_15",
        "genotype_16",
        "genotype_17",
        "genotype_18",
        "genotype_19",
        "genotype_20",
        "genotype_21",
        "genotype_22",
        "genotype_23",
        "genotype_24",
        "genotype_25",
    ]
    columns = copy.deepcopy(genotypes)
    columns.append("person")
    data_genotypes = data_samples.loc[
        :, data_samples.columns.isin(columns)
    ]
    data_genotypes.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_genotypes.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    # Return information.
    return data_genotypes


def extract_persons_genotype(
    persons=None,
    data_samples_tissues_persons=None,
):
    """
    Extracts identifiers of persons with valid genotypes.

    arguments:
        persons (list<str>): identifiers of persons for which to consider
            genotype
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (list<str>): identifiers of persons with valid genotypes

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    # Organize data.
    data_genotypes = extract_genotypes(
        data_samples_tissues_persons=data_samples,
    )
    data_genotypes.dropna(
        subset=None,
        axis="index",
        how="any",
        thresh=1,
        inplace=True,
    )
    # Extract identifiers of persons.
    persons_genotype = utility.collect_unique_elements(
        elements_original=data_genotypes.index.to_list()
    )
    # Return information.
    return persons_genotype


def extract_persons_selection_genotype(
    persons_selection=None,
    data_samples_tissues_persons=None,
    report=None,
):
    """
    Extracts identifiers of persons with valid genotypes.

    arguments:
        persons_selection (list<str>): identifiers of persons from selection
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of persons with valid genotypes

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    # Organize data.
    persons_gtex = utility.collect_unique_elements(
        elements_original=data_samples["person"].to_list()
    )
    # Extract identifiers of persons.
    persons_gtex_genotype = extract_persons_genotype(
        persons=persons_gtex,
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    persons_selection_genotype = extract_persons_genotype(
        persons=persons_selection,
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Selection of persons with valid genotypes.")
        utility.print_terminal_partition(level=2)
        print("Count of persons in GTEx cohort: " + str(len(persons_gtex)))
        print(
            "Count of GTEx persons with valid genotypes: " +
            str(len(persons_gtex_genotype))
        )
        utility.print_terminal_partition(level=3)
        print("Count of selection persons: " + str(len(persons_selection)))
        print(
            "Count of selection persons with valid genotypes: " +
            str(len(persons_selection_genotype))
        )
        utility.print_terminal_partition(level=2)
    # Return information.
    return persons_selection_genotype


def impute_persons_genotypes(
    data_samples_tissues_persons_selection=None,
):
    """
    Imputes missing genotypes.

    Calculate means of principal components on genotype across selection of
    persons.

    arguments:
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues across selection samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across selection
            samples

    """

    # Copy data.
    data_samples_selection = data_samples_tissues_persons_selection.copy(
        deep=True
    )
    # Organize data.
    data_genotypes = extract_genotypes(
        data_samples_tissues_persons=data_samples_selection,
    )
    # Calculate values for imputation.
    # Calculate mean values of principal components on all available genotypes.
    imputations = data_genotypes.aggregate(
        lambda x: x.mean()
    )
    #imputations = series_imputation.to_dict()

    # Insert imputations to selections of persons and tissues.
    # This step should only fill missing values in genotype columns.
    # "Values not in the dict/Series/DataFrame will not be filled."
    if False:
        data_copy.apply(
            lambda x: x.fillna(
                imputations[x.name],
                inplace=True,
            ),
            axis="index",
        )
    data_samples_selection.fillna(
        value=imputations,
        #axis="columns",
        inplace=True,
    )
    # Return information.
    return data_samples_selection


def define_variables():
    """
    Defines a list of variables' names for regression analysis.

    arguments:

    raises:

    returns:
        (dict<list<str>>): names of independent variables for regression

    """


    # Logical variables to convert to binary.
    binary = [
        "respiration",
        "inflammation",
        "infection",
        "mononucleosis",
        "steroid",
        "ventilation",
        "refrigeration",
    ]
    # Categorical variables for which each person has multiple values.
    # Reduce dimensionality for these variables.
    category = [
        "tissues",
        "facilities",
        "batches_extraction",
        "batches_sequence",
    ]
    # Variables of raw values that might require standardization to adjust
    # scale.
    scale = [
        "sex_y",
        "age",
        "body",
        "climate",
        "hardiness",
        "respiration_binary",
        "inflammation_binary",
        "infection_binary",
        "mononucleosis_binary",
        "steroid_binary",
        "delay",
        "ventilation_duration",
        "refrigeration_duration",
        "refrigeration_binary",
    ]
    # Variables that relate to hypotheses of interest.
    hypothesis = [
        "sex_y_scale",
        "age_scale",
        "body_scale",
        "climate_scale",
        "hardiness_scale",
        "respiration_binary_scale",
        "inflammation_binary_scale",
        #"infection_binary_scale", # omit because most persons are True
        "mononucleosis_binary_scale",
        #"steroid_binary_scale", # omit because very few persons are True
        "ventilation_duration_scale",
    ]
    # Variables that relate to genotype.
    genotype = [
        "genotype_1",
        "genotype_2",
        "genotype_3",
        "genotype_4",
        "genotype_5",
        #"genotype_6", # <- omit
        #"genotype_7", # <- omit
        #"genotype_8", # <- omit
        #"genotype_9", # <- omit
        #"genotype_10", # <- omit
    ]
    # Variables that relate to technical methods.
    technique = [
        "refrigeration_binary_scale", # Use the binary since duration had so many missing values
        "delay_scale",

        "tissues_1",
        "tissues_2",
        "tissues_3",
        "tissues_4",
        "tissues_5",
        #"tissues_6", # <- omit
        #"tissues_7", # <- omit
        #"tissues_8", # <- omit
        #"tissues_9", # <- omit
        #"tissues_10", # <- omit

        #"facilities_1", # <- omit due to multicollinearity
        #"facilities_2", # <- omit due to multicollinearity
        #"facilities_3", # <- omit

        #"batches_extraction_1",
        #"batches_extraction_2",
        #"batches_extraction_3",
        #"batches_extraction_4", # <- omit
        #"batches_extraction_5", # <- omit
        #"batches_extraction_6", # <- omit
        #"batches_extraction_7", # <- omit
        #"batches_extraction_8", # <- omit
        #"batches_extraction_9", # <- omit
        #"batches_extraction_10", # <- omit

        #"batches_sequence_1",
        #"batches_sequence_2",
        #"batches_sequence_3",
        #"batches_sequence_4",
        #"batches_sequence_5",
        #"batches_sequence_6", # <- omit
        #"batches_sequence_7", # <- omit
        #"batches_sequence_8", # <- omit
        #"batches_sequence_9", # <- omit
        #"batches_sequence_10", # <- omit
    ]

    # Compile variables.

    # Regression.
    regression = list()
    regression.extend(hypothesis)
    regression.extend(genotype)
    regression.extend(technique)

    # Heritability.
    heritability_simple = list()
    heritability_simple.extend(technique)
    heritability_complex = list()
    heritability_complex.extend(hypothesis)
    heritability_complex.extend(technique)

    # Quantitative Trait Loci (QTL).
    trait = list()
    trait.extend(hypothesis)
    trait.extend(technique)

    # Compile information.
    information = dict()
    information["binary"] = binary
    information["category"] = category
    information["scale"] = scale
    information["hypothesis"] = hypothesis
    information["genotype"] = genotype
    information["technique"] = technique
    information["regression"] = regression
    #information["heritability_simple"] = heritability_simple
    #information["heritability_complex"] = heritability_complex
    #information["trait"] = trait

    # Return information.
    return information


def impute_samples_duration(
    indicator=None,
    duration=None,
    replacement=None,
    report=None,
    data_samples_tissues_persons=None,
):
    """
    Impute duration variables across samples.

    Relevant variables include indicators of a person's ventilation or
    refrigeration.

    If indicator is false and duration is missing, then impute with a value of
    0.0 for duration.
    If indicator is true but corresponding duration is missing, then
    impute with the minimal nonzero value of duration.

    arguments:
        indicator (str): name of column for Boolean indicator variable
        duration (str): name of column for float duration variable
        replacement (bool): whether to replace original duration columns
        report (bool): whether to print reports
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues across selection samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across samples

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    # Report.
    if report:
        # Copy data.
        data_report = data_samples_tissues_persons.copy(deep=True)
        # Organize data.
        data_report.reset_index(
            level=None,
            inplace=True
        )
        columns = [
            "person",
            indicator,
            duration,
        ]
        data_report = data_report.loc[
            :, data_report.columns.isin(columns)
        ]
        data_report.drop_duplicates(
            subset=None,
            keep="first",
            inplace=True,
        )
        data_report.set_index(
            ["person"],
            append=False,
            drop=True,
            inplace=True
        )
        counter_false = 0
        counter_true = 0
        counter_null = 0
        counter_records = 0
        records = data_report.to_dict("records")
        for record in records:
            counter_records += 1
            if (
                (record[indicator] == False) and (math.isnan(record[duration]))
            ):
                counter_false += 1
            elif (
                (record[indicator] == True) and (math.isnan(record[duration]))
            ):
                counter_true += 1
            elif (
                math.isnan(record[indicator]) and math.isnan(record[duration])
            ):
                counter_null += 1
        utility.print_terminal_partition(level=2)
        print("Extent of imputation in duration variables.")
        print("Percentages are relative to counts of persons.")
        utility.print_terminal_partition(level=3)
        print(
            "False indicators with missing duration: " +
            str(round((counter_false / counter_records) * 100, 2)) + "%"
        )
        print("Impute with 0.0... obviously.")
        utility.print_terminal_partition(level=3)
        print(
            "True indicators with missing duration: " +
            str(round((counter_true / counter_records) * 100, 2)) + "%"
        )
        print("Impute with minimum value of duration.")
        utility.print_terminal_partition(level=3)
        print(
            "Missing values of both indicator and duration: " +
            str(round((counter_null / counter_records) * 100, 2)) + "%"
        )
        print("Impute with false and zero.")
        pass
    # Impute.
    duration_imputation = str(duration + "_imputation")
    durations_nonzero = data_samples[duration].loc[
        data_samples[duration] > 0
    ].to_numpy()
    minimum = numpy.nanmin(durations_nonzero)
    if report:
        print("Minimal value of duration: " + str(minimum))
    data_samples[duration_imputation] = data_samples.apply(
        lambda row:
            0.0 if (
                (row[indicator] == False) and (math.isnan(row[duration]))
            ) else
            (minimum if (
                (row[indicator] == True) and (math.isnan(row[duration]))
            ) else
            (0.0 if (
                math.isnan(row[indicator]) and math.isnan(row[duration])
            ) else
            (row[duration]))),
        axis="columns",
    )
    # Organize data.
    # Replace original duration column.
    if replacement:
        data_samples.drop(
            labels=[duration],
            axis="columns",
            inplace=True
        )
        data_samples.rename(
            columns={
                duration_imputation: duration,
            },
            inplace=True,
        )
    # Return information.
    return data_samples


def correct_samples_duration(
    indicator=None,
    duration=None,
    replacement=None,
    report=None,
    data_samples_tissues_persons=None,
):
    """
    Correct duration variables across samples.

    Relevant variables include indicators of a person's ventilation or
    refrigeration.

    If duration is 0.0 hours, then change indicator to false.

    arguments:
        indicator (str): name of column for Boolean indicator variable
        duration (str): name of column for float duration variable
        replacement (bool): whether to replace original duration columns
        report (bool): whether to print reports
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues across selection samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across samples

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    # Report.
    if report:
        # Copy data.
        data_report = data_samples_tissues_persons.copy(deep=True)
        # Organize data.
        data_report.reset_index(
            level=None,
            inplace=True
        )
        columns = [
            "person",
            indicator,
            duration,
        ]
        data_report = data_report.loc[
            :, data_report.columns.isin(columns)
        ]
        data_report.drop_duplicates(
            subset=None,
            keep="first",
            inplace=True,
        )
        data_report.set_index(
            ["person"],
            append=False,
            drop=True,
            inplace=True
        )
        counter_error = 0
        counter_null = 0
        counter_records = 0
        records = data_report.to_dict("records")
        for record in records:
            counter_records += 1
            if (
                (record[indicator] == True) and (record[duration] == 0)
            ):
                counter_error += 1
            elif (
                math.isnan(record[indicator]) and (record[duration] == 0)
            ):
                counter_null += 1
        utility.print_terminal_partition(level=2)
        print("Extent of correction of inconsistency in duration variables.")
        print("Percentages are relative to counts of persons.")
        utility.print_terminal_partition(level=3)
        print(
            "True indicators but zero durations: " +
            str(round((counter_error / counter_records) * 100, 2)) + "%"
        )
        print("Prioritize duration variable, change indicator variable.")
        utility.print_terminal_partition(level=3)
        print(
            "Missing indicators but zero durations: " +
            str(round((counter_null / counter_records) * 100, 2)) + "%"
        )
        print("Prioritize duration variable, change indicator false.")
        utility.print_terminal_partition(level=3)
        pass
    # Impute.
    indicator_correction = str(indicator + "_correction")
    data_samples[indicator_correction] = data_samples.apply(
        lambda row:
            False if (
                (row[indicator] == True) and (row[duration] == 0)
            ) else
            (False if (
                math.isnan(row[indicator]) and (row[duration] == 0)
            ) else
            (row[indicator])),
        axis="columns",
    )
    # Organize data.
    # Replace original duration column.
    if replacement:
        data_samples.drop(
            labels=[indicator],
            axis="columns",
            inplace=True
        )
        data_samples.rename(
            columns={
                indicator_correction: indicator,
            },
            inplace=True,
        )
    # Return information.
    return data_samples


def impute_correct_samples_duration(
    indicator=None,
    duration=None,
    replacement=None,
    report=None,
    data_samples_tissues_persons=None,
):
    """
    Impute duration variables across samples.

    Relevant variables include indicators of a person's ventilation or
    refrigeration.

    If indicator is false and duration is missing, then impute with a value of
    0.0 for duration.
    If indicator is true but corresponding duration is missing, then
    impute with the minimal nonzero value of duration.

    arguments:
        indicator (str): name of column for Boolean indicator variable
        duration (str): name of column for float duration variable
        replacement (bool): whether to replace original duration columns
        report (bool): whether to print reports
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues across selection samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across samples

    """

    # Impute values
    data_imputation = impute_samples_duration(
        indicator=indicator,
        duration=duration,
        replacement=replacement,
        report=report,
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    # Correct values.
    # Execute correction after imputation so as to change indicator variable to
    # correspond with imputation value of zero for missing indicator.
    if not replacement:
        # Base corrections on any imputations.
        duration_correction = str(duration + "_imputation")
    else:
        duration_correction = duration
    data_correction = correct_samples_duration(
        indicator=indicator,
        duration=duration_correction,
        replacement=replacement,
        report=report,
        data_samples_tissues_persons=data_imputation,
    )
    # Return information.
    return data_correction


def organize_samples_season_variables(
    data_samples_tissues_persons=None,
):
    """
    Organize variables that relate to the season of death.

    arguments:
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues across selection samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across samples

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    # Climate.
    # Seasons are cyclic.
    # The climate variable represents the general climate in which the person
    # died.
    # Winter and summer seasons are opposite extremes in climate, with spring
    # and fall in between.
    data_samples["climate"] = data_samples["season"].apply(
        lambda value:
            0 if (value == "winter") else
            (1 if ((value == "spring") or (value == "fall")) else
            (2 if (value == "summer") else
            (float("nan"))))
    )
    # Season sequence.
    data_samples["season_sequence"] = data_samples["season"].apply(
        lambda value:
            0 if (value == "winter") else
            (1 if (value == "spring") else
            (2 if (value == "summer") else
            (3 if (value == "fall") else
            (float("nan")))))
    )
    # Season categorical binary variables.
    data_samples["winter"] = data_samples["season"].apply(
        lambda value: 1 if (value == "winter") else 0
    )
    data_samples["spring"] = data_samples["season"].apply(
        lambda value: 1 if (value == "spring") else 0
    )
    data_samples["summer"] = data_samples["season"].apply(
        lambda value: 1 if (value == "summer") else 0
    )
    data_samples["fall"] = data_samples["season"].apply(
        lambda value: 1 if (value == "fall") else 0
    )
    # Return information.
    return data_samples


def organize_samples_sex_variables(
    data_samples_tissues_persons=None,
):
    """
    Organize variables that relate to biological sex.

    arguments:
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues across selection samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across samples

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    # Sex.
    data_samples["sex_x"] = data_samples["sex_text"].apply(
        lambda value:
            2 if (value == "female") else
            (1 if (value == "male") else
            float("nan"))
    )
    data_samples["sex_y"] = data_samples["sex_text"].apply(
        lambda value:
            0 if (value == "female") else
            (1 if (value == "male") else
            float("nan"))
    )
    # Return information.
    return data_samples


def convert_samples_logical_variables_binary(
    variables=None,
    data_samples_tissues_persons=None,
):
    """
    Organize variables that relate to biological sex.

    arguments:
        variables (list<str>): names of logical variables
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues across selection samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across samples

    """

    # Copy data.
    data = data_samples_tissues_persons.copy(deep=True)
    # Iterate on variables.
    for variable in variables:
        variable_binary = str(variable + ("_binary"))
        data[variable_binary] = data[variable].apply(
            lambda value:
                0 if (value == False) else
                (1 if (value == True) else
                float("nan"))
        )
        pass
    # Return information.
    return data


def organize_samples_properties(
    variables=None,
    data_samples_tissues_persons_selection=None,
):
    """
    Organize persons' properties across samples.

    arguments:
        variables (dict<list<str>>): names of independent variables for
            regression
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues across selection samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues across selection
            samples

    """

    # Copy data.
    data_samples_selection = data_samples_tissues_persons_selection.copy(
        deep=True
    )
    # Organize data.

    # Ventilation.
    #print("Ventilation")
    data_ventilation = impute_correct_samples_duration(
        indicator="ventilation",
        duration="ventilation_duration",
        data_samples_tissues_persons=data_samples_selection,
        replacement=True,
        report=False,
    )
    # Refrigeration.
    #print("Refrigeration")
    data_refrigeration = impute_correct_samples_duration(
        indicator="refrigeration",
        duration="refrigeration_duration",
        data_samples_tissues_persons=data_ventilation,
        replacement=True,
        report=False,
    )
    # Season.
    data_season = organize_samples_season_variables(
        data_samples_tissues_persons=data_refrigeration,
    )
    # Sex.
    data_sex = organize_samples_sex_variables(
        data_samples_tissues_persons=data_season,
    )
    # Convert logical variables to binary representation for regressions.
    data_binary = convert_samples_logical_variables_binary(
        variables=variables["binary"],
        data_samples_tissues_persons=data_sex,
    )
    # Return information.
    return data_binary


def collect_persons_unique_values(
    values=None,
    delimiter=None,
):
    """
    Collects unique values from a person's samples.

    arguments:
        values (list<str>): values from a person's samples
        delimiter (str): character for separation of values

    raises:

    returns:
        (str): values

    """

    #values_valid = list(filter(lambda value: not math.isnan(value), values))
    #values_type = list(map(lambda value: str(value), values_valid))
    values_combination = delimiter.join(values)
    values_split = values_combination.split(delimiter)
    values_strip = list(map(lambda value: value.strip(), values_split))
    values_unique = utility.collect_unique_elements(
        elements_original=values_strip
    )
    return values_unique


def extract_persons_properties(
    data_samples_tissues_persons=None,
    report=None,
):
    """
    Extracts information about persons.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about persons

    """

    # Organize data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_samples.rename_axis(
        "",
        axis="columns",
        inplace=True,
    )
    data_samples.reset_index(
        level=None,
        inplace=True
    )
    data_samples.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )

    # Split data by person.
    groups = data_samples.groupby(
        level=["person"],
    )
    data_collection = pandas.DataFrame()
    for name, data_original in groups:
        data_original.reset_index(
            level=None,
            inplace=True
        )
        data_novel = data_original.copy(deep=True)
        data_novel.drop(
            labels=[
                "sample",
                "tissue_major",
                "tissue_minor",
                "facilities",
                "batch_extraction",
                "batches_sequence",
            ],
            axis="columns",
            inplace=True
        )
        data_novel.drop_duplicates(
            subset=None,
            keep="first",
            inplace=True,
        )
        delimiter = ","
        data_novel["tissues"] = delimiter.join(collect_persons_unique_values(
            values=data_original["tissue_major"].dropna().to_list(),
            delimiter=",",
        ))
        data_novel["facilities"] = delimiter.join(
            collect_persons_unique_values(
                values=data_original["facilities"].dropna().to_list(),
                delimiter=",",
        ))
        data_novel["batches_extraction"] = delimiter.join(
            collect_persons_unique_values(
                values=data_original["batch_extraction"].dropna().to_list(),
                delimiter=",",
        ))
        data_novel["batches_sequence"] = delimiter.join(
            collect_persons_unique_values(
                values=data_original["batches_sequence"].dropna().to_list(),
                delimiter=",",
        ))
        data_collection = data_collection.append(
            data_novel,
            ignore_index=True,
        )

    # Organize data.
    #data_collection.reindex()
    data_collection.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("data_persons_properties_raw")
        print(data_collection)
        utility.print_terminal_partition(level=2)
    # Return information.
    return data_collection


def expand_persons_categories_matrix(
    category=None,
    values=None,
    delimiter=None,
    data=None,
):
    """
    Expands a person's categorical properties to a binary matrix.

    arguments:
        category (str): name of a categorical property
        values (list<str>): all unique values of the categorical property
            across persons
        delimiter (str): character for separation of values
        data (object): Pandas data frame of persons' properties

    raises:

    returns:
        (object): Pandas data frame of a binary matrix

    """

    # Copy data.
    data = data.copy(deep=True)
    # Iterate on persons.
    persons = data.index.to_list()
    data_collection = pandas.DataFrame()
    for person in persons:
        # Create an empty matrix with values of zero for all categories.
        data_empty = pandas.DataFrame(
            numpy.int(0),
            index=[0],
            columns=values
        )
        # Adjust matrix to match person's values of the categorical property.
        data_person = data_empty.copy(deep=True)
        data_person["person"] = person
        # Iterate on the person's values of the categorical property.
        values_person_raw = data.at[person, category]
        values_person = values_person_raw.split(delimiter)
        for value in values_person:
            data_person[value] = 1
        # Collect person's values.
        data_collection = data_collection.append(
            data_person,
            ignore_index=True,
        )
    # Organize data.
    data_collection.reset_index(
        level=None,
        inplace=True
    )
    data_collection.drop(
        labels="index",
        axis="columns",
        inplace=True
    )
    data_collection.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    # Return information.
    return data_collection


def reduce_dimensionality_categorical_variable(
    variable=None,
    data_persons_properties_raw=None,
    report=None,
):
    """
    Organizes information about persons.

    arguments:
        variable (dict<list<str>>): names of categorical variable
        data_persons_properties_raw (object): Pandas data frame of persons'
            properties
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information from dimensionality reduction

    """

    # Copy data.
    data_raw = data_persons_properties_raw.copy(deep=True)
    # Extract and collect all unique values of facilities and batches across
    # all persons.
    values = collect_persons_unique_values(
        values=data_raw[variable].dropna().to_list(),
        delimiter=",",
    )
    # Create two-dimensional matrices to represent persons' categories.
    # Persons have values of zero (False) or one (True) to indicate whether
    # they had a sample in each category.
    data_matrix = expand_persons_categories_matrix(
        category=variable,
        values=values,
        delimiter=",",
        data=data_raw,
    )
    # Reduce dimensionality.
    # Count of principal components to keep.
    if len(values) > 20:
        components = 20
    else:
        components = len(values)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Variable: " + str(variable))
        print("Count of unique values: " + str(len(values)))
        utility.print_terminal_partition(level=2)
    # Principal components.
    report_dimension = utility.calculate_principal_components(
        data=data_matrix,
        components=components,
        report=report,
    )
    # Return information.
    return report_dimension


def insert_components(
    prefix=None,
    data_components=None,
    data_collection=None,
):
    """
    Inserts principal components of a property into a data frame.

    arguments:
        prefix (str): name to use as a prefix for components
        data_components (object): Pandas data frame of a property's principal
            components across observations
        data_collection (object): Pandas data frame in which to insert the
            property's principal components

    raises:

    returns:
        (object): Pandas data frame including property's principal components

    """

    # Copy data.
    data_components = data_components.copy(deep=True)
    data_collection = data_collection.copy(deep=True)
    # Change names of columns for property's principal components.
    components = data_components.shape[1]
    translations = dict()
    for count in range(components):
        original = ("component_" + str(count + 1))
        novel = (prefix + "_" + str(count + 1))
        translations[original] = novel
    data_components.rename(
        columns=translations,
        inplace=True,
    )
    # Insert components.
    data_combination = data_collection.join(
        data_components,
        how="left",
    )
    # Return information.
    return data_combination


def organize_persons_category_dimensionality(
    variables=None,
    data_persons_properties_raw=None,
    report=None,
):
    """
    Organizes information about persons.

    arguments:
        variables (list<str>): names of categorical variables
        data_persons_properties_raw (object): Pandas data frame of persons'
            properties
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy data.
    data_raw = data_persons_properties_raw.copy(deep=True)
    data_collection = data_persons_properties_raw.copy(deep=True)
    # Organize dimensionality reduction on categorical variables.
    components_data = dict()
    variances_data = dict()
    for variable in variables:
        # Reduce dimensionality of categorical variable's values.
        report_dimension = reduce_dimensionality_categorical_variable(
            variable=variable,
            data_persons_properties_raw=data_raw,
            report=report,
        )
        components_data[variable] = (
            report_dimension["data_observations_components"]
        )
        variances_data[variable] = (
            report_dimension["data_components_variances"]
        )
        # Introduce variable's principal components to persons' properties.
        data_collection = insert_components(
            prefix=variable,
            data_components=report_dimension["data_observations_components"],
            data_collection=data_collection,
        )
    # Compile information.
    bin = dict()
    bin["data_persons_properties_dimension"] = data_collection
    bin["components_data"] = components_data
    bin["variances_data"] = variances_data
    # Return information.
    return bin


def standardize_scale_variables(
    variables=None,
    data_persons_properties=None,
):
    """
    Standardizes variables' values to z-score scale.

    arguments:
        variables (list<str>): names of numerical variables
        data_persons_properties_raw (object): Pandas data frame of persons'
            properties

    raises:

    returns:
        (object): Pandas data frame of information about persons

    """

    # Copy data.
    data = data_persons_properties.copy(deep=True)
    # Iterate on variables.
    for variable in variables:
        variable_scale = str(variable + ("_scale"))
        data[variable] = data[variable].astype(numpy.float32)
        data[variable_scale] = data[variable].pipe(
            lambda series: scipy.stats.zscore(
                series.to_numpy(),
                axis=0,
                ddof=1, # Sample standard deviation.
                nan_policy="omit",
            ),
        )
        pass
    # Return information.
    return data


def organize_persons_properties(
    persons=None,
    variables=None,
    data_persons_properties_raw=None,
    report=None,
):
    """
    Organizes information about persons.

    arguments:
        persons (list<str>): identifiers of persons
        variables (dict<list<str>>): names of independent variables for
            regression
        data_persons_properties_raw (object): Pandas data frame of persons'
            properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about persons

    """

    # Copy data.
    data_raw = data_persons_properties_raw.copy(deep=True)
    # Select data for persons.
    data_raw_persons = data_raw.loc[
        data_raw.index.isin(persons), :
    ]
    # Organize dimensionality reduction on categorical variables.
    bin_category = organize_persons_category_dimensionality(
        variables=variables["category"],
        data_persons_properties_raw=data_raw_persons,
        report=report,
    )
    # Standardize the scales of numerical variables for regression.
    # These variables are on ordinal or ratio (continuous) scales.
    data_persons_properties_scale = standardize_scale_variables(
        variables=variables["scale"],
        data_persons_properties=bin_category[
            "data_persons_properties_dimension"
        ],
    )
    # Compile information.
    information = dict()
    information["data_persons_properties"] = data_persons_properties_scale
    information["components_data"] = bin_category["components_data"]
    information["category_variances_data"] = bin_category["variances_data"]
    # Return information.
    return information


def organize_persons_properties_sets(
    persons_selection=None,
    data_samples_tissues_persons=None,
    data_persons_properties=None,
    report=None,
):
    """
    Extracts identifiers of persons with valid genotypes.

    arguments:
        persons_selection (list<str>): identifiers of persons from selection
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues across all samples
        data_persons_properties (object): Pandas data frame of persons'
            properties
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information for charts about persons' properties

    """

    # Copy data.
    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_persons = data_persons_properties.copy(deep=True)
    # Organize data.
    bin = dict()
    bin["gtex"] = utility.collect_unique_elements(
        elements_original=data_samples["person"].to_list()
    )
    bin["gtex_genotype"] = extract_persons_genotype(
        persons=bin["gtex"],
        data_samples_tissues_persons=data_samples,
    )
    bin["selection"] = persons_selection
    bin["selection_genotype"] = extract_persons_genotype(
        persons=persons_selection,
        data_samples_tissues_persons=data_samples,
    )
    bin["ventilation"] = data_persons_properties.loc[
        data_persons_properties["ventilation"] == True, :
    ].index.to_list()
    bin["non_ventilation"] = data_persons_properties.loc[
        data_persons_properties["ventilation"] == False, :
    ].index.to_list()
    bin["respiration"] = data_persons_properties.loc[
        data_persons_properties["respiration"] == True, :
    ].index.to_list()
    bin["inflammation"] = data_persons_properties.loc[
        data_persons_properties["inflammation"] == True, :
    ].index.to_list()
    bin["infection"] = data_persons_properties.loc[
        data_persons_properties["infection"] == True, :
    ].index.to_list()
    bin["steroid"] = data_persons_properties.loc[
        data_persons_properties["steroid"] == True, :
    ].index.to_list()
    bin["refrigeration"] = data_persons_properties.loc[
        data_persons_properties["refrigeration"] == True, :
    ].index.to_list()

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Organization of information for charts about persons' properties."
        )
        utility.print_terminal_partition(level=2)
        print(
            "Count of original persons in GTEx cohort: " +
            str(len(bin["gtex"]))
        )
        print(
            "Count of persons from selection by samples and tissues: " +
            str(len(bin["selection"]))
        )
        print(
            "Count of original persons with valid genotypes: " +
            str(len(bin["gtex_genotype"]))
        )
        print(
            "Count of selection persons with valid genotypes: " +
            str(len(bin["selection_genotype"]))
        )
        print(
            "Count of persons on ventilation: " +
            str(len(bin["ventilation"]))
        )
        print(
            "Count of persons with respiratory difficulties: " +
            str(len(bin["respiration"]))
        )
        print(
            "Count of persons on inflammation: " +
            str(len(bin["inflammation"]))
        )
        print(
            "Count of persons on infection: " +
            str(len(bin["infection"]))
        )
        print(
            "Count of persons on steroids: " +
            str(len(bin["steroid"]))
        )
        print(
            "Count of persons on refrigeration after death: " +
            str(len(bin["refrigeration"]))
        )
        utility.print_terminal_partition(level=2)
    # Return information.
    return bin



def organize_contingency_table_chi(
    persons_selection=None,
    variables_contingency=None,
    data_persons_properties=None,
    report=None,
):
    """
    Extracts identifiers of persons with valid genotypes.

    arguments:
        persons_selection (list<str>): identifiers of persons from selection
        variables_contingency (list<str>): names of variables for contingency
        data_persons_properties (object): Pandas data frame of persons'
            properties
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information for charts about persons' properties

    """

    # TODO: use pandas.crosstab()

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    # Organize data.
    # Select data for persons.
    data_persons = data_persons_properties.loc[
        data_persons_properties.index.isin(persons_selection), :
    ]
    data_properties = data_persons.loc[
        :, data_persons.columns.isin(variables_contingency)
    ]
    # Replace missing values with zero.
    data_properties.fillna(
        value=False,
        #axis="columns",
        inplace=True,
    )
    # Contingency table.
    data_contingency = pandas.crosstab(
        data_properties[variables_contingency[0]],
        data_properties[variables_contingency[1]],
        rownames=[variables_contingency[0]],
        colnames=[variables_contingency[1]],
    )
    (chi2, probability, freedom, expectation) = scipy.stats.chi2_contingency(
        data_contingency.to_numpy(),
        correction=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Contingency table and Chi2 test for independence."
        )
        print(
            str(variables_contingency[0]) + " versus " +
            str(variables_contingency[1])
        )
        utility.print_terminal_partition(level=2)
        print(data_contingency)
        print(data_contingency.to_numpy())
        utility.print_terminal_partition(level=3)
        print("chi2: " + str(chi2))
        print("probability: " + str(probability))

        if False:
            # CMV and EBV.
            persons_male = data_persons_properties.loc[
                data_persons_properties["sex_text"] == "male", :
            ].index.to_list()
            print("persons male: " + str(len(persons_male)))
            persons_female = data_persons_properties.loc[
                data_persons_properties["sex_text"] == "female", :
            ].index.to_list()
            print("persons female: " + str(len(persons_female)))
            persons_cmv_ebv = data_persons_properties.loc[
                data_persons_properties["cmv_ebv"] == True, :
            ].index.to_list()
            print("persons cmv_ebv: " + str(len(persons_cmv_ebv)))

    pass



###################################################
# TODO: Still need to organize code for heritability covariates and QTL covariates....

# Organize covariates for heritability analysis.


def organize_heritability_variables(
    variables=None,
    data_persons_properties=None,
):
    """
    Organizes information about families and persons for heritability analysis.

    As the data from GTEx do not include families, generate distinct families
    for each person for use in heritability analysis in GCTA.

    arguments:
        variables (list<str>): names of variables for regression
        data_persons_properties (object): Pandas data frame of persons and
            their properties

    raises:

    returns:
        (dict): information for heritability analysis
    """

    # Introduce empty identifiers for families for compatibility with Plink and
    # GCTA.
    # Organize data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    data_persons_properties["family"] = 0
    data_persons_properties.reset_index(
        level=None,
        inplace=True
    )
    # Extract families and persons.
    data_copy = data_persons_properties.copy(deep=True)
    data_families_persons = data_copy.loc[
        :, data_copy.columns.isin(["person", "family"])
    ]
    data_families_persons = data_families_persons[[
        "family",
        "person",
    ]]

    # Do not use categorical covariates in heritability analysis.
    # These would use GCTA's "--covar" flag.
    # Instead, use the quantitative adaptations.

    # Extract quantitative covariates.
    columns = copy.deepcopy(variables)
    columns.insert(0, "person")
    columns.insert(0, "family")
    data_copy = data_persons_properties.copy(deep=True)
    data_persons_variables = data_copy.loc[
        :, data_copy.columns.isin(columns)
    ]
    data_persons_variables = data_persons_variables[[*columns]]

    # Compile information.
    information = {
        "data_families_persons": data_families_persons,
        "data_persons_variables": data_persons_variables,
    }
    # Return information.
    return information


# Organize covariates for quantitative trait loci (QTL) analysis.


def organize_quantitative_trait_loci_variables(
    variables=None,
    data_persons_properties=None,
):
    """
    Organizes information about persons for quantitative trait loci (QTL)
    analysis.

    arguments:
        variables (list<str>): names of variables for regression
        data_persons_properties (object): Pandas data frame of persons and
            their properties

    raises:

    returns:
        (object): Pandas data frame of persons' variables for quantitative
            trait loci analysis
    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)

    # Extract quantitative covariates.
    columns = copy.deepcopy(variables)
    data_persons_variables = data_persons_properties.loc[
        :, data_persons_properties.columns.isin(columns)
    ]

    # Transpose columns to rows and rows to columns.
    # Organize variables across rows and persons across columns.
    data_transposition = data_persons_variables.transpose(copy=True)
    data_transposition.reset_index(
        level=None,
        inplace=True
    )
    data_transposition.rename_axis(
        index="",
        axis="index",
        copy=False,
        inplace=True
    )
    data_transposition.rename_axis(
        columns="",
        axis="columns",
        copy=False,
        inplace=True
    )
    data_transposition.rename(
        columns={
            "index": "id",
        },
        inplace=True,
    )

    # Return information.
    return data_transposition


###############################################################################33

# TODO: within this function...
# TODO: organize persons_selection and persons_genotype
# TODO: calculate covariate principal components using only the relevant persons
# TODO: 1. regression (all 631 persons with imputed genotypes)
# TODO: 2. heritability and trait with only the 555 or whatever persons with valid genotypes

# TODO: change "batch_isolation" to "batch_extraction"
#

def extract_organize_persons_properties(
    stringency=None,
    dock=None,
):
    """
    Extracts and organizes information about samples and genes for further
    analysis.
    # Organize ventilation and refrigeration variables.
    # If a person did not receive ventilation or refrigeration, then change
    # respective duration variables to 0.0.

    # TODO: translate season to climate...
    # TODO: do any other variable modifications...

    Selection of persons
    1. for modality analysis: persons_selection
    - persons with samples for at least 10 of 20 sex-neutral tissues
    2. for regression analysis: persons_selection
    - impute missing genotypes to maximize observations in regression
    - impute missing genotypes by mean across only persons_selection
    - principal components on genotypes only consider persons_selection anyway
    3. for heritability analysis: persons_selection_genotype
    - only include persons with valid genotypic (SNP) data
    4. for quantitative trait loci (QTL) analysis: persons_selection_genotype
    - only include persons with valid genotypic (SNP) data

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        dock (str): path to root or dock directory for source and product
            directories and files

        old arguments...

        persons_selection (list<str>): identifiers of persons from selection of
            samples, tissues, and persons for further analysis
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_samples_tissues_persons_selection (object): Pandas data frame of
            persons and tissues for all samples
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
    """

    # Remove previous files to avoid version or batch confusion.
    path_selection = os.path.join(dock, "selection", str(stringency))
    path_persons_properties = os.path.join(
        path_selection, "persons_properties"
    )
    utility.remove_directory(path=path_persons_properties)

    # Read source information from file.
    source = read_source_persons_properties(
        stringency=stringency,
        dock=dock,
    )

    ##########
    ##########
    ##########
    # Enhance person's properties within sample data before collapsing to
    # persons.

    # Extract identifiers of persons with valid genotypes before imputation.
    # These persons also meet the selection criteria from earlier.
    # This function takes the original data for all samples in GTEx.
    persons_selection_genotype = extract_persons_selection_genotype(
        persons_selection=source["persons_selection"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        report=True,
    )
    # Impute genotypes for selection persons for regression analysis.
    data_samples_imputation = impute_persons_genotypes(
        data_samples_tissues_persons_selection=(
            source["data_samples_tissues_persons_selection"]
        ),
    )
    # Define variables.
    variables = define_variables()
    # Organize persons' properties across samples.
    data_samples_organization = organize_samples_properties(
        variables=variables,
        data_samples_tissues_persons_selection=data_samples_imputation,
    )

    ##########
    ##########
    ##########
    # Collapse sample data to persons' properties.

    #####################
    # TODO: I need to accommodate different persons for different analyses...
    # 1. persons_selection for basic regression
    # 2. persons_genotype for heritability and QTL...
    # TODO: Plan
    # organize "extract_persons_properties" and "organize_persons_properties"
    # within a larger function that also accepts "persons_selection" or "persons_selection_genotype"
    # to subset sample data before organizing properties...

    # Extract information about persons' properties.
    data_persons_properties_raw = extract_persons_properties(
        data_samples_tissues_persons=data_samples_organization,
        report=True,
    )
    # Organize information about sets of persons.
    persons_sets = organize_persons_properties_sets(
        persons_selection=source["persons_selection"],
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_persons_properties=data_persons_properties_raw,
        report=True,
    )

    # Organize persons' properties.
    # Expand categorical variables and reduce dimensionality.
    bin_organization = organize_persons_properties(
        persons=persons_sets["ventilation"],
        variables=variables,
        data_persons_properties_raw=data_persons_properties_raw,
        report=True,
    )
    utility.print_terminal_partition(level=2)
    print("data_persons_properties")
    print(bin_organization["data_persons_properties"])

    bin_regression = dict()
    bin_regression["persons"] = source["persons_selection"]
    bin_regression["data_persons_properties"] = bin_organization["data_persons_properties"]

    bin_heritability = dict()
    bin_trait = dict()

    # Contingency tables and chi2.
    # Test on persons_selection and on persons_ventilation.
    organize_contingency_table_chi(
        persons_selection=source["persons_selection"],
        variables_contingency=["ventilation", "sex_text"],
        data_persons_properties=data_persons_properties_raw,
        report=True,
    )
    organize_contingency_table_chi(
        persons_selection=source["persons_selection"],
        variables_contingency=["ventilation", "mononucleosis"],
        data_persons_properties=data_persons_properties_raw,
        report=True,
    )
    organize_contingency_table_chi(
        persons_selection=persons_sets["ventilation"],
        variables_contingency=["sex_text", "mononucleosis"],
        data_persons_properties=data_persons_properties_raw,
        report=True,
    )


    # Compile information.
    information = dict()
    information["persons_selection_genotype"] = persons_selection_genotype
    information["data_samples_tissues_persons_imputation"] = (
        data_samples_imputation
    )
    information["category_variances_data"] = (
        bin_organization["category_variances_data"]
    )
    information["persons_sets"] = persons_sets
    information["bin_regression"] = bin_regression
    # Write product information to file.
    write_product_persons_properties(
        stringency=stringency,
        dock=dock,
        information=information
    )
    pass






    if False:
        # Copy data.
        data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
        data_samples_tissues_persons_selection = (
            data_samples_tissues_persons_selection.copy(deep=True)
        )
        data_gene_signal = data_gene_signal.copy(deep=True)

        if False:
            # Transpose data structure.
            # Organize genes across columns and samples across rows.
            data_transposition = data_gene_signal.transpose(copy=True)
            # Associate samples to persons and tissues.
            data_gene_signal_factor = assembly.associate_samples_persons_tissues(
                data_samples_tissues_persons=data_samples_tissues_persons,
                data_gene_sample=data_transposition,
            )

        # Select persons.


        # TODO: here is a critical point... for imputation, I want all the persons...

        # Split up this procedure by persons considered...


        ##########
        # Regression
        ##########
        #bin_regression =

        ##########
        # Heritability
        ##########
        #bin_heritability =

        ##########
        # Quantitative trait loci (QTL)
        ##########
        #bin_trait =


        # Organize information for heritability analysis.
        simple = organize_heritability_variables(
            variables=variables["heritability_simple"],
            data_persons_properties=bin["data_persons_properties"],
        )
        complex = organize_heritability_variables(
            variables=variables["heritability_complex"],
            data_persons_properties=bin["data_persons_properties"],
        )

        # Organize information for quantitative trait loci (QTL) analysis.
        data_persons_variables_trait = organize_quantitative_trait_loci_variables(
            variables=variables["trait"],
            data_persons_properties=bin["data_persons_properties"],
        )

        # Compile information.
        information = {
            "data_persons_properties_raw": data_persons_properties_raw,
            "data_persons_properties": bin["data_persons_properties"],

            "data_tissues_variance": bin["data_tissues_variance"],
            "data_facilities_variance": bin["data_facilities_variance"],
            "data_batches_extraction_variance": bin[
                "data_batches_extraction_variance"
            ],
            "data_batches_sequence_variance": bin[
                "data_batches_sequence_variance"
            ],

            "data_families_persons": simple["data_families_persons"],
            "data_persons_variables_simple": simple["data_persons_variables"],
            "data_persons_variables_complex": complex["data_persons_variables"],

            "data_persons_variables_trait": data_persons_variables_trait,

        }
        # Return information.
        return information

    pass


##########
##########
##########
# Summary of persons' properties.


def prepare_persons_properties_summary(
    persons=None,
    data_samples_tissues_persons=None,
    data_persons_properties=None,
):
    """
    Prepares summaries of persons' properties.

    arguments:
        persons (list<str>): identifiers of persons
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_persons_properties (object): Pandas data frame of persons'
            properties

    raises:

    returns:
        (dict): information
    """

    # Copy data.
    data_samples_tissues_persons = data_samples_tissues_persons.copy(deep=True)
    data_persons_properties = data_persons_properties.copy(deep=True)

    # Count persons in groups by multiple factors.
    utility.print_terminal_partition(level=1)
    print("Summary: counting persons in groups by factors.")
    data_persons_sex_decade_counts = (
        utility.count_data_factors_groups_elements(
            factors=["sex", "decade"],
            element="person",
            count="counts",
            data=data_persons_properties,
    ))
    utility.print_terminal_partition(level=2)
    print("Sex and decade...")
    utility.print_terminal_partition(level=3)
    print(data_persons_sex_decade_counts)
    utility.print_terminal_partition(level=3)
    sum_persons = data_persons_sex_decade_counts["counts"].sum()
    print("Sum persons: " + str(sum_persons))

    data_persons_sex_hardiness_counts = (
        utility.count_data_factors_groups_elements(
            factors=["sex", "hardiness"],
            element="person",
            count="counts",
            data=data_persons_properties,
    ))
    utility.print_terminal_partition(level=2)
    print("Sex and hardiness...")
    utility.print_terminal_partition(level=3)
    print(data_persons_sex_hardiness_counts)

    data_persons_hardiness_decade_counts = (
        utility.count_data_factors_groups_elements(
            factors=["hardiness", "decade"],
            element="person",
            count="counts",
            data=data_persons_properties,
    ))
    utility.print_terminal_partition(level=2)
    print("Hardiness and decade...")
    utility.print_terminal_partition(level=3)
    print(data_persons_hardiness_decade_counts)

    # Count tissues or persons in groups by single factor.
    # Count tissues per person and persons per tissue.
    data_tissues_per_person = utility.count_data_factors_groups_elements(
        factors=["person"],
        element="tissue_major",
        count="counts",
        data=data_samples_tissues_persons,
    )
    utility.print_terminal_partition(level=2)
    print("Tissues per person...")
    utility.print_terminal_partition(level=3)
    print(data_tissues_per_person)
    utility.print_terminal_partition(level=3)
    mean_tissues = data_tissues_per_person["counts"].mean()
    print("Mean tissues per person: " + str(mean_tissues))

    data_persons_per_tissue = utility.count_data_factors_groups_elements(
        factors=["tissue_major"],
        element="person",
        count="counts",
        data=data_samples_tissues_persons,
    )
    utility.print_terminal_partition(level=2)
    print("Persons per tissue...")
    utility.print_terminal_partition(level=3)
    print(data_persons_per_tissue)

    # Compile information.
    information = {

        # Sets of persons
        "data_persons_sex_decade_counts": data_persons_sex_decade_counts,
        "data_persons_sex_hardiness_counts": data_persons_sex_hardiness_counts,
        "data_persons_hardiness_decade_counts": data_persons_hardiness_decade_counts,

        # Sample coverage
        "data_tissues_per_person": data_tissues_per_person,
        "data_persons_per_tissue": data_persons_per_tissue,
    }
    # Return information.
    return information


##########
# Product.


def write_product_gene_annotation(
    stringency=None,
    dock=None,
    information=None
):
    """
    Writes product information to file.

    arguments:
        stringency (str): category, loose or tight, of selection criteria
        dock (str): path to root or dock directory for source and product
            directories and files.
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    path_gene_annotation = os.path.join(path_selection, "gene_annotation")
    utility.create_directories(path_gene_annotation)

    path_gene_annotation_gtex = os.path.join(
        path_gene_annotation, "data_gene_annotation_gtex.pickle"
    )
    path_gene_annotation_gtex_text = os.path.join(
        path_gene_annotation, "data_gene_annotation_gtex.tsv"
    )
    path_gene_annotation_gencode = os.path.join(
        path_gene_annotation, "data_gene_annotation_gencode.pickle"
    )
    path_gene_annotation_gencode_text = os.path.join(
        path_gene_annotation, "data_gene_annotation_gencode.tsv"
    )

    path_genes_gtex_annotation = os.path.join(
        path_gene_annotation, "genes_gtex_annotation.pickle"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_gene_annotation_gtex"],
        path_gene_annotation_gtex
    )
    pandas.to_pickle(
        information["data_gene_annotation_gencode"],
        path_gene_annotation_gencode
    )
    information["data_gene_annotation_gtex"].to_csv(
        path_or_buf=path_gene_annotation_gtex_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_gene_annotation_gencode"].to_csv(
        path_or_buf=path_gene_annotation_gencode_text,
        sep="\t",
        header=True,
        index=True,
    )

    with open(path_genes_gtex_annotation, "wb") as file_product:
        pickle.dump(information["genes_gtex_annotation"], file_product)

    pass


def write_product_samples_genes_signals_charts(
    dock=None,
    stringency=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    path_samples_genes_signals = os.path.join(
        path_selection, "samples_genes_signals"
    )
    utility.create_directories(path_samples_genes_signals)

    path_signals_initial = os.path.join(
        path_selection, "signals_initial.pickle"
    )
    path_signals_final = os.path.join(
        path_selection, "signals_final.pickle"
    )

    # Write information to file.
    with open(path_signals_initial, "wb") as file_product:
        pickle.dump(information["signals_initial"], file_product)
    with open(path_signals_final, "wb") as file_product:
        pickle.dump(information["signals_final"], file_product)

    pass


def write_product_samples_genes_signals(
    dock=None,
    stringency=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        information (object): information to write to file.

    raises:

    returns:

    """

    # Priority information for further analysis.

    # Information for plotting charts.
    write_product_samples_genes_signals_charts(
        dock=dock,
        stringency=stringency,
        information=information,
    )

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    path_samples_genes_signals = os.path.join(
        path_selection, "samples_genes_signals"
    )
    utility.create_directories(path_samples_genes_signals)

    path_gene_signal = os.path.join(
        path_samples_genes_signals, "data_gene_signal.pickle"
    )
    path_samples_tissues_persons = os.path.join(
        path_samples_genes_signals, "data_samples_tissues_persons.pickle"
    )
    path_samples_tissues_persons_text = os.path.join(
        path_samples_genes_signals, "data_samples_tissues_persons.tsv"
    )

    path_genes = os.path.join(
        path_samples_genes_signals, "genes.pickle"
    )
    path_genes_text = os.path.join(
        path_samples_genes_signals, "genes.txt"
    )

    path_samples = os.path.join(
        path_samples_genes_signals, "samples.pickle"
    )
    path_persons = os.path.join(
        path_samples_genes_signals, "persons.pickle"
    )
    path_persons_text = os.path.join(
        path_samples_genes_signals, "persons.txt"
    )
    path_tissues = os.path.join(
        path_samples_genes_signals, "tissues.pickle"
    )

    # Write information to file.

    pandas.to_pickle(
        information["data_gene_signal"],
        path_gene_signal
    )
    pandas.to_pickle(
        information["data_samples_tissues_persons"],
        path_samples_tissues_persons
    )
    information["data_samples_tissues_persons"].to_csv(
        path_or_buf=path_samples_tissues_persons_text,
        sep="\t",
        header=True,
        index=True,
    )
    with open(path_genes, "wb") as file_product:
        pickle.dump(information["genes"], file_product)
    utility.write_file_text_list(
        elements=information["genes"],
        delimiter="\n",
        path_file=path_genes_text
    )
    with open(path_samples, "wb") as file_product:
        pickle.dump(information["samples"], file_product)
    with open(path_persons, "wb") as file_product:
        pickle.dump(information["persons"], file_product)
    utility.write_file_text_list(
        elements=information["persons"],
        delimiter="\n",
        path_file=path_persons_text
    )
    with open(path_tissues, "wb") as file_product:
        pickle.dump(information["tissues"], file_product)

    pass


def write_product_persons_properties_charts(
    dock=None,
    stringency=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    path_persons_properties = os.path.join(
        path_selection, "persons_properties"
    )
    path_charts = os.path.join(
        path_persons_properties, "charts"
    )
    utility.create_directories(path_charts)

    path_persons_sets = os.path.join(
        path_charts, "persons_sets.pickle"
    )

    # Write information to file.
    with open(path_persons_sets, "wb") as file_product:
        pickle.dump(information["persons_sets"], file_product)

    pass


def write_product_persons_properties_analysis_regression(
    path_parent=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        path_parent (str): path to parent directory
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_directory = os.path.join(
        path_parent, "regression"
    )
    utility.create_directories(path_directory)

    path_data_persons_properties = os.path.join(
        path_directory, "data_persons_properties.pickle"
    )
    information["data_persons_properties"].to_pickle(
        path=path_data_persons_properties
    )

    path_data_persons_properties_text = os.path.join(
        path_directory, "data_persons_properties.tsv"
    )
    information["data_persons_properties"].to_csv(
        path_or_buf=path_data_persons_properties_text,
        sep="\t",
        header=True,
        index=True,
    )

    pass


def write_product_persons_properties(
    dock=None,
    stringency=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    path_persons_properties = os.path.join(
        path_selection, "persons_properties"
    )
    utility.create_directories(path_persons_properties)

    # Information for plotting charts.
    write_product_persons_properties_charts(
        stringency=stringency,
        dock=dock,
        information=information
    )
    # Information for subsequent analyses.
    write_product_persons_properties_analysis_regression(
        path_parent=path_persons_properties,
        information=information["bin_regression"]
    )


    if False:
        # Specify directories and files.
        path_selection = os.path.join(dock, "selection", str(stringency))
        path_persons_properties = os.path.join(
            path_selection, "persons_properties"
        )
        path_charts = os.path.join(
            path_persons_properties, "charts"
        )
        utility.create_directories(path_charts)

        path_persons_sets = os.path.join(
            path_charts, "persons_sets.pickle"
        )

        # Write information to file.
        with open(path_persons_sets, "wb") as file_product:
            pickle.dump(information["persons_sets"], file_product)

    pass






def write_product_heritability(
    dock=None,
    stringency=None,
    model=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        model (str): category, simple or complex, of regression model
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    utility.create_directories(path_selection)

    # Variables for heritability analysis.
    path_heritability = os.path.join(
        path_selection, "heritability", str(model)
    )
    utility.create_directories(path_heritability)
    path_families_persons = os.path.join(
        path_heritability, "data_families_persons.pickle"
    )
    path_families_persons_text = os.path.join(
        path_heritability, "families_persons.tsv"
    )
    path_persons_variables = os.path.join(
        path_heritability, "persons_variables.tsv"
    )

    # Write information to file.
    if model == "simple":
        pandas.to_pickle(
            information["data_families_persons"],
            path_families_persons
        )
        information["data_families_persons"].to_csv(
            path_or_buf=path_families_persons_text,
            sep="\t",
            na_rep="NA",
            header=False,
            index=False,
        )
        information["data_persons_variables_simple"].to_csv(
            path_or_buf=path_persons_variables,
            sep="\t",
            na_rep="NA",
            header=False,
            index=False,
        )
    elif model == "complex":
        pandas.to_pickle(
            information["data_families_persons"],
            path_families_persons
        )
        information["data_families_persons"].to_csv(
            path_or_buf=path_families_persons_text,
            sep="\t",
            na_rep="NA",
            header=False,
            index=False,
        )
        information["data_persons_variables_complex"].to_csv(
            path_or_buf=path_persons_variables,
            sep="\t",
            na_rep="NA",
            header=False,
            index=False,
        )
    pass


def write_product_trait(
    dock=None,
    stringency=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    utility.create_directories(path_selection)

    # Variables for quantitative trait loci (QTL) analysis.
    path_trait = os.path.join(
        path_selection, "trait"
    )
    utility.create_directories(path_trait)

    path_data_persons_variables = os.path.join(
        path_trait, "data_persons_variables.tsv"
    )

    # Write information to file.
    information["data_persons_variables_trait"].to_csv(
        path_or_buf=path_data_persons_variables,
        sep="\t",
        na_rep="NA",
        header=True,
        index=False,
    )
    pass


def write_product_charts(
    dock=None,
    stringency=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    utility.create_directories(path_selection)

    path_data_tissues_per_person = os.path.join(
        path_selection, "data_tissues_per_person.pickle"
    )
    path_data_persons_per_tissue = os.path.join(
        path_selection, "data_persons_per_tissue.pickle"
    )
    path_data_persons_sex_decade_counts = os.path.join(
        path_selection, "data_persons_sex_decade_counts.pickle"
    )
    path_data_persons_sex_hardiness_counts = os.path.join(
        path_selection, "data_persons_sex_hardiness_counts.pickle"
    )
    path_data_persons_hardiness_decade_counts = os.path.join(
        path_selection, "data_persons_hardiness_decade_counts.pickle"
    )

    path_data_tissues_variance = os.path.join(
        path_selection, "data_tissues_variance.pickle"
    )
    path_data_facilities_variance = os.path.join(
        path_selection, "data_facilities_variance.pickle"
    )
    path_data_batches_extraction_variance = os.path.join(
        path_selection, "data_batches_extraction_variance.pickle"
    )
    path_data_batches_sequence_variance = os.path.join(
        path_selection, "data_batches_sequence_variance.pickle"
    )

    # Write information to file.

    pandas.to_pickle(
        information["data_tissues_per_person"],
        path_data_tissues_per_person
    )
    pandas.to_pickle(
        information["data_persons_per_tissue"],
        path_data_persons_per_tissue
    )
    pandas.to_pickle(
        information["data_persons_sex_decade_counts"],
        path_data_persons_sex_decade_counts
    )
    pandas.to_pickle(
        information["data_persons_sex_hardiness_counts"],
        path_data_persons_sex_hardiness_counts
    )
    pandas.to_pickle(
        information["data_persons_hardiness_decade_counts"],
        path_data_persons_hardiness_decade_counts
    )

    pandas.to_pickle(
        information["data_tissues_variance"],
        path_data_tissues_variance,
    )
    pandas.to_pickle(
        information["data_facilities_variance"],
        path_data_facilities_variance,
    )
    pandas.to_pickle(
        information["data_batches_extraction_variance"],
        path_data_batches_extraction_variance,
    )
    pandas.to_pickle(
        information["data_batches_sequence_variance"],
        path_data_batches_sequence_variance,
    )

    pass


def write_product(
    dock=None,
    stringency=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        stringency (str): category, loose or tight, of selection criteria
        information (object): information to write to file.

    raises:

    returns:

    """

    # Priority information for further analysis.

    # Information for heritability analysis.
    write_product_heritability(
        dock=dock,
        stringency=stringency,
        model="simple",
        information=information,
    )
    write_product_heritability(
        dock=dock,
        stringency=stringency,
        model="complex",
        information=information,
    )

    # Information for quantitative trait loci (QTL) analysis.
    write_product_trait(
        dock=dock,
        stringency=stringency,
        information=information,
    )

    # Information for plotting charts.
    write_product_charts(
        dock=dock,
        stringency=stringency,
        information=information,
    )

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", str(stringency))
    utility.create_directories(path_selection)

    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    #path_gene_signal_factor = os.path.join(
    #    path_selection, "data_gene_signal_factor.pickle"
    #)
    path_genes_selection = os.path.join(
        path_selection, "genes_selection.pickle"
    )
    path_genes_selection_text = os.path.join(
        path_selection, "genes_selection.txt"
    )

    path_samples_tissues_persons = os.path.join(
        path_selection, "data_samples_tissues_persons.pickle"
    )
    path_samples_tissues_persons_text = os.path.join(
        path_selection, "data_samples_tissues_persons.tsv"
    )
    path_samples_selection = os.path.join(
        path_selection, "samples_selection.pickle"
    )
    path_persons_selection = os.path.join(
        path_selection, "persons_selection.pickle"
    )
    path_persons_selection_text = os.path.join(
        path_selection, "persons_selection.txt"
    )
    path_tissues = os.path.join(
        path_selection, "tissues_selection.pickle"
    )

    path_persons_properties_raw = os.path.join(
        path_selection, "data_persons_properties_raw.pickle"
    )
    path_persons_properties_raw_text = os.path.join(
        path_selection, "data_persons_properties_raw.tsv"
    )
    path_persons_properties = os.path.join(
        path_selection, "data_persons_properties.pickle"
    )
    path_persons_properties_text = os.path.join(
        path_selection, "data_persons_properties.tsv"
    )

    # Write information to file.

    pandas.to_pickle(
        information["data_gene_signal"],
        path_gene_signal
    )
    #pandas.to_pickle(
    #    information["data_gene_signal_factor"],
    #    path_gene_signal_factor
    #)
    with open(path_genes_selection, "wb") as file_product:
        pickle.dump(information["genes_selection"], file_product)
    utility.write_file_text_list(
        elements=information["genes_selection"],
        delimiter="\n",
        path_file=path_genes_selection_text
    )
    with open(path_samples_gtex, "wb") as file_product:
        pickle.dump(information["samples_gtex"], file_product)
    with open(path_samples_selection, "wb") as file_product:
        pickle.dump(information["samples_selection"], file_product)
    with open(path_persons, "wb") as file_product:
        pickle.dump(information["persons_selection"], file_product)
    utility.write_file_text_list(
        elements=information["persons_selection"],
        delimiter="\n",
        path_file=path_persons_selection_text
    )
    with open(path_tissues, "wb") as file_product:
        pickle.dump(information["tissues_selection"], file_product)

    pandas.to_pickle(
        information["data_samples_tissues_persons"],
        path_samples_tissues_persons
    )
    information["data_samples_tissues_persons"].to_csv(
        path_or_buf=path_samples_tissues_persons_text,
        sep="\t",
        header=True,
        index=True,
    )

    pandas.to_pickle(
        information["data_persons_properties_raw"],
        path_persons_properties_raw
    )
    information["data_persons_properties_raw"].to_csv(
        path_or_buf=path_persons_properties_raw_text,
        sep="\t",
        header=True,
        index=True,
    )
    pandas.to_pickle(
        information["data_persons_properties"],
        path_persons_properties
    )
    information["data_persons_properties"].to_csv(
        path_or_buf=path_persons_properties_text,
        sep="\t",
        header=True,
        index=True,
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

    # Enable automatic garbage collection to clear memory.
    gc.enable()

    ##################################################
    ##################################################
    ##################################################

    # Select genes' annotations.
    if False:
        select_organize_genes_annotations(
            stringency="tight",
            dock=dock
        )
        utility.print_terminal_partition(level=2)
        print("Completed selection of gene annotations.")
        # Pause procedure.
        #time.sleep(60.0)

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # Select samples, genes, and their signals.
    if False:
        select_organize_samples_genes_signals(
            stringency="tight",
            dock=dock
        )
        utility.print_terminal_partition(level=2)
        print("Completed selection of genes, samples, signals.")
        # Pause procedure.
        #time.sleep(60.0)

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # Organize persons and their properties for analyses.
    if True:
        extract_organize_persons_properties(
            stringency="tight",
            dock=dock,
        )

    if False:

        # Read source information from file.
        source = read_source(dock=dock)

        # Selection: tight.
        # Select data for pantissue aggregate analysis.
        # Criteria:
        # 1. select genes on autosomes that encode proteins
        # 2. select samples from whole tissues, not cell lines
        # 3. select samples from tissues with coverage of >=100 samples
        # 4. select samples from sex-neutral tissues
        # 5. select genes with signal (>= 0.1 TPM) in >=50% of samples
        # 6. select samples with signal (>= 0.1 TPM) in >=50% of genes
        bin_selection = select_organize_samples_genes_signals(
            data_gene_annotation_gtex=source["data_gene_annotation_gtex"],
            data_gene_annotation_gencode=source["data_gene_annotation_gencode"],
            data_samples_tissues_persons=source["data_samples_tissues_persons"],
            data_gene_signal=data_gene_signal,
            stringency="tight",
        )

        if False:
            # Organize information from selection of samples, genes, and signals.
            bin_organization = extract_organize_persons_properties(
                persons_selection=bin["persons_selection"],
                data_samples_tissues_persons=source["data_samples_tissues_persons"],
                data_samples_tissues_persons_selection=(
                    bin["data_samples_tissues_persons"]
                ),
                data_gene_signal=bin["data_gene_signal"],
            )
            # Prepare summary data about persons' properties.
            summary = prepare_persons_properties_summary(
                persons=bin["persons_selection"],
                data_samples_tissues_persons=bin["data_samples_tissues_persons"],
                data_persons_properties=organization["data_persons_properties"],
            )
            # Compile information.
            information = dict()
            information.update(bin_selection)
            information.update(bin_organization)
            information.update(summary)

            # Write product information to file.
            write_product(
                dock=dock,
                stringency="tight",
                information=information,
            )

    pass


if (__name__ == "__main__"):
    execute_procedure()
