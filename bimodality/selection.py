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

# Relevant

import numpy
import pandas

# Custom

import assembly
import measurement
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
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_gene_count = os.path.join(
        path_assembly, "data_gene_count.feather"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.feather"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(
        path_gene_annotation
    )
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_count = pandas.read_feather(
        path=path_gene_count
    )
    data_gene_signal = pandas.read_feather(
        path=path_gene_signal,
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_count": data_gene_count,
        "data_gene_signal": data_gene_signal,
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
    utility.print_terminal_partition(level=2)

    # Organize data.

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


# Summary.


def summarize_samples_genes(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Summarize selection of samples and genes.

    arguments:
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
    # Counts of samples and genes.
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
    print(data_adipose)

    pass


# Genes.


def select_gene_annotation(
    data_gene_annotation=None,
):
    """
    Selects genes that encode proteins.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' annotations

    """

    utility.print_terminal_partition(level=2)
    print(
        "Select annotations for genes that encode proteins and are on " +
        "autosomes."
    )

    # Select entries for genes.
    print(data_gene_annotation.shape)
    data_gene_annotation = (
        data_gene_annotation.loc[data_gene_annotation["feature"] == "gene", :]
    )
    data_gene_annotation.drop(
        labels="feature",
        axis="columns",
        inplace=True
    )
    print(data_gene_annotation.shape)

    # Select entries for genes on autosomes (non-sex chromosomes).
    # There are 2422 genes on the X chromosome, of which 850 encode proteins.
    # There are 567 genes on the Y chromosome, of which 66 encode proteins.
    # There are 37 genes on the mitochondrial chromosome, of which 13 encode
    # proteins.
    print(
        "count of genes in GENCODE annotations: " +
        str(data_gene_annotation.shape[0])
    )
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
    print(
        "count of GENCODE genes on nuclear (not mito) autosomes: " +
        str(data_gene_annotation.shape[0])
    )
    # Select entries for genes that encode proteins.
    data_gene_annotation = (
        data_gene_annotation.loc[
            data_gene_annotation["gene_type"] == "protein_coding", :
        ]
    )
    print(
        "count of GENCODE genes that encode proteins: " +
        str(data_gene_annotation.shape[0])
    )

    print(data_gene_annotation.iloc[0:10, 0:7])

    return data_gene_annotation


def select_genes(
    data_gene_annotation=None,
    data_gene_signal=None
):
    """
    Selects genes that encode proteins.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and patients

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients.

    """

    # Select protein-coding genes from signal data.
    # Signal matrix has genes on the index dimension.
    utility.print_terminal_partition(level=2)
    print(
        "Selection of signal genes that encode proteins."
    )
    # Describe original count of genes.
    print("signal total genes: " + str(data_gene_signal.shape[0]))
    genes_signal = data_gene_signal.index.to_list()
    print("signal total genes: " + str(len(genes_signal)))

    # Extract identifiers of protein-coding genes.
    genes_protein = data_gene_annotation.index.to_list()
    print(
        "count of GENCODE genes of type 'protein_coding': " +
        str(len(genes_protein))
    )

    # Filter gene signals.
    # Remove indices.
    data_gene_signal.reset_index(
        level=None,
        inplace=True
    )
    data_gene_signal_protein = data_gene_signal.loc[
        data_gene_signal["gene"].isin(genes_protein), :
    ]
    genes = data_gene_signal_protein.index.to_list()

    utility.print_terminal_partition(level=2)
    print(
        "signal genes that encode proteins and on autosomes: " +
        str(data_gene_signal_protein.shape[0])
    )
    print(
        "signal genes that encode proteins and on autosomes: " +
        str(len(genes))
    )
    print(
        "unique signal genes that encode proteins and on autosomes: " +
        str(len(utility.collect_unique_elements(elements_original=genes)))
    )

    utility.print_terminal_partition(level=2)

    # Organize data.
    data_gene_signal_protein = organize_data_axes_indices(
        data=data_gene_signal_protein
    )

    return data_gene_signal_protein


# Samples.


def select_samples_by_exclusion(
    data_samples_tissues_persons=None,
):
    """
    Selects samples by exclusion of specific tissues or persons.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (list<str>): identifiers of samples

    """

    # Define tissues for exclusion.
    # These tissues are extra-corporeal tissues, such as cell lines.
    # skin, fibroblast
    # blood, lymphocyte
    tissues_extracorp = []
    tissues_extracorp.append("fibroblast")
    tissues_extracorp.append("lymphocyte")
    data_samples_corp = data_samples_tissues_persons.loc[
        ~data_samples_tissues_persons["tissue_minor"].isin(tissues_extracorp),
    ]
    # These tissues have fewer than 100 total samples.
    tissues_coverage = []
    tissues_coverage.append("bladder")
    tissues_coverage.append("cervix")
    tissues_coverage.append("fallopius")
    tissues_coverage.append("kidney")
    data_samples_coverage = data_samples_corp.loc[
        ~data_samples_corp["tissue_major"].isin(tissues_coverage),
    ]
    # Extract identifiers of samples.
    samples = data_samples_coverage.index.to_list()
    # Return information.
    return samples


def select_samples_by_inclusion(
    data_samples_tissues_persons=None,
):
    """
    Selects samples by inclusion of specific tissues or persons.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (list<str>): identifiers of samples

    """

    # Define tissues for inclusion.
    # These tissues are not sex-specific.
    tissues_sex = [
        "adipose", # 797
        "adrenal", # 258
        "artery", # 786
        "blood", # 767
        "brain", # 382
        #"breast", # 459
        "colon", # 555
        "esophagus", # 710
        "heart", # 561
        "intestine", # 187
        "liver", # 226
        "lung", # 578
        "muscle", # 803
        "nerve", # 619
        #"ovary", # 180
        "pancreas", # 328
        "pituitary", # 283
        #"prostate", # 245
        "salivary", # 162
        "skin", # 912
        "spleen", # 241
        "stomach", # 359
        #"testis", # 361
        "thyroid", # 653
        #"uterus", # 142
        #"vagina", # 156
    ]

    samples = select_samples_by_tissue(
        tissues=tissues_sex,
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    # Return information.
    return samples


def select_samples_by_tissue(
    tissues=None,
    data_samples_tissues_persons=None,
):
    """
    Selects samples by specific tissues.

    arguments:
        tissues (list<str>): identifiers of tissues
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (list<str>): identifiers of samples

    """

    # Select samples by major tissues.
    data_samples = data_samples_tissues_persons.loc[
        data_samples_tissues_persons["tissue_major"].isin(tissues), :
    ]
    # Extract identifiers of samples.
    samples = data_samples.index.to_list()
    # Return information.
    return samples


def select_samples(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Selects genes' signals for samples of interest.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    # Select samples by exclusion.
    # Exclude samples for extra-corporeal tissues, such as cell lines.
    # Exclude samples for tissues with fewer than 50 total samples.
    samples_exclusion = select_samples_by_exclusion(
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    #print(len(samples_exclusion))
    # Select samples by inclusion.
    samples_inclusion = select_samples_by_inclusion(
        data_samples_tissues_persons=data_samples_tissues_persons,
    )
    #print(len(samples_inclusion))

    # Select samples that are in both exclusion and inclusion lists.
    if False:
        samples = utility.filter_common_elements(
            list_one=samples_exclusion,
            list_two=samples_inclusion,
        )
        utility.print_terminal_partition(level=2)
        print("Count of samples to select: " + str(len(samples)))
        # Select genes' signals for samples of interest.
        utility.print_terminal_partition(level=1)
        data_selection = data_gene_signal.loc[:, samples]

    data_exclusion = data_gene_signal.loc[
        :, ~data_gene_signal.columns.isin(samples_exclusion)
    ]
    data_inclusion = data_gene_signal.loc[
        :, data_gene_signal.columns.isin(samples_inclusion)
    ]

    # Return information.
    return data_inclusion


def select_samples_genes_by_few_tissues(
    tissues_major=None,
    data_gene_sample=None
):
    """
    Selects samples and genes by tissues of interest.

    arguments:
        tissues_major (list<str>): major tissues of interest
        data_gene_sample (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons

    """

    data_gene_sample.reset_index(
        level=["tissue_major", "tissue_minor"],
        inplace=True
    )
    data_selection = data_gene_sample.loc[
        data_gene_sample["tissue_major"].isin(tissues_major),
    ]
    data_selection.set_index(
        ["tissue_major", "tissue_minor"],
        append=False,
        drop=True,
        inplace=True
    )
    return data_selection


# Signal.


def filter_rows_columns_by_threshold_proportion(
    data=None,
    dimension=None,
    threshold=None,
    proportion=None,
):
    """
    Filters either rows or columns.

    Persistence of a row or column requires at least a specific count of values
    beyond a specific threshold.

    Filter rows by consideration of values across columns in each row.
    Filter columns by consideration of values across rows in each column.

    arguments:
        data (object): Pandas data frame of genes' signals across samples
        dimension (str): dimension to filter, either "row" or "column"
        threshold (float): minimal signal
        proportion (float): minimal proportion of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Determine count from proportion.
    if dimension == "row":
        # Filter rows by consideration of columns for each row.
        columns = data.shape[1]
        count = round(proportion * columns)
    elif dimension == "column":
        # Filter columns by consideration of rows for each columns.
        rows = data.shape[0]
        count = round(proportion * rows)

    # Determine whether values exceed threshold.
    data_threshold = (data >= threshold)
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    # This aggregation operation produces a series.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )

    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count, : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count]
    return data_pass


def select_samples_genes_signals(
    threshold=None,
    proportion_gene=None,
    proportion_sample=None,
    data_gene_signal=None
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

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples
    """

    utility.print_terminal_partition(level=2)
    print("shape of original data...")
    print(data_gene_signal.shape)

    # Filter genes by their signals across samples.
    # Filter to keep only genes with signals beyond threshold in at least one
    # sample.
    data_row = filter_rows_columns_by_threshold_proportion(
        data=data_gene_signal,
        dimension="row",
        threshold=threshold,
        proportion=proportion_gene
    )
    print("shape of data after signal filter on genes...")
    print(data_row.shape)
    # Filter samples by their signals across genes.
    # Filter to keep only samples with signals beyond threshold in at least one
    # gene.
    data_column = filter_rows_columns_by_threshold_proportion(
        data=data_row,
        dimension="column",
        threshold=threshold,
        proportion=proportion_sample
    )
    print("shape of data after signal filter on samples...")
    print(data_column.shape)

    # Return information.
    return data_column


# Families, persons, tissues, samples.

# TODO: I need a new function to organize category and quantity covariates.
# TODO: try to organize the whole "families_persons" thing
# TODO: make sure that both the "families" and the "persons" get written to the text file...

# TODO: rename this function
# TODO: within this function
# 1. extract samples, tissues, persons as simple lists
# 2. select relevant persons and their attributes
# 3. organize exports in proper formats

def extract_gene_signal_families_persons_tissues_samples(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Extracts families, persons, tissues, and samples from selection of genes'
    signals.

    As the data from GTEx do not include families, generate distinct families
    for each person for use in GCTA.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (dict<list<str>>): collection of persons, tissues, and samples
    """

    utility.print_terminal_partition(level=2)
    print("extract samples, persons, and tissues")

    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_signal.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_gene_signal_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_sample=data_transposition,
    )
    print(data_gene_signal_factor)
    # Remove indices.
    data_gene_signal_factor.reset_index(
        level=None,
        inplace=True
    )

    # Extract selection samples.
    samples = utility.collect_unique_elements(
        elements_original=data_gene_signal_factor["sample"].to_list()
    )
    # Extract selection tissues.
    tissues_major = utility.collect_unique_elements(
        elements_original=data_gene_signal_factor["tissue_major"].to_list()
    )
    # Extract selection persons.
    persons = utility.collect_unique_elements(
        elements_original=data_gene_signal_factor["person"].to_list()
    )

    # Organize families and persons.
    #data_gene_signal_factor["family"] = data_gene_signal_factor["person"]
    data_gene_signal_factor["family"] = 0
    data_families_persons = data_gene_signal_factor.loc[
        :, data_gene_signal_factor.columns.isin(["family", "person"])
    ]
    data_families_persons.rename_axis(
        "",
        axis="columns",
        inplace=True,
    )
    data_families_persons.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_families_persons.reindex()

    # TODO: Why do I set "person" as the index? This doesn't make sense...
    # TODO: Oh... I think I need "person" as index in the distribution procedure...

    data_families_persons.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )

    families_persons = utility.convert_dataframe_to_records(
        data=data_families_persons
    )

    utility.print_terminal_partition(level=1)
    print("data families persons...")
    print(data_families_persons)

    # Summary.
    utility.print_terminal_partition(level=2)
    print("count of unique samples: " + str(len(samples)))
    print("count of unique persons: " + str(len(persons)))
    print("count of unique tissues_major: " + str(len(tissues_major)))

    # Compile information.
    information = {
        "families_persons": families_persons,
        "data_families": data_families_persons,
        "persons": persons,
        "tissues_major": tissues_major,
        "samples": samples,
    }

    # Return information.
    return information


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
    utility.create_directory(path_selection)

    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation.pickle"
    )
    path_gene_count = os.path.join(
        path_selection, "data_gene_count.pickle"
    )
    path_gene_count_factor = os.path.join(
        path_selection, "data_gene_count_factor.pickle"
    )
    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    path_gene_signal_factor = os.path.join(
        path_selection, "data_gene_signal_factor.pickle"
    )

    path_families_persons = os.path.join(
        path_selection, "families_persons.tsv"
    )
    path_data_families = os.path.join(
        path_selection, "data_families.pickle"
    )
    path_persons = os.path.join(
        path_selection, "persons.txt"
    )
    path_tissues = os.path.join(
        path_selection, "tissues.txt"
    )
    path_samples = os.path.join(
        path_selection, "samples.txt"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_gene_annotation"],
        path_gene_annotation
    )
    pandas.to_pickle(
        information["data_gene_count"],
        path_gene_count
    )
    pandas.to_pickle(
        information["data_gene_count_factor"],
        path_gene_count_factor
    )
    pandas.to_pickle(
        information["data_gene_signal"],
        path_gene_signal
    )
    pandas.to_pickle(
        information["data_gene_signal_factor"],
        path_gene_signal_factor
    )

    pandas.to_pickle(
        information["data_families"],
        path_data_families
    )
    information["data_families"].to_csv(
        path_or_buf=path_families_persons,
        columns=["family", "person",],
        sep="\t",
        na_rep="NA",
        header=False,
        index=True,
    )
    pandas.to_pickle(
        information["data_families"],
        path_data_families
    )
    utility.write_file_text_list(
        elements=information["persons"],
        delimiter="\n",
        path_file=path_persons
    )
    utility.write_file_text_list(
        elements=information["tissues_major"],
        delimiter="\n",
        path_file=path_tissues
    )
    utility.write_file_text_list(
        elements=information["samples"],
        delimiter="\n",
        path_file=path_samples
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

    # Remove previous files to avoid version or batch confusion.
    path_selection = os.path.join(dock, "selection")
    utility.remove_directory(path=path_selection)

    # Read source information from file.
    source = read_source(dock=dock)

    # Organize data.
    data_gene_count = organize_data_axes_indices(
        data=source["data_gene_count"]
    )
    data_gene_signal = organize_data_axes_indices(
        data=source["data_gene_signal"]
    )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal,
    )
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_count,
    )

    # Select genes that encode proteins and are on autosomes.
    # Select annotations for protein-coding genes on autosomes.
    # Count of genes: 18366
    data_gene_annotation = select_gene_annotation(
        data_gene_annotation=source["data_gene_annotation"],
    )
    data_gene_signal_gene = select_genes(
        data_gene_annotation=data_gene_annotation,
        data_gene_signal=data_gene_signal,
    )
    data_gene_count_gene = select_genes(
        data_gene_annotation=data_gene_annotation,
        data_gene_signal=data_gene_count,
    )

    # Select samples by tissues.
    # Count of samples: 15701
    data_gene_signal_sample = select_samples(
        data_gene_signal=data_gene_signal_gene,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
    )
    data_gene_count_sample = select_samples(
        data_gene_signal=data_gene_count_gene,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
    )

    # Select genes and samples by signals.
    # threshold: 1.0
    # proportion   genes_count
    # 0.0          18366
    # 0.1          15445 <-- conservative but effective
    # 0.25         13971
    # 0.5          12435
    # 10% of samples must have signals beyond threshold for a gene to pass.
    # 50% of genes must have signals beyond threshold for a sample to pass.
    data_gene_signal_selection = select_samples_genes_signals(
        threshold=1.0,
        proportion_gene=0.1,
        proportion_sample=0.50,
        data_gene_signal=data_gene_signal_sample,
    )
    data_gene_count_selection = select_samples_genes_signals(
        threshold=1.0,
        proportion_gene=0.1,
        proportion_sample=0.50,
        data_gene_signal=data_gene_count_sample,
    )

    # Summarize original counts of samples and genes.
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_count_selection,
    )
    summarize_samples_genes(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal_selection,
    )

    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_count_selection.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_gene_count_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_sample=data_transposition,
    )
    print(data_gene_count_factor)

    # Transpose data structure.
    # Organize genes across columns and samples across rows.
    data_transposition = data_gene_signal_selection.transpose(copy=True)
    # Associate samples to persons and tissues.
    data_gene_signal_factor = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_sample=data_transposition,
    )
    print(data_gene_signal_factor)

    #########################
    # TODO:
    # export categorical and quantitative covariates of persons...
    # 1. select relevant persons
    # 2. export person properties in a convenient format (Pandas to pickle) for category procedure
    # 3. export person properties in separate files for GCTA...
    #########################

    # Extract gene signal families, persons, tissues, samples.
    collection = extract_gene_signal_families_persons_tissues_samples(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal_selection,
    )

    # Compile information.
    information = {
        "data_gene_annotation": data_gene_annotation,
        "data_gene_count": data_gene_count_selection,
        "data_gene_count_factor": data_gene_count_factor,
        "data_gene_signal": data_gene_signal_selection,
        "data_gene_signal_factor": data_gene_signal_factor,
        "families_persons": collection["families_persons"],
        "data_families": collection["data_families"],
        "persons": collection["persons"],
        "tissues_major": collection["tissues_major"],
        "samples": collection["samples"],
    }
    # Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
