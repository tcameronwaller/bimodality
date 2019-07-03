"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os
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
    path_split = os.path.join(dock, "split")
    path_genes = os.path.join(
        path_split, "genes.txt"
    )
    path_signal = os.path.join(
        path_split, "genes_samples_signals.pickle"
    )
    # Read information from file.
    genes = utility.read_file_text_list(path_genes)
    with open(path_signal, "rb") as file_source:
        genes_samples_signals = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes": genes,
        "genes_samples_signals": genes_samples_signals,
    }


##########
# Association of genes, samples, tissues, and persons


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

    data_samples_tissues_persons.reset_index(level="sample", inplace=True)
    samples_tissues_persons = utility.convert_dataframe_to_records(
        data=data_samples_tissues_persons
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
    data_gene_signal=None
):
    """
    Collects samples, persons, and tissues for each gene's signal.

    arguments:
        reference (dict<dict<str>>): Tissue and person for each sample.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples.

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
    data = data_gene_signal.transpose(copy=True)
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


def associate_genes_samples_persons_tissues(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Associates samples, persons, and tissues for each gene's signal.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, persons,
            and tissues

    """

    # Associate genes' signals to persons and tissues.
    utility.print_terminal_partition(level=2)
    print(
        "Association of genes' signals to samples, tissues, and persons."
    )
    print(data_samples_tissues_persons)
    print(data_gene_signal)

    # Create reference for tissues and persons of each sample.
    reference = collect_samples_persons_tissues_reference(
        data_samples_tissues_persons=data_samples_tissues_persons
    )

    # Include identifiers for person and tissue for each sample in data for
    # genes' signals.
    data_gene_signal = collect_genes_samples_persons_tissues(
        reference=reference,
        data_gene_signal=data_gene_signal
    )
    #print(data_gene_signal)
    print(data_gene_signal.iloc[0:10, 0:7])

    # Return information.
    return data_gene_signal


##########
#Aggregation.


def calculate_mean_median_signal_by_persons_tissues(
    data_gene_signal_assembly=None,
    data_gene_signal_selection=None,
):
    """
    Calculates the mean and median value of genes' signals in groups by
    persons and tissues.



    arguments:
        data_gene_signal_assembly (object): Pandas data frame of genes' signals
            for all samples, tissues, and persons.
        data_gene_signal_selection (object): Pandas data frame of genes'
            signals for samples, tissues, and persons of interest.

    raises:

    returns:
        (dict): Pandas data frames of genes' signals.

    """

    # Organization of signals.
    utility.print_terminal_partition(level=1)
    print("Aggregation of signals.")
    print(
        "Calculation of mean and median values of genes' signals in groups " +
        "by persons and tissues."
    )

    # Calculate mean signal for each gene in all samples for each person and
    # tissue.
    utility.print_terminal_partition(level=2)
    print(
        "Calculation of mean values of signals for each gene across " +
        "redundant samples for each person and tissue."
    )
    print(
        "Calculate mean values before imputation to avoid bias for persons " +
        "and tissues with extra samples."
    )
    data_gene_signal_mean = calculate_gene_signal_mean_by_persons_tissues(
        data_gene_signal=data_gene_signal_selection
    )
    print(data_gene_signal_mean.iloc[0:10, 0:7])

    # Calculate median signals.
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
    print("Calculation of median values of each gene's signal " +
        "within each tissue across all persons."
    )
    print("These median values are for imputation.")
    print(
        "Calculate imputation values before selection of samples, persons, " +
        "and tissues so as to use as many measurements as are available."
    )
    # The map of associations between persons, tissues, and samples does not
    # need filtration.
    # The lists of tissues and persons of interest control the selection of
    # samples.
    data_gene_signal_median = calculate_gene_signal_median_by_tissues(
        data_gene_signal=data_gene_signal_assembly
    )
    print(data_gene_signal_median.iloc[0:25, 0:10])

    # Compile information.
    information = {
        "mean": data_gene_signal_mean,
        "median": data_gene_signal_median
    }

    # Return information.
    return information


def calculate_gene_signal_mean_by_persons_tissues(
    data_gene_signal=None
):
    """
    Calculates the mean values of genes' signals across all samples from each
    tissue of each person.

    Redundant samples exist for many persons and tissues.
    Calculate the mean signal for each gene from all of these redundant
    samples.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals for
            samples, tissues, and persons.

    raises:

    returns:
        (object): Pandas data frame of mean values of genes' signals for all
            tissues of all persons.

    """

    # Drop the sample designations.
    data_gene_signal.reset_index(level=["sample"], inplace=True)
    data_gene_signal.drop(
        labels="sample",
        axis="columns",
        inplace=True
    )
    # Define groups.
    groups = data_gene_signal.groupby(level=["person", "tissue"])
    # Apply calculation to groups.
    #data_mean = groups.aggregate(lambda x: x.mean())
    data_mean = groups.aggregate(statistics.mean)
    return data_mean


def calculate_gene_signal_median_by_tissues(
    data_gene_signal=None
):
    """
    Calculates the median values of signals for all genes across specific
    persons and tissues.

    Redundant samples exist for many persons and tissues.
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
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons.

    raises:

    returns:
        (object): Pandas data frame of median values of genes' signals for all
            tissues.

    """

    # Calculate mean values by person and tissue.
    data_gene_signal_mean = calculate_gene_signal_mean_by_persons_tissues(
        data_gene_signal=data_gene_signal
    )
    data_gene_signal_mean.reset_index(level=["person"], inplace=True)
    data_gene_signal_mean.drop(
        labels="person",
        axis="columns",
        inplace=True
    )
    # Define groups.
    groups = data_gene_signal_mean.groupby(level=["tissue"])
    # Apply calculation to groups.
    #data_median = groups.aggregate(lambda x: x.median())
    data_median = groups.aggregate(statistics.median)
    return data_median


##################################################
# Scrap


def calculate_gene_signal_mean_by_persons_tissues_old(
    persons_tissues_samples=None,
    data_gene_signal=None
):
    """
    Calculates the mean values of genes' signals in each tissue of each
    person.

    Redundant samples exist for many persons and tissues.
    Calculate the mean signal for each gene from all of these redundant
    samples.

    Product data format.
    Indices by person and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # person tissue
    # bob     kidney  5.0    10.0   7.0
    # bob     skin    5.0    10.0   7.0
    # bob     liver   5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # april   liver   2.0    3.0    4.0
    # april   skin    2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        data_gene_signal (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of mean values of genes' signals for all
            tissues of all persons.

    """

    # Create index for convenient access.
    data_gene_signal_index = data_gene_signal.set_index("Name")
    # Determine list of all genes with signals.
    genes = data_gene_signal["Name"].to_list()
    # Collect mean signals of all genes for each person and tissue.
    summary = list()
    # Iterate on persons.
    for person in persons_tissues_samples:
        # Iterate on tissues for which person has samples.
        for tissue in persons_tissues_samples[person]:
            # Collect signals for all genes in each sample.
            signals_genes = dict()
            for sample in persons_tissues_samples[person][tissue]:
                for gene in genes:
                    signal = data_gene_signal_index.at[gene, sample]
                    # Determine whether an entry already exists for the gene.
                    if gene in signals_genes:
                        signals_genes[gene].append(signal)
                    else:
                        signals_genes[gene] = list([signal])
            # Compile record.
            record = dict()
            record["person"] = person
            record["tissue"] = tissue
            # Calculate mean signal for each gene in each person and tissue.
            for gene in genes:
                signal_mean = statistics.mean(signals_genes[gene])
                record[gene] = signal_mean
            # Include record.
            summary.append(record)
    # Organize data.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    #data_summary_index = data_summary.set_index(
    #    ["person", "tissue"], append=False, drop=False
    #)
    return data_summary


def calculate_gene_signal_median_by_tissues_old(
    tissues=None,
    persons=None,
    tissues_persons_samples=None,
    data_gene_signal_mean=None
):
    """
    Calculates the median values of signals for all genes across specific
    persons and tissues.

    Redundant samples exist for many persons and tissues.
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
        persons (list<str>): persons with signal for tissues of interest.
        tissues_persons_samples (dict<dict<list<str>>): Samples for each
            person of each tissue.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of median values of genes' signals for all
            tissues.

    """

    # Create index for convenient access.
    data_gene_signal_mean_index = data_gene_signal_mean.set_index(
        ["person", "tissue"], append=False, drop=True
    )
    # Extract genes from names of columns.
    genes = data_gene_signal_mean_index.columns.to_list()
    # Collect median signals of all genes for each tissue.
    summary = list()
    # Iterate on genes.
    for gene in genes:
        # Compile record.
        record = dict()
        record["gene"] = gene
        # Iterate on tissues.
        for tissue in tissues:
            # Collect gene's signals in the tissue across all persons.
            signals_tissue = list()
            # Iterate on persons.
            for person in tissues_persons_samples[tissue]:
                signal_tissue_person = data_gene_signal_mean_index.loc[
                    (person, tissue), gene
                ]
                signals_tissue.append(signal_tissue_person)
            # Calculate median value of gene's signal in tissue across all
            # persons.
            signal_tissue = statistics.median(signals_tissue)
            record[tissue] = signal_tissue
        # Include record.
        summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    return data_summary


def collect_gene_signal_by_persons_tissues_old(
    tissues=None,
    persons=None,
    persons_tissues_samples=None,
    data_gene_signal_median_tissue=None,
    data_gene_signal_mean=None
):
    """
    Collects the signals for all genes in each tissue of each person.

    Product data format.
    Indices by person and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # person tissue
    # bob     kidney  5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons with signal for tissues of interest.
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        data_gene_signal_median_tissue (object): Pandas data frame of median
            values of genes' signals for all tissues.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific persons and tissues.

    """

    # Create indices for convenient access.
    data_gene_signal_mean_index = data_gene_signal_mean.set_index(
        ["person", "tissue"], append=False, drop=True
    )
    data_gene_signal_median_tissue_index = (
        data_gene_signal_median_tissue.set_index("gene")
    )
    # Extract genes from names of columns.
    genes = data_gene_signal_mean_index.columns.to_list()
    # Collect signals of all genes for each person and tissue.
    summary = list()
    # Iterate on persons of interest.
    for person in persons:
        # Iterate on tissues of interest.
        for tissue in tissues:
            # Compile record.
            record = dict()
            record["person"] = person
            record["tissue"] = tissue
            # Determine whether the person has samples for the tissue.
            if tissue in persons_tissues_samples[person]:
                # Collect mean values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_mean_index.loc[
                        (person, tissue), gene
                    ]
                    record[gene] = signal
            else:
                # Collect imputation values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_median_tissue_index.loc[
                        gene, tissue
                    ]
                    record[gene] = signal
            # Include record.
            summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    #data_summary_index = data_summary.set_index(
    #    ["person", "tissue"], append=False, drop=False
    #)
    return data_summary

##################################################


##########
# Imputation.


def organize_genes_signals(
    tissues=None,
    persons=None,
    persons_tissues_samples=None,
    data_gene_signal_median=None,
    data_gene_signal_mean=None
):
    """
    Organizes the signals for all genes in each tissue of each person.

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons with signal for tissues of interest.
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        data_gene_signal_median_tissue (object): Pandas data frame of median
            values of genes' signals for all tissues.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific persons and tissues.

    """

    # Organization of signals.
    utility.print_terminal_partition(level=1)
    print("Organization of signals.")
    print(
        "Imputation of genes' signals for some persons and tissues will " +
        "enhance coverage."
    )

    # Collection of genes' aggregate signals for all persons and tissues.
    utility.print_terminal_partition(level=2)
    print(
        "Collection of genes' signals for all persons and tissues with " +
        "imputation of misses."
    )
    data_gene_signal_imputation = collect_gene_signal_real_imputation(
        tissues=tissues,
        persons=persons,
        persons_tissues_samples=persons_tissues_samples,
        data_gene_signal_median=data_gene_signal_median,
        data_gene_signal_mean=data_gene_signal_mean
    )

    # Create indices for convenient access.
    data_gene_signal_imputation.set_index(
        ["person", "tissue"],
        append=False,
        drop=True,
        inplace=True
    )
    data_gene_signal_imputation.rename_axis(
        columns="gene",
        axis="columns",
        copy=False,
        inplace=True
    )

    print(data_gene_signal_imputation.iloc[0:25, 0:5])
    print(data_gene_signal_imputation.shape)

    return data_gene_signal_imputation


def collect_gene_signal_real_imputation(
    tissues=None,
    persons=None,
    persons_tissues_samples=None,
    data_gene_signal_median=None,
    data_gene_signal_mean=None
):
    """
    Collects the signals for all genes in each tissue of each person.

    Product data format.
    Indices by person and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # person tissue
    # bob     kidney  5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons with signal for tissues of interest.
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        data_gene_signal_median_tissue (object): Pandas data frame of median
            values of genes' signals for all tissues.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific persons and tissues.

    """

    # Extract genes from names of columns.
    genes = data_gene_signal_mean.columns.to_list()
    # Collect signals of all genes for each person and tissue.
    summary = list()
    # Iterate on persons of interest.
    for person in persons:
        # Iterate on tissues of interest.
        for tissue in tissues:
            # Compile record.
            record = dict()
            record["person"] = person
            record["tissue"] = tissue
            # Determine whether the person has samples for the tissue.
            if tissue in persons_tissues_samples[person]:
                # Collect mean values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_mean.loc[
                        (person, tissue), gene
                    ]
                    record[gene] = signal
            else:
                # Collect imputation values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_median.loc[
                        tissue, gene
                    ]
                    record[gene] = signal
            # Include record.
            summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    return data_summary



##########
# Transformation to standard space.


def transform_gene_signal_standard(data_gene_signal=None):
    """
    Transforms values of genes' signals to standard or z-score space.

    arguments:
        data_gene_signal (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific persons and tissues.

    """

    # Transform signals to standard score space.
    utility.print_terminal_partition(level=2)
    print(
        "Transformation of genes' signals to standard score (z-score) " +
        "space."
    )
    data_gene_signal_standard = calculate_standard_score_gene_signal_by_tissue(
        data_gene_signal=data_gene_signal
    )
    print(data_gene_signal_standard.iloc[0:10, 0:10])
    # Compare summary statistics before and after transformation.
    utility.print_terminal_partition(level=3)
    print("Summary statistics for gene signals in base-2 logarithmic space.")
    groups = data_gene_signal.groupby(level="tissue")
    print(groups.describe())
    utility.print_terminal_partition(level=3)
    print("Summary statistics for gene signals in z-score space.")
    groups_standard = data_gene_signal_standard.groupby(level="tissue")
    print(groups_standard.describe())
    print("Mean...")
    print(groups_standard.mean())
    print("Standard deviation...")
    print(groups_standard.std())

    return data_gene_signal_standard


def calculate_standard_score_gene_signal_by_tissue(data_gene_signal=None):
    """
    Calculates the standard (z-score) of genes' signals for each person and
    tissue.

    The standard scores are relative to tissue.
    The values of mean and standard deviation are across all persons for each
    tissue.

    arguments:
        data_gene_signal (object): Pandas data frame of signals for all genes
            across specific persons and tissues.

    raises:

    returns:
        (object): Pandas data frame of standard score signals for all genes
            across specific persons and tissues.

    """

    # lambda x: math.log((x + 1), 2)

    groups = data_gene_signal.groupby(level="tissue")
    # For some reason, the scipy function throws an error about illegal values.
    #groups.transform(lambda value: print(value.to_list()))
    #data_standard_index = (
    #    groups.transform(lambda value: scipy.stats.zscore(value.to_list()))
    #)
    data_standard = groups.transform(lambda x: (x - x.mean()) / x.std())
    return data_standard


########################################################
# TODO: go ahead and calculate z-scores...
########################################################

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
    path_organization = os.path.join(dock, "organization")
    utility.confirm_path_directory(path_organization)
    path_mean = os.path.join(
        path_organization, "data_gene_signal_mean.pickle"
    )
    path_median = os.path.join(
        path_organization, "data_gene_signal_median.pickle"
    )
    path_imputation = os.path.join(
        path_organization, "data_gene_signal_imputation.pickle"
    )
    path_log = os.path.join(
        path_organization, "data_gene_signal_log.pickle"
    )
    path_standard = os.path.join(
        path_organization, "data_gene_signal_standard.pickle"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_gene_signal_mean"], path_mean
    )
    pandas.to_pickle(
        information["data_gene_signal_median"], path_median
    )
    pandas.to_pickle(
        information["data_gene_signal_imputation"], path_imputation
    )
    pandas.to_pickle(
        information["data_gene_signal_log"], path_log
    )
    pandas.to_pickle(
        information["data_gene_signal_standard"], path_standard
    )

    pass


###############################################################################
# Procedure


def execute_procedure(
    selection=None,
    count=None,
    data_gene_sample_signal=None
):
    """
    Function to execute module's main behavior.

    arguments:
        selection (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        data_gene_sample_signal (object): Pandas data frame of a gene's signals
            across samples

    raises:

    returns:
        (dict): report about gene, list of gene's signals across persons

    """

    # Aggregate genes' signals within groups by person and major tissue
    # category.
    data_gene_signal_tissue = aggregate_gene_signal_tissue(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=source["data_gene_signal"],
    )


    print(source["data_gene_signal_assembly"].iloc[0:25, 0:10])
    print(source["data_gene_signal_assembly"].shape)
    print(source["data_gene_signal_selection"].iloc[0:25, 0:10])
    print(source["data_gene_signal_selection"].shape)

    #data_gene_signal_assembly = source["data_gene_signal_assembly"].iloc[ : , 0:1000]
    #data_gene_signal_selection = source["data_gene_signal_selection"].iloc[ : , 0:1000]

    genes = source["data_gene_signal_selection"].columns.to_list()
    print(str(len(genes)))
    genes_new = genes[0:100]

    #data_gene_signal_assembly = source["data_gene_signal_assembly"].loc[ : , genes_new]
    #data_gene_signal_selection = source["data_gene_signal_selection"].loc[ : , genes_new]
    data_gene_signal_assembly = source["data_gene_signal_assembly"]
    data_gene_signal_selection = source["data_gene_signal_selection"]



    # Calculate mean and median values of genes' signals for each group by
    # person and tissue.
    aggregation = calculate_mean_median_signal_by_persons_tissues(
        data_gene_signal_assembly=data_gene_signal_assembly,
        data_gene_signal_selection=data_gene_signal_selection,
    )

    # Collect real and imputation values of genes' signals for persons and
    # tissues of interest.
    data_gene_signal_imputation = organize_genes_signals(
        tissues=source["tissues"],
        persons=source["persons"],
        persons_tissues_samples=source["persons_tissues_samples"],
        data_gene_signal_median=aggregation["median"],
        data_gene_signal_mean=aggregation["mean"]
    )

    # Transform values of genes' signals to base-2 logarithmic space.
    data_gene_signal_log = transform_gene_signal_log(
        data_gene_signal=data_gene_signal_imputation
    )

    # Transform values of genes' signals to standard or z-score space.
    data_gene_signal_standard = transform_gene_signal_standard(
        data_gene_signal=data_gene_signal_log
    )

    # Compile information.
    information = {
        "report": report,
        "gene_persons_signals": gene_persons_signals,
    }
    # Return information.
    return information





    ###################################################################

    # TODO: I think organize the entire mean, median, imputation stuff within its own function...


    if False:
        print(data_gene_signal_mean.shape)
        data_gene_signal_mean_index = data_gene_signal_mean.set_index(
            ["person", "tissue"], append=False, drop=True
        )
        print(data_gene_signal_mean_index.iloc[0:25, 0:5])
        # Demonstration of access to single values with multiple levels of indices.
        print("Access to values at degrees of reference.")
        print(data_gene_signal_mean_index.loc[
            ("GTEX-111CU", "Skin"), "ENSG00000078808.12"
        ])
        print(data_gene_signal_mean_index.loc[
            ("GTEX-111CU", "Muscle"), "ENSG00000078808.12"
        ])
        print(data_gene_signal_mean_index.loc[
            ("GTEX-111CU", "Nerve"), "ENSG00000078808.12"
        ])
        print(data_gene_signal_mean_index.loc["GTEX-111CU", ])
        # Extract genes from names of columns.
        #headers = data_gene_signal_mean_index.columns.to_list()
        #print(headers)
        #print("count of genes: " + str(len(headers)))

        # Reverse an index to a column.
        #dataframe["new_column"] = dataframe.index
        #dataframe.reset_index(level=["index_name"])

    pass


if (__name__ == "__main__"):
    execute_procedure()
