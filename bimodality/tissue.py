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
import sklearn
import sklearn.datasets
import sklearn.decomposition

# Custom

import assembly
import measurement
import selection
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
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_selection = os.path.join(dock, "selection")
    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_signal": data_gene_signal,
    }


##########
# Summary


def calculate_gene_signal_mean_by_minor_tissues(
    data=None
):
    """
    Calculates the mean values of genes' signals across all samples from each
    minor tissue.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples,
            persons, and tissues

    raises:

    returns:
        (object): Pandas data frame of mean values of genes' signals for all
            tissues of all persons.

    """

    # Drop the sample designations.
    data.reset_index(level=["sample", "person"], inplace=True)
    data.drop(
        labels=["sample", "person"],
        axis="columns",
        inplace=True
    )
    # Define groups.
    groups = data.groupby(level=["tissue_major", "tissue_minor"])
    # Apply calculation to groups.
    #data_mean = groups.aggregate(lambda x: x.mean())
    data_mean = groups.aggregate(statistics.mean)
    return data_mean


##########
# Normalization, Standardization


def standardize_gene_signal(data_gene_signal=None):
    """
    Transforms values of genes' signals to standard or z-score space.

    Data has genes across rows and samples across columns.

           sample_1 sample_2 sample_3 sample_4 sample_5
    gene_1 ...      ...      ...      ...      ...
    gene_2 ...      ...      ...      ...      ...
    gene_3 ...      ...      ...      ...      ...
    gene_4 ...      ...      ...      ...      ...
    gene_5 ...      ...      ...      ...      ...

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Transform signals to standard score space.
    data_gene_signal_standard = calculate_standard_score_gene_signal_by_gene(
        data_gene_signal=data_gene_signal
    )
    print(data_gene_signal_standard.iloc[0:10, 0:10])
    # Compare summary statistics before and after transformation.
    utility.print_terminal_partition(level=3)
    print("Summary statistics for gene signals before standardization.")
    data_mean = data_gene_signal.apply(
        lambda x: x.mean(),
        axis="columns"
    )
    print("Mean")
    print(data_mean.iloc[0:10])
    data_deviation = data_gene_signal.apply(
        lambda x: x.std(),
        axis="columns"
    )
    print("Standard deviation")
    print(data_deviation.iloc[0:10])
    utility.print_terminal_partition(level=3)
    print("Summary statistics for gene signals after standardization.")
    data_mean = data_gene_signal_standard.apply(
        lambda x: x.mean(),
        axis="columns"
    )
    print("Mean")
    print(data_mean.iloc[0:10])
    data_deviation = data_gene_signal_standard.apply(
        lambda x: x.std(),
        axis="columns"
    )
    print("Standard deviation")
    print(data_deviation.iloc[0:10])
    return data_gene_signal_standard


def calculate_standard_score_gene_signal_by_gene(data_gene_signal=None):
    """
    Calculates the standard (z-score) of genes' signals across each gene.

    The standard scores are relative to gene.
    The values of mean and standard deviation are across all samples for each
    gene.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    data_gene_standard = data_gene_signal.apply(
        lambda x: (x - x.mean()) / x.std(),
        axis="columns",
    )
    return data_gene_standard


def normalize_standardize_gene_signal(
    data_gene_signal=None
):
    """
    Normalizes and standardizes genes' signals across all samples.

    Data has genes across rows and samples across columns.

           sample_1 sample_2 sample_3 sample_4 sample_5
    gene_1 ...      ...      ...      ...      ...
    gene_2 ...      ...      ...      ...      ...
    gene_3 ...      ...      ...      ...      ...
    gene_4 ...      ...      ...      ...      ...
    gene_5 ...      ...      ...      ...      ...

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples

    """

    # Transform genes' signals to logarithmic space.
    # Transform values of genes' signals to base-2 logarithmic space.
    print(data_gene_signal)
    data_normal = measurement.transform_gene_signal_log(
        data_gene_signal=data_gene_signal,
        pseudo_count=1.0,
    )
    # Transform genes' signals to standard or z-score space.
    # Calculate standard score for each gene.
    # Standard score is undefined if the series has inadequate variance.
    # Consider filtering before calculation.
    # Alternatively, filter afterwards.
    data_normal_standard = standardize_gene_signal(
        data_gene_signal=data_normal
    )
    # Some genes have inadequate variance, such as due to having too few
    # nonzero measurements, to have real standard values.
    # Remove genes, not samples, with missing values.
    print("shape of original data frame: " + str(data_normal_standard.shape))
    data_normal_standard.dropna(
        axis="index",
        inplace=True,
    )
    print("shape of data frame: " + str(data_normal_standard.shape))
    return data_normal_standard


##########
# Principle components


def calculate_principal_components(
    data=None,
    components=None,
):
    """
    Calculates the principal components.

    arguments:
        data (object): Pandas data frame of signals with features across rows
            and observations across columns
        components (int): count of principle components

    raises:

    returns:
        (object): Pandas data frame of principle components for each factor

    """

    # Organize data for principle component analysis.
    # Organize features (genes) across columns and instances (samples) across
    # rows.
    data = data.transpose(copy=True)
    # Organize data as an array of arrays.
    matrix = data.to_numpy()
    # Execute principle component analysis.
    pca = sklearn.decomposition.PCA(n_components=components)
    #report = pca.fit_transform(matrix)
    pca.fit(matrix)
    # Report extent of variance that each principal component explains.
    variance_ratios = pca.explained_variance_ratio_
    component_numbers = range(len(variance_ratios) + 1)
    variance_series = {
        "component": component_numbers[1:],
        "variance": variance_ratios
    }
    data_component_variance = pandas.DataFrame(data=variance_series)
    # Transform data by principal components.
    matrix_component = pca.transform(matrix)
    # Match components to samples.
    samples = data.index.to_list()
    records = []
    for index, row in enumerate(matrix_component):
        record = {}
        record["sample"] = samples[index]
        for count in range(components):
            name = "component_" + str(count + 1)
            record[name] = row[count]
        records.append(record)
    data_component = (
        utility.convert_records_to_dataframe(records=records)
    )
    data_component.set_index(
        ["sample"],
        append=False,
        drop=True,
        inplace=True
    )
    # Compile information.
    information = {
        "data_sample_component": data_component,
        "data_component_variance": data_component_variance,
    }
    # Return information.
    return information


def calculate_report_gene_sample_principal_components(
    data=None,
    data_samples_factors=None,
    components=None,
):
    """
    Calculates the principal components for genes as features and samples as
    observations.

    arguments:
        data (object): Pandas data frame of signals with features across rows
            and observations across columns
        data_samples_factors (object): Pandas data frame of factors for each
            sample
        components (int): count of principle components

    raises:

    returns:
        (object): Pandas data frame of principle components for each factor

    """

    # Describe variance across categories of tissues.
    # Normalize and standardized gene's signals for principal components.
    data_normal_standard = normalize_standardize_gene_signal(
        data_gene_signal=data
    )
    print("Data after normalization and standardization.")
    print(data_normal_standard)

    report = calculate_principal_components(
        data=data_normal_standard,
        components=10
    )
    utility.print_terminal_partition(level=2)
    print("Report from principal component analysis...")
    print("Explained variance by each principal component...")
    print(report["data_component_variance"])
    utility.print_terminal_partition(level=3)
    print(report["data_sample_component"])
    # Associate samples to major and minor tissue types.
    data_factor_component = assembly.associate_samples_persons_tissues(
        data_samples_tissues_persons=data_samples_factors,
        data_gene_sample=report["data_sample_component"],
    )
    utility.print_terminal_partition(level=3)
    print(data_factor_component)

    # Compile information.
    information = {
        "data_component_variance": report["data_component_variance"],
        "data_sample_component": report["data_sample_component"],
        "data_factor_component": data_factor_component,
    }
    # Return information.
    return information




##########
# Differential expression


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


def organize_differential_expression(
    data_samples_tissues_patients=None,
    data_gene_count=None,
):
    """
    Organize data for differential expression analysis.

    arguments:
        data_samples_tissues_patients (object): Pandas data frame of patients
            and tissues for all samples
        data_gene_count (object): Pandas data frame of genes' counts for all
            samples

    raises:

    returns:

    """

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
    pass


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
    path_tissue = os.path.join(dock, "tissue")
    utility.confirm_path_directory(path_tissue)
    path_gene_factor = os.path.join(
        path_tissue, "data_gene_signal_factor.pickle"
    )
    path_component_variance_all = os.path.join(
        path_tissue, "data_component_variance_all.pickle"
    )
    path_factor_component_all = os.path.join(
        path_tissue, "data_factor_component_all.pickle"
    )
    path_component_variance_few = os.path.join(
        path_tissue, "data_component_variance_few.pickle"
    )
    path_factor_component_few = os.path.join(
        path_tissue, "data_factor_component_few.pickle"
    )

    # Write information to file.
    if False:
        pandas.to_pickle(
            information["data_gene_signal_factor"],
            path_gene_factor
        )
        pass

    pandas.to_pickle(
        information["data_component_variance_all"],
        path_component_variance_all
    )
    pandas.to_pickle(
        information["data_factor_component_all"],
        path_factor_component_all
    )

    pandas.to_pickle(
        information["data_component_variance_few"],
        path_component_variance_few
    )
    pandas.to_pickle(
        information["data_factor_component_few"],
        path_factor_component_few
    )

    if False:
        information["data_component_variance"].to_csv(
            path_or_buf=path_component_variance_text,
            sep="\t",
            header=True,
            index=False,
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

    # Calculate principal components.
    # 1. Compare all 19 non-sexual tissues with samples for more than 100
    # persons.
    utility.print_terminal_partition(level=1)
    print("Principal components across all 19 non-sexual tissues...")
    utility.print_terminal_partition(level=2)
    report_component_all = calculate_report_gene_sample_principal_components(
        data=source["data_gene_signal"],
        data_samples_factors=source["data_samples_tissues_persons"],
        components=10
    )

    # Calculate principal components.
    # 2. Compare remainder tissues after removal of blood, brain, and
    # pituitary, as these tissues seem to dominate variance.
    utility.print_terminal_partition(level=1)
    print("Principal components across fewer tissues...")
    utility.print_terminal_partition(level=2)
    # Define tissues for inclusion.
    # These tissues are not sex-specific.
    tissues = [
        "adipose", # 552
        "adrenal", # 190
        "artery", # 551
        #"blood", # 407
        #"brain", # 254
        "colon", # 371
        "esophagus", # 513
        "heart", # 399
        "liver", # 175
        "lung", # 427
        "muscle", # 564
        "nerve", # 414
        "pancreas", # 248
        #"pituitary", # 183
        "skin", # 583
        "intestine", # 137
        "spleen", # 162
        "stomach", # 262
        "thyroid", # 446
    ]
    samples = selection.select_samples_by_tissue(
        tissues=tissues,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
    )
    data_few = selection.select_samples_genes(
        samples=samples,
        threshold=True,
        data_gene_signal=source["data_gene_signal"],
    )
    print(data_few)
    report_component_few = calculate_report_gene_sample_principal_components(
        data=data_few,
        data_samples_factors=source["data_samples_tissues_persons"],
        components=10
    )

    # Calculate mean of genes' signals across all samples for each minor tissue
    # type.
    if False:
        data_gene_tissue = calculate_gene_signal_mean_by_minor_tissues(
            data=source["data_gene_sample"]
        )
        print(data_gene_tissue)

        # Prepare report of data after normalization and standardization.
        # Transpose data structure.
        # Organize genes across columns and samples across rows.
        data_transposition = data_gene_signal_normal_standard.transpose(copy=True)
        # Associate samples to persons and tissues.
        data_gene_signal_factor = assembly.associate_samples_persons_tissues(
            data_samples_tissues_persons=source["data_samples_tissues_persons"],
            data_gene_sample=data_transposition,
        )
        print(
            "Genes' signals after normalization and standardization with " +
            "samples and factors."
        )
        print(data_gene_signal_factor)

    # Compile information.
    information = {
        "data_component_variance_all": (
            report_component_all["data_component_variance"]
        ),
        "data_factor_component_all": (
            report_component_all["data_factor_component"]
        ),
        "data_component_variance_few": (
            report_component_few["data_component_variance"]
        ),
        "data_factor_component_few": (
            report_component_few["data_factor_component"]
        ),
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
