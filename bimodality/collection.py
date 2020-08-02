"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os
import copy
import pickle
import itertools
import statistics
import math
import numpy

# Relevant

import pandas
import gseapy

# Custom

import utility
import assembly
import integration

#dir()
#importlib.reload()

###############################################################################
# Functionality



##########
# Initialization


def initialize_directories(dock=None):
    """
    Initialize directories for procedure's product files.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = dock
    paths["collection"] = os.path.join(paths["dock"], "collection")
    paths["covid19"] = os.path.join(paths["collection"], "covid19")
    paths["cytokine"] = os.path.join(paths["collection"], "cytokine")
    paths["integrin"] = os.path.join(paths["collection"], "integrin")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=paths["collection"])
    utility.create_directories(path=paths["collection"])
    utility.create_directories(path=paths["covid19"])
    utility.create_directories(path=paths["cytokine"])
    utility.create_directories(path=paths["integrin"])
    # Return information.
    return paths


##########
# Collection of COVID-19 genes


def read_source_covid19_genes(
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    genes_candidacy = integration.read_source_genes_sets_candidacy(
        cohort="selection",
        dock=dock,
    )
    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    translations_genes = (
        assembly.read_source_gene_name_identifier_translations(
            dock=dock
    ))
    path_collection = os.path.join(
        dock, "annotation", "gene_sets", "genes_covid19",
        "collection_2020-07-20"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    # Read file names.
    files_collection = utility.extract_directory_file_names(
        path=path_collection,
    )
    # Read and organize data.
    bin = dict()
    utility.print_terminal_partition(level=2)
    print("reading files...")
    utility.print_terminal_partition(level=3)
    for file_name in files_collection:
        utility.print_terminal_partition(level=4)
        print(file_name)
        name = file_name.split(".")[0]
        study = str(name.replace("genes_", ""))
        identifier = str(name.replace("genes_pubmed_", ""))
        data_key = str("data_" + name)
        print(data_key)
        # Specify directories and files.
        path_file = os.path.join(path_collection, file_name)
        # Read information from file.
        bin[study] = dict()
        bin[study]["study"] = study
        bin[study]["identifier"] = identifier
        bin[study]["data"] = pandas.read_csv(
            path_file,
            sep="\t",
            header=0,
            low_memory=False,
        )
        pass
    # Return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "translations_genes": translations_genes,
        "genes_selection": genes_candidacy["any"],
        "bin_studies": bin,
    }


def translate_study_genes_identifiers(
    bin_studies=None,
    data_gene_annotation=None,
    translations_genes=None,
    report=None,
):
    """
    Translates genes' names from all studies to Ensembl identifiers.

    arguments:
        bin_studies (dict): collection of information about each study
        data_gene_annotation (object): Pandas data frame of genes' annotations
        translations_genes (dict<str>): pairwise custom translations of genes'
            names to Ensembl identifiers, see
            assembly.read_source_gene_name_identifier_translations()
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about each study with genes'
            identifiers

    """

    utility.print_terminal_partition(level=2)
    print("translating genes' names to Ensembl identifiers...")
    print("following genes' names do not match...")
    utility.print_terminal_partition(level=2)
    # Iterate on studies.
    bin_studies = copy.deepcopy(bin_studies)
    for study in bin_studies.keys():
        data_study = bin_studies[study]["data"]
        # Determine whether the study already includes genes' identifiers.
        if ("identifier" not in data_study.columns.to_list()):
            data_study["identifier"] = data_study["name"].apply(
                lambda gene_name:
                    assembly.translate_gene_name_to_identifier(
                        name=gene_name,
                        data_gene_annotation=data_gene_annotation,
                        translations_genes=translations_genes,
                    )
            )
            bin_studies[study]["data"] = data_study
            # Report.
            if report:
                print(data_study)
                pass
            pass
        pass
    utility.print_terminal_partition(level=2)
    print("end translation...")
    utility.print_terminal_partition(level=2)
    return bin_studies


def collect_studies_unique_gene_identifiers(
    bin_studies=None,
):
    """
    Collects unique identifiers of genes from all studies.

    arguments:
        bin_studies (dict): collection of information about each study

    raises:

    returns:
        (list<str>): unique identifiers of genes from all studies

    """

    genes_collection = []
    for study in bin_studies.keys():
        data_study = bin_studies[study]["data"]
        genes_study = data_study["identifier"].to_list()
        genes_collection.extend(genes_study)
    # Determine valid, non null values of the gene's fold change.
    genes_valid = list(filter(
        lambda identifier: ("ENSG" in str(identifier)),
        genes_collection
    ))
    genes_unique = utility.collect_unique_elements(
        elements_original=genes_valid
    )
    return genes_unique


def translate_study_comparisons_identifiers(
    bin_studies=None,
    report=None,
):
    """
    Collects and organizes unique designations of comparisons from all studies.

    arguments:
        bin_studies (dict): collection of information about each study
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about each study with unique
            comparisons in each study

    """

    bin_studies = copy.deepcopy(bin_studies)
    for study in bin_studies.keys():
        identifier_study = bin_studies[study]["identifier"]
        data_study = bin_studies[study]["data"]
        columns = data_study.columns.to_list()
        comparisons = list(filter(
            lambda value: (not value in ["identifier", "name"]),
            columns
        ))
        comparisons_unique = sorted(utility.collect_unique_elements(
            elements_original=comparisons
        ))
        bin_studies[study]["comparisons"] = comparisons_unique
        # Organize study's unique comparisons.
        bin_studies[study]["comparisons_translation"] = dict()
        for comparison in comparisons_unique:
            name = str(identifier_study + "_" + comparison)
            bin_studies[study]["comparisons_translation"][comparison] = name
            pass
        # Translate comparison columns.
        data_study.rename(
            columns=bin_studies[study]["comparisons_translation"],
            inplace=True,
        )
        bin_studies[study]["data"] = data_study
        # Report.
        if report:
            print(data_study)
            pass
        pass
    return bin_studies


def collect_studies_unique_comparisons(
    bin_studies=None,
):
    """
    Collects unique identifiers of comparisons from all studies.

    arguments:
        bin_studies (dict): collection of information about each study

    raises:

    returns:
        (list<str>): unique identifiers of comparisons from all studies

    """

    collection = []
    for study in bin_studies.keys():
        data_study = bin_studies[study]["data"]
        columns = data_study.columns.to_list()
        comparisons = list(filter(
            lambda value: (not value in ["identifier", "name"]),
            columns
        ))
        collection.extend(comparisons)
    comparisons_unique = sorted(utility.collect_unique_elements(
        elements_original=collection
    ))
    return comparisons_unique


def determine_gene_study_comparison_value(
    gene_identifier=None,
    comparison=None,
    data_study=None,
    warn=None,
):
    """
    Determines a gene's value of fold change for a comparison in a study.

    The function that calls this function already verifies that the study
    includes the gene and the comparison.

    arguments:
        gene_identifier (str): unique identifier of a gene
        comparison (str): name of a comparison
        data_study (object): Pandas data frame of information about comparisons
            across genes in a study
        warn (bool): whether to print warnings

    raises:

    returns:
        (float): gene's value of fold change for a comparison in a study

    """

    # Select study's information for gene.
    data_study = data_study.copy(deep=True)
    data_study_gene = data_study.loc[
        data_study["identifier"] == gene_identifier, :
    ].copy(deep=True)
    data_study_gene.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
        #ignore_index=True,
    )
    # Determine valid, non null values of the gene's fold change.
    values_valid = list(filter(
        lambda value: (not math.isnan(value)),
        data_study_gene[comparison].to_list()
    ))
    if len(values_valid) > 1:
        value = statistics.mean(values_valid)
        if warn:
            utility.print_terminal_partition(level=3)
            print(
                "warning: gene has multiple fold change values for a single " +
                "study and comparison." +
                gene_identifier
            )
    elif len(values_valid) == 1:
        value = values_valid[0]
    else:
        value = float("nan")
    # Return information.
    return value


def collect_genes_annotations_studies_comparisons_valid_change(
    genes_identifiers=None,
    bin_studies=None,
    data_gene_annotation=None,
    warn=None,
):
    """
    Collects unique studies in which each gene has a valid fold change.

    arguments:
        genes_identifiers (list<str>): unique identifiers of genes in all
            studies
        bin_studies (dict): collection of information about each study
        data_gene_annotation (object): Pandas data frame of genes' annotations
        warn (bool): whether to print warnings

    raises:

    returns:
        (object): Pandas data frame of studies for each gene

    """

    # Iterate across genes.
    records = list()
    for gene_identifier in genes_identifiers:
        # Determine whether gene has a valid identifier.
        if ("ENSG" in str(gene_identifier)):
            # Collect studies and comparisons for the gene.
            studies_gene = list()
            comparisons_gene = list()
            accumulations = list()
            depletions = list()
            for study in bin_studies.keys():
                data_study = bin_studies[study]["data"]
                # Determine whether the gene has valid fold change in study.
                if (gene_identifier in data_study["identifier"].to_list()):
                    # Study mentions current gene.
                    # Determine whether the gene has a valid fold change.
                    # Consider each study comparison for the gene.
                    comparisons_study = list(filter(
                        lambda value: (not value in ["identifier", "name"]),
                        data_study.columns.to_list()
                    ))
                    # Iterate across study's comparisons.
                    for comparison in comparisons_study:
                        value = determine_gene_study_comparison_value(
                            gene_identifier=gene_identifier,
                            comparison=comparison,
                            data_study=data_study,
                            warn=warn,
                        )
                        if not math.isnan(value):
                            studies_gene.append(study)
                            comparisons_gene.append(comparison)
                            if value >= 1:
                                accumulations.append(value)
                            elif value < 1:
                                depletions.append(value)
                            pass
                        pass
                    pass
                pass
            # Collect unique studies.
            studies_gene_unique = sorted(utility.collect_unique_elements(
                elements_original=studies_gene
            ))
            comparisons_gene_unique = sorted(utility.collect_unique_elements(
                elements_original=comparisons_gene
            ))
            # Organize record.
            record = dict()
            record["identifier"] = gene_identifier
            record["studies"] = len(studies_gene_unique)
            record["comparisons"] = len(comparisons_gene_unique)
            record["reference"] = ";".join(studies_gene_unique)
            record["accumulations"] = len(accumulations)
            record["depletions"] = len(depletions)
            annotations = assembly.access_gene_contextual_annotations(
                gene_identifier=gene_identifier,
                data_gene_annotation=data_gene_annotation,
            )
            record.update(annotations)
            records.append(record)
            pass
        pass
    # Organize data.
    data = pandas.DataFrame(data=records)
    data.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    return data


def collect_gene_studies_comparisons_folds(
    gene_identifier=None,
    entry_original=None,
    comparisons=None,
    bin_studies=None,
    warn=None,
):
    """
    Collects information about a gene's fold changes across multiple studies
    and comparisons.

    This function assumes that comparisons are not redundant across studies.
    Similar comparisons are distinct if they belong to different studies.

    arguments:
        gene_identifier (str): unique identifier of a gene
        entry_original (dict): entry of information about a gene
        comparisons (list<str>): unique designations of comparisons in all
            studies
        bin_studies (dict): collection of information about each study
        warn (bool): whether to print warnings

    raises:

    returns:
        (dict): novel entry for gene

    """

    # Copy original entry.
    entry = copy.deepcopy(entry_original)
    # Iterate across studies.
    for study in bin_studies.keys():
        data_study = bin_studies[study]["data"]
        # Determine whether study mentions the current gene.
        if (gene_identifier in data_study["identifier"].to_list()):
            # Consider each study comparison for the gene.
            comparisons_study = list(filter(
                lambda value: (not value in ["identifier", "name"]),
                data_study.columns.to_list()
            ))
            # Iterate across study's comparisons.
            for comparison in comparisons_study:
                # Determine whether gene's entry already has a value for the
                # comparison.
                if (comparison not in entry.keys()):
                    # Collect comparison's fold change value for the gene.
                    value = determine_gene_study_comparison_value(
                        gene_identifier=gene_identifier,
                        comparison=comparison,
                        data_study=data_study,
                    )
                    entry[comparison] = value
                else:
                    # Determine whether gene's entry has a null value for the
                    # comparison.
                    if math.isnan(entry[comparison]):
                        # Gene does not have a valid value for the study and
                        # comparison.
                        value = determine_gene_study_comparison_value(
                            gene_identifier=gene_identifier,
                            comparison=comparison,
                            data_study=data_study,
                            warn=warn,
                        )
                        entry[comparison] = value
                    else:
                        # Gene already has a valid value for the study and
                        # comparison.
                        value = determine_gene_study_comparison_value(
                            gene_identifier=gene_identifier,
                            comparison=comparison,
                            data_study=data_study,
                            warn=warn,
                        )
                        entry[comparison] = statistics.mean(
                            entry[comparison], value
                        )
                        # Warn.
                        if warn:
                            utility.print_terminal_partition(level=3)
                            print(
                                "warning: encountered multiple values for a " +
                                "gene, study, and comparison."
                            )
                    pass
                pass
            pass
        pass
    # Collect null values for any comparisons that do not apply to the gene.
    for comparison in comparisons:
        if comparison not in entry.keys():
            entry[comparison] = float("nan")
            pass
        pass
    # Return novel entry.
    return entry


def integrate_genes_studies_comparisons(
    genes_identifiers=None,
    comparisons=None,
    bin_studies=None,
    data_gene_annotation=None,
    warn=None,
):
    """
    Collects and organizes genes that show differential expression in multiple
    studies and comparisons.

    arguments:
        genes_identifiers (list<str>): unique identifiers of genes in all
            studies
        comparisons (list<str>): unique designations of comparisons in all
            studies
        bin_studies (dict): collection of information about each study
        data_gene_annotation (object): Pandas data frame of genes' annotations
        warn (bool): whether to print warnings

    raises:

    returns:
        (object): Pandas data frame of genes summaries across comparisons and
            studies

    """

    # Iterate across genes.
    entries = dict()
    for gene_identifier in genes_identifiers:
        # Determine whether gene has a valid identifier.
        if ("ENSG" in str(gene_identifier)):
            # Determine whether an entry already exists for the gene.
            # Multiple gene names might translate to a single identifier.
            if gene_identifier not in entries.keys():
                # An entry does not already exist for the gene.
                # Create new entry for gene.
                entry = dict()
                entry["identifier"] = gene_identifier
                # Collect information for gene.
                entry_novel = collect_gene_studies_comparisons_folds(
                    gene_identifier=gene_identifier,
                    entry_original=entry,
                    comparisons=comparisons,
                    bin_studies=bin_studies,
                    warn=warn,
                )
            else:
                # Collect information for gene.
                entry_novel = collect_gene_studies_comparisons_folds(
                    gene_identifier=gene_identifier,
                    entry_original=entries[identifier],
                    comparisons=comparisons,
                    bin_studies=bin_studies,
                    warn=warn,
                )
                pass
            # Collect entry.
            entries_novel = dict()
            entries_novel[entry_novel["identifier"]] = entry_novel
            entries.update(entries_novel)
            pass
        pass
    # Organize data.
    records = entries.values()
    data = pandas.DataFrame(data=records)
    return data


def calculate_mean_positive_negative_log_folds(
    data_studies_comparisons_folds=None,
    report=None,
):
    """
    Calculates the mean positive and mean negative logarithmic fold changes
    separately for each gene across studies and comparisons.

    arguments:
        data_studies_comparisons_folds (object): Pandas data frame of genes'
            logarithmic fold changes across studies and comparisons
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of genes' mean positive and negative
            logarithmic fold changes

    """

    def calculate_directional_mean(series=None, direction=None,):
        array = series.to_numpy()
        if direction == "positive":
            array_direction = numpy.where(array > 0, array, numpy.nan)
        elif direction == "negative":
            array_direction = numpy.where(array < 0, array, numpy.nan)
        array_valid = array_direction[numpy.where(
            ~numpy.isnan(array_direction)
        )]
        if array_valid.shape[0] > 0:
            mean = numpy.nanmean(array_valid)
        else:
            mean = numpy.nan
        return mean

    # Copy data.
    data = data_studies_comparisons_folds.copy(deep=True)
    # Calculate mean positive and negative logarithmic fold changes separately.
    data["log2_fold_positive"] = data.apply(
        lambda row: calculate_directional_mean(
            series=row,
            direction="positive",
        ), axis="columns", # Apply function to each row of data.
    )
    data["log2_fold_negative"] = data.apply(
        lambda row: calculate_directional_mean(
            series=row,
            direction="negative",
        ),
        axis="columns", # Apply function to each row of data.
    )
    # Organize data.
    data = data.loc[
        :, data.columns.isin(["log2_fold_positive", "log2_fold_negative"])
    ]
    # Return information.
    return data


def integrate_organize_genes_studies_comparisons(
    genes_identifiers=None,
    comparisons=None,
    bin_studies=None,
    data_genes_studies_valid=None,
    data_gene_annotation=None,
    report=None,
    warn=None,
):
    """
    Collects and organizes genes that show differential expression in multiple
    studies and comparisons.

    arguments:
        genes_identifiers (list<str>): unique identifiers of genes in all
            studies
        comparisons (list<str>): unique designations of comparisons in all
            studies
        bin_studies (dict): collection of information about each study
        data_genes_studies_valid (object): Pandas data frame of studies for
            each gene
        data_gene_annotation (object): Pandas data frame of genes' annotations
        report (bool): whether to print reports
        warn (bool): whether to print warnings

    raises:

    returns:
        (object): Pandas data frame of genes summaries across comparisons and
            studies

    """

    # Collect genes' summaries across studies and comparisons.
    data = integrate_genes_studies_comparisons(
        genes_identifiers=genes_identifiers,
        comparisons=comparisons,
        bin_studies=bin_studies,
        data_gene_annotation=data_gene_annotation,
        warn=warn,
    )
    data.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Calculate logarithmic fold changes.
    # Transform genes' signals to logarithmic space.
    # Transform values of genes' signals to base-2 logarithmic space.
    data_log = utility.calculate_pseudo_logarithm_signals(
        pseudo_count=0.0,
        base=2.0,
        data=data,
    )
    # Calculate mean values of logarithmic fold changes across studies and
    # comparisons.
    # Keep positive and negative logarithmic fold changes separate.
    data_directional_aggregate = calculate_mean_positive_negative_log_folds(
        data_studies_comparisons_folds=data_log,
        report=True,
    )
    # Combine data.
    # Logarithmic fold changes across studies and comparisons.
    # Aggregate positive and negative logarithmic fold changes.
    # Genes' contextual annotations and references to studies.
    data_folds_aggregates = data_log.join(
        data_directional_aggregate,
        how="left",
        on="identifier"
    )
    data_annotation = data_folds_aggregates.join(
        data_genes_studies_valid,
        how="left",
        on="identifier"
    )
    # Introduce counts of studies and comparisons for each gene.
    data_annotation = data_annotation[[
        #"identifier",
        "name",
        "chromosome",
        "start",
        "end",
        "strand",
        "studies",
        "comparisons",
        "accumulations",
        "depletions",
        "log2_fold_positive",
        "log2_fold_negative",
        *comparisons,
        "reference",
    ]]
    data_annotation.sort_values(
        by=["studies", "comparisons"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    return data_annotation


def collect_integrate_organize_genes_studies_comparisons(
    bin_studies=None,
    data_gene_annotation=None,
    translations_genes=None,
    report=None,
):
    """
    Collects and organizes genes that show differential expression in multiple
    studies and comparisons.

    arguments:
        bin_studies (dict): collection of information about each study
        data_gene_annotation (object): Pandas data frame of genes' annotations
        translations_genes (dict<str>): pairwise custom translations of genes'
            names to Ensembl identifiers, see
            assembly.read_source_gene_name_identifier_translations()
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of genes summaries across comparisons and
            studies

    """

    # Translate genes' names to Ensembl identifiers.
    bin_studies_genes = translate_study_genes_identifiers(
        bin_studies=bin_studies,
        data_gene_annotation=data_gene_annotation,
        translations_genes=translations_genes,
        report=False,
    )
    # Collect unique genes' identifiers across all studies.
    genes_identifiers = collect_studies_unique_gene_identifiers(
        bin_studies=bin_studies_genes,
    )
    # Report.
    if report:
        print("unique genes from all studies: " + str(len(genes_identifiers)))
        pass

    # Translate comparison identifiers from all studies.
    bin_studies_comparisons = translate_study_comparisons_identifiers(
        bin_studies=bin_studies_genes,
        report=False,
    )
    # Collect unique comparisons across all studies.
    comparisons = collect_studies_unique_comparisons(
        bin_studies=bin_studies_comparisons,
    )
    # Report.
    if report:
        print(
            "unique comparisons from all studies: " +
            str(len(comparisons))
        )
        pass
    # Collect studies that report valid fold changes for each gene.
    data_genes_studies_valid = (
        collect_genes_annotations_studies_comparisons_valid_change(
            genes_identifiers=genes_identifiers,
            bin_studies=bin_studies_comparisons,
            data_gene_annotation=data_gene_annotation,
            warn=True,
    ))
    # Collect and integrate genes' fold changes across comparisons in all
    # studies.
    data_genes_comparisons_studies = (
        integrate_organize_genes_studies_comparisons(
            genes_identifiers=genes_identifiers,
            comparisons=comparisons,
            bin_studies=bin_studies_comparisons,
            data_genes_studies_valid=data_genes_studies_valid,
            data_gene_annotation=data_gene_annotation,
            report=report,
            warn=True,
    ))
    # Return information.
    return data_genes_comparisons_studies


def select_covid19_genes_by_studies_fold_directions(
    data_genes_comparisons_studies=None,
    genes_selection=None,
    threshold_studies=None,
    report=None,
):
    """
    Collects and organizes genes that show differential expression in multiple
    studies and comparisons.

    arguments:
        data_genes_comparisons_studies (object): Pandas data frame of genes'
            differential expression in studies
        genes_selection (list<str>): identifiers of genes
        threshold_studies (int): minimal count of studies
        report (bool): whether to print reports

    raises:

    returns:
        (dict<list<str>>): sets of genes

    """

    # Copy data.
    data = data_genes_comparisons_studies.copy(deep=True)
    # Select data for genes that match selection for study.
    data = data.loc[data.index.isin(genes_selection), :]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Count of differential expression genes that match selection" +
            "of genes for study."
        )
        print("selection genes DE COVID-19: " + str(data.shape[0]))
    # Select data for genes that match threshold.
    data_studies = data.loc[
        data["studies"] >= threshold_studies, :
    ]
    genes_studies = utility.collect_unique_elements(
        elements_original=data_studies.index.to_list()
    )
    # Select data for genes that show accumulation in majority of studies.
    data_accumulation = data_studies.loc[
        data_studies["accumulations"] > (data_studies["depletions"] + 1), :
    ]
    genes_accumulation = utility.collect_unique_elements(
        elements_original=data_accumulation.index.to_list()
    )
    # Select data for genes that show depletion in majority of studies.
    data_depletion = data_studies.loc[
        data_studies["depletions"] > (data_studies["accumulations"] + 1), :
    ]
    genes_depletion = utility.collect_unique_elements(
        elements_original=data_depletion.index.to_list()
    )
    # Select data for genes that show depletion in majority of studies.
    data_mix = data_studies.loc[
        lambda datum:
            (datum["accumulations"] == datum["depletions"]) |
            (datum["accumulations"] == (datum["depletions"] + 1)) |
            (datum["depletions"] == (datum["accumulations"] + 1))
    ]
    genes_mix = utility.collect_unique_elements(
        elements_original=data_mix.index.to_list()
    )
    # Compile information.
    bin = dict()
    bin["studies"] = genes_studies
    bin["up"] = genes_accumulation
    bin["down"] = genes_depletion
    bin["mix"] = genes_mix
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Count of differential expression genes that match threshold."
        )
        print("threshold count of studies: " + str(threshold_studies))
        print("genes DE COVID-19 by studies: " + str(len(genes_studies)))
        print("genes by accumulation: " + str(len(genes_accumulation)))
        print("genes by depletion: " + str(len(genes_depletion)))
        print("genes by mix of folds: " + str(len(genes_mix)))
    # Return information.
    return bin


def write_product_covid_genes(
    paths=None,
    information=None,
):
    """
    Writes product information to file.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_data = os.path.join(
        paths["covid19"], "data_genes_comparisons_studies.pickle"
    )
    path_data_text = os.path.join(
        paths["covid19"], "data_genes_comparisons_studies.tsv"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_genes_comparisons_studies"],
        path_data
    )
    information["data_genes_comparisons_studies"].to_csv(
        path_or_buf=path_data_text,
        sep="\t",
        header=True,
        index=True,
    )

    # Specify directories and files.
    path_sets_genes = os.path.join(
        paths["covid19"], "sets_genes.pickle"
    )
    # Write information to file.
    with open(path_sets_genes, "wb") as file_product:
        pickle.dump(information["sets_genes"], file_product)
    # Write individual gene sets.
    for set in information["sets_genes"].keys():
        # Specify directories and files.
        path_set = os.path.join(paths["covid19"], str(set + ".txt"))
        # Write information to file.
        utility.write_file_text_list(
            elements=information["sets_genes"][set],
            delimiter="\n",
            path_file=path_set
        )
    pass


# Collection of cytokine genes




# Collection of integrin genes



###############################################################################
# Procedure


def read_collect_organize_report_write_covid19_genes(
    paths=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_covid19_genes(
        dock=paths["dock"],
    )

    # Collect and organize genes from each study and comparison.
    # Do not combine similar comparisons from different studies. For example,
    # do not combine information from "cd4_t_cell" in two different studies.
    # Such a combination of fold changes would complicate interpretation and
    # would be problematic as different studies use different tissues such as
    # PBMCs versus BALF.
    # Instead, organize comparisons that are distinct for each study.
    data_genes_comparisons_studies = (
        collect_integrate_organize_genes_studies_comparisons(
            bin_studies=source["bin_studies"],
            data_gene_annotation=source["data_gene_annotation"],
            translations_genes=source["translations_genes"],
            report=True,
    ))

    utility.print_terminal_partition(level=1)
    print(data_genes_comparisons_studies)

    # Select genes that show differential expression in multiple studies.
    sets_genes = select_covid19_genes_by_studies_fold_directions(
        data_genes_comparisons_studies=data_genes_comparisons_studies,
        genes_selection=source["genes_selection"],
        threshold_studies=1,
        report=True,
    )

    # Compile information.
    bin = dict()
    bin["data_genes_comparisons_studies"] = data_genes_comparisons_studies
    bin["sets_genes"] = sets_genes
    write_product_covid_genes(
        information=bin,
        paths=paths,
    )
    pass


def execute_procedure(
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Initialize directories.
    paths = initialize_directories(dock=dock)
    # Organize genes with differential expression in COVID-19.
    read_collect_organize_report_write_covid19_genes(
        paths=paths,
    )

    pass


if (__name__ == "__main__"):
    execute_procedure()
