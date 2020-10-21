"""
...

"""

###############################################################################
# Notes

# If the count of bootstrap shuffles is too small, then it is likely not to
# have any values equal to or greater than the real value, in which case the
# probability is zero.

###############################################################################
# Installation and importation

# Standard

import os
import math
import statistics
import pickle
import copy
import random
import itertools

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import assembly
import distribution
import prediction
#import plot
import utility

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
    paths["integration"] = os.path.join(paths["dock"], "integration")
    paths["groups"] = os.path.join(paths["dock"], "integration", "groups")
    paths["gene_set_summaries"] = os.path.join(
        paths["dock"], "integration", "gene_set_summaries"
    )
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=paths["integration"])
    utility.create_directory(path=paths["integration"])
    utility.create_directories(path=paths["groups"])
    utility.create_directories(path=paths["gene_set_summaries"])
    # Define paths.
    cohorts = list()
    cohorts.append("selection")
    cohorts.append("respiration")
    cohorts.append("ventilation")
    for cohort in cohorts:
        paths[cohort] = dict()
        paths[cohort]["set"] = dict()
        paths[cohort]["set"]["cardinality"] = os.path.join(
            paths["integration"], cohort, "set", "cardinality"
        )
        paths[cohort]["set"]["collection_candidacy"] = os.path.join(
            paths["integration"], cohort, "set", "collection_candidacy"
        )
        paths[cohort]["set"]["allocation"] = os.path.join(
            paths["integration"], cohort, "set", "allocation"
        )
        paths[cohort]["set"]["prediction_ontology"] = os.path.join(
            paths["integration"], cohort, "set", "prediction_ontology"
        )
        paths[cohort]["set"]["prediction_interaction"] = dict()
        paths[cohort]["set"]["prediction_interaction"]["any"] = os.path.join(
            paths["integration"], cohort, "set", "prediction_interaction",
            "any"
        )
        paths[cohort]["set"]["prediction_interaction"]["multimodal"] = (
            os.path.join(
                paths["integration"], cohort, "set", "prediction_interaction",
                "multimodal"
        ))
        paths[cohort]["population"] = os.path.join(
            paths["integration"], cohort, "population"
        )
        paths[cohort]["correlation"] = os.path.join(
            paths["integration"], cohort, "correlation"
        )
        # Initialize directories.
        utility.create_directories(path=paths[cohort]["set"]["cardinality"])
        utility.create_directories(
            path=paths[cohort]["set"]["collection_candidacy"]
        )
        utility.create_directories(path=paths[cohort]["set"]["allocation"])
        utility.create_directories(
            path=paths[cohort]["set"]["prediction_ontology"]
        )
        utility.create_directories(
            path=paths[cohort]["set"]["prediction_interaction"]["any"]
        )
        utility.create_directories(
            path=paths[cohort]["set"]["prediction_interaction"]["multimodal"]
        )
        utility.create_directories(path=paths[cohort]["population"])
        utility.create_directories(path=paths[cohort]["correlation"])
    # Return information.
    return paths


##########
# Read source


def read_source_annotation_sets(dock=None):
    """
    Reads and organizes source information from file

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    def match_format(name=None):
        return (".csv" in name)

    # Specify directories and files.
    # Read information from file.

    # Ontology.
    path_annotation = os.path.join(dock, "annotation")
    path_function = os.path.join(path_annotation, "function")
    path_structure = os.path.join(path_annotation, "structure")
    path_hypothesis = os.path.join(path_annotation, "hypothesis")

    # Select files.
    files_function = os.listdir(path=path_function)
    files_function_format = list(filter(match_format, files_function))
    files_structure = os.listdir(path=path_structure)
    files_structure_format = list(filter(match_format, files_structure))
    files_hypothesis = os.listdir(path=path_hypothesis)
    files_hypothesis_format = list(filter(match_format, files_hypothesis))

    # Collect information from files.
    collection_function = dict()
    for file in files_function_format:
        # Read information.
        path_data = os.path.join(
            path_function, file
        )
        data = pandas.read_csv(
            path_data,
            sep="\t",
            header=0,
        )
        # Extract name of set.
        name = file.replace(".csv", "")
        # Collect information.
        collection_function[name] = data
        pass

    # Collect information from files.
    collection_structure = dict()
    for file in files_structure_format:
        # Read information.
        path_data = os.path.join(
            path_structure, file
        )
        data = pandas.read_csv(
            path_data,
            sep="\t",
            header=0,
        )
        # Extract name of set.
        name = file.replace(".csv", "")
        # Collect information.
        collection_structure[name] = data
        pass

    # Collect information from files.
    collection_hypothesis = dict()
    for file in files_hypothesis_format:
        # Read information.
        path_data = os.path.join(
            path_hypothesis, file
        )
        data = pandas.read_csv(
            path_data,
            sep="\t",
            header=0,
        )
        # Extract name of set.
        name = file.replace(".csv", "")
        # Collect information.
        collection_hypothesis[name] = data
        pass

    # Compile and return information.
    return {
        "function": collection_function,
        "structure": collection_structure,
        "hypothesis": collection_hypothesis,
    }


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
    # Read information from file.

    # Annotation.
    sets_query = read_source_annotation_query_genes_sets(dock=dock)
    genes_query = read_source_annotation_query_genes_all(
        dock=dock,
    )

    # Selection.
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )
    path_persons_properties = os.path.join(
        path_selection, "data_persons_properties.pickle"
    )
    path_genes_selection = os.path.join(
        path_selection, "genes_selection.pickle"
    )
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)

    # Distribution.
    path_distribution = os.path.join(dock, "distribution")
    path_collection = os.path.join(path_distribution, "collection")
    path_distribution_report = os.path.join(
        path_collection, "data_gene_report.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        path_collection, "data_signals_genes_persons.pickle"
    )
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    data_gene_distribution_report = pandas.read_pickle(
        path_distribution_report
    )

    # Probability.
    path_probability = os.path.join(dock, "probability")
    path_data_genes_permutation_probabilities = os.path.join(
        path_probability, "data_genes_discoveries.pickle"
    )
    data_genes_permutation_probabilities = pandas.read_pickle(
        path_data_genes_permutation_probabilities
    )
    path_genes_probability = os.path.join(
        path_probability, "genes_probability.pickle"
    )
    with open(path_genes_probability, "rb") as file_source:
        genes_probability = pickle.load(file_source)
    path_sets_genes_multimodal_permutation = os.path.join(
        path_probability, "sets_genes_multimodal.pickle"
    )
    with open(path_sets_genes_multimodal_permutation, "rb") as file_source:
        sets_genes_multimodal_permutation = pickle.load(file_source)

    # Heritability.
    path_heritability = os.path.join(dock, "heritability")
    path_collection = os.path.join(path_heritability, "collection")
    path_data_genes_heritabilities_complex = os.path.join(
        path_collection, "data_genes_heritabilities_complex.pickle"
    )
    data_genes_heritabilities_complex = pandas.read_pickle(
        path_data_genes_heritabilities_complex
    )

    # Candidacy.
    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_unimodal = os.path.join(
        path_candidacy, "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)

    # Prediction.
    path_prediction = os.path.join(dock, "prediction")
    path_data_regression_genes = os.path.join(
        path_prediction, "data_regression_genes.pickle"
    )
    path_sets_genes_prediction = os.path.join(
        path_prediction, "sets.pickle"
    )
    path_genes_prediction = os.path.join(
        path_prediction, "genes_prediction.pickle"
    )
    path_genes_multimodal_prediction = os.path.join(
        path_prediction, "genes_multimodal_prediction.pickle"
    )
    data_regression_genes = pandas.read_pickle(
        path_data_regression_genes
    )
    with open(path_sets_genes_prediction, "rb") as file_source:
        sets_genes_prediction = pickle.load(file_source)
    with open(path_genes_prediction, "rb") as file_source:
        genes_prediction = pickle.load(file_source)
    with open(path_genes_multimodal_prediction, "rb") as file_source:
        genes_multimodal_prediction = pickle.load(file_source)


    # Compile and return information.
    return {
        "sets_query": sets_query,
        "genes_query": genes_query,

        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "genes_selection": genes_selection,

        "data_signals_genes_persons": data_signals_genes_persons,
        "data_gene_distribution_report": data_gene_distribution_report,

        "data_genes_permutation_probabilities": (
            data_genes_permutation_probabilities
        ),
        "genes_probability": genes_probability,
        "sets_genes_multimodal_permutation": sets_genes_multimodal_permutation,

        "data_genes_heritabilities_complex": data_genes_heritabilities_complex,

        "genes_unimodal": genes_unimodal,
        "genes_multimodal": genes_multimodal,

        "data_regression_genes": data_regression_genes,
        "sets_genes_prediction": sets_genes_prediction,
        "genes_prediction": genes_prediction,
        "genes_multimodal_prediction": genes_multimodal_prediction,
    }


##########
# Gene sets


# TODO: DAVID gene ontologies need to be updated for new prediction sets...


def read_source_genes_sets_collection(
    count=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        count (int): threshold by count of studies
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    if count == 1:
        set = "study_one"
    elif count == 2:
        set = "study_two"
    elif count == 3:
        set = "study_three"
        pass
    # Specify directories and files.
    path_sets_genes_covid19 = os.path.join(
        dock, "collection", "covid19", set, "sets_genes.pickle"
    )
    # Read information from file.
    with open(path_sets_genes_covid19, "rb") as file_source:
        sets_genes_covid19 = pickle.load(file_source)
    # Compile information.
    bin = dict()
    bin["covid19"] = sets_genes_covid19
    # Return information.
    return bin


def read_source_genes_sets_candidacy(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_genes_any = os.path.join(
        dock, "candidacy", cohort, "distribution", "any", "genes.pickle"
    )
    path_genes_multimodal = os.path.join(
        dock, "candidacy", cohort, "distribution", "multimodal", "genes.pickle"
    )
    path_genes_nonmultimodal = os.path.join(
        dock, "candidacy", cohort, "distribution", "nonmultimodal",
        "genes.pickle"
    )
    path_genes_unimodal = os.path.join(
        dock, "candidacy", cohort, "distribution", "unimodal", "genes.pickle"
    )
    # Read information from file.
    with open(path_genes_any, "rb") as file_source:
        genes_any = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)
    with open(path_genes_nonmultimodal, "rb") as file_source:
        genes_nonmultimodal = pickle.load(file_source)
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    # Compile information.
    bin = dict()
    bin["any"] = genes_any
    bin["multimodal"] = genes_multimodal
    bin["nonmultimodal"] = genes_nonmultimodal
    bin["unimodal"] = genes_unimodal
    # Return information.
    return bin


def read_source_genes_sets_heritability(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (list<str>): identifiers of genes

    """

    # Specify directories and files.
    path_genes_validity = os.path.join(
        dock, "heritability", cohort, "collection", "genes_validity.pickle"
    )
    path_genes_heritability = os.path.join(
        dock, "heritability", cohort, "collection", "genes.pickle"
    )
    # Read information from file.
    with open(path_genes_validity, "rb") as file_source:
        genes_validity = pickle.load(file_source)
    with open(path_genes_heritability, "rb") as file_source:
        genes_heritability = pickle.load(file_source)
    # Compile information.
    bin = dict()
    bin["validity"] = genes_validity
    bin["heritability"] = genes_heritability
    # Return information.
    return bin


def read_source_genes_sets_prediction(
    cohort=None,
    model=None,
    genes_query_scope=None,
    dock=None,
):
    """
    Reads and organizes source information from file.

    Priority parameters:
    cohort="selection"
    model="selection_main"
    genes_query_scope="covid19_multimodal"

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): name of regression model
        genes_query_scope (str): the set of genes by which to query
            associations, thereby reducing the scope of multiple hypotheses
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_sets_genes = os.path.join(
        dock, "prediction", cohort, model, "genes_associations",
        "sets_genes.pickle"
    )
    # Read information from file.
    with open(path_sets_genes, "rb") as file_source:
        sets_genes = pickle.load(file_source)
    # Compile information.
    sets = copy.deepcopy(sets_genes[genes_query_scope])
    # Return information.
    return sets


def read_source_annotation_query_genes_all(
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        set (str): name of set in query table
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (list<str>): identifiers of genes

    """

    # Specify directories and files.
    # Read information from file.

    # Annotation.
    path_data_genes_query = os.path.join(
        dock, "annotation", "gene_sets", "genes_query", "query_genes.tsv"
    )
    data_genes_query = pandas.read_csv(
        path_data_genes_query,
        sep="\t",
        header=0,
    )
    genes_unique = utility.collect_unique_elements(
        elements_original=data_genes_query["identifier"].to_list()
    )
    # Return information.
    return genes_unique


def read_source_annotation_query_genes_set(
    set=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        set (str): name of set in query table
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (list<str>): identifiers of genes

    """

    # Specify directories and files.
    # Read information from file.

    # Annotation.
    path_data_genes_query = os.path.join(
        dock, "annotation", "gene_sets", "genes_query", "query_genes.tsv"
    )
    data_genes_query = pandas.read_csv(
        path_data_genes_query,
        sep="\t",
        header=0,
    )
    data_genes_query_pass = data_genes_query.loc[
        data_genes_query[set] == 1, :
    ]
    genes = data_genes_query_pass["identifier"].to_list()
    # Return information.
    return genes


def read_source_annotation_query_genes_sets(dock=None):
    """
    Reads and organizes source information from file

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    feature = read_source_annotation_query_genes_set(
        set="feature",
        dock=dock,
    )
    heritability = read_source_annotation_query_genes_set(
        set="heritability",
        dock=dock,
    )
    distribution = read_source_annotation_query_genes_set(
        set="distribution",
        dock=dock,
    )
    correlation = read_source_annotation_query_genes_set(
        set="correlation",
        dock=dock,
    )
    covid19 = read_source_annotation_query_genes_set(
        set="covid19",
        dock=dock,
    )
    covid19_up_prediction = read_source_annotation_query_genes_set(
        set="covid19_up_prediction",
        dock=dock,
    )
    covid19_drug = read_source_annotation_query_genes_set(
        set="covid19_drug",
        dock=dock,
    )
    # Compile and return information.
    return {
        "feature": feature,
        "heritability": heritability,
        "distribution": distribution,
        "correlation": correlation,
        "covid19": covid19,
        "covid19_up_prediction": covid19_up_prediction,
        "covid19_drug": covid19_drug,
    }


def read_source_genes_sets_correlation(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_genes_correlation = os.path.join(
        dock, "integration", cohort, "correlation", "genes.pickle"
    )
    # Read information from file.
    with open(path_genes_correlation, "rb") as file_source:
        genes_correlation = pickle.load(file_source)
    # Return information.
    return genes_correlation


def read_source_genes_sets_combination(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Compile information.
    bin = dict()
    # Read genes sets.
    genes_candidacy = read_source_genes_sets_candidacy(
        cohort=cohort,
        dock=dock,
    )
    genes_collection = read_source_genes_sets_collection(
        count=1,
        dock=dock,
    )
    #genes_heritability = read_source_genes_sets_heritability(
    #    cohort="selection",
    #    dock=dock,
    #)
    # Organize genes sets.
    bin["any"] = genes_candidacy["any"]
    bin["multimodal"] = genes_candidacy["multimodal"]
    bin["nonmultimodal"] = genes_candidacy["nonmultimodal"]
    bin["unimodal"] = genes_candidacy["unimodal"]
    bin["covid19"] = genes_collection["covid19"]["any"]
    bin["covid19_multimodal"] = utility.filter_common_elements(
        list_one=genes_collection["covid19"]["any"],
        list_two=genes_candidacy["multimodal"],
    )
    bin["covid19_up"] = genes_collection["covid19"]["up"]
    bin["covid19_down"] = genes_collection["covid19"]["down"]
    bin["covid19_mix"] = genes_collection["covid19"]["mix"]
    bin["covid19_up_multimodal"] = utility.filter_common_elements(
        list_one=genes_collection["covid19"]["up"],
        list_two=genes_candidacy["multimodal"],
    )
    bin["covid19_down_multimodal"] = utility.filter_common_elements(
        list_one=genes_collection["covid19"]["down"],
        list_two=genes_candidacy["multimodal"],
    )
    bin["covid19_mix_multimodal"] = utility.filter_common_elements(
        list_one=genes_collection["covid19"]["mix"],
        list_two=genes_candidacy["multimodal"],
    )
    #bin["heritability"] = genes_heritability["heritability"]
    # Return information.
    return bin


# TODO: include ontology gene sets... DAVID enrichment sets on multimodal genes and COVID-19 genes

def read_source_genes_sets_combination_integration(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Compile information.
    bin = dict()
    # Read genes sets.
    bin["collection"] = read_source_genes_sets_collection(
        count=1,
        dock=dock,
    )
    bin["candidacy"] = read_source_genes_sets_candidacy(
        cohort=cohort,
        dock=dock,
    )
    bin["heritability"] = read_source_genes_sets_heritability(
        cohort=cohort,
        dock=dock,
    )
    bin["combination"] = read_source_genes_sets_combination(
        cohort=cohort,
        dock=dock,
    )
    bin["prediction"] = read_source_genes_sets_prediction(
        cohort=cohort,
        model=str(cohort + "_main"),
        genes_query_scope="covid19_multimodal",
        dock=dock,
    )
    bin["query"] = read_source_annotation_query_genes_sets(
        dock=dock,
    )
    # Return information.
    return bin


def write_product_collection_candidacy_gene_sets(
    cohort=None,
    set=None,
    genes=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        set (str): name of set of genes
        genes (list<str>): identifiers of genes in set
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_genes = os.path.join(
        paths[cohort]["set"]["collection_candidacy"], str(set + ".txt")
    )
    # Write information to file.
    utility.write_file_text_list(
        elements=genes,
        delimiter="\n",
        path_file=path_genes
    )
    pass


def read_organize_report_write_collection_candidacy_gene_sets(
    paths=None
):
    """
    Organizes sets of genes by integration of information from multiple
    procedures.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_genes_sets_collection_candidacy(
        cohort="selection",
        dock=paths["dock"],
    )
    # Write sets of genes to file.
    for set in source.keys():
        # Report.
        #utility.print_terminal_partition(level=2)
        #print(set + " genes: " + str(len(source[set])))
        write_product_collection_candidacy_gene_sets(
            cohort="selection",
            set=set,
            genes=source[set],
            paths=paths,
        )
        pass
    pass


##########
# Contingency table comparisons


def compile_gene_set_attribution_table(
    genes_master=None,
    genes_sets=None,
    report=None,
):
    """
    Organizes a table of genes' attribution to sets.

    arguments:
        genes_master (list<str>): inclusive list of genes' identifiers
        genes_sets (dict): collections of genes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of genes' sets

    """

    records = list()
    for gene in genes_master:
        record = dict()
        record["identifier"] = gene
        for set in genes_sets.keys():
            record[set] = (gene in genes_sets[set])
            pass
        records.append(record)
        pass
    # Organize data.
    data = utility.convert_records_to_dataframe(records=records)
    data.set_index(
        "identifier",
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
        columns="sets",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(data)
    return data


def organize_contingency_table_chi(
    observations=None,
    variables_contingency=None,
    data=None,
    report=None,
):
    """
    Extracts identifiers of persons with valid genotypes.

    arguments:
        observations (list<str>): identifiers of rows to include
        variables_contingency (list<str>): names of variables for contingency
        data (object): Pandas data frame of observations in categories
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy data.
    data = data.copy(deep=True)
    # Organize data.
    data_observations = data.loc[data.index.isin(observations), :]
    data_categories = data_observations.loc[
        :, data_observations.columns.isin(variables_contingency)
    ]
    # Replace missing values with zero.
    data_categories.fillna(
        value=False,
        #axis="columns",
        inplace=True,
    )
    # Contingency table.
    data_contingency = pandas.crosstab(
        data_categories[variables_contingency[0]],
        data_categories[variables_contingency[1]],
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
        utility.print_terminal_partition(level=3)
        print(data_contingency)
        print(data_contingency.to_numpy())
        utility.print_terminal_partition(level=4)
        print("chi2: " + str(chi2))
        print("probability: " + str(probability))
    pass




def read_organize_report_write_contingency_table_comparisons(
    paths=None,
    report=None,
):
    """
    Organizes contingency table comparisons across genes and persons.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Read source information from file.
    sets_candidacy = read_source_genes_sets_candidacy(
        cohort="selection",
        dock=paths["dock"],
    )
    sets_covid19_one = read_source_genes_sets_collection(
        count=1,
        dock=paths["dock"],
    )
    sets_covid19_two = read_source_genes_sets_collection(
        count=2,
        dock=paths["dock"],
    )
    sets_covid19_three = read_source_genes_sets_collection(
        count=3,
        dock=paths["dock"],
    )
    sets_heritability = read_source_genes_sets_heritability(
        cohort="selection",
        dock=paths["dock"],
    )
    # Compile gene attribution table.
    sets = dict()
    sets["any"] = sets_candidacy["any"]
    sets["multimodal"] = sets_candidacy["multimodal"]
    sets["unimodal"] = sets_candidacy["unimodal"]
    sets["validity"] = sets_heritability["validity"]
    sets["heritability"] = sets_heritability["heritability"]
    sets["covid19"] = sets_covid19_one["covid19"]["any"]
    data_gene_sets = compile_gene_set_attribution_table(
        genes_master=sets["any"],
        genes_sets=sets,
        report=report,
    )
    # Aggregate a contingency table.
    organize_contingency_table_chi(
        observations=sets_heritability["validity"],
        variables_contingency=["multimodal", "heritability"],
        data=data_gene_sets,
        report=True,
    )
    organize_contingency_table_chi(
        observations=sets_candidacy["any"],
        variables_contingency=["multimodal", "covid19"],
        data=data_gene_sets,
        report=True,
    )

    pass




##########
# Comparisons of gene pan-tissue signals between groups of persons


def read_source_genes_signals_persons_properties(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", cohort,
        "data_persons_properties.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    # Return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "data_signals_genes_persons": data_signals_genes_persons,
    }


def combine_persons_properties_sets(
    sets_first=None,
    sets_second=None,
    bin=None,
):
    """
    Combines sets of persons.

    arguments:
        sets_first (list<str>): names of sets to place first in name
        sets_second (list<str>): names of sets to place second in name
        bin (dict): identifiers of persons in groups by their properties

    raises:

    returns:
        (dict): identifiers of persons in groups by their properties

    """

    # Iterate.
    for set_first in sets_first:
        for set_second in sets_second:
            set_name = str(set_first + "_" + set_second)
            bin[set_name] = utility.filter_common_elements(
                list_one=bin[set_first],
                list_two=bin[set_second],
            )
            pass
        pass
    return bin


def organize_persons_properties_sets(
    data_persons_properties=None,
    report=None,
):
    """
    Extracts identifiers of persons.

    arguments:
        data_persons_properties (object): Pandas data frame of persons'
            properties
        report (bool): whether to print reports

    raises:

    returns:
        (dict<list<str>>): identifiers of persons in groups by their properties

    """

    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    # Organize data.
    bin = dict()
    # Sex.
    bin["male"] = data_persons_properties.loc[
        data_persons_properties["sex_text"] == "male", :
    ].index.to_list()
    bin["female"] = data_persons_properties.loc[
        data_persons_properties["sex_text"] == "female", :
    ].index.to_list()
    # Age.
    bin["young"] = data_persons_properties.loc[
        data_persons_properties["age_grade"] == 0, :
    ].index.to_list()
    bin["old"] = data_persons_properties.loc[
        data_persons_properties["age_grade"] == 2, :
    ].index.to_list()
    bin["younger"] = data_persons_properties.loc[
        data_persons_properties["age_halves"] == 0, :
    ].index.to_list()
    bin["older"] = data_persons_properties.loc[
        data_persons_properties["age_halves"] == 1, :
    ].index.to_list()
    # Race.
    bin["race_europe"] = data_persons_properties.loc[
        data_persons_properties["race"] == "europe", :
    ].index.to_list()
    bin["race_africa"] = data_persons_properties.loc[
        data_persons_properties["race"] == "africa", :
    ].index.to_list()
    bin["race_asia"] = data_persons_properties.loc[
        data_persons_properties["race"] == "asia", :
    ].index.to_list()
    bin["race_america"] = data_persons_properties.loc[
        data_persons_properties["race"] == "america", :
    ].index.to_list()
    bin["race_other"] = data_persons_properties.loc[
        data_persons_properties["race"] == "other", :
    ].index.to_list()
    bin["race_white"] = data_persons_properties.loc[
        data_persons_properties["race_white"] == 1, :
    ].index.to_list()
    bin["race_not_white"] = data_persons_properties.loc[
        data_persons_properties["race_white"] == 0, :
    ].index.to_list()
    # Ventilation.
    bin["breath"] = data_persons_properties.loc[
        data_persons_properties["ventilation"] == False, :
    ].index.to_list()
    bin["ventilation"] = data_persons_properties.loc[
        data_persons_properties["ventilation"] == True, :
    ].index.to_list()
    # Leukocytosis.
    bin["calm"] = data_persons_properties.loc[
        data_persons_properties["leukocyte"] == False, :
    ].index.to_list()
    bin["leukocytosis"] = data_persons_properties.loc[
        data_persons_properties["leukocyte"] == True, :
    ].index.to_list()
    # Chronic inflammation.
    bin["calm_chronic"] = data_persons_properties.loc[
        data_persons_properties["inflammation"] == False, :
    ].index.to_list()
    bin["inflammation"] = data_persons_properties.loc[
        data_persons_properties["inflammation"] == True, :
    ].index.to_list()
    # Combinations.
    bin = combine_persons_properties_sets(
        sets_first=[
            "male", "female", "young", "younger", "old", "older",
            "race_white", "race_not_white"
        ],
        sets_second=["breath", "ventilation", "calm", "leukocytosis"],
        bin=bin,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Count of persons in each group...")
        utility.print_terminal_partition(level=2)
        for group in bin.keys():
            #utility.print_terminal_partition(level=4)
            print(group + " persons: " + str(len(bin[group])))
            pass
        utility.print_terminal_partition(level=2)
        pass
    # Return information.
    return bin


def compare_genes_signals_persons_groups_two(
    gene_identifier=None,
    comparison=None,
    group_1_persons=None,
    group_2_persons=None,
    group_1_label=None,
    group_2_label=None,
    data_signals_genes_persons=None,
    data_gene_annotation=None,
    report=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        gene_identifier (str): identifier of a gene
        comparison (str): name for comparison
        group_1_persons (list<str>): identifiers of persons in first group
        group_2_persons (list<str>): identifiers of persons in second group
        group_1_label (str): name of first group of persons
        group_2_label (str): name of first group of persons
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        data_gene_annotation (object): Pandas data frame of genes' annotations
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information about comparison

    """

    # Compile information.
    bin = dict()
    bin["gene_identifier"] = gene_identifier
    bin["gene_name"] = assembly.access_gene_name(
        identifier=gene_identifier,
        data_gene_annotation=data_gene_annotation,
    )
    bin["comparison"] = comparison
    bin["title"] = str(bin["gene_name"] + "_" + comparison)
    bin["group_1_persons"] = len(group_1_persons)
    bin["group_2_persons"] = len(group_2_persons)
    bin["group_1_label"] = group_1_label
    bin["group_2_label"] = group_2_label
    # Copy data.
    data_signals = data_signals_genes_persons.copy(deep=True)
    # Select data for persons in groups.
    data_group_1 = data_signals.loc[
        data_signals.index.isin(group_1_persons), :
    ]
    data_group_2 = data_signals.loc[
        data_signals.index.isin(group_2_persons), :
    ]
    # Select gene's signals for each group.
    bin["group_1_values"] = data_group_1[gene_identifier].dropna().to_numpy()
    bin["group_2_values"] = data_group_2[gene_identifier].dropna().to_numpy()
    bin["group_1_valids"] = bin["group_1_values"].size
    bin["group_2_valids"] = bin["group_2_values"].size
    # Calculate probability by t test.
    t_statistic, p_value = scipy.stats.ttest_ind(
        bin["group_1_values"],
        bin["group_2_values"],
        equal_var=True,
        nan_policy="omit",
    )
    bin["probability"] = p_value
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("gene: " + bin["gene_name"])
        print("comparison: " + bin["comparison"])
        print(
            bin["group_1_label"] + " (" + str(bin["group_1_valids"]) + ")" +
            " versus " +
            bin["group_2_label"] + " (" + str(bin["group_2_valids"]) + ")"
        )
        print("t test p-value: " + str(bin["probability"]))
    # Return information.
    return bin


def compare_genes_signals_persons_groups_four(
    gene_identifier=None,
    comparison=None,
    group_1_persons=None,
    group_2_persons=None,
    group_3_persons=None,
    group_4_persons=None,
    group_1_label=None,
    group_2_label=None,
    group_3_label=None,
    group_4_label=None,
    data_signals_genes_persons=None,
    data_gene_annotation=None,
    report=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        gene_identifier (str): identifier of a gene
        comparison (str): name for comparison
        group_1_persons (list<str>): identifiers of persons in first group
        group_2_persons (list<str>): identifiers of persons in second group
        group_3_persons (list<str>): identifiers of persons in third group
        group_4_persons (list<str>): identifiers of persons in fourth group
        group_1_label (str): name of first group of persons
        group_2_label (str): name of second group of persons
        group_3_label (str): name of third group of persons
        group_4_label (str): name of fourth group of persons
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        data_gene_annotation (object): Pandas data frame of genes' annotations
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information about comparison

    """

    # Compile information.
    bin = dict()
    bin["gene_identifier"] = gene_identifier
    bin["gene_name"] = assembly.access_gene_name(
        identifier=gene_identifier,
        data_gene_annotation=data_gene_annotation,
    )
    bin["comparison"] = comparison
    bin["title"] = str(bin["gene_name"] + "_" + comparison)
    bin["group_1_persons"] = len(group_1_persons)
    bin["group_2_persons"] = len(group_2_persons)
    bin["group_3_persons"] = len(group_3_persons)
    bin["group_4_persons"] = len(group_4_persons)
    bin["group_1_label"] = group_1_label
    bin["group_2_label"] = group_2_label
    bin["group_3_label"] = group_3_label
    bin["group_4_label"] = group_4_label
    # Copy data.
    data_signals = data_signals_genes_persons.copy(deep=True)
    # Select data for persons in groups.
    data_group_1 = data_signals.loc[
        data_signals.index.isin(group_1_persons), :
    ]
    data_group_2 = data_signals.loc[
        data_signals.index.isin(group_2_persons), :
    ]
    data_group_3 = data_signals.loc[
        data_signals.index.isin(group_3_persons), :
    ]
    data_group_4 = data_signals.loc[
        data_signals.index.isin(group_4_persons), :
    ]
    # Select gene's signals for each group.
    bin["group_1_values"] = data_group_1[gene_identifier].dropna().to_numpy()
    bin["group_2_values"] = data_group_2[gene_identifier].dropna().to_numpy()
    bin["group_3_values"] = data_group_3[gene_identifier].dropna().to_numpy()
    bin["group_4_values"] = data_group_4[gene_identifier].dropna().to_numpy()
    bin["group_1_valids"] = bin["group_1_values"].size
    bin["group_2_valids"] = bin["group_2_values"].size
    bin["group_3_valids"] = bin["group_3_values"].size
    bin["group_4_valids"] = bin["group_4_values"].size
    # Calculate probability by t test.
    t_statistic, p_value = scipy.stats.ttest_ind(
        bin["group_1_values"],
        bin["group_2_values"],
        equal_var=True,
        nan_policy="omit",
    )
    bin["probability_1_2"] = p_value
    t_statistic, p_value = scipy.stats.ttest_ind(
        bin["group_3_values"],
        bin["group_4_values"],
        equal_var=True,
        nan_policy="omit",
    )
    bin["probability_3_4"] = p_value
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("gene: " + bin["gene_name"])
        print("comparison: " + bin["comparison"])
        print(
            bin["group_1_label"] + " (" + str(bin["group_1_valids"]) + ")" +
            " versus " +
            bin["group_2_label"] + " (" + str(bin["group_2_valids"]) + ")"
        )
        print("t test p-value: " + str(bin["probability_1_2"]))
        utility.print_terminal_partition(level=4)
        print(
            bin["group_3_label"] + " (" + str(bin["group_3_valids"]) + ")" +
            " versus " +
            bin["group_4_label"] + " (" + str(bin["group_4_valids"]) + ")"
        )
        print("t test p-value: " + str(bin["probability_3_4"]))
    # Return information.
    return bin


def organize_genes_person_two_groups_signal_comparisons(
    genes=None,
    sets_persons=None,
    data_signals_genes_persons=None,
    data_gene_annotation=None,
    report=None,
):
    """
    Organizes and comparisons between two groups across multiple genes.

    arguments:
        genes (list<str>): identifiers of genes
        sets_persons (dict<list<str>>): identifiers of persons in groups by
            their properties
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        data_gene_annotation (object): Pandas data frame of genes' annotations
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): information about comparisons

    """

    comparisons = list()
    for gene in genes:
        comparisons.append(compare_genes_signals_persons_groups_two(
            gene_identifier=gene,
            comparison="sex",
            group_1_persons=sets_persons["female"],
            group_2_persons=sets_persons["male"],
            group_1_label="female",
            group_2_label="male",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_two(
            gene_identifier=gene,
            comparison="age",
            group_1_persons=sets_persons["young"],
            group_2_persons=sets_persons["old"],
            group_1_label="young",
            group_2_label="old",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_two(
            gene_identifier=gene,
            comparison="ager",
            group_1_persons=sets_persons["younger"],
            group_2_persons=sets_persons["older"],
            group_1_label="younger",
            group_2_label="older",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_two(
            gene_identifier=gene,
            comparison="leukocytosis",
            group_1_persons=sets_persons["calm"],
            group_2_persons=sets_persons["leukocytosis"],
            group_1_label="calm",
            group_2_label="leukocytosis",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_two(
            gene_identifier=gene,
            comparison="inflammation",
            group_1_persons=sets_persons["calm_chronic"],
            group_2_persons=sets_persons["inflammation"],
            group_1_label="calm",
            group_2_label="inflammation",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_two(
            gene_identifier=gene,
            comparison="ventilation",
            group_1_persons=sets_persons["breath"],
            group_2_persons=sets_persons["ventilation"],
            group_1_label="breath",
            group_2_label="ventilation",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
    return comparisons


def organize_genes_person_four_groups_signal_comparisons(
    genes=None,
    sets_persons=None,
    data_signals_genes_persons=None,
    data_gene_annotation=None,
    report=None,
):
    """
    Organizes and comparisons between four groups across multiple genes.

    arguments:
        genes (list<str>): identifiers of genes
        sets_persons (dict<list<str>>): identifiers of persons in groups by
            their properties
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        data_gene_annotation (object): Pandas data frame of genes' annotations
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): information about comparisons

    """

    comparisons = list()
    for gene in genes:
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="sex_ventilation",
            group_1_persons=sets_persons["male_breath"],
            group_2_persons=sets_persons["male_ventilation"],
            group_3_persons=sets_persons["female_breath"],
            group_4_persons=sets_persons["female_ventilation"],
            group_1_label="male_breath",
            group_2_label="male_vent",
            group_3_label="female_breath",
            group_4_label="female_vent",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="sex_leukocytosis",
            group_1_persons=sets_persons["male_calm"],
            group_2_persons=sets_persons["male_leukocytosis"],
            group_3_persons=sets_persons["female_calm"],
            group_4_persons=sets_persons["female_leukocytosis"],
            group_1_label="male_calm",
            group_2_label="male_leuk",
            group_3_label="female_calm",
            group_4_label="female_leuk",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="age_ventilation",
            group_1_persons=sets_persons["young_breath"],
            group_2_persons=sets_persons["young_ventilation"],
            group_3_persons=sets_persons["old_breath"],
            group_4_persons=sets_persons["old_ventilation"],
            group_1_label="young_breath",
            group_2_label="young_vent",
            group_3_label="old_breath",
            group_4_label="old_vent",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="age_leukocytosis",
            group_1_persons=sets_persons["young_calm"],
            group_2_persons=sets_persons["young_leukocytosis"],
            group_3_persons=sets_persons["old_calm"],
            group_4_persons=sets_persons["old_leukocytosis"],
            group_1_label="young_calm",
            group_2_label="young_leuk",
            group_3_label="old_calm",
            group_4_label="old_leuk",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="ager_ventilation",
            group_1_persons=sets_persons["younger_breath"],
            group_2_persons=sets_persons["younger_ventilation"],
            group_3_persons=sets_persons["older_breath"],
            group_4_persons=sets_persons["older_ventilation"],
            group_1_label="younger_breath",
            group_2_label="younger_vent",
            group_3_label="older_breath",
            group_4_label="older_vent",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="ager_leukocytosis",
            group_1_persons=sets_persons["younger_calm"],
            group_2_persons=sets_persons["younger_leukocytosis"],
            group_3_persons=sets_persons["older_calm"],
            group_4_persons=sets_persons["older_leukocytosis"],
            group_1_label="younger_calm",
            group_2_label="younger_leuk",
            group_3_label="older_calm",
            group_4_label="older_leuk",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="race_ventilation",
            group_1_persons=sets_persons["race_white_breath"],
            group_2_persons=sets_persons["race_white_ventilation"],
            group_3_persons=sets_persons["race_not_white_breath"],
            group_4_persons=sets_persons["race_not_white_ventilation"],
            group_1_label="white_breath",
            group_2_label="white_vent",
            group_3_label="other_breath",
            group_4_label="other_vent",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="race_leukocytosis",
            group_1_persons=sets_persons["race_white_calm"],
            group_2_persons=sets_persons["race_white_leukocytosis"],
            group_3_persons=sets_persons["race_not_white_calm"],
            group_4_persons=sets_persons["race_not_white_leukocytosis"],
            group_1_label="white_calm",
            group_2_label="white_leuk",
            group_3_label="other_calm",
            group_4_label="other_leuk",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="ventilation_age",
            group_1_persons=sets_persons["young_breath"],
            group_2_persons=sets_persons["old_breath"],
            group_3_persons=sets_persons["young_ventilation"],
            group_4_persons=sets_persons["old_ventilation"],
            group_1_label="breath_young",
            group_2_label="breath_old",
            group_3_label="vent_young",
            group_4_label="vent_old",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
        comparisons.append(compare_genes_signals_persons_groups_four(
            gene_identifier=gene,
            comparison="leukocytosis_age",
            group_1_persons=sets_persons["young_calm"],
            group_2_persons=sets_persons["old_calm"],
            group_3_persons=sets_persons["young_leukocytosis"],
            group_4_persons=sets_persons["old_leukocytosis"],
            group_1_label="calm_young",
            group_2_label="calm_old",
            group_3_label="leuk_young",
            group_4_label="leuk_old",
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
            report=report,
        ))
    return comparisons


def write_product_genes_signals_persons_groups(
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_comparisons_two = os.path.join(
        paths["groups"], "genes_comparisons_two_groups.pickle"
    )
    path_comparisons_four = os.path.join(
        paths["groups"], "genes_comparisons_four_groups.pickle"
    )
    # Write information to file.
    with open(path_comparisons_two, "wb") as file_product:
        pickle.dump(information["comparisons_two"], file_product)
    with open(path_comparisons_four, "wb") as file_product:
        pickle.dump(information["comparisons_four"], file_product)
    pass


def read_organize_report_write_genes_signals_persons_groups(
    paths=None,
    report=None,
):
    """
    Organizes the comparisons of genes' signals between groups of persons.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_genes_signals_persons_properties(
        cohort="selection",
        dock=paths["dock"],
    )
    # Standardize scale of genes' signals.
    data_signals_standard = prediction.standardize_genes_signals(
        data_signals_genes_persons=source["data_signals_genes_persons"],
        report=report,
    )
    # Organize groups of persons by their properties.
    sets_persons = organize_persons_properties_sets(
        data_persons_properties=source["data_persons_properties"],
        report=report,
    )
    # Organize genes' signals in groups of persons.
    # BHLHE40: ENSG00000134107
    # SAMHD1: ENSG00000101347

    comparisons_two = organize_genes_person_two_groups_signal_comparisons(
        genes=[
            "ENSG00000147050", # KDM6A
            "ENSG00000177575", # CD163
            "ENSG00000169738", # DCXR
            "ENSG00000134986", # NREP
            "ENSG00000154134", # ROBO3
            "ENSG00000153234", # NR4A2
        ],
        sets_persons=sets_persons,
        data_signals_genes_persons=data_signals_standard,
        data_gene_annotation=source["data_gene_annotation"],
        report=report,
    )
    comparisons_four = organize_genes_person_four_groups_signal_comparisons(
        genes=[
            "ENSG00000158050", # DUSP2
            "ENSG00000146013", # GFRA3
            "ENSG00000154134", # ROBO3
            "ENSG00000153234", # NR4A2
            "ENSG00000134986", # NREP
        ],
        sets_persons=sets_persons,
        data_signals_genes_persons=data_signals_standard,
        data_gene_annotation=source["data_gene_annotation"],
        report=report,
    )
    # Write information to file.
    information = dict()
    information["comparisons_two"] = comparisons_two
    information["comparisons_four"] = comparisons_four
    write_product_genes_signals_persons_groups(
        information=information,
        paths=paths,
    )
    pass


##########
# Pairwise correlations


def read_source_pairwise_gene_correlations(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Read genes sets.
    if False:
        genes_sets = read_source_organize_genes_sets_collection_candidacy_query(
            cohort=cohort,
            dock=dock,
        )
    genes_sets = read_source_genes_sets_combination_integration(
        cohort="selection",
        dock=dock,
    )
    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    # Return information.
    return {
        "genes_sets": genes_sets,
        "data_gene_annotation": data_gene_annotation,
        "data_signals_genes_persons": data_signals_genes_persons,
    }


def extract_ontology_gene_sets(
    source_ontology=None,
):
    """
    Defines sets of genes' identifiers from gene ontology enrichment.

    arguments:
        source_ontology (dict<object>): collection of Pandas data frames with
            identifiers and names of genes

    raises:

    returns:
        (dict<list<str>>): identifiers of genes

    """

    collection = dict()
    for name in source_ontology.keys():
        data = source_ontology[name]
        identifiers = data["identifier"].to_numpy().tolist()
        # Collection information.
        collection[name] = identifiers
    # Return information.
    return collection


def organize_genes_pairwise_signals_correlations(
    method=None,
    threshold_high=None,
    threshold_low=None,
    discovery=None,
    genes=None,
    data_signals_genes_persons=None,
):
    """
    Calculates Pearson correlation coefficients between pairs of genes.

    arguments:
        method (str): method for correlation, pearson, spearman, or kendall
        threshold_high (float): value must be greater than this threshold
        threshold_low (float): value must be less than this threshold
        discovery (float): value of false discovery rate for which to include
            correlation coefficients
        genes (list<str>): identifiers of genes for which to calculate pairwise
            correlations
        data_signals_genes_persons (object): Pandas data frame of genes'
            pan-tissue signals across persons

    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' pairwise
            correlations

    """

    # Filter query genes to ensure that they have signals.
    genes_valid = utility.filter_common_elements(
        list_one=genes,
        list_two=data_signals_genes_persons.columns.tolist(),
    )
    # Calculate correlation matrix
    data_correlation = utility.organize_feature_signal_correlations(
        method=method, # pearson (normal distribution), spearman
        threshold_high=threshold_high,
        threshold_low=threshold_low,
        count=2, # accommodate the value 1.0 for self pairs (A, A)
        discovery=discovery,
        features=genes_valid,
        data_signal=data_signals_genes_persons,
    )
    # Return information.
    return data_correlation


def select_translate_gene_identifiers_data_columns(
    genes_query=None,
    data_gene_annotation=None,
    data_signals_genes_persons=None,
):
    """
    Organize information for chart.

    arguments:
        genes_query (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons

    raises:

    returns:
        (object): Pandas data frame of pan-tissue signals across genes and
            persons

    """

    # Copy data.
    data_signals_genes_persons = data_signals_genes_persons.copy(deep=True)
    # Select data for genes of interest.
    data_selection = data_signals_genes_persons.loc[
        :, data_signals_genes_persons.columns.isin(genes_query)
    ]
    # Organize signals.
    # Translate identifiers of genes.
    #identifiers = data_signals_genes_persons.columns.to_list()
    translations = dict()
    for identifier in genes_query:
        translations[identifier] = assembly.access_gene_name(
            identifier=identifier,
            data_gene_annotation=data_gene_annotation,
        )
        pass
    data_selection.rename(
        columns=translations,
        inplace=True,
    )
    # Return information.
    return data_selection


def write_product_pairwise_gene_correlations(
    cohort=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_data_genes_correlations = os.path.join(
        paths[cohort]["correlation"],
        "data_genes_correlations.pickle"
    )
    path_genes_correlation = os.path.join(
        paths[cohort]["correlation"], "genes.pickle"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_genes_correlations"],
        path_data_genes_correlations
    )
    with open(path_genes_correlation, "wb") as file_product:
        pickle.dump(
            information["genes_correlation"], file_product
        )
    pass


def read_organize_report_write_pairwise_gene_correlations(
    paths=None,
    report=None,
):
    """
    Organizes evaluation of subpopulation structure on the basis of pan-tissue
    expression of genes of interest.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_pairwise_gene_correlations(
        cohort="selection",
        dock=paths["dock"],
    )
    # Specify set of genes for which to calculate pairwise correlations.
    if True:
        genes_correlation = (
            source["genes_sets"]["combination"]["covid19_multimodal"]
        )
    if False:
        genes_correlation = (
            source["genes_sets"]["prediction"]["priority"]
        )
    if False:
        genes_correlation = (
            source["genes_sets"]["combination"]
            ["covid19_up_multimodal"]
        )
    if False:
        genes_correlation = (
            source["genes_sets"]["query"]["correlation"]
        )
    if False:
        genes_correlation = (
            source["genes_sets"]["query"]["covid19_up_prediction"]
        )

    # Calculate correlations between pairs of genes.
    # Use Spearman correlations.
    data_correlation = organize_genes_pairwise_signals_correlations(
        method="spearman", # pearson (normal distribution), spearman
        threshold_high=0.75, # 1.0, 0.75, 0.5, 0.25, 0.0
        threshold_low=-0.75, # -1.0, -0.75, -0.5, -0.25, -0.0
        discovery=0.05,
        genes=genes_correlation,
        data_signals_genes_persons=source["data_signals_genes_persons"],
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("correlation genes: " + str(data_correlation.shape[0]))
        print(data_correlation)
        utility.print_terminal_partition(level=2)
    # Compile information.
    information = dict()
    information["data_genes_correlations"] = data_correlation
    information["genes_correlation"] = data_correlation.index.to_list()
    # Write information to file.
    write_product_pairwise_gene_correlations(
        cohort="selection",
        information=information,
        paths=paths,
    )
    pass


##########
# Gene query summary tables

# TODO: need to finish this...

def read_source_gene_sets_covid19_prediction_heritability(
    dock=None,
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

    # TODO: read in COVID-19 DE gene table
    # TODO: read in prediction context summary tables
    # TODO: read in heritability summary

    # Specify directories and files.
    path_regressions_discoveries = os.path.join(
        dock, "prediction", cohort, model, "regressions_discoveries",
        "regressions_discoveries.pickle"
    )
    pass


def read_organize_report_write_gene_sets_covid19_prediction_heritability(
    paths=None,
    report=None,
):
    """
    Organizes summary tables for gene's COVID-19 expression, associations to
    regression variables, and heritabilities of pan-tissue signals.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_gene_sets_covid19_prediction_heritability(
        dock=paths["dock"],
    )
    # Write information to file.
    information = dict()
    pass




#####################################
##################################################
##### old stuff#####


def read_source_genes_sets_function(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): regression model, either technique or hypothesis
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_sets_genes = os.path.join(
        dock, "function", "sets_genes.pickle"
    )
    # Read information from file.
    with open(path_sets_genes, "rb") as file_source:
        sets_genes = pickle.load(file_source)
    # Return information.
    return sets_genes


def read_source_gene_sets(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    # Read genes sets.
    genes_candidacy = read_source_genes_sets_candidacy(
        cohort=cohort,
        dock=dock,
    )
    genes_heritability = read_source_genes_sets_heritability(
        cohort=cohort,
        dock=dock,
    )
    genes_prediction = read_source_genes_sets_prediction(
        cohort=cohort,
        dock=dock,
    )
    genes_function = read_source_genes_sets_function(
        cohort=cohort,
        dock=dock,
    )
    # Return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_candidacy": genes_candidacy,
        "genes_heritability": genes_heritability,
        "genes_prediction": genes_prediction,
        "genes_function": genes_function,
    }


def allocate_query_items_sets(
    variable=None,
    query_items=None,
    sets=None,
    report=None,
    set_report=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variable (str): name of variable
        query_items (list<str>): identifiers of interest
        sets (dict<list<str>>): sets of identifiers
        report (bool): whether to print reports
        set_report (str): name of set for report

    raises:

    returns:
        (dict): sets of identifiers that match and do not match query

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("query variable: " + str(variable))
        print("count query items: " + str(len(query_items)))
        #print(sets.keys())
        print(
            "count total items in example set: " + str(len(sets[set_report]))
        )
        pass
    # Collect items by inclusion and exclusion in sets.
    sets_genes = dict()
    sets_genes["inclusion"] = dict()
    sets_genes["exclusion"] = dict()
    records = list()
    for set in sets.keys():
        sets_genes["inclusion"][set] = utility.filter_unique_common_elements(
            list_one=query_items,
            list_two=copy.deepcopy(sets[set]),
        )
        sets_genes["exclusion"][set] = utility.filter_unique_exclusion_elements(
            elements_exclusion=query_items,
            elements_total=copy.deepcopy(sets[set]),
        )
        record = dict()
        record["groups"] = set
        record["total"] = len(sets[set])
        record["query"] = len(sets_genes["inclusion"][set])
        record["other"] = len(sets_genes["exclusion"][set])
        records.append(record)
    data = utility.convert_records_to_dataframe(records=records)
    # Compile information.
    bin = dict()
    bin["sets_genes"] = sets_genes
    bin["data"] = data
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        if False:
            print(
                "count inclusion: " +
                str(len(bin["sets_genes"]["inclusion"][set_report]))
            )
            print(
                "count exclusion: " +
                str(len(bin["sets_genes"]["exclusion"][set_report]))
            )
        utility.print_terminal_partition(level=4)
        print(data)
        pass
    # Return information.
    return bin


def write_product_integration_gene_sets(
    cohort=None,
    variable=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        variable (str): name of independent regression variable
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_sets_genes = os.path.join(
        paths[cohort]["set"]["allocation"], str(variable + ".pickle")
    )
    path_data = os.path.join(
        paths[cohort]["set"]["cardinality"], str(variable + ".pickle")
    )
    # Write information to file.
    with open(path_sets_genes, "wb") as file_product:
        pickle.dump(information["sets_genes"], file_product)
    pandas.to_pickle(
        information["data"],
        path_data
    )
    pass


def collect_unique_genes_prediction_variables_ontology_sets(
    allocation=None,
    variables=None,
    sets=None,
    report=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        allocation (dict): collections of genes
        variables (list<str>): names of variables from prediction procedure
        sets (list<str>): names of sets from functional ontologies
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of genes

    """

    # Collect union of genes from variables and sets.
    genes_collection = list()
    for variable in variables:
        for set in sets:
            genes_new = allocation[variable]["inclusion"][set]
            genes_collection.extend(genes_new)
            pass
        pass
    # Collect unique union genes.
    genes_unique = utility.collect_unique_elements(
        elements_original=genes_collection
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "count of union genes from variables and sets: " +
            str(len(genes_unique))
        )
        pass
    # Return information.
    return genes_unique


def write_product_integration_prediction_ontology_genes(
    cohort=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        variable (str): name of independent regression variable
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_genes = os.path.join(
        paths[cohort]["set"]["prediction_ontology"], "genes.pickle"
    )
    path_genes_text = os.path.join(
        paths[cohort]["set"]["prediction_ontology"], "genes.txt"
    )
    # Write information to file.
    with open(path_genes, "wb") as file_product:
        pickle.dump(information["genes"], file_product)
    utility.write_file_text_list(
        elements=information["genes"],
        delimiter="\n",
        path_file=path_genes_text
    )
    pass


def read_organize_report_write_integration_gene_sets(paths=None):
    """
    Organizes sets of genes by integration of information from multiple
    procedures.

    1. genes with multimodal distributions of pan-tissue signals
    2. genes with heritable pan-tissue signals
    3. genes with enrichment allocation to functional ontological groups

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_gene_sets(
        cohort="selection",
        dock=paths["dock"],
    )
    # Allocations of genes in sets from functional ontologies.
    allocation = dict()
    for variable in source["genes_prediction"]["multimodal"].keys():
        bin_allocation = allocate_query_items_sets(
            variable=variable,
            query_items=source["genes_prediction"]["multimodal"][variable],
            sets=source["genes_function"],
            report=False,
            set_report="defense",
        )
        allocation[variable] = bin_allocation["sets_genes"]
        write_product_integration_gene_sets(
            cohort="selection",
            variable=variable,
            information=bin_allocation,
            paths=paths,
        )

    # Collect genes from multiple functional ontology sets of interest.
    # 1. specify prediction variable: ventilation_binary_scale
    # 2. specify functional ontology sets: ["defense", "inflammation", "cytokine"]
    # ventilation_binary_scale
    # mononucleosis_binary_scale
    # leukocyte_binary_scale
    # inflammation_binary_scale

    genes_union = collect_unique_genes_prediction_variables_ontology_sets(
        allocation=allocation,
        variables=[
            "ventilation_binary_scale",
            #"mononucleosis_binary_scale",
            #"leukocyte_binary_scale",
            #"inflammation_binary_scale",
        ],
        sets=[
            "defense",
            "migration",
            "proliferation",
            "inflammation",
            "cytokine",
        ],
        report=False,
    )
    # Compile information.
    information = dict()
    information["genes"] = genes_union
    # Write information to file.
    write_product_integration_prediction_ontology_genes(
        cohort="selection",
        information=information,
        paths=paths,
    )
    pass


##########
# Interaction gene sets


def read_source_genes_sets_prediction_interaction_cohort_model(
    cohort=None,
    model=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): name of regression model
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<list<str>>): sets of genes' identifiers

    """

    # Specify directories and files.
    path_sets_genes = os.path.join(
        dock, "prediction", cohort, model, "genes",
        "sets_genes.pickle"
    )
    # Read information from file.
    with open(path_sets_genes, "rb") as file_source:
        sets_genes = pickle.load(file_source)
    # Return information.
    return sets_genes


def read_source_genes_sets_prediction_interaction(
    cohort=None,
    dock=None,
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<list<str>>): sets of genes' identifiers

    """

    # Read source information from file.
    sets_main = read_source_genes_sets_prediction_interaction_cohort_model(
        cohort=cohort,
        model="selection_main",
        dock=dock,
    )
    sets_sex_ventilation = (
        read_source_genes_sets_prediction_interaction_cohort_model(
            cohort=cohort,
            model="selection_sex_ventilation",
            dock=dock,
    ))
    sets_age_ventilation = (
        read_source_genes_sets_prediction_interaction_cohort_model(
            cohort=cohort,
            model="selection_age_ventilation",
            dock=dock,
    ))
    # Compile information.
    bin = dict()
    bin["main"] = sets_main
    bin["sex_ventilation"] = sets_sex_ventilation
    bin["age_ventilation"] = sets_age_ventilation
    # Return information.
    return bin


def collect_inclusion_exclusion_unique_elements(
    inclusion=None,
    exclusion_one=None,
    exclusion_two=None,
):
    """
    Organizes sets of genes by integration of information from multiple
    procedures.

    arguments:
        inclusion (list): set of elements to keep
        exclusion_one (list): set of elements to remove from inclusion
        exclusion_two (list): set of elements to remove from inclusion

    raises:

    returns:
        (list): set of unique elements

    """

    # Combine exclusion elements.
    exclusion = copy.deepcopy(exclusion_one)
    exclusion.extend(exclusion_two)
    # Exclude and collect unique elements.
    inclusion_unique = utility.filter_unique_exclusion_elements(
        elements_exclusion=exclusion,
        elements_total=inclusion,
    )
    # Return information.
    return inclusion_unique


def organize_prediction_interaction_gene_sets_distribution(
    sex_age_ventilation=None,
    sex_ventilation=None,
    age_ventilation=None,
):
    """
    Organizes sets of genes by integration of information from multiple
    procedures.

    arguments:
        sex_age_ventilation (dict): sets of genes from regression
        sex_ventilation (dict): sets of genes from regression
        age_ventilation (dict): sets of genes from regression

    raises:

    returns:
        (dict): sets of genes

    """

    # Collect unique genes that associate with sex-ventilation interaction.
    sex_ventilation_unique = collect_inclusion_exclusion_unique_elements(
        inclusion=sex_ventilation["sex_risk*ventilation_binary_scale"],
        exclusion_one=sex_age_ventilation["ventilation_binary_scale"],
        exclusion_two=list(),
    )
    # Collect unique genes that associate with age-ventilation interaction.
    age_ventilation_unique = collect_inclusion_exclusion_unique_elements(
        inclusion=age_ventilation["age*ventilation_binary_scale"],
        exclusion_one=sex_age_ventilation["ventilation_binary_scale"],
        exclusion_two=list(),
    )
    # Compile information.
    bin = dict()
    bin["sex_ventilation_unique"] = sex_ventilation_unique
    bin["age_ventilation_unique"] = age_ventilation_unique
    # Return information.
    return bin


def write_product_prediction_interaction_genes(
    cohort=None,
    distribution=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        distribution (str): group of genes by their distribution of pan-tissue
            signals
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_sex_ventilation = os.path.join(
        paths[cohort]["set"]["prediction_interaction"][distribution],
        "sex_ventilation.pickle"
    )
    path_sex_ventilation_text = os.path.join(
        paths[cohort]["set"]["prediction_interaction"][distribution],
        "sex_ventilation.txt"
    )
    path_age_ventilation = os.path.join(
        paths[cohort]["set"]["prediction_interaction"][distribution],
        "age_ventilation.pickle"
    )
    path_age_ventilation_text = os.path.join(
        paths[cohort]["set"]["prediction_interaction"][distribution],
        "age_ventilation.txt"
    )
    # Write information to file.
    with open(path_sex_ventilation, "wb") as file_product:
        pickle.dump(information["sex_ventilation_unique"], file_product)
    utility.write_file_text_list(
        elements=information["sex_ventilation_unique"],
        delimiter="\n",
        path_file=path_sex_ventilation_text
    )
    with open(path_age_ventilation, "wb") as file_product:
        pickle.dump(information["age_ventilation_unique"], file_product)
    utility.write_file_text_list(
        elements=information["age_ventilation_unique"],
        delimiter="\n",
        path_file=path_age_ventilation_text
    )
    pass

# TODO: I think this sub-procedure is now obsolete...
# TODO: I now know that it is better to include the main terms (sex, age, ventilation) in the model
# TODO: with the interaction terms (sex*ventilation, age*ventilation)...
# TODO: this sub-procedure combines genes from separate models... not ideal...
def read_organize_report_write_prediction_interaction_gene_sets(paths=None):
    """
    Organizes sets of genes by integration of information from multiple
    procedures.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_genes_sets_prediction_interaction(
        cohort="selection",
        dock=paths["dock"],
    )
    # Organize sets for genes with any distribution of pan-tissue signals.
    sets_any = organize_prediction_interaction_gene_sets_distribution(
        sex_age_ventilation=source["main"]["any"],
        sex_ventilation=source["sex_ventilation"]["any"],
        age_ventilation=source["age_ventilation"]["any"],
    )
    # Organize sets for genes with multimodal distribution of pan-tissue
    # signals.
    sets_multimodal = organize_prediction_interaction_gene_sets_distribution(
        sex_age_ventilation=source["main"]["multimodal"],
        sex_ventilation=source["sex_ventilation"]["multimodal"],
        age_ventilation=source["age_ventilation"]["multimodal"],
    )
    # Write information to file.
    write_product_prediction_interaction_genes(
        cohort="selection",
        distribution="any",
        information=sets_any,
        paths=paths,
    )
    write_product_prediction_interaction_genes(
        cohort="selection",
        distribution="multimodal",
        information=sets_multimodal,
        paths=paths,
    )
    pass


# Summary report on integration of genes


def organize_genes_integration_report(
    genes=None,
    data_gene_annotation=None,
    data_gene_distribution_report=None,
    data_genes_permutation_probabilities=None,
    data_genes_heritabilities=None,
    data_regression_genes=None,
):
    """
    Organizes information about genes.

    Data structure
    identifier        name    coefficient  p_coefficient  dip     p_dip   ...
     ENSG00000078808   TK421   0.75         0.003          0.12    0.005  ...
     ENSG00000029363   R2D2    0.55         0.975          0.23    0.732  ...

    arguments:
        genes (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_gene_distribution_report (object): Pandas data frame of
            information about genes' distributions
        data_genes_permutation_probabilities (object): Pandas data frame of
            genes' probabilities from permutations
        data_genes_heritabilities (object): Pandas data frame of genes'
            heritabilities
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across genes

    raises:

    returns:
        (object): Pandas data frame of genes' properties

    """

    # Organize data.

    # Collect information about genes.
    records = list()
    # Iterate on genes.
    for gene in genes:
        # Create record for gene.
        record = dict()
        # Access basic information about gene.
        record["identifier"] = gene
        record["name"] = assembly.access_gene_name(
            identifier=gene,
            data_gene_annotation=data_gene_annotation,
        )
        record["chromosome"] = data_gene_annotation.loc[gene, "seqname"]
        record["start"] = data_gene_annotation.loc[gene, "start"]
        record["end"] = data_gene_annotation.loc[gene, "end"]
        # Permutation.
        record["permutation_coefficient"] = (
            data_genes_permutation_probabilities.loc[
                gene, "discovery_coefficient"
            ]
        )
        record["permutation_mixture"] = (
            data_genes_permutation_probabilities.loc[
                gene, "discovery_mixture"
            ]
        )
        record["permutation_dip"] = data_genes_permutation_probabilities.loc[
            gene, "discovery_dip"
        ]
        # Heritability.
        record["heritability_phenotype"] = (
            data_genes_heritabilities.loc[gene, "phenotype"]
        )
        record["heritability_genotype"] = (
            data_genes_heritabilities.loc[gene, "genotype"]
        )
        record["heritability_proportion"] = (
            data_genes_heritabilities.loc[gene, "proportion"]
        )
        record["heritability_discovery"] = (
            data_genes_heritabilities.loc[gene, "discovery"]
        )
        # Prediction.
        record["prediction_r_square"] = (
            data_regression_genes.loc[gene, "r_square"]
        )
        record["prediction_sex_discovery"] = (
            data_regression_genes.loc[gene, "female_scale_discovery"]
        )
        record["prediction_age_discovery"] = (
            data_regression_genes.loc[gene, "age_scale_discovery"]
        )
        record["prediction_body_discovery"] = (
            data_regression_genes.loc[gene, "body_scale_discovery"]
        )
        record["prediction_hardiness_discovery"] = (
            data_regression_genes.loc[gene, "hardiness_scale_discovery"]
        )


        # Compile information.
        records.append(record)
        pass

    # Organize data.
    data = utility.convert_records_to_dataframe(
        records=records
    )
    data = data[[
        "identifier",
        "name",
        "chromosome",
        "start",
        "end",
        "permutation_coefficient",
        "permutation_mixture",
        "permutation_dip",
        "heritability_phenotype",
        "heritability_genotype",
        "heritability_proportion",
        "heritability_discovery",
        "prediction_r_square",
        "prediction_sex_discovery",
        "prediction_age_discovery",
        "prediction_body_discovery",
        "prediction_hardiness_discovery"
    ]]

    data.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    data.rename_axis(
        index="identifier",
        axis="index",
        copy=False,
        inplace=True,
    )
    data.rename_axis(
        columns="properties",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Return information.
    return data



###############scrap???############################


# Rank

def select_rank_genes(
    data_genes_integration=None,
    selection=None,
    rank=None,
):
    """
    Prepares ranks of genes for subsequent functional analyses.

    arguments:
        data_genes_integration (object): Pandas data frame of genes'
            properties, bimodality scores, and probabilities
        selection (list<str>): genes for selection
        rank (str): property to use for ranks

    raises:

    returns:
        (object): Pandas data frame of genes' properties

    """

    # Copy data.
    data_copy = data_genes_integration.copy(deep=True)
    # Select data for genes of interest.
    data = data_copy.loc[
        data_copy.index.isin(selection), :
    ]
    # Rank genes.
    data.sort_values(
        by=[rank],
        axis="index",
        ascending=False,
        inplace=True,
    )
    # Organize information.
    # Return information.
    return data


# Exportation


def export_genes_function(
    data_genes=None,
    identifier=None,
    rank=None,
):
    """
    Prepares table of genes for export to functional analysis.

    arguments:
        data_genes (object): Pandas data frame of genes' properties, bimodality
            scores, and probabilities
        identifier (str): identifier to use, "ensembl" or "hugo"
        rank (bool): whether to include gene rankings in the table

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality scores,
            and probabilities

    """

    # Copy data.
    data_copy = data_genes.copy(deep=True)
    # Organize data.
    data_copy.reset_index(
        level=None,
        inplace=True
    )
    # Select columns.
    columns = list()
    if identifier == "ensembl":
        columns.append("identifier")
    elif identifier == "hugo":
        columns.append("name")
    if rank:
        columns.append("combination")
    data = data_copy.loc[
        :, data_copy.columns.isin(columns)
    ]

    # Return information.
    return data


# Write.


def write_product_query(dock=None, information=None):
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
    path_integration = os.path.join(dock, "integration")
    path_query = os.path.join(
        path_integration, "query_correlations"
    )
    utility.create_directories(path_query)

    # Iterate on sets.
    collection = information["bin_query_correlation"]
    for name in collection.keys():
        # Specify path to file.
        file_name = ("data_" + name + ".pickle")
        path_data = os.path.join(
            path_query, file_name
        )
        # Write information to file.
        collection[name].to_pickle(
            path_data
        )
        pass

    pass


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

    # Ontology.
    write_product_query(
        dock=dock,
        information=information,
    )

    # Specify directories and files.
    path_integration = os.path.join(dock, "integration")
    utility.create_directory(path_integration)

    path_data_correlation_unimodal = os.path.join(
        path_integration, "data_correlation_genes_unimodal.pickle"
    )
    path_data_correlation_multimodal = os.path.join(
        path_integration, "data_correlation_genes_multimodal.pickle"
    )
    path_data_correlation_multimodal_hardiness = os.path.join(
        path_integration, "data_correlation_multimodal_hardiness.pickle"
    )
    path_data_correlation_multimodal_sex_age_body = os.path.join(
        path_integration, "data_correlation_multimodal_sex_age_body.pickle"
    )
    path_data_correlation_multimodal_union = os.path.join(
        path_integration, "data_correlation_multimodal_union.pickle"
    )

    path_data_persons_genes_components = os.path.join(
        path_integration, "data_persons_genes_components.pickle"
    )
    path_data_persons_genes_variances = os.path.join(
        path_integration, "data_persons_genes_variances.pickle"
    )

    path_data_genes_integration = os.path.join(
        path_integration, "data_genes_integration.pickle"
    )
    path_data_genes_integration_text = os.path.join(
        path_integration, "data_genes_integration.tsv"
    )


    if False:

        path_export_genes_selection = os.path.join(
            path_integration, "export_genes_selection.tsv"
        )
        path_export_genes_total = os.path.join(
            path_integration, "export_genes_total.tsv"
        )

    # Write information to file.
    information["data_correlation_genes_unimodal"].to_pickle(
        path_data_correlation_unimodal
    )
    information["data_correlation_genes_multimodal"].to_pickle(
        path_data_correlation_multimodal
    )
    information["data_correlation_multimodal_hardiness"].to_pickle(
        path_data_correlation_multimodal_hardiness
    )
    information["data_correlation_multimodal_sex_age_body"].to_pickle(
        path_data_correlation_multimodal_sex_age_body
    )
    information["data_correlation_multimodal_union"].to_pickle(
        path_data_correlation_multimodal_union
    )

    information["data_persons_genes_components"].to_pickle(
        path_data_persons_genes_components
    )
    information["data_persons_genes_variances"].to_pickle(
        path_data_persons_genes_variances
    )

    information["data_genes_integration"].to_pickle(
        path_data_genes_integration
    )
    information["data_genes_integration"].to_csv(
        path_or_buf=path_data_genes_integration_text,
        columns=None,
        sep="\t",
        na_rep="",
        header=True,
        index=True,
    )

    if False:
        with open(path_genes_integration, "wb") as file_product:
            pickle.dump(
                information["genes_integration"], file_product
            )

        # Exportation
        information["data_export_genes_selection"].to_csv(
            path_or_buf=path_export_genes_selection,
            columns=None,
            sep="\t",
            na_rep="",
            header=False,
            index=False,
        )
        information["data_export_genes_total"].to_csv(
            path_or_buf=path_export_genes_total,
            columns=None,
            sep="\t",
            na_rep="",
            header=False,
            index=False,
        )

    pass


###############################################################################
# Procedure

# TODO: split this procedure up into subroutines
# 1. organize sets of multimodal genes in functional groups from gene ontologies
# 1.1. all multimodal genes in ontology functional groups without any query genes
# 1.2. all multimodal genes in ontology functional groups with query genes mapped onto these groups
# 2. organize pairwise correlation matrices for sets of genes
# 3. maybe also organize code here for the sort and cluster asymmetric matrices of genes' pan-tissue signals across persons...


def execute_procedure(dock=None):
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

    # TODO: introduce "report" flag for each sub-routine...

    # Organize sets of genes from collection and candidacy.
    # These sets are top-down, so do not include with the bottom-up query
    # genes.
    read_organize_report_write_collection_candidacy_gene_sets(
        paths=paths
    )

    # Organize contingency table comparisons.
    read_organize_report_write_contingency_table_comparisons(
        paths=paths,
        report=True,
    )

    # Calculate and organize correlations in pan-tissue signals between pairs
    # of relevant genes.
    read_organize_report_write_pairwise_gene_correlations(
        paths=paths,
        report=False,
    )

    # Organize comparisons of genes' signals between groups of persons.
    read_organize_report_write_genes_signals_persons_groups(
        paths=paths,
        report=False,
    )


    if False:
        # Integrate and organize summary tables for query sets of genes with their
        # COVID-19 differential expression, variable associations from regressions,
        # and heritabilities.
        read_organize_report_write_gene_sets_covid19_prediction_heritability(
            paths=paths,
            report=True,
        )

    # TODO: this subprocedure is now obsolete... ?
    # Organize sets of genes that are specific to regression interaction
    # variables for sex, age, and ventilation.
    # Subtract from these sets genes that associate with each variable alone.
    #read_organize_report_write_prediction_interaction_gene_sets(paths=paths)

    if False:
        # Organize sets of genes by integration of distribution modality,
        # functional gene ontologies, regression associations, and queries.
        read_organize_report_write_integration_gene_sets(paths=paths)

    if False:
        ##########
        # Correlations between pairs of genes

        # Calculate and organize correlations for groups of genes from regression
        # analysis.
        bin_prediction = organize_gene_correlations_multimodal_prediction(
            sets=source["sets_genes_prediction"],
            data_signals_genes_persons=source["data_signals_genes_persons"],
        )

        # Calculate and organize correlations for groups of genes from Gene
        # Ontology enrichment analysis.
        bin_query_correlation = organize_gene_correlations_multimodal_query(
            data_signals_genes_persons=source["data_signals_genes_persons"],
            genes_selection=source["genes_selection"],
            sets_query=source["sets_query"],
        )

        ##########
        # Integration

        # Select genes that are both multimodal and significant from permutation.
        collection = dict()
        collection["multimodal"] = source["genes_multimodal"]
        collection["probability"] = source["genes_probability"]
        collection["prediction"] = source["genes_prediction"]
        genes_multimodal_permutation_prediction = utility.select_elements_by_sets(
            names=["multimodal", "probability", "prediction"],
            sets=collection,
            count=3,
        )
        utility.print_terminal_partition(level=2)
        print(
            "genes significant, multimodal, and predictable: " +
            str(len(genes_multimodal_permutation_prediction))
        )

        # Integrate and organize information about all genes.
        data_genes_integration = organize_genes_integration_report(
            genes=genes_multimodal_permutation_prediction,
            data_gene_annotation=source["data_gene_annotation"],
            data_gene_distribution_report=source["data_gene_distribution_report"],
            data_genes_permutation_probabilities=(
                source["data_genes_permutation_probabilities"]
            ),
            data_genes_heritabilities=(
                source["data_genes_heritabilities_complex"]
            ),
            data_regression_genes=source["data_regression_genes"],
        )

        utility.print_terminal_partition(level=1)
        print(data_genes_integration)

    if False:

        # Select genes integration of methods from candidacy, probability, and
        # heritability procedures.
        sets_genes_integration = dict()
        sets_genes_integration["candidacy"] = source["genes_candidacy"]
        sets_genes_integration["probability"] = source["genes_probability"]
        sets_genes_integration["heritability"] = source["genes_heritability"]
        genes_integration = utility.select_elements_by_sets(
            names=["candidacy", "probability", "heritability"],
            sets=sets_genes_integration,
            count=3,
        )
        utility.print_terminal_partition(level=2)
        print(
            "selection of genes by candidacy, probability, and heritability: " +
            str(len(genes_integration))
        )
        utility.print_terminal_partition(level=2)

        # Select genes of interest.
        data_genes_selection = select_rank_genes(
            data_genes_integration=data_genes_integration,
            selection=genes_integration,
            rank="combination",
        )
        print(data_genes_selection)

        # Exportation of genes for functional analysis.
        utility.print_terminal_partition(level=1)
        data_export_genes_selection = export_genes_function(
            data_genes=data_genes_selection,
            identifier="ensembl",
            rank=False,
        )
        print(data_export_genes_selection)
        data_export_genes_total = export_genes_function(
            data_genes=data_genes_integration,
            identifier="ensembl",
            rank=False,
        )
        print(data_export_genes_total)

        # Define genes of interest.
        # Genes of interest are few for thorough summary.
        genes_report = define_report_genes(
            data_summary_genes=data_summary_genes,
            rank="combination"
        )

        # Organize thorough summaries for a few genes of interest.
        reports = prepare_reports_genes(
            genes=genes_report,
            data_gene_annotation=source["data_gene_annotation"],
            genes_signals_patients_tissues=(
                source["genes_signals_patients_tissues"]
            ),
            shuffles=source["shuffles"][0:500],
            genes_scores_distributions=source["genes_scores_distributions"],
        )

        # Compile information.
        information = {
            "data_genes_integration": data_genes_integration,
            "data_correlation_genes_unimodal": data_correlation_genes_unimodal,
            "data_correlation_genes_multimodal": data_correlation_genes_multimodal,
            "data_correlation_multimodal_hardiness": bin_prediction["hardiness"],
            "data_correlation_multimodal_sex_age_body": (
                bin_prediction["sex_age_body"]
            ),
            "data_correlation_multimodal_union": bin_prediction["union"],

            "bin_query_correlation": bin_query_correlation,

            #"data_genes_selection": data_genes_selection,
            #"data_export_genes_selection": data_export_genes_selection,
            #"data_export_genes_total": data_export_genes_total,
        }
        #Write product information to file.
        write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
