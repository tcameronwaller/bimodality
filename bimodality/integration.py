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
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=paths["integration"])
    utility.create_directory(path=paths["integration"])
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
        paths[cohort]["set"]["allocation"] = os.path.join(
            paths["integration"], cohort, "set", "allocation"
        )
        paths[cohort]["set"]["prediction_ontology"] = os.path.join(
            paths["integration"], cohort, "set", "prediction_ontology"
        )
        paths[cohort]["population"] = os.path.join(
            paths["integration"], cohort, "population"
        )
        paths[cohort]["correlation"] = os.path.join(
            paths["integration"], cohort, "correlation"
        )
        # Initialize directories.
        utility.create_directories(path=paths[cohort]["set"]["cardinality"])
        utility.create_directories(path=paths[cohort]["set"]["allocation"])
        utility.create_directories(
            path=paths[cohort]["set"]["prediction_ontology"]
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
    path_annotation = os.path.join(dock, "annotation_2020-04-28")
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

    genes_cytokine = read_source_annotation_query_genes_set(
        set="cytokine",
        dock=dock,
    )
    genes_antigen = read_source_annotation_query_genes_set(
        set="antigen",
        dock=dock,
    )
    genes_infection = read_source_annotation_query_genes_set(
        set="infection",
        dock=dock,
    )
    genes_inflammation = read_source_annotation_query_genes_set(
        set="inflammation",
        dock=dock,
    )
    genes_heme = read_source_annotation_query_genes_set(
        set="heme",
        dock=dock,
    )
    genes_oxidation = read_source_annotation_query_genes_set(
        set="oxidation",
        dock=dock,
    )
    genes_population = read_source_annotation_query_genes_set(
        set="population",
        dock=dock,
    )


    # Compile and return information.
    return {
        "cytokine": genes_cytokine,
        "antigen": genes_antigen,
        "infection": genes_infection,
        "inflammation": genes_inflammation,
        "heme": genes_heme,
        "oxidation": genes_oxidation,
        "population": genes_population,
    }


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
    path_annotation = os.path.join(dock, "annotation_2020-04-28")
    path_data_genes_query = os.path.join(path_annotation, "genes_query.csv")
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
    path_annotation = os.path.join(dock, "annotation_2020-04-28")
    path_data_genes_query = os.path.join(path_annotation, "genes_query.csv")
    data_genes_query = pandas.read_csv(
        path_data_genes_query,
        sep="\t",
        header=0,
    )
    genes = data_genes_query["identifier"].to_list()
    # Return information.
    return genes


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
    path_genes_selection = os.path.join(
        dock, "selection", "tight", "samples_genes_signals",
        "genes.pickle"
    )
    path_genes_unimodal = os.path.join(
        dock, "candidacy", cohort, "unimodal", "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        dock, "candidacy", cohort, "multimodal", "genes_multimodal.pickle"
    )
    # Read information from file.
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)
    # Compile information.
    bin = dict()
    bin["selection"] = genes_selection
    bin["unimodal"] = genes_unimodal
    bin["multimodal"] = genes_multimodal
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
        (object): source information

    """

    # Specify directories and files.
    path_genes_selection = os.path.join(
        dock, "heritability", cohort, "collection", "genes_selection.pickle"
    )
    path_genes_unimodal = os.path.join(
        dock, "heritability", cohort, "collection", "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        dock, "heritability", cohort, "collection", "genes_multimodal.pickle"
    )
    # Read information from file.
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)
    # Compile information.
    bin = dict()
    bin["selection"] = genes_selection
    bin["unimodal"] = genes_unimodal
    bin["multimodal"] = genes_multimodal
    # Return information.
    return bin


def read_source_genes_sets_prediction(
    cohort=None,
    model=None,
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
        dock, "prediction", cohort, model, "genes", "sets_genes.pickle"
    )
    # Read information from file.
    with open(path_sets_genes, "rb") as file_source:
        sets_genes = pickle.load(file_source)
    # Return information.
    return sets_genes


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
        model="hypothesis",
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
        report=True,
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
# Populations by gene expression
# 1. selection cohort
# 2. respiration cohort
# 3. ventilation cohort


def read_source_gene_person_populations(
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
    print(source["genes_function"].keys())

    """

    # Specify directories and files.
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )
    # Read information from file.
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    # Read genes sets.
    genes_candidacy = read_source_genes_sets_candidacy(
        cohort=cohort,
        dock=dock,
    )
    # Return information.
    return {
        "data_signals_genes_persons": data_signals_genes_persons,
        "genes_candidacy": genes_candidacy,
    }


def organize_persons_genes_components(
    genes=None,
    data_signals_genes_persons=None,
    report=None,
):
    """
    Organizes a principal components analysis on genes' pan-tissue signals as
    features across persons as instances.

    arguments:
        genes (list<str>): identifiers of genes
        data_signals_genes_persons (object): Pandas data frame of genes'
            pan-tissue signals across persons
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' pairwise
            correlations

    """

    # Copy data.
    data_signals = data_signals_genes_persons.copy(deep=True)
    # Select genes of interest.
    data_selection = data_signals.loc[
        :, data_signals.columns.isin(genes)
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Selection of genes with pan-tissue signals across persons.")
        utility.print_terminal_partition(level=3)
        print(data_selection)
    # Reduce dimensionality.
    components = min(int(len(genes)), int(data_selection.shape[0]))
    result = utility.calculate_principal_components(
        data=data_selection,
        components=components,
        report=report,
    )
    # Return information.
    return result


def write_product_gene_person_populations(
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
    path_data_components = os.path.join(
        paths[cohort]["population"], "data_persons_genes_components.pickle"
    )
    path_data_variances = os.path.join(
        paths[cohort]["population"], "data_persons_genes_variances.pickle"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_persons_genes_components"],
        path_data_components
    )
    pandas.to_pickle(
        information["data_persons_genes_variances"],
        path_data_variances
    )
    pass


def read_organize_report_write_integration_gene_person_populations(
    paths=None
):
    """
    Organizes evaluation of subpopulation structure on the basis of pan-tissue
    expression of genes of interest.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_gene_person_populations(
        cohort="selection",
        dock=paths["dock"],
    )
    # Calculate principal components on genes across persons.
    bin = organize_persons_genes_components(
        genes=source["genes_candidacy"]["multimodal"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
        report=False,
    )
    # Compile information.
    information = dict()
    information["data_persons_genes_components"] = bin[
        "data_observations_components"
    ]
    information["data_persons_genes_variances"] = bin[
        "data_components_variances"
    ]
    # Write information to file.
    write_product_gene_person_populations(
        cohort="selection",
        information=information,
        paths=paths,
    )
    pass


##########
# Pairwise correlations


def organize_gene_correlations_multimodal_prediction(
    sets=None,
    data_signals_genes_persons=None,
):
    """
    Calculates Pearson correlation coefficients between pairs of genes.

    arguments:
        sets (dict<dict<list<str>>>): sets of genes that associate
            significantly to variables in regression
        data_signals_genes_persons (object): Pandas data frame of genes'
            pan-tissue signals across persons


    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' pairwise
            correlations

    """

    # Calculate correlations between genes that associate with hypothetical
    # variables in regression.
    # Hardiness.
    data_hardiness = utility.organize_feature_signal_correlations(
        method="spearman", # pearson (normal distribution), spearman
        threshold_high=0.0, # 1.0, 0.75, 0.5, 0.0
        threshold_low=-0.0, # -1.0, -0.75, -0.5, -0.0
        count=2, # accommodate the value 1.0 for self pairs (A, A)
        discovery=0.05,
        features=sets["multimodal"]["hardiness_scale"],
        data_signal=data_signals_genes_persons,
    )
    # Sex, age, body.
    genes_sex_age_body = list()
    genes_sex_age_body.extend(sets["multimodal"]["female_scale"])
    genes_sex_age_body.extend(sets["multimodal"]["age_scale"])
    genes_sex_age_body.extend(sets["multimodal"]["body_scale"])
    genes_sex_age_body_unique = utility.collect_unique_elements(
        elements_original=genes_sex_age_body
    )
    data_sex_age_body = utility.organize_feature_signal_correlations(
        method="spearman", # pearson (normal distribution), spearman
        threshold_high=0.0, # 1.0, 0.75, 0.5, 0.0
        threshold_low=-0.0, # -1.0, -0.75, -0.5, -0.0
        count=2, # accommodate the value 1.0 for self pairs (A, A)
        discovery=0.05,
        features=genes_sex_age_body_unique,
        data_signal=data_signals_genes_persons,
    )
    #utility.print_terminal_partition(level=1)
    #print("data_sex_age_body")
    #utility.print_terminal_partition(level=2)
    #print(data_sex_age_body)
    # Sex, age, body.
    genes_union = list()
    genes_union.extend(sets["multimodal"]["female_scale"])
    genes_union.extend(sets["multimodal"]["age_scale"])
    genes_union.extend(sets["multimodal"]["body_scale"])
    genes_union.extend(sets["multimodal"]["hardiness_scale"])
    genes_union_unique = utility.collect_unique_elements(
        elements_original=genes_union
    )
    data_union = utility.organize_feature_signal_correlations(
        method="spearman", # pearson (normal distribution), spearman
        threshold_high=0.5, # 1.0, 0.75, 0.5, 0.0
        threshold_low=-0.5, # -1.0, -0.75, -0.5, -0.0
        count=2, # accommodate the value 1.0 for self pairs (A, A)
        discovery=0.05,
        features=genes_union_unique,
        data_signal=data_signals_genes_persons,
    )
    # Compile information.
    collection = dict()
    collection["hardiness"] = data_hardiness
    collection["sex_age_body"] = data_sex_age_body
    collection["union"] = data_union

    # Return information.
    return collection


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


def organize_gene_correlations_multimodal_query(
    data_signals_genes_persons=None,
    genes_selection=None,
    sets_query=None,
):
    """
    Calculates Pearson correlation coefficients between pairs of genes.

    arguments:
        data_signals_genes_persons (object): Pandas data frame of genes'
            pan-tissue signals across persons
        genes_selection (list<str>): identifiers of genes
        sets_query (dict<list<str>>): identifiers of genes in sets

    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' pairwise
            correlations

    """

    # Function sets.
    collection = dict()
    for name in sets_query:
        # Filter genes.
        set_genes = utility.filter_common_elements(
            list_one=sets_query[name],
            list_two=genes_selection,
        )
        # Calculate correlation matrix
        data_correlation = utility.organize_feature_signal_correlations(
            method="spearman", # pearson (normal distribution), spearman
            threshold_high=0.0, # 1.0, 0.75, 0.5, 0.0
            threshold_low=-0.0, # -1.0, -0.75, -0.5, -0.0
            count=2, # accommodate the value 1.0 for self pairs (A, A)
            discovery=0.05,
            features=set_genes,
            data_signal=data_signals_genes_persons,
        )
        # Compile information.
        collection[name] = data_correlation

    # Return information.
    return collection


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

    # Organize sets of genes by integration of distribution modality,
    # functional gene ontologies, regression associations, and queries.
    read_organize_report_write_integration_gene_sets(paths=paths)
    # Organize recognition of groups of persons from pan-tissue expression of
    # genes of interest.
    read_organize_report_write_integration_gene_person_populations(paths=paths)

    if False:
        ##########
        # Correlations between pairs of genes

        # Calculate correlations between pairs of genes.
        # Use Spearman correlations for both unimodal and multimodal.
        data_correlation_genes_unimodal = (
            utility.organize_feature_signal_correlations(
                method="spearman", # pearson (normal distribution), spearman
                threshold_high=0.5, # 1.0, 0.75, 0.5, 0.0
                threshold_low=-0.5, # -1.0, -0.75, -0.5, -0.0
                count=2, # accommodate the value 1.0 for self pairs (A, A)
                discovery=0.05,
                features=source["genes_unimodal"],
                data_signal=source["data_signals_genes_persons"],
        ))
        utility.print_terminal_partition(level=2)
        print("Unimodal genes...")
        utility.print_terminal_partition(level=3)
        print(data_correlation_genes_unimodal.shape)
        data_correlation_genes_multimodal = (
            utility.organize_feature_signal_correlations(
                method="spearman", # pearson (normal distribution), spearman
                threshold_high=0.5, # 1.0, 0.75, 0.5, 0.0
                threshold_low=-0.5, # -1.0, -0.75, -0.5, -0.0
                count=2, # accommodate the value 1.0 for self pairs (A, A)
                discovery=0.05,
                features=source["genes_multimodal"],
                data_signal=source["data_signals_genes_persons"],
        ))
        utility.print_terminal_partition(level=2)
        print("Multimodal genes...")
        utility.print_terminal_partition(level=3)
        print(data_correlation_genes_multimodal.shape)

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
