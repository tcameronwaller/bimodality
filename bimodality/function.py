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

# Relevant

import pandas
import gseapy

# Custom

import utility
import prediction

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
    paths["function"] = os.path.join(paths["dock"], "function")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=paths["function"])
    utility.create_directory(path=paths["function"])
    # Return information.
    return paths


def read_source(
    group=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        group (str): group of persons, either selection or ventilation
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
    path_data_enrichment_report = os.path.join(
        dock, "template", "gene_ontology_sets",
        "multimodal_genes_1011_biological_process", "webgestalt_2020-05-05",
        "enrichment_results_wg_result1588718242.txt"
    )
    path_data_child_parent = os.path.join(
        dock, "template", "annotation_2020-05-05", "ontology",
        "data_gene_ontology_child_parent.csv"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_enrichment_report = pandas.read_csv(
        path_data_enrichment_report,
        sep="\t",
        header=0,
        low_memory=False,
    )
    data_child_parent = pandas.read_csv(
        path_data_child_parent,
        sep="\t",
        header=0,
        low_memory=False,
    )
    # Read genes sets.
    genes_selection = prediction.read_source_genes_sets_selection(
        group=group,
        dock=dock,
    )
    # Return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_selection": genes_selection,
        "data_enrichment_report": data_enrichment_report,
        "data_child_parent": data_child_parent,
    }



def collect_parent_child_enrichment_gene_sets(
    data_enrichment_report=None,
    data_child_parent=None,
    report=None,
):
    """
    Extracts information about persons.

    arguments:
        data_enrichment_report (object): Pandas data frame of genes'
            identifiers in each child set
        data_child_parent (object): Pandas data frame of child sets that belong
            to each parent set
        report (bool): whether to print reports

    raises:

    returns:
        (dict<list<str>>): identifiers of genes in each parent set

    """

    # Organize data.
    data_child_parent.reset_index(
        level=None,
        inplace=True
    )
    data_child_parent.set_index(
        ["identifier_parent"],
        append=False,
        drop=True,
        inplace=True
    )
    data_enrichment_report.reset_index(
        level=None,
        inplace=True
    )
    data_enrichment_report.set_index(
        ["geneSet"],
        append=False,
        drop=True,
        inplace=True
    )
    # Split data by parent.
    groups = data_child_parent.groupby(
        level=["identifier_parent"],
    )
    # Collect unique genes of all children from each parent.
    parents_genes = dict()
    for name, data_parent in groups:
        # Extract identifiers of children from the parent.
        identifiers_children = data_parent["identifier_child"].to_list()
        # Collect genes from each child.
        children_genes = list()
        for child in identifiers_children:
            if child in data_enrichment_report.index.to_list():
                # Extract identifiers of genes.
                genes_child_raw = data_enrichment_report.at[child, "userId"]
                genes_child = genes_child_raw.split(";")
                # Collect genes.
                children_genes.extend(genes_child)
        # Collect unique genes from parent.
        children_genes_unique = utility.collect_unique_elements(
            elements_original=children_genes,
        )
        parents_genes[name] = children_genes_unique
    # Organize data.
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("unique genes from each parent set")
        for name, data_parent in groups:
            utility.print_terminal_partition(level=3)
            print("parent: " + name)
            print("count of children sets: " + str(data_parent.shape[0]))
            print("count of children genes: " + str(len(parents_genes[name])))
            print(data_parent)
        utility.print_terminal_partition(level=2)
    # Return information.
    return parents_genes


def collect_union_gene_set(
    sets=None,
):
    """
    Collects the union of elements from multiple sets.

    arguments:
        sets (dict<dict<list<str>>>): sets of genes

    raises:

    returns:
        (dict<dict<list<str>>>): sets of genes

    """

    sets = copy.deepcopy(sets)
    union = list()
    for set in sets.keys():
        union.extend(sets[set])
    sets["union"] = utility.collect_unique_elements(
        elements_original=union
    )
    # Return information.
    return sets


def collect_orphan_gene_set(
    sets=None,
    genes_query=None,
):
    """
    Collects the union of elements from multiple sets.

    arguments:
        sets (dict<dict<list<str>>>): sets of genes
        genes_query (list<str>): identifiers of genes in original enrichment
            query

    raises:

    returns:
        (dict<dict<list<str>>>): sets of genes

    """

    sets = copy.deepcopy(sets)
    orphans_raw = list(filter(
        lambda gene: not gene in sets["union"],
        genes_query
    ))
    sets["orphan"] = utility.collect_unique_elements(
        elements_original=orphans_raw
    )
    # Return information.
    return sets



###############################################################################
# Procedure


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
    # Read source information from file.
    source = read_source(
        group="selection",
        dock=paths["dock"],
    )

    print(source["data_enrichment_report"])
    print(source["data_child_parent"])
    # Collect genes' identifiers in each parent set.
    sets_enrichment = collect_parent_child_enrichment_gene_sets(
        data_enrichment_report=source["data_enrichment_report"],
        data_child_parent=source["data_child_parent"],
        report=True,
    )
    # Determine orphan genes from query set.
    sets_union = collect_union_gene_set(
        sets=sets_enrichment,
    )
    print("union genes: " + str(len(sets_union["union"])))
    sets_orphan = collect_orphan_gene_set(
        sets=sets_union,
        genes_query=source["genes_selection"]["multimodal"],
    )
    print("orphan genes: " + str(len(sets_orphan["orphan"])))

    # TODO: now create "sets_genes" by consulting genes_multimodal to determine orphan genes without any parents


    pass


if (__name__ == "__main__"):
    execute_procedure()
