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

# Relevant

import pandas
import gseapy

# Custom

import utility
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
        group (str): group of sets, either parentage or orphan
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # TODO: change group to "parentage" or "orphan", and read in the appropriate enrichment report...

    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    # Read genes sets.
    genes_candidacy = integration.read_source_genes_sets_candidacy(
        cohort="selection",
        dock=dock,
    )
    # Read tables from clusters of gene ontology sets.
    bin = dict()
    clusters = list()
    clusters.append(dict(name="defense", file="data_cluster_one.csv"))
    clusters.append(dict(name="migration", file="data_cluster_two.csv"))
    clusters.append(dict(name="proliferation", file="data_cluster_three.csv"))
    clusters.append(dict(name="inflammation", file="data_cluster_four.csv"))
    clusters.append(dict(name="cytokine", file="data_cluster_five.csv"))
    clusters.append(dict(name="apoptosis", file="data_cluster_six.csv"))
    clusters.append(dict(name="lipid", file="data_cluster_seven.csv"))
    clusters.append(dict(name="collagen", file="data_cluster_eight.csv"))
    clusters.append(dict(name="calcium", file="data_cluster_nine.csv"))
    clusters.append(dict(name="reproduction", file="data_cluster_ten.csv"))
    for cluster in clusters:
        # Specify directories and files.
        path_report = os.path.join(
            dock, "template", "ontology_gene_sets", "david",
            "multimodal_genes_755_biological_process_2020-05-07",
            "split_cluster_tables", cluster["file"]
        )
        # Read information from file.
        entry = dict()
        entry["name"] = cluster["name"]
        entry["report"] = pandas.read_csv(
            path_report,
            sep="\t",
            header=0,
            low_memory=False,
        )
        bin[cluster["name"]] = entry
    # Return information.
    return {
        "genes_candidacy": genes_candidacy,
        "cluster_reports": bin,
    }


def collect_report_ontology_parentage_orphan_genes(
    cluster_reports=None,
    genes_query=None,
    report=None,
):
    """
    Extracts information about persons.

    arguments:
        cluster_reports (dict): reports for each cluster
        genes_query (list<str>): identifiers of genes in original enrichment
            query
        report (bool): whether to print reports

    raises:

    returns:
        (dict<list<str>>): identifiers of genes in each parent set

    """

    # Collect genes.
    genes_collection = list()
    # Iterate on cluster reports.
    for key in cluster_reports.keys():
        # Organize data.
        data_report = cluster_reports[key]["report"]
        #print(cluster_reports[key]["name"])
        #print(data_report)
        data_report.rename_axis(
            index="set",
            axis="index",
            copy=False,
            inplace=True,
        )
        records = utility.convert_dataframe_to_records(data=data_report)
        # Iterate on sets within cluster.
        for record in records:
            # Extract identifiers of genes.
            genes_set_raw = record["Genes"]
            genes_set = genes_set_raw.split(", ")
            # Collect genes.
            genes_collection.extend(genes_set)
    # Collect unique genes from parent.
    genes_parentage = utility.collect_unique_elements(
        elements_original=genes_collection,
    )
    # Collect orphan genes.
    genes_orphan = list(filter(
        lambda gene: not gene in genes_parentage,
        genes_query
    ))
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("unique parentage and orphan genes")
        print("parentage genes: " + str(len(genes_parentage)))
        print("orphan genes: " + str(len(genes_orphan)))
        utility.print_terminal_partition(level=2)
    pass


def collect_ontology_enrichment_cluster_gene_sets(
    cluster_reports=None,
    report=None,
):
    """
    Extracts information about persons.

    arguments:
        cluster_reports (dict): reports for each cluster
        report (bool): whether to print reports

    raises:

    returns:
        (dict<list<str>>): identifiers of genes in each parent set

    """

    # Collect unique genes of all children sets from each parent cluster.
    parents_genes = dict()
    # Iterate on cluster reports.
    for key in cluster_reports.keys():
        # Organize data.
        name = cluster_reports[key]["name"]
        data_report = cluster_reports[key]["report"]
        data_report.rename_axis(
            index="set",
            axis="index",
            copy=False,
            inplace=True,
        )
        records = utility.convert_dataframe_to_records(data=data_report)
        # Collect genes from each children set in parent cluster.
        genes_child = list()
        # Iterate on children sets within cluster.
        for record in records:
            # Extract identifiers of genes.
            genes_set_raw = record["Genes"]
            genes_set = genes_set_raw.split(", ")
            # Collect genes.
            genes_child.extend(genes_set)
        # Collect unique genes from parent.
        genes_child_unique = utility.collect_unique_elements(
            elements_original=genes_child,
        )
        parents_genes[name] = genes_child_unique
    # Organize data.
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("unique genes from each parent set")
        for key in cluster_reports.keys():
            utility.print_terminal_partition(level=3)
            print("parent: " + cluster_reports[key]["name"])
            print(
                "count of children sets: " +
                str(cluster_reports[key]["report"].shape[0])
            )
            print(
                "count of children genes: " +
                str(len(parents_genes[cluster_reports[key]["name"]]))
            )
            #print(data_parent)
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


def write_product_genes(
    group=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        group (str): group of persons, either selection or ventilation
        information (object): information to write to file.
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_sets_genes = os.path.join(
        paths["function"], "sets_genes.pickle"
    )
    path_genes_orphan_text = os.path.join(
        paths["function"], "genes_orphan.txt"
    )
    # Write information to file.
    with open(path_sets_genes, "wb") as file_product:
        pickle.dump(information["sets_genes"], file_product)
    utility.write_file_text_list(
        elements=information["sets_genes"]["orphan"],
        delimiter="\n",
        path_file=path_genes_orphan_text
    )
    pass


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
    # Report percentage of query genes that are orphans.
    # Orphan genes do not have assignment to any set.
    collect_report_ontology_parentage_orphan_genes(
        cluster_reports=source["cluster_reports"],
        genes_query=source["genes_candidacy"]["multimodal"],
        report=True,
    )
    # Collect genes' identifiers in each parent set.
    sets_clusters = collect_ontology_enrichment_cluster_gene_sets(
        cluster_reports=source["cluster_reports"],
        report=True,
    )
    # Determine orphan genes from query set.
    sets_union = collect_union_gene_set(
        sets=sets_clusters,
    )
    print("union genes: " + str(len(sets_union["union"])))
    sets_orphan = collect_orphan_gene_set(
        sets=sets_union,
        genes_query=source["genes_candidacy"]["multimodal"],
    )
    print("orphan genes: " + str(len(sets_orphan["orphan"])))

    # Compile information.
    bin = dict()
    bin["sets_genes"] = sets_orphan
    write_product_genes(
        group="selection",
        information=bin,
        paths=paths,
    )

    pass


if (__name__ == "__main__"):
    execute_procedure()
