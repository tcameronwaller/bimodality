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


# Collection of COVID-19 genes


def read_source_covid_genes(
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
    path_collection = os.path.join(
        dock, "annotation", "gene_sets", "genes_covid19", "collection"
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
        designation = file_name.split(".")[0]
        identifier = str(designation.replace("genes_", ""))
        data_key = str("data_" + designation)
        print(data_key)
        # Specify directories and files.
        path_file = os.path.join(path_collection, file_name)
        # Read information from file.
        bin[identifier] = dict()
        bin[identifier]["identifier"] = identifier
        bin[identifier]["data"] = pandas.read_csv(
            path_file,
            sep="\t",
            header=0,
            low_memory=False,
        )
        pass
    # Return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_selection": genes_candidacy["any"],
        "bin_studies": bin,
    }


def collect_sort_unique_study_comparisons(
    bin_studies=None,
):
    """
    Collects and organizes unique designations of comparisons from all studies.

    arguments:
        bin_studies (dict): collection of information about each study

    raises:

    returns:
        (list<str>): unique designations of comparisons from all studies

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
    comparisons_unique = utility.collect_unique_elements(
        elements_original=collection
    )
    return sorted(comparisons_unique)


def collect_unique_gene_names(
    bin_studies=None,
):
    """
    Collects unique names of genes from all studies.

    arguments:
        bin_studies (dict): collection of information about each study

    raises:

    returns:
        (list<str>): unique names of genes from all studies

    """

    collection = []
    for study in bin_studies.keys():
        data_study = bin_studies[study]["data"]
        genes_names_study = data_study["name"].to_list()
        collection.extend(genes_names_study)
    genes_names_unique = utility.collect_unique_elements(
        elements_original=collection
    )
    return genes_names_unique


def determine_gene_study_comparison_value(
    gene=None,
    comparison=None,
    data_study_gene=None,
):
    """
    Determines whether a gene belongs to a comparison set in a study.

    arguments:
        gene (str): name of a gene
        comparison (str): name of a comparison
        data_study_gene (object): Pandas data frame of a gene's inclusion in
           comparison sets in a study

    raises:

    returns:
        (bool): whether gene belongs to the comparison set in the study

    """

    # Determine whether current study reports comparison.
    if comparison in data_study_gene.columns.to_list():
        values = data_study_gene[comparison].to_list()
        match = any(list(map(
            lambda value: (value > 0.5),
            values
        )))
    else:
        match = False
    # Return information.
    return match


def collect_gene_studies_comparisons(
    gene_name=None,
    entry_original=None,
    comparisons=None,
    bin_studies=None,
):
    """
    Collects and organizes genes that show differential expression in multiple
    studies and comparisons.

    arguments:
        gene_name (str): name of a gene
        entry_original (dict): entry of information about a gene
        comparisons (list<str>): unique designations of comparisons in all
            studies
        bin_studies (dict): collection of information about each study

    raises:

    returns:
        (dict): novel entry for gene

    """

    # Copy original entry.
    entry = copy.deepcopy(entry_original)
    entry["names"].append(gene_name)
    # Iterate across studies.
    for study in bin_studies.keys():
        data_study = bin_studies[study]["data"]
        genes_names_study = data_study["name"].to_list()
        # Determine whether study mentions the current name of the gene.
        # If a gene has multiple names, then each will have its own separate
        # consideration in this function.
        if (gene_name in genes_names_study):
            # Study mentions current gene.
            if study not in entry["studies"]:
                entry["studies"].append(study)
            # Select study's information for gene.
            data_study_gene = data_study.loc[
                data_study["name"] == gene_name, :
            ].copy(deep=True)
            data_study_gene.drop_duplicates(
                subset=None,
                keep="first",
                inplace=True,
                #ignore_index=True,
            )
            # Need to consider each comparison for the gene.
            # Iterate across comparisons.
            for comparison in comparisons:
                # Determine whether entry already has a value for the
                # comparison.
                if (comparison in entry.keys()):
                    # Determine whether record has a true value for the
                    # comparison.
                    # If value is already true, then ignore.
                    # If value is not already true, then consider.
                    if not entry[comparison]:
                        # While the gene already has a false value for the
                        # comparison, reconsider to include any true
                        # values from novel names of the gene or novel studies.
                        value = determine_gene_study_comparison_value(
                            gene=gene_name,
                            comparison=comparison,
                            data_study_gene=data_study_gene,
                        )
                        entry[comparison] = value
                        pass
                else:
                    # As the gene does not have a value for the comparison,
                    # consider.
                    value = determine_gene_study_comparison_value(
                        gene=gene_name,
                        comparison=comparison,
                        data_study_gene=data_study_gene,
                    )
                    entry[comparison] = value
                    pass
                pass
            pass
        pass
    entry["reference"] = ";".join(entry["studies"])
    # Return novel entry.
    return entry


def collect_genes_summaries(
    genes_names=None,
    comparisons=None,
    bin_studies=None,
    data_gene_annotation=None,
    dock=None,
):
    """
    Collects and organizes genes that show differential expression in multiple
    studies and comparisons.

    arguments:
        genes_names (list<str>): unique names of genes in all studies
        comparisons (list<str>): unique designations of comparisons in all
            studies
        bin_studies (dict): collection of information about each study
        data_gene_annotation (object): Pandas data frame of genes' annotations
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of genes summaries across comparisons and
            studies

    """

    # Iterate across genes.
    entries = dict()
    for gene_name in genes_names:
        # Access information about gene.
        identifier = assembly.translate_gene_name_to_identifier(
            name=gene_name,
            data_gene_annotation=data_gene_annotation,
            dock=dock,
        )
        # Determine whether an entry already exists for the gene.
        # Multiple gene names might translate to a single identifier.
        if identifier not in entries.keys():
            # An entry does not already exist for the gene.
            # Create new entry for gene.
            entry = dict()
            entry["identifier"] = identifier
            entry["name"] = gene_name
            entry["names"] = list()
            entry["studies"] = list()
            # Collect information for gene.
            entry_novel = collect_gene_studies_comparisons(
                gene_name=gene_name,
                entry_original=entry,
                comparisons=comparisons,
                bin_studies=bin_studies,
            )
        else:
            # Collect information for gene.
            entry_novel = collect_gene_studies_comparisons(
                gene_name=gene_name,
                entry_original=entries[identifier],
                comparisons=comparisons,
                bin_studies=bin_studies,
            )
        # Collect entry.
        entries_novel = dict()
        entries_novel[entry_novel["identifier"]] = entry_novel
        entries.update(entries_novel)
    # Organize data.
    records = entries.values()
    data = pandas.DataFrame(data=records)
    return data


def collect_organize_genes_summaries(
    genes_names=None,
    comparisons=None,
    bin_studies=None,
    data_gene_annotation=None,
    dock=None,
):
    """
    Collects and organizes genes that show differential expression in multiple
    studies and comparisons.

    arguments:
        genes_names (list<str>): unique names of genes in all studies
        comparisons (list<str>): unique designations of comparisons in all
            studies
        bin_studies (dict): collection of information about each study
        data_gene_annotation (object): Pandas data frame of genes' annotations
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of genes summaries across comparisons and
            studies

    """

    # Collect genes' summaries across studies and comparisons.
    data = collect_genes_summaries(
        genes_names=genes_names,
        comparisons=comparisons,
        bin_studies=bin_studies,
        data_gene_annotation=data_gene_annotation,
        dock=dock,
    )
    # Organize data.
    data["count_studies"] = data.apply(
        lambda row: len(row["studies"]),
        axis="columns",
    )
    data.drop(
        labels=["studies"],
        axis="columns",
        inplace=True
    )
    data["count_comparisons"] = data.apply(
        lambda row: len(list(
            itertools.compress(
                row.loc[comparisons].to_list(), row.loc[comparisons].to_list()
            )
        )),
        axis="columns",
    )
    data = data[[
        "identifier",
        "name",
        "names",
        *comparisons,
        "reference",
        "count_studies",
        "count_comparisons",
    ]]
    data.sort_values(
        by=["count_studies"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    data.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    return data


def collect_organize_genes_studies_comparisons(
    bin_studies=None,
    data_gene_annotation=None,
    dock=None,
    report=None,
):
    """
    Collects and organizes genes that show differential expression in multiple
    studies and comparisons.

    arguments:
        bin_studies (dict): collection of information about each study
        data_gene_annotation (object): Pandas data frame of genes' annotations
        dock (str): path to root or dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of genes summaries across comparisons and
            studies

    """

    # Collect unique comparisons from all studies.
    comparisons = collect_sort_unique_study_comparisons(
        bin_studies=bin_studies,
    )
    # Collect unique gene names.
    genes_names = collect_unique_gene_names(
        bin_studies=bin_studies,
    )
    if report:
        utility.print_terminal_partition(level=2)
        print("unique genes: " + str(len(genes_names)))
    # Collect relevance of comparisons and references to studies for each gene.
    data_genes_comparisons_studies = collect_organize_genes_summaries(
        genes_names=genes_names,
        comparisons=comparisons,
        bin_studies=bin_studies,
        data_gene_annotation=data_gene_annotation,
        dock=dock,
    )
    # Return information.
    return data_genes_comparisons_studies


def select_covid19_genes_by_studies(
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
        (list<str>): identifiers of genes

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
    data = data.loc[
        data["count_studies"] >= threshold_studies, :
    ]
    genes_covid19 = utility.collect_unique_elements(
        elements_original=data.index.to_list()
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Count of differential expression genes that match threshold."
        )
        print("threshold count of studies: " + str(threshold_studies))
        print("threshold genes DE COVID-19: " + str(data.shape[0]))
        print("DE COVID-19 genes: " + str(len(genes_covid19)))
    # Return information.
    return genes_covid19


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
    path_genes = os.path.join(
        paths["covid19"], "genes.pickle"
    )
    path_genes_text = os.path.join(
        paths["covid19"], "genes.txt"
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
    with open(path_genes, "wb") as file_product:
        pickle.dump(
            information["genes"], file_product
        )
    utility.write_file_text_list(
        elements=information["genes"],
        delimiter="\n",
        path_file=path_genes_text
    )

    pass


# Collection of cytokine genes




# Collection of integrin genes



###############################################################################
# Procedure


def read_collect_organize_report_write_covid_genes(
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
    source = read_source_covid_genes(
        dock=paths["dock"],
    )
    # Collect and organize genes from each study and comparison.
    data_genes_comparisons_studies = (
        collect_organize_genes_studies_comparisons(
            bin_studies=source["bin_studies"],
            data_gene_annotation=source["data_gene_annotation"],
            dock=paths["dock"],
            report=True,
    ))
    # Select genes that show differential expression in multiple studies.
    genes_covid19 = select_covid19_genes_by_studies(
        data_genes_comparisons_studies=data_genes_comparisons_studies,
        genes_selection=source["genes_selection"],
        threshold_studies=2,
        report=True,
    )
    # Compile information.
    bin = dict()
    bin["data_genes_comparisons_studies"] = data_genes_comparisons_studies
    bin["genes"] = genes_covid19
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
    read_collect_organize_report_write_covid_genes(
        paths=paths,
    )

    pass


if (__name__ == "__main__"):
    execute_procedure()
