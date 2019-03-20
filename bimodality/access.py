"""
...
"""

###############################################################################
# Notes

# Specifications are for the Genotype-Tissue Expression (GTEx) analysis version
# 7 (GTEx V7).

# Sample notes of interest
# SAMPID
# SMATSSCR
# SMTS
# SMTSD
# SMGTC




###############################################################################
# Installation and importation

# Standard

import os
import wget
import gzip

# Relevant

# Custom

import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def define_source():
    """
    Defines and translates names of files for source data

    arguments:

    raises:

    returns:
        (dict): definitions of paths and translations of file names

    """

    # Compile and return information.
    return {
        "description_sample": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v7/annotations/"
            ),
            "name": ("GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx"),
            "suffix": "xlsx",
            "compression": False,
        },
        "description_patient": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v7/annotations/"
            ),
            "name": ("GTEx_Analysis_v7_Annotations_SubjectPhenotypesDD.xlsx"),
            "suffix": "xlsx",
            "compression": False,
        },
        "attribute_sample": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v7/annotations/"
            ),
            "name": ("GTEx_v7_Annotations_SampleAttributesDS.txt"),
            "suffix": "txt",
            "compression": False,
        },
        "attribute_patient": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v7/annotations/"
            ),
            "name": ("GTEx_v7_Annotations_SubjectPhenotypesDS.txt"),
            "suffix": "txt",
            "compression": False,
        },
        "count_gene": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/"
            ),
            "name": (
                "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz"
            ),
            "suffix": "gct.gz",
            "compression": True,
        },
        "signal_gene": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/"
            ),
            "name": (
                "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
            ),
            "suffix": "gct.gz",
            "compression": True,
        },
    }


def download_files(reference=None, path_local=None):
    """
    Downloads data files from remote server to local files.

    https://gtexportal.org/home/datasets

    arguments:
        reference (dict): paths and names of files for specific information
        path_local (str): path to file on drive

    returns:

    raises:

    """

    wrap_download_file(
        reference=reference, key="description_sample", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="description_patient", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="attribute_sample", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="attribute_patient", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="count_gene", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="signal_gene", path_local=path_local
    )
    pass


def wrap_download_file(reference=None, key=None, path_local=None):
    """
    Wraps the download specifications.

    arguments:
        reference (dict): paths and names of files for specific information
        key (str): key for reference
        path_local (str): path to file on drive

    returns:

    raises:

    """

    download_file(
        path_remote=reference[key]["path"],
        name_remote=reference[key]["name"],
        path_local=path_local,
        name_local=(key + "." + reference[key]["suffix"]),
    )
    pass


def download_file(
    path_remote=None, name_remote=None, path_local=None, name_local=None
):
    """
    Downloads a data file from remote server and saves to local file.

    https://gtexportal.org/home/datasets

    arguments:
        path_remote (str): path to file on server
        name_remote (str): name of file on server
        path_local (str): path to file on drive
        name_local (str): name of file on local drive

    returns:

    raises:

    """

    path_name_remote = path_remote + name_remote
    path_name_local = os.path.join(path_local, name_local)
    utility.remove_file(path_name_local)
    wget.download(path_name_remote, path_name_local)


def extract_files(reference=None, path_local=None):
    """
    Extracts data files from compression.

    arguments:
        reference (dict): paths and names of files for specific information
        path_local (str): path to file on drive

    returns:

    raises:

    """

    wrap_extract_file(
        reference=reference, key="description_sample", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="description_patient", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="attribute_sample", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="attribute_patient", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="count_gene", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="signal_gene", path_local=path_local
    )
    pass


def wrap_extract_file(reference=None, key=None, path_local=None):
    """
    Wraps the extraction specifications.

    arguments:
        reference (dict): paths and names of files for specific information
        key (str): key for reference
        path_local (str): path to file on drive

    returns:

    raises:

    """

    if reference[key]["compression"]:
        name_file = (key + "." + reference[key]["suffix"])
        path_file = os.path.join(path_local, name_file)
        utility.decompress_file_gzip(path_file)
    pass




def read_source(directory=None):
    """
    Reads and organizes source information from file

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_conversion = os.path.join(directory, "conversion")
    path_compartments = os.path.join(path_conversion, "compartments.pickle")
    path_processes = os.path.join(path_conversion, "processes.pickle")
    path_reactions = os.path.join(path_conversion, "reactions.pickle")
    path_metabolites = os.path.join(path_conversion, "metabolites.pickle")
    path_source = os.path.join(directory, "source")
    path_customization = os.path.join(path_source, "customization")
    path_networkx = os.path.join(
        path_conversion, "network_elements_networkx.pickle"
    )
    path_simplification_metabolites = os.path.join(
        path_customization, "simplification_metabolites.tsv"
    )
    #path_network = os.path.join(directory, "network")
    #path_nodes_reactions = os.path.join(path_network, "nodes_reactions.pickle")
    #path_nodes_metabolites = os.path.join(
    #    path_network, "nodes_metabolites.pickle"
    #)
    # Read information from file.
    with open(path_compartments, "rb") as file_source:
        compartments = pickle.load(file_source)
    with open(path_processes, "rb") as file_source:
        processes = pickle.load(file_source)
    with open(path_reactions, "rb") as file_source:
        reactions = pickle.load(file_source)
    with open(path_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    with open(path_networkx, "rb") as file_source:
        information = pickle.load(file_source)
    simplification_metabolites = utility.read_file_table(
        path_file=path_simplification_metabolites,
        names=None,
        delimiter="\t"
    )
    # Compile and return information.
    return {
        "compartments": compartments,
        "processes": processes,
        "reactions": reactions,
        "metabolites": metabolites,
        "nodes": information["nodes"],
        "links": information["links"],
        "nodes_reactions_identifiers": (
            information["nodes_reactions_identifiers"]
        ),
        "nodes_reactions": information["nodes_reactions"],
        "nodes_metabolites_identifiers": (
            information["nodes_metabolites_identifiers"]
        ),
        "nodes_metabolites": information["nodes_metabolites"],
        "simplification_metabolites": simplification_metabolites
    }


def write_product(directory=None, information=None):
    """
    Writes product information to file

    arguments:
        directory (str): directory for product files
        information (object): information to write to file

    raises:

    returns:

    """

    # Specify directories and files.
    path = os.path.join(directory, "analysis")
    utility.confirm_path_directory(path)
    path_nodes_metabolites = os.path.join(path, "nodes_metabolites.pickle")
    path_nodes_reactions_text = os.path.join(path, "nodes_reactions.tsv")
    path_nodes_metabolites_text = os.path.join(path, "nodes_metabolites.tsv")
    path_network_reactions = os.path.join(path, "network_reactions.tsv")
    path_network_metabolites = os.path.join(path, "network_metabolites.tsv")
    path_simplification_metabolites = os.path.join(
        path, "simplification_metabolites.tsv"
    )
    # Write information to file.
    with open(path_nodes_metabolites, "wb") as file_product:
        pickle.dump(information["nodes_metabolites"], file_product)
    utility.write_file_table(
        information=information["nodes_reactions_text"],
        path_file=path_nodes_reactions_text,
        names=information["nodes_reactions_text"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["nodes_metabolites_text"],
        path_file=path_nodes_metabolites_text,
        names=information["nodes_metabolites_text"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["network_reactions"],
        path_file=path_network_reactions,
        names=information["network_reactions"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["network_metabolites"],
        path_file=path_network_metabolites,
        names=information["network_metabolites"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["simplification_metabolites"],
        path_file=path_simplification_metabolites,
        #names=information["simplification_metabolites"][0].keys(),
        names=[
            "metabolite", "name", "compartment", "omission", "replication",
            "default", "category",
            "degree_total", "degree_general", "degree_compartmental",
            "note"
        ],
        delimiter="\t"
    )


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

    # Define path.

    directories = ["access"]
    #path_local = os.path.join(os.sep, *directories)
    path_local = os.path.join(dock, *directories)
    utility.create_directory(path_local)
    # Define reference.
    reference = define_source()
    # Download files.
    download_files(reference=reference, path_local=path_local)
    # Extract files.
    extract_files(reference=reference, path_local=path_local)

    pass


if (__name__ == "__main__"):
    execute_procedure()
