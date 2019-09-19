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


# TODO: update to version 8 of GTEx


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
                "https://storage.googleapis.com/gtex_analysis_v8/annotations/"
            ),
            "name": ("GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx"),
            "suffix": "xlsx",
            "compression": False,
        },
        "description_person": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v8/annotations/"
            ),
            "name": ("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx"),
            "suffix": "xlsx",
            "compression": False,
        },
        "attribute_sample": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v8/annotations/"
            ),
            "name": ("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
            "suffix": "txt",
            "compression": False,
        },
        "attribute_person": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v8/annotations/"
            ),
            "name": ("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"),
            "suffix": "txt",
            "compression": False,
        },
        "count_gene": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
            ),
            "name": (
                "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
            ),
            "suffix": "gct.gz",
            "compression": True,
        },
        "signal_gene": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
            ),
            "name": (
                "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
            ),
            "suffix": "gct.gz",
            "compression": True,
        },
        "annotation_gene_gtex": {
            "path": (
                "https://storage.googleapis.com/gtex_analysis_v8/reference/"
            ),
            "name": (
                "gencode.v26.GRCh38.genes.gtf"
            ),
            "suffix": "gtf",
            "compression": False,
        },
        # Update on 19 September 2019
        # Gencode Release 31 (GRCh38.p12)
        "annotation_gene_gencode": {
            "path": (
                "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/" +
                "release_31/"
            ),
            "name": (
                "gencode.v31.annotation.gtf.gz"
            ),
            "suffix": "gtf.gz",
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
        reference=reference, key="description_person", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="attribute_sample", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="attribute_person", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="count_gene", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="signal_gene", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="annotation_gene_gtex", path_local=path_local
    )
    wrap_download_file(
        reference=reference, key="annotation_gene_gencode",
        path_local=path_local
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
        reference=reference, key="description_person", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="attribute_sample", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="attribute_person", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="count_gene", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="signal_gene", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="annotation_gene_gtex", path_local=path_local
    )
    wrap_extract_file(
        reference=reference, key="annotation_gene_gencode",
        path_local=path_local
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
