"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os

# Relevant

import pandas
import gseapy

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
    path_analysis = os.path.join(dock, "analysis")
    path_summary = os.path.join(
        path_analysis, "data_summary_genes.pickle"
    )
    path_rank_genes = os.path.join(
        path_analysis, "data_rank_genes.txt"
    )
    path_genes = os.path.join(
        path_analysis, "genes.txt"
    )

    # Read information from file.
    data_summary_genes = pandas.read_pickle(path_summary)
    data_rank_genes = pandas.read_csv(
        path_rank_genes,
        sep="\t",
        header=None,
        #nrows=1000,
    )
    genes = utility.read_file_text_list(path_genes)
    # Compile and return information.
    return {
        "data_summary_genes": data_summary_genes,
        "data_rank_genes": data_rank_genes,
        "genes": genes,
    }


# Write.


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
    path_analysis = os.path.join(dock, "analysis")
    utility.create_directory(path_analysis)
    path_probabilities = os.path.join(
        path_analysis, "genes_probabilities.pickle"
    )
    path_summary = os.path.join(
        path_analysis, "data_summary_genes.pickle"
    )
    path_summary_text = os.path.join(
        path_analysis, "data_summary_genes.tsv"
    )
    path_summary_text_alternative = os.path.join(
        path_analysis, "data_summary_genes_alternative.tsv"
    )
    path_rank_genes = os.path.join(
        path_analysis, "data_rank_genes.rnk"
    )
    path_genes = os.path.join(
        path_analysis, "genes.txt"
    )
    # Write information to file.
    with open(path_probabilities, "wb") as file_product:
        pickle.dump(information["genes_probabilities"], file_product)
    information["data_summary_genes"].to_pickle(path_summary)
    information["data_summary_genes"].to_csv(
        path_or_buf=path_summary_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_summary_genes"].reset_index(
        level="identifier", inplace=True
    )
    summary_genes = utility.convert_dataframe_to_records(
        data=information["data_summary_genes"]
    )
    utility.write_file_text_table(
        information=summary_genes,
        path_file=path_summary_text_alternative,
        names=summary_genes[0].keys(),
        delimiter="\t",
        header=True,
    )
    information["data_rank_genes"].to_csv(
        path_or_buf=path_rank_genes,
        sep="\t",
        header=False,
        index=False,
    )
    utility.write_file_text_list(
        elements=information["genes"],
        delimiter="\n",
        path_file=path_genes
    )
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
    #source["data_summary_genes"]
    #source["data_rank_genes"]
    #source["genes"]

    result = gseapy.prerank(
        rnk=source["data_rank_genes"],
        gene_sets="MSigDB_Oncogenic_Signatures",
        processes=7,
        permutation_num=100,
        outdir=dock,
        format="png"
    )
    print(result)
    print(result.res2d)
    names = gseapy.get_library_name()
    for name in names:
        print(name)



    if False:
        # Compile information.
        information = {
        }
        #Write product information to file.
        write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
