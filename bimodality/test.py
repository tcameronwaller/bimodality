
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
#import cmapPy.pandasGEXpress.parse
#import cmapPy.pandasGEXpress.gct2gctx

import gtfparse

# Custom

import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


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

    print("Hello beautiful world!")

    path_access = os.path.join(dock, "access")
    path_map = os.path.join(path_access, "annotation_gene_gencode.gtf")
    utility.print_file_lines(path_file=path_map, start=0, stop=10)


    # TODO: need to remove version number from Ensembl ids...


    pass


if (__name__ == "__main__"):
    execute_procedure()
