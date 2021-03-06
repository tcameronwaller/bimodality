"""
...
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard.
import argparse
import textwrap

# Relevant.

# Custom.

import access
import assembly
import selection

import split
import distribution
import candidacy
import stratification

import shuffle
import permutation
import probability

import heritability
import prediction
import function
import collection
import integration


import modality
import plot
#import test
import utility
#import expecto

#dir()
#importlib.reload()

###############################################################################
# Functionality


def define_interface_parsers():
    """
    Defines and parses arguments from terminal's interface.

    arguments:

    raises:

    returns:
        (object): arguments from terminal

    """

    # Define description.
    description = define_general_description()
    # Define epilog.
    epilog = define_general_epilog()
    # Define arguments.
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(title="procedures")
    parser_main = define_main_subparser(subparsers=subparsers)
    # Parse arguments.
    return parser.parse_args()


def define_general_description():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        Access data from Genotype-Tissue Expression (GTEx) and do other stuff.

        --------------------------------------------------
    """)
    return description


def define_general_epilog():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def define_main_subparser(subparsers=None):
    """
    Defines subparser for procedures that adapt a model of human metabolism.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to parser

    """

    if False:
        parser_main.add_argument(
            "-measurement", "--measurement", dest="measurement",
            action="store_true",
            help=(
                "Selection of genes and measurements of interest for further " +
                "analysis."
            )
        )
        parser_main.add_argument(
            "-sample", "--sample", dest="sample", action="store_true",
            help=(
                "Selection of samples, patients, and tissues of interest for " +
                "further analysis."
            )
        )
        parser_main.add_argument(
            "-tissue", "--tissue", dest="tissue", action="store_true",
            help=(
                "Comparison of minor and major categories of tissues."
            )
        )



    # Define description.
    description = define_main_description()
    # Define epilog.
    epilog = define_main_epilog()
    # Define parser.
    parser_main = subparsers.add_parser(
        name="main",
        description=description,
        epilog=epilog,
        help="Help for main routine.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser_main.add_argument(
        "-dock", "--dock", dest="dock", type=str, required=True,
        help=(
            "Path to root or dock directory for source and product " +
            "directories and files."
        )
    )
    parser_main.add_argument(
        "-access", "--access", dest="access", action="store_true",
        help="Access raw information from GTEx Portal."
    )
    parser_main.add_argument(
        "-assembly", "--assembly", dest="assembly", action="store_true",
        help=(
            "Preliminary assembly of relevant information."
        )
    )
    parser_main.add_argument(
        "-selection", "--selection", dest="selection", action="store_true",
        help=(
            "Coordination of selection of samples and genes of interest for " +
            "further analysis."
        )
    )
    parser_main.add_argument(
        "-split", "--split", dest="split", action="store_true",
        help="Split genes' signals across samples by genes."
    )

    parser_main.add_argument(
        "-distribution", "--distribution", dest="distribution",
        action="store_true",
        help=(
            "Analyze real and shuffle signals for a single gene."
        )
    )
    parser_main.add_argument(
        "-candidacy", "--candidacy", dest="candidacy", action="store_true",
        help=(
            "Evaluate genes' candidacy by bimodal distribution."
        )
    )
    parser_main.add_argument(
        "-stratification", "--stratification", dest="stratification",
        action="store_true",
        help=(
            "Evaluate cohort stratification by bimodal genes."
        )
    )

    parser_main.add_argument(
        "-shuffle", "--shuffle", dest="shuffle",
        action="store_true",
        help=(
            "Prepare indices of iterative random shuffles, permutations " +
            "for use in pipe procedure."
        )
    )
    parser_main.add_argument(
        "-permutation", "--permutation", dest="permutation",
        action="store_true",
        help=(
            "Applies shuffles to gene's signals."
        )
    )
    parser_main.add_argument(
        "-count", "--count", dest="count", type=int, required=False,
        help=(
            "Count of shuffles to generate in shuffle procedure."
        )
    )
    parser_main.add_argument(
        "-local", "--local", dest="local", action="store_true",
        help=(
            "Execute pipe procedure with local resources using " +
            "multiprocessing package."
        )
    )
    parser_main.add_argument(
        "-remote", "--remote", dest="remote", action="store_true",
        help=(
            "Execute pipe procedure with remote resources using Sun Grid " +
            "Engine."
        )
    )
    parser_main.add_argument(
        "-gene", "--gene", dest="gene", type=str, required=False,
        help=(
            "Identifier of a single gene for pipe procedure."
        )
    )

    parser_main.add_argument(
        "-probability", "--probability", dest="probability",
        action="store_true",
        help=(
            "Collection from pipe procedure of scores and random " +
            "distributions for each gene."
        )
    )
    parser_main.add_argument(
        "-combination", "--combination", dest="combination",
        action="store_true",
        help=(
            "Combination of multiple scores for distribution modality."
        )
    )

    parser_main.add_argument(
        "-heritability", "--heritability", dest="heritability",
        action="store_true",
        help=(
            "Collection of results from heritability analysis in GCTA."
        )
    )
    parser_main.add_argument(
        "-prediction", "--prediction", dest="prediction", action="store_true",
        help=(
            "Evaluate persons' properties as covariates of gene distributions."
        )
    )
    parser_main.add_argument(
        "-function", "--function", dest="function", action="store_true",
        help=(
            "Collect and organize sets of genes from set over " +
            "representation in DAVID."
        )
    )
    parser_main.add_argument(
        "-collection", "--collection", dest="collection", action="store_true",
        help=(
            "Collect genes that have differential expression in persons " +
            "with COVID-19."
        )
    )
    parser_main.add_argument(
        "-integration", "--integration", dest="integration",
        action="store_true",
        help=(
            "Integration of information about genes."
        )
    )

    parser_main.add_argument(
        "-category", "--category", dest="category",
        action="store_true",
        help=(
            "Analysis of groups of persons for individual genes."
        )
    )

    parser_main.add_argument(
        "-structure", "--structure", dest="structure", action="store_true",
        help="Data organization for genomic structure."
    )


    parser_main.add_argument(
        "-analysis", "--analysis", dest="analysis", action="store_true",
        help="Analysis of real, shuffle, and simulation data sets."
    )
    parser_main.add_argument(
        "-modality", "--modality", dest="modality", action="store_true",
        help="Definition and test of measures for modality."
    )
    parser_main.add_argument(
        "-plot", "--plot", dest="plot", action="store_true",
        help="Plot charts."
    )
    parser_main.add_argument(
        "-test", "--test", dest="test", action="store_true",
        help="Temporary code to test functionality."
    )

    parser_main.add_argument(
        "-expecto", "--expecto", dest="expecto", action="store_true",
        help="Access, read, and convert Expecto data."
    )

    # Define behavior.
    parser_main.set_defaults(func=evaluate_main_parameters)
    # Return parser.
    return parser_main


def define_main_description():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        Bimodality's main procedure

        Do stuff.

        --------------------------------------------------
    """)
    return description


def define_main_epilog():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        main routine

        --------------------------------------------------
        additional notes...


        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def evaluate_main_parameters(arguments):
    """
    Evaluates parameters for model procedure.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    if False:
        if arguments.measurement:
            # Report status.
            print("... executing measurement procedure ...")
            # Execute procedure.
            measurement.execute_procedure(dock=arguments.dock)
        if arguments.sample:
            # Report status.
            print("... executing sample procedure ...")
            # Execute procedure.
            sample.execute_procedure(dock=arguments.dock)
        if arguments.tissue:
            # Report status.
            print("... executing tissue procedure ...")
            # Execute procedure.
            tissue.execute_procedure(dock=arguments.dock)
    pass

    print("--------------------------------------------------")
    print("... call to main routine ...")
    # Execute procedure.
    if arguments.access:
        # Report status.
        print("... executing access procedure ...")
        # Execute procedure.
        access.execute_procedure(dock=arguments.dock)
    if arguments.assembly:
        # Report status.
        print("... executing assembly procedure ...")
        # Execute procedure.
        assembly.execute_procedure(dock=arguments.dock)
    if arguments.selection:
        # Report status.
        print("... executing selection procedure ...")
        # Execute procedure.
        selection.execute_procedure(dock=arguments.dock)
    if arguments.split:
        # Report status.
        print("... executing split procedure ...")
        # Execute procedure.
        split.execute_procedure(dock=arguments.dock)

    if arguments.distribution:
        if arguments.local:
            # Report status.
            print("... executing distribution procedure ...")
            # Execute procedure.
            distribution.execute_procedure_local(dock=arguments.dock)
        elif arguments.remote:
            # Execute procedure.
            distribution.execute_procedure_remote(
                dock=arguments.dock, gene=arguments.gene
            )
    if arguments.candidacy:
        # Report status.
        print("... executing candidacy procedure ...")
        # Execute procedure.
        candidacy.execute_procedure(dock=arguments.dock)
    if arguments.stratification:
        # Report status.
        print("... executing stratification procedure ...")
        # Execute procedure.
        stratification.execute_procedure(dock=arguments.dock)

    if arguments.shuffle:
        # Report status.
        print("... executing permutation procedure ...")
        # Execute procedure.
        shuffle.execute_procedure(dock=arguments.dock, count=arguments.count)
    if arguments.permutation:
        if arguments.local:
            # Report status.
            print("... executing permutation procedure ...")
            # Execute procedure.
            permutation.execute_procedure_local(dock=arguments.dock)
        elif arguments.remote:
            # Execute procedure.
            permutation.execute_procedure_remote(
                dock=arguments.dock, gene=arguments.gene
            )

    if arguments.probability:
        # Report status.
        print("... executing probability procedure ...")
        # Execute procedure.
        probability.execute_procedure(dock=arguments.dock)
    if arguments.combination:
        # Report status.
        print("... executing combination procedure ...")
        # Execute procedure.
        combination.execute_procedure(dock=arguments.dock)
    if arguments.integration:
        # Report status.
        print("... executing integration procedure ...")
        # Execute procedure.
        integration.execute_procedure(dock=arguments.dock)
    if arguments.category:
        # Report status.
        print("... executing category procedure ...")
        # Execute procedure.
        category.execute_procedure_local(dock=arguments.dock)
    if arguments.prediction:
        # Report status.
        print("... executing prediction procedure ...")
        # Execute procedure.
        prediction.execute_procedure(dock=arguments.dock)
    if arguments.heritability:
        # Report status.
        print("... executing heritability procedure ...")
        # Execute procedure.
        heritability.execute_procedure(dock=arguments.dock)
    if arguments.structure:
        # Report status.
        print("... executing structure procedure ...")
        # Execute procedure.
        structure.execute_procedure(dock=arguments.dock)

    if arguments.analysis:
        # Report status.
        print("... executing analysis procedure ...")
        # Execute procedure.
        analysis.execute_procedure(dock=arguments.dock)
    if arguments.function:
        # Report status.
        print("... executing function procedure ...")
        # Execute procedure.
        function.execute_procedure(dock=arguments.dock)
    if arguments.collection:
        # Report status.
        print("... executing collection procedure ...")
        # Execute procedure.
        collection.execute_procedure(dock=arguments.dock)
    if arguments.modality:
        # Report status.
        print("... executing modality procedure ...")
        # Execute procedure.
        modality.execute_procedure(dock=arguments.dock)
    if arguments.plot:
        # Report status.
        print("... executing plot procedure ...")
        # Execute procedure.
        plot.execute_procedure(dock=arguments.dock)
    if arguments.test:
        # Report status.
        print("... executing test procedure ...")
        # Execute procedure.
        test.execute_procedure(dock=arguments.dock)
    if arguments.expecto:
        # Report status.
        print("... executing expecto procedure ...")
        # Execute procedure.
        expecto.execute_procedure(dock=arguments.dock)


###############################################################################
# Procedure


def execute_procedure():
    """
    Function to execute module's main behavior.

    arguments:

    returns:

    raises:

    """

    # TODO: I want 2 separate procedures: 1. definition, 2. analysis

    # Parse arguments from terminal.
    arguments = define_interface_parsers()
    # Call the appropriate function.
    arguments.func(arguments)


if (__name__ == "__main__"):
    execute_procedure()
