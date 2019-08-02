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
import measurement
import sample
import selection
import tissue

import split

import candidacy
import category

import distribution
import organization
import restriction
import aggregation
import permutation
import collection
import combination

#import integration

if False:
    import analysis
    import function
import metric
import plot
#import test
import utility

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
        "-selection", "--selection", dest="selection", action="store_true",
        help=(
            "Coordination of selection of samples and genes of interest for " +
            "further analysis."
        )
    )
    parser_main.add_argument(
        "-tissue", "--tissue", dest="tissue", action="store_true",
        help=(
            "Comparison of minor and major categories of tissues."
        )
    )

    parser_main.add_argument(
        "-split", "--split", dest="split", action="store_true",
        help="Split genes' signals across samples by genes."
    )

    parser_main.add_argument(
        "-permutation", "--permutation", dest="permutation",
        action="store_true",
        help=(
            "Prepare indices of iterative random shuffles, permutations " +
            "for use in pipe procedure."
        )
    )
    parser_main.add_argument(
        "-count", "--count", dest="count", type=int, required=False,
        help=(
            "Count of shuffles to generate in shuffle procedure."
        )
    )

    parser_main.add_argument(
        "-candidacy", "--candidacy", dest="candidacy", action="store_true",
        help=(
            "Evaluate genes' candidacy prior to shuffles."
        )
    )

    parser_main.add_argument(
        "-distribution", "--distribution", dest="distribution",
        action="store_true",
        help=(
            "Analyze real and shuffle signals for a single gene."
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
        "-organization", "--organization", dest="organization",
        action="store_true",
        help="Organization of information for each gene about signals " +
        "across persons and tissues."
    )
    parser_main.add_argument(
        "-restriction", "--restriction", dest="restriction",
        action="store_true",
        help="Restriction of analysis for each gene to specific persons " +
        "and tissues."
    )
    parser_main.add_argument(
        "-aggregation", "--aggregation", dest="aggregation",
        action="store_true",
        help="Aggregation of gene's signals across tissues within persons."
    )

    parser_main.add_argument(
        "-collection", "--collection", dest="collection",
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
        "-category", "--category", dest="category",
        action="store_true",
        help=(
            "Analysis of groups of persons for individual genes."
        )
    )

    parser_main.add_argument(
        "-analysis", "--analysis", dest="analysis", action="store_true",
        help="Analysis of real, shuffle, and simulation data sets."
    )
    parser_main.add_argument(
        "-function", "--function", dest="function", action="store_true",
        help="Functional analyses on ranks of genes."
    )
    parser_main.add_argument(
        "-metric", "--metric", dest="metric", action="store_true",
        help="Definition and test of metrics for modality."
    )
    parser_main.add_argument(
        "-plot", "--plot", dest="plot", action="store_true",
        help="Plot charts."
    )
    parser_main.add_argument(
        "-test", "--test", dest="test", action="store_true",
        help="Temporary code to test functionality."
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
    if arguments.selection:
        # Report status.
        print("... executing selection procedure ...")
        # Execute procedure.
        selection.execute_procedure(dock=arguments.dock)
    if arguments.tissue:
        # Report status.
        print("... executing tissue procedure ...")
        # Execute procedure.
        tissue.execute_procedure(dock=arguments.dock)
    if arguments.split:
        # Report status.
        print("... executing split procedure ...")
        # Execute procedure.
        split.execute_procedure(dock=arguments.dock)

    if arguments.candidacy:
        # Report status.
        print("... executing candidacy procedure ...")
        # Execute procedure.
        candidacy.execute_procedure_local(dock=arguments.dock)
    if arguments.permutation:
        # Report status.
        print("... executing permutation procedure ...")
        # Execute procedure.
        permutation.execute_procedure(dock=arguments.dock, count=arguments.count)

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

    if arguments.organization:
        # Report status.
        print("... executing organization procedure ...")
        # Execute procedure.
        organization.test(dock=arguments.dock)
    if arguments.restriction:
        # Report status.
        print("... executing restriction procedure ...")
        # Execute procedure.
        restriction.test(dock=arguments.dock)
    if arguments.aggregation:
        # Report status.
        print("... executing aggregation procedure ...")
        # Execute procedure.
        aggregation.test(dock=arguments.dock)

    if arguments.collection:
        # Report status.
        print("... executing collection procedure ...")
        # Execute procedure.
        collection.execute_procedure(dock=arguments.dock)
    if arguments.combination:
        # Report status.
        print("... executing combination procedure ...")
        # Execute procedure.
        combination.execute_procedure(dock=arguments.dock)
    if arguments.category:
        # Report status.
        print("... executing category procedure ...")
        # Execute procedure.
        category.execute_procedure(dock=arguments.dock)

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
    if arguments.metric:
        # Report status.
        print("... executing metric procedure ...")
        # Execute procedure.
        metric.execute_procedure(dock=arguments.dock)
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
