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
import selection
import assembly
import organization
import analysis
import metric
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
        "-d", "--dock", dest="dock", type=str, required=True,
        help=(
            "Path to root or dock directory for source and product " +
            "directories and files."
        )
    )
    parser_main.add_argument(
        "-a", "--access", dest="access", action="store_true",
        help="Access raw information from GTEx Portal."
    )
    parser_main.add_argument(
        "-s", "--selection", dest="selection", action="store_true",
        help="Selection of tissues and patients of interest."
    )
    parser_main.add_argument(
        "-b", "--assembly", dest="assembly", action="store_true",
        help=(
            "Collection of tissues, patients, and genes' signals for each " +
            "sample."
        )
    )
    parser_main.add_argument(
        "-o", "--organization", dest="organization", action="store_true",
        help="Organization of information for patients, tissues, samples, " +
        "and genes."
    )
    parser_main.add_argument(
        "-n", "--analysis", dest="analysis", action="store_true",
        help="Analysis of real, shuffle, and simulation data sets."
    )
    parser_main.add_argument(
        "-m", "--metric", dest="metric", action="store_true",
        help="Definition and test of metrics for modality."
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
    if arguments.selection:
        # Report status.
        print("... executing selection procedure ...")
        # Execute procedure.
        selection.execute_procedure(dock=arguments.dock)
    if arguments.assembly:
        # Report status.
        print("... executing assembly procedure ...")
        # Execute procedure.
        assembly.execute_procedure(dock=arguments.dock)
    if arguments.organization:
        # Report status.
        print("... executing organization procedure ...")
        # Execute procedure.
        organization.execute_procedure(dock=arguments.dock)
    if arguments.analysis:
        # Report status.
        print("... executing analysis procedure ...")
        # Execute procedure.
        analysis.execute_procedure(dock=arguments.dock)
    if arguments.metric:
        # Report status.
        print("... executing metric procedure ...")
        # Execute procedure.
        metric.execute_procedure(dock=arguments.dock)




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
