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
import copy
import random
import itertools
import functools
import multiprocessing
import datetime
import time


# Relevant

import numpy
import pandas
import scipy.stats
#from sklearn.linear_model import LinearRegression
import statsmodels.api
import statsmodels.stats.outliers_influence

# Custom

import integration
import selection
import distribution
import modality
import utility
import assembly

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
    paths["prediction"] = os.path.join(paths["dock"], "prediction")
    # Remove previous files to avoid version or batch confusion.
    #utility.remove_directory(path=paths["prediction"])
    #utility.create_directory(path=paths["prediction"])
    # Define paths for cohorts of persons.
    cohorts = list()
    cohorts.append("selection")
    cohorts.append("respiration")
    cohorts.append("ventilation")
    for cohort in cohorts:
        paths[cohort] = dict()
        # Define path for comprehensive information from regression on model.
        paths[cohort]["regression"] = os.path.join(
            paths["prediction"], cohort, "regression"
        )
        # Define paths for groups of genes by their distributions.
        paths[cohort]["genes"] = os.path.join(
            paths["prediction"], cohort, "genes"
        )
        paths[cohort]["distribution"] = dict()
        groups = list()
        groups.append("any")
        groups.append("multimodal")
        groups.append("nonmultimodal")
        groups.append("unimodal")
        for group in groups:
            paths[cohort]["distribution"][group] = os.path.join(
                paths["prediction"], cohort, "genes", "distribution", group
            )
            # Initialize directories.
            utility.create_directories(path=paths[cohort]["regression"])
            utility.create_directories(
                path=paths[cohort]["distribution"][group]
            )
    # Return information.
    return paths


def read_source(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
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
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", cohort,
        "data_persons_properties.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    # Read genes sets.
    genes_candidacy = integration.read_source_genes_sets_candidacy(
        cohort=cohort,
        dock=dock,
    )
    genes_heritability = integration.read_source_genes_sets_heritability(
        cohort="selection",
        dock=dock,
    )
    # Return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "genes_candidacy": genes_candidacy,
        "genes_heritability": genes_heritability,
        "data_signals_genes_persons": data_signals_genes_persons,
    }


def read_source_regression(
    cohort=None,
    paths=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_regression_genes = os.path.join(
        paths[cohort]["regression"], "data_regression_genes.pickle"
    )
    # Read information from file.
    data_regression_genes = pandas.read_pickle(
        path_data_regression_genes
    )
    # Compile and return information.
    return {
        "data_regression_genes": data_regression_genes,
    }


##########
# Organization


def calculate_report_independent_variable_correlations(
    pairs=None,
    data_persons_properties=None,
    method=None,
):
    """
    Reports the correlation coefficients between pairs of independent
    variables.

    arguments:
        pairs (list<tuple<str>>): pairs of names of independent variables
        data_persons_properties (object): Pandas data frame of persons and
            their relevant properties
        method (str): method for correlation, pearson, spearman, or kendall

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)

    # Organize data.
    data = data_persons_properties.copy(deep=True)

    # Iterate on pairs of independent variables.
    for pair in pairs:
        # Select data for pair of features.
        data_pair = data.loc[:, list(pair)]
        # Remove observations with missing values for either feature.
        data_pair.dropna(
            axis="index",
            how="any",
            inplace=True,
        )
        #print(data_pair)
        # Determine observations for which pair of features has matching
        # values.
        count = data_pair.shape[0]
        # Calculate correlation.
        if count > 1:
            if method == "pearson":
                correlation, probability = scipy.stats.pearsonr(
                    data_pair[pair[0]].to_numpy(),
                    data_pair[pair[1]].to_numpy(),
                )
            elif method == "spearman":
                correlation, probability = scipy.stats.spearmanr(
                    data_pair[pair[0]].to_numpy(),
                    data_pair[pair[1]].to_numpy(),
                )
            elif method == "kendall":
                correlation, probability = scipy.stats.kendalltau(
                    data_pair[pair[0]].to_numpy(),
                    data_pair[pair[1]].to_numpy(),
                )
            pass
        else:
            correlation = float("nan")
            probability = float("nan")
            pass
        # Report correlation.
        utility.print_terminal_partition(level=3)
        print(pair[0] + " ... versus ... " + pair[1])
        print("correlation: " + str(round(correlation, 6)))
        print("probability: " + str(probability))

    utility.print_terminal_partition(level=2)
    pass


def organize_dependent_independent_variables(
    variables=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
    report=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        variables (list<str>): names of independent variables for regression
        data_persons_properties (object): Pandas data frame of persons and
            their relevant properties
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of dependent and independent variables

    """

    # Remove all columns from persons properties except the covariates
    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    data_signals_genes_persons = data_signals_genes_persons.copy(
        deep=True
    )
    # Organize data.
    data_persons_properties = data_persons_properties.loc[
        :, data_persons_properties.columns.isin(variables)
    ]

    # Join signals on persons' properties.
    data_variables = data_persons_properties.join(
        data_signals_genes_persons,
        how="left",
        on="person"
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=1)
        print("data_variables")
        utility.print_terminal_partition(level=2)
        print(data_variables)
        pass
    # Return information.
    return data_variables


##########
# Regression


def regress_signal_ordinary_residuals(
    dependence=None,
    independence=None,
    proportion=None,
    data=None,
):
    """
    Regresses a quantitative continuous dependent variable against multiple
    independent variables and returns relevant parameters and statistics.

    Data format must have observations across rows and dependent and
    independent variables (features) across columns.

    Description of formats for StatsModels...

    Format of dependent variable is a vector of scalar values.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of independent variable(s) is a matrix: a first dimension vector of
    observations and for each observation a second dimension vector of scalar
    values.
    StatsModels also requires a constant for the intercept.
    [
        [1.3, 5.2, 1.0],
        [1.5, 5.1, 1.0],
        [1.2, 5.5, 1.0],
        ...
    ]

    Description of formats for SKLearn...

    Format of dependent variable is an array of observations, and each
    observation is an array of features.
    [
        [1.3],
        [1.5],
        [1.2],
        ...
    ]

    arguments:
        dependence (str): name of dependent variable
        independence (list<str>): names of independent variables
        proportion (float): minimal proportion of total observations that must
            have nonmissing values
        data (object): Pandas data frame of values and associations between
            dependent and independent variables

    raises:

    returns:
        (dict): collection of residuals for regression summary of information
            about regression model


    """

    if False:
        # Regress by SciKit Learn.
        regression = LinearRegression(
            fit_intercept=True,
            normalize=False,
            copy_X=True,
            n_jobs=None
        ).fit(
            values_independence,
            values_dependence,
            sample_weight=None
        )
        score = regression.score(
            values_independence,
            values_dependence,
            sample_weight=None
        )

    # Organize data.
    data = data.copy(deep=True)
    # Determine minimal count of observations.
    count = data.shape[0]
    threshold = proportion * count
    # Remove observations with any missing values.
    columns = copy.deepcopy(independence)
    columns.insert(0, dependence)
    data = data.loc[
        :, data.columns.isin(columns)
    ]
    data.dropna(
        axis="index",
        how="any",
        #subset=[dependence],
        inplace=True,
    )
    data = data[[*columns]]

    # Note
    # It is very important to keep track of the order of independent variables
    # in order to match parameters and probabilities to the correct variables.

    # Determine whether data have sufficient observations for regression.
    if data.shape[0] > threshold:
        # Extract values of dependent and independent variables.
        values_dependence = data[dependence].to_numpy()
        values_independence = data.loc[ :, independence].to_numpy()
        # Introduce constant value for intercept.
        values_independence_intercept = statsmodels.api.add_constant(
            values_independence,
            prepend=True,
        )
        # Define model.
        model = statsmodels.api.OLS(
            values_dependence,
            values_independence_intercept,
            missing="drop",
        )
        report = model.fit()
        #utility.print_terminal_partition(level=2)
        #print(data)
        #print(independence)
        #print(values_independence)
        #print(report.summary())
        #utility.print_terminal_partition(level=3)
        #print(dir(report))
        #print(report.params)
        #print(report.pvalues)

        # Organize residuals.
        residuals = report.resid

        # Compile summary information.
        counter = 1
        parameters = dict()
        probabilities = dict()
        inflations = dict()
        parameters["intercept_parameter"] = report.params[0]
        probabilities["intercept_probability"] = report.pvalues[0]
        inflations["intercept_inflation"] = float("nan")
        # Iterate on each independent variable.
        # Accommodate index for intercept.
        for variable in independence:
            # Coefficient or parameter.
            parameter = str(variable + ("_parameter"))
            parameters[parameter] = report.params[counter]
            # Probability.
            probability = str(variable + ("_probability"))
            probabilities[probability] = report.pvalues[counter]
            # Variance Inflation Factor (VIF).
            inflation = str(variable + ("_inflation"))
            inflations[inflation] = (
                statsmodels.stats.outliers_influence.variance_inflation_factor(
                    values_independence_intercept,
                    counter
                )
            )
            # Increment index.
            counter += 1
            pass
        summary = {
            "freedom": report.df_model,
            "observations": report.nobs,
            "r_square": report.rsquared,
            "r_square_adjust": report.rsquared_adj,
            "log_likelihood": report.llf,
            "akaike": report.aic,
            "bayes": report.bic,
            "condition": report.condition_number,
        }
        summary.update(parameters)
        summary.update(probabilities)
        summary.update(inflations)
    else:
        # Compile information.
        #probabilities = list(
        #    itertools.repeat(float("nan"), len(values_independence_intercept))
        #)
        parameters = dict()
        probabilities = dict()
        inflations = dict()
        parameters["intercept_parameter"] = float("nan")
        probabilities["intercept_probability"] = float("nan")
        inflations["intercept_inflation"] = float("nan")
        for variable in independence:
            parameter = str(variable + ("_parameter"))
            parameters[parameter] = float("nan")
            probability = str(variable + ("_probability"))
            probabilities[probability] = float("nan")
            inflation = str(variable + ("_inflation"))
            inflations[inflation] = float("nan")
            pass
        summary = {
            "freedom": float("nan"),
            "observations": float("nan"),
            "r_square": float("nan"),
            "r_square_adjust": float("nan"),
            "log_likelihood": float("nan"),
            "akaike": float("nan"),
            "bayes": float("nan"),
            "condition": float("nan"),
        }
        summary.update(parameters)
        summary.update(probabilities)
        summary.update(inflations)
        residuals = numpy.empty(0)
    # Compile information.
    bin = dict()
    bin["report"] = summary
    bin["residuals"] = residuals
    # Return information.
    return bin


def regress_cases(
    cases=None,
    variables=None,
    data_variables=None,
    data_gene_annotation=None,
    report=None,
):
    """
    Drives iterative regression across cases.

    arguments:
        cases (list<str>): names of independent cases for which to regress
            variables
        variables (list<str>): names of independent variables for regression
        data_variables (object): Pandas data frame of dependent and independent
            variables
        data_gene_annotation (object): Pandas data frame of genes' annotations
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of residuals for regression on each case and data
            summary on regressions

    """

    # TODO: I need to collect residuals...

    # Regress.
    count_total = len(cases)
    records = list()
    residuals_genes = dict()
    counter = 0
    for case in cases:
        # Report.
        counter += 1
        if report:
            #print("count " + str(counter) + ": " + str(case))
            percentage = (counter / count_total) * 100
            if (percentage % 10 == 0):
                print("complete cases: " + str(round(percentage)) + "%")
        bin_case = regress_signal_ordinary_residuals(
            dependence=case,
            independence=variables,
            proportion=0.5,
            data=data_variables,
        )
        # Collect residuals.
        residuals_genes[case] = bin_case["residuals"]
        # Collect reports.
        report = bin_case["report"]
        report["case"] = case
        report["name"] = assembly.access_gene_name(
            identifier=case,
            data_gene_annotation=data_gene_annotation,
        )
        records.append(report)
        pass
    # Organize data.
    data_regression = utility.convert_records_to_dataframe(
        records=records
    )
    data_regression.set_index(
        ["case"],
        append=False,
        drop=True,
        inplace=True
    )
    columns = list()
    columns.append("name")
    columns.append("observations")
    columns.append("freedom")
    columns.append("r_square")
    columns.append("r_square_adjust")
    columns.append("log_likelihood")
    columns.append("akaike")
    columns.append("bayes")
    columns.append("condition")
    columns.append("intercept_parameter")
    columns.append("intercept_probability")
    for variable in variables:
        columns.append(str(variable + ("_parameter")))
        columns.append(str(variable + ("_probability")))
        columns.append(str(variable + ("_inflation")))
        pass
    data_regression = data_regression[[*columns]]
    # Compile information.
    bin = dict()
    bin["residuals_genes"] = residuals_genes
    bin["data_regression_genes"] = data_regression
    # Return information.
    return bin


def report_regression_models_quality(
    variables=None,
    data_regression_models=None,
):
    """
    Drives iterative regression across cases.

    arguments:
        cases (list<str>): names of independent cases for which to regress
            variables
        variables (list<str>): names of independent variables for regression
        data_regression_models (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:

    """

    # Copy data.
    data = data_regression_models.copy(deep=True)

    # Define columns of interest.
    columns = list()
    columns.append("observations")
    columns.append("freedom")
    columns.append("r_square")
    columns.append("r_square_adjust")
    columns.append("log_likelihood")
    columns.append("akaike")
    columns.append("bayes")
    columns.append("condition")
    #columns.append("intercept_inflation")
    for variable in variables:
        columns.append(str(variable + ("_inflation")))
        pass
    # Select and sort columns of interest.
    data = data.loc[
        :, data.columns.isin(columns)
    ]
    data = data[[*columns]]
    # Aggregate data by mean.
    series = data.aggregate(
        lambda x: x.mean()
    )
    # Report information.
    utility.print_terminal_partition(level=2)
    print("Mean statistics for regression models.")
    print(series)
    utility.print_terminal_partition(level=2)
    pass


def calculate_regression_discoveries(
    variables=None,
    threshold=None,
    data=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variables (list<str>): names of independent regression variables for
            which to calculate false discovery rates
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries
        data (object): Pandas data frame of parameters and statistics from
            regressions across cases

    raises:

    returns:
        (object): Pandas data frame of parameters and statistics from
            regressions across cases

    """

    # Iterate on relevant variables.
    for variable in variables:
        # Calculate false discovery rate from probability.
        data = utility.calculate_false_discovery_rate(
            threshold=threshold,
            probability=str(variable + ("_probability")),
            discovery=str(variable + ("_discovery")),
            significance=str(variable + ("_significance")),
            data_probabilities=data,
        )
        pass
    # Return information.
    return data


def organize_data_regress_cases_report(
    genes=None,
    variables_regression=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
    data_gene_annotation=None,
    threshold_discovery=None,
    discovery=None,
    report=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        genes (list<str>): identifiers of genes for regression
        variables_regression (list<str>): names of independent variables for
            regression
        data_persons_properties (object): Pandas data frame of persons and
            their relevant properties
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        data_gene_annotation (object): Pandas data frame of genes' annotations
        threshold_discovery (float): threshold by false discovery rate
        discovery (bool): whether to calculate false discovery rates across all
            genes
        report (bool): whether to print reports

    raises:

    returns:

    """

    if report:
        # Calculate correlation coefficients between pairs of independent
        # variables.
        pairs = list()
        pairs.append(("tissues_1", "facilities_1"))
        pairs.append(("tissues_1", "facilities_2"))
        pairs.append(("delay_scale", "facilities_1"))
        pairs.append(("delay_scale", "facilities_2"))
        calculate_report_independent_variable_correlations(
            pairs=pairs,
            data_persons_properties=data_persons_properties,
            method="spearman",
        )

    if report:
        utility.print_terminal_partition(level=3)
        print("genes for regression: " + str(len(genes)))
        utility.print_terminal_partition(level=3)

    # Organize dependent and independent variables for regression analysis.
    data_variables = organize_dependent_independent_variables(
        variables=variables_regression,
        data_persons_properties=data_persons_properties,
        data_signals_genes_persons=data_signals_genes_persons,
        report=report,
    )

    # Regress each gene's signal across persons.
    # Iterate on genes.
    bin_regression = regress_cases(
        cases=genes,
        variables=variables_regression,
        data_variables=data_variables,
        data_gene_annotation=data_gene_annotation,
        report=True,
    )
    # Summarize regression quality.
    # Report means of statistics across independent models.
    if report:
        report_regression_models_quality(
            variables=variables_regression,
            data_regression_models=bin_regression["data_regression_genes"],
        )
        pass

    # Review.
    # 4 February 2020
    # Model statistics match in the individual regression's summary and in the
    # compilation across cases.
    # Independent variables' parameters and probabilities match in the
    # individual regression's summary and in the compilation across cases.

    if discovery:
        # Adjust probabilities for multiple hypothesis tests.
        # Calculate false discovery rate by Benjamini-Hochberg's method.
        data_regression_genes = calculate_regression_discoveries(
            variables=variables_regression,
            threshold=threshold_discovery,
            data=bin_regression["data_regression_genes"],
        )
    else:
        data_regression_genes = bin_regression["data_regression_genes"]
        pass
    # Compile information.
    bin = dict()
    bin["residuals_genes"] = bin_regression["residuals_genes"]
    bin["data_regression_genes"] = data_regression_genes
    # Return information.
    return bin


##########
# Sets


def select_genes_by_association_variables_separate(
    variables=None,
    coefficient_sign=None,
    count_selection=None,
    data_regression_genes=None,
):
    """
    Selects genes that associate to multiple variables separately.

    arguments:
        variables (list<str>): names of independent regression
            variables for which to select genes by significant association
        coefficient_sign (str): option to select only genes with negative,
            positive, or any sign of their coefficients for the variable
        count_selection (int): count of genes with greatest and least values
            of variable's regression parameter to keep
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:
        (dict<list<str>>): identifiers of genes that associate significantly to
            each variable in regression

    """

    # Collect sets of genes for each variables separately.
    sets = dict()
    for variable in variables:
        significance = str(variable + ("_significance"))
        parameter = str(variable + ("_parameter"))
        data_copy = data_regression_genes.copy(deep=True)
        data_copy = data_copy.loc[
            :, data_copy.columns.isin([significance, parameter])
        ]
        # Select regression data by significance of the variable.
        data_significance = data_copy.loc[data_copy[significance], :]
        # Select regression data by sign of the variable's regression
        # coefficient or parameter.
        if coefficient_sign == "negative":
            data_sign = data_significance.loc[
                (data_significance[parameter] < 0), :
            ]
        elif coefficient_sign == "positive":
            data_sign = data_significance.loc[
                (data_significance[parameter] > 0), :
            ]
        elif coefficient_sign == "any":
            data_sign = data_significance.copy(deep=True)
            pass
        # Select regression data by rankings on the variable's parameter.
        data_descent = data_sign.sort_values(
            by=[parameter],
            axis="index",
            ascending=False,
            inplace=False,
        )
        genes_greatest = data_descent.iloc[
            0:count_selection, :
        ].index.to_list()
        data_ascent = data_sign.sort_values(
            by=[parameter],
            axis="index",
            ascending=True,
            inplace=False,
        )
        genes_least = data_ascent.iloc[
            0:count_selection, :
        ].index.to_list()
        # Compile unique identifiers of genes.
        if coefficient_sign == "any":
            genes = copy.deepcopy(genes_greatest)
            genes.extend(genes_least)
        elif coefficient_sign == "negative":
            genes = copy.deepcopy(genes_least)
        elif coefficient_sign == "positive":
            genes = copy.deepcopy(genes_greatest)
            pass
        sets[variable] = utility.collect_unique_elements(
            elements_original=genes
        )
    # Return information.
    return sets


def select_genes_by_association_variables_together(
    variables=None,
    data_regression_genes=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variables (list<str>): names of independent regression
            variables for which to select genes by significant association
        count_selection (int): count of genes with greatest and least values
            of variable's regression parameter to keep
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:
        (list<str>): identifiers of genes that associate simultaneously to
            multiple variables

    """

    # Copy data.
    data = data_regression_genes.copy(deep=True)
    # Collect sets of genes for each variables separately.
    columns = list()
    for variable in variables:
        significance = str(variable + ("_significance"))
        columns.append(significance)
        pass
    # Select columns.
    data_columns = data.loc[
        :, data.columns.isin(columns)
    ]
    # Select genes with true significances across all variables.
    series_all = data_columns.all(
        axis="columns",
        bool_only=True,
    )
    series_true = series_all.loc[series_all]
    genes = series_true.index.to_list()
    # Return information.
    return genes


def select_genes_by_association_regression_variables(
    variables_separate=None,
    variables_together=None,
    threshold_r_square=None,
    threshold_discovery=None,
    coefficient_sign=None,
    count_selection=None,
    genes=None,
    data_regression_genes=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variables_separate (list<str>): names of independent regression
            variables for which to select genes by significant association
        variables_together (list<str>): names of multiple variables for which
            to select genes that associate significantly with all together
        threshold_r_square (float): threshold by r square
        threshold_discovery (float): threshold by false discovery rate
        coefficient_sign (str): option to select only genes with negative,
            positive, or any sign of their coefficients for the variable
        count_selection (int): count of genes with greatest and least values
            of variable's regression parameter to keep
        genes (list<str>): identifiers of genes
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:
        (dict<list<str>>): identifiers of genes that associate significantly to
            each variable in regression

    """

    # Copy data.
    data_regression = data_regression_genes.copy(deep=True)
    # Select regression data for genes of interest.
    data_genes = data_regression.loc[data_regression.index.isin(genes), :]
    # Select regression data by threshold on r square.
    data_threshold_r = data_genes.loc[
        data_genes["r_square"] >= threshold_r_square, :
    ]
    # Adjust probabilities for multiple hypothesis tests.
    # Calculate false discovery rate by Benjamini-Hochberg's method.
    data_discovery = calculate_regression_discoveries(
        variables=variables_separate,
        threshold=threshold_discovery,
        data=data_threshold_r,
    )
    # Select genes that associate significantly with each variable of interest
    # separately.
    # Any, or logic.
    sets = select_genes_by_association_variables_separate(
        variables=variables_separate,
        coefficient_sign=coefficient_sign,
        count_selection=count_selection,
        data_regression_genes=data_discovery,
    )
    # Select genes that associate significantly with multiple query variables
    # together.
    # All, and logic.
    set_together = select_genes_by_association_variables_together(
        variables=variables_together,
        data_regression_genes=data_discovery,
    )
    # Compile information.
    sets.update(dict(query=set_together))
    # Return information.
    return sets


def select_genes_by_gene_distributions_variable_associations(
    variables_separate=None,
    variables_together=None,
    threshold_r_square=None,
    threshold_discovery=None,
    coefficient_sign=None,
    count_selection=None,
    genes_sets=None,
    data_regression_genes=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variables_separate (list<str>): names of independent regression
            variables for which to select genes by significant association
        variables_together (list<str>): names of multiple variables for which
            to select genes that associate significantly with all together
        threshold_r_square (float): threshold by r square
        threshold_discovery (float): threshold by false discovery rate
        coefficient_sign (str): option to select only genes with negative,
            positive, or any sign of their coefficients for the variable
        count_selection (int): count of genes with greatest and least values
            of variable's regression parameter to keep
        genes_sets (dict<list<str>>): identifiers of all genes in groups by
            distributions of their pan-tissue signals
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:
        (dict<dict<list<str>>>): sets of genes that associate significantly to
            variables in regression

    """

    # Select genes from each cohort that associate significantly with each
    # variable.
    sets = dict()
    for set in genes_sets.keys():
        sets[set] = select_genes_by_association_regression_variables(
            variables_separate=variables_separate,
            variables_together=variables_together,
            threshold_r_square=threshold_r_square,
            threshold_discovery=threshold_discovery,
            coefficient_sign=coefficient_sign,
            count_selection=count_selection,
            genes=genes_sets[set],
            data_regression_genes=data_regression_genes,
        )
    # Return information.
    return sets


def collect_union_sets_genes(
    sets=None,
    union_variables=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        sets (dict<dict<list<str>>>): sets of genes that associate
            significantly to variables in regression
        union_variables (list<str>): names of variables for which to collect
            union

    raises:

    returns:
        (dict<dict<list<str>>>): sets of genes that associate significantly to
            variables in regression

    """

    sets = copy.deepcopy(sets)
    for distribution in sets.keys():
        union = list()
        for variable in union_variables:
            union.extend(sets[distribution][variable])
        sets[distribution]["union"] = utility.collect_unique_elements(
            elements_original=union
        )
    # Return information.
    return sets


def report_association_variables_sets_genes(
    variables=None,
    sets=None,
    genes_sets=None,
    distribution_master=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variables (list<str>): names of independent regression variables
        sets (dict<dict<list<str>>>): sets of genes that associate
            significantly to variables in regression
        genes_sets (dict<list<str>>): identifiers of all genes in groups by
            distributions of their pan-tissue signals
        distribution_master (str): distribution group of genes for which to
            report allocation to variables

    raises:

    returns:
        (dict<dict<list<str>>>): sets of genes that associate significantly to
            variables in regression

    """

    sets = copy.deepcopy(sets)
    utility.print_terminal_partition(level=2)
    print("Counts of genes.")
    for distribution in genes_sets.keys():
        print(distribution + " genes: " + str(len(genes_sets[distribution])))
    utility.print_terminal_partition(level=3)
    print("Counts of genes that associate significantly with each variable.")
    print("master distribution genes: " + distribution_master)
    for variable in variables:
        utility.print_terminal_partition(level=4)
        print(variable)
        genes_master = genes_sets[distribution_master]
        count_master = len(genes_master)
        count_variable = len(sets[distribution_master][variable])
        percentage = (count_variable / count_master) * 100
        print("count variable genes: " + str(count_variable))
        print("master genes: " + str(round(percentage, 3)) + "%")
        pass
    pass


##########
# Tissue-tissue correlation


def calculate_pairwise_correlations(
    data=None,
    key=None,
):
    """
    Calculates the pairwise correlation coefficients for multiple pairs of
    features.

    Data has observations across rows and features across columns.

    index key           feature_1 feature_2 feature_3 feature_4 feature_5
    0     observation_1 ...       ...       ...       ...       ...
    1     observation_2 ...       ...       ...       ...       ...
    2     observation_3 ...       ...       ...       ...       ...
    3     observation_4 ...       ...       ...       ...       ...
    4     observation_5 ...       ...       ...       ...       ...

    arguments:
        key (str): name of column to match observations
        data (object): Pandas data frame of features across observations

    raises:

    returns:
        (dict<dict<dict<float>>>): matrix of counts of common observations and
            correlation coefficients from pairwise comparison of features

    """

    # Organize data.
    data_copy = data.copy(deep=True)
    data_copy.set_index(
        keys=[key],
        drop=True,
        append=False,
        inplace=True,
    )
    features = data_copy.columns.values.tolist()
    # Determine pairs of features available for each observation.
    # Use sequence-specific permutations to simplify assembly of matrix.
    pairs = list(itertools.permutations(features, 2))

    # Collect counts and correlations across pairs of features.
    correlations = dict()
    # Iterate on pairs of features.
    for pair in pairs:
        # Select data for pair of features.
        data_pair = data_copy.loc[:, list(pair)]
        # Remove observations with missing values for either feature.
        data_pair.dropna(
            axis="index",
            how="any",
            inplace=True,
        )
        # Determine observations for which pair of features has matching
        # values.
        count = data_pair.shape[0]
        # Calculate correlation.
        if count > 1:
            correlation = scipy.stats.pearsonr(
                data_pair[pair[0]].values,
                data_pair[pair[1]].values,
            )[0]
            pass
        else:
            correlation = float("nan")
            pass

        # Collect information.
        if not pair[0] in correlations:
            correlations[pair[0]] = dict()
            pass
        correlations[pair[0]][pair[1]] = dict()
        correlations[pair[0]][pair[1]]["count"] = count
        correlations[pair[0]][pair[1]]["correlation"] = correlation

    # Return information.
    return correlations


def calculate_mean_correlation(
    observation,
    correlations=None,
):
    """
    Calculates the mean of pairwise correlation coefficients for multiple pairs
    of features from a single observation.

    arguments:
        observation (object): Pandas series of features from a single
            observation
        correlations (dict<dict<dict<float>>>): matrix of counts of common
            observations and correlation coefficients from pairwise comparison
            of features

    raises:

    returns:
        (float): mean of correlation coefficients

    """

    # Determine pairs of features for which observation has valid values.
    # Remove features with missing values for the observation.
    observation.dropna(
        inplace=True,
    )
    # Extract names of valid features.
    features = observation.index.values.tolist()
    # Determine pairs of features.
    pairs = list(itertools.combinations(features, 2))
    #print(features)
    #print(len(features))
    #print(pairs)
    # Collect counts and correlations across pairs of features.
    records = list()
    # Iterate on pairs of features.
    for pair in pairs:
        count = correlations[pair[0]][pair[1]]["count"]
        correlation = correlations[pair[0]][pair[1]]["correlation"]
        # Collection information.
        record = dict()
        record["count"] = count
        record["correlation"] = correlation
        records.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(
        records=records
    )
    # Apply Fisher's z transformation to correlation coefficients.
    # Fisher's z transformation is equivalent to the inverse hyperbolic
    # tangent.
    if False:
        data["transformation"] = data["correlation"].map(
            lambda x: math.atanh(x),
        )
    data["transformation"] = data["correlation"].apply(
        lambda x: numpy.arctanh(x),
    )

    # Calculate mean across transformations of correlations.
    size = data.shape[0]
    data["numerator"] = data.aggregate(
        lambda row: ((row["count"] - 3) * row["transformation"]),
        axis="columns",
    )
    data["denominator"] = data.aggregate(
        #lambda row: (row["count"] - (3 * size)), # incorrect
        #lambda row: ((row["count"] - 3) * size), # incorrect
        lambda row: ((row["count"] - 3)),
        axis="columns",
    )
    numerator = data["numerator"].aggregate("sum")
    denominator = data["denominator"].aggregate("sum")
    mean_transformation_weight = (numerator / denominator)
    mean_transformation_simple = data["transformation"].aggregate("mean")

    # Apply inverse Fisher's z transformation to mean.
    # Inverse Fisher's z transformation is equivalent to the hyperbolic
    # tangent.
    mean_weight = math.tanh(mean_transformation_weight)
    mean_simple = math.tanh(mean_transformation_simple)

    return mean_weight


def calculate_mean_pairwise_correlations(
    correlations=None,
    data=None,
    key=None,
):
    """
    Calculates the mean pairwise correlation coefficients for multiple pairs of
    features across observations.

    Data has observations across rows and features across columns.

    index key           feature_1 feature_2 feature_3 feature_4 feature_5
    0     observation_1 ...       ...       ...       ...       ...
    1     observation_2 ...       ...       ...       ...       ...
    2     observation_3 ...       ...       ...       ...       ...
    3     observation_4 ...       ...       ...       ...       ...
    4     observation_5 ...       ...       ...       ...       ...

    arguments:
        correlations (dict<dict<dict<float>>>): matrix of counts of common
            observations and correlation coefficients from pairwise comparison
            of features
        key (str): name of column to match observations
        data (object): Pandas data frame of features across observations

    raises:

    returns:
        (object): Pandas data frame of mean correlations across observations

    """

    # Organize data.
    data_copy = data.copy(deep=True)
    data_copy.set_index(
        keys=[key],
        drop=True,
        append=False,
        inplace=True,
    )
    features = data_copy.columns.values.tolist()

    # Iterate on observations.
    data_copy["correlation"] = data_copy.apply(
        calculate_mean_correlation,
        axis="columns",
        raw=False,
        correlations=correlations,
    )

    # Organize data.
    data_copy.drop(
        labels=features,
        axis="columns",
        inplace=True
    )

    # Return data.
    return data_copy


def calculate_mean_tissue_pairwise_correlations(
    data_gene_persons_tissues_signals=None
):
    """
    Calculates the mean pairwise correlation coefficients between available
    tissues across persons.

    Source data has persons across rows and tissues across columns.

    index    tissue_1 tissue_2 tissue_3 tissue_4 tissue_5
    person_1 ...      ...      ...      ...      ...
    person_2 ...      ...      ...      ...      ...
    person_3 ...      ...      ...      ...      ...
    person_4 ...      ...      ...      ...      ...
    person_5 ...      ...      ...      ...      ...

    Product data.

    index    correlation_mean
    person_1 ...
    person_2 ...
    person_3 ...
    person_4 ...
    person_5 ...

    Method
    1. determine pairs of available tissues
    2. determine counts of usable persons for each pairwise tissue comparison
    3. calculate correlation coefficients between pairs of tissues across
        persons
    4. organize counts and correlation coefficients within a tissue-tissue
        adjacency matrix
    5. calculate mean correlation coefficient for each person
    5.1. determine pairs of tissues available for each person
    5.2. collect counts and correlation coefficients for each pair of tissues
    5.3. apply Fisher's z transformation to correlation coefficients
    5.4. calculate mean Fisher's z transformation with consideration of counts
    5.5. apply inverse Fisher's z transformation

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of mean tissue-tissue correlation
            coefficients across persons

    """

    # Organize data.
    tissues = data_gene_persons_tissues_signals.columns.values.tolist()
    data_index = data_gene_persons_tissues_signals.reset_index(
        level="person",
        drop=False,
        inplace=False
    )
    # Calculate pairwise correlations between pairs of tissues.
    correlations = calculate_pairwise_correlations(
        key="person",
        data=data_index,
    )

    #utility.print_terminal_partition(level=2)
    #print(correlations["adipose"]["artery"])
    #print(correlations["artery"]["adipose"])
    #print(correlations["lung"]["heart"])
    #print(correlations["heart"]["lung"])
    #utility.print_terminal_partition(level=2)

    # Calculate mean tissue-tissue correlations across persons.
    data_gene_mean_correlations = calculate_mean_pairwise_correlations(
        correlations=correlations,
        key="person",
        data=data_index,
    )

    # Return information.
    return data_gene_mean_correlations


########## Keep these functions#
# Analysis of Variance (ANOVA)


def select_values_by_binary_category(
    category=None,
    value=None,
    data=None
):
    """
    Selects a gene's aggregate scores by binary category.

    The dependent variable of interest must have column name "value".

    arguments:
        category (str): name of a category
        value (int): value of a category
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (list<float>): values for a group

    """

    # Select values of gene's scores by attribute category.
    data_selection = data.loc[
        data[category] == value, :
    ]
    return data_selection["value"].values.tolist()


def evaluate_variance_by_binary_groups(
    groups=None,
    threshold=None,
    data=None
):
    """
    Evaluates the variance between values in binary groups by one-way analysis
    of variance.

    The probability is sometimes a missing value.
    This missing value does not seem to indicate missing values in the groups.
    Potential problems
        1. a specific group lacks values
        2. multiple groups have the same values

    arguments:
        groups (list<str>): names of binary groups
        threshold (int): minimal count of values in a group for inclusion
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (dict<float>): values of results from the analysis of variance

    """

    # Collect values in each group.
    groups_values = list()
    for group in groups:
        values = select_values_by_binary_category(
            category=group,
            value=1,
            data=data,
        )
        if len(values) > threshold:
            groups_values.append(values)
    # Analyze variance.
    # Use splat operator to extract arguments from list.
    # Alternative: scipy.stats.f_oneway(*[df[col] for col in df])
    # The probability is sometimes a missing value.
    # Even when this probability is missing, there do not seem to be
    # irregularities in the values.
    if (len(groups_values) > 1):
        report = scipy.stats.f_oneway(*groups_values)
        # Compile information.
        information = {
            "statistic": report[0],
            "probability": report[1],
        }
    else:
        # Compile information.
        information = {
            "statistic": float("nan"),
            "probability": float("nan"),
        }
    # Return information.
    return information


# Write.


def write_product_regression(
    cohort=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_residuals_genes = os.path.join(
        paths[cohort]["regression"], "residuals_genes.pickle"
    )
    path_data_regression_genes = os.path.join(
        paths[cohort]["regression"], "data_regression_genes.pickle"
    )
    path_data_regression_genes_text = os.path.join(
        paths[cohort]["regression"], "data_regression_genes.tsv"
    )

    # Write information to file.
    with open(path_residuals_genes, "wb") as file_product:
        pickle.dump(information["residuals_genes"], file_product)
    pandas.to_pickle(
        information["data_regression_genes"],
        path_data_regression_genes
    )
    information["data_regression_genes"].to_csv(
        path_or_buf=path_data_regression_genes_text,
        sep="\t",
        header=True,
        index=True,
    )

    pass


def write_product_genes_set_variable(
    cohort=None,
    distribution=None,
    variable=None,
    sets_genes=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        distribution (str): group of genes by distribution of pan-tissue
            signals
        variable (str): name of independent regression variable
        sets_genes (dict<list<str>>): sets of genes that associate
            significantly to variables in regression
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_variable = os.path.join(
        paths[cohort]["distribution"][distribution],
        str(variable + ".txt")
    )
    # Write information to file.
    utility.write_file_text_list(
        elements=sets_genes[distribution][variable],
        delimiter="\n",
        path_file=path_variable
    )
    pass


def write_product_genes(
    cohort=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        information (object): information to write to file.
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_sets_genes = os.path.join(
        paths[cohort]["genes"], "sets_genes.pickle"
    )
    # Write information to file.
    with open(path_sets_genes, "wb") as file_product:
        pickle.dump(information["sets_genes"], file_product)

    # Write individual gene sets.
    distributions = list()
    distributions.append("any")
    distributions.append("multimodal")
    distributions.append("nonmultimodal")
    distributions.append("unimodal")
    for distribution in distributions:
        for variable in information["sets_genes"][distribution].keys():
            write_product_genes_set_variable(
                cohort=cohort,
                distribution=distribution,
                variable=variable,
                sets_genes=information["sets_genes"],
                paths=paths,
            )
    pass


###############################################################################
# Procedure


def organize_data_regress_cases_report_write(
    cohort=None,
    paths=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Report.
    utility.print_terminal_partition(level=1)
    print("... Prediction procedure for: " + str(cohort) + " persons...")
    utility.print_terminal_partition(level=2)

    # Read source information from file.
    source = read_source(
        cohort=cohort,
        dock=paths["dock"],
    )
    # Specify genes sets.
    genes_sets = source["genes_candidacy"]
    #genes_sets = source["genes_heritability"]
    # Define variables for regression.
    variables = selection.define_variables()

    # Regression on hypothesis variables.
    if False:
        #genes_regression = random.sample(genes_sets["any"], 100)
        genes_regression = genes_sets["any"]#[0:100]
        bin_regression = organize_data_regress_cases_report(
            genes=genes_regression,
            variables_regression=variables[cohort]["model_hypothesis"],
            data_persons_properties=source["data_persons_properties"],
            data_signals_genes_persons=source["data_signals_genes_persons"],
            data_gene_annotation=source["data_gene_annotation"],
            threshold_discovery=0.05,
            discovery=True,
            report=False,
        )
        write_product_regression(
            cohort=cohort,
            information=bin_regression,
            paths=paths,
        )

    # Select and organize sets of genes by their associations and overlap
    # with modality sets.
    # Select genes with significant association with each hypothetical
    # variable of interest.
    if True:
        variables_interest = variables[cohort]["model_hypothesis"]
        variables_query = variables["query"]["six"] # one, two,
        #model = "technique"
        #variables_collection = variables["batch"]
        source_regression = read_source_regression(
            cohort=cohort,
            paths=paths,
        )
        # Select genes by their association to variables of hypothetical
        # interest.
        sets_genes = select_genes_by_gene_distributions_variable_associations(
            variables_separate=variables_interest,
            variables_together=variables_query,
            genes_sets=genes_sets,
            threshold_r_square=0.25, # 0.1, 0.25
            threshold_discovery=0.05,
            coefficient_sign="any", # "negative", "positive", or "any"
            count_selection=1000, # select count of genes with greatest absolute values of regression parameters (beta coefficients)
            data_regression_genes=source_regression["data_regression_genes"],
        )
        # Include union sets.
        sets_genes = collect_union_sets_genes(
            sets=sets_genes,
            union_variables=variables_interest,
        )
        # Include query set in variables of interest.
        variables_summary = copy.deepcopy(variables_interest)
        variables_summary.append("query")
        variables_summary.append("union")
        report_association_variables_sets_genes(
            variables=variables_summary,
            sets=sets_genes,
            genes_sets=genes_sets,
            distribution_master="multimodal", # "any" or "multimodal"
        )
        # Compile information.
        bin_sets = dict()
        bin_sets["sets_genes"] = sets_genes
        write_product_genes(
            cohort=cohort,
            information=bin_sets,
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
    # Organize data and regress cases.
    if True:
        organize_data_regress_cases_report_write(
            cohort="selection",
            paths=paths,
        )
    if False:
        organize_data_regress_cases_report_write(
            cohort="respiration",
            paths=paths,
        )
    if False:
        organize_data_regress_cases_report_write(
            cohort="ventilation",
            paths=paths,
        )
    pass


if (__name__ == "__main__"):
    execute_procedure()
