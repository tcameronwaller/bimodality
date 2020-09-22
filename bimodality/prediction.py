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


def initialize_directories(
    cohorts_models=None,
    restore=None,
    dock=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        cohorts_models (dict<list<str>>): models in each cohort
        restore (bool): whether to remove previous versions of data
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
    paths["context_summaries"] = os.path.join(
        paths["prediction"], "context_summaries"
    )
    # Initialize directories.
    utility.create_directories(
        path=paths["context_summaries"]
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["prediction"])
        utility.create_directory(path=paths["prediction"])
    # Define paths for cohorts of persons.
    for cohort in cohorts_models.keys():
        paths[cohort] = dict()
        # Define paths for regression models.
        for model in cohorts_models[cohort]:
            paths[cohort][model] = dict()
            # Define path for comprehensive information from regression on
            # model.
            paths[cohort][model]["regression"] = os.path.join(
                paths["prediction"], cohort, model, "regression"
            )
            # Define paths for sets of genes.
            paths[cohort][model]["regressions_discoveries"] = os.path.join(
                paths["prediction"], cohort, model, "regressions_discoveries"
            )
            paths[cohort][model]["genes_associations"] = os.path.join(
                paths["prediction"], cohort, model, "genes_associations"
            )
            # Initialize directories.
            utility.create_directories(
                path=paths[cohort][model]["regression"]
            )
            utility.create_directories(
                path=paths[cohort][model]["regressions_discoveries"]
            )
            utility.create_directories(
                path=paths[cohort][model]["genes_associations"]
            )
            pass
        pass
    # Return information.
    return paths


def read_source_regression(
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
    if False:
        genes_heritability = integration.read_source_genes_sets_heritability(
            cohort="selection",
            dock=dock,
        )
    # Return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "genes_candidacy": genes_candidacy,
        #"genes_heritability": genes_heritability,
        "data_signals_genes_persons": data_signals_genes_persons,
    }


def read_source_association(
    cohort=None,
    model=None,
    cohort_multimodality=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): name of regression model
        cohort_multimodality (str): cohort from which to draw multimodal genes
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_regression_genes = os.path.join(
        dock, "prediction", cohort, model, "regression",
        "data_regression_genes.pickle"
    )
    # Read information from file.
    data_regression_genes = pandas.read_pickle(
        path_data_regression_genes
    )
    # Read genes sets.
    # Specify the same set of multimodal genes for all cohorts.
    sets_genes = integration.read_source_genes_sets_collection_candidacy(
        cohort=cohort_multimodality,
        dock=dock,
    )
    if False:
        genes_heritability = integration.read_source_genes_sets_heritability(
            cohort="selection",
            dock=dock,
        )
    # Compile and return information.
    return {
        "data_regression_genes": data_regression_genes,
        "sets_genes": sets_genes,
        #"genes_heritability": genes_heritability,
    }


def read_source_association_summary(
    cohorts_models=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohorts_models (dict<list<str>>): models in each cohort
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Read regression summaries for each cohort and model.
    cohorts_models_regressions = dict()
    for cohort in cohorts_models.keys():
        cohorts_models_regressions[cohort] = dict()
        for model in cohorts_models[cohort]:
            # Specify directories and files.
            path_regressions_discoveries = os.path.join(
                dock, "prediction", cohort, model, "regressions_discoveries",
                "regressions_discoveries.pickle"
            )
            # Read information from file.
            with open(path_regressions_discoveries, "rb") as file_source:
                regressions = pickle.load(file_source)
            cohorts_models_regressions[cohort][model] = regressions
            pass
        pass
    # Read sets of genes for each cohort and model.
    cohorts_models_sets_genes = dict()
    for cohort in cohorts_models.keys():
        cohorts_models_sets_genes[cohort] = dict()
        for model in cohorts_models[cohort]:
            # Specify directories and files.
            path_sets_genes = os.path.join(
                dock, "prediction", cohort, model, "genes_associations",
                "sets_genes.pickle"
            )
            # Read information from file.
            with open(path_sets_genes, "rb") as file_source:
                sets_genes = pickle.load(file_source)
            cohorts_models_sets_genes[cohort][model] = sets_genes
            pass
        pass
    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    # Read information from file.
    path_data_regression_models = os.path.join(
        dock, "annotation", "regression",
        "table_regression_models.tsv"
    )
    data_regression_models = pandas.read_csv(
        path_data_regression_models,
        sep="\t",
        header=0,
    )
    # Select relevant information.
    data_regression_models = data_regression_models.loc[
        (data_regression_models["type"].isin(["hypothesis", "genotype"])), :
    ]
    variables = data_regression_models["variable"].to_list()
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "cohorts_models_regressions": cohorts_models_regressions,
        "cohorts_models_sets_genes": cohorts_models_sets_genes,
        "variables": variables
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


def standardize_genes_signals(
    data_signals_genes_persons=None,
    report=None,
):
    """
    Calculates the standard (z-score) of genes' pan-tissue signals across
    persons.

    arguments:
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of pan-tissue signals across genes and
            persons

    """

    # Calculate standard score.
    # This method inserts missing values if the standard deviation is zero.
    data_standard = data_signals_genes_persons.apply(
        lambda series: scipy.stats.zscore(
            series.to_numpy(),
            axis=0,
            ddof=1, # sample standard deviation
            nan_policy="omit", # ignore missing values
        ),
        axis="index", # Apply function to each column of data.
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary statistics for signals after standardization.")
        print(data_standard.iloc[0:10, 0:10])
        data_mean = data_standard.aggregate(
            lambda x: x.mean(),
            axis="index"
        )
        print("Mean")
        print(data_mean.iloc[0:10])
        utility.print_terminal_partition(level=3)
        data_deviation = data_standard.aggregate(
            lambda x: x.std(),
            axis="index"
        )
        print("Standard deviation")
        print(data_deviation.iloc[0:10])
        utility.print_terminal_partition(level=2)
    return data_standard


def organize_dependent_independent_variables(
    variables=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
    standardization=None,
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
        standardization (bool): whether to standardize the dependent variable
            by z score
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
    if standardization:
        data_signals_standard = standardize_genes_signals(
            data_signals_genes_persons=data_signals_genes_persons,
            report=report,
        )
    else:
        data_signals_standard = data_signals_genes_persons
    # Join signals on persons' properties.
    data_variables = data_persons_properties.join(
        data_signals_standard,
        how="left",
        on="person"
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=1)
        print("independent variables")
        utility.print_terminal_partition(level=3)
        print(data_persons_properties)
        utility.print_terminal_partition(level=2)
        print("data_variables")
        utility.print_terminal_partition(level=3)
        print(data_variables)
        utility.print_terminal_partition(level=2)
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
        #values_independence = data.loc[ :, independence].to_numpy()
        data_independence = data.loc[ :, independence]
        # Introduce constant value for intercept.
        # If any column in the independent variables already has constant
        # values, then the function skips it by default.
        # It is necessary to change parameter "has_constant" to avoid this
        # conditional behavior.
        data_independence_intercept = statsmodels.api.add_constant(
            data_independence,
            prepend=True, # insert intercept constant first
            has_constant="add", # Introduce new intercept constant regardless
        )
        # Define model.
        model = statsmodels.api.OLS(
            values_dependence,
            data_independence_intercept,
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
        # Collect parameters, probabilities, and statistics.
        report_parameters = pandas.Series(data=report.params)
        report_probabilities = pandas.Series(data=report.pvalues)
        parameters = dict()
        probabilities = dict()
        inflations = dict()
        if ("const" in report_parameters.index):
            #parameters["intercept_parameter"] = report.params[0]
            parameters["intercept_parameter"] = report_parameters["const"]
        else:
            parameters["intercept_parameter"] = float("nan")
            utility.print_terminal_partition(level=4)
            print("Warning: regression data does not have constant intercept.")
            print(dependence)
            utility.print_terminal_partition(level=4)
        if ("const" in report_probabilities.index):
            #probabilities["intercept_probability"] = report.pvalues[0]
            probabilities["intercept_probability"] = (
                report_probabilities["const"]
            )
        else:
            probabilities["intercept_probability"] = float("nan")
            utility.print_terminal_partition(level=4)
            print("Warning: regression data does not have constant intercept.")
            print(dependence)
            utility.print_terminal_partition(level=4)
        inflations["intercept_inflation"] = float("nan")
        # Iterate on each independent variable.
        # Initiate counter at 1 to assume that intercept is at index 0.
        counter = 1
        # Accommodate index for intercept.
        for variable in independence:
            # Coefficient or parameter.
            parameter = str(variable + ("_parameter"))
            #parameters[parameter] = report.params[counter]
            parameters[parameter] = report_parameters[variable]
            # Probability.
            probability = str(variable + ("_probability"))
            #probabilities[probability] = report.pvalues[counter]
            probabilities[probability] = report_probabilities[variable]
            # Variance Inflation Factor (VIF).
            inflation = str(variable + ("_inflation"))
            inflation_value = (
                statsmodels.stats.outliers_influence.variance_inflation_factor(
                    data_independence_intercept.to_numpy(),
                    counter
                )
            )
            inflations[inflation] = round(inflation_value, 2)
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
        # Monitor progress.
        counter += 1
        #print("count " + str(counter) + ": " + str(case))
        percentage = (counter / count_total) * 100
        # Report.
        if ((percentage % 10 == 0) and report):
            print("complete cases: " + str(round(percentage)) + "%")
            pass
        bin_case = regress_signal_ordinary_residuals(
            dependence=case,
            independence=variables,
            proportion=0.5,
            data=data_variables,
        )
        # Collect residuals.
        residuals_genes[case] = bin_case["residuals"]
        # Collect reports.
        record = bin_case["report"]
        record["case"] = case
        record["name"] = assembly.access_gene_name(
            identifier=case,
            data_gene_annotation=data_gene_annotation,
        )
        records.append(record)
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


def prepare_regressions_quality_summary(
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
        (object): Pandas series of mean measures of regression quality

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
        lambda column: numpy.nanmean(column.to_numpy())
    )
    # Return information
    return series


def organize_data_regress_cases_report(
    genes=None,
    variables_regression=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
    data_gene_annotation=None,
    standardization=None,
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
        standardization (bool): whether to standardize the dependent variable
            by z score
        report (bool): whether to print reports

    raises:

    returns:

    """

    if False:
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

    if False:
        utility.print_terminal_partition(level=3)
        print("genes for regression: " + str(len(genes)))
        utility.print_terminal_partition(level=3)

    # Organize dependent and independent variables for regression analysis.
    data_variables = organize_dependent_independent_variables(
        variables=variables_regression,
        data_persons_properties=data_persons_properties,
        data_signals_genes_persons=data_signals_genes_persons,
        standardization=standardization,
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
    series_quality = prepare_regressions_quality_summary(
        variables=variables_regression,
        data_regression_models=bin_regression["data_regression_genes"],
    )
    # Report means of statistics across independent models.
    if report:
        # Report information.
        utility.print_terminal_partition(level=2)
        print("Mean statistics for regression models.")
        print(series_quality)
        utility.print_terminal_partition(level=2)
        pass

    # Review.
    # 4 February 2020
    # Model statistics match in the individual regression's summary and in the
    # compilation across cases.
    # Independent variables' parameters and probabilities match in the
    # individual regression's summary and in the compilation across cases.

    # Compile information.
    bin = dict()
    bin["residuals_genes"] = bin_regression["residuals_genes"]
    bin["data_regression_genes"] = bin_regression["data_regression_genes"]
    # bin["series_quality"] =
    # Return information.
    return bin


##########
# Sets


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
    variables = copy.deepcopy(variables)
    variables.insert(0, "intercept")
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


def select_regression_cases_thresholds_discoveries(
    threshold_r_square_adjust=None,
    threshold_discovery=None,
    variables_model=None,
    sets_cases=None,
    data_regression_genes=None,
):
    """
    Select regression summaries for subsets of cases and calculate false
    discovery rates.

    This function accepts a summary of regressions (R-squares, coefficients,
    etc) across cases (genes) for a single cohort and model. This function also
    accepts multiple subsets of those cases (genes) of interest and a single
    threshold for the value of each regression's R-square adjust statistic.

    This function filters the summary of regressions by the threshold on
    R-square adjust. It then filters the summary by each subset of cases.

    For each of these subsets, this function then calculates the false
    discovery rates for each variable in the regression model.

    arguments:
        threshold_r_square_adjust (float): threshold by R-square adjust
        threshold_discovery (float): threshold by false discovery rate
        variables_model (list<str>): names of variables in regression model
        sets_cases (dict<list<str>>): identifiers of cases in sets of interest
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:
        (dict<object>): collection of Pandas data frames with statistics from
            regressions across cases

    """

    # Copy data.
    data_regressions = data_regression_genes.copy(deep=True)
    # Filter regression cases by threshold on R-square adjust.
    data_regressions_threshold = data_regressions.loc[
        data_regressions["r_square_adjust"] >= threshold_r_square_adjust, :
    ]
    # Iterate across sets of genes of interest.
    # Each set of genes is a query against the total regression cases.
    regressions_discoveries = dict()
    for set in sets_cases.keys():
        # Select regression data for genes of interest.
        data_cases = data_regressions_threshold.loc[
            data_regressions_threshold.index.isin(sets_cases[set]), :
        ]
        # Adjust parameter probabilities for multiple hypothesis tests.
        # Calculate false discovery rate by Benjamini-Hochberg's method.
        data_discoveries = calculate_regression_discoveries(
            variables=variables_model,
            threshold=threshold_discovery,
            data=data_cases,
        )
        # Organize data.
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
        columns.append("intercept_discovery")
        columns.append("intercept_significance")
        for variable in variables_model:
            columns.append(str(variable + ("_parameter")))
            columns.append(str(variable + ("_probability")))
            columns.append(str(variable + ("_discovery")))
            columns.append(str(variable + ("_significance")))
            columns.append(str(variable + ("_inflation")))
            pass
        data_discoveries = data_discoveries[[*columns]]
        # Collect information.
        regressions_discoveries[set] = data_discoveries
    # Return information.
    return regressions_discoveries


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
    # Select genes that associate significantly with each variable of interest
    # separately.
    # Any, or logic.
    sets = select_genes_by_association_variables_separate(
        variables=variables_separate,
        coefficient_sign=coefficient_sign,
        count_selection=count_selection,
        data_regression_genes=data_genes,
    )
    # Select genes that associate significantly with multiple query variables
    # together.
    # All, and logic.
    set_together = select_genes_by_association_variables_together(
        variables=variables_together,
        data_regression_genes=data_genes,
    )
    # Compile information.
    sets.update(dict(query=set_together))
    # Return information.
    return sets


def select_genes_by_gene_distributions_variable_associations(
    variables_separate=None,
    variables_together=None,
    coefficient_sign=None,
    count_selection=None,
    sets_genes=None,
    regressions_discoveries=None,
):
    """
    Selects and scales regression parameters.

    1. Subset regression data by a set of query genes of interest.
    2. Correct for multiple hypotheses within the scope of the query genes of
    interest.
    3. Identify genes that associate significantly with each variable.

    arguments:
        variables_separate (list<str>): names of independent regression
            variables for which to select genes by significant association
        variables_together (list<str>): names of multiple variables for which
            to select genes that associate significantly with all together
        coefficient_sign (str): option to select only genes with negative,
            positive, or any sign of their coefficients for the variable
        count_selection (int): count of genes with greatest and least values
            of variable's regression parameter to keep
        sets_genes (dict<list<str>>): identifiers of genes in sets of interest
        regressions_discoveries (dict<object>): collection of Pandas data
            frames with statistics from regressions across cases

    raises:

    returns:
        (dict<dict<list<str>>>): sets of genes that associate significantly to
            variables in regression

    """

    # Copy data.
    regressions_discoveries = copy.deepcopy(regressions_discoveries)
    # Select genes from each cohort that associate significantly with each
    # variable.
    sets = dict()
    for set in sets_genes.keys():
        sets[set] = select_genes_by_association_regression_variables(
            variables_separate=variables_separate,
            variables_together=variables_together,
            coefficient_sign=coefficient_sign,
            count_selection=count_selection,
            genes=sets_genes[set],
            data_regression_genes=regressions_discoveries[set],
        )
    # Return information.
    return sets


def collect_union_priority_sets_genes(
    sets=None,
    union_variables=None,
    priority_variables=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        sets (dict<dict<list<str>>>): sets of genes that associate
            significantly to variables in regression
        union_variables (list<str>): names of variables for which to collect
            union
        priority_variables (list<str>): names of variables for which to collect
            priority genes

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
        priority = list()
        for variable in priority_variables:
            if variable in sets[distribution]:
                priority.extend(sets[distribution][variable])
        sets[distribution]["priority"] = utility.collect_unique_elements(
            elements_original=priority
        )
    # Return information.
    return sets


def report_association_variables_sets_genes(
    variables=None,
    sets_genes_association=None,
    sets_genes_original=None,
    set_master=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variables (list<str>): names of independent regression variables
        sets_genes_association (dict<list<str>>): identifiers of genes of
            interest that associate with variables of interest
        sets_genes_original (dict<list<str>>): identifiers of genes in sets of
            interest
        set_master (str): set of genes for which to report allocation to
            variables

    raises:

    returns:
        (dict<dict<list<str>>>): sets of genes that associate significantly to
            variables in regression

    """

    sets_genes_association = copy.deepcopy(sets_genes_association)
    utility.print_terminal_partition(level=2)
    print("Counts of genes.")
    for set in sets_genes_original.keys():
        print(set + " genes: " + str(len(sets_genes_original[set])))
    utility.print_terminal_partition(level=3)
    print("Counts of genes that associate significantly with each variable.")
    print("master set genes: " + set_master)
    for variable in variables:
        utility.print_terminal_partition(level=4)
        print(variable)
        genes_master = sets_genes_original[set_master]
        count_master = len(genes_master)
        count_variable = len(sets_genes_association[set_master][variable])
        percentage = (count_variable / count_master) * 100
        print("count variable genes: " + str(count_variable))
        print("master genes: " + str(round(percentage, 3)) + "%")
        pass
    pass


##########
# Summary


def collect_unique_union_genes_cohorts_models_variables(
    set_query=None,
    variables=None,
    cohorts_models=None,
    cohorts_models_sets_genes=None,
    report=None,
):
    """
    Collects and summarizes genes' associations across multiple cohorts,
    regression models, and regression independent variables.

    arguments:
        set_query (str): name of an original query of set of genes by which to
            collect significant associations to regression variables
        variables (list<str>): names of independent regression variables for
            which to summarize genes' associations
        cohorts_models (dict<list<str>>): models in each cohort
        cohorts_models_sets_genes (dict<dict<dict<list<str>>>): collection of
            sets of genes by cohort, model, and variables of interest
        report (bool): whether to print reports
    raises:

    returns:
        (object): Pandas data frame of genes' associations across multiple
            cohorts, regression models, and regression independent variables

    """

    # Collect genes.
    genes = list()
    # Iterate on cohorts, models, and variables.
    # Cohorts.
    for cohort in cohorts_models_sets_genes.keys():
        if cohort in cohorts_models.keys():
            #utility.print_terminal_partition(level=2)
            #print("cohort: " + cohort)
            # Models.
            for model in cohorts_models_sets_genes[cohort].keys():
                if model in cohorts_models[cohort]:
                    #utility.print_terminal_partition(level=3)
                    #print("model: " + model)
                    # Define sets by variables.
                    sets_variables = (
                        cohorts_models_sets_genes[cohort][model][set_query]
                    )
                    # Variables.
                    for variable in sets_variables.keys():
                        if variable in variables:
                            #utility.print_terminal_partition(level=4)
                            #print("variable: " + variable)
                            # Collect genes for cohort, model, and variable.
                            genes.extend(sets_variables[variable])
                        pass
                    pass
                pass
            pass
        pass
    # Collect unique genes.
    genes_unique = utility.collect_unique_elements(
        elements_original=genes
    )
    # Return information.
    return genes_unique


def assign_binary_associations_across_genes(
    genes_true=None,
    genes_total=None,
    record_original=None,
):
    """
    Collects and summarizes genes' associations across multiple cohorts,
    regression models, and regression independent variables.

    arguments:
        genes_true (list<str>): identifiers of genes
        genes_total (list<str>): identifiers of genes
        record_original (dict): information for row in summary table
    raises:

    returns:
        (dict): novel record

    """

    # Copy record.
    record_novel = copy.deepcopy(record_original)
    # Assign binary values for each gene.
    for gene in genes_total:
        if gene in genes_true:
            record_novel[gene] = 1
        else:
            record_novel[gene] = float("nan")
            pass
        pass
    # Return information.
    return record_novel


def collect_regressions_discoveries_across_genes(
    genes_union=None,
    variable=None,
    data_regression_genes=None,
    record_original=None,
):
    """
    Collects and summarizes genes' associations across multiple cohorts,
    regression models, and regression independent variables.

    arguments:
        genes_union (list<str>): identifiers of genes that associate with any
            of the variables of interest
        variable (str): name of variable in regression model
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases
        record_original (dict): information for row in summary table
    raises:

    returns:
        (dict): novel record

    """

    # Copy record.
    record_novel = copy.deepcopy(record_original)
    # Organize information.
    parameter_name = str(variable + ("_parameter"))
    discovery_name = str(variable + ("_discovery"))
    significance_name = str(variable + ("_significance"))
    # Assign binary values for each gene.
    for gene in genes_union:
        # Determine whether gene is in table.
        if gene in data_regression_genes.index.to_list():
            # Access regression information for gene.
            parameter = data_regression_genes.at[gene, parameter_name]
            discovery = data_regression_genes.at[gene, discovery_name]
            significance = data_regression_genes.at[gene, significance_name]
            # Determine whether gene associates significantly with the
            # variable.
            if significance:
                # Organize information.
                parameter_string =str(round(parameter, 3))
                discovery_string =str(numpy.format_float_scientific(
                    discovery, precision=3
                ))
                entry = str(
                    str(parameter_string) +
                    " (fdr: " + str(discovery_string) + ")"
                )
                record_novel[gene] = entry
            else:
                record_novel[gene] = ""
                pass
        else:
            record_novel[gene] = ""
            pass
        pass
    # Return information.
    return record_novel


def collect_genes_cohorts_models_variables(
    genes_union=None,
    set_query=None,
    variables=None,
    cohorts_models=None,
    cohorts_models_regressions=None,
    cohorts_models_sets_genes=None,
    report=None,
):
    """
    Collects and summarizes genes' associations across multiple cohorts,
    regression models, and regression independent variables.

    arguments:
        genes_union (list<str>): identifiers of genes that associate with any
            of the variables of interest
        set_query (str): name of an original query of set of genes by which to
            collect significant associations to regression variables
        variables (list<str>): names of independent regression variables for
            which to summarize genes' associations
        cohorts_models (dict<list<str>>): models in each cohort
        cohorts_models_regressions (dict<object>): collection of statistics
            from regressions for each cohort and model
        cohorts_models_sets_genes (dict<dict<dict<list<str>>>): collection of
            sets of genes by cohort, model, and variables of interest
        report (bool): whether to print reports
    raises:

    returns:
        (list<dict>): records of genes' associations across multiple
            cohorts, regression models, and regression independent variables

    """

    # Collect information about genes across cohorts, models, and variables.
    records = list()
    # Iterate on cohorts, models, and variables.
    # Cohorts.
    for cohort in cohorts_models_sets_genes.keys():
        if cohort in cohorts_models.keys():
            #utility.print_terminal_partition(level=2)
            #print("cohort: " + cohort)
            # Models.
            for model in cohorts_models_sets_genes[cohort].keys():
                if model in cohorts_models[cohort]:
                    #utility.print_terminal_partition(level=3)
                    #print("model: " + model)
                    # Define sets by variables.
                    sets_variables = (
                        cohorts_models_sets_genes[cohort][model][set_query]
                    )
                    data_regression_genes = (
                        cohorts_models_regressions[cohort][model][set_query]
                    )
                    # Variables.
                    for variable in sets_variables.keys():
                        if variable in variables:
                            #utility.print_terminal_partition(level=4)
                            #print("variable: " + variable)
                            # Create record for cohort, model, and variable.
                            record = dict()
                            record["cohort"] = cohort
                            record["model"] = model
                            record["variable"] = variable
                            record = (
                                collect_regressions_discoveries_across_genes(
                                    genes_union=genes_union,
                                    variable=variable,
                                    data_regression_genes=data_regression_genes,
                                    record_original=record,
                            ))
                            # Collect record.
                            records.append(record)
                            pass
                        pass
                    pass
                pass
            pass
        pass
    # Return information.
    return records


def collect_organize_genes_cohorts_models_variables(
    genes_union=None,
    set_query=None,
    variables=None,
    cohorts_models=None,
    cohorts_models_regressions=None,
    cohorts_models_sets_genes=None,
    report=None,
):
    """
    Collects and summarizes genes' associations across multiple cohorts,
    regression models, and regression independent variables.

    arguments:
        genes_union (list<str>): identifiers of genes that associate with any
            of the variables of interest
        set_query (str): name of an original query of set of genes by which to
            collect significant associations to regression variables
        variables (list<str>): names of independent regression variables for
            which to summarize genes' associations
        cohorts_models (dict<list<str>>): models in each cohort
        cohorts_models_regressions (dict<object>): collection of statistics
            from regressions for each cohort and model
        cohorts_models_sets_genes (dict<dict<dict<list<str>>>): collection of
            sets of genes by cohort, model, and variables of interest
        report (bool): whether to print reports
    raises:

    returns:
        (object): Pandas data frame of genes' associations across multiple
            cohorts, regression models, and regression independent variables

    """

    # Collect information about genes across cohorts, models, and variables.
    records = collect_genes_cohorts_models_variables(
        genes_union=genes_union,
        set_query=set_query,
        variables=variables,
        cohorts_models=cohorts_models,
        cohorts_models_regressions=cohorts_models_regressions,
        cohorts_models_sets_genes=cohorts_models_sets_genes,
        report=report,
    )
    # Organize data.
    data = pandas.DataFrame(data=records)
    data.set_index(
        ["cohort", "model", "variable"],
        append=False,
        drop=True,
        inplace=True
    )
    data.rename_axis(
        columns="identifier",
        axis="columns",
        copy=False,
        inplace=True
    )
    data_transposition = data.transpose(copy=True)
    data_transposition.reset_index(
        level=None,
        inplace=True
    )
    data_transposition.set_index(
        ["identifier"],
        append=False,
        drop=True,
        inplace=True
    )
    # Sort genes by their counts of valid variable associations.
    data_transposition["count_valid"] = data_transposition.apply(
        lambda series:
            len(list(filter(
                lambda value: (len(str(value)) > 0),
                series.to_list()
            ))),
        axis="columns", # Apply function to each row of data.
    )
    data_transposition.sort_values(
        by=["count_valid"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    # Return information.
    return data_transposition


def organize_genes_association_summary_query_set(
    set_query=None,
    variables=None,
    union_variables=None,
    priority_variables=None,
    cohorts_models=None,
    cohorts_models_regressions=None,
    cohorts_models_sets_genes=None,
    data_gene_annotation=None,
    report=None,
):
    """
    Collects and summarizes genes' associations across multiple cohorts,
    regression models, and regression independent variables.

    arguments:
        set_query (str): name of an original query of set of genes by which to
            collect significant associations to regression variables
        variables (list<str>): names of independent regression variables for
            which to summarize genes' associations
        union_variables (list<str>): names of variables for which to collect
            union
        priority_variables (list<str>): names of variables for which to collect
            priority genes
        cohorts_models (dict<list<str>>): models in each cohort
        cohorts_models_regressions (dict<object>): collection of statistics
            from regressions for each cohort and model
        cohorts_models_sets_genes (dict<dict<dict<list<str>>>): collection of
            sets of genes by cohort, model, and variables of interest
        data_gene_annotation (object): Pandas data frame of genes' annotations
        report (bool): whether to print reports
    raises:

    returns:
        (object): Pandas data frame of genes' associations across multiple
            cohorts, regression models, and regression independent variables

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=1)
        print("Summary of associations for query gene set: " + str(set_query))
        pass
    # Collect unique identifiers of genes that association with any variables
    # across cohorts and models.
    genes_union = collect_unique_union_genes_cohorts_models_variables(
        set_query=set_query,
        variables=union_variables,
        cohorts_models=cohorts_models,
        cohorts_models_sets_genes=cohorts_models_sets_genes,
        report=report,
    )
    genes_priority = collect_unique_union_genes_cohorts_models_variables(
        set_query=set_query,
        variables=priority_variables,
        cohorts_models=cohorts_models,
        cohorts_models_sets_genes=cohorts_models_sets_genes,
        report=report,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("unique union genes: " + str(len(genes_union)))
        print("unique priority genes: " + str(len(genes_priority)))
        utility.print_terminal_partition(level=2)
    # Collect summary of genes' associations across cohorts, models, and
    # variables.
    data_summary = collect_organize_genes_cohorts_models_variables(
        genes_union=genes_union,
        set_query=set_query,
        variables=variables,
        cohorts_models=cohorts_models,
        cohorts_models_regressions=cohorts_models_regressions,
        cohorts_models_sets_genes=cohorts_models_sets_genes,
        report=report,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(data_summary)
        print("unique union genes: " + str(len(data_summary.index.to_list())))
    # Return information.
    return data_summary


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
    model=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): name of regression model
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_residuals_genes = os.path.join(
        paths[cohort][model]["regression"], "residuals_genes.pickle"
    )
    path_data_regression_genes = os.path.join(
        paths[cohort][model]["regression"], "data_regression_genes.pickle"
    )
    path_data_regression_genes_text = os.path.join(
        paths[cohort][model]["regression"], "data_regression_genes.tsv"
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
    model=None,
    set=None,
    variable=None,
    sets_genes=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): name of regression model
        set (str): name of set of genes for association to variables
        variable (str): name of independent regression variable
        sets_genes (dict<list<str>>): sets of genes that associate
            significantly to variables in regression
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Initialize directories.
    path_set = os.path.join(
        paths[cohort][model]["genes_associations"], "sets", set
    )
    utility.create_directories(
        path=path_set
    )

    # Specify directories and files.
    path_variable = os.path.join(
        path_set,
        str(variable + ".txt")
    )
    # Write information to file.
    utility.write_file_text_list(
        elements=sets_genes[set][variable],
        delimiter="\n",
        path_file=path_variable
    )
    pass


def write_product_genes_associations(
    cohort=None,
    model=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): name of regression model
        information (object): information to write to file.
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_regressions_discoveries = os.path.join(
        paths[cohort][model]["regressions_discoveries"],
        "regressions_discoveries.pickle"
    )
    path_sets_genes = os.path.join(
        paths[cohort][model]["genes_associations"], "sets_genes.pickle"
    )
    # Write information to file.
    with open(path_regressions_discoveries, "wb") as file_product:
        pickle.dump(information["regressions_discoveries"], file_product)
    with open(path_sets_genes, "wb") as file_product:
        pickle.dump(information["sets_genes"], file_product)

    # Write individual gene sets.
    for set in information["sets_genes"].keys():
        for variable in information["sets_genes"][set].keys():
            write_product_genes_set_variable(
                cohort=cohort,
                model=model,
                set=set,
                variable=variable,
                sets_genes=information["sets_genes"],
                paths=paths,
            )
    pass


def write_product_genes_associations_summary(
    set_query=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        set_query (str): name of an original query of set of genes by which to
            collect significant associations to regression variables
        information (object): information to write to file.
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Initialize directories.
    path_set_query = os.path.join(
        paths["context_summaries"], set_query
    )
    utility.create_directories(
        path=path_set_query
    )
    # Specify directories and files.
    path_data = os.path.join(
        path_set_query, "data_summary.tsv"
    )
    path_genes_identifiers = os.path.join(
        path_set_query, "genes_identifiers.pickle"
    )
    path_genes_names = os.path.join(
        path_set_query, "genes_names.txt"
    )
    # Write information to file.
    information["data_summary"].to_csv(
        path_or_buf=path_data,
        sep="\t",
        header=True,
        index=True,
    )
    with open(path_genes_identifiers, "wb") as file_product:
        pickle.dump(information["genes_identifiers"], file_product)
    utility.write_file_text_list(
        elements=information["genes_names"],
        delimiter="\n",
        path_file=path_genes_names
    )
    pass


###############################################################################
# Procedure


def organize_data_regress_cases_report_write(
    cohort=None,
    model=None,
    paths=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): name of regression model
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source_regression(
        cohort=cohort,
        dock=paths["dock"],
    )

    # Specify genes sets.
    genes_sets = source["genes_candidacy"]
    #genes_sets = source["genes_heritability"]

    # Define variables for regression.
    variables = selection.read_source_organize_regression_model_variables(
        model=model,
        dock=paths["dock"],
    )
    # Regression on hypothesis variables.
    #genes_regression = random.sample(genes_sets["any"], 1000)
    genes_regression = genes_sets["any"]#[0:10]
    utility.print_terminal_partition(level=2)
    print("regression genes: " + str(len(genes_regression)))
    bin_regression = organize_data_regress_cases_report(
        genes=genes_regression,
        variables_regression=variables["model"],
        data_persons_properties=source["data_persons_properties"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
        data_gene_annotation=source["data_gene_annotation"],
        standardization=True,
        report=True,
    )
    write_product_regression(
        cohort=cohort,
        model=model,
        information=bin_regression,
        paths=paths,
    )
    pass


def organize_regression_gene_associations_report_write(
    cohort=None,
    model=None,
    paths=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        model (str): name of regression model
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source.
    source = read_source_association(
        cohort=cohort,
        model=model,
        cohort_multimodality="selection", # respiration, ventilation
        dock=paths["dock"],
    )

    # Define variables for regression associations.
    variables = selection.read_source_organize_regression_model_variables(
        model=model,
        dock=paths["dock"],
    )
    variables_query = selection.define_regression_variable_queries()
    variables_interest = variables["model"]
    variables_query_interest = variables_query["six"]

    # Select regression summaries for genes of interest.
    # Impose thresholds by R-square adjust.
    # Calculate false discovery rates relative to each query set of genes.
    regressions_discoveries = select_regression_cases_thresholds_discoveries(
        threshold_r_square_adjust=0.25,
        threshold_discovery=0.05,
        variables_model=variables["model"],
        sets_cases=source["sets_genes"],
        data_regression_genes=source["data_regression_genes"],
    )
    # Report.
    data_example = regressions_discoveries["covid19_up_multimodal"]
    utility.print_terminal_partition(level=2)
    print("regression data for genes in 'covid19_up_multimodal'...")
    print(data_example.iloc[0:10, 0:10])
    utility.print_terminal_partition(level=2)

    # Select genes by their association to variables of hypothetical
    # interest.
    sets_genes_association = (
        select_genes_by_gene_distributions_variable_associations(
            variables_separate=variables_interest,
            variables_together=variables_query_interest,
            sets_genes=source["sets_genes"],
            coefficient_sign="any", # "negative", "positive", or "any"
            count_selection=1000, # select count of genes with greatest absolute values of regression parameters (beta coefficients)
            regressions_discoveries=regressions_discoveries,
    ))

    # Include union sets.
    if False:
        variables_priority = [
            "sex_y_scale", "age_scale", "ventilation_duration_scale",
            "ventilation_binary_scale", "sex_y*ventilation_binary_scale",
            "age*ventilation_binary_scale", "sex_y*age_scale",
        ]
    sets_genes_association = collect_union_priority_sets_genes(
        sets=sets_genes_association,
        union_variables=variables_interest,
        priority_variables=variables["hypothesis"],
    )
    # Include query set in variables of interest.
    variables_summary = copy.deepcopy(variables_interest)
    variables_summary.append("query")
    variables_summary.append("union")
    variables_summary.append("priority")

    report_association_variables_sets_genes(
        variables=variables_summary,
        sets_genes_association=sets_genes_association,
        sets_genes_original=source["sets_genes"],
        set_master="covid19_multimodal", # "any", "multimodal", "covid19_multimodal"
    )
    # Compile information.
    bin_sets = dict()
    bin_sets["regressions_discoveries"] = regressions_discoveries
    bin_sets["sets_genes"] = sets_genes_association
    write_product_genes_associations(
        cohort=cohort,
        model=model,
        information=bin_sets,
        paths=paths,
    )
    pass


def organize_summary_gene_set_associations_report_write(
    cohorts_models=None,
    paths=None,
):
    """
    Prepares a summary of genes and their associations to variables across all
    cohorts and models.

    There will be a summary for each original set of genes, such as genes with
    multimodal distributions or genes with differential expression in COVID-19.

    arguments:
        cohorts_models (dict<list<str>>): models in each cohort
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source.
    source = read_source_association_summary(
        cohorts_models=cohorts_models,
        dock=paths["dock"],
    )
    # Define an example of the original query gene sets.
    # These are the same sets of genes for selection and calculation of false
    # discovery rates on the regression statistics.
    # As the genes of interest and their false discovery rates are
    # interdependent, each original set of genes will have a separate summary
    # table.
    example_sets_genes = (
        source["cohorts_models_sets_genes"]["selection"]["selection_main"]
    )
    # Iterate on original query sets of genes.
    # Examples of original query sets are genes with multimodal distributions
    # or genes that have differential expression in COVID-19.
    for set_query in example_sets_genes.keys():
        # Organize summary for the query set.
        data_summary = organize_genes_association_summary_query_set(
            set_query=set_query,
            variables=source["variables"],
            union_variables=source["variables"],
            priority_variables=[
                "sex_y_scale", "age_scale", "ventilation_duration_scale",
                "ventilation_binary_scale", "sex_y*ventilation_binary_scale",
                "age*ventilation_binary_scale", "sex_y*age_scale",
            ],
            cohorts_models=cohorts_models,
            cohorts_models_regressions=source["cohorts_models_regressions"],
            cohorts_models_sets_genes=source["cohorts_models_sets_genes"],
            data_gene_annotation=source["data_gene_annotation"],
            report=True,
        )
        # Collect names of genes.
        genes_identifiers = data_summary.index.to_list()
        genes_names = list()
        for identifier in genes_identifiers:
            name = assembly.access_gene_name(
                identifier=identifier,
                data_gene_annotation=source["data_gene_annotation"],
            )
            genes_names.append(name)
            pass
        # Compile information.
        bin = dict()
        bin["data_summary"] = data_summary
        bin["genes_identifiers"] = genes_identifiers
        bin["genes_names"] = genes_names
        write_product_genes_associations_summary(
            set_query=set_query,
            information=bin,
            paths=paths,
        )
        pass
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

    # Define cohorts of persons.
    cohorts = [
        "selection",
        "respiration",
        "ventilation",
    ]
    # Define regression models for each cohort.
    cohorts_models = dict()
    cohorts_models["selection"] = [
        "selection_main",
        "selection_lite",
        "selection_race",
        "selection_sex_age",
        "selection_sex_ventilation",
        "selection_age_ventilation",
        "selection_race_ventilation",
        "selection_sex_leukocyte",
        "selection_age_leukocyte",
        "selection_race_leukocyte",
    ]
    cohorts_models["respiration"] = [
        "respiration_main",
        "respiration_lite",
        "respiration_race",
        "respiration_sex_age",
    ]
    cohorts_models["ventilation"] = [
        "ventilation_main",
        "ventilation_lite",
        "ventilation_race",
        "ventilation_sex_age",
        "ventilation_sex_leukocyte",
        "ventilation_age_leukocyte",
        "ventilation_race_leukocyte",
    ]
    # Initialize directories.
    paths = initialize_directories(
        cohorts_models=cohorts_models,
        restore=False,
        dock=dock,
    )
    # Execute procedure for each cohort.
    for cohort in cohorts:
        # Execute procedure for each regression model.
        for model in cohorts_models[cohort]:
            # Organize data, regress across genes, and write information to file.
            if False:
                utility.print_terminal_partition(level=2)
                print("cohort: " + str(cohort))
                print("model: " + str(model))
                utility.print_terminal_partition(level=2)
                organize_data_regress_cases_report_write(
                    cohort=cohort,
                    model=model,
                    paths=paths,
                )
                pass
            # Select and organize sets of genes by their associations and overlap
            # with modality sets.
            # Select genes with significant association with each hypothetical
            # variable of interest.
            if True:
                utility.print_terminal_partition(level=2)
                print("cohort: " + str(cohort))
                print("model: " + str(model))
                utility.print_terminal_partition(level=2)
                organize_regression_gene_associations_report_write(
                    cohort=cohort,
                    model=model,
                    paths=paths,
                )
                pass
            pass
        pass
    pass
    # Collect summaries of regression quality for each cohort and model.
    # TODO: I need to implement this...
    # TODO: use "prepare_regressions_quality_summary"
    if False:
        organize_summary_gene_set_associations_report_write(
            cohorts_models=cohorts_models,
            paths=paths,
        )
    # Collect summaries of genes' associations with variables of interest
    # across cohorts and models.
    if False:
        organize_summary_gene_set_associations_report_write(
            cohorts_models=cohorts_models,
            paths=paths,
        )
    pass


if (__name__ == "__main__"):
    execute_procedure()
