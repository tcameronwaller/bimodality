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
import gc
import functools

# Relevant

import numpy
import pandas
import gtfparse


# Custom

import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality



##########
# Organization of persons' genotypes.

# TODO: do not introduce genotype PCs to sample data yet... wait until selection procedure...
# TODO: organize information from genotype PCs...


def read_source_persons_genotypes(
    cohort=None,
    method=None,
    dock=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        method (str): method to use for genotypes, "plink" or "gcta"
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    if method == "plink":
        # Specify directories and files.
        path_genotype_component = os.path.join(
            dock, "access_private", "relation", cohort, "plink",
            "components.eigenvec"
        )
        path_genotype_variance = os.path.join(
            dock, "access_private", "relation", cohort, "plink",
            "components.eigenval"
        )
        # Read information from file.
        data_genotype_components = pandas.read_csv(
            path_genotype_component,
            sep="\s+",
            header=0,
            names=[
                "person",
                "component_1", "component_2", "component_3", "component_4",
                "component_5", "component_6", "component_7", "component_8",
                "component_9", "component_10", "component_11", "component_12",
                "component_13", "component_14", "component_15", "component_16",
                "component_17", "component_18", "component_19", "component_20",
                "component_21", "component_22", "component_23", "component_24",
                "component_25",
            ],
            low_memory=False,
        )
        genotype_variances = utility.read_file_text_list(
            delimiter="\n",
            path_file=path_genotype_variance
        )
    elif method == "gcta":
        # Specify directories and files.
        path_genotype_component = os.path.join(
            dock, "access_private", "relation", cohort, "gcta", "bed_bim_fam",
            "components.eigenvec",
        )
        path_genotype_variance = os.path.join(
            dock, "access_private", "relation", cohort, "gcta", "bed_bim_fam",
            "components.eigenval"
        )
        # Read information from file.
        data_genotype_components = pandas.read_csv(
            path_genotype_component,
            sep="\s+",
            header=None,
            names=[
                "family", "person",
                "component_1", "component_2", "component_3", "component_4",
                "component_5", "component_6", "component_7", "component_8",
                "component_9", "component_10", "component_11", "component_12",
                "component_13", "component_14", "component_15", "component_16",
                "component_17", "component_18", "component_19", "component_20",
                "component_21", "component_22", "component_23", "component_24",
                "component_25",
            ],
            low_memory=False,
        )
        genotype_variances = utility.read_file_text_list(
            delimiter="\n",
            path_file=path_genotype_variance
        )
        pass
    # Compile and return information.
    return {
        "data_genotype_components": data_genotype_components,
        "genotype_variances": genotype_variances,
    }


def organize_person_genotype_variance(
    genotype_variances=None,
):
    """
    Associates samples, tissues, and persons and collects their attributes.

    arguments:
        genotype_variances (list<float>): variances from principal components
            on genotypes

    raises:

    returns:
        (object): Pandas data frame of variances of each principal component

    """

    count = 1
    records = list()
    for variance in genotype_variances:
        if (len(variance) > 1) and (not math.isnan(float(variance))):
            record = dict()
            record["variance"] = float(variance)
            record["component"] = count
            count += 1
            records.append(record)
            pass
        pass
    data = utility.convert_records_to_dataframe(
        records=records
    )
    return data


def organize_person_genotype_source_data(
    data_genotype_components=None,
    genotype_variances=None,
):
    """
    Organizes source data for sample attributes.

    arguments:
        data_genotype_components (object): Pandas data frame of principal
            components from genotypes for persons, as calculated in GCTA or
            PLINK2
        genotype_variances (list<float>): variances from principal components
            on genotypes

    raises:

    returns:
        (dict): collection of data

    """

    # Organize data.
    bin = dict()
    bin["data_genotype_components"] = data_genotype_components.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=False
    )
    bin["data_genotype_variances"] = organize_person_genotype_variance(
        genotype_variances=genotype_variances,
    )
    # Return information.
    return bin


def write_product_persons_genotypes(
    cohort=None,
    information=None,
    dock=None,
):
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
    path_directory = os.path.join(dock, "assembly", "genotype", cohort)
    utility.create_directories(path_directory)
    path_data_genotype_components = os.path.join(
        path_directory, "data_genotype_component.pickle"
    )
    path_data_genotype_variances = os.path.join(
        path_directory, "data_genotype_variance.pickle"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_genotype_components"],
        path_data_genotype_components
    )
    pandas.to_pickle(
        information["data_genotype_variances"],
        path_data_genotype_variances
    )

    pass


def organize_persons_genotypes(
    dock=None,
):
    """
    Organize information about samples.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    # Remove previous files to avoid version or batch confusion.
    path_genotype = os.path.join(dock, "assembly", "genotype")
    utility.remove_directory(path=path_genotype)
    # Iterate on cohorts.
    cohorts = ["selection", "respiration", "ventilation"]
    for cohort in cohorts:
        # Read source information from file.
        source = read_source_persons_genotypes(
            cohort=cohort,
            method="plink",
            dock=dock,
        )
        # Organize source data.
        bin_data = organize_person_genotype_source_data(
            data_genotype_components=source["data_genotype_components"],
            genotype_variances=source["genotype_variances"],
        )
        # Compile information.
        information = dict()
        information["data_genotype_components"] = bin_data["data_genotype_components"]
        information["data_genotype_variances"] = bin_data["data_genotype_variances"]
        #Write product information to file.
        write_product_persons_genotypes(
            cohort=cohort,
            information=information,
            dock=dock,
        )
    pass


##########
# Organization of samples' attributes, specifically associations to persons
# and tissues.


def read_source_sample(dock=None):
    """
    Reads and organizes source information from file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.

    # Access
    path_access = os.path.join(dock, "access")
    path_attribute_sample = os.path.join(path_access, "attribute_sample.txt")
    path_attribute_person = os.path.join(path_access, "attribute_person.txt")

    # Private access
    path_access_private = os.path.join(dock, "access_private")
    path_attribute_sample_private = os.path.join(
        path_access_private, "sample_attribute.txt"
    )
    path_attribute_person_private = os.path.join(
        path_access_private, "person_attribute.txt"
    )
    path_person_ancestry = os.path.join(
        path_access_private, "persons_ancestry.txt"
    )

    # Customization
    path_customization = os.path.join(dock, "customization")
    path_tissues_major = os.path.join(
        path_customization, "translation_tissues_major.tsv"
    )
    path_tissues_minor = os.path.join(
        path_customization, "translation_tissues_minor.tsv"
    )

    path_data_person_smoke = os.path.join(
        dock, "annotation_2020-06-18", "smoke",
        "data_person_smoke_threshold-1-year.csv",
    )

    # Signal
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")

    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_sample_attribute = pandas.read_csv(
        path_attribute_sample,
        sep="\t",
        header=0,
    )
    data_sample_attribute_private = pandas.read_csv(
        path_attribute_sample_private,
        sep="\t",
        header=9,
        low_memory=False,
    )
    data_person_attribute = pandas.read_csv(
        path_attribute_person,
        sep="\t",
        header=0,
    )
    data_person_attribute_private = pandas.read_csv(
        path_attribute_person_private,
        sep="\t",
        header=9,
        low_memory=False,
    )
    data_person_smoke = pandas.read_csv(
        path_data_person_smoke,
        sep="\t",
        header=0,
    )
    data_person_ancestry = pandas.read_csv(
        path_person_ancestry,
        sep="\t",
        header=0,
        low_memory=False,
    )
    data_tissues_major = pandas.read_csv(
        path_tissues_major,
        sep="\t",
        header=0,
    )
    data_tissues_minor = pandas.read_csv(
        path_tissues_minor,
        sep="\t",
        header=0,
    )

    # Extract identifiers of samples with measurements for genes.
    samples_gtex = read_gene_signal_samples(dock=dock)

    # Compile and return information.
    return {
        "data_sample_attribute": data_sample_attribute,
        "data_sample_attribute_private": data_sample_attribute_private,
        "data_person_attribute": data_person_attribute,
        "data_person_attribute_private": data_person_attribute_private,
        "data_person_smoke": data_person_smoke,
        "data_person_ancestry": data_person_ancestry,
        "data_tissues_major": data_tissues_major,
        "data_tissues_minor": data_tissues_minor,
        "samples_gtex": samples_gtex,
    }


def summarize_raw_data_sample(
    data_person_attribute=None,
    data_sample_attribute=None,
):
    """
    Optimizes data types.

    arguments:
        data_person_attribute (object): Pandas data frame of attributes for all
            samples.
        data_sample_attribute (object): Pandas data frame of attributes for all
            samples.

    raises:

    returns:

    """

    # Samples beginning with code person "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of persons' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of persons' attributes.")
    print(data_person_attribute)
    print(data_person_attribute.iloc[0:10, 0:10])

    # Summarize table of samples' attributes.
    utility.print_terminal_partition(level=2)
    print("Summary of table of samples' attributes.")
    print(data_sample_attribute)
    print(data_sample_attribute.iloc[0:10, 0:7])

    pass


def organize_sample_attribute_source_data(
    data_sample_attribute=None,
    data_sample_attribute_private=None,
    data_person_attribute=None,
    data_person_attribute_private=None,
    data_person_smoke=None,
    data_person_ancestry=None,
    data_tissues_major=None,
    data_tissues_minor=None,
):
    """
    Organizes source data for sample attributes.

    arguments:
        data_sample_attribute (object): Pandas data frame of public attributes
            for all samples
        data_sample_attribute_private (object): Pandas data frame of private
            attributes for all samples
        data_person_attribute (object): Pandas data frame of public attributes
            for persons
        data_person_attribute_private (object): Pandas data frame of private
            attributes for persons
        data_person_smoke (object): Pandas data frame of curation annotation of
            person's status of smoke usage
        data_person_ancestry (object): Pandas data frame of ancestry for
            persons
        data_tissues_major (object): Pandas data frame of translations for
            names of major tissues
        data_tissues_minor (object): Pandas data frame of translations for
            names of minor tissues

    raises:

    returns:
        (dict): collection of data

    """

    # Organize data.
    bin = dict()
    bin["data_sample_attribute"] = data_sample_attribute.set_index(
        ["SAMPID"],
        append=False,
        drop=True,
        inplace=False
    )
    bin["data_sample_attribute_private"] = (
        data_sample_attribute_private.set_index(
            ["SAMPID"],
            append=False,
            drop=True,
            inplace=False
    ))
    bin["data_person_attribute"] = data_person_attribute.set_index(
        ["SUBJID"],
        append=False,
        drop=True,
        inplace=False
    )
    bin["data_person_attribute_private"] = (
        data_person_attribute_private.set_index(
            ["SUBJID"],
            append=False,
            drop=True,
            inplace=False
    ))
    bin["data_person_smoke"] = (
        data_person_smoke.set_index(
            ["SUBJID"],
            append=False,
            drop=True,
            inplace=False
    ))
    bin["data_person_ancestry"] = data_person_ancestry.set_index(
        ["IID"],
        append=False,
        drop=True,
        inplace=False
    )
    bin["data_tissues_major"] = data_tissues_major.set_index(
        ["source"],
        append=False,
        drop=True,
        inplace=False
    )
    bin["data_tissues_minor"] = data_tissues_minor.set_index(
        ["source"],
        append=False,
        drop=True,
        inplace=False
    )
    # Return information.
    return bin


def organize_person_health_variables_collections(
    variables_health=None,
    data_person_attribute=None,
    data_person_attribute_private=None,
    report=None,
):
    """
    Associates samples, tissues, and persons and collects their attributes.

    arguments:
        variables_health (dict<list<str>>): names of variables
        data_person_attribute (object): Pandas data frame of public attributes
            for persons
        data_person_attribute_private (object): Pandas data frame of private
            attributes for persons
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of data

    """

    # Copy data.
    data_private = data_person_attribute_private.copy(deep=True)
    # Select data for each collection of variables.
    bin = dict()
    for collection in variables_health.keys():
        # Select relevant variables.
        data_selection = data_private.loc[
            :, data_private.columns.isin(variables_health[collection])
        ]
        # Rename index.
        data_selection.rename_axis(
            index="person",
            axis="index",
            copy=False,
            inplace=True,
        )
        # Translate values.
        bin[collection] = data_selection.applymap(
            lambda value: translate_binary_boolean(
                value=value,
                translation="binary",
            )
        )
        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print("Collections of persons' health variables.")
            print("Collection: " + str(collection))
            utility.print_terminal_partition(level=2)
            print("before translation...")
            print(data_selection)
            utility.print_terminal_partition(level=3)
            print("after translation...")
            print(bin[collection])
    # Return information.
    return bin


def extract_gtex_sample_person_identifier(sample=None):
    """
    Extracts the person's identifier from a sample's identifier.

    arguments:
        sample (str): identifier of a sample

    raises:

    returns:
        (str): identifier of a person

    """

    split_strings = sample.split("-")
    person = "-".join(split_strings[0:2])
    return person


def translate_sex(value=None):
    """
    Translates annotations of sex.

    arguments:
        value (int): Annotation of sex

    raises:

    returns:
        (str): Name of sex

    """

    if value == 1:
        sex = "male"
    elif value == 2:
        sex = "female"
    return sex


def translate_race(value=None):
    """
    Translates annotations of race.

    arguments:
        value (int): Annotation of race

    raises:

    returns:
        (str): Name of race

    """

    if value == 1:
        race = "asia"
    elif value == 2:
        race = "africa"
    elif value == 3:
        race = "europe"
    elif value == 4:
        race = "america"
    elif (value == 98 or value == 99):
        race = "other"

    return race


def translate_ethnicity(value=None):
    """
    Translates annotations of ethnicity.

    arguments:
        value (int): annotation of ethnicity

    raises:

    returns:
        (str): name of ethnicity

    """

    if value == 1:
        ethnicity = "hispanic"
    elif (value == 0 or value == 97 or value == 98 or value == 99):
        ethnicity = "other"

    return ethnicity


def translate_hardiness(value=None):
    """
    Translates annotations of hardiness.

    arguments:
        value (int): Annotation of hardiness

    raises:

    returns:
        (str): Name of hardiness

    """

    # The original encoding represents ventilator cases as 0.
    # If persons were on a ventilator, then this treatment implies that they
    # had a long, slow death.
    # Recode the ventilator cases as 5, the most exteme value.

    if value == 0:
        return 5
    else:
        return value


def translate_season(value=None):
    """
    Translates annotations of season of death.

    arguments:
        value (int): Annotation of season

    raises:

    returns:
        (str): name of season

    """

    if isinstance(value, str):
        if str(value).strip().lower() == "fall":
            season = "fall"
        elif str(value).strip().lower() == "spring":
            season = "spring"
        elif str(value).strip().lower() == "summer":
            season = "summer"
        elif str(value).strip().lower() == "winter":
            season = "winter"
        else:
            season = float("nan")
    else:
        season = float("nan")
    return season


def translate_ancestry(value=None):
    """
    Translates annotations of ancestry.

    arguments:
        value (int): Annotation of ancestry

    raises:

    returns:
        (str): name of ancestry

    """

    if value == "Caucasian":
        ancestry = "europe"
    elif value == "African":
        ancestry = "africa"
    elif value == "Asian":
        ancestry = "asia"
    elif value == "Mexican":
        ancestry = "america"
    elif value == "Indian":
        ancestry = "india"

    return ancestry


def determine_person_ancestry(
    person=None,
    data_person_attribute_private=None,
    data_person_ancestry=None,
):
    """
    Determines a person's categorical ancestry.

    arguments:
        person (str): identifier of a person
        data_person_attribute_private (object): Pandas data frame of private
            attributes for persons
        data_person_ancestry (object): Pandas data frame of ancestry for
            persons

    raises:

    returns:
        (str): person's categorical ancestry

    """

    # Extract person's race attribute.
    race_raw = data_person_attribute_private.at[person, "RACE"]
    race_translation = translate_race(value=race_raw)

    # Extract person's ethnicity attribute.
    ethnicity_raw = data_person_attribute_private.at[person, "ETHNCTY"]
    ethnicity_translation = translate_ethnicity(value=ethnicity_raw)

    # Determine consensus between race and ethnicity.
    if (ethnicity_translation == "hispanic") and (race_translation == "other"):
        race = "america"
    else:
        race = race_translation

    # Determine ancestry from genotype ancestral analysis by Meghana Pagadala.
    if person in data_person_ancestry.index.values.tolist():
        ancestry_raw = data_person_ancestry.at[person, "ethnicity"]
        ancestry_translation = translate_ancestry(value=ancestry_raw)
    else:
        ancestry_translation = "other"
    # Determine consensus between race and ancestry.
    if race == ancestry_translation:
        ancestry = ancestry_translation
    elif race == "other" and ancestry_translation != "other":
        ancestry = ancestry_translation
    elif ancestry_translation == "other" and race != "other":
        ancestry = race
    else:
        print("person: " + person)
        print(
            "potential discrepancy... race: " + race +
            " ancestry: " + ancestry_translation
        )
        ancestry = ancestry_translation

    return ancestry


def define_person_binary_health_variables():
    """
    Defines a list of variables' names.

    arguments:

    raises:

    returns:
        (dict<list<str>>): names of variables

    """

    # Respiratory disorders.
    respiration = [
        "MHCOPD", "MHCLRD", "MHCOUGHU", "MHASTHMA", "MHSRCDSS",
        "MHPNMNIA", "MHPNMIAB", "MHTBHX",
    ]
    # Inflammatory disorders.
    inflammation = [
        "MHFVRU", "MHTEMPU", "MHARTHTS", "MHRA", "MHLAPTHU",
        "MHASCITES", "MHLUPUS", "MHSCLRDRM", "MHLVRDIS", "MHNEPH", "MHPRKNSN",
        "MHALZHMR", "MHALZDMT", "MHDMNTIA", "MHALS", "MHMS",
    ]
    # Abnormal white blood cell count.
    leukocyte = [
        "MHABNWBC",
    ]
    # Infections.
    infection = [
        "MHFNGINF", "MHBCTINF", "MHSEPSIS", "MHPSBLDCLT", "MHOPPINF",
        "MHINFLNE", "MHOSTMYLTS", "MHGNRR12M", "MHSYPH12M", "LBPRRVDRL",
        "LBRPR", "MHSTD", "MHENCEPHA", "MHFLU", "MHREYES", "MHSARS",
        "MHMENINA", "MHRBSANML", "MHSMLPXCT", "MHHIVCT", "LBHIV1NT", "LBHIVAB",
        "LBHIVO", "MHHEPCCT", "LBHCV1NT", "LBHBHCVAB", "MHHEPBCT", "LBHBCABM",
        "LBHBCABT", "LBHBSAB", "LBHBSAG", "MHWNVCT", "MHWNVHX",
    ]
    # Infectious mononucleosis.
    # Reactivation of CMV or EBV.
    mononucleosis = [
        "LBCMVTAB", "LBEBVGAB", "LBEBVMAB",
    ]
    # Steroids.
    steroid = [
        "MHHGH", "MHSTRDLT",
    ]
    # Cardiovascular disease.
    heart = [
        "MHCVD", "MHHMPHLIA", "MHHMPHLIAB", "MHHTN", "MHHRTATT", "MHHRTDIS",
        "MHHRTDISB",
    ]
    # Diabetes.
    diabetes = [
        "MHT1D", "MHT2D",
    ]
    # Compile information.
    bin = dict()
    bin["respiration"] = respiration
    bin["inflammation"] = inflammation
    bin["leukocyte"] = leukocyte
    bin["infection"] = infection
    bin["mononucleosis"] = mononucleosis
    bin["steroid"] = steroid
    bin["heart"] = heart
    bin["diabetes"] = diabetes
    # Return information.
    return bin


def translate_binary_boolean(
    value=None,
    translation=None,
):
    """
    Translates values of a binary boolean variable.

    arguments:
        value (int): value of variable from GTEx
        translation (str): translation, either binary or boolean

    raises:

    returns:
        (bool): whether value was true

    """

    # Translate values.
    if isinstance(value, str):
        if (str(value).strip().lower() == "yes"):
            if translation == "binary":
                value_translation = 1
            elif translation == "boolean":
                value_translation = True
        elif (str(value).strip().lower() == "no"):
            if translation == "binary":
                value_translation = 0
            elif translation == "boolean":
                value_translation = False
        else:
            value_translation = float("nan")
    elif (isinstance(value, int) or isinstance(value, float)):
        # Accommodate floats or inexact entries.
        if (0.5 < value and value < 1.5):
            if translation == "binary":
                value_translation = 1
            elif translation == "boolean":
                value_translation = True
        elif (0 <= value and value <= 0.5):
            if translation == "binary":
                value_translation = 0
            elif translation == "boolean":
                value_translation = False
        else:
            value_translation = float("nan")
    else:
        value_translation = float("nan")
    return value_translation


def determine_person_boolean_binary_any(
    person=None,
    variables=None,
    data_person_attribute_private=None,
):
    """
    Determines a person's categorical annotation from or logic on multiple
    binary variables.

    arguments:
        person (str): identifier of a person
        variables (list<str>): names of variables, if any of which is 1, then
        the function returns 1
        data_person_attribute_private (object): Pandas data frame of private
            attributes for persons

    raises:

    returns:
        (bool): whether person has a value of 1 for any of the variables

    """

    # Collect raw values of variables.
    values_raw = list()
    for variable in variables:
        value = data_person_attribute_private.at[person, variable]
        values_raw.append(value)
    # Translate raw values.
    values_translation = list(map(
        lambda value: translate_binary_boolean(
            value=value, translation="boolean"), values_raw
    ))
    # Remove missing values.
    # Missing values evaluate to true.
    values_valid_true = list(filter(
        lambda value: not math.isnan(value),
        values_translation
    ))
    # Invert values.
    values_valid_false = list(map(
        lambda value: not value, values_valid_true
    ))
    # Determine whether any of the values represent true.
    if len(values_valid_true) > 0:
        if any(values_valid_true):
            match = True
        elif any(values_valid_false):
            match = False
    else:
        match = float("nan")
    return match


def determine_sample_associations_attributes(
    sample=None,
    variables_health=None,
    data_sample_attribute=None,
    data_sample_attribute_private=None,
    data_person_attribute=None,
    data_person_attribute_private=None,
    data_person_smoke=None,
    data_person_ancestry=None,
    data_tissues_major=None,
    data_tissues_minor=None,
):
    """
    Associates samples, tissues, and persons and collects their attributes.

    arguments:
        sample (str): identifier of a sample
        variables_health (dict<list<str>>): names of variables
        data_sample_attribute (object): Pandas data frame of public attributes
            for all samples
        data_sample_attribute_private (object): Pandas data frame of private
            attributes for all samples
        data_person_attribute (object): Pandas data frame of public attributes
            for persons
        data_person_attribute_private (object): Pandas data frame of private
            attributes for persons
        data_person_smoke (object): Pandas data frame of curation annotation of
            person's status of smoke usage
        data_person_ancestry (object): Pandas data frame of ancestry for
            persons
        data_tissues_major (object): Pandas data frame of translations for
            names of major tissues
        data_tissues_minor (object): Pandas data frame of translations for
            names of minor tissues

    raises:

    returns:
        (dict): information about a sample's associations and attributes

    """

    # Access sample attributes.
    #removal = data_sample_attribute_private.at[sample, "SMTORMVE"]
    batch_extraction = data_sample_attribute_private.at[sample, "SMNABTCH"]
    batches_sequence = data_sample_attribute_private.at[sample, "SMGEBTCH"]
    facilities = data_sample_attribute_private.at[sample, "SMCENTER"]
    #autolysis = data_sample_attribute_private.at[sample, "SMATSSCR"]
    major = data_sample_attribute_private.at[sample, "SMTS"]
    minor = data_sample_attribute_private.at[sample, "SMTSD"]
    tissue_major = data_tissues_major.at[major, "product"]
    tissue_minor = data_tissues_minor.at[minor, "product"]

    # Access person attributes.
    person = extract_gtex_sample_person_identifier(sample=sample)
    sex_raw = data_person_attribute_private.at[person, "SEX"]
    sex_text = translate_sex(value=sex_raw)
    decade = data_person_attribute.at[person, "AGE"]
    age = data_person_attribute_private.at[person, "AGE"]
    body = data_person_attribute_private.at[person, "BMI"]
    hardiness_raw = data_person_attribute_private.at[person, "DTHHRDY"]
    hardiness = translate_hardiness(value=hardiness_raw)
    season_raw = data_person_attribute_private.at[person, "DTHSEASON"]
    season = translate_season(value=season_raw)
    delay = data_person_attribute_private.at[person, "TRDNISCH"]
    delay_start = data_person_attribute_private.at[person, "TRISCHD"]
    delay_incision = data_person_attribute_private.at[person, "TRCHSTIND"]
    refrigeration = translate_binary_boolean(
        value=data_person_attribute_private.at[person, "DTHRFG"],
        translation="boolean",
    )
    refrigeration_duration = (
        data_person_attribute_private.at[person, "DTHRFGD"]
    )
    refrigeration_unit = (
        data_person_attribute_private.at[person, "DTHRFGDU"]
    )
    ventilation = translate_binary_boolean(
        value=data_person_attribute_private.at[person, "DTHVNT"],
        translation="boolean",
    )
    ventilation_verification = translate_binary_boolean(
        value=data_person_attribute_private.at[person, "TRVNTSR"],
        translation="boolean",
    )
    ventilation_duration = (
        data_person_attribute_private.at[person, "DTHVNTD"]
    )
    ventilation_unit = (
        data_person_attribute_private.at[person, "DTHVNTDU"]
    )
    # Determine person's categorical ancestry.
    ancestry = determine_person_ancestry(
        person=person,
        data_person_attribute_private=data_person_attribute_private,
        data_person_ancestry=data_person_ancestry,
    )
    # Determine persons' history of respiratory conditions.
    respiration = determine_person_boolean_binary_any(
        person=person,
        variables=variables_health["respiration"],
        data_person_attribute_private=data_person_attribute_private,
    )
    # Determine persons' history of respiratory conditions.
    smoke = data_person_smoke.at[person, "smoke"]
    # Determine persons' history of chronic inflammation.
    inflammation = determine_person_boolean_binary_any(
        person=person,
        variables=variables_health["inflammation"],
        data_person_attribute_private=data_person_attribute_private,
    )
    leukocyte = determine_person_boolean_binary_any(
        person=person,
        variables=variables_health["leukocyte"],
        data_person_attribute_private=data_person_attribute_private,
    )
    # Determine persons' history of infection.
    infection = determine_person_boolean_binary_any(
        person=person,
        variables=variables_health["infection"],
        data_person_attribute_private=data_person_attribute_private,
    )
    mononucleosis = determine_person_boolean_binary_any(
        person=person,
        variables=variables_health["mononucleosis"],
        data_person_attribute_private=data_person_attribute_private,
    )
    # Determine persons' history of steroid use.
    steroid = determine_person_boolean_binary_any(
        person=person,
        variables=variables_health["steroid"],
        data_person_attribute_private=data_person_attribute_private,
    )
    # Determine persons' history of cardiovascular disease.
    heart = determine_person_boolean_binary_any(
        person=person,
        variables=variables_health["heart"],
        data_person_attribute_private=data_person_attribute_private,
    )
    # Determine persons' history of cardiovascular disease.
    diabetes = determine_person_boolean_binary_any(
        person=person,
        variables=variables_health["diabetes"],
        data_person_attribute_private=data_person_attribute_private,
    )
    # Compile and return information.
    information = {
        "sample": sample,
        "facilities": facilities,
        "batch_extraction": batch_extraction,
        "batches_sequence": batches_sequence,
        "tissue_major": tissue_major,
        "tissue_minor": tissue_minor,
        "person": person,
        "ancestry": ancestry,
        "sex_text": sex_text,
        "decade": decade,
        "age": age,
        "body": body,
        "hardiness": hardiness,
        "season": season,
        "delay": delay,
        "delay_start": delay_start,
        "delay_incision": delay_incision,
        "refrigeration": refrigeration,
        "refrigeration_duration": refrigeration_duration,
        "refrigeration_unit": refrigeration_unit,
        "ventilation": ventilation,
        "ventilation_verification": ventilation_verification,
        "ventilation_duration": ventilation_duration,
        "ventilation_unit": ventilation_unit,
        "respiration": respiration,
        "smoke": smoke,
        "inflammation": inflammation,
        "leukocyte": leukocyte,
        "infection": infection,
        "mononucleosis": mononucleosis,
        "steroid": steroid,
        "heart": heart,
        "diabetes": diabetes,
    }
    return information


def collect_samples_tissues_persons(
    samples=None,
    variables_health=None,
    data_sample_attribute=None,
    data_sample_attribute_private=None,
    data_person_attribute=None,
    data_person_attribute_private=None,
    data_person_smoke=None,
    data_person_ancestry=None,
    data_tissues_major=None,
    data_tissues_minor=None,
):
    """
    Associates samples, tissues, and persons and collects their attributes.

    arguments:
        samples (list<str>): identifiers of samples with signals for genes
        variables_health (dict<list<str>>): names of variables
        data_sample_attribute (object): Pandas data frame of public attributes
            for all samples
        data_sample_attribute_private (object): Pandas data frame of private
            attributes for all samples
        data_person_attribute (object): Pandas data frame of public attributes
            for persons
        data_person_attribute_private (object): Pandas data frame of private
            attributes for persons
        data_person_smoke (object): Pandas data frame of curation annotation of
            person's status of smoke usage
        data_person_ancestry (object): Pandas data frame of ancestry for
            persons
        data_tissues_major (object): Pandas data frame of translations for
            names of major tissues
        data_tissues_minor (object): Pandas data frame of translations for
            names of minor tissues

    raises:

    returns:
        (list<dict<str>>): information about persons and tissues for samples

    """

    # Collect tissues and persons for each sample.
    samples_tissues_persons = list()
    for sample in samples:
        record = determine_sample_associations_attributes(
            sample=sample,
            variables_health=variables_health,
            data_sample_attribute=data_sample_attribute,
            data_sample_attribute_private=data_sample_attribute_private,
            data_person_attribute=data_person_attribute,
            data_person_attribute_private=data_person_attribute_private,
            data_person_smoke=data_person_smoke,
            data_person_ancestry=data_person_ancestry,
            data_tissues_major=data_tissues_major,
            data_tissues_minor=data_tissues_minor,
        )
        samples_tissues_persons.append(record)
    # Return information.
    return samples_tissues_persons


def organize_samples_tissues_persons(
    dock=None,
):
    """
    Organize information about samples.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    # Remove previous files to avoid version or batch confusion.
    path_assembly_sample = os.path.join(dock, "assembly", "sample")
    utility.remove_directory(path=path_assembly_sample)

    # Read source information from file.
    source = read_source_sample(dock=dock)


    utility.print_terminal_partition(level=1)
    print("data_person_attribute_private")
    print(source["data_person_attribute_private"])

    # Summarize the structures of the raw data.
    #summarize_raw_data_sample(
    #    data_person_attribute=source["data_person_attribute"],
    #    data_sample_attribute=source["data_sample_attribute"],
    #)

    # Extract association of samples to persons and tissues.
    utility.print_terminal_partition(level=1)
    print("Association of samples to persons and tissues.")

    # Organize samples by hierarchy of persons and tissues.
    # Collect tissues and persons for each sample.
    utility.print_terminal_partition(level=2)
    print("Collection of tissues and persons for each sample.")
    print(
        "Extract sample identifiers from the original data for genes' " +
        "signals."
    )
    print("These samples are comprehensive before any filters or selection.")
    print("Extract person identifiers from sample identifiers.")
    # Extract identifiers of samples with measurements for genes.
    # Extract names of columns.
    # Pandas series
    #headers = data_gene_signal.columns
    # NumPy array
    #headers = data_gene_signal.columns.values
    # List
    #headers = source["data_gene_signal"].columns.to_list()
    # Exclude name and description.
    #samples = headers[2:]

    # Organize source data.
    bin_data = organize_sample_attribute_source_data(
        data_person_attribute=source["data_person_attribute"],
        data_person_attribute_private=source["data_person_attribute_private"],
        data_person_smoke=source["data_person_smoke"],
        data_person_ancestry=source["data_person_ancestry"],
        data_sample_attribute=source["data_sample_attribute"],
        data_sample_attribute_private=source["data_sample_attribute_private"],
        data_tissues_major=source["data_tissues_major"],
        data_tissues_minor=source["data_tissues_minor"],
    )
    # Define collections of health variables for persons.
    variables_health = define_person_binary_health_variables()
    # Organize collections of health variables across persons.
    bin_health = organize_person_health_variables_collections(
        variables_health=variables_health,
        data_person_attribute=bin_data["data_person_attribute"],
        data_person_attribute_private=(
            bin_data["data_person_attribute_private"]
        ),
        report=False,
    )
    # Collect information about samples.
    samples_tissues_persons = collect_samples_tissues_persons(
        samples=source["samples_gtex"],
        variables_health=variables_health,
        data_person_attribute=bin_data["data_person_attribute"],
        data_person_attribute_private=(
            bin_data["data_person_attribute_private"]
        ),
        data_person_smoke=(bin_data["data_person_smoke"]),
        data_person_ancestry=bin_data["data_person_ancestry"],
        data_sample_attribute=bin_data["data_sample_attribute"],
        data_sample_attribute_private=(
            bin_data["data_sample_attribute_private"]
        ),
        data_tissues_major=bin_data["data_tissues_major"],
        data_tissues_minor=bin_data["data_tissues_minor"],
    )
    data_samples_tissues_persons = utility.convert_records_to_dataframe(
        records=samples_tissues_persons
    )
    data_samples_tissues_persons.set_index(
        ["sample"],
        append=False,
        drop=True,
        inplace=True
    )
    data_samples_tissues_persons.rename_axis(
        index="sample",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_samples_tissues_persons.rename_axis(
        columns="properties",
        axis="columns",
        copy=False,
        inplace=True
    )
    print(data_samples_tissues_persons)
    print(data_samples_tissues_persons.iloc[0:10, :])
    print(data_samples_tissues_persons.shape)

    # Compile information.
    information = {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "collections_health_variables": bin_health,
    }

    #Write product information to file.
    write_product_sample(dock=dock, information=information)


##########
# Organization of genes' annotations.


def read_source_gene_annotation(dock=None):
    """
    Reads and organizes source information from file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_access = os.path.join(dock, "access")
    path_gene_annotation_gtex = os.path.join(
        path_access, "annotation_gene_gtex.gtf"
    )
    path_gene_annotation_gencode = os.path.join(
        path_access, "annotation_gene_gencode.gtf"
    )
    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_gene_annotation_gtex = gtfparse.read_gtf(path_gene_annotation_gtex)
    data_gene_annotation_gencode = gtfparse.read_gtf(path_gene_annotation_gencode)
    # Compile and return information.
    return {
        "data_gene_annotation_gtex": data_gene_annotation_gtex,
        "data_gene_annotation_gencode": data_gene_annotation_gencode,
    }


def summarize_raw_data_gene_annotation(
    data_gene_annotation=None,
):
    """
    Optimizes data types.

    arguments:
        data_gene_annotation (object): Pandas data frame of annotation of
            genes.

    raises:

    returns:

    """

    # Samples beginning with code person "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of genes' annotations.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' annotations from GENCODE version 31.")
    print(data_gene_annotation)
    print(data_gene_annotation.iloc[0:10, 0:10])

    pass


def extract_gene_identifier(value=None):
    """
    Removes designation of version from Ensembl gene identifier.

    arguments:
        value (str): raw identifier of gene including version

    raises:

    returns:
        (str): identifier of gene without designation of version

    """

    value_split = value.split(".")
    identifier = value_split[0]
    return identifier


def organize_genes_annotations_set(
    data_gene_annotation=None,
):
    """
    Organizes a single set of genes' annotations.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' annotations

    """

    # Summarize the structures of the raw data.
    summarize_raw_data_gene_annotation(
        data_gene_annotation=data_gene_annotation,
    )

    # Organize annotations of genes.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' annotations.")

    # Define and select relevant columns.
    print(data_gene_annotation.shape)
    columns = [
        "gene_id",
        "gene_name",
        "feature",
        "gene_type",
        "seqname",
        "start",
        "end",
        "strand",
    ]
    data_gene_annotation = data_gene_annotation.loc[
        : , data_gene_annotation.columns.isin(columns)
    ]

    print(data_gene_annotation.shape)
    print(data_gene_annotation.iloc[0:10, 0:15])
    # Remove version designations from genes' identifiers.
    data_gene_annotation["identifier"] = data_gene_annotation["gene_id"].apply(
        lambda value: extract_gene_identifier(
            value=value,
        )
    )
    data_gene_annotation.drop(
        labels="gene_id",
        axis="columns",
        inplace=True
    )

    # Remove redundant records.
    # Some genes' records are redundant.
    data_gene_annotation.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    # Organize axes.
    data_gene_annotation.set_index(
        ["identifier"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data_gene_annotation.iloc[0:10, 0:15])
    utility.print_terminal_partition(level=1)

    return data_gene_annotation


def organize_genes_annotations(
    dock=None
):
    """
    Organizes genes' annotations.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of genes' annotations.

    """

    # Remove previous files to avoid version or batch confusion.
    path_assembly_annotation = os.path.join(dock, "assembly", "annotation")
    utility.remove_directory(path=path_assembly_annotation)

    # Read source information from file.
    source = read_source_gene_annotation(dock=dock)

    data_gene_annotation_gtex = organize_genes_annotations_set(
        data_gene_annotation=source["data_gene_annotation_gtex"],
    )
    data_gene_annotation_gencode = organize_genes_annotations_set(
        data_gene_annotation=source["data_gene_annotation_gencode"],
    )

    # Compile information.
    information = {
        "data_gene_annotation_gtex": data_gene_annotation_gtex,
        "data_gene_annotation_gencode": data_gene_annotation_gencode,
    }

    # Collect garbage to clear memory.
    gc.collect()

    #Write product information to file.
    write_product_gene_annotation(dock=dock, information=information)


def access_gene_name(
    identifier=None,
    data_gene_annotation=None,
):
    """
    Combines elements in ordered pairs.

    arguments:
        identifier (str): identifier of gene
        data_gene_annotation (object): Pandas data frame of genes' annotations


    returns:
        (str): name of gene

    raises:

    """

    name = data_gene_annotation.at[identifier, "gene_name"].replace(".", "-")
    return name


##########
# Organization of genes' counts.


def read_source_gene_count(dock=None):
    """
    Reads and organizes source information from file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_access = os.path.join(dock, "access")
    path_gene_count = os.path.join(path_access, "count_gene.gct")
    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_gene_count = pandas.read_csv(
        path_gene_count,
        sep="\t",
        header=2,
        #nrows=1000,
    )
    # Compile and return information.
    return {
        "data_gene_count": data_gene_count,
    }


def summarize_raw_data_gene_count(
    data_gene_count=None,
):
    """
    Optimizes data types.

    arguments:
        data_gene_count (object): Pandas data frame of genes' counts across
            samples.

    raises:

    returns:

    """

    # Samples beginning with code person "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of genes' signals in all samples.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' counts in all samples.")
    print("Genes' counts represent reads mapped to genes.")
    print(data_gene_count)
    print(data_gene_count.iloc[0:10, 0:10])
    print("Count of genes: " + str(data_gene_count.shape[0]))
    print("Count of samples: " + str(data_gene_count.shape[1]))

    pass


def organize_genes_counts(
    dock=None,
):
    """
    Organizes counts of genes.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of genes' counts across samples

    """

    # Read source information from file.
    source = read_source_gene_count(dock=dock)

    # Summarize the structures of the raw data.
    summarize_raw_data_gene_count(
        data_gene_count=source["data_gene_count"],
    )

    # Organize genes' signals.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' counts.")

    # Organize data with names and indices.
    data_gene_count_raw = organize_data_axes_indices(
        data=source["data_gene_count"]
    )

    # Optimize data types.
    data_gene_count = convert_data_types(
        data=data_gene_count_raw,
        type="int32"
    )

    # Remove indices for compatibility with feather file format.
    data_gene_count.reset_index(
        level=None,
        inplace=True
    )

    # Compile information.
    information = {
        "data_gene_count": data_gene_count,
    }

    # Collect garbage to clear memory.
    gc.collect()

    #Write product information to file.
    write_product_gene_count(dock=dock, information=information)


##########
# Organization of genes' signals.


def read_gene_signal_samples(dock=None):
    """
    Reads and extracts sample identifiers from data for signals across samples
    and genes.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (list<str>): list of sample identifiers

    """

    # Specify directories and files.
    path_access = os.path.join(dock, "access")
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")

    # Read lines for samples with genes' signals.
    lines = utility.read_file_text_lines(
        path_file=path_gene_signal,
        start=2,
        stop=3,
    )
    line = lines[0]
    # Split line's content by delimiter.
    headers = line.split("\t")
    headers.remove("Name")
    headers.remove("Description")
    samples = headers
    return samples


def read_gene_signal_genes(dock=None):
    """
    Reads and extracts gene identifiers from data for signals across samples
    and genes.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (list<str>): list of gene identifiers

    """

    # Report.
    utility.print_terminal_partition(level=1)
    print("reading genes' identifiers from GTEx signals...")

    # Specify directories and files.
    path_access = os.path.join(dock, "access")
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")

    # Read information from file.
    identifiers_raw = utility.read_file_text_lines_elements(
        path_file=path_gene_signal,
        delimiter="\t",
        index=0,
        start=3,
        stop=1000000, # 1000000
        report=True,
    )
    # Convert identifiers.
    identifiers = list(map(
        lambda value: extract_gene_identifier(
            value=value,
        ),
        identifiers_raw
    ))
    # Report.
    utility.print_terminal_partition(level=2)
    print("count of genes: " + str(len(identifiers)))
    # Return information.
    return identifiers


def compile_variables_types(
    variables=None,
    type=None,
):
    """
    Compiles types of variables.

    arguments:
        variables (list<str>): identifiers of variables
        type (str): variable type

    raises:

    returns:
        (dict<str>): variable types

    """

    types = dict.fromkeys(variables, type)
    return types


def read_source_gene_signal(dock=None):
    """
    Reads and organizes source information from file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_access = os.path.join(dock, "access")
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")

    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)

    # Extract identifiers of samples with measurements for genes.
    samples = read_gene_signal_samples(dock=dock)
    # Designate variable types to conserve memory.
    types = compile_variables_types(
        variables=samples,
        type="float32",
    )
    types_new = {
        "Name": "str",
        "Description": "str",
    }
    types.update(types_new)

    data_gene_signal = pandas.read_csv(
        path_gene_signal,
        sep="\t",
        header=2,
        dtype=types,
        #nrows=50000,
    )
    # Compile and return information.
    return data_gene_signal


def organize_data_axes_indices(data=None):
    """
    Organizes data with names and indices.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples.

    """

    # The Pandas function to rename axis copies data deeply by default and
    # requires a lot of memory.
    utility.print_terminal_partition(level=2)

    # Remove version designations from genes' identifiers.
    print("Remove version designations from genes' identifiers.")
    data["gene"] = data["Name"].apply(
        lambda value: extract_gene_identifier(
            value=value,
        )
    )
    data.drop(
        labels="Name",
        axis="columns",
        inplace=True
    )
    # Drop column for genes' descriptions.
    print("Drop column for genes' description.")
    data.drop(
        labels="Description",
        axis="columns",
        inplace=True,
    )
    print(data.iloc[0:10, 0:10])
    print("Organize data with names and indices.")
    data.set_index(
        "gene",
        drop=True,
        inplace=True,
    )
    data.rename_axis(
        index="gene",
        axis="index",
        copy=False,
        inplace=True,
    )
    data.rename_axis(
        columns="sample",
        axis="columns",
        copy=False,
        inplace=True
    )
    print(data.iloc[0:10, 0:10])
    return data


def convert_data_types(data=None, type=None):
    """
    Optimizes data types.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples
        type (str): variable type to which to convert values in data

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples.

    """

    utility.print_terminal_partition(level=2)
    print("Optimize efficiency of storage of data.")

    # Summarize data type and size.
    print("Print technical information about data.")
    print(data.info())

    # All values are of the same type.

    if False:
        # data.select_dtypes(include=["float"])
        data_type = data.apply(pandas.to_numeric, downcast="float")

    # Determine minimal and maximal values.
    utility.print_terminal_partition(level=3)
    print("Before down cast type.")
    maximum = data.values.max()
    minimum = data.values.min()
    print("Maximum: " + str(maximum))
    print("Minimum: " + str(minimum))
    utility.print_terminal_partition(level=3)

    # Optimize data types.
    # Store genes' signals as numeric values of type float in 32 bits.
    data_type = data.astype(type)

    # Determine minimal and maximal values.
    utility.print_terminal_partition(level=3)
    print("After down cast type.")
    maximum = data_type.values.max()
    minimum = data_type.values.min()
    print("Maximum: " + str(maximum))
    print("Minimum: " + str(minimum))
    utility.print_terminal_partition(level=3)

    # Summarize data type and size.
    print("Print technical information about data.")
    print(data_type.info())

    return data_type


def summarize_raw_data_gene_signal(
    data_gene_signal=None,
):
    """
    Optimizes data types.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples.

    raises:

    returns:

    """

    # Samples beginning with code person "K-562" seem to be exceptions.
    # These samples do not seem to have any attributes or measurements.
    # These samples all seem to be for tissue "Bone Marrow".

    # Summarize the structures of the raw data.
    utility.print_terminal_partition(level=1)
    print("Summary of structures of raw data tables.")

    # Summarize table of genes' signals in all samples.
    utility.print_terminal_partition(level=2)
    print("Summary of table of genes' signals in all samples.")
    print("Genes' signals are in Transcripts Per Million (TPM).")
    print(data_gene_signal)
    print(data_gene_signal.iloc[0:10, 0:10])
    print("Count of genes: " + str(data_gene_signal.shape[0]))
    print("Count of samples: " + str(data_gene_signal.shape[1]))
    #print(source["data_gene_signal"].loc[:,"GTEX-1117F-0226-SM-5GZZ7"])

    pass


def organize_genes_signals(
    dock=None,
):
    """
    Collects tissues, persons, and genes' signals for each sample.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Remove previous files to avoid version or batch confusion.
    path_assembly_signal = os.path.join(dock, "assembly", "signal")
    utility.remove_directory(path=path_assembly_signal)

    # Read source information from file.
    samples_gtex = read_gene_signal_samples(dock=dock)
    genes_gtex = read_gene_signal_genes(dock=dock)
    data_gene_signal = read_source_gene_signal(dock=dock)

    # Report successful read.
    utility.print_terminal_partition(level=1)
    print("Read signals successfully.")

    # Summarize the structures of the raw data.
    summarize_raw_data_gene_signal(
        data_gene_signal=data_gene_signal,
    )

    # Organize genes' signals.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' signals.")

    # Organize data with names and indices.
    data_gene_signal = organize_data_axes_indices(
        data=data_gene_signal,
    )

    # Optimize data types.
    # Use this function to determine whether lower memory types are adequate
    # for values.
    # Then read in the data with the lower memory type to conserve memory.
    if False:
        data_gene_signal = convert_data_types(
            data=data_gene_signal,
            type="float32"
        )

    # Remove indices for compatibility with feather file format.
    data_gene_signal.reset_index(
        level=None,
        inplace=True
    )

    # Compile information.
    information = {
        "samples_gtex": samples_gtex,
        "genes_gtex": genes_gtex,
        "data_gene_signal": data_gene_signal
    }

    # Collect garbage to clear memory.
    gc.collect()

    # Write product information to file.
    write_product_gene_signal(dock=dock, information=information)


##########
# Association of samples, tissues, and persons


def collect_samples_persons_tissues_reference(
    data_samples_tissues_persons=None,
):
    """
    Collects person and tissue for each sample.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (dict<dict<str>>): person and tissue for each sample

    """

    data_samples = data_samples_tissues_persons.copy(deep=True)
    data_samples.reset_index(level="sample", inplace=True)
    samples_tissues_persons = utility.convert_dataframe_to_records(
        data=data_samples
    )
    reference = dict()
    for record in samples_tissues_persons:
        sample = record["sample"]
        person = record["person"]
        tissue_major = record["tissue_major"]
        tissue_minor = record["tissue_minor"]
        if sample not in reference:
            reference[sample] = dict()
            reference[sample]["person"] = person
            reference[sample]["tissue_major"] = tissue_major
            reference[sample]["tissue_minor"] = tissue_minor
    return reference


def collect_genes_samples_persons_tissues(
    reference=None,
    data_gene_sample=None
):
    """
    Collects samples, persons, and tissues for each gene's signal.

    arguments:
        reference (dict<dict<str>>): Tissue and person for each sample.
        data_gene_sample (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons.

    """

    def match_person(sample):
        return reference[sample]["person"]
    def match_tissue_major(sample):
        return reference[sample]["tissue_major"]
    def match_tissue_minor(sample):
        return reference[sample]["tissue_minor"]
    # Designate person and tissue of each sample.
    data = data_gene_sample.copy(deep=True)
    data.reset_index(level="sample", inplace=True)
    data["person"] = (
        data["sample"].apply(match_person)
    )
    data["tissue_major"] = (
        data["sample"].apply(match_tissue_major)
    )
    data["tissue_minor"] = (
        data["sample"].apply(match_tissue_minor)
    )
    data.set_index(
        ["sample", "person", "tissue_major", "tissue_minor"],
        append=False,
        drop=True,
        inplace=True
    )
    return data


def associate_samples_persons_tissues(
    data_samples_tissues_persons=None,
    data_gene_sample=None
):
    """
    Associates samples, persons, and tissues for each gene's signal.

    The table data_gene_sample needs to have samples across rows and genes
    across columns.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_sample (object): Pandas data frame of genes' signals across
            samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Associate genes' signals to persons and tissues.
    print("Associate samples to factors.")
    #print(data_samples_tissues_persons)
    #print(data_gene_sample)
    # Create reference for tissues and persons of each sample.
    reference = collect_samples_persons_tissues_reference(
        data_samples_tissues_persons=data_samples_tissues_persons
    )
    # Include identifiers for person and tissue for each sample in data for
    # genes' signals.
    data_gene_sample_person_tissue = collect_genes_samples_persons_tissues(
        reference=reference,
        data_gene_sample=data_gene_sample
    )
    #print(data_gene_signal)
    print(data_gene_sample_person_tissue.iloc[0:10, 0:7])
    # Return information.
    return data_gene_sample_person_tissue


##########
# Product.


def write_product_sample(dock=None, information=None):
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
    path_assembly = os.path.join(dock, "assembly", "sample")
    utility.create_directories(path_assembly)

    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_samples_tissues_persons_text = os.path.join(
        path_assembly, "data_samples_tissues_persons.tsv"
    )
    path_collections_health = os.path.join(
        path_assembly, "collections_health_variables.pickle"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_samples_tissues_persons"],
        path_samples_tissues_persons
    )
    information["data_samples_tissues_persons"].to_csv(
        path_or_buf=path_samples_tissues_persons_text,
        sep="\t",
        header=True,
        index=True,
    )
    with open(path_collections_health, "wb") as file_product:
        pickle.dump(information["collections_health_variables"], file_product)

    pass


def write_product_gene_annotation(dock=None, information=None):
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
    path_assembly = os.path.join(dock, "assembly", "annotation")
    utility.create_directories(path_assembly)
    path_gene_annotation_gtex = os.path.join(
        path_assembly, "data_gene_annotation_gtex.pickle"
    )
    path_gene_annotation_gtex_text = os.path.join(
        path_assembly, "data_gene_annotation_gtex.tsv"
    )
    path_gene_annotation_gencode = os.path.join(
        path_assembly, "data_gene_annotation_gencode.pickle"
    )
    path_gene_annotation_gencode_text = os.path.join(
        path_assembly, "data_gene_annotation_gencode.tsv"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_gene_annotation_gtex"],
        path_gene_annotation_gtex
    )
    pandas.to_pickle(
        information["data_gene_annotation_gencode"],
        path_gene_annotation_gencode
    )
    information["data_gene_annotation_gtex"].to_csv(
        path_or_buf=path_gene_annotation_gtex_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_gene_annotation_gencode"].to_csv(
        path_or_buf=path_gene_annotation_gencode_text,
        sep="\t",
        header=True,
        index=True,
    )

    pass


def write_product_gene_count(dock=None, information=None):
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
    path_assembly = os.path.join(dock, "assembly")
    utility.create_directory(path_assembly)
    path_gene_count = os.path.join(
        path_assembly, "data_gene_count.feather"
    )

    # Write information to file.
    information["data_gene_count"].to_feather(
        fname=path_gene_count,
    )
    pass


def write_product_gene_signal(dock=None, information=None):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files.
        information (object): information to write to file.

    raises:

    returns:

    """

    # Collect garbage to clear memory.
    gc.collect()

    # Specify directories and files.
    path_assembly = os.path.join(dock, "assembly", "signal")
    utility.create_directories(path_assembly)
    path_samples_gtex = os.path.join(
        path_assembly, "samples_gtex.pickle"
    )
    path_genes_gtex = os.path.join(
        path_assembly, "genes_gtex.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.feather"
    )

    # Write information to file.
    with open(path_samples_gtex, "wb") as file_product:
        pickle.dump(information["samples_gtex"], file_product)
    with open(path_genes_gtex, "wb") as file_product:
        pickle.dump(information["genes_gtex"], file_product)
    information["data_gene_signal"].to_feather(
        fname=path_gene_signal,
    )
    #pandas.to_pickle(
    #    information["data_gene_signal"],
    #    path=path_gene_signal,
    #    compression="gzip",
    #)
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

    # Memory conservation
    # Data for genes' signals is extensive.
    # Conserve memory.
    # Avoid unnecessary copies of the data.
    # Containerize portions of script within separate functions.
    # Collect garbage frequently.
    # Optimize data types of genes' signals.
    # Organize independent portions of data separately to avoid storing data
    # within memory unnecessarily.

    # Enable automatic garbage collection to clear memory.
    gc.enable()

    ##################################################
    ##################################################
    ##################################################

    if False:

        # Organize persons' genotypes.
        organize_persons_genotypes(dock=dock)

        # Organize associations of samples to persons and tissues.
        organize_samples_tissues_persons(dock=dock)

        # Collect garbage to clear memory.
        gc.collect()

        # Organize genes' annotations.
        organize_genes_annotations(dock=dock)

        # Collect garbage to clear memory.
        gc.collect()

        pass

    ##################################################
    ##################################################
    ##################################################

    # Organize genes' counts.
    #organize_genes_counts(dock=dock)

    # Collect garbage to clear memory.
    #gc.collect()

    ##################################################
    ##################################################
    ##################################################

    if True:

        # Organize genes' signals.
        organize_genes_signals(dock=dock)

        # Collect garbage to clear memory.
        gc.collect()

        pass

    ##################################################
    ##################################################
    ##################################################

    pass


if (__name__ == "__main__"):
    execute_procedure()
