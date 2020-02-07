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
# Source.


def read_source(dock=None):
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
    path_attribute_sample = os.path.join(path_access, "attribute_sample.txt")
    path_attribute_person = os.path.join(path_access, "attribute_person.txt")
    path_gene_annotation = os.path.join(
        path_access, "annotation_gene_gencode.gtf"
    )
    path_gene_count = os.path.join(path_access, "count_gene.gct")
    path_gene_signal = os.path.join(path_access, "signal_gene.gct")

    path_access_private = os.path.join(dock, "access_private")
    path_attribute_person_private = os.path.join(
        path_access_private, "attribute_person.txt"
    )

    path_customization = os.path.join(dock, "customization")
    path_tissues_major = os.path.join(
        path_customization, "translation_tissues_major.tsv"
    )
    path_tissues_minor = os.path.join(
        path_customization, "translation_tissues_minor.tsv"
    )

    # TODO: include translation for person attributes



    # Read information from file.
    #utility.print_file_lines(path_file=path_annotation_gene, start=0, stop=10)
    data_person_attribute = pandas.read_csv(
        path_attribute_person,
        sep="\t",
        header=0,
    )
    data_sample_attribute = pandas.read_csv(
        path_attribute_sample,
        sep="\t",
        header=0,
    )
    data_gene_annotation = gtfparse.read_gtf(path_gene_annotation)
    data_gene_count = pandas.read_csv(
        path_gene_count,
        sep="\t",
        header=2,
        #nrows=1000,
    )
    data_gene_signal = pandas.read_csv(
        path_gene_signal,
        sep="\t",
        header=2,
        #nrows=1000,
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
    # Compile and return information.
    return {
        "data_person_attribute": data_person_attribute,
        "data_sample_attribute": data_sample_attribute,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_count": data_gene_count,
        "data_gene_signal": data_gene_signal,
        "data_tissues_major": data_tissues_major,
        "data_tissues_minor": data_tissues_minor,
    }


##########
# Organization of samples' attributes, specifically associations to persons
# and tissues.


def read_gene_signal_samples(path=None):
    """
    Reads and extracts sample identifiers from genes' signals.

    arguments:
        path (str): path to file of genes' signals across samples

    raises:

    returns:
        (list<str>): list of sample identifiers

    """

    # Read lines for samples with genes' signals.
    lines = utility.read_file_text_lines(
        path_file=path,
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
    path_genotype_component_gcta = os.path.join(
        path_access_private, "relation", "gcta", "components.eigenvec"
    )
    path_genotype_component_plink = os.path.join(
        path_access_private, "relation", "plink", "components.eigenvec"
    )
    path_genotype_variance_plink = os.path.join(
        path_access_private, "relation", "plink", "components.eigenval"
    )

    # Customization
    path_customization = os.path.join(dock, "customization")
    path_tissues_major = os.path.join(
        path_customization, "translation_tissues_major.tsv"
    )
    path_tissues_minor = os.path.join(
        path_customization, "translation_tissues_minor.tsv"
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
    data_person_ancestry = pandas.read_csv(
        path_person_ancestry,
        sep="\t",
        header=0,
        low_memory=False,
    )
    data_genotype_component_gcta = pandas.read_csv(
        path_genotype_component_gcta,
        sep="\s+",
        header=None,
        names=[
            "family", "person",
            "component_1", "component_2", "component_3", "component_4",
            "component_5", "component_6", "component_7", "component_8",
            "component_9", "component_10",
        ],
        low_memory=False,
    )
    data_genotype_component_plink = pandas.read_csv(
        path_genotype_component_plink,
        sep="\s+",
        header=0,
        names=[
            "person",
            "component_1", "component_2", "component_3", "component_4",
            "component_5", "component_6", "component_7", "component_8",
            "component_9", "component_10",
        ],
        low_memory=False,
    )
    genotype_variances_plink = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_genotype_variance_plink
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
    samples = read_gene_signal_samples(path=path_gene_signal)

    # Compile and return information.
    return {
        "data_sample_attribute": data_sample_attribute,
        "data_sample_attribute_private": data_sample_attribute_private,
        "data_person_attribute": data_person_attribute,
        "data_person_attribute_private": data_person_attribute_private,
        "data_person_ancestry": data_person_ancestry,
        "data_genotype_component_gcta": data_genotype_component_gcta,
        "data_genotype_component_plink": data_genotype_component_plink,
        "genotype_variances_plink": genotype_variances_plink,
        "data_tissues_major": data_tissues_major,
        "data_tissues_minor": data_tissues_minor,
        "samples": samples,
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

    if value == "Fall":
        season = "fall"
    elif value == "Spring":
        season = "spring"
    elif value == "Summer":
        season = "summer"
    elif value == "Winter":
        season = "winter"
    else:
        season = "other"

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


def extract_person_genotypes(
    person=None,
    source=None,
    data_genotype_gcta=None,
    data_genotype_plink=None,
):
    """
    Extracts principal components from a person's genotype.

    arguments:
        person (str): identifier of a person
        source (str): source of genotype principal components, either gcta or
            plink
        data_genotype_gcta (object): Pandas data frame of principal components
            from genotypes for persons, as calculated in GCTA
        data_genotype_plink (object): Pandas data frame of principal components
            from genotypes for persons, as calculated in PLINK2

    raises:

    returns:
        (dict): principal components from a person's genotype

    """

    # Determine source of genotype principal components.
    if source == "plink":
        # Use principal components from genotype analysis in PLINK2.
        if person in data_genotype_plink.index.values.tolist():
            genotype_1 = data_genotype_plink.at[person, "component_1"]
            genotype_2 = data_genotype_plink.at[person, "component_2"]
            genotype_3 = data_genotype_plink.at[person, "component_3"]
            genotype_4 = data_genotype_plink.at[person, "component_4"]
            genotype_5 = data_genotype_plink.at[person, "component_5"]
            genotype_6 = data_genotype_plink.at[person, "component_6"]
            genotype_7 = data_genotype_plink.at[person, "component_7"]
            genotype_8 = data_genotype_plink.at[person, "component_8"]
            genotype_9 = data_genotype_plink.at[person, "component_9"]
            genotype_10 = data_genotype_plink.at[person, "component_10"]
            genotype_11 = data_genotype_plink.at[person, "component_11"]
            genotype_12 = data_genotype_plink.at[person, "component_12"]
            genotype_13 = data_genotype_plink.at[person, "component_13"]
            genotype_14 = data_genotype_plink.at[person, "component_14"]
            genotype_15 = data_genotype_plink.at[person, "component_15"]
            genotype_16 = data_genotype_plink.at[person, "component_16"]
            genotype_17 = data_genotype_plink.at[person, "component_17"]
            genotype_18 = data_genotype_plink.at[person, "component_18"]
            genotype_19 = data_genotype_plink.at[person, "component_19"]
            genotype_20 = data_genotype_plink.at[person, "component_20"]
        else:
            genotype_1 = float("nan")
            genotype_2 = float("nan")
            genotype_3 = float("nan")
            genotype_4 = float("nan")
            genotype_5 = float("nan")
            genotype_6 = float("nan")
            genotype_7 = float("nan")
            genotype_8 = float("nan")
            genotype_9 = float("nan")
            genotype_10 = float("nan")
            genotype_11 = float("nan")
            genotype_12 = float("nan")
            genotype_13 = float("nan")
            genotype_14 = float("nan")
            genotype_15 = float("nan")
            genotype_16 = float("nan")
            genotype_17 = float("nan")
            genotype_18 = float("nan")
            genotype_19 = float("nan")
            genotype_20 = float("nan")
        pass
    elif source == "gcta":
        # Use principal components from genotype analysis in GCTA.
        if person in data_genotype_gcta.index.values.tolist():
            genotype_1 = data_genotype_gcta.at[person, "component_1"]
            genotype_2 = data_genotype_gcta.at[person, "component_2"]
            genotype_3 = data_genotype_gcta.at[person, "component_3"]
            genotype_4 = data_genotype_gcta.at[person, "component_4"]
            genotype_5 = data_genotype_gcta.at[person, "component_5"]
            genotype_6 = data_genotype_gcta.at[person, "component_6"]
            genotype_7 = data_genotype_gcta.at[person, "component_7"]
            genotype_8 = data_genotype_gcta.at[person, "component_8"]
            genotype_9 = data_genotype_gcta.at[person, "component_9"]
            genotype_10 = data_genotype_gcta.at[person, "component_10"]
        else:
            genotype_1 = float("nan")
            genotype_2 = float("nan")
            genotype_3 = float("nan")
            genotype_4 = float("nan")
            genotype_5 = float("nan")
            genotype_6 = float("nan")
            genotype_7 = float("nan")
            genotype_8 = float("nan")
            genotype_9 = float("nan")
            genotype_10 = float("nan")
        pass

    # Compile and return information.
    information = {
        "genotype_1": genotype_1,
        "genotype_2": genotype_2,
        "genotype_3": genotype_3,
        "genotype_4": genotype_4,
        "genotype_5": genotype_5,
        "genotype_6": genotype_6,
        "genotype_7": genotype_7,
        "genotype_8": genotype_8,
        "genotype_9": genotype_9,
        "genotype_10": genotype_10,
        "genotype_11": genotype_11,
        "genotype_12": genotype_12,
        "genotype_13": genotype_13,
        "genotype_14": genotype_14,
        "genotype_15": genotype_15,
        "genotype_16": genotype_16,
        "genotype_17": genotype_17,
        "genotype_18": genotype_18,
        "genotype_19": genotype_19,
        "genotype_20": genotype_20,
    }
    return information


def determine_sample_associations_attributes(
    sample=None,
    data_sample_attribute=None,
    data_sample_attribute_private=None,
    data_person_attribute=None,
    data_person_attribute_private=None,
    data_person_ancestry=None,
    data_genotype_component_gcta=None,
    data_genotype_component_plink=None,
    data_tissues_major=None,
    data_tissues_minor=None,
):
    """
    Associates samples, tissues, and persons and collects their attributes.

    arguments:
        sample (str): identifier of a sample
        data_sample_attribute (object): Pandas data frame of public attributes
            for all samples
        data_sample_attribute_private (object): Pandas data frame of private
            attributes for all samples
        data_person_attribute (object): Pandas data frame of public attributes
            for persons
        data_person_attribute_private (object): Pandas data frame of private
            attributes for persons
        data_person_ancestry (object): Pandas data frame of ancestry for
            persons
        data_genotype_component_gcta (object): Pandas data frame of principal
            components from genotypes for persons, as calculated in GCTA
        data_genotype_component_plink (object): Pandas data frame of principal
            components from genotypes for persons, as calculated in PLINK2
        data_tissues_major (object): Pandas data frame of translations for
            names of major tissues
        data_tissues_minor (object): Pandas data frame of translations for
            names of minor tissues

    raises:

    returns:
        (dict): information about a sample's associations and attributes

    """

    # Access tissue attributes.
    #removal = data_sample_attribute_private.at[sample, "SMTORMVE"]
    batch_isolation = data_sample_attribute_private.at[sample, "SMNABTCH"]
    batches_analysis = data_sample_attribute_private.at[sample, "SMGEBTCH"]
    facilities = data_sample_attribute_private.at[sample, "SMCENTER"]
    #autolysis = data_sample_attribute_private.at[sample, "SMATSSCR"]
    major = data_sample_attribute_private.at[sample, "SMTS"]
    minor = data_sample_attribute_private.at[sample, "SMTSD"]
    tissue_major = data_tissues_major.at[major, "product"]
    tissue_minor = data_tissues_minor.at[minor, "product"]

    # Access person attributes.
    person = extract_gtex_sample_person_identifier(sample=sample)
    sex_raw = data_person_attribute_private.at[person, "SEX"]
    sex = translate_sex(value=sex_raw)
    decade = data_person_attribute.at[person, "AGE"]
    age = data_person_attribute_private.at[person, "AGE"]
    body = data_person_attribute_private.at[person, "BMI"]
    hardiness_raw = data_person_attribute_private.at[person, "DTHHRDY"]
    hardiness = translate_hardiness(value=hardiness_raw)
    season_raw = data_person_attribute_private.at[person, "DTHSEASON"]
    season = translate_season(value=season_raw)
    delay = data_person_attribute_private.at[person, "TRDNISCH"]
    #place = data_person_attribute_private.at[person, "DTHPLCE"]

    # Determine person's categorical ancestry.
    ancestry = determine_person_ancestry(
        person=person,
        data_person_attribute_private=data_person_attribute_private,
        data_person_ancestry=data_person_ancestry,
    )

    # Include principal components from persons' genotypes.
    genotypes = extract_person_genotypes(
        person=person,
        source="plink", # "plink" or "gcta"
        data_genotype_gcta=data_genotype_component_gcta,
        data_genotype_plink=data_genotype_component_plink,
    )

    # Compile and return information.
    information = {
        "sample": sample,
        "facilities": facilities,
        "batch_isolation": batch_isolation,
        "batches_analysis": batches_analysis,
        "tissue_major": tissue_major,
        "tissue_minor": tissue_minor,
        "person": person,
        "sex": sex,
        "decade": decade,
        "age": age,
        "ancestry": ancestry,
        "body": body,
        "hardiness": hardiness,
        "season": season,
        "delay": delay,
    }
    information.update(genotypes)
    return information


def collect_samples_tissues_persons(
    samples=None,
    data_sample_attribute=None,
    data_sample_attribute_private=None,
    data_person_attribute=None,
    data_person_attribute_private=None,
    data_person_ancestry=None,
    data_genotype_component_gcta=None,
    data_genotype_component_plink=None,
    data_tissues_major=None,
    data_tissues_minor=None,
):
    """
    Associates samples, tissues, and persons and collects their attributes.

    arguments:
        samples (list<str>): identifiers of samples with signals for genes
        data_sample_attribute (object): Pandas data frame of public attributes
            for all samples
        data_sample_attribute_private (object): Pandas data frame of private
            attributes for all samples
        data_person_attribute (object): Pandas data frame of public attributes
            for persons
        data_person_attribute_private (object): Pandas data frame of private
            attributes for persons
        data_person_ancestry (object): Pandas data frame of ancestry for
            persons
        data_genotype_component_gcta (object): Pandas data frame of principal
            components from genotypes for persons, as calculated in GCTA
        data_genotype_component_plink (object): Pandas data frame of principal
            components from genotypes for persons, as calculated in PLINK2
        data_tissues_major (object): Pandas data frame of translations for
            names of major tissues
        data_tissues_minor (object): Pandas data frame of translations for
            names of minor tissues

    raises:

    returns:
        (list<dict<str>>): information about persons and tissues for samples

    """

    # Organize data.
    data_sample_attribute.set_index(
        ["SAMPID"],
        append=False,
        drop=True,
        inplace=True
    )
    data_sample_attribute_private.set_index(
        ["SAMPID"],
        append=False,
        drop=True,
        inplace=True
    )
    data_person_attribute.set_index(
        ["SUBJID"],
        append=False,
        drop=True,
        inplace=True
    )
    data_person_attribute_private.set_index(
        ["SUBJID"],
        append=False,
        drop=True,
        inplace=True
    )
    data_person_ancestry.set_index(
        ["IID"],
        append=False,
        drop=True,
        inplace=True
    )
    data_genotype_component_gcta.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    data_genotype_component_plink.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    data_tissues_major.set_index(
        ["source"],
        append=False,
        drop=True,
        inplace=True
    )
    data_tissues_minor.set_index(
        ["source"],
        append=False,
        drop=True,
        inplace=True
    )
    # Collect tissues and persons for each sample.
    samples_tissues_persons = list()
    for sample in samples:
        record = determine_sample_associations_attributes(
            sample=sample,
            data_sample_attribute=data_sample_attribute,
            data_sample_attribute_private=data_sample_attribute_private,
            data_person_attribute=data_person_attribute,
            data_person_attribute_private=data_person_attribute_private,
            data_person_ancestry=data_person_ancestry,
            data_genotype_component_gcta=data_genotype_component_gcta,
            data_genotype_component_plink=data_genotype_component_plink,
            data_tissues_major=data_tissues_major,
            data_tissues_minor=data_tissues_minor,
        )
        samples_tissues_persons.append(record)
    # Return information.
    return samples_tissues_persons


def organize_genotype_variance(
    genotype_variances_plink=None,
):
    """
    Associates samples, tissues, and persons and collects their attributes.

    arguments:
        genotype_variances_plink (object): variances of genotypes' principal
            components, as calculated in PLINK2

    raises:

    returns:
        (object): Pandas data frame of variance of each principal component

    """

    count = 1
    records = list()
    for variance in genotype_variances_plink:
        record = dict()
        record["variance"] = variance
        record["component"] = count
        count += 1
        pass
    data = utility.convert_records_to_dataframe(
        records=records
    )
    return data


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

    # Read source information from file.
    source = read_source_sample(dock=dock)

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
    # Collect information about samples.
    samples_tissues_persons = collect_samples_tissues_persons(
        samples=source["samples"],
        data_person_attribute=source["data_person_attribute"],
        data_person_attribute_private=source["data_person_attribute_private"],
        data_person_ancestry=source["data_person_ancestry"],
        data_genotype_component_gcta=source["data_genotype_component_gcta"],
        data_genotype_component_plink=source["data_genotype_component_plink"],
        data_sample_attribute=source["data_sample_attribute"],
        data_sample_attribute_private=source["data_sample_attribute_private"],
        data_tissues_major=source["data_tissues_major"],
        data_tissues_minor=source["data_tissues_minor"],
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

    # Organize variance of principal components on genotype.
    data_genotype_variance_plink = organize_genotype_variance(
        genotype_variances_plink=source["genotype_variances_plink"],
    )

    # Compile information.
    information = {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_genotype_variance_plink": data_genotype_variance_plink,
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


def extract_gene_identifier(string):
    """
    Removes designation of version from Ensembl gene identifier.

    arguments:
        string (str): Ensembl identifier of gene.

    raises:

    returns:
        (str): Ensemble identifier of gene without designation of version.

    """

    string_split = string.split(".")
    identifier_gene = string_split[0]
    return identifier_gene


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
    ]
    data_gene_annotation = data_gene_annotation.loc[
        : , data_gene_annotation.columns.isin(columns)
    ]

    print(data_gene_annotation.shape)
    print(data_gene_annotation.iloc[0:10, 0:15])
    # Remove version designations from genes' identifiers.
    data_gene_annotation["identifier"] = (
        data_gene_annotation["gene_id"].apply(extract_gene_identifier)
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
    data["gene"] = data["Name"].apply(extract_gene_identifier)
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
    samples = read_gene_signal_samples(path=path_gene_signal)
    # Designate variable types to conserve memory.
    types = compile_variables_types(
        variables=samples,
        type="float64",
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
        #dtype=types,
        #nrows=1000,
    )
    # Compile and return information.
    return {
        "data_gene_signal": data_gene_signal,
    }


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

    # Read source information from file.
    source = read_source_gene_signal(dock=dock)

    # Summarize the structures of the raw data.
    summarize_raw_data_gene_signal(
        data_gene_signal=source["data_gene_signal"],
    )

    # Organize genes' signals.
    utility.print_terminal_partition(level=1)
    print("Organization of genes' signals.")

    # Organize data with names and indices.
    data_gene_signal_raw = organize_data_axes_indices(
        data=source["data_gene_signal"]
    )

    # Optimize data types.
    data_gene_signal = convert_data_types(
        data=data_gene_signal_raw,
        type="float32"
    )

    # Remove indices for compatibility with feather file format.
    data_gene_signal.reset_index(
        level=None,
        inplace=True
    )

    # Compile information.
    information = {
        "data_gene_signal": data_gene_signal
    }

    # Collect garbage to clear memory.
    gc.collect()

    #Write product information to file.
    write_product_gene_signal(dock=dock, information=information)


def transform_gene_signal_log(
    data_gene_signal=None,
    pseudo_count=None,
):
    """
    Transforms values of genes' signals to base-two logarithmic space.

    arguments:
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.
        pseudo_count (float): value to add to signal before transformation

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific persons and tissues.

    """

    # Transform genes' signals to base-2 logarithmic space.
    # Transform before calculation of any median or mean values. <- maybe false
    # Logarithms are not distributive.
    # To accommodate gene signals of 0.0, add a pseudo count of 2.0 to all
    # counts before calculation of the base-2 logarithm.
    # log-2 signal = log2(TPM + pseudo-count)
    # pseudo-count = 1.0
    if False:
        utility.print_terminal_partition(level=2)
        print("Transformation of genes' signals to base-2 logarithmic space.")
        print(
            "To accommodate gene signals of 0.0, add a pseudo count of " +
            "1.0-5.0 to all counts before calculation of the base-2 logarithm."
        )
    data_gene_signal_log = calculate_logarithm_gene_signal(
        pseudo_count=pseudo_count,
        data_gene_signal=data_gene_signal
    )
    #print(data_gene_signal_log.iloc[0:10, 0:10])

    return data_gene_signal_log


def calculate_logarithm_gene_signal(pseudo_count=None, data_gene_signal=None):
    """
    Calculates the base-2 logarithm of genes' signals in each sample.

    Original gene signals are in transcript counts per million (TPM).

    To accommodate gene signals of 0.0, add a pseudo count of 1.0 to all counts
    before calculation of the base-2 logarithm.

    arguments:
        pseudo_count (float): Pseudo count to add to gene signal before
            transformation to avoid values of zero.
        data_gene_signal (object): Pandas data frame of signals for all genes
            across specific persons and tissues.

    raises:

    returns:
        (object): Pandas data frame of base-2 logarithmic signals for all genes
            across specific persons and tissues.

    """

    # lambda x: math.log((x + 1), 2)

    # An alternative approach would be to set the label columns to indices.
    if False:
        data_gene_signal_index = data_gene_signal.set_index(
            ["Name", "Description"], append=True, drop=True
        )
        data_gene_signal_log = data_gene_signal_index.apply(
            lambda value: 2 * value
        )
        data_log = data_signal_index.copy()
        data_log.iloc[0:, 2:] = data_log.iloc[0:, 2:].applymap(
            lambda value: math.log((value + 1.0), 2)
        )

    data_log = data_gene_signal.applymap(
        lambda value: math.log((value + pseudo_count), 2)
    )
    # Reverse an index to a column.
    return data_log


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
    path_assembly = os.path.join(dock, "assembly")
    utility.create_directory(path_assembly)
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_samples_tissues_persons_text = os.path.join(
        path_assembly, "data_samples_tissues_persons.tsv"
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
    path_assembly = os.path.join(dock, "assembly")
    utility.create_directory(path_assembly)
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
    path_assembly = os.path.join(dock, "assembly")
    utility.create_directory(path_assembly)
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.feather"
    )

    # Write information to file.
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

    # Remove previous files to avoid version or batch confusion.
    path_assembly = os.path.join(dock, "assembly")
    utility.remove_directory(path=path_assembly)

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

    # Organize associations of samples to persons and tissues.
    organize_samples_tissues_persons(dock=dock)

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    # Organize genes' annotations.
    organize_genes_annotations(dock=dock)

    # Collect garbage to clear memory.
    gc.collect()

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

    # Organize genes' signals.
    organize_genes_signals(dock=dock)

    # Collect garbage to clear memory.
    gc.collect()

    ##################################################
    ##################################################
    ##################################################

    pass


if (__name__ == "__main__"):
    execute_procedure()
