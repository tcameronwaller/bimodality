"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard
import os
import csv
import copy
import textwrap
import string
import gzip
import shutil
import textwrap
import itertools
import math

# Relevant

import pandas
import sklearn
import scipy
import numpy
import statsmodels.api

# Custom

#dir()
#importlib.reload()

###############################################################################
# Functionality


# General.


def convert_string_low_alpha_num(characters):
    """
    Converts string of characters to lower case with only alphabetical or
    numerical characters.

    arguments:
        characters (str): characters in a string

    raises:

    returns:
        (str): characters in a string

    """

    # Convert all characters to lower case.
    characters_lower = characters.lower()
    # Remove all characters other than alphabetical or numerical characters.
    characters_novel = characters_lower
    for character in characters_lower:
        if (
            (character not in string.ascii_letters) and
            (character not in string.digits)
        ):
            characters_novel = characters_novel.replace(character, "")
    return characters_novel


def extract_subdirectory_names(path=None):
    """
    Extracts names of subdirectories.

    Reads all contents of a directory and returns only those of subdirectories.

    arguments:
        path (str): path to directory

    raises:

    returns:
        (list<str>): names of subdirectories

    """

    contents = os.listdir(path=path)
    names = list(filter(
        lambda content: os.path.isdir(os.path.join(path, content)),
        contents
    ))
    return names


def create_directories(path=None):
    """
    Creates a directory if it does not already exist.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if not os.path.exists(path):
        os.makedirs(path)


def remove_file(path=None):
    """
    Removes a file if it exists.

    arguments:
        path (str): path to file

    raises:

    returns:

    """

    if os.path.exists(path):
        os.remove(path)


def remove_directory(path=None):
    """
    Removes a directory if it exists.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if os.path.exists(path):
        # Only for empty directory.
        #os.rmdir(path)
        # For non empty directory.
        shutil.rmtree(path)


def remove_empty_directory(path=None):
    """
    Removes a directory if it is empty.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if (os.path.exists(path)) and (len(os.listdir(path)) < 1):
        os.rmdir(path)


def create_directory(path=None):
    """
    Creates a directory if it is does not already exist.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if not (os.path.exists(path)):
        os.mkdir(path)


def compress_file_gzip(path=None):
    """
    Copies and compresses a file by gzip.

    arguments:
        path (str): path to file

    raises:

    returns:

    """

    with open(path, "rb") as file_source:
        with gzip.open((path + ".gz"), "wb") as file_product:
            shutil.copyfileobj(file_source, file_product)
    pass


def decompress_file_gzip(path=None):
    """
    Copies and decompresses a file by gzip.

    Removes any ".gz" suffix from the file name.

    arguments:
        path (str): path to file

    raises:

    returns:

    """

    split_strings = path.split(".")
    if split_strings[-1] == "gz":
        path_novel = ".".join(split_strings[0:-1])
    else:
        path_novel = path
    with gzip.open(path, "rb") as file_source:
        with open(path_novel, "wb") as file_product:
            shutil.copyfileobj(file_source, file_product)
    pass


def print_terminal_partition(level=None):
    """
    Prints string to terminal to partition reports for clarity.

    arguments:
        level (int): level or type of partition to create and print.

    raises:

    returns:

    """

    if level == 1:
        partition = textwrap.dedent("""\



            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



        """)
    elif level == 2:
        partition = textwrap.dedent("""\



            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



        """)
    elif level == 3:
        partition = textwrap.dedent("""\



            ----------
            ----------
            ----------



        """)
    elif level == 4:
        partition = textwrap.dedent("""\

            ----------

        """)
    else:
        partition = ""
    print(partition)
    pass


def read_file_text_table(path_file=None, names=None, delimiter=None):
    """
    Reads and organizes source information from file

    This function reads and organizes relevant information from file.

    arguments:
        path_file (str): path to directory and file
        names (list<str>): names for values in each row of table
        delimiter (str): delimiter between values in the table

    returns:
        (list<dict>): tabular information from file

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    with open(path_file, "r") as file_source:
        reader = csv.DictReader(
            file_source, fieldnames=names, delimiter=delimiter
        )
        information = list(map(lambda row: dict(row), list(reader)))
    # Return information
    return information


def read_file_text(
    path_file=None
):
    """
    Reads information from file as a text string.

    arguments:
        path_file (str): path to directory and file

    returns:
        (str): information from file

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    with open(path_file, "r") as file_source:
        content = file_source.read()
    # Return information
    return content


def read_file_text_list(
    delimiter=None,
    path_file=None
):
    """
    Reads and organizes source information from file.

    Delimiters include "\n", "\t", ",", " ".

    arguments:
        delimiter (str): delimiter between elements in list
        path_file (str): path to directory and file

    returns:
        (list<str>): information from file

    raises:

    """

    # Read information from file
    content = read_file_text(path_file=path_file)
    # Split content by line delimiters.
    values_split = content.split(delimiter)
    values_strip = list(map(lambda value: value.strip(), values_split))
    # Return information
    return values_strip


def read_file_text_lines(
    path_file=None,
    start=None,
    stop=None,
):
    """
    Reads a range of lines from a text file.

    arguments:
        path_file (str): path to directory and file
        start (int): index of line to start reading
        stop (int): index of line to stop reading

    returns:
        (list<str>): information from file

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    lines = list()
    with open(path_file, "r") as file_source:
        line = file_source.readline()
        count = 0
        while line:
            if (start <= count and count < stop):
                lines.append(line.strip())
            line = file_source.readline()
            count += 1
    # Return information.
    return lines


def read_file_text_lines_elements(
    path_file=None,
    delimiter=None,
    index=None,
    start=None,
    stop=None,
    report=None,
):
    """
    Reads all lines from a text file and collects elements from a single
    index within each line.

    arguments:
        path_file (str): path to directory and file
        delimiter (str): delimiter between elements in list
        index (int): index of element to collect from each line
        start (int): index of line to start reading
        stop (int): index of line to stop reading
        report (bool): whether to print element from each line

    returns:
        (list<str>): information from file

    raises:

    """

    # Read information from file.
    elements = list()
    with open(path_file, "r") as file_source:
        line = file_source.readline()
        count = 0
        while line:
            if (start <= count and count < stop):
                row = line.strip().split(delimiter)
                element = row[index]
                elements.append(element)
                if report:
                    print("line " + str(count) + ": " + str(element))
            line = file_source.readline()
            count += 1
    # Return information.
    return elements


def write_file_text_table(
    information=None,
    path_file=None,
    names=None,
    delimiter=None,
    header=None,
):
    """
    Writes to file in text format information from a list of dictionaries.

    arguments:
        information (list<str>): information
        path_file (str): path to directory and file
        names (list<str>): names for values in each row of table
        delimiter (str): delimiter between values in the table
        header (bool): whether to write column headers in file

    returns:

    raises:

    """

    # Write information to file
    #with open(out_file_path_model, "w") as out_file:
    #    out_file.write(content_identifier)
    with open(path_file, "w") as file_product:
        writer = csv.DictWriter(
            file_product, fieldnames=names, delimiter=delimiter
        )
        if header:
            writer.writeheader()
        writer.writerows(information)


def write_file_text_list(
    elements=None,
    delimiter=None,
    path_file=None
):
    """
    Writes to file in text format information from an array of strings.

    Delimiters include "\n", "\t", ",", " ".

    arguments:
        elements (list<str>): character elements
        delimiter (str): delimiter between elements in list
        path_file (str): path to directory and file

    returns:

    raises:

    """

    # Write information to file
    with open(path_file, "w") as file_product:
        string = delimiter.join(elements)
        file_product.write(string)
    pass


def print_file_lines(path_file=None, start=None, stop=None):
    """
    Reads and prints a specific range of lines from a file.

    arguments:
        path_file (str): path to directory and file
        start (int): index of line in file at which to start reading
        stop (int): index of line in file at which to stop reading

    returns:

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    with open(path_file, "r") as file_source:
        line = file_source.readline()
        count = 0
        while line and (start <= count and count < stop):
            print("***line {}***: {}".format(count, line.strip()))
            line = file_source.readline()
            count += 1


def find(match=None, sequence=None):
    """
    Finds the first element in a sequence to match a condition, otherwise none

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (object | NoneType): first element from sequence to match condition or
            none

    raises:

    """

    for element in sequence:
        if match(element):
            return element
    return None


def find_index(match=None, sequence=None):
    """
    Finds index of first element in sequence to match a condition, otherwise -1

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (int): index of element if it exists or -1 if it does not exist

    raises:

    """

    for index, element in enumerate(sequence):
        if match(element):
            # Element matches condition
            # Return element's index
            return index
    # Not any elements match condition
    # Return -1
    return -1


def find_all(match=None, sequence=None):
    """
    Finds all elements in a sequence to match a condition, otherwise none.

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (list<dict> | NoneType): elements from sequence to match condition or
            none

    raises:

    """

    matches = []
    for element in sequence:
        if match(element):
            matches.append(element)
    if len(matches) > 0:
        return matches
    else:
        return None


def select_elements_by_sets(
    names=None,
    sets=None,
    count=None,
):
    """
    Selects unique elements that belong to at a minimal count of specific sets.

    This altorithm assumes that each element can belong to multiple sets but
    that each set's elements are unique or nonredundant.

    Data format
    - elements (dict)
    -- element (list)
    --- set (str)

    arguments:
        names (list<str>): names of sets to consider
        sets (dict<list<str>>): sets of elements
        count (int): minimal count of sets to which an element must belong

    returns:
        (list): unique elements that belong to minimal count of specific sets

    raises:

    """

    # Collect sets to which each element belongs.
    elements = dict()
    for name in names:
        for element in collect_unique_elements(elements_original=sets[name]):
            if element not in elements:
                elements[element] = list()
                elements[element].append(name)
            else:
                elements[element].append(name)
        pass
    pass

    # Filter elements by count of sets.
    passes = list()
    for element in elements:
        if len(elements[element]) >= count:
            passes.append(element)

    # Return elements that pass filter.
    return passes


def collect_unique_elements(elements_original=None):
    """
    Collects unique elements

    arguments:
        elements_original (list): sequence of elements

    returns:
        (list): unique elements

    raises:

    """

    elements_novel = []
    for element in elements_original:
        if element not in elements_novel:
            elements_novel.append(element)
    return elements_novel


def collect_value_from_records(key=None, records=None):
    """
    Collects a single value from multiple records

    arguments:
        key (str): key of value in each record
        records (list<dict>): sequence of records

    raises:

    returns:
        (list): values from records

    """

    def access(record):
        return record[key]
    return list(map(access, records))


def collect_values_from_records(key=None, records=None):
    """
    Collects values from multiple records.

    arguments:
        key (str): key of value in each record
        records (list<dict>): sequence of records

    raises:

    returns:
        (list): values from records

    """

    collection = []
    for record in records:
        collection.extend(record[key])
    return collection


def compare_lists_by_inclusion(list_one=None, list_two=None):
    """
    Compares lists by inclusion.

    Returns True if all elements in list_two are in list_one.

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (bool): whether first list includes all elements from second

    raises:

    """

    def match(element_two=None):
        return element_two in list_one
    matches = list(map(match, list_two))
    return all(matches)


def compare_lists_by_mutual_inclusion(list_one=None, list_two=None):
    """
    Compares lists by mutual inclusion

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (bool): whether each list includes all elements from the other

    raises:

    """

    forward = compare_lists_by_inclusion(
        list_one=list_one,
        list_two=list_two
    )
    reverse = compare_lists_by_inclusion(
        list_one=list_two,
        list_two=list_one
    )
    return forward and reverse


def filter_common_elements(list_one=None, list_two=None):
    """
    Filters elements by whether both of two lists include them

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (list): elements that both of two lists include

    raises:

    """

    def match(element_two=None):
        return element_two in list_one
    return list(filter(match, list_two))


def filter_unique_union_elements(list_one=None, list_two=None):
    """
    Filters unique elements from union of two lists.

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (list): elements that both of two lists include

    raises:

    """

    union = list_one + list_two
    unique = collect_unique_elements(elements_original=union)
    return unique


def filter_unique_exclusion_elements(
    elements_exclusion=None,
    elements_total=None
):
    """
    Filters unique elements by exclusion.

    arguments:
        elements_exclusion (list): list of elements
        elements_total (list): list of elements

    returns:
        (list): elements from total list not in exclusion list

    raises:

    """

    elements_novel = []
    for element in elements_total:
        if element not in elements_exclusion:
            elements_novel.append(element)
    elements_unique = collect_unique_elements(elements_original=elements_novel)
    return elements_unique


def collect_records_targets_by_categories(
    target=None,
    category=None,
    records=None
):
    """
    Collects values of a target attribute that occur together in records with
    each value of another category attribute.
    Each record has a single value of the target attribute.
    Each record can have either a single value or multiple values of the
    category attribute.
    These collections do not necessarily include only unique values of the
    target attribute.

    arguments:
        target (str): name of attribute in records to collect for each category
        category (str): name of attribute in records to define categories
        records (list<dict>): records with target and category attributes

    raises:

    returns:
        (dict<list<str>>): values of the target attribute that occur together
            in records with each value of the category attribute

    """

    def collect_record_target_by_category(
        target_value=None,
        category_value=None,
        collection_original=None,
    ):
        collection_novel = copy.deepcopy(collection_original)
        # Determine whether collection includes the category's value.
        if category_value in collection_novel.keys():
            # Collection includes the category's value.
            target_values = collection_novel[category_value]
            target_values.append(target_value)
            # Include target's value in collection.
            collection_novel[category_value] = target_values
        else:
            # Collection does not include the category's value.
            # Include category's value and target's value in collection.
            collection_novel[category_value] = [target_value]
        return collection_novel

    collection = {}
    for record in records:
        target_value = record[target]
        category_values = record[category]
        if isinstance(category_values, list):
            for category_value in category_values:
                collection = collect_record_target_by_category(
                    target_value=target_value,
                    category_value=category_value,
                    collection_original=collection
                )
        else:
            category_value = category_values
            collection = collect_record_target_by_category(
                target_value=target_value,
                category_value=category_value,
                collection_original=collection
            )
    return collection


def collect_values_from_records_in_reference(
    key=None, identifiers=None, reference=None
):
    """
    Collects a single value from a specific record in a reference.

    arguments:
        key (str): key of value in record
        identifiers (list<str>): identifiers of records in reference
        reference (dict<dict<str>>): reference of records

    raises:

    returns:
        (list<str>): values from records

    """

    values = []
    for identifier in identifiers:
        record = reference[identifier]
        value = record[key]
        values.append(value)
    return values


def filter_nonempty_elements(elements_original=None):
    """
    Filters nonempty elements.

    arguments:
        elements_original (list<str>): sequence of elements

    returns:
        (list<str>): non-empty elements

    raises:

    """

    elements_novel = []
    for element in elements_original:
        if len(str(element)) > 0:
            elements_novel.append(element)
    return elements_novel


def filter_entries_identifiers(
    identifiers=None,
    entries_original=None
):
    """
    Filters nodes and links by identifiers.

    arguments:
        identifiers (list<str>): identifiers of elements to keep
        entries_original (dict<dict>): entries

    raises:

    returns:
        (dict<dict>): entries

    """

    entries_novel = {}
    for entry in entries_original.values():
        if entry["identifier"] in identifiers:
            entries_novel[entry["identifier"]] = entry
    return entries_novel


def calculate_standard_score(
    value=None,
    mean=None,
    deviation=None,
):
    """
    Calculates the standard score, z-score, of a value.

    arguments:
        value (float): value to transform to standard score space
        mean (float): mean of distribution from which value comes
        deviation (float): standard deviation of distribution from which value
            comes

    raises:

    returns:
        (float): value in standard score space

    """

    if deviation != 0:
        return ((value - mean) / deviation)
    else:
        return float("nan")


def calculate_standard_scores(
    values=None,
    mean=None,
    deviation=None,
):
    """
    Calculates the standard scores, z-scores, of values.

    arguments:
        value (list<float>): values to transform to standard score space
        mean (float): mean of distribution from which value comes
        deviation (float): standard deviation of distribution from which value
            comes

    raises:

    returns:
        (list<float>): value in standard score space

    """

    values_standard = list()
    for value in values:
        value_standard = calculate_standard_score(
            value=value,
            mean=mean,
            deviation=deviation,
        )
        values_standard.append(value_standard)
    return values_standard


def combine_unique_elements_pairwise_orderless(
    elements=None,
):
    """
    Combines elements in orderless pairs.

    ABCD: AB AC AD BC BD CD

    arguments:
        elements (list): elements of interest

    returns:
        (list<tuple>): unique pairs of elements

    raises:

    """

    # Select unique elements.
    elements_unique = collect_unique_elements(elements_original=elements)
    # Combine elements in pairs.
    pairs = list(itertools.combinations(elements_unique, 2))
    # Return information.
    return pairs


def combine_unique_elements_pairwise_order(
    elements=None,
):
    """
    Combines elements in ordered pairs.

    ABCD: AB AC AD BA BC BD CA CB CD DA DB DC

    arguments:
        elements (list): elements of interest

    returns:
        (list<tuple>): unique pairs of elements

    raises:

    """

    # Select unique elements.
    elements_unique = collect_unique_elements(elements_original=elements)
    # Combine elements in pairs.
    pairs = list(itertools.permutations(elements_unique, 2))
    # Return information.
    return pairs


# Pandas

def convert_records_to_dataframe(records=None):
    """
    Converts information from list of dictionaries to Pandas data frame.

    arguments:
        records (list<dict>): records.

    raises:

    returns:
        (object): Pandas data frame.

    """

    return pandas.DataFrame(data=records)


def convert_dataframe_to_records(data=None):
    """
    Converts information from Pandas data frame to list of dictionaries.

    arguments:
        data (object): Pandas data frame.

    raises:

    returns:
        (list<dict>): records.

    """

    data_copy = data.copy(deep=True)
    data_copy.reset_index(
        level=None,
        inplace=True,
    )
    return data_copy.to_dict(orient="records")


def segregate_data_two_thresholds(
    data=None,
    abscissa=None,
    ordinate=None,
    threshold_abscissa=None,
    selection_abscissa=None,
    threshold_ordinate=None,
    selection_ordinate=None,
):
    """
    Segregates data by values against thresholds on two dimensions.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with independent variable
        ordinate (str): name of data column with dependent variable
        threshold_abscissa (float): threshold for abscissa
        selection_abscissa (str): selection criterion for abscissa's values
            against threshold
        threshold_ordinate (float): threshold for ordinate
        selection_ordinate (str): selection criterion for ordinate's values
            against threshold

    raises:

    returns:
        (dict): collection of data that pass and fail against thresholds

    """

    # Copy data.
    data_pass = data.copy(deep=True)
    data_fail = data.copy(deep=True)
    # Abscissa.
    if selection_abscissa == ">=":
        data_pass = data_pass.loc[
            data_pass[abscissa] >= threshold_abscissa,
            :
        ]
    elif selection_abscissa == "<=":
        data_pass = data_pass.loc[
            data_pass[abscissa] <= threshold_abscissa,
            :
        ]
    # Ordinate.
    if selection_ordinate == ">=":
        data_pass = data_pass.loc[
            data_pass[ordinate] >= threshold_ordinate,
            :
        ]
    elif selection_ordinate == "<=":
        data_pass = data_pass.loc[
            data_pass[ordinate] <= threshold_ordinate,
            :
        ]

    # Fail.
    data_fail = data_fail.loc[
        ~data_fail.index.isin(data_pass.index),
        :
    ]

    # Return information.
    collection = dict()
    collection["pass"] = data_pass
    collection["fail"] = data_fail
    return collection


def filter_rows_columns_by_threshold_proportion(
    data=None,
    dimension=None,
    threshold=None,
    proportion=None,
):
    """
    Filters either rows or columns by whether a certain proportion of values
    are greater than a minimal thresholds.

    Persistence of a row or column requires at least a specific proportion of
    values beyond a specific threshold.

    Filter rows by consideration of values across columns in each row.
    Filter columns by consideration of values across rows in each column.

    arguments:
        data (object): Pandas data frame of values
        dimension (str): dimension to filter, either "row" or "column"
        threshold (float): minimal signal
        proportion (float): minimal proportion of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of values

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Determine count from proportion.
    if dimension == "row":
        # Filter rows by consideration of columns for each row.
        columns = data.shape[1]
        count = round(proportion * columns)
    elif dimension == "column":
        # Filter columns by consideration of rows for each columns.
        rows = data.shape[0]
        count = round(proportion * rows)

    # Determine whether values exceed threshold.
    data_threshold = (data >= threshold)
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    # This aggregation operation produces a series.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )

    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count, : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count]
    return data_pass


def filter_rows_columns_by_threshold_outer_proportion(
    data=None,
    dimension=None,
    threshold_high=None,
    threshold_low=None,
    proportion=None,
):
    """
    Filters either rows or columns by whether a certain proportion of values
    are outside of low and high thresholds.

    Persistence of a row or column requires at least a specific proportion of
    values beyond a specific threshold.

    Filter rows by consideration of values across columns in each row.
    Filter columns by consideration of values across rows in each column.

    arguments:
        data (object): Pandas data frame of values
        dimension (str): dimension to filter, either "row" or "column"
        threshold_high (float): value must be greater than this threshold
        threshold_low (float): value must be less than this threshold
        proportion (float): minimal proportion of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of values

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Copy data.
    data = data.copy(deep=True)

    # Determine count from proportion.
    if dimension == "row":
        # Filter rows by consideration of columns for each row.
        columns = data.shape[1]
        count = round(proportion * columns)
    elif dimension == "column":
        # Filter columns by consideration of rows for each columns.
        rows = data.shape[0]
        count = round(proportion * rows)

    # Determine whether values exceed threshold.
    data_threshold = data.applymap(
        lambda value: ((value <= threshold_low) or (threshold_high <= value))
    )
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    # This aggregation operation produces a series.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )

    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count, : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count]
    return data_pass


def calculate_false_discovery_rate(
    threshold=None,
    probability=None,
    discovery=None,
    significance=None,
    data_probabilities=None,
):
    """
    Calculates false discovery rates from probabilities.

    arguments:
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries
        probability (str): name of column of probabilities
        discovery (str): name of column of false discovery rates
        significance (str): name of column of significance
        data_probabilities (object): Pandas data frame of probabilities

    raises:

    returns:
        (object): Pandas data frame of probabilities and false discovery rates

    """

    # Copy data.
    data_copy = data_probabilities.copy(deep=True)

    # False discovery rate method cannot accommodate missing values.
    # Remove null values.
    data_null_boolean = pandas.isna(data_copy[probability])
    data_null = data_copy.loc[data_null_boolean]
    data_valid = data_copy.dropna(
        axis="index",
        how="any",
        inplace=False,
    )
    # Calculate false discovery rates from probabilities.
    probabilities = data_valid[probability].to_numpy()
    report = statsmodels.stats.multitest.multipletests(
        probabilities,
        alpha=threshold,
        method="fdr_bh",
        is_sorted=False,
    )
    significances = report[0]
    discoveries = report[1]
    data_valid[significance] = significances
    data_valid[discovery] = discoveries
    data_null[significance] = False
    data_null[discovery] = float("nan")

    # Combine null and valid portions of data.
    data_discoveries = data_valid.append(
        data_null,
        ignore_index=False,
    )

    # Return information.
    return data_discoveries


def calculate_pseudo_logarithm_signals(
    pseudo_count=None,
    base=None,
    data=None
):
    """
    Adds a pseudo count to signals and then calculates their logarithm.

    arguments:
        pseudo_count (float): Pseudo count to add to gene signal before
            transformation to avoid values of zero
        base (float): logarithmic base
        data (object): Pandas data frame of signals

    raises:

    returns:
        (object): Pandas data frame of base-2 logarithmic signals for all genes
            across specific persons and tissues.

    """

    data_log = data.applymap(
        lambda value: math.log((value + pseudo_count), base)
    )
    return data_log


def count_data_factors_groups_elements(
    factors=None,
    element=None,
    count=None,
    data=None,
):
    """
    Counts elements in groups by factor.

    arguments:
        factors (list<str>): names of factor columns in data
        element (str): name of column in data for elements to count
        count (str): name for column of counts in count data
        data (object): Pandas data frame of elements in groups by factors

    raises:

    returns:
        (object): Pandas data frame of counts of elements in each factor group
    """

    # Copy data.
    data_copy = data.copy(deep=True)

    # Organize data.
    data_copy.reset_index(
        level=None,
        inplace=True,
    )
    columns = factors
    columns = copy.deepcopy(factors)
    columns.append(element)
    data_columns = data_copy.loc[
        :, data_copy.columns.isin(columns)
    ]
    data_columns.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_columns.set_index(
        factors,
        append=False,
        drop=True,
        inplace=True
    )
    # Count rows (elements) in each factor group of the data.
    # This process is similar to iteration on groups and collection of the
    # counts of elements in each group.
    # Notice the use of ".size()" instead of ".shape[0]" to count rows.
    # As the groups are series, ".size()" counts the elements properly.
    data_counts = data_columns.groupby(
        level=factors,
        sort=True,
        as_index=False
    ).size().to_frame(
        name=count
    )
    data_counts.reset_index(
        level=None,
        inplace=True,
    )
    return data_counts


# Principal components

# TODO: switch to statsmodels PCA
# https://www.statsmodels.org/devel/generated/statsmodels.multivariate.pca.pca.html

def calculate_principal_components(
    data=None,
    components=None,
    report=None,
):
    """
    Calculates the principal components.

    Format of data should have features across columns and observations across
    rows.

    This function then transforms the data frame to a matrix.

    This matrix has observations across dimension zero and features across
    dimension one.

    arguments:
        data (object): Pandas data frame of signals with features across rows
            and observations across columns
        components (int): count of principle components
        report (bool): whether to print reports to terminal

    raises:

    returns:
        (dict): information about data's principal components

    """

    # Organize data.
    data_copy = data.copy(deep=True)

    # Organize data for principle component analysis.
    # Organize data as an array of arrays.
    matrix = data_copy.to_numpy()
    # Execute principle component analysis.
    pca = sklearn.decomposition.PCA(n_components=components)
    #report = pca.fit_transform(matrix)
    pca.fit(matrix)
    # Report extent of variance that each principal component explains.
    variance_ratios = pca.explained_variance_ratio_
    component_numbers = range(len(variance_ratios) + 1)
    variance_series = {
        "component": component_numbers[1:],
        "variance": variance_ratios
    }
    data_component_variance = pandas.DataFrame(data=variance_series)
    # Transform data by principal components.
    matrix_component = pca.transform(matrix)
    # Match components to samples.
    observations = data_copy.index.to_list()
    records = []
    for index, row in enumerate(matrix_component):
        record = {}
        record["observation"] = observations[index]
        for count in range(components):
            name = "component_" + str(count + 1)
            record[name] = row[count]
        records.append(record)
    data_component = (
        convert_records_to_dataframe(records=records)
    )
    data_component.set_index(
        ["observation"],
        append=False,
        drop=True,
        inplace=True
    )
    # Report
    if report:
        print_terminal_partition(level=3)
        print("dimensions of original data: " + str(data.shape))
        print("dimensions of novel data: " + str(data_component.shape))
        print("proportional variance:")
        print(data_component_variance)
        print_terminal_partition(level=3)
    # Compile information.
    information = {
        "data_observations_components": data_component,
        "data_components_variances": data_component_variance,
    }
    # Return information.
    return information


# Pairwise correlations in symmetric adjacency matrix


def collect_pairwise_correlations_pobabilities(
    method=None,
    features=None,
    data=None,
):
    """
    Calculates and collects correlation coefficients and probabilities across
    observations between pairs of features. These values are symmetrical and
    not specific to order in pairs of features.

    arguments:
        method (str): method for correlation, pearson, spearman, or kendall
        features (list<str>): names of features
        data (object): Pandas data frame of features (columns) across
            observations (rows)

    raises:

    returns:
        (dict): correlation coefficients and probabilities for pairs of
            features

    """

    # Organize data.
    data_copy = data.copy(deep=True)

    # Determine ordered, pairwise combinations of genes.
    # ABCD: AB AC AD BC BD CD
    pairs = combine_unique_elements_pairwise_orderless(
        elements=features,
    )

    # Collect counts and correlations across pairs of features.
    collection = dict()
    # Iterate on features.
    # Generate appropriate values for self-pairs.
    for feature in features:
        collection[feature] = dict()
        collection[feature][feature] = dict()
        collection[feature][feature]["count"] = float("nan")
        collection[feature][feature]["correlation"] = 1.0
        collection[feature][feature]["probability"] = 0.0
        collection[feature][feature]["discovery"] = float("nan")
        collection[feature][feature]["significance"] = True

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

        # Collect information.
        if not pair[0] in collection:
            collection[pair[0]] = dict()
            pass
        collection[pair[0]][pair[1]] = dict()
        collection[pair[0]][pair[1]]["count"] = count
        collection[pair[0]][pair[1]]["correlation"] = correlation
        collection[pair[0]][pair[1]]["probability"] = probability

        pass

    # Return information.
    return collection


def calculate_pairwise_correlations_discoveries(
    collection=None,
    threshold=None,
):
    """
    Calculates and collects Benjamini-Hochberg false discovery rates from
    probabilities of correlations across observations between pairs of
    features.

    arguments:
        collection (dict): correlation coefficients and probabilities for pairs
            of features
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries

    raises:

    returns:
        (dict): correlation coefficients, probabilities, and discoveries for
            pairs of features

    """

    # Extract probabilities.
    data_probability = extract_pairs_values_long(
        collection=collection,
        key="probability",
        name="probability",
    )

    # Calculate false discovery rates.
    data_discovery = calculate_false_discovery_rate(
        threshold=threshold,
        probability="probability",
        discovery="discovery",
        significance="significance",
        data_probabilities=data_probability,
    )

    # Introduce false discovery rates to collection.
    collection_discovery = copy.deepcopy(collection)
    records = convert_dataframe_to_records(data=data_discovery)
    for record in records:
        one = record["feature_one"]
        two = record["feature_two"]
        collection_discovery[one][two]["discovery"] = record["discovery"]
        collection_discovery[one][two]["significance"] = record["significance"]
        pass
    # Return information.
    return collection_discovery


def extract_pairs_values_long(
    collection=None,
    key=None,
    name=None,
):
    """
    Extracts  a matrix of correlation coefficients.

    arguments:
        collection (dict): collection of correlation coefficients and
            probabilities across observations between pairs of features
        key (str): key to extract values, either correlation or probability
        name (str): name of column for values

    raises:

    returns:
        (object): Pandas data frame of correlation coefficients

    """

    # Copy.
    collection = copy.deepcopy(collection)

    # Collect counts and correlations across pairs of features.
    records = list()
    # Iterate on dimension one of entities.
    for feature_one in collection:
        # Iterate on dimension two of entities.
        for feature_two in collection[feature_one]:
            # Collect information.
            record = dict()
            record["feature_one"] = feature_one
            record["feature_two"] = feature_two
            record[name] = collection[feature_one][feature_two][key]
            # Collect information.
            records.append(record)
            pass
        pass

    # Organize data.
    data = convert_records_to_dataframe(records=records)
    #data.set_index(
    #    ["feature_one", "feature_two"],
    #    append=False,
    #    drop=True,
    #    inplace=True
    #)
    # Return information.
    return data


def organize_symmetric_adjacency_matrix(
    features=None,
    collection=None,
    key=None,
    threshold=None,
    fill=None,
):
    """
    Organizes a symmetric adjacency matrix of correlation coefficients.

    arguments:
        features (list<str>): names of features
        collection (dict): collection of correlation coefficients and
            probabilities across observations between pairs of features
        key (str): key to extract values, either correlation or probability
        threshold (bool): whether to include only values that are significant
            by the false discovery rate
        fill (float): value with which to fill missing values or those that are
            insignificant

    raises:

    returns:
        (object): Pandas data frame of correlation coefficients

    """

    # Collect counts and correlations across pairs of features.
    records = list()
    # Iterate on dimension one of entities.
    for feature_one in features:
        # Collect information.
        record = dict()
        record["feature"] = feature_one

        # Iterate on dimension two of entities.
        for feature_two in features:
            # Access value.
            if feature_two in collection[feature_one]:
                value = collection[feature_one][feature_two][key]
                significance = (
                    collection[feature_one][feature_two]["significance"]
                )
            elif feature_one in collection[feature_two]:
                value = (collection[feature_two][feature_one][key])
                significance = (
                    collection[feature_two][feature_one]["significance"]
                )
            else:
                significance = False
            pass
            # Collect information.
            if threshold:
                if significance:
                    record[feature_two] = value
                else:
                    record[feature_two] = fill
            else:
                record[feature_two] = value
        # Collect information.
        records.append(record)
        pass

    # Organize data.
    data = convert_records_to_dataframe(records=records)
    data.sort_values(
        by=["feature"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    data.set_index(
        ["feature"],
        append=False,
        drop=True,
        inplace=True
    )
    data.rename_axis(
        columns="features",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Return information.
    return data


def cluster_adjacency_matrix(
    data=None,
):
    """
    Organizes a matrix of correlation coefficients.

    This function is currently specific to a symmetric adjacency matrix.

    arguments:
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    # Cluster.
    columns = data.columns.to_numpy()#.tolist()
    rows = data.index.to_numpy()#.tolist()
    matrix = numpy.transpose(data.to_numpy())
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "single", "complete", "average"
        metric="euclidean",
        optimal_ordering=True,
    )
    dendrogram = scipy.cluster.hierarchy.dendrogram(
        linkage,
    )
    dimension = len(matrix)
    # Access seriation from dendrogram leaves.
    leaves = dendrogram["leaves"]
    # Generate new matrix with values of 1.0 that will persist on the diagonal.
    #matrix_order = numpy.zeros((dimension, dimension))
    matrix_cluster = numpy.full(
        (dimension, dimension),
        fill_value=1.0,
    )
    # Sort matrix values.
    a,b = numpy.triu_indices(dimension, k=1)
    matrix_cluster[a,b] = (
        matrix[[leaves[i] for i in a], [leaves[j] for j in b]]
    )
    matrix_cluster[b,a] = matrix_cluster[a,b]
    # Sort matrix row and column labels.
    indices = range(0, dimension)
    rows_sort = list(map(
        lambda index: rows[leaves[index]],
        indices
    ))
    columns_sort = list(map(
        lambda index: columns[leaves[index]],
        indices
    ))
    # Organize data.
    data_cluster = pandas.DataFrame(
        data=matrix_cluster,
        index=rows_sort,
        columns=columns_sort,
    )
    # Return information.
    return data_cluster


def cluster_data_columns(
    data=None,
):
    """
    Clusters features on columns by their similarities across instances on
    rows.

    arguments:
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    # Cluster.
    columns = data.columns.to_numpy()#.tolist()
    rows = data.index.to_numpy()#.tolist()
    index_name = data.index.name
    # Plan to cluster across columns.
    # Organize columns across dimension zero.
    matrix = numpy.transpose(data.to_numpy())
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "single", "complete", "average"
        metric="euclidean",
        optimal_ordering=True,
    )
    dendrogram = scipy.cluster.hierarchy.dendrogram(
        linkage,
    )
    # Access seriation from dendrogram leaves.
    leaves = dendrogram["leaves"]
    # Sort matrix row and column labels.
    indices = range(0, len(matrix))
    matrix_cluster = list(map(
        lambda index: matrix[leaves[index]],
        indices
    ))
    columns_sort = list(map(
        lambda index: columns[leaves[index]],
        indices
    ))
    # Organize data.
    data_cluster = pandas.DataFrame(
        data=numpy.transpose(matrix_cluster),
        index=rows,
        columns=columns_sort,
    )
    data_cluster.rename_axis(
        index_name,
        axis="index",
        inplace=True,
    )
    # Return information.
    return data_cluster


def cluster_data_rows(
    data=None,
):
    """
    Clusters instances on rows by their similarities across features on
    columns.

    arguments:
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    # Cluster.
    columns = data.columns.to_numpy()#.tolist()
    rows = data.index.to_numpy()#.tolist()
    index_name = data.index.name
    # Plan to cluster across columns.
    # Organize rows across dimension zero.
    matrix = data.to_numpy()
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "single", "complete", "average"
        metric="euclidean",
        optimal_ordering=True,
    )
    dendrogram = scipy.cluster.hierarchy.dendrogram(
        linkage,
    )
    # Access seriation from dendrogram leaves.
    leaves = dendrogram["leaves"]
    # Sort matrix row and column labels.
    indices = range(0, len(matrix))
    matrix_cluster = list(map(
        lambda index: matrix[leaves[index]],
        indices
    ))
    rows_sort = list(map(
        lambda index: rows[leaves[index]],
        indices
    ))
    # Organize data.
    data_cluster = pandas.DataFrame(
        data=matrix_cluster,
        index=rows_sort,
        columns=columns,
    )
    data_cluster.rename_axis(
        index_name,
        axis="index",
        inplace=True,
    )
    # Return information.
    return data_cluster


def filter_features_by_threshold_outer_count(
    data=None,
    threshold_high=None,
    threshold_low=None,
    count=None,
):
    """
    Filters features in a symmetric adjacency matrix by whether a certain count
    of values are outside of low and high thresholds.

    arguments:
        data (object): Pandas data frame of values
        threshold_high (float): value must be greater than this threshold
        threshold_low (float): value must be less than this threshold
        count (int): minimal count of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of values

    """

    def count_true(slice=None, count=None):
        values = slice.to_numpy().tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Copy data.
    data = data.copy(deep=True)
    # Determine whether values exceed threshold.
    data_threshold = data.applymap(
        lambda value: ((value <= threshold_low) or (threshold_high <= value))
    )
    # This aggregation operation produces a series.
    data_count_row = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis="columns",
    )
    # Select rows with appropriate values.
    data_row = data.loc[data_count_row, : ]
    # Extract identifiers of features.
    features = collect_unique_elements(
        elements_original=data_row.index.tolist()
    )
    # Select identical features across rows and columns.
    data_feature_row = data.loc[
        data.index.isin(features),
        :
    ]
    data_feature_column = data_feature_row.loc[
        :, data_feature_row.columns.isin(features)
    ]
    # Return information.
    return data_feature_column


def organize_feature_signal_correlations(
    method=None,
    threshold_high=None,
    threshold_low=None,
    count=None,
    discovery=None,
    features=None,
    data_signal=None,
):
    """
    Calculates and organizes correlation coefficients between pairs of features
    across observations.

    arguments:
        method (str): method for correlation, pearson, spearman, or kendall
        threshold_high (float): value must be greater than this threshold
        threshold_low (float): value must be less than this threshold
        count (int): minimal count of rows or columns that must
            pass threshold
        discovery (float): value of false discovery rate for which to include
            correlation coefficients
        features (list<str>): identifiers of features
        data_signal (object): Pandas data frame of features' signals across
            observations

    raises:

    returns:
        (object): Pandas data frame of genes' pairwise correlations

    """

    # Organize data.
    data = data_signal.copy(deep=True)

    # Calculate correlations between gene pairs of their pantissue signals
    # across persons.
    # Only collect genes with correlations beyond thresholds.
    collection = collect_pairwise_correlations_pobabilities(
        method=method,
        features=features,
        data=data,
    )
    # Calculate false discovery rates for each primary feature.
    collection_discovery = calculate_pairwise_correlations_discoveries(
        collection=collection,
        threshold=discovery,
    )

    # Organize data matrix for correlation coefficients.
    # Filter correlation coefficients by value of false discovery rate.
    # Clustering algorithm requires that matrix not include missing values.
    data_correlation = organize_symmetric_adjacency_matrix(
        features=features,
        collection=collection_discovery,
        key="correlation",
        threshold=True,
        fill=0.0,
    )
    # Filter genes by their pairwise correlations.
    data_pass = filter_features_by_threshold_outer_count(
        data=data_correlation,
        threshold_high=threshold_high,
        threshold_low=threshold_low,
        count=count,
    )
    # Cluster data.
    data_cluster = cluster_adjacency_matrix(
        data=data_pass,
    )
    #data_cluster = data_pass
    # Return information.
    return data_cluster


# Metabolic information.


def prepare_curation_report(
    compartments=None,
    processes=None,
    reactions=None,
    metabolites=None
):
    """
    Prepares a summary report on curation of metabolic sets and entities.

    arguments:
        compartments (dict<dict>): information about compartments
        processes (dict<dict>): information about processes
        reactions (dict<dict>): information about reactions
        metabolites (dict<dict>): information about metabolites

    returns:
        (str): report of summary information

    raises:

    """

    # Count compartments.
    count_compartments = len(compartments)
    # Count processes.
    count_processes = len(processes)
    # Count reactions.
    count_reactions = len(reactions)
    # Count metabolites.
    count_metabolites = len(metabolites)
    # Count reactions with references to MetaNetX.
    count_one = count_entities_with_references(
        references=["metanetx"],
        entities=reactions
    )
    proportion_one = count_one / count_reactions
    percentage_one = round((proportion_one * 100), 2)
    # Count reactions with references either to genes or enzyme commission.
    count_two = count_entities_with_references(
        references=["gene", "enzyme"],
        entities=reactions
    )
    proportion_two = count_two / count_reactions
    percentage_two = round((proportion_two * 100), 2)
    # Count metabolites with references to MetaNetX.
    count_three = count_entities_with_references(
        references=["metanetx"],
        entities=metabolites
    )
    proportion_three = count_three / count_metabolites
    percentage_three = round((proportion_three * 100), 2)
    # Count metabolites with references to Human Metabolome Database (HMDB) and
    # PubChem.
    count_four = count_entities_with_references(
        references=["hmdb", "pubchem"],
        entities=metabolites
    )
    proportion_four = count_four / count_metabolites
    percentage_four = round((proportion_four * 100), 2)
    # Compile information.
    report = textwrap.dedent("""\

        --------------------------------------------------
        curation report

        compartments: {count_compartments}
        processes: {count_processes}
        reactions: {count_reactions}
        metabolites: {count_metabolites}

        reactions in MetaNetX: {count_one} ({percentage_one} %)
        reactions with gene or enzyme: {count_two} ({percentage_two} %)
        metabolites in MetaNetX: {count_three} ({percentage_three} %)
        metabolites with HMDB or PubChem: {count_four} ({percentage_four} %)

        --------------------------------------------------
    """).format(
        count_compartments=count_compartments,
        count_processes=count_processes,
        count_reactions=count_reactions,
        count_metabolites=count_metabolites,
        count_one=count_one,
        percentage_one=percentage_one,
        count_two=count_two,
        percentage_two=percentage_two,
        count_three=count_three,
        percentage_three=percentage_three,
        count_four=count_four,
        percentage_four=percentage_four
    )
    # Return information.
    return report


def count_entities_with_references(references=None, entities=None):
    """
    Counts entities with any of specific references.

    arguments:
        references (list<str>): identifiers of references
        entities (dict<dict>): information about entities

    returns:
        (int): count of entities with specific reference

    raises:

    """

    count = 0
    for entity in entities.values():
        matches = []
        for reference in references:
            if reference in entity["references"].keys():
                if len(entity["references"][reference]) > 0:
                    matches.append(True)
        if any(matches):
            count += 1
    return count


def collect_reaction_participants_value(
    key=None, criteria=None, participants=None
):
    """
    Collects a value from a reaction's specific participants

    arguments:
        key (str): key of value to collect from each participant
        criteria (dict<list>): criteria by which to select participants
        participants (list<dict>): information about a reaction's participants

    returns:
        (list<str>): values from a reaction's participants

    raises:

    """

    participants_match = filter_reaction_participants(
        criteria=criteria, participants=participants
    )
    return collect_value_from_records(key=key, records=participants_match)


def filter_reaction_participants(criteria=None, participants=None):
    """
    Filters a reaction's participants by multiple criteria

    arguments:
        criteria (dict<list>): criteria by which to select participants
        participants (list<dict>): information about a reaction's participants

    returns:
        (list<dict>): information about a reaction's participants

    raises:

    """

    def match(participant):
        if "metabolites" in criteria:
            match_metabolite = (
                participant["metabolite"] in criteria["metabolites"]
            )
        else:
            match_metabolite = True
        if "compartments" in criteria:
            match_compartment = (
                participant["compartment"] in criteria["compartments"]
            )
        else:
            match_compartment = True
        if "roles" in criteria:
            match_role = participant["role"] in criteria["roles"]
        else:
            match_role = True
        return match_metabolite and match_compartment and match_role
    return list(filter(match, participants))






###############################################################################
# Procedure
