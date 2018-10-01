#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxx-------------------MODULE FOR FILE HANDLING-----------------xxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import sys
import glob
import shutil
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling
# ------------------------------------------------------------------------------------------------------------------- #

def remove_file(file_name):
    """
    Removes the file "file_name" in the constituent directory.
    Args:
         file_name  : Name of the file to be removed from the current directory
    Returns:
        None
    """
    try:
        os.remove(file_name)
    except OSError:
        pass


def remove_similar_files(common_text):
    """
    Removes similar files based on the string "common_text".
    Args:
        common_text : String containing partial name of the files to be deleted
    Returns:
        None
    """
    for residual_file in glob.glob(common_text):
        remove_file(residual_file)


def group_similar_files(text_list, common_text, exceptions=''):
    """
    Groups similar files based on the string "common_text". Writes the similar files
    onto the list 'text_list' (only if this string is not empty) and appends the similar
    files to a list 'python_list'.
    Args:
        text_list   : Name of the output text file with names grouped based on the 'common_text'
        common_text : String containing partial name of the files to be grouped
        exceptions  : String containing the partial name of the files that need to be excluded
    Returns:
        list_files  : Python list containing the names of the grouped files
    """
    list_files = glob.glob(common_text)
    if exceptions != '':
        list_exception = exceptions.split(',')
        for file_name in glob.glob(common_text):
            for text in list_exception:
                bool_test = re.search(str(text), file_name)
                if bool_test:
                    try:
                        list_files.remove(file_name)
                    except ValueError:
                        pass

    list_files.sort()
    if len(text_list) != 0:
        with open(str(text_list), "w") as f:
            for index in range(0, len(list_files)):
                f.write(str(list_files[index]) + "\n")

    return list_files


def copy_files(in_path, out_path, common_text, exceptions=''):
    """
    Copies similar files based on the string "common_text" from the directory specified by "in_path"
    onto the directory specified by "out_path".
    Args:
        in_path     : Path of the directory from which files are to be copied
        out_path    : Path of the directory to which files are to be copied
        common_text : String containing partial name of the files to be copied
        exceptions  : String containing the partial name of the files that need to be excluded
    Returns:
        None
    """
    owd = os.getcwd()
    os.chdir(str(in_path))

    list_copy = group_similar_files("", common_text=str(common_text), exceptions=str(exceptions))
    for file_name in list_copy:
        shutil.copy(os.path.join(str(in_path), str(file_name)), str(out_path))

    os.chdir(owd)


def text_list_to_python_list(text_list):
    """
    Returns data in the file 'text_list' as a python_list.
    Args:
        text_list   : Input file containing filenames
    Returns:
        python_list : List of all the elements in the file 'text_list'
    Raises:
        Error : File 'text_list 'Not Found
    """
    if os.path.isfile(text_list):
        with open(text_list, "r+") as f:
            python_list = f.read().split()
            return python_list
    else:
        print "Error : File '{0}' Not Found".format(str(text_list))
        sys.exit(1)


def python_list_to_text_list(python_list, text_list):
    """
    Put the data from the input 'python_list' to a file 'text_list' line-wise.
    Args:
        python_list : Python_list from which data has to be read
        text_list   : Name of the text file onto which data has to be appended
    Returns:
        None
    """
    with open(str(text_list), "w") as f:
        for element in python_list:
            f.write(str(element) + "\n")


def list_lists_to_list(list_lists, text_list):
    """
    Groups filenames from a list 'list_lists' onto a single file 'text_list'.
    Args:
        list_lists  : List containing the names of different lists
        text_list   : Name of the file onto which all the filenames from the 'list_lists' have to be appended
    Returns:
        list_name   : Python list containing the names of all the constituent files
    """
    list_name = []
    for file_name in list_lists:
        with open(str(file_name), 'r') as f:
            file_list = f.read().split()
            for elements in file_list:
                list_name.append(str(elements))
    python_list_to_text_list(list_name, str(text_list))

    return list_name


def check_ifexists(path):
    """
    Checks if a file or directory exists.
    Args:
        path        : Path whose existence is to be checked
    Returns:
        True        : Returns True only when the path exists
    """
    if path[-1] == '/':
        if os.path.isdir(str(path)):
            return True
            pass
        else:
            print "Error : Directory '{0}' Does Not Exist".format(str(path))

    else:
        if os.path.isfile(str(path)):
            return True
            pass
        else:
            print "Error : Path '{0}' Does Not Exist".format(str(path))


def display_text(text):
    """
    Displays text mentioned in the string 'text_to_display'
    Args:
        text    : Text to be displayed
    Returns:
        None
    """
    print "\n{0} {1} {0}".format("#", "-" * (12 + len(text)))
    print "{0} {1} {2} {1} {0}".format("#", "-" * 5, str(text))
    print "{0} {1} {0}\n".format("#", "-" * (12 + len(text)))

# ------------------------------------------------------------------------------------------------------------------- #
