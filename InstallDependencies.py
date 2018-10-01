#!/usr/bin/env python
"""
Installs Necessary Modules To Get The Reduction Pipeline Going.
"""

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import shutil
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Checks For Dependencies Required For Pipeline
# ------------------------------------------------------------------------------------------------------------------- #
MissingLibraries = []

try:
    import numpy
except ImportError:
    numpy = None
    MissingLibraries.append("numpy")

try:
    import scipy
except ImportError:
    scipy = None
    MissingLibraries.append("scipy")

try:
    import pyraf
except ImportError:
    pyraf = None
    MissingLibraries.append("pyraf")

try:
    import astropy
except ImportError:
    astropy = None
    MissingLibraries.append("astropy")

try:
    import ephem
except ImportError:
    ephem = None
    MissingLibraries.append("ephem")

try:
    import matplotlib
except ImportError:
    matplotlib = None
    MissingLibraries.append("matplotlib")

try:
    import imreg_dft
except ImportError:
    imreg_dft = None
    MissingLibraries.append("imreg_dft")

try:
    import specutils
except ImportError:
    specutils = None
    MissingLibraries.append("specutils")


if len(MissingLibraries) != 0:
    print "#" + "-" * 40 + "#"
    print "ERROR: Following Dependencies Not Found"
    print "#" + "-" * 40 + "#"

    for index, library in enumerate(MissingLibraries):
        print "{0} : {1}".format(index + 1, library)

    print "#" + "-" * 40 + "#"
    sys.exit(2)
else:
    print("Success: All Dependencies Satisfied")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Copying The Script In The Directory Specified
# ------------------------------------------------------------------------------------------------------------------- #
DIR_HOME = os.getenv('HOME')
DIR_CURNT = os.getcwd()
DIR_INSTALL = input("Enter Installation Path: ").strip()

if os.path.isdir(DIR_INSTALL):
    query_delete = input("Installation Directory Already Exists. Replace It? (Y/N)") or "N"
    if query_delete.upper() == 'Y':
        try:
            shutil.rmtree(DIR_INSTALL, ignore_errors=True)
        except:
            print "Could Not Delete The Installation Directory Specified. Check Permissions"

print("#-----Copying Files To The Installation Directory-----#")
shutil.copytree(DIR_CURNT + "/BasicReduction", DIR_INSTALL + "/")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Make Modules Executable
# ------------------------------------------------------------------------------------------------------------------- #
ListModules = ['Photometry', 'Spectroscopy', 'ObsPlan', 'PhotPreProcess', 'FluxCalib', 'SpecPreProcess', 'Align']

print("#-----Converting Scripts To Executables-----#")

for module in ListModules:
    PathModule = DIR_INSTALL + "/" + str(module)
    os.chmod(PathModule, 0744)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Making A Backup Of IRAF And Performing Necessary Changes
# ------------------------------------------------------------------------------------------------------------------- #
print("#-----Backing Up Current IRAF Directory-----#")

if os.path.exists(DIR_HOME + "/iraf"):
    shutil.copytree(DIR_HOME + "/iraf", DIR_HOME + "/iraf.bak")
    shutil.rmtree(DIR_HOME + "/iraf")
    print("Current IRAF Directory Backed Up Into '{0}'".format(DIR_HOME + "/iraf.bak"))
    shutil.copytree(DIR_INSTALL + "/iraf", DIR_HOME)
else:
    shutil.copytree(DIR_INSTALL + "/iraf", DIR_HOME + "/iraf")

print("#-----Running 'mkiraf' In The New Installation Directory-----#")

os.chdir(DIR_HOME + "/iraf")
os.system("mkiraf")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Modifying '.bashrc' File To Specify Installation Directory To PATH & Specifying IRAF Home To File 'IRAFHome'
# ------------------------------------------------------------------------------------------------------------------- #
print("#-----Specifying Installation Directory To The PATH Environment Variable-----#")
modify_bash = raw_input("Modify The '.bashrc' To Reflect The Change? (Y/N): ").upper()
if modify_bash == "Y":
    with open(DIR_HOME + "/.bashrc", 'a') as bash:
        IRAF_PATH = "\nexport PATH=$PATH:{0}".format(DIR_INSTALL)
        bash.write(IRAF_PATH)
    os.system("bash")
else:
    print """
    Change PATH Variable In '.bashrc' To Point To The Installation Directory:
    Add 'export PATH=$PATH:{0}'""".format(DIR_INSTALL)

print("#-----Providing Location Of Home IRAF Directory To The File IRAFHome-----#")
with open(DIR_INSTALL + "/IRAFHome.py", 'w') as firaf:
    firaf.write("irafhome = '%s'".format(DIR_INSTALL + "/iraf"))
# ------------------------------------------------------------------------------------------------------------------- #
