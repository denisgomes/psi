"""Project level settings"""

import os

NAME = "PSI"

VERSION = "0.0.1"

PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT_DIRECTORY = os.path.abspath(os.path.join(PACKAGE_ROOT_DIRECTORY,
                                                      ".."))
LICENSE = open(os.path.join(PROJECT_ROOT_DIRECTORY, "LICENSE")).read()

README = open(os.path.join(PROJECT_ROOT_DIRECTORY, "README.rst")).read()

COPYRIGHT = "Copyright (c) 2019 Denis Gomes"

DESCRIPTION = "An engineering pipe stress design and analysis software."

DATA_DIRECTORY = os.path.join(PACKAGE_ROOT_DIRECTORY, "data")
BEAM_DATA_FILE = os.path.join(DATA_DIRECTORY, "wshapes.csv")
PIPE_DATA_FILE = os.path.join(DATA_DIRECTORY, "pipes.csv")
MATERIAL_DATA_FILE = os.path.join(DATA_DIRECTORY, "materials.csv")
INSULATION_DATA_FILE = os.path.join(DATA_DIRECTORY, "insulation.csv")
UNITS_DIRECTORY = os.path.join(DATA_DIRECTORY, "units")
TEMPLATE_DIRECTORY = os.path.join(PACKAGE_ROOT_DIRECTORY, "templates")
