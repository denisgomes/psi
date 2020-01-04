# Copyright (c) 2019 Denis Gomes
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
#    * Neither the name of Pipe Stress Infinity (PSI) nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

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

DESCRIPTION = "The pipe stress design and analysis software."

DATA_DIRECTORY = os.path.join(PACKAGE_ROOT_DIRECTORY, "data")
BEAM_DATA_FILE = os.path.join(DATA_DIRECTORY, "wshapes.csv")
PIPE_DATA_FILE = os.path.join(DATA_DIRECTORY, "pipes.csv")
MATERIAL_DATA_FILE = os.path.join(DATA_DIRECTORY, "materials.csv")
INSULATION_DATA_FILE = os.path.join(DATA_DIRECTORY, "insulation.csv")
UNITS_DIRECTORY = os.path.join(DATA_DIRECTORY, "units")
TEMPLATE_DIRECTORY = os.path.join(PACKAGE_ROOT_DIRECTORY, "templates")