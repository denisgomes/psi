# Pipe Stress Infinity (PSI) - The pipe stress analysis and design software.
# Copyright (c) 2021 Denis Gomes <denisgomes@consultant.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Project level settings"""

import os

NAME = "PSI"
VERSION = "0.0.1"
COPYRIGHT = "Copyright (c) 2021 Denis Gomes"
DESCRIPTION = "The pipe stress design and analysis software."
PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
DATA_DIRECTORY = os.path.join(PACKAGE_ROOT_DIRECTORY, "data")
BEAM_DATA_FILE = os.path.join(DATA_DIRECTORY, "wshapes.csv")
PIPE_DATA_FILE = os.path.join(DATA_DIRECTORY, "pipes.csv")
MATERIAL_DATA_FILE = os.path.join(DATA_DIRECTORY, "materials.csv")
INSULATION_DATA_FILE = os.path.join(DATA_DIRECTORY, "insulation.csv")
FLUID_DATA_FILE = os.path.join(DATA_DIRECTORY, "fluids.csv")
UNITS_DIRECTORY = os.path.join(DATA_DIRECTORY, "units")
TEMPLATE_DIRECTORY = os.path.join(PACKAGE_ROOT_DIRECTORY, "templates")
SPRING_DIRECTORY = os.path.join(DATA_DIRECTORY, "springs")
