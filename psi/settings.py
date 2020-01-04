# Pipe Stress Infinity (PSI) - The pipe stress design and analysis software.
# Copyright (c) 2019 Denis Gomes

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

"""Default application settings.

When a model is first saved or when model options are saved, the application
settings are saved as model settings. Application settings are settings that
reside on the application side whereas model settings are transfered via the
saved file. The model settings, if one exists are merged with the application
settings when a model is opened. The model settings are then used for the
remainder of the session.

Applications settings and model settings may diverge over the course of time
due to different versions. When the software is update new application settings
will be merged onto the model settings to ensure, backwards compatibility
between version.
"""

from psi import VERSION


class Options:
    # Options should be encapsulated to associate units and default values
    pass


options = {
    # units systems, options: ("english", "si")
    "core.units": "english",

    # vertical direction, options: ("y", "z")
    "core.vertical": "y",

    # for corrosion allowance and mill tolerance, reduced thickness
    # is considered for stress calculations
    "core.stress_cases_corroded": False,

    # activate bourdon pressure effect
    "core.bourdon_effect": False,

    # include pressure thrust
    "core.pressure_thrust": False,

    # program version
    "core.version": VERSION,
}

options_types = {}


def validate(settings):
    """Validate all settings"""
    pass


def merge(settings):
    """Merge application settings with model settings"""
    pass
