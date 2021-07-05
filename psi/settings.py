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

"""Default application settings.

A configuration object is created when a model is first instantiated. The
program will use the model settings defined in the configuration instance.

Application and model settings may diverge over time due to different versions.
When the software is updated new application settings will be merged with the
model settings to ensure backwards compatibility between version. Note that
application settings in this context are the settings in this file with respect
to the latest version of the software.
"""

from psi import VERSION
from psi import units
from psi.units import DEFAULT_UNITS


@units.define(translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness",
              tref="temperature")
class Configuration:
    """Default model settings.

    Parameters
    ----------
    units : str
        Define the model units. 'english' by default.

    vertical : str
        Define the vertical direction for the model. 'y' by default.

    stress_case_corroded : bool
        Use reduced thickness to evaluate pipe stresses. Reduced pipe wall is
        not used for sections property calculations.

    bourdon_effect : bool
        Include bourdon effects.

    pressure_thrust : bool
        Include pressure thrust effects.

    liberal_stress : bool
        Use liberal stress allowable for expansion stresses.

    weak_springs : bool
        Include weak springs for numerical stability.

    translation_stiffness : float
        Support stiffness used in the translation directions.

    rotation_stiffness : float
        Support stiffness used in the rotation directions.

    axial_force : bool
        Include axial force due to structural loading for code stress
        calculations.

    tref : float
        Reference temperture used for thermal expansion calculations.

    timoshenko : bool
        Use Timoshenko beam formulation accounting for shear deformation.

    version : str
        Latest software version.
    """

    _app = None

    @property
    def app(self):
        return self._app

    def __init__(self, default_units=DEFAULT_UNITS):
        with units.Units(user_units=default_units):
            # core options
            self._units = default_units
            self.vertical = "y"
            self.stress_case_corroded = True
            self.bourdon_effect = False
            self.pressure_thrust = False
            self.liberal_stress = False
            self.weak_springs = False
            self.translation_stiffness = 10**10
            self.rotation_stiffness = 10**12
            self.axial_force = False
            self.tref = 70.0
            self.timoshenko = True
            self.version = VERSION

            # gui options
            self.debug = False
            self.background_top_color = (0, 0, 0)
            self.background_bottom_color = (0, 0, 0)
            self.pipe_color = (45, 45, 45)
            self.valve_color = (100, 100, 100)
            self.reducer_color = (30, 111, 121)
            self.flange_color = (33, 100, 244)
            self.anchor_color = (11, 23, 113)
            self.support_color = (111, 111, 111)
            self.spring_color = (133, 0, 123)
            self.snubber_color = (123, 111, 222)
            self.point_color = (10, 100, 10)
            self.line_color = (34, 123, 11)
            self.face_color = (22, 11, 11)

            self.epsilon = 0.01
            self.pick_color = (200, 200, 200)
            self.mouse_rotate_factor = 1.5
            self.mouse_translate_factor = 0.3

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, value):
        self._units = value
        self.app.units.set_user_units(value)

    def export(self, fname):
        """Export the current model settings"""
        pass

    def import_(self, fname):
        """Import settings from an outside source"""
        pass

    def merge(self, other):
        """Merge configuration settings from another configuration"""
        pass

    def validate(self, other):
        """Determine that the configuation types are consistant"""
        pass
