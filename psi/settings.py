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
        # Core Options #

        # units systems, options: ("english", "si")
        "core.units": "english",

        # vertical direction, options: ("y", "z")
        "core.vertical": "y",

        # ambient temperature
        "core.tref": 70,

        # for corrosion allowance and mill tolerance, reduced thickness
        # is considered for stress calculations
        "core.stress_cases_corroded": False,

        # activate bourdon pressure effect
        "code.bourdon_effect": False,

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
