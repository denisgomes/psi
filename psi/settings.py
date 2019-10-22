"""Default application settings.

When a model is first saved or when model options are saved, the application
settings are saved as model settings. Application settings are settings that
reside on the application side whereas model settings are transfered via the
saved file. The model settings, if one exists are merged with the application
settings when a model is opened. The application always uses the application
settings internally.
"""

from psi import VERSION


options = {
        ### Core Options ###

        # units systems, options: ("english", "si")
        "core.units": "english",

        # vertical direction, options: ("y", "z")
        "core.vertical": "y",

        # ambient temperature
        "core.tref": 70,

        # rigid element stiffness
        "core.rigid_stiffness": 1e12,

        # for corrosion allowance and mill tolerance, reduced thickness
        # is considered for stress calculations
        "core.stress_cases_corroded": False,

        # program version
        "core.version": VERSION,
        }

options_types = {}


def validate(settings):
    """Validate all settings"""

    pass


def merge(settings):
    """Merger model settings with application settings"""

    pass
