"""Default application settings.

When a model is first saved or when model options are saved, the application
settings are saved as model settings. Application settings are settings that
reside on the application side whereas model settings are transfered via the
saved file. The model settings, if one exists are merged with the application
settings when a model is opened. The application always uses the application
settings internally.
"""

from openpipe import VERSION


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

        ### Graphics Options ###
        "graphics.debug": False,
        "graphics.culling": True,
        "graphics.background_top_color": (255, 255, 255),
        "graphics.background_bottom_color": (255, 255, 255),
        "graphics.pick_color": (153, 153, 255),
        "graphics.face_color": (108, 126, 148),
        "graphics.line_color": (0, 0, 255),
        "graphics.point_color": (0, 0, 0),
        "graphics.pipe_color": (108, 126, 148),
        "graphics.reducer_color": (108, 126, 148),
        "graphics.bend_color": (108, 126, 148),
        "graphics.valve_color": (29, 229, 30),
        "graphics.flange_color": (200, 123, 0),
        "graphics.anchor_color": (29, 100, 30),
        "graphics.support_color": (29, 200, 30),
        "graphics.spring_color": (23, 88, 234),
        "graphics.phong": True,
        "graphics.edge": True,
        "graphics.epsilon": 0.1,

        ### GUI Options ###
        "gui.mouse_translate_factor": 0.05,
        "gui.mouse_rotate_factor": 0.75,
        "gui.mouse_coarse_zoom_factor": 8,
        "gui.mouse_fine_zoom_factor": 0.25,
        }


options_types = {}


def validate(settings):
    """Validate all settings"""

    pass


def merge(settings):
    """Merger model settings with application settings"""

    pass
