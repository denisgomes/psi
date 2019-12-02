# Copyright (c) 2019 Denis Gomes
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.

#    * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.

#    * Neither the name of Pipe Stress Infinity (PSI) nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.

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
