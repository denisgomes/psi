# Copyright (c) 2019 Denis Gomes
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
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

"""Display various result reports generated from the loadcases solved."""

from datetime import datetime
import sys
from contextlib import redirect_stdout

from jinja2 import Environment, FileSystemLoader

from psi import TEMPLATE_DIRECTORY
from psi.settings import options
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi.units import Quantity


class Report(Entity, ActiveEntityMixin):

    def __init__(self, name, loadcases):
        super(Report, self).__init__(name)
        self.loadcases = loadcases
        self.env = Environment(loader=FileSystemLoader(TEMPLATE_DIRECTORY))

    @property
    def parent(self):
        return self.app.reports

    def to_screen(self):
        raise NotImplementedError("implement")

    def to_excel(self):
        raise NotImplementedError("implement")

    def to_word(self):
        raise NotImplementedError("implement")


class Movements(Report):
    """Nodal displacement results."""

    def __init__(self, name, loadcases):
        """Create a movements report instance.

        Parameters
        ----------
        name : str
            Unique name for report object.

        loadcases : list of loadcases
            Loadcases for which results are displayed.
        """
        super(Movements, self).__init__(name, loadcases)

        if len(loadcases) == 1:
            self.template = self.env.get_template("single_case_movements")
        else:
            self.template = self.env.get_template("multiple_case_movements")

    def to_screen(self):
        """Print movement report results to screen."""
        version = options["core.version"]
        date = datetime.now().date()
        jobname = self.app.models.active_object.jobname
        time = datetime.strftime(datetime.now(), "%I:%M %p")

        with redirect_stdout(sys.__stdout__):
            print(self.template.render(version=version,
                                       date=date,
                                       time=time,
                                       jobname=jobname,
                                       licensed_to="Humanity",
                                       report_type=self.__class__.__name__,
                                       report_desc="Displacement Report",
                                       units=Quantity.user_units,
                                       loadcases=self.loadcases,
                                       zip=zip,     # pass zip
                                       ))


class Reactions(Report):
    """Support reaction results."""

    def __init__(self, name, loadcases):
        """Create a reactions report instance.

        Parameters
        ----------
        name : str
            Unique name for report object.

        loadcases : list of loadcases
            Loadcases for which results are displayed.
        """
        super(Reactions, self).__init__(name, loadcases)

        if len(loadcases) == 1:
            self.template = self.env.get_template("single_case_reactions")
        else:
            self.template = self.env.get_template("multiple_case_reactions")

    def to_screen(self):
        """Print reaction report results to screen."""
        version = options["core.version"]
        date = datetime.now().date()
        jobname = self.app.models.active_object.jobname
        time = datetime.strftime(datetime.now(), "%I:%M %p")

        with redirect_stdout(sys.__stdout__):
            print(self.template.render(version=version,
                                       date=date,
                                       time=time,
                                       jobname=jobname,
                                       licensed_to="Humanity",
                                       report_type=self.__class__.__name__,
                                       report_desc="Support Reaction Report",
                                       units=Quantity.user_units,
                                       loadcases=self.loadcases,
                                       zip=zip,     # pass zip
                                       ))


class Forces(Report):
    """Internal forces results."""

    def __init__(self, name, loadcases):
        """Create a forces report instance.

        Parameters
        ----------
        name : str
            Unique name for report object.

        loadcases : list of loadcases
            Loadcases for which results are displayed.
        """
        super(Forces, self).__init__(name, loadcases)

        if len(loadcases) == 1:
            self.template = self.env.get_template("single_case_forces")
        else:
            self.template = self.env.get_template("multiple_case_forces")

    def to_screen(self):
        """Print forces report results to screen."""
        version = options["core.version"]
        date = datetime.now().date()
        jobname = self.app.models.active_object.jobname
        time = datetime.strftime(datetime.now(), "%I:%M %p")

        with redirect_stdout(sys.__stdout__):
            print(self.template.render(version=version,
                                       date=date,
                                       time=time,
                                       jobname=jobname,
                                       licensed_to="Humanity",
                                       report_type=self.__class__.__name__,
                                       report_desc="Internal Forces Report",
                                       units=Quantity.user_units,
                                       loadcases=self.loadcases,
                                       zip=zip,     # pass zip
                                       ))


class Stresses(Report):
    """Code stress output results"""
    pass


class ReportContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(ReportContainer, self).__init__()
        self.Movements = Movements
        self.Reactions = Reactions
        self.Forces = Forces
        self.Stresses = Stresses
