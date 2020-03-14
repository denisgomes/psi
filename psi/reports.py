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

"""Display various result reports generated from the loadcases solved."""

from datetime import datetime
import sys
from contextlib import redirect_stdout

from tqdm import tqdm
from jinja2 import Environment, FileSystemLoader

from psi import TEMPLATE_DIRECTORY
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
        """Create a movement report instance.

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
        version = self.app.models.active_object.settings.version
        date = datetime.now().date()
        jobname = self.app.models.active_object.jobname
        time = datetime.strftime(datetime.now(), "%I:%M %p")

        with redirect_stdout(sys.__stdout__):
            tqdm.write(self.template.render(version=version,
                                        date=date,
                                        time=time,
                                        jobname=jobname,
                                        licensed_to="Community",
                                        report_type=self.__class__.__name__,
                                        report_desc="Displacement Report",
                                        units=Quantity.user_units,
                                        loadcases=self.loadcases,
                                        zip=zip,     # pass zip
                                        enumerate=enumerate,
                                        abs=abs,
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
        version = self.app.models.active_object.settings.version
        date = datetime.now().date()
        jobname = self.app.models.active_object.jobname
        time = datetime.strftime(datetime.now(), "%I:%M %p")

        with redirect_stdout(sys.__stdout__):
            tqdm.write(self.template.render(version=version,
                                        date=date,
                                        time=time,
                                        jobname=jobname,
                                        licensed_to="Community",
                                        report_type=self.__class__.__name__,
                                        report_desc="Reactions Report",
                                        units=Quantity.user_units,
                                        loadcases=self.loadcases,
                                        zip=zip,     # pass zip
                                        enumerate=enumerate,
                                        abs=abs,
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
        # version = options["core.version"]
        version = self.app.models.active_object.settings.version
        date = datetime.now().date()
        jobname = self.app.models.active_object.jobname
        time = datetime.strftime(datetime.now(), "%I:%M %p")

        with redirect_stdout(sys.__stdout__):
            tqdm.write(self.template.render(version=version,
                                        date=date,
                                        time=time,
                                        jobname=jobname,
                                        licensed_to="Community",
                                        report_type=self.__class__.__name__,
                                        report_desc="Forces Report",
                                        units=Quantity.user_units,
                                        loadcases=self.loadcases,
                                        zip=zip,     # pass zip
                                        enumerate=enumerate,
                                        abs=abs,
                                        ))


class Stresses(Report):
    """Code stress output results"""

    def __init__(self, name, loadcases):
        """Create a forces report instance.

        Parameters
        ----------
        name : str
            Unique name for report object.

        loadcases : list of loadcases
            Loadcases for which results are displayed.
        """
        super(Stresses, self).__init__(name, loadcases)

        if len(loadcases) == 1:
            self.template = self.env.get_template("single_case_stresses")
        else:
            self.template = self.env.get_template("multiple_case_stresses")

    def to_screen(self):
        """Print forces report results to screen."""
        # version = options["core.version"]
        version = self.app.models.active_object.settings.version
        date = datetime.now().date()
        jobname = self.app.models.active_object.jobname
        time = datetime.strftime(datetime.now(), "%I:%M %p")

        with redirect_stdout(sys.__stdout__):
            tqdm.write(self.template.render(version=version,
                                        date=date,
                                        time=time,
                                        jobname=jobname,
                                        licensed_to="Community",
                                        report_type=self.__class__.__name__,
                                        report_desc="Stresses Report",
                                        units=Quantity.user_units,
                                        loadcases=self.loadcases,
                                        zip=zip,     # pass zip
                                        enumerate=enumerate,
                                        abs=abs,
                                        ))


class ReportContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(ReportContainer, self).__init__()
        self.Movements = Movements
        self.Reactions = Reactions
        self.Forces = Forces
        self.Stresses = Stresses
