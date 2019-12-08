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

"""Display various result reports generated from the loadcases solved"""

import datetime
import sys
from contextlib import redirect_stdout

from jinja2 import Environment, FileSystemLoader

from psi import TEMPLATE_DIRECTORY
from psi.settings import options
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi.units import Quantity


class Result(Entity, ActiveEntityMixin):

    def __init__(self, name, loadcases):
        super(Result, self).__init__(name)
        self.loadcases = loadcases
        self.env = Environment(loader=FileSystemLoader(TEMPLATE_DIRECTORY))

    @property
    def parent(self):
        return self.app.results

    def to_screen(self):
        raise NotImplementedError("implement")

    def to_excel(self):
        raise NotImplementedError("implement")

    def to_word(self):
        raise NotImplementedError("implement")


class Stresses(Result):
    """Code stress output results"""
    pass


class Movements(Result):
    """Nodal displacement results"""

    def __init__(self, name, loadcases):
        super(Movements, self).__init__(name, loadcases)

        if len(loadcases) == 1:
            self.template = self.env.get_template("single_case_displacement")
        else:
            self.template = self.env.get_template("multiple_case_displacement")

    def to_screen(self):
        version = options["core.version"]
        date = datetime.datetime.now().date()
        jobname = self.app.models.active_object.jobname
        time = datetime.datetime.now().time()

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


class Reactions(Result):
    """Support reaction results"""
    pass


class ResultContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(ResultContainer, self).__init__()
        self.Stresses = Stresses
        self.Movements = Movements
        self.Reactions = Reactions
