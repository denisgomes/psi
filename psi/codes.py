"""Implementation of different piping codes"""


from openpipe.core.entity import (Entity, EntityContainer, ActiveEntityMixin,
                                  ActiveEntityContainerMixin)


class Code(Entity, ActiveEntityMixin):

    def __init__(self, name):
        super(Code, self).__init__(name)

        self.activate()

    @property
    def parent(self):
        return self.app.codes

    def activate(self):
        self.parent.activate(self)

    def apply(self, elements=None):
        pass


class B311(Code):

    def __init__(self, name):
        super(B311, self).__init__(name)

    def thke(self):
        """Effective thickness used for stress calculations"""
        pass

    def Shoop(self, section, pressure):
        """Hoop stress due to pressure"""
        pass

    def Slp(self, section, pressure):
        """Longitudinal stress due to pressure.

        Exact formulation for the pressure stress per code section 102.3.2.
        The pressure stress is a primary stress and by definition can result
        in large gross deformations and overall failure, i.e., pipe rupture.
        """
        do = section.od
        di = section.id
        p = pressure.pres   # the pressure load

        return p*di**2 / (do**2-di**2)

    def Slb(self, section):
        pass

    def Sh(self, material, loadcase):
        """Material hot allowable"""

        # extract temperature of loadcase
        exp = loadcase.loads
        temp = exp.temp

        if loadcase.stype == "SUS":
            return material.sh[temp]
        elif loadcase.stype == "OCC":
            # k factor is 1.15 for events lasting 8hrs or less and
            # 1.20 for 1hr or less per code para. 102.3.3, use 1.15 to
            # be conservative
            return 1.15 * material.sh[temp]
        elif loadcase.stype == "EXP":
            pass

    def Sa(self):
        """Thermal stress allowable using left over sustained stress.

        Allowed only if the user selects to use liberal stress allowables.
        """
        pass


class CodeContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(CodeContainer, self).__init__()
        self.B311 = B311

    def apply(self, elements=None):
        pass
