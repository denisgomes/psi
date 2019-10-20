"""Implementation of code based stress intensification factors.

B31.1
-----
Stress intensification and flexibility factors are based on Mandatory Appendix
D. Per code, the calculated SIFs must always be greater than or eqaul to 1.0.

SIFs provided are valid for a D/t ratio of less than 100. Beyond this limit,
the pipe behaves like large duct piping and must be modeling using shell
elements.
"""

from __future__ import division


class Fitting(object):
    """A piping component such as a reducer or an elbow.

    The same sif is applied to the end points of the element. Note that the B31
    codes do not specify how to handle a situation where a point has multiple
    sifs defined.

    The fitting sif is baked into the respective element and has a flexibility
    factor.
    """

    def __init__(self, element):
        self.element = element

    def kfac(self):
        """Flexibility factor of fitting"""
        return 1.0

    def sif(self):
        """In and out-of-plane stress intensification factor"""
        sifi = sifo = 1.0

        return (sifi, sifo)


class Reducer(Fitting):
    """For concentric reducers"""

    def __init__(self, element):
        super(Reducer, self).__init__(element)

    def sif(self):
        assert self.element.type == "Bend", "element must be a bend"

        code = self.element.code
        section = self.element.section

        if isinstance(code, "B311"):
            pass


class Bend(Fitting):

    def __init__(self, element):
        super(Bend, self).__init__(element)

    def sif(self):
        if isinstance(self.element.code, "B311"):
            pass


class Point(object):
    """A tee type intersection or a welding connection"""

    def __init__(self, element, point):
        self.element = element
        self.point = point

    def is_intersection(self):
        """A tee type intersection"""
        return self.point.vertex.edges == 3

    def is_connection(self):
        """A welding connection"""
        return self.point.vertex.edges == 2

    def sif(self):
        """In and out-of-plane stress intensification factor"""
        sifi = sifo = 1.0

        return (sifi, sifo)


class Intersection(Point):
    """A tee type intersection"""

    def __init__(self, element, point):
        super(Intersection, self).__init__(element, point)

        assert self.intersection(), "invalid intersection point"


class Connection(Point):
    """A joint between two pipe spools"""

    def __init__(self, element, point):
        super(Connection, self).__init__(element, point)

        assert self.connection(), "invalid connection point"


class Welding(Intersection):
    """Welding tees per ASME B16.9"""

    def __init__(self, element, point):
        super(Welding, self).__init__(element, point)

    def sif(self, code, do, tn, dob, rx, tc):
        """
        Parameters
        ----------
        code : str
            The code used to calculate the parameters

        do : float
            Outer diameter of run (i.e. header) pipe

        tn : float
            Nomimal thickness of run (i.e. header) pipe

        dob : float
            Outer diameter of branch pipe

        rx : float
            External crotch radius of welding tees

        tc : float
            Crotch thickness
        """
        if code == "B31.1":
            # per mandatory appendix D
            r = (do-tn) / 2

            if rx >= dob/8 and tc >= 1.5*tn:
                h = 4.4*tn / r
            else:
                h = 3.1*tn / r

            # in-plane and out-of-plane sifs are the same
            # for B31.1 the higher is used
            sif = 0.9 / h**(2/3)

            # must be larger than or equal to 1.0
            if sif < 1.0:
                sif = 1.0

            sifi = sifo = sif

            return (sifi, sifo)


class Unreinforced(Intersection):
    """Unreinforced fabricated tee"""

    def __init__(self, element, point):
        super(Unreinforced, self).__init__(element, point)

    def sif(self, code, do, tn):
        """
        Parameters
        ----------
        code : str
            The code used to calculate the parameters

        do : float
            Outer diameter of run (i.e. header) pipe

        tn : float
            Nomimal thickness of run (i.e. header) pipe
        """
        if code == "B31.1":

            r = (do-tn) / 2

            h = tn / r

            # in-plane and out-of-plane sifs are the same
            # for B31.1 the higher is used
            sif = 0.9 / h**(2/3)

            if sif < 1.0:
                sif = 1.0

            sifi = sifo = sif

            return (sifi, sifo)


class Reinforced(Intersection):
    """Reinforced fabricated tee"""

    def __init__(self, element, point):
        super(Reinforced, self).__init__(element, point)

    def sif(self, code, do, tn, tr):
        """
        Parameters
        ----------
        code : str
            The code used to calculate the parameters

        do : float
            Outer diameter of run (i.e. header) pipe

        tn : float
            Nomimal thickness of run (i.e. header) pipe

        tr : float
            Pad thickness
        """
        if code == "B31.1":

            r = (do-tn) / 2

            if tr > 1.5*tn:
                h = 4.05*tn / r
            else:
                h = (tn+tr/2)**(5/2) / (r*tn**(3/2))

            # in-plane and out-of-plane sifs are the same
            # for B31.1 the higher is used
            sif = 0.9 / h**(2/3)

            if sif < 1.0:
                sif = 1.0

            sifi = sifo = sif

            return (sifi, sifo)


class Weldolet(Intersection):
    """Olet fitting with welded outlet branch"""

    def __init__(self, element, point):
        super(Weldolet, self).__init__(element, point)

    def sif(self, code, do, tn):
        """
        Parameters
        ----------
        code : str
            The code used to calculate the parameters

        do : float
            Outer diameter of run (i.e. header) pipe

        tn : float
            Nomimal thickness of run (i.e. header) pipe
        """
        if code == "B31.1":

            r = (do-tn) / 2

            h = 3.3*tn / r

            # in-plane and out-of-plane sifs are the same
            # for B31.1 the higher is used
            sif = 0.9 / h**(2/3)

            if sif < 1.0:
                sif = 1.0

            sifi = sifo = sif

            return (sifi, sifo)


class Sockolet(Intersection):
    """Olet fitting with socket welded outlet branch.

    A sockolet is similar to a weldolet with the exception that the branch
    pipe goes into a socket type connection and is welded. A sockolet can be
    modeled by defining a weldolot 'intersection' sif at the header pipe and
    a socket 'connection' sif at the branch connection.

    Note that ASME codes do not describe any way to combine sifs and so in
    this method is way to get around that. Others just stick with the higher
    of the two sifs when their is no other choice since adding them together
    is too conservative.
    """

    def __init__(self, element, point):
        super(Sockolet, self).__init__(element, point)

    def sif(self, code, do, tn):
        """
        Parameters
        ----------
        code : str
            The code used to calculate the parameters

        do : float
            Outer diameter of run (i.e. header) pipe

        tn : float
            Nomimal thickness of run (i.e. header) pipe
        """
        if code == "B31.1":

            r = (do-tn) / 2

            h = 3.3*tn / r

            # in-plane and out-of-plane sifs are the same
            # for B31.1 the higher is used
            sif = 0.9 / h**(2/3)

            if sif < 1.0:
                sif = 1.0

            sifi = sifo = sif

            return (sifi, sifo)


class Sweepolet(Intersection):
    """Contoured integrally reinforced insert with butt-welded branch"""

    def __init__(self, element, point):
        super(Sweepolet, self).__init__(element, point)

    def sif(self, code, do, tn, dob, rx, tc):
        """
        Parameters
        ----------
        code : str
            The code used to calculate the parameters

        do : float
            Outer diameter of run (i.e. header) pipe

        tn : float
            Nomimal thickness of run (i.e. header) pipe

        dob : float
            Outer diameter of branch pipe

        rx : float
            External crotch radius of welded-in contour insert

        tc : float
            Crotch thickness
        """
        if code == "B31.1":

            r = (do-tn) / 2

            if rx >= dob/8 and tc >= 1.5*tn:
                h = 4.4*tn / r
            else:
                h = 3.1*tn / r

            # in-plane and out-of-plane sifs are the same
            # for B31.1 the higher is used
            sif = 0.9 / h**(2/3)

            if sif < 1.0:
                sif = 1.0

            sifi = sifo = sif

            return (sifi, sifo)


class ButtWeld(Connection):
    """Buttwelded piping connections"""

    def __init__(self, element, point):
        super(Welding, self).__init__(element, point)

    def sif(self, code):
        if code == "B31.1":
            return 1.0


class SIFContainer(object):

    def __init__(self):
        self.Welding = Welding
        self.Unreinforced = Unreinforced
        self.Reinforced = Reinforced
        self.Weldolet = Weldolet
        self.Sockolet = Sockolet
        self.Sweepolet = Sweepolet
        self.Weldolet = Weldolet
        self.ButtWeld = ButtWeld
