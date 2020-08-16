"""Primitive bounding volumes"""

from psi.utils.euclid import Vector3


class AABB(object):
    """Axis aligned bounding box"""

    __slots__ = ("_min", "_max", "_is_init")

    def __init__(self):
        self._min = Vector3()
        self._max = Vector3()
        self._is_init = False

    def _get_min(self):
        return self._min[:]

    def _set_min(self, xmin, ymin, zmin):
        self._min[:] = xmin, ymin, zmin

    min = property(_get_min, lambda self, args: self._set_min(*args))

    def _get_max(self):
        return self._max[:]

    def _set_max(self, xmax, ymax, zmax):
        self._max[:] = xmax, ymax, zmax

    max = property(_get_max, lambda self, args: self._set_max(*args))

    def _get_extent(self):
        return self._min[:] + self._max[:]

    def _set_extent(self, xmin, ymin, zmin, xmax, ymax, zmax):
        """Set min and max at the same time"""
        self._min[:] = xmin, ymin, zmin
        self._max[:] = xmax, ymax, zmax

    extent = property(_get_extent, lambda self, args: self._set_extent(*args))

    def size(self):
        diag = self._max - self._min
        return diag[:]

    def center(self):
        ctr = 0.5 * (self._max+self._min)
        return ctr[:]

    def area(self):
        """Calculate the surface area of the box"""
        l, w, h = (self._max - self._min)[:]
        return 2 * (l*w + w*h + l*h)

    def contains(self, other, percent=0):
        """Check to determine if the other bound is inside self

        The percent parameter allows for a loose check against a bigger
        volume.

        Values for percent varies from 0 to 1 (ie 10% is 0.1).
        """
        loose = 0.5 * percent * (self._max-self._min)
        sminx, sminy, sminz = (self._min-loose)[:]  # loose extent
        smaxx, smaxy, smaxz = (self._max+loose)[:]
        ominx, ominy, ominz, omaxx, omaxy, omaxz = other._get_extent()
        return (sminx <= ominx and sminy <= ominy and sminz <= ominz and
                smaxx >= omaxx and smaxy >= omaxy and smaxz >= omaxz)

    def intersects(self, other):
        """Checks for intersection of two aabbs.  Edge cases are counted as
        valid intersections.
        """
        sxmin, symin, szmin, sxmax, symax, szmax = self._get_extent()
        oxmin, oymin, ozmin, oxmax, oymax, ozmax = other._get_extent()

        return not (oxmin >= sxmax or oxmax <= sxmin or
                    oymin >= symax or oymax <= symin or
                    ozmin >= szmax or ozmax <= szmin)

    def clear(self):
        """Set the min and max positions to default values"""
        self._set_extent(0, 0, 0, 0, 0, 0)
        self._is_init = False

    def update(self, vertices):
        """Update the AABB by passing in a new set of vertices"""
        assert len(vertices) % 3 == 0, "len must be a multiple of 3"
        xcomps = vertices[0::3]
        ycomps = vertices[1::3]
        zcomps = vertices[2::3]
        xmin = min(xcomps)
        ymin = min(ycomps)
        zmin = min(zcomps)
        xmax = max(xcomps)
        ymax = max(ycomps)
        zmax = max(zcomps)
        self._min.x = xmin
        self._min.y = ymin
        self._min.z = zmin
        self._max.x = xmax
        self._max.y = ymax
        self._max.z = zmax

    def merge(self, other):
        """Merge other's bounding volume with self's volume in-place"""
        if self._is_init:
            extent = (map(min, self._get_min(), other._get_min()) +
                      map(max, self._get_max(), other._get_max()))
            self._set_extent(*extent)
        else:   # prevents merging with default 0 volume
            extent = other._get_min() + other._get_max()
            self._set_extent(*extent)
            self._is_init = True
