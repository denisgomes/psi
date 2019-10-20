"""An implementation of 3 dimensional vectors for basic geometric usage."""


import math


class Vec3:

    __slots__ = ("x", "y", "z")

    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        """Add two vectors together"""
        return Vec3(self.x+other.x, self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        """Subtract one vector from another"""
        return Vec3(self.x-other.x, self.y-other.y, self.z-other.z)

    def __mul__(self, other):
        """Multiply two vectors or scale a single vector"""
        if isinstance(other, (int, float)):
            return Vec3(other*self.x, other*self.y, other*self.z)
        elif isinstance(other, Vec3):
            return Vec3(self.x*other.x, self.y*other.y, self.z*other.z)
        else:
            raise NotImplementedError("unsupported operation")

    __rmul__ = __mul__

    def __truediv__(self, other):
        """Divide two vectors or scale a single vector"""
        if isinstance(other, (int, float)):
            return Vec3(other/self.x, other/self.y, other/self.z)
        elif isinstance(other, Vec3):
            return Vec3(self.x/other.x, self.y/other.y, self.z/other.z)
        else:
            raise NotImplementedError("unsupported operation")

    __rtruediv__ = __truediv__

    def norm(self):
        """Get the magnitude or length of a vector"""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def norm2(self):
        """Get the magnitude or length of a vector squared"""
        return self.x**2 + self.y**2 + self.z**2

    def unit(self):
        """Get the unit directional vector with length equals 1"""
        norm = self.norm()
        return Vec3(self.x/norm, self.y/norm, self.z/norm)

    def __repr__(self):
        return f"{self.__class__.__name__}(x={self.x}, y={self.y}, z={self.z})"
