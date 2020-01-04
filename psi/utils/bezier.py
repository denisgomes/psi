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

"""Utility functions used throughout the project"""


# http://stackoverflow.com/questions/7267226/range-for-floats
def frange(start, stop, step=1.0):
    while start < stop:
        yield start
        start += step


# http://stackoverflow.com/questions/246525/how-can-i-draw-a-bezier-curve-using-pythons-pil
def bezier(xys, weights):
    '''Each control point must have a weight greater than zero.  If all the
    weights are 1, then it is a regular bezier curve.  The points locating the
    start and the end of the curve also have control points but they have a
    weight of 1.

    Modified the example taken from startoverflow to generate a rational bezier
    curve which can exactly represent conics sections such as arcs.  When
    working with CAD models accuracy is very important for viable results.
    '''
    n = len(xys)
    combinations = pascal_row(n-1)

    def bezier(ts):
        result = []
        for t in ts:
            tpowers = (t**i for i in range(n))
            upowers = reversed([(1-t)**i for i in range(n)])
            coefs = [c*a*b for c, a, b in zip(combinations, tpowers, upowers)]
            result.append(
                tuple(sum([coef*p*w for coef, p, w in zip(coefs, ps, weights)])/
                    sum([coef*w for coef, w in zip(coefs, weights)]) for ps in zip(*xys)))
        return result
    return bezier


def pascal_row(n):
    # This returns the nth row of Pascal's Triangle
    result = [1]
    x, numerator = 1, n
    for denominator in range(1, n//2+1):
        # print(numerator,denominator,x)
        x *= numerator
        x /= denominator
        result.append(x)
        numerator -= 1
    if n & 1 == 0:
        # n is even
        result.extend(reversed(result[:-1]))
    else:
        result.extend(reversed(result))
    return result
