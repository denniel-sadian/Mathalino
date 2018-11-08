"""
Python 3
July 2, 2017
Denniel Luis Saway Sadian
"""

from fractions import Fraction as Fraction
from math import sin, cos, tan, radians, pi


class ArithmeticSequence:
    """ Arithmetic sequence for fractions and integers or floats.
        common difference: An = A1 + (n - 1)d
                      sum: Sn = n/2 [2(a1) + (n-1)d] """

    def __init__(self, sequence: list):
        """ The sequence argument is a list, so the user can give repeating
            values for testing. Fraction terms, integer or float terms
            can be combined. """
        # passing the list to self.sequence variable
        self.sequence = list(sequence)
        # special list used in checking whether the given
        # list is an arithmetic sequence
        # all items from this list must be 'yes' in order to be
        # an arithmetic sequence,
        # otherwise it is not an arithmetic sequence
        self.correct = list()
        # getting the common difference from self.sequence
        if len(self.sequence) > 1:
            self.cd = self.sequence[1] - self.sequence[0]
            if self.is_it():  # getting the sum if an arithmetic sequence
                self.sum = (len(self.sequence) / 2) * \
                           ((2 * self.sequence[0]) +
                            (len(self.sequence) - 1) * self.cd)
        else:
            self.cd = self.sequence[0]

    def is_it(self):
        """ Checking whether the sequence is an arithmetic sequence. """
        for i in self.sequence:
            try:
                if i + self.cd == self.sequence[self.sequence.index(i) + 1]:
                    self.correct.append('yes')
                else:
                    self.correct.append('no')
            except IndexError:
                pass
        return 'no' not in self.correct

    def expand(self, span):
        """ Expanding the arithmetic sequence """
        if self.is_it():
            for i in range(span):
                self.sequence.append(self.sequence[-1] + self.cd)
            return ArithmeticSequence(self.sequence)
        else:
            raise ValueError(f'{self.sequence} is not an arithmetic sequence.')

    def between_terms(self, numbers):
        """ Finding terms between an arithmetic sequence.
            example: (3, _, _, _, 15)
            If you have like this (3, 2, 1, _, _, _, -3),
            then only consider this (1, _, _, _, -3)
            equivalent to: x = ArithmeticSequence([1, -3])
                           x.between_terms(3)
                           print(x.sequence)
                       ... [1, 0.0, -1.0, -2.0, -3.0] """
        if len(self.sequence) == 2:
            self.cd = (self.sequence[1] - self.sequence[0]) / (numbers + 1)
            self.sequence = [self.sequence[0]]
            for i in range(numbers + 1):
                self.sequence.append(self.sequence[0] + (i + 2 - 1) * self.cd)
        else:
            raise ValueError("sequence contained more than two terms, "
                             "but two were only allowed.")


class GeometricSequence:
    """ Geometric Sequence for fractions, integers or floats
        common ratio: An = A1 r ** n-1 """

    def __init__(self, sequence):
        self.sequence = list(sequence)
        self.correct = list()
        if len(self.sequence) > 1:
            if self.sequence[0] < self.sequence[1]:
                self.cr = Fraction(self.sequence[0] / self.sequence[1])
            else:
                self.cr = self.sequence[0] / self.sequence[1]
        else:
            self.cr = self.sequence[0]

    def is_it(self):
        for i in self.sequence:
            try:
                if i / self.cr == self.sequence[self.sequence.index(i) + 1]:
                    self.correct.append('yes')
                else:
                    self.correct.append('no')
            except IndexError:
                pass
        return 'no' not in self.correct

    def expand(self, span):
        if self.is_it():
            for i in range(span):
                self.sequence.append(self.sequence[-1] / self.cr)
            return GeometricSequence(self.sequence)
        else:
            raise ValueError(f'{self.sequence} is not a geometric sequence.')


def square(side):
    """ Square area: A = side ** 2 """
    return round(side ** 2, 2)


def rectangle(base, height):
    """ Rectangle area: A = base * height """
    return round(base * height, 2)


def circle(radius):
    """ Circle area: A = radius * radius * pi """
    return round(radius * radius * pi, 2)


def triangle(base, height):
    """ Triangle area: A = (base * height) / 2 """
    return round((base * height) / 2, 2)


def trapezoid(base1, base2, height):
    """ Trapezoid area: A = ((base1 + base2) * height) / 2 """
    return round(((base1 + base2) * height) / 2, 2)


def cube(edge):
    """ Cube volume: V = edge ** 3 """
    return round(edge ** 3, 2)


def rectangular_prism(width, height, length):
    """ Rectangular prism volume: V = width * height * length """
    return round(width * height * length, 2)


def irregular_prism(base_area, height):
    """ Irregular prism volume: V = base_area * height """
    return round(base_area * height, 2)


def cylinder(radius, height):
    """ Cylinder volume: V = pi * radius ** 2 """
    return round(pi * radius ** 2 * height, 2)


def pyramid(length, width, height):
    """ Pyramid volume: V = (length * width * height) / 3 """
    return round((length * width * height) / 3, 2)


def cone(radius, height):
    """ Cone volume: V = (pi * radius ** 2 * height) / 3 """
    return round((pi * radius ** 2 * height) / 3, 2)


def circumference(diameter):
    """ Circumference: diameter * pi """
    return round(diameter * pi, 2)


def perimeter_s(side):
    """ Square perimeter: side * 4 """
    return round(side * 4)


def perimeter_r(length, height):
    """ Rectangle perimeter: (length + length) + (height + height) """
    return round((length + length) + (height + height), 2)


def speed(d, t):
    """ Scalar speed: S = d / t """
    return round(float(d) / float(t), 2)


def distance(s, t):
    """ Scalar distance: D = s * t """
    return round(float(s) * float(t), 2)


def time_s(d, s):
    """ Scalar time: T = d / s """
    return round(float(d) / float(s), 2)


def velocity(d, t):
    """ Vector velocity: V = d / t """
    return round(float(d) / float(t), 2)


def displacement(v, t):
    """ Vector displacement: D = v / t """
    return round(float(v) * float(t), 2)


def time_v(d, v):
    """ Vector time: T = d / v """
    return round(float(d) / float(v), 2)


def trigo_ratios(theta, hypotenuse, adjacent, opposite):
    """ Finds trigonometric ratios and returns a dictionary """
    # solving for adjacent, opposite and angle
    if opposite == 0 and adjacent == 0:
        return {'adjacent': round(hypotenuse * cos(radians(theta)), 2),
                'opposite': round(hypotenuse * sin(theta), 2),
                'hypotenuse': hypotenuse, 'angle': 90 - theta}
    # solving for opposite, hypotenuse and angle
    elif opposite == 0 and hypotenuse == 0:
        return {'adjacent': adjacent,
                'opposite': round(adjacent * tan(radians(theta)), 2),
                'hypotenuse': round(adjacent / cos(radians(theta)), 2),
                'angle': 90 - theta}
    # solving for adjacent, hypotenuse and angle
    elif adjacent == 0 and hypotenuse == 0:
        return {'adjacent': round(opposite * tan(radians(theta)), 2),
                'opposite': opposite,
                'hypotenuse': round(opposite * sin(radians(theta)), 2),
                'angle': 90 - theta}


def average(*args):
    """ Finds average """
    return round(sum(args) / len(args), 2)


def molar_mass(convert, **kw):
    """ Conversion of molar mass """
    if convert is 'gm' and kw['g'] and kw['a']:  # gram to mole
        return round((kw['g'] * 1) / kw['a'], 2)
    elif convert is 'mg' and kw['g'] and kw['a']:  # mole to gram
        return round((kw['g'] * kw['a']) / 1, 2)
    elif convert is 'am' and kw['g']:  # atom to mole
        return round((kw['g'] * 1) / (6.022 * 10 ** 23), 2)
    elif convert is 'ma' and kw['g']:  # mole to atom
        return round((kw['g'] * (6.022 * 10 ** 23)) / 1, 2)
    elif convert is 'ag' and kw['g'] and kw['a']:  # atom to gram
        return round((((kw['g'] * 1) / (6.022 * 10 ** 23)) * kw['a']) / 1, 2)
    elif convert is 'ga' and kw['g'] and kw['a']:  # gram to atom
        return round((((kw['g'] * 1) / kw['a']) * (6.022 * 10 ** 23)) / 1, 2)


if __name__ == '__main__':  # tests
    s0 = ArithmeticSequence([Fraction('1/6'), Fraction('1/3'),
                             Fraction('1/2'), Fraction('2/3')])
    s1 = ArithmeticSequence([3, 6, 9, 12])
    s2 = ArithmeticSequence([12, 9, 6, 3])
    s3 = ArithmeticSequence([7, 16, 25, 34, 43, Fraction('52/1')])
    s4 = ArithmeticSequence([Fraction('1/6'), Fraction('1/3')])
    x = s4.expand(2)
    s5 = ArithmeticSequence([3, 2])
    z = s5.expand(5)
    print(f'{s0.sequence} = {s0.cd} {s0.is_it()} {s0.sum}')
    print(f'{s1.sequence} = {s1.cd} {s1.is_it()} {s1.sum}')
    print(f'{s2.sequence} = {s2.cd} {s2.is_it()} {s2.sum}')
    print(f'{s3.sequence} = {s3.cd} {s3.is_it()} {s3.sum}')
    print(x.sequence, x.cd, x.is_it(), x.sum)
    print(z.sequence, z.cd, z.is_it(), z.sum)
    x = ArithmeticSequence([Fraction('1/6'), Fraction('2/3')])
    x.between_terms(2)
    z = ArithmeticSequence([1, -3])
    z.between_terms(3)
    print(x.cd, x.sequence, x.is_it(), x.sum)
    print(z.cd, z.sequence, z.is_it(), z.sum)
    print(square(12), 'square')
    print(rectangle(12, 12), 'rectangle')
    print(circle(12), 'circle')
    print(triangle(12, 12), 'triangle')
    print(trapezoid(12, 12, 12), 'trapezoid')
    print(cube(12), 'cube')
    print(rectangular_prism(12, 12, 12), 'rectangular prism')
    print(irregular_prism(12, 12), 'irregular prism')
    print(cylinder(12, 12), 'cylinder')
    print(pyramid(12, 12, 12), 'pyramid')
    print(cone(12, 12), 'cone')
    print(perimeter_s(12), 'square perimeter')
    print(perimeter_r(12, 12), 'rectangle perimeter')
    print(speed(12, 2), 'speed')
    print(distance(speed(12, 2), 2), 'distance')
    print(time_s(distance(speed(12, 2), 2), speed(12, 2)), 'time')
    print(velocity(12, 2), 'velocity')
    print(displacement(velocity(12, 2), 2), 'displacement')
    print(time_v(displacement(velocity(12, 2), 2), velocity(12, 2)), 'time')
    print(trigo_ratios(12, 2, 0, 0), 'trigonometry')
    print(average(98, 100, 99, 100), 'average')
    print(molar_mass('gm', g=12, a=14))
    print(molar_mass('mg', g=12, a=14))
    print(molar_mass('am', g=12))
    print(molar_mass('ma', g=12))
    print(molar_mass('ag', g=12, a=14))
    print(molar_mass('ga', g=12, a=14))
    print(molar_mass('ga', g=12, a=14))
    gs = GeometricSequence([5.2, 2.6, 1.3, 0.65])
    gs1 = gs.expand(3)
    print(gs1.cr, gs1.is_it(), gs1.sequence)
