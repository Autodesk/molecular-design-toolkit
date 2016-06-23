import math

import numpy as np

from moldesign import units as u
from moldesign.forcefields.forcefield import FFTerm
from moldesign.geom import coords as geo


class HarmonicBondTerm(FFTerm):
    def __init__(self, a1, a2, k, d0):
        self.atoms = [a1, a2]
        self.k = k
        self.d0 = d0
        self.bond = a1.bond_graph.get(a2, None)

    def __str__(self):
        return 'Harmonic bond ({atoms}), k={k.magnitude:6.3f}{k.units}, equil={d0.magnitude:6.2}{d0.units}'.format(
            atoms=','.join(map(str, self.atoms)),
            k=self.k,
            d0=self.d0)

    def coord(self):
        return geo.distance(*self.atoms)

    def energy(self):
        return u.default.convert(self.k * (self.coord() - self.d0)**2)


class HarmonicAngleTerm(FFTerm):
    def __init__(self, a1, a2, a3, k, theta0):
        self.atoms = [a1, a2, a3]
        self.k = k
        self.theta0 = theta0
        self.bonds = [a1.bond_graph.get(a2,None),
                      a2.bond_graph.get(a3,None)]

    def __str__(self):
        return 'Harmonic angle ({atoms}), k={k.magnitude:6.3f}{k.units}, equil={d0.magnitude:6.2}{d0.units}'.format(
            atoms=','.join(map(str, self.atoms)),
            k=self.k.defunits(),
            d0=self.theta0.defunits())

    def coord(self):
        return geo.angle(*self.atoms)

    def energy(self):
        return u.default.convert(self.k * (self.coord() - self.theta0)**2)


class PeriodicTorsionTerm(FFTerm):
    def __init__(self, a1, a2, a3, a4, n, v_n, gamma):
        self.atoms = [a1, a2, a3, a4]
        self.n = n
        self.v_n = v_n
        self.gamma = gamma
        self.bonds = [a1.bond_graph.get(a2, None),
                      a2.bond_graph.get(a3, None),
                      a3.bond_graph.get(a4, None)]

    def __str__(self):
        return 'Periodic torsion ({atoms}),' \
               ' n={n}, k={k.magnitude:6.3f}{k.units}, equil={d0.magnitude:6.2}{d0.units}'.format(
            atoms=','.join(map(str, self.atoms)),
            k=self.v_n.defunits(),
            n=self.n,
            d0=self.gamma)

    def coord(self):
        return geo.dihedral(*self.atoms)

    def energy(self):
        return u.default.convert(0.5 * self.v_n * (1.0 + np.cos(self.n * self.coord() - self.gamma)))


class LennardJonesSigmaEps(FFTerm):
    def __init__(self, atom, sigma, epsilon):
        self.atom = atom
        self.sigma = sigma
        self.epsilon = epsilon

    def get_sigma(self, other):
        return (self.sigma + other.sigma) / 2.0

    def get_epsilon(self, other):
        return math.sqrt(self.sigma * other.sigma)