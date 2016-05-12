import math
import numpy as np
import moldesign as mdt
import moldesign.geometry as geo
import moldesign.units as u

from moldesign.interfaces.ambertools import assign_forcefield


class FFTerm(object):
    pass


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


class FFParameters(object):
    """
    This object contains assigned force field parameters for a specific system
    The base focuses on the AMBER / CHARMM - type force
    """

    # TODO: this needs to describe attenuation for things close together
    # TODO: deal with nonbonded exceptions
    def __init__(self, bonds, angles, dihedrals, partial_charges, lennard_jones):
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals

        self.partial_charges = partial_charges  # maps atoms to their partial charges
        self.lennard_jones = lennard_jones  # maps atoms to LJ terms

        # lookups
        self.bond_term = {term.bond: term for term in self.bonds}
        self.angle_term = {tuple(term.atoms): term for term in self.angles}
        for term in self.angles: self.angle_term[tuple(reversed(term.atoms))] = term
        self.dihedral_term = {}
        for term in self.dihedrals:
            self.dihedral_term.setdefault(tuple(term.atoms), []).append(term)
            self.dihedral_term.setdefault(tuple(reversed(term.atoms)), []).append(term)



class ResidueTemplate(object):
    def __init__(self, mol, charges, ffparams=None, minimized_mol=None):
        self.mol = mol
        self.charges = {atom.pdbname:charge for atom,charge in charges.iteritems()}

        if ffparams is None:
            self.ffparams = self.get_ffparams()
        else:
            self.ffparams = ffparams

        self.minimized_mol = minimized_mol

    def get_ffparams(self):
        chargemap = {atom: self.charges[atom.pdbname] for atom in self.mol.atoms}
        return mdt.interfaces.ambertools.get_gaff_parameters(self.mol, chargemap)




__all__ = ['assign_forcefield']