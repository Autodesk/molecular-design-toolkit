"""
This module mostly contains potential energy classes for users to import.
TODO: look into http://molmod.github.io/
TODO: log errors rather than just printing them
"""
import moldesign
import numpy as np
from moldesign.basemethods import EnergyModelBase
from moldesign import units as u

try:
    from moldesign.interfaces.openmm import OpenMMPotential
except ImportError:
    print 'WARNING: OpenMM interface could not be imported'

try:
    from moldesign.interfaces.pyscf_interface import PySCFPotential
except (ImportError, OSError):
    print 'WARNING: PySCF interface could not be imported'
    # TODO: convert to logging


"""Gives the user our default RHF implementation"""
def RHF(**kwargs):
    return PySCFPotential(theory='rhf', **kwargs)


class HarmonicOscillator(EnergyModelBase):
    def __init__(self, k):
        self.k = k
        assert k.dimensionality == {'[mass]':1,'[time]':-2}
        super(EnergyModelBase,self).__init__()

    def calculate(self, requests):
        x = self.mol.atoms[0].x
        energy = 0.5 * self.k * (x**2)
        forces = np.zeros(self.mol.ndims) * u.hartree / u.bohr
        forces[0] = - self.k * x
        return dict(potential_energy=energy, forces=forces)


class Spring(EnergyModelBase):
    """Two atoms attached by a spring"""
    def __init__(self, d0, k):
        self.d0 = d0
        self.k = k
        super(Spring, self).__init__()

    def prep(self):
        assert self.mol.natoms == 2
        self._prepped = True

    def calculate(self, requests):
        dvec = self.mol.atoms[1].position - self.mol.atoms[0].position
        dist = np.sqrt(dvec.dot(dvec))
        pe = 0.5 * self.k * (dist - self.d0)**2
        f = self.k * dvec * (dist - self.d0) / dist
        forcearray = u.to_units_array(f.tolist() + (-f).tolist())
        return {'potential_energy': pe,
                'forces': forcearray}
