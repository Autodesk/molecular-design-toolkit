from types import FunctionType

import numpy as np
import sympy
from moldesign import DotDict
from sympy.abc import x

from moldesign.basemethods import EnergyModelBase, Parameter
from moldesign.keywords import ElectronicWfn
import moldesign.units as u


class DiabatizedWfn(ElectronicWfn):
    def __init__(self,statevec=None):
        if statevec:
            self.wfn = [ VectorOrbital(s) for s in statevec ]


class VectorOrbital(object):
    def __init__(self,vector):
        self.vector=vector

    def overlap(self,other):
        return self.vector.dot(other.vector)


class Analytical2StatePotential(EnergyModelBase):

    PARAMETERS = DotDict(v11=Parameter(types=[FunctionType]),
                         v12=Parameter(types=[FunctionType]),
                         v22=Parameter(types=[FunctionType]))
    DEFAULT_PROPERTIES = ['state_energies','state_forces',
                          'ci_vecs','nacv']
    ALL_PROPERTIES = DEFAULT_PROPERTIES

    def __init__(self,v11=None,v22=None,v12=None):
        """
        Generates analytical adiabatic potential based on diabatic 2-state problems.
        Works with unit-less quantitites, but will output in atomic units.
        :param v11: sympy expression for first diabat
        :param v22: sympy expression for second diabat
        :param v12: sympy expression for diabatic coupling
        :param defaults: values of free parameters
        :return: ElectronicInterface object that generates the adiabatic curves
        """
        super(Analytical2StatePotential,self).__init__()
        hamiltonian = sympy.Matrix([[v11,v12],[v12,v22]])
        eigen0,eigen1 = hamiltonian.eigenvects()
        self.s0_energy,mult0,ev0 = eigen0
        self.s1_energy,mult1,ev1 = eigen1
        ev0 = ev0[0]
        ev1 = ev1[0]
        self.ev0_expr = ev0
        self.ev1_expr = ev1
        norm0 = ev0.norm()
        norm1 = ev1.norm()
        dev1dx = ev1.diff(x)
        self.nacv01_expr = (ev0.dot(dev1dx)) /(norm0 * norm1)
        self.lambdify()

    def lambdify(self):
        # Turn expressions into functions
        self.S0 = sympy.lambdify(x, self.s0_energy)
        self.S1 = sympy.lambdify(x, self.s1_energy)
        self.forceS0 = sympy.lambdify(x, -self.s0_energy.diff(x))
        self.forceS1 = sympy.lambdify(x, -self.s1_energy.diff(x))
        self.nacv01 = sympy.lambdify(x, self.nacv01_expr)
        self.ev0 = sympy.lambdify(x, self.ev0_expr)
        self.ev1 = sympy.lambdify(x, self.ev1_expr)

    def calculate(self,requests):
        geom = self.mol
        xpos = geom.positions[0].value_in(u.a0)
        energies = np.array([self.S0(xpos),self.S1(xpos)])*u.hartree

        #Forces - just set the first element
        forces = [ np.zeros(geom.ndim)*u.hartree/u.a0  for i in xrange(2)]
        forces[0][0] = self.forceS0(xpos) * u.hartree/u.a0
        forces[1][0] = self.forceS1(xpos) * u.hartree/u.a0

        #non-adiabatic coupling - just set the first element
        nacv = {}
        for pair in ((0,1),(1,0)):
            nacv[pair] = np.zeros(geom.ndim) / u.a0
        try:
            nacv[0,1][0] = self.nacv01(xpos) / u.a0
        except ZeroDivisionError:
            pass
        nacv[1,0] = -nacv[0,1]
        statevec = [ self.ev0_expr.evalf(subs={x:xpos}).normalized(),
                     self.ev1_expr.evalf(subs={x:xpos}).normalized() ]
        return dict(state_energies=energies,
                    state_forces=forces,
                    nacv=nacv,
                    civecs=statevec)

    def __getstate__(self):
        odict = self.__dict__.copy()
        for name in 'S0 S1 forceS0 forceS1 nacv01 ev0 ev1'.split():
            del odict[name]
        return odict

    def __setstate__(self,idict):
        self.__dict__.update(idict)
        self.lambdify()



class Piecewise1DPotential(object):
    def __init__(self,p_left,p_right):
        self.p_left = p_left
        self.p_right = p_right

    def compute(self,geom):
        xpos = geom.positions[0]
        if xpos < 0.0 * u.angstrom:
            return self.p_left.compute(geom)
        else:
            return self.p_right.compute(geom)

#Note that Tully problems specify 2000 m_e mass particle
def Tully1(A=0.01,B=1.6,C=0.005,D=1.0):
    v11_right = A*(1.0-sympy.exp(-B*x))
    v11_left = -A*(1.0-sympy.exp(B*x))
    v12 = C*sympy.exp(-D*x**2)
    left_pot = Analytical2StatePotential(v11_left, -v11_left, v12)
    right_pot = Analytical2StatePotential(v11_right, -v11_right, v12)
    return Piecewise1DPotential(left_pot,right_pot)

def Tully2(A=0.1, B=0.28, C=0.015, D=0.06, E0 = 0.05):
    v11 = 0.0
    v22 = -A * sympy.exp( -B * x**2 ) + E0
    v12 = C * sympy.exp( -D*x**2 )
    return Analytical2StatePotential(v11, v22, v12)

def Tully3(A=6e-4,B=0.1,C=0.9):
    v11 = A
    v22=-v11
    v12_left = B*sympy.exp(C*x)
    v12_right = B*(2-sympy.exp(-C*x) )
    left_pot = Analytical2StatePotential(v11, v22, v12_left)
    right_pot = Analytical2StatePotential(v11, v22, v12_right)
    return Piecewise1DPotential(left_pot,right_pot)
