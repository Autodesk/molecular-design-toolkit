# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import copy
import unittest

from moldesign.base import geometry

from moldesign.units import *




#Default widths taken from Thompson, Punwong, and Martinez, Chem Phys 370 (2010), 70-77.
#doi:10.1016/j.chemphys.2010.03.020
#These are reported in Bohr, and are therefore related to their alpha by (2*alpha)^-1/2
from base.atom import Atom

martinez_alpha={1:4.7/a0**2,
                6:22.7/a0**2,
                7:19.0/a0**2,
                8:12.2/a0**2,
                9:8.5/a0**2,
                16:16.7/a0**2,
                17:7.4/a0**2}
def_width = {ATOMIC_NUMBERS: 2.0*alpha**(-0.5)
             for atnum,alpha in martinez_alpha.items()}

class FrozenGaussianAtom(Atom):
    """
    Atom augmented with wavepacket/basis set information
    """
    def __init__(self,*args,**kwargs):
        super(FrozenGaussianAtom,self).__init__(*args,**kwargs)
        if 'width' not in kwargs:
            self.width = def_width[self.atnum]

class GaussianWavepacket(geometry.MDGeometry):
    """
    Associates a gaussian wavepacket with the geometry. Wavepacket has the form
    $\Psi(r;r_0,p_0,a) = \prod \left( 2 a \pi \right)^{-\frac{1}{4}} \exp( \frac{(r-r_0)^2}{2a} + \frac{ i p_0}{\hbar} (r-r_0)]$
    where a = width**2
    """
    def __init__(self,*args,**kwargs):
        super(GaussianWavepacket,self).__init__(*args,**kwargs)
        self.widths = self._init_dim_array('width',a0)
        self.prefactor = 1.0 / np.sqrt( np.product( self.widths * sqrtpi ) )


    def __call__(self,R):
        """
        Returns wavepacket amplitude at position R
        :param R: self.ndims dimensional vector of positions
        :return: complex wavepacket amplitude
        """
        oscill = np.exp( (1.0j*self.momenta/hbar)*(R-self.positions) )
        envelope = np.exp( -(R-self.positions)**2/(2.0*self.widths**2))
        return self.prefactor * np.product( oscill * envelope )

    def __mul__(self,othr):
        """
        Returns product of two gaussian wavepackets, (i.e. centroid)
        using self's complex conjugate
        :param othr: other gaussian wavepacket
        :return: product wavepacket
        """
        result=copy.deepcopy(self)
        #convert widths to prefactor form
        a = 1/(2.0*self.widths**2)
        b = 1/(2.0*othr.widths**2)

        result.widths = (2*(a+b))**(-0.5)
        result.positions = (a*self.positions + b*othr.positions)/(a+b)
        result.momenta = -self.momenta + othr.momenta

        prefactor = np.conj(self.prefactor) * othr.prefactor
        for idim in xrange(self.ndim):
            prefactor *= self.prefac_1d(self.positions[idim], self.momenta[idim], a[idim],
                                      othr.positions[idim], othr.momenta[idim], b[idim])
        result.prefactor = prefactor
        return result


    def prefac_1d( self, x1, p1, a1, x2, p2, a2):
        newa = -(a1*a2)/(a1+a2)
        newx = (a1*x1 + a2*x2) / (a1+a2)
        newp = p2-p1
        result = np.exp( (newa * ((x1-x2)**2)) +
                         ((1.0j/hbar) * (newp*newx + x1*p1 - x2*p2)) )
        return result

    def olap_1d(self, x1, p1, a1, x2, pw, a2):
        new_width = (a1+a2)
        prefac = self.prefac_1d( x1, p1, a1, x2, p2, a2)
        return new_width * prefac * sqrtpi * sqrt2

    def overlap(self,other):
        prod = self*other
        return prod.integrate()

    def kinetic_matrix_element(self, othr):
        """ <CG1 | T | CG2>
        """
        #shamelessly copied from FMS repository
        a = 1/(2.0*self.widths**2)
        b = 1/(2.0*othr.widths**2)
        dR = self.positions-other.positions
        pij = b*self.momenta + a*othr.momenta
        total = 0.0j
        prefacs = ( 4.0j * a * b * dR * pij +
                    2.0 * a * b *(a+b) +
                    4* (R**2) * (a**2) * (b**2) +
                    pij**2/(a+b)**2 )
        for idim in xrange(self.ndim):
            total += prefacs[idim] *\
                self.olap_1d(self.positions[idim], self.momenta[idim], a[idim],
                othr.positions[idim], othr.momenta[idim], b[idim])

        return total

    def integrate(self):
        """
        Returns integral over all space
        """
        return self.prefactor * np.product( self.widths * sqrtpi * sqrt2 )




def _test_setup():
    """
    Create a 3D test wavepacket
    :return: wavepacket object
    """
    atom = FrozenGaussianAtom(atnum = 6, mass=12.0*amu, ndim=3)
    atom.pos = np.array([0.0, 0.5,0.25]) * a0
    atom.mom = np.array([0.0, 0.0, -0.05]) * amu * a0 / fs
    wavepacket1 = GaussianWavepacket([atom],None)

    atom = FrozenGaussianAtom(atnum = 6, mass=12.0*amu, ndim=3)
    atom.pos = np.array([0.0, 0.0, 0.0]) * a0
    atom.mom = np.array([0.0, 0.0, -0.05]) * amu * a0 / fs
    wavepacket2 = GaussianWavepacket([atom],None)
    return wavepacket1,wavepacket2


class TestWavepacketMath(unittest.TestCase):
    def test_multiplication(self):
        """Tests multiplication of wavepackets by comparing analytical and numerical results"""
        wp1,wp2 = _test_setup()
        for iii in range(100):
            testpos = np.random.rand(3)*a0
            numprod = np.conj(wp1(testpos)) * wp2(testpos)
            anaprod = (wp1*wp2)(testpos)
            self.assertAlmostEqual(numprod/anaprod,1.0,7,
                                   msg="numprod=%s, anaprod=%s"%(numprod,anaprod))

if __name__=='__main__':
    unittest.main()