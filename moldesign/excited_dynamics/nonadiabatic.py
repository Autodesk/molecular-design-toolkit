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
"""This module contains mixins that are required for a variety of nonadiabatic calculations."""
from moldesign.base.geometry import MultiStateGeometry

from moldesign.units import *


class NPITDC(MultiStateGeometry):
    """
    Multi-state dynamics mix-in class that calculates the non-adiabatic coupling
    between electronic states using Levine's NPI method.
    """
    def calc_tdc(self,wfns_i):
        """
        Returns the norm-preserving interpolation (NPI) time-derivative coupling (TDC)
        of Meek and Levine (dx.doi.org/10.1021/jz5009449)

        The TDC between all pairs of states is returned as a matrix;
        the element [i,j] is the TDC between states j and k,
        $ TDC_{jk} = < \phi_j(t) | \frac{d}{dt} \phi_k(t) > $

        :param wfns_i: Electronic wavefunctions (i.e., self.el_wfns) from the previous timestep. They must implement the .overlap method
        :return: $TDC_{jk}$ [1/TIME]
        """
        wfns_f = self.el_wfns

        tdolap = np.zeros((self.nstates,self.nstates))
        for j in xrange(self.nstates):
            for k in xrange(self.nstates):
                olap = wfns_i[j].overlap( wfns_f[k] )
                if olap>1.0:
                    if olap-1.0 < 1e-15: olap=1.0
                    else: raise ValueError('Overlap greater than 1: olap = 1+%f'%(olap-1.0))
                tdolap[j,k] = olap

        dt = self.timestep
        sn1 = np.arcsin(tdolap)
        cn1 = np.arccos(tdolap)

        TDC = np.zeros(tdolap.shape)/dt

        for j in xrange(tdolap.shape[0]):
            for k in xrange(tdolap.shape[1]):
                if j==k: continue #TODO: VERIFY THIS!!!!!!!!
                A = -sinx_div_x( cn1[j,j] - sn1[j,k] )
                B = sinx_div_x( cn1[j,j] + sn1[j,k] )
                C = sinx_div_x( cn1[k,k] - sn1[k,j] )
                D = sinx_div_x( cn1[k,k] + sn1[k,j] )

                if len(tdolap) == 2:
                    E=0.0
                else:
                    wlj = np.sqrt(1.0-tdolap[j,j]**2-tdolap[k,j]**2)
                    if wlj == 0.0:
                        E = 0.0
                    else:
                        wlk = (-tdolap[j,k]*tdolap[j,j] -tdolap[k,k]*tdolap[k,j]) / wlj
                        E = 2.0 * np.arcsin( wlj * (wlj * wlk * (np.sqrt(
                            (1.0-wlj**2)*(1.0-wlk**2) ) * np.arcsin(wlk) ) )) / (
                            np.arcsin(wlj)**2 - np.arcsin(wlk)**2 )

                TDC[j,k] = (cn1[j,j]*(A+B) + sn1[k,j]*(C+D) + E) / (2.0*dt)

        return TDC

def sinx_div_x(x):
    try:
        result = np.sin(x) / x
    except FloatingPointError:
        result = 1.0
    return result


class NACV_TDC(MultiStateGeometry):
    """
    Multi-state dynamics mix-in class that calculates the non-adiabatic coupling
    between electronic states using a non-adiabatic coupling vector.
    """

    def calc_tdc(self,wfns_i):
        """
        Returns nonadiabatic coupling
        $\sigma_{ij} = \mathbf{\dot{R}} \cdot < \phi_i (\mathbf{R}) | \nabla_{\mathbf{R}} | \phi_j (\mathbf{R}) >$
        :param wfns_i: Not used, passed only for compatibility
        :return: $TDC_{jk}$ [1/TIME]
        """
        TDC = np.zeros((self.nstates,self.nstates))/self.timestep

        for i in xrange(self.nstates):
            for j in xrange(self.nstates):
                if i==j: continue
                TDC[i,j] = ( self.nacv(i,j) * self.momenta/self.masses ).sum()

        return TDC

class MomentumAdjust(MultiStateGeometry):
    def adjust_energy(self, energy, start_state):
        """
        Adjust the energy of the geometry by scaling the momentum in the "vec" direction.
        self.electronic_state, and, therefore, self.potential_energy, must already be updated.
        :param energy: Adjust total energy to this value [ENERGY]
        :return: momentum scaling factor
        """

        vec = self.get_adjustment_direction(start_state)

        #Normalize it for numerical stability; gives A units of momentum
        vec = vec.magnitude / np.linalg.norm(vec)
        kinetic_target = energy - self.potential_energy
        kinetic_initial = self.kinetic_energy
        if kinetic_target.magnitude < 0.0:
            raise FrustratedHop(kinetic_target,self.kinetic_energy)

        #Solve quadratic equation: kinetic_target = (p-Av)^2/2m

        a = (vec * vec / self.masses).sum()
        b = -2.0*(vec * self.momenta / self.masses).sum()
        c = (self.momenta * self.momenta / self.masses).sum() - 2.0*kinetic_target

        #Make sure answers are real
        disc = b*b - 4.0*a*c
        if disc.magnitude < 0.0: raise FrustratedHop()

        sqdisc = np.sqrt(disc)
        r1 = (-b + sqdisc)/(2.0*a)
        r2 = (-b - sqdisc)/(2.0*a)

        A = r1 if abs(r1) < abs(r2) else r2

        self.momenta -= A * vec
        np.testing.assert_almost_equal(self.kinetic_energy.to(hartree).magnitude,
                                       kinetic_target.to(hartree).magnitude,
                                       decimal=10)
        self.log.info('Adjusted kinetic energy from %s to %s'%(str(kinetic_initial),
                                                               str(kinetic_target)))
        return A


class NACVAdjustment(MultiStateGeometry):
    def get_adjustment_direction(self,start_state):
        return self.nacv(start_state,self.electronic_state)
