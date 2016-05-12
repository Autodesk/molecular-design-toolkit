from moldesign.excited_dynamics import wavepackets
from moldesign.base import propagators
from moldesign.units import *
from base.atom import Atom


class AIMS_Atom(Atom):
    def __init__(self,*args,**kwargs):
        super(AIMS_Atom,self).__init__(*args,**kwargs)
        self.width = 0.0*def_length


class AIMS_BasisFn(propagators.VelocityVerletMD, wavepackets.GaussianWavepacket):
    """
    Holds a time-dependent nuclear/electronic product basis function.
    In practical terms, this is a nuclear configuration with
    specified gaussian widths, an associated electronic wavefunction.
    """
    def __init__(self,*args,**kwargs):
        super(AIMS_BasisFn,self).__init__(*args,**kwargs)
        self.state_id=-1
        self.semiclassical_phase = 0.0 * radians
        self.amplitude = 0.0j





class AIMS_Wavefn(object):
    """
    Holds a time-dependent nuclear/electronic product basis function.
    In practical terms, this is a nuclear configuration with
    specified gaussian widths, an associated electronic wavefunction.
    """
    def __init__(self):
        pass
