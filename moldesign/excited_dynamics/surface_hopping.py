from moldesign.base.geometry import MultiStateGeometry,MDGeometry

from moldesign.units import *


#TODO: Make a hop its own object?


class SurfaceHoppingGeometry(MultiStateGeometry,MDGeometry):
    def __init__(self,atoms,nstates=2,**kwargs):
        super(SurfaceHoppingGeometry,self).__init__(atoms,**kwargs)
        self.nstates = nstates
        self.amps = np.array( [ (1.0+0.j) if self.electronic_state==i else (0.0j)
                                for i in xrange(self.nstates)] )

    def ampdot(self,t,amps,energies,tdc):
        """
        Calculate surface hopping electron state time derivatives
        $ i \hbar \dot{c}_j = \sum_n c_n \left( V_{n}\delta{jn} - i \hbar \sigma_{jn} $
        :return: time derivatives of complex amplitudes for each state [1/TIME]
        These assume the amplitude is of the form A(t)*exp( V*t/(i hbar) ),
        and return the time derivative of A(t)
        :return type:
        """
        nstates=len(amps)
        adot = (1.0j)*np.zeros( len(energies) )*tdc[0,0]
        for i in xrange(nstates):
            for j in xrange(nstates):
                adot[i] += amps[j] * tdc[i,j] * np.exp(
                    -imi*(energies[j]-energies[i])*t / hbar )
        return adot


    def propagate_hop(self):
        start_state = self.electronic_state
        start_el_energies = self.el_energies
        start_energy = self.total_energy
        start_wfn = self.el_wfns

        self.propagate(timestep = self.timestep)

        self.tdc = self.calc_tdc(start_wfn)
        self.propagate_amplitudes(self.timestep, start_el_energies)

        hopped = self.hop()
        if hopped:
            try:
                self.adjust_energy(start_energy,start_state)
            except FrustratedHop:
                self.electronic_state = start_state
                self.log.info('Hop failed, remaining on state %d'%start_state)



    def propagate_amplitudes(self,timestep,initial_energies,thresh=1.e-6):
        """
        Uses an adaptive 4th-order Runge Kutta integrator.
        Integrates from starttime to self.time
        Energies are interpolated from the start to the end of the timestep
        Uses c(t)=A(t)exp(-i Vi t/hbar) for numerical stability (with t=0 @ starttime)
        This is probably the interaction picture
        :param timestep: [TIME]
        :param initial_energies: state energies at starttime [ENERGY]
        :param thresh: Threshold for norm convergence.
        """
        energy_interpolator = interpolator(0.0*timestep,initial_energies,
                                           timestep,self.el_energies)

        for islice in xrange(8):
            nsteps = 2**islice
            dt = timestep/(1.0*nsteps)
            time = 0.0
            amps = self.amps.copy()

            #At each step, integrate A(t) from t to t+dt
            for istep in xrange(nsteps):
                k1 = self.ampdot( dt*0.0, amps, energy_interpolator(0.0*dt) , self.tdc)
                k2 = self.ampdot( dt/2.0,amps+k1*dt/2.0, energy_interpolator(dt/2.0), self.tdc)
                k3 = self.ampdot( dt/2.0,amps+k2*dt/2.0, energy_interpolator(dt/2.0), self.tdc)
                k4 = self.ampdot( dt, amps+k3*dt, energy_interpolator(dt), self.tdc)

                amps += dt*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
                #transform from A(t) to c(t)
                amps = amps * np.exp(-imi*self.el_energies*dt/hbar)
                time += dt

            #Check that runge kutta has converged
            assert time == timestep
            thisnorm = np.linalg.norm(amps)
            if islice>0:
                if np.abs( (thisnorm - lastnorm)/lastnorm) < thresh:
                    self.log.info('Runge Kutta converged, %d steps'%nsteps)
                    break
            lastnorm = thisnorm

        else:
            raise ValueError('Runge Kutta did not converge')

        self.amps = amps


class FSSH(SurfaceHoppingGeometry):
    """
    Mixin that adds the fewest-switches surface hopping prescription.
    See DOI: 10.1021/jp512024w for a straightforward description
    :return: True if hopped, otherwise false
    """
    def hop(self):
        i_st = self.electronic_state
        a = np.outer(self.amps , self.amps.conj())
        probs = np.zeros(self.nstates)
        for f_st in xrange(self.nstates):
            if i_st == f_st: continue

            numer = self.timestep * self.tdc[i_st,f_st]
            if np.abs(numer*a[f_st,i_st]) < 1.0e-12: #Call this numerical 0
                probs[f_st] = 0.0
            else:
                probs[f_st] = (2.0 * self.timestep * self.tdc[i_st,f_st] ) * (
                    a[f_st,i_st]/a[i_st,i_st]).real

            if probs[f_st] < 0.0: probs[f_st] = 0.0

        #Switch state
        rand = self.random.random()
        cumprob = 0.0
        self.sh_probabilities = probs
        self.sh_rand = rand
        for istate,prob in enumerate(probs):
            cumprob += prob
            if cumprob > rand:
                self.log.info('FSSH hop from state %d to %d (prob=%f)'%
                              (self.electronic_state, istate, prob))
                assert istate != self.electronic_state
                self.electronic_state = istate
                return True
        return False

class NoHopping(SurfaceHoppingGeometry):
    """
    Mixin to make a "surface hopping" trajectory that never actually hops.
    Useful to examine 1st order scattering onto other electronic states
    """
    def hop(self): return False


def interpolator( xi, yi, xf, yf):
    dx = xf-xi
    def interp(x):
        y = yi*(x-xi)/dx + yf*(xf-x)/dx
        return y
    return interp


def test_amplitude_prop():
    from moldesign.util import test_support as TS

    mol = SurfaceHoppingGeometry( TS.ethyl_atoms(),
        electronic_state=1,
        potential_model=None )
    mol.el_energies = [0.010,0.020] * hartree
    mol.tdc = np.array([[0.0,0.1],[-0.1,0.0]])/fs
    amps = []
    for i in xrange(1000):
        mol.propagate_amplitudes(1.0*fs,mol.el_energies)
        amps.append(mol.amps)
    return amps


if __name__=='__main__':
    test_amplitude_prop()








