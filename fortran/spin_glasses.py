import _spin_glasses
import f90wrap.runtime
import logging

class Precision(f90wrap.runtime.FortranModule):
    """
    Module precision
    
    
    Defined at precision.f90 lines 1-7
    
    """
    @property
    def sl(self):
        """
        Element sl ftype=integer pytype=int
        
        
        Defined at precision.f90 line 5
        
        """
        return _spin_glasses.f90wrap_precision__get__sl()
    
    @property
    def dl(self):
        """
        Element dl ftype=integer pytype=int
        
        
        Defined at precision.f90 line 6
        
        """
        return _spin_glasses.f90wrap_precision__get__dl()
    
    def __str__(self):
        ret = ['<precision>{\n']
        ret.append('    sl : ')
        ret.append(repr(self.sl))
        ret.append(',\n    dl : ')
        ret.append(repr(self.dl))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

precision = Precision()

class Montecarlo(f90wrap.runtime.FortranModule):
    """
    Module montecarlo
    
    
    Defined at montecarlo.f90 lines 1-267
    
    """
    @staticmethod
    def initializespins(n, s):
        """
        initializespins(n, s)
        
        
        Defined at montecarlo.f90 lines 18-35
        
        Parameters
        ----------
        n : int
        s : int array
        
        """
        _spin_glasses.f90wrap_initializespins(n=n, s=s)
    
    @staticmethod
    def metropolisstep(n, s, j, h, beta, deltaenergy):
        """
        metropolisstep(n, s, j, h, beta, deltaenergy)
        
        
        Defined at montecarlo.f90 lines 40-92
        
        Parameters
        ----------
        n : int
        s : int array
        j : float array
        h : float array
        beta : float
        deltaenergy : float
        
        -------------------------------------------------
           the deltaEnergy subroutine is structured in this way
           subroutine deltaEnergy(n, s1, s2, J, h, dE)
               implicit none
               integer, intent(in) :: n
               integer, intent(in) :: s1(n), s2(n)
               real(dl), intent(in) :: J(n,n), h(n)
               real(dl), intent(out) :: dE
                calculate dE
           end subroutine deltaEnergy
        -------------------------------------------------
        > perform n updates for each Metropolis-Hastings step
        """
        _spin_glasses.f90wrap_metropolisstep(n=n, s=s, j=j, h=h, beta=beta, \
            deltaenergy=deltaenergy)
    
    @staticmethod
    def swapconfigurations(n, s1, s2, beta1, beta2, j, h, deltaenergy):
        """
        swapconfigurations(n, s1, s2, beta1, beta2, j, h, deltaenergy)
        
        
        Defined at montecarlo.f90 lines 96-122
        
        Parameters
        ----------
        n : int
        s1 : int array
        s2 : int array
        beta1 : float
        beta2 : float
        j : float array
        h : float array
        deltaenergy : float
        
        """
        _spin_glasses.f90wrap_swapconfigurations(n=n, s1=s1, s2=s2, beta1=beta1, \
            beta2=beta2, j=j, h=h, deltaenergy=deltaenergy)
    
    @staticmethod
    def simulatedannealing(n_beta, betas, n_equil, j, h, n, s, deltaenergy):
        """
        simulatedannealing(n_beta, betas, n_equil, j, h, n, s, deltaenergy)
        
        
        Defined at montecarlo.f90 lines 126-165
        
        Parameters
        ----------
        n_beta : int
        betas : float array
        n_equil : int
        j : float array
        h : float array
        n : int
        s : int array
        deltaenergy : float
        
        """
        _spin_glasses.f90wrap_simulatedannealing(n_beta=n_beta, betas=betas, \
            n_equil=n_equil, j=j, h=h, n=n, s=s, deltaenergy=deltaenergy)
    
    @staticmethod
    def simulatedannealingsampling(n_beta, betas, n_equil, j, h, n, n_samples, \
        spins, deltaenergy):
        """
        simulatedannealingsampling(n_beta, betas, n_equil, j, h, n, n_samples, spins, \
            deltaenergy)
        
        
        Defined at montecarlo.f90 lines 170-194
        
        Parameters
        ----------
        n_beta : int
        betas : float array
        n_equil : int
        j : float array
        h : float array
        n : int
        n_samples : int
        spins : int array
        deltaenergy : float
        
        """
        _spin_glasses.f90wrap_simulatedannealingsampling(n_beta=n_beta, betas=betas, \
            n_equil=n_equil, j=j, h=h, n=n, n_samples=n_samples, spins=spins, \
            deltaenergy=deltaenergy)
    
    @staticmethod
    def paralleltempering(n, n_replicas, betas, n_sweeps, n_burn_in, n_equil, j, h, \
        s, deltaenergy):
        """
        paralleltempering(n, n_replicas, betas, n_sweeps, n_burn_in, n_equil, j, h, s, \
            deltaenergy)
        
        
        Defined at montecarlo.f90 lines 197-267
        
        Parameters
        ----------
        n : int
        n_replicas : int
        betas : float array
        n_sweeps : int
        n_burn_in : int
        n_equil : int
        j : float array
        h : float array
        s : int array
        deltaenergy : float
        
        """
        _spin_glasses.f90wrap_paralleltempering(n=n, n_replicas=n_replicas, betas=betas, \
            n_sweeps=n_sweeps, n_burn_in=n_burn_in, n_equil=n_equil, j=j, h=h, s=s, \
            deltaenergy=deltaenergy)
    
    _dt_array_initialisers = []
    

montecarlo = Montecarlo()

class Graphs(f90wrap.runtime.FortranModule):
    """
    Module graphs
    
    
    Defined at graphs.f90 lines 1-172
    
    """
    @staticmethod
    def generatesquare2dlattice(l, jconn):
        """
        generatesquare2dlattice(l, jconn)
        
        
        Defined at graphs.f90 lines 17-84
        
        Parameters
        ----------
        l : int
        jconn : int array
        
        """
        _spin_glasses.f90wrap_generatesquare2dlattice(l=l, jconn=jconn)
    
    @staticmethod
    def generatesquarechimeragraph(l, jconn):
        """
        generatesquarechimeragraph(l, jconn)
        
        
        Defined at graphs.f90 lines 88-172
        
        Parameters
        ----------
        l : int
        jconn : int array
        
        """
        _spin_glasses.f90wrap_generatesquarechimeragraph(l=l, jconn=jconn)
    
    _dt_array_initialisers = []
    

graphs = Graphs()

class Problems(f90wrap.runtime.FortranModule):
    """
    Module problems
    
    
    Defined at problems.f90 lines 1-245
    
    """
    @staticmethod
    def henfrustratedclusterloop(n, alpha, r, jconn, jproblem):
        """
        henfrustratedclusterloop(n, alpha, r, jconn, jproblem)
        
        
        Defined at problems.f90 lines 21-167
        
        Parameters
        ----------
        n : int
        alpha : float
        r : float
        jconn : int array
        jproblem : float array
        
        """
        _spin_glasses.f90wrap_henfrustratedclusterloop(n=n, alpha=alpha, r=r, \
            jconn=jconn, jproblem=jproblem)
    
    @staticmethod
    def kingfrustratedclusterloop():
        """
        kingfrustratedclusterloop()
        
        
        Defined at problems.f90 lines 174-177
        
        
        """
        _spin_glasses.f90wrap_kingfrustratedclusterloop()
    
    @staticmethod
    def deceptiveclusterloop():
        """
        deceptiveclusterloop()
        
        
        Defined at problems.f90 lines 183-185
        
        
        """
        _spin_glasses.f90wrap_deceptiveclusterloop()
    
    _dt_array_initialisers = []
    

problems = Problems()

class Model(f90wrap.runtime.FortranModule):
    """
    Module model
    
    
    Defined at model.f90 lines 1-126
    
    """
    @staticmethod
    def energy(n, s, h, j):
        """
        energy = energy(n, s, h, j)
        
        
        Defined at model.f90 lines 18-37
        
        Parameters
        ----------
        n : int
        s : int array
        h : float array
        j : float array
        
        Returns
        -------
        energy : float
        
        """
        energy = _spin_glasses.f90wrap_energy(n=n, s=s, h=h, j=j)
        return energy
    
    @staticmethod
    def delta_energy(n, s1, s2, j, h):
        """
        de = delta_energy(n, s1, s2, j, h)
        
        
        Defined at model.f90 lines 42-54
        
        Parameters
        ----------
        n : int
        s1 : int array
        s2 : int array
        j : float array
        h : float array
        
        Returns
        -------
        de : float
        
        """
        de = _spin_glasses.f90wrap_delta_energy(n=n, s1=s1, s2=s2, j=j, h=h)
        return de
    
    @staticmethod
    def makelineartemperatureschedule(n_beta, betas, t_min, t_max):
        """
        makelineartemperatureschedule(n_beta, betas, t_min, t_max)
        
        
        Defined at model.f90 lines 58-71
        
        Parameters
        ----------
        n_beta : int
        betas : float array
        t_min : float
        t_max : float
        
        """
        _spin_glasses.f90wrap_makelineartemperatureschedule(n_beta=n_beta, betas=betas, \
            t_min=t_min, t_max=t_max)
    
    @staticmethod
    def makeexponentialtemperatureshedulemaxmin(n_beta, betas, exponent, t_max):
        """
        makeexponentialtemperatureshedulemaxmin(n_beta, betas, exponent, t_max)
        
        
        Defined at model.f90 lines 75-98
        
        Parameters
        ----------
        n_beta : int
        betas : float array
        exponent : float
        t_max : float
        
        """
        _spin_glasses.f90wrap_makeexponentialtemperatureshedulemaxmin(n_beta=n_beta, \
            betas=betas, exponent=exponent, t_max=t_max)
    
    @staticmethod
    def makeexponentialtemperaturesheduleminmax(n_beta, betas, exponent, t_min):
        """
        makeexponentialtemperaturesheduleminmax(n_beta, betas, exponent, t_min)
        
        
        Defined at model.f90 lines 102-126
        
        Parameters
        ----------
        n_beta : int
        betas : float array
        exponent : float
        t_min : float
        
        """
        _spin_glasses.f90wrap_makeexponentialtemperaturesheduleminmax(n_beta=n_beta, \
            betas=betas, exponent=exponent, t_min=t_min)
    
    _dt_array_initialisers = []
    

model = Model()

