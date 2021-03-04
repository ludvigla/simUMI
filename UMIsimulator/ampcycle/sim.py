import collections
from collections import Counter
import numpy as np

class Pool(object):
    """
    A class to represent a pool of Unique Molecular Identifiers (UMIs).

    ...

    Attributes
    ----------
    umis : dict
        dictionary of input UMIs as keys and UMI counts as values
    umis_amplified : dict
        dictionary of amplified UMIs as keys and UMI counts as values
    efficiency_min, efficiency_max : float
        minimum and maximum efficiency of amplification
    efficiency_dict : dict
        dictionary of amplification efficiencies associated with each UMI

    Methods
    -------
    show()
        Prints information about the Pool class object
    createUMIsWithBias(n, efficiency_min=0.1, efficiency_max=0.9, error_rate=0.00001, UMI_length=7)
        Creates a new pool of efficiency biased UMIs
    SimulateCycleWithErrors(verbose=False)
        Runs a PCR cycle
    PCRcyclesWithErrorsAndBias(PCR_cycles)
        Runs a PCR amplification of UMIs for `PCR_cycles` number of cycles
    downsampleUMIs(efficiency=1e-2)
        Downsample an amplified UMI pool
    addSequencingErrors(seq_error_rate=0.01)
        Add sequencing errors to an amplififed UMI pool

    Examples
    ----------
    >>> pool = UMIsimulator.ampcycle.Pool()
    >>> pool.show()
    An uninitiated 'Pool' class object

    """

    def __init__(self):
        self.umis = None
        self.umis_amplified = None
        self.efficiency_min = None
        self.efficiency_max = None
        self.efficiency_dict = None

    def show(self):
        if (self.umis == None):
            print("An uninitiated 'Pool' class object")
        else:
            print("An object of class 'Pool' with a total of {} amplified UMIs ({} unique)".format(sum(self.umis_amplified.values()), len(self.umis_amplified.keys())))

    def createUMIsWithBias(self, n, efficiency_min=0.65, efficiency_max=0.9, error_rate=0.00001, UMI_length=7):
        """
        Create a pool of UMIs. 
        
        Initiates an object of class 'Pool' by creating a random 
        pool of UMIs to be used for simulation. An efficiency dictionary 
        is also created based on the minimum and maximum efficiencies
        provided.

        Parameters
        ----------
        n : int
            number of UMIs in the pool
        efficiency_min, efficiency_max : float
            minimum and maximum amplification efficiency
        error_rate : float
            mutation error rate
        UMI_length : int
            UMI length

        Returns
        ----------
        Pool
            an object of class 'Pool' with a randomized UMI pool 

        Examples
        ----------
        >>> pool = UMIsimulator.ampcycle.Pool()
        >>> pool.createUMIsWithBias(10, 0.2, 0.8)
        >>> pool.show()
        An object of class 'Pool' with a total of 10 amplified UMIs (10 unique)

        """
        assert efficiency_min > 0 and efficiency_max <= 1, "PCR efficiencies must be between 0 and 1"

        UMIs = []
        efficiency_dict = collections.defaultdict(float)
        for umi in range(0, n):
            UMI = ""
            for base in range(0, UMI_length):
                UMI += np.random.choice(["A", "C", "G", "T"])

            UMIs.append(UMI)
            efficiency_dict[UMI] = np.random.uniform(efficiency_min, efficiency_max)

        self.umis = Counter(UMIs)
        self.umis_amplified = Counter(UMIs)
        self.efficiency_min = efficiency_min
        self.efficiency_max = efficiency_max
        self.efficiency_dict = efficiency_dict
        self.error_rate = error_rate
        self.UMI_length = UMI_length

        return self

    def SimulateCycleWithErrors(self, verbose=False):
        """
        Simulate a PCR cycle

        Runs one round of PCR amplification where each UMI is amplified 
        according to it's associated amplification efficiency. The amplified
        UMIs then goes through a mutation process to emulate the errors 
        introduced by the DNA polymerase. DNA bases are randomly mutated at a 
        probability defined by the `error_rate` set when the `Pool` object was created.
        Both the input UMIs and amplified + mutated UMIs are returned.

        Returns
        ----------
        Pool
            an object of class 'Pool' with amplified UMIs

        Examples
        ----------
        >>> pool = UMIsimulator.ampcycle.Pool()
        >>> pool.createUMIsWithBias(10, 0.2, 0.8)
        >>> pool.SimulateCycleWithErrors()
        >>> pool.show()
        An object of class 'Pool' with a total of __ amplified UMIs (__ unique)

        """

        verboseprint = print if verbose else lambda *a, **k: None

        errors_dict = {"A": ["C", "G", "T"],
                       "C": ["A", "G", "T"],
                       "G": ["C", "A", "T"],
                       "T": ["C", "G", "A"]}

        ints = np.random.binomial(list(self.umis_amplified.values()), list(self.efficiency_dict.values()))
        c = []
        for x, y in zip(list(self.umis_amplified.keys()), ints):
            c.extend([x] * y)
        amplified_UMIs = Counter(c)

        if self.error_rate == 0:
            self.umis_amplified += amplified_UMIs
            return self
        else:
            # Sample from poisson distribution with mu = error_rate*len(UMI)*UMI_length
            sample_mutated = np.random.poisson(self.error_rate * sum(amplified_UMIs.values()) * self.UMI_length, 1)

            # mutated_UMIs = np.random.choice(UMI, sample_mutated, replace = False)
            p_vals_UMI = np.divide(list(amplified_UMIs.values()), sum(amplified_UMIs.values()))
            mutated_UMIs = np.random.choice(list(amplified_UMIs.keys()), size=sample_mutated, p=p_vals_UMI)

            mutated_UMIs_list = []
            for umi in mutated_UMIs:
                random_base = np.random.choice(range(len(umi)))
                umi = list(umi)
                umi[random_base] = np.random.choice(errors_dict[umi[random_base]])
                umi = "".join(umi)
                mutated_UMIs_list.append(umi)

            # Add mutated UMIs (no need to replace if len(pool) >> len(mutated))
            verboseprint("\tAdded %s UMIs to the pool\n\t%s were mutated" % (sum(amplified_UMIs.values()), len(mutated_UMIs_list)))
            self.umis_amplified += amplified_UMIs
            self.umis_amplified += Counter(mutated_UMIs_list)

            for umi in self.umis_amplified.keys():
                if not self.efficiency_dict[umi]:
                    self.efficiency_dict[umi] = np.random.uniform(self.efficiency_min, self.efficiency_max)

            return self


    def PCRcyclesWithErrorsAndBias(self, PCR_cycles, verbose=False):
        """
        Runs a PCR amplification of UMIs for `PCR_cycles` number of cycles on an object of class 'Pool'.
        
        The 'Pool' object needs to be initiated first by running `createUMIsWithBias`. 
        The input pool of UMIs is then amplified for each cycle based on the 
        amplification efficiency dictionary. For each cycle, UMIs can be mutated 
        by a rate specified as the 'error_rate' defined when the 'Pool' object was initiated.
        For each new UMI, an efficiency is randomly selected based on the 'efficiency_min' and
        'efficiency_max' parameters defined when the 'Pool' object was initiated.

        Parameters
        ----------
        PCR_cycles : int
            Number of PCR cycles

        Returns
        ----------
        Pool
            an object of class 'Pool' with amplified UMIs

        Examples
        ----------
        >>> pool = UMIsimulator.ampcycle.Pool()
        >>> pool.createUMIsWithBias(10, 0.2, 0.8)
        >>> pool.PCRcyclesWithErrorsAndBias(10)
        >>> pool.show()
        An object of class 'Pool' with a total of __ amplified UMIs (__ unique)

        """
        assert self.umis != None, "The UMI pool has not been initiated"

        verboseprint = print if verbose else lambda *a, **k: None

        cls = self.__class__
        for cycle in range(0, PCR_cycles):
            verboseprint("Running cycle %s:" % (cycle + 1))
            self = cls.SimulateCycleWithErrors(self, verbose)

        return self


    def downsampleUMIs(self, efficiency=1e-2):
        """
        subsample from UMI pool using an efficiency factor
        
        Parameters
        ----------
        efficiency : float
            proportion of UMIs to be kept from the pool

        Returns
        ----------
        Pool
            an object of class 'Pool' with a downsampled UMI pool

        Examples
        ----------
        >>> pool = UMIsimulator.ampcycle.Pool()
        >>> pool.createUMIsWithBias(10, 0.2, 0.8)
        >>> pool.PCRcyclesWithErrorsAndBias(10)
        >>> pool.downsampleUMIs()
        >>> pool.show()
        An object of class 'Pool' with a total of __ amplified UMIs (__ unique)
        
        """

        assert self.umis != None, "The UMI pool object has not been initiated"
        assert self.umis_amplified != None, "There are no amplified UMIs present in the pool"

        diluted_size = int(round(sum(self.umis_amplified.values()) * efficiency))
        p_vals_UMI = np.divide(list(self.umis_amplified.values()), sum(self.umis_amplified.values()))
        sampled_UMIs = np.random.choice(list(self.umis_amplified.keys()), size=diluted_size, p=p_vals_UMI)
        self.umis_amplified = Counter(sampled_UMIs)
        return self


    def addSequencingErrors(self, seq_error_rate=0.01):
        """

        Add sequencing errors to UMI pool
        
        Adds sequencing errors to the UMI pool as defined by the 
        `seq_error_rate` parameter. 

        Parameters
        ----------
        seq_error_rate : float
            sequencing error rate

        Returns
        ----------
        Pool
            an object of class 'Pool' with mutated UMIs

        Examples
        ----------
        >>> pool = UMIsimulator.ampcycle.Pool()
        >>> pool.createUMIsWithBias(10, 0.2, 0.8)
        >>> pool.PCRcyclesWithErrorsAndBias(10)
        >>> pool.addSequencingErrors()
        >>> pool.show()
        An object of class 'Pool' with a total of __ amplified UMIs (__ unique)

        """

        assert self.umis != None, "The UMI pool object has not been initiated"
        assert self.umis_amplified != None, "There are no amplified UMIs present in the pool"

        errors_dict = {"A": ["C", "G", "T"],
                       "C": ["A", "G", "T"],
                       "G": ["C", "A", "T"],
                       "T": ["C", "G", "A"]}
        umi_list_mutated = []
        for umi in self.umis_amplified.keys():
            umi_list = umi * self.umis_amplified[umi]
            mutated = [np.random.choice(errors_dict[L]) if np.random.random() < seq_error_rate else L for L in umi_list]
            for ls in zip(*[iter(mutated)]*self.UMI_length):
                umi_list_mutated += ["".join(ls)]

        self.umis_amplified = Counter(umi_list_mutated)

        return self

