import multiprocessing as mp
import numpy as np
import pandas as pd
from UMIsimulator.dedup import dedup_percentile,dedup_naive,dedup_unique,dedup_hierarchical,dedup_graph,dedup_adj,dedup_dir_adj
import UMIsimulator.ampcycle

class SimUMI(object):
    """
    A class to represent an object to be used for simulations for
    PCR amplification of UMIs.

    ...

    Attributes
    ----------
    iterations : int
        number of iterations to be run in the simulation
    umi_length : int
        length of UMI
    pcr_cycles : int
        number of pcr cycles to run
    dna_pol_error_rate : float
        error rate for the DNA polymerase used in the PCR reaction
    seq_pol_error_rate : float
        error rate for sequencing
    eff_min, eff_max : float
        minimum and maximum efficiency of amplification
    number_of_umis : int
        number of UMIs to start with in the pool

    Methods
    -------
    updateVariable(variable=None)
        Update variable to vary in simulation, all other variables are fixed
    iterator(variable, iter_list, dedup_methods, verbose=False)
        Create an iterator for umi simulation

    Examples
    ----------
    >>> pool = UMIsimulator.ampcycle.Pool()
    >>> pool.show()
    An uninitiated 'Pool' class object

    """

    def __init__(self,
                 iterations=100,
                 umi_length=7,
                 pcr_cycles=10,
                 dna_pol_error_rate=1e-5,
                 seq_error_rate=1e-3,
                 eff_min=0.65,
                 eff_max=0.9,
                 number_of_umis=20):

        self.iterations = iterations
        self.umi_length = umi_length
        self.pcr_cycles = pcr_cycles
        self.dna_pol_error_rate = dna_pol_error_rate
        self.seq_error_rate = seq_error_rate
        self.eff_min = eff_min
        self.eff_max = eff_max
        self.number_of_umis = number_of_umis
        self.variable = None

    def updateVariable(self, variable=None):

        if variable == "umi_length":
            self.umi_length = self.variable
        elif variable == "rt_error_rate":
            self.rt_error_rate = self.variable
        elif variable == "pcr_cycles":
            self.pcr_cycles = self.variable
        elif variable == "dna_pol_error_rate":
            self.dna_pol_error_rate = self.variable
        elif variable == "seq_error_rate":
            self.seq_error_rate = self.variable
        elif variable == "eff_min":
            self.eff_min = self.variable
        elif variable == "eff_max":
            self.eff_max = self.variable
        elif variable == "number_of_umis":
            self.number_of_umis = self.variable
        else:
            raise ValueError("Not a valid input parameter")

    def iterator(self, variable, iter_list, dedup_methods, verbose=False):

        verboseprint = print if verbose else lambda *a, **k: None

        df = pd.DataFrame()

        for value in iter_list:
            self.variable = value
            self.updateVariable(variable)

            tmp_df = simulationCycle(self.iterations,
                                     self.number_of_umis,
                                     self.umi_length,
                                     self.pcr_cycles,
                                     self.dna_pol_error_rate,
                                     self.seq_error_rate,
                                     self.eff_min,
                                     self.eff_max,
                                     dedup_methods)
            tmp_df[variable] = value

            verboseprint("Finished %s for iterations variable '%s' set to %s" % (self.iterations, variable, value))

            df = df.append(tmp_df)

        return df


def mp_iter(
    number_of_UMIs, 
    length_of_UMIs, 
    PCR_cycles,
    error_rate, 
    seq_error_rate, 
    efficiency_min,
    efficiency_max,
    efficiency = None,
    seed = None):
    """
    
    Create a pool of amplified UMIs

    Parameters
    ----------
    number_of_UMIs : int
        number of UMIs to start with
    length_of_UMIs : int
        UMI length
    PCR_cycles : int
        number of PCR cycles
    error_rate : float
        error rate for DNA polymerase
    seq_error_rate : float
        sequencing error rate
    efficiency_min, efficiency_max : float
        minimum and maximum amplification efficiency
    efficiency : float
        downsample proportion

    Returns
    ----------
    Pool
        an object of class 'Pool' with an amplified UMI pool 

    """

    if seed != None:
        np.random.seed(seed)
    pool = UMIsimulator.ampcycle.Pool()
    pool.createUMIsWithBias(number_of_UMIs, efficiency_min, efficiency_max, error_rate, length_of_UMIs)
    pool.PCRcyclesWithErrorsAndBias(PCR_cycles)
    if efficiency != None:
        pool.downsampleUMIs(efficiency)
    pool.addSequencingErrors(seq_error_rate)
    return pool


def CountsPerMethod(
    iterations, 
    number_of_UMIs, 
    length_of_UMIs, 
    PCR_cycles,
    error_rate, 
    seq_error_rate, 
    efficiency_min, 
    efficiency_max,
    dedup_methods,
    efficiency = None,
    seed = None):
    """
    
    Estimate UMI counts 

    Runs a simulation process and using multiprocessing and to 
    generate a set of amplified UMI pools. Then estimates UMI counts 
    from and amplified UMI pool using a set of predefined duplicated 
    removal methods.

    Parameters
    ----------
    iterations : int
        number of iterations to run
    number_of_UMIs : int
        number of UMIs to start with
    length_of_UMIs : int
        UMI length
    PCR_cycles : int
        number of PCR cycles
    error_rate : float
        error rate for DNA polymerase
    seq_error_rate : float
        sequencing error rate
    efficiency_min, efficiency_max : float
        minimum and maximum amplification efficiency
    dedup_methods : list
        list of duplicate removal methods to test
    efficiency : float
        downsample proportion

    Returns
    ----------
    dict
        a dictionary with siulation results for each dupicate removal method

    """

    AssertDedups(dedup_methods)

    process_pool = mp.Pool(processes=7)
    results = [process_pool.apply_async(mp_iter, args = (
    number_of_UMIs, length_of_UMIs, PCR_cycles, error_rate, 
    seq_error_rate, efficiency_min, efficiency_max, efficiency, seed)) for i in range(iterations)]

    naive = []; simple = []; complete = []; ward = []
    unique = []; perc = []; graph = []; adj = []; dir_adj = []

    for p in results:

        UMIs_pool = p.get()

        for dedup in dedup_methods:
            if dedup == "naive":
                naive.append(dedup_naive(UMIs_pool.umis_amplified))
            elif dedup == "hier_simple":
                simple.append(dedup_hierarchical(UMIs_pool.umis_amplified, method = "single"))
            elif dedup == "hier_complete":
                complete.append(dedup_hierarchical(UMIs_pool.umis_amplified, method = "complete"))
            elif dedup == "hier_ward":
                ward.append(dedup_hierarchical(UMIs_pool.umis_amplified, method = "ward"))
            elif dedup == "unique":
                unique.append(dedup_unique(UMIs_pool.umis_amplified))
            elif dedup == "perc":
                perc.append(dedup_percentile(UMIs_pool.umis_amplified))
            elif dedup == "graph":
                graph.append(dedup_graph(UMIs_pool.umis_amplified))
            elif dedup == "adj":
                adj.append(dedup_adj(UMIs_pool.umis_amplified))
            elif dedup == "dir_adj":
                dir_adj.append(dedup_dir_adj(UMIs_pool.umis_amplified))

    process_pool.terminate()

    d = dict(zip(["naive", "hier_simple", "hier_complete", "hier_ward", "unique", "perc", "graph", "adj", "dir_adj"], 
        [naive, simple, complete, ward, unique, perc, graph, adj, dir_adj]))
    d = {k: v for k, v in d.items() if len(v) > 0}

    return d


def simulationCycle(
    iterations, 
    number_of_UMIs, 
    length_of_UMIs, 
    PCR_cycles, 
    error_rate, 
    seq_error_rate, 
    efficiency_min, 
    efficiency_max,
    dedup_methods,
    efficiency = None,
    seed=None):
    """
    
    Summarize results for a simulation cycle 

    Runs simulation process and converts results into a DataFrame.

    Parameters
    ----------
    iterations : int
        number of iterations to run
    number_of_UMIs : int
        number of UMIs to start with
    length_of_UMIs : int
        UMI length
    PCR_cycles : int
        number of PCR cycles
    error_rate : float
        error rate for DNA polymerase
    seq_error_rate : float
        sequencing error rate
    efficiency_min, efficiency_max : float
        minimum and maximum amplification efficiency
    dedup_methods : list
        list of duplicate removal methods to test
    efficiency : float
        downsample proportion

    Returns
    ----------
    DataFrame
        a pandas DataFrame with four columns
        "number_of_molecules" = Number of UMIs input to the simulation
        "dedup" = duplicate removal method
        "count" = estimated UMI count
        "CV" = coefficient of variation calculated for estimated counts across all iterations

    """

    dedups = ["naive", "hier_simple", "hier_complete", "hier_ward", 
                     "unique", "perc", "graph", "adj", "dir_adj"]
    for dedup in dedup_methods:
        assert dedup in dedups, "'%s' is not a valid duplicate removal method" % dedup

    counts = CountsPerMethod(iterations, number_of_UMIs, length_of_UMIs, PCR_cycles, 
                            error_rate, seq_error_rate, efficiency_min, efficiency_max, dedup_methods, efficiency, seed)

    tmp_df = pd.DataFrame({'number_of_molecules': (number_of_UMIs,) * len(dedup_methods),
                           'dedup': dedup_methods,
                           'count': [np.mean(counts[x]) for x in dedup_methods],
                           'CV': [CV(counts[x]) for x in dedup_methods]})
    return tmp_df


def CV(x):
    """

    calculate coefficient of variation

    Parameters
    ----------
    x : array
        an array like object

    Returns
    ----------
    float
        coefficient of variation

    """
    return np.std(x) / np.mean(x)


def AssertDedups(dedup_methods):
    dedups = ["naive", "hier_simple", "hier_complete", "hier_ward", 
                     "unique", "perc", "graph", "adj", "dir_adj"]
    for dedup in dedup_methods:
        assert dedup in dedups, "'%s' is not a valid duplicate removal method" % dedup
