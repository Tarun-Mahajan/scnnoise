import numpy as np
import pandas as pd
from plotnine import *
import multiprocessing as mp
import itertools
import random
from scipy.interpolate import interp1d
class LangevinOneGene:
    """
    Chemical Langevin Equation (CLE) implementation for one gene using model M5 [1]

    Attributes
    ----------
    number_species : int
        Number of distinct molecular species to simulate
    params : dict of str: float
        A dictionary containing all the parameters for the gene expression model
    meanVarFunc : function(dict of str: float, float, ndarray)
        A function to calculate mean and variance for the chemical langevin equation
    time_step : float
        Time step to progress Euler-Maruyama iterations
    degradation_names : list of str
        List containing names of degradation parameters for the chemical species in the system
    replication_time : float
        Time at which replication happens
    division_time : float
        Time at which cell division happens
    dosage_compensation_factor : float
        Factor by which burst frequency decreases at replication
    dosageCompensationFunc : function(cls, bool)
        A function to update parameters at replication
    params_stochastic : list of bool
        List which states whether parameters are stochastic or not
    params_stochastic_strength : ndarray
        1D array giving noise strength for the system parameters if these are stochastic

    Notes
    -----
    This class implements the M5 model for gene expression introduced in [1].
    Brielfy, under this model gene expression is modeled by a two-state process,
    where the two states represent the active and inactive states for the
    promoter. Further, cell-cycle variation is also included. The model accounts
    for both nascent and mature mRNA. Here, we simulate this model using the
    CLE formulation. We adapt the bursty-langevin approximation developed in [2]
    to the M5 model.

    .. [1] Cao, Zhixing, and Ramon Grima. "Analytical distributions for detailed
    models of stochastic gene expression in eukaryotic cells." Proceedings of
    the National Academy of Sciences 117.9 (2020): 4682-4692.
    .. [2] Yan, Ching-Cher Sanders, et al. "Efficient and flexible implementation
    of Langevin simulation for gene burst production." Scientific reports 7.1
    (2017): 1-16.
    """
    def __init__(self, number_species, params, meanVarFunc, time_step, degradation_names, \
                 replication_time, division_time,\
                 dosage_compensation_factor = None, dosageCompensationFunc = None,\
                 params_stochastic = None, params_stochastic_strength = None):
        self.params = params
        self.meanVarFunc = meanVarFunc
        self.params_stochastic = params_stochastic
        self.params_stochastic_strength = params_stochastic_strength
        self.time_step = time_step
        self.degradation_names = degradation_names
        self.dosage_compensation_factor = dosage_compensation_factor
        self.dosageCompensationFunc = dosageCompensationFunc
        self.division_time = division_time
        self.replication_time = replication_time
        self.production_accumulated = np.zeros(number_species)
        self.accumulation_flag = [0, 0]

    def CLEComputeMeanVar(self, molecule_count, time_step):
        mean_, var_ = self.meanVarFunc(self.params, time_step, molecule_count)
        return mean_, var_

    def getProduction(self, mean_, var_):
        production = np.zeros(len(mean_))
        for _ in range(len(mean_)):
            brownian_noise = np.random.normal(size = 1)
            production[_] = mean_[_] + \
            1.0 * np.sqrt(var_[_]) * brownian_noise[0]
        return production

    def updateCount(self, molecule_count_cur, production, degradation_name, \
    idx, time_step):
        _ = idx
        molecule_count = molecule_count_cur
        brownian_noise = np.random.normal(size = 1)
        degradation_term = self.params[degradation_name] * \
        molecule_count_cur * time_step
        degradation_noise = np.sqrt(self.params[degradation_name] * \
        molecule_count_cur * time_step) * brownian_noise[0]
        molecule_count += - degradation_term - degradation_noise
        if self.accumulation_flag[_]:
            self.production_accumulated[_] += production
            if self.production_accumulated[_] > 0:
                molecule_count += self.production_accumulated[_]
                self.accumulation_flag[_] = 0
                self.production_accumulated[_] = 0
        else:
            if production > 0:
                molecule_count += production
            else:
                self.accumulation_flag[_] = 1
                self.production_accumulated[_] = production
        return molecule_count

    def updateMoleculeCount(self, mean_, var_, molecule_count_cur, \
    degradation_names, time_step):
        molecule_count = np.zeros(len(molecule_count_cur))
        production = self.getProduction(mean_, var_)
        for _ in range(len(molecule_count_cur)):
            molecule_count[_] = self.updateCount(molecule_count_cur[_], \
            production[_], degradation_names[_], _, time_step)

        return molecule_count

    def checkWhereCellCycle(self, time_, cycle_num):
        current_time_point = time_
        next_time_point = time_ + self.time_step
        cycle_begin_time = cycle_num * self.division_time
        cycle_end_time = (cycle_num + 1) * self.division_time
        cycle_replication_time = cycle_num * self.division_time + \
        self.replication_time
        where_cell_cycle = 3
        if (next_time_point < cycle_replication_time and \
        current_time_point >= cycle_begin_time):
            where_cell_cycle = 0
        elif (next_time_point < cycle_end_time and \
        current_time_point >= cycle_replication_time):
            where_cell_cycle = 1
        elif next_time_point < cycle_end_time and \
        next_time_point > cycle_replication_time and \
        current_time_point < cycle_replication_time:
            where_cell_cycle = 2
        else:
            where_cell_cycle = 3
        return where_cell_cycle

    def updateState(self, time_list, molecule_count, molecule_count_cur):
        time_list.append(self.time_step)
        time_step = 2 * self.time_step
        count_negative_flag = 0
        molecule_count_cur_update = 0 * molecule_count_cur
        count = 0
        while not count_negative_flag:
            print('Count = ' + str(count) + ' Product = ' + str(np.prod(molecule_count_cur_update)))
            time_step /= 2
            mean_, var_ = self.CLEComputeMeanVar(molecule_count_cur, time_step)
            molecule_count_cur_update = self.updateMoleculeCount(mean_, var_, \
            molecule_count_cur, self.degradation_names, time_step)
            if np.prod(molecule_count_cur_update) > 0:
                count_negative_flag = 1
            count += 1
        for _ in range(len(molecule_count)):
                (molecule_count[_]).append(molecule_count_cur_update[_])

    def dosageCompensate(self, compensation_flag):
        if self.dosageCompensationFunc is not None:
            self.dosageCompensationFunc(self, compensation_flag)

    def extractCurrentMoleculeCount(self, molecule_count, molecule_count_cur):
        for _ in range(len(molecule_count)):
            molecule_count_cur[_] = molecule_count[_][-1]

    def divideMolecules(self, molecule_count_cur):
        molecules_before_division = np.copy(molecule_count_cur)
        molecule_count_cur[:] = 0
        for _ in range(len(molecule_count_cur)):
            if _ == 0:
                molecule_count_cur[_] = molecules_before_division[_]/2
            else:
                for i_ in range(int(molecules_before_division[_])):
                    if np.random.randint(0, 2) == 1:
                        molecule_count_cur[_] += 1
        return molecule_count_cur

    def simulateCellCycle(self, time_, molecule_count):
        time_list = [0]
        molecule_count_cur = np.zeros(len(molecule_count))
        cycle_num = 0
        while sum(time_list) < time_:
#             self.sample_time_step()
            where_cell_cycle = self.checkWhereCellCycle(sum(time_list), \
            cycle_num)
#             print('Time = ' + str(sum(time_list)))
            print('Time = ' + str(sum(time_list)) + ' mRNA_n = ' + str(molecule_count[0][-1]) +\
                 ' mRNA_m = ' + str(molecule_count[1][-1]) + ' where = ' + str(where_cell_cycle))
            if where_cell_cycle == 0 or where_cell_cycle == 1:
                if where_cell_cycle == 0:
                    self.params["gene_copy"] = 2
                else:
                    self.params["gene_copy"] = 4
                self.extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
                self.updateState(time_list, molecule_count, molecule_count_cur)
            elif where_cell_cycle == 2:
                self.params["gene_copy"] = 4
                self.dosageCompensate(compensation_flag = 1)
                self.extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
                self.updateState(time_list, molecule_count, molecule_count_cur)
            else:
                self.params["gene_copy"] = 2
                self.dosageCompensate(compensation_flag = 0)
                self.extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
                molecule_count_cur = self.divideMolecules(molecule_count_cur)
                self.updateState(time_list, molecule_count, molecule_count_cur)
                cycle_num += 1
        self.time_list = time_list
        self.molecule_count = molecule_count
