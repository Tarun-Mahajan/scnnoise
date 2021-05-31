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
                 cell_cycle_simulate = 0, replication_time = None, division_time = None,\
                 dosage_compensation_factor = None, dosageCompensationFunc = None,\
                 params_stochastic = None, params_stochastic_strength = None, \
                 time_step_reduce_lim = 10):
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
        self.number_species = number_species
        self.cell_cycle_simulate = cell_cycle_simulate
        self.time_step_reduce_lim = time_step_reduce_lim

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

    def getNoiseParam(self, var_, time_step):
        production = np.zeros(len(var_))
        for _ in range(len(var_)):
            brownian_noise = np.random.normal(size = 1)
            production[_] = 1.0 * np.sqrt(var_[_] * time_step) * brownian_noise[0]
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

    def updateParam(self, param_noise_cur, production):
        param_noise = param_noise_cur + production
        return param_noise

    def updateMoleculeCount(self, mean_, var_, molecule_count_cur, \
    degradation_names, time_step):
        molecule_count = np.zeros(len(molecule_count_cur))
        production = self.getProduction(mean_, var_)
        for _ in range(len(molecule_count_cur)):
            molecule_count[_] = self.updateCount(molecule_count_cur[_], \
            production[_], degradation_names[_], _, time_step)

        return molecule_count

    def updateParamStochastic(self, param_noise_cur, time_step):
        param_noise = np.zeros(len(param_noise_cur))
        production = self.getNoiseParam(self.params_stochastic_strength, time_step)
        for _ in range(len(param_noise_cur)):
            param_noise[_] = self.updateParam(param_noise_cur[_], \
            production[_])

        return param_noise

    def checkWhereCellCycle(self, time_, time_step, cycle_num):
        where_cell_cycle = 4
        if self.cell_cycle_simulate:
            current_time_point = time_
            next_time_point = time_ + time_step
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

    def updateState(self, time_list, molecule_count, molecule_count_cur, check, where_cell_cycle, \
                   param_noise_list, param_noise_cur):
        def checkRestart(count, time_step_reduce_lim, molecule_count_cur, \
                        molecule_count_cur_update, param_noise_cur_update, param_noise_cur, \
                         time_step, time_step_orig):
            if count == time_step_reduce_lim:
                for _ in range(len(molecule_count_cur)):
                    if molecule_count_cur_update[_] <= 0:
                        molecule_count_cur[_] *= 1.5
                if param_noise_list is not None:
                    for _ in range(len(param_noise_cur)):
                        if param_noise_cur_update[_] <= 0:
                            param_noise_cur[_] *= 1.5
                count = 0
                time_step = 2 * time_step_orig
        time_step = 2 * self.time_step
        count_negative_flag = 0
        molecule_count_cur_update = 0 * molecule_count_cur
        if param_noise_list is not None:
            param_noise_cur_update = 0 * param_noise_cur
        else:
            param_noise_cur_update = param_noise_cur
        count = 0
        while not count_negative_flag:
            checkRestart(count, self.time_step_reduce_lim, molecule_count_cur, \
                         molecule_count_cur_update, param_noise_cur_update, \
                         param_noise_cur, time_step, self.time_step)
#             print('Count = ' + str(count) + ' Product1 = ' + str(np.prod(molecule_count_cur_update)) + \
#                  ' Product2 = ' + str(np.prod(param_noise_cur_update)))

            time_step /= 2
            if param_noise_list is not None:
                for _ in range(len(param_noise_cur)):
                    self.params[self.params_stochastic[_]] = param_noise_cur[_]
#                     print('Params = ' + str(self.params[self.params_stochastic[_]]))
            mean_, var_ = self.CLEComputeMeanVar(molecule_count_cur, time_step)
            molecule_count_cur_update = self.updateMoleculeCount(mean_, var_, \
            molecule_count_cur, self.degradation_names, time_step)
            if param_noise_list is not None:
                param_noise_cur_update = self.updateParamStochastic(param_noise_cur, time_step)
            if np.prod(molecule_count_cur_update) > 0:
                if param_noise_list is None:
                    count_negative_flag = 1
                else:
                    if np.prod(param_noise_cur_update) > 0:
                        count_negative_flag = 1
            count += 1
        if not check:
            for _ in range(len(molecule_count)):
                (molecule_count[_]).append(molecule_count_cur_update[_])
            if param_noise_list is not None:
                for _ in range(len(param_noise_list)):
                    (param_noise_list[_]).append(param_noise_cur_update[_])
            if where_cell_cycle != 2 and where_cell_cycle != 3:
                time_list.append(time_step)
            else:
                time_list.append(self.time_step)
        else:
            return time_step

    def dosageCompensate(self, compensation_flag):
        if self.dosageCompensationFunc is not None:
            self.dosageCompensationFunc(self, compensation_flag)

    def divideMolecules(self, molecule_count_cur, where_cell_cycle):
        if where_cell_cycle == 3:
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

    def checkBoundary(self, where_cell_cycle, cycle_num, time_list, molecule_count):
        molecule_count_cur = np.zeros(self.number_species)
        if where_cell_cycle == 2:
            self.params["gene_copy"] *= 2
            self.dosageCompensate(compensation_flag = 1)
            self.extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
            time_step = self.updateState(time_list, molecule_count, molecule_count_cur, 1, \
                                         where_cell_cycle)
            where_cell_cycle = self.checkWhereCellCycle(sum(time_list), time_step, \
                                                        cycle_num)
            self.params["gene_copy"] /= 2
            self.dosageCompensate(compensation_flag = 0)
            return where_cell_cycle
        elif where_cell_cycle == 3:
            self.params["gene_copy"] /= 2
            self.dosageCompensate(compensation_flag = 0)
            self.extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
            molecule_count_cur = self.divideMolecules(molecule_count_cur)
            time_step = self.updateState(time_list, molecule_count, molecule_count_cur, 1, \
                                        where_cell_cycle)
            where_cell_cycle = self.checkWhereCellCycle(sum(time_list), time_step, \
                                                        cycle_num)
            self.params["gene_copy"] *= 2
            self.dosageCompensate(compensation_flag = 1)
            return where_cell_cycle
        else:
            return where_cell_cycle
    def setGeneCopy(self, gene_copy, where_cell_cycle):
        if where_cell_cycle == 0 or where_cell_cycle == 3 or where_cell_cycle == 4:
            self.params["gene_copy"]  = gene_copy
        else:
            self.params["gene_copy"] = 2 * gene_copy

    def replicationChanges(self, where_cell_cycle):
        if where_cell_cycle == 2:
            self.dosageCompensate(compensation_flag = 1)
        elif where_cell_cycle == 3:
            self.dosageCompensate(compensation_flag = 0)
        else:
            dummy_ = None

    def simulateCellCycle(self, time_, molecule_count):
        def updateCellCycleNum(where_cell_cycle, cycle_num):
            if where_cell_cycle == 3:
                cycle_num[0] += 1

        def extractCurrentMoleculeCount(molecule_count, molecule_count_cur):
            for _ in range(len(molecule_count)):
                molecule_count_cur[_] = molecule_count[_][-1]

        def extractCurrentParamStochastic(param_noise_list, param_noise_cur):
            for _ in range(len(param_noise_list)):
                param_noise_cur[_] = param_noise_list[_][-1]

        time_list = [0]
        molecule_count_cur = np.zeros(len(molecule_count))

        if self.params_stochastic is not None:
            param_noise_cur = np.zeros(len(self.params_stochastic))
            param_noise_list = [None] * len(self.params_stochastic)
            for _ in range(len(self.params_stochastic)):
                param_noise_list[_] = [self.params[self.params_stochastic[_]]]
        else:
            param_noise_cur = None
            param_noise_list = None
        cycle_num = [0]
        gene_copy = np.copy(self.params["gene_copy"])
#         print('gene copy = ' + str(gene_copy))
        if not self.cell_cycle_simulate:
            where_cell_cycle = 4
        else:
            where_cell_cycle = 0
        while sum(time_list) < time_:
            #             print('Time = ' + str(sum(time_list)))

#             self.sample_time_step()
#             print('time step = ' + str(self.time_step))
            where_cell_cycle = self.checkWhereCellCycle(sum(time_list), self.time_step, \
                                                        cycle_num[0])
#             print('Time = ' + str(sum(time_list)) + ' mRNA_n = ' + str(molecule_count[0][-1]) +\
#                   ' mRNA_m = ' + str(molecule_count[1][-1]) + ' where = ' + str(where_cell_cycle) + \
#                   ' cycle = ' + str(cycle_num))
            self.setGeneCopy(gene_copy, where_cell_cycle)
            self.replicationChanges(where_cell_cycle)
            extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
            if self.params_stochastic is not None:
                extractCurrentParamStochastic(param_noise_list, param_noise_cur)
            molecule_count_cur = self.divideMolecules(molecule_count_cur, where_cell_cycle)
            self.updateState(time_list, molecule_count, molecule_count_cur, 0, where_cell_cycle, \
                            param_noise_list, param_noise_cur)
            updateCellCycleNum(where_cell_cycle, cycle_num)
        self.time_list = time_list
        self.molecule_count = molecule_count
        self.param_noise_list = param_noise_list
