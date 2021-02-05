class LangevinOneGene:
    del __init__(number_variables, params, meanVarFuncs, varFuncs, params_stochastic, \
    params_stochastic_strength, time_step, params_dosage_compensated = None, \
    dosage_compensation_factor = None):
        self.number_variables = number_variables
        self.params = params
        self.meanVarFunc = meanVarFunc
        self.params_stochastic = params_stochastic
        self.params_stochastic_strength = params_stochastic_strength
        self.time_step = time_step
        self.degradation_names = degradation_names
        self.params_dosage_compensated = params_dosage_compensated
        self.dosage_compensation_factor = dosage_compensation_factor

    def CLEComputeMeanVar(self, molecule_count):
        mean_, var_ = self.meanVarFunc(self.params, molecule_count)
        return mean_, var_

    def getProduction(mean_, var_):
        production = np.zeros(len(mean_))
        for _ in range(len(mean_)):
            brownian_noise = np.random.normal(size = 1)
            production[_] = mean_[_] + \
            np.sqrt(var_[_]) * brownian_noise[0]
        return production

    def updateCount(self, molecule_count_cur, production, degradation_name, \
    idx):
        molecule_count = 0
        brownian_noise = np.random.normal(size = 1)
        degradation_term = self.params[degradation_name] * \
        molecule_count_cur * self.time_step
        degradation_noise = np.sqrt(self.params[degradation_name] * \
        molecule_count_cur * time_step) * brownian_noise
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
    degradation_names):
        molecule_count = np.zeros(len(molecule_count_cur))
        production = self.getProduction(mean_, var_)
        for _ in range(len(molecule_count_cur)):
            molecule_count[_] = self.updateCount(self, molecule_count_cur[_], \
            production[_], degradation_names[_], _)

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

    def updateState(self, time_list_, molecule_count, molecule_count_cur, \
    mean_, var_):
        time_list.append(self.time_step)
        mean_, var_ = self.CLEComputeMeanVar(molecule_count_cur)
        molecule_count_cur = self.updateMoleculeCount(mean_, var_, \
        molecule_count_cur, self.degradation_names)
        for _ in range(len(molecule_count)):
                (molecule_count[_]).append(molecule_count_cur[_])

    def dosageCompensate(self, compensation_flag):
        if compensation_flag:
            for _ in range(len(self.params_dosage_compensated)):
                self.params[self.params_dosage_compensated[_]] *= \
                self.dosage_compensation_factor
        else:
            for _ in range(len(self.params_dosage_compensated)):
                self.params[self.params_dosage_compensated[_]] /= \
                self.dosage_compensation_factor

    def extractCurrentMoleculeCount(molecule_count, molecule_count_cur):
        for _ in range(len(molecule_count)):
            molecule_count_cur[_] = molecule_count[_][-1]

    def divideMolecules(molecule_count_cur):
        molecules_before_division = molecule_count_cur
        molecule_count_cur[:] = 0
        for _ in range(len(molecule_count_cur)):
            for i_ in range(int(molecules_before_division[_])):
                if np.random.randint(0, 2) == 1:
                    molecule_count_cur[_] += 1

    def simulateCellCycle(self, time_, molecule_count):
        time_list = [0]
        molecule_count_cur = np.zeros(len(molecule_count))
        cycle_num = 0
        while sum(time_list) < time_:
#             print('Time = ' + str(self.accumulation_flag))
#             self.sample_time_step()
            where_cell_cycle = self.checkWhereCellCycle(self, sum(time_list), \
            cycle_num)
            if where_cell_cycle == 0 or where_cell_cycle == 1:
                if where_cell_cycle == 0:
                    self.params["gene_copy"] = 2
                else:
                    self.params["gene_copy"] = 4
                self.extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
                self.updateState(time_list_, molecule_count, molecule_count_cur, \
                mean_, var_)
            elif where_cell_cycle == 2:
                self.params["gene_copy"] = 4
                self.dosageCompensate(compensation_flag = 1)
                self.extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
                self.updateState(time_list_, molecule_count, molecule_count_cur, \
                mean_, var_)
            else:
                self.params["gene_copy"] = 2
                self.dosageCompensate(compensation_flag = 0)
                self.extractCurrentMoleculeCount(molecule_count, molecule_count_cur)
                self.divideMolecules(molecule_count_cur)
                self.updateState(time_list_, molecule_count, molecule_count_cur, \
                mean_, var_)
                cycle_num += 1
        self.time_list = time_list
        self.molecule_count = molecule_count

class LangevinGRN:
    
