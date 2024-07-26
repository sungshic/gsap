import random
from bisect import bisect
from operator import itemgetter

import numpy as np
import numpy.ma as ma
import itertools

__all__ = [
    "SNPSetManager",
    "PopulationManager",
    "PermPatternGenerator",
]

class SNPSetManager:
    _snp_set = None
    _snp_lut = None
    _max_snp_idx = None

    def __init__(self, max_snp_idx):
        self._snp_set = set() # initialize
        self._snp_lut = {} # dict initialize
        self._max_snp_idx = max_snp_idx

    # snp_list is a list of tuples (snp_idx, snp_details)
    def addSNPs(self, snp_list):
        for snp_idx, snp_details in snp_list:
            self._snp_set.add(snp_idx)
            if (snp_idx not in self._snp_lut):
                self._snp_lut[snp_idx] = [] # initialize
            self._snp_lut[snp_idx].append(snp_details) # add the details to the list for snp_idx
        
    def getSNPSetAsASortedList(self):
        return sorted(list(self._snp_set))
    
class PopulationManager:
    _ind_lut = None
    _gen_lut = None
    _p_gen = None
    _cur_gen_no = 0
    _evo_history = None
    _best_ind = None
    _snp_m = None
    _final_pop = None

    def __init__(self, state_list, unit_len, snp_m):
        self._ind_lut = [] # a list of generations of populations
        self._gen_lut = {} # a lut of UUID to a generation to which the individual belongs
        self._evo_history = []
        self._p_gen = PermPatternGenerator(state_list, unit_len)
        self._snp_m = snp_m

    def initializeGeneration(self, generation):
        if (len(self._ind_lut) <= generation): # if the generation is never seen before
            self._ind_lut.append({}) # initialize the dictionary for that generation

    def incrementGeneration(self):
        self._cur_gen_no += 1

    def generateIndividual(self, snp_idx_list): #, generation=None):
#         if generation:
#             gen_no = generation
#         else:
#             gen_no = self._cur_gen_no
        perm_len = len(snp_idx_list)
        ind_pattern = self._p_gen.generatePermutation(perm_len) # get a random individual
        return zip(snp_idx_list, ind_pattern)
        #return self.assignIndividualToGeneration(ind, snp_idx_list, gen_no)

#     def generateOffspring(self, individual): #, generation=None):
#         #import ipdb; ipdb.set_trace()
#         ind_hexkey = individual[-1] # to retrieve the last hex key assigned for the individual 
#         gen_no = self._gen_lut[ind_hexkey] + 1 # next generation w.r.t. the current one
#         ind_pattern = self._ind_lut[gen_no-1][ind_hexkey]['pattern']
#         snp_idx_list = self._ind_lut[gen_no-1][ind_hexkey]['snp_idx_list']
#         
#         offspring_pattern = self.assignIndividualToGeneration(ind_pattern, snp_idx_list, gen_no)
#         individual += offspring_pattern
#         if individual.fitness.valid:
#             self._ind_lut[gen_no][offspring_key[0]]['fitness'] = individual.fitness.values[0] # copy fitness
#         
#         return individual

#     def assignIndividualToGeneration(self, ind_pattern, snp_idx_list, generation=None):
#         if generation:
#             gen_no = generation
#         else:
#             gen_no = self._cur_gen_no
#         self.initializeGeneration(gen_no) # initialize dict if need be
#         ind_key = uuid.uuid4()
#         self._ind_lut[gen_no][ind_key.hex] = {'pattern':None, 'snp_idx_list':None, 'fitness':None}
#         self._ind_lut[gen_no][ind_key.hex]['pattern'] = zip(snp_idx_list, ind_pattern) # assign an UUID look-up to that individual
#         #self._ind_lut[gen_no][ind_key.hex]['snp_idx_list'] = snp_idx_list
#         self._gen_lut[ind_key.hex] = gen_no
#         #self._ind_lut[generation][ind_key.hex]['fitness'] = None
#         return zip(snp_idx_list, ind_pattern)
# 
    def getIndividual(self, ind_hexkey):
        gen_no = self._gen_lut[ind_hexkey]
        return self._ind_lut[gen_no][ind_hexkey]

    def getPopulation(self, generation=None):
        if (generation):
            gen_no = generation
        else:
            gen_no = self._cur_gen_no
        if (len(self._ind_lut) > gen_no):
            return self._ind_lut[gen_no]
        else:
            return None

    def mutateIndividual(self, ind, mut_prob):
        #ind = self.getIndividual(ind_hexkey)
        # do nothing for now.
#         for idx in range(len(ind)):
#             if random.random() > mut_prob:
#                 continue
#             ind[idx] = (ind[idx][0], self._p_gen.generatePermutation(1)[0])
        c_state_idx_list = list(next(zip(*ind)))
        for idx in range(self._snp_m._max_snp_idx+1):
            if (random.random() > mut_prob):
                continue
            if (idx in c_state_idx_list):
                c_idx_idx = bisect(c_state_idx_list, idx) - 1 
                ind[c_idx_idx] = (ind[c_idx_idx][0], self._p_gen.generatePermutation(1)[0])
            else:
                ind.append((idx, self._p_gen.generatePermutation(1)[0]))
        ind.sort(key=itemgetter(0))
        for idx in sorted(range(len(ind)), reverse=True):
            if (ind[idx][1] == 3): # unconstrained
                ind.pop(idx) # remove from the pattern


    def matePopulation(self, generation=None):
        if (generation):
            gen_no = generation
        else:
            gen_no = self._cur_gen_no
        pop = self.getPopulation(gen_no)


    # fitness function
    # this function evaluates individuals against a score board lookup in the simulation study
    # in the Dual-EA trial, a call to FBA will need to be used instead.
    def evaluateFitness(self, individual):
        #import ipdb; ipdb.set_trace()
        #ind_hexkey = individual[-1] # to retrieve the last hex key assigned for the individual
        #ind_data = pm.getIndividual(ind_hexkey)
        ind_pattern = [pattern for idx, pattern in individual[:]]
        #ind_data = pm.getIndividual(individual)
        # now do something with the individual's data to evaluate its fitness
        #fitness = sum(ind_pattern) # interim mock fitness calculation
        totalscore = 0
        ind_set = set(individual)
        for answerkeys in self._snp_m._answersheet:
            for answerkey, score in answerkeys:
                answerset = set(answerkey)
                if (answerset.intersection(ind_set) == answerset):
                    totalscore += score
                    break

        #ind_data['fitness'] = fitness
        individual.fitness.values = (totalscore,)
        return totalscore,

    def printEvolvedAnswers(self, individual):
        #import ipdb; ipdb.set_trace()
        #ind_hexkey = individual[-1] # to retrieve the last hex key assigned for the individual
        #ind_data = pm.getIndividual(ind_hexkey)
        #ind_data = pm.getIndividual(individual)
        # now do something with the individual's data to evaluate its fitness
        #fitness = sum(ind_pattern) # interim mock fitness calculation
        totalscore = 0
        ind_set = set(individual)
        for answerkeys in self._snp_m._answersheet:
            for answerkey, score in answerkeys:
                answerset = set(answerkey)
                if (answerset.intersection(ind_set) == answerset):
                    print(answerset)
                    totalscore += score
                    break

        #ind_data['fitness'] = fitness
        #individual.fitness.values = (totalscore,)

    def getIntersectionOfAnswers(self, population):
        final_set = None
        for ind in population:
            if (not final_set):
                final_set = set(ind)
            else:
                final_set = final_set.intersection(set(ind))

        return final_set

    def getCumulutiveFitnessOfAnswers(self, population):
        cum_fitness = {}
        for ind in population:
            cur_fitness_val = ind.fitness.values[0]
            for idx in range(len(ind)):
                answer_tuple = ind[idx]
                if (cum_fitness.has_key(answer_tuple)):
                    cum_fitness[answer_tuple] += cur_fitness_val
                else:
                    cum_fitness[answer_tuple] = cur_fitness_val

        ordered_answer_list = cum_fitness.items()
        ordered_answer_list.sort(key=itemgetter(1), reverse=True)
        return ordered_answer_list



    def getWeightedCDF(self, population): #gen_no):
        #choices = self.getPopulation(gen_no).items()
        #ind_uuids, ind_details = zip(*choices)
        total = 0
        cum_weights = []
        for ind in population: #ind_details['fitness']:
            total += ind.fitness.values[0]
            cum_weights.append(total)
        return cum_weights
        
    def getWeightedRandomChoice(self, population, cum_weights=None): #choices):
        if (not cum_weights):
            cum_weights = self.getWeightedCDF(population)

        #choices = self.getPopulation(gen_no).items()
        #ind_uuids, ind_details = zip(*choices)
        total = cum_weights[-1]
        num_trials = 0
        i = None
        while (True):
            num_trials += 1
            x = random.random() * total
            i = bisect(cum_weights[:-1], x)
            if (population[i].fitness.valid):
                break
            if (num_trials == 100):
                return None, None
        return i, population[i]

    def getRandomChoice(self, population, except_idx=None):
        #choices = self.getPopulation(gen_no).items()
        #ind_uuids, ind_details = zip(*choices)
        num_trials = 0
        i = None
        while (True):
            num_trials += 1
            i = random.randint(0, len(population)-1)
            if (except_idx):
                while (i == except_idx):
                    i = random.randint(0, len(population)-1)
            if (population[i].fitness.valid):
                break
            if (num_trials == 1000):
                return None, None
        return i, population[i]

    # calling genRanges(100000000,600000000,100)
    # would return a generator of the following tuples
    # (100000000,100000100)
    # (100000100,100000200)
    # (100000200,100000300)
    # ...
    # (599999900,600000000)
    def genRanges(self, start, stop, step):
        current = start
        while (current < stop):
            next_current = current + step
            if (next_current < stop):
                yield (current, next_current)
            else:
                yield (current, stop)
            current = next_current

    # calling genRangesFromList([0,2,5,10,14], 5)
    # would return a generator of the following tuples
    # (0,2)
    # (2,3)
    # (3,4)
    def genRangesFromList(self, idx_list, step):
        cur_mark = (0, idx_list[0])
        for idx_idx, idx in enumerate(idx_list):
            if (idx - cur_mark[1] >= step):
                yield (cur_mark[0], idx_idx)
                cur_mark = (idx_idx, idx)
            elif (idx_idx == len(idx_list) - 1): # last idx
                yield (cur_mark[0], idx_idx)

    def determineCrossoverRange(self, uc_state_idx_list, idx_tuple):
        uc_state_idx_range = uc_state_idx_list[idx_tuple[0]:idx_tuple[1]]
        range_len = len(uc_state_idx_range)
        range_minmax_diff = uc_state_idx_list[idx_tuple[1]] - uc_state_idx_list[idx_tuple[0]]
        if (not (range_minmax_diff == range_len)): # if constraints exist in this range
            # perform crossover
            cxpoint = random.randint(0, range_len - 1)
#             cxpoint2 = random.randint(1, range_len - 2)
#             if cxpoint2 >= cxpoint1: 
#                     cxpoint2 += 1
#             else: # Swap the two cx points if order is not right
#                 cxpoint1, cxpoint2 = cxpoint2, cxpoint1
            # now find the corresponding crossover idxs for constrained sets
            if (range_len - cxpoint > cxpoint - 0):
                cx_idx_min = uc_state_idx_range[cxpoint]
                cx_idx_max = uc_state_idx_range[range_len-1]
            else:
                cx_idx_min = uc_state_idx_range[0]
                cx_idx_max = uc_state_idx_range[cxpoint]

            return (cx_idx_min, cx_idx_max)

        return None

#             ind1_cx_idxs = [idx for idx in range(cx_idx_min, cx_idx_max+1) if idx in ind1_state_idx_set]
#             ind2_cx_idxs = [idx for idx in range(cx_idx_min, cx_idx_max+1) if idx in ind2_state_idx_set]
            
    # crossover function
    # ind1 and ind2 are of type numpy.ndarray
    # carrying an array of N genes as part of the metabolic pathway model.
    def crossoverSNPAlleles(self, ind1, ind2, cx_prob):
        """ 
        ind1 and ind2 are of type IDEA_MPModel
        IDEA_MPModel consists of a numpy.array of the following dictionary items:
        is_constrained is a flux constraint condition of a boolean type;
        if a gene and its associated metabolic enzyme's flux is_constrained,  
        the lowerbound and upperbound of that flux is adjusted to account for the constraint. 
        If the constraint is to silence the flux, lowerbound = upperbound = 0
        If the constraint is to have a positive flux, lowerbound = 0, and upperbound > 0
        If the constraint is to have a negative flux, lowerbound < 0, and upperbound = 0
        
        Execute a two points crossover with copy on the input individuals. The
        copy is required because the slicing in numpy returns a view of the data,
        which leads to a self overwritting in the swap operation. It prevents
        ::
        
            >>> import numpy
            >>> a = numpy.array((1,2,3,4))
            >>> b = numpy.array((5.6.7.8))
            >>> a[1:3], b[1:3] = b[1:3], a[1:3]
            >>> print(a)
            [1 6 7 4]
            >>> print(b)
            [5 6 7 8]
        """

        # sort individual state patterns by the state_idx
        ind1.sort(key=itemgetter(0))
        ind2.sort(key=itemgetter(0))

        # get the two state_idx sets
        ind1_state_idx_set = set(next(zip(*ind1)))
        ind2_state_idx_set = set(next(zip(*ind2)))
#         common_state_idx = ind1_state_idx_set.intersection(ind2_state_idx_set)

        # using the homology of smallest unconstrained regions enclosing the constrained sets
        # this will solve the potential problem of the common constrained set 
        # being too small (size<2) for crossover homology
        ind1_unconstrained_idx_set = set(range(min(ind1_state_idx_set)-1, max(ind1_state_idx_set)+2)).difference(ind1_state_idx_set)
        ind2_unconstrained_idx_set = set(range(min(ind2_state_idx_set)-1, max(ind2_state_idx_set)+2)).difference(ind2_state_idx_set)
        uc_state_idx_union = ind1_unconstrained_idx_set.union(ind2_unconstrained_idx_set)
        size = len(uc_state_idx_union)

#         cxpoint1 = random.randint(1, size - 1)
#         cxpoint2 = random.randint(1, size - 2)
#         if cxpoint2 >= cxpoint1: 
#                 cxpoint2 += 1
#         else: # Swap the two cx points if order is not right
#             cxpoint1, cxpoint2 = cxpoint2, cxpoint1
        
        # now find the corresponding crossover idxs for constrained sets
        uc_state_idx_list = list(uc_state_idx_union)

        idx_tuple_list = list(self.genRangesFromList(uc_state_idx_list, step=20))

        for idx_tuple in idx_tuple_list:
            #print 'idx_tuple: ' + str(idx_tuple)
            #import ipdb; ipdb.set_trace()
            if (random.random() > cx_prob):
                continue

            cx_range = self.determineCrossoverRange(uc_state_idx_list, idx_tuple)
            if (cx_range):
                cx_idx_min, cx_idx_max = cx_range
                ind1_c_state_idxs = [idx for idx in range(cx_idx_min, cx_idx_max+1) if idx in ind1_state_idx_set]
                ind2_c_state_idxs = [idx for idx in range(cx_idx_min, cx_idx_max+1) if idx in ind2_state_idx_set]
                # example c_state_idxs given ind1 of [(3, 2.0), (5, 1.0), (1000, 3.0)]
                # a c_state_idxs would be [3, 5]
                # do the crossover
                ind1_cx_idxs = []
                ind2_cx_idxs = []
                if (ind1_c_state_idxs):
#                     crossover_range1 = ind1[min(ind1_c_state_idxs):max(ind1_c_state_idxs)+1] 
                    ind1_cx_idxs = [idx_idx for idx_idx, entry in enumerate(ind1) if entry[0] in ind1_c_state_idxs]
                if (ind2_c_state_idxs):
#                     crossover_range2 = ind2[min(ind2_c_state_idxs):max(ind2_c_state_idxs)+1]
                    ind2_cx_idxs = [idx_idx for idx_idx, entry in enumerate(ind2) if entry[0] in ind2_c_state_idxs]
                ind1_cx_idxs.sort(reverse=True)
                for cx_candidate_idx in ind1_cx_idxs:
                    cx_candidate = ind1.pop(cx_candidate_idx)
                    ind2.append(cx_candidate)

                ind2_cx_idxs.sort(reverse=True)
                for cx_candidate_idx in ind2_cx_idxs:
                    cx_candidate = ind2.pop(cx_candidate_idx)
                    ind1.append(cx_candidate)
                
                ind1.sort(key=itemgetter(0))
                ind2.sort(key=itemgetter(0))
#         cx_idx_min = uc_state_idx_list[cxpoint1]
#         cx_idx_max = uc_state_idx_list[cxpoint2]
#         ind1_cx_idxs = [idx for idx in range(cx_idx_min, cx_idx_max+1) if idx in ind1_state_idx_set]
#         ind2_cx_idxs = [idx for idx in range(cx_idx_min, cx_idx_max+1) if idx in ind2_state_idx_set]
        # do the crossover
#         crossover_range1 = []
#         crossover_range2 = []
#         if ind1_cx_idxs:
#             crossover_range1 = ind1[min(ind1_cx_idxs):max(ind1_cx_idxs)+1] 
#         if ind2_cx_idxs:
#             crossover_range2 = ind2[min(ind2_cx_idxs):max(ind2_cx_idxs)+1]
#         for cx_candidate_idx in crossover_range1:
#             cx_candidate = ind1.pop(cx_candidate_idx)
#             ind2.append(cx_candidate)
# 
#         for cx_candidate_idx in crossover_range2:
#             cx_candidate = ind2.pop(cx_candidate_idx)
#             ind1.append(cx_candidate)
#         
#         ind1.sort(key=itemgetter(0))
#         ind2.sort(key=itemgetter(0))

#         size = len(ind1)
#         cxpoint1 = random.randint(1, size)
#         cxpoint2 = random.randint(1, size - 1)
#         if cxpoint2 >= cxpoint1:
#             cxpoint2 += 1
#         else: # Swap the two cx points
#             cxpoint1, cxpoint2 = cxpoint2, cxpoint1
# 
#         ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
#             = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()
        del ind1.fitness.values
        del ind2.fitness.values
            
        return ind1, ind2


class PermPatternGenerator:

    _state_list = None
    _unit_len = None
    _lut = {}

    def __init__(self, state_list, unit_len):
        self._state_list = state_list
        self._unit_len = unit_len

        # initialize permutation luts
        self._lut[unit_len] = np.array(list(itertools.product(self._state_list, repeat=unit_len)))
        for short_split in range(1,unit_len):
            self._lut[short_split] = np.array(list(itertools.product(self._state_list, repeat=short_split)))
        

    # a function to randomly generate an individual pattern to form a population
    def generatePermutation(self, perm_len):
        perm_pattern = np.array([]) # initialize    
        split_idx_list = range(perm_len)[::self._unit_len] # splitting a perm_len long array into _unit_len long chunks
        last_idx = split_idx_list[-1] # get the last idx in the list
        unit_lut = self._lut[self._unit_len]
        unit_lut_len = len(unit_lut)
        for start_idx in split_idx_list:
            if (start_idx != last_idx):
                rand_idx = random.randint(0, unit_lut_len-1)
                unit_pattern = unit_lut[rand_idx].copy()
                perm_pattern = np.append(perm_pattern, unit_pattern)
            else:
                last_lut = self._lut[perm_len - last_idx]
                last_lut_len = len(last_lut)
                rand_idx = random.randint(0, last_lut_len-1)
                last_pattern = last_lut[rand_idx].copy()
                perm_pattern = np.append(perm_pattern, last_pattern)
            
        return perm_pattern

    
