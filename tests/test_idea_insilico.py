import pytest
from copy import deepcopy

from deap import base, creator, tools

from gsap.sim.idea_insilico import PopulationManager, SNPSetManager

# arrange fixtures
@pytest.fixture
def gsap_sim_common(tmp_path):
    snp_m = SNPSetManager(2000)
    # snp_m.addSNPs([(1940, 'mutS'), (1954, 'ribD'), (1000, 'trpC'), (1999, 'trpA'), (1995, 'yheD')])
    snp_m.addSNPs(
        [(1940, "mutS"), (1954, "ribD"), (1000, "trpC"), (1999, "trpA"), (1995, "yheD")]
    )
    snp_m._answersheet = [
        [
            ([(18, 0), (200, 1), (360, 1), (942, 2), (1000, 0)], 10),
            ([(18, 0), (294, 1)], 5),
            ([(294, 1)], 3),
            ([(1000, 0)], 3),
        ],
        [
            ([(1900, 2), (1940, 0), (1954, 0), (1994, 0), (1995, 1), (1998, 1)], 20),
            ([(1900, 2), (1954, 0), (1995, 1)], 15),
            ([(1900, 2), (1954, 0)], 10),
            ([(1940, 0), (1954, 0)], 9),
            ([(1900, 2)], 4),
            ([(1940, 0)], 4),
            ([(1954, 0)], 3),
        ],
    ]
    # state space: 0 for silenced, 1 for pos flux, 2 for neg flux, 3 for unconstrained.
    pm = PopulationManager([0, 1, 2, 3], 5, snp_m)

    creator.create("FitnessMax", base.Fitness, weights=(1000.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    # creator.create("Individual", str, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()

    # toolbox.register("attr_state", random.randint, 0, 3) # is an element out of a set {0, 1, 2, 3}
    toolbox.register(
        "attr_state", pm.generateIndividual, snp_m.getSNPSetAsASortedList()
    )  # a 3 var long pattern at 0th gen
    # 0 denotes the unconstrained, 1 the muted, 2 the positive flux, and 3 the negative flux states
    toolbox.register(
        "individual", tools.initIterate, creator.Individual, toolbox.attr_state
    )
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)


    toolbox.register("evaluate", pm.evaluateFitness)
    toolbox.register("mate", pm.crossoverSNPAlleles)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)
    toolbox.register("clone", deepcopy)

    return (toolbox, snp_m, pm)


@pytest.mark.skip(reason="test performed manually, as it takes too long.")
def test_idea_insilico(gsap_sim_common):  # {
    toolbox, snp_m, pm = gsap_sim_common
    # random.seed(64)

    pop = toolbox.population(n=3000)

    # CXPB: prob of crossover
    # MUTPB: prob of mutation
    # NGEN: number of generations evo runs
    CXPB, MUTPB, NGEN = 0.5, 0.001, 40

    # evaluate the population
    fitnesses = [toolbox.evaluate(ind) for ind in pop]
    print("  Evaluated %i individuals" % len(pop))

    # evolution begins
    for g in range(NGEN):  # {
        print("#### Generation %i ####" % g)

        # selection
        offspring_candidates = toolbox.select(pop, len(pop))
        # clone the candidates and generate offsprings out of the clones
        # clones = [toolbox.clone(candidate) for candidate in offspring_candidates]
        offsprings = list(map(toolbox.clone, offspring_candidates))
        # import ipdb; ipdb.set_trace()
        # offsprings = [pm.generateOffspring(clone) for clone in clones]

        # crossover in the current generation
        cum_weights = pm.getWeightedCDF(offsprings)
        for i in range(int(len(offsprings) / 2)):
            # pick a random offspring weighted by fitness
            # pick a second random offspring unweighted
            # crossover
            ind1_idx, ind1 = pm.getWeightedRandomChoice(offsprings, cum_weights)
            ind2_idx, ind2 = pm.getRandomChoice(offsprings, except_idx=ind1_idx)
            if ind1 and ind2:
                pm.crossoverSNPAlleles(ind1, ind2, CXPB)
        # mutagenesis
        for offspring in offsprings:
            pm.mutateIndividual(offspring, MUTPB)

        # evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offsprings if not ind.fitness.valid]
        fitnesses = [toolbox.evaluate(ind) for ind in invalid_ind]

        print("evaluated %i offsprings" % len(invalid_ind))

        # generation hand over
        pop[:] = offsprings

        # # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x * x for x in fits)
        std = abs(sum2 / length - mean**2) ** 0.5
        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)
    # }

    print("-- End of (successful) evolution --")
    best_ind = tools.selBest(pop, 1)[0]
    for history in pop:
        pm._evo_history.append(history)
    pm._best_ind = best_ind
    pm._final_pop = pop

    return True


#     print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
#     # Numpy equality function (operators.eq) between two arrays returns the
#     # equality element wise, which raises an exception in the if similar()
#     # check of the hall of fame. Using a different equality function like
#     # numpy.array_equal or numpy.allclose solve this issue.
#     hof = tools.HallOfFame(1, similar=np.array_equal)
#
#     stats = tools.Statistics(lambda ind: ind.fitness.values)
#     stats.register("avg", np.mean)
#     stats.register("std", np.std)
#     stats.register("min", np.min)
#     stats.register("max", np.max)
#
#     algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=40, stats=stats,
#                         halloffame=hof)
#
#     return pop, stats, hof
# }

