import numpy as np
from copy import deepcopy as dc

n_spins_curr = None
random_ints = None
random_ints_used = 0
random_floats = None
random_floats_used = 0

def _ising_neighbour_helper(n_spins):
    global n_spins_curr
    global random_ints
    global random_ints_used
    if n_spins_curr == n_spins:
        if random_ints_used < len(random_ints):
            res = random_ints[random_ints_used]
            random_ints_used += 1
        else:
            random_ints = np.random.default_rng().choice(n_spins, size=16384)
            res = random_ints[0]
            random_ints_used = 1
    else:
        n_spins_curr = n_spins
        random_ints = np.random.default_rng().choice(n_spins, size=16384)
        res = random_ints[0]
        random_ints_used = 1
    return res

def _acceptance_helper():
    global random_floats
    global random_floats_used
    if (not (random_floats is None)) and (random_floats_used < len(random_floats)):
        res = random_floats[random_floats_used]
        random_floats_used += 1
    else:
        random_floats = np.random.default_rng().uniform(size=16384)
        res = random_floats[0]
        random_floats_used = 1
    return res

def _reset_helpers():
    global n_spins_curr
    global random_ints
    global random_ints_used
    global random_floats
    global random_floats_used

    n_spins_curr = None
    random_ints = None
    random_ints_used = 0
    random_floats
    random_floats_used = 0

def initial_state(annealing_inputs):
    n_spins = annealing_inputs[0].shape[0]
    return np.random.default_rng().choice([-1, 1], size=n_spins)

def ising_cost(annealing_inputs, state):
    J, h, c = tuple(annealing_inputs)[:3]
    return np.dot(state, np.dot(J, state)) + np.dot(h, state) + c

def ising_neighbour(annealing_inputs, current_state):
    new_state = dc(current_state)
    n_spins = len(new_state)
    #spin = np.random.default_rng().choice(n_spins)
    spin = _ising_neighbour_helper(n_spins)
    new_state[spin] *= -1
    return new_state

def boltzmann_acceptance_rule(current_cost, candidate_cost, temperature):
    if candidate_cost <= current_cost:
        accept = True
    else:

        acceptance_probability = \
            np.exp(-(candidate_cost-current_cost)/temperature)

        #test_probability = np.random.default_rng().uniform()
        test_probability = _acceptance_helper()
        accept = test_probability <= acceptance_probability
    return accept


def temperature_schedule(annealing_inputs, acceptance_parameter_generator_inputs, iterations, iteration):
    #T_max = 100.0

    T_max = acceptance_parameter_generator_inputs[0]
    scale = np.log((iterations+1)/(iteration+1))/np.log(iterations+1)

    T = T_max*scale
    return T


def simulated_annealing_run(extra_inputs, iterations,
                           initial_state_generator,
                           cost_function,
                           candidate_generator,
                           acceptance_rule,
                           acceptance_parameter_generator,
                           acceptance_parameter_generator_inputs,
                           verbose=False):

    best_cost = float('inf')
    state = initial_state_generator(extra_inputs)
    cost = cost_function(extra_inputs, state)
    #costs = np.zeros(iterations+1, dtype=float)
    costs = []
    costs.append(cost)
    if verbose: last_print = -np.float('inf')
    last_costsave = 0.0
    for iteration in range(iterations):
        if verbose:
            pc = 100*(iteration)/iterations
            if pc - last_print >= 0.999:
                last_print = pc
                print(f"{pc:.2f}% complete. Starting_cost={costs[0]}. Current cost={cost}.", end="\r")

        acceptance_parameter = acceptance_parameter_generator(extra_inputs, acceptance_parameter_generator_inputs, \
            iterations, iteration)
            
        candidate_state = candidate_generator(extra_inputs, state)
        candidate_cost = cost_function(extra_inputs, candidate_state)
        accept = acceptance_rule(cost, candidate_cost, acceptance_parameter)
        if accept:
            state, cost = candidate_state, candidate_cost
            if cost < best_cost:
                best_state, best_cost = state, cost
        pm = (iteration/iterations)*1000
        if pm - last_costsave >= 1.0:
            last_costsave = pm
            costs.append(cost)
    costs = np.array(costs)
    return best_state, best_cost, costs

def simulated_annealing(extra_inputs, iterations, runs,
                        initial_state_generator,
                        cost_function,
                        candidate_generator,
                        acceptance_rule,
                        acceptance_parameter_generator,
                        acceptance_parameter_generator_inputs,
                        end_cost=None,
                        store_cost=None,
                        store=None,
                        verbose=False):

    _reset_helpers()

    if verbose:
        if runs == 1:
            outer_verbose, inner_verbose = False, True
        else:
            outer_verbose, inner_verbose = True, False
    else:
        outer_verbose, inner_verbose = False, False

    if (store_cost is None) and (not (store is None)):
        raise ValueError

    if (not (store_cost is None)) and (store is None):
        raise ValueError

    best_cost = float('inf')
    for run in range(runs):
        if outer_verbose:
            pc = 100*(run)/runs
            print(f"{pc:.2f}% complete. (best cost: {best_cost})   ", end="\r")

        state, cost, costs = simulated_annealing_run(extra_inputs, iterations,
                             initial_state_generator,
                             cost_function,
                             candidate_generator,
                             acceptance_rule,
                             acceptance_parameter_generator,
                             acceptance_parameter_generator_inputs,
                             verbose=inner_verbose)
        if (not (store_cost is None)) and cost <= store_cost:
            store.append([cost, state])
        if cost < best_cost:
            best_state, best_cost, best_costs = state, cost, costs
        if (not (end_cost is None)) and (best_cost <= end_cost):
            if verbose: print(f"breaking early with cost {best_cost}...")
            break
    return best_state, best_cost, best_costs

#def anneal(J, h, c, n_refls, refl_size, n_triplets, iterations, runs, max_temperature=10.0, end_cost=None, store_cost=None, \
def anneal(J, h, c, iterations, runs, max_temperature=10.0, end_cost=None, store_cost=None, \
    input_state=None, verbose=False):

    if (not (store_cost is None)):
        store = []
    else:
        store = None

    if input_state is None:
        isfunc = initial_state
    else:
        isfunc = lambda annealing_inputs: input_state

    #state, cost, costs = simulated_annealing((J, h, c, n_refls, refl_size, n_triplets), iterations, runs, \
    state, cost, costs = simulated_annealing((J, h, c), iterations, runs, \
        isfunc, ising_cost, ising_neighbour, boltzmann_acceptance_rule, \
        temperature_schedule, (max_temperature,), end_cost=end_cost, store_cost=store_cost, \
        store=store, verbose=verbose)

    if (not (store_cost is None)):
        return state, cost, costs, store
    else:
        return state, cost, costs
    