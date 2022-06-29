import numpy as np
from copy import deepcopy as dc

def initial_state(annealing_inputs):
    n_spins = annealing_inputs[0].shape[0]
    return np.random.default_rng().choice([-1, 1], size=n_spins)

def ising_cost(annealing_inputs, state):
    J, h, c = tuple(annealing_inputs)
    return np.dot(state, np.dot(J, state)) + np.dot(h, state) + c

def ising_neighbour(annealing_inputs, current_state):
    new_state = dc(current_state)
    n_spins = len(new_state)
    spin = np.random.default_rng().choice(n_spins)
    new_state[spin] *= -1
    return new_state

def boltzmann_acceptance_rule(current_cost, candidate_cost, temperature):
    if candidate_cost <= current_cost:
        accept = True
    else:

        acceptance_probability = \
            np.exp(-(candidate_cost-current_cost)/temperature)

        test_probability = np.random.default_rng().uniform()
        accept = test_probability <= acceptance_probability
    return accept

def temperature_schedule(annealing_inputs, iterations, iteration):
    T_max = 10.0
    scale = np.log((iterations+1)/(iteration+1))/np.log(iterations+1)
    T = T_max*scale
    return T


def simulated_annealing_run(extra_inputs, iterations,
                           initial_state_generator,
                           cost_function,
                           candidate_generator,
                           acceptance_rule,
                           acceptance_parameter_generator,
                           verbose=False):

    best_cost = float('inf')
    state = initial_state_generator(extra_inputs)
    cost = cost_function(extra_inputs, state)
    costs = np.zeros(iterations+1, dtype=float)
    costs[0] = cost
    for iteration in range(iterations):
        if verbose:
            pc = 100*(iteration)/iterations
            print(f"{pc:.2f}% complete.", end="\r")

        acceptance_parameter = acceptance_parameter_generator(extra_inputs, \
            iterations, iteration)
            
        candidate_state = candidate_generator(extra_inputs, state)
        candidate_cost = cost_function(extra_inputs, candidate_state)
        accept = acceptance_rule(cost, candidate_cost, acceptance_parameter)
        if accept:
            state, cost = candidate_state, candidate_cost
            if cost < best_cost:
                best_state, best_cost = state, cost
        costs[iteration+1] = cost
    return best_state, best_cost, costs

def simulated_annealing(extra_inputs, iterations, runs,
                        initial_state_generator,
                        cost_function,
                        candidate_generator,
                        acceptance_rule,
                        acceptance_parameter_generator,
                        end_cost=None,
                        store_cost=None,
                        store=None,
                        verbose=False):
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
                             verbose=inner_verbose)
        if (not (store_cost is None)) and cost <= store_cost:
            store.append([cost, state])
        if cost < best_cost:
            best_state, best_cost, best_costs = state, cost, costs
        if (not (end_cost is None)) and (best_cost <= end_cost):
            if verbose: print(f"breaking early with cost {best_cost}...")
            break
    return best_state, best_cost, best_costs

def anneal(J, h, c, iterations, runs, end_cost=None, store_cost=None, \
    verbose=False):

    if (not (store_cost is None)):
        store = []
    else:
        store = None

    state, cost, costs = simulated_annealing((J, h, c), iterations, runs, \
        initial_state, ising_cost, ising_neighbour, boltzmann_acceptance_rule, \
        temperature_schedule, end_cost=end_cost, store_cost=store_cost, \
        store=store, verbose=verbose)

    if (not (store_cost is None)):
        return state, cost, costs, store
    else:
        return state, cost, costs
    