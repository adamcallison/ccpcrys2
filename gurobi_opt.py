import numpy as np
import gurobipy as gp

def ising_to_qubo(J, h, c):
    M = 4*J
    v = np.zeros_like(h)
    nspins = len(h)
    for j in range(nspins):
        val = 2*(h[j] - ( np.sum(J[j]) + np.sum(J[:,j]) ))
        v[j] = val
    k = c + np.sum(J) + np.sum(h)
    return M, v, k

def optimize(J, h, c, time_limit=None, heuristic_frac=None, verbose=False):
    if (not (heuristic_frac is None)) and (time_limit is None):
        raise ValueError
    m=gp.Model("ccpcrys")
    if not verbose:
        m.Params.LogToConsole = 0
    nspins = len(h)
    b = m.addMVar(nspins, vtype=gp.GRB.BINARY, name='b')
    M, v, k = ising_to_qubo(J, h, c)
    obj = (b@M@b) + (v@b) + k
    m.setObjective(obj, gp.GRB.MINIMIZE)
    if not (time_limit is None):
        m.params.TimeLimit = time_limit
    if not (heuristic_frac is None):
        m.params.Heuristics = heuristic_frac


    m.params.PoolSolutions = 1000
    m.params.PoolSearchMode = 2

    m.optimize()

    sols = []
    costs = []
    solcount = m.SolCount
    for j in range(solcount):
        m.params.SolutionNumber = j
        state = 1 - (2*np.array(m.Xn))
        cost = np.dot(state, np.dot(J, state)) + np.dot(h, state) + c
        sols.append(state)
        costs.append(cost)
    costs = np.array(costs)

    return sols, costs