import numpy as np

def total_spins_old(variable_count, variable_bit_count, triplet_count):
    return (variable_bit_count*variable_count) + triplet_count

def total_spins(variable_count, variable_bit_count, triplet_count):
    return (variable_bit_count*variable_count) + (2*triplet_count)

def spin_index_old(variable_idx, bit_idx, variable_count, variable_bit_count, \
    triplet_count):
    if (type(bit_idx) is str) and (bit_idx.lower() == 'triplet'):
        num_spins = total_spins(variable_count, variable_bit_count, \
            triplet_count)
        spin_idx = num_spins - (variable_idx+1)
    elif (type(bit_idx) is int) and (bit_idx >= 0) and \
        (bit_idx < variable_bit_count):
        spin_idx = (variable_bit_count*variable_idx) + bit_idx
    else:
        raise ValueError
    return spin_idx

def spin_index(variable_idx, bit_idx, variable_count, variable_bit_count, \
    triplet_count):
    if (type(bit_idx) is tuple) and (bit_idx[0].lower() == 'triplet'):
        num_spins = total_spins(variable_count, variable_bit_count, \
            triplet_count)
        if bit_idx[1] == 0:
            spin_idx = num_spins - ((2*variable_idx)+2)
        elif bit_idx[1] == 1:
            spin_idx = num_spins - ((2*variable_idx)+1)
        else:
            raise ValueError
    elif (type(bit_idx) is int) and (bit_idx >= 0) and \
        (bit_idx < variable_bit_count):
        spin_idx = (variable_bit_count*variable_idx) + bit_idx
    else:
        raise ValueError
    return spin_idx

def variable_times_variable(variable1_idx, variable2_idx, variable1_sign, \
    variable2_sign, variable_bit_count, variable_count, triplet_count):

    v1_spins = tuple(
        spin_index(variable1_idx, bit_idx, variable_count, variable_bit_count, \
            triplet_count) for bit_idx in range(variable_bit_count)
    )

    v2_spins = tuple(
        spin_index(variable2_idx, bit_idx, variable_count, variable_bit_count, \
            triplet_count) for bit_idx in range(variable_bit_count)
    )

    J, h, c = {}, {}, 0.0

    c += 1.0

    for j in range(variable_bit_count):
        h[v1_spins[j]] = h.get(v1_spins[j], 0.0) - variable1_sign*(2**j)/(2**variable_bit_count)
        h[v2_spins[j]] = h.get(v2_spins[j], 0.0) - variable2_sign*(2**j)/(2**variable_bit_count)

    for j in range(variable_bit_count):
        for k in range(variable_bit_count):
            coeff = variable1_sign*variable2_sign*\
                (2**(j+k))/(4**variable_bit_count)
            tmp1, tmp2 = v1_spins[j], v2_spins[k]
            if tmp1 == tmp2:
                c += coeff
            else:
                J[(tmp1, tmp2)] = J.get((tmp1, tmp2), 0.0) + 0.5*coeff
                J[(tmp2, tmp1)] = J.get((tmp2, tmp1), 0.0) + 0.5*coeff

    J = {key: (np.pi**2)*val for key, val in J.items()}
    h = {key: (np.pi**2)*val for key, val in h.items()}
    c = (np.pi**2)*c
    return J, h, c

def variable_times_sum(variable_idx, variable_sign, variable_bit_count, \
    variable_count, triplet_idx, triplet_count):

    t_spin0 = spin_index(triplet_idx, ('triplet', 0), variable_count, \
        variable_bit_count, triplet_count)
    t_spin1 = spin_index(triplet_idx, ('triplet', 1), variable_count, \
        variable_bit_count, triplet_count)


    v_spins = tuple(
        spin_index(variable_idx, bit_idx, variable_count, variable_bit_count, \
            triplet_count) for bit_idx in range(variable_bit_count)
    )

    J, h, c = {}, {}, 0.0

    c += 3.0

    h[t_spin0] = h.get(t_spin0, 0.0) - 1
    h[t_spin1] = h.get(t_spin1, 0.0) - 2
    for j in range(variable_bit_count):
        h[v_spins[j]] = h.get(v_spins[j], 0.0) - 3*variable_sign*(2**j)/(2**variable_bit_count)
        coeff = variable_sign*(2**j)/(2**variable_bit_count)
        J[(v_spins[j], t_spin0)] = J.get((v_spins[j], t_spin0), 0.0) + 0.5*coeff
        J[(t_spin0, v_spins[j])] = J.get((t_spin0, v_spins[j]), 0.0) + 0.5*coeff
        coeff = variable_sign*(2**j)/(2**(variable_bit_count-1))
        J[(v_spins[j], t_spin1)] = J.get((v_spins[j], t_spin1), 0.0) + 0.5*coeff
        J[(t_spin1, v_spins[j])] = J.get((t_spin1, v_spins[j]), 0.0) + 0.5*coeff

    J = {key: (np.pi**2)*val for key, val in J.items()}
    h = {key: (np.pi**2)*val for key, val in h.items()}
    c = (np.pi**2)*c
    return J, h, c

def sum_times_sum(variable_bit_count, variable_count, triplet_idx, \
    triplet_count):

    t_spin0 = spin_index(triplet_idx, ('triplet', 0), variable_count, \
        variable_bit_count, triplet_count)
    t_spin1 = spin_index(triplet_idx, ('triplet', 1), variable_count, \
        variable_bit_count, triplet_count)

    J, h, c = {}, {}, 0.0

    c += 14.0
    h[t_spin0] = h.get(t_spin0, 0.0) - 6
    h[t_spin1] = h.get(t_spin1, 0.0) - 12

    J[(t_spin0, t_spin1)] = J.get((t_spin0, t_spin1), 0.0) + 0.5*4
    J[(t_spin1, t_spin0)] = J.get((t_spin1, t_spin0), 0.0) + 0.5*4

    J = {key: (np.pi**2)*val for key, val in J.items()}
    h = {key: (np.pi**2)*val for key, val in h.items()}
    c = (np.pi**2)*c
    return J, h, c

def merge_dict_by_sum(d1, d2):
    d = dict(d1)
    for key, val in d2.items():
        try:
            d[key] += val
        except KeyError:
            d[key] = val
    return d

def merge_dict_by_diff(d1, d2):
    d = dict(d1)
    for key, val in d2.items():
        try:
            d[key] -= val
        except KeyError:
            d[key] = -val
    return d

def triplet_hamiltonian(variable1_idx, variable2_idx, variable3_idx, \
    variable1_sign, variable2_sign, variable3_sign, variable_count, \
    variable_bit_count, triplet_idx, triplet_count):
    
    J, h, c = {}, {}, 0.0

    # t1t1
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable1_idx, variable1_sign, variable1_sign, variable_bit_count, variable_count, triplet_count)
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += cterm

    # t2t2
    Jterm, hterm, cterm = variable_times_variable(variable2_idx, \
        variable2_idx, variable2_sign, variable2_sign, variable_bit_count, variable_count, triplet_count)
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += cterm

    # t3t3
    Jterm, hterm, cterm = variable_times_variable(variable3_idx, \
        variable3_idx, variable3_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += cterm

    # t1t2
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable2_idx, variable1_sign, variable2_sign, variable_bit_count, variable_count, triplet_count)
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += 2*cterm

    # t1t3
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable3_idx, variable1_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += 2*cterm

    # t2t3
    Jterm, hterm, cterm = variable_times_variable(variable2_idx, \
        variable3_idx, variable2_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += 2*cterm

    #t1t
    Jterm, hterm, cterm = variable_times_sum(variable1_idx, variable1_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_diff(J, Jterm)
    h = merge_dict_by_diff(h, hterm)
    c -= 2*cterm


    #t2t
    Jterm, hterm, cterm = variable_times_sum(variable2_idx,  variable2_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_diff(J, Jterm)
    h = merge_dict_by_diff(h, hterm)
    c -= 2*cterm

    #t3t
    Jterm, hterm, cterm = variable_times_sum(variable3_idx,  variable3_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_diff(J, Jterm)
    h = merge_dict_by_diff(h, hterm)
    c -= 2*cterm

    #tt
    Jterm, hterm, cterm = sum_times_sum(variable_bit_count, variable_count, \
        triplet_idx, triplet_count)
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += cterm

    return J, h, c

def ham_dicts_to_ham(nqubits, Jdict, hdict):
    J = np.zeros((nqubits, nqubits), dtype=float)
    h = np.zeros(nqubits, dtype=float)
    for key, val in Jdict.items():
        J[key[0], key[1]] = val
    for key, val in hdict.items():
        h[key[0]] = val
    return J, h

def hamiltonian(variable_count, variable_bit_count, triplet_refls, \
    triplet_signs, triplet_weights=None, verbose=False):
    if triplet_weights is None:
        triplet_weights = [1.0 for t in triplet_refls]

    triplet_count = len(triplet_refls)

    J, h, c = {}, {}, 0.0

    for j, triplet in enumerate(triplet_refls):
        if verbose:
            print(f"Processing triplet {j+1} of {triplet_count}...   ", end="\r")
        triplet_idx = j
        variable1_idx, variable2_idx, variable3_idx = tuple(triplet)
        variable1_sign, variable2_sign, variable3_sign = tuple(triplet_signs[j])

        Jt, ht, ct = triplet_hamiltonian(variable1_idx, variable2_idx, \
            variable3_idx, variable1_sign, variable2_sign, variable3_sign, \
            variable_count, variable_bit_count, triplet_idx, triplet_count)

        Jt = {key:triplet_weights[j]*val for key, val in Jt.items()}
        ht = {key:triplet_weights[j]*val for key, val in ht.items()}
        J = merge_dict_by_sum(J, Jt)
        h = merge_dict_by_sum(h, ht)
        c += triplet_weights[j]*ct

    J = {key:val for key, val in J.items() if not (val == 0.0) }
    h = {key:val for key, val in h.items() if not (val == 0.0) }

    nqubits = total_spins(variable_count, variable_bit_count, triplet_count)

    J, h = ham_dicts_to_ham(nqubits, J, h)

    return J, h, c

