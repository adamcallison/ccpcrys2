import numpy as np

def total_spins(variable_count, variable_bit_count, triplet_count):
    return (variable_bit_count*variable_count) + triplet_count

def spin_index(variable_idx, bit_idx, variable_count, variable_bit_count, \
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

def zero_hamiltonian_old(variable_count, variable_bit_count, triplet_count):
    num_spins = total_spins(variable_count, variable_bit_count, triplet_count)
    J = np.zeros((num_spins, num_spins), dtype=float)
    h = np.zeros(num_spins, dtype=float)
    c = 0.0
    return J, h, c

def zero_hamiltonian(variable_count, variable_bit_count, triplet_count):
    # has unnecessary inputs
    J, h, c = {}, {}, 0.0
    return J, h, c

def variable_times_variable_old(variable1_idx, variable2_idx, variable1_sign, \
    variable2_sign, variable_bit_count, variable_count, triplet_count):

    v1_spins = tuple(
        spin_index(variable1_idx, bit_idx, variable_count, variable_bit_count, \
            triplet_count) for bit_idx in range(variable_bit_count)
    )

    v2_spins = tuple(
        spin_index(variable2_idx, bit_idx, variable_count, variable_bit_count, \
            triplet_count) for bit_idx in range(variable_bit_count)
    )

    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    c += 1.0

    for j in range(variable_bit_count):
        h[v1_spins[j]] -= variable1_sign*(2**j)/(2**variable_bit_count)
        h[v2_spins[j]] -= variable2_sign*(2**j)/(2**variable_bit_count)

    for j in range(variable_bit_count):
        for k in range(variable_bit_count):
            coeff = variable1_sign*variable2_sign*\
                (2**(j+k))/(4**variable_bit_count)
            tmp1, tmp2 = v1_spins[j], v2_spins[k]
            if tmp1 == tmp2:
                c += coeff
            else:
                J[tmp1, tmp2] += 0.5*coeff
                J[tmp2, tmp1] += 0.5*coeff

    return (np.pi**2)*J, (np.pi**2)*h, (np.pi**2)*c

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

    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    c += 1.0

    for j in range(variable_bit_count):
        #h[v1_spins[j]] -= variable1_sign*(2**j)/(2**variable_bit_count)
        #h[v2_spins[j]] -= variable2_sign*(2**j)/(2**variable_bit_count)
        h[v1_spins[j]] = h.get(v1_spins[j], 0.0) - variable1_sign*(2**j)/(2**variable_bit_count)
        h[v2_spins[j]] = h.get(v2_spins[j], 0.0) - variable2_sign*(2**j)/(2**variable_bit_count)
        # this may be flawed even in the "old" case, not a problem due to h being 0?? need to check

    for j in range(variable_bit_count):
        for k in range(variable_bit_count):
            coeff = variable1_sign*variable2_sign*\
                (2**(j+k))/(4**variable_bit_count)
            tmp1, tmp2 = v1_spins[j], v2_spins[k]
            if tmp1 == tmp2:
                c += coeff
            else:
                #J[tmp1, tmp2] += 0.5*coeff
                #J[tmp2, tmp1] += 0.5*coeff
                J[(tmp1, tmp2)] = J.get((tmp1, tmp2), 0.0) + 0.5*coeff
                J[(tmp2, tmp1)] = J.get((tmp2, tmp1), 0.0) + 0.5*coeff

    #return (np.pi**2)*J, (np.pi**2)*h, (np.pi**2)*c
    J = {key: (np.pi**2)*val for key, val in J.items()}
    h = {key: (np.pi**2)*val for key, val in h.items()}
    c = (np.pi**2)*c
    return J, h, c
    
def variable_times_sum_old(variable_idx, variable_sign, variable_bit_count, \
    variable_count, triplet_idx, triplet_count):

    t_spin = spin_index(triplet_idx, 'triplet', variable_count, \
        variable_bit_count, triplet_count)

    v_spins = tuple(
        spin_index(variable_idx, bit_idx, variable_count, variable_bit_count, \
            triplet_count) for bit_idx in range(variable_bit_count)
    )

    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    c += 3.0

    h[t_spin] -= 1
    for j in range(variable_bit_count):
        h[v_spins[j]] -= 3*variable_sign*(2**j)/(2**variable_bit_count)
        coeff = variable_sign*(2**j)/(2**variable_bit_count)
        J[v_spins[j], t_spin] += 0.5*coeff
        J[t_spin, v_spins[j]] += 0.5*coeff

    return (np.pi**2)*J, (np.pi**2)*h, (np.pi**2)*c

def variable_times_sum(variable_idx, variable_sign, variable_bit_count, \
    variable_count, triplet_idx, triplet_count):

    t_spin = spin_index(triplet_idx, 'triplet', variable_count, \
        variable_bit_count, triplet_count)

    v_spins = tuple(
        spin_index(variable_idx, bit_idx, variable_count, variable_bit_count, \
            triplet_count) for bit_idx in range(variable_bit_count)
    )

    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    c += 3.0

    #h[t_spin] -= 1
    h[t_spin] = h.get(t_spin, 0.0) - 1
    for j in range(variable_bit_count):
        #h[v_spins[j]] -= 3*variable_sign*(2**j)/(2**variable_bit_count)
        h[v_spins[j]] = h.get(v_spins[j], 0.0) - 3*variable_sign*(2**j)/(2**variable_bit_count)
        coeff = variable_sign*(2**j)/(2**variable_bit_count)
        #J[v_spins[j], t_spin] += 0.5*coeff
        #J[t_spin, v_spins[j]] += 0.5*coeff
        J[(v_spins[j], t_spin)] = J.get((v_spins[j], t_spin), 0.0) + 0.5*coeff
        J[(t_spin, v_spins[j])] = J.get((t_spin, v_spins[j]), 0.0) + 0.5*coeff

    #return (np.pi**2)*J, (np.pi**2)*h, (np.pi**2)*c
    J = {key: (np.pi**2)*val for key, val in J.items()}
    h = {key: (np.pi**2)*val for key, val in h.items()}
    c = (np.pi**2)*c
    return J, h, c

def sum_times_sum_old(variable_bit_count, variable_count, triplet_idx, \
    triplet_count):

    t_spin = spin_index(triplet_idx, 'triplet', variable_count, \
        variable_bit_count, triplet_count)

    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    c += 10.0
    h[t_spin] -= 6

    return (np.pi**2)*J, (np.pi**2)*h, (np.pi**2)*c

def sum_times_sum(variable_bit_count, variable_count, triplet_idx, \
    triplet_count):

    t_spin = spin_index(triplet_idx, 'triplet', variable_count, \
        variable_bit_count, triplet_count)

    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    c += 10.0
    #h[t_spin] -= 6
    h[t_spin] = h.get(t_spin, 0.0) - 6

    #return (np.pi**2)*J, (np.pi**2)*h, (np.pi**2)*c
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

def triplet_hamiltonian_old(variable1_idx, variable2_idx, variable3_idx, \
    variable1_sign, variable2_sign, variable3_sign, variable_count, \
    variable_bit_count, triplet_idx, triplet_count):
    
    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    # t1t1
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable1_idx, variable1_sign, variable1_sign, variable_bit_count, variable_count, triplet_count)
    J += Jterm
    h += hterm
    c += cterm

    # t2t2
    Jterm, hterm, cterm = variable_times_variable(variable2_idx, \
        variable2_idx, variable2_sign, variable2_sign, variable_bit_count, variable_count, triplet_count)
    J += Jterm
    h += hterm
    c += cterm

    # t3t3
    Jterm, hterm, cterm = variable_times_variable(variable3_idx, \
        variable3_idx, variable3_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    J += Jterm
    h += hterm
    c += cterm

    # t1t2
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable2_idx, variable1_sign, variable2_sign, variable_bit_count, variable_count, triplet_count)
    J += 2*Jterm
    h += 2*hterm
    c += 2*cterm

    # t1t3
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable3_idx, variable1_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    J += 2*Jterm
    h += 2*hterm
    c += 2*cterm

    # t2t3
    Jterm, hterm, cterm = variable_times_variable(variable2_idx, \
        variable3_idx, variable2_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    J += 2*Jterm
    h += 2*hterm
    c += 2*cterm

    #t1t
    Jterm, hterm, cterm = variable_times_sum(variable1_idx, variable1_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    J -= 2*Jterm
    h -= 2*hterm
    c -= 2*cterm

    #t2t
    Jterm, hterm, cterm = variable_times_sum(variable2_idx,  variable2_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    J -= 2*Jterm
    h -= 2*hterm
    c -= 2*cterm

    #t3t
    Jterm, hterm, cterm = variable_times_sum(variable3_idx,  variable3_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    J -= 2*Jterm
    h -= 2*hterm
    c -= 2*cterm

    #tt
    Jterm, hterm, cterm = sum_times_sum(variable_bit_count, variable_count, \
        triplet_idx, triplet_count)
    J += Jterm
    h += hterm
    c += cterm

    return J, h, c

def triplet_hamiltonian(variable1_idx, variable2_idx, variable3_idx, \
    variable1_sign, variable2_sign, variable3_sign, variable_count, \
    variable_bit_count, triplet_idx, triplet_count):
    
    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    # t1t1
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable1_idx, variable1_sign, variable1_sign, variable_bit_count, variable_count, triplet_count)
    #J += Jterm
    #h += hterm
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += cterm

    # t2t2
    Jterm, hterm, cterm = variable_times_variable(variable2_idx, \
        variable2_idx, variable2_sign, variable2_sign, variable_bit_count, variable_count, triplet_count)
    #J += Jterm
    #h += hterm
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += cterm

    # t3t3
    Jterm, hterm, cterm = variable_times_variable(variable3_idx, \
        variable3_idx, variable3_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    #J += Jterm
    #h += hterm
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += cterm

    # t1t2
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable2_idx, variable1_sign, variable2_sign, variable_bit_count, variable_count, triplet_count)
    #J += 2*Jterm
    #h += 2*hterm
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += 2*cterm

    # t1t3
    Jterm, hterm, cterm = variable_times_variable(variable1_idx, \
        variable3_idx, variable1_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    #J += 2*Jterm
    #h += 2*hterm
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += 2*cterm

    # t2t3
    Jterm, hterm, cterm = variable_times_variable(variable2_idx, \
        variable3_idx, variable2_sign, variable3_sign, variable_bit_count, variable_count, triplet_count)
    #J += 2*Jterm
    #h += 2*hterm
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += 2*cterm

    #t1t
    Jterm, hterm, cterm = variable_times_sum(variable1_idx, variable1_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    #J -= 2*Jterm
    #h -= 2*hterm
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_diff(J, Jterm)
    h = merge_dict_by_diff(h, hterm)
    c -= 2*cterm


    #t2t
    Jterm, hterm, cterm = variable_times_sum(variable2_idx,  variable2_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    #J -= 2*Jterm
    #h -= 2*hterm
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_diff(J, Jterm)
    h = merge_dict_by_diff(h, hterm)
    c -= 2*cterm

    #t3t
    Jterm, hterm, cterm = variable_times_sum(variable3_idx,  variable3_sign, \
        variable_bit_count, variable_count, triplet_idx, triplet_count)
    #J -= 2*Jterm
    #h -= 2*hterm
    Jterm = {key:2*val for key, val in Jterm.items()}
    hterm = {key:2*val for key, val in hterm.items()}
    J = merge_dict_by_diff(J, Jterm)
    h = merge_dict_by_diff(h, hterm)
    c -= 2*cterm

    #tt
    Jterm, hterm, cterm = sum_times_sum(variable_bit_count, variable_count, \
        triplet_idx, triplet_count)
    #J += Jterm
    #h += hterm
    J = merge_dict_by_sum(J, Jterm)
    h = merge_dict_by_sum(h, hterm)
    c += cterm

    return J, h, c

def hamiltonian_old(variable_count, variable_bit_count, triplet_refls, \
    triplet_signs, triplet_weights=None, verbose=False):
    if triplet_weights is None:
        triplet_weights = [1.0 for t in triplet_refls]

    triplet_count = len(triplet_refls)

    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    for j, triplet in enumerate(triplet_refls):
        if verbose:
            print(f"Processing triplet {j+1} of {triplet_count}...   ", end="\r")
        triplet_idx = j
        variable1_idx, variable2_idx, variable3_idx = tuple(triplet)
        variable1_sign, variable2_sign, variable3_sign = tuple(triplet_signs[j])

        Jt, ht, ct = triplet_hamiltonian(variable1_idx, variable2_idx, \
            variable3_idx, variable1_sign, variable2_sign, variable3_sign, \
            variable_count, variable_bit_count, triplet_idx, triplet_count)

        J += triplet_weights[j]*Jt
        h += triplet_weights[j]*ht
        c += triplet_weights[j]*ct

    return J, h, c

def hamiltonian(variable_count, variable_bit_count, triplet_refls, \
    triplet_signs, triplet_weights=None, verbose=False):
    if triplet_weights is None:
        triplet_weights = [1.0 for t in triplet_refls]

    triplet_count = len(triplet_refls)

    J, h, c = zero_hamiltonian(variable_count, variable_bit_count, \
        triplet_count)

    for j, triplet in enumerate(triplet_refls):
        if verbose:
            print(f"Processing triplet {j+1} of {triplet_count}...   ", end="\r")
        triplet_idx = j
        variable1_idx, variable2_idx, variable3_idx = tuple(triplet)
        variable1_sign, variable2_sign, variable3_sign = tuple(triplet_signs[j])

        Jt, ht, ct = triplet_hamiltonian(variable1_idx, variable2_idx, \
            variable3_idx, variable1_sign, variable2_sign, variable3_sign, \
            variable_count, variable_bit_count, triplet_idx, triplet_count)

        #J += triplet_weights[j]*Jt
        #h += triplet_weights[j]*ht
        Jt = {key:triplet_weights[j]*val for key, val in Jt.items()}
        ht = {key:triplet_weights[j]*val for key, val in ht.items()}
        J = merge_dict_by_sum(J, Jt)
        h = merge_dict_by_sum(h, ht)
        c += triplet_weights[j]*ct

    J = {key:val for key, val in J.items() if not (val == 0.0) }
    h = {key:val for key, val in h.items() if not (val == 0.0) }

    return J, h, c

