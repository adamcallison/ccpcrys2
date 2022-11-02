import numpy as np
import scipy.optimize as spo

def state_to_ints(state, variable_bit_count, triplet_count):
    sum_piece = state[-2*triplet_count:][::-1]
    sum_piece_processed = []
    for j in range(triplet_count):
        if (sum_piece[(2*j)], sum_piece[(2*j)+1]) == (1, 1):
            sum_piece_processed.append(0)
        elif (sum_piece[(2*j)], sum_piece[(2*j)+1]) == (1, -1):
            sum_piece_processed.append(1)
        elif (sum_piece[(2*j)], sum_piece[(2*j)+1]) == (-1, 1):
            sum_piece_processed.append(2)
        elif (sum_piece[(2*j)], sum_piece[(2*j)+1]) == (-1, -1):
            sum_piece_processed.append(3)
        else:
            raise ValueError
    sum_piece = sum_piece_processed
    state = state[:-2*triplet_count]

    assert len(state) % variable_bit_count == 0

    variable_count = len(state) // variable_bit_count

    variable_ints = []

    for variable_idx in range(variable_count):
        start = variable_idx*variable_bit_count
        end = (variable_idx+1)*variable_bit_count
        variable_state = state[start:end]
        variable_int = 0
        for j, spin in enumerate(variable_state):
            bit = (0 if spin == 1 else 1)
            variable_int += ((2**j)*bit)
        variable_ints.append(variable_int)
    
    return np.array(variable_ints), np.array(sum_piece)

def state_to_phases(state, variable_bit_count, triplet_count):
    variable_ints, sum_piece = state_to_ints(state, variable_bit_count, \
        triplet_count)

    N = 2**variable_bit_count

    phases = 2*np.pi*( (variable_ints/N) + (1/(2*N)))

    sum_targets = 2*np.pi*sum_piece

    return phases, sum_targets

def phase_shift(int_to_refl, translation):
    phase_shifts = np.zeros(len(int_to_refl))
    for j in range(len(int_to_refl)):
        refl = np.array([int(x.strip()) for x in ( int_to_refl[j][1:-1].split(',') ) ])
        phase_shift_val = np.dot(refl, translation) 
        phase_shifts[j] = phase_shift_val
    return phase_shifts

def standardize_phases(phases):
    phase_fracs = phases / (2*np.pi)
    phase_ints = np.floor(phase_fracs)
    new_phases_fracs = phase_fracs - phase_ints
    new_phases = new_phases_fracs*(2*np.pi)
    return new_phases

def phase_match(source_phases, target_phases, int_to_refl):
    def func(translation):
        phase_shifts = phase_shift(int_to_refl, translation)
        shifted_phases = standardize_phases(source_phases + phase_shifts)
        #shifted_phases = (source_phases + phase_shifts)
        diff = np.sum(np.abs(shifted_phases - target_phases))/len(phase_shifts)*(180/np.pi)
        return diff
    #x0 = np.zeros(3)
    #res = spo.minimize(func, x0)
    bounds = ((-1000000.0,1000000.0),)*3
    res = spo.shgo(func, bounds, iters=6, sampling_method='simplicial')
    translation = res.x
    diff = res.fun
    success = res.success
    phase_shifts = phase_shift(int_to_refl, translation)
    shifted_phases = source_phases + phase_shifts
    #print(res)
    return shifted_phases, phase_shifts, translation, diff, success

