import numpy as np

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

