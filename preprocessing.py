def friedel_standard(r):
    rl = [int(x) for x in r.strip('[').strip(']').split(',')]
    if (rl[0] == 0) and (rl[2] == 0):
        return rl[1] >= 0
    elif (rl[2] == 0):
        return rl[0] >= 0
    else:
        return rl[2] >= 0
    raise Exception("Shouldn't have reached here")

def friedel_invert(r):
    rl = [int(x) for x in r.strip('[').strip(']').split(',')]
    rl = [-1*x for x in rl]
    r = ''.join(str(rl).split())
    return r

def friedel_standardise(r):
    if not friedel_standard(r):
        return friedel_invert(r)
    return r

def process_triplets(triplets):
    triplet_refls = []
    triplet_signs = []
    refl_to_int = {}
    curr_int = -1
    triplets_used = []
    for t_idx, t in enumerate(triplets):
        triplet = []
        signs = []
        for r in t:
            rstand = friedel_standardise(r)
            if rstand == r:
                sign = 1
            else:
                sign = -1
            signs.append(sign)
            try:
                refl_int = refl_to_int[rstand]
            except KeyError:
                curr_int += 1
                refl_to_int[rstand] = curr_int
                refl_int = curr_int
            triplet.append(refl_int)
        if not triplet in triplet_refls:
            triplet_refls.append(triplet)
            triplet_signs.append(signs)
            triplets_used.append(t_idx)
    int_to_refl = {value:key for key, value in refl_to_int.items()}
    return triplet_refls, triplet_signs, refl_to_int, int_to_refl, triplets_used

    