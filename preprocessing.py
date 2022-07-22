import numpy as np

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

def process_triplets(triplets, prune=False):
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

    counts_dict = {}
    for j, trs in enumerate(triplet_refls):
        for k, r in enumerate(trs):
            counts_dict[r] = counts_dict.get(r, 0) + 1

    if prune:
        triplet_refls_pruned, counts_dict_pruned, int_to_prunedint, prunedint_to_int = prune_triplets(triplet_refls, counts_dict)
        triplets_used_pruned, triplet_signs_pruned = [], []
        for j, x in enumerate(triplets_used):
            #print(x)
            if triplet_refls[j] in triplet_refls_pruned:
                triplets_used_pruned.append(x)
                triplet_signs_pruned.append(triplet_signs[j])
        triplets_used, triplet_signs = triplets_used_pruned, triplet_signs_pruned
        triplet_refls = [[int_to_prunedint[r] for r in t] for t in triplet_refls_pruned]
        counts_dict = {int_to_prunedint[key]:val for key, val in counts_dict.items() if key in int_to_prunedint.keys()}
        int_to_refl = {int_to_prunedint[key]:val for key, val in int_to_refl.items() if key in int_to_prunedint.keys()}
        refl_to_int = {key:int_to_prunedint[val] for key, val in refl_to_int.items() if val in int_to_prunedint.keys()}

    return triplet_refls, triplet_signs, counts_dict, refl_to_int, int_to_refl, triplets_used

def prune_triplets(triplet_refls, counts_dict, verbose=False):
    # need to deal with triplets_used somehow
    total_triplets = len(triplet_refls)
    loopidx = -1
    while True:
        loopidx += 1
        if verbose:
            print(f"{len(triplet_refls)} of {total_triplets} triplets remain on loop iteration {loopidx}   ", end="\r" )
        if np.all( ~(np.array(list(counts_dict.values())) == 1) ):
            break
        triplet_refls_pruned = []
        counts_dict_pruned = dict(counts_dict)
        for tr in triplet_refls:
            r0, r1, r2 = tuple(tr)
            if (counts_dict[r0] <= 1) or (counts_dict[r1] <= 1) or (counts_dict[r2] <= 1):
                counts_dict_pruned[r0] -= 1
                counts_dict_pruned[r1] -= 1
                counts_dict_pruned[r2] -= 1
            else:
                triplet_refls_pruned.append(tr)
        counts_dict = counts_dict_pruned
        triplet_refls = triplet_refls_pruned
    counts_dict = {key:val for key, val in counts_dict.items() if not val == 0}

    newint_to_oldint = {}
    for newint, oldint in enumerate(list(counts_dict.keys())):
        newint_to_oldint[newint] = oldint
    oldint_to_newint = {val:key for key, val in newint_to_oldint.items()}

    return triplet_refls, counts_dict, oldint_to_newint, newint_to_oldint