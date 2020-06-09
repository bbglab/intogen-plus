import operator

import numpy as np

# Cython version
from intogen_combination.schulze_strongest_path_cython import strongest_path

# Python version
#from intogen_combination.schulze_strongest_path import strongest_path


def combination_ranking(ballot_dict, weights):

    # INIT
    all_candidates = list(set.union(*[set(ballot_dict[voter].keys()) for voter in ballot_dict]))
    all_candidates_to_idx = {k: i for i, k in enumerate(all_candidates)}
    size = len(all_candidates)

    pref = np.zeros(size**2, dtype=np.float64)  #[0.0] * (size**2)
    spath = np.zeros(size**2, dtype=np.float64)  #[0.0] * (size**2)

    # PREPARE
    if weights is None:
        weights = {}
        for voter in ballot_dict:
            weights[voter] = 1.
        weights = dict(weights)

    for voter in ballot_dict:
        d = ballot_dict[voter]
        for i in all_candidates:
            if i not in d.keys():
                for j in d:
                    pref[all_candidates_to_idx[j]*size + all_candidates_to_idx[i]] += weights[voter]
            else:
                r = d[i]
                for j in d:
                    if d[j] < r:
                        pref[all_candidates_to_idx[j]*size + all_candidates_to_idx[i]] += weights[voter]

    # STRONGEST PATH
    strongest_path(size, pref, spath)

    # COMBINATION RANKING
    scores_dict = {}
    for i in range(size):
        score = 0
        for j in range(size):
            if spath[i*size + j] < spath[j*size + i]:
                score += 1
        scores_dict[all_candidates[i]] = score
    sorted_scores = sorted(scores_dict.items(), key=operator.itemgetter(1), reverse=True)

    ranking = {}
    prev_score = None
    prev_rank = None
    counter = 1
    while len(sorted_scores) > 0:
        c = sorted_scores.pop()
        if prev_score is None:
            ranking[c[0]] = 1
            prev_score = c[1]
            prev_rank = 1
        elif prev_score == c[1]:
            ranking[c[0]] = prev_rank
        elif prev_score < c[1]:
            ranking[c[0]] = counter
            prev_score = c[1]
            prev_rank = counter
        counter += 1

    return ranking
