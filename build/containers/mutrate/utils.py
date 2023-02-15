from itertools import product
from collections import defaultdict
import numpy as np

CB = dict(zip(list('ACGT'), list('TGCA')))


def mut_key_generator():
    """
    Returns:
        Generates all possible lexicographic pairs
            1st component: substitution;
            2nd component: flanks
    """
    subs = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
    for s in sorted(subs):
        for c in sorted(product({'A', 'C', 'G', 'T'}, repeat=2)):
            yield tuple([s, ''.join(c)])


def shortkey_to_lex(key):
    """
    Args:
        key: signature key in short key format
    Returns:
        tuple representing a lexicographic pair key
    """
    t1 = ''.join([key[1], key[-1]])
    t2 = ''.join([key[0], key[2]])
    return tuple([t1, t2])


def complementary(key):
    """
    Args:
        key: short key format
    Returns:
        DNA complementary in short key format
    """
    ref = key[1]
    alt = key[-1]
    ctxt = (key[0], key[2])
    return ''.join([CB[ctxt[1]], CB[ref], CB[ctxt[0]], '>', CB[alt]])


def purine_to_pyrimidine(m):
    """
    Args:
        m: tuple of str: mutation key: substitution and flanks
    Returns:
        reverse complement of mutation profile m
    """
    t1 = CB[m[0][0]] + CB[m[0][1]]
    t2 = CB[m[1][1]] + CB[m[1][0]]
    return tuple([t1, t2])


def maf_to_makeup(annotmuts):
    """
    Args:
        annotmuts: dataframe (maf)
    Returns:
        From maf it creates dict {sample --> {context_key --> count}}
    """
    makeup = defaultdict(lambda: defaultdict(int))
    for i, row in annotmuts.iterrows():
        ref1 = row['ref']
        ref3 = row['ref3_cod']
        alt1 = row['mut']
        sample = row['sampleID']
        if (alt1 in set('ACGT')) and (ref1 in set('ACGT')):
            key = (ref1 + alt1, ref3[0] + ref3[-1])
            if ref1 not in set('CT'):
                key = purine_to_pyrimidine(key)
            makeup[sample][key] += 1
    return makeup


def lex_to_sigfit(mk):
    """
    Args:
        mk: tuple of str: sub and flanks
    Returns:
        str: context in sigfit notation
    """
    return mk[1][0] + mk[0][0] + mk[1][1] + '>' + mk[1][0] + mk[0][1] + mk[1][1]


def reverse_complement(seq):
    """
    Args:
        seq: str
    Returns:
        reverse complement of seq
    """
    res = ''
    for a in seq[::-1]:
        res += CB[a]
    return res


def normalize_profile(profile, triplets):
    """
    Args:
        profile: count for each context
        triplets: count for each possible tnt
    Returns:
        profile that would have resulted from equally abundant triplets
    """
    abundance = {}
    for b in {'C', 'T'}:
        all_pairs = iter([''.join(l) for l in product({'A', 'C', 'G', 'T'}, repeat=2)])
        for p in all_pairs:
            tnt = p[0] + b + p[1]
            try:
                abundance[tnt] = triplets[tnt] + triplets[reverse_complement(tnt)]
            except KeyError as err:
                abundance[tnt] = 1
    cond_prob = defaultdict(float)
    for mut_type in mut_key_generator():
        ref_tnt = mut_type[1][0] + mut_type[0][0] + mut_type[1][1]
        if abundance[ref_tnt] != 0:
            cond_prob[mut_type] = profile[mut_type] / abundance[ref_tnt]
        else:
            cond_prob[mut_type] = 0
    total_cond_prob = sum(cond_prob.values())
    norm_signature = defaultdict(int)
    for mut_type in mut_key_generator():
        norm_signature[mut_type] = cond_prob[mut_type] / total_cond_prob
    return norm_signature


def denorm_profile(profile, triplets):
    """
    Given a normalised makeup, it produces the expected relative frequencies in a
    new genomic element given its trinucleotide relative abundance. Both inputs
    are expected to contain relative frequencies -- in particular, their values
    must add up to 1.
    """
    coef_matrix = np.zeros((96+1, 96))
    profile_vector = np.array([profile[k] for k in mut_key_generator()])
    total = sum(triplets.values())
    tnt_relative_abundance = {k: v / total for k, v in triplets.items()}
    abundance_dict = defaultdict(int)
    for mtype in mut_key_generator():
        abundance_dict[mtype] = tnt_relative_abundance[mtype[1][0] + mtype[0][0] + mtype[1][1]]
    abundance_vector = np.array([abundance_dict[k] for k in mut_key_generator()])
    for i in range(96):
        for j in range(96):
            if i == j:
                if abundance_vector[j] != 0:
                    coef_matrix[i][j] = (profile_vector[j] - 1) / abundance_vector[j]
                else:
                    coef_matrix[i][j] = 1
            else:
                coef_matrix[i][j] = profile_vector[i] / abundance_vector[j]
    coef_matrix[96] = np.ones(96)  # this is the (96+1)-th row
    coef_matrix = np.delete(coef_matrix, (0), axis=0)
    b = np.zeros(96)
    b[95] = 1
    x = np.linalg.solve(coef_matrix, b)
    return dict(zip(list(mut_key_generator()), x))


def denorm_subs(profile, subs_abundance):
    """
    Given a normalized profile, it produces the expected relative frequencies
    given a new genomic scope defined by site abundances.
    """
    coef_matrix = np.zeros((96 + 1, 96))
    profile_vector = np.array([profile[k] for k in mut_key_generator()])
    total_subs = sum(subs_abundance.values())
    subs_rel_abundance = {k: v / total_subs for k, v in subs_abundance.items()}
    abundance_vector = np.array([subs_rel_abundance[k] for k in mut_key_generator()])
    for i, j in product(range(96), repeat=2):
        if i == j:
            if abundance_vector[j] != 0:
                coef_matrix[i][j] = (profile_vector[j] - 1) / abundance_vector[j]
            else:
                coef_matrix[i][j] = 1
        else:
            coef_matrix[i][j] = profile_vector[i] / abundance_vector[j]
    coef_matrix[96] = np.ones((96))  # this is the (96 + 1)-th row
    coef_matrix = np.delete(coef_matrix, (0), axis=0)
    b = np.zeros(96)
    b[95] = 1
    x = np.linalg.solve(coef_matrix, b)
    return dict(zip(list(mut_key_generator()), x))

