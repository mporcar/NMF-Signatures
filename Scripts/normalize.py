import itertools as itt
from collections import defaultdict


CB = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


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

def mut_key_generator():
    """

    Returns:
        Yields all tuples representing lexicographic sortable contexts
        1st component: substitution
        2nd component: flanks

    """

    subs = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
    for s in sorted(subs):
        for c in itt.product({'A', 'C', 'G', 'T'}, repeat=2):
            yield tuple([s, ''.join(c)])



def normalize_profile(profile, tnt_abundance):
    """

    Args:
        profile: dict: count for each context
        tnt_abundance: dict: count for each possible tnt

    Returns:
        profile that would have resulted from uniform tri-nucleotide content

    """

    abundance = {}
    for b in {'C', 'T'}:
        all_pairs = iter([''.join(l) for l in itt.product({'A', 'C', 'G', 'T'}, repeat=2)])
        for p in all_pairs:
            tnt = p[0] + b + p[1]
            abundance[tnt] = tnt_abundance[tnt] + tnt_abundance[reverse_complement(tnt)]
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
