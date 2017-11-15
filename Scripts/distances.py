from scipy.stats import multinomial
import numpy as np
import math
from scipy.spatial.distance import euclidean
 

def hellinger(p, q):
    return np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)) / np.sqrt(2)

def kullback_leibler(profile_1, profile_2):
    """
    Args:
        profile_1: profile
        profile:2: profile

    Returns:
        Kullback-Leibler divergence D_KL(p || q)
    """
    p_array = np.array([profile_1[k] for k in sorted(profile_1)])
    q_array = np.array([profile_2[k] for k in sorted(profile_2)])
    kl = 0.0
    for i in range(len(p_array)):
        if p_array[i] < 0.0002:
            continue
        else:
            if q_array[i] < 0.0002:
                continue
            else:
                kl += p_array[i] * (math.log(p_array[i]/q_array[i]))
    return kl


def battacharyya(profile_1, profile_2):
    """

    Args:
        profile_1: profile
        profile_2: profile

    Returns:
        battacharya divergence

    """

    #p_array = np.array([np.sqrt(profile_1[k]) for k in sorted(profile_1)])
    #q_array = np.array([np.sqrt(profile_2[k]) for k in sorted(profile_2)])
    p_array = np.array([np.sqrt(k) for k in profile_2])
    q_array = np.array([np.sqrt(k) for k in profile_1])
    return(1 - np.dot(p_array, q_array))

def kullback_leibler_distance(profile_1, profile_2):
    """

    Args:
        profile_1: profile
        profile:2: profile

    Returns:
        Kullback-Leibler distance (D_KL(p || q)+D_KL(q || p))/2
    """
    return((kullback_leibler(profile_1,profile_2)+kullback_leibler(profile_2,profile_1))/2)