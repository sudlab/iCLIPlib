''' This module contains a varety of diffrent functions for calculating
distance between two profiles'''

from .utils import spread
import numpy as np


##################################################
def calcAverageDistance(profile1, profile2):
    ''' This function calculates the average distance of all
    pairwise distances in two profiles
    
    Parameters
    ----------
    profile1, profile2 : pandas.Series
        Two series to calculate average distance between
    
    Returns
    -------
    int
        Average distance between each pairwise combination of
        cross-links from the two profiles.
        
    '''

    def _cartesian(x, y):
        return np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])

    positions = _cartesian(profile1.index.values, profile2.index.values)
    counts = _cartesian(profile1.values, profile2.values)
    counts = np.prod(counts, axis=1)
    distances = np.abs(positions[:, 0] - positions[:, 1])
    mean_distance = (distances.astype("float64") *
                     counts).sum() / np.sum(counts)

    return mean_distance


##################################################
def findMinDistance(profile1, profile2):
    '''Finds mean distance between each cross-link in profile1
    and the closest cross-link in profile2
    
    Parameters
    ----------
    profile1, profile2 : pandas.Series of int
        Profiles to calculate distances on. Index is expected to be genome
        or transcriptome coordinates.
        
    Returns
    -------
    float
        Mean of the distance between each cross-link in profile1 and the
        closet cross-link in profile2_ready
    
    Note
    ----
    In order that `findMinDistance(profile1, profile2) == findMinDistance(
    profile2, profile1)`, the number of cross-links on each site is
    ignored.
    
    '''

    locations1 = profile1.index.values

    locations2 = profile2.index.values  # .astype("int16")

    mat1 = np.repeat(locations1, locations2.size).reshape(
        (locations1.size, locations2.size))
    mat2 = np.tile(locations2, locations1.size).reshape(
        (locations1.size, locations2.size))

    distances = np.abs(mat1-mat2).min(axis=1)

    return distances.mean()


##################################################
def corr_profile(profile1, profile2, nspread, profile2_ready=False):
    """Calculate the spearmans correlation coefficient btween two
    (possibly extended) cross-link profiles.
    
    Parameters
    ----------
    profile1, profile2 : pandas.Series of int
        Cross-link profiles. Index is positions, value is counts.
    nspread : int
        Number of bases to extend each profile in each direction
    profile2_ready : bool, optional
        Has profile2 already been reindexed and spread (see below)
    
    Return
    ------
    float
        Spearmans correlation coefficient between the two profiles
        
    Notes
    -----
    During each call profile is reindexed so that bases with no cross-links
    are included in the index, and cross-links sites are extended.
    
    As this is a slow process, and it is imagined that if randomisaitons
    are applied to profile1, but profile2 is held constant, it would be
    more efficent to apply this only once to profile2 and compare it to 
    many randomisaitons of profile1. If `profile2_read=True`, it is assumed
    that `profile2` is already supplied with 0 count positions indexed and
    cross-links extended. The can be achieved with:
        
        profile2 = profile2.reindex(range(start, end))
        profile2 = profile2.fillna(x)
        profile2 = profile2.spread(profile2, nspread)
        
    A known flaw in this measure of the relationship between two profiles
    is that if profiles contain many 0s, even if the 0s are in the same 
    positions, the correlation will be low. 
    """
    
    profile1 = profile1.reindex(
                    range(int(profile1.index.values.min())-1,
                          int(profile1.index.values.max())+1)).fillna(0)
    profile1 = spread(profile1, nspread, False)
        
    if not profile2_ready:
        profile2 = profile2.reindex(
                       range(int(profile2.index.values.min()),
                             int(profile2.index.values.max()))).fillna(0)
        profile2 = spread(profile2, nspread, False)
        
    return profile1.corr(profile2, method="spearman")
