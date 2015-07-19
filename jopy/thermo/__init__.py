# -*- coding: utf-8 -*-
from __future__ import print_function, division

import matplotlib.pyplot as plt

from jopy.thermo.utils import check_T

import numpy as np
from numpy import empty_like
from jopy.utils import transition_factor





def eta_carnot(T_cold,T_hot):
    check_T(T_cold)
    check_T(T_hot)
    return 1. - T_cold / T_hot

def __lmtd(Delta_T1, Delta_T2):
    return (Delta_T1 - Delta_T2) / np.log(Delta_T1 / Delta_T2)

def _lmtd_libpf(Delta_T1, Delta_T2, deadBand=0.1):
    """Logarithmic mean temperature difference
    
    This function returns the logarithmic mean temperature 
    difference as calculated from temperature differences at
    both end. It is a relatively robust solution that also
    is used in libpf for process simulations.
    
    Parameters:
    -----------
        Delta_T1 : float or numpy.array
            Temperature difference on one side of the heat exchanger        
        Delta_T2 : float or numpy.array
            Temperature difference on the other side of the heat exchanger 
        deadBand : float or numpy.array
            Absolute temperature difference that is considered to be 0.0 K
        
    Returns:
    --------
        float or numpy.array
            The calculated logarithmic mean temperature difference
    
    """    
    same = np.abs(Delta_T1-Delta_T2) <= 0.1
    if (np.abs(Delta_T1) <= deadBand) or same:             # small first value
        if (np.abs(Delta_T2) <= deadBand) or same:         # small second value
            LMTD = (Delta_T1 + Delta_T2) / 2.0
        else:                                              # normal second value
            Delta_T1clip = deadBand*Delta_T2 / np.abs(Delta_T2) #{absolute value is set to deadBand, sign is the same as Delta_T2}
            LMTD = ((Delta_T1clip - Delta_T2) / (Delta_T1 - Delta_T2)) * __lmtd(Delta_T1clip,Delta_T2)
    else:
        if (np.abs(Delta_T2) <= deadBand):
            Delta_T2clip = deadBand*Delta_T1/np.abs(Delta_T1)
            LMTD = ((Delta_T1 - Delta_T2clip) / (Delta_T1 - Delta_T2)) * __lmtd(Delta_T1,Delta_T2clip)
        else:
            if ((Delta_T1 * Delta_T2) <= 0.0):             # The deltas have different signs 
                if (np.abs(Delta_T1) <= np.abs(Delta_T2)): # smaller first value
                    Delta_T1clip = deadBand*Delta_T2/np.abs(Delta_T2)
                    LMTD = ((Delta_T1clip - Delta_T2) / (Delta_T1 - Delta_T2)) * __lmtd(Delta_T1clip,Delta_T2)
                else:                                      # smaller second value
                    Delta_T2clip = deadBand*Delta_T1/np.abs(Delta_T1)
                    LMTD = ((Delta_T1 - Delta_T2clip) / (Delta_T1 - Delta_T2)) * __lmtd(Delta_T1,Delta_T2clip)
            else:
                LMTD = __lmtd(Delta_T1,Delta_T2)
                
    return LMTD


def _lmtd_jopy(Delta_T1, Delta_T2, dead_band=0.1):
    """Logarithmic mean temperature difference
    
    This function returns the logarithmic mean temperature 
    difference as calculated from temperature differences at
    both end. 
    
    Parameters:
    -----------
        Delta_T1 : float or numpy.array
            Temperature difference on one side of the heat exchanger        
        Delta_T2 : float or numpy.array
            Temperature difference on the other side of the heat exchanger 
        dead_band : float or numpy.array
            Absolute temperature difference that is considered to be 0.0 K
        
    Returns:
    --------
        float or numpy.array
            The calculated logarithmic mean temperature difference
    
    """
    
    res            = empty_like(Delta_T1)
    mask           = np.abs(Delta_T1) < np.abs(Delta_T2) 
    res[mask]      = np.sign(Delta_T2[mask])    
    Delta_T1[mask] = res[mask] * np.maximum(res[mask] * Delta_T1[mask], dead_band)
    mask           = np.logical_not(mask)
    res[mask]      = np.sign(Delta_T1[mask])
    Delta_T2[mask] = res[mask] * np.maximum(res[mask] * Delta_T2[mask], dead_band)
    
    # Now the inputs are sanitised
    mask = np.abs(Delta_T1-Delta_T2) <= dead_band
    res[mask] = (Delta_T1[mask] + Delta_T2[mask]) / 2.0
    #res[mask] = transition_factor(-dead_band, dead_band, Delta_T1-Delta_T2[mask])
    mask = np.logical_not(mask)
    res[mask] = __lmtd(Delta_T1[mask], Delta_T2[mask])
    return res
    

lmtd = _lmtd_jopy 
"""Please see :py:func:`._lmtd_jopy` for documentation"""

if __name__ == "__main__":
    
    var = np.linspace(-20, 20)
    con = np.zeros_like(var)+10
    plt.plot(var,lmtd(var,con),'o')
    plt.show()
    
    #print(np.vectorize(_lmtd_libpf)([0,1,5,10,10,10,10],[10,10,10,10,5,1,0]))
    #print(np.vectorize(_lmtd_libpf)([-0,-1,-5,-10,-10,-10,-10],[10,10,10,10,5,1,0]))
    
    #res = _lmtd_libpf(10.0, 10.0)
    #print(res,10.0)
