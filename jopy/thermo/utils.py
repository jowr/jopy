# -*- coding: utf-8 -*-
from __future__ import print_function, division

import numpy as np
import warnings

class UnitError(ValueError):
    pass

def check_T(T):
    """Check for valid temperature
    
    Parameters
    ----------
        T : float or numpy.array in Kelvin [K]
        
    Returns
    -------
    bool
        True if valid, else False
        
    """
    if np.any(T<0.0): raise UnitError("negative temperature: "+str(T))
    if np.any(T>1e4): raise UnitError("too high temperature: "+str(T))
    
    if np.any(T<250.0): warnings.warn("very low temperature: "+str(T))
    if np.any(T>1.0e3): warnings.warn("very high temperature: "+str(T))
    
    return True


def check_p(p):
    """Check for valid pressure
    
    Parameters
    ----------
        p : float or numpy.array in Pascal [Pa]
        
    Returns
    -------
    bool
        True if valid, else False
        
    """
    if np.any(p<0.0): raise UnitError("negative pressure: "+str(p))
    if np.any(p>1e8): raise UnitError("too high pressure: "+str(p))
    
    if np.any(p<100.0): warnings.warn("very low pressure: "+str(p))
    if np.any(p>1.0e7): warnings.warn("very high pressure: "+str(p))
    
    return True 

    