# -*- coding: utf-8 -*-
from __future__ import print_function, division

import numpy as np
from numpy import pi,e

def module_class_dict(mod):
    """Get all classes from a module in a dict with the names
     
    Parameters
    ----------
    mod : Python module
        The module to extract the classes from.
         
    """
    mod_name = str(mod).split("'")[1]
    ret = {}
    for name, cls in mod.__dict__.items():
        if isinstance(cls, type):
            if mod_name in str(cls):
                ret[name] = cls
    return ret
    
# def module_class_dict(mod):
#     """
#     Returns a list of names of the abstract base
#     classes that should not be instantiated. Can
#     be used to build an ignore list.
#     """
#     mod_name = str(mod).split("'")[1]
#     print(mod_name)
#     ignList = {}
#     for i in inspect.getmembers(mod):
#         if inspect.isclass(i[1]):
#             ignList[i[0]] = i[1]
#     return ignList

    #dict([(name, cls) for name, cls in mod.__dict__.items() if isinstance(cls, type)])

def transition_factor(start=0.25, stop=0.75, position=0.5, order=2):
    """Weighting factor for smooth transition (from 0 to 1)
    
    This function returns a value between 0 and 1. A smooth transition
    is achieved by means of defining the position and the transition
    interval from start to stop parameter. Outside this interval,
    the 0 and the 1 are returned, respectively. This transition function
    with up to two smooth derivatives was taken from [1]. If you provide
    an order higher than 2, the generalised logistic function [2] will be
    used to calculated the transition curve.

    Parameters:
    -----------
    start : float
        start of transition interval; default 0.25
    stop : float 
        end of transition interval; default 0.75
    position : float 
        current position; default 0.5
    order : integer
        Smooth up to which derivative?; default 2

    Returns:
    --------
    float : 
        smooth transition between 0 and 1 from start to stop [-]
        
    
    Use tFactor in an equation like this:
    tFactor = transition_factor(start=start,stop=stop,position=position);
    smoothed = tFactor*1stValue + (1 - tFactor)*2ndValue;

    [1] Christoph C Richter, Proposal of New Object-Oriented Equation-Based Model
        Libraries for Thermodynamic Systems, PhD thesis, Technical University
        Carolo-Wilhelmina Braunschweig, 2008
    [2] Generalised logistic function on http://en.wikipedia.org/wiki/Generalised_logistic_function
    """
    a_map = [-1./2., -2./pi, -3./4., -8./pi] #First parameters
    b_map = [ 1./2.,  1./2.,  1./2.,  1./2.] #Second parameters

    #Rename variables to match with Richter2008, p.68ff
    phi = 0.0 #"current phase";
    a = 0.0 # "multiplier";
    b = 0.0 # "addition";
    x_t = 0.0 # "Start of transition";
    x = 0.0 # "Current position";
    DELTAx = 0.0 # "Length of transition";

    #Parameters for generalised logistic function
    A = 0. #"Lower asymptote";
    K = 1. #"Upper asymptote";
    B = 8. #"Growth rate";
    nu= 1. #"Symmetry changes";
    Q = nu #"Zero correction";
    M = nu #"Maximum growth for Q = nu";
    X = 0.
    END =     0.
    START =   0.
    factor =  0.

    order = int(order)
    if order < 0:
        raise ValueError("This function only supports positive values for the order of smooth derivatives.")

    swapper = None
    if start>stop:
        swapper = start
        start = stop
        stop = swapper
        #raise Exception("There is only support for positive differences, please provide start < stop.")

    position = np.array(position)
    res = np.zeros_like(position)
    res[position < start] = 1.0
    res[position > stop ] = 0.0
    theMask  = (position >= start) & (position <= stop)
    position = position[theMask]

    #0th to 2nd order
    if order <= 2:
        a      = a_map[order];
        b      = b_map[order];
        x      = position;
        DELTAx = stop - start;
        x_t    = start + 0.5*DELTAx;
        phi    = (x - x_t) / DELTAx * pi;
    else: #higher order
        #We need to do some arbitrary scaling:
        END   =  4.0
        START = -2.0
        factor= (END-START) / (stop-start)
        X     = START + (position - start) * factor


    resTMP = np.zeros_like(position)

    if   (order == 0):
        for i in range(len(position)):
            resTMP[i] = a                                * np.sin(phi[i])                           + b
    elif (order == 1):
        for i in range(len(position)):
            resTMP[i] = a * ( 1./2. * np.cos(phi[i])     * np.sin(phi[i]) + 1./2.*phi[i])           + b
    elif (order == 2):
        for i in range(len(position)):
            resTMP[i] = a * ( 1./3. * np.cos(phi[i])**2. * np.sin(phi[i]) + 2./3. * np.sin(phi[i])) + b
    else:
        for i in range(len(position)):
            resTMP[i] = 1. - (A + (K-A)  / np.power( 1. + Q * np.power(e,(-B*(X[i] - M))),1./nu))

    res[theMask] = resTMP
    
    if swapper is None: return 1-res
    else: return res


def transition_factor_alt(center = 0.5, length = 1.0, position = 0.0, order = 2):
    """Please see :py:func:`.transition_factor` for documentation"""    
    return transition_factor(start=center-0.5*length, stop=center+0.5*length, position=position, order=order)
