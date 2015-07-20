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
    Delta_T1 : scalar or array_like
        Temperature difference on one side of the heat exchanger        
    Delta_T2 : scalar or array_like
        Temperature difference on the other side of the heat exchanger 
    deadBand : scalar or array_like
        Absolute temperature difference that is considered to be 0.0 K
        
    Returns:
    --------
    scalar or array_like
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
        Delta_T1 : scalar or array_like
            Temperature difference on one side of the heat exchanger        
        Delta_T2 : scalar or array_like
            Temperature difference on the other side of the heat exchanger 
        dead_band : scalar or array_like
            Absolute temperature difference that is considered to be 0.0 K
        
    Returns:
    --------
        scalar or array_like
            The calculated logarithmic mean temperature difference
    
    """
    Delta_T1       = np.asarray(Delta_T1)
    Delta_T2       = np.asarray(Delta_T2)    
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
    

def lmtd(Delta_T1, Delta_T2, dead_band=0.1):
    """Logarithmic mean temperature difference
    
    This function returns the logarithmic mean temperature 
    difference as calculated from temperature differences at
    both end. 
    
    Parameters:
    -----------
        Delta_T1 : scalar or array_like
            Temperature difference on one side of the heat exchanger        
        Delta_T2 : scalar or array_like
            Temperature difference on the other side of the heat exchanger 
        dead_band : scalar or array_like
            Absolute temperature difference that is considered to be the minimum
        
    Returns:
    --------
        array_like
            The calculated logarithmic mean temperature difference
            
    Please see :py:func:`._lmtd_jopy` for more documentation.
    
    """
    return _lmtd_jopy(Delta_T1, Delta_T2, dead_band)   
    #return np.vectorize(_lmtd_libpf)(Delta_T1, Delta_T2, dead_band)
    
    def SimpleRankineCycle(states=[],:
        """Logarithmic mean temperature difference
        
        This function returns the logarithmic mean temperature 
        difference as calculated from temperature differences at
        both end. 
        
        Parameters:
        -----------
            Delta_T1 : scalar or array_like
                Temperature difference on one side of the heat exchanger        
            Delta_T2 : scalar or array_like
                Temperature difference on the other side of the heat exchanger 
            dead_band : scalar or array_like
                Absolute temperature difference that is considered to be the minimum
            
        Returns:
        --------
            array_like
                The calculated logarithmic mean temperature difference
                
        Please see :py:func:`._lmtd_jopy` for more documentation.
        
        """

   
    def SimpleRankineCycle(T3, p3, T1, p1, epsilon_e, epsilon_p, fluid, axObj, plotterObj,label=False):

        state.update(CoolProp.PT_INPUTS,p1,T1)
        h1 = state.hmass()
        s1 = state.smass()

        p2 = p3
        state.update(CoolProp.PSmass_INPUTS,p2,s1)
        h2 = h1 + (state.hmass() - h1) / epsilon_p
        state.update(CoolProp.HmassP_INPUTS,h2,p2)
        s2 = state.smass()
        T2 = state.T()

        state.update(CoolProp.PT_INPUTS,p3,T3)
        h3 = state.hmass()
        s3 = state.smass()

        p4 = p1
        state.update(CoolProp.PSmass_INPUTS,p4,s3)
        h4 = h3 - epsilon_e * (h3 - state.hmass())
        state.update(CoolProp.HmassP_INPUTS,h4,p4)
        s4 = state.smass()
        T4 = state.T()

        w_net = h3 - h4
        q_boiler = h3 - h2
        eta_c = w_net / q_boiler

        #Ts = PropsPlot(fluid, 'Ts')
        #Ts.draw_isolines('P', [p1, p3], num=10)
        #Ts.set_axis_limits([0., 12., 200., 900.])

        #axObj.plot(s_tp/1e3,T_tp-273.15 , color=plotterObj._black, ls='-', alpha=1.0)

        isoObj  = IsoLines(fluid, "Ts", "Q")
        isoqual = isoObj.get_isolines([0.0,1.0], num=2)

        x = np.append(isoqual[ 0]['x'],isoqual[-1]['x'][::-1])/1e3
        y = np.append(isoqual[ 0]['y'],isoqual[-1]['y'][::-1])-273.15
        axObj.plot(x,y, color=plotterObj._black, ls='-', alpha=1.0)

        isoObj  = IsoLines(fluid, "Ts", "P")
        prange  = [p1,2e5,5e5,10e5,p3]
        isobars = isoObj.get_isolines(prange, num=len(prange))

        p = -1
        for c,i in enumerate(isobars):
            x = i['x']/1e3
            y = i['y']-273.15
            dp = prange[c]/1e5 - p
            p = prange[c]/1e5
            s = PropsSI('S','P',p*1e5,'Q',0.5,fluid)/1e3
            #print "Delta p: {0}".format(dp)
            if abs(dp)>0.8: #c%2==0 :
                axObj.plot(  x,    y, color=plotterObj._black, ls='-', alpha=0.50)
                if label:
                    putXLabel(xv=x, yv=y, x=s, text="{0:3.1f} bar".format(p), axis=axObj)

        #for i in range(len(x)):
        #    axObj.plot(  x[i]/1e3,   y[i]-273.15, color=plotterObj._black, ls='-', alpha=0.5)
        #    putXLabel(xv=x[i]/1e3,yv=y[i]-273.15, x=0, text="", axis=axObj)



        # Create the process lines
        A = []
        A.append({'H':h1,'P':p1,'S':s1,'T':T1})
        A.append({'H':h2,'P':p2,'S':s2,'T':T2})
        A.append({'H':h3,'P':p3,'S':s3,'T':T3})
        A.append({'H':h4,'P':p4,'S':s4,'T':T4})

        A.append(A[0].copy())

        processes = []

        for i in range(len(A)-1):
            s = np.linspace(         A[i]['S'] ,          A[i+1]['S'] ,num=points)
            p = np.logspace(np.log10(A[i]['P']), np.log10(A[i+1]['P']),num=points)
            dic = {}
            dic['P'] = p
            dic['S'] = s
            dic['T'] = PropsSI('T','P',p,'S',s,fluid)
            processes.append(dic)

        x = []
        y = []
        for lin in processes:
            #axObj.plot(lin['S']/1e3,lin['T']-273.15,color=plotterObj._black, linestyle='--')
            x.extend(lin['S']/1e3)
            y.extend(lin['T']-273.15)

        plotterObj.plotData([x],[y],ax=axObj,legend=False)

        x = np.array([s1,s2,s3,s4])
        y = np.array([T1,T2,T3,T4])

        #print x
        #print y
        #print " "

        plotterObj.plotData([x/1e3],[y-273.15],ax=axObj,legend=False)

        #axObj.plot(x/1e3,y-273.15,'o',color=plotterObj._black)

        #plotterObj.drawLegend(ax=axObj,loc=0) # the relative size of legend markers vs. original
        axObj.set_xlabel(ur"Specific entropy $s$ / \si{\kilo\joule\per\kilo\gram\per\kelvin}")
        axObj.set_ylabel(ur"Temperature $T$ / \si{\celsius}")
        axObj.set_xlim([-0.25,1.60])
        axObj.set_ylim([-25,325])

        #plotterObj.plotData([x], [y], ax=axObj)

        #ax = Ts.axis
        #ax.text(s1/1000., T1,' 1', fontsize=10, rotation=0, color='r')
        #ax.text(s2/1000., T2,' 2', fontsize=10, rotation=0, color='r')
        #ax.text(s3/1000., T3,' 3', fontsize=10, rotation=0, color='r')
        #ax.text(s4/1000., T4,' 4', fontsize=10, rotation=0, color='r')
        #ax.text(8., 850., "Efficiency: %.1f%%" %(eta_c*100.))
        #ax.text(8., 800., "Net work: %d kJ/kg" %(w_net/1000))
        #ax.text(8., 750., "Heat input: %d kJ/kg" %(q_boiler/1000))

    simPlotterObj = BasePlotter()
    figPV = simPlotterObj.getFigure(**sixupProps)
    simPlotterObj.ccycle = simPlotterObj.multiplyCycle(simPlotterObj.getColorCycle(length=3),doubles=2)
    simPlotterObj.scycle = cycle(['-'])
    simPlotterObj.mcycle = cycle(['None'])





if __name__ == "__main__":
    
    var = np.linspace(-20, 20)
    con = np.zeros_like(var)+10
    plt.plot(var,lmtd(var,con),'o')
    plt.show()
    
    #print(np.vectorize(_lmtd_libpf)([0,1,5,10,10,10,10],[10,10,10,10,5,1,0]))
    #print(np.vectorize(_lmtd_libpf)([-0,-1,-5,-10,-10,-10,-10],[10,10,10,10,5,1,0]))
    
    #res = _lmtd_libpf(10.0, 10.0)
    #print(res,10.0)
