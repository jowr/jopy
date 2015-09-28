# -*- coding: utf-8 -*-
from __future__ import print_function, division

if __name__ == "__main__":
    from jopy.recip.mechanisms import RecipExplicit, RecipImplicit
    #from math import pi
    import numpy as np
    import matplotlib.pyplot as plt
        
    me = RecipImplicit()
    metoo = RecipExplicit()
    cr = 0.05
    cl = 0.1
    bo = 0.1
    pp = 0.75*cr
    cv = 20e-6
    me.set_geometry(cr, cl, bo, pp, cv)
    metoo.set_geometry(cr, cl, bo, pp, cv)
    
    full = me.revolution(1000)
    
    pos = me.l(full)
    postoo = metoo.l(full)
    #fultoo = metoo._calc_theta_from_distance(postoo)
    
    #plt.figure()
    #plt.plot(full,pos)
    #plt.plot(full,postoo)
    #plot(fultoo,postoo,':')
    
    print("  TDC :  BDC ")
    print("{0:8.2f} : {1:8.2f}".format(np.degrees(me.theta_0_TDC),    np.degrees(me.theta_0_BDC)))
    print("{0:8.2f} : {1:8.2f}".format(np.degrees(metoo.theta_0_TDC), np.degrees(metoo.theta_0_BDC)))
    print("{0:8.2f} : {1:8.2f}".format(np.degrees(me.TDC()),    np.degrees(me.BDC())))
    print("{0:8.2f} : {1:8.2f}".format(np.degrees(metoo.TDC()), np.degrees(metoo.BDC())))
    print("{0:8.4f} : {1:8.4f}".format(me.l_cr_max,    me.l_cr_min))
    print("{0:8.4f} : {1:8.4f}".format(metoo.l_cr_max, metoo.l_cr_min))
    print("{0:8.4f} : {1:8.4f}".format(np.min(pos),   np.max(pos)))
    print("{0:8.4f} : {1:8.4f}".format(np.min(postoo),np.max(postoo)))
    #print(pi)    
    
    plt.figure()
    for p in np.linspace(0, 0.8, 5):
        pp = p*cr
        me.set_geometry(cr, cl, bo, pp, cv)
        metoo.set_geometry(cr, cl, bo, pp, cv)
        diff = (me.V(full)-metoo.V(full))/metoo.V(full)*100.0
        plt.plot(full,diff,label=str(p)+": (Dubbel-Bjarne)/Bjarne")
    plt.legend(loc=3)
    
    plt.figure()
    plt.plot(full,me.V(full)*1e6,label='Dubbel')
    plt.plot(full,metoo.V(full)*1e6,label='Bjarne')
    plt.legend(loc=3)
    
    plt.figure()
    #plt.plot(full,me.V(full)*1e6,label='Dubbel')
    plt.plot(full,metoo.dVdtheta(full)*1e6,label='Bjarne')
    plt.legend(loc=3)
    
    plt.figure()
    #plt.plot(full,me.V(full)*1e6,label='Dubbel')
    plt.plot(full,metoo.d2Vdtheta2(full)*1e6,label='Bjarne')
    plt.legend(loc=3)
    
    plt.show()
    
    me.info()
    metoo.info()
    