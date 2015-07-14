

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    xdata = np.linspace(0,100)
    ydata = []
    for _ in range(4):
        ydata.append(np.random.normal(size=len(xdata)))

    plt.figure()
    for y in ydata: plt.plot(xdata,y)
    #plt.show()
    
    import matplotlib as mpl
    
    from jopy.styles.mplib import BaseStyle
    plt.figure()
    bs = BaseStyle()
    bs.update_rc_params()    
    for y in ydata: plt.plot(xdata,y)
    
    from jopy.styles.mplib import DtuStyle
    plt.figure()
    bs = DtuStyle()
    bs.update_rc_params()    
    for y in ydata: plt.plot(xdata,y)
    
    from jopy.styles.mplib import IpuStyle
    plt.figure()
    bs = IpuStyle()
    bs.update_rc_params()    
    for y in ydata: plt.plot(xdata,y)
    
    plt.show()
#     
#     
#     
#     
#     from jopy.recip.mechanisms import RecipExplicit, RecipImplicit
#     #from math import pi
#     import numpy as np
#     
#         
#     me = RecipImplicit()
#     metoo = RecipExplicit()
#     cr = 0.05
#     cl = 0.1
#     bo = 0.1
#     pp = 0.8*cr
#     cv = 20e-6
#     me.set_geometry(cr, cl, bo, pp, cv)
#     metoo.set_geometry(cr, cl, bo, pp, cv)
#     
#     full = me.revolution(1000)
#     
#     pos = me.calc_distance_to_head(full)
#     postoo = metoo.calc_distance_to_head(full)
#     #fultoo = metoo._calc_theta_from_distance(postoo)
#     
#     #plt.figure()
#     #plt.plot(full,pos)
#     #plt.plot(full,postoo)
#     #plot(fultoo,postoo,':')
#     
#     print("  TDC :  BDC ")
#     print("{0:8.2f} : {1:8.2f}".format(np.degrees(me._theta_TDC),    np.degrees(me._theta_BDC)))
#     print("{0:8.2f} : {1:8.2f}".format(np.degrees(metoo._theta_TDC), np.degrees(metoo._theta_BDC)))
#     print("{0:8.2f} : {1:8.2f}".format(np.degrees(me.TDC()),    np.degrees(me.BDC())))
#     print("{0:8.2f} : {1:8.2f}".format(np.degrees(metoo.TDC()), np.degrees(metoo.BDC())))
#     print("{0:8.2f} : {1:8.2f}".format(np.degrees(me._l_max),    np.degrees(me._l_min)))
#     print("{0:8.2f} : {1:8.2f}".format(np.degrees(metoo._l_max), np.degrees(metoo._l_min)))
#     print("{0:8.4f} : {1:8.4f}".format(np.min(pos),   np.max(pos)))
#     print("{0:8.4f} : {1:8.4f}".format(np.min(postoo),np.max(postoo)))
#     #print(pi)    
#     
#     plt.figure()
#     for p in np.linspace(0, 0.8, 5):
#         pp = p*cr
#         me.set_geometry(cr, cl, bo, pp, cv)
#         metoo.set_geometry(cr, cl, bo, pp, cv)
#         diff = (me.volume(full)-metoo.volume(full))/metoo.volume(full)*100.0
#         plt.plot(full,diff,label=str(p)+": (Dubbel-Bjarne)/Bjarne")
#         
#     plt.legend(loc=3)
#     
#     plt.figure()
#     plt.plot(full,me.volume(full)*1e6,label='Dubbel')
#     plt.plot(full,metoo.volume(full)*1e6,label='Bjarne')
#     plt.legend(loc=3)
#     
#     plt.show()
#     
#     me.info()
#     metoo.info()
    