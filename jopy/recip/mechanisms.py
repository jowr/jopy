# -*- coding: utf-8 -*-
from __future__ import print_function, division

from ..base import JopyBaseClass

import numpy as np
from numpy import pi
from texttable import Texttable
from scipy.optimize import minimize, minimize_scalar
from abc import ABCMeta, abstractmethod

class RecipBase(JopyBaseClass):
    """
    The basic object that handles many tasks for the reciprocating 
    machine. It exploits as many interelations as possible and the 
    object can be optimised heavily by using explicit equations
    instead of the solver calls.
    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self):
        super(RecipBase, self).__init__()
        # direct geometry 
        self._cr = None # crankshaft radius
        self._cl = None # conrod length
        self._bo = None # bore
        self._pp = None # piston pin offset
        # Angle values used for offset calculations
        self._theta_0_TDC = None
        self._theta_0_BDC = None
        # Determines whether to use offset or not
        self._of = None
        # Other variables used to save calculation time
        self._l_min    = None # clearance from dead volume
        self._l_max    = None 
        self._l_cr_max = None
        self._l_cr_min = None
        
        
        
    @property
    def cr(self):
        """Crankshaft radius [m]"""
        return self._cr

    #@cr.setter
    #def cr(self, value):
    #    self._cr = value
    #
    #@cr.deleter
    #def cr(self):
    #    del self._cr


    @property
    def cl(self):
        """Conrod length [m]"""
        return self._cl

    #@cl.setter
    #def cl(self, value):
    #    self._cl = value
    #
    #@cl.deleter
    #def cl(self):
    #    del self._cl


    @property
    def bo(self):
        """Cylinder bore (diameter) [m]"""
        return self._bo

    #@bo.setter
    #def bo(self, value):
    #    self._bo = value
    #
    #@bo.deleter
    #def bo(self):
    #    del self._bo


    @property
    def pp(self):
        """Piston pin offset from crankshaft centre [m]"""
        return self._pp

    #@pp.setter
    #def pp(self, value):
    #    self._pp = value
    #
    #@pp.deleter
    #def pp(self):
    #    del self._pp
    
    @property
    def l_min(self):
        """Clearance height at TDC/minimum volume [m]"""
        if self._l_min is None: self._l_min = self._calc_distance_to_head(self.theta_0_TDC)
        return self._l_min
    
    @property
    def l_max(self):
        """Clearance height at BDC/maximum volume [m]"""
        if self._l_max is None: self._l_max = self._calc_distance_to_head(self.theta_0_BDC)
        return self._l_max
    
    @property
    def l_cr_min(self):
        """Clearance height at TDC/minimum volume [m]"""
        if self._l_cr_min is None: self._l_cr_min = self._calc_distance_to_shaft(self.theta_0_BDC)
        return self._l_cr_min
    
    @property
    def l_cr_max(self):
        """Clearance height at TDC/minimum volume [m]"""
        if self._l_cr_max is None: self._l_cr_max = self._calc_distance_to_shaft(self.theta_0_TDC)
        return self._l_cr_max


    @property
    def theta_0_TDC(self):
        """Crankshaft angle at TDC without offset [rad] - 0 = crankshaft top, not TDC"""
        if self._theta_0_TDC is None: self._theta_0_TDC = self._calc_theta_0_TDC()
        return self._theta_0_TDC

    #@theta_0_TDC.setter
    #def theta_0_TDC(self, value):
    #    self._theta_0_TDC = value
    #
    #@theta_0_TDC.deleter
    #def theta_0_TDC(self):
    #    del self._theta_0_TDC


    @property
    def theta_0_BDC(self):
        """Crankshaft angle at BDC without offset [rad] - pi = crankshaft bottom, not BDC"""
        if self._theta_0_BDC is None: self._theta_0_BDC = self._calc_theta_0_BDC()
        return self._theta_0_BDC

    #@theta_0_BDC.setter
    #def theta_0_BDC(self, value):
    #    self._theta_0_BDC = value
    #
    #@theta_0_BDC.deleter
    #def theta_0_BDC(self):
    #    del self._theta_0_BDC
    
    

    
    def set_geometry(self,cr,cl,bo,pp=0,cv=0,of=False):
        """ 
        Update the geometric variables and perform some calculations.
        
        Parameters
        ----------
        cr : float
            crankshaft radius [m]
        cl : float
            conrod length [m]
        bo : float
            bore [m]
        pp : float
            piston pin offset [m]
        cv : float
            clearance volume at TDC in [m³]
        of : float
            offset for the crank angle position [boolean]
            Allow for theta_TDC other than 0. If false (theta_TDC = 0)
            is assured by precalculating the angular offset.
            
        """
        self._cr = cr
        self._cl = cl
        self._bo = bo
        self._pp = pp
        # Angle values used for offset calculations
        self._theta_0_TDC = None
        self._theta_0_BDC = None
        # Determines whether to use offset or not
        self._of = of
        # Other variables used to save calculation time
        self._l_min    = cv / self.A() # distance from piston to head at TDC
        self._l_max    = None 
        self._l_cr_max = None
        self._l_cr_min = None
        
    def _calc_theta_0(self,theta):
        """Process the crankshaft angle input to yield a value without offset."""
        if not self._of: return theta + self.theta_0_TDC
        else: return theta
        
    def _calc_theta(self,theta_0):
        """Process the crankshaft angle input to yield a value with offset."""
        if not self._of: return theta_0 - self.theta_0_TDC
        else: return theta_0
    
    def _head_to_crank(self,position):
        return self.l_cr_max - position + self.l_min
        
    def _crank_to_head(self,position):
        return self.l_cr_max - position + self.l_min        
    
    @staticmethod
    def _calc_dthetadt(rpm):
        """Return the radians per time"""
        return 2.0*pi*rpm/60.0
    
    def _calc_theta_bounds(self,TDCtoBDC=None,BDCtoTDC=None,dldtheta=None,dldt=None):
        """Determine whether we travel from TDC to BDC or the other way."""
        small = 1e-12
        downstroke = None
        if TDCtoBDC is not None:
            downstroke = np.asanyarray(TDCtoBDC)*1. 
        elif BDCtoTDC is not None:
            downstroke = np.asanyarray(BDCtoTDC)*-1.
        elif dldtheta is not None:
            downstroke = np.asanyarray(dldtheta)*-1.
        elif dldt is not None:
            downstroke = np.asanyarray(dldt)*-1.
        else:
            self.autolog(str(TDCtoBDC))
            self.autolog(str(BDCtoTDC))
            self.autolog(str(dldtheta))
            self.autolog(str(dldt))
            raise ValueError("Unable to detect downstroke or upstroke.")
        
        goingdo = downstroke >  small
        goingup = downstroke < -small
        undefin = np.logical_not(np.logical_or(goingdo,goingup))
        if np.sum(undefin)>0: 
            self.autolog(str(TDCtoBDC))
            self.autolog(str(BDCtoTDC))
            self.autolog(str(dldtheta))
            self.autolog(str(dldt))
            raise ValueError("Unable to detect downstroke or upstroke, the derivative is close to zero.")
                
        ret = np.empty(downstroke.shape, dtype=(float,2))
        ret[goingdo] = (        self.TDC(), self.BDC())
        ret[goingup] = (-2.0*pi+self.BDC(), self.TDC())
        return ret,goingdo,goingup

    
    def stroke(self):
        """Stroke length of the current geometry"""
        return self.l(self.BDC())-self.l(self.TDC())
       
    def TDC(self):
        """Crankshaft angle at TDC"""
        return self._calc_theta(self.theta_0_TDC)
    
    def BDC(self):
        """Crankshaft angle at BDC"""
        return self._calc_theta(self.theta_0_BDC)
    
    # Mandatory functions needed to fill the geometry variables,
    # the methods are abstract even though they have an implementation
    # to force the user to actively reimplement the calls. This is
    # done to draw attention to the fact the solvers are highly 
    # inefficient. 
    @abstractmethod
    def _calc_theta_0_TDC(self):
        """Calculate the crankshaft angle at TDC""" 
        def f(x): return self._calc_distance_to_head(x)
        res = minimize_scalar(f, bounds=(-0.5*pi, 0.5*pi), method='bounded')
        self.autolog(str(res))
        return res.x
    
    @abstractmethod
    def _calc_theta_0_BDC(self): 
        """Calculate the crankshaft angle at BDC"""
        def f(x): return -self._calc_distance_to_head(x)
        res = minimize_scalar(f, bounds=(0.5*pi, 1.5*pi), method='bounded')
        self.autolog(str(res))
        return res.x
    
    # Some of the functions are written in a way that causes 
    # indefinite recursion - two functions calling each other 
    # require a redefinition of at least one of them.      
    @abstractmethod
    def _calc_distance_to_head(self,theta_0): 
        """Calculate the distance from cylinder head to piston pin"""
        return self._crank_to_head(self._calc_distance_to_shaft(theta_0))
    
    @abstractmethod
    def _calc_distance_to_shaft(self,theta_0): 
        """Calculate the distance from crankshaft centre to piston pin"""
        return self._head_to_crank(self._calc_distance_to_head(theta_0))
    
    @abstractmethod
    def _calc_theta_0_from_head(self,distance,TDCtoBDC=None,BDCtoTDC=None,dldtheta=None,dldt=None): 
        """Calculate the crankshaft angle from the distance between piston pin and head"""
        return self._calc_theta_0_from_crank(self._head_to_crank(distance),TDCtoBDC=TDCtoBDC,BDCtoTDC=BDCtoTDC,dldtheta=dldtheta,dldt=dldt)
    
    @abstractmethod
    def _calc_theta_0_from_crank(self,distance,TDCtoBDC=None,BDCtoTDC=None,dldtheta=None,dldt=None): 
        """Calculate the crankshaft angle from the distance between piston pin and crankshaft centre"""
        return self._calc_theta_0_from_head(self._crank_to_head(distance),TDCtoBDC=TDCtoBDC,BDCtoTDC=BDCtoTDC,dldtheta=dldtheta,dldt=dldt)
    
    # Some shorthand notations for further calculations
    def A(self):
        """Piston surface facing the control volume [m²]"""
        return pi * np.power(self.bo,2.0) / 4.0 
    
    def V(self,theta):
        """Volume in cylinder"""
        return self._calc_distance_to_head(self._calc_theta_0(theta)) * self.A()

    def l(self,theta):
        """Piston position from TDC + clearance height"""
        return self._calc_distance_to_head(self._calc_theta_0(theta))
    
    def theta(self,position,TDCtoBDC=None,BDCtoTDC=None,dldtheta=None,dldt=None): 
        return self._calc_theta_0_from_head
    
    def dldtheta(self,theta): 
        """dl/dtheta - derivative of piston position w.r.t. crankshaft angle"""
        raise NotImplementedError("Missing function")
    
    def dldt(self,theta,rpm): 
        """dl/dt - derivative of piston position w.r.t. time a.k.a. the piston velocity"""
        return self.dldtheta(theta)*self._calc_dthetadt(rpm)
    
    def d2ldtheta2(self,theta): 
        """d2l/dtheta2 - 2nd derivative of piston position w.r.t. crankshaft angle"""
        raise NotImplementedError("Missing function")
    
    def d2ldt2(self,theta,rpm): 
        """dl/dt - 2nd derivative of piston position w.r.t. time a.k.a. the piston acceleration"""
        return self.d2ldtheta2(theta)*np.power(self._calc_dthetadt(rpm),2.0)
    
    def dVdtheta(self,theta):
        """dV/dtheta - derivative of volume w.r.t. crankshaft angle"""
        return self.dldtheta(theta)*self.A()
    
    def d2Vdtheta2(self,theta):
        """d2V/dtheta2 - 2nd derivative of volume w.r.t. crankshaft angle"""
        return self.d2ldtheta2(theta)*self.A()
    

    def info(self):        
        table = Texttable()
        table.set_deco(Texttable.HEADER)
        # table.set_chars(['-', '|', '+', '='])
        table.set_cols_dtype(['t',  # text
                      'f',  # float (decimal) 
                      't',  # text  
                      't']) # text 
        table.set_cols_align(["l", "r", "l", "l"])
        table.set_cols_valign(["m", "m", "m", "m"])
        # table.set_cols_width([10, 12, 13, 13, 13])
        table.header(["Variable","Value","Unit","Description"])
        table.add_row(["cr", self._cr,"m","crankshaft radius"])
        table.add_row(["cl", self._cl,"m","conrod length"])
        table.add_row(["bo", self._bo,"m","bore"])
        table.add_row(["pp", self._pp,"m","piston pin offset"])
        table.add_row(["cv", self.V(self.TDC())*1e6,"cm3","clearance volume at TDC"])
        table.add_row(["TDC", np.degrees(self.TDC()),"deg","angle of piston TDC"])
        table.add_row(["BDC", np.degrees(self.BDC()),"deg","angle of piston BDC"])
        table.add_row(["Max V", self.V(self.BDC())*1e6,"cm3","volume at BDC"])
        table.add_row(["Min V", self.V(self.TDC())*1e6,"cm3","volume at TDC"])
        print("\n"+table.draw()+"\n")
        
    
    def revolution(self,num):
        # num: steps you get
        full = np.linspace(-pi,pi,num+1)
        return full[:-1]
    
    def plotVolume(self):
        full = self.revolution(1000)
        volu = self.volume(full)
        
        import matplotlib as mpl
        #mpl.use('Qt4Agg')
        from matplotlib.pyplot import plot, show
        import matplotlib.pylab as plt
        
        if self.DEBUG:
            fig, axs = plt.subplots(2, 1, sharex=True)
            rev = self.revolution(2000)
            pos  = self._position(rev)
            axs[0].plot(rev*180/pi,pos*100)
            axs[0].plot(self.TDC()*180/pi,self._position(self.TDC())*100,'o')
            iMin  = np.where(pos==pos.min())
            axs[0].plot(rev[iMin]*180/pi,self._position(rev[iMin])*100,'o')
            iMax  = np.where(pos==pos.max())
            axs[0].plot(rev[iMax]*180/pi,self._position(rev[iMax])*100,'o')
            axs[0].set_ylabel(r'Piston position (cm)')
            ax = axs[1]
            self.autolog("Position: ", str(pos.min()), str(pos.max()), str(self.stroke()))
            self.autolog("Volume:    ",str(volu.min()), str(volu.max()))
        else:
            fig, ax = plt.subplots(1, 1)
            
        ax.plot(full*180/pi,volu*1e6)
        ax.plot(self.TDC()*180/pi,self.volume(self.TDC())*1e6,'o')
        iMin  = np.where(volu==volu.min())
        ax.plot(full[iMin]*180/pi,self.volume(full[iMin])*1e6,'o')
        iMax  = np.where(volu==volu.max())
        ax.plot(full[iMax]*180/pi,self.volume(full[iMax])*1e6,'o')
        ax.set_xlabel(r'Crankshaft angle $\theta$ (deg)')
        ax.set_ylabel(r'Cylinder volume $V$ (cm$^3$)')
        show()
        
        

class RecipExplicit(RecipBase):
    """Class that implements the formulations from Bjarne Dindler Rasmussen's PhD 
    thesis and adds an explicit equation for the position to crankshaft angle 
    conversion."""
    
    def __init__(self):
        super(RecipExplicit, self).__init__()
        # Dimensionless properties
        self._sigma = None
        self._lambda = None
        
    def set_geometry(self,cr,cl,bo,pp=0,cv=0,of=False):
        """ 
        Update the geometric variables and perform some calculations.
        
        cr : float
            crankshaft radius [m]
        cl : float
            conrod length [m]
        bo : float
            bore [m]
        pp : float
            piston pin offset [m]
        cv : float
            clearance volume at TDC in [m3]
        of : float
            offset for the crank angle position [boolean]
            Allow for theta_TDC other than 0. If false (theta_TDC = 0)
            is assured by precalculating the angular offset.
        """
        super(RecipExplicit, self).set_geometry(cr,cl,bo,pp,cv,of)
        # Dimensionless properties
        self._sigma = pp / cl
        self._lambda = cr / cl
        
        #Carry out some basic comparisons
        _l_cr_max = np.sqrt(np.power(self.cl+self.cr,2.0)-np.power(self.pp,2.0))
        _l_cr_min = np.sqrt(np.power(self.cl-self.cr,2.0)-np.power(self.pp,2.0))
        
        #if self.l_cr_max <> _l_cr_max:
        #    raise ValueError("l_cr_max is not equal: {1:7.5f} vs. {1:7.5f}".format(self.l_cr_max,_l_cr_max))
        #if self.l_cr_min <> _l_cr_min:
        #    raise ValueError("l_cr_min is not equal: {1:7.5f} vs. {1:7.5f}".format(self.l_cr_min,_l_cr_min))       
        
        
    # Mandatory functions needed to fill the geometry variables
    def _calc_theta_0_TDC(self):
        """Calculate the crankshaft angle at TDC""" 
        return -np.arcsin(self.pp/(self.cl+self.cr))

    def _calc_theta_0_BDC(self): 
        """Calculate the crankshaft angle at BDC"""
        return -np.arcsin(self.pp/(self.cl-self.cr)) + pi
    
    # Some of the functions are written in a way that causes 
    # indefinite recursion - two functions calling each other 
    # require a redefinition of at least one of them.
    def _calc_distance_to_head(self,theta_0): 
        """Calculate the distance from cylinder head to piston pin"""
        return super(RecipExplicit, self)._calc_distance_to_head(theta_0)
    
    def _calc_distance_to_shaft(self,theta_0): 
        """Calculate the distance from crankshaft centre to piston pin"""
        def p(x): return np.power(x,2)
        return self.cl*np.sqrt(1.0 - p(self.cr*np.sin(theta_0) + self.pp)/p(self.cl)) + self.cr*np.cos(theta_0)
    
    def _calc_theta_0_from_head(self,distance,TDCtoBDC=None,BDCtoTDC=None,dldtheta=None,dldt=None): 
        """Calculate the crankshaft angle from the distance between piston pin and head"""
        return super(RecipExplicit, self)._calc_theta_0_from_head(distance,TDCtoBDC=TDCtoBDC,BDCtoTDC=BDCtoTDC,dldtheta=dldtheta,dldt=dldt)
    
    def _calc_theta_0_from_crank(self,distance,TDCtoBDC=None,BDCtoTDC=None,dldtheta=None,dldt=None): 
        """Calculate the crankshaft angle from the distance between piston pin and crankshaft centre"""
        def p(x): return np.power(x,2)
        bounds,do,up = self._calc_theta_bounds(TDCtoBDC, BDCtoTDC, dldtheta, dldt)
        x = np.empty_like(distance)
        x[do] = - np.sqrt(self._calc_theta_of_pos_root(distance))
        x[up] = + np.sqrt(self._calc_theta_of_pos_root(distance))    
        return -2.0*np.arctan((2.0*self.cr*self.pp + x)/(-p(self.cl) + p(self.cr) + 2.0*self.cr*distance + p(distance) + p(self.pp)))
    
    def _calc_theta_of_pos_root(self, position):
        def p(x): return np.power(x,2)
        return -1.0*(p(self.cl) - 2.0*self.cl*self.cr + p(self.cr) - p(position) - p(self.pp))*(p(self.cl) + 2.0*self.cl*self.cr + p(self.cr) - p(position) - p(self.pp))

    
    def dldtheta(self,theta): 
        """dl/dtheta - derivative of piston position w.r.t. crankshaft angle"""
        def p(x): return np.power(x,2)
        theta_0 = self._calc_theta_0(theta)
        dcrankdtheta = -self.cr*np.sin(theta_0) - self.cr*(self.cr*np.sin(theta_0) + self.pp)*np.cos(theta_0)/(self.cl*np.sqrt(1.0 - p(self.cr*np.sin(theta_0) + self.pp)/p(self.cl)))
        return dcrankdtheta * -1.0 
    
    def d2ldtheta2(self,theta): 
        """d2l/dtheta2 - 2nd derivative of piston position w.r.t. crankshaft angle"""
        def p(x,y=2.0): return np.power(x,y)
        theta_0 = self._calc_theta_0(theta)
        d2crankdtheta2 = -self.cr*np.cos(theta_0) - p(self.cr)*p(np.cos(theta_0))/(self.cl*np.sqrt(1.0 - p(self.cr*np.sin(theta_0) + self.pp)/p(self.cl))) + self.cr*(self.cr*np.sin(theta_0) + self.pp)*np.sin(theta_0)/(self.cl*np.sqrt(1.0 - p(self.cr*np.sin(theta_0) + self.pp)/p(self.cl))) - p(self.cr)*p(self.cr*np.sin(theta_0) + self.pp)*p(np.cos(theta_0))/(p(self.cl,3)*p(1.0 - p(self.cr*np.sin(theta_0) + self.pp)/p(self.cl),3./2.))
        #d2crankdtheta2 = -cr*cos(theta) - cr**2*cos(theta)**2/(cl*sqrt(1.0 - (cr*sin(theta) + pp)**2/cl**2)) + cr*(cr*sin(theta) + pp)*sin(theta)/(cl*sqrt(1.0 - (cr*sin(theta) + pp)**2/cl**2)) - cr**2*(cr*sin(theta) + pp)**2*cos(theta)**2/(cl**3*(1.0 - (cr*sin(theta) + pp)**2/cl**2)**(3/2))
        return d2crankdtheta2 * -1.0 

         

    
    
class RecipImplicit(RecipBase):
    """Alternative definition of the reciprocating machine that uses solvers to determine 
    minimun and maximun values. Mostly used for testing of the analytical solutions. 
    Definition of piston movement and formulae taken from Dubbel, pages P5 to P7"""
    
    def set_geometry(self,cr,cl,bo,pp=0,cv=0,of=False):
        """ 
        Update the geometric variables and perform some calculations.
        
        cr : float
            crankshaft radius [m]
        cl : float
            conrod length [m]
        bo : float
            bore [m]
        pp : float
            piston pin offset [m]
        cv : float
            clearance volume at TDC in [m3]
        of : float
            offset for the crank angle position [boolean]
            Allow for theta_TDC other than 0. If false (theta_TDC = 0)
            is assured by precalculating the angular offset.
        """
        super(RecipImplicit, self).set_geometry(cr,cl,bo,pp,cv,of)
        
        #Carry out some basic comparisons
        _l_max = np.sin(   self.theta_0_TDC+0.5*pi) * (self.cl + self.cr)
        _l_min = np.cos(pi-self.theta_0_BDC       ) * (self.cl - self.cr)
        
        #if self.l_max <> _l_max:
        #    raise ValueError("l_cr_max is not equal: {1:7.5f} vs. {1:7.5f}".format(self.l_max,_l_max))
        #if self.l_min <> _l_min:
        #    raise ValueError("l_cr_min is not equal: {1:7.5f} vs. {1:7.5f}".format(self.l_min,_l_min))   
    
    # Custom helper functions
    def _beta(self,_theta_0):
        """Conrod angle"""
        return (self._pp+self._cr*np.sin(_theta_0))/self._cl
    
    # Mandatory functions needed to fill the geometry variables
    def _calc_theta_0_TDC(self):
        """Calculate the crankshaft angle at TDC""" 
        return super(RecipImplicit, self)._calc_theta_0_TDC()

    def _calc_theta_0_BDC(self): 
        """Calculate the crankshaft angle at BDC"""
        return super(RecipImplicit, self)._calc_theta_0_BDC()
    
    # Some of the functions are written in a way that causes 
    # indefinite recursion - two functions calling each other 
    # require a redefinition of at least one of them.
    def _calc_distance_to_head(self,theta_0): 
        """Calculate the distance from cylinder head to piston pin"""
        return np.sqrt(np.power(self.cl+self.cr,2.0)-np.power(self.pp,2.0))-self.cl*np.cos(self._beta(-theta_0 + pi))+self.cr*np.cos(-theta_0 + pi) + self.l_min
    
    def _calc_distance_to_shaft(self,theta_0): 
        """Calculate the distance from crankshaft centre to piston pin"""
        def p(x): return np.power(x,2)
        return self.cl*np.sqrt(1.0 - p(self.cr*np.sin(theta_0) + self.pp)/p(self.cl)) + self.cr*np.cos(theta_0)
    
    def _calc_theta_0_from_head(self,distance,TDCtoBDC=None,BDCtoTDC=None,dldtheta=None,dldt=None): 
        """Calculate the crankshaft angle from the distance between piston pin and head"""
        bounds,do,up = self._calc_theta_bounds(TDCtoBDC, BDCtoTDC, dldtheta, dldt)
        pos    = np.asarray(distance)
        guess  = np.ones_like(pos) * np.mean(bounds)
        def f(x): return np.sum(np.power(pos-self.l(x),2.0)) 
        res = minimize(f, guess, bounds=bounds, method='L-BFGS-B')
        self.autolog(str(res))
        return res.x
    
    def _calc_theta_0_from_crank(self,distance,TDCtoBDC=None,BDCtoTDC=None,dldtheta=None,dldt=None): 
        """Calculate the crankshaft angle from the distance between piston pin and crankshaft centre"""
        return super(RecipImplicit, self)._calc_theta_0_from_head(distance,TDCtoBDC=TDCtoBDC,BDCtoTDC=BDCtoTDC,dldtheta=dldtheta,dldt=dldt)
    
   

# class CylinderHetaTransfer(object):
#     
#     def basicRePrCorrelation(self,fluid,T_f,rho_f,char_vel,char_len,a,b,c,DEBUG=False):
#         
#         # sanitise two-phase situations        
#         Q = PropsSI('Q', 'T', T_f, 'D', rho_f, fluid)
#         
#         mu_f     = numpy.array([])
#         lambda_f = numpy.array([])
#         cp_f     = numpy.array([])
#         for i,qVal in enumerate(Q):
#             mu_f     = numpy.append(mu_f,    [PropsSI('V', 'T', T_f[i], 'D', rho_f[i], fluid)])
#             lambda_f = numpy.append(lambda_f,[PropsSI('L', 'T', T_f[i], 'D', rho_f[i], fluid)])
#             if 0>Q[i] or 1<Q[i]:        
#                 #mu_f     = numpy.append(mu_f,    [PropsSI('V', 'T', T_f[i], 'D', rho_f[i], fluid)])
#                 #lambda_f = numpy.append(lambda_f,[PropsSI('L', 'T', T_f[i], 'D', rho_f[i], fluid)])
#                 cp_f     = numpy.append(cp_f,    [PropsSI('C', 'T', T_f[i], 'D', rho_f[i], fluid)])
#             else:
#                 #mu_f     = numpy.append(mu_f,    [-1])
#                 #lambda_f = numpy.append(lambda_f,[-1])
#                 cp_f     = numpy.append(cp_f,    [1e5])
#             
#         
#         
#         if mu_f.any     <= 0: print "Invalid viscosity, make sure transport properties are calculated."
#         if lambda_f.any <= 0: print "Invalid thermal conductivity, make sure transport properties are calculated."
#         if cp_f.any     <= 0: print "Invalid heat capacity, make sure properties are calculated correctly."
#         Pr       = cp_f * mu_f / lambda_f 
#         Re       = (rho_f * char_vel * char_len) / mu_f
#         Nu       = a * numpy.power(Re,b) * numpy.power(Pr,c) 
#         h        = Nu * lambda_f / char_len
#         if DEBUG:
#             printDebug(lambda_f,namespace=locals())
#             printDebug(cp_f,namespace=locals())
#             printDebug(mu_f,namespace=locals())
#             printDebug(Pr,namespace=locals())
#             printDebug(Re,namespace=locals())
#             printDebug(Nu,namespace=locals())
#             printDebug(h,namespace=locals())
#         return h,Nu
#        
#     def Kornhauser1994(self,fluid,T_f,rho_f,char_vel,char_len,a=0.56,b=0.69,DEBUG=False):
#         lambda_f = CP.PropsU('L', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         cp_f     = CP.PropsU('C', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         alpha_f  = lambda_f / (rho_f * cp_f)
#         Pe       = (char_vel * char_len**2) / (4*alpha_f)
#         Nu       = a * Pe**b
#         h        = Nu * lambda_f / char_len
#         if DEBUG:
#             printDebug(lambda_f,namespace=locals())
#             printDebug(cp_f,namespace=locals())
#             printDebug(alpha_f,namespace=locals())
#             printDebug(Pe,namespace=locals())
#             printDebug(Nu,namespace=locals())
#             printDebug(h,namespace=locals())
#         return h,Nu
#     
#         
#     def Destoop1986(self,fluid,T_f,rho_f,Lambda,Gamma,DEBUG=False):
#         """
#         Destoop's simplified correlation for ammonia compressors
#         """
#         h,Nu = self.basicRePrCorrelation(fluid,T_f,rho_f,Lambda,Gamma,0.6,0.8,0.6,DEBUG=DEBUG)
#         return h,Nu 
#     
#     
#     def Annand1963(self,fluid,T_f,rho_f,Lambda,Gamma,DEBUG=False):
#         """
#         Annand's IC engine correlation
#         """
#         h,Nu = self.basicRePrCorrelation(fluid,T_f,rho_f,Lambda,Gamma,0.575,0.7,0.0,DEBUG=DEBUG)
#         return h,Nu 
#     
#     
#     def Annand1963b(self,fluid,T_f,rho_f,char_vel,char_len,T_w,DEBUG=False):
#         """
#         Annand's SI engine correlation
#         """
#         c1 = 0.575
#         #c1 = 0.35 # for comparison with Irimescu
#         c2 = 4.3e-9 
#         
#         mu_f     = CP.PropsU('V', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         lambda_f = CP.PropsU('L', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         #cp_f     = CP.PropsU('C', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         if mu_f.any     <= 0: print "Invalid viscosity, make sure transport properties are calculated."
#         if lambda_f.any <= 0: print "Invalid thermal conductivity, make sure transport properties are calculated."
#         #if cp_f.any     <= 0: print "Invalid heat capacity, make sure properties are calculated correctly."
#         #Pr       = cp_f * mu_f / lambda_f 
#         Re       = (rho_f * char_vel * char_len) / mu_f
#         
#         h  = c1 * lambda_f / char_len * numpy.power(Re,0.7) + c2 * (numpy.power(T_f,4)-numpy.power(T_w,4))/(T_f-T_w)
#         Nu = h / lambda_f * char_len
# 
#         return h,Nu 
#     
#         
#     def Woschni1967(self,fluid,T_f,rho_f,Lambda,Gamma,DEBUG=False):
#         """
#         Annand's IC engine correlation
#         """
#         h,Nu = self.basicRePrCorrelation(fluid,T_f,rho_f,Lambda,Gamma,0.035,0.7,0.0,DEBUG=DEBUG)
#         return h,Nu 
#     
#     
#     def Adair1972(self,theta,omega,T_f,rho_f,fluid,bore,position,DEBUG=False):
#         """
#         Adair's heat transfer correlation
#         theta: crankshaft angle in rad with TDC at 0.
#         omega: angular velocity RPS (revolutions per second), get converted into rad/s
#         T_f: bulk fluid temperature
#         rho_f: bulk fluid density
#         fluid: fluid string
#         bore: cylinder bore
#         position: distance piston-head
#         """
# 
#         theta = numpy.mod(theta-math.pi,2*math.pi) # Fix TDC-BDC problem 
#         #omega = der(crankshaftAngle) 
#         # Equation 15 from Adair et al.
#         omega    = 2.*math.pi*omega 
#         # TODO: Check conversion from 1/s to rad/s
#         omega_g1 = 2.*omega *       (1.04+numpy.cos(2.*theta))
#         omega_g2 = 2.*omega * 0.5 * (1.04+numpy.cos(2.*theta))
#           
#         omega_g = 0. * omega_g1
#         
#         for i in range(len(theta)):
#             if theta[i]>0.5*math.pi and theta[i]<1.5*math.pi:
#                 omega_g[i] = omega_g2[i]
#             else:
#                 omega_g[i] = omega_g1[i]
#            
#         #surfaceArea = pistonCrossArea + 2. * math.sqrt(pistonCrossArea*math.pi)*position
#         volume      = math.pi * numpy.power(bore,2.) / 4. * position #Get volumes"
#         surfaceArea = math.pi * numpy.power(bore,2.) / 4. * 2. + math.pi * bore * position 
#         
#         d_e = 6. / surfaceArea * volume
#         
#         Gamma  = d_e
#         Lambda = 0.5 * d_e * omega_g 
#         
#         h,Nu = self.basicRePrCorrelation(fluid,T_f,rho_f,Lambda,Gamma,0.053,0.8,0.6,DEBUG=DEBUG)
# 
#         return h,Nu
#     #    # There is a small mistake in equation 19 of the paper, DeltaT goes in the numerator.
#     #    Ts = PropsSI('T', 'H', h_f, 'P', p_f, fluid) 
#     #    q_w = -1. * h * (Ts - T_wall)
#     #    Q_flows[i] = surfaceAreas[i]*q_w[i]
#     
#     
#     def BasicGnielinski(self,fluid,T_f,rho_f,char_vel,char_len,L,zeta,xtra,K,DEBUG=False):
#         
#         D = char_len 
#         
#         mu_f        = CP.PropsU('V', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         lambda_f    = CP.PropsU('L', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         cp_f        = CP.PropsU('C', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         if mu_f.any     <= 0: print "Invalid viscosity, make sure transport properties are calculated."
#         if lambda_f.any <= 0: print "Invalid thermal conductivity, make sure transport properties are calculated."
#         if cp_f.any     <= 0: print "Invalid heat capacity, make sure properties are calculated correctly."
#         Pr          = cp_f * mu_f / lambda_f 
#         Re          = (rho_f * char_vel * char_len) / mu_f        
#         numerator   = (zeta/8.) * (Re-xtra) * Pr 
#         denominator = 1 + 12.7 * numpy.sqrt(zeta/8.) * (numpy.power(Pr,2./3.)-1.)
#         Nu          = numerator / denominator * (1 + numpy.power(D/L,2./3.)) * K  
#         h           = Nu * lambda_f / char_len
#         if DEBUG:
#             printDebug(lambda_f,namespace=locals())
#             printDebug(cp_f,namespace=locals())
#             printDebug(mu_f,namespace=locals())
#             printDebug(Pr,namespace=locals())
#             printDebug(Re,namespace=locals())
#             printDebug(Nu,namespace=locals())
#             printDebug(h,namespace=locals())
#         return h,Nu
# 
# 
#     def Gnielinski1976(self,fluid,T_f,rho_f,char_vel,char_len,L,T_w,rho_w,DEBUG=False):        
#         mu_f        = CP.PropsU('V', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         mu_w        = CP.PropsU('V', 'T', T_w, 'D', rho_w, fluid, 'SI') 
#         if mu_f.any     <= 0: print "Invalid viscosity, make sure transport properties are calculated."
#         Re          = (rho_f * char_vel * char_len) / mu_f       
#         
#         zeta        = numpy.power((1.82 * numpy.log10(Re)-1.64),-2)
#         xtra        = 1000.
#         K           = numpy.power(mu_f/mu_w,0.14)
#         #for xi in numpy.array(numpy.power(mu_f/mu_w,0.14)):
#         #    print xi 
#         
#         h,Nu = self.BasicGnielinski(fluid,T_f,rho_f,char_vel,char_len,L,zeta,xtra,K,DEBUG=DEBUG)
#         return h,Nu
#         
#         
#     def Gnielinski2010(self,fluid,T_f,rho_f,char_vel,char_len,L,T_w,DEBUG=False):
#         mu_f        = CP.PropsU('V', 'T', T_f, 'D', rho_f, fluid, 'SI') 
#         if mu_f.any     <= 0: print "Invalid viscosity, make sure transport properties are calculated."
#         Re          = (rho_f * char_vel * char_len) / mu_f       
#         
#         zeta        = numpy.power((1.80 * numpy.log10(Re)-1.50),-2)
#         xtra        = 0.
#         K           = 1. #numpy.power(T_f/T_w,0.25)
#         
#         h,Nu = self.BasicGnielinski(fluid,T_f,rho_f,char_vel,char_len,L,zeta,xtra,K,DEBUG=DEBUG)
#         return h,Nu
# 
# class SteadyState(object):
#     
#     def idealNozzle(self,fluid,h_up,p_up,p_down,DEBUG=False):
#         """Model of a nozzle flow using perfect gas parameter gamma. 
#         
#         """        
#         from pyrp.refpropClasses import RefpropSI
#         RP = RefpropSI()
#         RP.SETUPFLEX(FluidNames=fluid)
#         T_up,p_up,rho_up,Dl_up,Dv_up,q_up,e_up,h_up,s_up,cv_up,cp_up,w_up = RP.PHFLSH(p_up,h_up)
#         v_up = 1/rho_up 
#         gamma0    = cp_up / cv_up
#         
#         def function(gamma):
#             p_thr_crit = p_up*numpy.power(2/(gamma+1),gamma/(gamma-1))
#             p_thr_a    = max(p_thr_crit,p_down)
#             v_thr_a    = ((p_up*v_up**gamma)/p_thr_a)**(1./gamma)
#             T_thr,p_thr,rho_thr,Dl_thr,Dv_thr,q_thr,e_thr,h_thr,s_thr,cv_thr,cp_thr,w_thr = RP.PSFLSH(p_thr_a,s_up)
#             return (v_thr_a-1/rho_thr)**2
#         
#         res = minimize(function, gamma0)
#         #if DEBUG:
#         #    print res, "\n"
#         gamma = res.x[0]
#         p_thr_crit = p_up*(2/(gamma+1))**(gamma/(gamma-1))
#         p_thr    = max(p_thr_crit,p_down)
#         T_thr,p_thr,rho_thr,Dl_thr,Dv_thr,q_thr,e_thr,h_thr,s_thr,cv_thr,cp_thr,w_thr = RP.PSFLSH(p_thr,s_up)
#         vel_thr = ((h_up - h_thr)*2.)**(1/2.)
#         if DEBUG: 
#             print p_thr, vel_thr, w_thr, "\n"  
#         return T_thr,p_thr,rho_thr,Dl_thr,Dv_thr,q_thr,e_thr,h_thr,s_thr,cv_thr,cp_thr,w_thr,vel_thr,gamma
#     
#     
#     def simplyfiedModel(self,fluid,r,l,q,b=0.1,c=20e-6):
#         mechanism = MechanicsAlt()
#         mechanism.setGeometry(r,l,q,b=b,c=c)
#         
#         fluid = "pentane"
#         
#         from pyrp.refpropClasses import RefpropSI
#         RP = RefpropSI()
#         RP.SETUPFLEX(FluidNames=fluid)
#         
#         d_thr_su   = 0.022 
#         d_leak = 0.00022
#         d_thr_ex   = 0.022
#         
#         N_exp = 1000.
#         N_exp_s = N_exp / 60.
#         
#         V_s = 750e-6
#         V_0 =  36e-6
#         
#         p_supply = 15e5
#         T_supply = 155+273.15
#         T_su,p_su,rho_su,Dl_su,Dv_su,q_su,e_su,h_su,s_su,cv_su,cp_su,w_su = RP.TPFLSH(T_supply,p_supply)
#         
#         p_exhaust = 1.5e5
#         T_exhaust = 90+273.15
#         T_ex,p_ex,rho_ex,Dl_ex,Dv_ex,q_ex,e_ex,h_ex,s_ex,cv_ex,cp_ex,w_ex = RP.TPFLSH(T_exhaust,p_exhaust)
#         
#         StSt = SteadyState()
#         
#         def cycle(r_p_su1,r_M_dot_leak):
#             p_su1 = r_p_su1 * p_su 
#             T_thr_su,p_thr_su,rho_thr_su,Dl_thr_su,Dv_thr_su,q_thr_su,e_thr_su,h_thr_su,s_thr_su,cv_thr_su,cp_thr_su,w_thr_su,vel_thr_su,gamma_thr_su = StSt.idealNozzle(fluid, h_su, p_su, p_su1,DEBUG=True)
#             A_thr_su     = (pi*d_thr_su^2)/4
#             V_dot_thr_su = vel_thr_su*A_thr_su
#             M_dot  = V_dot_thr_su * rho_thr_su
#             T_su1,p_su1,rho_su1,Dl_su1,Dv_su1,q_su1,e_su1,h_su1,s_su1,cv_su1,cp_su1,w_su1 = RP.PHFLSH(p_su1,h_su)
#             # No heat transfer
#             T_su2,p_su2,rho_su2,Dl_su2,Dv_su2,q_su2,e_su2,h_su2,s_su2,cv_su2,cp_su2,w_su2 = RP.PHFLSH(p_su1,h_su)
#             # Leakage
#             T_leak,p_leak,rho_leak,Dl_leak,Dv_leak,q_leak,e_leak,h_leak,s_leak,cv_leak,cp_leak,w_leak,vel_leak,gamma_leak = StSt.idealNozzle(fluid, h_su2, p_su2, p_ex,DEBUG=True)
#             A_leak     = (pi*d_leak^2)/4
#             V_dot_leak = vel_leak*A_leak
#             M_dot_leak  = V_dot_leak * rho_leak
#             
#             r_M_dot_leak = M_dot_leak / M_dot
#             M_dot_in   = M_dot - M_dot_leak
#             
#             W_dot_suc = N_exp_s*(V_s-V_0)*p_su2
#             



if __name__ == "__main__":
    raise ValueError("Do not call this package directly.")
#    StSt = SteadyState()
#    p = linspace(5e5,12.5e5,10)
#    for pi in p:
#        StSt.idealNozzle("pentane", 6e5, 15e5, pi,DEBUG=True) 
#    
#    import matplotlib as mpl
#    mpl.use('Qt4Agg')
#    
#    from matplotlib.pyplot import plot, show
#    
#    me = Mechanics()
#    r = 0.05
#    l = 0.11
#    q = 0.03
#    me.setGeometry(r,l,q,b=0.1,c=20e-6)
#    
#    full = me.revolution(1000)
#    posi = me.position(full)
#    
#    
#    print posi.max(), posi[0], posi.min()
#    print me.stroke()
#    print posi.max() - posi.min()
#    
#    plot(full,posi)
#    plot(full[0],posi[0],'o')
#    
#    print me.position(me.TDC())
#    print me.position(me.BDC())
#    
#    print me.theta(0.3*me.position(me.BDC()), afterTDC=True, result=True)
#    print me.theta(0.3*me.position(me.BDC()), afterTDC=True) * 180 / pi
#    print ""
#    
#    print me.theta(0.3*me.position(me.BDC()), afterTDC=False, result=True)
#    print me.theta(0.3*me.position(me.BDC()), afterTDC=False) * 180 / pi
#    print ""
#    
#    me.info()
#    
#    show()

#print ""
#print me.stroke(), me.position(0) - me.position(me.offsetBDC()-me.offsetTDC())


#class Mechanics2(object):
#    def __init__(self):
#        self.DEBUG = False
#        
#        self.stroke = 0.1
#        self.bore   = 0.1
#        self.V_min  = 20e-6
#        self.steps  = 400
#        self.theta  = self.getAngles(self.steps)
#        self.l_con  = 0.15
#        
#        self.A_pis  = pi * self.bore**2 / 4
#        self.x_min  = self.V_min / self.A_pis
#        self.V_max  = self.V_min + self.stroke * self.A_pis
#        self.x      = self.positionNoOffset(self.theta, None)
#        self.V      = self.x * self.A_pis
#        
#    def info(self):
#        print "stroke:         " + str(self.stroke)
#        print "bore:           " + str(self.bore)
#        print "piston position:" + str(self.x_min) + " to " + str(self.x_min+self.stroke)
#        print "maximum volume: " + str(self.V_max)
#        print "minimum volume: " + str(self.V_min)
#    
#    def validate(self):
#        if self.DEBUG: print self.V_min, self.V_max
#        if self.DEBUG: print self.x[0]*self.A_pis, self.x[round(self.steps/2)]*self.A_pis
#        if self.DEBUG: print "\nThis should be a sequence from 0 to 180:"
#        if self.DEBUG: print self.thetaNoOffset(numpy.linspace(me.x[0],me.x[round(me.steps/2)],10))*180/pi
#
#    def positionNoOffset(self,theta,pos):
#        # r: crankshaft radius
#        r = self.stroke / 2.
#        # l: conrod length
#        l = self.l_con
#        # theta: crankshaft angle from TDC
#        theta = numpy.array(theta)
#        position = numpy.array(r - r*cos(theta)+l-sqrt(l**2-(r*sin(theta))**2)) + self.x_min
#        if pos is None: 
#            return position 
#        else: 
#            residuals = position - pos
#            return residuals
#
#    def thetaNoOffset(self,pos):
#        pos = numpy.array(pos)
#        x0     = [(0.25*pi) for s in pos] # initial guess for cut-off
#        #bounds = (0.1*pi, 0.9*pi)
#        #res = minimize(self.positionNoOffset, x0, bounds=bounds, args=(pos))
#        res = fsolve(self.positionNoOffset, x0, args=(pos), xtol=1e-06)
#        #print "Found best x of ",res.x," with an error of ",residuals," after ",iterationCounter," iterations."
#        return res
#
#    def getAngles(self,num):
#        # num: steps you get
#        full = numpy.linspace(0,2*pi,num+1)
#        return full[:-1]
