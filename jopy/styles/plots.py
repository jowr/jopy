# -*- coding: utf-8 -*-
from __future__ import print_function, division
import matplotlib
import matplotlib.pyplot as plt 
import copy

class Figure(matplotlib.figure.Figure):
    
    def _get_axis(self,**kwargs):
        ax = kwargs.pop('ax', self._get_axes()[0])
        return ax
    
    
    def _get_axes(self,**kwargs):        
        ax = kwargs.pop('ax', [])
        ax = kwargs.pop('axs', ax)       
        if ax is None or len(ax)<1:
            try: ax = super(Figure, self)._get_axes()
            except: ax = [plt.gca()]; pass
        return ax
        
        
    def get_legend_handles_labels_axis(self,ax=None,axs=None):
        """Extracts the handles and labels from an axis or from a list of axes. 
        
        Useful for manual legend processing and customisation.
        
        """
        ax = self._get_axes(ax=ax,axs=axs)
        handles = []; labels  = []
        for a in ax:
            handlestmp, labelstmp = a.get_legend_handles_labels()
            handles.extend(handlestmp)
            labels.extend(labelstmp)
        return handles, labels, ax[0]
    
    
    def draw_legend(self, **kwargs):
        """Puts a legend on the provided axis. 
        
        Can be used with kwargs like ncol=2 and alike, which are passed 
        on to the corrresponding pyplot routines.
                
        """
        tc = kwargs.pop('textcolour', matplotlib.rcParams["text.color"])
        tc = kwargs.pop('textcolor', tc)

        #kwargs.setdefault('loc', 0)
        #kwargs.setdefault('frameon', True)

        h, l, a = self.get_legend_handles_labels_axis(ax=kwargs.pop('ax', None),axs=kwargs.pop('axs', None))
        #handles = copy.copy(kwargs.pop('handles', handles))
        handles = []
        for h in kwargs.pop('handles', h):
            handles.append(copy.copy(h))
            handles[-1].set_alpha(1.0)
        labels = []
        for l in kwargs.pop('labels', l):
            labels.append(copy.copy(l))
        legend = a.legend(handles,labels,**kwargs)
        try:
            rect = legend.get_frame()
            rect.set_facecolor(matplotlib.rcParams["grid.color"])
            rect.set_linewidth(0)
            rect.set_edgecolor(tc)
            # Change the alpha value, make sure it is visible
            def set_alpha(objList):
                for o in objList:
                    try: o.set_alpha(1.0)
                    except: matplotlib.artist.setp(o, alpha=1.0); pass
                    #mpl.artist.setp(o, markersize=6)

            #mpl.artist.setp(o, alpha=np.max([1.0,o.get_alpha()]))
            #    h.set_alpha(np.max([1.0,h.get_alpha()]))
            #    #mpl.artist.setp(h, alpha=np.max([1.0,h.get_alpha()]))
            #    mpl.artist.setp(h, markersize=6)
            set_alpha(legend.legendHandles)
            set_alpha(legend.get_lines())
            set_alpha(legend.get_patches())
            #
            #for h in legend.legendHandles:
            #    h.set_alpha(np.max([1.0,h.get_alpha()]))
            #    #mpl.artist.setp(h, alpha=np.max([1.0,h.get_alpha()]))
            #    mpl.artist.setp(h, markersize=6)
            # Change the legend label colors to almost black, too
            for t in legend.texts:
                t.set_color(tc)
        except AttributeError:
            # There are no labled objects
            pass
        return legend
    
    def to_file(self, name, **kwargs):
        dic = dict(bbox_inches='tight')
        dic.update(**kwargs)
        self.savefig(name, **dic)
    
    def to_raster(self, name, **kwargs):
        dic = dict(dpi=300)
        dic.update(**kwargs)
        if name.endswith(".png") or name.endswith(".jpg"):
            self.to_file(name, **dic)
        else:
            raise ValueError("You can only save jpg and png images as raster images.")
        
    def to_power_point(self, name, **kwargs):
        dic = dict(dpi=600, transparent=True)
        dic.update(**kwargs)
        if name.endswith(".png"):
            self.to_raster(name, **dic)
        else:
            raise ValueError("You should use png images with MS PowerPoint.")