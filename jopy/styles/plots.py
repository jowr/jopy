# -*- coding: utf-8 -*-
from __future__ import print_function, division
import matplotlib

class Figure(matplotlib.figure.Figure):
    
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