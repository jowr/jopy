# -*- coding: utf-8 -*-
from __future__ import print_function, division

from ..base import JopyBaseClass


import matplotlib as mpl
import matplotlib.cm as mplcm
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import brewer2mpl
from itertools import cycle
import sys
import platform

class BaseStyle(JopyBaseClass):
    default_map = "cubehelix_kindl"
    default_lst = 4
    
    def __init__(self):
        JopyBaseClass.__init__(self)
        
        self._figure = None
        self._lines  = {}
        self._ccycle = None # Colour cycle
        self._scycle = None # Style cycle
        self._mcycle = None # Marker cycle
        self._black = "#222222"
        #self._black = "#000000"
        #self._black = 'green'
        self._lgrey = "#FBFBFB" #"GhostWhite"
        self._linewi = 0.75
        self._color_maps = {}
        self._register_color_maps()
        self._color_lists = {}
        self._register_color_lists()
    
    def update_rc_params(self):
        #mpl.rcParams['legend.fontsize'] = 'medium'
        mpl.rcParams['font.size']       = 11.0
        mpl.rcParams['font.serif']      = "Bitstream Vera Serif, New Century Schoolbook, Century Schoolbook L, Utopia, ITC Bookman, Bookman, Nimbus Roman No9 L, Times New Roman, Times, Palatino, Charter, serif"
        mpl.rcParams['font.sans-serif'] = "Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica, Avant Garde, sans-serif"
        mpl.rcParams['font.cursive']    = "Apple Chancery, Textile, Zapf Chancery, Sand, cursive"
        mpl.rcParams['font.fantasy']    = "Comic Sans MS, Chicago, Charcoal, Impact, Western, fantasy"
        mpl.rcParams['font.monospace']  = "Bitstream Vera Sans Mono, Andale Mono, Nimbus Mono L, Courier New, Courier, Fixed, Terminal, monospace"
                
        #mpl.rcParams['text.usetex'] = True
        mpl.rcParams['text.latex.unicode']=True
        mpl.rcParams['text.latex.preamble'] = self._mplpreamble()
        
        mpl.rcParams['mathtext.fontset'] = "cm" # Should be 'cm' (Computer Modern), 'stix', 'stixsans'
        
        if platform.system() == 'Windows':
            mpl.rcParams['pdf.fonttype'] = 42
            mpl.rcParams['ps.fonttype']  = 42

        
        # ######################
        mpl.rcParams["patch.edgecolor"]      = self._black
        mpl.rcParams["text.color"]           = self._black
        mpl.rcParams["axes.edgecolor"]       = self._black   # axes edge color
        mpl.rcParams["axes.linewidth"]       = self._linewi    # edge linewidth
        mpl.rcParams["axes.labelcolor"]      = self._black
        mpl.rcParams["xtick.major.pad"]      = 6      # distance to major tick label in points
        mpl.rcParams["xtick.minor.pad"]      = 6      # distance to the minor tick label in points
        mpl.rcParams["xtick.color"]          = self._black      # color of the tick labels
        mpl.rcParams["ytick.major.pad"]      = 6      # distance to major tick label in points
        mpl.rcParams["ytick.minor.pad"]      = 6      # distance to the minor tick label in points
        mpl.rcParams["ytick.color"]          = self._black      # color of the tick labels
        mpl.rcParams["grid.color"]           = self._lgrey   # grid color
        mpl.rcParams["legend.numpoints"]     = 1      # the number of points in the legend line
        mpl.rcParams["legend.fontsize"]      = "medium"
        mpl.rcParams["legend.scatterpoints"] = 1 # number of scatter points
        mpl.rcParams["path.simplify"]        = True   # When True, simplify paths by removing "invisible" points
        mpl.rcParams["savefig.dpi"]          = 300      # figure dots per inch
        #mpl.rcParams["savefig.format"]       = "pdf"      # png, ps, pdf, svg
        mpl.rcParams["savefig.bbox"]         = "tight" # 'tight' or 'standard'.
        #
        # Update the colours to be black again
        #mpl.rcParams["patch.edgecolor"]      = self._black
        #mpl.rcParams["text.color"]           = self._black
        #mpl.rcParams["axes.edgecolor"]       = self._black   # axes edge color
        #mpl.rcParams["axes.labelcolor"]      = self._black
        #mpl.rcParams["xtick.color"]          = "000000"      # color of the tick labels
        #mpl.rcParams["ytick.color"]          = "000000"      # color of the tick labels
        #
        mpl.rcParams['contour.negative_linestyle'] = 'solid'
        #
        #mpl.use("pgf")
        #mpl.rcParams['text.usetex'] = True
        mpl.rcParams["pgf.texsystem"] = "pdflatex"
        mpl.rcParams["pgf.preamble"]  = self._mplpreamble()
        mpl.rcParams["legend.handlelength"] = 1.75 # the length of the legend lines in fraction of fontsize
        #: [
        #          r"\usepackage[utf8x]{inputenc}",
        #          r"\usepackage[T1]{fontenc}",
        #          r"\usepackage{cmbright}",
        #          ]
        # }
        #mpl.rcParams["axes.color_cycle"]  = ', '.join([self.rgb_to_hex(col) for col in self.color_cycle()])
        #print(self.cycle_to_list(self.color_cycle()))
        mpl.rcParams["axes.color_cycle"]  = self.cycle_to_list(self.color_cycle())
        #axes.color_cycle    : b, g, r, c, m, y, k  # color cycle for plot lines
                                            # as list of string colorspecs:
                                            # single letter, long name, or
                                            # web-style hex
        
    def _mplpreamble(self):
        preamble = []
        preamble.append(r'\usepackage[binary-units=true,abbreviations=true]{siunitx}')
        preamble.append(r'\usepackage{amsmath}')
        #preamble.append(r'\usepackage{amssymb}')
        preamble.append(r'\usepackage[version=3]{mhchem}')
        return preamble 
    
    def _add_to_list(self,objs,lst):
        """Takes a list of object, adds them to lst and returns the new list. If any of the objects was in the list before, it will be removed."""
        # https://docs.python.org/2/library/itertools.html#recipes
        #seen = set()
        #seen_add = seen.add
        #for element in ifilterfalse(seen.__contains__, objs+lst):
        #    seen_add(element)
        #    yield element
        res = []
        for i in objs+lst:
            try: iy = i.strip()
            except: iy = i; pass 
            if iy not in res:
                res.append(iy)
        return res
            
    def _add_to_rc_font_list(self,objs,lst):
        mpl.rcParams[lst] = self._add_to_list(objs.split(','), mpl.rcParams[lst])


    def cycle_to_list(self, cyc):
        """Convert a cycle to a list of unique entries
        
        Takes a cycle object and extract elements until
        it receives and object that has been extracted before
        
        """
        lst = []
        while True:
            obj = next(cyc)
            if obj not in lst: lst.append(obj)
            else: return lst

    def multiply_list(self, lst, doubles=1):
        out = []
        for i in lst:
            for _ in range(doubles):
                out.append(i)
        return out

    def multiply_cycle(self, cyc, doubles=1):
        out = self.cycleToList(cyc)
        out = self.multiplyList(out, doubles=doubles)
        return cycle(out)

    # http://stackoverflow.com/questions/214359/converting-hex-color-to-rgb-and-vice-versa
    def hex_to_rgb(self,hex_string):
        rgb = mpl.colors.hex2color(hex_string)
        return tuple([int(255*x) for x in rgb])

    def rgb_to_hex(self,rgb_tuple):
        return mpl.colors.rgb2hex([1.0*x/255 for x in rgb_tuple])
    
    def _to_rgb(self,inp):
        cc = mpl.colors.ColorConverter()
        return cc.to_rgb(inp) 


    def get_color_map(self,name=default_map,reverse=False):
        """A function to get a matplotlib colour map by name"""
        if name.endswith('_r'):
            name = name[:-2] # remove "_r"
            reverse = not reverse
        # Use the standard maps
        if reverse:
            cm = mplcm.get_cmap(name+'_r')
        else:
            cm = mplcm.get_cmap(name)
        return cm #LinearSegmentedColormap.from_list(cm)

    
    def _register_color_maps(self):
        """A function to create and register the custom colour map objects 
        in a way matplotlib can digest. The cubehelix (including Kindl et al., 
        the Brewer3 colour maps (YlOrRd, PuBuGn, YlGnBu) all provide proper 
        desaturation in grey-scale.
        
        """
        specs = {}
        # We start out with the custom cubehelix maps
        #========= =======================================================
        #Keyword Description
        #========= =======================================================
        #gamma     gamma factor to emphasise either low intensity values
        #          (gamma < 1), or high intensity values (gamma > 1);
        #          defaults to 1.0.
        #s         the start color; defaults to 0.5 (i.e. purple).
        #r         the number of r,g,b rotations in color that are made
        #          from the start to the end of the color scheme; defaults
        #          to -1.5 (i.e. -> B -> G -> R -> B).
        #h         the hue parameter which controls how saturated the
        #          colors are. If this parameter is zero then the color
        #          scheme is purely a greyscale; defaults to 1.0.
        #========= =======================================================
        # 0 = blue, 1 = red, 2 = green
        specs['cubehelix_alt']   = mpl._cm.cubehelix(h=1.5) # standard colours but more intensity
        specs['cubehelix_blue']  = mpl._cm.cubehelix(s=0.3,r=-0.5,h=1.5) # blue colours and higher intensity
        specs['cubehelix_red']   = mpl._cm.cubehelix(s=1.3,r=-0.5,h=1.5) # blue colours and higher intensity
        specs['cubehelix_green'] = mpl._cm.cubehelix(s=2.3,r=-0.5,h=1.5) # blue colours and higher intensity
        specs['cubehelix_kindl'] = mpl._cm.cubehelix(gamma=1.4,s=0.4,r=-0.8,h=2.0) # blue colours and higher intensity

        # A Python version of Matteo Niccoli's colourmaps
        # http://nbviewer.ipython.org/github/kwinkunks/notebooks/blob/master/Matteo_colourmaps.ipynb
        #

        ## The original data is approximated with polynomials
        #p0 = np.array([ -3.00669779e-36, 6.26525719e-33, -5.87240910e-30, 3.25751282e-27, -1.18087586e-24, 2.89863878e-22, -4.75604889e-20, 4.67614077e-18, -1.13599364e-16, -4.16063333e-14, 7.27326802e-12, -6.41161566e-10, 3.52560300e-08, -1.23850903e-06, 2.67527478e-05, -3.20955377e-04, 1.51205802e-03, 8.78290363e-03, 2.40494252e-02])
        #p1 = np.array([  8.13719543e-37, -1.77388587e-33, 1.75833582e-30, -1.04750030e-27, 4.17412745e-25, -1.17054612e-22, 2.36196641e-20, -3.43234035e-18, 3.50200516e-16, -2.30787699e-14, 6.79825384e-13, 3.17731005e-11, -4.59193023e-09, 2.36050767e-07, -6.49812536e-06, 9.42940406e-05, -6.24155259e-04, 4.04636648e-03, 1.11088863e-02])
        #p2 = np.array([ -1.85874163e-35, 4.32740491e-32, -4.60790627e-29, 2.97271126e-26, -1.29744258e-23, 4.05138291e-21, -9.33419995e-19, 1.61295373e-16, -2.10418623e-14, 2.06972791e-12, -1.52201298e-10, 8.23377786e-09, -3.19603306e-07, 8.58912760e-06, -1.52305419e-04, 1.67708019e-03, -1.05536314e-02, 3.80789592e-02, 5.82194596e-03])
        #x  = range(256)
        #LinL = np.empty((256,3))
        #LinL[:,0] = np.polyval(p0,x)
        #LinL[:,1] = np.polyval(p1,x)
        #LinL[:,2] = np.polyval(p2,x)
        #LinL = np.clip(LinL,0,1)
        LinL = np.array([[  1.43000000e-02, 1.43000000e-02, 1.43000000e-02],
 [  4.04000000e-02, 1.25000000e-02, 3.25000000e-02],
 [  5.16000000e-02, 1.54000000e-02, 4.43000000e-02],
 [  6.16000000e-02, 1.84000000e-02, 5.30000000e-02],
 [  6.99000000e-02, 2.15000000e-02, 6.15000000e-02],
 [  8.14000000e-02, 2.29000000e-02, 6.87000000e-02],
 [  8.57000000e-02, 2.73000000e-02, 7.63000000e-02],
 [  9.28000000e-02, 3.05000000e-02, 8.05000000e-02],
 [  1.00800000e-01, 3.30000000e-02, 8.46000000e-02],
 [  1.06400000e-01, 3.56000000e-02, 9.39000000e-02],
 [  1.11000000e-01, 3.96000000e-02, 9.79000000e-02],
 [  1.18400000e-01, 4.15000000e-02, 1.02000000e-01],
 [  1.22400000e-01, 4.53000000e-02, 1.06200000e-01],
 [  1.26300000e-01, 4.89000000e-02, 1.10500000e-01],
 [  1.30800000e-01, 5.18000000e-02, 1.15000000e-01],
 [  1.35700000e-01, 5.41000000e-02, 1.20000000e-01],
 [  1.41300000e-01, 5.55000000e-02, 1.25600000e-01],
 [  1.45300000e-01, 5.90000000e-02, 1.25600000e-01],
 [  1.50600000e-01, 6.03000000e-02, 1.30900000e-01],
 [  1.53300000e-01, 6.37000000e-02, 1.34400000e-01],
 [  1.56900000e-01, 6.64000000e-02, 1.38500000e-01],
 [  1.62400000e-01, 6.78000000e-02, 1.42500000e-01],
 [  1.65400000e-01, 7.08000000e-02, 1.47100000e-01],
 [  1.70300000e-01, 7.29000000e-02, 1.50400000e-01],
 [  1.74200000e-01, 7.55000000e-02, 1.54200000e-01],
 [  1.79600000e-01, 7.68000000e-02, 1.59500000e-01],
 [  1.80900000e-01, 7.83000000e-02, 1.77500000e-01],
 [  1.79800000e-01, 8.07000000e-02, 1.96700000e-01],
 [  1.78300000e-01, 8.29000000e-02, 2.15900000e-01],
 [  1.78200000e-01, 8.43000000e-02, 2.34100000e-01],
 [  1.76500000e-01, 8.66000000e-02, 2.51400000e-01],
 [  1.77000000e-01, 8.86000000e-02, 2.64600000e-01],
 [  1.76100000e-01, 9.11000000e-02, 2.78200000e-01],
 [  1.75400000e-01, 9.33000000e-02, 2.92200000e-01],
 [  1.77300000e-01, 9.51000000e-02, 3.02600000e-01],
 [  1.75800000e-01, 9.72000000e-02, 3.17400000e-01],
 [  1.75100000e-01, 9.99000000e-02, 3.29000000e-01],
 [  1.74400000e-01, 1.02500000e-01, 3.40500000e-01],
 [  1.73500000e-01, 1.04700000e-01, 3.53400000e-01],
 [  1.74800000e-01, 1.07100000e-01, 3.62700000e-01],
 [  1.74700000e-01, 1.09400000e-01, 3.73900000e-01],
 [  1.72600000e-01, 1.12200000e-01, 3.85800000e-01],
 [  1.73100000e-01, 1.15300000e-01, 3.94000000e-01],
 [  1.73100000e-01, 1.17500000e-01, 4.05100000e-01],
 [  1.73100000e-01, 1.19700000e-01, 4.16100000e-01],
 [  1.72000000e-01, 1.22400000e-01, 4.26800000e-01],
 [  1.73000000e-01, 1.26100000e-01, 4.33000000e-01],
 [  1.71600000e-01, 1.28500000e-01, 4.44500000e-01],
 [  1.71000000e-01, 1.31400000e-01, 4.54000000e-01],
 [  1.70600000e-01, 1.34000000e-01, 4.64200000e-01],
 [  1.66400000e-01, 1.38100000e-01, 4.73900000e-01],
 [  1.58200000e-01, 1.48200000e-01, 4.71700000e-01],
 [  1.48700000e-01, 1.58300000e-01, 4.68300000e-01],
 [  1.42200000e-01, 1.65300000e-01, 4.69900000e-01],
 [  1.35200000e-01, 1.72900000e-01, 4.69400000e-01],
 [  1.28100000e-01, 1.79600000e-01, 4.70800000e-01],
 [  1.25400000e-01, 1.85400000e-01, 4.71900000e-01],
 [  1.20700000e-01, 1.90400000e-01, 4.76200000e-01],
 [  1.16700000e-01, 1.96200000e-01, 4.77300000e-01],
 [  1.16600000e-01, 2.00400000e-01, 4.81400000e-01],
 [  1.14700000e-01, 2.05700000e-01, 4.82300000e-01],
 [  1.13700000e-01, 2.09300000e-01, 4.88800000e-01],
 [  1.09100000e-01, 2.14900000e-01, 4.90400000e-01],
 [  1.08900000e-01, 2.19000000e-01, 4.94400000e-01],
 [  1.07400000e-01, 2.23400000e-01, 4.98400000e-01],
 [  1.10800000e-01, 2.27000000e-01, 5.02200000e-01],
 [  1.09600000e-01, 2.31500000e-01, 5.05600000e-01],
 [  1.05800000e-01, 2.36700000e-01, 5.07700000e-01],
 [  1.04800000e-01, 2.40200000e-01, 5.14500000e-01],
 [  1.04700000e-01, 2.44400000e-01, 5.18200000e-01],
 [  1.06800000e-01, 2.48200000e-01, 5.22300000e-01],
 [  1.08600000e-01, 2.52100000e-01, 5.26400000e-01],
 [  1.06700000e-01, 2.56800000e-01, 5.29000000e-01],
 [  1.06100000e-01, 2.60700000e-01, 5.34600000e-01],
 [  1.05900000e-01, 2.64800000e-01, 5.38600000e-01],
 [  1.05600000e-01, 2.69000000e-01, 5.42700000e-01],
 [  9.69000000e-02, 2.76900000e-01, 5.34300000e-01],
 [  8.79000000e-02, 2.84700000e-01, 5.25100000e-01],
 [  8.32000000e-02, 2.90700000e-01, 5.21800000e-01],
 [  7.93000000e-02, 2.97200000e-01, 5.15300000e-01],
 [  6.86000000e-02, 3.04400000e-01, 5.06800000e-01],
 [  6.39000000e-02, 3.10600000e-01, 5.00600000e-01],
 [  5.86000000e-02, 3.16600000e-01, 4.95500000e-01],
 [  5.36000000e-02, 3.22700000e-01, 4.88800000e-01],
 [  3.88000000e-02, 3.29200000e-01, 4.81700000e-01],
 [  4.09000000e-02, 3.34300000e-01, 4.78600000e-01],
 [  3.45000000e-02, 3.40100000e-01, 4.72200000e-01],
 [  3.00000000e-02, 3.45800000e-01, 4.66500000e-01],
 [  2.90000000e-02, 3.50600000e-01, 4.64700000e-01],
 [  2.26000000e-02, 3.56400000e-01, 4.57800000e-01],
 [  1.54000000e-02, 3.61900000e-01, 4.52900000e-01],
 [  1.46000000e-02, 3.67000000e-01, 4.48700000e-01],
 [  1.69000000e-02, 3.71800000e-01, 4.46400000e-01],
 [  1.17000000e-02, 3.77300000e-01, 4.40000000e-01],
 [  5.50000000e-03, 3.82800000e-01, 4.33400000e-01],
 [  5.20000000e-03, 3.87300000e-01, 4.32700000e-01],
 [  8.00000000e-04, 3.92700000e-01, 4.26700000e-01],
 [  0.00000000e+00, 3.97700000e-01, 4.22000000e-01],
 [  3.00000000e-04, 4.02400000e-01, 4.20000000e-01],
 [  1.30000000e-03, 4.07200000e-01, 4.16600000e-01],
 [  0.00000000e+00, 4.13200000e-01, 4.10700000e-01],
 [  0.00000000e+00, 4.17900000e-01, 4.07100000e-01],
 [  0.00000000e+00, 4.24200000e-01, 3.97700000e-01],
 [  0.00000000e+00, 4.30200000e-01, 3.91900000e-01],
 [  0.00000000e+00, 4.34500000e-01, 3.89000000e-01],
 [  0.00000000e+00, 4.39500000e-01, 3.84900000e-01],
 [  0.00000000e+00, 4.45600000e-01, 3.77600000e-01],
 [  0.00000000e+00, 4.50800000e-01, 3.72800000e-01],
 [  0.00000000e+00, 4.56300000e-01, 3.66600000e-01],
 [  0.00000000e+00, 4.61300000e-01, 3.59700000e-01],
 [  0.00000000e+00, 4.67000000e-01, 3.54200000e-01],
 [  0.00000000e+00, 4.71600000e-01, 3.50400000e-01],
 [  0.00000000e+00, 4.76300000e-01, 3.46400000e-01],
 [  0.00000000e+00, 4.82000000e-01, 3.37500000e-01],
 [  0.00000000e+00, 4.87000000e-01, 3.33100000e-01],
 [  0.00000000e+00, 4.91800000e-01, 3.25600000e-01],
 [  0.00000000e+00, 4.96900000e-01, 3.19800000e-01],
 [  0.00000000e+00, 5.02400000e-01, 3.12600000e-01],
 [  0.00000000e+00, 5.06000000e-01, 3.10100000e-01],
 [  0.00000000e+00, 5.11800000e-01, 3.01200000e-01],
 [  0.00000000e+00, 5.16600000e-01, 2.93800000e-01],
 [  0.00000000e+00, 5.21400000e-01, 2.87100000e-01],
 [  0.00000000e+00, 5.26800000e-01, 2.81600000e-01],
 [  0.00000000e+00, 5.31500000e-01, 2.72600000e-01],
 [  0.00000000e+00, 5.35600000e-01, 2.67500000e-01],
 [  0.00000000e+00, 5.40700000e-01, 2.59700000e-01],
 [  0.00000000e+00, 5.46100000e-01, 2.55200000e-01],
 [  0.00000000e+00, 5.51900000e-01, 2.50600000e-01],
 [  0.00000000e+00, 5.56500000e-01, 2.46900000e-01],
 [  0.00000000e+00, 5.62400000e-01, 2.39600000e-01],
 [  0.00000000e+00, 5.67800000e-01, 2.36000000e-01],
 [  0.00000000e+00, 5.72700000e-01, 2.33800000e-01],
 [  0.00000000e+00, 5.77800000e-01, 2.28700000e-01],
 [  0.00000000e+00, 5.82900000e-01, 2.25000000e-01],
 [  0.00000000e+00, 5.88300000e-01, 2.18000000e-01],
 [  0.00000000e+00, 5.93100000e-01, 2.14600000e-01],
 [  0.00000000e+00, 5.99100000e-01, 2.08900000e-01],
 [  0.00000000e+00, 6.03600000e-01, 2.05600000e-01],
 [  0.00000000e+00, 6.08400000e-01, 1.99900000e-01],
 [  0.00000000e+00, 6.13800000e-01, 1.96100000e-01],
 [  0.00000000e+00, 6.18900000e-01, 1.89900000e-01],
 [  0.00000000e+00, 6.23900000e-01, 1.84800000e-01],
 [  0.00000000e+00, 6.29000000e-01, 1.75900000e-01],
 [  0.00000000e+00, 6.34500000e-01, 1.70700000e-01],
 [  0.00000000e+00, 6.38100000e-01, 1.63800000e-01],
 [  0.00000000e+00, 6.43800000e-01, 1.59200000e-01],
 [  0.00000000e+00, 6.48900000e-01, 1.51900000e-01],
 [  0.00000000e+00, 6.53600000e-01, 1.41000000e-01],
 [  0.00000000e+00, 6.59000000e-01, 1.32200000e-01],
 [  0.00000000e+00, 6.64500000e-01, 1.22200000e-01],
 [  0.00000000e+00, 6.65600000e-01, 9.09000000e-02],
 [  0.00000000e+00, 6.64400000e-01, 3.22000000e-02],
 [  3.51000000e-02, 6.66000000e-01, 0.00000000e+00],
 [  7.97000000e-02, 6.70300000e-01, 0.00000000e+00],
 [  1.12900000e-01, 6.73900000e-01, 0.00000000e+00],
 [  1.39200000e-01, 6.77600000e-01, 0.00000000e+00],
 [  1.56600000e-01, 6.81400000e-01, 0.00000000e+00],
 [  1.76500000e-01, 6.84700000e-01, 0.00000000e+00],
 [  1.89000000e-01, 6.89000000e-01, 0.00000000e+00],
 [  2.03000000e-01, 6.92800000e-01, 0.00000000e+00],
 [  2.16700000e-01, 6.96600000e-01, 0.00000000e+00],
 [  2.29900000e-01, 7.00300000e-01, 0.00000000e+00],
 [  2.39100000e-01, 7.04400000e-01, 0.00000000e+00],
 [  2.51700000e-01, 7.08100000e-01, 0.00000000e+00],
 [  2.57400000e-01, 7.12400000e-01, 0.00000000e+00],
 [  2.67900000e-01, 7.16200000e-01, 0.00000000e+00],
 [  2.79000000e-01, 7.20000000e-01, 0.00000000e+00],
 [  2.87800000e-01, 7.24000000e-01, 0.00000000e+00],
 [  2.96500000e-01, 7.28000000e-01, 0.00000000e+00],
 [  3.05200000e-01, 7.31900000e-01, 0.00000000e+00],
 [  3.10100000e-01, 7.36200000e-01, 0.00000000e+00],
 [  3.18700000e-01, 7.40200000e-01, 0.00000000e+00],
 [  3.27200000e-01, 7.44100000e-01, 0.00000000e+00],
 [  3.34500000e-01, 7.48200000e-01, 0.00000000e+00],
 [  3.40600000e-01, 7.52300000e-01, 0.00000000e+00],
 [  3.60400000e-01, 7.54900000e-01, 0.00000000e+00],
 [  3.89800000e-01, 7.56300000e-01, 0.00000000e+00],
 [  4.16900000e-01, 7.57400000e-01, 0.00000000e+00],
 [  4.46100000e-01, 7.58000000e-01, 0.00000000e+00],
 [  4.68100000e-01, 7.59400000e-01, 0.00000000e+00],
 [  4.90000000e-01, 7.61200000e-01, 0.00000000e+00],
 [  5.08900000e-01, 7.62700000e-01, 0.00000000e+00],
 [  5.30400000e-01, 7.63700000e-01, 0.00000000e+00],
 [  5.50000000e-01, 7.64900000e-01, 0.00000000e+00],
 [  5.69800000e-01, 7.66000000e-01, 0.00000000e+00],
 [  5.82500000e-01, 7.68800000e-01, 0.00000000e+00],
 [  5.99900000e-01, 7.70100000e-01, 0.00000000e+00],
 [  6.17300000e-01, 7.71300000e-01, 0.00000000e+00],
 [  6.31400000e-01, 7.73000000e-01, 0.00000000e+00],
 [  6.48700000e-01, 7.74100000e-01, 0.00000000e+00],
 [  6.63200000e-01, 7.76300000e-01, 0.00000000e+00],
 [  6.75700000e-01, 7.78200000e-01, 0.00000000e+00],
 [  6.91200000e-01, 7.79500000e-01, 0.00000000e+00],
 [  7.06100000e-01, 7.80800000e-01, 0.00000000e+00],
 [  7.22200000e-01, 7.81800000e-01, 0.00000000e+00],
 [  7.30500000e-01, 7.85200000e-01, 0.00000000e+00],
 [  7.44200000e-01, 7.86600000e-01, 0.00000000e+00],
 [  7.58000000e-01, 7.88000000e-01, 0.00000000e+00],
 [  7.70900000e-01, 7.89600000e-01, 0.00000000e+00],
 [  7.83300000e-01, 7.91500000e-01, 0.00000000e+00],
 [  7.87200000e-01, 7.89100000e-01, 9.51000000e-02],
 [  7.97200000e-01, 7.90300000e-01, 1.98800000e-01],
 [  8.07200000e-01, 7.91700000e-01, 2.56000000e-01],
 [  8.11600000e-01, 7.94900000e-01, 3.00100000e-01],
 [  8.21100000e-01, 7.96400000e-01, 3.39700000e-01],
 [  8.30800000e-01, 7.98000000e-01, 3.71000000e-01],
 [  8.35000000e-01, 8.01100000e-01, 4.02800000e-01],
 [  8.45000000e-01, 8.02600000e-01, 4.29200000e-01],
 [  8.54800000e-01, 8.04100000e-01, 4.55500000e-01],
 [  8.60200000e-01, 8.07300000e-01, 4.73500000e-01],
 [  8.65800000e-01, 8.10000000e-01, 4.99300000e-01],
 [  8.75800000e-01, 8.11600000e-01, 5.18400000e-01],
 [  8.85600000e-01, 8.13000000e-01, 5.40200000e-01],
 [  8.89500000e-01, 8.16400000e-01, 5.60200000e-01],
 [  8.99400000e-01, 8.18000000e-01, 5.77500000e-01],
 [  9.07700000e-01, 8.20200000e-01, 5.91800000e-01],
 [  9.10600000e-01, 8.24100000e-01, 6.09400000e-01],
 [  9.20500000e-01, 8.25700000e-01, 6.25300000e-01],
 [  9.28400000e-01, 8.27800000e-01, 6.42000000e-01],
 [  9.34300000e-01, 8.30700000e-01, 6.57600000e-01],
 [  9.41500000e-01, 8.32900000e-01, 6.76200000e-01],
 [  9.51200000e-01, 8.34800000e-01, 6.86800000e-01],
 [  9.54900000e-01, 8.38400000e-01, 7.02600000e-01],
 [  9.62200000e-01, 8.40800000e-01, 7.17000000e-01],
 [  9.71200000e-01, 8.42900000e-01, 7.28700000e-01],
 [  9.70800000e-01, 8.48200000e-01, 7.40900000e-01],
 [  9.71300000e-01, 8.53000000e-01, 7.55500000e-01],
 [  9.69100000e-01, 8.59100000e-01, 7.65500000e-01],
 [  9.69900000e-01, 8.64200000e-01, 7.74600000e-01],
 [  9.70300000e-01, 8.69100000e-01, 7.87100000e-01],
 [  9.71000000e-01, 8.74000000e-01, 7.99900000e-01],
 [  9.69500000e-01, 8.80000000e-01, 8.06700000e-01],
 [  9.69600000e-01, 8.85100000e-01, 8.18800000e-01],
 [  9.68600000e-01, 8.90800000e-01, 8.27800000e-01],
 [  9.68100000e-01, 8.96200000e-01, 8.37800000e-01],
 [  9.68800000e-01, 9.01300000e-01, 8.46700000e-01],
 [  9.69600000e-01, 9.06400000e-01, 8.55700000e-01],
 [  9.70300000e-01, 9.11500000e-01, 8.64700000e-01],
 [  9.70800000e-01, 9.16300000e-01, 8.77300000e-01],
 [  9.69100000e-01, 9.22400000e-01, 8.83800000e-01],
 [  9.69200000e-01, 9.27300000e-01, 8.96100000e-01],
 [  9.69900000e-01, 9.32300000e-01, 9.05100000e-01],
 [  9.69300000e-01, 9.38100000e-01, 9.10800000e-01],
 [  9.71400000e-01, 9.42500000e-01, 9.23000000e-01],
 [  9.71200000e-01, 9.47800000e-01, 9.31100000e-01],
 [  9.70000000e-01, 9.53700000e-01, 9.38100000e-01],
 [  9.70700000e-01, 9.58700000e-01, 9.47000000e-01],
 [  9.71300000e-01, 9.63800000e-01, 9.56000000e-01],
 [  9.72600000e-01, 9.68700000e-01, 9.64800000e-01],
 [  9.73800000e-01, 9.73800000e-01, 9.71100000e-01],
 [  9.78000000e-01, 9.78000000e-01, 9.78000000e-01],
 [  9.82400000e-01, 9.82400000e-01, 9.82400000e-01],
 [  9.86800000e-01, 9.86800000e-01, 9.86800000e-01],
 [  9.91200000e-01, 9.91200000e-01, 9.91200000e-01],
 [  9.95600000e-01, 9.95600000e-01, 9.95600000e-01],
 [  1.00000000e+00, 1.00000000e+00, 1.00000000e+00]])
        b3 = LinL[:,2] # value of blue at sample n
        b2 = LinL[:,2] # value of blue at sample n
        b1 = np.linspace(0, 1, len(b2)) # position of sample n - ranges from 0 to 1

        # Setting up columns for tuples
        g3 = LinL[:,1]
        g2 = LinL[:,1]
        g1 = np.linspace(0,1,len(g2))

        r3 = LinL[:,0]
        r2 = LinL[:,0]
        r1 = np.linspace(0,1,len(r2))

        # Creating tuples
        R = zip(r1,r2,r3)
        G = zip(g1,g2,g3)
        B = zip(b1,b2,b3)

        # Transposing
        RGB = zip(R,G,B)
        rgb = zip(*RGB)

        # Creating dictionary
        k = ['red', 'green', 'blue']
        specs['matteoniccoli'] = dict(zip(k,rgb)) 
        
        for name in specs:
            mplcm.register_cmap(name=name, data=specs[name])
            mplcm.register_cmap(name=name+"_r", data=mplcm._reverse_cmap_spec(specs[name]))
            self._color_maps[name] = self.get_color_map(name)
            self._color_maps[name+"_r"] = self.get_color_map(name+"_r")
            
    def _register_color_lists(self, length=default_lst):
        cc = mpl.colors.ColorConverter()
        self._color_lists['matplotlib'] = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        self._color_lists['simple']     = ['#FF0000', '#FFD300', '#3914AF', '#00CC00']
        self._color_lists['brewer']     = ['#e66101', '#fdb863', '#b2abd2', '#5e3c99']
        for i in self._color_lists:
            self._color_lists[i] = [cc.to_rgb(item) for item in self._color_lists[i]]
        self._color_lists['cblind']     = [(0.9,0.6,0.0), (0.35,0.70,0.90), (0.00,0.60,0.50), (0.00,0.45,0.70), (0.80,0.40,0.00), (0.80,0.60,0.70)]
        self._color_lists['brewer1']    = brewer2mpl.get_map('Set1', 'qualitative', length).mpl_colors
        self._color_lists['brewer2']    = brewer2mpl.get_map('Set2', 'qualitative', length).mpl_colors
        

    def color_cycle(self, name='map', cmap=None, length=None):
        '''Returns the current colour cycle, creates a new one if necessary.
        
        Parameters
        ----------
        name : str
            selector for colour cycle:
            'cmap': use a colourmap to generate cycle, see http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
            'DTU': uses the colours from the DTU design guide
            'DTU_dark': uses the colours from the DTU design guide, darker colours first, good for presentations
            'cblind': A scheme for colourblind people
            'matplotlib': The standard matplotlib scheme
            'simple': A simple four-colour scheme that work for greyscale and colourblind people
            'brewer1' and 'brewer2': See http://colorbrewer2.org, works for colourblind and greyscale
            
        '''
        if name=='map':
            if length is None: length = self.default_lst
            if cmap is None: cmap = self.default_map
            cmap  = self.get_color_map(cmap)
            clist = [cmap(i) for i in np.linspace(0.25, 0.75, length)]
        else:
            clist = self._color_lists[name]
            
        if length is not None:
            if length<1: 
                return cycle(['none'])
            elif length<=len(clist):
                return cycle(list(clist)[0:length])
            elif length>len(clist):
                self.autolog("Colour cycle is too short, cannot extend it.")
        return cycle(clist)
    
    def _show_info(self,show=True):
        self.update_rc_params()
        lsts = self._color_lists.keys() 
        l = len(lsts) 
        import matplotlib.pyplot as plt
        line_fig = plt.figure()
        xdata = np.linspace(0,6)
        for i, m in enumerate(lsts):
            plt.subplot(1,l,i+1)
            plt.axis("off")
            for j in self._color_lists[m]:
                plt.plot(xdata,np.random.normal(size=len(xdata)),lw=1.5,color=j)
            plt.title(m,rotation=45)
        plt.tight_layout()

        xdata=np.outer(np.arange(0,1,0.01),np.ones(10))
        maps = [m for m in self._color_maps.keys() if not m.endswith("_r")]
        l=len(maps)
        map_fig = plt.figure()
        for i, m in enumerate(maps):
            plt.subplot(1,l,i+1)
            plt.axis("off")
            plt.imshow(xdata,aspect='auto',cmap=plt.get_cmap(m),origin="lower")
            plt.title(m,rotation=45)
        plt.tight_layout()
        
        if show: plt.show()
        #else: plt.close('all')
        return line_fig,map_fig


class DtuStyle(BaseStyle):
    def __init__(self):
        BaseStyle.__init__(self)
        
    def update_rc_params(self):
        BaseStyle.update_rc_params(self)
        self._add_to_rc_font_list('Myriad Pro, Arial, Verdana','font.sans-serif')
        self._add_to_rc_font_list('Minion Pro, Times New Roman, Times','font.serif') 
        #mpl.rcParams['font.monospace'] = 'CMU Typewriter Text'
         
        #mpl.rcParams['font.family'] = 'sans-serif'
        #mpl.rcParams['mathtext.it'] = 'Myriad Pro:italic'
        #mpl.rcParams['mathtext.bf'] = 'Myriad Pro:bold'
         
        mpl.rcParams['font.family'] = 'serif'
        #mpl.rcParams['mathtext.it'] = 'Minion Pro:italic'
        #mpl.rcParams['mathtext.bf'] = 'Minion Pro:bold'
         
        mpl.rcParams['mathtext.fontset'] = "custom" # Should be 'cm' (Computer Modern), 'stix', 'stixsans'
        mpl.rcParams['mathtext.cal'] = 'serif:italic'
        mpl.rcParams['mathtext.rm']  = 'serif'
        mpl.rcParams['mathtext.tt']  = 'monospace'
        mpl.rcParams['mathtext.it']  = 'serif:italic'
        mpl.rcParams['mathtext.bf']  = 'serif:bold'
        mpl.rcParams['mathtext.sf']  = 'sans'
         
         
    def _mplpreamble(self):
        preamble = BaseStyle._mplpreamble(self)
        preamble.append(r'\usepackage{MnSymbol}')
        preamble.append(r'\usepackage[lf]{MyriadPro} % Sans font')
        preamble.append(r'\usepackage[lf]{MinionPro} % Serif font')
        return preamble 
    
    def _register_color_maps(self):
        BaseStyle._register_color_maps(self)
        
        cc    = mpl.colors.ColorConverter()
        redwhite    = [cc.to_rgba('#FF0000',alpha=1),cc.to_rgba('#FFFFFF',alpha=1)]
        greenwhite  = [cc.to_rgba('#99CC33',alpha=1),cc.to_rgba('#FFFFFF',alpha=1)]
        purplewhite = [cc.to_rgba('#660099',alpha=1),cc.to_rgba('#FFFFFF',alpha=1)]
        yellowwhite = [cc.to_rgba('#FFCC00',alpha=1),cc.to_rgba('#FFFFFF',alpha=1)]
        # create map and register it together with reversed colours
        maps = []
        maps.append(LinearSegmentedColormap.from_list('DTU1', redwhite))
        maps.append(LinearSegmentedColormap.from_list('DTU1_r', redwhite[::-1]))
        maps.append(LinearSegmentedColormap.from_list('DTU2', greenwhite))
        maps.append(LinearSegmentedColormap.from_list('DTU2_r', greenwhite[::-1]))
        maps.append(LinearSegmentedColormap.from_list('DTU3', purplewhite))
        maps.append(LinearSegmentedColormap.from_list('DTU3_r', purplewhite[::-1]))
        maps.append(LinearSegmentedColormap.from_list('DTU4', yellowwhite))
        maps.append(LinearSegmentedColormap.from_list('DTU4_r', yellowwhite[::-1]))
        
        for cmap in maps:
            mplcm.register_cmap(cmap=cmap)
            self._color_maps[cmap.name] = cmap

    def _register_color_lists(self, length=BaseStyle.default_lst):
        BaseStyle._register_color_lists(self)
        cc = mpl.colors.ColorConverter()
        self._color_lists['DTU']      = [ cc.to_rgb(item) for item in ['#FF0000', '#99CC33', '#660099', '#FFCC00', '#999999', '#000000', '#33CCFF', '#3366CC', '#FF9900', '#CC3399', '#66CC00'] ]
        self._color_lists['DTU_dark'] = [ cc.to_rgb(item) for item in ['#FF0000', '#660099', '#99CC33', '#3366CC', '#999999', '#FFCC00', '#000000', '#33CCFF', '#FF9900', '#CC3399', '#66CC00'] ]
        

class IpuStyle(BaseStyle):
    
    default_map = "IPU"
    
    def __init__(self):
        BaseStyle.__init__(self)
        
    def update_rc_params(self):
        BaseStyle.update_rc_params(self)
        self._add_to_rc_font_list('Helvetica, Arial, cmbright, Verdana','font.sans-serif')
        self._add_to_rc_font_list('Times New Roman, Times','font.serif') 
        self._add_to_rc_font_list('sourcecodepro','font.monospace')
        mpl.rcParams['font.family'] = 'sans-serif'
        #mpl.rcParams['font.monospace'] = 'CMU Typewriter Text'
        mpl.rcParams['mathtext.fontset'] = "custom" # Should be 'cm' (Computer Modern), 'stix', 'stixsans'
        mpl.rcParams['mathtext.cal'] = 'sans:italic'
        mpl.rcParams['mathtext.rm']  = 'sans'
        mpl.rcParams['mathtext.tt']  = 'monospace'
        mpl.rcParams['mathtext.it']  = 'sans:italic'
        mpl.rcParams['mathtext.bf']  = 'sans:bold'
        mpl.rcParams['mathtext.sf']  = 'sans'
                
         
    def _mplpreamble(self):
        preamble = BaseStyle._mplpreamble(self)
        preamble.append(r'\usepackage{cmbright}')
        preamble.append(r'\usepackage{helvet}')
        #preamble.append(r'\usepackage{sansmath}')
        #preamble.append(r'\sansmath')
        #preamble.append(r'\renewcommand*{\familydefault}{\sfdefault}')
        return preamble 
    
    def _register_color_maps(self):
        BaseStyle._register_color_maps(self)
        rgb = [
          (  0./255. ,   0./255. ,   0./255.),     
          (  0./255. , 102./255. ,  51./255.),
          (114./255. , 121./255. , 126./255.),
          ( 91./255. , 172./255. ,  38./255.),
          (217./255. , 220./255. , 222./255.),
          (255./255. , 255./255. , 255./255.)]

        # create map and register it together with reversed colours
        maps = []
        maps.append(LinearSegmentedColormap.from_list('IPU'  , rgb))
        maps.append(LinearSegmentedColormap.from_list('IPU_r', rgb[::-1]))
    
        for cmap in maps:
            mplcm.register_cmap(cmap=cmap)
            self._color_maps[cmap.name] = cmap
    
    def _register_color_lists(self, length=BaseStyle.default_lst):
        BaseStyle._register_color_lists(self)
        self._color_lists['IPU']      = [
          (  0./255. ,   0./255. ,   0./255.),     
          (  0./255. , 102./255. ,  51./255.),
          (114./255. , 121./255. , 126./255.),
          ( 91./255. , 172./255. ,  38./255.),
          (217./255. , 220./255. , 222./255.),
          #(255./255. , 255./255. , 255./255.)
        ]
        
    def color_cycle(self, name='IPU', cmap=None, length=None):
        return BaseStyle.color_cycle(self, name, map, length)


class VdgStyle(BaseStyle):
    def __init__(self):
        BaseStyle.__init__(self)
        mpl.rcParams['lines.linewidth'] = 1.25     # line width in points
        
    def update_rc_params(self):
        BaseStyle.update_rc_params(self)
        self._add_to_rc_font_list('Helvetica, Arial, cmbright, Verdana','font.sans-serif')
        self._add_to_rc_font_list('Times New Roman, Times','font.serif')
        self._add_to_rc_font_list('sourcecodepro','font.monospace')
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['mathtext.fontset'] = "custom" # Should be 'cm' (Computer Modern), 'stix', 'stixsans'
        mpl.rcParams['mathtext.cal'] = 'sans:italic'
        mpl.rcParams['mathtext.rm']  = 'sans'
        mpl.rcParams['mathtext.tt']  = 'monospace'
        mpl.rcParams['mathtext.it']  = 'sans:italic'
        mpl.rcParams['mathtext.bf']  = 'sans:bold'
        mpl.rcParams['mathtext.sf']  = 'sans'
                
         
    def _mplpreamble(self):
        preamble = BaseStyle._mplpreamble(self)
        preamble.append(r'\usepackage{cmbright}')
        preamble.append(r'\usepackage{helvet}')
        #preamble.append(r'\usepackage{sansmath}')
        #preamble.append(r'\sansmath')
        #preamble.append(r'\renewcommand*{\familydefault}{\sfdefault}')
        return preamble 

    def _register_color_lists(self, length=BaseStyle.default_lst):
        BaseStyle._register_color_lists(self)
        self._color_lists['VDG']      = [
          (  0./255. ,  51./255. , 102./255.),     
          (205./255. ,  51./255. ,  51./255.),
          (  0./255. , 109./255. , 148./255.),
          (127./255. , 186./255. ,  50./255.),
          #(255./255. , 255./255. , 255./255.)
        ]
        
    def color_cycle(self, name='VDG', cmap=None, length=None):
        return BaseStyle.color_cycle(self, name, map, length)
