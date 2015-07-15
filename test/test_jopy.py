"""
Tests for `jopy` module.
"""
import matplotlib; matplotlib.use('Agg')
from matplotlib.figure import Figure

from jopy.recip.mechanisms import RecipExplicit, RecipImplicit, RecipBase
import numpy as np
from jopy.styles.mplib import BaseStyle, DtuStyle, IpuStyle

LOCAL=True


class TestJopy(object):

    @classmethod
    def setup_class(cls):
        pass

    def test_something(self):
        pass

    @classmethod
    def teardown_class(cls):
        pass

class TestJopyRecip(object):
    exp = RecipExplicit()
    imp = RecipImplicit()

    @classmethod
    def setup_class(cls):
        cr = 0.05
        cl = 0.15
        bo = 0.09
        pp = 0.00
        cv = 20e-6
        TestJopyRecip.exp.set_geometry(cr,cl,bo,pp,cv)
        TestJopyRecip.imp.set_geometry(cr,cl,bo,pp,cv)
        pass

    def test_recip_objects(self):
        assert not isinstance(TestJopyRecip.exp, RecipImplicit)
        assert isinstance(TestJopyRecip.imp, RecipBase)
        
    def test_recip_functions(self):
        rev = TestJopyRecip.exp.revolution(100)
        diff = np.abs((TestJopyRecip.exp.volume(rev)-TestJopyRecip.imp.volume(rev))/TestJopyRecip.exp.volume(rev))
        #print(np.max(diff),np.mean(diff))
        assert np.mean(diff)<0.005 # less than 0.5%

    @classmethod
    def teardown_class(cls):
        pass

class TestJopyStyle(object):
    objs = []

    @classmethod
    def setup_class(cls):
        TestJopyStyle.objs.append(BaseStyle())
        TestJopyStyle.objs.append(DtuStyle())
        TestJopyStyle.objs.append(IpuStyle())
        pass

    def test_style_show(self):
        for obj in TestJopyStyle.objs:
            line_fig,map_fig = obj._show_info(show=False)
            assert isinstance(line_fig, Figure)
            assert isinstance(map_fig, Figure)

    @classmethod
    def teardown_class(cls):
        pass
