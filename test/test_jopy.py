"""
Tests for `jopy` module.
"""
from jopy.recip import RecipExplicit, RecipImplicit
import numpy as np


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
        bo = 0.05
        pp = 0.005
        cv = 10e-6
        exp.set_geometry(cr,cl,bo,pp,cv)
        imp.set_geometry(cr,cl,bo,pp,cv)
        pass

    def test_recip_objects(self):
        assert not isinstance(exp, RecipImplicit)
        assert isinstance(imp, RecipBase)
        
    def test_recip_functions(self):
        rev = exp.revolution(100)
        assert np.max(np.abs(exp.volume(rev)-imp.volume(rev)))<1e-10

    @classmethod
    def teardown_class(cls):
        pass
