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
        TestJopyRecip.exp.set_geometry(cr,cl,bo,pp,cv)
        TestJopyRecip.imp.set_geometry(cr,cl,bo,pp,cv)
        pass

    def test_recip_objects(self):
        assert not isinstance(TestJopyRecip.exp, RecipImplicit)
        assert isinstance(TestJopyRecip.imp, RecipBase)
        
    def test_recip_functions(self):
        rev = TestJopyRecip.exp.revolution(100)
        assert np.max(np.abs(TestJopyRecip.exp.volume(rev)-TestJopyRecip.imp.volume(rev)))<1e-10

    @classmethod
    def teardown_class(cls):
        pass
