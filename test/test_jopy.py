"""
Tests for `jopy` module.
"""
from jopy.recip import RecipMachine, RecipMachineAlt


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

    @classmethod
    def setup_class(cls):
        pass

    def test_something(self):
        obj = RecipMachine()
        alt = RecipMachineAlt()
        assert not isinstance(obj, RecipMachineAlt)
        assert isinstance(alt, RecipMachine)

    @classmethod
    def teardown_class(cls):
        pass
