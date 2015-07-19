# -*- coding: utf-8 -*-
from __future__ import print_function, division

import jopy.thermo

class TestUtils(object):

    @classmethod
    def setup_class(cls):
        pass

    def test_lmtd(self):
        res = jopy.thermo.lmtd(0.0, 0.0)
        assert res == 0.0
        res = jopy.thermo.lmtd(10.0, 10.0)
        print(res,10.0)

    @classmethod
    def teardown_class(cls):
        pass