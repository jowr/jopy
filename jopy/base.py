# -*- coding: utf-8 -*-
from __future__ import print_function, division

class JopyBaseClass(object):
    """The base class for all objects
    
    The mother of all classes in the jopy module. Implements 
    basic functionality for debugging and exception handling.
    """     
    def __init__(self):
        """Summary line.

        Extended description of function, just a usage example 
        for the NumPy style docstrings. See also: 
        http://sphinx-doc.org/ext/example_numpy.html#example-numpy
        
        The mother of all classes in the jopy module. Implements 
        basic functionality for debugging and exception handling.
        
        Parameters
        ----------
        arg1 : int
            Description of arg1
        arg2 : str
            Description of arg2
        
        Returns
        -------
        bool
            Description of return value
        
        """
        self.DEBUG = False
         
    @property
    def DEBUG(self):
        return self._DEBUG
    @DEBUG.setter
    def DEBUG(self, value):
        self._DEBUG = value
    @DEBUG.deleter
    def DEBUG(self):
        del self._DEBUG 
        
    def autolog(self, message):
        """Centralised logging facility

        Use this function in your code to write to the log files. It can
        also be extended to perform some more sophisticated actions 
        for advanced error detection.
        
        Function name and line number get prepended automatically.
    
        Parameters
        ----------
        message : str
            message to log
    
        """
        import inspect, logging
        # Get the previous frame in the stack, otherwise it would
        # be this function!!!
        func = inspect.currentframe().f_back.f_code
        msg = "%s: %s in %s:%i" % (
            message, 
            func.co_name, 
            func.co_filename, 
            func.co_firstlineno
        )
        # Dump the message + the name of this function to the log.
        logging.debug(msg)
        if self.DEBUG:
            print(msg)