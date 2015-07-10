'''
Created on 3 Apr 2014

@author: jowr
'''
class JopyBaseClass(object):     
    def __init__(self):
        """The mother of all classes in the jopy 
        module. Implements basic functionality 
        for debugging and exception handling."""
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
        "Automatically log the current function details."
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
            print msg