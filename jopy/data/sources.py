
import os
from blaze.interactive import Data

def get_sqlite_handle(path):
    """Gets a blaze object for an sqlite database

    Does not simplify things, but help me remeber to use blaze more
    
    Parameters
    ----------
    path : str
        The actual path to the sqlite file
    
    Returns
    -------
    blaze.interactive.Data
        The database handle or None
    
    """
    if os.path.isfile(os.path.abspath(path)):
        tt = os.path.abspath(path)
        #tt = pathname2url(tt)
        sql = Data('sqlite:///'+tt) # A SQLite database
        return sql
    else: 
        print("Could not find "+path+".")
    return None 