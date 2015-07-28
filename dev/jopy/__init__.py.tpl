from __future__ import print_function

__author__ = '{{ author }}'
__email__ = '{{ email }}'
__version__ = '{{ version }}'

# Fix matplotlib backend on X-less machines
try: 
	import matplotlib
	USE_MATPLOTLIB=True
except Exception as e:
	print(str(e))
	print("Disabling matplotlib functions.")
	USE_MATPLOTLIB=False
	
if USE_MATPLOTLIB:
	import os
	if "DISPLAY" not in os.environ: 
		matplotlib.use('Agg')
		USE_DISPLAY=False
	else:
		USE_DISPLAY=True
else:
	USE_DISPLAY=False

# Add more imports here is you like 