__author__ = '{{ author }}'
__email__ = '{{ email }}'
__version__ = '{{ version }}'

# Fix matplotlib backend on X-less machines 
import os, matplotlib
if "DISPLAY" not in os.environ: matplotlib.use('Agg')

# Add more import here is you like 