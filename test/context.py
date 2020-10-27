"""Define the path so tests can import the L1A reading routines.

Sample useage::
    from .context import oci_l1a

File taken from https://docs.python-guide.org/writing/structure/#setup-py
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# import oci_l1a
