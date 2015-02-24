"""
Example of running external tests with python tests.

20111116 rsz Created.
"""

# Set up path.
import sys
from os.path import normpath, dirname, join
local_module = normpath(join(dirname(__file__), '..'))
sys.path.insert(0, local_module)

from ezTest import *
from ezTestDemo import ezTestDemo

#######################################################################
if __name__=='__main__':
  suites = [
    ezTestDemo(),
    ezConsoleSuite('./ezTestDemo', 4),
  ]
  
  ezTestRunner(suites, opt=ezTestOpt()).run()
