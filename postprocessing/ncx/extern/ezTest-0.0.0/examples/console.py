"""
Example of running external test (script, executable) within python wrapper.
ezTest options like "--color" are propagated to the console test. See ezConsoleSuite for interface criteria.

20111116 rsz Created.
"""

# Set up path.
import sys
from os.path import normpath, dirname, join
local_module = normpath(join(dirname(__file__), '..'))
sys.path.insert(0, local_module)

from ezTest import *

#######################################################################
if __name__=='__main__':
  suites = [
    ezConsoleSuite('./ezTestDemo', 4), # The command and number of possible tests.
  ]
  
  ezTestRunner(suites, opt=ezTestOpt()).run()
