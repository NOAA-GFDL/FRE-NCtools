"""
Simple demo of ezTest python interface.

python ezTestDemo.py --color --verbose

2011116 rsz Created.
"""

# Set up path.
import sys
from os.path import normpath, dirname, join
local_module = normpath(join(dirname(__file__), '..'))
sys.path.insert(0, local_module)

from ezTest import *

#######################################################################
class ezTestDemo(ezTestSuite):
  """ A test suite with single setup/teardown and many tests. """
  
  def __init__(self):
    self.name = 'ezTest Example Suite'
    
    # Ordered list of functions to run and descriptive name.
    # The runner will assign a test number.
    # Each test function should return True or False.
    self.tests = [
      [self.testExample1, 'First example test.'],
      [self.testExample2, 'Second example test that fails.'],
    ]
    
  def testExample1(self): return True
  def testExample2(self): return False

#######################################################################
if __name__=='__main__':
  # Create list of suites.
  suites = [ezTestDemo(), ezTestDemo()]
  # Get command-line options.
  cliOptions = ezTestOpt()
  
  runner = ezTestRunner(suites, opt=cliOptions)
  runner.run()