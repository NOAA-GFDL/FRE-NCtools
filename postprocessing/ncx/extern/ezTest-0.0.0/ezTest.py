"""
Template unit test class and runner.
To use, do 'from ezTest import ezTestRunner' and use ezTest or ezConsoleTest.
For demo, run 'python ezTest.py --verbose --color'.

20111021 rsz Created.
20111116 rsz Added color.

"""
# NB: time.clock() is best on windows, lousy on linux, so fall back on time.time().
import os, sys, time, traceback, shlex, subprocess, tempfile, hashlib
import numpy as np

#######################################################################
def ezDeleteFiles(names):
  """ Delete file(s). """
  if type(names) <> type([]):
    names = [names]
    
  for n in names:
    if os.path.exists(n):
      os.unlink(n)
        
#######################################################################
def ezHash(data):
  """ Return sha512 digest hash for an input array or string. """
  h = hashlib.sha512()
  if hasattr(data, 'shape'):
    # Loop over outer dimension in case data is very large or read from disk.
    for i in xrange(data.shape[0]):
      h.update(data[i,...])
  else:
    h.update(data)
    
  return h.hexdigest()

#######################################################################
def ezShell(c, flat=False):
  """
  Execute shell command and return output lines as list or flat string.
  Inlines stderr with stdout to preserve print order and ease parsing.
  """
  p = subprocess.Popen(shlex.split(c), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  if p.stdout:
    if flat:
      return p.stdout.read()
    else:
      return p.stdout.readlines()
  else:
    return ""

#######################################################################
def ezTmpFilenames(prefix='', suffix='.tmp', n=1):
  """ Return single or many filenames. """
  def New(): 
    f = tempfile.NamedTemporaryFile(prefix=prefix, delete=True)
    f.close()
    return f.name + suffix
  
  if n == 1:
    return New()
  else:
    fs = []
    while n > 0:
      fs.append(New())
      n = n - 1
      
    return fs

#######################################################################
class ezCommand(object):
  """
  Run a command and capture it's output string, error string and exit status.
  http://www.daniweb.com/software-development/python/code/257449
  
  Modified to keep stdout and stderr together so they are shown correctly in order upon testing.
  """

  def __init__(self, command):
    self.command = command
 
  def run(self, shell=True):
    import subprocess as sp
    # Keep stdout and stderr together inlined for messages to be useful.
    process = sp.Popen(self.command, shell = shell, stdout = sp.PIPE, stderr = sp.STDOUT)
    self.pid = process.pid
    self.output, self.error = process.communicate()
    self.failed = process.returncode
    return self
   
  @property
  def returncode(self):
    return self.failed
  
#######################################################################
class ezMsg:
  """ Simple error message wrapper for test functions. """
  def __init__(self):
    self.msg = ''

#######################################################################
class ezTestSuite(object):
  def __init__(self):
    self.name = ''
    self.opt = {}
    self.tests = []
    
  def echo(self, msg):
    if self.opt.get('verbose'): print msg
    
  def GetNumTests(self): 
    return len(self.tests)

  def rm(self, filenames):
    if self.opt.get('clean'): ezDeleteFiles(filenames)

  def setUp(self): 
    pass
    
  def tearDown(self): 
    pass

  def tmpfiles(self, suffix='.tmp', n=1):
    """ Return single or many filenames. """
    prefix = self.opt['tmpdir']
    if prefix:
      # Ensure trailing slash.
      #if sys.platform in ['win32']:
      #  prefix += '\\'
      #else:
      prefix += '/'
    
    def New(): 
      f = tempfile.NamedTemporaryFile(prefix=prefix, delete=True)
      f.close()
      return f.name + suffix
    
    if n == 1:
      return New()
    else:
      fs = []
      while n > 0:
        fs.append(New())
        n = n - 1
        
    return fs

#######################################################################
class ezConsoleSuite(object):
  """
  Wraps command-line (console) invoked test program suite to allow the Python test runner to combine heterogenous test suites and report them nicely together.
  
  Console programs must conform to these ezTest criteria:
  - exit code is 0 for success, otherwise the number of failed unit tests.
  - must accept last arguments as test case ids to invoke. This would run just two cases, numbers 2 and 4: "./test_myExample 2 4"
  
  The test runner will report overall runtime for this suite, the console stdout and stderr, and will tally the complete total of passed tests with other suites.
  """
  def __init__(self, cmd, n):
    """
    cmd: Commandline to execute test.
    n: Total number of possible tests for suite.
    """
    self.cmd = cmd
    self.n = n
    
  def GetNumTests(self): return self.n
  
#######################################################################
class ezTestRunner(object):
  """ 
  Supply list of ezTest suites and they will be run automatically.
  Select specific tests of suites or ranges to run with 'tests' option.
  """
  def __init__(self, suites, opt={}):
    self.suites = suites
    self.opt = opt
    self.tests = self.CreateTestLists()
    
  def CreateTestLists(self):
    """
    Create list of test id lists where each list is for a test suite.
    Parses self.suites, which may look like '1,2,3/1,2' or 'none/1:3/1,2/1,4:6'.
    The test ids returned are 1-based, so remember to offset them  by -1 if using ids directly to index arrays.
    """
    nsuites = len(self.suites)
    list = [[]]*nsuites
    ctr = 0
    tests = self.opt.get('tests')
    if not tests: 
      tests = '/'.join(':'*nsuites) # All tests.
      
    # Parse it. Create list of lists, where test ids are 1-based.
    for s in tests.split('/'):
      sublist = []
      if s == ':':
        sublist = range(1, self.suites[ctr].GetNumTests()+1)
      elif s <> 'none':
        split = s.split(',')
        for t in split:
          if t.count(':'):
            u,v = t.split(':')
            sublist.extend(range(int(float(u)), int(float(v))+1))
          else:
            sublist.append(int(float(t)))
      
      list[ctr] = sublist
      ctr += 1

    return list
    
  def InitColors(self):
    if self.opt['color']:
      self.R = '\x1b[91m'
      self.G = '\x1b[92m'
      self.Y = '\x1b[93m'
      self.B = '\x1b[94m'
      self.E = '\x1b[0m'
    else:
      self.R = ''
      self.G = ''
      self.Y = ''
      self.B = ''
      self.E = ''
      
  def RunConsoleSuite(self, suite, list):
    """
    suite: An ezConsoleSuite.
    list: A list of test ids to run (1-based).
    """
    n = len(list)
    print '%s%s (Running %d of %d Possible Tests)%s' % (self.B, suite.cmd, n, suite.GetNumTests(), self.E)
    ids = str(list).replace(',','')[1:-1]
    clean = '--clean %d' % (self.opt['clean'])
    tmpdir = ''
    if self.opt['tmpdir']:
      tmpdir = '--tmpdir %s' % (self.opt['tmpdir'])
    
    color = ''
    if self.opt['color']:
      color = '--color'

    verbose = ''
    if self.opt['verbose']:
      verbose = '--verbose'
      
    cmd = '%s %s %s %s %s %s' % (suite.cmd, clean, tmpdir, color, verbose, ids)
    
    #ticSuite = time.clock()
    com = ezCommand(cmd).run()
    #tocSuite = time.clock()
    
    # If fatal exit, return code may be large like 127, which is not num failed tests.
    if com.returncode > n: 
      nPass = 0
    else: 
      nPass = n - com.returncode
      
    print com.output
    if com.error: print com.error
    sys.stdout.flush()
    
    return nPass

  def RunPySuite(self, suite, list):
    """
    suite: A Python suite like ezTestDemo.
    list: A list of test ids to run (1-based).
    """
    
    n = len(list)
    print '%s%s (Running %d of %d Possible Tests)%s' % (self.B, suite.name, n, len(suite.tests), self.E)
    nPass = 0
    ticSuite = time.time()
    suite.setUp()
    for tid in list:
      fun,desc = suite.tests[tid-1]
      print '%s%d. %s %s' % (self.B, tid, desc, self.E)    
      tic = time.time()
      passed = False
      
      try:
        passed = fun()
      except Exception,err:
        print self.R,
        traceback.print_exc(file=sys.stdout)
        print self.E,

      toc = time.time()          
      if not passed: 
        print self.R + 'FAIL',
      else: 
        if self.opt['verbose']: print self.G + 'PASS',
        nPass += 1
        
      if self.opt['verbose'] or not passed:
        print '(%.5f sec)%s' % ((toc-tic), self.E)
        
      if self.opt['verbose']:
        print '%s-----------------------------------------------------------------%s' % (self.B, self.E)

    suite.tearDown()
    tocSuite = time.time()
    if nPass == n:
      finalColor = self.G
    else:
      finalColor = self.Y
      
    print '%s%d of %d tests passed in %.5f seconds.%s' % (finalColor, nPass, n, tocSuite-ticSuite, self.E)
    sys.stdout.flush()
    
    return nPass
    
  def run(self):
    self.InitColors()
    
    # Start time.
    ticTotal = time.time()
    nSuites = 0

    # Count total number of tests for all suites.
    nAllTests = 0

    # Total number of passed tests.
    nAllPass = 0
    suiteId = -1
    for list in self.tests:
      suiteId += 1
      n = len(list)
      if n < 1:
        continue
        
      s = self.suites[suiteId]
      nSuites += 1

      if isinstance(s, ezConsoleSuite):
        nPass = self.RunConsoleSuite(s, list)
      else:
        nPass = self.RunPySuite(s, list)
        print ''
        
      nAllTests += n
      nAllPass += nPass
      
    tocTotal = time.time()
    if nSuites > 1:
      if nAllPass == nAllTests:
        finalColor = self.G
      else:
        finalColor = self.Y
        
      print '%s%d of %d total tests passed in %.5f seconds.%s' % (finalColor, nAllPass, nAllTests, tocTotal-ticTotal, self.E)
    
#######################################################################
def ezTestOptUsage():
  u = 'SYNOPSIS\n'
  u += '\t%s -- Run unit test framework.\n\n' % sys.argv[0]
  u += 'USAGE\n'
  u += '\t%s [options]\n\n' % sys.argv[0]
  u += 'OPTIONS\n'
  u += '\t-c, --clean FLAG   Enable or disable temporary file removal.\n'
  u += '\t-C, --color        Enable colored messages.\n'
  u += '\t-h, --help         Show this message.\n'
  u += '\t-t, --tests ARG    Select which test cases instead of all.\n'
  u += '\t                   Suites are delimited by /.\n'
  u += '\t                   Example for 3 suites for some/all/no tests is\n'
  u += '\t                   1,2,3:10,15/:/none\n'
  u += '\t-T, --tmpdir T     Specify temporary directory for suites.\n'
  u += '\t-v, --verbose      Print verbose messages.\n'

  return u
 
#######################################################################
def ezTestOpt():
  """
  Returns dictionary of parsed command-line options.
  """
  import getopt
  shortopt = "c:Cht:T:v"
  longopt = ['clean=', 'color', "help", 'tests=', 'tmpdir=', 'verbose']
  try:
    opts, args = getopt.getopt(sys.argv[1:], shortopt,longopt)
  except getopt.GetoptError:
    print ezTestOptUsage()
    sys.exit(2)

  options = {}
  for k in longopt: options[k.replace('=','')] = None
  options['clean'] = 1
  options['tmpdir'] = ''
  
  for op, arg in opts:
    if op in ("-h", "--help"):
      print ezTestOptUsage()
      sys.exit(1)
    elif op in ("-c", "--clean"):
      options['clean'] = int(arg)
    elif op in ("-C", "--color"):
      options['color'] = 1
    elif op in ("-t", "--tests"):
      options['tests'] = arg
    elif op in ("-T", "--tmpdir"):
      options['tmpdir'] = arg
    elif op in ("-v", "--verbose"):
      options['verbose'] = 1

  return options
  
