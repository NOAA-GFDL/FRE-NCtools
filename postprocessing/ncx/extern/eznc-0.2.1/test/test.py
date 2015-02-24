"""
20111114 rsz Created.
20120107 rsz Added test_ezNcDateTimeSlicer, test_ezNcFiles.
"""
from ezTest import *

###########################################################
if __name__=='__main__':
  suites = [
    ezConsoleSuite('./test_ezNcBuffers', 13),
    ezConsoleSuite('./test_ezNcCF', 4),
    ezConsoleSuite('./test_ezNcCompare', 2),
    ezConsoleSuite('./test_ezNcCopySlabTask', 1),
    ezConsoleSuite('./test_ezNc', 49),
    ezConsoleSuite('./test_ezNcDateTime', 8),
    ezConsoleSuite('./test_ezNcDateTimeSlicer', 1),
    ezConsoleSuite('./test_ezNcFiles', 2),
    ezConsoleSuite('./test_ezNcSlab', 3),
    ezConsoleSuite('./test_ezNcSlabIterator', 4),
    ezConsoleSuite('./test_ezNcSlicer', 1),
    ezConsoleSuite('./test_ezNcUtil', 15),
    ezConsoleSuite('./test_ezNcVariant', 3),
    ezConsoleSuite('./test_ezPointerVector', 1),
  ]
  
  ezTestRunner(suites, opt=ezTestOpt()).run()