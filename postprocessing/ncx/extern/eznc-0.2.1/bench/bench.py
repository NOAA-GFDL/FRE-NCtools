"""
Run different combinations of netcdf file settings.

rm -f /tmp/*.nc
python bench.py > /tmp/tmp.log
tail -f /tmp/tmp.log
grep -v ncbench /tmp/tmp.log > /tmp/tmp.csv

20111019 rsz Created.
"""
import sys, os, string
import subprocess

#######################################################################
def run(c):
  p =subprocess.Popen(c.split(), stdout=subprocess.PIPE)
  if p.stdout:
    return p.stdout.readlines()
  else:
    return ""

#######################################################################
def main():
  nt = 100
  nz = 50
  ny = 100
  nx = 200
  hpad = 0
  salign = 4
  spad = 0
  ralign = 4
  TMPFILE = '/tmp/tmp.nc'
  
  args = ' nt nx ny nz format initsz blksz hpad salign spad ralign prefill chsize chelems chpreempt'
  
  # Results. Make args comma delimited, remove first comma, add stats labels.
  print args.replace(' ',',')[1:] + ',Total(s),Meta(s),Static(s),Rec(s),StaticGB,RecGB,MBps,VirtMemMax,RssMemMax'
  
  for format in '32 64 hdf'.split():
    for initsz in [0, 1.5e9]:
      for blksz in [0, 2**14, 2**16, 2**18, 2**20]:
        for prefill in '0 1'.split():
          for chsize in [0, 2**13, 2**15, 2**16]:
            for chelems in [0, 100, 1000]:
              for chpreempt in [0, .5, 1]:
                vars = {}
                for v in args.split():
                  vars[v] = eval(v)
                  
                f = run('mktemp -u /tmp/tmp.XXX.nc')[0].strip()
                pat = string.Template('./ncbench ' + args.replace(' ', ' $') + ' ' + f)

                cmd = pat.substitute(vars)
                print cmd
                stdout = run(cmd)
                os.remove(f)
                print cmd[10:-15].replace(' ',',')[1:] + stdout[1].strip()
                sys.stdout.flush()

#######################################################################
if __name__=='__main__':
  main()