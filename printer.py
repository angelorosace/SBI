from datetime import datetime

_is_verbose = False

# Flags the program to allow debugging print statements
def set_is_verbose():
  global _is_verbose
  _is_verbose = True

# When the program is verbose, prints the current timestamp
# and everything to be printed to standard output
def prnt(*args):
  if _is_verbose:
    print('[%s] -' % datetime.now(), *args)