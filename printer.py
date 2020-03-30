from datetime import datetime

_is_verbose = False

def set_is_verbose():
  global _is_verbose
  _is_verbose = True

def prnt(*args):
  if _is_verbose:
    print('[%s] -' % datetime.now(), *args)