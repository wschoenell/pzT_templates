__author__ = 'william'

from pickle import Pickler
import os

class AtomicPickler(Pickler):
  def __init__(self, protocol):
    # You may want to replace this with a fake file object that just
    # discards writes.
    blackhole = open(os.devnull, 'w')

    Pickler.__init__(self, blackhole, protocol)
    self.depth = 0

  def save(self, o):
    self.depth += 1
    if self.depth == 1:
      return Pickler.save(self, o)
    self.depth -= 1
    return

def is_atomically_pickleable(o, protocol=None):
  pickler = AtomicPickler(protocol)
  try:
    pickler.dump(o)
    return True
  except:
    # Hopefully this exception was actually caused by dump(), and not
    # something like a KeyboardInterrupt
    return False