"""Python implementation of selected C++ parallelUtil classes.
"""

# $Source: /usr/data0/leipzig_work/tmat_cvs/src/parallel_util.py,v $

# /* *********************************************************************** */
# /*                                                                         */
# /* Copyright (C) 2010  Kort Travis                                         */
# /*                                                                         */
# /*                                                                         */
# /* This program is free software; you can redistribute it and/or modify    */
# /* it under the terms of the GNU Lesser General Public License as          */
# /* published by the Free Software Foundation; version 2.1 of the License.  */
# /*                                                                         */
# /* This program is distributed in the hope that it will be useful,         */
# /* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
# /* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
# /* GNU Lesser General Public License for more details.                     */
# /*                                                                         */
# /* You should have received a copy of the GNU Lesser General Public        */
# /* License along with this program; if not, write to the Free Software     */
# /* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
# /* USA.                                                                    */
# /*                                                                         */
# /* *********************************************************************** */


import sys, os, threading, copy

# define parallel python methods _only_ if the module is installed:
global pp_loaded
try:
  import pp
  pp_loaded = True
except ImportError:
  pp_loaded = False


__all__ = ['worker_thread', 'multi_mutex', 'threaded_iterator', 'threaded_block']
 
# ------------------ GLOBAL variables used by this module ----------------------------
global N_threads_, NUM_THREADS_

# -------- init globals: -------------
N_threads_ = 0
NUM_THREADS_ = int(os.getenv('OMP_NUM_THREADS') or 1)
  
# ------------------------------------------------------------------------------------

class worker_thread(threading.Thread):
  #  Do _not_ use "__slots__" here: to allow attachment of arbitrary attributes.

  def __init__(self, group=None, target=None, name=None, args=(), kwargs={}):
    global N_threads_, NUM_THREADS_    
    threading.Thread.__init__(self, group=group, target=target, name=name, args=args, kwargs=kwargs)
    
    self.__thread_offset = N_threads_
    N_threads_ = (N_threads_ + 1) % NUM_THREADS_
    
    self.__return_value = None
    
       
  @staticmethod
  def __new__(cls, *args, **kwargs):
    return super(worker_thread, cls).__new__(cls, *args, **kwargs)

  @staticmethod
  def NUM_THREADS():
    """
    number of parallel _worker_ threads as suggested by the environment
    (often this will be OMP_NUM_THREADS)
    """
    global NUM_THREADS_
    return NUM_THREADS_

  def thread_offset(self):
    return self.__thread_offset

  def run(self):
    ' run thread-function '
    if self._Thread__target:
      self.__return_value = self._Thread__target(*self._Thread__args, **self._Thread__kwargs)
    
  def collect(self):
    ' join the thread, returning the thread-function\'s return value '    
    self.join()
    return self.__return_value
    
  
  def return_value(self):
    ' the return value from the thread-function (to allow non-threaded usage) '
    return self.__return_value  


class multi_mutex:
  ' set of id-associated mutexes '
  __slots__ = ['locks_']

  @staticmethod 
  def __new__(cls, *args, **kwargs):
    return super(multi_mutex, cls).__new__(cls, *args, **kwargs)
  
  def __init__(self):
    self.locks_ = {}
     
  def __del__(self):
    for id in self.locks_:
      self.locks_.pop(id)
      
  def acquire(self, id, blocking=True):
    ' acquire lock associated with id for current thread (or block) '
    test = None
    if id in self.locks_:
      test = self.locks_[id].acquire(blocking)
    else:
      l = threading.Lock()
      self.locks_[id] = l  
      test = l.acquire(blocking)
    
    return test
    
  def unlock(self, id):
    ' release lock associated with id '
    self.locks_.pop(id)
  
  def debug_print(self):
    print 'multi_mutex: \n  %s' % self.locks_


class threaded_iterator:
  ' round-robin partitioning of any iterable '
  
  # keep local copies of thread bookkeeping info:
  __slots__ = ['__TLS', '__iterable', 'NUM_THREADS']

  @staticmethod 
  def __new__(cls, *args, **kwargs):
    return super(threaded_iterator, cls).__new__(cls, *args, **kwargs)
  
  def __init__(self, iterable, **kwargs):
    """
    threaded_iterator(iterable, **kwargs)
      kwargs: attributes _supercede_ any currentThread().thread_offset or other default attribute values 
    """
    self.__TLS = threading.local()
    self.__iterable = iterable

    # transfer any attributes from kwargs to TLS:
    # (note: these attributes supercede those initialized below)
    for attr, val in kwargs.iteritems():
      self.__TLS.__dict__[attr] = val

    if (not hasattr(self, 'NUM_THREADS')):
      self.NUM_THREADS = worker_thread.NUM_THREADS()
        
    if (not hasattr(self.__TLS, 'thread_offset')):
      if hasattr(threading.currentThread(), 'thread_offset'):
        self.__TLS.thread_offset = threading.currentThread().thread_offset() # shall be the derived "worker_thread" class
      elif self.NUM_THREADS > 1:
        self.__TLS.thread_offset = self.NUM_THREADS               # main thread gets _empty_ partition    
      else:
        self.__TLS.thread_offset = 0                              # main thread gets _entire_ partition    
  
  def __del__(self):
    pass  
  
  def __thread_offset(self):
    return self.__TLS.thread_offset
    
  def __own_index(self, n):
    ' index position in partition of iterable associated with current thread '
    return (n % self.NUM_THREADS) == self.__thread_offset()

  def __iter__(self):
    for n, item  in enumerate(self.__iterable):
      if not self.__own_index(n):
        continue
      else:
        yield item

  def debug_print(self):
    print 'parallel_util::threaded_iterator:'
    print '  thread_offset: %d' % self.__TLS.thread_offset
    print '  __iterable: %s' % self.__iterable
    print '  NUM_THREADS: %d' % self.NUM_THREADS

          
def threaded_block(target=None, args=(), kwargs={}):
  """
  <result list> = threaded_block(target=None, args=(), kwargs={}) 
    execute target method in worker threads, combine results into <result list>
  """
  thread_list = [worker_thread(target=target, args=args, kwargs=kwargs) for n in range(worker_thread.NUM_THREADS())]
  result_list = []
  
  for th in thread_list:
    th.start()
  for th in thread_list:
    result_list.append(th.collect())

  return result_list
  
if pp_loaded:
  def pp_threaded_block(pp_server, target, args=(), kwargs={}, depfuncs=(), modules=(), callback=None, callbackargs=(), group='default', globals=None, NUM_THREADS=None):
    """
    <result list> = pp_threaded_block(pp_server, target, args=(), kwargs={}, depfuncs=(), modules=(), callback=None, callbackargs=(), group='default', globals=None) 
      execute target method in worker threads, combine results into <result list>
      (note: target method must accept kwargs "thread_offset", and "NUM_THREADS")
    """
    if (not NUM_THREADS):
      NUM_THREADS = worker_thread.NUM_THREADS()
    
    """
    # allow thread-function to have keyword arguments:
    def thread_foo(args_, kwargs_):
      return target(*args_, **kwargs_)
    """
      
    result_holding_list = []
    result_list = []
       
    kwargs_dup = copy.deepcopy(kwargs) # make a copy in order to add "thread_offset" attribute
    kwargs_dup['NUM_THREADS'] = NUM_THREADS
    for thread_offset in range(NUM_THREADS):
      kwargs_dup['thread_offset'] = thread_offset
      args_ = list(args)
      args_.append(kwargs_dup) # workaround: no kwargs in parallel python
      result_holding_list.append(pp_server.submit(target, tuple(args_), (), modules, callback, callbackargs, group, globals=globals))
    for held_result in result_holding_list:
      result_list.append(held_result())

    return result_list



  
