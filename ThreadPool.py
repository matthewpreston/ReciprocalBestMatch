# ThreadPool.py - A simple implementation of a threading pool. Gratuitously
# plagiarized this: https://www.metachris.com/2016/04/python-threadpool/
# but with a few modifications for my sanity
#
# Written By: Matt Preston (and Chris Hager)
# Written On: May 31, 2017
# Revised On: Never

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import six

from threading import Thread
import logging

Queue = six.moves.queue

class Worker(Thread):
    """Thread executing tasks from a given tasks queue"""
    
    def __init__(self, tasks, name=None):
        Thread.__init__(self, name=name)
        self.tasks = tasks
        self.daemon = True
        self.start()

    def run(self):
        while True:
            try:
                func, args, kargs = self.tasks.get()
                func(*args, **kargs)
            except Exception as e:
                logging.error(e, exc_info=True)
            finally:
                self.tasks.task_done()

class ThreadPool(object):
    """Pool of threads consuming tasks from a queue"""
    
    def __init__(self, numThreads=0, names=None):
        self.tasks = Queue.Queue()
        self.numThreads = 0
        self.logger = logging.getLogger("%s.%s" % (__name__,
                                                   self.__class__.__name__))
        self.AddThreads(numThreads, names)    

    def AddThreads(self, numThreads, names=None):
        """Adds threads to the pool"""
        if isinstance(names, list): # Names given
            assert numThreads == len(names), \
                "Unequal pairing of names to number of threads"
            for name in names:
                Worker(self.tasks, name)
        elif names is None:         # No names
            for _ in range(numThreads):
                Worker(self.tasks)
        else:                       # Invalid type for 'names'
            self.logger.warning("ThreadPool.__init__(): 'names' is not a list "
                                "or not None, ignoring")
            for _ in range(numThreads):
                Worker(self.tasks)
        self.numThreads += numThreads
            
    def AddTask(self, func, *args, **kargs):
        """Add a task to the queue"""
        self.tasks.put((func, args, kargs))

    def HaveMainThreadWork(self):
        """
        Let main thread work through the tasks as well; automatically waits for
        threads to finish as well
        """
        while not self.tasks.empty():
            try:
                func, args, kargs = self.tasks.get()
                func(*args, **kargs)
                self.tasks.task_done()
            except Queue.Empty as e:
                pass
            except Exception as e:
                self.logger.error(e, exc_info=True)
                self.tasks.task_done()
        self.Join()
        
    def Join(self):
        """Wait for completion of all the tasks in the queue"""
        self.tasks.join()
